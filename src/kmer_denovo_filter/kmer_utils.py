"""K-mer utility functions."""

import logging
import os
import subprocess
import tempfile
import threading
import time

import ahocorasick

logger = logging.getLogger(__name__)

_COMP = str.maketrans("ACGTacgt", "TGCAtgca")


def _is_symbolic(allele):
    """Return True if *allele* is a symbolic VCF allele with no literal sequence.

    Symbolic alleles include ``<DEL>``, ``<INS>``, ``<DUP>``, breakend
    notation containing ``[`` or ``]``, and the overlapping-deletion
    marker ``*``.
    """
    if not allele:
        return True
    return allele[0] == "<" or allele == "*" or "[" in allele or "]" in allele


def reverse_complement(seq):
    """Return the reverse complement of a DNA sequence."""
    return seq.translate(_COMP)[::-1]


def canonicalize(kmer):
    """Return the canonical (lexicographically smaller) form of a k-mer."""
    rc = kmer.translate(_COMP)[::-1]
    return kmer if kmer < rc else rc


def build_kmer_automaton(canonical_kmers):
    """Build an Aho-Corasick automaton from canonical k-mers.

    Adds both the forward and reverse-complement of each canonical
    k-mer so that reads can be scanned without per-position
    canonicalization.  The stored value for every pattern is the
    original canonical k-mer.

    Args:
        canonical_kmers: Iterable of canonical k-mer strings.

    Returns:
        An :class:`ahocorasick.Automaton` ready for ``iter()``,
        or ``None`` if *canonical_kmers* is empty.
    """
    A = ahocorasick.Automaton()
    for kmer in canonical_kmers:
        A.add_word(kmer, kmer)
        rc = reverse_complement(kmer)
        if rc != kmer:
            A.add_word(rc, kmer)
    if len(A) == 0:
        return None
    A.make_automaton()
    return A


def estimate_automaton_memory_gb(n_kmers):
    """Estimate the in-memory size of an Aho-Corasick automaton in GB.

    Empirically measured at ~928 bytes per pattern for 31-mer DNA
    sequences (including trie node, failure links, and Python string
    value overhead).  With forward + reverse-complement, the pattern
    count is roughly ``2 × n_kmers``.

    Args:
        n_kmers: Number of canonical k-mers.

    Returns:
        Estimated memory in GB (float).
    """
    n_patterns = n_kmers * 2  # fwd + RC
    bytes_per_pattern = 928  # empirically measured
    total_bytes = n_patterns * bytes_per_pattern
    return total_bytes / (1024**3)


# ── Jellyfish-backed k-mer query (disk-backed, low memory) ─────────


def _extract_read_kmers(seq, kmer_size):
    """Extract canonicalized k-mers from a read sequence.

    Args:
        seq: Read sequence string.
        kmer_size: Length of k-mers.

    Returns:
        Tuple of (canon_at_pos, unique_candidates):
            canon_at_pos: dict mapping query position to canonical k-mer.
            unique_candidates: deduplicated list of canonical k-mers
                preserving first-seen order.
    """
    seq_len = len(seq)
    if seq_len < kmer_size:
        return {}, []

    seq_upper = seq.upper()
    canon_at_pos = {}
    candidates = []

    for i in range(seq_len - kmer_size + 1):
        kmer = seq_upper[i:i + kmer_size]
        if "N" in kmer:
            continue
        canon = canonicalize(kmer)
        canon_at_pos[i] = canon
        candidates.append(canon)

    unique_candidates = list(dict.fromkeys(candidates))
    return canon_at_pos, unique_candidates


class JellyfishKmerQuery:
    """Query k-mers against a Jellyfish index file without loading into memory.

    Uses ``jellyfish query -s /dev/stdin`` subprocesses so that the
    k-mer hash table lives in jellyfish's memory-mapped address space
    (shared via OS page cache across workers) rather than in Python.

    For WGS-scale discovery with hundreds of millions of proband-unique
    k-mers, this uses ~2–10 GB (the jellyfish hash file, memory-mapped)
    instead of ~20–460 GB (Python set or Aho-Corasick automaton).

    A result cache avoids re-querying k-mers that have already been
    seen.  Adjacent reads in a sorted BAM share most of their k-mers,
    so the cache dramatically reduces the number of subprocess calls.

    Usage::

        q = JellyfishKmerQuery("/path/to/proband_unique.jf")
        for read in bam:
            hits = q.scan_read(read.query_sequence, kmer_size=31)
            # hits: set of (canonical_kmer, query_start_index)
        q.close()
    """

    def __init__(self, jf_path):
        self.jf_path = jf_path
        self._cache = {}  # canonical_kmer -> bool (present in index)

    def _subprocess_query(self, kmers):
        """Query k-mers via a jellyfish subprocess.

        Uses ``communicate()`` to safely handle stdin/stdout without
        risk of deadlock, even for large batches.

        Args:
            kmers: Iterable of canonical k-mer strings.

        Returns:
            Set of canonical k-mers that are present (count > 0).
        """
        fasta_block = "".join(
            f">{i}\n{kmer}\n" for i, kmer in enumerate(kmers)
        )
        proc = subprocess.Popen(
            ["jellyfish", "query", self.jf_path, "-s", "/dev/stdin"],
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        stdout, _ = proc.communicate(fasta_block.encode())

        hits = set()
        for line in stdout.decode().split("\n"):
            line = line.strip()
            if not line:
                continue
            parts = line.split()
            if len(parts) >= 2 and parts[1] != "0":
                hits.add(parts[0])
        return hits

    def query_batch(self, canonical_kmers):
        """Query a batch of canonical k-mers against the index.

        Results are cached so that repeated queries for the same k-mer
        (common with overlapping reads) skip the subprocess entirely.

        Args:
            canonical_kmers: List of canonical k-mer strings.

        Returns:
            Set of canonical k-mers that are present (count > 0) in the
            index.
        """
        if not canonical_kmers:
            return set()

        uncached = [k for k in canonical_kmers if k not in self._cache]
        if uncached:
            hits = self._subprocess_query(uncached)
            for k in uncached:
                self._cache[k] = k in hits

        return {k for k in canonical_kmers if self._cache.get(k, False)}

    def scan_read(self, seq, kmer_size):
        """Scan a read sequence for proband-unique k-mers.

        Extracts all k-mers via a sliding window, canonicalizes each,
        and batch-queries the jellyfish index.

        Args:
            seq: Read sequence string (uppercase recommended).
            kmer_size: Length of k-mers.

        Returns:
            Tuple of (unique_in_read, kmer_hit_indices):
                unique_in_read: set of canonical k-mers found.
                kmer_hit_indices: set of query start indices.
        """
        canon_at_pos, unique_candidates = _extract_read_kmers(seq, kmer_size)

        if not unique_candidates:
            return set(), set()

        hits = self.query_batch(unique_candidates)

        unique_in_read = set()
        kmer_hit_indices = set()
        for pos, canon in canon_at_pos.items():
            if canon in hits:
                unique_in_read.add(canon)
                kmer_hit_indices.add(pos)

        return unique_in_read, kmer_hit_indices

    def close(self):
        """Clear the result cache."""
        self._cache.clear()

    def __del__(self):
        self.close()


# ── Kraken2-backed bacterial content classification ────────────────


# Standard NCBI taxonomy ID for the Bacteria domain
_BACTERIA_TAXID = 2
_HUMAN_TAXID = 9606

# Interval between Kraken2 memory heartbeat log messages (seconds)
_KRAKEN2_HEARTBEAT_INTERVAL = 30
# Timeout when joining the heartbeat thread after Kraken2 completes (seconds)
_KRAKEN2_HEARTBEAT_JOIN_TIMEOUT = 2


def _read_proc_rss_kb(pid):
    """Read RSS memory in kB for *pid* from ``/proc/{pid}/status``.

    Returns ``None`` when the file is unavailable (non-Linux, or the
    process has already exited).
    """
    try:
        with open(f"/proc/{pid}/status") as fh:
            for line in fh:
                if line.startswith("VmRSS:"):
                    return int(line.split()[1])
    except OSError:
        pass
    return None


class Kraken2Runner:
    """Classify reads with kraken2 and tally bacterial content.

    Wraps the ``kraken2`` binary in a subprocess-based interface
    analogous to :class:`JellyfishKmerQuery`.  Reads are written to a
    temporary FASTQ, classified, and the kraken2 per-read output is
    parsed to count how many reads are classified as bacterial
    (``taxid 2`` or any descendant in the LCA hierarchy).

    The ``--confidence`` threshold (default 0.0) controls how strict
    the LCA classification must be.  A value of 0.2 requires at least
    20 % of k-mers in the read to agree on the assigned clade.

    Usage::

        kr = Kraken2Runner("/path/to/kraken2_db")
        result = kr.classify_sequences({"read1": "ACGT...", "read2": ...})
        print(result.bacterial_reads)  # set of read names
        print(result.summary())        # human-readable counts
    """

    class Result:
        """Container for a kraken2 classification run.

        Attributes:
            total: Total reads classified.
            classified: Number of reads that received a taxonomic assignment.
            unclassified: Number of reads with no assignment.
            bacterial_read_names: Set of read names assigned to Bacteria
                (taxid 2) or descendant taxa.
            bacterial_count: Number of bacterial reads.
            human_count: Number of reads assigned to Homo sapiens (taxid 9606)
                or descendants.
            root_count: Number of reads assigned to root (taxid 1) with no
                more specific classification.
        """

        def __init__(self):
            self.total = 0
            self.classified = 0
            self.unclassified = 0
            self.bacterial_read_names = set()
            self.bacterial_count = 0
            self.human_count = 0
            self.root_count = 0

        def summary(self):
            """Return a human-readable summary string."""
            pct = (
                f"{100 * self.bacterial_count / self.total:.1f}"
                if self.total > 0
                else "0.0"
            )
            return (
                f"kraken2: {self.total} reads, "
                f"{self.classified} classified, "
                f"{self.bacterial_count} bacterial ({pct}%), "
                f"{self.human_count} human, "
                f"{self.root_count} root"
            )

        @property
        def bacterial_fraction(self):
            """Fraction of reads classified as bacterial (0.0–1.0)."""
            if self.total == 0:
                return 0.0
            return round(self.bacterial_count / self.total, 4)

    def __init__(self, db_path, *, confidence=0.0, threads=1, memory_mapping=False):
        self.db_path = db_path
        self.confidence = confidence
        self.threads = threads
        self.memory_mapping = memory_mapping

    # ── taxonomy helpers ───────────────────────────────────────────

    @staticmethod
    def _load_bacterial_taxids(db_path):
        """Load the set of taxonomy IDs that descend from Bacteria.

        Parses ``taxonomy/nodes.dmp`` within the kraken2 database
        directory.  Any taxid whose lineage passes through taxid 2
        (Bacteria) is included.

        If the taxonomy directory is not present (e.g. a pre-built
        database without nodes.dmp), returns ``None`` so that the
        caller can emit a warning and fall back to direct-taxid-only
        matching.
        """
        nodes_path = os.path.join(db_path, "taxonomy", "nodes.dmp")
        if not os.path.isfile(nodes_path):
            return None

        parent_map = {}
        try:
            with open(nodes_path) as fh:
                for line in fh:
                    parts = line.split("\t|\t")
                    if len(parts) < 3:
                        continue
                    child_id = int(parts[0].strip())
                    parent_id = int(parts[1].strip())
                    parent_map[child_id] = parent_id
        except (OSError, ValueError):
            return None

        # Walk up from every node; cache membership
        bacterial = set()
        non_bacterial = set()

        def _is_bacterial(taxid):
            if taxid in bacterial:
                return True
            if taxid in non_bacterial:
                return False
            path = []
            cur = taxid
            while cur not in bacterial and cur not in non_bacterial:
                if cur == _BACTERIA_TAXID:
                    bacterial.update(path)
                    bacterial.add(cur)
                    return True
                if cur == 1 or cur == 0 or cur not in parent_map:
                    non_bacterial.update(path)
                    non_bacterial.add(cur)
                    return False
                path.append(cur)
                cur = parent_map[cur]
            if cur in bacterial:
                bacterial.update(path)
                return True
            non_bacterial.update(path)
            return False

        for taxid in parent_map:
            _is_bacterial(taxid)

        return bacterial

    @staticmethod
    def _extract_taxids_from_kmer_string(kmer_string):
        """Extract integer taxonomy IDs from kraken2 k-mer output field."""
        if not kmer_string:
            return set()

        taxids = set()
        # Paired-end output can include the '|:|' delimiter between mates.
        for token in kmer_string.replace("|:|", " ").split():
            taxid, _, _ = token.partition(":")
            if not taxid:
                continue
            try:
                taxids.add(int(taxid))
            except ValueError:
                continue
        return taxids

    @staticmethod
    def _estimate_mmap_db_size_bytes(db_path):
        """Estimate Kraken2 DB bytes potentially mapped by --memory-mapping.

        Returns:
            Tuple ``(total_bytes, file_sizes)``, where ``file_sizes`` is a
            mapping of basename → file size in bytes for ``*.k2d`` files.
        """
        file_sizes = {}
        try:
            for name in os.listdir(db_path):
                if not name.endswith(".k2d"):
                    continue
                path = os.path.join(db_path, name)
                if os.path.isfile(path):
                    file_sizes[name] = os.path.getsize(path)
        except OSError:
            return None, {}
        return sum(file_sizes.values()), file_sizes

    # ── classification ─────────────────────────────────────────────

    def classify_sequences(self, sequences, tmpdir=None):
        """Classify named sequences and return a :class:`Result`.

        Args:
            sequences: Dict mapping read name → sequence string,
                **or** a list of ``(name, sequence)`` tuples.
            tmpdir: Optional directory for temporary FASTQ file.

        Returns:
            A :class:`Kraken2Runner.Result` with tallied counts.
        """
        if isinstance(sequences, dict):
            items = sequences.items()
        else:
            items = sequences

        result = self.Result()

        # Materialize into a list so we can count without consuming
        items = list(items)
        if not items:
            return result

        result.total = len(items)

        # Write a temporary FASTQ
        fd, fastq_path = tempfile.mkstemp(
            suffix=".fq", prefix="kraken2_", dir=tmpdir,
        )
        try:
            with os.fdopen(fd, "w") as fh:
                for name, seq in items:
                    qual = "I" * len(seq)  # dummy Phred 40
                    fh.write(f"@{name}\n{seq}\n+\n{qual}\n")

            # Build kraken2 command
            cmd = [
                "kraken2",
                "--db", self.db_path,
                "--threads", str(self.threads),
                "--confidence", str(self.confidence),
                "--output", "-",          # per-read output to stdout
                "--report", "/dev/null",  # suppress summary report
            ]
            if self.memory_mapping:
                cmd.append("--memory-mapping")
            cmd.append(fastq_path)

            proc = subprocess.Popen(
                cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
            )

            # Background thread: log RSS memory heartbeats while kraken2 runs
            kraken2_start = time.monotonic()
            stop_heartbeat = threading.Event()
            peak_rss_kb = [None]
            previous_rss_kb = [None]

            if self.memory_mapping:
                total_db_bytes, db_file_sizes = self._estimate_mmap_db_size_bytes(
                    self.db_path
                )
                logger.info(
                    "[Kraken2] memory-mapping enabled; Kraken2 avoids preloading "
                    "the full DB into process-local RAM and faults DB pages "
                    "from disk on demand",
                )
                if total_db_bytes is not None:
                    files_desc = ", ".join(
                        f"{name}={size / 1_073_741_824:.2f} GB"
                        for name, size in sorted(db_file_sizes.items())
                    )
                    logger.info(
                        "[Kraken2] memory-mapping DB footprint: %.2f GB across "
                        "%d .k2d files (%s)",
                        total_db_bytes / 1_073_741_824,
                        len(db_file_sizes),
                        files_desc,
                    )
                    logger.info(
                        "[Kraken2] memory-mapping expectation: RSS may grow as "
                        "pages are touched; practical upper bound is near the "
                        "touched DB footprint (up to ~%.2f GB if all DB pages "
                        "are touched)",
                        total_db_bytes / 1_073_741_824,
                    )

            def _heartbeat():
                while not stop_heartbeat.wait(_KRAKEN2_HEARTBEAT_INTERVAL):
                    rss = _read_proc_rss_kb(proc.pid)
                    elapsed = time.monotonic() - kraken2_start
                    if rss is not None:
                        if peak_rss_kb[0] is None or rss > peak_rss_kb[0]:
                            peak_rss_kb[0] = rss
                        logger.info(
                            "[Kraken2] heartbeat — %.0f s elapsed, "
                            "RSS: %.1f GB (peak: %.1f GB, delta: %+0.1f GB)",
                            elapsed, rss / 1_048_576,
                            peak_rss_kb[0] / 1_048_576,
                            (
                                0.0 if previous_rss_kb[0] is None
                                else (rss - previous_rss_kb[0]) / 1_048_576
                            ),
                        )
                        previous_rss_kb[0] = rss
                    else:
                        logger.info(
                            "[Kraken2] heartbeat — %.0f s elapsed "
                            "(memory info unavailable)",
                            elapsed,
                        )

            heartbeat_thread = threading.Thread(
                target=_heartbeat, daemon=True, name="kraken2-heartbeat",
            )
            heartbeat_thread.start()
            try:
                stdout, stderr = proc.communicate()
            finally:
                stop_heartbeat.set()
                heartbeat_thread.join(timeout=_KRAKEN2_HEARTBEAT_JOIN_TIMEOUT)

            elapsed = time.monotonic() - kraken2_start
            if proc.returncode != 0:
                logger.warning(
                    "kraken2 exited with code %d after %.0f s: %s",
                    proc.returncode, elapsed,
                    stderr.decode(errors="replace").strip()[:500],
                )
                return result

            if peak_rss_kb[0] is not None:
                logger.info(
                    "[Kraken2] classification complete — %d reads in %.0f s "
                    "(peak RSS: %.1f GB)",
                    result.total, elapsed, peak_rss_kb[0] / 1_048_576,
                )
            else:
                logger.info(
                    "[Kraken2] classification complete — %d reads in %.0f s "
                    "(peak RSS: unavailable)",
                    result.total, elapsed,
                )

            # Load bacterial taxid set for lineage-aware matching
            bacterial_taxids = self._load_bacterial_taxids(self.db_path)
            if bacterial_taxids is None:
                logger.warning(
                    "Kraken2 taxonomy lineage matching is unavailable "
                    "(missing/unreadable taxonomy/nodes.dmp under DB: %s). "
                    "Falling back to exact taxid==2 matching only; "
                    "bacterial fractions may be severely undercounted.",
                    self.db_path,
                )

            # Parse per-read output
            # Format: C/U\tread_name\ttaxid\tlength\tkmers_string
            for line in stdout.decode(errors="replace").split("\n"):
                line = line.strip()
                if not line:
                    continue
                parts = line.split("\t")
                if len(parts) < 3:
                    continue

                status = parts[0]
                read_name = parts[1]
                try:
                    taxid = int(parts[2])
                except ValueError:
                    continue
                kmer_taxids = self._extract_taxids_from_kmer_string(
                    parts[4] if len(parts) >= 5 else "",
                )

                if status == "U":
                    result.unclassified += 1
                    continue

                result.classified += 1

                # Check bacterial assignment
                is_bacterial = False
                if bacterial_taxids is not None:
                    is_bacterial = taxid in bacterial_taxids
                else:
                    # Fallback: only exact taxid 2
                    is_bacterial = taxid == _BACTERIA_TAXID
                # Be conservative: when Kraken2 exposes explicit human
                # k-mer evidence for a read, avoid flagging it as
                # bacterial even if the assigned LCA is in Bacteria.
                if is_bacterial and _HUMAN_TAXID in kmer_taxids:
                    is_bacterial = False

                if is_bacterial:
                    result.bacterial_count += 1
                    result.bacterial_read_names.add(read_name)
                elif taxid == _HUMAN_TAXID:
                    result.human_count += 1
                elif taxid == 1:
                    result.root_count += 1

        finally:
            try:
                os.unlink(fastq_path)
            except OSError:
                pass

        return result


def read_supports_alt(
    read, variant_pos, ref, alt, min_baseq=0, *,
    aligned_pairs=None, seq=None, quals=None,
):
    """Return True if *read* carries the alternate allele at *variant_pos*.

    Extracts the exact read sequence aligned to the reference span of the
    variant and compares it strictly to the candidate alternate allele.
    Handles SNPs, MNPs, insertions, deletions, and complex indels natively.

    Returns ``False`` for symbolic alleles or when *alt* is ``None``.

    Args:
        min_baseq: Minimum base quality threshold for bases considered as
            alt support.
        aligned_pairs: Optional pre-computed result of
            ``read.get_aligned_pairs(matches_only=False)``.  Computed from
            *read* when not provided.
        seq: Optional pre-decoded ``read.query_sequence``.  Decoded from
            *read* when not provided.
        quals: Optional pre-decoded ``read.query_qualities``. Decoded from
            *read* only when ``min_baseq > 0`` and not provided.
    """
    if alt is None or _is_symbolic(alt):
        return False

    if seq is None:
        seq = read.query_sequence
    if seq is None:
        return False
    if min_baseq > 0 and quals is None:
        quals = read.query_qualities

    if aligned_pairs is None:
        aligned_pairs = read.get_aligned_pairs(matches_only=False)

    extracted_seq = []
    in_variant_region = False

    for qpos, rpos in aligned_pairs:
        # Stop collecting once we reach or pass the end of the reference allele span
        if rpos is not None and rpos >= variant_pos + len(ref):
            break

        # Start collecting when we hit the exact start of the variant
        if rpos == variant_pos:
            in_variant_region = True

        if in_variant_region:
            # qpos is None for deleted bases (skip), otherwise append the read base
            if qpos is not None:
                if (
                    min_baseq > 0 and quals is not None
                    and quals[qpos] < min_baseq
                ):
                    return False
                extracted_seq.append(seq[qpos])

    # If the variant region was skipped entirely due to read boundaries
    if not in_variant_region:
        return False

    return "".join(extracted_seq).upper() == alt.upper()


def extract_variant_spanning_kmers(
    read, variant_pos, k, min_baseq=0, ref=None, alt=None,
    *, aligned_pairs=None, seq=None, quals=None,
):
    """Extract canonical k-mers from a read that span the variant position.

    Args:
        read: pysam AlignedSegment
        variant_pos: 0-based reference position of the variant
        k: k-mer size
        min_baseq: Minimum base quality threshold
        ref: Reference allele string (for INDEL handling)
        alt: Alternate allele string (for INDEL handling)
        aligned_pairs: Ignored; kept for API compatibility.
        seq: Optional pre-decoded ``read.query_sequence``.  Decoded from
            *read* when not provided.
        quals: Optional pre-decoded ``read.query_qualities``.  Decoded from
            *read* when not provided.

    Returns:
        Set of canonical k-mer strings spanning the variant position.
    """
    try:
        read_pos_at_variant = read.get_reference_positions(full_length=True).index(variant_pos)
    except ValueError:
        return set()

    if seq is None:
        seq = read.query_sequence
    if seq is None:
        return set()
    if quals is None:
        quals = read.query_qualities

    # For insertions the variant occupies len(alt) bases in the read.
    # Extend the window so k-mers spanning the right junction are captured.
    alt_len = len(alt) if alt and not _is_symbolic(alt) else 1
    variant_end_in_read = read_pos_at_variant + alt_len - 1

    kmers = set()
    start_min = max(0, read_pos_at_variant - k + 1)
    start_max = min(len(seq) - k, variant_end_in_read)

    # Pre-compute a boolean array marking bad positions (N or low quality)
    # so that the inner loop can use a sliding-window counter instead of
    # re-scanning every k-mer from scratch: O(window) instead of O(k × window).
    window_end = start_max + k  # exclusive upper bound of bases touched
    seq_upper = seq[start_min:window_end].upper()
    bad = bytearray(window_end - start_min)
    for i, ch in enumerate(seq_upper):
        if ch == 'N':
            bad[i] = 1
    if quals is not None and min_baseq > 0:
        for i in range(window_end - start_min):
            if quals[start_min + i] < min_baseq:
                bad[i] = 1

    # Initialise sliding-window bad-base count for the first k-mer
    bad_count = sum(bad[:min(k, len(bad))])

    for s in range(start_min, start_max + 1):
        offset = s - start_min
        if offset > 0:
            # Slide: remove leftmost base of previous window, add new rightmost
            bad_count -= bad[offset - 1]
            bad_count += bad[offset + k - 1]
        if bad_count:
            continue
        kmers.add(canonicalize(seq[s:s + k]))

    return kmers
