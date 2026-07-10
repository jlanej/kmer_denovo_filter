"""Kraken2-backed non-human content classification engine.

This module wraps the external ``kraken2`` binary and reduces its per-read
LCA output to per-domain tallies and a non-human fraction.  It depends only
on the Python standard library and a ``kraken2`` executable on ``PATH``; it
never imports pysam, so it can classify any batch of named sequences
regardless of where they came from.
"""

import logging
import os
import struct
import subprocess
import tempfile
import threading
import time

from nonhuman_screen.result import TaxonomicFractions

logger = logging.getLogger(__name__)


# Standard NCBI taxonomy IDs for major clades
_BACTERIA_TAXID = 2
_ARCHAEA_TAXID = 2157
_FUNGI_TAXID = 4751
_EUKARYOTA_TAXID = 2759
_METAZOA_TAXID = 33208
_VIRIDIPLANTAE_TAXID = 33090
_VIRUSES_TAXID = 10239
_HUMAN_TAXID = 9606
# UniVec Core (NCBI taxid 81077): synthetic sequencing-vector and adapter
# sequences included in PrackenDB.  These are artificial constructs — not
# real biological organisms — and must be excluded from non-human counts so
# that reads misclassified here (due to vector k-mer contamination in library
# preparation or shared k-mers with human DNA) are not errantly flagged as
# evidence of microbial contamination.
_UNIVEC_CORE_TAXID = 81077

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
    """Classify reads with kraken2 and tally non-human content.

    Wraps the ``kraken2`` binary in a subprocess-based interface
    analogous to :class:`JellyfishKmerQuery`.  Reads are written to a
    temporary FASTQ, classified, and the kraken2 per-read output is
    parsed to count how many reads are classified as bacterial
    (``taxid 2``), archaeal (``taxid 2157``), fungal (``taxid 4751``),
    protist (eukaryotic but not metazoan, fungal, or plant), viral
    (``taxid 10239``), or non-human (any classified read definitively
    outside the human lineage).  All non-human tallies apply a
    conservative **human homology guard**: if a read's k-mer detail
    string includes any k-mer that voted for human (taxid 9606), the
    read is excluded from every non-human numerator.

    Viral reads receive the same human homology guard as all other
    non-human categories.  This is particularly important for viruses
    that can integrate into the human genome (e.g. endogenous
    retroviruses, HBV, HPV), whose integrated copies may share k-mers
    with the human reference.  A read carrying both viral and human
    k-mer evidence is conservatively excluded from the viral count.

    **UniVec Core handling**: PrackenDB includes UniVec Core (taxid
    81077), a curated set of synthetic sequencing-vector and adapter
    sequences.  These are artificial constructs and do not represent
    genuine biological contamination.  Reads classified under UniVec
    Core are tracked independently (``univec_core_read_names`` /
    ``univec_core_count``) and explicitly excluded from the consolidated
    non-human fraction (NHF) so that vector library-preparation
    artefacts — or human reads whose k-mers happen to match vector
    sequences — are not errantly counted as microbial contamination.
    The human homology guard is applied first; UniVec Core exclusion
    provides a second, unconditional safety net for any
    UniVec-classified read regardless of k-mer voting.

    The ``--confidence`` threshold (default 0.0) controls how strict
    the LCA classification must be.  A value of 0.2 requires at least
    20 % of k-mers in the read to agree on the assigned clade.

    Usage::

        kr = Kraken2Runner("/path/to/kraken2_db")
        result = kr.classify_sequences({"read1": "ACGT...", "read2": ...})
        print(result.bacterial_read_names)  # set of read names
        print(result.summary())             # human-readable counts
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
            archaeal_read_names: Set of read names assigned to Archaea
                (taxid 2157) or descendant taxa.
            archaeal_count: Number of archaeal reads.
            fungal_read_names: Set of read names assigned to Fungi
                (taxid 4751) or descendant taxa.
            fungal_count: Number of fungal reads.
            protist_read_names: Set of read names assigned to protist taxa
                (eukaryotic but not metazoan, fungal, or plant).
            protist_count: Number of protist reads.
            viral_read_names: Set of read names assigned to Viruses
                (taxid 10239) or descendant taxa.  Reads with any human
                k-mer evidence are excluded (human homology guard), which
                is particularly relevant for integrating viruses such as
                endogenous retroviruses, HBV, and HPV.
            viral_count: Number of viral reads (after human homology guard).
            univec_core_read_names: Set of read names assigned to UniVec
                Core (taxid 81077) or descendant taxa.  These are synthetic
                sequencing-vector and adapter sequences, not biological
                organisms.  The human homology guard is applied: reads with
                any human k-mer evidence are excluded.
            univec_core_count: Number of UniVec Core reads (after human
                homology guard).
            nonhuman_read_names: Set of read names definitively classified
                as non-human (any clade outside the human lineage).
                UniVec Core reads are excluded from this set.
            nonhuman_count: Number of non-human reads.
            unclassified_read_names: Set of read names that received no
                taxonomic assignment from Kraken2 (status ``"U"``).
            human_lineage_read_names: Set of read names that are classified
                but are **not** in ``nonhuman_read_names`` and **not** in
                ``univec_core_read_names``.  This covers reads directly
                assigned to the human clade, reads cleared by the human
                homology guard (HHG), and reads assigned to broad
                taxonomic ranks on the human-to-root lineage path
                (e.g. Root, Eukaryota).  Together with
                ``nonhuman_read_names``, ``univec_core_read_names``, and
                ``unclassified_read_names`` these four sets form a
                partition of all classified+unclassified reads, so the
                corresponding per-variant fractions sum to 1.
            human_lineage_count: Number of human-lineage reads (see
                ``human_lineage_read_names``).
            human_count: Number of reads assigned to Homo sapiens (taxid 9606)
                or descendants.
            root_count: Number of reads assigned to root (taxid 1) with no
                more specific classification.
            per_read_detail: Dict mapping read name to a dict with keys
                ``status`` (``"C"`` or ``"U"``), ``taxid`` (int),
                ``domain`` (str), ``guard_status`` (str),
                ``is_nonhuman`` (bool), and ``kmer_string`` (str).
                Populated during :meth:`classify_sequences` for every
                parsed read so the per-read audit trail can be written
                to the companion BED file.
        """

        def __init__(self):
            self.total = 0
            self.classified = 0
            self.unclassified = 0
            self.bacterial_read_names = set()
            self.bacterial_count = 0
            self.archaeal_read_names = set()
            self.archaeal_count = 0
            self.fungal_read_names = set()
            self.fungal_count = 0
            self.protist_read_names = set()
            self.protist_count = 0
            self.viral_read_names = set()
            self.viral_count = 0
            self.univec_core_read_names = set()
            self.univec_core_count = 0
            self.nonhuman_read_names = set()
            self.nonhuman_count = 0
            self.unclassified_read_names = set()
            self.human_lineage_read_names = set()
            self.human_lineage_count = 0
            self.human_count = 0
            self.root_count = 0
            self.per_read_detail = {}
            # False when the taxonomy dumps (nodes.dmp/names.dmp) could not be
            # loaded and classification fell back to exact-taxid matching, in
            # which case non-human fractions are undercounted.  Consumers that
            # gate on non-human content should treat a False here as "unknown"
            # rather than "clean".
            self.taxonomy_available = True

        def summary(self):
            """Return a human-readable summary string."""
            pct = (
                f"{100 * self.bacterial_count / self.total:.1f}"
                if self.total > 0
                else "0.0"
            )
            nh_pct = (
                f"{100 * self.nonhuman_count / self.total:.1f}"
                if self.total > 0
                else "0.0"
            )
            return (
                f"kraken2: {self.total} reads, "
                f"{self.classified} classified, "
                f"{self.bacterial_count} bacterial ({pct}%), "
                f"{self.archaeal_count} archaeal, "
                f"{self.fungal_count} fungal, "
                f"{self.protist_count} protist, "
                f"{self.viral_count} viral, "
                f"{self.univec_core_count} univec_core, "
                f"{self.nonhuman_count} non-human ({nh_pct}%), "
                f"{self.human_count} human, "
                f"{self.root_count} root"
            )

        @property
        def bacterial_fraction(self):
            """Fraction of reads classified as bacterial (0.0–1.0)."""
            if self.total == 0:
                return 0.0
            return round(self.bacterial_count / self.total, 4)

        @property
        def nonhuman_fraction(self):
            """Fraction of reads classified as non-human (0.0–1.0)."""
            if self.total == 0:
                return 0.0
            return round(self.nonhuman_count / self.total, 4)

        def fractions(self):
            """Return per-domain :class:`~nonhuman_screen.result.TaxonomicFractions`."""
            return TaxonomicFractions.from_result(self)

    def __init__(self, db_path, *, confidence=0.0, threads=1, memory_mapping=False):
        self.db_path = db_path
        self.confidence = confidence
        self.threads = threads
        self.memory_mapping = memory_mapping

    # ── database introspection ─────────────────────────────────────

    @staticmethod
    def read_kmer_length(db_path):
        """Return the k-mer length stored in a Kraken2 database.

        Kraken2 databases store build-time options (including the k-mer
        length ``k``) in ``opts.k2d`` as a binary ``IndexOptions`` struct.
        The first field of that struct is ``k`` (a ``size_t``, 8 bytes on
        64-bit platforms), so reading the first 8 bytes as a little-endian
        unsigned 64-bit integer gives the k-mer length used by the database.

        The PrackenDB pre-built database
        (``k2_NCBI_reference_20251007.tar.gz``) is built with Kraken2's
        default k-mer length of **35**.

        Args:
            db_path: Path to the Kraken2 database directory (must contain
                ``opts.k2d``, or a versioned subdirectory that does).

        Returns:
            The k-mer length as an integer, or ``None`` when ``opts.k2d``
            cannot be found or parsed.
        """
        # Look for opts.k2d directly in db_path, then one level deeper
        # (PrackenDB extracts into a versioned subdirectory).
        candidate_dirs = [db_path]
        try:
            for entry in os.scandir(db_path):
                if entry.is_dir():
                    candidate_dirs.append(entry.path)
        except OSError:
            pass

        for d in candidate_dirs:
            opts_path = os.path.join(d, "opts.k2d")
            if not os.path.isfile(opts_path):
                continue
            try:
                with open(opts_path, "rb") as fh:
                    data = fh.read(8)
                if len(data) == 8:
                    (k,) = struct.unpack("<Q", data)
                    if 1 <= k <= 256:  # sanity-check the value
                        return k
            except OSError:
                pass
        return None

    # ── taxonomy helpers ───────────────────────────────────────────

    @staticmethod
    def _load_parent_map(db_path):
        """Parse ``nodes.dmp`` from *db_path* and return a parent map.

        Tries ``taxonomy/nodes.dmp`` first, then ``nodes.dmp`` at the
        database root (PrackenDB layout).

        Returns:
            ``{child_taxid: parent_taxid}`` dict, or ``None`` if no
            readable ``nodes.dmp`` is found.
        """
        nodes_path = os.path.join(db_path, "taxonomy", "nodes.dmp")
        if not os.path.isfile(nodes_path):
            nodes_path = os.path.join(db_path, "nodes.dmp")
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
        return parent_map

    @staticmethod
    def _load_name_map(db_path):
        """Parse ``names.dmp`` and return a taxid → scientific name map.

        Tries ``taxonomy/names.dmp`` first, then ``names.dmp`` at the
        database root (PrackenDB layout).  Only rows whose name class
        (4th ``\\t|\\t``-delimited field) equals ``scientific name`` are
        retained.  Spaces in names are replaced with underscores.

        Args:
            db_path: Path to the Kraken2 database directory.

        Returns:
            ``{taxid: name_string}`` dict, or ``None`` if no readable
            ``names.dmp`` is found.
        """
        names_path = os.path.join(db_path, "taxonomy", "names.dmp")
        if not os.path.isfile(names_path):
            names_path = os.path.join(db_path, "names.dmp")
            if not os.path.isfile(names_path):
                logger.warning(
                    "names.dmp not found under %s; taxon names will be "
                    "unavailable in the per-read detail file.",
                    db_path,
                )
                return None

        name_map = {}
        try:
            with open(names_path) as fh:
                for line in fh:
                    parts = line.split("\t|\t")
                    if len(parts) < 4:
                        continue
                    # 4th field ends with "\t|\n"; strip to get name class
                    name_class = parts[3].replace("\t|", "").strip()
                    if name_class != "scientific name":
                        continue
                    try:
                        taxid = int(parts[0].strip())
                    except ValueError:
                        continue
                    name = parts[1].strip().replace(" ", "_")
                    name_map[taxid] = name
        except OSError:
            return None
        return name_map

    @staticmethod
    def _descendants_of(parent_map, root_taxid):
        """Return the set of taxids that descend from *root_taxid*.

        Walks the full taxonomy tree cached in *parent_map* (a
        ``{child: parent}`` dict) and returns all taxids whose lineage
        passes through *root_taxid* (inclusive).
        """
        members = set()
        non_members = set()

        def _is_member(taxid):
            if taxid in members:
                return True
            if taxid in non_members:
                return False
            path = []
            cur = taxid
            while cur not in members and cur not in non_members:
                if cur == root_taxid:
                    members.update(path)
                    members.add(cur)
                    return True
                if cur == 1 or cur == 0 or cur not in parent_map:
                    non_members.update(path)
                    non_members.add(cur)
                    return False
                path.append(cur)
                cur = parent_map[cur]
            if cur in members:
                members.update(path)
                return True
            non_members.update(path)
            return False

        for taxid in parent_map:
            _is_member(taxid)
        return members

    @staticmethod
    def _ancestors_of(parent_map, taxid):
        """Return the set of taxids on the lineage from *taxid* to root.

        Walks upward through *parent_map* and collects every node
        between *taxid* and root (inclusive).
        """
        ancestors = set()
        cur = taxid
        while cur in parent_map:
            ancestors.add(cur)
            parent = parent_map[cur]
            if parent == cur:  # root
                break
            cur = parent
        return ancestors

    @staticmethod
    def _load_bacterial_taxids(db_path):
        """Load the set of taxonomy IDs that descend from Bacteria.

        Parses ``taxonomy/nodes.dmp`` (or ``nodes.dmp`` at the database
        root for pre-built archives like PrackenDB) within the kraken2
        database directory.  Any taxid whose lineage passes through
        taxid 2 (Bacteria) is included.

        If neither location contains a readable ``nodes.dmp``, returns
        ``None`` so that the caller can emit a warning and fall back to
        direct-taxid-only matching.
        """
        parent_map = Kraken2Runner._load_parent_map(db_path)
        if parent_map is None:
            return None
        return Kraken2Runner._descendants_of(parent_map, _BACTERIA_TAXID)

    @staticmethod
    def _load_all_taxid_sets(db_path):
        """Load taxonomy and return descendant sets for all domains.

        Returns a dict with keys ``bacterial``, ``archaeal``, ``fungal``,
        ``protist``, ``viral``, ``univec_core``, ``human_lineage``, and
        ``human_clade``.  Each value is a set of NCBI taxonomy IDs.

        ``protist`` is defined as eukaryotic taxa that are **not**
        Metazoa, Fungi, or Viridiplantae.

        ``viral`` contains all descendants of Viruses (taxid 10239).
        Reads classified as viral are treated with particular care because
        some viruses (e.g. endogenous retroviruses, HBV, HPV) can integrate
        into the human genome. The human homology guard (checking for human
        k-mer evidence in the per-read detail string) conservatively
        excludes any read with both viral and human k-mer evidence from
        the viral numerator.

        ``univec_core`` contains all descendants of UniVec Core (taxid
        81077), a set of synthetic sequencing-vector and adapter sequences
        included in PrackenDB.  These artificial constructs are excluded
        from the consolidated non-human fraction unconditionally: they do
        not represent biological contamination and may share k-mers with
        human DNA, so classifying them as non-human would produce false
        positives.

        ``human_lineage`` contains every taxid on the path from human
        (9606) to root — these are taxonomic ranks too broad to be
        confidently assigned as non-human.

        ``human_clade`` contains human (9606) and any descendant
        subspecies or populations.

        Returns ``None`` when ``nodes.dmp`` is unavailable.
        """
        parent_map = Kraken2Runner._load_parent_map(db_path)
        if parent_map is None:
            return None

        bacterial = Kraken2Runner._descendants_of(parent_map, _BACTERIA_TAXID)
        archaeal = Kraken2Runner._descendants_of(parent_map, _ARCHAEA_TAXID)
        fungal = Kraken2Runner._descendants_of(parent_map, _FUNGI_TAXID)
        eukaryota = Kraken2Runner._descendants_of(parent_map, _EUKARYOTA_TAXID)
        metazoa = Kraken2Runner._descendants_of(parent_map, _METAZOA_TAXID)
        viridiplantae = Kraken2Runner._descendants_of(
            parent_map, _VIRIDIPLANTAE_TAXID,
        )
        protist = eukaryota - metazoa - fungal - viridiplantae
        viral = Kraken2Runner._descendants_of(parent_map, _VIRUSES_TAXID)
        univec_core = Kraken2Runner._descendants_of(
            parent_map, _UNIVEC_CORE_TAXID,
        )

        human_lineage = Kraken2Runner._ancestors_of(parent_map, _HUMAN_TAXID)
        human_clade = Kraken2Runner._descendants_of(parent_map, _HUMAN_TAXID)

        return {
            "bacterial": bacterial,
            "archaeal": archaeal,
            "fungal": fungal,
            "protist": protist,
            "viral": viral,
            "univec_core": univec_core,
            "human_lineage": human_lineage,
            "human_clade": human_clade,
        }

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

        # Log the database k-mer length so users can see which k-mer size
        # Kraken2 will use for LCA classification.
        kmer_len = self.read_kmer_length(self.db_path)
        if kmer_len is not None:
            logger.info("[Kraken2] database k-mer length: %d", kmer_len)
        else:
            logger.debug(
                "[Kraken2] could not read k-mer length from opts.k2d "
                "(db_path: %s)", self.db_path,
            )

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
                "--output", "/dev/stdout", # per-read output to stdout
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

            def _heartbeat():
                while not stop_heartbeat.wait(_KRAKEN2_HEARTBEAT_INTERVAL):
                    rss = _read_proc_rss_kb(proc.pid)
                    elapsed = time.monotonic() - kraken2_start
                    if rss is not None:
                        logger.info(
                            "[Kraken2] heartbeat — %.0f s elapsed, "
                            "RSS: %.1f GB",
                            elapsed, rss / 1_048_576,
                        )
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

            logger.info(
                "[Kraken2] classification complete — %d reads in %.0f s",
                result.total, elapsed,
            )

            # Load taxonomy sets for lineage-aware matching
            taxid_sets = self._load_all_taxid_sets(self.db_path)
            result.taxonomy_available = taxid_sets is not None
            if taxid_sets is None:
                logger.warning(
                    "Kraken2 taxonomy lineage matching is unavailable "
                    "(missing/unreadable taxonomy/nodes.dmp under DB: %s). "
                    "Falling back to exact taxid matching only; "
                    "non-human fractions may be severely undercounted.",
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
                    result.unclassified_read_names.add(read_name)
                    result.per_read_detail[read_name] = {
                        "status": "U",
                        "taxid": 0,
                        "domain": "Unclassified",
                        "guard_status": "UNCLASSIFIED",
                        "is_nonhuman": False,
                        "kmer_string": "",
                    }
                    continue

                result.classified += 1

                # Human homology guard: if any k-mer voted for human,
                # conservatively exclude this read from ALL non-human
                # numerators to avoid over-flagging human reads.
                has_human_kmer = _HUMAN_TAXID in kmer_taxids

                if taxid_sets is not None:
                    is_bacterial = taxid in taxid_sets["bacterial"]
                    is_archaeal = taxid in taxid_sets["archaeal"]
                    is_fungal = taxid in taxid_sets["fungal"]
                    is_protist = taxid in taxid_sets["protist"]
                    is_viral = taxid in taxid_sets["viral"]
                    is_univec_core = taxid in taxid_sets["univec_core"]
                    is_human = taxid in taxid_sets["human_clade"]
                    is_nonhuman = (
                        taxid not in taxid_sets["human_lineage"]
                        and taxid not in taxid_sets["human_clade"]
                        and taxid not in taxid_sets["univec_core"]
                    )
                else:
                    # Fallback: only exact taxid matching
                    is_bacterial = taxid == _BACTERIA_TAXID
                    is_archaeal = taxid == _ARCHAEA_TAXID
                    is_fungal = taxid == _FUNGI_TAXID
                    is_protist = False  # cannot determine without tree
                    is_viral = taxid == _VIRUSES_TAXID
                    is_univec_core = taxid == _UNIVEC_CORE_TAXID
                    is_human = taxid == _HUMAN_TAXID
                    is_nonhuman = taxid not in (_HUMAN_TAXID, 1, _UNIVEC_CORE_TAXID)

                # Determine domain label *before* the guard clears flags
                if is_bacterial:
                    _domain = "Bacteria"
                elif is_archaeal:
                    _domain = "Archaea"
                elif is_fungal:
                    _domain = "Fungi"
                elif is_protist:
                    _domain = "Protist"
                elif is_viral:
                    _domain = "Viruses"
                elif is_univec_core:
                    _domain = "UniVec_Core"
                elif is_human:
                    _domain = "Human"
                elif taxid == 1:
                    _domain = "Root"
                elif (taxid_sets is not None
                      and taxid in taxid_sets["human_lineage"]):
                    _domain = "Ambiguous_Ancestor"
                else:
                    _domain = "Root"

                # Apply human homology guard to all non-human categories.
                if has_human_kmer:
                    is_bacterial = False
                    is_archaeal = False
                    is_fungal = False
                    is_protist = False
                    is_viral = False
                    is_univec_core = False
                    is_nonhuman = False

                # Determine guard status string
                if is_human:
                    _guard = "HUMAN"
                elif has_human_kmer:
                    _guard = "HHG"
                elif _domain == "UniVec_Core":
                    _guard = "UVC"
                elif is_nonhuman:
                    _guard = "PASS"
                else:
                    _guard = "PASS"

                if is_bacterial:
                    result.bacterial_count += 1
                    result.bacterial_read_names.add(read_name)
                if is_archaeal:
                    result.archaeal_count += 1
                    result.archaeal_read_names.add(read_name)
                if is_fungal:
                    result.fungal_count += 1
                    result.fungal_read_names.add(read_name)
                if is_protist:
                    result.protist_count += 1
                    result.protist_read_names.add(read_name)
                if is_viral:
                    result.viral_count += 1
                    result.viral_read_names.add(read_name)
                if is_univec_core:
                    result.univec_core_count += 1
                    result.univec_core_read_names.add(read_name)
                if is_nonhuman:
                    result.nonhuman_count += 1
                    result.nonhuman_read_names.add(read_name)
                if not is_nonhuman and not is_univec_core:
                    # Classified but neither definitively non-human nor
                    # UniVec Core.  Covers: human clade, HHG-guarded,
                    # Root, and Ambiguous_Ancestor reads.
                    result.human_lineage_read_names.add(read_name)
                    result.human_lineage_count += 1
                if is_human:
                    result.human_count += 1
                elif taxid == 1:
                    result.root_count += 1

                result.per_read_detail[read_name] = {
                    "status": status,
                    "taxid": taxid,
                    "domain": _domain,
                    "guard_status": _guard,
                    "is_nonhuman": is_nonhuman,
                    "kmer_string": parts[4] if len(parts) >= 5 else "",
                }

        finally:
            try:
                os.unlink(fastq_path)
            except OSError:
                pass

        return result


# Public, import-friendly alias for the nested result container.
ClassificationResult = Kraken2Runner.Result
