"""BAM/CRAM scanning engine and alignment helpers.

Provides the multi-backend k-mer hit scanner used by both VCF-mode and
discovery-mode pipelines.  Two backends are supported:

* **Aho-Corasick automaton** — fast C-level multi-pattern matching;
  used when the proband-unique k-mer set fits in Python memory.
* **JellyfishKmerQuery** — disk-backed queries via ``jellyfish query``;
  used for large k-mer sets where the Aho-Corasick automaton would
  exceed available RAM.

No function in this module accepts an ``argparse.Namespace`` object.
"""

import collections
import logging

import pysam

from kmer_denovo_filter.kmer_utils import (
    JellyfishKmerQuery,
    _extract_read_kmers,
    build_kmer_automaton,
)

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# FASTA I/O (local helper to avoid circular import with utils.py)
# ---------------------------------------------------------------------------


def _load_kmers_from_fasta(fasta_path):
    """Load k-mer strings from a FASTA file.

    Reads the file line-by-line, yielding only sequence lines (not
    header lines starting with ``>``).
    """
    kmers = set()
    with open(fasta_path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line and not line.startswith(">"):
                kmers.add(line)
    return kmers


# ---------------------------------------------------------------------------
# CIGAR / alignment helpers
# ---------------------------------------------------------------------------


def _extract_softclips(cigartuples):
    """Extract left and right soft-clip lengths from CIGAR tuples.

    Args:
        cigartuples: List of ``(operation, length)`` tuples from
            ``pysam.AlignedSegment.cigartuples``.  May be ``None``
            for unmapped reads.

    Returns:
        ``(softclip_left, softclip_right)`` tuple of integers.
    """
    if not cigartuples:
        return (0, 0)
    # CIGAR ops: 4 = CSOFT_CLIP, 5 = CHARD_CLIP
    # Hard clips may appear outside soft clips: e.g. 5H10S80M5S3H
    left = 0
    for op, length in cigartuples:
        if op == 4:  # soft clip
            left = length
            break
        elif op == 5:  # hard clip — skip and keep looking
            continue
        else:
            break

    right = 0
    for op, length in reversed(cigartuples):
        if op == 4:
            right = length
            break
        elif op == 5:
            continue
        else:
            break

    # Avoid double-counting when there is only one non-hard-clip CIGAR op
    non_hard = [t for t in cigartuples if t[0] != 5]
    if len(non_hard) == 1 and non_hard[0][0] == 4:
        right = 0

    return (left, right)


def _collect_kmer_ref_positions(read, kmer_hit_indices, kmer_size):
    """Map query-level k-mer hit positions to reference coordinates.

    Args:
        read: A pysam.AlignedSegment (must be mapped).
        kmer_hit_indices: Set of query start indices where novel k-mers
            were found.
        kmer_size: Length of k-mers.

    Returns:
        Counter keyed by reference position with coverage counts.
    """
    cov = collections.Counter()
    aligned_pairs = read.get_aligned_pairs(matches_only=True)
    query_to_ref = {qpos: rpos for qpos, rpos in aligned_pairs}
    for start_idx in kmer_hit_indices:
        for qpos in range(start_idx, start_idx + kmer_size):
            rpos = query_to_ref.get(qpos)
            if rpos is not None:
                cov[rpos] += 1
    return cov


def _infer_sv_type(region_a, region_b):
    """Infer SV type from two linked regions.

    Returns one of: INTRA (same chromosome) or BND (translocation).
    Distinguishing DEL/DUP/INV would require SA tag strand information
    which is not carried in the region data structure.
    """
    if region_a[0] != region_b[0]:
        return "BND"
    return "INTRA"


# ---------------------------------------------------------------------------
# Read alignment metadata collection
# ---------------------------------------------------------------------------


def _collect_read_alignment_metadata(
    child_bam, ref_fasta, read_names,
    informative_reads_by_variant=None,
):
    """Collect alignment metadata for informative reads from a BAM/CRAM.

    For each read in *read_names*, collects **all** alignment records
    (primary + supplementary) and stores per-alignment metadata needed
    for the genomic span BED file.

    Args:
        child_bam: Path to child BAM/CRAM.
        ref_fasta: Path to reference FASTA (may be ``None`` for BAM).
        read_names: Set of read names to collect metadata for.
        informative_reads_by_variant: Optional dict mapping variant keys
            to read-name sets for targeted BAM fetching.

    Returns:
        Dict mapping ``read_name`` → list of alignment record dicts,
        each with keys:

        - ``chrom`` (str): Reference contig name.
        - ``start`` (int): 0-based aligned start (``reference_start``).
        - ``end`` (int): 0-based exclusive aligned end (``reference_end``).
        - ``mapq`` (int): Mapping quality.
        - ``softclip_left`` (int): Soft-clipped bases at left end.
        - ``softclip_right`` (int): Soft-clipped bases at right end.
        - ``has_sa`` (bool): Whether the read has an SA tag.
        - ``is_supplementary`` (bool): Whether this record is supplementary.
    """
    if not read_names:
        return {}

    alignment_meta = {}  # read_name -> list of record dicts

    bam = pysam.AlignmentFile(
        child_bam, reference_filename=ref_fasta if ref_fasta else None,
    )

    def _process_read(read):
        if read.query_name not in read_names:
            return
        if read.is_unmapped:
            return
        sc_left, sc_right = _extract_softclips(read.cigartuples)
        has_sa = read.has_tag("SA")
        rec = {
            "chrom": read.reference_name,
            "start": read.reference_start,
            "end": read.reference_end,
            "mapq": read.mapping_quality,
            "softclip_left": sc_left,
            "softclip_right": sc_right,
            "has_sa": has_sa,
            "is_supplementary": read.is_supplementary,
        }
        alignment_meta.setdefault(read.query_name, []).append(rec)

    used_targeted_fetch = False
    if informative_reads_by_variant:
        loci_to_names = {}
        for var_key, names in informative_reads_by_variant.items():
            if not names:
                continue
            parts = var_key.split(":")
            if len(parts) < 2:
                continue
            chrom = parts[0]
            try:
                pos = int(parts[1])
            except ValueError:
                continue
            target_names = set(names).intersection(read_names)
            if not target_names:
                continue
            loci_to_names.setdefault((chrom, pos), set()).update(target_names)

        if loci_to_names:
            used_targeted_fetch = True
            seen = set()  # (read_name, is_supplementary, start) dedup key
            for (chrom, pos), target_names in sorted(loci_to_names.items()):
                for read in bam.fetch(chrom, pos, pos + 1):
                    key = (read.query_name, read.is_supplementary,
                           read.reference_start)
                    if key not in seen:
                        seen.add(key)
                        _process_read(read)

    if not used_targeted_fetch:
        for read in bam.fetch(until_eof=True):
            _process_read(read)

    bam.close()
    return alignment_meta


# ---------------------------------------------------------------------------
# Multi-process k-mer scanning engine
# ---------------------------------------------------------------------------

# Module-level worker state — set by _init_scan_worker via
# multiprocessing pool initialiser.
_worker_automaton = None       # Aho-Corasick automaton (small k-mer sets)
_worker_jf_query = None        # JellyfishKmerQuery (large k-mer sets)
_worker_kmer_size = None
_worker_min_distinct_kmers_per_read = 1

# Number of reads to accumulate before issuing a single jellyfish
# subprocess call.  Larger batches amortize subprocess overhead
# but temporarily hold more read objects in memory.
_JF_READ_BATCH_SIZE = 5000


def _init_scan_worker(proband_data, kmer_size,
                      min_distinct_kmers_per_read=1):
    """Initializer for per-contig scan workers.

    *proband_data* may be:
    - A path ending in ``.jf`` — opens a :class:`JellyfishKmerQuery`
      that queries k-mers against a memory-mapped jellyfish hash.
      This is the low-memory path used for WGS discovery mode with
      hundreds of millions of proband-unique k-mers.
    - A FASTA file path — loads k-mers into a Python set and builds
      an Aho-Corasick automaton (fast but memory-intensive).
    - A Python set — builds an Aho-Corasick automaton directly.
    """
    global _worker_automaton, _worker_jf_query, _worker_kmer_size
    global _worker_min_distinct_kmers_per_read

    _worker_automaton = None
    _worker_jf_query = None

    if isinstance(proband_data, str) and proband_data.endswith(".jf"):
        # Jellyfish-backed mode: each query subprocess memory-maps
        # the same .jf file; the OS page cache is shared across workers.
        _worker_jf_query = JellyfishKmerQuery(proband_data)
    elif isinstance(proband_data, str):
        kmers = _load_kmers_from_fasta(proband_data)
        _worker_automaton = build_kmer_automaton(kmers)
        del kmers
    else:
        _worker_automaton = build_kmer_automaton(proband_data)

    _worker_kmer_size = kmer_size
    _worker_min_distinct_kmers_per_read = min_distinct_kmers_per_read


def _process_informative_read(read, unique_in_read, kmer_hit_indices,
                              kmer_size, reads_seen, read_hits,
                              read_sv_meta, kmer_coverage, read_coverage):
    """Record an informative read's hits, coverage, and SV metadata.

    Returns 1 if the read is unmapped-informative, 0 otherwise.
    Mutates *reads_seen*, *read_hits*, *read_sv_meta*, *kmer_coverage*,
    and *read_coverage* in place.
    """
    dedup_key = (read.query_name, read.is_supplementary)
    if dedup_key in reads_seen:
        return 0

    reads_seen.add(dedup_key)
    if read.is_unmapped:
        return 1

    read_hits.append((
        read.reference_name,
        read.reference_start,
        read.reference_end,
        read.query_name,
        unique_in_read,
        read.is_supplementary,
    ))
    # Map novel k-mer query positions to reference coords
    chrom = read.reference_name
    cov = _collect_kmer_ref_positions(
        read, kmer_hit_indices, kmer_size,
    )
    kmer_coverage[chrom] += cov
    # Count one read per touched position
    for pos in cov:
        read_coverage[chrom][pos] += 1

    # Collect SV metadata for this informative read
    max_clip = 0
    if read.cigartuples:
        for op, length in read.cigartuples:
            if op == 4 and length > max_clip:  # soft clip
                max_clip = length
    read_sv_meta[dedup_key] = {
        "has_sa": read.has_tag("SA"),
        "sa_str": read.get_tag("SA") if (
            read.has_tag("SA") and not read.is_supplementary
        ) else None,
        "is_paired": read.is_paired,
        "is_proper_pair": read.is_proper_pair,
        "mate_is_unmapped": (
            read.mate_is_unmapped if read.is_paired else False
        ),
        "max_clip": max_clip,
    }
    return 0


def _scan_contig_for_hits(child_bam, ref_fasta, contig):
    """Scan reads mapped to *contig* for proband-unique k-mers.

    When *contig* is ``None``, unmapped reads are scanned instead.

    Supports two scanning backends:
    - **Aho-Corasick automaton** (``_worker_automaton``) — fast C-level
      multi-pattern matching; used when the k-mer set fits in memory.
    - **JellyfishKmerQuery** (``_worker_jf_query``) — disk-backed
      queries via ``jellyfish query``; used for large k-mer sets in
      discovery mode where the set is too large for Aho-Corasick.

    Returns:
        (read_hits, reads_seen, unmapped_informative, total_reads_scanned,
         read_sv_meta, kmer_coverage, read_coverage)

    ``read_sv_meta`` is a dict keyed by ``(query_name, is_supplementary)``
    with per-read SV metadata (has_sa, sa_str, is_paired, is_proper_pair,
    mate_is_unmapped, max_clip) collected for each informative read so
    that annotation and linking can be done without re-scanning the BAM.

    ``kmer_coverage`` is a dict mapping chrom to a Counter of reference
    positions overlapped by novel k-mers (counts total k-mer base
    overlaps across all reads).

    ``read_coverage`` is a dict mapping chrom to a Counter of reference
    positions where at least one novel k-mer was found, counting the
    number of distinct reads touching each position.
    """
    automaton = _worker_automaton
    jf_query = _worker_jf_query
    kmer_size = _worker_kmer_size
    min_dk_per_read = _worker_min_distinct_kmers_per_read
    bam = pysam.AlignmentFile(
        child_bam, reference_filename=ref_fasta if ref_fasta else None,
    )

    read_hits = []
    reads_seen = set()
    read_sv_meta = {}
    kmer_coverage = collections.defaultdict(collections.Counter)
    read_coverage = collections.defaultdict(collections.Counter)
    unmapped_informative = 0
    total_reads_scanned = 0

    if contig is None:
        try:
            iterator = bam.fetch("*")
        except (ValueError, KeyError):
            bam.close()
            return (read_hits, reads_seen, unmapped_informative,
                    total_reads_scanned, read_sv_meta, kmer_coverage,
                    read_coverage)
    else:
        iterator = bam.fetch(contig=contig)

    if jf_query is not None:
        # ── Batched jellyfish path ─────────────────────────────────
        # Reads are accumulated in batches so that k-mers from many
        # reads are queried in a single jellyfish subprocess call.
        # This reduces subprocess overhead from O(n_reads) to
        # O(n_reads / batch_size).
        pending = []   # (read, canon_at_pos)
        pending_kmers = set()

        for read in iterator:
            if read.is_secondary:
                continue
            if read.is_duplicate:
                continue

            total_reads_scanned += 1
            seq = read.query_sequence
            if seq is None:
                continue

            canon_at_pos, unique_candidates = _extract_read_kmers(
                seq, kmer_size,
            )
            pending_kmers.update(unique_candidates)
            pending.append((read, canon_at_pos))

            if len(pending) < _JF_READ_BATCH_SIZE:
                continue

            # Query all unique k-mers from this batch in one subprocess.
            # Avoid unbounded per-worker cache growth by processing each
            # batch against this local hit set, then clearing cache.
            batch_hits = set()
            if pending_kmers:
                batch_hits = jf_query.query_batch(list(pending_kmers))
                pending_kmers = set()

            # Process each read against batch-level hits
            for read_obj, c_at_pos in pending:
                unique_in_read = set()
                kmer_hit_indices = set()
                for pos, canon in c_at_pos.items():
                    if canon in batch_hits:
                        unique_in_read.add(canon)
                        kmer_hit_indices.add(pos)

                if len(unique_in_read) < min_dk_per_read:
                    continue

                unmapped_informative += _process_informative_read(
                    read_obj, unique_in_read, kmer_hit_indices,
                    kmer_size, reads_seen, read_hits,
                    read_sv_meta, kmer_coverage, read_coverage,
                )
            pending = []
            jf_query.close()

        # Flush remaining reads
        if pending:
            batch_hits = set()
            if pending_kmers:
                batch_hits = jf_query.query_batch(list(pending_kmers))
            for read_obj, c_at_pos in pending:
                unique_in_read = set()
                kmer_hit_indices = set()
                for pos, canon in c_at_pos.items():
                    if canon in batch_hits:
                        unique_in_read.add(canon)
                        kmer_hit_indices.add(pos)

                if len(unique_in_read) < min_dk_per_read:
                    continue

                unmapped_informative += _process_informative_read(
                    read_obj, unique_in_read, kmer_hit_indices,
                    kmer_size, reads_seen, read_hits,
                    read_sv_meta, kmer_coverage, read_coverage,
                )
            jf_query.close()
    else:
        # ── Aho-Corasick path (or no backend) ──────────────────────
        for read in iterator:
            if read.is_secondary:
                continue
            if read.is_duplicate:
                continue

            total_reads_scanned += 1
            seq = read.query_sequence
            if seq is None:
                continue

            unique_in_read = set()
            kmer_hit_indices = set()
            if automaton is not None:
                for _end_idx, canonical_kmer in automaton.iter(seq):
                    unique_in_read.add(canonical_kmer)
                    kmer_hit_indices.add(_end_idx - kmer_size + 1)

            if len(unique_in_read) < min_dk_per_read:
                continue

            unmapped_informative += _process_informative_read(
                read, unique_in_read, kmer_hit_indices,
                kmer_size, reads_seen, read_hits,
                read_sv_meta, kmer_coverage, read_coverage,
            )

    bam.close()
    return (read_hits, reads_seen, unmapped_informative,
            total_reads_scanned, read_sv_meta, kmer_coverage,
            read_coverage)
