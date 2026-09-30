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
import re

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


_CIGAR_OPS = {op: code for code, op in enumerate("MIDNSHP=X")}


def _parse_cigar(cigar):
    """Parse a CIGAR string into (op, length) tuples, like cigartuples."""
    return [
        (_CIGAR_OPS[op], int(length))
        for length, op in re.findall(r"(\d+)([MIDNSHP=X])", cigar)
    ]


def _query_span(cigartuples):
    """Return (start, end) of the aligned read bases, in reference order.

    Leading soft and hard clips come before *start*; matches, mismatches
    and insertions consume the read.  CIGARs are in reference
    orientation, so this indexes the read as stored for its strand.
    """
    start = 0
    for op, length in cigartuples or ():
        if op not in (4, 5):
            break
        start += length
    aligned = sum(length for op, length in cigartuples or ()
                  if op in (0, 1, 7, 8))
    return start, start + aligned


def _clip_side(cigartuples):
    """Return which end of an alignment faces a breakpoint.

    ``"+"`` when the alignment is clipped (soft or hard) more on its
    right, so the breakpoint follows it on the reference; ``"-"`` when
    more on its left; ``None`` when unclipped or equally clipped.  CIGARs
    are in reference orientation, so this holds on either strand.
    """
    if not cigartuples:
        return None
    left = right = 0
    for op, length in cigartuples:
        if op not in (4, 5):
            break
        left += length
    for op, length in reversed(cigartuples):
        if op not in (4, 5):
            break
        right += length
    if right > left:
        return "+"
    if left > right:
        return "-"
    return None


_SV_TYPE_BY_STRANDS = {
    ("+", "-"): "DEL", ("-", "+"): "DUP",
    ("+", "+"): "INV", ("-", "-"): "INV",
}


def _infer_sv_type(region_a, region_b, strands=None):
    """Infer the SV type joining two regions (or two breakends in one).

    *strands* is the breakpoint orientation at (region_a, region_b), with
    region_a the earlier one: ``"+"`` when the sequence joined at the
    breakpoint lies left of it, ``"-"`` when it lies right of it.  Returns
    ``BND`` across chromosomes; on one chromosome ``DEL`` for (+, -),
    ``DUP`` for (-, +), ``INV`` for (+, +) or (-, -), and ``INTRA`` when
    the orientation is unknown.
    """
    if region_a[0] != region_b[0]:
        return "BND"
    return _SV_TYPE_BY_STRANDS.get(strands, "INTRA")


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

# Structural-variant evidence recorded for each informative read.
_SV_MIN_INDEL = 50  # bp: CIGAR insertions/deletions at least this long
_SV_MIN_CLIP = 20   # bp: soft clips at least this long mark a breakpoint


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


def _reset_scan_worker():
    """Drop the state set by :func:`_init_scan_worker`.

    Needed when the scan runs in-process rather than in a pool worker, so
    the automaton or query cache does not outlive the scan.
    """
    global _worker_automaton, _worker_jf_query, _worker_kmer_size
    global _worker_min_distinct_kmers_per_read

    _worker_automaton = None
    _worker_jf_query = None
    _worker_kmer_size = None
    _worker_min_distinct_kmers_per_read = 1


class _InformativeReadWriter:
    """Write informative reads to a BAM file, tagged ``dk:i:1``.

    The file is created when the first read arrives, so a contig without
    informative reads leaves no file behind.
    """

    def __init__(self, path, header):
        self.path = path
        self._header = header
        self._bam = None

    def write(self, read):
        if self._bam is None:
            self._bam = pysam.AlignmentFile(
                self.path, "wb", header=self._header,
            )
        read.set_tag("dk", 1, value_type="i")
        self._bam.write(read)

    def close(self):
        if self._bam is not None:
            self._bam.close()


def _process_informative_read(read, unique_in_read, kmer_hit_indices,
                              kmer_size, reads_seen, read_hits,
                              read_sv_meta, kmer_coverage, read_coverage,
                              bam_out=None):
    """Record an informative read's hits, coverage, and SV metadata.

    When *bam_out* is given (an :class:`_InformativeReadWriter`), the read
    is also written to it, whether mapped or not.

    An unmapped read placed next to its mapped mate is recorded in
    *read_sv_meta* only by that position (``placed_unmapped``), and its
    mate's mapping quality when the MQ tag gives it: evidence of
    unmappable sequence at the locus.

    Returns 1 if the read is unmapped-informative, 0 otherwise.
    Mutates *reads_seen*, *read_hits*, *read_sv_meta*, *kmer_coverage*,
    and *read_coverage* in place.
    """
    dedup_key = (read.query_name, read.is_supplementary)
    if dedup_key in reads_seen:
        return 0

    reads_seen.add(dedup_key)
    if bam_out is not None:
        bam_out.write(read)
    if read.is_unmapped:
        if read.reference_id >= 0:
            read_sv_meta[dedup_key] = {
                "placed_unmapped": (read.reference_name, read.reference_start),
                "mate_mapq": (read.get_tag("MQ") if read.has_tag("MQ")
                              else None),
            }
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
    # update(), not +=: Counter.__iadd__ rescans the whole Counter on
    # every call, which makes a contig scan quadratic.
    kmer_coverage[chrom].update(cov)
    # Count one read per touched position
    for pos in cov:
        read_coverage[chrom][pos] += 1

    # Collect SV metadata for this informative read
    sc_left, sc_right = _extract_softclips(read.cigartuples)
    # Where long soft clips start, candidate breakpoints, and on which side
    # of the alignment: "-" before its start, "+" after its end
    clips = []
    if sc_left >= _SV_MIN_CLIP:
        clips.append((read.reference_start, "-"))
    if sc_right >= _SV_MIN_CLIP:
        clips.append((read.reference_end, "+"))
    # A discordant pair's mate position and strand can link two breakpoint
    # regions; its mapping quality is known only when the MQ tag is present.
    mate = None
    if (read.is_paired and not read.mate_is_unmapped
            and not read.is_proper_pair):
        mate = (
            read.next_reference_name, read.next_reference_start,
            read.get_tag("MQ") if read.has_tag("MQ") else None,
            read.mate_is_reverse,
        )
    # The longest CIGAR insertion or deletion of at least _SV_MIN_INDEL bp
    large_indel, longest = None, 0
    for op, length in read.cigartuples or ():
        if op in (1, 2) and length >= _SV_MIN_INDEL and length > longest:
            large_indel, longest = ("INS" if op == 1 else "DEL"), length
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
        "max_clip": max(sc_left, sc_right),
        "mapq": read.mapping_quality,
        "pos": (chrom, read.reference_start),
        "end": read.reference_end,
        "is_reverse": read.is_reverse,
        "clip_side": _clip_side(read.cigartuples),
        "qspan": _query_span(read.cigartuples),
        "clips": clips,
        "large_indel": large_indel,  # "DEL", "INS" or None
        "mate": mate,
    }
    return 0


def _scan_contig_for_hits(child_bam, ref_fasta, contig, informative_bam=None):
    """Scan reads mapped to *contig* for proband-unique k-mers.

    When *contig* is ``None``, unmapped reads are scanned instead.  When
    *informative_bam* is given, every informative read recorded by this
    scan is also written to that path (created only if there is at least
    one such read), tagged ``dk:i:1``.

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

    bam_out = None
    if informative_bam is not None:
        bam_out = _InformativeReadWriter(informative_bam, bam.header)

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
                    bam_out=bam_out,
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
                    bam_out=bam_out,
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
                bam_out=bam_out,
            )

    if bam_out is not None:
        bam_out.close()
    bam.close()
    return (read_hits, reads_seen, unmapped_informative,
            total_reads_scanned, read_sv_meta, kmer_coverage,
            read_coverage)
