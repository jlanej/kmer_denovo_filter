"""BAM/CRAM scanning engine and alignment helpers.

Provides the k-mer hit scanner of discovery mode's Module 3.  The
proband-unique k-mers live in a
:class:`~kmer_denovo_filter.core.kmer_set.KmerSet` (a sorted array of
2-bit codes, loaded once by the parent and shared copy-on-write by the
forked workers), and every batch of reads is tested against it with
vectorized ``numpy`` operations: no per-k-mer Python loop and no
subprocess per batch.

No function in this module accepts an ``argparse.Namespace`` object.
"""

import collections
import logging
import re

import pysam

from kmer_denovo_filter.core.kmer_set import KmerSet

logger = logging.getLogger(__name__)


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
# multiprocessing pool initialiser.  With the fork start method the
# parent builds the k-mer set before the pool starts and every worker
# shares its array copy-on-write; with spawn each worker gets a copy.
_worker_kmer_set = None        # KmerSet of the proband-unique k-mers
_worker_kmer_size = None
_worker_min_distinct_kmers_per_read = 1

# Reads are scanned in batches: a batch ends at this many reads or this
# many bases, whichever comes first.  The batch's transient arrays are a
# few tens of bytes per base, so 1.5 million bases cost about 100 MB.
_SCAN_BATCH_READS = 10_000
_SCAN_BATCH_BASES = 1_500_000

# Structural-variant evidence recorded for each informative read.
_SV_MIN_INDEL = 50  # bp: CIGAR insertions/deletions at least this long
_SV_MIN_CLIP = 20   # bp: soft clips at least this long mark a breakpoint


def _init_scan_worker(kmer_set, kmer_size, min_distinct_kmers_per_read=1):
    """Initializer for the scan workers.

    *kmer_set* may be:
    - A :class:`KmerSet`, as the discovery pipeline builds it once in the
      parent (shared copy-on-write by forked workers).
    - A FASTA file path — loaded into a :class:`KmerSet`.
    - An iterable of k-mer strings — encoded into a :class:`KmerSet`.
    """
    global _worker_kmer_set, _worker_kmer_size
    global _worker_min_distinct_kmers_per_read

    if isinstance(kmer_set, KmerSet):
        _worker_kmer_set = kmer_set
    elif isinstance(kmer_set, str):
        _worker_kmer_set = KmerSet.from_fasta(kmer_set, kmer_size)
    else:
        _worker_kmer_set = KmerSet.from_kmers(kmer_set or (), kmer_size)

    _worker_kmer_size = kmer_size
    _worker_min_distinct_kmers_per_read = min_distinct_kmers_per_read


def _reset_scan_worker():
    """Drop the state set by :func:`_init_scan_worker`.

    Needed when the scan runs in-process rather than in a pool worker, so
    the k-mer set does not outlive the scan.
    """
    global _worker_kmer_set, _worker_kmer_size
    global _worker_min_distinct_kmers_per_read

    _worker_kmer_set = None
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
                              bam_out=None, deferred_coverage=None,
                              window=None):
    """Record an informative read's hits, coverage, and SV metadata.

    When *bam_out* is given (an :class:`_InformativeReadWriter`), the read
    is also written to it, whether mapped or not.

    An unmapped read placed next to its mapped mate is recorded in
    *read_sv_meta* only by that position (``placed_unmapped``), and its
    mate's mapping quality when the MQ tag gives it: evidence of
    unmappable sequence at the locus.

    A read whose de-duplication key may also be reported by another
    window of the same contig - a supplementary alignment, or a read
    whose mate maps to this contig outside *window* ``(start, end)`` -
    has its coverage put in *deferred_coverage* (as ``(key, coverage)``)
    instead of *kmer_coverage* and *read_coverage*, when that list is
    given: the caller adds it only if the key was not seen by an earlier
    window, exactly as a scan of the whole contig would have skipped the
    read.

    Returns 1 if the read is unmapped-informative, 0 otherwise.
    Mutates *reads_seen*, *read_hits*, *read_sv_meta*, *kmer_coverage*,
    *read_coverage* and *deferred_coverage* in place.
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
    if deferred_coverage is not None and _key_may_recur_in_contig(read, window):
        deferred_coverage.append((dedup_key, cov))
    else:
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


def _key_may_recur_in_contig(read, window):
    """Whether another record of this contig, outside *window*, may carry
    the read's de-duplication key: its mate (same name, both primary) when
    it maps to this contig outside the window - or when the record does
    not say where the mate is - or, for a supplementary alignment, another
    supplementary part of the read anywhere."""
    if window is None:
        return False
    if read.is_supplementary:
        return True
    if not read.is_paired:
        return False
    if read.next_reference_id < 0:
        return True
    if read.next_reference_id != read.reference_id:
        return False
    start, end = window
    return read.next_reference_start < start or read.next_reference_start >= end


def _scan_contig_for_hits(child_bam, ref_fasta, contig, informative_bam=None,
                          start=None, end=None):
    """Scan the reads of *contig* (or of a window of it) for proband-unique
    k-mers.

    When *contig* is ``None``, unmapped reads are scanned instead.  With
    *start* and *end*, only the reads whose alignment starts in
    ``[start, end)`` are scanned (a read overlapping the window from before
    it belongs to the window of its start), so the windows of a contig
    together scan each read once, in the order of the file.  When
    *informative_bam* is given, every informative read recorded by this
    scan is also written to that path (created only if there is at least
    one such read), tagged ``dk:i:1``.

    Reads are taken in batches (:data:`_SCAN_BATCH_READS` reads or
    :data:`_SCAN_BATCH_BASES` bases) and each batch is scanned at once
    with :meth:`KmerSet.scan`; the informative reads of a batch are then
    recorded in file order, so the result is that of scanning the reads
    one by one.

    Returns:
        (read_hits, reads_seen, unmapped_informative, total_reads_scanned,
         read_sv_meta, kmer_coverage, read_coverage, deferred_coverage)

    ``read_hits`` holds, per informative mapped read, its contig, start,
    end, name, the set of proband-unique k-mers it carries (as 2-bit
    codes: only their number is used downstream) and whether it is
    supplementary.

    ``deferred_coverage`` is empty for a whole-contig scan.  For a window
    it holds the coverage of the reads whose key another window of the
    contig may also report (see :func:`_process_informative_read`), as
    ``(key, Counter)`` in file order, for the caller to add once the
    windows are folded in order.

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
    kmer_set = _worker_kmer_set
    kmer_size = _worker_kmer_size
    min_dk_per_read = _worker_min_distinct_kmers_per_read
    if kmer_set is None:
        kmer_set = KmerSet(kmer_size or 1)
    bam = pysam.AlignmentFile(
        child_bam, reference_filename=ref_fasta if ref_fasta else None,
    )

    read_hits = []
    reads_seen = set()
    read_sv_meta = {}
    kmer_coverage = collections.defaultdict(collections.Counter)
    read_coverage = collections.defaultdict(collections.Counter)
    deferred_coverage = []
    unmapped_informative = 0
    total_reads_scanned = 0
    window = None

    if contig is None:
        try:
            iterator = bam.fetch("*")
        except (ValueError, KeyError):
            bam.close()
            return (read_hits, reads_seen, unmapped_informative,
                    total_reads_scanned, read_sv_meta, kmer_coverage,
                    read_coverage, deferred_coverage)
    elif start is None and end is None:
        iterator = bam.fetch(contig=contig)
    else:
        iterator = bam.fetch(contig=contig, start=start, stop=end)
        window = (start, end)

    bam_out = None
    if informative_bam is not None:
        bam_out = _InformativeReadWriter(informative_bam, bam.header)

    pending = []        # the reads of the batch, in file order
    pending_seqs = []   # their sequences, as bytes
    pending_bases = 0

    def _scan_batch():
        nonlocal unmapped_informative
        for index, positions, codes in kmer_set.scan(pending_seqs):
            unique_in_read = set(
                codes.tolist() if hasattr(codes, "tolist") else codes)
            if len(unique_in_read) < min_dk_per_read:
                continue
            unmapped_informative += _process_informative_read(
                pending[index], unique_in_read, set(positions.tolist()),
                kmer_size, reads_seen, read_hits, read_sv_meta,
                kmer_coverage, read_coverage, bam_out=bam_out,
                deferred_coverage=deferred_coverage, window=window,
            )

    for read in iterator:
        if read.is_secondary:
            continue
        if read.is_duplicate:
            continue
        if start is not None and read.reference_start < start:
            continue  # scanned by the window its alignment starts in

        total_reads_scanned += 1
        seq = read.query_sequence
        if seq is None or len(seq) < kmer_size:
            continue

        pending.append(read)
        pending_seqs.append(seq.encode("ascii"))
        pending_bases += len(seq)
        if (len(pending) >= _SCAN_BATCH_READS
                or pending_bases >= _SCAN_BATCH_BASES):
            _scan_batch()
            pending = []
            pending_seqs = []
            pending_bases = 0

    if pending:
        _scan_batch()

    if bam_out is not None:
        bam_out.close()
    bam.close()
    return (read_hits, reads_seen, unmapped_informative,
            total_reads_scanned, read_sv_meta, kmer_coverage,
            read_coverage, deferred_coverage)
