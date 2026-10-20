"""K-mer utility functions."""

import logging

# The non-human contamination engine (Kraken2Runner) and allele-support
# helpers (read_supports_alt, _is_symbolic) now live in the standalone
# ``nonhuman-screen`` package. They are re-exported here so existing imports
# — and tests — that reference ``kmer_denovo_filter.kmer_utils`` keep working
# unchanged.
from nonhuman_screen.engine import (  # noqa: F401
    Kraken2Runner,
    _ARCHAEA_TAXID,
    _BACTERIA_TAXID,
    _EUKARYOTA_TAXID,
    _FUNGI_TAXID,
    _HUMAN_TAXID,
    _METAZOA_TAXID,
    _UNIVEC_CORE_TAXID,
    _VIRIDIPLANTAE_TAXID,
    _VIRUSES_TAXID,
    _read_proc_rss_kb,
)
from nonhuman_screen.alleles import _is_symbolic, read_supports_alt  # noqa: F401

logger = logging.getLogger(__name__)

_COMP = str.maketrans("ACGTacgt", "TGCAtgca")

def reverse_complement(seq):
    """Return the reverse complement of a DNA sequence."""
    return seq.translate(_COMP)[::-1]


def canonicalize(kmer):
    """Return the canonical (lexicographically smaller) form of a k-mer."""
    rc = kmer.translate(_COMP)[::-1]
    return kmer if kmer < rc else rc


# ── Per-read k-mer extraction (the reference implementation) ───────


def _extract_read_kmers(seq, kmer_size):
    """Extract canonicalized k-mers from a read sequence.

    This is the per-k-mer reference that discovery mode's vectorized scan
    (:mod:`kmer_denovo_filter.core.kmer_set`) is checked against: the
    windows of the upper-cased read, those holding an N skipped, each
    canonicalized.

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
