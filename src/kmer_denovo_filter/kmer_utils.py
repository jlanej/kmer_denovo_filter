"""K-mer utility functions."""

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
    rc = reverse_complement(kmer)
    return min(kmer, rc)


def extract_variant_spanning_kmers(
    read, variant_pos, k, min_baseq=0, ref=None, alt=None,
):
    """Extract canonical k-mers from a read that span the variant position.

    Args:
        read: pysam AlignedSegment
        variant_pos: 0-based reference position of the variant
        k: k-mer size
        min_baseq: Minimum base quality threshold
        ref: Reference allele string (for INDEL handling)
        alt: Alternate allele string (for INDEL handling)

    Returns:
        Set of canonical k-mer strings spanning the variant position.
    """
    aligned_pairs = read.get_aligned_pairs(matches_only=True)
    read_pos_at_variant = None
    for qpos, rpos in aligned_pairs:
        if rpos == variant_pos:
            read_pos_at_variant = qpos
            break

    if read_pos_at_variant is None:
        return set()

    seq = read.query_sequence
    quals = read.query_qualities
    if seq is None:
        return set()

    # For insertions the variant occupies len(alt) bases in the read.
    # Extend the window so k-mers spanning the right junction are captured.
    alt_len = len(alt) if alt and not _is_symbolic(alt) else 1
    variant_end_in_read = read_pos_at_variant + alt_len - 1

    kmers = set()
    start_min = max(0, read_pos_at_variant - k + 1)
    start_max = min(len(seq) - k, variant_end_in_read)

    for s in range(start_min, start_max + 1):
        kmer = seq[s:s + k]
        if "N" in kmer or "n" in kmer:
            continue
        if quals is not None and min_baseq > 0:
            if any(q < min_baseq for q in quals[s:s + k]):
                continue
        kmers.add(canonicalize(kmer))

    return kmers
