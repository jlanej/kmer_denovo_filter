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


def read_supports_alt(read, variant_pos, ref, alt):
    """Return True if *read* carries the alternate allele at *variant_pos*.

    For SNPs the read base must match *alt*.  For deletions the reference
    bases beyond the anchor must be absent from the read.  For insertions
    the inserted sequence must be present following the anchor.

    Returns ``False`` for symbolic alleles or when *alt* is ``None``.
    """
    if alt is None or _is_symbolic(alt):
        return False

    seq = read.query_sequence
    if seq is None:
        return False

    aligned_pairs = read.get_aligned_pairs(matches_only=False)

    # Build ref-pos → read-pos map (only for aligned positions)
    ref_to_qpos = {}
    for qpos, rpos in aligned_pairs:
        if rpos is not None and qpos is not None:
            ref_to_qpos[rpos] = qpos

    anchor_qpos = ref_to_qpos.get(variant_pos)
    if anchor_qpos is None:
        return False

    if len(ref) == 1 and len(alt) == 1:
        # SNP
        return seq[anchor_qpos].upper() == alt.upper()

    if len(alt) < len(ref):
        # Deletion – anchor base must match and deleted ref positions must
        # be absent from the read alignment.
        if seq[anchor_qpos].upper() != alt[0].upper():
            return False
        for i in range(len(alt), len(ref)):
            if ref_to_qpos.get(variant_pos + i) is not None:
                return False
        return True

    if len(alt) > len(ref):
        # Insertion – the inserted bases must follow the anchor in the read.
        ins_bases = alt[len(ref):]
        for i, base in enumerate(ins_bases):
            rp = anchor_qpos + len(ref) + i
            if rp >= len(seq):
                return False
            if seq[rp].upper() != base.upper():
                return False
        return True

    return False


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
