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

    Extracts the exact read sequence aligned to the reference span of the
    variant and compares it strictly to the candidate alternate allele.
    Handles SNPs, MNPs, insertions, deletions, and complex indels natively.

    Returns ``False`` for symbolic alleles or when *alt* is ``None``.
    """
    if alt is None or _is_symbolic(alt):
        return False

    seq = read.query_sequence
    if seq is None:
        return False

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
                extracted_seq.append(seq[qpos])

    # If the variant region was skipped entirely due to read boundaries
    if not in_variant_region:
        return False

    return "".join(extracted_seq).upper() == alt.upper()


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
