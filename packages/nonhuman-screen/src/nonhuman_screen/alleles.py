"""Allele-support determination for reads at a variant locus.

Pure-Python helpers that decide whether an aligned read carries a given
alternate allele.  They operate on a pysam ``AlignedSegment`` passed in by
the caller but do **not** import pysam themselves, so this module has no
third-party dependency.
"""


def _is_symbolic(allele):
    """Return True if *allele* is a symbolic VCF allele with no literal sequence.

    Symbolic alleles include ``<DEL>``, ``<INS>``, ``<DUP>``, breakend
    notation containing ``[`` or ``]``, and the overlapping-deletion
    marker ``*``.
    """
    if not allele:
        return True
    return allele[0] == "<" or allele == "*" or "[" in allele or "]" in allele


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

