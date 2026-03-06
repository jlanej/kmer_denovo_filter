"""K-mer utility functions."""

import logging

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


def estimate_automaton_memory_gb(n_kmers, kmer_size=31):
    """Estimate the in-memory size of an Aho-Corasick automaton in GB.

    Empirically measured at ~928 bytes per pattern for 31-mer DNA
    sequences (including trie node, failure links, and Python string
    value overhead).  With forward + reverse-complement, the pattern
    count is roughly ``2 × n_kmers``.

    Args:
        n_kmers: Number of canonical k-mers.
        kmer_size: Length of each k-mer (default 31).

    Returns:
        Estimated memory in GB (float).
    """
    n_patterns = n_kmers * 2  # fwd + RC
    bytes_per_pattern = 928  # empirically measured
    total_bytes = n_patterns * bytes_per_pattern
    return total_bytes / (1024**3)


def _format_elapsed_simple(seconds):
    """Format elapsed seconds as a human-readable string."""
    if seconds < 60:
        return f"{seconds:.1f}s"
    if seconds < 3600:
        minutes = int(seconds // 60)
        secs = seconds % 60
        return f"{minutes}m {secs:.1f}s"
    hours = int(seconds // 3600)
    minutes = int((seconds % 3600) // 60)
    return f"{hours}h {minutes}m"


def read_supports_alt(read, variant_pos, ref, alt, *, aligned_pairs=None, seq=None):
    """Return True if *read* carries the alternate allele at *variant_pos*.

    Extracts the exact read sequence aligned to the reference span of the
    variant and compares it strictly to the candidate alternate allele.
    Handles SNPs, MNPs, insertions, deletions, and complex indels natively.

    Returns ``False`` for symbolic alleles or when *alt* is ``None``.

    Args:
        aligned_pairs: Optional pre-computed result of
            ``read.get_aligned_pairs(matches_only=False)``.  Computed from
            *read* when not provided.
        seq: Optional pre-decoded ``read.query_sequence``.  Decoded from
            *read* when not provided.
    """
    if alt is None or _is_symbolic(alt):
        return False

    if seq is None:
        seq = read.query_sequence
    if seq is None:
        return False

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
