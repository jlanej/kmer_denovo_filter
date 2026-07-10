"""Unit tests for allele-support determination.

``read_supports_alt`` accepts pre-computed ``seq`` and ``aligned_pairs``, so it
can be exercised without a pysam read object.
"""

from nonhuman_screen.alleles import _is_symbolic, read_supports_alt


def _pairs(start, length):
    """Fully-aligned (qpos, rpos) pairs for a length-*length* read at *start*."""
    return [(i, start + i) for i in range(length)]


class TestIsSymbolic:
    def test_symbolic_and_literal(self):
        assert _is_symbolic("<DEL>")
        assert _is_symbolic("*")
        assert _is_symbolic("N[chr2:321]")
        assert _is_symbolic("")
        assert _is_symbolic(None)
        assert not _is_symbolic("A")
        assert not _is_symbolic("ACGT")


class TestReadSupportsAlt:
    def test_snp_support(self):
        # read of 10 bp aligned at ref 100..109; base at ref 105 is 'T'
        seq = "AAAAA" "T" "AAAA"
        assert read_supports_alt(
            None, 105, "A", "T", seq=seq, aligned_pairs=_pairs(100, 10),
        )

    def test_snp_no_support(self):
        seq = "AAAAA" "G" "AAAA"  # base at 105 is 'G', not the 'T' alt
        assert not read_supports_alt(
            None, 105, "A", "T", seq=seq, aligned_pairs=_pairs(100, 10),
        )

    def test_symbolic_alt_never_supports(self):
        seq = "AAAAAAAAAA"
        assert not read_supports_alt(
            None, 105, "A", "<DEL>", seq=seq, aligned_pairs=_pairs(100, 10),
        )

    def test_none_alt_never_supports(self):
        seq = "AAAAAAAAAA"
        assert not read_supports_alt(
            None, 105, "A", None, seq=seq, aligned_pairs=_pairs(100, 10),
        )

    def test_variant_region_beyond_read_returns_false(self):
        # read covers ref 100..109; variant at 200 is not spanned
        seq = "AAAAAAAAAA"
        assert not read_supports_alt(
            None, 200, "A", "T", seq=seq, aligned_pairs=_pairs(100, 10),
        )

    def test_deletion_support(self):
        # ref allele "AGG" (pos 105..107) deleted -> alt "A".
        # aligned_pairs: read base present at rpos 105, deleted (qpos None) at 106,107
        pairs = [(i, 100 + i) for i in range(6)]  # rpos 100..105 -> qpos 0..5
        pairs += [(None, 106), (None, 107)]        # deleted ref bases
        pairs += [(6, 108), (7, 109)]              # resume
        seq = "AAAAAA" "AA"  # qpos 0..7
        assert read_supports_alt(
            None, 105, "AGG", "A", seq=seq, aligned_pairs=pairs,
        )
