"""Tests for k-mer utility functions."""

from kmer_denovo_filter.kmer_utils import (
    canonicalize,
    extract_variant_spanning_kmers,
    reverse_complement,
)


class TestReverseComplement:
    def test_palindrome(self):
        assert reverse_complement("ACGT") == "ACGT"

    def test_poly_a(self):
        assert reverse_complement("AAAA") == "TTTT"

    def test_mixed(self):
        assert reverse_complement("AACG") == "CGTT"

    def test_single_base(self):
        assert reverse_complement("A") == "T"
        assert reverse_complement("C") == "G"

    def test_lowercase(self):
        assert reverse_complement("acgt") == "acgt"


class TestCanonicalize:
    def test_palindrome(self):
        assert canonicalize("ACGT") == "ACGT"

    def test_returns_smaller(self):
        assert canonicalize("TTTT") == "AAAA"
        assert canonicalize("AAAA") == "AAAA"

    def test_forward_is_smaller(self):
        assert canonicalize("ACCC") == "ACCC"

    def test_revcomp_is_smaller(self):
        assert canonicalize("GGGT") == "ACCC"


class TestExtractVariantSpanningKmers:
    """Tests using a mock read object."""

    class MockRead:
        """Minimal mock of pysam.AlignedSegment."""

        def __init__(self, seq, aligned_pairs, quals=None):
            self.query_sequence = seq
            self._aligned_pairs = aligned_pairs
            self.query_qualities = quals

        def get_aligned_pairs(self):
            return self._aligned_pairs

    def test_simple_spanning(self):
        # Read: ACGTACGT aligned at ref pos 100-107
        seq = "ACGTACGT"
        pairs = [(i, 100 + i) for i in range(8)]
        read = self.MockRead(seq, pairs)

        # k=4, variant at ref pos 102 (read pos 2)
        kmers = extract_variant_spanning_kmers(read, 102, 4, min_baseq=0)
        # K-mers spanning pos 2: start 0 (pos 0-3), start 1 (pos 1-4), start 2 (pos 2-5)
        # But we also need start_max = min(8-4, 2) = 2
        # start_min = max(0, 2-4+1) = 0
        # So starts: 0, 1, 2
        assert len(kmers) == 3

    def test_no_overlap(self):
        seq = "ACGT"
        pairs = [(i, 100 + i) for i in range(4)]
        read = self.MockRead(seq, pairs)

        # Variant at pos 200, not overlapping
        kmers = extract_variant_spanning_kmers(read, 200, 3, min_baseq=0)
        assert len(kmers) == 0

    def test_quality_filter(self):
        seq = "ACGTACGT"
        pairs = [(i, 100 + i) for i in range(8)]
        # Low quality at position 2
        quals = [30, 30, 5, 30, 30, 30, 30, 30]
        read = self.MockRead(seq, pairs, quals)

        kmers = extract_variant_spanning_kmers(read, 102, 4, min_baseq=20)
        # All k-mers that include position 2 should be filtered
        assert len(kmers) == 0

    def test_n_in_kmer(self):
        seq = "ACNTACGT"
        pairs = [(i, 100 + i) for i in range(8)]
        read = self.MockRead(seq, pairs)

        kmers = extract_variant_spanning_kmers(read, 102, 4, min_baseq=0)
        # K-mers containing N should be excluded
        # start 0: ACNT - has N
        # start 1: CNTA - has N
        # start 2: NTAC - has N
        assert len(kmers) == 0

    def test_canonical_output(self):
        seq = "TTTTAAAA"
        pairs = [(i, 100 + i) for i in range(8)]
        read = self.MockRead(seq, pairs)

        kmers = extract_variant_spanning_kmers(read, 103, 4, min_baseq=0)
        # All k-mers should be canonicalized
        for kmer in kmers:
            rc = reverse_complement(kmer)
            assert kmer <= rc

    def test_none_sequence(self):
        pairs = [(i, 100 + i) for i in range(8)]
        read = self.MockRead(None, pairs)
        kmers = extract_variant_spanning_kmers(read, 102, 4, min_baseq=0)
        assert len(kmers) == 0

    def test_deletion(self):
        """K-mers should span the deletion junction."""
        # Read: ACGTCGT (7 bases), deletion of base at ref pos 104
        # Reference would be: A C G T X C G T aligned at 100-107
        # Read has deletion at 104: aligned at 100-103, then 105-107
        seq = "ACGTCGT"
        pairs = [
            (0, 100), (1, 101), (2, 102), (3, 103),
            (None, 104),  # deleted base
            (4, 105), (5, 106), (6, 107),
        ]
        read = self.MockRead(seq, pairs)

        # Variant anchor at ref pos 103 (the base before deletion),
        # REF=TX, ALT=T
        kmers = extract_variant_spanning_kmers(
            read, 103, 4, min_baseq=0, ref="TX", alt="T",
        )
        # read_pos_at_variant = 3 (qpos for rpos 103)
        # alt_len = 1, variant_end_in_read = 3
        # start_min = max(0, 3-4+1) = 0
        # start_max = min(7-4, 3) = 3
        # 4 k-mers: starts 0,1,2,3
        assert len(kmers) == 4
        # K-mers at starts 1,2,3 span the junction (include pos 3 and 4+)
        # All capture the novel deletion sequence

    def test_insertion(self):
        """K-mers should span the full insertion including right junction."""
        # Read: 12 bases with 4-base insertion in the middle
        # Ref positions: 100,101,102,103, then insertion, 104,105,106,107
        seq = "ACTGCATATCGA"
        pairs = [
            (0, 100), (1, 101), (2, 102), (3, 103),
            (4, None), (5, None), (6, None), (7, None),  # inserted bases
            (8, 104), (9, 105), (10, 106), (11, 107),
        ]
        read = self.MockRead(seq, pairs)

        # Variant anchor at ref pos 103, REF=G, ALT=GCATA
        kmers = extract_variant_spanning_kmers(
            read, 103, 4, min_baseq=0, ref="G", alt="GCATA",
        )
        # read_pos_at_variant = 3
        # alt_len = 5, variant_end_in_read = 3 + 5 - 1 = 7
        # start_min = max(0, 3-4+1) = 0
        # start_max = min(12-4, 7) = 7
        # 8 k-mers: starts 0..7 (= k + alt_len - 1 = 4 + 5 - 1 = 8)
        assert len(kmers) == 8

    def test_insertion_without_ref_alt(self):
        """Without ref/alt, only k k-mers are extracted (old behavior)."""
        seq = "ACTGCATATCGA"
        pairs = [
            (0, 100), (1, 101), (2, 102), (3, 103),
            (4, None), (5, None), (6, None), (7, None),
            (8, 104), (9, 105), (10, 106), (11, 107),
        ]
        read = self.MockRead(seq, pairs)

        # No ref/alt → falls back to alt_len=1 (old behaviour)
        kmers = extract_variant_spanning_kmers(read, 103, 4, min_baseq=0)
        assert len(kmers) == 4  # only k k-mers
