"""Tests for k-mer utility functions."""

from kmer_denovo_filter.kmer_utils import (
    _extract_read_kmers,
    _is_symbolic,
    build_kmer_automaton,
    canonicalize,
    extract_variant_spanning_kmers,
    read_supports_alt,
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

        def get_aligned_pairs(self, matches_only=False):
            if matches_only:
                return [
                    (q, r) for q, r in self._aligned_pairs
                    if q is not None and r is not None
                ]
            return self._aligned_pairs

        def get_reference_positions(self, full_length=False):
            if full_length:
                return [r for q, r in self._aligned_pairs if q is not None]
            return [r for q, r in self._aligned_pairs
                    if q is not None and r is not None]

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

    def test_symbolic_del_treated_as_snp(self):
        """Symbolic <DEL> allele should not crash; treated like a SNP."""
        seq = "ACGTACGT"
        pairs = [(i, 100 + i) for i in range(8)]
        read = self.MockRead(seq, pairs)

        kmers = extract_variant_spanning_kmers(
            read, 102, 4, min_baseq=0, ref="A", alt="<DEL>",
        )
        # alt_len defaults to 1 for symbolic alleles → same as SNP
        assert len(kmers) == 3

    def test_symbolic_ins_treated_as_snp(self):
        """Symbolic <INS> allele should not crash; treated like a SNP."""
        seq = "ACGTACGT"
        pairs = [(i, 100 + i) for i in range(8)]
        read = self.MockRead(seq, pairs)

        kmers = extract_variant_spanning_kmers(
            read, 102, 4, min_baseq=0, ref="A", alt="<INS>",
        )
        assert len(kmers) == 3

    def test_star_allele_treated_as_snp(self):
        """Star (*) allele should not crash; treated like a SNP."""
        seq = "ACGTACGT"
        pairs = [(i, 100 + i) for i in range(8)]
        read = self.MockRead(seq, pairs)

        kmers = extract_variant_spanning_kmers(
            read, 102, 4, min_baseq=0, ref="A", alt="*",
        )
        assert len(kmers) == 3

    def test_breakend_allele_treated_as_snp(self):
        """Breakend notation should not crash; treated like a SNP."""
        seq = "ACGTACGT"
        pairs = [(i, 100 + i) for i in range(8)]
        read = self.MockRead(seq, pairs)

        kmers = extract_variant_spanning_kmers(
            read, 102, 4, min_baseq=0, ref="A", alt="A]chr2:100]",
        )
        assert len(kmers) == 3


class TestIsSymbolic:
    def test_symbolic_del(self):
        assert _is_symbolic("<DEL>") is True

    def test_symbolic_ins(self):
        assert _is_symbolic("<INS>") is True

    def test_symbolic_dup(self):
        assert _is_symbolic("<DUP>") is True

    def test_symbolic_inv(self):
        assert _is_symbolic("<INV>") is True

    def test_symbolic_cnv(self):
        assert _is_symbolic("<CNV>") is True

    def test_star_allele(self):
        assert _is_symbolic("*") is True

    def test_breakend_right(self):
        assert _is_symbolic("A]chr2:100]") is True

    def test_breakend_left(self):
        assert _is_symbolic("[chr2:100[A") is True

    def test_normal_snp(self):
        assert _is_symbolic("A") is False

    def test_normal_insertion(self):
        assert _is_symbolic("ACGT") is False

    def test_none(self):
        assert _is_symbolic(None) is True

    def test_empty_string(self):
        assert _is_symbolic("") is True


class TestReadSupportsAlt:
    """Tests for read_supports_alt using mock read objects."""

    class MockRead:
        """Minimal mock of pysam.AlignedSegment."""

        def __init__(self, seq, aligned_pairs, quals=None):
            self.query_sequence = seq
            self._aligned_pairs = aligned_pairs
            self.query_qualities = quals

        def get_aligned_pairs(self, matches_only=False):
            if matches_only:
                return [
                    (q, r) for q, r in self._aligned_pairs
                    if q is not None and r is not None
                ]
            return self._aligned_pairs

    def test_snp_supports_alt(self):
        # Read has G at ref pos 102 (index 2 in "ACGTACGT"), alt is G
        seq = "ACGTACGT"
        pairs = [(i, 100 + i) for i in range(8)]
        read = self.MockRead(seq, pairs)
        assert read_supports_alt(read, 102, "A", "G") is True

    def test_snp_does_not_support_alt(self):
        seq = "ACGTACGT"
        pairs = [(i, 100 + i) for i in range(8)]
        read = self.MockRead(seq, pairs)
        # Read has G at pos 102, alt is T
        assert read_supports_alt(read, 102, "A", "T") is False

    def test_deletion_supports_alt(self):
        # Read: ACGTCGT (7 bases), deletion of base at ref pos 104
        seq = "ACGTCGT"
        pairs = [
            (0, 100), (1, 101), (2, 102), (3, 103),
            (None, 104),  # deleted base
            (4, 105), (5, 106), (6, 107),
        ]
        read = self.MockRead(seq, pairs)
        # REF=TX, ALT=T at anchor pos 103
        assert read_supports_alt(read, 103, "TX", "T") is True

    def test_deletion_does_not_support_alt(self):
        # Read has all bases aligned (no deletion)
        seq = "ACGTACGT"
        pairs = [(i, 100 + i) for i in range(8)]
        read = self.MockRead(seq, pairs)
        # REF=TA, ALT=T (deletion) but read has the base at 104
        assert read_supports_alt(read, 103, "TA", "T") is False

    def test_insertion_supports_alt(self):
        # Read: ACTGCATATCGA with 4-base insertion
        seq = "ACTGCATATCGA"
        pairs = [
            (0, 100), (1, 101), (2, 102), (3, 103),
            (4, None), (5, None), (6, None), (7, None),  # inserted
            (8, 104), (9, 105), (10, 106), (11, 107),
        ]
        read = self.MockRead(seq, pairs)
        # REF=G, ALT=GCATA at anchor pos 103
        assert read_supports_alt(read, 103, "G", "GCATA") is True

    def test_insertion_does_not_support_alt(self):
        # Read matches reference (no insertion)
        seq = "ACGTACGT"
        pairs = [(i, 100 + i) for i in range(8)]
        read = self.MockRead(seq, pairs)
        # REF=T, ALT=TGCA but read doesn't have those inserted bases
        assert read_supports_alt(read, 103, "T", "TGCA") is False

    def test_symbolic_allele_returns_false(self):
        seq = "ACGTACGT"
        pairs = [(i, 100 + i) for i in range(8)]
        read = self.MockRead(seq, pairs)
        assert read_supports_alt(read, 102, "A", "<DEL>") is False

    def test_none_alt_returns_false(self):
        seq = "ACGTACGT"
        pairs = [(i, 100 + i) for i in range(8)]
        read = self.MockRead(seq, pairs)
        assert read_supports_alt(read, 102, "A", None) is False

    def test_none_sequence_returns_false(self):
        pairs = [(i, 100 + i) for i in range(8)]
        read = self.MockRead(None, pairs)
        assert read_supports_alt(read, 102, "A", "T") is False

    def test_variant_pos_not_in_read(self):
        seq = "ACGT"
        pairs = [(i, 100 + i) for i in range(4)]
        read = self.MockRead(seq, pairs)
        assert read_supports_alt(read, 200, "A", "T") is False

    def test_mnp_supports_alt(self):
        # Read has CG at ref pos 100-101, variant REF=AT ALT=CG
        seq = "CGAAA"
        pairs = [(i, 100 + i) for i in range(5)]
        read = self.MockRead(seq, pairs)
        assert read_supports_alt(read, 100, "AT", "CG") is True

    def test_mnp_does_not_support_alt(self):
        # Read has AT (reference) at ref pos 100-101
        seq = "ATAAA"
        pairs = [(i, 100 + i) for i in range(5)]
        read = self.MockRead(seq, pairs)
        assert read_supports_alt(read, 100, "AT", "CG") is False

    def test_mnp_partial_match_returns_false(self):
        # Read has CA at ref pos 100-101 (only first base matches alt)
        seq = "CAAAA"
        pairs = [(i, 100 + i) for i in range(5)]
        read = self.MockRead(seq, pairs)
        assert read_supports_alt(read, 100, "AT", "CG") is False

    def test_complex_deletion_supports_alt(self):
        # REF=ACGT ALT=AT — anchor + T remain, CG are deleted
        # Read: ATXXX with deletion of ref positions 102 and 103
        seq = "ATXXX"
        pairs = [
            (0, 100), (1, 101),
            (None, 102), (None, 103),  # deleted CG
            (2, 104), (3, 105), (4, 106),
        ]
        read = self.MockRead(seq, pairs)
        assert read_supports_alt(read, 100, "ACGT", "AT") is True

    def test_complex_deletion_wrong_tail_returns_false(self):
        # Same deletion but read has AC instead of AT after anchor
        seq = "ACXXX"
        pairs = [
            (0, 100), (1, 101),
            (None, 102), (None, 103),
            (2, 104), (3, 105), (4, 106),
        ]
        read = self.MockRead(seq, pairs)
        assert read_supports_alt(read, 100, "ACGT", "AT") is False

    def test_insertion_anchor_mismatch_returns_false(self):
        # Anchor base in read does not match alt[0]; should return False
        # seq[3] = 'X', which maps to ref pos 103 (the anchor) — 'X' != alt[0]='G'
        seq = "ACTXCATATCGA"
        pairs = [
            (0, 100), (1, 101), (2, 102), (3, 103),
            (4, None), (5, None), (6, None), (7, None),  # inserted
            (8, 104), (9, 105), (10, 106), (11, 107),
        ]
        read = self.MockRead(seq, pairs)
        # anchor at pos 103 maps to seq[3]='X', which does not match alt[0]='G'
        assert read_supports_alt(read, 103, "G", "GCATA") is False

    def test_insertion_tandem_repeat_no_false_positive(self):
        # Read perfectly matches reference (no insertion),
        # downstream ref resembles inserted bases — must NOT be a false positive.
        seq = "ACGTCATA"
        pairs = [(i, 100 + i) for i in range(8)]
        read = self.MockRead(seq, pairs)
        # REF=T, ALT=TCATA — but read has no insertion
        assert read_supports_alt(read, 103, "T", "TCATA") is False

    def test_low_quality_base_at_variant_position_rejects_alt_support(self):
        seq = "ACGTACGT"
        pairs = [(i, 100 + i) for i in range(8)]
        quals = [30, 30, 5, 30, 30, 30, 30, 30]
        read = self.MockRead(seq, pairs, quals)
        assert read_supports_alt(read, 102, "A", "G", min_baseq=20) is False

    def test_high_quality_alt_base_returns_true(self):
        seq = "ACGTACGT"
        pairs = [(i, 100 + i) for i in range(8)]
        quals = [30] * 8
        read = self.MockRead(seq, pairs, quals)
        assert read_supports_alt(read, 102, "A", "G", min_baseq=20) is True

    def test_min_baseq_zero_does_not_fetch_query_qualities(self):
        class LazyQualRead:
            query_sequence = "ACGTACGT"
            _aligned_pairs = [(i, 100 + i) for i in range(8)]

            @property
            def query_qualities(self):
                raise AssertionError("query_qualities should not be accessed")

            def get_aligned_pairs(self, matches_only=False):
                if matches_only:
                    return [
                        (q, r) for q, r in self._aligned_pairs
                        if q is not None and r is not None
                    ]
                return self._aligned_pairs

        read = LazyQualRead()
        assert read_supports_alt(read, 102, "A", "G", min_baseq=0) is True


class TestBuildKmerAutomaton:
    """Tests for the Aho-Corasick automaton builder."""

    def test_empty_set(self):
        """An empty k-mer set should produce None (no automaton)."""
        A = build_kmer_automaton(set())
        assert A is None

    def test_single_canonical_kmer(self):
        """A single canonical k-mer should be found in a matching sequence."""
        kmer = canonicalize("ACGTAC")
        A = build_kmer_automaton({kmer})
        hits = {v for _, v in A.iter("ACGTAC")}
        assert kmer in hits

    def test_finds_reverse_complement(self):
        """The automaton should match the reverse-complement of a canonical k-mer."""
        kmer = canonicalize("AACGT")   # canonical = "AACGT"
        rc = reverse_complement("AACGT")    # "ACGTT"
        A = build_kmer_automaton({kmer})
        # Searching for the reverse complement should return the canonical form
        hits = {v for _, v in A.iter(rc)}
        assert kmer in hits

    def test_palindrome_kmer(self):
        """A palindromic k-mer is its own reverse complement."""
        kmer = "ACGT"  # palindrome
        assert reverse_complement(kmer) == kmer
        A = build_kmer_automaton({kmer})
        hits = {v for _, v in A.iter("ACGT")}
        assert kmer in hits

    def test_multiple_kmers_in_sequence(self):
        """Multiple k-mers found in a single sequence."""
        kmers = {canonicalize("AAACCC"), canonicalize("CCCGGG")}
        A = build_kmer_automaton(kmers)
        seq = "AAACCCGGG"
        hits = {v for _, v in A.iter(seq)}
        assert hits == kmers

    def test_no_match(self):
        """No matches when k-mers are absent from the sequence."""
        A = build_kmer_automaton({canonicalize("TTTTT")})
        hits = list(A.iter("AACCC"))
        assert hits == []

    def test_returns_canonical_value(self):
        """The value stored for every match must be the canonical k-mer."""
        kmer = canonicalize("GGGTA")  # canonical = "GGGTA" or its rc
        A = build_kmer_automaton({kmer})
        rc = reverse_complement(kmer)
        # Search with whichever orientation matches
        for seq in (kmer, rc):
            for _, val in A.iter(seq):
                assert val == kmer


class TestExtractReadKmers:
    """Tests for _extract_read_kmers helper."""

    def test_basic_extraction(self):
        """Extracts canonical k-mers with correct positions."""
        seq = "ACGTACGT"
        canon_at_pos, unique_cands = _extract_read_kmers(seq, 5)
        assert len(canon_at_pos) == 4  # positions 0..3
        assert all(isinstance(k, str) for k in canon_at_pos.values())
        assert len(unique_cands) > 0

    def test_short_sequence(self):
        """Sequence shorter than kmer_size returns empty results."""
        canon_at_pos, unique_cands = _extract_read_kmers("ACG", 5)
        assert canon_at_pos == {}
        assert unique_cands == []

    def test_skips_n_containing_kmers(self):
        """K-mers containing N should be skipped."""
        seq = "ACNGTACGT"
        canon_at_pos, unique_cands = _extract_read_kmers(seq, 5)
        # positions 0-4 span the N, so many k-mers are skipped
        for pos, kmer in canon_at_pos.items():
            assert "N" not in kmer

    def test_deduplicates_candidates(self):
        """unique_candidates should contain no duplicates."""
        seq = "AAAAAAAAAA"  # many identical k-mers
        canon_at_pos, unique_cands = _extract_read_kmers(seq, 5)
        assert len(unique_cands) == len(set(unique_cands))

    def test_canonicalization(self):
        """Extracted k-mers are canonical."""
        seq = "ACGTACGT"
        canon_at_pos, unique_cands = _extract_read_kmers(seq, 5)
        for kmer in unique_cands:
            assert kmer == canonicalize(kmer)

    def test_consistency_with_manual(self):
        """Results match manual canonicalization."""
        seq = "ACGTAC"
        kmer_size = 4
        canon_at_pos, unique_cands = _extract_read_kmers(seq, kmer_size)
        expected_pos = {}
        for i in range(len(seq) - kmer_size + 1):
            kmer = seq[i:i + kmer_size]
            expected_pos[i] = canonicalize(kmer)
        assert canon_at_pos == expected_pos


import os
import subprocess
import pytest

from kmer_denovo_filter.kmer_utils import JellyfishKmerQuery


@pytest.fixture
def jf_index(tmp_path):
    """Create a small jellyfish index for testing."""
    fa = tmp_path / "test.fa"
    fa.write_text(">seq\nACGTACGTACGTACGTACGT\n")
    jf = str(tmp_path / "test.jf")
    subprocess.run(
        ["jellyfish", "count", "-m", "5", "-s", "10M", "-t", "1",
         "-C", str(fa), "-o", jf],
        check=True,
    )
    return jf


class TestJellyfishKmerQueryCache:
    """Tests for JellyfishKmerQuery caching behaviour."""

    def test_query_batch_returns_hits(self, jf_index):
        """query_batch returns present k-mers."""
        q = JellyfishKmerQuery(jf_index)
        # ACGTA is in the index (from "ACGTACGTACGTACGTACGT")
        hits = q.query_batch([canonicalize("ACGTA")])
        assert len(hits) > 0
        q.close()

    def test_query_batch_no_hits(self, jf_index):
        """query_batch returns empty set for absent k-mers."""
        q = JellyfishKmerQuery(jf_index)
        hits = q.query_batch([canonicalize("TTTTT")])
        assert len(hits) == 0
        q.close()

    def test_cache_avoids_subprocess(self, jf_index):
        """Second call to query_batch should use cache (no subprocess)."""
        q = JellyfishKmerQuery(jf_index)
        kmer = canonicalize("ACGTA")
        hits1 = q.query_batch([kmer])
        # Second call — same k-mer should come from cache
        hits2 = q.query_batch([kmer])
        assert hits1 == hits2
        q.close()

    def test_cache_populated(self, jf_index):
        """Cache is populated after query_batch."""
        q = JellyfishKmerQuery(jf_index)
        kmer = canonicalize("ACGTA")
        q.query_batch([kmer])
        assert kmer in q._cache
        q.close()

    def test_mixed_cached_uncached(self, jf_index):
        """Querying mix of cached and uncached k-mers works."""
        q = JellyfishKmerQuery(jf_index)
        kmer1 = canonicalize("ACGTA")
        q.query_batch([kmer1])  # cache kmer1

        kmer2 = canonicalize("TTTTT")
        # kmer1 is cached, kmer2 is not
        hits = q.query_batch([kmer1, kmer2])
        assert kmer1 in hits
        assert kmer2 not in hits
        q.close()

    def test_scan_read_basic(self, jf_index):
        """scan_read finds k-mers present in index."""
        q = JellyfishKmerQuery(jf_index)
        unique, indices = q.scan_read("ACGTACGTACGTACGTACGT", 5)
        assert len(unique) > 0
        assert len(indices) > 0
        q.close()

    def test_scan_read_no_match(self, jf_index):
        """scan_read returns empty for non-matching sequence."""
        q = JellyfishKmerQuery(jf_index)
        unique, indices = q.scan_read("TTTTTTTTTTTTTTT", 5)
        assert len(unique) == 0
        assert len(indices) == 0
        q.close()

    def test_scan_read_caches_results(self, jf_index):
        """scan_read populates the cache for future lookups."""
        q = JellyfishKmerQuery(jf_index)
        q.scan_read("ACGTACGTACGTACGTACGT", 5)
        # Cache should contain k-mers from the read
        assert len(q._cache) > 0
        q.close()

    def test_batch_then_scan_uses_cache(self, jf_index):
        """Pre-populating cache via query_batch speeds up scan_read."""
        q = JellyfishKmerQuery(jf_index)
        seq = "ACGTACGTACGTACGTACGT"
        # Pre-populate cache
        _, unique_cands = _extract_read_kmers(seq, 5)
        q.query_batch(unique_cands)
        cache_size_after_batch = len(q._cache)

        # scan_read should use cache, not add new entries
        unique, indices = q.scan_read(seq, 5)
        assert len(q._cache) == cache_size_after_batch
        assert len(unique) > 0
        q.close()

    def test_empty_batch(self, jf_index):
        """Empty batch returns empty set."""
        q = JellyfishKmerQuery(jf_index)
        hits = q.query_batch([])
        assert hits == set()
        q.close()

    def test_close_clears_cache(self, jf_index):
        """close() clears the internal cache."""
        q = JellyfishKmerQuery(jf_index)
        q.query_batch([canonicalize("ACGTA")])
        assert len(q._cache) > 0
        q.close()
        assert len(q._cache) == 0
