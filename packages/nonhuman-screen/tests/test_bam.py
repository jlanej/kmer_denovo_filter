"""Tests for the BAM/allele layer (classify_variants_alt_reads etc.).

Uses a tiny synthetic BAM and stubs kraken2 by replacing
``Kraken2Runner.classify_sequences``, so no binary or database is needed.
"""

from unittest import mock

import pytest

pysam = pytest.importorskip("pysam")

from nonhuman_screen.bam import (  # noqa: E402
    classify_reads_from_bam,
    classify_variant_alt_reads,
    classify_variants_alt_reads,
    reads_supporting_alt,
)
from nonhuman_screen.engine import Kraken2Runner  # noqa: E402


def _make_bam(path, reads):
    """Write+index a coordinate-sorted BAM. *reads* = list of (name, start, seq)."""
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": "chr1", "LN": 1000}],
    }
    with pysam.AlignmentFile(str(path), "wb", header=header) as bam:
        for name, start, seq in sorted(reads, key=lambda r: r[1]):
            a = pysam.AlignedSegment(bam.header)
            a.query_name = name
            a.query_sequence = seq
            a.flag = 0
            a.reference_id = 0
            a.reference_start = start
            a.mapping_quality = 60
            a.cigartuples = [(0, len(seq))]  # all match/mismatch
            a.query_qualities = pysam.qualitystring_to_array("I" * len(seq))
            bam.write(a)
    pysam.index(str(path))


@pytest.fixture
def snp_bam(tmp_path):
    # SNP at 0-based ref pos 105, ref 'A', alt 'T'. Reads are 10bp @100..109;
    # base at offset 5 == ref pos 105.
    reads = [
        ("alt1", 100, "AAAAA" "T" "AAAA"),   # supports T
        ("alt2", 100, "AAAAA" "T" "AAAA"),   # supports T
        ("ref1", 100, "AAAAA" "A" "AAAA"),   # ref allele, no support
    ]
    bam = tmp_path / "snp.bam"
    _make_bam(bam, reads)
    return str(bam)


def _fake_classify(self, sequences, tmpdir=None):
    """Stub: classify 'alt1' as non-human (bacterial), everything else human."""
    r = Kraken2Runner.Result()
    names = list(sequences)
    r.total = len(names)
    r.classified = len(names)
    for n in names:
        if n == "alt1":
            r.nonhuman_read_names.add(n)
            r.nonhuman_count += 1
            r.bacterial_read_names.add(n)
            r.bacterial_count += 1
        else:
            r.human_lineage_read_names.add(n)
            r.human_lineage_count += 1
    return r


class TestReadsSupportingAlt:
    def test_only_alt_reads_returned(self, snp_bam):
        names = reads_supporting_alt(snp_bam, "chr1", 105, "A", "T")
        assert names == {"alt1", "alt2"}


class TestClassifyVariantAltReads:
    def test_single_variant_nhf(self, snp_bam):
        with mock.patch.object(
            Kraken2Runner, "classify_sequences",
            autospec=True, side_effect=_fake_classify,
        ):
            v = classify_variant_alt_reads(snp_bam, "db", "chr1", 105, "A", "T")
        assert v.supporting_reads == 2               # alt1 + alt2, not ref1
        assert v.nonhuman_fraction == 0.5            # alt1 non-human of {alt1, alt2}
        assert v.fractions.bacterial == 0.5
        assert v.variant_key == "chr1:105:A:T"

    def test_no_supporting_reads(self, snp_bam):
        # alt 'G' is carried by no read -> empty support, zero fraction, no kraken2
        with mock.patch.object(
            Kraken2Runner, "classify_sequences",
            autospec=True, side_effect=_fake_classify,
        ) as m:
            v = classify_variant_alt_reads(snp_bam, "db", "chr1", 105, "A", "G")
        assert v.supporting_reads == 0
        assert v.nonhuman_fraction == 0.0
        assert m.call_count == 0  # empty union -> engine not invoked


class TestBatching:
    def test_single_kraken2_call_across_variants(self, snp_bam):
        variants = [
            ("chr1", 105, "A", "T"),
            ("chr1", 105, "A", "T"),  # same locus/allele, shared reads
        ]
        with mock.patch.object(
            Kraken2Runner, "classify_sequences",
            autospec=True, side_effect=_fake_classify,
        ) as m:
            results = classify_variants_alt_reads(snp_bam, "db", variants)
        assert len(results) == 2
        assert all(v.nonhuman_fraction == 0.5 for v in results)
        assert m.call_count == 1  # batched: one classify call for all variants


class TestClassifyReadsFromBam:
    def test_loci_targeted(self, snp_bam):
        loci = {("chr1", 105): {"alt1", "ref1"}}
        with mock.patch.object(
            Kraken2Runner, "classify_sequences",
            autospec=True, side_effect=_fake_classify,
        ):
            res = classify_reads_from_bam(snp_bam, "db", loci=loci)
        assert res.total == 2  # alt1 + ref1 fetched at the locus
        assert res.nonhuman_read_names == {"alt1"}
