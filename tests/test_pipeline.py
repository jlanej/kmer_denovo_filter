"""Integration test for the full pipeline using synthetic data."""

import json
import logging
import os
import tempfile

import pysam
import pytest

from kmer_denovo_filter.cli import parse_args
from kmer_denovo_filter.pipeline import (
    _format_elapsed,
    _format_file_size,
    _validate_inputs,
    run_discovery_pipeline,
    run_pipeline,
)

GIAB_DIR = os.path.join(os.path.dirname(__file__), "data", "giab")
GIAB_DATA_EXISTS = os.path.isfile(os.path.join(GIAB_DIR, "HG002_child.bam"))


def _create_ref_fasta(path, chrom="chr1", length=200):
    """Create a small reference FASTA with non-repetitive sequence.

    Uses MD5 solely to produce a deterministic pseudo-random base sequence
    for test data (not for any cryptographic purpose).
    """
    import hashlib
    seq = []
    bases = "ACGT"
    for i in range(length):
        h = hashlib.md5(str(i).encode()).hexdigest()
        seq.append(bases[int(h, 16) % 4])
    seq = "".join(seq)
    with open(path, "w") as fh:
        fh.write(f">{chrom}\n{seq}\n")
    pysam.faidx(path)
    return seq


def _create_bam(path, ref_fasta, chrom, reads):
    """Create a BAM file from a list of (name, pos, seq, quals) tuples.

    ``pos`` is 0-based.  An optional fifth element may supply a CIGAR
    tuple list; otherwise a simple all-M alignment is used.
    """
    header = pysam.AlignmentHeader.from_references(
        [chrom], [300],
    )
    with pysam.AlignmentFile(path, "wb", header=header) as bam:
        for entry in reads:
            name, pos, seq, *rest = entry
            quals = rest[0] if rest else None
            cigar = rest[1] if len(rest) > 1 else [(0, len(seq))]
            seg = pysam.AlignedSegment()
            seg.query_name = name
            seg.query_sequence = seq
            seg.flag = 0
            seg.reference_id = 0
            seg.reference_start = pos
            seg.mapping_quality = 60
            seg.cigar = cigar
            if quals is not None:
                seg.query_qualities = pysam.qualitystring_to_array(quals)
            else:
                seg.query_qualities = pysam.qualitystring_to_array(
                    "I" * len(seq)
                )
            bam.write(seg)
    pysam.sort("-o", path, path)
    pysam.index(path)


def _create_vcf(path, chrom, variants, sample="HG002"):
    """Create a VCF with given variants.

    ``variants`` is a list of (pos_1based, ref, alt) tuples.
    """
    header = pysam.VariantHeader()
    header.add_sample(sample)
    header.add_line(f"##contig=<ID={chrom},length=300>")
    with pysam.VariantFile(path, "w", header=header) as vcf:
        for pos, ref, alt in variants:
            rec = vcf.new_record(
                contig=chrom, start=pos - 1, stop=pos - 1 + len(ref),
                alleles=(ref, alt),
            )
            vcf.write(rec)


class TestPipelineIntegration:
    """End-to-end pipeline tests with synthetic data."""

    @pytest.fixture()
    def tmpdir(self, tmp_path):
        return str(tmp_path)

    def test_denovo_variant_detected(self, tmpdir):
        """A variant supported only in the child should yield DKU > 0."""
        chrom = "chr1"
        ref_fa = os.path.join(tmpdir, "ref.fa")
        ref_seq = _create_ref_fasta(ref_fa, chrom, 200)

        # Variant at 1-based position 51 (0-based 50)
        var_pos_0 = 50

        # Build child read that contains a mutation at position 50
        child_seq = list(ref_seq[40:80])
        child_seq[10] = "G" if ref_seq[50] != "G" else "T"
        child_seq = "".join(child_seq)

        child_bam = os.path.join(tmpdir, "child.bam")
        _create_bam(
            child_bam, ref_fa, chrom,
            [("read1", 40, child_seq, None)],
        )

        # Parent reads match reference (no mutation)
        parent_seq = ref_seq[40:80]
        mother_bam = os.path.join(tmpdir, "mother.bam")
        _create_bam(
            mother_bam, ref_fa, chrom,
            [("mread1", 40, parent_seq, None)],
        )

        father_bam = os.path.join(tmpdir, "father.bam")
        _create_bam(
            father_bam, ref_fa, chrom,
            [("fread1", 40, parent_seq, None)],
        )

        in_vcf = os.path.join(tmpdir, "input.vcf")
        alt_base = child_seq[10]
        _create_vcf(in_vcf, chrom, [(51, ref_seq[50], alt_base)])

        out_vcf = os.path.join(tmpdir, "output.vcf.gz")
        metrics_json = os.path.join(tmpdir, "metrics.json")
        info_bam = os.path.join(tmpdir, "informative.bam")

        args = parse_args([
            "--child", child_bam,
            "--mother", mother_bam,
            "--father", father_bam,
            "--ref-fasta", ref_fa,
            "--vcf", in_vcf,
            "--output", out_vcf,
            "--metrics", metrics_json,
            "--kmer-size", "5",
            "--informative-reads", info_bam,
            "--proband-id", "HG002",
        ])
        run_pipeline(args)

        # Check output VCF is bgzipped and tabix-indexed
        assert os.path.exists(out_vcf)
        assert os.path.exists(out_vcf + ".tbi")

        # Check output VCF has DKU annotation as FORMAT (proband matched)
        vcf_out = pysam.VariantFile(out_vcf)
        records = list(vcf_out)
        assert len(records) == 1
        assert records[0].samples["HG002"]["DKU"] > 0
        assert records[0].samples["HG002"]["DKT"] is not None
        assert records[0].samples["HG002"]["DKA"] > 0
        # Proportion fields
        assert records[0].samples["HG002"]["DKU_DKT"] > 0
        assert records[0].samples["HG002"]["DKA_DKT"] > 0
        # De novo variant: child-unique k-mers are not in parents
        assert records[0].samples["HG002"]["MAX_PKC"] is not None
        assert records[0].samples["HG002"]["AVG_PKC"] is not None
        assert records[0].samples["HG002"]["MIN_PKC"] is not None
        assert records[0].samples["HG002"]["MAX_PKC_ALT"] is not None
        assert records[0].samples["HG002"]["AVG_PKC_ALT"] is not None
        assert records[0].samples["HG002"]["MIN_PKC_ALT"] is not None
        vcf_out.close()

        # Check metrics file
        with open(metrics_json) as fh:
            metrics = json.load(fh)
        assert metrics["total_variants"] == 1
        assert metrics["child_unique_kmers"] > 0

        # Check informative reads BAM
        assert os.path.exists(info_bam)
        assert os.path.exists(info_bam + ".bai")
        bam_info = pysam.AlignmentFile(info_bam)
        info_reads = list(bam_info)
        assert len(info_reads) >= 1
        # Each read should have a DV tag
        for read in info_reads:
            assert read.has_tag("DV")
            assert read.query_name == "read1"
        bam_info.close()

    def test_inherited_variant_no_unique(self, tmpdir):
        """When the variant is also in a parent, DKU should be 0."""
        chrom = "chr1"
        ref_fa = os.path.join(tmpdir, "ref.fa")
        ref_seq = _create_ref_fasta(ref_fa, chrom, 200)

        # Same mutation in child and mother
        mut_seq = list(ref_seq[40:80])
        mut_seq[10] = "G" if ref_seq[50] != "G" else "T"
        mut_seq = "".join(mut_seq)

        child_bam = os.path.join(tmpdir, "child.bam")
        _create_bam(
            child_bam, ref_fa, chrom,
            [("read1", 40, mut_seq, None)],
        )

        mother_bam = os.path.join(tmpdir, "mother.bam")
        _create_bam(
            mother_bam, ref_fa, chrom,
            [("mread1", 40, mut_seq, None)],
        )

        father_bam = os.path.join(tmpdir, "father.bam")
        _create_bam(
            father_bam, ref_fa, chrom,
            [("fread1", 40, ref_seq[40:80], None)],
        )

        in_vcf = os.path.join(tmpdir, "input.vcf")
        alt_base = mut_seq[10]
        _create_vcf(in_vcf, chrom, [(51, ref_seq[50], alt_base)])

        out_vcf = os.path.join(tmpdir, "output.vcf.gz")
        info_bam = os.path.join(tmpdir, "informative.bam")

        args = parse_args([
            "--child", child_bam,
            "--mother", mother_bam,
            "--father", father_bam,
            "--ref-fasta", ref_fa,
            "--vcf", in_vcf,
            "--output", out_vcf,
            "--kmer-size", "5",
            "--informative-reads", info_bam,
            "--proband-id", "HG002",
        ])
        run_pipeline(args)

        vcf_out = pysam.VariantFile(out_vcf)
        records = list(vcf_out)
        assert len(records) == 1
        assert records[0].samples["HG002"]["DKU"] == 0
        assert records[0].samples["HG002"]["DKA"] == 0
        # Proportion fields should be 0.0 for inherited variants
        assert records[0].samples["HG002"]["DKU_DKT"] == pytest.approx(0.0)
        assert records[0].samples["HG002"]["DKA_DKT"] == pytest.approx(0.0)
        # Inherited variant: child k-mers are shared with parent, so max_pkc >= 1
        assert records[0].samples["HG002"]["MAX_PKC"] >= 1
        assert records[0].samples["HG002"]["MIN_PKC"] >= 1
        vcf_out.close()

        # No informative reads for inherited variant
        bam_info = pysam.AlignmentFile(info_bam)
        info_reads = list(bam_info)
        assert len(info_reads) == 0
        bam_info.close()

    def test_empty_vcf(self, tmpdir):
        """Pipeline should handle an empty VCF without error."""
        chrom = "chr1"
        ref_fa = os.path.join(tmpdir, "ref.fa")
        _create_ref_fasta(ref_fa, chrom, 200)

        for label in ("child", "mother", "father"):
            bam = os.path.join(tmpdir, f"{label}.bam")
            _create_bam(bam, ref_fa, chrom, [])

        in_vcf = os.path.join(tmpdir, "input.vcf")
        _create_vcf(in_vcf, chrom, [])

        out_vcf = os.path.join(tmpdir, "output.vcf.gz")

        args = parse_args([
            "--child", os.path.join(tmpdir, "child.bam"),
            "--mother", os.path.join(tmpdir, "mother.bam"),
            "--father", os.path.join(tmpdir, "father.bam"),
            "--ref-fasta", ref_fa,
            "--vcf", in_vcf,
            "--output", out_vcf,
        ])
        run_pipeline(args)

        vcf_out = pysam.VariantFile(out_vcf)
        records = list(vcf_out)
        assert len(records) == 0
        vcf_out.close()

    def test_deletion_variant_detected(self, tmpdir):
        """A de novo deletion should yield DKU > 0."""
        chrom = "chr1"
        ref_fa = os.path.join(tmpdir, "ref.fa")
        ref_seq = _create_ref_fasta(ref_fa, chrom, 200)

        # Deletion: remove 2 bases at 0-based positions 51-52
        # VCF anchor at 1-based pos 51 (0-based 50)
        anchor_0 = 50
        del_len = 2

        # Child read: ref_seq[40:51] + ref_seq[53:80] with a 2-base deletion
        child_seq = ref_seq[40:anchor_0 + 1] + ref_seq[anchor_0 + 1 + del_len:80]
        left_match = anchor_0 + 1 - 40   # 11
        right_match = 80 - (anchor_0 + 1 + del_len)  # 27
        child_cigar = [(0, left_match), (2, del_len), (0, right_match)]

        child_bam = os.path.join(tmpdir, "child.bam")
        _create_bam(
            child_bam, ref_fa, chrom,
            [("read1", 40, child_seq, None, child_cigar)],
        )

        # Parents match reference (no deletion)
        parent_seq = ref_seq[40:80]
        mother_bam = os.path.join(tmpdir, "mother.bam")
        _create_bam(
            mother_bam, ref_fa, chrom,
            [("mread1", 40, parent_seq, None)],
        )
        father_bam = os.path.join(tmpdir, "father.bam")
        _create_bam(
            father_bam, ref_fa, chrom,
            [("fread1", 40, parent_seq, None)],
        )

        ref_allele = ref_seq[anchor_0:anchor_0 + 1 + del_len]
        alt_allele = ref_seq[anchor_0]
        in_vcf = os.path.join(tmpdir, "input.vcf")
        _create_vcf(in_vcf, chrom, [(anchor_0 + 1, ref_allele, alt_allele)])

        out_vcf = os.path.join(tmpdir, "output.vcf.gz")
        args = parse_args([
            "--child", child_bam,
            "--mother", mother_bam,
            "--father", father_bam,
            "--ref-fasta", ref_fa,
            "--vcf", in_vcf,
            "--output", out_vcf,
            "--kmer-size", "5",
            "--proband-id", "HG002",
        ])
        run_pipeline(args)

        vcf_out = pysam.VariantFile(out_vcf)
        records = list(vcf_out)
        assert len(records) == 1
        assert records[0].samples["HG002"]["DKU"] > 0, \
            "Deletion should produce DKU > 0"
        vcf_out.close()

    def test_insertion_variant_detected(self, tmpdir):
        """A de novo insertion should yield DKU > 0."""
        chrom = "chr1"
        ref_fa = os.path.join(tmpdir, "ref.fa")
        ref_seq = _create_ref_fasta(ref_fa, chrom, 200)

        # Insertion: insert 3 bases after 0-based position 50
        anchor_0 = 50
        ins_bases = "GCA"

        # Child read: ref_seq[40:51] + "GCA" + ref_seq[51:80]
        child_seq = ref_seq[40:anchor_0 + 1] + ins_bases + ref_seq[anchor_0 + 1:80]
        left_match = anchor_0 + 1 - 40   # 11
        right_match = 80 - (anchor_0 + 1)  # 29
        child_cigar = [
            (0, left_match), (1, len(ins_bases)), (0, right_match),
        ]

        child_bam = os.path.join(tmpdir, "child.bam")
        _create_bam(
            child_bam, ref_fa, chrom,
            [("read1", 40, child_seq, None, child_cigar)],
        )

        # Parents match reference (no insertion)
        parent_seq = ref_seq[40:80]
        mother_bam = os.path.join(tmpdir, "mother.bam")
        _create_bam(
            mother_bam, ref_fa, chrom,
            [("mread1", 40, parent_seq, None)],
        )
        father_bam = os.path.join(tmpdir, "father.bam")
        _create_bam(
            father_bam, ref_fa, chrom,
            [("fread1", 40, parent_seq, None)],
        )

        ref_allele = ref_seq[anchor_0]
        alt_allele = ref_seq[anchor_0] + ins_bases
        in_vcf = os.path.join(tmpdir, "input.vcf")
        _create_vcf(in_vcf, chrom, [(anchor_0 + 1, ref_allele, alt_allele)])

        out_vcf = os.path.join(tmpdir, "output.vcf.gz")
        args = parse_args([
            "--child", child_bam,
            "--mother", mother_bam,
            "--father", father_bam,
            "--ref-fasta", ref_fa,
            "--vcf", in_vcf,
            "--output", out_vcf,
            "--kmer-size", "5",
            "--proband-id", "HG002",
        ])
        run_pipeline(args)

        vcf_out = pysam.VariantFile(out_vcf)
        records = list(vcf_out)
        assert len(records) == 1
        assert records[0].samples["HG002"]["DKU"] > 0, \
            "Insertion should produce DKU > 0"
        vcf_out.close()

    def test_summary_includes_pkc_fields(self, tmpdir):
        """Summary output should include MAX_PKC and AVG_PKC columns."""
        chrom = "chr1"
        ref_fa = os.path.join(tmpdir, "ref.fa")
        ref_seq = _create_ref_fasta(ref_fa, chrom, 200)

        child_seq = list(ref_seq[40:80])
        child_seq[10] = "G" if ref_seq[50] != "G" else "T"
        child_seq = "".join(child_seq)

        child_bam = os.path.join(tmpdir, "child.bam")
        _create_bam(
            child_bam, ref_fa, chrom,
            [("read1", 40, child_seq, None)],
        )

        parent_seq = ref_seq[40:80]
        mother_bam = os.path.join(tmpdir, "mother.bam")
        _create_bam(mother_bam, ref_fa, chrom, [("mread1", 40, parent_seq, None)])
        father_bam = os.path.join(tmpdir, "father.bam")
        _create_bam(father_bam, ref_fa, chrom, [("fread1", 40, parent_seq, None)])

        in_vcf = os.path.join(tmpdir, "input.vcf")
        alt_base = child_seq[10]
        _create_vcf(in_vcf, chrom, [(51, ref_seq[50], alt_base)])

        out_vcf = os.path.join(tmpdir, "output.vcf.gz")
        summary_txt = os.path.join(tmpdir, "summary.txt")

        args = parse_args([
            "--child", child_bam,
            "--mother", mother_bam,
            "--father", father_bam,
            "--ref-fasta", ref_fa,
            "--vcf", in_vcf,
            "--output", out_vcf,
            "--summary", summary_txt,
            "--kmer-size", "5",
            "--proband-id", "HG002",
        ])
        run_pipeline(args)

        assert os.path.exists(summary_txt)
        with open(summary_txt) as fh:
            summary = fh.read()

        # Summary should contain MAX_PKC and AVG_PKC in the header and stats
        assert "MAX_PKC" in summary
        assert "AVG_PKC" in summary
        assert "MIN_PKC" in summary
        assert "MAX_PKC_ALT" in summary
        assert "AVG_PKC_ALT" in summary
        assert "MIN_PKC_ALT" in summary
        # Summary should contain proportion fields
        assert "DKU_DKT" in summary
        assert "DKA_DKT" in summary

    def test_info_annotation_when_proband_unmatched(self, tmpdir):
        """When --proband-id does not match a VCF sample, use INFO fields."""
        chrom = "chr1"
        ref_fa = os.path.join(tmpdir, "ref.fa")
        ref_seq = _create_ref_fasta(ref_fa, chrom, 200)

        child_seq = list(ref_seq[40:80])
        child_seq[10] = "G" if ref_seq[50] != "G" else "T"
        child_seq = "".join(child_seq)

        child_bam = os.path.join(tmpdir, "child.bam")
        _create_bam(
            child_bam, ref_fa, chrom,
            [("read1", 40, child_seq, None)],
        )

        parent_seq = ref_seq[40:80]
        mother_bam = os.path.join(tmpdir, "mother.bam")
        _create_bam(mother_bam, ref_fa, chrom, [("mread1", 40, parent_seq, None)])
        father_bam = os.path.join(tmpdir, "father.bam")
        _create_bam(father_bam, ref_fa, chrom, [("fread1", 40, parent_seq, None)])

        in_vcf = os.path.join(tmpdir, "input.vcf")
        alt_base = child_seq[10]
        _create_vcf(in_vcf, chrom, [(51, ref_seq[50], alt_base)])

        out_vcf = os.path.join(tmpdir, "output.vcf.gz")

        # Use a proband ID that doesn't match any VCF sample
        args = parse_args([
            "--child", child_bam,
            "--mother", mother_bam,
            "--father", father_bam,
            "--ref-fasta", ref_fa,
            "--vcf", in_vcf,
            "--output", out_vcf,
            "--kmer-size", "5",
            "--proband-id", "NO_MATCH",
        ])
        run_pipeline(args)

        vcf_out = pysam.VariantFile(out_vcf)
        records = list(vcf_out)
        assert len(records) == 1
        assert "DKU" in records[0].info
        assert records[0].info["DKU"] > 0
        assert "DKT" in records[0].info
        assert "DKA" in records[0].info
        assert records[0].info["DKA"] > 0
        assert "DKU_DKT" in records[0].info
        assert records[0].info["DKU_DKT"] > 0
        assert "DKA_DKT" in records[0].info
        assert records[0].info["DKA_DKT"] > 0
        assert "MAX_PKC" in records[0].info
        assert "AVG_PKC" in records[0].info
        assert "MIN_PKC" in records[0].info
        assert "MAX_PKC_ALT" in records[0].info
        assert "AVG_PKC_ALT" in records[0].info
        assert "MIN_PKC_ALT" in records[0].info
        vcf_out.close()

    def test_info_annotation_when_no_proband_id(self, tmpdir):
        """When --proband-id is not supplied, use INFO fields."""
        chrom = "chr1"
        ref_fa = os.path.join(tmpdir, "ref.fa")
        ref_seq = _create_ref_fasta(ref_fa, chrom, 200)

        child_seq = list(ref_seq[40:80])
        child_seq[10] = "G" if ref_seq[50] != "G" else "T"
        child_seq = "".join(child_seq)

        child_bam = os.path.join(tmpdir, "child.bam")
        _create_bam(
            child_bam, ref_fa, chrom,
            [("read1", 40, child_seq, None)],
        )

        parent_seq = ref_seq[40:80]
        mother_bam = os.path.join(tmpdir, "mother.bam")
        _create_bam(mother_bam, ref_fa, chrom, [("mread1", 40, parent_seq, None)])
        father_bam = os.path.join(tmpdir, "father.bam")
        _create_bam(father_bam, ref_fa, chrom, [("fread1", 40, parent_seq, None)])

        in_vcf = os.path.join(tmpdir, "input.vcf")
        alt_base = child_seq[10]
        _create_vcf(in_vcf, chrom, [(51, ref_seq[50], alt_base)])

        out_vcf = os.path.join(tmpdir, "output.vcf.gz")

        # No --proband-id at all
        args = parse_args([
            "--child", child_bam,
            "--mother", mother_bam,
            "--father", father_bam,
            "--ref-fasta", ref_fa,
            "--vcf", in_vcf,
            "--output", out_vcf,
            "--kmer-size", "5",
        ])
        run_pipeline(args)

        vcf_out = pysam.VariantFile(out_vcf)
        records = list(vcf_out)
        assert len(records) == 1
        assert "DKU" in records[0].info
        assert records[0].info["DKU"] > 0
        assert "DKT" in records[0].info
        assert "DKA" in records[0].info
        assert records[0].info["DKA"] > 0
        assert "DKU_DKT" in records[0].info
        assert records[0].info["DKU_DKT"] > 0
        assert "DKA_DKT" in records[0].info
        assert records[0].info["DKA_DKT"] > 0
        assert "MAX_PKC" in records[0].info
        assert "AVG_PKC" in records[0].info
        assert "MIN_PKC" in records[0].info
        assert "MAX_PKC_ALT" in records[0].info
        assert "AVG_PKC_ALT" in records[0].info
        assert "MIN_PKC_ALT" in records[0].info
        vcf_out.close()


class TestFormatElapsed:
    """Unit tests for the _format_elapsed helper."""

    def test_seconds_only(self):
        assert _format_elapsed(5.123) == "5.1s"

    def test_minutes_and_seconds(self):
        result = _format_elapsed(125.5)
        assert result == "2m 5.5s"

    def test_hours(self):
        result = _format_elapsed(3725)
        assert result == "1h 2m 5s"

    def test_zero(self):
        assert _format_elapsed(0) == "0.0s"


class TestFormatFileSize:
    """Unit tests for the _format_file_size helper."""

    def test_bytes(self, tmp_path):
        f = tmp_path / "small.txt"
        f.write_text("hello")
        assert "B" in _format_file_size(str(f))

    def test_missing_file(self):
        assert _format_file_size("/no/such/file") == "?"

    def test_empty_file(self, tmp_path):
        f = tmp_path / "empty.txt"
        f.write_bytes(b"")
        assert _format_file_size(str(f)) == "0.0 B"


class TestValidateInputs:
    """Tests for input validation in _validate_inputs."""

    @pytest.fixture()
    def tmpdir(self, tmp_path):
        return str(tmp_path)

    def _make_args(self, tmpdir, **overrides):
        """Build a minimal args namespace with valid dummy files."""
        # Create dummy files so validation passes by default
        for name in (
            "child.bam", "child.bam.bai",
            "mother.bam", "mother.bam.bai",
            "father.bam", "father.bam.bai",
            "input.vcf",
        ):
            path = os.path.join(tmpdir, name)
            if not os.path.exists(path):
                open(path, "w").close()

        defaults = {
            "child": os.path.join(tmpdir, "child.bam"),
            "mother": os.path.join(tmpdir, "mother.bam"),
            "father": os.path.join(tmpdir, "father.bam"),
            "vcf": os.path.join(tmpdir, "input.vcf"),
            "output": os.path.join(tmpdir, "output.vcf.gz"),
            "ref_fasta": None,
            "kmer_size": 31,
            "min_baseq": 20,
            "min_mapq": 20,
            "threads": 4,
            "debug_kmers": False,
        }
        defaults.update(overrides)
        import argparse
        return argparse.Namespace(**defaults)

    def test_valid_inputs_pass(self, tmpdir):
        args = self._make_args(tmpdir)
        # Should not raise
        _validate_inputs(args)

    def test_missing_child_bam(self, tmpdir):
        args = self._make_args(tmpdir, child="/no/such/file.bam")
        with pytest.raises(SystemExit):
            _validate_inputs(args)

    def test_missing_vcf(self, tmpdir):
        args = self._make_args(tmpdir, vcf="/no/such/input.vcf")
        with pytest.raises(SystemExit):
            _validate_inputs(args)

    def test_missing_ref_fasta(self, tmpdir):
        args = self._make_args(tmpdir, ref_fasta="/no/such/ref.fa")
        with pytest.raises(SystemExit):
            _validate_inputs(args)

    def test_cram_without_ref_fasta(self, tmpdir):
        cram = os.path.join(tmpdir, "child.cram")
        open(cram, "w").close()
        open(cram + ".crai", "w").close()
        args = self._make_args(tmpdir, child=cram, ref_fasta=None)
        with pytest.raises(SystemExit):
            _validate_inputs(args)

    def test_missing_bam_index(self, tmpdir):
        bam = os.path.join(tmpdir, "noindex.bam")
        open(bam, "w").close()
        args = self._make_args(tmpdir, child=bam)
        with pytest.raises(SystemExit):
            _validate_inputs(args)

    def test_even_kmer_size(self, tmpdir):
        args = self._make_args(tmpdir, kmer_size=30)
        with pytest.raises(SystemExit):
            _validate_inputs(args)

    def test_kmer_size_too_small(self, tmpdir):
        args = self._make_args(tmpdir, kmer_size=1)
        with pytest.raises(SystemExit):
            _validate_inputs(args)

    def test_kmer_size_too_large(self, tmpdir):
        args = self._make_args(tmpdir, kmer_size=300)
        with pytest.raises(SystemExit):
            _validate_inputs(args)

    def test_negative_min_baseq(self, tmpdir):
        args = self._make_args(tmpdir, min_baseq=-1)
        with pytest.raises(SystemExit):
            _validate_inputs(args)

    def test_negative_min_mapq(self, tmpdir):
        args = self._make_args(tmpdir, min_mapq=-5)
        with pytest.raises(SystemExit):
            _validate_inputs(args)

    def test_zero_threads(self, tmpdir):
        args = self._make_args(tmpdir, threads=0)
        with pytest.raises(SystemExit):
            _validate_inputs(args)


class TestProgressLogging:
    """Tests verifying structured progress messages are logged."""

    @pytest.fixture()
    def tmpdir(self, tmp_path):
        return str(tmp_path)

    def test_pipeline_logs_step_markers(self, tmpdir, caplog):
        """The pipeline should log step markers and elapsed times."""
        chrom = "chr1"
        ref_fa = os.path.join(tmpdir, "ref.fa")
        ref_seq = _create_ref_fasta(ref_fa, chrom, 200)

        child_seq = list(ref_seq[40:80])
        child_seq[10] = "G" if ref_seq[50] != "G" else "T"
        child_seq = "".join(child_seq)

        child_bam = os.path.join(tmpdir, "child.bam")
        _create_bam(child_bam, ref_fa, chrom,
                     [("read1", 40, child_seq, None)])

        parent_seq = ref_seq[40:80]
        mother_bam = os.path.join(tmpdir, "mother.bam")
        _create_bam(mother_bam, ref_fa, chrom,
                     [("mread1", 40, parent_seq, None)])
        father_bam = os.path.join(tmpdir, "father.bam")
        _create_bam(father_bam, ref_fa, chrom,
                     [("fread1", 40, parent_seq, None)])

        in_vcf = os.path.join(tmpdir, "input.vcf")
        alt_base = child_seq[10]
        _create_vcf(in_vcf, chrom, [(51, ref_seq[50], alt_base)])

        out_vcf = os.path.join(tmpdir, "output.vcf.gz")

        args = parse_args([
            "--child", child_bam,
            "--mother", mother_bam,
            "--father", father_bam,
            "--ref-fasta", ref_fa,
            "--vcf", in_vcf,
            "--output", out_vcf,
            "--kmer-size", "5",
        ])

        with caplog.at_level(logging.INFO):
            run_pipeline(args)

        log_text = caplog.text
        # Check step markers appear
        assert "[Step 1/5]" in log_text
        assert "[Step 2/5]" in log_text
        assert "[Step 3/5]" in log_text
        assert "[Step 4/5]" in log_text
        assert "[Step 5/5]" in log_text
        # Check configuration summary
        assert "pipeline starting" in log_text
        assert "k-mer size:" in log_text
        # Check file size appears in config summary
        assert "B)" in log_text  # file size like "(123.0 B)"
        # Check child k-mer extraction progress
        assert "reads scanned" in log_text
        assert "k-mers collected" in log_text
        # Check parent scan progress
        assert "Mother done" in log_text or "Mother scan" in log_text
        assert "Father done" in log_text or "Father scan" in log_text
        assert "child k-mers found" in log_text
        # Check annotation progress with running tallies
        assert "de novo so far" in log_text
        assert "total reads" in log_text
        assert "var/s" in log_text
        assert "remaining" in log_text
        # Check unique k-mer percentage
        assert "% unique" in log_text
        # Check completion message
        assert "Pipeline finished successfully" in log_text


@pytest.mark.skipif(
    not GIAB_DATA_EXISTS,
    reason="GIAB test data not available",
)
class TestGIABIntegration:
    """Integration tests using real GIAB HG002 trio data."""

    @pytest.fixture()
    def tmpdir(self, tmp_path):
        return str(tmp_path)

    def test_giab_denovo_pipeline(self, tmpdir):
        """Run the full pipeline on GIAB child-private SNVs without a ref."""
        out_vcf = os.path.join(tmpdir, "annotated.vcf.gz")
        metrics_json = os.path.join(tmpdir, "metrics.json")
        summary_txt = os.path.join(tmpdir, "summary.txt")

        args = parse_args([
            "--child", os.path.join(GIAB_DIR, "HG002_child.bam"),
            "--mother", os.path.join(GIAB_DIR, "HG004_mother.bam"),
            "--father", os.path.join(GIAB_DIR, "HG003_father.bam"),
            "--vcf", os.path.join(GIAB_DIR, "candidates.vcf.gz"),
            "--output", out_vcf,
            "--metrics", metrics_json,
            "--summary", summary_txt,
            "--kmer-size", "31",
            "--proband-id", "HG002",
        ])
        run_pipeline(args)

        # Annotated VCF should exist with all variants annotated
        assert os.path.exists(out_vcf)
        assert os.path.exists(out_vcf + ".tbi")
        vcf_out = pysam.VariantFile(out_vcf)
        records = list(vcf_out)
        assert len(records) == 20
        for rec in records:
            assert rec.samples["HG002"]["DKU"] is not None
            assert rec.samples["HG002"]["DKT"] is not None
            assert rec.samples["HG002"]["DKA"] is not None
        vcf_out.close()

        # At least some variants should be flagged as likely de novo
        dnm_count = sum(
            1 for r in records if r.samples["HG002"]["DKU"] > 0
        )
        assert dnm_count > 0, "Expected at least some likely DNMs"

        # Metrics file should exist
        assert os.path.exists(metrics_json)
        with open(metrics_json) as fh:
            metrics = json.load(fh)
        assert metrics["total_variants"] == 20
        assert metrics["variants_with_unique_reads"] > 0

        # Summary file should exist and contain key sections
        assert os.path.exists(summary_txt)
        with open(summary_txt) as fh:
            summary = fh.read()
        assert "De Novo Variant Summary" in summary
        assert "Likely de novo" in summary
        assert "DE_NOVO" in summary


class TestDiscoveryPipeline:
    """End-to-end tests for VCF-free discovery mode."""

    @pytest.fixture()
    def tmpdir(self, tmp_path):
        return str(tmp_path)

    def test_discovery_denovo_detected(self, tmpdir):
        """A de novo variant should produce a BED region and informative reads."""
        chrom = "chr1"
        ref_fa = os.path.join(tmpdir, "ref.fa")
        ref_seq = _create_ref_fasta(ref_fa, chrom, 200)

        # Child read with a mutation at position 50
        child_seq = list(ref_seq[30:90])
        child_seq[20] = "G" if ref_seq[50] != "G" else "T"
        child_seq = "".join(child_seq)

        child_bam = os.path.join(tmpdir, "child.bam")
        # Create multiple copies to exceed min_child_count
        _create_bam(
            child_bam, ref_fa, chrom,
            [
                ("read1", 30, child_seq, None),
                ("read2", 30, child_seq, None),
                ("read3", 30, child_seq, None),
                ("read4", 30, child_seq, None),
            ],
        )

        # Parent reads match reference (no mutation)
        parent_seq = ref_seq[30:90]
        mother_bam = os.path.join(tmpdir, "mother.bam")
        _create_bam(
            mother_bam, ref_fa, chrom,
            [
                ("mread1", 30, parent_seq, None),
                ("mread2", 30, parent_seq, None),
                ("mread3", 30, parent_seq, None),
            ],
        )
        father_bam = os.path.join(tmpdir, "father.bam")
        _create_bam(
            father_bam, ref_fa, chrom,
            [
                ("fread1", 30, parent_seq, None),
                ("fread2", 30, parent_seq, None),
                ("fread3", 30, parent_seq, None),
            ],
        )

        out_prefix = os.path.join(tmpdir, "disc_out")

        args = parse_args([
            "--child", child_bam,
            "--mother", mother_bam,
            "--father", father_bam,
            "--ref-fasta", ref_fa,
            "--out-prefix", out_prefix,
            "--min-child-count", "3",
            "--kmer-size", "5",
        ])
        run_discovery_pipeline(args)

        # Check BED file with per-region read/k-mer counts
        bed_path = f"{out_prefix}.bed"
        assert os.path.exists(bed_path)
        with open(bed_path) as fh:
            bed_lines = [l.strip() for l in fh if l.strip()]
        assert len(bed_lines) >= 1, "Expected at least one candidate region"
        # BED lines should have 5 columns: chrom, start, end, reads, kmers
        parts = bed_lines[0].split("\t")
        assert len(parts) == 5
        assert int(parts[3]) >= 1  # at least 1 read
        assert int(parts[4]) >= 1  # at least 1 k-mer

        # Check informative reads BAM
        info_bam = f"{out_prefix}.informative.bam"
        assert os.path.exists(info_bam)
        assert os.path.exists(info_bam + ".bai")
        bam_info = pysam.AlignmentFile(info_bam)
        info_reads = list(bam_info)
        assert len(info_reads) >= 1
        for read in info_reads:
            assert read.has_tag("dk")
            assert read.get_tag("dk") == 1
        bam_info.close()

        # Check metrics JSON with per-region detail
        metrics_path = f"{out_prefix}.metrics.json"
        assert os.path.exists(metrics_path)
        with open(metrics_path) as fh:
            metrics = json.load(fh)
        assert metrics["mode"] == "discovery"
        assert metrics["proband_unique_kmers"] > 0
        assert metrics["informative_reads"] > 0
        assert metrics["candidate_regions"] >= 1
        # Per-region detail should be present
        assert "regions" in metrics
        assert len(metrics["regions"]) >= 1
        for region in metrics["regions"]:
            assert "chrom" in region
            assert "start" in region
            assert "end" in region
            assert "size" in region
            assert region["reads"] >= 1
            assert region["unique_kmers"] >= 1

        # Check summary text file
        summary_path = f"{out_prefix}.summary.txt"
        assert os.path.exists(summary_path)
        with open(summary_path) as fh:
            summary = fh.read()
        assert "Discovery Mode Summary" in summary
        assert "K-mer Filtering" in summary
        assert "Proband-unique k-mers" in summary
        assert "Candidate regions" in summary
        assert "Per-Region Results" in summary

    def test_discovery_inherited_no_regions(self, tmpdir):
        """When child shares k-mers with a parent, no regions should appear."""
        chrom = "chr1"
        ref_fa = os.path.join(tmpdir, "ref.fa")
        ref_seq = _create_ref_fasta(ref_fa, chrom, 200)

        # Same mutation in child and mother
        mut_seq = list(ref_seq[30:90])
        mut_seq[20] = "G" if ref_seq[50] != "G" else "T"
        mut_seq = "".join(mut_seq)

        child_bam = os.path.join(tmpdir, "child.bam")
        _create_bam(
            child_bam, ref_fa, chrom,
            [
                ("read1", 30, mut_seq, None),
                ("read2", 30, mut_seq, None),
                ("read3", 30, mut_seq, None),
                ("read4", 30, mut_seq, None),
            ],
        )

        mother_bam = os.path.join(tmpdir, "mother.bam")
        _create_bam(
            mother_bam, ref_fa, chrom,
            [
                ("mread1", 30, mut_seq, None),
                ("mread2", 30, mut_seq, None),
                ("mread3", 30, mut_seq, None),
            ],
        )

        father_bam = os.path.join(tmpdir, "father.bam")
        _create_bam(
            father_bam, ref_fa, chrom,
            [
                ("fread1", 30, ref_seq[30:90], None),
                ("fread2", 30, ref_seq[30:90], None),
                ("fread3", 30, ref_seq[30:90], None),
            ],
        )

        out_prefix = os.path.join(tmpdir, "disc_inh")

        args = parse_args([
            "--child", child_bam,
            "--mother", mother_bam,
            "--father", father_bam,
            "--ref-fasta", ref_fa,
            "--out-prefix", out_prefix,
            "--min-child-count", "3",
            "--kmer-size", "5",
        ])
        run_discovery_pipeline(args)

        bed_path = f"{out_prefix}.bed"
        assert os.path.exists(bed_path)
        with open(bed_path) as fh:
            bed_lines = [l.strip() for l in fh if l.strip()]
        # No candidate regions expected for inherited variation
        assert len(bed_lines) == 0

        metrics_path = f"{out_prefix}.metrics.json"
        with open(metrics_path) as fh:
            metrics = json.load(fh)
        assert metrics["proband_unique_kmers"] == 0

        # Summary should still be written for inherited/empty case
        summary_path = f"{out_prefix}.summary.txt"
        assert os.path.exists(summary_path)
        with open(summary_path) as fh:
            summary = fh.read()
        assert "Discovery Mode Summary" in summary

    def test_discovery_with_prebuilt_ref_jf(self, tmpdir):
        """Discovery mode should accept a prebuilt --ref-jf."""
        chrom = "chr1"
        ref_fa = os.path.join(tmpdir, "ref.fa")
        ref_seq = _create_ref_fasta(ref_fa, chrom, 200)

        child_seq = list(ref_seq[30:90])
        child_seq[20] = "G" if ref_seq[50] != "G" else "T"
        child_seq = "".join(child_seq)

        child_bam = os.path.join(tmpdir, "child.bam")
        _create_bam(
            child_bam, ref_fa, chrom,
            [
                ("read1", 30, child_seq, None),
                ("read2", 30, child_seq, None),
                ("read3", 30, child_seq, None),
                ("read4", 30, child_seq, None),
            ],
        )

        parent_seq = ref_seq[30:90]
        mother_bam = os.path.join(tmpdir, "mother.bam")
        _create_bam(
            mother_bam, ref_fa, chrom,
            [
                ("mread1", 30, parent_seq, None),
                ("mread2", 30, parent_seq, None),
                ("mread3", 30, parent_seq, None),
            ],
        )
        father_bam = os.path.join(tmpdir, "father.bam")
        _create_bam(
            father_bam, ref_fa, chrom,
            [
                ("fread1", 30, parent_seq, None),
                ("fread2", 30, parent_seq, None),
                ("fread3", 30, parent_seq, None),
            ],
        )

        # Pre-build the reference jellyfish index
        import subprocess
        ref_jf = os.path.join(tmpdir, "ref.k5.jf")
        subprocess.run([
            "jellyfish", "count", "-m", "5", "-s", "10M",
            "-t", "1", "-C", ref_fa, "-o", ref_jf,
        ], check=True)
        assert os.path.exists(ref_jf)

        out_prefix = os.path.join(tmpdir, "disc_prebuilt")

        args = parse_args([
            "--child", child_bam,
            "--mother", mother_bam,
            "--father", father_bam,
            "--ref-fasta", ref_fa,
            "--ref-jf", ref_jf,
            "--out-prefix", out_prefix,
            "--min-child-count", "3",
            "--kmer-size", "5",
        ])
        run_discovery_pipeline(args)

        bed_path = f"{out_prefix}.bed"
        assert os.path.exists(bed_path)
        metrics_path = f"{out_prefix}.metrics.json"
        assert os.path.exists(metrics_path)
        summary_path = f"{out_prefix}.summary.txt"
        assert os.path.exists(summary_path)

    def test_discovery_empty_child(self, tmpdir):
        """Discovery should handle an empty child BAM without error."""
        chrom = "chr1"
        ref_fa = os.path.join(tmpdir, "ref.fa")
        _create_ref_fasta(ref_fa, chrom, 200)

        for label in ("child", "mother", "father"):
            bam = os.path.join(tmpdir, f"{label}.bam")
            _create_bam(bam, ref_fa, chrom, [])

        out_prefix = os.path.join(tmpdir, "disc_empty")

        args = parse_args([
            "--child", os.path.join(tmpdir, "child.bam"),
            "--mother", os.path.join(tmpdir, "mother.bam"),
            "--father", os.path.join(tmpdir, "father.bam"),
            "--ref-fasta", ref_fa,
            "--out-prefix", out_prefix,
            "--min-child-count", "3",
            "--kmer-size", "5",
        ])
        run_discovery_pipeline(args)

        bed_path = f"{out_prefix}.bed"
        assert os.path.exists(bed_path)
        with open(bed_path) as fh:
            assert fh.read().strip() == ""

        metrics_path = f"{out_prefix}.metrics.json"
        with open(metrics_path) as fh:
            metrics = json.load(fh)
        assert metrics["child_candidate_kmers"] == 0

        # Summary should be written even for empty case
        summary_path = f"{out_prefix}.summary.txt"
        assert os.path.exists(summary_path)
        with open(summary_path) as fh:
            summary = fh.read()
        assert "Discovery Mode Summary" in summary
        assert "Candidate regions:" in summary


class TestDiscoveryValidation:
    """Tests for input validation in discovery mode."""

    @pytest.fixture()
    def tmpdir(self, tmp_path):
        return str(tmp_path)

    def _make_discovery_args(self, tmpdir, **overrides):
        """Build a minimal args namespace for discovery mode."""
        for name in (
            "child.bam", "child.bam.bai",
            "mother.bam", "mother.bam.bai",
            "father.bam", "father.bam.bai",
            "ref.fa",
        ):
            path = os.path.join(tmpdir, name)
            if not os.path.exists(path):
                open(path, "w").close()

        defaults = {
            "child": os.path.join(tmpdir, "child.bam"),
            "mother": os.path.join(tmpdir, "mother.bam"),
            "father": os.path.join(tmpdir, "father.bam"),
            "vcf": None,
            "output": None,
            "ref_fasta": os.path.join(tmpdir, "ref.fa"),
            "ref_jf": None,
            "out_prefix": os.path.join(tmpdir, "output"),
            "min_child_count": 3,
            "kmer_size": 31,
            "min_baseq": 20,
            "min_mapq": 20,
            "threads": 4,
            "debug_kmers": False,
        }
        defaults.update(overrides)
        import argparse
        return argparse.Namespace(**defaults)

    def test_valid_discovery_inputs(self, tmpdir):
        args = self._make_discovery_args(tmpdir)
        _validate_inputs(args)

    def test_discovery_missing_ref_fasta(self, tmpdir):
        args = self._make_discovery_args(
            tmpdir, ref_fasta=None, ref_jf=None,
        )
        with pytest.raises(SystemExit):
            _validate_inputs(args)

    def test_discovery_ref_jf_not_found(self, tmpdir):
        args = self._make_discovery_args(
            tmpdir, ref_jf="/no/such/ref.jf",
        )
        with pytest.raises(SystemExit):
            _validate_inputs(args)

    def test_discovery_min_child_count_zero(self, tmpdir):
        args = self._make_discovery_args(tmpdir, min_child_count=0)
        with pytest.raises(SystemExit):
            _validate_inputs(args)

    def test_discovery_with_ref_jf_only(self, tmpdir):
        """When --ref-jf is provided, --ref-fasta can be None."""
        ref_jf = os.path.join(tmpdir, "ref.k31.jf")
        open(ref_jf, "w").close()
        args = self._make_discovery_args(
            tmpdir, ref_fasta=None, ref_jf=ref_jf,
        )
        _validate_inputs(args)
