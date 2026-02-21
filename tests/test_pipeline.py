"""Integration test for the full pipeline using synthetic data."""

import json
import os
import tempfile

import pysam
import pytest

from kmer_denovo_filter.cli import parse_args
from kmer_denovo_filter.pipeline import run_pipeline

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
        vcf_out.close()


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
