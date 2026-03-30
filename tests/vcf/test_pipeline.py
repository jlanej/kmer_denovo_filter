"""VCF-mode pipeline tests (split from tests/test_pipeline.py)."""

import json
import logging
import os
import tempfile

import pysam
import pytest

import kmer_denovo_filter.pipeline as pipeline_mod
import kmer_denovo_filter.vcf.pipeline as vcf_pipeline_mod
from kmer_denovo_filter.cli import parse_args
from kmer_denovo_filter.pipeline import (
    _validate_inputs,
    run_pipeline,
    run_discovery_pipeline,
)
from kmer_denovo_filter.utils import (
    _estimate_fasta_sequence_count,
    _format_elapsed,
    _format_file_size,
)

from tests.helpers import _create_bam, _create_ref_fasta, _create_vcf

GIAB_DIR = os.path.join(os.path.dirname(__file__), "..", "data", "giab")
GIAB_DATA_EXISTS = os.path.isfile(os.path.join(GIAB_DIR, "HG002_child.bam"))

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

    def test_vcf_kraken2_bacterial_fraction_annotations(self, tmpdir, monkeypatch):
        """VCF mode writes DKU_BF/DKA_BF and all non-human fraction fields."""
        chrom = "chr1"
        ref_fa = os.path.join(tmpdir, "ref.fa")
        ref_seq = _create_ref_fasta(ref_fa, chrom, 200)
        var_pos_0 = 50

        child_seq = list(ref_seq[40:80])
        child_seq[10] = "G" if ref_seq[var_pos_0] != "G" else "T"
        child_seq = "".join(child_seq)

        child_bam = os.path.join(tmpdir, "child.bam")
        _create_bam(
            child_bam, ref_fa, chrom,
            [("read1", 40, child_seq, None)],
        )

        parent_seq = ref_seq[40:80]
        mother_bam = os.path.join(tmpdir, "mother.bam")
        father_bam = os.path.join(tmpdir, "father.bam")
        _create_bam(mother_bam, ref_fa, chrom, [("mread1", 40, parent_seq, None)])
        _create_bam(father_bam, ref_fa, chrom, [("fread1", 40, parent_seq, None)])

        in_vcf = os.path.join(tmpdir, "input.vcf")
        _create_vcf(in_vcf, chrom, [(51, ref_seq[var_pos_0], child_seq[10])])
        out_vcf = os.path.join(tmpdir, "output.vcf.gz")
        kraken_db = os.path.join(tmpdir, "kraken_db")
        os.makedirs(kraken_db)

        def _mock_check_tool(_name):
            return True

        def _mock_run_kraken2(*_args, **_kwargs):
            result = pipeline_mod.Kraken2Runner.Result()
            result.total = 1
            result.classified = 1
            result.bacterial_count = 1
            result.bacterial_read_names.add("read1")
            result.nonhuman_count = 1
            result.nonhuman_read_names.add("read1")
            return result

        monkeypatch.setattr(vcf_pipeline_mod, "_check_tool", _mock_check_tool)
        monkeypatch.setattr(
            vcf_pipeline_mod, "_run_kraken2_on_reads", _mock_run_kraken2,
        )

        args = parse_args([
            "--child", child_bam,
            "--mother", mother_bam,
            "--father", father_bam,
            "--ref-fasta", ref_fa,
            "--vcf", in_vcf,
            "--output", out_vcf,
            "--kmer-size", "5",
            "--proband-id", "HG002",
            "--kraken2-db", kraken_db,
        ])
        run_pipeline(args)

        with pysam.VariantFile(out_vcf) as vcf_out:
            rec = next(iter(vcf_out))
            sample = rec.samples["HG002"]
            # Bacterial fraction
            assert "DKU_BF" in sample
            assert "DKA_BF" in sample
            assert sample["DKU_BF"] == pytest.approx(1.0)
            assert sample["DKA_BF"] == pytest.approx(1.0)
            # Archaeal fraction (0 — no archaeal reads)
            assert "DKU_AF" in sample
            assert "DKA_AF" in sample
            assert sample["DKU_AF"] == pytest.approx(0.0)
            # Fungal fraction (0)
            assert "DKU_FF" in sample
            assert "DKA_FF" in sample
            # Protist fraction (0)
            assert "DKU_PF" in sample
            assert "DKA_PF" in sample
            # Viral fraction (0)
            assert "DKU_VF" in sample
            assert "DKA_VF" in sample
            assert sample["DKU_VF"] == pytest.approx(0.0)
            assert sample["DKA_VF"] == pytest.approx(0.0)
            # UniVec Core fraction (0)
            assert "DKU_UCF" in sample
            assert "DKA_UCF" in sample
            assert sample["DKU_UCF"] == pytest.approx(0.0)
            assert sample["DKA_UCF"] == pytest.approx(0.0)
            # Non-human fraction (1.0 — same as bacterial here)
            assert "DKU_NHF" in sample
            assert "DKA_NHF" in sample
            assert sample["DKU_NHF"] == pytest.approx(1.0)
            assert sample["DKA_NHF"] == pytest.approx(1.0)
            # Unclassified fraction (0 — no unclassified reads)
            assert "DKU_UF" in sample
            assert "DKA_UF" in sample
            assert sample["DKU_UF"] == pytest.approx(0.0)
            assert sample["DKA_UF"] == pytest.approx(0.0)
            # Human lineage fraction (0 — no human lineage reads)
            assert "DKU_HLF" in sample
            assert "DKA_HLF" in sample
            assert sample["DKU_HLF"] == pytest.approx(0.0)
            assert sample["DKA_HLF"] == pytest.approx(0.0)
            # Sum-to-1 property
            dku_sum = (
                sample["DKU_NHF"] + sample["DKU_UCF"]
                + sample["DKU_HLF"] + sample["DKU_UF"]
            )
            assert dku_sum == pytest.approx(1.0)
            dka_sum = (
                sample["DKA_NHF"] + sample["DKA_UCF"]
                + sample["DKA_HLF"] + sample["DKA_UF"]
            )
            assert dka_sum == pytest.approx(1.0)

    def test_vcf_kraken2_fractions_sum_to_one_mixed(self, tmpdir, monkeypatch):
        """NHF + UCF + HLF + UF = 1.0 with a mix of read categories."""
        chrom = "chr1"
        ref_fa = os.path.join(tmpdir, "ref.fa")
        ref_seq = _create_ref_fasta(ref_fa, chrom, 200)
        var_pos_0 = 50

        child_seq = list(ref_seq[40:80])
        child_seq[10] = "G" if ref_seq[var_pos_0] != "G" else "T"
        child_seq = "".join(child_seq)

        child_bam = os.path.join(tmpdir, "child.bam")
        # Create 4 reads to cover all categories
        _create_bam(
            child_bam, ref_fa, chrom,
            [
                ("read_bact", 40, child_seq, None),
                ("read_human", 40, child_seq, None),
                ("read_uncl", 40, child_seq, None),
                ("read_ucf", 40, child_seq, None),
            ],
        )

        parent_seq = ref_seq[40:80]
        mother_bam = os.path.join(tmpdir, "mother.bam")
        father_bam = os.path.join(tmpdir, "father.bam")
        _create_bam(mother_bam, ref_fa, chrom, [("mread1", 40, parent_seq, None)])
        _create_bam(father_bam, ref_fa, chrom, [("fread1", 40, parent_seq, None)])

        in_vcf = os.path.join(tmpdir, "input.vcf")
        _create_vcf(in_vcf, chrom, [(51, ref_seq[var_pos_0], child_seq[10])])
        out_vcf = os.path.join(tmpdir, "output.vcf.gz")
        kraken_db = os.path.join(tmpdir, "kraken_db")
        os.makedirs(kraken_db)

        def _mock_check_tool(_name):
            return True

        def _mock_run_kraken2(*_args, **_kwargs):
            result = pipeline_mod.Kraken2Runner.Result()
            result.total = 4
            result.classified = 3
            result.unclassified = 1
            # Bacterial (non-human)
            result.bacterial_count = 1
            result.bacterial_read_names.add("read_bact")
            result.nonhuman_count = 1
            result.nonhuman_read_names.add("read_bact")
            # Human lineage (classified, not NHF, not UCF)
            result.human_lineage_count = 1
            result.human_lineage_read_names.add("read_human")
            result.human_count = 1
            # Unclassified
            result.unclassified_read_names.add("read_uncl")
            # UniVec Core
            result.univec_core_count = 1
            result.univec_core_read_names.add("read_ucf")
            return result

        monkeypatch.setattr(vcf_pipeline_mod, "_check_tool", _mock_check_tool)
        monkeypatch.setattr(
            vcf_pipeline_mod, "_run_kraken2_on_reads", _mock_run_kraken2,
        )

        args = parse_args([
            "--child", child_bam,
            "--mother", mother_bam,
            "--father", father_bam,
            "--ref-fasta", ref_fa,
            "--vcf", in_vcf,
            "--output", out_vcf,
            "--kmer-size", "5",
            "--proband-id", "HG002",
            "--kraken2-db", kraken_db,
        ])
        run_pipeline(args)

        with pysam.VariantFile(out_vcf) as vcf_out:
            rec = next(iter(vcf_out))
            sample = rec.samples["HG002"]
            # Each of the 4 reads should contribute 0.25 to its category
            nhf = sample["DKU_NHF"]
            ucf = sample["DKU_UCF"]
            hlf = sample["DKU_HLF"]
            uf = sample["DKU_UF"]
            assert nhf == pytest.approx(0.25)
            assert ucf == pytest.approx(0.25)
            assert hlf == pytest.approx(0.25)
            assert uf == pytest.approx(0.25)
            # Sum-to-1 property
            assert nhf + ucf + hlf + uf == pytest.approx(1.0)
            # DKA versions
            dka_nhf = sample["DKA_NHF"]
            dka_ucf = sample["DKA_UCF"]
            dka_hlf = sample["DKA_HLF"]
            dka_uf = sample["DKA_UF"]
            assert dka_nhf + dka_ucf + dka_hlf + dka_uf == pytest.approx(1.0)

    def test_kraken2_read_extraction_uses_variant_locus(self, tmpdir, monkeypatch):
        """Kraken2 read extraction uses locus-targeted fetches when available."""
        chrom = "chr1"
        ref_fa = os.path.join(tmpdir, "ref.fa")
        _create_ref_fasta(ref_fa, chrom, 200)

        read_name = "pair1"
        non_variant_seq = "A" * 40
        variant_seq = "T" * 40
        child_bam = os.path.join(tmpdir, "child.bam")
        _create_bam(
            child_bam, ref_fa, chrom,
            [
                (read_name, 10, non_variant_seq, None),
                (read_name, 50, variant_seq, None),
            ],
        )

        captured = {}

        def _mock_classify(self, sequences, tmpdir=None):
            captured.update(dict(sequences))
            return pipeline_mod.Kraken2Runner.Result()

        monkeypatch.setattr(
            pipeline_mod.Kraken2Runner, "classify_sequences", _mock_classify,
        )

        pipeline_mod._run_kraken2_on_reads(
            child_bam=child_bam,
            ref_fasta=ref_fa,
            read_names={read_name},
            kraken2_db=os.path.join(tmpdir, "kraken_db"),
            informative_reads_by_variant={"chr1:50": {read_name}},
        )

        assert set(captured) == {read_name}
        assert captured[read_name] == variant_seq

    def test_kraken2_memory_mapping_is_forwarded(self, tmpdir, monkeypatch):
        """Kraken2 memory-mapping option is forwarded to the runner."""
        chrom = "chr1"
        ref_fa = os.path.join(tmpdir, "ref.fa")
        _create_ref_fasta(ref_fa, chrom, 200)
        child_bam = os.path.join(tmpdir, "child.bam")
        _create_bam(child_bam, ref_fa, chrom, [("read1", 50, "A" * 40, None)])

        captured = {"memory_mapping": None}
        real_result_cls = pipeline_mod.Kraken2Runner.Result

        class _FakeKraken2Runner:
            Result = real_result_cls

            def __init__(
                self, _db_path, *, confidence=0.0, threads=1,
                memory_mapping=False,
            ):
                assert confidence == 0.0
                assert threads == 1
                captured["memory_mapping"] = memory_mapping

            def classify_sequences(self, _sequences, tmpdir=None):
                return self.Result()

        monkeypatch.setattr(vcf_pipeline_mod, "Kraken2Runner", _FakeKraken2Runner)

        vcf_pipeline_mod._run_kraken2_on_reads(
            child_bam=child_bam,
            ref_fasta=ref_fa,
            read_names={"read1"},
            kraken2_db=os.path.join(tmpdir, "kraken_db"),
            memory_mapping=True,
        )

        assert captured["memory_mapping"] is True

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

    def test_deletion_dka_positive(self, tmpdir):
        """A de novo deletion must yield DKA > 0 and DKA_DKT > 0.

        DKA counts fragments where at least one k-mer is absent from both
        parents *and* the read sequence exactly matches the deletion alt
        allele.  A positive DKA is the rigorous proof that the allele-
        matching step correctly attributes the deletion read to the ALT.
        """
        chrom = "chr1"
        ref_fa = os.path.join(tmpdir, "ref.fa")
        ref_seq = _create_ref_fasta(ref_fa, chrom, 200)

        anchor_0 = 50
        del_len = 2

        child_seq = ref_seq[40:anchor_0 + 1] + ref_seq[anchor_0 + 1 + del_len:80]
        child_cigar = [(0, anchor_0 + 1 - 40), (2, del_len),
                       (0, 80 - (anchor_0 + 1 + del_len))]
        child_bam = os.path.join(tmpdir, "child.bam")
        _create_bam(child_bam, ref_fa, chrom,
                    [("read1", 40, child_seq, None, child_cigar)])

        parent_seq = ref_seq[40:80]
        mother_bam = os.path.join(tmpdir, "mother.bam")
        _create_bam(mother_bam, ref_fa, chrom,
                    [("mread1", 40, parent_seq, None)])
        father_bam = os.path.join(tmpdir, "father.bam")
        _create_bam(father_bam, ref_fa, chrom,
                    [("fread1", 40, parent_seq, None)])

        ref_allele = ref_seq[anchor_0:anchor_0 + 1 + del_len]
        alt_allele = ref_seq[anchor_0]
        in_vcf = os.path.join(tmpdir, "input.vcf")
        _create_vcf(in_vcf, chrom, [(anchor_0 + 1, ref_allele, alt_allele)])

        out_vcf = os.path.join(tmpdir, "output.vcf.gz")
        args = parse_args([
            "--child", child_bam, "--mother", mother_bam,
            "--father", father_bam, "--ref-fasta", ref_fa,
            "--vcf", in_vcf, "--output", out_vcf,
            "--kmer-size", "5", "--proband-id", "HG002",
        ])
        run_pipeline(args)

        vcf_out = pysam.VariantFile(out_vcf)
        records = list(vcf_out)
        assert len(records) == 1
        rec = records[0].samples["HG002"]
        assert rec["DKU"] > 0, "Deletion read must produce DKU > 0"
        assert rec["DKA"] > 0, \
            "Deletion read must produce DKA > 0: allele-matching must " \
            "confirm that the deletion-carrying read exactly supports ALT"
        assert rec["DKA_DKT"] > 0, "DKA_DKT ratio must be positive for deletion"
        vcf_out.close()

    def test_insertion_dka_positive(self, tmpdir):
        """A de novo insertion must yield DKA > 0 and DKA_DKT > 0.

        Proves that the allele-matching step (read_supports_alt) correctly
        recognises insertion-carrying reads as supporting the ALT allele,
        producing a positive DKA score.
        """
        chrom = "chr1"
        ref_fa = os.path.join(tmpdir, "ref.fa")
        ref_seq = _create_ref_fasta(ref_fa, chrom, 200)

        anchor_0 = 50
        ins_bases = "GCA"

        child_seq = (ref_seq[40:anchor_0 + 1] + ins_bases
                     + ref_seq[anchor_0 + 1:80])
        child_cigar = [(0, anchor_0 + 1 - 40), (1, len(ins_bases)),
                       (0, 80 - (anchor_0 + 1))]
        child_bam = os.path.join(tmpdir, "child.bam")
        _create_bam(child_bam, ref_fa, chrom,
                    [("read1", 40, child_seq, None, child_cigar)])

        parent_seq = ref_seq[40:80]
        mother_bam = os.path.join(tmpdir, "mother.bam")
        _create_bam(mother_bam, ref_fa, chrom,
                    [("mread1", 40, parent_seq, None)])
        father_bam = os.path.join(tmpdir, "father.bam")
        _create_bam(father_bam, ref_fa, chrom,
                    [("fread1", 40, parent_seq, None)])

        ref_allele = ref_seq[anchor_0]
        alt_allele = ref_seq[anchor_0] + ins_bases
        in_vcf = os.path.join(tmpdir, "input.vcf")
        _create_vcf(in_vcf, chrom, [(anchor_0 + 1, ref_allele, alt_allele)])

        out_vcf = os.path.join(tmpdir, "output.vcf.gz")
        args = parse_args([
            "--child", child_bam, "--mother", mother_bam,
            "--father", father_bam, "--ref-fasta", ref_fa,
            "--vcf", in_vcf, "--output", out_vcf,
            "--kmer-size", "5", "--proband-id", "HG002",
        ])
        run_pipeline(args)

        vcf_out = pysam.VariantFile(out_vcf)
        records = list(vcf_out)
        assert len(records) == 1
        rec = records[0].samples["HG002"]
        assert rec["DKU"] > 0, "Insertion read must produce DKU > 0"
        assert rec["DKA"] > 0, \
            "Insertion read must produce DKA > 0: allele-matching must " \
            "confirm that the insertion-carrying read exactly supports ALT"
        assert rec["DKA_DKT"] > 0, "DKA_DKT ratio must be positive for insertion"
        vcf_out.close()

    def test_deletion_dka_is_allele_specific_not_all_spanning(self, tmpdir):
        """DKA counts only reads that carry the deletion, not every read
        spanning the site.

        Setup:
        - Parents have *no* reads at all (empty BAMs), making every child
          k-mer technically unique so that DKU equals DKT.
        - Child has two reads at the deletion site: one carrying the
          2-base deletion (del_read) and one matching reference (ref_read).
        - Both reads produce unique k-mers (parents are empty), so
          DKU = DKT = 2.
        - Only del_read passes read_supports_alt, so DKA = 1 < DKU.

        This rigorously demonstrates that DKA is allele-specific: even when
        a read has unique k-mers at a deletion site, it is NOT counted in
        DKA unless it actually carries the alt allele.
        """
        chrom = "chr1"
        ref_fa = os.path.join(tmpdir, "ref.fa")
        ref_seq = _create_ref_fasta(ref_fa, chrom, 200)

        anchor_0 = 50
        del_len = 2

        # Read 1: carries the 2-base deletion
        del_seq = ref_seq[40:anchor_0 + 1] + ref_seq[anchor_0 + 1 + del_len:80]
        del_cigar = [(0, anchor_0 + 1 - 40), (2, del_len),
                     (0, 80 - (anchor_0 + 1 + del_len))]

        # Read 2: reference allele (no deletion) — same region
        ref_read_seq = ref_seq[40:80]

        child_bam = os.path.join(tmpdir, "child.bam")
        _create_bam(child_bam, ref_fa, chrom, [
            ("del_read", 40, del_seq, None, del_cigar),
            ("ref_read", 40, ref_read_seq, None),
        ])

        # Parents: empty BAMs so every child k-mer is unique
        mother_bam = os.path.join(tmpdir, "mother.bam")
        _create_bam(mother_bam, ref_fa, chrom, [])
        father_bam = os.path.join(tmpdir, "father.bam")
        _create_bam(father_bam, ref_fa, chrom, [])

        ref_allele = ref_seq[anchor_0:anchor_0 + 1 + del_len]
        alt_allele = ref_seq[anchor_0]
        in_vcf = os.path.join(tmpdir, "input.vcf")
        _create_vcf(in_vcf, chrom, [(anchor_0 + 1, ref_allele, alt_allele)])

        out_vcf = os.path.join(tmpdir, "output.vcf.gz")
        args = parse_args([
            "--child", child_bam, "--mother", mother_bam,
            "--father", father_bam, "--ref-fasta", ref_fa,
            "--vcf", in_vcf, "--output", out_vcf,
            "--kmer-size", "5", "--proband-id", "HG002",
        ])
        run_pipeline(args)

        vcf_out = pysam.VariantFile(out_vcf)
        records = list(vcf_out)
        assert len(records) == 1
        rec = records[0].samples["HG002"]
        dkt = rec["DKT"]
        dku = rec["DKU"]
        dka = rec["DKA"]
        assert dkt == 2, \
            f"Both reads span the deletion site; expected DKT=2, got {dkt}"
        assert dku == 2, \
            f"Both reads have unique k-mers (empty parents); expected DKU=2, got {dku}"
        assert dka == 1, \
            f"Only the deletion read supports ALT; expected DKA=1, got {dka}"
        vcf_out.close()

    def test_insertion_dka_is_allele_specific_not_all_spanning(self, tmpdir):
        """DKA counts only reads that carry the insertion, not every read
        spanning the site.

        Setup mirrors test_deletion_dka_is_allele_specific_not_all_spanning
        but for an insertion:
        - Empty parent BAMs → all child k-mers are unique (DKU = DKT).
        - Child has an insertion read and a ref-allele read.
        - Only the insertion read passes read_supports_alt → DKA = 1 < DKU.
        """
        chrom = "chr1"
        ref_fa = os.path.join(tmpdir, "ref.fa")
        ref_seq = _create_ref_fasta(ref_fa, chrom, 200)

        anchor_0 = 50
        ins_bases = "GCA"

        # Read 1: carries the insertion
        ins_seq = (ref_seq[40:anchor_0 + 1] + ins_bases
                   + ref_seq[anchor_0 + 1:80])
        ins_cigar = [(0, anchor_0 + 1 - 40), (1, len(ins_bases)),
                     (0, 80 - (anchor_0 + 1))]

        # Read 2: reference allele (no insertion)
        ref_read_seq = ref_seq[40:80]

        child_bam = os.path.join(tmpdir, "child.bam")
        _create_bam(child_bam, ref_fa, chrom, [
            ("ins_read", 40, ins_seq, None, ins_cigar),
            ("ref_read", 40, ref_read_seq, None),
        ])

        # Parents: empty BAMs so every child k-mer is unique
        mother_bam = os.path.join(tmpdir, "mother.bam")
        _create_bam(mother_bam, ref_fa, chrom, [])
        father_bam = os.path.join(tmpdir, "father.bam")
        _create_bam(father_bam, ref_fa, chrom, [])

        ref_allele = ref_seq[anchor_0]
        alt_allele = ref_seq[anchor_0] + ins_bases
        in_vcf = os.path.join(tmpdir, "input.vcf")
        _create_vcf(in_vcf, chrom, [(anchor_0 + 1, ref_allele, alt_allele)])

        out_vcf = os.path.join(tmpdir, "output.vcf.gz")
        args = parse_args([
            "--child", child_bam, "--mother", mother_bam,
            "--father", father_bam, "--ref-fasta", ref_fa,
            "--vcf", in_vcf, "--output", out_vcf,
            "--kmer-size", "5", "--proband-id", "HG002",
        ])
        run_pipeline(args)

        vcf_out = pysam.VariantFile(out_vcf)
        records = list(vcf_out)
        assert len(records) == 1
        rec = records[0].samples["HG002"]
        dkt = rec["DKT"]
        dku = rec["DKU"]
        dka = rec["DKA"]
        assert dkt == 2, \
            f"Both reads span the insertion site; expected DKT=2, got {dkt}"
        assert dku == 2, \
            f"Both reads have unique k-mers (empty parents); expected DKU=2, got {dku}"
        assert dka == 1, \
            f"Only the insertion read supports ALT; expected DKA=1, got {dka}"
        vcf_out.close()

    def test_decomposed_indel_allele_specific_dka(self, tmpdir):
        """Two different deletion alleles at the same anchor position receive
        independent DKA scores.

        When the child carries a 1-base deletion but NOT a 3-base deletion,
        the record for the 1-base deletion must have DKA > 0 while the
        record for the 3-base deletion must have DKA = 0.

        This is the indel equivalent of test_decomposed_multiallelic_unique_
        annotations and proves that DKA allele-matching discriminates between
        different indel lengths at the same position.
        """
        chrom = "chr1"
        ref_fa = os.path.join(tmpdir, "ref.fa")
        ref_seq = _create_ref_fasta(ref_fa, chrom, 200)

        anchor_0 = 50

        # Child carries a 1-base deletion (del1)
        del1_len = 1
        del1_seq = ref_seq[40:anchor_0 + 1] + ref_seq[anchor_0 + 1 + del1_len:80]
        del1_cigar = [(0, anchor_0 + 1 - 40), (2, del1_len),
                      (0, 80 - (anchor_0 + 1 + del1_len))]

        child_bam = os.path.join(tmpdir, "child.bam")
        _create_bam(child_bam, ref_fa, chrom,
                    [("read1", 40, del1_seq, None, del1_cigar)])

        # Parents match reference (neither deletion present)
        parent_seq = ref_seq[40:80]
        mother_bam = os.path.join(tmpdir, "mother.bam")
        _create_bam(mother_bam, ref_fa, chrom,
                    [("mread1", 40, parent_seq, None)])
        father_bam = os.path.join(tmpdir, "father.bam")
        _create_bam(father_bam, ref_fa, chrom,
                    [("fread1", 40, parent_seq, None)])

        # VCF: two decomposed records at the same anchor position
        del1_ref = ref_seq[anchor_0:anchor_0 + 1 + del1_len]  # 2-base REF
        del1_alt = ref_seq[anchor_0]                           # 1-base ALT

        del2_len = 3
        del2_ref = ref_seq[anchor_0:anchor_0 + 1 + del2_len]  # 4-base REF
        del2_alt = ref_seq[anchor_0]                           # 1-base ALT

        in_vcf = os.path.join(tmpdir, "input.vcf")
        _create_vcf(in_vcf, chrom, [
            (anchor_0 + 1, del1_ref, del1_alt),
            (anchor_0 + 1, del2_ref, del2_alt),
        ])

        out_vcf = os.path.join(tmpdir, "output.vcf.gz")
        args = parse_args([
            "--child", child_bam, "--mother", mother_bam,
            "--father", father_bam, "--ref-fasta", ref_fa,
            "--vcf", in_vcf, "--output", out_vcf,
            "--kmer-size", "5", "--proband-id", "HG002",
        ])
        run_pipeline(args)

        vcf_out = pysam.VariantFile(out_vcf)
        records = list(vcf_out)
        assert len(records) == 2, "Both decomposed records should be present"

        rec_del1 = [r for r in records if r.alts and r.alts[0] == del1_alt
                    and len(r.ref) == len(del1_ref)]
        rec_del2 = [r for r in records if r.alts and r.alts[0] == del2_alt
                    and len(r.ref) == len(del2_ref)]
        assert len(rec_del1) == 1, "Expected one record for the 1-base deletion"
        assert len(rec_del2) == 1, "Expected one record for the 3-base deletion"

        # 1-base deletion: child carries it → DKA > 0
        dka_del1 = rec_del1[0].samples["HG002"]["DKA"]
        assert dka_del1 > 0, \
            f"Child carries 1-base deletion; expected DKA > 0, got {dka_del1}"

        # 3-base deletion: child does NOT carry it → DKA = 0
        dka_del2 = rec_del2[0].samples["HG002"]["DKA"]
        assert dka_del2 == 0, \
            f"Child does NOT carry 3-base deletion; expected DKA = 0, got {dka_del2}"
        vcf_out.close()

    def test_decomposed_insertion_allele_specific_dka(self, tmpdir):
        """Two insertion alleles at one position receive independent DKA.

        Child carries only the 3-base insertion allele; a second 1-base
        insertion record at the same anchor must remain DKA=0.
        """
        chrom = "chr1"
        ref_fa = os.path.join(tmpdir, "ref.fa")
        ref_seq = _create_ref_fasta(ref_fa, chrom, 200)

        anchor_0 = 50
        ins3 = "GCA"
        child_seq = (ref_seq[40:anchor_0 + 1] + ins3
                     + ref_seq[anchor_0 + 1:80])
        child_cigar = [(0, anchor_0 + 1 - 40), (1, len(ins3)),
                       (0, 80 - (anchor_0 + 1))]

        child_bam = os.path.join(tmpdir, "child.bam")
        _create_bam(child_bam, ref_fa, chrom,
                    [("read1", 40, child_seq, None, child_cigar)])

        parent_seq = ref_seq[40:80]
        mother_bam = os.path.join(tmpdir, "mother.bam")
        _create_bam(mother_bam, ref_fa, chrom,
                    [("mread1", 40, parent_seq, None)])
        father_bam = os.path.join(tmpdir, "father.bam")
        _create_bam(father_bam, ref_fa, chrom,
                    [("fread1", 40, parent_seq, None)])

        ref_allele = ref_seq[anchor_0]
        alt_ins3 = ref_allele + ins3
        alt_ins1 = ref_allele + "T"
        in_vcf = os.path.join(tmpdir, "input.vcf")
        _create_vcf(in_vcf, chrom, [
            (anchor_0 + 1, ref_allele, alt_ins3),
            (anchor_0 + 1, ref_allele, alt_ins1),
        ])

        out_vcf = os.path.join(tmpdir, "output.vcf.gz")
        args = parse_args([
            "--child", child_bam, "--mother", mother_bam,
            "--father", father_bam, "--ref-fasta", ref_fa,
            "--vcf", in_vcf, "--output", out_vcf,
            "--kmer-size", "5", "--proband-id", "HG002",
        ])
        run_pipeline(args)

        vcf_out = pysam.VariantFile(out_vcf)
        records = list(vcf_out)
        assert len(records) == 2, "Both decomposed insertion records should be present"

        rec_ins3 = [r for r in records if r.alts and r.alts[0] == alt_ins3]
        rec_ins1 = [r for r in records if r.alts and r.alts[0] == alt_ins1]
        assert len(rec_ins3) == 1, "Expected one record for the 3-base insertion"
        assert len(rec_ins1) == 1, "Expected one record for the 1-base insertion"

        dka_ins3 = rec_ins3[0].samples["HG002"]["DKA"]
        assert dka_ins3 > 0, \
            f"Child carries 3-base insertion; expected DKA > 0, got {dka_ins3}"

        dka_ins1 = rec_ins1[0].samples["HG002"]["DKA"]
        assert dka_ins1 == 0, \
            f"Child does NOT carry 1-base insertion; expected DKA = 0, got {dka_ins1}"
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

    def test_decomposed_multiallelic_unique_annotations(self, tmpdir):
        """Two VCF records at the same position with different ALTs should
        each receive independent annotations (variant key includes ref:alt).
        """
        chrom = "chr1"
        ref_fa = os.path.join(tmpdir, "ref.fa")
        ref_seq = _create_ref_fasta(ref_fa, chrom, 200)

        # Position 51 (1-based), ref_base at 0-based 50
        ref_base = ref_seq[50]
        # Pick two different ALT bases
        alt_bases = [b for b in "ACGT" if b != ref_base]
        alt1 = alt_bases[0]
        alt2 = alt_bases[1]

        # Child has alt1 at position 50 (de novo for alt1)
        child_seq = list(ref_seq[40:80])
        child_seq[10] = alt1
        child_seq = "".join(child_seq)

        child_bam = os.path.join(tmpdir, "child.bam")
        _create_bam(
            child_bam, ref_fa, chrom,
            [("read1", 40, child_seq, None)],
        )

        # Parents match reference
        parent_seq = ref_seq[40:80]
        mother_bam = os.path.join(tmpdir, "mother.bam")
        _create_bam(mother_bam, ref_fa, chrom,
                     [("mread1", 40, parent_seq, None)])
        father_bam = os.path.join(tmpdir, "father.bam")
        _create_bam(father_bam, ref_fa, chrom,
                     [("fread1", 40, parent_seq, None)])

        # Decomposed multiallelic: two records at the same position
        in_vcf = os.path.join(tmpdir, "input.vcf")
        _create_vcf(in_vcf, chrom, [
            (51, ref_base, alt1),   # child carries this allele
            (51, ref_base, alt2),   # child does NOT carry this allele
        ])

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
        assert len(records) == 2, "Both decomposed records should be present"

        rec_alt1 = [r for r in records if r.alts and r.alts[0] == alt1]
        rec_alt2 = [r for r in records if r.alts and r.alts[0] == alt2]
        assert len(rec_alt1) == 1, f"Expected one record for ALT={alt1}"
        assert len(rec_alt2) == 1, f"Expected one record for ALT={alt2}"

        # alt1 record: child carries this allele → DKA > 0
        assert rec_alt1[0].samples["HG002"]["DKA"] > 0
        # alt2 record: child does NOT carry this allele → DKA == 0
        assert rec_alt2[0].samples["HG002"]["DKA"] == 0
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





class TestEstimateFastaSequenceCount:
    """Unit tests for sampled FASTA entry count estimation."""

    def test_empty_file(self, tmp_path):
        fa = tmp_path / "empty.fa"
        fa.write_bytes(b"")
        assert _estimate_fasta_sequence_count(str(fa)) == (0, False)

    def test_small_file_exact_count(self, tmp_path):
        fa = tmp_path / "small.fa"
        fa.write_text(">0\nAAAA\n>1\nCCCC\n>2\nGGGG\n")
        assert _estimate_fasta_sequence_count(str(fa)) == (3, False)

    def test_large_file_extrapolates(self, tmp_path):
        fa = tmp_path / "large.fa"
        n_entries = 1500
        with open(fa, "w") as fh:
            for i in range(n_entries):
                fh.write(f">{i:05d}\n")
                fh.write("ACGTACGTACGTACGTACGT\n")

        estimate, extrapolated = _estimate_fasta_sequence_count(
            str(fa), sample_lines=1000,
        )

        assert extrapolated is True
        assert estimate == n_entries


class TestScanParentHashSizing:
    """Verify that _scan_parent_jellyfish sizes hash dynamically."""

    def test_hash_size_from_n_filter_kmers(self, tmp_path, monkeypatch):
        """When n_filter_kmers is provided, hash_size = max(2*n, 10M)."""
        from unittest.mock import MagicMock, patch, call
        from kmer_denovo_filter.core import jellyfish_wrappers as jfw

        # Create dummy files
        kmer_fa = tmp_path / "kmers.fa"
        kmer_fa.write_text(">0\nACGT\n>1\nTTTT\n")
        parent_bam = tmp_path / "parent.bam"
        parent_bam.write_bytes(b"dummy")
        parent_dir = str(tmp_path / "pdir")

        # Build a mock that simulates jellyfish producing a single jf file
        captured_jf_cmd = []

        class FakeProc:
            def __init__(self, cmd, **kw):
                if cmd[0] == "jellyfish" and cmd[1] == "count":
                    captured_jf_cmd.clear()
                    captured_jf_cmd.extend(cmd)
                    # Create the output file so the dump path is exercised
                    jf_out = cmd[cmd.index("-o") + 1]
                    open(jf_out, "wb").close()
                self.stdout = MagicMock()
                self.stdout.close = MagicMock()
                self.stderr = MagicMock()
                self.stderr.read = MagicMock(return_value=b"")
                self.returncode = 0
                self.pid = 1234

            def wait(self, timeout=None):
                pass

            def communicate(self):
                return (b"", b"")

        # Patch subprocess.Popen and the dump Popen
        with patch.object(jfw.subprocess, "Popen", side_effect=FakeProc):
            # Also mock the dump to avoid real jellyfish dump
            with patch.object(jfw.os.path, "exists", return_value=False):
                result = jfw._scan_parent_jellyfish(
                    str(parent_bam), None, str(kmer_fa), 31,
                    parent_dir, threads=4,
                    n_filter_kmers=15_000_000,
                )

        # Verify hash size: max(15M * 2, 10M) = 30M
        assert "-s" in captured_jf_cmd
        s_idx = captured_jf_cmd.index("-s")
        assert captured_jf_cmd[s_idx + 1] == "30000000"

    def test_hash_size_defaults_to_10m_minimum(self, tmp_path, monkeypatch):
        """When n_filter_kmers is small, hash stays at 10M minimum."""
        from unittest.mock import MagicMock, patch
        from kmer_denovo_filter.core import jellyfish_wrappers as jfw

        kmer_fa = tmp_path / "kmers.fa"
        kmer_fa.write_text(">0\nACGT\n>1\nTTTT\n")
        parent_bam = tmp_path / "parent.bam"
        parent_bam.write_bytes(b"dummy")
        parent_dir = str(tmp_path / "pdir")

        captured_jf_cmd = []

        class FakeProc:
            def __init__(self, cmd, **kw):
                if cmd[0] == "jellyfish" and cmd[1] == "count":
                    captured_jf_cmd.clear()
                    captured_jf_cmd.extend(cmd)
                    jf_out = cmd[cmd.index("-o") + 1]
                    open(jf_out, "wb").close()
                self.stdout = MagicMock()
                self.stdout.close = MagicMock()
                self.stderr = MagicMock()
                self.stderr.read = MagicMock(return_value=b"")
                self.returncode = 0
                self.pid = 1234

            def wait(self, timeout=None):
                pass

            def communicate(self):
                return (b"", b"")

        with patch.object(jfw.subprocess, "Popen", side_effect=FakeProc):
            with patch.object(jfw.os.path, "exists", return_value=False):
                result = jfw._scan_parent_jellyfish(
                    str(parent_bam), None, str(kmer_fa), 31,
                    parent_dir, threads=4,
                    n_filter_kmers=100,
                )

        # max(100*2, 10M) = 10M
        s_idx = captured_jf_cmd.index("-s")
        assert captured_jf_cmd[s_idx + 1] == "10000000"

    def test_counts_fasta_entries_when_n_filter_kmers_is_none(
        self, tmp_path, monkeypatch,
    ):
        """When n_filter_kmers is None, k-mer count is read from FASTA."""
        from unittest.mock import MagicMock, patch
        from kmer_denovo_filter.core import jellyfish_wrappers as jfw

        # Write a FASTA with 5 entries
        kmer_fa = tmp_path / "kmers.fa"
        entries = "".join(f">{i}\nACGTACGT\n" for i in range(5))
        kmer_fa.write_text(entries)
        parent_bam = tmp_path / "parent.bam"
        parent_bam.write_bytes(b"dummy")
        parent_dir = str(tmp_path / "pdir")

        captured_jf_cmd = []

        class FakeProc:
            def __init__(self, cmd, **kw):
                if cmd[0] == "jellyfish" and cmd[1] == "count":
                    captured_jf_cmd.clear()
                    captured_jf_cmd.extend(cmd)
                    jf_out = cmd[cmd.index("-o") + 1]
                    open(jf_out, "wb").close()
                self.stdout = MagicMock()
                self.stdout.close = MagicMock()
                self.stderr = MagicMock()
                self.stderr.read = MagicMock(return_value=b"")
                self.returncode = 0
                self.pid = 1234

            def wait(self, timeout=None):
                pass

            def communicate(self):
                return (b"", b"")

        with patch.object(jfw.subprocess, "Popen", side_effect=FakeProc):
            with patch.object(jfw.os.path, "exists", return_value=False):
                result = jfw._scan_parent_jellyfish(
                    str(parent_bam), None, str(kmer_fa), 31,
                    parent_dir, threads=4,
                )

        # 5 entries → max(5*2, 10M) = 10M
        s_idx = captured_jf_cmd.index("-s")
        assert captured_jf_cmd[s_idx + 1] == "10000000"



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
        """Run the full pipeline on GIAB candidates (SNVs + SV-like DNMs)."""
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
        assert len(records) == 22
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
        assert metrics["total_variants"] == 22
        assert metrics["variants_with_unique_reads"] > 0

        # Summary file should exist and contain key sections
        assert os.path.exists(summary_txt)
        with open(summary_txt) as fh:
            summary = fh.read()
        assert "De Novo Variant Summary" in summary
        assert "Likely de novo" in summary
        assert "DE_NOVO" in summary





class TestModuleSeparation:
    """Verify that vcf/ and discovery/ sub-packages are properly isolated."""

    def test_discovery_does_not_import_vcf(self):
        """The discovery sub-package must never import from vcf/."""
        import ast

        discovery_pipeline = os.path.join(
            os.path.dirname(__file__), "..", "..", "src",
            "kmer_denovo_filter", "discovery", "pipeline.py",
        )
        with open(discovery_pipeline) as f:
            tree = ast.parse(f.read(), filename=discovery_pipeline)

        for node in ast.walk(tree):
            if isinstance(node, ast.ImportFrom) and node.module:
                assert "kmer_denovo_filter.vcf" not in node.module, (
                    f"discovery/pipeline.py imports from vcf/: {node.module}"
                )
            if isinstance(node, ast.Import):
                for alias in node.names:
                    assert "kmer_denovo_filter.vcf" not in alias.name, (
                        f"discovery/pipeline.py imports from vcf/: {alias.name}"
                    )

    def test_subpackages_exist(self):
        """Both vcf/ and discovery/ sub-packages should be importable."""
        import kmer_denovo_filter.vcf
        import kmer_denovo_filter.discovery
        assert hasattr(kmer_denovo_filter.vcf, "run_pipeline")
        assert hasattr(kmer_denovo_filter.discovery, "run_discovery_pipeline")

    def test_backward_compat_pipeline_reexports(self):
        """The pipeline.py shim re-exports all public names."""
        from kmer_denovo_filter import pipeline
        # VCF-mode functions
        assert callable(pipeline.run_pipeline)
        assert callable(pipeline._write_annotated_vcf)
        assert callable(pipeline._parse_vcf_variants)
        # Discovery-mode functions
        assert callable(pipeline.run_discovery_pipeline)
        assert callable(pipeline._anchor_and_cluster)
        assert callable(pipeline._parse_candidate_summary)
        assert hasattr(pipeline, "SULOVARI_DNM_REGIONS")
        # Shared
        assert callable(pipeline._validate_inputs)
        assert callable(pipeline._check_tool)
