"""Discovery-mode pipeline tests (split from tests/test_pipeline.py)."""

import collections
import json
import os
import tempfile

import pysam
import pytest

import kmer_denovo_filter.pipeline as pipeline_mod
import kmer_denovo_filter.discovery.pipeline as discovery_pipeline_mod
import kmer_denovo_filter.core.bam_scanner as bam_scanner_mod
from kmer_denovo_filter.cli import parse_args
from kmer_denovo_filter.pipeline import (
    _classify_regions,
    _validate_inputs,
    _write_bedgraph,
    _write_bedpe,
    run_discovery_pipeline,
)
from kmer_denovo_filter.utils import (
    _collect_kmer_ref_positions,
)

from tests.helpers import (
    _create_bam,
    _create_bam_with_supplementary,
    _create_ref_fasta,
)

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
            bed_lines = [l.strip() for l in fh if l.strip() and not l.startswith("#")]
        assert len(bed_lines) >= 1, "Expected at least one candidate region"
        # BED lines should have 10 columns: chrom, start, end, reads, kmers,
        # split_reads, discordant_pairs, max_clip_len, unmapped_mates, class
        parts = bed_lines[0].split("\t")
        assert len(parts) == 10
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
            # SV annotation fields
            assert "split_reads" in region
            assert "discordant_pairs" in region
            assert "max_clip_len" in region
            assert "unmapped_mates" in region
            assert "class" in region
            assert region["class"] in ("SV", "SMALL", "AMBIGUOUS")

        # Check BEDPE file exists
        bedpe_path = f"{out_prefix}.sv.bedpe"
        assert os.path.exists(bedpe_path)

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
        assert "Split" in summary
        assert "Class" in summary

    def test_discovery_ignores_kraken2_for_now(self, tmpdir, monkeypatch):
        """Discovery mode should not invoke kraken2 even if args are provided."""
        chrom = "chr1"
        ref_fa = os.path.join(tmpdir, "ref.fa")
        ref_seq = _create_ref_fasta(ref_fa, chrom, 200)

        mut_seq = list(ref_seq[30:90])
        mut_seq[20] = "G" if ref_seq[50] != "G" else "T"
        mut_seq = "".join(mut_seq)

        child_bam = os.path.join(tmpdir, "child.bam")
        _create_bam(
            child_bam, ref_fa, chrom,
            [("read1", 30, mut_seq, None), ("read2", 30, mut_seq, None)],
        )
        mother_bam = os.path.join(tmpdir, "mother.bam")
        father_bam = os.path.join(tmpdir, "father.bam")
        _create_bam(mother_bam, ref_fa, chrom, [("mread1", 30, ref_seq[30:90], None)])
        _create_bam(father_bam, ref_fa, chrom, [("fread1", 30, ref_seq[30:90], None)])

        out_prefix = os.path.join(tmpdir, "disc_no_kraken")
        kraken_db = os.path.join(tmpdir, "kraken_db")
        os.makedirs(kraken_db)

        original_check_tool = pipeline_mod._check_tool

        called = {
            "checked_kraken2": False,
            "ran_kraken2": False,
        }

        def _check_tool_no_kraken(name):
            if name == "kraken2":
                called["checked_kraken2"] = True
                return False
            return original_check_tool(name)

        def _track_if_called(*_args, **_kwargs):
            called["ran_kraken2"] = True
            result = pipeline_mod.Kraken2Runner.Result()
            return result

        monkeypatch.setattr(discovery_pipeline_mod, "_check_tool", _check_tool_no_kraken)

        args = parse_args([
            "--child", child_bam,
            "--mother", mother_bam,
            "--father", father_bam,
            "--ref-fasta", ref_fa,
            "--out-prefix", out_prefix,
            "--min-child-count", "2",
            "--kmer-size", "5",
            "--kraken2-db", kraken_db,
        ])
        run_discovery_pipeline(args)

        metrics_path = f"{out_prefix}.metrics.json"
        with open(metrics_path) as fh:
            metrics = json.load(fh)
        assert "kraken2" not in metrics
        assert called["checked_kraken2"] is False
        assert called["ran_kraken2"] is False

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
            bed_lines = [l.strip() for l in fh if l.strip() and not l.startswith("#")]
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
            data_lines = [
                l for l in fh if l.strip() and not l.startswith("#")
            ]
        assert len(data_lines) == 0

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

    def test_discovery_min_supporting_reads_filters(self, tmpdir):
        """Regions with fewer reads than --min-supporting-reads are removed."""
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
            [("mread1", 30, parent_seq, None)],
        )
        father_bam = os.path.join(tmpdir, "father.bam")
        _create_bam(
            father_bam, ref_fa, chrom,
            [("fread1", 30, parent_seq, None)],
        )

        # First run without the filter — should produce regions
        out1 = os.path.join(tmpdir, "disc_nofilt")
        args1 = parse_args([
            "--child", child_bam, "--mother", mother_bam,
            "--father", father_bam, "--ref-fasta", ref_fa,
            "--out-prefix", out1, "--min-child-count", "3",
            "--kmer-size", "5",
        ])
        run_discovery_pipeline(args1)
        with open(f"{out1}.bed") as fh:
            baseline = [l.strip() for l in fh if l.strip() and not l.startswith("#")]
        assert len(baseline) >= 1

        # Run with very high --min-supporting-reads — should filter all
        out2 = os.path.join(tmpdir, "disc_hifilt")
        args2 = parse_args([
            "--child", child_bam, "--mother", mother_bam,
            "--father", father_bam, "--ref-fasta", ref_fa,
            "--out-prefix", out2, "--min-child-count", "3",
            "--kmer-size", "5",
            "--min-supporting-reads", "100",
        ])
        run_discovery_pipeline(args2)
        with open(f"{out2}.bed") as fh:
            filtered = [l.strip() for l in fh if l.strip() and not l.startswith("#")]
        assert len(filtered) == 0

    def test_discovery_min_distinct_kmers_filters(self, tmpdir):
        """Regions with fewer k-mers than --min-distinct-kmers are removed."""
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
            [("mread1", 30, parent_seq, None)],
        )
        father_bam = os.path.join(tmpdir, "father.bam")
        _create_bam(
            father_bam, ref_fa, chrom,
            [("fread1", 30, parent_seq, None)],
        )

        out = os.path.join(tmpdir, "disc_kfilt")
        args = parse_args([
            "--child", child_bam, "--mother", mother_bam,
            "--father", father_bam, "--ref-fasta", ref_fa,
            "--out-prefix", out, "--min-child-count", "3",
            "--kmer-size", "5",
            "--min-distinct-kmers", "10000",
        ])
        run_discovery_pipeline(args)
        with open(f"{out}.bed") as fh:
            filtered = [l.strip() for l in fh if l.strip() and not l.startswith("#")]
        assert len(filtered) == 0

    def test_discovery_min_distinct_kmers_per_read_filters(self, tmpdir):
        """Reads with fewer distinct kmers than --min-distinct-kmers-per-read
        are excluded before region clustering and bedGraph generation."""
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
            [("mread1", 30, parent_seq, None)],
        )
        father_bam = os.path.join(tmpdir, "father.bam")
        _create_bam(
            father_bam, ref_fa, chrom,
            [("fread1", 30, parent_seq, None)],
        )

        # Run with a very high per-read threshold — should filter all reads
        out = os.path.join(tmpdir, "disc_pread")
        args = parse_args([
            "--child", child_bam, "--mother", mother_bam,
            "--father", father_bam, "--ref-fasta", ref_fa,
            "--out-prefix", out, "--min-child-count", "3",
            "--kmer-size", "5",
            "--min-distinct-kmers-per-read", "10000",
        ])
        run_discovery_pipeline(args)
        with open(f"{out}.bed") as fh:
            filtered = [l.strip() for l in fh
                        if l.strip() and not l.startswith("#")]
        assert len(filtered) == 0

    def test_discovery_min_distinct_kmers_per_read_default(self, tmpdir):
        """Default --min-distinct-kmers-per-read is k/4."""
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
            [("mread1", 30, parent_seq, None)],
        )
        father_bam = os.path.join(tmpdir, "father.bam")
        _create_bam(
            father_bam, ref_fa, chrom,
            [("fread1", 30, parent_seq, None)],
        )

        # k=5, so default min_distinct_kmers_per_read = 5//4 = 1
        out = os.path.join(tmpdir, "disc_default")
        args = parse_args([
            "--child", child_bam, "--mother", mother_bam,
            "--father", father_bam, "--ref-fasta", ref_fa,
            "--out-prefix", out, "--min-child-count", "3",
            "--kmer-size", "5",
        ])
        run_discovery_pipeline(args)

        with open(f"{out}.metrics.json") as fh:
            metrics = json.load(fh)
        assert metrics["filters"]["min_distinct_kmers_per_read"] == 1

    def test_discovery_bed_includes_filter_header(self, tmpdir):
        """BED file should include a #filters header comment."""
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
            [("mread1", 30, parent_seq, None)],
        )
        father_bam = os.path.join(tmpdir, "father.bam")
        _create_bam(
            father_bam, ref_fa, chrom,
            [("fread1", 30, parent_seq, None)],
        )

        out = os.path.join(tmpdir, "disc_fhdr")
        args = parse_args([
            "--child", child_bam, "--mother", mother_bam,
            "--father", father_bam, "--ref-fasta", ref_fa,
            "--out-prefix", out, "--min-child-count", "3",
            "--kmer-size", "5",
        ])
        run_discovery_pipeline(args)

        with open(f"{out}.bed") as fh:
            header_lines = [l for l in fh if l.startswith("#filters:")]
        assert len(header_lines) == 1
        hdr = header_lines[0]
        assert "min_distinct_kmers_per_read=" in hdr
        assert "min_supporting_reads=" in hdr
        assert "min_distinct_kmers=" in hdr

    def test_discovery_cluster_distance(self, tmpdir):
        """--cluster-distance is accepted and passed to the pipeline."""
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
            [("mread1", 30, parent_seq, None)],
        )
        father_bam = os.path.join(tmpdir, "father.bam")
        _create_bam(
            father_bam, ref_fa, chrom,
            [("fread1", 30, parent_seq, None)],
        )

        out = os.path.join(tmpdir, "disc_cd")
        args = parse_args([
            "--child", child_bam, "--mother", mother_bam,
            "--father", father_bam, "--ref-fasta", ref_fa,
            "--out-prefix", out, "--min-child-count", "3",
            "--kmer-size", "5",
            "--cluster-distance", "100",
        ])
        run_discovery_pipeline(args)
        bed_path = f"{out}.bed"
        assert os.path.exists(bed_path)
        with open(bed_path) as fh:
            bed_lines = [l.strip() for l in fh if l.strip() and not l.startswith("#")]
        assert len(bed_lines) >= 1

    def test_discovery_parent_max_count(self, tmpdir):
        """--parent-max-count allows low-count parental k-mers through."""
        chrom = "chr1"
        ref_fa = os.path.join(tmpdir, "ref.fa")
        ref_seq = _create_ref_fasta(ref_fa, chrom, 200)

        # Child mutation at position 50
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

        # Mother has one read with the mutation (count=1 in parent)
        mother_bam = os.path.join(tmpdir, "mother.bam")
        _create_bam(
            mother_bam, ref_fa, chrom,
            [("mread1", 30, child_seq, None)],
        )
        father_bam = os.path.join(tmpdir, "father.bam")
        _create_bam(
            father_bam, ref_fa, chrom,
            [("fread1", 30, ref_seq[30:90], None)],
        )

        # Default parent-max-count=0: mutation k-mers found in mother (count=1)
        # should be removed → no regions
        out_strict = os.path.join(tmpdir, "disc_strict")
        args_strict = parse_args([
            "--child", child_bam, "--mother", mother_bam,
            "--father", father_bam, "--ref-fasta", ref_fa,
            "--out-prefix", out_strict, "--min-child-count", "3",
            "--kmer-size", "5", "--parent-max-count", "0",
        ])
        run_discovery_pipeline(args_strict)
        with open(f"{out_strict}.bed") as fh:
            strict_lines = [l.strip() for l in fh if l.strip() and not l.startswith("#")]
        assert len(strict_lines) == 0, (
            "With parent-max-count=0, k-mers present in mother should be removed"
        )

        # Relaxed parent-max-count=1: k-mers with count=1 in mother are
        # tolerated → regions should appear
        out_relaxed = os.path.join(tmpdir, "disc_relaxed")
        args_relaxed = parse_args([
            "--child", child_bam, "--mother", mother_bam,
            "--father", father_bam, "--ref-fasta", ref_fa,
            "--out-prefix", out_relaxed, "--min-child-count", "3",
            "--kmer-size", "5", "--parent-max-count", "1",
        ])
        run_discovery_pipeline(args_relaxed)
        with open(f"{out_relaxed}.bed") as fh:
            relaxed_lines = [l.strip() for l in fh if l.strip() and not l.startswith("#")]
        assert len(relaxed_lines) >= 1, (
            "With parent-max-count=1, k-mers with count=1 in mother "
            "should be tolerated"
        )

    def test_discovery_multiple_regions(self, tmpdir):
        """Two distant mutations should produce two separate BED regions."""
        chrom = "chr1"
        ref_fa = os.path.join(tmpdir, "ref.fa")
        ref_seq = _create_ref_fasta(ref_fa, chrom, 200)

        # Mutation 1 at position 50 (read covering 30-90)
        child_seq1 = list(ref_seq[30:90])
        child_seq1[20] = "G" if ref_seq[50] != "G" else "T"
        child_seq1 = "".join(child_seq1)

        # Mutation 2 at position 150 (read covering 130-190)
        child_seq2 = list(ref_seq[130:190])
        child_seq2[20] = "G" if ref_seq[150] != "G" else "T"
        child_seq2 = "".join(child_seq2)

        child_bam = os.path.join(tmpdir, "child.bam")
        _create_bam(
            child_bam, ref_fa, chrom,
            [
                ("read1a", 30, child_seq1, None),
                ("read2a", 30, child_seq1, None),
                ("read3a", 30, child_seq1, None),
                ("read4a", 30, child_seq1, None),
                ("read1b", 130, child_seq2, None),
                ("read2b", 130, child_seq2, None),
                ("read3b", 130, child_seq2, None),
                ("read4b", 130, child_seq2, None),
            ],
        )

        parent_seq1 = ref_seq[30:90]
        parent_seq2 = ref_seq[130:190]
        mother_bam = os.path.join(tmpdir, "mother.bam")
        _create_bam(
            mother_bam, ref_fa, chrom,
            [
                ("mread1", 30, parent_seq1, None),
                ("mread2", 130, parent_seq2, None),
            ],
        )
        father_bam = os.path.join(tmpdir, "father.bam")
        _create_bam(
            father_bam, ref_fa, chrom,
            [
                ("fread1", 30, parent_seq1, None),
                ("fread2", 130, parent_seq2, None),
            ],
        )

        # Use cluster-distance=0 so the two mutations stay separate
        out = os.path.join(tmpdir, "disc_multi")
        args = parse_args([
            "--child", child_bam, "--mother", mother_bam,
            "--father", father_bam, "--ref-fasta", ref_fa,
            "--out-prefix", out, "--min-child-count", "3",
            "--kmer-size", "5", "--cluster-distance", "0",
        ])
        run_discovery_pipeline(args)

        bed_path = f"{out}.bed"
        with open(bed_path) as fh:
            bed_lines = [l.strip() for l in fh if l.strip() and not l.startswith("#")]
        assert len(bed_lines) == 2, (
            f"Expected 2 regions for 2 distant mutations, got {len(bed_lines)}"
        )

        # Verify metrics match
        metrics_path = f"{out}.metrics.json"
        with open(metrics_path) as fh:
            metrics = json.load(fh)
        assert metrics["candidate_regions"] == 2
        assert len(metrics["regions"]) == 2





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


class TestDiscoverySV:
    """Tests for SV support in discovery mode."""

    @pytest.fixture()
    def tmpdir(self, tmp_path):
        return str(tmp_path)

    def test_supplementary_read_in_informative_bam(self, tmpdir):
        """A supplementary read carrying proband-unique k-mers should appear
        in the informative BAM."""
        chrom = "chr1"
        ref_fa = os.path.join(tmpdir, "ref.fa")
        ref_seq = _create_ref_fasta(ref_fa, chrom, 200)

        # Create a mutation at position 50
        child_seq = list(ref_seq[30:90])
        child_seq[20] = "G" if ref_seq[50] != "G" else "T"
        child_seq = "".join(child_seq)

        # Create child BAM with primary + supplementary for same read name
        child_bam = os.path.join(tmpdir, "child.bam")
        _create_bam_with_supplementary(
            child_bam, ref_fa, [chrom], [300],
            [
                # Primary reads
                {"name": "read1", "pos": 30, "seq": child_seq, "flag": 0},
                {"name": "read2", "pos": 30, "seq": child_seq, "flag": 0},
                {"name": "read3", "pos": 30, "seq": child_seq, "flag": 0},
                {"name": "read4", "pos": 30, "seq": child_seq, "flag": 0},
                # Supplementary read carrying same mutation sequence
                {"name": "read1", "pos": 30, "seq": child_seq,
                 "flag": 0x800},  # supplementary flag
            ],
        )

        parent_seq = ref_seq[30:90]
        mother_bam = os.path.join(tmpdir, "mother.bam")
        _create_bam(mother_bam, ref_fa, chrom,
                     [("mread1", 30, parent_seq, None)])
        father_bam = os.path.join(tmpdir, "father.bam")
        _create_bam(father_bam, ref_fa, chrom,
                     [("fread1", 30, parent_seq, None)])

        out_prefix = os.path.join(tmpdir, "sv_supp")
        args = parse_args([
            "--child", child_bam, "--mother", mother_bam,
            "--father", father_bam, "--ref-fasta", ref_fa,
            "--out-prefix", out_prefix, "--min-child-count", "3",
            "--kmer-size", "5",
        ])
        run_discovery_pipeline(args)

        # Check informative BAM includes supplementary
        info_bam = f"{out_prefix}.informative.bam"
        bam_info = pysam.AlignmentFile(info_bam)
        info_reads = list(bam_info)
        bam_info.close()

        # Should have reads including the supplementary
        supp_reads = [r for r in info_reads if r.is_supplementary]
        primary_reads = [r for r in info_reads if not r.is_supplementary]
        assert len(primary_reads) >= 1
        assert len(supp_reads) >= 1
        # All reads should have dk tag
        for read in info_reads:
            assert read.has_tag("dk")
            assert read.get_tag("dk") == 1

    def test_unlinked_snp_region_is_small(self, tmpdir):
        """A simple point mutation (no split reads) should be classified SMALL."""
        chrom = "chr1"
        ref_fa = os.path.join(tmpdir, "ref.fa")
        ref_seq = _create_ref_fasta(ref_fa, chrom, 200)

        child_seq = list(ref_seq[30:90])
        child_seq[20] = "G" if ref_seq[50] != "G" else "T"
        child_seq = "".join(child_seq)

        child_bam = os.path.join(tmpdir, "child.bam")
        _create_bam(child_bam, ref_fa, chrom, [
            ("read1", 30, child_seq, None),
            ("read2", 30, child_seq, None),
            ("read3", 30, child_seq, None),
            ("read4", 30, child_seq, None),
        ])

        parent_seq = ref_seq[30:90]
        mother_bam = os.path.join(tmpdir, "mother.bam")
        _create_bam(mother_bam, ref_fa, chrom,
                     [("mread1", 30, parent_seq, None)])
        father_bam = os.path.join(tmpdir, "father.bam")
        _create_bam(father_bam, ref_fa, chrom,
                     [("fread1", 30, parent_seq, None)])

        out_prefix = os.path.join(tmpdir, "sv_snp")
        args = parse_args([
            "--child", child_bam, "--mother", mother_bam,
            "--father", father_bam, "--ref-fasta", ref_fa,
            "--out-prefix", out_prefix, "--min-child-count", "3",
            "--kmer-size", "5",
        ])
        run_discovery_pipeline(args)

        # Check BED: class should be SMALL
        bed_path = f"{out_prefix}.bed"
        with open(bed_path) as fh:
            bed_lines = [l.strip() for l in fh if l.strip() and not l.startswith("#")]
        assert len(bed_lines) >= 1
        parts = bed_lines[0].split("\t")
        assert parts[5] == "0"   # split_reads
        assert parts[9] == "SMALL"

        # Check metrics: class should be SMALL
        metrics_path = f"{out_prefix}.metrics.json"
        with open(metrics_path) as fh:
            metrics = json.load(fh)
        for region in metrics["regions"]:
            assert region["split_reads"] == 0
            assert region["class"] == "SMALL"

        # BEDPE should be empty (header only)
        bedpe_path = f"{out_prefix}.sv.bedpe"
        assert os.path.exists(bedpe_path)
        with open(bedpe_path) as fh:
            lines = [l for l in fh if not l.startswith("#")]
        assert len(lines) == 0

    def test_split_read_sv_detection(self, tmpdir):
        """A read with SA tag bridging two regions should produce SV
        classification and BEDPE entry."""
        chrom = "chr1"
        ref_fa = os.path.join(tmpdir, "ref.fa")
        ref_seq = _create_ref_fasta(ref_fa, chrom, 200)

        # Create two different mutation sequences at distant positions
        # Mutation 1: position 50 (read at 30-90)
        child_seq1 = list(ref_seq[30:90])
        child_seq1[20] = "G" if ref_seq[50] != "G" else "T"
        child_seq1 = "".join(child_seq1)

        # Mutation 2: position 150 (read at 130-190)
        child_seq2 = list(ref_seq[130:190])
        child_seq2[20] = "G" if ref_seq[150] != "G" else "T"
        child_seq2 = "".join(child_seq2)

        # SA tag pointing from first region to second (1-based pos)
        sa_tag = f"{chrom},131,+,30M30S,60,0;"

        child_bam = os.path.join(tmpdir, "child.bam")
        _create_bam_with_supplementary(
            child_bam, ref_fa, [chrom], [300],
            [
                # Primary reads with mutation 1, with SA tags
                {"name": "sv_read1", "pos": 30, "seq": child_seq1,
                 "flag": 0, "sa_tag": sa_tag},
                {"name": "sv_read2", "pos": 30, "seq": child_seq1,
                 "flag": 0, "sa_tag": sa_tag},
                {"name": "read3", "pos": 30, "seq": child_seq1, "flag": 0},
                {"name": "read4", "pos": 30, "seq": child_seq1, "flag": 0},
                # Supplementary reads with mutation 2 at position 130
                {"name": "sv_read1", "pos": 130, "seq": child_seq2,
                 "flag": 0x800},
                {"name": "sv_read2", "pos": 130, "seq": child_seq2,
                 "flag": 0x800},
                # Primary reads with mutation 2
                {"name": "read5", "pos": 130, "seq": child_seq2, "flag": 0},
                {"name": "read6", "pos": 130, "seq": child_seq2, "flag": 0},
            ],
        )

        parent_seq1 = ref_seq[30:90]
        parent_seq2 = ref_seq[130:190]
        mother_bam = os.path.join(tmpdir, "mother.bam")
        _create_bam(mother_bam, ref_fa, chrom, [
            ("mread1", 30, parent_seq1, None),
            ("mread2", 130, parent_seq2, None),
        ])
        father_bam = os.path.join(tmpdir, "father.bam")
        _create_bam(father_bam, ref_fa, chrom, [
            ("fread1", 30, parent_seq1, None),
            ("fread2", 130, parent_seq2, None),
        ])

        out_prefix = os.path.join(tmpdir, "sv_split")
        args = parse_args([
            "--child", child_bam, "--mother", mother_bam,
            "--father", father_bam, "--ref-fasta", ref_fa,
            "--out-prefix", out_prefix, "--min-child-count", "3",
            "--kmer-size", "5", "--cluster-distance", "0",
        ])
        run_discovery_pipeline(args)

        # Check that split_reads > 0 in at least one region
        metrics_path = f"{out_prefix}.metrics.json"
        with open(metrics_path) as fh:
            metrics = json.load(fh)
        assert metrics["candidate_regions"] >= 1
        has_split = any(r["split_reads"] > 0 for r in metrics["regions"])
        assert has_split, "Expected at least one region with split_reads > 0"
        # Region with split reads should be classified SV
        for region in metrics["regions"]:
            if region["split_reads"] >= 2:
                assert region["class"] == "SV"

    def test_bedpe_written(self, tmpdir):
        """BEDPE output file should always be written."""
        chrom = "chr1"
        ref_fa = os.path.join(tmpdir, "ref.fa")
        ref_seq = _create_ref_fasta(ref_fa, chrom, 200)

        child_seq = list(ref_seq[30:90])
        child_seq[20] = "G" if ref_seq[50] != "G" else "T"
        child_seq = "".join(child_seq)

        child_bam = os.path.join(tmpdir, "child.bam")
        _create_bam(child_bam, ref_fa, chrom, [
            ("read1", 30, child_seq, None),
            ("read2", 30, child_seq, None),
            ("read3", 30, child_seq, None),
            ("read4", 30, child_seq, None),
        ])

        parent_seq = ref_seq[30:90]
        mother_bam = os.path.join(tmpdir, "mother.bam")
        _create_bam(mother_bam, ref_fa, chrom,
                     [("mread1", 30, parent_seq, None)])
        father_bam = os.path.join(tmpdir, "father.bam")
        _create_bam(father_bam, ref_fa, chrom,
                     [("fread1", 30, parent_seq, None)])

        out_prefix = os.path.join(tmpdir, "sv_bedpe")
        args = parse_args([
            "--child", child_bam, "--mother", mother_bam,
            "--father", father_bam, "--ref-fasta", ref_fa,
            "--out-prefix", out_prefix, "--min-child-count", "3",
            "--kmer-size", "5",
        ])
        run_discovery_pipeline(args)

        bedpe_path = f"{out_prefix}.sv.bedpe"
        assert os.path.exists(bedpe_path)
        with open(bedpe_path) as fh:
            header = fh.readline()
        assert header.startswith("#chrom1")

    def test_classify_regions_unit(self):
        """Unit test for _classify_regions logic."""
        regions = [
            ("chr1", 100, 200),
            ("chr1", 500, 600),
            ("chr2", 100, 200),
            ("chr3", 100, 200),
            ("chr4", 100, 200),
        ]
        annotations = {
            ("chr1", 100, 200): {"split_reads": 3, "discordant_pairs": 0,
                                 "max_clip_len": 0, "unmapped_mates": 0},
            ("chr1", 500, 600): {"split_reads": 0, "discordant_pairs": 0,
                                 "max_clip_len": 0, "unmapped_mates": 0},
            ("chr2", 100, 200): {"split_reads": 1, "discordant_pairs": 0,
                                 "max_clip_len": 0, "unmapped_mates": 0},
            # High discordant pairs → SV even with 0 split reads
            ("chr3", 100, 200): {"split_reads": 0, "discordant_pairs": 3,
                                 "max_clip_len": 0, "unmapped_mates": 0},
            # High unmapped mates → SV even with 0 split reads
            ("chr4", 100, 200): {"split_reads": 0, "discordant_pairs": 0,
                                 "max_clip_len": 0, "unmapped_mates": 2},
        }
        # Region A is linked to region B
        sv_links = [
            {"region_a": ("chr1", 100, 200), "region_b": ("chr1", 500, 600),
             "supporting_reads": {"r1"}, "sv_type_hint": "INTRA"},
        ]
        _classify_regions(regions, annotations, sv_links)

        # Region with >= 2 split reads → SV
        assert annotations[("chr1", 100, 200)]["class"] == "SV"
        # Region linked via sv_links → SV (even with 0 split reads)
        assert annotations[("chr1", 500, 600)]["class"] == "SV"
        # Region with 1 split read but not linked → AMBIGUOUS
        assert annotations[("chr2", 100, 200)]["class"] == "AMBIGUOUS"
        # Region with >= 2 discordant pairs → SV
        assert annotations[("chr3", 100, 200)]["class"] == "SV"
        # Region with >= 2 unmapped mates → SV
        assert annotations[("chr4", 100, 200)]["class"] == "SV"

    def test_write_bedpe_format(self, tmpdir):
        """Unit test for _write_bedpe output format."""
        links = [
            {"region_a": ("chr1", 100, 200), "region_b": ("chr1", 5000, 5100),
             "supporting_reads": {"r1", "r2"}, "sv_type_hint": "INTRA"},
            {"region_a": ("chr1", 100, 200), "region_b": ("chr2", 100, 200),
             "supporting_reads": {"r3"}, "sv_type_hint": "BND"},
        ]
        bedpe_path = os.path.join(tmpdir, "test.bedpe")
        _write_bedpe(links, bedpe_path)

        with open(bedpe_path) as fh:
            lines = fh.readlines()

        assert lines[0].startswith("#chrom1")
        parts1 = lines[1].strip().split("\t")
        assert len(parts1) == 9
        assert parts1[0] == "chr1"
        assert parts1[6] == "SV_1"
        assert parts1[7] == "2"  # 2 supporting reads
        assert parts1[8] == "INTRA"

        parts2 = lines[2].strip().split("\t")
        assert parts2[3] == "chr2"
        assert parts2[8] == "BND"

    def test_sv_bedpe_cli_argument(self, tmpdir):
        """--sv-bedpe argument should be accepted."""
        chrom = "chr1"
        ref_fa = os.path.join(tmpdir, "ref.fa")
        ref_seq = _create_ref_fasta(ref_fa, chrom, 200)

        child_seq = list(ref_seq[30:90])
        child_seq[20] = "G" if ref_seq[50] != "G" else "T"
        child_seq = "".join(child_seq)

        child_bam = os.path.join(tmpdir, "child.bam")
        _create_bam(child_bam, ref_fa, chrom, [
            ("read1", 30, child_seq, None),
            ("read2", 30, child_seq, None),
            ("read3", 30, child_seq, None),
            ("read4", 30, child_seq, None),
        ])

        parent_seq = ref_seq[30:90]
        mother_bam = os.path.join(tmpdir, "mother.bam")
        _create_bam(mother_bam, ref_fa, chrom,
                     [("mread1", 30, parent_seq, None)])
        father_bam = os.path.join(tmpdir, "father.bam")
        _create_bam(father_bam, ref_fa, chrom,
                     [("fread1", 30, parent_seq, None)])

        custom_bedpe = os.path.join(tmpdir, "custom.bedpe")
        out_prefix = os.path.join(tmpdir, "sv_cli")
        args = parse_args([
            "--child", child_bam, "--mother", mother_bam,
            "--father", father_bam, "--ref-fasta", ref_fa,
            "--out-prefix", out_prefix, "--min-child-count", "3",
            "--kmer-size", "5",
            "--sv-bedpe", custom_bedpe,
        ])
        run_discovery_pipeline(args)
        assert os.path.exists(custom_bedpe)


# ── bedGraph output tests ──────────────────────────────────────────────





# ── bedGraph output tests ──────────────────────────────────────────────


class TestWriteBedgraph:
    """Unit tests for _write_bedgraph."""

    @staticmethod
    def _data_lines(path):
        """Read non-header data lines from a bedgraph file."""
        with open(path) as fh:
            return [
                l for l in fh.read().strip().split("\n")
                if l and not l.startswith("#") and not l.startswith("track")
            ]

    def test_basic_bedgraph(self, tmp_path):
        """Adjacent positions with the same value are merged."""
        cov = {
            "chr1": collections.Counter({10: 2, 11: 2, 12: 2, 20: 1}),
        }
        out = str(tmp_path / "out.bedgraph")
        _write_bedgraph(cov, out)
        lines = self._data_lines(out)
        assert lines == ["chr1\t10\t13\t2", "chr1\t20\t21\t1"]

    def test_different_values_not_merged(self, tmp_path):
        """Adjacent positions with different values stay separate."""
        cov = {
            "chr1": collections.Counter({5: 3, 6: 1}),
        }
        out = str(tmp_path / "out.bedgraph")
        _write_bedgraph(cov, out)
        lines = self._data_lines(out)
        assert lines == ["chr1\t5\t6\t3", "chr1\t6\t7\t1"]

    def test_multi_chrom_sorted(self, tmp_path):
        """Chromosomes are written in sorted order."""
        cov = {
            "chr2": collections.Counter({100: 1}),
            "chr1": collections.Counter({50: 2}),
        }
        out = str(tmp_path / "out.bedgraph")
        _write_bedgraph(cov, out)
        lines = self._data_lines(out)
        assert lines[0].startswith("chr1")
        assert lines[1].startswith("chr2")

    def test_empty_coverage(self, tmp_path):
        """Empty coverage produces only the header."""
        out = str(tmp_path / "out.bedgraph")
        _write_bedgraph({}, out)
        lines = self._data_lines(out)
        assert lines == []

    def test_single_position(self, tmp_path):
        """Single position produces a 1-bp interval."""
        cov = {"chrX": collections.Counter({42: 5})}
        out = str(tmp_path / "out.bedgraph")
        _write_bedgraph(cov, out)
        lines = self._data_lines(out)
        assert lines == ["chrX\t42\t43\t5"]

    def test_header_present(self, tmp_path):
        """bedGraph file should have a descriptive header line."""
        cov = {"chr1": collections.Counter({10: 1})}
        out = str(tmp_path / "out.bedgraph")
        _write_bedgraph(cov, out)
        with open(out) as fh:
            first_line = fh.readline()
        assert first_line.startswith("#track type=bedGraph")
        assert "min_reads" in first_line

    def test_min_reads_filter(self, tmp_path):
        """Positions below min_reads threshold are filtered out."""
        cov = {
            "chr1": collections.Counter({10: 5, 20: 2, 30: 8}),
        }
        rc = {
            "chr1": collections.Counter({10: 4, 20: 1, 30: 3}),
        }
        out = str(tmp_path / "out.bedgraph")
        _write_bedgraph(cov, out, read_coverage=rc, min_reads=3)
        lines = self._data_lines(out)
        # Only pos 10 (4 reads) and 30 (3 reads) pass min_reads=3
        assert len(lines) == 2
        assert "chr1\t10\t11\t5" in lines
        assert "chr1\t30\t31\t8" in lines





class TestCollectKmerRefPositions:
    """Unit tests for _collect_kmer_ref_positions."""

    def test_simple_alignment(self, tmp_path):
        """All-match alignment maps k-mer positions directly."""
        bam_path = str(tmp_path / "test.bam")
        header = pysam.AlignmentHeader.from_references(["chr1"], [300])
        with pysam.AlignmentFile(bam_path, "wb", header=header) as bam:
            seg = pysam.AlignedSegment()
            seg.query_name = "read1"
            seg.query_sequence = "ACGTACGT"
            seg.flag = 0
            seg.reference_id = 0
            seg.reference_start = 100
            seg.mapping_quality = 60
            seg.cigar = [(0, 8)]  # 8M
            seg.query_qualities = pysam.qualitystring_to_array("I" * 8)
            bam.write(seg)

        with pysam.AlignmentFile(bam_path, "rb",
                                 check_sq=False) as bam:
            read = next(bam.fetch(until_eof=True))
            # k-mer of size 3 starting at query pos 2 → covers qpos 2,3,4
            cov = _collect_kmer_ref_positions(read, {2}, 3)
            # ref positions should be 102, 103, 104
            assert cov[102] == 1
            assert cov[103] == 1
            assert cov[104] == 1
            assert len(cov) == 3

    def test_insertion_skips_ref(self, tmp_path):
        """Query positions in insertions have no ref mapping."""
        bam_path = str(tmp_path / "test.bam")
        header = pysam.AlignmentHeader.from_references(["chr1"], [300])
        with pysam.AlignmentFile(bam_path, "wb", header=header) as bam:
            seg = pysam.AlignedSegment()
            seg.query_name = "read1"
            seg.query_sequence = "ACGTACGT"
            seg.flag = 0
            seg.reference_id = 0
            seg.reference_start = 100
            seg.mapping_quality = 60
            # 3M 2I 3M → qpos 3,4 are insertions (no ref pos)
            seg.cigar = [(0, 3), (1, 2), (0, 3)]
            seg.query_qualities = pysam.qualitystring_to_array("I" * 8)
            bam.write(seg)

        with pysam.AlignmentFile(bam_path, "rb",
                                 check_sq=False) as bam:
            read = next(bam.fetch(until_eof=True))
            # k-mer of size 4 starting at query pos 2 → qpos 2,3,4,5
            # qpos 2→ref102, 3→None(ins), 4→None(ins), 5→ref103
            cov = _collect_kmer_ref_positions(read, {2}, 4)
            assert cov[102] == 1
            assert cov[103] == 1
            assert 3 not in cov  # ins positions skipped
            assert len(cov) == 2





class TestBedgraphDiscoveryIntegration:
    """Integration test: bedGraph is produced by discovery pipeline."""

    def test_bedgraph_written(self, tmp_path):
        chrom = "chr1"
        ref_fa = str(tmp_path / "ref.fa")
        ref_seq = _create_ref_fasta(ref_fa, chrom, length=200)

        # Child has a de novo mutation at position 50
        child_seq = list(ref_seq[30:90])
        child_seq[20] = "G" if ref_seq[50] != "G" else "T"
        child_seq = "".join(child_seq)

        child_bam = str(tmp_path / "child.bam")
        _create_bam(child_bam, ref_fa, chrom, [
            ("read1", 30, child_seq, None),
            ("read2", 30, child_seq, None),
            ("read3", 30, child_seq, None),
            ("read4", 30, child_seq, None),
        ])

        parent_seq = ref_seq[30:90]
        mother_bam = str(tmp_path / "mother.bam")
        _create_bam(mother_bam, ref_fa, chrom, [
            ("mread1", 30, parent_seq, None),
            ("mread2", 30, parent_seq, None),
            ("mread3", 30, parent_seq, None),
        ])

        father_bam = str(tmp_path / "father.bam")
        _create_bam(father_bam, ref_fa, chrom, [
            ("fread1", 30, parent_seq, None),
            ("fread2", 30, parent_seq, None),
            ("fread3", 30, parent_seq, None),
        ])

        out_prefix = str(tmp_path / "disco")
        args = parse_args([
            "--child", child_bam, "--mother", mother_bam,
            "--father", father_bam, "--ref-fasta", ref_fa,
            "--out-prefix", out_prefix, "--min-child-count", "3",
            "--kmer-size", "5",
        ])
        run_discovery_pipeline(args)

        bedgraph = f"{out_prefix}.kmer_coverage.bedgraph"
        assert os.path.exists(bedgraph), "bedGraph file not created"
        with open(bedgraph) as fh:
            content = fh.read()
        assert content.strip(), "bedGraph should have content"
        for line in content.strip().split("\n"):
            if line.startswith("#") or line.startswith("track"):
                continue
            cols = line.split("\t")
            assert len(cols) == 4, f"Expected 4 columns, got {len(cols)}"
            assert cols[0] == chrom
            assert int(cols[1]) < int(cols[2])  # start < end
            assert int(cols[3]) > 0  # positive coverage





class TestJellyfishBatchScanMemory:
    """Unit tests for Jellyfish-backed batch scanning memory behavior."""

    def test_scan_contig_clears_cache_per_batch(self, monkeypatch):
        class _FakeRead:
            def __init__(self, seq):
                self.query_sequence = seq
                self.is_secondary = False
                self.is_duplicate = False

        class _FakeBam:
            def __init__(self, *args, **kwargs):
                self._reads = [
                    _FakeRead("AAAAA"),
                    _FakeRead("CCCCC"),
                    _FakeRead("GGGGG"),
                ]

            def fetch(self, contig=None):
                return iter(self._reads)

            def close(self):
                return None

        class _FakeJFQuery:
            def __init__(self):
                self.query_calls = []
                self.close_calls = 0

            def query_batch(self, kmers):
                self.query_calls.append(set(kmers))
                return set(kmers)

            def close(self):
                self.close_calls += 1

        jf_query = _FakeJFQuery()
        seen_reads = []

        def _fake_extract_read_kmers(seq, kmer_size):
            return ({0: seq}, [seq])

        def _fake_process_informative_read(
            read,
            unique_in_read,
            kmer_hit_indices,
            kmer_size,
            reads_seen,
            read_hits,
            read_sv_meta,
            kmer_coverage,
            read_coverage,
        ):
            seen_reads.append((read.query_sequence, unique_in_read))
            return 0

        class _FakePysam:
            AlignmentFile = _FakeBam

        monkeypatch.setattr(bam_scanner_mod, "pysam", _FakePysam)
        monkeypatch.setattr(
            bam_scanner_mod, "_extract_read_kmers", _fake_extract_read_kmers,
        )
        monkeypatch.setattr(
            bam_scanner_mod, "_process_informative_read",
            _fake_process_informative_read,
        )
        monkeypatch.setattr(bam_scanner_mod, "_JF_READ_BATCH_SIZE", 2)
        monkeypatch.setattr(bam_scanner_mod, "_worker_automaton", None)
        monkeypatch.setattr(bam_scanner_mod, "_worker_jf_query", jf_query)
        monkeypatch.setattr(bam_scanner_mod, "_worker_kmer_size", 5)
        monkeypatch.setattr(
            bam_scanner_mod, "_worker_min_distinct_kmers_per_read", 1,
        )

        (
            _read_hits,
            _reads_seen,
            _unmapped_informative,
            total_reads_scanned,
            _read_sv_meta,
            _kmer_coverage,
            _read_coverage,
        ) = bam_scanner_mod._scan_contig_for_hits("child.bam", "ref.fa", "chr1")

        assert total_reads_scanned == 3
        assert len(jf_query.query_calls) == 2
        assert jf_query.query_calls[0] == {"AAAAA", "CCCCC"}
        assert jf_query.query_calls[1] == {"GGGGG"}
        assert jf_query.close_calls == 2
        assert len(seen_reads) == 3
