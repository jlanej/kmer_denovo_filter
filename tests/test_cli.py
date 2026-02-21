"""Tests for CLI argument parsing."""

import pytest

from kmer_denovo_filter.cli import parse_args


class TestParseArgs:
    REQUIRED_ARGS = [
        "--child", "child.bam",
        "--mother", "mother.bam",
        "--father", "father.bam",
        "--vcf", "input.vcf",
        "--output", "output.vcf",
    ]

    def test_required_args(self):
        args = parse_args(self.REQUIRED_ARGS)
        assert args.child == "child.bam"
        assert args.mother == "mother.bam"
        assert args.father == "father.bam"
        assert args.ref_fasta is None
        assert args.vcf == "input.vcf"
        assert args.output == "output.vcf"

    def test_defaults(self):
        args = parse_args(self.REQUIRED_ARGS)
        assert args.kmer_size == 31
        assert args.min_baseq == 20
        assert args.min_mapq == 20
        assert args.threads == 4
        assert args.debug_kmers is False
        assert args.metrics is None
        assert args.informative_reads is None
        assert args.summary is None
        assert args.ref_fasta is None

    def test_custom_kmer_size(self):
        args = parse_args(self.REQUIRED_ARGS + ["--kmer-size", "25"])
        assert args.kmer_size == 25

    def test_short_kmer_size(self):
        args = parse_args(self.REQUIRED_ARGS + ["-k", "21"])
        assert args.kmer_size == 21

    def test_short_ref_fasta(self):
        extra = [
            "--child", "c.bam", "--mother", "m.bam", "--father", "f.bam",
            "-r", "ref.fa", "--vcf", "in.vcf", "-o", "out.vcf",
        ]
        args = parse_args(extra)
        assert args.ref_fasta == "ref.fa"
        assert args.output == "out.vcf"

    def test_metrics_flag(self):
        args = parse_args(self.REQUIRED_ARGS + ["--metrics", "metrics.json"])
        assert args.metrics == "metrics.json"

    def test_debug_kmers(self):
        args = parse_args(self.REQUIRED_ARGS + ["--debug-kmers"])
        assert args.debug_kmers is True

    def test_min_baseq(self):
        args = parse_args(self.REQUIRED_ARGS + ["--min-baseq", "30"])
        assert args.min_baseq == 30

    def test_min_mapq(self):
        args = parse_args(self.REQUIRED_ARGS + ["--min-mapq", "10"])
        assert args.min_mapq == 10

    def test_threads(self):
        args = parse_args(self.REQUIRED_ARGS + ["--threads", "8"])
        assert args.threads == 8

    def test_short_threads(self):
        args = parse_args(self.REQUIRED_ARGS + ["-t", "2"])
        assert args.threads == 2

    def test_informative_reads(self):
        args = parse_args(
            self.REQUIRED_ARGS + ["--informative-reads", "info.bam"]
        )
        assert args.informative_reads == "info.bam"

    def test_summary_flag(self):
        args = parse_args(self.REQUIRED_ARGS + ["--summary", "summary.txt"])
        assert args.summary == "summary.txt"

    def test_ref_fasta_optional(self):
        args = parse_args(self.REQUIRED_ARGS)
        assert args.ref_fasta is None

    def test_ref_fasta_provided(self):
        args = parse_args(self.REQUIRED_ARGS + ["--ref-fasta", "ref.fa"])
        assert args.ref_fasta == "ref.fa"

    def test_missing_required(self):
        with pytest.raises(SystemExit):
            parse_args(["--child", "child.bam"])
