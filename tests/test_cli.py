"""Tests for CLI argument parsing."""

import pytest

from kmer_denovo_filter.cli import parse_args, parse_vcf_args, parse_discovery_args


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
        assert args.memory is None
        assert args.metrics is None
        assert args.informative_reads is None
        assert args.summary is None
        assert args.ref_fasta is None
        assert args.ref_jf is None
        assert args.min_child_count == 3
        assert args.out_prefix is None
        assert args.kraken2_memory_mapping is False

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

    def test_memory_default(self):
        args = parse_args(self.REQUIRED_ARGS)
        assert args.memory is None

    def test_memory_custom(self):
        args = parse_args(self.REQUIRED_ARGS + ["--memory", "64"])
        assert args.memory == 64.0

    def test_memory_fractional(self):
        args = parse_args(self.REQUIRED_ARGS + ["--memory", "128.5"])
        assert args.memory == 128.5

    def test_informative_reads(self):
        args = parse_args(
            self.REQUIRED_ARGS + ["--informative-reads", "info.bam"]
        )
        assert args.informative_reads == "info.bam"

    def test_summary_flag(self):
        args = parse_args(self.REQUIRED_ARGS + ["--summary", "summary.txt"])
        assert args.summary == "summary.txt"

    def test_proband_id(self):
        args = parse_args(self.REQUIRED_ARGS + ["--proband-id", "HG002"])
        assert args.proband_id == "HG002"

    def test_proband_id_default(self):
        args = parse_args(self.REQUIRED_ARGS)
        assert args.proband_id is None

    def test_ref_fasta_optional(self):
        args = parse_args(self.REQUIRED_ARGS)
        assert args.ref_fasta is None

    def test_ref_fasta_provided(self):
        args = parse_args(self.REQUIRED_ARGS + ["--ref-fasta", "ref.fa"])
        assert args.ref_fasta == "ref.fa"

    def test_missing_required(self):
        with pytest.raises(SystemExit):
            parse_args(["--child", "child.bam"])

    def test_vcf_optional(self):
        """--vcf is no longer required at the parser level."""
        args = parse_args([
            "--child", "child.bam",
            "--mother", "mother.bam",
            "--father", "father.bam",
        ])
        assert args.vcf is None
        assert args.output is None

    def test_discovery_mode_args(self):
        """Discovery mode arguments are parsed correctly."""
        args = parse_args([
            "--child", "child.bam",
            "--mother", "mother.bam",
            "--father", "father.bam",
            "--ref-fasta", "ref.fa",
            "--out-prefix", "trio1",
            "--min-child-count", "5",
            "--ref-jf", "ref.k31.jf",
        ])
        assert args.vcf is None
        assert args.out_prefix == "trio1"
        assert args.min_child_count == 5
        assert args.ref_jf == "ref.k31.jf"

    def test_ref_jf_default(self):
        args = parse_args(self.REQUIRED_ARGS)
        assert args.ref_jf is None

    def test_min_child_count_default(self):
        args = parse_args(self.REQUIRED_ARGS)
        assert args.min_child_count == 3

    def test_out_prefix(self):
        args = parse_args(self.REQUIRED_ARGS + ["--out-prefix", "myprefix"])
        assert args.out_prefix == "myprefix"

    def test_cluster_distance_default(self):
        args = parse_args(self.REQUIRED_ARGS)
        assert args.cluster_distance == 500

    def test_cluster_distance_custom(self):
        args = parse_args(self.REQUIRED_ARGS + ["--cluster-distance", "1000"])
        assert args.cluster_distance == 1000

    def test_min_supporting_reads_default(self):
        args = parse_args(self.REQUIRED_ARGS)
        assert args.min_supporting_reads == 1

    def test_min_supporting_reads_custom(self):
        args = parse_args(
            self.REQUIRED_ARGS + ["--min-supporting-reads", "3"]
        )
        assert args.min_supporting_reads == 3

    def test_min_distinct_kmers_default(self):
        args = parse_args(self.REQUIRED_ARGS)
        assert args.min_distinct_kmers == 1

    def test_min_distinct_kmers_custom(self):
        args = parse_args(
            self.REQUIRED_ARGS + ["--min-distinct-kmers", "5"]
        )
        assert args.min_distinct_kmers == 5

    def test_parent_max_count_default(self):
        args = parse_args(self.REQUIRED_ARGS)
        assert args.parent_max_count == 0

    def test_parent_max_count_custom(self):
        args = parse_args(self.REQUIRED_ARGS + ["--parent-max-count", "2"])
        assert args.parent_max_count == 2

    def test_discovery_mode_new_args(self):
        """Discovery mode new arguments are parsed correctly."""
        args = parse_args([
            "--child", "child.bam",
            "--mother", "mother.bam",
            "--father", "father.bam",
            "--ref-fasta", "ref.fa",
            "--out-prefix", "trio1",
            "--cluster-distance", "250",
            "--min-supporting-reads", "4",
            "--min-distinct-kmers", "3",
            "--parent-max-count", "1",
        ])
        assert args.cluster_distance == 250
        assert args.min_supporting_reads == 4
        assert args.min_distinct_kmers == 3
        assert args.parent_max_count == 1

    def test_min_distinct_kmers_per_read_default(self):
        args = parse_args(self.REQUIRED_ARGS)
        assert args.min_distinct_kmers_per_read is None

    def test_min_distinct_kmers_per_read_custom(self):
        args = parse_args(
            self.REQUIRED_ARGS + ["--min-distinct-kmers-per-read", "5"]
        )
        assert args.min_distinct_kmers_per_read == 5

    def test_tmp_dir_default(self):
        args = parse_args(self.REQUIRED_ARGS)
        assert args.tmp_dir is None

    def test_tmp_dir_custom(self):
        args = parse_args(self.REQUIRED_ARGS + ["--tmp-dir", "/scratch/tmp"])
        assert args.tmp_dir == "/scratch/tmp"

    def test_kraken2_memory_mapping_flag(self):
        args = parse_args(self.REQUIRED_ARGS + ["--kraken2-memory-mapping"])
        assert args.kraken2_memory_mapping is True

    def test_kraken2_read_detail_default(self):
        args = parse_args(self.REQUIRED_ARGS)
        assert args.kraken2_read_detail is None

    def test_kraken2_read_detail_custom(self):
        args = parse_args(
            self.REQUIRED_ARGS + ["--kraken2-read-detail", "/tmp/detail.bed.gz"],
        )
        assert args.kraken2_read_detail == "/tmp/detail.bed.gz"


class TestParseVcfArgs:
    """Tests for the VCF-mode parser (kmer-denovo)."""

    REQUIRED_ARGS = [
        "--child", "child.bam",
        "--mother", "mother.bam",
        "--father", "father.bam",
        "--vcf", "input.vcf",
        "--output", "output.vcf",
    ]

    def test_required_args(self):
        args = parse_vcf_args(self.REQUIRED_ARGS)
        assert args.child == "child.bam"
        assert args.mother == "mother.bam"
        assert args.father == "father.bam"
        assert args.vcf == "input.vcf"
        assert args.output == "output.vcf"

    def test_vcf_is_required(self):
        """--vcf is required in the VCF parser."""
        with pytest.raises(SystemExit):
            parse_vcf_args([
                "--child", "child.bam",
                "--mother", "mother.bam",
                "--father", "father.bam",
                "--output", "output.vcf",
            ])

    def test_output_is_required(self):
        """--output is required in the VCF parser."""
        with pytest.raises(SystemExit):
            parse_vcf_args([
                "--child", "child.bam",
                "--mother", "mother.bam",
                "--father", "father.bam",
                "--vcf", "input.vcf",
            ])

    def test_defaults(self):
        args = parse_vcf_args(self.REQUIRED_ARGS)
        assert args.kmer_size == 31
        assert args.min_baseq == 20
        assert args.min_mapq == 20
        assert args.threads == 4
        assert args.debug_kmers is False
        assert args.memory is None
        assert args.metrics is None
        assert args.informative_reads is None
        assert args.summary is None
        assert args.ref_fasta is None
        assert args.kraken2_memory_mapping is False
        assert args.proband_id is None

    def test_no_discovery_flags(self):
        """VCF parser must not accept discovery-only flags."""
        for flag in ("--out-prefix", "--ref-jf", "--min-child-count",
                     "--cluster-distance", "--min-supporting-reads",
                     "--min-distinct-kmers", "--min-bedgraph-reads",
                     "--min-distinct-kmers-per-read", "--parent-max-count",
                     "--sv-bedpe", "--candidate-summary"):
            with pytest.raises(SystemExit):
                parse_vcf_args(self.REQUIRED_ARGS + [flag, "value"])

    def test_vcf_accepts_kraken2_flags(self):
        args = parse_vcf_args(
            self.REQUIRED_ARGS + [
                "--kraken2-db", "/db",
                "--kraken2-confidence", "0.5",
                "--kraken2-memory-mapping",
                "--kraken2-read-detail", "/tmp/detail.bed.gz",
                "--kraken2-span-bed", "/tmp/span.bed.gz",
                "--no-expanded-bed",
            ]
        )
        assert args.kraken2_db == "/db"
        assert args.kraken2_confidence == 0.5
        assert args.kraken2_memory_mapping is True
        assert args.kraken2_read_detail == "/tmp/detail.bed.gz"
        assert args.kraken2_span_bed == "/tmp/span.bed.gz"
        assert args.no_expanded_bed is True

    def test_shared_args(self):
        args = parse_vcf_args(self.REQUIRED_ARGS + [
            "--kmer-size", "25",
            "--min-baseq", "30",
            "--threads", "8",
            "--memory", "64",
            "--debug-kmers",
            "--jf-hash-size", "2G",
            "--tmp-dir", "/scratch",
            "--ref-fasta", "ref.fa",
        ])
        assert args.kmer_size == 25
        assert args.min_baseq == 30
        assert args.threads == 8
        assert args.memory == 64.0
        assert args.debug_kmers is True
        assert args.jf_hash_size == "2G"
        assert args.tmp_dir == "/scratch"
        assert args.ref_fasta == "ref.fa"


class TestParseDiscoveryArgs:
    """Tests for the discovery-mode parser (kmer-discovery)."""

    REQUIRED_ARGS = [
        "--child", "child.bam",
        "--mother", "mother.bam",
        "--father", "father.bam",
        "--out-prefix", "trio1",
    ]

    def test_required_args(self):
        args = parse_discovery_args(self.REQUIRED_ARGS)
        assert args.child == "child.bam"
        assert args.mother == "mother.bam"
        assert args.father == "father.bam"
        assert args.out_prefix == "trio1"

    def test_out_prefix_is_required(self):
        """--out-prefix is required in the discovery parser."""
        with pytest.raises(SystemExit):
            parse_discovery_args([
                "--child", "child.bam",
                "--mother", "mother.bam",
                "--father", "father.bam",
            ])

    def test_defaults(self):
        args = parse_discovery_args(self.REQUIRED_ARGS)
        assert args.kmer_size == 31
        assert args.min_baseq == 20
        assert args.threads == 4
        assert args.debug_kmers is False
        assert args.memory is None
        assert args.ref_fasta is None
        assert args.ref_jf is None
        assert args.min_child_count == 3
        assert args.cluster_distance == 500
        assert args.min_supporting_reads == 1
        assert args.min_distinct_kmers == 1
        assert args.min_bedgraph_reads == 3
        assert args.min_distinct_kmers_per_read is None
        assert args.parent_max_count == 0
        assert args.sv_bedpe is None
        assert args.candidate_summary is None

    def test_no_vcf_flags(self):
        """Discovery parser must not accept VCF-only flags."""
        for flag in ("--vcf", "--output", "--min-mapq", "--proband-id",
                     "--informative-reads", "--metrics", "--summary",
                     "--kraken2-db", "--kraken2-confidence",
                     "--kraken2-memory-mapping", "--kraken2-read-detail",
                     "--kraken2-span-bed", "--no-expanded-bed"):
            with pytest.raises(SystemExit):
                parse_discovery_args(self.REQUIRED_ARGS + [flag, "value"])

    def test_discovery_specific_args(self):
        args = parse_discovery_args(self.REQUIRED_ARGS + [
            "--ref-fasta", "ref.fa",
            "--ref-jf", "ref.k31.jf",
            "--min-child-count", "5",
            "--cluster-distance", "250",
            "--min-supporting-reads", "4",
            "--min-distinct-kmers", "3",
            "--min-bedgraph-reads", "5",
            "--min-distinct-kmers-per-read", "7",
            "--parent-max-count", "1",
            "--sv-bedpe", "out.bedpe",
            "--candidate-summary", "summary.txt",
        ])
        assert args.ref_fasta == "ref.fa"
        assert args.ref_jf == "ref.k31.jf"
        assert args.min_child_count == 5
        assert args.cluster_distance == 250
        assert args.min_supporting_reads == 4
        assert args.min_distinct_kmers == 3
        assert args.min_bedgraph_reads == 5
        assert args.min_distinct_kmers_per_read == 7
        assert args.parent_max_count == 1
        assert args.sv_bedpe == "out.bedpe"
        assert args.candidate_summary == "summary.txt"

    def test_shared_args(self):
        args = parse_discovery_args(self.REQUIRED_ARGS + [
            "--kmer-size", "25",
            "--min-baseq", "30",
            "--threads", "8",
            "--memory", "64",
            "--debug-kmers",
            "--jf-hash-size", "2G",
            "--tmp-dir", "/scratch",
        ])
        assert args.kmer_size == 25
        assert args.min_baseq == 30
        assert args.threads == 8
        assert args.memory == 64.0
        assert args.debug_kmers is True
        assert args.jf_hash_size == "2G"
        assert args.tmp_dir == "/scratch"


class TestEntryPoints:
    """Tests that entry point functions are importable and callable."""

    def test_vcf_main_importable(self):
        from kmer_denovo_filter.cli import vcf_main
        assert callable(vcf_main)

    def test_discovery_main_importable(self):
        from kmer_denovo_filter.cli import discovery_main
        assert callable(discovery_main)

    def test_main_still_importable(self):
        from kmer_denovo_filter.cli import main
        assert callable(main)

    def test_vcf_main_missing_args_exits(self):
        """vcf_main exits when required args are missing."""
        from kmer_denovo_filter.cli import vcf_main
        with pytest.raises(SystemExit):
            vcf_main([])

    def test_discovery_main_missing_args_exits(self):
        """discovery_main exits when required args are missing."""
        from kmer_denovo_filter.cli import discovery_main
        with pytest.raises(SystemExit):
            discovery_main([])
