"""Command-line interface for kmer-denovo."""

import argparse
import sys

from kmer_denovo_filter.pipeline import run_pipeline


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        prog="kmer-denovo",
        description="De novo variant curation using k-mer analysis",
    )
    parser.add_argument(
        "--child", required=True, help="Child BAM/CRAM file (indexed)"
    )
    parser.add_argument(
        "--mother", required=True, help="Mother BAM/CRAM file (indexed)"
    )
    parser.add_argument(
        "--father", required=True, help="Father BAM/CRAM file (indexed)"
    )
    parser.add_argument(
        "--ref-fasta", "-r", default=None,
        help="Reference FASTA with .fai index (required for CRAM input)",
    )
    parser.add_argument(
        "--vcf", default=None,
        help="Input VCF with candidate variants. When omitted, the tool "
             "runs in VCF-free discovery mode (requires --out-prefix)",
    )
    parser.add_argument(
        "--output", "-o", default=None, help="Output annotated VCF"
    )
    parser.add_argument(
        "--ref-jf", default=None,
        help="Path to a precomputed Jellyfish reference index. "
             "Defaults to [ref-fasta].k[kmer-size].jf",
    )
    parser.add_argument(
        "--min-child-count", type=int, default=3,
        help="Minimum child k-mer occurrences in discovery mode (default: 3)",
    )
    parser.add_argument(
        "--out-prefix", default=None,
        help="Output prefix for discovery mode files "
             "([prefix].bed, [prefix].informative.bam, [prefix].metrics.json)",
    )
    parser.add_argument(
        "--metrics", default=None, help="Output summary metrics JSON file"
    )
    parser.add_argument(
        "--summary", default=None,
        help="Output human-readable summary of variant stats and likely DNMs",
    )
    parser.add_argument(
        "--kmer-size", "-k", type=int, default=31,
        help="K-mer size (default: 31)",
    )
    parser.add_argument(
        "--min-baseq", type=int, default=20,
        help="Minimum base quality for read k-mers (default: 20)",
    )
    parser.add_argument(
        "--min-mapq", type=int, default=20,
        help="Minimum mapping quality for child reads (default: 20)",
    )
    parser.add_argument(
        "--informative-reads", default=None,
        help="Output BAM with reads carrying informative (child-unique) "
             "k-mers for IGV visualization",
    )
    parser.add_argument(
        "--proband-id", default=None,
        help="Sample ID of the proband in the VCF. When provided and "
             "matching a VCF sample, DKU/DKT/DKA are written as FORMAT "
             "fields on that sample; otherwise they are written as INFO "
             "fields.",
    )
    parser.add_argument(
        "--candidate-summary", default=None,
        help="Path to a VCF-mode summary.txt for candidate comparison "
             "in discovery mode. High-quality de novos (DKA_DKT > 0.25, "
             "DKA > 10) are checked against discovered regions.",
    )
    parser.add_argument(
        "--threads", "-t", type=int, default=4,
        help="Number of threads for jellyfish (default: 4)",
    )
    parser.add_argument(
        "--debug-kmers", action="store_true", default=False,
        help="Enable per-variant debug output",
    )
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)

    # Validate mode: either VCF mode or discovery mode
    if args.vcf is not None:
        # VCF mode requires --output
        if args.output is None:
            print("error: --output is required when --vcf is provided",
                  file=sys.stderr)
            sys.exit(2)
        run_pipeline(args)
    else:
        # Discovery mode requires --out-prefix
        if args.out_prefix is None:
            print("error: either --vcf (with --output) or --out-prefix "
                  "(for discovery mode) must be provided",
                  file=sys.stderr)
            sys.exit(2)
        from kmer_denovo_filter.pipeline import run_discovery_pipeline
        run_discovery_pipeline(args)
