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
        "--ref-fasta", "-r", required=True,
        help="Reference FASTA with .fai index",
    )
    parser.add_argument(
        "--vcf", required=True, help="Input VCF with candidate variants"
    )
    parser.add_argument(
        "--output", "-o", required=True, help="Output annotated VCF"
    )
    parser.add_argument(
        "--metrics", default=None, help="Output summary metrics JSON file"
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
    run_pipeline(args)
