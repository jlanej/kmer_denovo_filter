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
             "([prefix].bed, [prefix].informative.bam, "
             "[prefix].sv.bedpe, [prefix].metrics.json)",
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
        "--informative-reads", default=None,
        help="Output BAM with reads carrying informative (child-unique) "
             "k-mers for IGV visualization",
    )
    parser.add_argument(
        "--min-mapq", type=int, default=20,
        help="Minimum mapping quality for child reads in VCF mode (default: 20). "
             "Not applied in discovery mode, which scans all primary reads "
             "regardless of mapping quality.",
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
        "--cluster-distance", type=int, default=500,
        help="Maximum gap (bp) for merging adjacent regions in discovery "
             "mode (default: 500)",
    )
    parser.add_argument(
        "--min-supporting-reads", type=int, default=1,
        help="Minimum number of supporting reads per region in discovery "
             "mode (default: 1)",
    )
    parser.add_argument(
        "--min-distinct-kmers", type=int, default=1,
        help="Minimum number of distinct proband-unique k-mers per region "
             "in discovery mode (default: 1)",
    )
    parser.add_argument(
        "--min-bedgraph-reads", type=int, default=3,
        help="Minimum number of distinct reads with at least one de novo "
             "k-mer at a position for it to be included in the bedGraph "
             "and read coverage BED (default: 3)",
    )
    parser.add_argument(
        "--min-distinct-kmers-per-read", type=int, default=None,
        help="Minimum number of distinct proband-unique k-mers a read "
             "must carry to be retained for region and bedGraph output. "
             "Applied before --min-supporting-reads and "
             "--min-bedgraph-reads filters. (default: k/4)",
    )
    parser.add_argument(
        "--parent-max-count", type=int, default=0,
        help="Maximum k-mer count in a parent before the k-mer is "
             "considered parental. K-mers with count > this value in "
             "either parent are removed (default: 0)",
    )
    parser.add_argument(
        "--sv-bedpe", default=None,
        help="Output BEDPE file for linked SV breakpoint pairs in "
             "discovery mode (default: [out-prefix].sv.bedpe)",
    )
    parser.add_argument(
        "--jf-hash-size", default=None,
        help="Initial hash size for Jellyfish count (e.g. '2G', '500M'). "
             "By default, estimated from the child BAM file size. "
             "A larger value reduces the chance of hash overflow "
             "(which creates multi-file indexes that require more memory "
             "to merge/dump), but uses more upfront RAM. "
             "Recommended: set to ~50%% of available memory.",
    )
    parser.add_argument(
        "--tmp-dir", default=None,
        help="Directory for temporary files (jellyfish indexes, intermediate "
             "FASTA files). Defaults to a subdirectory next to the output "
             "files. IMPORTANT: avoid RAM-backed filesystems like tmpfs "
             "(/tmp on many HPC systems), as intermediate files can exceed "
             "100 GB for WGS data.",
    )
    parser.add_argument(
        "--threads", "-t", type=int, default=4,
        help="Number of threads for jellyfish (default: 4)",
    )
    parser.add_argument(
        "--memory", type=float, default=None,
        help="Available memory in GB. On HPC systems (e.g. SLURM), set "
             "this to the allocated memory (e.g. --memory 64 for a 64 GB "
             "allocation) so that worker counts and hash sizes are tuned "
             "correctly. When omitted, auto-detected from the system.",
    )
    parser.add_argument(
        "--debug-kmers", action="store_true", default=False,
        help="Enable per-variant debug output",
    )
    parser.add_argument(
        "--kraken2-db", default=None,
        help="Path to a Kraken2 database for bacterial content flagging. "
             "In VCF mode, informative reads are classified with kraken2 and "
             "bacterial fraction annotations are added to the output VCF. "
             "Currently ignored in discovery mode. Requires kraken2 to be on "
             "PATH.",
    )
    parser.add_argument(
        "--kraken2-confidence", type=float, default=0.0,
        help="Kraken2 confidence threshold (0.0–1.0) for LCA "
             "classification (default: 0.0)",
    )
    parser.add_argument(
        "--kraken2-memory-mapping", action="store_true", default=False,
        help="Enable Kraken2 --memory-mapping to reduce RAM usage by mapping "
             "database files from disk (may run slower)",
    )
    parser.add_argument(
        "--kraken2-read-detail", default=None,
        help="Output path for the per-read Kraken2 classification detail "
             "BED file (bgzipped + tabix-indexed). When --kraken2-db is "
             "provided and this is omitted, auto-derived from --output "
             "(e.g. my_trio.annotated.kraken2_reads.bed.gz).",
    )
    parser.add_argument(
        "--kraken2-span-bed", default=None,
        help="Output path for the species-annotated genomic span BED file "
             "(bgzipped + tabix-indexed). Maps each classified read's aligned "
             "reference span to its Kraken2-assigned species, including "
             "soft-clip lengths and split-read indicators. When --kraken2-db "
             "is provided and this is omitted, auto-derived from --output "
             "(e.g. my_trio.annotated.kraken2_spans.bed.gz).",
    )
    parser.add_argument(
        "--no-expanded-bed", action="store_true", default=False,
        help="Disable generation of the soft-clip-expanded span BED file. "
             "By default, when --kraken2-db is provided, both the standard "
             "span BED and an expanded span BED (coordinates extended by "
             "soft-clip lengths) are written.",
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
