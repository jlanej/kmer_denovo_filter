"""Command-line interface for kmer-denovo and kmer-discovery."""

import argparse
import sys


# ── Shared argument helpers ──────────────────────────────────────────


def _add_shared_args(parser):
    """Add arguments common to both VCF and discovery modes."""
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
        help="Reference FASTA with .fai index (required for CRAM input; "
             "also required for kmer-discovery unless --ref-jf is provided)",
    )
    parser.add_argument(
        "--kmer-size", "-k", type=int, default=31,
        help="K-mer size (default: 31)",
    )
    parser.add_argument(
        "--min-baseq", type=int, default=20,
        help="Minimum base quality for read k-mers in kmer-denovo "
             "(default: 20). Not used by kmer-discovery, which counts "
             "k-mers from all bases.",
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
        "--jf-hash-size", default=None,
        help="Hash size for the child's jellyfish count (e.g. '4G', '500M'). "
             "By default sized from --memory and from what the child's reads "
             "hold; the count runs with --disk, so a table that fills is "
             "dumped and merged, never doubled in memory. Give this only to "
             "override the automatic size.",
    )
    parser.add_argument(
        "--no-child-bloom-filter", action="store_true", default=False,
        help="Count the child's k-mers without the Bloom filter. By default "
             "(with --min-child-count 2 or more) a Bloom filter keeps k-mers "
             "seen once - most of a WGS child's distinct k-mers, sequencing "
             "errors - out of the hash, which shrinks it several-fold; the "
             "filter admits a k-mer from its second occurrence, and the dump "
             "bound is lowered by one to compensate.",
    )
    parser.add_argument(
        "--tmp-dir", default=None,
        help="Directory for temporary files (jellyfish indexes, intermediate "
             "FASTA files). Defaults to a subdirectory next to the output "
             "files. IMPORTANT: avoid RAM-backed filesystems like tmpfs "
             "(/tmp on many HPC systems), as intermediate files can exceed "
             "100 GB for WGS data.",
    )


def _add_nhf_screen_args(parser, db_required=False, compat=False):
    """The non-human fraction screen's options (kmer-discovery and kmer-nhf).

    The option names are VCF mode's (``--kraken2-db``,
    ``--kraken2-confidence``, ``--kraken2-memory-mapping``); the compat
    parser already defines those three, so *compat* adds only the one
    that is discovery's own.
    """
    if not compat:
        parser.add_argument(
            "--kraken2-db", default=None, required=db_required,
            help="Path to a Kraken2 database for the non-human fraction "
                 "screen: the informative reads of every SV and AMBIGUOUS "
                 "region (and of every junction end), and a sample of all "
                 "informative reads, are classified with kraken2, and each "
                 "screened region gets its non-human fraction in the BED "
                 "columns nhf_reads and nhf ('.' = not screened: a blank is "
                 "not 0). Annotated, never filtered. Without this option the "
                 "screen is recorded as disabled; kraken2 or the database "
                 "failing is recorded as failed and the run still succeeds "
                 "(kmer-nhf adds the screen to a finished run). Requires "
                 "kraken2 on PATH.",
        )
        parser.add_argument(
            "--kraken2-confidence", type=float, default=0.0,
            help="Kraken2 confidence threshold (0.0–1.0) for LCA "
                 "classification (default: 0.0)",
        )
        parser.add_argument(
            "--kraken2-memory-mapping", action="store_true", default=False,
            help="Enable Kraken2 --memory-mapping to reduce RAM usage by "
                 "mapping database files from disk (may run slower)",
        )
    parser.add_argument(
        "--nhf-sample-reads", type=int, default=200000, metavar="N",
        help="With --kraken2-db: at most this many informative reads (read "
             "ends, chosen by a hash of the read name) for the run-level "
             "non-human fraction (default: 200000)",
    )


# ── VCF-mode parser ─────────────────────────────────────────────────


def parse_vcf_args(argv=None):
    """Parse arguments for the VCF annotation pipeline (kmer-denovo)."""
    parser = argparse.ArgumentParser(
        prog="kmer-denovo",
        description="De novo variant curation using k-mer analysis (VCF mode)",
    )
    _add_shared_args(parser)
    parser.add_argument(
        "--vcf", required=True,
        help="Input VCF with candidate variants",
    )
    parser.add_argument(
        "--output", "-o", required=True, help="Output annotated VCF"
    )
    parser.add_argument(
        "--metrics", default=None, help="Output summary metrics JSON file"
    )
    parser.add_argument(
        "--summary", default=None,
        help="Output human-readable summary of variant stats and likely DNMs",
    )
    parser.add_argument(
        "--informative-reads", default=None,
        help="Output BAM with reads carrying informative (child-unique) "
             "k-mers for IGV visualization",
    )
    parser.add_argument(
        "--min-mapq", type=int, default=20,
        help="Minimum mapping quality for child reads (default: 20)",
    )
    parser.add_argument(
        "--proband-id", default=None,
        help="Sample ID of the proband in the VCF. When provided and "
             "matching a VCF sample, DKU/DKT/DKA are written as FORMAT "
             "fields on that sample; otherwise they are written as INFO "
             "fields.",
    )
    parser.add_argument(
        "--kraken2-db", default=None,
        help="Path to a Kraken2 database for non-human content classification. "
             "Informative reads are classified with kraken2 and "
             "non-human fraction annotations are added to the output VCF. "
             "Requires kraken2 to be on PATH.",
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
    parser.add_argument(
        "--report", default=None,
        help="Output path for an interactive HTML report summarizing "
             "the k-mer filtering results with Plotly visualizations. "
             "When omitted, no report is generated.",
    )
    return parser.parse_args(argv)


def parse_discovery_args(argv=None):
    """Parse arguments for the VCF-free discovery pipeline (kmer-discovery)."""
    parser = argparse.ArgumentParser(
        prog="kmer-discovery",
        description="VCF-free de novo k-mer discovery pipeline",
    )
    _add_shared_args(parser)
    # Discovery mode has no base-quality filter; None records that the
    # flag was not given, so discovery_main can warn when it is.
    parser.set_defaults(min_baseq=None)
    parser.add_argument(
        "--out-prefix", required=True,
        help="Output prefix for discovery mode files "
             "([prefix].bed, [prefix].informative.bam, "
             "[prefix].sv.bedpe, [prefix].kmer_coverage.bedgraph, "
             "[prefix].read_coverage.bed, [prefix].metrics.json, "
             "[prefix].summary.txt)",
    )
    parser.add_argument(
        "--ref-jf", default=None,
        help="Path to a precomputed Jellyfish reference index. "
             "Defaults to [ref-fasta].k[kmer-size].jf",
    )
    parser.add_argument(
        "--min-child-count", type=int, default=3,
        help="Minimum child k-mer occurrences (default: 3)",
    )
    parser.add_argument(
        "--candidate-summary", default=None,
        help="Path to a VCF-mode summary.txt for candidate comparison. "
             "High-quality de novos (DKA_DKT > 0.25, "
             "DKA > 10) are checked against discovered regions.",
    )
    parser.add_argument(
        "--dnm-regions", default=None, metavar="TSV",
        help="Known de novo events to check the discovered regions "
             "against: a tab-separated file with columns chrom, 1-based "
             "pos, size in bp ('.' if unknown) and event type. Adds a "
             "per-event table to the metrics JSON, summary and report: an "
             "event is supported when split, discordant or "
             "breakpoint-clipped reads lie within --cluster-distance of "
             "one of its ends; regions anywhere in its span are reported "
             "separately. See "
             "examples/HG002_trio/sulovari2023_dnm_regions.tsv.",
    )
    parser.add_argument(
        "--dnm-negative-regions", default=None, metavar="TSV",
        help="Events in the same format that should NOT be found - "
             "inherited variants, whose k-mers are in a parent. Evaluated "
             "like --dnm-regions and reported beside them as the "
             "false-positive control.",
    )
    parser.add_argument(
        "--dnm-null", type=int, default=20, metavar="N",
        help="With --dnm-regions: also evaluate N random loci per event, "
             "of the same size on the same chromosome, and report how often "
             "chance alone supports them (default: 20; 0 disables)",
    )
    parser.add_argument(
        "--discard-proband-kmers", action="store_true",
        help="Do not keep the proband-unique k-mer index "
             "([out-prefix].proband_unique.jf); by default it is copied "
             "out of the temporary directory for later use",
    )
    parser.add_argument(
        "--cluster-distance", type=int, default=500,
        help="Maximum gap (bp) for merging adjacent regions; also how far "
             "a split alignment or mate may land outside a region and "
             "still link to it in the BEDPE, and how far a region may lie "
             "from a --dnm-regions event and still count toward it "
             "(default: 500)",
    )
    parser.add_argument(
        "--min-supporting-reads", type=int, default=1,
        help="Minimum number of supporting reads per region (default: 1)",
    )
    parser.add_argument(
        "--min-distinct-kmers", type=int, default=1,
        help="Minimum number of distinct proband-unique k-mers per region "
             "(default: 1)",
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
             "must carry to be retained for region, bedGraph and "
             "informative BAM output. "
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
        help="Output BEDPE file for linked SV breakpoint pairs "
             "(default: [out-prefix].sv.bedpe)",
    )
    parser.add_argument(
        "--report", default=None,
        help="Output path for an interactive HTML report summarizing "
             "the discovery results with Plotly visualizations. "
             "When omitted, no report is generated.",
    )
    _add_nhf_screen_args(parser)
    return parser.parse_args(argv)


def parse_nhf_args(argv=None):
    """Parse arguments for the after-the-fact screen (``kmer-nhf``)."""
    parser = argparse.ArgumentParser(
        prog="kmer-nhf",
        description="The non-human fraction screen of an existing "
                    "kmer-discovery run: classify its informative reads with "
                    "kraken2 exactly as kmer-discovery --kraken2-db does, and "
                    "rewrite the BED columns nhf_reads and nhf, the metrics "
                    "block nonhuman_fraction and the summary section in "
                    "place. For runs made before the option existed, or "
                    "with a new database.",
    )
    parser.add_argument(
        "--out-prefix", required=True,
        help="The run's output prefix: [prefix].informative.bam, [prefix].bed, "
             "[prefix].sv.bedpe and [prefix].metrics.json are read; the BED, "
             "the metrics and [prefix].summary.txt are rewritten in place",
    )
    _add_nhf_screen_args(parser, db_required=True)
    parser.add_argument(
        "--cluster-distance", type=int, default=500,
        help="The run's --cluster-distance: how far outside a region an "
             "unmapped read may be placed and still belong to it "
             "(default: 500)",
    )
    parser.add_argument(
        "--threads", "-t", type=int, default=4,
        help="Number of threads for kraken2 (default: 4)",
    )
    parser.add_argument(
        "--tmp-dir", default=None,
        help="Directory for kraken2's temporary FASTQ (default: a "
             "subdirectory next to the output files)",
    )
    parser.add_argument(
        "--sv-bedpe", default=None,
        help="The run's BEDPE when it is not [prefix].sv.bedpe",
    )
    parser.add_argument(
        "--report", default=None,
        help="Also regenerate the interactive HTML report from the "
             "rewritten metrics and summary, at this path",
    )
    return parser.parse_args(argv)


def parse_args(argv=None):
    """Parse arguments accepting both VCF and discovery flags.

    Retained for backward compatibility and existing tests.
    New code should use :func:`parse_vcf_args` or
    :func:`parse_discovery_args` directly.
    """
    parser = argparse.ArgumentParser(
        prog="kmer-denovo",
        description="De novo variant curation using k-mer analysis",
    )
    _add_shared_args(parser)
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
             "[prefix].sv.bedpe, [prefix].kmer_coverage.bedgraph, "
             "[prefix].read_coverage.bed, [prefix].metrics.json, "
             "[prefix].summary.txt)",
    )
    parser.add_argument(
        "--metrics", default=None, help="Output summary metrics JSON file"
    )
    parser.add_argument(
        "--summary", default=None,
        help="Output human-readable summary of variant stats and likely DNMs",
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
        "--dnm-regions", default=None, metavar="TSV",
        help="Known de novo events to check the discovered regions "
             "against in discovery mode: a tab-separated file with columns "
             "chrom, 1-based pos, size in bp ('.' if unknown) and event "
             "type. Adds a per-event table to the metrics JSON, summary and "
             "report: an event is supported when split, discordant or "
             "breakpoint-clipped reads lie within --cluster-distance of one "
             "of its ends; regions anywhere in its span are reported "
             "separately.",
    )
    parser.add_argument(
        "--dnm-negative-regions", default=None, metavar="TSV",
        help="Events in the same format that should NOT be found - "
             "inherited variants, whose k-mers are in a parent. Evaluated "
             "like --dnm-regions and reported beside them as the "
             "false-positive control.",
    )
    parser.add_argument(
        "--dnm-null", type=int, default=20, metavar="N",
        help="With --dnm-regions: also evaluate N random loci per event, "
             "of the same size on the same chromosome, and report how often "
             "chance alone supports them (default: 20; 0 disables)",
    )
    parser.add_argument(
        "--discard-proband-kmers", action="store_true",
        help="Do not keep the proband-unique k-mer index "
             "([out-prefix].proband_unique.jf); by default it is copied "
             "out of the temporary directory for later use",
    )
    parser.add_argument(
        "--cluster-distance", type=int, default=500,
        help="Maximum gap (bp) for merging adjacent regions in discovery "
             "mode; also how far a split alignment or mate may land "
             "outside a region and still link to it, and how far a region "
             "may lie from a --dnm-regions event and still count toward it "
             "(default: 500)",
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
             "must carry to be retained for region, bedGraph and "
             "informative BAM output. "
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
        "--kraken2-db", default=None,
        help="Path to a Kraken2 database for non-human content classification. "
             "In VCF mode, informative reads are classified with kraken2 and "
             "non-human fraction annotations are added to the output VCF. "
             "In discovery mode, the informative reads of every SV and "
             "AMBIGUOUS region and a sample of all informative reads are "
             "classified, and the BED gets the columns nhf_reads and nhf. "
             "Requires kraken2 to be on PATH.",
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
    parser.add_argument(
        "--report", default=None,
        help="Output path for an interactive HTML report summarizing "
             "the k-mer filtering results with Plotly visualizations. "
             "When omitted, no report is generated.",
    )
    _add_nhf_screen_args(parser, compat=True)
    return parser.parse_args(argv)


def vcf_main(argv=None):
    """Entry point for the VCF annotation pipeline (``kmer-denovo``)."""
    from kmer_denovo_filter.vcf.pipeline import run_pipeline

    args = parse_vcf_args(argv)
    run_pipeline(args)


def discovery_main(argv=None):
    """Entry point for the VCF-free discovery pipeline (``kmer-discovery``)."""
    from kmer_denovo_filter.discovery.pipeline import run_discovery_pipeline

    args = parse_discovery_args(argv)
    if args.min_baseq is not None:
        print("warning: --min-baseq has no effect in kmer-discovery; "
              "k-mers are counted from all bases",
              file=sys.stderr)
    run_discovery_pipeline(args)


def nhf_main(argv=None):
    """Entry point for the after-the-fact non-human fraction screen (``kmer-nhf``).

    Exit status 0 when the run was screened, 1 when the screen failed
    (the failure is recorded in the run's files all the same) and 2 when
    the prefix is not a kmer-discovery run.
    """
    import logging
    import os

    from kmer_denovo_filter.discovery import nhf
    from kmer_denovo_filter.utils import _resolve_tmp_dir

    args = parse_nhf_args(argv)
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s %(levelname)s %(message)s")
    out_dir = os.path.dirname(os.path.abspath(args.out_prefix)) or "."
    tmp_root = _resolve_tmp_dir(args.tmp_dir, out_dir)
    try:
        block = nhf.rescreen_run(
            args.out_prefix, args.kraken2_db, bedpe_path=args.sv_bedpe,
            cluster_distance=args.cluster_distance,
            sample_reads=args.nhf_sample_reads,
            confidence=args.kraken2_confidence,
            memory_mapping=args.kraken2_memory_mapping,
            threads=args.threads, tmpdir=tmp_root,
        )
    except FileNotFoundError as exc:
        print(f"error: {exc} not found: is {args.out_prefix} the output "
              f"prefix of a kmer-discovery run?", file=sys.stderr)
        sys.exit(2)
    finally:
        try:
            if not args.tmp_dir:
                os.rmdir(tmp_root)  # only succeeds if empty
        except OSError:
            pass
    if args.report:
        from kmer_denovo_filter.report import generate_report
        generate_report(
            output_path=args.report,
            discovery_metrics_path=f"{args.out_prefix}.metrics.json",
            discovery_summary_path=f"{args.out_prefix}.summary.txt",
        )
    print(nhf.describe(block))
    if block["status"] == nhf.STATUS_FAILED:
        sys.exit(1)


def parse_report_args(argv=None):
    """Parse arguments for the standalone report generator (``kmer-report``)."""
    parser = argparse.ArgumentParser(
        prog="kmer-report",
        description=(
            "Generate an interactive HTML report from kmer-denovo / "
            "kmer-discovery output files.  Useful for regenerating a "
            "report from existing pipeline outputs without re-running "
            "the full pipeline."
        ),
    )
    parser.add_argument(
        "--output", "-o", required=True,
        help="Output path for the HTML report.",
    )
    parser.add_argument(
        "--vcf-metrics", default=None,
        help="Path to VCF-mode metrics.json produced by kmer-denovo.",
    )
    parser.add_argument(
        "--vcf-summary", default=None,
        help="Path to VCF-mode summary.txt produced by kmer-denovo.",
    )
    parser.add_argument(
        "--vcf", default=None,
        help="Path to annotated VCF produced by kmer-denovo "
             "(used for Kraken2 annotations if present).",
    )
    parser.add_argument(
        "--discovery-metrics", default=None,
        help="Path to discovery-mode metrics.json produced by kmer-discovery.",
    )
    parser.add_argument(
        "--discovery-summary", default=None,
        help="Path to discovery-mode summary.txt produced by kmer-discovery.",
    )
    return parser.parse_args(argv)


def report_main(argv=None):
    """Entry point for the standalone report generator (``kmer-report``)."""
    from kmer_denovo_filter.report import generate_report

    args = parse_report_args(argv)
    result = generate_report(
        output_path=args.output,
        vcf_metrics_path=args.vcf_metrics,
        vcf_summary_path=args.vcf_summary,
        vcf_path=args.vcf,
        discovery_metrics_path=args.discovery_metrics,
        discovery_summary_path=args.discovery_summary,
    )
    print(f"Report written to: {result}")


def main(argv=None):
    """Backward-compatible entry point that dispatches by mode."""
    from kmer_denovo_filter.vcf.pipeline import run_pipeline

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
        from kmer_denovo_filter.discovery.pipeline import run_discovery_pipeline
        run_discovery_pipeline(args)
