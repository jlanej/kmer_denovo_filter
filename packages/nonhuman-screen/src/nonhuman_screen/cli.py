"""Command-line interface: ``nonhuman-screen``.

Subcommands
-----------
``classify`` — screen reads from a BAM/CRAM for non-human content.  With
``--variants`` it computes the allele-based non-human fraction for every ALT
allele (the reads supporting each ALT are classified); otherwise it classifies
all mapped reads and reports a single batch summary.
"""

from __future__ import annotations

import argparse
import json
import logging
import shutil
import sys

from nonhuman_screen import __version__
from nonhuman_screen.alleles import _is_symbolic


def _build_parser():
    parser = argparse.ArgumentParser(
        prog="nonhuman-screen",
        description="Classify sequencing reads by non-human taxonomic content.",
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    sub = parser.add_subparsers(dest="command", required=True)

    c = sub.add_parser(
        "classify",
        help="Screen reads from a BAM/CRAM for non-human content.",
    )
    c.add_argument("--bam", required=True, help="Input BAM/CRAM.")
    c.add_argument(
        "--kraken2-db", required=True,
        help="Path to a kraken2 database directory.",
    )
    c.add_argument(
        "--ref-fasta", default=None,
        help="Reference FASTA (required for CRAM).",
    )
    c.add_argument(
        "--variants", default=None,
        help="VCF/BCF of variants. When given, compute the allele-based "
             "non-human fraction for every ALT allele instead of a whole-BAM "
             "summary.",
    )
    c.add_argument(
        "--min-baseq", type=int, default=0,
        help="Minimum base quality for a read to count as ALT support "
             "(default: 0).",
    )
    c.add_argument(
        "--out-prefix", default=None,
        help="Write outputs to <prefix>.variant_nhf.tsv / <prefix>.summary.json. "
             "When omitted, results are printed to stdout.",
    )
    c.add_argument("--threads", type=int, default=1, help="kraken2 threads.")
    c.add_argument(
        "--confidence", type=float, default=0.0,
        help="kraken2 confidence threshold 0.0-1.0 (default: 0.0).",
    )
    c.add_argument(
        "--memory-mapping", action="store_true", default=False,
        help="Pass --memory-mapping to kraken2 to reduce RAM.",
    )
    return parser


def _read_vcf_variants(path):
    """Yield ``(chrom, pos0, ref, alt)`` for every concrete ALT allele."""
    import pysam

    with pysam.VariantFile(path) as vcf:
        for rec in vcf:
            if rec.alts is None:
                continue
            for alt in rec.alts:
                # Skip symbolic (<DEL>), breakend (N[chr2:321[), and the
                # spanning-deletion (*) alleles — none has a literal sequence to
                # match reads against.
                if _is_symbolic(alt):
                    continue
                yield (rec.chrom, rec.start, rec.ref, alt)


_TSV_COLUMNS = (
    "variant_key", "supporting_reads", "nonhuman_fraction",
    "bacterial", "archaeal", "fungal", "protist", "viral", "univec_core",
    "human_lineage", "unclassified",
)


def _classify_variants(args):
    from nonhuman_screen.bam import classify_variants_alt_reads

    variants = list(_read_vcf_variants(args.variants))
    logging.info("Classifying ALT-supporting reads for %d alleles", len(variants))
    results = classify_variants_alt_reads(
        args.bam, args.kraken2_db, variants,
        ref_fasta=args.ref_fasta, min_baseq=args.min_baseq,
        confidence=args.confidence, threads=args.threads,
        memory_mapping=args.memory_mapping,
    )

    lines = ["\t".join(_TSV_COLUMNS)]
    for v in results:
        f = v.fractions
        lines.append("\t".join(str(x) for x in (
            v.variant_key, v.supporting_reads, v.nonhuman_fraction,
            f.bacterial, f.archaeal, f.fungal, f.protist, f.viral,
            f.univec_core, f.human_lineage, f.unclassified,
        )))
    tsv = "\n".join(lines) + "\n"

    if args.out_prefix:
        tsv_path = f"{args.out_prefix}.variant_nhf.tsv"
        with open(tsv_path, "w") as fh:
            fh.write(tsv)
        with open(f"{args.out_prefix}.summary.json", "w") as fh:
            json.dump([v.to_dict() for v in results], fh, indent=2)
        logging.info("Wrote %s", tsv_path)
    else:
        sys.stdout.write(tsv)
    return 0


def _classify_all(args):
    import pysam

    from nonhuman_screen.engine import Kraken2Runner

    sequences = {}
    with pysam.AlignmentFile(
        args.bam, reference_filename=args.ref_fasta or None
    ) as bam:
        for read in bam.fetch(until_eof=True):
            n = read.query_name
            if read.query_sequence and n not in sequences:
                sequences[n] = read.query_sequence

    runner = Kraken2Runner(
        args.kraken2_db, confidence=args.confidence, threads=args.threads,
        memory_mapping=args.memory_mapping,
    )
    result = runner.classify_sequences(sequences)
    summary = {
        "reads": result.total,
        "classified": result.classified,
        "taxonomy_available": result.taxonomy_available,
        "nonhuman_fraction": result.nonhuman_fraction,
        "fractions": result.fractions().to_dict(),
    }
    out = json.dumps(summary, indent=2)
    if args.out_prefix:
        with open(f"{args.out_prefix}.summary.json", "w") as fh:
            fh.write(out + "\n")
    else:
        sys.stdout.write(out + "\n")
    logging.info("%s", result.summary())
    return 0


def main(argv=None):
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s",
    )
    args = _build_parser().parse_args(argv)

    if args.command == "classify":
        if shutil.which("kraken2") is None:
            sys.stderr.write(
                "error: kraken2 not found on PATH. Install kraken2 "
                "(see docs/database.md) and try again.\n"
            )
            return 2
        if args.variants:
            return _classify_variants(args)
        return _classify_all(args)

    return 1  # pragma: no cover - argparse requires a subcommand


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
