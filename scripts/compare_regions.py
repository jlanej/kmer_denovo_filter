#!/usr/bin/env python3
"""Compare discovery regions and VCF variants using bedGraph k-mer coverage.

For each VCF variant this script checks whether its position overlaps with
a discovery BED region and whether the bedGraph shows k-mer signal at that
location.  Regions are then classified as:

  CONCORDANT       – VCF variant has k-mer signal in the bedGraph AND
                     overlaps a discovery BED region.
  VCF_ONLY         – VCF variant has k-mer signal but no overlapping
                     discovery region.
  DISCOVERY_ONLY   – Discovery region has no overlapping VCF variant
                     but does have bedGraph k-mer signal.
  DISCORDANT       – VCF variant (with k-mer signal) has no discovery
                     region overlap.

A final summary is written to stdout (and optionally to a file).

Usage
-----
    python scripts/compare_regions.py \\
        --bedgraph  tests/example_output_discovery/giab_discovery.kmer_coverage.bedgraph \\
        --discovery tests/example_output_discovery/giab_discovery.bed \\
        --vcf       tests/example_output/annotated.vcf.gz \\
        [--output   comparison_summary.txt]
"""

import argparse
import collections
import sys

import pysam


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------

def load_bedgraph(path):
    """Return a dict mapping chrom -> list of (start, end, count) tuples.

    Positions are 0-based, half-open [start, end).
    """
    intervals = collections.defaultdict(list)
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#") or line.startswith("track"):
                continue
            parts = line.split("\t")
            chrom, start, end, count = parts[0], int(parts[1]), int(parts[2]), int(parts[3])
            intervals[chrom].append((start, end, count))
    return dict(intervals)


def load_discovery_bed(path):
    """Return a dict mapping chrom -> list of region dicts.

    Each region dict has keys: start, end, reads, unique_kmers, split_reads,
    discordant_pairs, max_clip_len, unmapped_mates, class.
    Positions are 0-based, half-open [start, end).
    """

    def _int(val, default=0):
        try:
            return int(val)
        except (ValueError, TypeError):
            return default

    regions = collections.defaultdict(list)
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            regions[parts[0]].append({
                "start": int(parts[1]),
                "end": int(parts[2]),
                "reads": _int(parts[3] if len(parts) > 3 else None),
                "unique_kmers": _int(parts[4] if len(parts) > 4 else None),
                "split_reads": _int(parts[5] if len(parts) > 5 else None),
                "discordant_pairs": _int(parts[6] if len(parts) > 6 else None),
                "max_clip_len": _int(parts[7] if len(parts) > 7 else None),
                "unmapped_mates": _int(parts[8] if len(parts) > 8 else None),
                "class": parts[9].strip() if len(parts) > 9 else "UNKNOWN",
            })
    return dict(regions)


def load_vcf_variants(path):
    """Return a list of dicts with per-variant metadata.

    VCF positions are converted to 0-based coordinates for comparison.
    """
    variants = []
    with pysam.VariantFile(path) as vcf:
        sample_names = list(vcf.header.samples)
        for rec in vcf:
            # Convert 1-based VCF POS to 0-based
            pos0 = rec.pos - 1
            dku, dka = None, None
            if sample_names:
                sample = rec.samples[sample_names[0]]
                dku = sample.get("DKU")
                dka = sample.get("DKA")
            else:
                dku = rec.info.get("DKU")
                dka = rec.info.get("DKA")
            variants.append({
                "chrom": rec.chrom,
                "pos0": pos0,  # 0-based
                "pos1": rec.pos,  # 1-based (for display)
                "ref": rec.ref,
                "alt": ",".join(str(a) for a in rec.alts),
                "dku": dku,
                "dka": dka,
            })
    return variants


# ---------------------------------------------------------------------------
# Overlap helpers
# ---------------------------------------------------------------------------

def _has_bedgraph_signal(chrom, pos0, bedgraph, window=0):
    """Return True if any bedGraph interval covers [pos0-window, pos0+window+1)."""
    intervals = bedgraph.get(chrom, [])
    q_start = pos0 - window
    q_end = pos0 + window + 1
    for start, end, count in intervals:
        if count > 0 and start < q_end and end > q_start:
            return True
    return False


def _overlapping_discovery_regions(chrom, pos0, discovery, window=0):
    """Return list of discovery region dicts that overlap pos0 ± window."""
    regions = discovery.get(chrom, [])
    q_start = pos0 - window
    q_end = pos0 + window + 1
    return [r for r in regions if r["start"] < q_end and r["end"] > q_start]


def _vcf_variants_in_region(chrom, reg_start, reg_end, variants_by_chrom):
    """Return VCF variants whose position falls in [reg_start, reg_end)."""
    return [
        v for v in variants_by_chrom.get(chrom, [])
        if reg_start <= v["pos0"] < reg_end
    ]


# ---------------------------------------------------------------------------
# Comparison logic
# ---------------------------------------------------------------------------

def compare(bedgraph, discovery, variants, window=0):
    """Classify VCF variants and discovery regions.

    Args:
        bedgraph:  dict chrom -> [(start, end, count), ...]
        discovery: dict chrom -> [(start, end), ...]
        variants:  list of variant dicts from load_vcf_variants()
        window:    extra bases around each VCF position to search (default 0)

    Returns:
        A dict with keys:
            "concordant"      – variants with signal AND discovery overlap
            "vcf_only"        – variants with signal but NO discovery overlap
            "discordant_vcf"  – variants with signal, no discovery overlap
            "discovery_only"  – discovery regions with no VCF variant inside
            "no_signal"       – variants with no bedGraph signal at all
    """
    # Index variants by chromosome for fast region lookups
    variants_by_chrom = collections.defaultdict(list)
    for v in variants:
        variants_by_chrom[v["chrom"]].append(v)

    concordant = []
    vcf_only = []
    no_signal = []

    for v in variants:
        has_signal = _has_bedgraph_signal(v["chrom"], v["pos0"], bedgraph, window)
        overlaps = _overlapping_discovery_regions(v["chrom"], v["pos0"], discovery, window)

        if has_signal and overlaps:
            concordant.append({"variant": v, "regions": overlaps})
        elif has_signal and not overlaps:
            vcf_only.append({"variant": v})
        else:
            no_signal.append({
                "variant": v,
                "has_discovery": bool(overlaps),
                "discovery_regions": overlaps,
            })

    # Find discovery regions with no VCF variant inside
    discovery_only = []
    for chrom, regions in sorted(discovery.items()):
        for region in regions:
            inside = _vcf_variants_in_region(
                chrom, region["start"], region["end"], variants_by_chrom,
            )
            if not inside:
                discovery_only.append({"chrom": chrom, **region})

    return {
        "concordant": concordant,
        "vcf_only": vcf_only,
        "no_signal": no_signal,
        "discovery_only": discovery_only,
    }


# ---------------------------------------------------------------------------
# Summary formatting
# ---------------------------------------------------------------------------

def _fmt_variant(v):
    return f"{v['chrom']}:{v['pos1']} {v['ref']}>{v['alt']}"


def _fmt_region_stats(region):
    """Return a compact stats string for a discovery region."""
    return (
        f"reads={region['reads']}"
        f"  unique_kmers={region['unique_kmers']}"
        f"  split_reads={region['split_reads']}"
        f"  class={region['class']}"
    )


def format_summary(result, window=0):
    """Return a human-readable summary string."""
    lines = []
    lines.append("=" * 60)
    lines.append("  bedGraph / Discovery / VCF Region Comparison")
    lines.append("=" * 60)

    if window:
        lines.append(f"  Search window: ±{window} bp around each VCF position")
    else:
        lines.append("  Search window: exact position overlap")
    lines.append("")

    # ── Concordant ────────────────────────────────────────────────
    concordant = result["concordant"]
    lines.append(f"CONCORDANT  ({len(concordant)} variants)")
    lines.append("  VCF variant has k-mer signal AND overlaps a discovery region")
    lines.append("-" * 60)
    for item in concordant:
        v = item["variant"]
        for region in item["regions"]:
            region_coord = f"{v['chrom']}:{region['start']}-{region['end']}"
            lines.append(
                f"  {_fmt_variant(v)}"
                f"  DKU={v['dku']}  DKA={v['dka']}"
                f"  region={region_coord}"
                f"  {_fmt_region_stats(region)}"
            )
    if not concordant:
        lines.append("  (none)")
    lines.append("")

    # ── VCF-only (signal but no discovery region) ─────────────────
    vcf_only = result["vcf_only"]
    lines.append(f"VCF_ONLY  ({len(vcf_only)} variants)")
    lines.append("  VCF variant has k-mer signal but no overlapping discovery region")
    lines.append("-" * 60)
    for item in vcf_only:
        v = item["variant"]
        lines.append(
            f"  {_fmt_variant(v)}"
            f"  DKU={v['dku']}  DKA={v['dka']}"
        )
    if not vcf_only:
        lines.append("  (none)")
    lines.append("")

    # ── No signal ─────────────────────────────────────────────────
    no_signal = result["no_signal"]
    lines.append(f"NO_SIGNAL  ({len(no_signal)} variants)")
    lines.append("  VCF variant has no bedGraph k-mer signal at its position")
    lines.append("-" * 60)
    for item in no_signal:
        v = item["variant"]
        if item["has_discovery"]:
            for region in item["discovery_regions"]:
                region_coord = f"{v['chrom']}:{region['start']}-{region['end']}"
                lines.append(
                    f"  {_fmt_variant(v)}"
                    f"  DKU={v['dku']}  DKA={v['dka']}"
                    f"  +discovery={region_coord}"
                    f"  {_fmt_region_stats(region)}"
                )
        else:
            lines.append(
                f"  {_fmt_variant(v)}"
                f"  DKU={v['dku']}  DKA={v['dka']}"
            )
    if not no_signal:
        lines.append("  (none)")
    lines.append("")

    # ── Discovery-only (no VCF variant inside) ───────────────────
    disc_only = result["discovery_only"]
    lines.append(f"DISCOVERY_ONLY  ({len(disc_only)} regions)")
    lines.append("  Discovery region has no overlapping VCF variant")
    lines.append("-" * 60)
    for item in disc_only:
        lines.append(
            f"  {item['chrom']}:{item['start']}-{item['end']}"
            f"  ({item['end'] - item['start']} bp)"
            f"  {_fmt_region_stats(item)}"
        )
    if not disc_only:
        lines.append("  (none)")
    lines.append("")

    # ── Overall counts ────────────────────────────────────────────
    total_vcf = len(concordant) + len(vcf_only) + len(no_signal)
    lines.append("=" * 60)
    lines.append("  Summary")
    lines.append("=" * 60)
    lines.append(f"  Total VCF variants:            {total_vcf}")
    lines.append(f"  Concordant (signal + region):  {len(concordant)}")
    lines.append(f"  VCF-only (signal, no region):  {len(vcf_only)}")
    lines.append(f"  No k-mer signal:               {len(no_signal)}")
    lines.append(f"  Discovery-only regions:        {len(disc_only)}")
    lines.append("=" * 60)

    return "\n".join(lines)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        prog="compare-regions",
        description=(
            "Compare discovery BED regions and VCF variants using "
            "bedGraph k-mer coverage.  Produces a concordance summary."
        ),
    )
    parser.add_argument(
        "--bedgraph", "-b", required=True,
        help="bedGraph file from the discovery pipeline "
             "([prefix].kmer_coverage.bedgraph)",
    )
    parser.add_argument(
        "--discovery", "-d", required=True,
        help="Discovery BED file ([prefix].bed)",
    )
    parser.add_argument(
        "--vcf", "-v", required=True,
        help="Annotated VCF (or VCF.gz) from the VCF-mode pipeline",
    )
    parser.add_argument(
        "--output", "-o", default=None,
        help="Write summary to this file in addition to stdout",
    )
    parser.add_argument(
        "--window", "-w", type=int, default=0,
        help="Extra bases around each VCF position to search for signal "
             "and discovery regions (default: 0, exact position overlap)",
    )
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)

    bedgraph = load_bedgraph(args.bedgraph)
    discovery = load_discovery_bed(args.discovery)
    variants = load_vcf_variants(args.vcf)

    result = compare(bedgraph, discovery, variants, window=args.window)
    summary = format_summary(result, window=args.window)

    print(summary)
    if args.output:
        with open(args.output, "w") as fh:
            fh.write(summary + "\n")


if __name__ == "__main__":
    main()
