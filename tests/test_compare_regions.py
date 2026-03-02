"""Tests for scripts/compare_regions.py.

Tests cover:
- Loading each input format (bedGraph, discovery BED, VCF)
- The compare() classification logic
- The format_summary() output
- An integration smoke-test using the GIAB example files (when available)
"""

import os
import sys
import tempfile

import pysam
import pytest

# Make the scripts directory importable
_SCRIPTS_DIR = os.path.join(os.path.dirname(__file__), "..", "scripts")
sys.path.insert(0, os.path.abspath(_SCRIPTS_DIR))

from compare_regions import (  # noqa: E402
    compare,
    format_summary,
    load_bedgraph,
    load_discovery_bed,
    load_vcf_variants,
    main,
    parse_args,
)

GIAB_DIR = os.path.join(os.path.dirname(__file__), "data", "giab")
EXAMPLE_OUTPUT_DIR = os.path.join(os.path.dirname(__file__), "example_output")
EXAMPLE_DISCOVERY_DIR = os.path.join(
    os.path.dirname(__file__), "example_output_discovery",
)

GIAB_DATA_EXISTS = (
    os.path.isfile(os.path.join(GIAB_DIR, "HG002_child.bam"))
    and os.path.isfile(
        os.path.join(EXAMPLE_DISCOVERY_DIR, "giab_discovery.kmer_coverage.bedgraph")
    )
)


# ---------------------------------------------------------------------------
# Helpers for synthetic test data
# ---------------------------------------------------------------------------

def _write_bedgraph_file(path, intervals):
    """Write a bedGraph file from [(chrom, start, end, count)] tuples."""
    with open(path, "w") as fh:
        for chrom, start, end, count in intervals:
            fh.write(f"{chrom}\t{start}\t{end}\t{count}\n")


def _write_bed_file(path, regions):
    """Write a minimal BED file from [(chrom, start, end)] tuples."""
    with open(path, "w") as fh:
        for chrom, start, end in regions:
            fh.write(f"{chrom}\t{start}\t{end}\t.\t.\t.\t.\t.\t.\t.\n")


def _write_vcf_gz(path, variants):
    """Create a bgzipped VCF with sample 'PROBAND' from a list of dicts.

    Each dict should have: chrom, pos (1-based), ref, alt.
    Optional: dku (int), dka (int).
    """
    header = pysam.VariantHeader()
    header.add_sample("PROBAND")
    header.add_meta(
        "FORMAT",
        items=[
            ("ID", "DKU"), ("Number", "1"), ("Type", "Integer"),
            ("Description", "De novo k-mer unique count"),
        ],
    )
    header.add_meta(
        "FORMAT",
        items=[
            ("ID", "DKA"), ("Number", "1"), ("Type", "Integer"),
            ("Description", "De novo k-mer alt count"),
        ],
    )
    # Add reference contigs using add_line (required by pysam for chrom assignment)
    chroms = sorted({v["chrom"] for v in variants})
    for chrom in chroms:
        header.add_line(f"##contig=<ID={chrom},length=300000000>")

    with pysam.VariantFile(path, "wz", header=header) as vcf:
        for v in variants:
            rec = vcf.new_record()
            rec.chrom = v["chrom"]
            rec.pos = v["pos"]  # 1-based
            rec.ref = v["ref"]
            rec.alts = (v["alt"],)
            rec.samples["PROBAND"]["DKU"] = v.get("dku", 0)
            rec.samples["PROBAND"]["DKA"] = v.get("dka", 0)
            vcf.write(rec)
    pysam.tabix_index(path, preset="vcf", force=True)


# ---------------------------------------------------------------------------
# Unit tests: load_bedgraph
# ---------------------------------------------------------------------------

class TestLoadBedgraph:
    def test_basic_load(self, tmp_path):
        bg = tmp_path / "test.bedgraph"
        _write_bedgraph_file(str(bg), [
            ("chr1", 100, 200, 5),
            ("chr1", 200, 300, 3),
            ("chr2", 50, 60, 10),
        ])
        result = load_bedgraph(str(bg))
        assert "chr1" in result
        assert "chr2" in result
        assert (100, 200, 5) in result["chr1"]
        assert (200, 300, 3) in result["chr1"]
        assert (50, 60, 10) in result["chr2"]

    def test_empty_file(self, tmp_path):
        bg = tmp_path / "empty.bedgraph"
        bg.write_text("")
        result = load_bedgraph(str(bg))
        assert result == {}

    def test_comments_and_track_lines_skipped(self, tmp_path):
        bg = tmp_path / "track.bedgraph"
        bg.write_text(
            "track type=bedGraph\n"
            "#comment\n"
            "chr1\t100\t200\t7\n"
        )
        result = load_bedgraph(str(bg))
        assert len(result["chr1"]) == 1
        assert result["chr1"][0] == (100, 200, 7)


# ---------------------------------------------------------------------------
# Unit tests: load_discovery_bed
# ---------------------------------------------------------------------------

class TestLoadDiscoveryBed:
    def test_basic_load(self, tmp_path):
        bed = tmp_path / "discovery.bed"
        _write_bed_file(str(bed), [
            ("chr1", 100, 200),
            ("chr2", 500, 600),
        ])
        result = load_discovery_bed(str(bed))
        assert (100, 200) in result["chr1"]
        assert (500, 600) in result["chr2"]

    def test_empty_file(self, tmp_path):
        bed = tmp_path / "empty.bed"
        bed.write_text("")
        result = load_discovery_bed(str(bed))
        assert result == {}

    def test_comment_lines_skipped(self, tmp_path):
        bed = tmp_path / "commented.bed"
        bed.write_text("#header\nchr1\t0\t100\t.\t.\t.\t.\t.\t.\t.\n")
        result = load_discovery_bed(str(bed))
        assert result["chr1"] == [(0, 100)]


# ---------------------------------------------------------------------------
# Unit tests: load_vcf_variants
# ---------------------------------------------------------------------------

class TestLoadVcfVariants:
    def test_basic_load(self, tmp_path):
        vcf_path = str(tmp_path / "test.vcf.gz")
        _write_vcf_gz(vcf_path, [
            {"chrom": "chr1", "pos": 101, "ref": "A", "alt": "T",
             "dku": 5, "dka": 5},
            {"chrom": "chr2", "pos": 501, "ref": "G", "alt": "C",
             "dku": 0, "dka": 0},
        ])
        variants = load_vcf_variants(vcf_path)
        assert len(variants) == 2
        assert variants[0]["chrom"] == "chr1"
        assert variants[0]["pos1"] == 101
        assert variants[0]["pos0"] == 100  # 0-based
        assert variants[0]["dku"] == 5
        assert variants[1]["chrom"] == "chr2"

    def test_empty_vcf(self, tmp_path):
        vcf_path = str(tmp_path / "empty.vcf.gz")
        _write_vcf_gz(vcf_path, [])
        variants = load_vcf_variants(vcf_path)
        assert variants == []


# ---------------------------------------------------------------------------
# Unit tests: compare()
# ---------------------------------------------------------------------------

class TestCompare:
    """Unit tests for the compare() classification function."""

    def _make_bedgraph(self, intervals):
        """Dict from [(chrom, start, end, count)]."""
        bg = {}
        for chrom, start, end, count in intervals:
            bg.setdefault(chrom, []).append((start, end, count))
        return bg

    def _make_discovery(self, regions):
        """Dict from [(chrom, start, end)]."""
        disc = {}
        for chrom, start, end in regions:
            disc.setdefault(chrom, []).append((start, end))
        return disc

    def _make_variants(self, records):
        """List of variant dicts from [(chrom, pos1, ref, alt, dku, dka)]."""
        return [
            {"chrom": c, "pos0": p - 1, "pos1": p,
             "ref": r, "alt": a, "dku": dku, "dka": dka}
            for c, p, r, a, dku, dka in records
        ]

    def test_concordant(self):
        """Variant with signal AND discovery region is concordant."""
        # bedgraph covers pos 100-200, discovery covers 50-300, variant at pos 150
        bg = self._make_bedgraph([("chr1", 100, 200, 5)])
        disc = self._make_discovery([("chr1", 50, 300)])
        variants = self._make_variants([("chr1", 151, "A", "T", 3, 3)])

        result = compare(bg, disc, variants)
        assert len(result["concordant"]) == 1
        assert len(result["vcf_only"]) == 0
        assert len(result["no_signal"]) == 0
        assert len(result["discovery_only"]) == 0

    def test_no_signal(self):
        """Variant with no bedGraph signal goes to no_signal."""
        bg = {}  # no signal
        disc = self._make_discovery([("chr1", 50, 300)])
        variants = self._make_variants([("chr1", 151, "A", "T", 0, 0)])

        result = compare(bg, disc, variants)
        assert len(result["no_signal"]) == 1
        assert result["no_signal"][0]["has_discovery"] is True
        assert len(result["concordant"]) == 0

    def test_vcf_only(self):
        """Variant with signal but no discovery region is vcf_only."""
        bg = self._make_bedgraph([("chr1", 100, 200, 5)])
        disc = {}  # no discovery regions
        variants = self._make_variants([("chr1", 151, "A", "T", 3, 3)])

        result = compare(bg, disc, variants)
        assert len(result["vcf_only"]) == 1
        assert len(result["concordant"]) == 0

    def test_discovery_only(self):
        """Discovery region with no VCF variant inside is discovery_only."""
        bg = self._make_bedgraph([("chr1", 500, 600, 3)])
        disc = self._make_discovery([("chr1", 500, 600)])
        variants = self._make_variants([("chr1", 151, "A", "T", 0, 0)])

        result = compare(bg, disc, variants)
        assert len(result["discovery_only"]) == 1
        assert result["discovery_only"][0] == {
            "chrom": "chr1", "start": 500, "end": 600,
        }

    def test_window_expands_overlap(self):
        """A window parameter extends the overlap check around the VCF position."""
        # signal at 200-300, variant at 150 — exact check fails
        bg = self._make_bedgraph([("chr1", 200, 300, 5)])
        disc = self._make_discovery([("chr1", 100, 400)])
        variants = self._make_variants([("chr1", 151, "A", "T", 3, 3)])

        # Without window: no signal
        result_no_win = compare(bg, disc, variants, window=0)
        assert len(result_no_win["no_signal"]) == 1

        # With window=60: pos0=150, search 90-211 → hits [200,300)
        result_win = compare(bg, disc, variants, window=60)
        assert len(result_win["concordant"]) == 1

    def test_empty_inputs(self):
        """All-empty inputs produce empty results."""
        result = compare({}, {}, [])
        assert result["concordant"] == []
        assert result["vcf_only"] == []
        assert result["no_signal"] == []
        assert result["discovery_only"] == []

    def test_multiple_variants_and_regions(self):
        """Multiple variants across chromosomes are classified correctly."""
        bg = self._make_bedgraph([
            ("chr1", 100, 200, 5),  # covers variant A
            ("chr2", 300, 400, 2),  # covers variant C
        ])
        disc = self._make_discovery([
            ("chr1", 50, 250),    # covers variant A
            ("chr2", 280, 420),   # covers variant C
            ("chr3", 1000, 2000), # no VCF variant
        ])
        variants = self._make_variants([
            ("chr1", 151, "A", "T", 3, 3),  # A: concordant
            ("chr1", 500, "G", "C", 0, 0),  # B: no signal
            ("chr2", 351, "T", "A", 2, 2),  # C: concordant
        ])
        result = compare(bg, disc, variants)
        assert len(result["concordant"]) == 2
        assert len(result["no_signal"]) == 1
        assert len(result["discovery_only"]) == 1
        assert result["discovery_only"][0]["chrom"] == "chr3"


# ---------------------------------------------------------------------------
# Unit tests: format_summary
# ---------------------------------------------------------------------------

class TestFormatSummary:
    def _minimal_result(self):
        return {
            "concordant": [],
            "vcf_only": [],
            "no_signal": [],
            "discovery_only": [],
        }

    def test_section_headers_present(self):
        summary = format_summary(self._minimal_result())
        assert "CONCORDANT" in summary
        assert "VCF_ONLY" in summary
        assert "NO_SIGNAL" in summary
        assert "DISCOVERY_ONLY" in summary
        assert "Summary" in summary

    def test_concordant_variant_appears(self):
        result = self._minimal_result()
        result["concordant"].append({
            "variant": {
                "chrom": "chr1", "pos1": 101, "ref": "A", "alt": "T",
                "dku": 5, "dka": 4,
            },
            "regions": [(50, 200)],
        })
        summary = format_summary(result)
        assert "chr1:101 A>T" in summary
        assert "DKU=5" in summary

    def test_discovery_only_region_appears(self):
        result = self._minimal_result()
        result["discovery_only"].append(
            {"chrom": "chr2", "start": 500, "end": 800},
        )
        summary = format_summary(result)
        assert "chr2:500-800" in summary

    def test_window_annotation(self):
        summary = format_summary(self._minimal_result(), window=50)
        assert "±50 bp" in summary

    def test_counts_in_summary(self):
        result = self._minimal_result()
        result["concordant"].append({
            "variant": {
                "chrom": "chr1", "pos1": 100, "ref": "A", "alt": "T",
                "dku": 1, "dka": 1,
            },
            "regions": [(50, 200)],
        })
        result["discovery_only"].append(
            {"chrom": "chr2", "start": 0, "end": 100},
        )
        summary = format_summary(result)
        assert "Total VCF variants:            1" in summary
        assert "Concordant (signal + region):  1" in summary
        assert "Discovery-only regions:        1" in summary


# ---------------------------------------------------------------------------
# Unit tests: parse_args
# ---------------------------------------------------------------------------

class TestParseArgs:
    BASE = [
        "--bedgraph", "bg.bedgraph",
        "--discovery", "disc.bed",
        "--vcf", "in.vcf.gz",
    ]

    def test_required_args(self):
        args = parse_args(self.BASE)
        assert args.bedgraph == "bg.bedgraph"
        assert args.discovery == "disc.bed"
        assert args.vcf == "in.vcf.gz"

    def test_defaults(self):
        args = parse_args(self.BASE)
        assert args.output is None
        assert args.window == 0

    def test_output_flag(self):
        args = parse_args(self.BASE + ["--output", "out.txt"])
        assert args.output == "out.txt"

    def test_window_flag(self):
        args = parse_args(self.BASE + ["--window", "100"])
        assert args.window == 100

    def test_short_flags(self):
        args = parse_args([
            "-b", "bg.bedgraph", "-d", "disc.bed", "-v", "in.vcf.gz",
            "-w", "50",
        ])
        assert args.bedgraph == "bg.bedgraph"
        assert args.window == 50


# ---------------------------------------------------------------------------
# Unit tests: main() (CLI entry point)
# ---------------------------------------------------------------------------

class TestMain:
    def test_main_runs_and_writes_output(self, tmp_path):
        """main() should write a summary file when --output is given."""
        # Bedgraph
        bg_path = str(tmp_path / "cov.bedgraph")
        _write_bedgraph_file(bg_path, [("chr1", 100, 200, 5)])

        # Discovery BED
        bed_path = str(tmp_path / "disc.bed")
        _write_bed_file(bed_path, [("chr1", 50, 300)])

        # VCF
        vcf_path = str(tmp_path / "in.vcf.gz")
        _write_vcf_gz(vcf_path, [
            {"chrom": "chr1", "pos": 151, "ref": "A", "alt": "T",
             "dku": 3, "dka": 3},
        ])

        out_path = str(tmp_path / "summary.txt")
        main([
            "--bedgraph", bg_path,
            "--discovery", bed_path,
            "--vcf", vcf_path,
            "--output", out_path,
        ])

        assert os.path.isfile(out_path)
        with open(out_path) as fh:
            content = fh.read()
        assert "CONCORDANT" in content
        assert "chr1:151 A>T" in content

    def test_main_no_output_file(self, tmp_path, capsys):
        """main() without --output should print to stdout."""
        bg_path = str(tmp_path / "cov.bedgraph")
        _write_bedgraph_file(bg_path, [])

        bed_path = str(tmp_path / "disc.bed")
        _write_bed_file(bed_path, [])

        vcf_path = str(tmp_path / "in.vcf.gz")
        _write_vcf_gz(vcf_path, [])

        main(["--bedgraph", bg_path, "--discovery", bed_path, "--vcf", vcf_path])
        captured = capsys.readouterr()
        assert "Summary" in captured.out


# ---------------------------------------------------------------------------
# Integration test: GIAB example files
# ---------------------------------------------------------------------------

@pytest.mark.skipif(
    not GIAB_DATA_EXISTS,
    reason="GIAB example output not available",
)
class TestGiabIntegration:
    """Smoke-test compare_regions against the committed GIAB example files."""

    def test_giab_compare_produces_results(self):
        bedgraph = os.path.join(
            EXAMPLE_DISCOVERY_DIR, "giab_discovery.kmer_coverage.bedgraph",
        )
        discovery = os.path.join(
            EXAMPLE_DISCOVERY_DIR, "giab_discovery.bed",
        )
        vcf = os.path.join(EXAMPLE_OUTPUT_DIR, "annotated.vcf.gz")

        bg = load_bedgraph(bedgraph)
        disc = load_discovery_bed(discovery)
        variants = load_vcf_variants(vcf)
        result = compare(bg, disc, variants)

        total = (
            len(result["concordant"])
            + len(result["vcf_only"])
            + len(result["no_signal"])
        )
        assert total == len(variants), "Every VCF variant must be classified"
        # At least some concordant regions expected for GIAB de novo calls
        assert len(result["concordant"]) >= 1

    def test_giab_summary_format(self):
        bedgraph = os.path.join(
            EXAMPLE_DISCOVERY_DIR, "giab_discovery.kmer_coverage.bedgraph",
        )
        discovery = os.path.join(
            EXAMPLE_DISCOVERY_DIR, "giab_discovery.bed",
        )
        vcf = os.path.join(EXAMPLE_OUTPUT_DIR, "annotated.vcf.gz")

        bg = load_bedgraph(bedgraph)
        disc = load_discovery_bed(discovery)
        variants = load_vcf_variants(vcf)
        result = compare(bg, disc, variants)
        summary = format_summary(result)

        assert "CONCORDANT" in summary
        assert "DISCOVERY_ONLY" in summary
        assert "Total VCF variants:" in summary

    def test_giab_main_writes_output(self, tmp_path):
        out_path = str(tmp_path / "giab_comparison.txt")
        main([
            "--bedgraph", os.path.join(
                EXAMPLE_DISCOVERY_DIR, "giab_discovery.kmer_coverage.bedgraph",
            ),
            "--discovery", os.path.join(
                EXAMPLE_DISCOVERY_DIR, "giab_discovery.bed",
            ),
            "--vcf", os.path.join(EXAMPLE_OUTPUT_DIR, "annotated.vcf.gz"),
            "--output", out_path,
        ])
        assert os.path.isfile(out_path)
        with open(out_path) as fh:
            content = fh.read()
        assert "Summary" in content
        assert "Total VCF variants:" in content
