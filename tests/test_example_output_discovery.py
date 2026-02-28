"""Regression tests for discovery-mode example output.

These tests compare freshly-generated discovery pipeline output (via the
``generated_discovery_output`` session fixture) against the committed
reference files in ``tests/example_output_discovery/``.  When a test
fails it prints a unified diff so the exact change is immediately
visible.

Mirrors the pattern established in ``test_example_output.py`` for
VCF-mode output.
"""

import difflib
import json
import os

import pysam
import pytest

EXAMPLE_OUTPUT_DISCOVERY_DIR = os.path.join(
    os.path.dirname(__file__), "example_output_discovery",
)
GIAB_DIR = os.path.join(os.path.dirname(__file__), "data", "giab")
GIAB_DISCOVERY_DATA_EXISTS = (
    os.path.isfile(os.path.join(GIAB_DIR, "HG002_child.bam"))
    and os.path.isfile(os.path.join(GIAB_DIR, "mini_ref.fa"))
    and os.path.isfile(os.path.join(GIAB_DIR, "mini_ref.fa.k31.jf"))
)


def _unified_diff(expected_lines, actual_lines, label):
    """Return a unified diff string between two line lists."""
    return "\n".join(difflib.unified_diff(
        expected_lines, actual_lines,
        fromfile=f"expected/{label}",
        tofile=f"generated/{label}",
        lineterm="",
    ))


@pytest.mark.skipif(
    not GIAB_DISCOVERY_DATA_EXISTS,
    reason="GIAB discovery test data not available",
)
class TestDiscoveryExampleOutput:
    """Regression tests for discovery-mode output files."""

    def test_bed_matches(self, generated_discovery_output):
        """BED file must match the committed example exactly."""
        expected_path = os.path.join(
            EXAMPLE_OUTPUT_DISCOVERY_DIR, "giab_discovery.bed",
        )
        generated_path = generated_discovery_output["bed"]

        with open(expected_path) as fh:
            expected_lines = fh.read().splitlines()
        with open(generated_path) as fh:
            generated_lines = fh.read().splitlines()

        if expected_lines != generated_lines:
            diff = _unified_diff(
                expected_lines, generated_lines, "giab_discovery.bed",
            )
            pytest.fail(
                f"BED file differs from expected:\n{diff}"
            )

    def test_metrics_json_matches(self, generated_discovery_output):
        """metrics.json must match the committed example exactly."""
        expected_path = os.path.join(
            EXAMPLE_OUTPUT_DISCOVERY_DIR, "giab_discovery.metrics.json",
        )
        generated_path = generated_discovery_output["metrics"]

        with open(expected_path) as fh:
            expected = json.load(fh)
        with open(generated_path) as fh:
            generated = json.load(fh)

        if expected != generated:
            exp_lines = json.dumps(
                expected, indent=2, sort_keys=True,
            ).splitlines()
            gen_lines = json.dumps(
                generated, indent=2, sort_keys=True,
            ).splitlines()
            diff = _unified_diff(exp_lines, gen_lines, "metrics.json")
            pytest.fail(
                f"metrics.json differs from expected:\n{diff}"
            )

    def test_summary_text_matches(self, generated_discovery_output):
        """summary.txt must match the committed example exactly."""
        expected_path = os.path.join(
            EXAMPLE_OUTPUT_DISCOVERY_DIR, "giab_discovery.summary.txt",
        )
        generated_path = generated_discovery_output["summary"]

        with open(expected_path) as fh:
            expected_lines = fh.read().splitlines()
        with open(generated_path) as fh:
            generated_lines = fh.read().splitlines()

        if expected_lines != generated_lines:
            diff = _unified_diff(
                expected_lines, generated_lines, "summary.txt",
            )
            pytest.fail(
                f"summary.txt differs from expected:\n{diff}"
            )

    def test_informative_bam_exists(self, generated_discovery_output):
        """Informative BAM and index must be generated."""
        assert os.path.isfile(generated_discovery_output["bam"])
        assert os.path.isfile(generated_discovery_output["bam_bai"])

    def test_informative_bam_has_reads(self, generated_discovery_output):
        """Informative BAM should contain reads with dk:i:1 tag."""
        bam = pysam.AlignmentFile(generated_discovery_output["bam"])
        reads = list(bam)
        bam.close()

        assert len(reads) > 0, "Informative BAM should contain reads"
        for read in reads:
            assert read.has_tag("dk"), (
                f"Read {read.query_name} missing dk tag"
            )
            assert read.get_tag("dk") == 1

    def test_region_count_matches(self, generated_discovery_output):
        """BED region count must match metrics.json candidate_regions."""
        with open(generated_discovery_output["bed"]) as fh:
            bed_count = sum(1 for line in fh if line.strip())

        with open(generated_discovery_output["metrics"]) as fh:
            metrics = json.load(fh)

        assert bed_count == metrics["candidate_regions"], (
            f"BED has {bed_count} regions but metrics says "
            f"{metrics['candidate_regions']}"
        )

    def test_metrics_regions_detail(self, generated_discovery_output):
        """metrics.json must contain per-region detail array."""
        with open(generated_discovery_output["metrics"]) as fh:
            metrics = json.load(fh)

        assert "regions" in metrics, "metrics.json missing 'regions' array"
        assert len(metrics["regions"]) == metrics["candidate_regions"]
        for region in metrics["regions"]:
            for key in ("chrom", "start", "end", "size", "reads",
                        "unique_kmers"):
                assert key in region, f"Region missing key: {key}"
