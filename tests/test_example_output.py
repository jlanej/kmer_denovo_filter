"""Regression tests that fail when example output changes.

These tests compare freshly-generated pipeline output (via the
``generated_example_output`` session fixture) against the committed
reference files in ``tests/example_output/``.  When a test fails it
prints a unified diff (or structured comparison) so the exact change
is immediately visible.
"""

import difflib
import json
import os

import pysam
import pytest

EXAMPLE_OUTPUT_DIR = os.path.join(os.path.dirname(__file__), "example_output")
GIAB_DIR = os.path.join(os.path.dirname(__file__), "data", "giab")
GIAB_DATA_EXISTS = os.path.isfile(os.path.join(GIAB_DIR, "HG002_child.bam"))


ANNOTATION_IDS = {
    "DKU", "DKT", "DKA", "DKU_DKT", "DKA_DKT",
    "MAX_PKC", "AVG_PKC", "MIN_PKC",
    "MAX_PKC_ALT", "AVG_PKC_ALT", "MIN_PKC_ALT",
}


def _unified_diff(expected_lines, actual_lines, label):
    """Return a unified diff string between two line lists."""
    return "\n".join(difflib.unified_diff(
        expected_lines, actual_lines,
        fromfile=f"expected/{label}",
        tofile=f"generated/{label}",
        lineterm="",
    ))


def _vcf_data_lines(path):
    """Extract non-header data lines from a bgzipped VCF."""
    lines = []
    with pysam.VariantFile(path) as vcf:
        for rec in vcf:
            lines.append(str(rec).rstrip())
    return lines


def _vcf_header_definitions(path):
    """Extract FORMAT/INFO definition lines from a VCF header.

    Returns only the lines added by our tool (DKU, DKT, DKA, etc.),
    sorted for stable comparison.
    """
    defs = []
    with pysam.VariantFile(path) as vcf:
        for rec in vcf.header.records:
            if rec.type == "FORMAT" or rec.type == "INFO":
                if rec.get("ID", "") in ANNOTATION_IDS:
                    defs.append(str(rec))
    return sorted(defs)


def _vcf_sample_fields(path):
    """Extract per-sample annotation values for stable comparison.

    Returns a list of dicts, one per variant, keyed by annotation field.
    """
    results = []
    with pysam.VariantFile(path) as vcf:
        for rec in vcf:
            variant_id = f"{rec.chrom}:{rec.pos} {rec.ref}>{','.join(str(a) for a in rec.alts)}"
            row = {"variant": variant_id}
            sample = rec.samples["HG002"]
            for f in sorted(ANNOTATION_IDS):
                val = sample.get(f)
                row[f] = val
            results.append(row)
    return results


@pytest.mark.skipif(
    not GIAB_DATA_EXISTS,
    reason="GIAB test data not available",
)
class TestExampleOutput:
    """Regression tests comparing generated output against committed examples."""

    def test_metrics_json_matches(self, generated_example_output):
        """metrics.json must match the committed example exactly."""
        expected_path = os.path.join(EXAMPLE_OUTPUT_DIR, "metrics.json")
        generated_path = generated_example_output["metrics"]

        with open(expected_path) as fh:
            expected = json.load(fh)
        with open(generated_path) as fh:
            generated = json.load(fh)

        if expected != generated:
            exp_lines = json.dumps(expected, indent=2, sort_keys=True).splitlines()
            gen_lines = json.dumps(generated, indent=2, sort_keys=True).splitlines()
            diff = _unified_diff(exp_lines, gen_lines, "metrics.json")
            pytest.fail(
                f"metrics.json differs from expected:\n{diff}"
            )

    def test_summary_text_matches(self, generated_example_output):
        """summary.txt must match the committed example exactly."""
        expected_path = os.path.join(EXAMPLE_OUTPUT_DIR, "summary.txt")
        generated_path = generated_example_output["summary"]

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

    def test_vcf_variant_count_matches(self, generated_example_output):
        """The number of variant records must match."""
        expected_path = os.path.join(EXAMPLE_OUTPUT_DIR, "annotated.vcf.gz")
        generated_path = generated_example_output["vcf"]

        expected_lines = _vcf_data_lines(expected_path)
        generated_lines = _vcf_data_lines(generated_path)

        if len(expected_lines) != len(generated_lines):
            pytest.fail(
                f"Variant count differs: expected {len(expected_lines)}, "
                f"got {len(generated_lines)}"
            )

    def test_vcf_annotations_match(self, generated_example_output):
        """Per-variant annotation values must match the committed VCF."""
        expected_path = os.path.join(EXAMPLE_OUTPUT_DIR, "annotated.vcf.gz")
        generated_path = generated_example_output["vcf"]

        expected = _vcf_sample_fields(expected_path)
        generated = _vcf_sample_fields(generated_path)

        mismatches = []
        for exp, gen in zip(expected, generated):
            for key in exp:
                if exp[key] != gen.get(key):
                    mismatches.append(
                        f"  Variant {exp['variant']}: "
                        f"{key} expected={exp[key]!r} got={gen.get(key)!r}"
                    )

        if mismatches:
            detail = "\n".join(mismatches)
            pytest.fail(
                f"VCF annotation values differ from expected "
                f"({len(mismatches)} field(s)):\n{detail}"
            )

    def test_vcf_header_definitions_match(self, generated_example_output):
        """FORMAT/INFO definitions added by our tool must match."""
        expected_path = os.path.join(EXAMPLE_OUTPUT_DIR, "annotated.vcf.gz")
        generated_path = generated_example_output["vcf"]

        expected_defs = _vcf_header_definitions(expected_path)
        generated_defs = _vcf_header_definitions(generated_path)

        if expected_defs != generated_defs:
            diff = _unified_diff(
                expected_defs, generated_defs, "vcf_header_definitions",
            )
            pytest.fail(
                f"VCF header definitions differ:\n{diff}"
            )

    def test_vcf_data_lines_match(self, generated_example_output):
        """Full VCF data lines (tab-separated) must match exactly."""
        expected_path = os.path.join(EXAMPLE_OUTPUT_DIR, "annotated.vcf.gz")
        generated_path = generated_example_output["vcf"]

        expected_lines = _vcf_data_lines(expected_path)
        generated_lines = _vcf_data_lines(generated_path)

        if expected_lines != generated_lines:
            diff = _unified_diff(
                expected_lines, generated_lines, "annotated.vcf.gz",
            )
            pytest.fail(
                f"VCF data lines differ from expected:\n{diff}"
            )
