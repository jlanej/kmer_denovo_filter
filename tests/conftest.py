"""Shared fixtures for kmer_denovo_filter tests."""

import os
import tempfile

import pytest

from kmer_denovo_filter.cli import parse_args
from kmer_denovo_filter.pipeline import run_pipeline

GIAB_DIR = os.path.join(os.path.dirname(__file__), "data", "giab")
GIAB_DATA_EXISTS = os.path.isfile(os.path.join(GIAB_DIR, "HG002_child.bam"))
EXAMPLE_OUTPUT_DIR = os.path.join(os.path.dirname(__file__), "example_output")


@pytest.fixture(scope="session")
def generated_example_output():
    """Run the GIAB pipeline once and return the output directory.

    This session-scoped fixture generates the full set of example output
    files (annotated VCF, metrics JSON, summary text) by running the
    pipeline against the GIAB HG002 trio test data.  The output is
    written to a temporary directory that is shared across all tests in
    the session, so the pipeline only runs once.

    Returns a dict mapping output names to their file paths::

        {
            "vcf":     "<tmpdir>/annotated.vcf.gz",
            "vcf_tbi": "<tmpdir>/annotated.vcf.gz.tbi",
            "metrics": "<tmpdir>/metrics.json",
            "summary": "<tmpdir>/summary.txt",
        }
    """
    if not GIAB_DATA_EXISTS:
        pytest.skip("GIAB test data not available")

    tmpdir = tempfile.mkdtemp(prefix="kmer_example_output_")

    out_vcf = os.path.join(tmpdir, "annotated.vcf.gz")
    metrics_json = os.path.join(tmpdir, "metrics.json")
    summary_txt = os.path.join(tmpdir, "summary.txt")

    args = parse_args([
        "--child", os.path.join(GIAB_DIR, "HG002_child.bam"),
        "--mother", os.path.join(GIAB_DIR, "HG004_mother.bam"),
        "--father", os.path.join(GIAB_DIR, "HG003_father.bam"),
        "--vcf", os.path.join(GIAB_DIR, "candidates.vcf.gz"),
        "--output", out_vcf,
        "--metrics", metrics_json,
        "--summary", summary_txt,
        "--proband-id", "HG002",
    ])
    run_pipeline(args)

    return {
        "vcf": out_vcf,
        "vcf_tbi": out_vcf + ".tbi",
        "metrics": metrics_json,
        "summary": summary_txt,
    }
