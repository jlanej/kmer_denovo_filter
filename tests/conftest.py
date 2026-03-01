"""Shared fixtures for kmer_denovo_filter tests."""

import os
import tempfile

import pytest

from kmer_denovo_filter.cli import parse_args
from kmer_denovo_filter.pipeline import run_pipeline, run_discovery_pipeline

GIAB_DIR = os.path.join(os.path.dirname(__file__), "data", "giab")
GIAB_DATA_EXISTS = os.path.isfile(os.path.join(GIAB_DIR, "HG002_child.bam"))
GIAB_DISCOVERY_DATA_EXISTS = (
    GIAB_DATA_EXISTS
    and os.path.isfile(os.path.join(GIAB_DIR, "mini_ref.fa"))
    and os.path.isfile(os.path.join(GIAB_DIR, "mini_ref.fa.k31.jf"))
)
EXAMPLE_OUTPUT_DIR = os.path.join(os.path.dirname(__file__), "example_output")
EXAMPLE_OUTPUT_DISCOVERY_DIR = os.path.join(
    os.path.dirname(__file__), "example_output_discovery",
)


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


@pytest.fixture(scope="session")
def generated_discovery_output():
    """Run the GIAB discovery pipeline once and return output paths.

    This session-scoped fixture generates discovery-mode output files
    (BED, metrics JSON, summary text, informative BAM) by running the
    pipeline against the GIAB HG002 trio test data with a mini reference
    built from perfect-match reads.

    Returns a dict mapping output names to their file paths::

        {
            "bed":     "<tmpdir>/giab_discovery.bed",
            "metrics": "<tmpdir>/giab_discovery.metrics.json",
            "summary": "<tmpdir>/giab_discovery.summary.txt",
            "bam":     "<tmpdir>/giab_discovery.informative.bam",
            "bam_bai": "<tmpdir>/giab_discovery.informative.bam.bai",
        }
    """
    if not GIAB_DISCOVERY_DATA_EXISTS:
        pytest.skip("GIAB discovery test data not available")

    tmpdir = tempfile.mkdtemp(prefix="kmer_discovery_output_")
    out_prefix = os.path.join(tmpdir, "giab_discovery")

    args = parse_args([
        "--child", os.path.join(GIAB_DIR, "HG002_child.bam"),
        "--mother", os.path.join(GIAB_DIR, "HG004_mother.bam"),
        "--father", os.path.join(GIAB_DIR, "HG003_father.bam"),
        "--ref-fasta", os.path.join(GIAB_DIR, "mini_ref.fa"),
        "--ref-jf", os.path.join(GIAB_DIR, "mini_ref.fa.k31.jf"),
        "--out-prefix", out_prefix,
        "--min-child-count", "3",
        "--kmer-size", "31",
        "--candidate-summary", os.path.join(
            EXAMPLE_OUTPUT_DIR, "summary.txt",
        ),
    ])
    run_discovery_pipeline(args)

    return {
        "bed": f"{out_prefix}.bed",
        "metrics": f"{out_prefix}.metrics.json",
        "summary": f"{out_prefix}.summary.txt",
        "bam": f"{out_prefix}.informative.bam",
        "bam_bai": f"{out_prefix}.informative.bam.bai",
        "bedpe": f"{out_prefix}.sv.bedpe",
    }
