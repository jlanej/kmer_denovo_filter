"""Integration tests comparing VCF-mode candidates to discovery-mode regions.

These tests verify that high-quality de novo candidates identified via
the VCF-based pipeline (DKA_DKT > 0.25, DKA > 10) are captured within
the genomic regions discovered by the VCF-free discovery pipeline.

As variant genotyping is refined for the discovery pipeline, these
comparisons can be extended to verify that true variation is captured.
"""

import json
import os

import pytest

from kmer_denovo_filter.pipeline import (
    _compare_candidates_to_regions,
    _parse_candidate_summary,
)

EXAMPLE_OUTPUT_DIR = os.path.join(os.path.dirname(__file__), "example_output")
EXAMPLE_OUTPUT_DISCOVERY_DIR = os.path.join(
    os.path.dirname(__file__), "example_output_discovery",
)
GIAB_DIR = os.path.join(os.path.dirname(__file__), "data", "giab")
GIAB_DISCOVERY_DATA_EXISTS = (
    os.path.isfile(os.path.join(GIAB_DIR, "HG002_child.bam"))
    and os.path.isfile(os.path.join(GIAB_DIR, "mini_ref.fa"))
    and os.path.isfile(os.path.join(GIAB_DIR, "mini_ref.fa.k31.jf"))
)


class TestParseCandidateSummary:
    """Unit tests for _parse_candidate_summary()."""

    def test_parses_high_quality_candidates(self):
        """Should extract candidates with DKA_DKT > 0.25 and DKA > 10."""
        summary_path = os.path.join(EXAMPLE_OUTPUT_DIR, "summary.txt")
        candidates = _parse_candidate_summary(summary_path)

        # From the summary, only two variants meet both criteria:
        #   chr11:55003995 T>C  DKA=24 DKA_DKT=0.4528
        #   chr11:55008577 C>T  DKA=13 DKA_DKT=0.3824
        assert len(candidates) == 2

        chroms = {c["chrom"] for c in candidates}
        assert chroms == {"chr11"}

        positions = {c["pos"] for c in candidates}
        assert positions == {55003995, 55008577}

    def test_candidate_fields(self):
        """Each candidate should have all expected fields."""
        summary_path = os.path.join(EXAMPLE_OUTPUT_DIR, "summary.txt")
        candidates = _parse_candidate_summary(summary_path)
        for cand in candidates:
            assert "chrom" in cand
            assert "pos" in cand
            assert "ref" in cand
            assert "alt" in cand
            assert "dka" in cand
            assert "dka_dkt" in cand
            assert cand["dka_dkt"] > 0.25
            assert cand["dka"] > 10

    def test_custom_thresholds(self):
        """Should respect custom DKA_DKT and DKA thresholds."""
        summary_path = os.path.join(EXAMPLE_OUTPUT_DIR, "summary.txt")
        # Very strict thresholds — only chr11:55003995 with DKA=24
        candidates = _parse_candidate_summary(
            summary_path, dka_dkt_min=0.40, dka_min=20,
        )
        assert len(candidates) == 1
        assert candidates[0]["pos"] == 55003995


class TestCompareCandidatesToRegions:
    """Unit tests for _compare_candidates_to_regions()."""

    def test_candidate_inside_region(self):
        """A candidate within a region should be marked as captured."""
        candidates = [{"chrom": "chr1", "pos": 100, "ref": "A", "alt": "T",
                        "dka": 15, "dka_dkt": 0.5, "call": "DE_NOVO"}]
        regions = [("chr1", 50, 200)]
        results = _compare_candidates_to_regions(candidates, regions)
        assert results[0]["captured"] is True
        assert results[0]["region"] == "chr1:51-200"

    def test_candidate_outside_region(self):
        """A candidate outside all regions should not be captured."""
        candidates = [{"chrom": "chr1", "pos": 300, "ref": "A", "alt": "T",
                        "dka": 15, "dka_dkt": 0.5, "call": "DE_NOVO"}]
        regions = [("chr1", 50, 200)]
        results = _compare_candidates_to_regions(candidates, regions)
        assert results[0]["captured"] is False
        assert results[0]["region"] is None

    def test_candidate_wrong_chrom(self):
        """A candidate on a different chromosome should not be captured."""
        candidates = [{"chrom": "chr2", "pos": 100, "ref": "A", "alt": "T",
                        "dka": 15, "dka_dkt": 0.5, "call": "DE_NOVO"}]
        regions = [("chr1", 50, 200)]
        results = _compare_candidates_to_regions(candidates, regions)
        assert results[0]["captured"] is False

    def test_real_data_overlap(self):
        """High-quality candidates should overlap with discovery regions."""
        summary_path = os.path.join(EXAMPLE_OUTPUT_DIR, "summary.txt")
        bed_path = os.path.join(
            EXAMPLE_OUTPUT_DISCOVERY_DIR, "giab_discovery.bed",
        )

        candidates = _parse_candidate_summary(summary_path)
        regions = []
        with open(bed_path) as fh:
            for line in fh:
                parts = line.strip().split("\t")
                regions.append((parts[0], int(parts[1]), int(parts[2])))

        results = _compare_candidates_to_regions(candidates, regions)

        # Both high-quality candidates should be captured
        for r in results:
            assert r["captured"], (
                f"High-quality candidate {r['chrom']}:{r['pos']} "
                f"{r['ref']}>{r['alt']} not captured by any discovery region"
            )


@pytest.mark.skipif(
    not GIAB_DISCOVERY_DATA_EXISTS,
    reason="GIAB discovery test data not available",
)
class TestIntegrationComparison:
    """Integration tests comparing VCF candidates to discovery regions.

    These tests use the full generated discovery output (which now
    includes --candidate-summary) to verify the comparison is present
    and correct.
    """

    def test_discovery_summary_contains_comparison(
        self, generated_discovery_output,
    ):
        """Discovery summary should include the candidate comparison."""
        with open(generated_discovery_output["summary"]) as fh:
            text = fh.read()

        assert "Candidate Comparison" in text
        assert "DKA_DKT > 0.25" in text
        assert "High-quality candidates" in text

    def test_discovery_metrics_contains_comparison(
        self, generated_discovery_output,
    ):
        """Discovery metrics should include candidate_comparison."""
        with open(generated_discovery_output["metrics"]) as fh:
            metrics = json.load(fh)

        assert "candidate_comparison" in metrics
        comp = metrics["candidate_comparison"]
        assert comp["hq_candidates"] == 2
        assert comp["captured"] == 2
        assert comp["capture_rate"] == 1.0
        assert len(comp["candidates"]) == 2

    def test_all_hq_candidates_captured(
        self, generated_discovery_output,
    ):
        """All high-quality candidates must fall within discovery regions."""
        with open(generated_discovery_output["metrics"]) as fh:
            metrics = json.load(fh)

        comp = metrics["candidate_comparison"]
        for cand in comp["candidates"]:
            assert cand["captured"], (
                f"Candidate {cand['variant']} not captured by any "
                f"discovery region"
            )
            assert cand["region"] is not None

    def test_captured_regions_are_valid(
        self, generated_discovery_output,
    ):
        """Each captured region label should correspond to a real region."""
        with open(generated_discovery_output["metrics"]) as fh:
            metrics = json.load(fh)

        region_labels = set()
        for r in metrics["regions"]:
            label = f"{r['chrom']}:{r['start'] + 1}-{r['end']}"
            region_labels.add(label)

        comp = metrics["candidate_comparison"]
        for cand in comp["candidates"]:
            if cand["captured"]:
                assert cand["region"] in region_labels, (
                    f"Region {cand['region']} not found in BED regions"
                )
