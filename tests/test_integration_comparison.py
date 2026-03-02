"""Integration tests comparing VCF-mode candidates to discovery-mode regions.

These tests verify that high-quality de novo candidates identified via
the VCF-based pipeline (DKA_DKT > 0.25, DKA > 10) are captured within
the genomic regions discovered by the VCF-free discovery pipeline.

Additionally, the curated de novo mutation (DNM) regions from Sulovari
et al. 2023 (PMC10006329) are evaluated against the discovery output to
verify that all known SV-like DNMs are nominated and characterised by
the VCF-free pipeline.
"""

import json
import os

import pytest

from kmer_denovo_filter.pipeline import (
    SULOVARI_DNM_REGIONS,
    _compare_candidates_to_regions,
    _evaluate_dnm_regions,
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

        # From the summary, three variants meet both criteria:
        #   chr8:125785997 A>AGAG…  DKA=24 DKA_DKT=0.3000  (SV-like insertion)
        #   chr11:55003995 T>C      DKA=24 DKA_DKT=0.4528
        #   chr11:55008577 C>T      DKA=13 DKA_DKT=0.3824
        assert len(candidates) == 3

        chroms = {c["chrom"] for c in candidates}
        assert chroms == {"chr8", "chr11"}

        positions = {c["pos"] for c in candidates}
        assert positions == {125785997, 55003995, 55008577}

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
                if line.startswith("#"):
                    continue
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
        assert comp["hq_candidates"] == 3
        assert comp["captured"] == 3
        assert comp["capture_rate"] == 1.0
        assert len(comp["candidates"]) == 3

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


class TestEvaluateDNMRegions:
    """Unit tests for _evaluate_dnm_regions()."""

    def test_all_detected_with_real_data(self):
        """All curated DNM regions should be detected in discovery output."""
        metrics_path = os.path.join(
            EXAMPLE_OUTPUT_DISCOVERY_DIR, "giab_discovery.metrics.json",
        )
        with open(metrics_path) as fh:
            metrics = json.load(fh)

        bed_path = os.path.join(
            EXAMPLE_OUTPUT_DISCOVERY_DIR, "giab_discovery.bed",
        )
        regions = []
        with open(bed_path) as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                parts = line.strip().split("\t")
                regions.append((parts[0], int(parts[1]), int(parts[2])))

        results = _evaluate_dnm_regions(regions, metrics["regions"])
        assert len(results) == 7
        for r in results:
            assert r["detected"], (
                f"DNM locus {r['locus']} ({r['event_type']}) "
                f"not detected by discovery"
            )
            assert r["assessment"] == "DETECTED"
            assert r["total_reads"] > 0
            assert r["total_unique_kmers"] > 0

    def test_region_overlap_point_event(self):
        """A point-event DNM inside a region should be detected."""
        dnm = [("chr1", 100, None, "sv_like")]
        regions = [("chr1", 50, 200)]
        detail = [{"chrom": "chr1", "start": 50, "end": 200, "size": 150,
                    "reads": 5, "unique_kmers": 10, "split_reads": 0,
                    "discordant_pairs": 0, "max_clip_len": 20,
                    "unmapped_mates": 0, "class": "SMALL"}]
        results = _evaluate_dnm_regions(regions, detail, dnm_regions=dnm)
        assert results[0]["detected"] is True
        assert results[0]["kmer_signal"] > 0

    def test_no_overlap(self):
        """A DNM with no overlapping discovery region → NOT_DETECTED."""
        dnm = [("chr1", 5000, 50, "deletion")]
        regions = [("chr1", 50, 200)]
        detail = [{"chrom": "chr1", "start": 50, "end": 200, "size": 150,
                    "reads": 5, "unique_kmers": 10, "split_reads": 0,
                    "discordant_pairs": 0, "max_clip_len": 20,
                    "unmapped_mates": 0, "class": "SMALL"}]
        results = _evaluate_dnm_regions(regions, detail, dnm_regions=dnm)
        assert results[0]["detected"] is False
        assert results[0]["assessment"] == "NOT_DETECTED"
        assert results[0]["kmer_signal"] == 0.0

    def test_adjacent_not_overlapping(self):
        """A DNM at region end boundary (half-open) → NOT_DETECTED."""
        # Region [50, 200) and DNM at 200 → adjacent, not overlapping
        dnm = [("chr1", 200, None, "sv_like")]
        regions = [("chr1", 50, 200)]
        detail = [{"chrom": "chr1", "start": 50, "end": 200, "size": 150,
                    "reads": 5, "unique_kmers": 10, "split_reads": 0,
                    "discordant_pairs": 0, "max_clip_len": 20,
                    "unmapped_mates": 0, "class": "SMALL"}]
        results = _evaluate_dnm_regions(regions, detail, dnm_regions=dnm)
        assert results[0]["detected"] is False

    def test_multi_region_overlap(self):
        """A large deletion spanning multiple regions aggregates evidence."""
        dnm = [("chr1", 100, 500, "deletion")]
        regions = [("chr1", 50, 200), ("chr1", 300, 500)]
        detail = [
            {"chrom": "chr1", "start": 50, "end": 200, "size": 150,
             "reads": 5, "unique_kmers": 10, "split_reads": 1,
             "discordant_pairs": 0, "max_clip_len": 40,
             "unmapped_mates": 2, "class": "SV"},
            {"chrom": "chr1", "start": 300, "end": 500, "size": 200,
             "reads": 3, "unique_kmers": 8, "split_reads": 0,
             "discordant_pairs": 1, "max_clip_len": 30,
             "unmapped_mates": 1, "class": "AMBIGUOUS"},
        ]
        results = _evaluate_dnm_regions(regions, detail, dnm_regions=dnm)
        assert results[0]["detected"] is True
        assert len(results[0]["discovery_regions"]) == 2
        assert results[0]["total_reads"] == 8
        assert results[0]["total_unique_kmers"] == 18
        assert results[0]["sv_class"] == "SV"

    def test_sv_class_priority(self):
        """Most severe SV class should be reported."""
        dnm = [("chr1", 100, 400, "deletion")]
        regions = [("chr1", 50, 200), ("chr1", 300, 450)]
        detail = [
            {"chrom": "chr1", "start": 50, "end": 200, "size": 150,
             "reads": 5, "unique_kmers": 10, "split_reads": 0,
             "discordant_pairs": 0, "max_clip_len": 0,
             "unmapped_mates": 0, "class": "SMALL"},
            {"chrom": "chr1", "start": 300, "end": 450, "size": 150,
             "reads": 3, "unique_kmers": 8, "split_reads": 0,
             "discordant_pairs": 0, "max_clip_len": 0,
             "unmapped_mates": 1, "class": "AMBIGUOUS"},
        ]
        results = _evaluate_dnm_regions(regions, detail, dnm_regions=dnm)
        assert results[0]["sv_class"] == "AMBIGUOUS"

    def test_result_fields(self):
        """Each evaluation result should have all expected fields."""
        metrics_path = os.path.join(
            EXAMPLE_OUTPUT_DISCOVERY_DIR, "giab_discovery.metrics.json",
        )
        with open(metrics_path) as fh:
            metrics = json.load(fh)

        bed_path = os.path.join(
            EXAMPLE_OUTPUT_DISCOVERY_DIR, "giab_discovery.bed",
        )
        regions = []
        with open(bed_path) as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                parts = line.strip().split("\t")
                regions.append((parts[0], int(parts[1]), int(parts[2])))

        results = _evaluate_dnm_regions(regions, metrics["regions"])
        expected_fields = {
            "locus", "event_type", "event_size", "detected",
            "discovery_regions", "total_reads", "total_unique_kmers",
            "max_clip_len", "unmapped_mates", "discordant_pairs",
            "split_reads", "sv_class", "kmer_signal", "assessment",
        }
        for r in results:
            assert set(r.keys()) == expected_fields

    def test_chr7_large_deletion_multi_region(self):
        """chr7 TRB locus 10.6kb deletion should span multiple regions."""
        metrics_path = os.path.join(
            EXAMPLE_OUTPUT_DISCOVERY_DIR, "giab_discovery.metrics.json",
        )
        with open(metrics_path) as fh:
            metrics = json.load(fh)

        bed_path = os.path.join(
            EXAMPLE_OUTPUT_DISCOVERY_DIR, "giab_discovery.bed",
        )
        regions = []
        with open(bed_path) as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                parts = line.strip().split("\t")
                regions.append((parts[0], int(parts[1]), int(parts[2])))

        results = _evaluate_dnm_regions(regions, metrics["regions"])
        chr7_del = [r for r in results if r["locus"] == "chr7:142786222"]
        assert len(chr7_del) == 1
        result = chr7_del[0]
        assert result["detected"] is True
        assert len(result["discovery_regions"]) >= 2
        assert result["sv_class"] == "SV"
        assert result["total_unique_kmers"] > 50

    def test_kmer_signal_per_locus(self):
        """Each detected locus should have positive k-mer signal density."""
        metrics_path = os.path.join(
            EXAMPLE_OUTPUT_DISCOVERY_DIR, "giab_discovery.metrics.json",
        )
        with open(metrics_path) as fh:
            metrics = json.load(fh)

        bed_path = os.path.join(
            EXAMPLE_OUTPUT_DISCOVERY_DIR, "giab_discovery.bed",
        )
        regions = []
        with open(bed_path) as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                parts = line.strip().split("\t")
                regions.append((parts[0], int(parts[1]), int(parts[2])))

        results = _evaluate_dnm_regions(regions, metrics["regions"])
        for r in results:
            if r["detected"]:
                assert r["kmer_signal"] > 0, (
                    f"{r['locus']} detected but kmer_signal is 0"
                )


@pytest.mark.skipif(
    not GIAB_DISCOVERY_DATA_EXISTS,
    reason="GIAB discovery test data not available",
)
class TestDNMRegionIntegration:
    """Integration tests for curated DNM region evaluation.

    These tests use the full generated discovery output to verify that
    all seven Sulovari et al. 2023 curated DNM regions are nominated
    by the VCF-free discovery pipeline.
    """

    def test_metrics_contains_dnm_evaluation(
        self, generated_discovery_output,
    ):
        """Discovery metrics should include dnm_evaluation."""
        with open(generated_discovery_output["metrics"]) as fh:
            metrics = json.load(fh)

        assert "dnm_evaluation" in metrics
        dnm_eval = metrics["dnm_evaluation"]
        assert dnm_eval["total_loci"] == 7
        assert dnm_eval["detected"] == 7
        assert dnm_eval["detection_rate"] == 1.0
        assert len(dnm_eval["loci"]) == 7

    def test_summary_contains_dnm_evaluation(
        self, generated_discovery_output,
    ):
        """Discovery summary should include the DNM evaluation section."""
        with open(generated_discovery_output["summary"]) as fh:
            text = fh.read()

        assert "Curated DNM Region Evaluation" in text
        assert "Sulovari et al. 2023" in text
        assert "Detected by discovery" in text
        assert "7 / 7 (100.0%)" in text

    def test_all_dnm_loci_detected(
        self, generated_discovery_output,
    ):
        """All 7 curated DNM loci must be detected."""
        with open(generated_discovery_output["metrics"]) as fh:
            metrics = json.load(fh)

        for locus in metrics["dnm_evaluation"]["loci"]:
            assert locus["detected"], (
                f"Curated DNM {locus['locus']} ({locus['event_type']}) "
                f"not detected"
            )
            assert locus["total_reads"] > 0
            assert locus["total_unique_kmers"] > 0

    def test_dnm_discovery_regions_are_valid(
        self, generated_discovery_output,
    ):
        """Each DNM region label should correspond to a real BED region."""
        with open(generated_discovery_output["metrics"]) as fh:
            metrics = json.load(fh)

        region_labels = set()
        for r in metrics["regions"]:
            label = f"{r['chrom']}:{r['start'] + 1}-{r['end']}"
            region_labels.add(label)

        for locus in metrics["dnm_evaluation"]["loci"]:
            for dr_label in locus["discovery_regions"]:
                assert dr_label in region_labels, (
                    f"DNM region {dr_label} for {locus['locus']} "
                    f"not found in BED regions"
                )

    def test_chr17_deletion_evidence(
        self, generated_discovery_output,
    ):
        """chr17:53340465 107bp deletion should have strong k-mer signal."""
        with open(generated_discovery_output["metrics"]) as fh:
            metrics = json.load(fh)

        loci_by_name = {
            l["locus"]: l for l in metrics["dnm_evaluation"]["loci"]
        }
        locus = loci_by_name["chr17:53340465"]
        assert locus["event_type"] == "deletion"
        assert locus["total_reads"] >= 10
        assert locus["total_unique_kmers"] >= 20
        assert locus["max_clip_len"] >= 100

    def test_chr8_sv_like_strong_signal(
        self, generated_discovery_output,
    ):
        """chr8:125785998 43bp SV-like event has strongest k-mer support."""
        with open(generated_discovery_output["metrics"]) as fh:
            metrics = json.load(fh)

        loci_by_name = {
            l["locus"]: l for l in metrics["dnm_evaluation"]["loci"]
        }
        locus = loci_by_name["chr8:125785998"]
        assert locus["total_reads"] >= 30
        assert locus["total_unique_kmers"] >= 40
        assert locus["kmer_signal"] >= 0.04
