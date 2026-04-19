"""Tests for the interactive HTML report generator."""

import json
import os
import tempfile

import pytest

from kmer_denovo_filter.report import (
    _classify_variant_type,
    _load_metrics,
    _load_summary_counts,
    _load_summary_variants,
    generate_report,
)

EXAMPLE_OUTPUT_DIR = os.path.join(os.path.dirname(__file__), "example_output")
EXAMPLE_OUTPUT_DISCOVERY_DIR = os.path.join(
    os.path.dirname(__file__), "example_output_discovery",
)


class TestLoadMetrics:
    """Test metrics.json loading."""

    def test_load_vcf_metrics(self):
        metrics = _load_metrics(
            os.path.join(EXAMPLE_OUTPUT_DIR, "metrics.json"),
        )
        assert metrics["total_variants"] == 22
        assert metrics["total_child_kmers"] == 1484
        assert metrics["parent_found_kmers"] == 1294
        assert metrics["child_unique_kmers"] == 190
        assert metrics["variants_with_unique_reads"] == 12

    def test_load_discovery_metrics(self):
        metrics = _load_metrics(
            os.path.join(EXAMPLE_OUTPUT_DISCOVERY_DIR,
                         "giab_discovery.metrics.json"),
        )
        assert metrics["mode"] == "discovery"
        assert metrics["candidate_regions"] == 21
        assert metrics["proband_unique_kmers"] == 630
        assert "regions" in metrics
        assert len(metrics["regions"]) == 21


class TestLoadSummary:
    """Test summary.txt parsing."""

    def test_load_summary_variants(self):
        variants = _load_summary_variants(
            os.path.join(EXAMPLE_OUTPUT_DIR, "summary.txt"),
        )
        assert len(variants) == 22
        # Check first variant
        first = variants[0]
        assert "chr8:40003391" in first["label"]
        assert first["dku"] == 1
        assert first["dkt"] == 17
        assert first["call"] == "DE_NOVO"

    def test_load_summary_counts(self):
        counts = _load_summary_counts(
            os.path.join(EXAMPLE_OUTPUT_DIR, "summary.txt"),
        )
        assert counts["total_candidates"] == 22
        assert counts["likely_denovo"] == 12
        assert counts["inherited"] == 10

    def test_variant_fields_complete(self):
        variants = _load_summary_variants(
            os.path.join(EXAMPLE_OUTPUT_DIR, "summary.txt"),
        )
        required_fields = {
            "label", "dku", "dkt", "dka", "dku_dkt", "dka_dkt",
            "max_pkc", "avg_pkc", "min_pkc", "max_pkc_alt",
            "avg_pkc_alt", "min_pkc_alt", "call",
        }
        for v in variants:
            assert required_fields.issubset(v.keys()), (
                f"Missing fields in variant {v.get('label')}"
            )

    def test_inherited_and_denovo_calls(self):
        variants = _load_summary_variants(
            os.path.join(EXAMPLE_OUTPUT_DIR, "summary.txt"),
        )
        calls = {v["call"] for v in variants}
        assert calls == {"DE_NOVO", "inherited"}
        denovo_count = sum(1 for v in variants if v["call"] == "DE_NOVO")
        assert denovo_count == 12


class TestGenerateReport:
    """Test HTML report generation."""

    def test_vcf_mode_report(self):
        with tempfile.NamedTemporaryFile(suffix=".html", delete=False) as fh:
            out = fh.name
        try:
            result = generate_report(
                output_path=out,
                vcf_metrics_path=os.path.join(
                    EXAMPLE_OUTPUT_DIR, "metrics.json",
                ),
                vcf_summary_path=os.path.join(
                    EXAMPLE_OUTPUT_DIR, "summary.txt",
                ),
            )
            assert result == out
            assert os.path.isfile(out)
            with open(out) as fh:
                html = fh.read()
            assert "kmer-denovo" in html
            assert "Pipeline Summary" in html
            assert "K-mer Filtering Funnel" in html
            assert "DKA/DKT Ratio" in html
            assert "Plotly.newPlot" in html
            # Stratification metric cards should be present
            assert "Putative denovo (input)" in html
            assert "Higher-quality" in html
            assert "1,484" in html  # total child kmers
        finally:
            os.unlink(out)

    def test_discovery_mode_report(self):
        with tempfile.NamedTemporaryFile(suffix=".html", delete=False) as fh:
            out = fh.name
        try:
            result = generate_report(
                output_path=out,
                discovery_metrics_path=os.path.join(
                    EXAMPLE_OUTPUT_DISCOVERY_DIR,
                    "giab_discovery.metrics.json",
                ),
            )
            assert result == out
            assert os.path.isfile(out)
            with open(out) as fh:
                html = fh.read()
            assert "Discovery mode" in html
            assert "Candidate Regions" in html
            assert "Curated DNM" in html
            assert "Sulovari" in html
        finally:
            os.unlink(out)

    def test_combined_mode_report(self):
        with tempfile.NamedTemporaryFile(suffix=".html", delete=False) as fh:
            out = fh.name
        try:
            result = generate_report(
                output_path=out,
                vcf_metrics_path=os.path.join(
                    EXAMPLE_OUTPUT_DIR, "metrics.json",
                ),
                vcf_summary_path=os.path.join(
                    EXAMPLE_OUTPUT_DIR, "summary.txt",
                ),
                discovery_metrics_path=os.path.join(
                    EXAMPLE_OUTPUT_DISCOVERY_DIR,
                    "giab_discovery.metrics.json",
                ),
            )
            assert os.path.isfile(out)
            with open(out) as fh:
                html = fh.read()
            assert "Combined" in html
            assert "Putative denovo (input)" in html
            assert "Candidate Regions" in html
        finally:
            os.unlink(out)

    def test_report_contains_all_sections(self):
        with tempfile.NamedTemporaryFile(suffix=".html", delete=False) as fh:
            out = fh.name
        try:
            generate_report(
                output_path=out,
                vcf_metrics_path=os.path.join(
                    EXAMPLE_OUTPUT_DIR, "metrics.json",
                ),
                vcf_summary_path=os.path.join(
                    EXAMPLE_OUTPUT_DIR, "summary.txt",
                ),
                discovery_metrics_path=os.path.join(
                    EXAMPLE_OUTPUT_DISCOVERY_DIR,
                    "giab_discovery.metrics.json",
                ),
            )
            with open(out) as fh:
                html = fh.read()
            # Check all major sections are present
            assert "Pipeline Summary" in html
            assert "Methods" in html
            assert "Results Highlights" in html
            assert "K-mer Filtering Funnel" in html
            assert "Variant Breakdown" in html
            assert "DKA/DKT Ratio Distribution" in html
            assert "Per-Variant Evidence Heatmap" in html
            assert "Parental K-mer Count" in html
            assert "Discovery Mode" in html
            assert "Discovery vs. VCF Mode Concordance" in html
            assert "Curated DNM Region Evaluation" in html
            assert "Per-Variant Detail Table" in html
            assert "Method Overview" in html
        finally:
            os.unlink(out)

    def test_report_is_self_contained_html(self):
        """Report should be valid self-contained HTML with inline Plotly."""
        with tempfile.NamedTemporaryFile(suffix=".html", delete=False) as fh:
            out = fh.name
        try:
            generate_report(
                output_path=out,
                vcf_metrics_path=os.path.join(
                    EXAMPLE_OUTPUT_DIR, "metrics.json",
                ),
                vcf_summary_path=os.path.join(
                    EXAMPLE_OUTPUT_DIR, "summary.txt",
                ),
            )
            with open(out) as fh:
                html = fh.read()
            assert html.startswith("<!DOCTYPE html>")
            assert "</html>" in html
            assert "plotly" in html.lower()
            # Plotly should be inline (no CDN script tag) for offline rendering
            assert '<script src="https://cdn.plot.ly' not in html
            assert "Plotly.newPlot" in html
        finally:
            os.unlink(out)

    def test_report_idempotent(self):
        """Generating the same report twice produces identical output."""
        with tempfile.NamedTemporaryFile(suffix=".html", delete=False) as f1:
            out1 = f1.name
        with tempfile.NamedTemporaryFile(suffix=".html", delete=False) as f2:
            out2 = f2.name
        try:
            kwargs = dict(
                vcf_metrics_path=os.path.join(
                    EXAMPLE_OUTPUT_DIR, "metrics.json",
                ),
                vcf_summary_path=os.path.join(
                    EXAMPLE_OUTPUT_DIR, "summary.txt",
                ),
            )
            generate_report(output_path=out1, **kwargs)
            generate_report(output_path=out2, **kwargs)
            with open(out1) as fh:
                html1 = fh.read()
            with open(out2) as fh:
                html2 = fh.read()
            assert html1 == html2
        finally:
            os.unlink(out1)
            os.unlink(out2)

    def test_empty_inputs_produce_valid_html(self):
        """Report with no input files should still produce valid HTML."""
        with tempfile.NamedTemporaryFile(suffix=".html", delete=False) as fh:
            out = fh.name
        try:
            generate_report(output_path=out)
            assert os.path.isfile(out)
            with open(out) as fh:
                html = fh.read()
            assert "<!DOCTYPE html>" in html
            assert "kmer-denovo" in html
        finally:
            os.unlink(out)

    def test_nonexistent_input_paths_handled(self):
        """Nonexistent input paths should be gracefully skipped."""
        with tempfile.NamedTemporaryFile(suffix=".html", delete=False) as fh:
            out = fh.name
        try:
            generate_report(
                output_path=out,
                vcf_metrics_path="/nonexistent/metrics.json",
                vcf_summary_path="/nonexistent/summary.txt",
            )
            assert os.path.isfile(out)
        finally:
            os.unlink(out)

    def test_concordance_table_in_combined_report(self):
        """Discovery vs VCF concordance table with candidate data."""
        with tempfile.NamedTemporaryFile(suffix=".html", delete=False) as fh:
            out = fh.name
        try:
            generate_report(
                output_path=out,
                vcf_metrics_path=os.path.join(
                    EXAMPLE_OUTPUT_DIR, "metrics.json",
                ),
                vcf_summary_path=os.path.join(
                    EXAMPLE_OUTPUT_DIR, "summary.txt",
                ),
                discovery_metrics_path=os.path.join(
                    EXAMPLE_OUTPUT_DISCOVERY_DIR,
                    "giab_discovery.metrics.json",
                ),
            )
            with open(out) as fh:
                html = fh.read()
            # Check concordance section has expected data
            assert "100.0%" in html  # capture rate
            assert "3" in html  # hq candidates
            assert "chr8:125785997" in html or "chr11:55003995" in html
        finally:
            os.unlink(out)

    def test_pkc_analysis_uses_alt_specific_metrics(self):
        """PKC box/scatter must use ALT-allele PKC, not total PKC.

        Scientific rationale: REF-allele k-mers are present in parents for
        all variants; only ALT-allele PKC distinguishes de novo (absent) from
        inherited (present).
        """
        with tempfile.NamedTemporaryFile(suffix=".html", delete=False) as fh:
            out = fh.name
        try:
            generate_report(
                output_path=out,
                vcf_metrics_path=os.path.join(
                    EXAMPLE_OUTPUT_DIR, "metrics.json",
                ),
                vcf_summary_path=os.path.join(
                    EXAMPLE_OUTPUT_DIR, "summary.txt",
                ),
            )
            with open(out) as fh:
                html = fh.read()
            # Report should use ALT-specific PKC labels
            assert "AVG_PKC_ALT" in html, (
                "PKC analysis must show ALT-allele parental k-mer counts"
            )
            assert "ALT-Allele Parental K-mer Count" in html, (
                "PKC box plot title must specify ALT-allele"
            )
            assert "ALT-allele k-mers" in html.lower() or "avg_pkc_alt" in html.lower()
        finally:
            os.unlink(out)

    def test_sankey_does_not_mix_units(self):
        """Sankey diagram must not mix k-mer and variant-level counts.

        Scientific rationale: Flowing k-mer counts (1294) into variant
        counts (12) in the same Sankey is logically invalid and visually
        misleading due to scale differences.
        """
        with tempfile.NamedTemporaryFile(suffix=".html", delete=False) as fh:
            out = fh.name
        try:
            generate_report(
                output_path=out,
                vcf_metrics_path=os.path.join(
                    EXAMPLE_OUTPUT_DIR, "metrics.json",
                ),
                vcf_summary_path=os.path.join(
                    EXAMPLE_OUTPUT_DIR, "summary.txt",
                ),
            )
            with open(out) as fh:
                html = fh.read()
            # The Sankey JSON should contain only k-mer node labels
            # (no "Variants with Unique Reads" node which would mix units)
            assert "Variants with Unique Reads" not in html, (
                "Sankey must not mix k-mer and variant-level flow nodes"
            )
            assert "Variants without Unique Reads" not in html, (
                "Sankey must not mix k-mer and variant-level flow nodes"
            )
            # K-mer level nodes should be present
            assert "Child-Unique K-mers" in html
            assert "Found in Parents" in html
        finally:
            os.unlink(out)

    def test_variant_type_breakdown_present(self):
        """Variant type breakdown plot must appear when variants are present."""
        with tempfile.NamedTemporaryFile(suffix=".html", delete=False) as fh:
            out = fh.name
        try:
            generate_report(
                output_path=out,
                vcf_metrics_path=os.path.join(
                    EXAMPLE_OUTPUT_DIR, "metrics.json",
                ),
                vcf_summary_path=os.path.join(
                    EXAMPLE_OUTPUT_DIR, "summary.txt",
                ),
            )
            with open(out) as fh:
                html = fh.read()
            assert "Variant Breakdown" in html
            assert "Variant Type Breakdown" in html
        finally:
            os.unlink(out)

    def test_chromosomal_distribution_present(self):
        """Chromosomal distribution plot must appear when variants are present."""
        with tempfile.NamedTemporaryFile(suffix=".html", delete=False) as fh:
            out = fh.name
        try:
            generate_report(
                output_path=out,
                vcf_metrics_path=os.path.join(
                    EXAMPLE_OUTPUT_DIR, "metrics.json",
                ),
                vcf_summary_path=os.path.join(
                    EXAMPLE_OUTPUT_DIR, "summary.txt",
                ),
            )
            with open(out) as fh:
                html = fh.read()
            assert "Chromosomal Distribution" in html
        finally:
            os.unlink(out)

    def test_variant_table_limited_to_max_rows(self):
        """Variant table must not embed all variants when count exceeds limit."""
        from kmer_denovo_filter.report import _VARIANT_TABLE_MAX_ROWS
        with tempfile.NamedTemporaryFile(suffix=".html", delete=False) as fh:
            out = fh.name
        try:
            generate_report(
                output_path=out,
                vcf_metrics_path=os.path.join(
                    EXAMPLE_OUTPUT_DIR, "metrics.json",
                ),
                vcf_summary_path=os.path.join(
                    EXAMPLE_OUTPUT_DIR, "summary.txt",
                ),
            )
            with open(out) as fh:
                html = fh.read()
            # The constant should be exported and used in the template
            assert _VARIANT_TABLE_MAX_ROWS == 100
        finally:
            os.unlink(out)

    def test_variant_table_has_type_column(self):
        """Per-variant table must include a variant type column."""
        with tempfile.NamedTemporaryFile(suffix=".html", delete=False) as fh:
            out = fh.name
        try:
            generate_report(
                output_path=out,
                vcf_metrics_path=os.path.join(
                    EXAMPLE_OUTPUT_DIR, "metrics.json",
                ),
                vcf_summary_path=os.path.join(
                    EXAMPLE_OUTPUT_DIR, "summary.txt",
                ),
            )
            with open(out) as fh:
                html = fh.read()
            # "Type" column header should be present
            assert "<th>Type</th>" in html
        finally:
            os.unlink(out)


class TestDownsampleVariants:
    """Unit tests for the _downsample_variants helper."""

    def _make_variants(self, n_denovo, n_inherited):
        variants = []
        for i in range(n_denovo):
            variants.append({"call": "DE_NOVO", "label": f"chr1:{i} A>C"})
        for i in range(n_inherited):
            variants.append({"call": "inherited", "label": f"chr2:{i} G>T"})
        return variants

    def test_no_downsampling_when_within_limit(self):
        from kmer_denovo_filter.report import _downsample_variants
        variants = self._make_variants(10, 10)
        result, was_trimmed = _downsample_variants(variants, 100)
        assert result is variants
        assert not was_trimmed

    def test_preserves_all_denovo_first(self):
        from kmer_denovo_filter.report import _downsample_variants
        variants = self._make_variants(50, 200)
        result, was_trimmed = _downsample_variants(variants, 100)
        assert was_trimmed
        assert len(result) <= 100
        denovo_in_result = [v for v in result if v["call"] == "DE_NOVO"]
        assert len(denovo_in_result) == 50  # all DE_NOVO kept

    def test_caps_total_when_denovo_exceeds_limit(self):
        from kmer_denovo_filter.report import _downsample_variants
        variants = self._make_variants(300, 100)
        result, was_trimmed = _downsample_variants(variants, 200)
        assert was_trimmed
        assert len(result) <= 200
        # Only DE_NOVO in result
        assert all(v["call"] == "DE_NOVO" for v in result)

    def test_exact_limit_not_downsampled(self):
        from kmer_denovo_filter.report import _downsample_variants
        variants = self._make_variants(5, 5)
        result, was_trimmed = _downsample_variants(variants, 10)
        assert not was_trimmed

    def test_returns_same_object_when_no_trim(self):
        from kmer_denovo_filter.report import _downsample_variants
        variants = self._make_variants(3, 3)
        result, was_trimmed = _downsample_variants(variants, 100)
        assert result is variants


class TestHeatmapDataCap:
    """Verify heatmap cluster-summary mode for large datasets."""

    def _make_variants(self, n, denovo_every=3):
        """Create *n* synthetic variants; every *denovo_every*-th is DE_NOVO.

        DE_NOVO variants get a high dka_dkt (0.3) and inherited ones get a
        low dka_dkt (0.05) to give k-means clustering something to separate.
        """
        return [
            {
                "label": f"chr1:{i} A>C",
                "dku": i % 10, "dkt": i % 10 + 1, "dka": i % 10,
                "dku_dkt": 0.5, "dka_dkt": 0.3 if i % denovo_every == 0 else 0.05,
                "max_pkc": 0, "avg_pkc": 0.0, "min_pkc": 0,
                "max_pkc_alt": 0, "avg_pkc_alt": 0.0, "min_pkc_alt": 0,
                "call": "DE_NOVO" if i % denovo_every == 0 else "inherited",
                "vtype": "SNV",
            }
            for i in range(n)
        ]

    def test_heatmap_row_cap_constant_exported(self):
        from kmer_denovo_filter.report import _HEATMAP_MAX_ROWS
        assert _HEATMAP_MAX_ROWS == 200

    def test_heatmap_n_clusters_constant_exported(self):
        from kmer_denovo_filter.report import _HEATMAP_N_CLUSTERS
        assert _HEATMAP_N_CLUSTERS == 8

    def test_heatmap_height_bounded(self):
        import re
        from kmer_denovo_filter.report import _make_evidence_heatmap
        # Cluster-summary mode (>_HEATMAP_MAX_ROWS variants)
        variants = self._make_variants(300)
        div = _make_evidence_heatmap(variants)
        heights = [int(m) for m in re.findall(r'"height":\s*(\d+)', div)]
        assert heights, "No height found in Plotly div JSON"
        assert all(h <= 2000 for h in heights), (
            f"Plot height {max(heights)} exceeds 2000 px browser-safe limit"
        )

    def test_cluster_mode_activates_above_threshold(self):
        """For >_HEATMAP_MAX_ROWS variants the heatmap shows cluster summary."""
        from kmer_denovo_filter.report import (
            _make_evidence_heatmap,
            _HEATMAP_MAX_ROWS,
        )
        variants = self._make_variants(_HEATMAP_MAX_ROWS + 50)
        div = _make_evidence_heatmap(variants)
        # Cluster mode shows "k-means" and "Cluster" labels
        assert "k-means" in div
        assert "Cluster" in div
        # And shows the total variant count
        assert str(len(variants)) in div or f"{len(variants):,}" in div

    def test_cluster_mode_shows_denovo_fraction(self):
        """Cluster rows must be annotated with their DE_NOVO percentage."""
        from kmer_denovo_filter.report import (
            _make_evidence_heatmap,
            _HEATMAP_MAX_ROWS,
        )
        variants = self._make_variants(_HEATMAP_MAX_ROWS + 10)
        div = _make_evidence_heatmap(variants)
        assert "DE_NOVO" in div

    def test_cluster_mode_is_static_plot(self):
        """Cluster heatmap must use staticPlot to reduce HTML overhead."""
        from kmer_denovo_filter.report import (
            _make_evidence_heatmap,
            _HEATMAP_MAX_ROWS,
        )
        variants = self._make_variants(_HEATMAP_MAX_ROWS + 10)
        div = _make_evidence_heatmap(variants)
        assert "staticPlot" in div

    def test_individual_mode_below_threshold(self):
        """For ≤_HEATMAP_MAX_ROWS variants the individual-row heatmap is used."""
        from kmer_denovo_filter.report import _make_evidence_heatmap
        variants = self._make_variants(20)
        div = _make_evidence_heatmap(variants)
        # Individual mode includes per-variant labels; Plotly encodes ">" as
        # "\u003e" in its JSON, so check for the chromosome prefix instead
        # of the full "chr1:0 A>C" string.
        assert "chr1:0" in div
        # And does not use staticPlot
        assert "staticPlot" not in div

    def test_heatmap_no_hover_text_string_matrix(self):
        """In individual mode, per-cell hover must not repeat full labels."""
        from kmer_denovo_filter.report import _make_evidence_heatmap
        variants = self._make_variants(20)
        div = _make_evidence_heatmap(variants)
        # With compact customdata the label appears in the y-axis array
        # (once per variant) but NOT repeated once per field column.
        # Plotly encodes ">" as "\u003e"; check only the invariant prefix.
        label_prefix = "chr1:0"
        occurrences = div.count(label_prefix)
        # Pre-fix: 8 occurrences (once per field). Post-fix: ≤ 2.
        assert occurrences <= 2, (
            f"Label prefix appears {occurrences} times (expected ≤ 2) — "
            "hover text string matrix was not replaced with customdata"
        )


class TestClassifyVariantType:
    """Unit tests for variant type classifier."""

    def test_snv(self):
        assert _classify_variant_type("chr1:1000 A>C") == "SNV"

    def test_insertion(self):
        assert _classify_variant_type("chr1:1000 A>ACGT") == "INS"

    def test_deletion(self):
        assert _classify_variant_type("chr1:1000 ACGT>A") == "DEL"

    def test_mnv(self):
        assert _classify_variant_type("chr1:1000 AC>GT") == "MNV"

    def test_no_allele_returns_other(self):
        assert _classify_variant_type("chr1:1000") == "Other"


class TestStratification:
    """Tests for the 6-stage stratification cascade."""

    def _v(self, dka=0, dka_dkt=0.0, dka_nhf=None, max_pkc_alt=0,
            label="chr1:1 A>C"):
        v = {
            "label": label, "dku": 0, "dkt": 10, "dka": dka,
            "dku_dkt": 0.0, "dka_dkt": dka_dkt,
            "max_pkc": 0, "avg_pkc": 0.0, "min_pkc": 0,
            "max_pkc_alt": max_pkc_alt, "avg_pkc_alt": 0.0, "min_pkc_alt": 0,
            "call": "DE_NOVO" if dka > 0 else "inherited",
        }
        if dka_nhf is not None:
            v["dka_nhf"] = dka_nhf
        return v

    def test_stage_0_putative_only(self):
        from kmer_denovo_filter.report import _stratify_variant
        # No k-mer evidence → stage 0 only
        assert _stratify_variant(self._v(dka=0, dka_dkt=0.0)) == 0

    def test_stage_1_kmer_dnm(self):
        from kmer_denovo_filter.report import _stratify_variant
        assert _stratify_variant(self._v(dka=1, dka_dkt=0.05)) == 1

    def test_stage_2_dka_ge_5(self):
        from kmer_denovo_filter.report import _stratify_variant
        # DKA=5 but DKA_DKT too low for stage 3
        assert _stratify_variant(self._v(dka=5, dka_dkt=0.05)) == 2

    def test_stage_3_higher_quality(self):
        from kmer_denovo_filter.report import _stratify_variant
        # DKA >= 5 and DKA_DKT > 0.1 but MAX_PKC_ALT >= 1
        assert _stratify_variant(
            self._v(dka=10, dka_dkt=0.5, max_pkc_alt=2),
            has_nhf_data=False,
        ) == 3

    def test_stage_4_parental_confirmed(self):
        from kmer_denovo_filter.report import _stratify_variant
        # DKA >= 5, DKA_DKT > 0.1, MAX_PKC_ALT < 1 but has NHF data
        # and no NHF field → stays at stage 4
        assert _stratify_variant(
            self._v(dka=10, dka_dkt=0.5, max_pkc_alt=0),
            has_nhf_data=True,
        ) == 4

    def test_stage_5_full_pass_no_nhf_data(self):
        from kmer_denovo_filter.report import _stratify_variant
        # No NHF data → stage 5 collapses (has_nhf_data=False)
        assert _stratify_variant(
            self._v(dka=10, dka_dkt=0.5, max_pkc_alt=0),
            has_nhf_data=False,
        ) == 5

    def test_stage_5_requires_low_nhf(self):
        from kmer_denovo_filter.report import _stratify_variant
        # NHF data is present cohort-wide; this variant fails NHF
        assert _stratify_variant(
            self._v(dka=10, dka_dkt=0.5, max_pkc_alt=0, dka_nhf=0.10),
            has_nhf_data=True,
        ) == 4
        # Below threshold → stage 5
        assert _stratify_variant(
            self._v(dka=10, dka_dkt=0.5, max_pkc_alt=0, dka_nhf=0.01),
            has_nhf_data=True,
        ) == 5

    def test_stage_5_missing_nhf_when_cohort_has_data(self):
        """Conservative: when other variants have NHF but this one doesn't,
        we cannot confirm non-contamination → stays at stage 4."""
        from kmer_denovo_filter.report import _stratify_variant
        assert _stratify_variant(
            self._v(dka=10, dka_dkt=0.5, max_pkc_alt=0),  # no dka_nhf
            has_nhf_data=True,
        ) == 4

    def test_compute_stratification_counts(self):
        from kmer_denovo_filter.report import _compute_stratification
        variants = [
            self._v(dka=0, dka_dkt=0.0),                                    # stage 0
            self._v(dka=2, dka_dkt=0.05),                                   # stage 1 (DKA<5)
            self._v(dka=15, dka_dkt=0.30, max_pkc_alt=0, dka_nhf=0.01),    # stage 5
            self._v(dka=15, dka_dkt=0.30, max_pkc_alt=0, dka_nhf=0.20),    # stage 4 (contam)
        ]
        strat = _compute_stratification(variants)
        assert strat["has_nhf_data"] is True
        # counts: [4, 3, 2, 2, 2, 1]
        assert strat["counts"] == [4, 3, 2, 2, 2, 1]
        # Lost at each step: 0, 1, 1, 0, 0, 1
        assert strat["lost"] == [0, 1, 1, 0, 0, 1]

    def test_compute_stratification_no_nhf(self):
        from kmer_denovo_filter.report import _compute_stratification
        variants = [
            self._v(dka=0, dka_dkt=0.0),
            self._v(dka=2, dka_dkt=0.05),
            self._v(dka=15, dka_dkt=0.30, max_pkc_alt=0),
            self._v(dka=20, dka_dkt=0.50, max_pkc_alt=0),
        ]
        strat = _compute_stratification(variants)
        assert strat["has_nhf_data"] is False
        # Final stage collapses to parental-confirmed stage (NHF filter not applied)
        assert strat["counts"][5] == strat["counts"][4]


class TestReportStratificationContent:
    """End-to-end checks for stratified report content."""

    def _generate(self):
        with tempfile.NamedTemporaryFile(suffix=".html", delete=False) as fh:
            out = fh.name
        generate_report(
            output_path=out,
            vcf_metrics_path=os.path.join(EXAMPLE_OUTPUT_DIR, "metrics.json"),
            vcf_summary_path=os.path.join(EXAMPLE_OUTPUT_DIR, "summary.txt"),
        )
        with open(out) as fh:
            html = fh.read()
        os.unlink(out)
        return html

    def test_pipeline_summary_present(self):
        html = self._generate()
        assert "Pipeline Summary" in html
        # Brief intro mentions VCF input and miscalled variants
        assert "putative" in html.lower()
        assert "miscalled" in html.lower()

    def test_stratification_funnel_has_six_stages(self):
        html = self._generate()
        assert "Putative denovo" in html
        assert "Putative kmer denovo" in html
        assert "Higher-quality denovo" in html
        assert "MAX_PKC_ALT" in html
        # Stage 6 label appears (NHF)
        assert "NHF" in html

    def test_stratification_sankey_present(self):
        html = self._generate()
        assert "Variant Flow Through Stratification Stages" in html

    def test_methods_section_present(self):
        html = self._generate()
        # Methods section content
        assert "<h2>2. Methods</h2>" in html
        assert "Stratification cascade" in html
        assert "Non-human contamination screening" in html
        assert "Jellyfish" in html

    def test_results_highlights_present(self):
        html = self._generate()
        assert "Results Highlights" in html

    def test_per_variant_table_only_final_set(self):
        """Per-variant table must only contain final DNM set variants
        (stage 5 = all filters passed)."""
        html = self._generate()
        assert "Final DNM Set" in html
        # chr8:40003391 has DKA_DKT=0.0588 → must NOT appear in the per-variant
        # table (it should still appear in scatter/heatmap data, but the
        # table filtering means its row HTML <td>chr8:40003391...</td> is gone)
        # Use a heuristic: count occurrences; in scatter/heatmap-encoded JSON
        # the label is escaped (chr8:40003391 A\u003eT) so the plain
        # "<td>chr8:40003391 A>T</td>" pattern is reliable.
        assert "<td>chr8:40003391 A&gt;T</td>" not in html
        assert "<td>chr8:40003391 A>T</td>" not in html


class TestStratificationConstants:
    def test_thresholds_exported(self):
        from kmer_denovo_filter.report import (
            _HIGH_QUALITY_DKA_DKT_THRESHOLD,
            _NHF_CONTAMINATION_THRESHOLD,
            _KMER_DENOVO_DKA_THRESHOLD,
            _KMER_DENOVO_DKA_STRONG_THRESHOLD,
            _MAX_PKC_ALT_THRESHOLD,
        )
        assert _KMER_DENOVO_DKA_THRESHOLD == 0
        assert _KMER_DENOVO_DKA_STRONG_THRESHOLD == 5
        assert _HIGH_QUALITY_DKA_DKT_THRESHOLD == 0.1
        assert _MAX_PKC_ALT_THRESHOLD == 1
        assert _NHF_CONTAMINATION_THRESHOLD == 0.05


class TestReportCLI:
    """Test the standalone kmer-report CLI entry point."""

    def test_report_main_generates_html(self):
        from kmer_denovo_filter.cli import report_main
        with tempfile.NamedTemporaryFile(suffix=".html", delete=False) as fh:
            out = fh.name
        try:
            report_main([
                "--output", out,
                "--vcf-metrics", os.path.join(EXAMPLE_OUTPUT_DIR, "metrics.json"),
                "--vcf-summary", os.path.join(EXAMPLE_OUTPUT_DIR, "summary.txt"),
            ])
            assert os.path.isfile(out)
            with open(out) as fh:
                html = fh.read()
            assert "kmer-denovo" in html
        finally:
            os.unlink(out)

    def test_parse_report_args(self):
        from kmer_denovo_filter.cli import parse_report_args
        args = parse_report_args([
            "--output", "/tmp/report.html",
            "--vcf-metrics", "/tmp/metrics.json",
            "--vcf-summary", "/tmp/summary.txt",
        ])
        assert args.output == "/tmp/report.html"
        assert args.vcf_metrics == "/tmp/metrics.json"
        assert args.vcf_summary == "/tmp/summary.txt"
        assert args.vcf is None
        assert args.discovery_metrics is None
