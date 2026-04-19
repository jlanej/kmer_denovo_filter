"""Interactive HTML report generation for kmer-denovo and kmer-discovery.

Generates a self-contained Plotly + Jinja2 HTML report from pipeline
output files (metrics.json, summary.txt, annotated VCF / discovery BED).
The report is designed as a publication-quality landing page that
communicates the scientific rationale and filtering results of the
k-mer de novo filtering strategy.
"""

import json
import logging
import os
import re
import statistics as stats

logger = logging.getLogger(__name__)

# Maximum rows shown in the per-variant detail table.  Only "higher quality"
# de novo variants (stage 3+: DKA_DKT > _HIGH_QUALITY_DKA_DKT_THRESHOLD) are
# included in the table; inherited and low-evidence calls are excluded because
# the unfiltered list (often 30 000+ rows) is not informative for a reviewer.
_VARIANT_TABLE_MAX_ROWS = 100

# ---------------------------------------------------------------------------
# Stratification thresholds
# ---------------------------------------------------------------------------
# Six progressively stricter filtering stages are applied to every input
# variant.  These thresholds are referenced by the funnel chart, the Sankey
# diagram, the per-variant table, and the contamination plots so the report
# tells a single, coherent story.
#
#   Stage 0 — Putative denovo (input)            : every variant in the VCF
#   Stage 1 — Putative kmer denovo               : DKA > 0
#   Stage 2 — Putative kmer denovo (stronger)    : DKA >= 5
#   Stage 3 — Higher-quality denovo              : DKA_DKT > 0.1
#   Stage 4 — Higher-quality denovo (parental)   : MAX_PKC_ALT < 1
#   Stage 5 — HQ, non-contamination             : DKA_NHF < 0.05
#
_KMER_DENOVO_DKA_THRESHOLD = 0          # Stage 1: DKA > 0
_KMER_DENOVO_DKA_STRONG_THRESHOLD = 5   # Stage 2: DKA >= 5
_HIGH_QUALITY_DKA_DKT_THRESHOLD = 0.1   # Stage 3: DKA_DKT > 0.1
_MAX_PKC_ALT_THRESHOLD = 1              # Stage 4: MAX_PKC_ALT < 1
_NHF_CONTAMINATION_THRESHOLD = 0.05     # Stage 5: DKA_NHF < 0.05

# Human-readable labels for the six stages — used by funnel, Sankey, and
# legends so all visualisations are consistent.
_STAGE_LABELS = [
    "Putative denovo<br>(input VCF)",
    "Putative kmer denovo<br>(DKA &gt; 0)",
    "Putative kmer denovo<br>(DKA &ge; 5)",
    "Higher-quality denovo<br>(DKA_DKT &gt; 0.1)",
    "Higher-quality denovo<br>(MAX_PKC_ALT &lt; 1)",
    "HQ, not contamination<br>(NHF &lt; 0.05)",
]
_STAGE_SHORT_LABELS = [
    "Putative",
    "Kmer DNM",
    "Kmer DNM (DKA\u22655)",
    "Higher quality",
    "HQ (PKC\u22640)",
    "HQ + non-contam",
]
# Visually distinct colour palette: each stage moves toward green to convey
# increasing confidence, while still being colour-blind safe.
_STAGE_COLORS = [
    "#4C78A8", "#F58518", "#E45756", "#72B7B2", "#EECA3B", "#54A24B",
]

# Maximum rows in the evidence heatmap.  Above this threshold the heatmap
# switches to cluster-summary mode: k-means is run on all variants and one row
# per cluster is shown (centroid values).  This scales to 100k+ variants.
_HEATMAP_MAX_ROWS = 200

# Number of k-means clusters used in cluster-summary heatmap mode.
_HEATMAP_N_CLUSTERS = 8

# Maximum data points in continuous scatter plots.  Beyond this the serialised
# Plotly JSON becomes multi-megabyte and slows or prevents rendering.
# DE_NOVO variants are always kept; inherited variants are uniformly sampled.
_SCATTER_MAX_POINTS = 2000


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _downsample_variants(variants, max_points):
    """Return at most *max_points* variants, keeping all DE_NOVO calls.

    When the full list fits within *max_points* the original list is returned
    unchanged (no copy).  Otherwise all DE_NOVO variants are kept and the
    inherited variants are uniformly sub-sampled to fill the remaining quota.

    Returns ``(sampled_variants, was_downsampled)`` so callers can annotate
    plots when the dataset has been reduced.
    """
    if len(variants) <= max_points:
        return variants, False

    denovo = [v for v in variants if v["call"] == "DE_NOVO"]
    inherited = [v for v in variants if v["call"] != "DE_NOVO"]

    if len(denovo) >= max_points:
        return denovo[:max_points], True

    remaining = max_points - len(denovo)
    step = max(1, len(inherited) // remaining)
    sampled_inherited = inherited[::step][:remaining]
    return denovo + sampled_inherited, True


def _kmeans_cluster(z_matrix, n_clusters, max_iter=100):
    """Lightweight k-means clustering using numpy (required by plotly).

    Uses k-means++ seeding for reproducible, well-separated initial centres.
    Returns a list of integer cluster labels with length ``len(z_matrix)``.
    When numpy is unavailable falls back to evenly-spaced bucket assignment
    so the heatmap always renders.
    """
    try:
        import numpy as np
    except ImportError:
        # Fallback: assign rows to buckets by their position (no real clustering)
        n = len(z_matrix)
        step = max(1, n // n_clusters)
        return [min(i // step, n_clusters - 1) for i in range(n)]

    X = np.array(z_matrix, dtype=np.float64)
    n_samples = len(X)

    if n_samples <= n_clusters:
        return list(range(n_samples))

    # K-means++ seeding (fixed seed=42 for reproducible report output;
    # identical input always produces identical cluster assignments so the
    # report is deterministic across regenerations)
    rng = np.random.RandomState(42)
    center_indices = [int(rng.randint(n_samples))]
    for _ in range(n_clusters - 1):
        # Distance from each point to its nearest existing centre
        sq_dists = np.min(
            [np.sum((X - X[ci]) ** 2, axis=1) for ci in center_indices],
            axis=0,
        )
        sq_dists = np.maximum(sq_dists, 0.0)
        total = sq_dists.sum()
        if total == 0:
            center_indices.append(int(rng.randint(n_samples)))
        else:
            center_indices.append(int(rng.choice(n_samples, p=sq_dists / total)))

    centers = X[center_indices].copy()
    labels = np.zeros(n_samples, dtype=np.int32)

    for _ in range(max_iter):
        # Vectorised squared-distance matrix: shape (n_samples, n_clusters).
        # Broadcast rather than np.einsum to keep the peak working memory to
        # one (n, k, f) tensor (~50 MB for 100k variants, k=8, f=8 float64).
        sq_dists = np.sum((X[:, None, :] - centers[None, :, :]) ** 2, axis=2)
        new_labels = np.argmin(sq_dists, axis=1).astype(np.int32)

        if np.array_equal(new_labels, labels):
            break
        labels = new_labels

        # Update cluster centres
        for k in range(n_clusters):
            mask = labels == k
            if mask.any():
                centers[k] = X[mask].mean(axis=0)

    return labels.tolist()


# ---------------------------------------------------------------------------
# Data loaders
# ---------------------------------------------------------------------------


def _load_metrics(metrics_path):
    """Load and return a metrics dict from a JSON file."""
    with open(metrics_path) as fh:
        return json.load(fh)


def _load_summary_variants(summary_path):
    """Parse per-variant rows from the summary.txt file.

    Returns a list of dicts with keys matching the summary table columns.
    Works for both VCF-mode and discovery-mode summary files.
    """
    variants = []
    in_table = False
    header_seen = False

    with open(summary_path) as fh:
        for line in fh:
            stripped = line.strip()
            if stripped.startswith("Per-Variant Results"):
                in_table = True
                continue
            if in_table and not header_seen:
                # Skip the column header line and the separator line
                if stripped.startswith("Variant") or stripped.startswith("---"):
                    continue
                if stripped.startswith("=") or not stripped:
                    continue
                # Once we see a data line (starts with chr or position),
                # we know the header has been consumed
                header_seen = True
                # Fall through to parse this line
            if in_table and header_seen:
                if not stripped or stripped.startswith("="):
                    break
                parts = stripped.split()
                if len(parts) < 13:
                    continue
                # Validate this is a data row by checking last column
                call = parts[-1]
                if call not in ("DE_NOVO", "inherited"):
                    continue
                # Variant label may contain spaces (e.g. allele description)
                # Format: label  DKU DKT DKA DKU_DKT DKA_DKT MAX_PKC AVG_PKC
                # MIN_PKC MAX_PKC_ALT AVG_PKC_ALT MIN_PKC_ALT Call
                # Parse from the right since variant label can have spaces
                call = parts[-1]
                min_pkc_alt = int(parts[-2])
                avg_pkc_alt = float(parts[-3])
                max_pkc_alt = int(parts[-4])
                min_pkc = int(parts[-5])
                avg_pkc = float(parts[-6])
                max_pkc = int(parts[-7])
                dka_dkt = float(parts[-8])
                dku_dkt = float(parts[-9])
                dka = int(parts[-10])
                dkt = int(parts[-11])
                dku = int(parts[-12])
                label = " ".join(parts[:-12])

                variants.append({
                    "label": label,
                    "dku": dku,
                    "dkt": dkt,
                    "dka": dka,
                    "dku_dkt": dku_dkt,
                    "dka_dkt": dka_dkt,
                    "max_pkc": max_pkc,
                    "avg_pkc": avg_pkc,
                    "min_pkc": min_pkc,
                    "max_pkc_alt": max_pkc_alt,
                    "avg_pkc_alt": avg_pkc_alt,
                    "min_pkc_alt": min_pkc_alt,
                    "call": call,
                })
    return variants


def _load_summary_counts(summary_path):
    """Parse the header counts section of summary.txt.

    Returns a dict with keys like 'total_candidates', 'likely_denovo',
    'inherited'.
    """
    counts = {}
    with open(summary_path) as fh:
        for line in fh:
            stripped = line.strip()
            if "Total candidates analyzed:" in stripped:
                counts["total_candidates"] = int(stripped.split(":")[-1].strip())
            elif "Likely de novo" in stripped:
                counts["likely_denovo"] = int(stripped.split(":")[-1].strip())
            elif "Inherited / unclear" in stripped:
                counts["inherited"] = int(stripped.split(":")[-1].strip())
    return counts


def _load_vcf_kraken2_annotations(vcf_path):
    """Load Kraken2 fraction annotations from an annotated VCF.

    Returns a list of dicts with variant label and fraction fields.
    Only returns data when Kraken2 annotations are present.
    """
    try:
        import pysam
    except ImportError:
        logger.warning("pysam not available; skipping VCF Kraken2 annotations")
        return []

    kraken2_fields = [
        "DKA_BF", "DKA_AF", "DKA_FF", "DKA_PF", "DKA_VF",
        "DKA_UCF", "DKA_NHF", "DKA_UF", "DKA_HLF",
    ]

    results = []
    try:
        with pysam.VariantFile(vcf_path) as vcf:
            # Check if any Kraken2 annotations are in the header
            header_ids = set()
            for rec in vcf.header.records:
                rid = rec.get("ID", "")
                if rid:
                    header_ids.add(rid)
            has_kraken2 = any(f in header_ids for f in kraken2_fields)
            if not has_kraken2:
                return []

            for rec in vcf:
                alt = ",".join(str(a) for a in rec.alts) if rec.alts else "."
                label = f"{rec.chrom}:{rec.pos} {rec.ref}>{alt}"
                row = {"label": label}
                for sample in rec.samples.values():
                    for field in kraken2_fields:
                        val = sample.get(field)
                        if val is not None:
                            row[field] = float(val)
                    break  # first sample only
                if len(row) > 1:
                    results.append(row)
    except Exception:
        logger.debug("Could not parse VCF for Kraken2 annotations", exc_info=True)
        return []

    return results


# ---------------------------------------------------------------------------
# Stratification: 4-stage filtering cascade
# ---------------------------------------------------------------------------


def _merge_kraken2_into_variants(variants, kraken2_data):
    """Attach Kraken2 fraction fields to each variant in place.

    Looks up each variant's label in *kraken2_data* and copies the per-variant
    fraction fields (``DKA_NHF``, ``DKA_HLF``, etc.) onto the variant dict
    using lowercase keys (``dka_nhf``, ``dka_hlf``, …) for downstream Python
    code.  Variants with no matching Kraken2 record are left unchanged so
    callers can detect "no contamination data" via ``"dka_nhf" not in v``.
    """
    if not kraken2_data:
        return variants
    kraken_map = {r["label"]: r for r in kraken2_data}
    fields = [
        "DKA_BF", "DKA_AF", "DKA_FF", "DKA_PF", "DKA_VF",
        "DKA_UCF", "DKA_NHF", "DKA_UF", "DKA_HLF",
    ]
    for v in variants:
        k = kraken_map.get(v["label"])
        if not k:
            continue
        for f in fields:
            if f in k:
                v[f.lower()] = float(k[f])
    return variants


def _stratify_variant(v, has_nhf_data=None):
    """Return the highest stratification stage (0–5) reached by *v*.

    Stages (cumulative; see module docstring for definitions):

    * 0 — Putative denovo (always reached for any variant in the input)
    * 1 — Putative kmer denovo (``dka > 0``)
    * 2 — Putative kmer denovo, stronger (``dka >= 5``)
    * 3 — Higher-quality denovo (``dka_dkt > 0.1``)
    * 4 — Higher-quality denovo, parental (``max_pkc_alt < 1``)
    * 5 — Higher-quality, non-contamination (stage 4 + ``dka_nhf < 0.05``)

    If *has_nhf_data* is False (no Kraken2 annotations available anywhere in
    the cohort) stage 5 collapses into stage 4 — i.e. any stage-4 variant
    is also counted as stage 5 — because we cannot meaningfully separate
    contamination without NHF data.  When *has_nhf_data* is True but the
    individual variant lacks a ``dka_nhf`` value the variant is treated as
    NOT passing stage 5 (conservative).
    """
    stage = 0
    if v.get("dka", 0) > _KMER_DENOVO_DKA_THRESHOLD:
        stage = 1
    if stage >= 1 and v.get("dka", 0) >= _KMER_DENOVO_DKA_STRONG_THRESHOLD:
        stage = 2
    if stage >= 2 and v.get("dka_dkt", 0.0) > _HIGH_QUALITY_DKA_DKT_THRESHOLD:
        stage = 3
    if stage >= 3 and v.get("max_pkc_alt", 999) < _MAX_PKC_ALT_THRESHOLD:
        stage = 4
    if stage >= 4:
        if not has_nhf_data:
            # No Kraken2 anywhere → cannot apply contamination filter, so
            # collapse stage 5 into stage 4 for reporting.
            stage = 5
        else:
            nhf = v.get("dka_nhf")
            if nhf is not None and nhf < _NHF_CONTAMINATION_THRESHOLD:
                stage = 5
    return stage


def _compute_stratification(variants, has_nhf_data=None):
    """Annotate variants with their stage and return summary counts.

    Adds an integer ``stage`` field (0–5) to every variant and returns a
    dict with both the cumulative counts at each stage and convenience
    sub-lists, e.g.::

        {
            "has_nhf_data": True,
            "counts": [22, 12, 8, 6, 4, 3],   # cumulative pass counts
            "labels": [...],                    # human-readable labels
            "lost":   [0,  10, 4, 2, 2, 1],    # variants dropped at each step
            "stage_groups": {
                0: [...all variants...],
                1: [...stage>=1...],
                2: [...stage>=2...],
                3: [...stage>=3...],
                4: [...stage>=4...],
                5: [...stage>=5...],
            }
        }
    """
    if has_nhf_data is None:
        has_nhf_data = any("dka_nhf" in v for v in variants)

    for v in variants:
        v["stage"] = _stratify_variant(v, has_nhf_data=has_nhf_data)

    counts = [
        len(variants),
        sum(1 for v in variants if v["stage"] >= 1),
        sum(1 for v in variants if v["stage"] >= 2),
        sum(1 for v in variants if v["stage"] >= 3),
        sum(1 for v in variants if v["stage"] >= 4),
        sum(1 for v in variants if v["stage"] >= 5),
    ]
    lost = [0] + [counts[i - 1] - counts[i] for i in range(1, 6)]
    stage_groups = {
        i: [v for v in variants if v["stage"] >= i] for i in range(6)
    }
    return {
        "has_nhf_data": has_nhf_data,
        "counts": counts,
        "labels": list(_STAGE_LABELS),
        "short_labels": list(_STAGE_SHORT_LABELS),
        "colors": list(_STAGE_COLORS),
        "lost": lost,
        "stage_groups": stage_groups,
    }


def _load_discovery_regions(metrics_path):
    """Load discovery region details from discovery metrics.json."""
    metrics = _load_metrics(metrics_path)
    return metrics.get("regions", [])


def _load_discovery_candidate_comparison(metrics_path):
    """Load candidate comparison from discovery metrics.json."""
    metrics = _load_metrics(metrics_path)
    return metrics.get("candidate_comparison", {})


def _load_discovery_dnm_evaluation(metrics_path):
    """Load DNM region evaluation from discovery metrics.json."""
    metrics = _load_metrics(metrics_path)
    return metrics.get("dnm_evaluation", {})


# ---------------------------------------------------------------------------
# Helper: self-contained inline Plotly bundle
# ---------------------------------------------------------------------------

def _get_plotly_bundle():
    """Return the minified Plotly.js bundle for inline embedding.

    Embedding Plotly inline makes the HTML report self-contained so it
    renders correctly in offline/HPC environments and always uses the JS
    version that matches the installed Python plotly package.
    """
    import plotly.offline as pyo
    return pyo.get_plotlyjs()


# ---------------------------------------------------------------------------
# Plot generators (return self-contained div HTML strings)
# ---------------------------------------------------------------------------

def _plotly_div(fig, div_id, config=None):
    """Return a self-contained ``<div>+<script>`` HTML snippet for a figure.

    Uses Plotly's own ``to_html(full_html=False, include_plotlyjs=False)``
    so that data is embedded directly as JS object literals (no JSON.parse
    step) and the correct Plotly.js API is used regardless of version.
    A fixed ``div_id`` is required for deterministic (idempotent) output.

    Pass ``config={"staticPlot": True}`` to disable all browser-side
    interactivity (hover, zoom, pan) for static summary plots.
    """
    import plotly.io as pio
    return pio.to_html(
        fig,
        full_html=False,
        include_plotlyjs=False,
        div_id=div_id,
        config=config if config is not None else {"responsive": True},
    )


    return _plotly_div(fig, div_id)


def _make_stratification_funnel(stratification, div_id="strat-funnel-plot"):
    """Bar chart showing the variant-level 6-stage filtering cascade.

    The six stages (defined in the module header) describe how many input
    variants survive at each progressively stricter filter.  This is the
    primary "results at a glance" plot for the report — a reviewer should be
    able to read off how many putative DNMs were called from the input VCF
    and how many remained after k-mer evidence and contamination filtering.
    """
    import plotly.graph_objects as go

    counts = stratification["counts"]
    labels = stratification["labels"]
    colors = stratification["colors"]
    has_nhf = stratification["has_nhf_data"]

    # If we have no contamination data, fade stage 6 so the reviewer sees
    # explicitly that the contamination filter could not be applied.
    bar_colors = list(colors)
    if not has_nhf:
        bar_colors[5] = "#cccccc"

    text = [f"{c:,}" for c in counts]
    if not has_nhf:
        # Stage 6 collapses to stage 5 when no NHF; annotate clearly.
        text[5] = f"{counts[5]:,}*"

    fig = go.Figure(data=[go.Bar(
        x=labels,
        y=counts,
        marker_color=bar_colors,
        text=text,
        textposition="outside",
        hovertemplate="<b>%{x}</b><br>Variants: %{y:,}<extra></extra>",
    )])

    # Annotate retained-percentage relative to stage 1
    if counts[0] > 0:
        for i in range(1, len(counts)):
            pct = 100.0 * counts[i] / counts[0]
            fig.add_annotation(
                x=labels[i], y=counts[i],
                text=f"{pct:.1f}% of input",
                showarrow=False, yshift=-18,
                font=dict(size=10, color="#666"),
            )

    title = "Variant Stratification Funnel"
    if not has_nhf:
        title += (
            "<br><sup>* No Kraken2 contamination data available — "
            "stage 6 collapses to stage 5 (NHF filter not applied)</sup>"
        )

    fig.update_layout(
        title=dict(text=title, font=dict(size=18)),
        yaxis_title="Variant Count",
        template="plotly_white",
        height=420,
        margin=dict(t=80, b=80),
        showlegend=False,
    )
    return _plotly_div(fig, div_id)


def _make_stratification_sankey(stratification, div_id="strat-sankey-plot"):
    """Sankey diagram showing how variants flow through the 6-stage cascade.

    Each stage emits two flows: variants that PASS to the next stage, and
    variants that are FILTERED OUT for failing the stage's criterion.  The
    diagram makes the attrition at each step explicit so the reviewer can see
    at a glance which filter is most/least selective for this dataset.
    """
    import plotly.graph_objects as go

    counts = stratification["counts"]
    short_labels = stratification["short_labels"]
    colors = stratification["colors"]
    has_nhf = stratification["has_nhf_data"]

    drop_reasons = [
        "Filtered: DKA = 0",
        "Filtered: DKA &lt; 5",
        "Filtered: DKA_DKT &le; 0.1",
        "Filtered: MAX_PKC_ALT &ge; 1",
        "Filtered: NHF &ge; 0.05 (contamination)",
    ]

    node_labels = [
        f"{short_labels[0]} ({counts[0]:,})",
        f"{short_labels[1]} ({counts[1]:,})",
        f"{short_labels[2]} ({counts[2]:,})",
        f"{short_labels[3]} ({counts[3]:,})",
        f"{short_labels[4]} ({counts[4]:,})",
        f"{short_labels[5]} ({counts[5]:,})"
        + ("" if has_nhf else " *no NHF data*"),
        # Drop nodes (indices 6–10)
        f"{drop_reasons[0]} ({counts[0] - counts[1]:,})",
        f"{drop_reasons[1]} ({counts[1] - counts[2]:,})",
        f"{drop_reasons[2]} ({counts[2] - counts[3]:,})",
        f"{drop_reasons[3]} ({counts[3] - counts[4]:,})",
        f"{drop_reasons[4]} ({counts[4] - counts[5]:,})"
        + ("" if has_nhf else " — N/A"),
    ]
    node_colors = list(colors) + ["#bbbbbb"] * 5

    # Edges: pass flows + drop flows
    sources = [0, 0,  1, 1,  2, 2,  3, 3,  4, 4]
    targets = [1, 6,  2, 7,  3, 8,  4, 9,  5, 10]
    values = [
        max(1, counts[1]),                   counts[0] - counts[1],
        max(1, counts[2]),                   counts[1] - counts[2],
        max(1, counts[3]),                   counts[2] - counts[3],
        max(1, counts[4]),                   counts[3] - counts[4],
        max(1, counts[5]),                   counts[4] - counts[5],
    ]
    # Replace any 0 with 1 so a missing flow is still drawn but tiny.
    values = [max(1, v) for v in values]

    fig = go.Figure(data=[go.Sankey(
        node=dict(
            pad=18, thickness=18,
            label=node_labels,
            color=node_colors,
        ),
        link=dict(source=sources, target=targets, value=values),
    )])
    fig.update_layout(
        title=dict(
            text="Variant Flow Through Stratification Stages",
            font=dict(size=18),
        ),
        template="plotly_white",
        height=420,
        margin=dict(t=60, b=20),
    )
    return _plotly_div(fig, div_id)


def _make_nhf_distribution_plot(variants, div_id="nhf-dist-plot"):
    """Histogram of per-variant non-human (NHF) fraction for contaminated variants.

    Returns ``None`` when no Kraken2 NHF data is available for any variant.
    Otherwise produces a histogram of NHF for variants with NHF >= 0.05
    (i.e. those that would be flagged as putative contamination) among kmer-DNM
    candidates (stage >= 1), so the reviewer sees the contamination profile
    of variants that actually trigger the contamination filter.
    """
    import plotly.graph_objects as go

    nhf_values = [v["dka_nhf"] for v in variants
                  if v.get("stage", 0) >= 1 and "dka_nhf" in v
                  and v["dka_nhf"] >= _NHF_CONTAMINATION_THRESHOLD]
    if not nhf_values:
        return None

    fig = go.Figure()
    fig.add_trace(go.Histogram(
        x=nhf_values, nbinsx=30,
        marker_color="#E45756", opacity=0.85,
        hovertemplate="NHF: %{x:.3f}<br>Count: %{y}<extra></extra>",
    ))
    fig.add_vline(
        x=_NHF_CONTAMINATION_THRESHOLD,
        line_dash="dash", line_color="#E45756", line_width=2,
        annotation_text=f"Contamination threshold ({_NHF_CONTAMINATION_THRESHOLD})",
        annotation_position="top right",
        annotation_font=dict(size=11, color="#E45756"),
    )
    fig.update_layout(
        title=dict(
            text=(f"Non-Human Fraction (NHF) Distribution — Putative Contamination "
                  f"(NHF &ge; {_NHF_CONTAMINATION_THRESHOLD}, n={len(nhf_values)})"),
            font=dict(size=16),
        ),
        xaxis_title="DKA_NHF (fraction of DKA reads classified as non-human)",
        yaxis_title="Variant Count",
        template="plotly_white",
        height=380,
        margin=dict(t=60, b=40),
    )
    return _plotly_div(fig, div_id)


def _make_kmer_funnel_chart(metrics, mode="vcf", div_id="funnel-plot"):
    """Create a waterfall/bar chart showing the k-mer filtering cascade."""
    import plotly.graph_objects as go

    if mode == "vcf":
        labels = [
            "Total Child<br>K-mers",
            "Found in<br>Parents",
            "Child-Unique<br>K-mers",
        ]
        total = metrics.get("total_child_kmers", 0)
        parent_found = metrics.get("parent_found_kmers", 0)
        unique = metrics.get("child_unique_kmers", 0)
        values = [total, parent_found, unique]
        colors = ["#4C78A8", "#E45756", "#54A24B"]
    else:
        labels = [
            "Child Candidate<br>K-mers",
            "Non-Reference<br>K-mers",
            "Proband-Unique<br>K-mers",
        ]
        values = [
            metrics.get("child_candidate_kmers", 0),
            metrics.get("non_ref_kmers", 0),
            metrics.get("proband_unique_kmers", 0),
        ]
        colors = ["#4C78A8", "#F58518", "#54A24B"]

    fig = go.Figure(data=[go.Bar(
        x=labels, y=values,
        marker_color=colors,
        text=[f"{v:,}" for v in values],
        textposition="outside",
        hovertemplate="%{x}<br>Count: %{y:,}<extra></extra>",
    )])
    fig.update_layout(
        title=dict(
            text="K-mer Filtering Funnel",
            font=dict(size=18),
        ),
        yaxis_title="K-mer Count",
        template="plotly_white",
        height=400,
        margin=dict(t=60, b=40),
    )

    # Add percentage annotations
    if values[0] > 0:
        for i in range(1, len(values)):
            pct = 100 * values[i] / values[0]
            fig.add_annotation(
                x=labels[i], y=values[i],
                text=f"{pct:.1f}%",
                showarrow=False, yshift=-15,
                font=dict(size=11, color="#666"),
            )

    return _plotly_div(fig, div_id)


def _make_sankey_diagram(metrics, mode="vcf", div_id="sankey-plot"):
    """Create a Sankey diagram showing the filtering flow."""
    import plotly.graph_objects as go

    if mode == "vcf":
        total = metrics.get("total_child_kmers", 0)
        parent_found = metrics.get("parent_found_kmers", 0)
        unique = metrics.get("child_unique_kmers", 0)

        node_labels = [
            f"Total Child K-mers<br>({total:,})",
            f"Found in Parents<br>({parent_found:,})",
            f"Child-Unique K-mers<br>({unique:,})",
        ]
        node_colors = ["#4C78A8", "#E45756", "#54A24B"]

        source = [0, 0]
        target = [1, 2]
        value = [max(1, parent_found), max(1, unique)]
    else:
        child_cand = metrics.get("child_candidate_kmers", 0)
        non_ref = metrics.get("non_ref_kmers", 0)
        proband_unique = metrics.get("proband_unique_kmers", 0)
        ref_kmers = child_cand - non_ref
        parent_kmers = non_ref - proband_unique

        node_labels = [
            f"Child Candidate K-mers ({child_cand:,})",
            f"Reference K-mers ({ref_kmers:,})",
            f"Non-Reference ({non_ref:,})",
            f"Parental K-mers ({parent_kmers:,})",
            f"Proband-Unique ({proband_unique:,})",
        ]
        node_colors = ["#4C78A8", "#BAB0AC", "#F58518", "#E45756", "#54A24B"]

        source = [0, 0, 2, 2]
        target = [1, 2, 3, 4]
        value = [max(1, ref_kmers), max(1, non_ref),
                 max(1, parent_kmers), max(1, proband_unique)]

    fig = go.Figure(data=[go.Sankey(
        node=dict(
            pad=20, thickness=20,
            label=node_labels,
            color=node_colors,
        ),
        link=dict(source=source, target=target, value=value),
    )])
    fig.update_layout(
        title=dict(
            text="K-mer Filtering Flow",
            font=dict(size=18),
        ),
        template="plotly_white",
        height=350,
        margin=dict(t=60, b=20),
    )
    return _plotly_div(fig, div_id)


def _make_dka_dkt_histogram(variants, div_id="histogram-plot"):
    """Create a histogram of DKA_DKT ratios with threshold marker.

    Variants are split into stratification groups for visual context:
    inherited / DKA = 0 (greyed) vs. kmer-DNM candidates (DKA > 0, blue).
    The dashed line marks the higher-quality threshold (DKA_DKT > 0.1).
    """
    import plotly.graph_objects as go

    low = [v["dka_dkt"] for v in variants if v.get("dka", 0) == 0]
    kmer_dnm = [v["dka_dkt"] for v in variants if v.get("dka", 0) > 0]

    fig = go.Figure()
    if low:
        fig.add_trace(go.Histogram(
            x=low, nbinsx=30,
            marker_color="#BAB0AC", opacity=0.75,
            name="DKA = 0 (no kmer evidence)",
            hovertemplate="DKA_DKT: %{x:.3f}<br>Count: %{y}<extra></extra>",
        ))
    if kmer_dnm:
        fig.add_trace(go.Histogram(
            x=kmer_dnm, nbinsx=30,
            marker_color="#4C78A8", opacity=0.85,
            name="Kmer DNM (DKA &gt; 0)",
            hovertemplate="DKA_DKT: %{x:.3f}<br>Count: %{y}<extra></extra>",
        ))
    # Threshold line
    fig.add_vline(
        x=_HIGH_QUALITY_DKA_DKT_THRESHOLD,
        line_dash="dash", line_color="#E45756", line_width=2,
        annotation_text=f"Higher-quality threshold ({_HIGH_QUALITY_DKA_DKT_THRESHOLD})",
        annotation_position="top right",
        annotation_font=dict(size=11, color="#E45756"),
    )
    fig.update_layout(
        barmode="overlay",
        title=dict(
            text="DKA/DKT Ratio Distribution",
            font=dict(size=18),
        ),
        xaxis_title="DKA_DKT Ratio",
        yaxis_title="Number of Variants",
        template="plotly_white",
        height=400,
        margin=dict(t=60, b=40),
        legend=dict(orientation="h", yanchor="bottom", y=1.02,
                    xanchor="right", x=1),
    )
    return _plotly_div(fig, div_id)


def _make_dka_vs_dkt_scatter(variants, div_id="scatter-plot"):
    """Create a scatter plot of DKA vs DKT colored by DKA_DKT ratio.

    The variant list is capped at ``_SCATTER_MAX_POINTS`` (DE_NOVO first) to
    keep the serialised Plotly JSON compact and ensure reliable rendering.
    """
    import plotly.graph_objects as go

    used, was_trimmed = _downsample_variants(variants, _SCATTER_MAX_POINTS)

    fig = go.Figure()
    fig.add_trace(go.Scatter(
        x=[v["dkt"] for v in used],
        y=[v["dka"] for v in used],
        mode="markers",
        marker=dict(
            size=[max(6, min(30, v["dku"] * 3)) for v in used],
            color=[v["dka_dkt"] for v in used],
            colorscale="Viridis",
            colorbar=dict(title="DKA_DKT"),
            showscale=True,
            line=dict(width=1, color="#333"),
        ),
        text=[v["label"] for v in used],
        hovertemplate=(
            "<b>%{text}</b><br>"
            "DKT: %{x}<br>DKA: %{y}<br>"
            "DKU: %{customdata[0]}<br>"
            "DKA_DKT: %{customdata[1]:.4f}<br>"
            "Call: %{customdata[2]}"
            "<extra></extra>"
        ),
        customdata=[[v["dku"], v["dka_dkt"], v["call"]] for v in used],
    ))

    # Add DKA=10 threshold line (minimum high-quality evidence requirement)
    fig.add_hline(y=10, line_dash="dot", line_color="#E45756", line_width=1.5,
                  annotation_text="DKA ≥ 10 threshold",
                  annotation_position="right",
                  annotation_font=dict(size=10, color="#E45756"))

    title_text = "DKA vs. DKT (size = DKU, color = DKA_DKT ratio)"
    if was_trimmed:
        title_text += (
            f"<br><sup>Showing {len(used)} of {len(variants)} variants "
            "(DE_NOVO first)</sup>"
        )

    fig.update_layout(
        title=dict(
            text=title_text,
            font=dict(size=18),
        ),
        xaxis_title="DKT (Total Spanning Fragments)",
        yaxis_title="DKA (Alt-Supporting Fragments with Unique K-mers)",
        template="plotly_white",
        height=500,
        margin=dict(t=60, b=40),
    )
    return _plotly_div(fig, div_id)


def _make_evidence_heatmap(variants, div_id="heatmap-plot"):
    """Evidence heatmap with automatic cluster-summary mode for large datasets.

    **Individual mode** (≤ ``_HEATMAP_MAX_ROWS`` variants)
        Shows one row per variant with compact ``customdata`` hover.

    **Cluster-summary mode** (> ``_HEATMAP_MAX_ROWS`` variants)
        Runs k-means (k = ``_HEATMAP_N_CLUSTERS``) on all variants using
        the 8 z-scored evidence features.  One row per cluster is displayed
        showing the centroid z-scores, sorted by DE_NOVO fraction (descending)
        so the most de-novo-enriched cluster appears at the top.  Each row is
        labelled with the cluster size and its DE_NOVO percentage.  The plot
        is rendered as a static image (no hover/zoom) to minimise HTML size —
        at 100k variants the entire cluster-summary section is < 50 KB.
    """
    import plotly.graph_objects as go

    fields = [
        "dku", "dkt", "dka", "dku_dkt", "dka_dkt",
        "max_pkc", "avg_pkc", "min_pkc",
    ]
    display_fields = [
        "DKU", "DKT", "DKA", "DKU_DKT", "DKA_DKT",
        "MAX_PKC", "AVG_PKC", "MIN_PKC",
    ]
    n_cols = len(fields)
    n_total = len(variants)

    # --- Z-score normalise per column over the full variant list -------------
    raw_all = [[v[f] for f in fields] for v in variants]
    z_all = [[0.0] * n_cols for _ in range(n_total)]
    for c in range(n_cols):
        col_vals = [raw_all[r][c] for r in range(n_total)]
        mean_val = stats.mean(col_vals) if col_vals else 0.0
        std_val = stats.pstdev(col_vals) if col_vals else 1.0
        if std_val == 0.0:
            std_val = 1.0
        for r in range(n_total):
            z_all[r][c] = (raw_all[r][c] - mean_val) / std_val

    if n_total > _HEATMAP_MAX_ROWS:
        # ── Cluster-summary mode ─────────────────────────────────────────────
        k = min(_HEATMAP_N_CLUSTERS, n_total)
        cluster_ids = _kmeans_cluster(z_all, k)

        # Aggregate per cluster
        cluster_data = {}
        for i, cl in enumerate(cluster_ids):
            if cl not in cluster_data:
                cluster_data[cl] = {"indices": [], "denovo": 0}
            cluster_data[cl]["indices"].append(i)
            if variants[i]["call"] == "DE_NOVO":
                cluster_data[cl]["denovo"] += 1

        # Sort clusters by DE_NOVO fraction descending
        sorted_clusters = sorted(
            cluster_data.items(),
            key=lambda kv: kv[1]["denovo"] / len(kv[1]["indices"]),
            reverse=True,
        )

        z_data = []
        y_labels = []
        hover_text = []
        raw_centroids = []
        for rank, (cl_id, info) in enumerate(sorted_clusters, start=1):
            idx = info["indices"]
            # idx is non-empty by construction: _kmeans_cluster assigns every
            # sample to exactly one cluster and empty clusters are skipped.
            centroid_z = [
                sum(z_all[i][c] for i in idx) / len(idx) for c in range(n_cols)
            ]
            centroid_raw = [
                sum(raw_all[i][c] for i in idx) / len(idx) for c in range(n_cols)
            ]
            pct = 100.0 * info["denovo"] / len(idx)
            z_data.append(centroid_z)
            raw_centroids.append(centroid_raw)
            y_labels.append(
                f"Cluster {rank}  (n={len(idx):,}, {pct:.0f}% DE_NOVO)"
            )
            hover_text.append([
                f"<b>Cluster {rank}</b><br>"
                f"n={len(idx):,} variants, {pct:.0f}% DE_NOVO<br>"
                f"{display_fields[c]}: {centroid_raw[c]:.3g} (centroid)<br>"
                f"Z-score: {centroid_z[c]:.2f}"
                for c in range(n_cols)
            ])

        title_text = (
            f"Per-Variant Evidence Heatmap — {n_total:,} variants summarised in "
            f"{len(sorted_clusters)} k-means clusters (k={k})<br>"
            f"<sup>Each row = cluster centroid; sorted by DE_NOVO fraction ↓</sup>"
        )
        is_cluster_mode = True

    else:
        # ── Individual-row mode ───────────────────────────────────────────────
        z_data = z_all
        raw_centroids = raw_all
        y_labels = [v["label"] for v in variants]
        hover_text = None
        title_text = "Per-Variant Evidence Heatmap (Z-score normalized)"
        is_cluster_mode = False

    # Build heatmap trace
    if is_cluster_mode:
        # Cluster rows: compact hover strings (few rows → string cost is tiny)
        heatmap = go.Heatmap(
            z=z_data,
            x=display_fields,
            y=y_labels,
            colorscale="RdBu_r",
            zmid=0,
            text=hover_text,
            hoverinfo="text",
            colorbar=dict(title="Z-score<br>(centroid)"),
        )
    else:
        # Individual rows: float customdata avoids per-cell label duplication
        heatmap = go.Heatmap(
            z=z_data,
            x=display_fields,
            y=y_labels,
            colorscale="RdBu_r",
            zmid=0,
            customdata=raw_centroids,
            hovertemplate=(
                "<b>%{y}</b><br>%{x}: %{customdata:.4g}"
                "<br>Z-score: %{z:.2f}<extra></extra>"
            ),
            colorbar=dict(title="Z-score"),
        )

    n_rows_shown = len(z_data)
    # Row height: 40 px max (readability), 20 px min (avoid over-compression),
    # targeting ~400 px for small cluster counts; cap total at 2000 px to stay
    # within browser canvas limits for the individual-row mode.
    row_px = min(40, max(20, 400 // max(1, n_rows_shown)))
    height = min(2000, max(400, row_px * n_rows_shown + 120))

    fig = go.Figure(data=heatmap)
    fig.update_layout(
        title=dict(text=title_text, font=dict(size=16 if is_cluster_mode else 18)),
        template="plotly_white",
        height=height,
        margin=dict(t=100 if is_cluster_mode else 80, b=40, l=300),
        yaxis=dict(autorange="reversed"),
    )

    # Cluster-summary plots are static: disabling interactivity removes all
    # browser-side event handler JS and keeps the HTML footprint small.
    cfg = {"staticPlot": True} if is_cluster_mode else {"responsive": True}
    return _plotly_div(fig, div_id, config=cfg)


def _make_pkc_boxplot(variants, div_id="pkc-box-plot"):
    """Create box plots of ALT-specific PKC metrics by call type.

    Uses ALT-allele parental k-mer counts (PKC_ALT), not total PKC, because
    reference-allele k-mers are present in parents for all variants.  Only
    the ALT-allele k-mer abundance in parents distinguishes de novo (absent)
    from inherited (present) variants.
    """
    import plotly.graph_objects as go

    denovo = [v for v in variants if v["call"] == "DE_NOVO"]
    inherited = [v for v in variants if v["call"] != "DE_NOVO"]

    fig = go.Figure()
    for label_group, group, color in [
        ("De Novo", denovo, "#54A24B"),
        ("Inherited", inherited, "#E45756"),
    ]:
        if not group:
            continue
        for metric, name in [
            ("max_pkc_alt", "MAX_PKC_ALT"),
            ("avg_pkc_alt", "AVG_PKC_ALT"),
            ("min_pkc_alt", "MIN_PKC_ALT"),
        ]:
            fig.add_trace(go.Box(
                y=[v[metric] for v in group],
                name=f"{name}<br>({label_group})",
                marker_color=color,
                boxmean=True,
            ))

    fig.update_layout(
        title=dict(
            text="ALT-Allele Parental K-mer Count (PKC_ALT) by Call Type",
            font=dict(size=18),
        ),
        yaxis_title="ALT-Allele K-mer Count in Parents",
        template="plotly_white",
        height=450,
        margin=dict(t=60, b=40),
        showlegend=False,
    )
    return _plotly_div(fig, div_id)


def _make_pkc_vs_dka_dkt_scatter(variants, div_id="pkc-scatter-plot"):
    """Create AVG_PKC_ALT vs DKA_DKT scatter plot.

    Uses ALT-allele parental k-mer count (avg_pkc_alt) on the y-axis.
    For genuine de novos, avg_pkc_alt should be near zero because the
    ALT-allele k-mers are absent from both parents.  For inherited variants,
    avg_pkc_alt is non-zero because the ALT allele is present in at least
    one parent.  This demonstrates the null hypothesis (parental coverage
    gap) can be rejected when inherited variants show high avg_pkc_alt.

    The variant list is capped at ``_SCATTER_MAX_POINTS`` to bound the
    serialised data size.
    """
    import plotly.graph_objects as go

    used, was_trimmed = _downsample_variants(variants, _SCATTER_MAX_POINTS)

    colors = ["#54A24B" if v["call"] == "DE_NOVO" else "#E45756"
              for v in used]

    fig = go.Figure()
    fig.add_trace(go.Scatter(
        x=[v["dka_dkt"] for v in used],
        y=[v["avg_pkc_alt"] for v in used],
        mode="markers",
        marker=dict(size=10, color=colors, line=dict(width=1, color="#333")),
        text=[v["label"] for v in used],
        customdata=[v["call"] for v in used],
        hovertemplate=(
            "<b>%{text}</b><br>"
            "Call: %{customdata}<br>"
            "DKA_DKT: %{x:.4f}<br>AVG_PKC_ALT: %{y:.1f}"
            "<extra></extra>"
        ),
    ))
    # DKA_DKT threshold line
    fig.add_vline(
        x=_HIGH_QUALITY_DKA_DKT_THRESHOLD,
        line_dash="dash", line_color="#ccc", line_width=1,
        annotation_text=f"DKA_DKT \u2265 {_HIGH_QUALITY_DKA_DKT_THRESHOLD}",
        annotation_position="top right",
        annotation_font=dict(size=10, color="#666"),
    )

    subtitle = (
        "De novo (green) should cluster at low AVG_PKC_ALT; "
        "inherited (red) at high AVG_PKC_ALT"
    )
    if was_trimmed:
        subtitle += (
            f" — showing {len(used)} of {len(variants)} variants (DE_NOVO first)"
        )

    fig.update_layout(
        title=dict(
            text=f"AVG_PKC_ALT vs. DKA_DKT Ratio<br><sup>{subtitle}</sup>",
            font=dict(size=16),
        ),
        xaxis_title="DKA_DKT Ratio",
        yaxis_title="AVG_PKC_ALT (Avg. ALT-Allele K-mer Count in Parents)",
        template="plotly_white",
        height=450,
        margin=dict(t=80, b=40),
    )
    return _plotly_div(fig, div_id)


def _make_contamination_bar(variants, kraken2_data, div_id="contamination-plot"):
    """Create stacked bar chart of Kraken2 read classification fractions.

    Filtered to variants with **NHF >= 0.05** (putative contamination) among
    kmer-DNM candidates (stage >= 1) so the plot focuses on variants where
    contamination is actually detected.  This provides an informative view
    of the contamination profile rather than showing all variants (most of
    which have negligible NHF).
    """
    import plotly.graph_objects as go

    if not kraken2_data:
        return None

    # Match kraken2_data to variants by label
    kraken_map = {r["label"]: r for r in kraken2_data}

    labels_with_data = []
    hlf_vals = []
    nhf_vals = []
    ucf_vals = []
    uf_vals = []

    for v in variants:
        # Only include kmer-DNM candidates (stage >= 1) with NHF >= 0.05
        if v.get("stage", 0) < 1:
            continue
        nhf_val = v.get("dka_nhf")
        if nhf_val is None or nhf_val < _NHF_CONTAMINATION_THRESHOLD:
            continue
        k = kraken_map.get(v["label"])
        if k:
            labels_with_data.append(v["label"])
            hlf_vals.append(k.get("DKA_HLF", 0))
            nhf_vals.append(k.get("DKA_NHF", 0))
            ucf_vals.append(k.get("DKA_UCF", 0))
            uf_vals.append(k.get("DKA_UF", 0))

    if not labels_with_data:
        return None

    fig = go.Figure()
    fig.add_trace(go.Bar(
        x=labels_with_data, y=hlf_vals,
        name="Human Lineage (DKA_HLF)", marker_color="#4C78A8",
    ))
    fig.add_trace(go.Bar(
        x=labels_with_data, y=nhf_vals,
        name="Non-Human (DKA_NHF)", marker_color="#E45756",
    ))
    fig.add_trace(go.Bar(
        x=labels_with_data, y=ucf_vals,
        name="UniVec Core (DKA_UCF)", marker_color="#F58518",
    ))
    fig.add_trace(go.Bar(
        x=labels_with_data, y=uf_vals,
        name="Unclassified (DKA_UF)", marker_color="#BAB0AC",
    ))

    fig.update_layout(
        barmode="stack",
        title=dict(
            text=(f"Kraken2 Read Classification — Putative Contamination Only "
                  f"(NHF &ge; {_NHF_CONTAMINATION_THRESHOLD}, n={len(labels_with_data)})"),
            font=dict(size=16),
        ),
        yaxis_title="Fraction of DKA Reads",
        template="plotly_white",
        height=450,
        margin=dict(t=60, b=120),
        xaxis_tickangle=-45,
        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1),
    )
    return _plotly_div(fig, div_id)


def _make_contamination_funnel(stratification, variants,
                               div_id="contam-funnel-plot"):
    """Bar chart showing the proportion of variants with NHF >= 0.05 at each stage.

    This shows how contamination (as measured by DKA_NHF) is distributed at
    each filtering stage, allowing the reviewer to see:
    - Baseline contamination rate in the full input set
    - How contamination is progressively removed by other filters
    - The residual contamination in the final passing set

    Returns ``None`` when no Kraken2 NHF data is available.
    """
    import plotly.graph_objects as go

    if not stratification["has_nhf_data"]:
        return None

    counts = stratification["counts"]
    short_labels = stratification["short_labels"]
    colors = stratification["colors"]

    # For each stage, compute the number and proportion of variants with
    # NHF >= 0.05 among those that pass that stage
    contam_counts = []
    contam_pcts = []
    for stage_idx in range(6):
        stage_variants = [v for v in variants if v.get("stage", 0) >= stage_idx]
        n_contam = sum(
            1 for v in stage_variants
            if v.get("dka_nhf") is not None
            and v["dka_nhf"] >= _NHF_CONTAMINATION_THRESHOLD
        )
        total = len(stage_variants)
        contam_counts.append(n_contam)
        contam_pcts.append(100.0 * n_contam / total if total > 0 else 0.0)

    if all(c == 0 for c in contam_counts):
        return None

    fig = go.Figure()
    fig.add_trace(go.Bar(
        x=short_labels,
        y=contam_pcts,
        marker_color=colors,
        text=[f"{c} ({p:.1f}%)" for c, p in zip(contam_counts, contam_pcts)],
        textposition="outside",
        hovertemplate=(
            "<b>%{x}</b><br>"
            "Contaminated: %{customdata[0]}<br>"
            "Total at stage: %{customdata[1]}<br>"
            "Proportion: %{y:.1f}%"
            "<extra></extra>"
        ),
        customdata=list(zip(contam_counts, counts)),
    ))

    fig.update_layout(
        title=dict(
            text=("Contamination Prevalence by Stage<br>"
                  f"<sup>Proportion of variants with DKA_NHF &ge; "
                  f"{_NHF_CONTAMINATION_THRESHOLD} at each filter level</sup>"),
            font=dict(size=16),
        ),
        yaxis_title="% Variants with NHF ≥ 0.05",
        template="plotly_white",
        height=400,
        margin=dict(t=80, b=60),
        showlegend=False,
    )
    return _plotly_div(fig, div_id)


def _make_discovery_region_scatter(regions, div_id="disc-scatter-plot"):
    """Create scatter plot of discovery regions: reads vs k-mers.

    Capped at ``_SCATTER_MAX_POINTS`` total regions to bound serialised size.
    """
    import plotly.graph_objects as go

    class_colors = {"SMALL": "#4C78A8", "AMBIGUOUS": "#F58518", "SV": "#E45756"}

    # Cap total regions; SV regions are highest priority, then AMBIGUOUS, SMALL
    used_regions = regions
    was_trimmed = False
    if len(regions) > _SCATTER_MAX_POINTS:
        sv = [r for r in regions if r.get("class") == "SV"]
        amb = [r for r in regions if r.get("class") == "AMBIGUOUS"]
        small = [r for r in regions if r.get("class") == "SMALL"]
        remaining = _SCATTER_MAX_POINTS
        keep_sv = sv[:remaining]
        remaining -= len(keep_sv)
        keep_amb = amb[:remaining]
        remaining -= len(keep_amb)
        keep_small = small[:remaining]
        used_regions = keep_sv + keep_amb + keep_small
        was_trimmed = True

    fig = go.Figure()
    for cls in ["SMALL", "AMBIGUOUS", "SV"]:
        cls_regions = [r for r in used_regions if r.get("class") == cls]
        if not cls_regions:
            continue
        fig.add_trace(go.Scatter(
            x=[r["reads"] for r in cls_regions],
            y=[r["unique_kmers"] for r in cls_regions],
            mode="markers",
            name=cls,
            marker=dict(
                size=[max(6, min(25, r.get("max_clip_len", 0) / 4 + 6))
                      for r in cls_regions],
                color=class_colors.get(cls, "#999"),
                line=dict(width=1, color="#333"),
            ),
            text=[
                f"{r['chrom']}:{r['start']+1}-{r['end']}<br>"
                f"Size: {r['size']}bp<br>"
                f"MaxClip: {r.get('max_clip_len', 0)}<br>"
                f"Class: {r.get('class', 'N/A')}"
                for r in cls_regions
            ],
            hovertemplate="%{text}<extra></extra>",
        ))

    title_text = "Discovery Regions: Reads vs. Distinct K-mers"
    if was_trimmed:
        title_text += (
            f"<br><sup>Showing {len(used_regions)} of {len(regions)} regions</sup>"
        )

    fig.update_layout(
        title=dict(text=title_text, font=dict(size=18)),
        xaxis_title="Supporting Reads",
        yaxis_title="Distinct Proband-Unique K-mers",
        template="plotly_white",
        height=450,
        margin=dict(t=60, b=40),
    )
    return _plotly_div(fig, div_id)


def _make_discovery_size_histogram(regions, div_id="disc-size-plot"):
    """Create histogram of discovery region sizes by class."""
    import plotly.graph_objects as go

    class_colors = {"SMALL": "#4C78A8", "AMBIGUOUS": "#F58518", "SV": "#E45756"}

    fig = go.Figure()
    for cls in ["SMALL", "AMBIGUOUS", "SV"]:
        cls_regions = [r for r in regions if r.get("class") == cls]
        if not cls_regions:
            continue
        fig.add_trace(go.Histogram(
            x=[r["size"] for r in cls_regions],
            name=cls,
            marker_color=class_colors.get(cls, "#999"),
            opacity=0.75,
        ))

    fig.update_layout(
        barmode="overlay",
        title=dict(
            text="Discovery Region Size Distribution",
            font=dict(size=18),
        ),
        xaxis_title="Region Size (bp)",
        yaxis_title="Count",
        template="plotly_white",
        height=400,
        margin=dict(t=60, b=40),
    )
    return _plotly_div(fig, div_id)


def _make_sv_evidence_chart(regions, div_id="sv-evidence-plot"):
    """Create grouped bar chart of SV evidence per region."""
    import plotly.graph_objects as go

    # Only include regions with some SV evidence
    sv_regions = [r for r in regions
                  if r.get("split_reads", 0) > 0
                  or r.get("discordant_pairs", 0) > 0
                  or r.get("unmapped_mates", 0) > 0
                  or r.get("max_clip_len", 0) > 50]
    if not sv_regions:
        return None

    labels = [f"{r['chrom']}:{r['start']+1}-{r['end']}" for r in sv_regions]

    fig = go.Figure()
    for field, name, color in [
        ("split_reads", "Split Reads", "#4C78A8"),
        ("discordant_pairs", "Discordant Pairs", "#E45756"),
        ("unmapped_mates", "Unmapped Mates", "#F58518"),
    ]:
        fig.add_trace(go.Bar(
            x=labels, y=[r.get(field, 0) for r in sv_regions],
            name=name, marker_color=color,
        ))

    fig.update_layout(
        barmode="group",
        title=dict(
            text="Structural Variant Evidence by Region",
            font=dict(size=18),
        ),
        yaxis_title="Count",
        template="plotly_white",
        height=400,
        margin=dict(t=60, b=120),
        xaxis_tickangle=-45,
    )
    return _plotly_div(fig, div_id)


def _make_threshold_sensitivity(variants, div_id="threshold-plot"):
    """Create threshold sensitivity plot: variants passing at each DKA_DKT."""
    import plotly.graph_objects as go

    dka_dkt_values = sorted([v["dka_dkt"] for v in variants])
    thresholds = [i * 0.01 for i in range(0, 101)]
    passing = []
    for t in thresholds:
        passing.append(sum(1 for v in dka_dkt_values if v >= t))

    fig = go.Figure()
    fig.add_trace(go.Scatter(
        x=thresholds, y=passing,
        mode="lines+markers",
        marker=dict(size=3, color="#4C78A8"),
        line=dict(color="#4C78A8", width=2),
        hovertemplate="Threshold: %{x:.2f}<br>Passing: %{y}<extra></extra>",
    ))
    fig.add_vline(
        x=_HIGH_QUALITY_DKA_DKT_THRESHOLD,
        line_dash="dash", line_color="#E45756", line_width=2,
        annotation_text=f"{_HIGH_QUALITY_DKA_DKT_THRESHOLD}",
        annotation_position="top right",
        annotation_font=dict(size=11, color="#E45756"),
    )
    fig.update_layout(
        title=dict(
            text="DKA_DKT Threshold Sensitivity",
            font=dict(size=18),
        ),
        xaxis_title="DKA_DKT Threshold",
        yaxis_title="Variants Passing",
        template="plotly_white",
        height=400,
        margin=dict(t=60, b=40),
    )
    return _plotly_div(fig, div_id)


# Pre-compiled regex for parsing alleles from variant labels.
# Format: "<chrom>:<pos> <ref>><alt>"
_VARIANT_LABEL_PATTERN = re.compile(r"\s+(\S+)>(\S+)$")


# ---------------------------------------------------------------------------
# New genomics-specific plot generators
# ---------------------------------------------------------------------------


def _classify_variant_type(label):
    """Classify a variant as SNV, INS, DEL, or MNV from its label string.

    Label format: ``<chrom>:<pos> <ref>><alt>``
    Returns one of "SNV", "INS", "DEL", "MNV", or "Other".
    """
    m = _VARIANT_LABEL_PATTERN.search(label)
    if not m:
        return "Other"
    ref, alt = m.group(1), m.group(2)
    if len(ref) == 1 and len(alt) == 1:
        return "SNV"
    if len(ref) < len(alt):
        return "INS"
    if len(ref) > len(alt):
        return "DEL"
    return "MNV"


def _make_variant_type_breakdown(variants, div_id="variant-type-plot"):
    """Create a grouped bar chart of variant types by call status.

    Shows the count of SNVs, insertions, deletions, and MNVs for
    DE_NOVO vs. inherited variants.  A genomics reviewer expects the
    de novo SNV rate (~38/genome for Illumina WGS) to be consistent
    with published trio studies.
    """
    import plotly.graph_objects as go

    type_order = ["SNV", "INS", "DEL", "MNV", "Other"]
    denovo_counts = {t: 0 for t in type_order}
    inherited_counts = {t: 0 for t in type_order}

    for v in variants:
        vtype = _classify_variant_type(v["label"])
        if v["call"] == "DE_NOVO":
            denovo_counts[vtype] += 1
        else:
            inherited_counts[vtype] += 1

    # Only include types that appear in the data
    present_types = [t for t in type_order
                     if denovo_counts[t] > 0 or inherited_counts[t] > 0]
    if not present_types:
        return None

    fig = go.Figure()
    fig.add_trace(go.Bar(
        x=present_types,
        y=[denovo_counts[t] for t in present_types],
        name="De Novo",
        marker_color="#54A24B",
        text=[denovo_counts[t] for t in present_types],
        textposition="outside",
    ))
    fig.add_trace(go.Bar(
        x=present_types,
        y=[inherited_counts[t] for t in present_types],
        name="Inherited / Unclear",
        marker_color="#E45756",
        text=[inherited_counts[t] for t in present_types],
        textposition="outside",
    ))
    fig.update_layout(
        barmode="group",
        title=dict(
            text="Variant Type Breakdown by Call Status",
            font=dict(size=18),
        ),
        xaxis_title="Variant Type",
        yaxis_title="Count",
        template="plotly_white",
        height=400,
        margin=dict(t=60, b=40),
        legend=dict(orientation="h", yanchor="bottom", y=1.02,
                    xanchor="right", x=1),
    )
    return _plotly_div(fig, div_id)


def _make_chromosomal_distribution(variants, div_id="chrom-dist-plot"):
    """Create a stacked bar chart of variant counts per chromosome.

    A genomics reviewer will check that de novo candidates are
    distributed across chromosomes as expected rather than clustered
    on a single chromosome (which may indicate a systematic error).
    """
    import plotly.graph_objects as go

    # Build chromosome order: numeric 1-22, then X, Y, M, then rest
    def _chrom_sort_key(chrom):
        c = chrom.replace("chr", "").upper()
        order = {"X": 23, "Y": 24, "M": 25, "MT": 25}
        try:
            return int(c)
        except ValueError:
            return order.get(c, 99)

    chrom_denovo = {}
    chrom_inherited = {}
    for v in variants:
        # Label: "chr8:40003391 A>C"
        chrom_part = v["label"].split(":")[0]
        if v["call"] == "DE_NOVO":
            chrom_denovo[chrom_part] = chrom_denovo.get(chrom_part, 0) + 1
        else:
            chrom_inherited[chrom_part] = chrom_inherited.get(chrom_part, 0) + 1

    all_chroms = sorted(
        set(chrom_denovo) | set(chrom_inherited),
        key=_chrom_sort_key,
    )
    if not all_chroms:
        return None

    fig = go.Figure()
    fig.add_trace(go.Bar(
        x=all_chroms,
        y=[chrom_denovo.get(c, 0) for c in all_chroms],
        name="De Novo",
        marker_color="#54A24B",
    ))
    fig.add_trace(go.Bar(
        x=all_chroms,
        y=[chrom_inherited.get(c, 0) for c in all_chroms],
        name="Inherited / Unclear",
        marker_color="#E45756",
    ))
    fig.update_layout(
        barmode="stack",
        title=dict(
            text="Chromosomal Distribution of Candidate Variants",
            font=dict(size=18),
        ),
        xaxis_title="Chromosome",
        yaxis_title="Variant Count",
        template="plotly_white",
        height=400,
        margin=dict(t=60, b=40),
        legend=dict(orientation="h", yanchor="bottom", y=1.02,
                    xanchor="right", x=1),
    )
    return _plotly_div(fig, div_id)


# ---------------------------------------------------------------------------
# HTML template
# ---------------------------------------------------------------------------

_HTML_TEMPLATE = """\
<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>kmer-denovo — Interactive Report</title>
<script>{{ plotly_bundle | safe }}</script>
<style>
  * { box-sizing: border-box; margin: 0; padding: 0; }
  body {
    font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Helvetica, Arial, sans-serif;
    line-height: 1.6; color: #333; background: #fafafa;
  }
  .container { max-width: 1200px; margin: 0 auto; padding: 20px; }
  h1 { font-size: 2em; margin-bottom: 10px; color: #1a1a2e; }
  h2 {
    font-size: 1.5em; margin: 40px 0 15px; padding-bottom: 8px;
    border-bottom: 2px solid #4C78A8; color: #1a1a2e;
  }
  h3 { font-size: 1.2em; margin: 20px 0 10px; color: #333; }
  p, .description { margin-bottom: 15px; color: #555; font-size: 0.95em; }
  .subtitle { color: #666; font-size: 1.1em; margin-bottom: 25px; }
  .metric-cards {
    display: grid;
    grid-template-columns: repeat(auto-fit, minmax(180px, 1fr));
    gap: 15px; margin: 20px 0;
  }
  .metric-card {
    background: white; border-radius: 8px; padding: 20px;
    box-shadow: 0 2px 8px rgba(0,0,0,0.08);
    text-align: center; border-top: 3px solid #4C78A8;
  }
  .metric-card .value {
    font-size: 2em; font-weight: 700; color: #1a1a2e;
    display: block; margin-bottom: 5px;
  }
  .metric-card .label { font-size: 0.85em; color: #666; }
  .metric-card.green { border-top-color: #54A24B; }
  .metric-card.red { border-top-color: #E45756; }
  .metric-card.orange { border-top-color: #F58518; }
  .plot-container {
    background: white; border-radius: 8px; padding: 20px; margin: 20px 0;
    box-shadow: 0 2px 8px rgba(0,0,0,0.08);
  }
  .plot-div { width: 100%; }
  table {
    width: 100%; border-collapse: collapse; margin: 15px 0;
    font-size: 0.9em; background: white;
    box-shadow: 0 2px 8px rgba(0,0,0,0.08); border-radius: 8px;
    overflow: hidden;
  }
  th {
    background: #4C78A8; color: white; padding: 12px 10px;
    text-align: left; font-weight: 600;
  }
  td { padding: 10px; border-bottom: 1px solid #eee; }
  tr:hover td { background: #f5f8ff; }
  .badge {
    display: inline-block; padding: 2px 8px; border-radius: 4px;
    font-size: 0.8em; font-weight: 600;
  }
  .badge-green { background: #d4edda; color: #155724; }
  .badge-red { background: #f8d7da; color: #721c24; }
  .badge-orange { background: #fff3cd; color: #856404; }
  .method-box {
    background: #f0f4f8; border-left: 4px solid #4C78A8;
    padding: 20px; margin: 20px 0; border-radius: 0 8px 8px 0;
  }
  .method-box h3 { margin-top: 0; }
  .section-rationale {
    background: #f8f9fa; padding: 12px 16px; margin: 10px 0;
    border-radius: 6px; font-size: 0.9em; color: #555;
    border-left: 3px solid #72B7B2;
  }
  .plot-caption {
    font-size: 0.85em; color: #666; margin-top: -8px;
    margin-bottom: 12px; padding: 0 8px; font-style: italic;
  }
  .summary-card {
    background: linear-gradient(135deg, #f0f4f8 0%, #e8efff 100%);
    border-radius: 10px; padding: 24px 28px; margin: 20px 0;
    box-shadow: 0 2px 12px rgba(0,0,0,0.06);
    border-left: 5px solid #4C78A8;
  }
  .summary-card h3 { margin-top: 0; color: #1a1a2e; }
  .summary-card p { color: #333; margin-bottom: 8px; }
  .summary-card ol, .summary-card ul {
    margin: 8px 0 8px 24px; color: #333; font-size: 0.95em;
  }
  .summary-card li { margin: 4px 0; }
  .stage-list { list-style: none; padding-left: 0; }
  .stage-list li {
    margin: 6px 0; padding-left: 28px; position: relative;
  }
  .stage-list li::before {
    content: ""; position: absolute; left: 0; top: 8px;
    width: 14px; height: 14px; border-radius: 3px;
  }
  .stage-list li.s0::before { background: #4C78A8; }
  .stage-list li.s1::before { background: #F58518; }
  .stage-list li.s2::before { background: #E45756; }
  .stage-list li.s3::before { background: #72B7B2; }
  .stage-list li.s4::before { background: #EECA3B; }
  .stage-list li.s5::before { background: #54A24B; }
  footer {
    margin-top: 40px; padding: 20px 0; border-top: 1px solid #ddd;
    text-align: center; color: #999; font-size: 0.85em;
  }
</style>
</head>
<body>
<div class="container">

<!-- ═══════════════════════════════════════════════════════════════════
     Header
     ═══════════════════════════════════════════════════════════════════ -->
<h1>🧬 kmer-denovo — Interactive Report</h1>
<p class="subtitle">
  De novo variant curation using alignment-independent k-mer analysis
  {% if mode == "vcf" %}(VCF mode){% elif mode == "discovery" %}(Discovery mode){% else %}(Combined){% endif %}
</p>

<!-- ═══════════════════════════════════════════════════════════════════
     Section 1: Pipeline Summary & Executive Snapshot
     ═══════════════════════════════════════════════════════════════════ -->
<h2>1. Pipeline Summary</h2>

{% if vcf_metrics %}
<div class="summary-card">
  <h3>What this pipeline does</h3>
  <p>
    This report summarises the results of <strong>kmer-denovo-filter</strong>
    applied to <strong>{{ "{:,}".format(vcf_metrics.get("total_variants", 0)) }}
    putative <em>de novo</em> variants</strong> supplied as input in a VCF
    (e.g. from a trio variant caller such as DeepVariant&nbsp;+&nbsp;GLnexus
    or DeNovoGear).  Without orthogonal evidence, the majority of these
    putative DNMs are <strong>simple miscalled variants</strong> — typically
    inherited alleles whose parental support was missed by the original
    caller, alignment artefacts, or sequence drawn from non-human
    contamination.  k-mer analysis offers an alignment-independent way to
    re-examine each candidate using the raw read sequence.
  </p>
  <p style="margin-top:14px;"><strong>Stratification used throughout this report:</strong></p>
  <ul class="stage-list">
    <li class="s0"><strong>Putative denovo</strong> &mdash; every variant from the input VCF.</li>
    <li class="s1"><strong>Putative kmer denovo</strong> &mdash; &ge; 1 child fragment carries an
      ALT-supporting child-unique k-mer (DKA &gt; 0).</li>
    <li class="s2"><strong>Putative kmer denovo (stronger)</strong> &mdash; DKA &ge; 5
      (at least 5 ALT-supporting fragments with unique k-mers).</li>
    <li class="s3"><strong>Higher-quality denovo</strong> &mdash; DKA_DKT &gt; 0.1 (more than
      10% of spanning fragments carry the unique-k-mer ALT signal).</li>
    <li class="s4"><strong>Higher-quality, parental confirmed</strong> &mdash;
      MAX_PKC_ALT &lt; 1 (ALT-allele k-mers absent from parents).</li>
    <li class="s5"><strong>HQ, not contamination</strong> &mdash; all prior filters
      <em>plus</em> the non-human read fraction NHF&nbsp;&lt;&nbsp;0.05.</li>
  </ul>
</div>
{% elif disc_metrics %}
<div class="summary-card">
  <h3>What this pipeline does</h3>
  <p>
    This report summarises a <strong>VCF-free discovery</strong> run of
    kmer-denovo-filter.  Instead of starting from putative DNMs in a VCF,
    the pipeline scans the proband BAM for genomic regions enriched for
    proband-unique k-mers (k-mers absent from both parents) and reports
    these as candidate <em>de novo</em> events.  The intent is to
    re-discover known DNMs and surface novel ones, including structural
    variants that VCF-based callers may miss.
  </p>
</div>
{% endif %}

<h3>At a glance</h3>
<div class="metric-cards">
  {% if stratification %}
  <div class="metric-card">
    <span class="value">{{ "{:,}".format(stratification.counts[0]) }}</span>
    <span class="label">Putative denovo (input)</span>
  </div>
  <div class="metric-card orange">
    <span class="value">{{ "{:,}".format(stratification.counts[1]) }}</span>
    <span class="label">Putative kmer denovo (DKA &gt; 0)</span>
  </div>
  <div class="metric-card red">
    <span class="value">{{ "{:,}".format(stratification.counts[2]) }}</span>
    <span class="label">Kmer denovo (DKA &ge; 5)</span>
  </div>
  <div class="metric-card">
    <span class="value">{{ "{:,}".format(stratification.counts[3]) }}</span>
    <span class="label">Higher-quality (DKA_DKT &gt; 0.1)</span>
  </div>
  <div class="metric-card green">
    <span class="value">{{ "{:,}".format(stratification.counts[4]) }}</span>
    <span class="label">HQ parental (MAX_PKC_ALT &lt; 1)</span>
  </div>
  <div class="metric-card green">
    <span class="value">{{ "{:,}".format(stratification.counts[5]) }}{% if not stratification.has_nhf_data %}*{% endif %}</span>
    <span class="label">HQ + non-contamination
      {% if not stratification.has_nhf_data %}<br><em>(no NHF data — same as parental)</em>{% endif %}
    </span>
  </div>
  {% endif %}
  {% if vcf_metrics %}
  <div class="metric-card">
    <span class="value">{{ "{:,}".format(vcf_metrics.get("total_child_kmers", 0)) }}</span>
    <span class="label">Total Child K-mers</span>
  </div>
  <div class="metric-card green">
    <span class="value">{{ "{:,}".format(vcf_metrics.get("child_unique_kmers", 0)) }}</span>
    <span class="label">Child-Unique K-mers</span>
  </div>
  {% endif %}

  {% if disc_metrics %}
  <div class="metric-card orange">
    <span class="value">{{ disc_metrics.get("candidate_regions", 0) }}</span>
    <span class="label">Candidate Regions</span>
  </div>
  <div class="metric-card">
    <span class="value">{{ "{:,}".format(disc_metrics.get("proband_unique_kmers", 0)) }}</span>
    <span class="label">Proband-Unique K-mers</span>
  </div>
  <div class="metric-card">
    <span class="value">{{ disc_metrics.get("informative_reads", 0) }}</span>
    <span class="label">Informative Reads</span>
  </div>
  {% endif %}
</div>

{% if strat_funnel_div %}
<div class="plot-container">
  {{ strat_funnel_div | safe }}
  <p class="plot-caption">
    Variant attrition through the four stratification stages.  Reading
    left-to-right shows how many of the input putative DNMs survive each
    progressively stricter filter.
  </p>
</div>
{% endif %}

{% if strat_sankey_div %}
<div class="plot-container">
  {{ strat_sankey_div | safe }}
  <p class="plot-caption">
    Variant flow diagram: each grey "Filtered" node shows how many variants
    were removed by the corresponding criterion.  This makes the
    selectivity of each filter explicit.
  </p>
</div>
{% endif %}

{% if sankey_div %}
<div class="plot-container">
  {{ sankey_div | safe }}
  <p class="plot-caption">
    K-mer-level flow: how the child's raw k-mer pool is reduced after
    parental subtraction.  This is k-mer-level (not variant-level) and is
    complementary to the variant funnel above.
  </p>
</div>
{% endif %}

<!-- ═══════════════════════════════════════════════════════════════════
     Section 2: Methods
     ═══════════════════════════════════════════════════════════════════ -->
<h2>2. Methods</h2>
<div class="method-box">
  <h3>Input &amp; rationale</h3>
  <p>
    {% if vcf_metrics %}
    The pipeline ingests a trio VCF of putative de novo variants (proband
    plus parental BAM/CRAM files).  Standard variant callers report many
    putative DNMs that are in fact inherited or artefactual; the goal of
    kmer-denovo-filter is to re-evaluate every candidate with
    alignment-independent k-mer evidence and to flag candidates whose
    supporting reads contain non-human sequence.
    {% else %}
    The pipeline scans the proband BAM directly (no input VCF), identifying
    contiguous genomic regions enriched for proband-unique k-mers as
    candidate de novo events.
    {% endif %}
  </p>
</div>

<div class="method-box">
  <h3>K-mer extraction &amp; parental subtraction</h3>
  <p>
    For each candidate variant we extract the set of k-mers present on
    reads spanning the locus in the proband, then query <a
    href="https://github.com/gmarcais/Jellyfish">Jellyfish</a> indexes
    built from each parent's sequencing reads.  K-mers absent from both
    parents are termed <strong>child-unique</strong>.  Per variant we
    record:
  </p>
  <ul style="margin-left:24px; color:#444; font-size:0.95em;">
    <li><strong>DKU</strong> &mdash; distinct child-unique k-mers at the locus.</li>
    <li><strong>DKT</strong> &mdash; total fragments spanning the locus.</li>
    <li><strong>DKA</strong> &mdash; fragments carrying ≥ 1 child-unique k-mer that supports the ALT allele.</li>
    <li><strong>DKA_DKT</strong> &mdash; DKA / DKT, analogous to a k-mer-based VAF.</li>
    <li><strong>PKC_ALT</strong> &mdash; ALT-allele k-mer count in each parent (MAX/AVG/MIN).</li>
  </ul>
</div>

<div class="method-box">
  <h3>Stratification cascade</h3>
  <p>
    Every putative DNM is placed into one of six cumulative stages.  These
    same stages are referenced by the funnel, the Sankey, the per-variant
    table, and the contamination plots so the report tells a single
    consistent story.
  </p>
  <ol style="margin-left:24px; color:#444; font-size:0.95em;">
    <li><strong>Putative denovo</strong> &mdash; every input variant (no filter).</li>
    <li><strong>Putative kmer denovo</strong> &mdash; <code>DKA &gt; 0</code>.  Removes the
      large majority of inherited / miscalled variants for which no proband
      fragment carries a child-unique ALT-supporting k-mer.</li>
    <li><strong>Putative kmer denovo (stronger)</strong> &mdash; <code>DKA &ge; 5</code>.
      Requires at least 5 ALT-supporting fragments with unique k-mers,
      reducing noise from single-read or low-depth signals.</li>
    <li><strong>Higher-quality denovo</strong> &mdash; <code>DKA_DKT &gt; 0.1</code>.  Requires
      that more than 10% of spanning fragments carry the unique-k-mer
      ALT signal &mdash; analogous to a minimum VAF requirement.</li>
    <li><strong>Higher-quality, parental confirmed</strong> &mdash;
      <code>MAX_PKC_ALT &lt; 1</code>.  Confirms that ALT-allele k-mers are
      absent from parental sequence, ruling out inherited variants missed
      by k-mer-level subtraction.</li>
    <li><strong>HQ, not contamination</strong> &mdash; all prior filters plus
      <code>DKA_NHF &lt; 0.05</code>.  Removes variants whose supporting
      reads are dominated (&ge; 5%) by Kraken2-classified non-human sequence.</li>
  </ol>
</div>

<div class="method-box">
  <h3>Non-human contamination screening</h3>
  <p>
    Microbial reads that misalign to the human reference are a documented
    source of false-positive DNMs.  We classify the DKA-supporting reads at
    each candidate site with <a
    href="https://ccb.jhu.edu/software/kraken2/">Kraken2</a> against a
    standard database (human + bacterial + viral + UniVec_Core) and
    annotate every variant with per-class read fractions:
  </p>
  <ul style="margin-left:24px; color:#444; font-size:0.95em;">
    <li><strong>DKA_HLF</strong> &mdash; human-lineage fraction.</li>
    <li><strong>DKA_NHF</strong> &mdash; non-human fraction (driver of stage-4 filter).</li>
    <li><strong>DKA_UCF</strong> &mdash; UniVec Core (vector / adapter) fraction.</li>
    <li><strong>DKA_UF</strong> &mdash; unclassified fraction.</li>
  </ul>
  <p style="margin-top:8px;">
    Variants with <code>DKA_NHF ≥ 0.05</code> are conservatively flagged
    as likely contamination and excluded from the final DNM set.
  </p>
</div>

<!-- ═══════════════════════════════════════════════════════════════════
     Section 3: Results highlights (auto-generated narrative)
     ═══════════════════════════════════════════════════════════════════ -->
{% if stratification and stratification.counts[0] > 0 %}
<h2>3. Results Highlights</h2>
<div class="summary-card">
  <p>
    Of the
    <strong>{{ "{:,}".format(stratification.counts[0]) }}</strong>
    putative DNMs supplied in the input VCF,
    <strong>{{ "{:,}".format(stratification.counts[1]) }}</strong>
    ({{ "%.1f"|format(100 * stratification.counts[1] / stratification.counts[0]) }}%)
    have any k-mer evidence (DKA &gt; 0).  Requiring DKA &ge; 5 narrows to
    <strong>{{ "{:,}".format(stratification.counts[2]) }}</strong>
    ({{ "%.1f"|format(100 * stratification.counts[2] / stratification.counts[0]) }}%
     of input).
    Applying the higher-quality cut-off (DKA_DKT &gt; 0.1) further reduces to
    <strong>{{ "{:,}".format(stratification.counts[3]) }}</strong>
    ({{ "%.1f"|format(100 * stratification.counts[3] / stratification.counts[0]) }}%
     of input).
    Confirming parental absence (MAX_PKC_ALT &lt; 1) yields
    <strong>{{ "{:,}".format(stratification.counts[4]) }}</strong>
    ({{ "%.1f"|format(100 * stratification.counts[4] / stratification.counts[0]) }}%
     of input).
    {% if stratification.has_nhf_data %}
    After excluding likely contamination (NHF &ge; 0.05),
    <strong>{{ "{:,}".format(stratification.counts[5]) }}</strong>
    ({{ "%.1f"|format(100 * stratification.counts[5] / stratification.counts[0]) }}%
     of input) candidates remain &mdash; these constitute the final
    higher-confidence DNM set reported in the table below.
    {% else %}
    <em>No Kraken2 contamination annotations were detected on the input
    VCF, so the contamination filter could not be applied; the
    parental-confirmed set (MAX_PKC_ALT &lt; 1) is reported as the final
    higher-confidence set.</em>
    {% endif %}
    The vast majority of putative DNMs that fail at stage 2 are
    <em>simple miscalled variants</em> &mdash; i.e. positions for which no
    proband read carries a child-unique ALT-supporting k-mer, indicating
    the variant is either inherited (parental support exists) or
    artefactual.
  </p>
</div>
{% endif %}

<!-- ═══════════════════════════════════════════════════════════════════
     Section 4: Variant Type &amp; Chromosomal Breakdown
     ═══════════════════════════════════════════════════════════════════ -->
{% if variant_type_div or chrom_dist_div %}
<h2>4. Variant Breakdown</h2>
<div class="section-rationale">
  <strong>What this shows:</strong> A genomics reviewer will
  immediately check whether the call set has the expected
  composition: ~38 SNVs and ~3–5 indels per WGS trio, distributed
  across all autosomes.  Excess calls on one chromosome or an
  unexpected type ratio are red flags for a systematic artefact.
</div>

{% if variant_type_div %}
<div class="plot-container">
  {{ variant_type_div | safe }}
  <p class="plot-caption">
    SNV / insertion / deletion / MNV counts split by call status.
  </p>
</div>
{% endif %}

{% if chrom_dist_div %}
<div class="plot-container">
  {{ chrom_dist_div | safe }}
  <p class="plot-caption">
    Per-chromosome variant counts.  Clustering on a single chromosome
    can indicate a systematic artefact.
  </p>
</div>
{% endif %}
{% endif %}

<!-- ═══════════════════════════════════════════════════════════════════
     Section 5: K-mer Filtering Funnel (k-mer level)
     ═══════════════════════════════════════════════════════════════════ -->
{% if funnel_div %}
<h2>5. K-mer Filtering Funnel</h2>
<div class="section-rationale">
  <strong>What this shows:</strong> Step-wise attrition of the
  <em>k-mer</em> pool (not the variant pool).  Together with the
  variant funnel above, this answers the reviewer's question
  "how aggressive is the parental subtraction?".
</div>
<div class="plot-container">
  {{ funnel_div | safe }}
</div>
{% endif %}

<!-- ═══════════════════════════════════════════════════════════════════
     Section 6: DKA/DKT Ratio Distribution
     ═══════════════════════════════════════════════════════════════════ -->
{% if variants %}
<h2>6. DKA/DKT Ratio Distribution &amp; Threshold Sensitivity</h2>
<div class="section-rationale">
  <strong>What this shows:</strong> The DKA_DKT ratio captures what
  fraction of spanning fragments carry child-unique allele-supporting
  k-mers &mdash; a direct measure of de novo evidence strength, analogous
  to variant allele fraction (VAF) in somatic callers.  Variants with
  DKA_DKT &gt; 0.1 are higher-quality candidates.
</div>

{% if histogram_div %}
<div class="plot-container">
  {{ histogram_div | safe }}
  <p class="plot-caption">
    Histogram of DKA_DKT ratio.  Variants with no k-mer evidence
    (DKA = 0) are shown in grey; kmer-DNM candidates (DKA &gt; 0)
    in blue.  The dashed red line marks the higher-quality threshold.
  </p>
</div>
{% endif %}

{% if threshold_div %}
<div class="plot-container">
  {{ threshold_div | safe }}
  <p class="plot-caption">
    Sensitivity analysis: number of variants passing for any choice of
    DKA_DKT threshold.  Use this to gauge the impact of moving the cut-off.
  </p>
</div>
{% endif %}

{% if scatter_div %}
<div class="plot-container">
  {{ scatter_div | safe }}
  <p class="plot-caption">
    DKA versus DKT for each variant; marker colour encodes DKA_DKT and
    size encodes DKU.  Genuine DNMs cluster along the diagonal;
    artefacts and inherited variants cluster near the x-axis.
  </p>
</div>
{% endif %}
{% endif %}

<!-- ═══════════════════════════════════════════════════════════════════
     Section 7: Per-Variant Evidence Heatmap
     ═══════════════════════════════════════════════════════════════════ -->
{% if heatmap_div %}
<h2>7. Per-Variant Evidence Heatmap</h2>
<div class="section-rationale">
  <strong>What this shows:</strong> Genomics equivalent of a gene
  expression heatmap &mdash; reveals variant groupings and outlier
  patterns that are invisible in tabular format.  Z-score normalisation
  per column ensures visual comparability across fields with different
  scales.  For very large variant lists this is rendered as a clustered
  summary (one row per k-means cluster) rather than per-variant rows.
</div>
<div class="plot-container">
  {{ heatmap_div | safe }}
</div>
{% endif %}

<!-- ═══════════════════════════════════════════════════════════════════
     Section 8: Parental K-mer Count (PKC) Analysis
     ═══════════════════════════════════════════════════════════════════ -->
{% if pkc_box_div %}
<h2>8. Parental K-mer Count (PKC) Analysis</h2>
<div class="section-rationale">
  <strong>What this shows:</strong> If k-mers are absent from parents,
  is that a coverage gap or true absence?  We use <strong>ALT-allele</strong>
  parental k-mer counts (PKC_ALT).  For genuine de novo variants ALT-allele
  k-mers should be absent from both parents (PKC_ALT ≈ 0); for inherited
  variants PKC_ALT is non-zero because the ALT allele is present in at
  least one parent.  High PKC_ALT at inherited sites confirms the parents
  are well-sequenced; zero PKC_ALT at de novo sites is therefore meaningful
  and not an artifact of parental under-sequencing.
</div>
<div class="plot-container">
  {{ pkc_box_div | safe }}
</div>

{% if pkc_scatter_div %}
<div class="plot-container">
  {{ pkc_scatter_div | safe }}
</div>
{% endif %}
{% endif %}

<!-- ═══════════════════════════════════════════════════════════════════
     Section 9: Non-Human Contamination Profile
     ═══════════════════════════════════════════════════════════════════ -->
{% if variants %}
<h2>9. Non-Human Contamination Profile (Kraken2)</h2>
<div class="section-rationale">
  <strong>Why this matters:</strong> Microbial sequence that
  cross-aligns to the human reference is a well-documented source of
  false-positive de novo calls.  The kmer-denovo-filter pipeline runs
  <a href="https://ccb.jhu.edu/software/kraken2/">Kraken2</a> on the
  DKA-supporting reads at each candidate variant and emits the
  per-class read-fraction tags listed in the methods (DKA_HLF,
  DKA_NHF, DKA_UCF, DKA_UF).  The final contamination filter excludes
  any variant with <code>DKA_NHF &ge; 0.05</code>.  The plots below
  focus on variants that actually exhibit contamination (NHF &ge; 0.05)
  rather than the full set, which is predominantly clean.
</div>

{% if contam_funnel_div %}
<div class="plot-container">
  {{ contam_funnel_div | safe }}
  <p class="plot-caption">
    Contamination prevalence at each stratification stage.  Shows the
    proportion of variants with DKA_NHF &ge; 0.05 at each filter level,
    revealing how contamination is progressively removed and what
    baseline contamination rate exists in the input set.
  </p>
</div>
{% endif %}

{% if nhf_dist_div %}
<div class="plot-container">
  {{ nhf_dist_div | safe }}
  <p class="plot-caption">
    Distribution of the non-human fraction (NHF) among variants with
    putative contamination (NHF &ge; 0.05).  These are the variants
    excluded by the contamination filter.
  </p>
</div>
{% endif %}

{% if contamination_div %}
<div class="plot-container">
  {{ contamination_div | safe }}
  <p class="plot-caption">
    Per-variant Kraken2 read classification breakdown, restricted to
    putative contamination variants (NHF &ge; 0.05) among kmer-DNM
    candidates.
  </p>
</div>
{% endif %}

{% if not stratification or not stratification.has_nhf_data %}
<div class="section-rationale" style="border-left-color:#F58518;">
  <strong>No Kraken2 contamination annotations were detected on the input VCF.</strong>
  The contamination filter has therefore not been applied, and
  any &ldquo;HQ + non-contamination&rdquo; counts above are equivalent to
  the parental-confirmed (MAX_PKC_ALT &lt; 1) counts.  Re-run the pipeline
  with the Kraken2 annotation step enabled (set <code>--kraken2-db</code>
  and related options) to populate this section.
</div>
{% endif %}
{% endif %}

<!-- ═══════════════════════════════════════════════════════════════════
     Section 10: Discovery Mode Results
     ═══════════════════════════════════════════════════════════════════ -->
{% if disc_regions %}
<h2>10. Discovery Mode: Region Landscape</h2>
<div class="section-rationale">
  <strong>What this shows:</strong> The VCF-free discovery mode
  identifies genomic regions enriched for proband-unique k-mers without
  prior variant knowledge.  The distribution and characteristics of
  discovered regions provide independent confirmation of de novo signal.
</div>

{% if disc_scatter_div %}
<div class="plot-container">
  {{ disc_scatter_div | safe }}
</div>
{% endif %}

{% if disc_size_div %}
<div class="plot-container">
  {{ disc_size_div | safe }}
</div>
{% endif %}

{% if sv_evidence_div %}
<div class="plot-container">
  {{ sv_evidence_div | safe }}
</div>
{% endif %}
{% endif %}

<!-- ═══════════════════════════════════════════════════════════════════
     Section 11: Discovery vs VCF Concordance
     ═══════════════════════════════════════════════════════════════════ -->
{% if candidate_comparison and candidate_comparison.get("candidates") %}
<h2>11. Discovery vs. VCF Mode Concordance</h2>
<div class="section-rationale">
  <strong>Scientific rationale:</strong> Cross-validation between two
  independent analytical modes is one of the strongest arguments for
  method validity.  High capture rates of VCF-mode candidates by
  discovery mode demonstrate that the k-mer signal is robust and not
  dependent on prior variant knowledge.
</div>

<div class="metric-cards">
  <div class="metric-card">
    <span class="value">{{ candidate_comparison.get("hq_candidates", 0) }}</span>
    <span class="label">High-Quality Candidates</span>
  </div>
  <div class="metric-card green">
    <span class="value">{{ candidate_comparison.get("captured", 0) }} / {{ candidate_comparison.get("hq_candidates", 0) }}</span>
    <span class="label">Captured by Discovery</span>
  </div>
  <div class="metric-card green">
    <span class="value">{{ "%.1f"|format(100 * candidate_comparison.get("capture_rate", 0)) }}%</span>
    <span class="label">Capture Rate</span>
  </div>
</div>

<table>
  <thead>
    <tr><th>Variant</th><th>DKA</th><th>DKA_DKT</th><th>Captured</th><th>Discovery Region</th></tr>
  </thead>
  <tbody>
    {% for c in candidate_comparison.get("candidates", []) %}
    <tr>
      <td>{{ c.variant }}</td>
      <td>{{ c.dka }}</td>
      <td>{{ "%.4f"|format(c.dka_dkt) }}</td>
      <td><span class="badge {% if c.captured %}badge-green{% else %}badge-red{% endif %}">
        {{ "✓" if c.captured else "✗" }}</span></td>
      <td>{{ c.region or "—" }}</td>
    </tr>
    {% endfor %}
  </tbody>
</table>
{% endif %}

<!-- ═══════════════════════════════════════════════════════════════════
     Section 12: Curated DNM Evaluation
     ═══════════════════════════════════════════════════════════════════ -->
{% if dnm_evaluation and dnm_evaluation.get("loci") %}
<h2>12. Curated DNM Region Evaluation (Sulovari et al. 2023)</h2>
<p class="description">
  Evaluation of discovery regions against curated de novo mutation loci
  from Sulovari et al. 2023.  This gold-standard validation demonstrates
  the method's ability to detect known structural de novo events.
</p>

<div class="metric-cards">
  <div class="metric-card">
    <span class="value">{{ dnm_evaluation.get("total_loci", 0) }}</span>
    <span class="label">Curated Loci</span>
  </div>
  <div class="metric-card green">
    <span class="value">{{ dnm_evaluation.get("detected", 0) }} / {{ dnm_evaluation.get("total_loci", 0) }}</span>
    <span class="label">Detected</span>
  </div>
  <div class="metric-card {% if dnm_evaluation.get("detection_rate", 0) >= 0.7 %}green{% else %}orange{% endif %}">
    <span class="value">{{ "%.1f"|format(100 * dnm_evaluation.get("detection_rate", 0)) }}%</span>
    <span class="label">Detection Rate</span>
  </div>
</div>

<table>
  <thead>
    <tr>
      <th>Locus</th><th>Event</th><th>Size</th><th>Reads</th>
      <th>K-mers</th><th>Signal</th><th>Class</th><th>Status</th>
    </tr>
  </thead>
  <tbody>
    {% for locus in dnm_evaluation.get("loci", []) %}
    <tr>
      <td>{{ locus.locus }}</td>
      <td>{{ locus.event_type }}</td>
      <td>{{ (locus.event_size|string ~ "bp") if locus.event_size else "—" }}</td>
      <td>{{ locus.total_reads }}</td>
      <td>{{ locus.total_unique_kmers }}</td>
      <td>{{ "%.4f"|format(locus.kmer_signal) }}</td>
      <td>{{ locus.sv_class }}</td>
      <td><span class="badge {% if locus.detected %}badge-green{% else %}badge-red{% endif %}">
        {{ locus.assessment }}</span></td>
    </tr>
    {% endfor %}
  </tbody>
</table>
{% endif %}

<!-- ═══════════════════════════════════════════════════════════════════
     Section 13: Per-Variant Detail Table  (final set only)
     ═══════════════════════════════════════════════════════════════════ -->
{% if variants %}
<h2>13. Per-Variant Detail Table — Final DNM Set</h2>
<p class="description">
  This table is restricted to the <strong>final DNM set</strong> only
  (all six filters passed).  The
  full per-variant list (typically tens of thousands of rows on a real
  trio) is not shown here — see <code>summary.txt</code> or the annotated
  VCF for the complete output.
  {% set total_hq = hq_variants | length %}
  {% set shown = variant_table_rows | length %}
  {% if total_hq > shown %}
  Showing {{ shown }} of {{ total_hq }} final-set variants
  (capped at {{ variant_table_max_rows }} rows).
  {% else %}
  {{ shown }} final-set variant{{ "" if shown == 1 else "s" }} shown.
  {% endif %}
</p>
{% if variant_table_rows %}
<table>
  <thead>
    <tr>
      <th>Variant</th><th>Type</th><th>DKU</th><th>DKT</th><th>DKA</th>
      <th>DKU_DKT</th><th>DKA_DKT</th><th>MAX_PKC_ALT</th>
      <th>AVG_PKC_ALT</th>
      {% if stratification and stratification.has_nhf_data %}<th>NHF</th>{% endif %}
      <th>Stage</th>
    </tr>
  </thead>
  <tbody>
    {% for v in variant_table_rows %}
    <tr>
      <td>{{ v.label }}</td>
      <td>{{ v.vtype }}</td>
      <td>{{ v.dku }}</td><td>{{ v.dkt }}</td><td>{{ v.dka }}</td>
      <td>{{ "%.4f"|format(v.dku_dkt) }}</td>
      <td>{{ "%.4f"|format(v.dka_dkt) }}</td>
      <td>{{ v.max_pkc_alt }}</td>
      <td>{{ "%.2f"|format(v.avg_pkc_alt) }}</td>
      {% if stratification and stratification.has_nhf_data %}
      <td>
        {% if "dka_nhf" in v %}{{ "%.3f"|format(v.dka_nhf) }}{% else %}—{% endif %}
      </td>
      {% endif %}
      <td><span class="badge badge-green">Final set</span></td>
    </tr>
    {% endfor %}
  </tbody>
</table>
{% else %}
<p class="description"><em>No variants passed all filters to reach the
  final DNM set (stage 6).</em></p>
{% endif %}
{% endif %}

<!-- ═══════════════════════════════════════════════════════════════════
     Section 14: Method Overview &amp; Interpretation Guide
     ═══════════════════════════════════════════════════════════════════ -->
<h2>14. Method Overview &amp; Interpretation Guide</h2>

<div class="method-box">
  <h3>Why K-mers?</h3>
  <p>
    K-mers provide alignment-independent evidence of novel sequence.
    Unlike read-level methods, k-mer analysis is robust to mapping
    artifacts, reference bias, and complex variants.  By decomposing
    reads into short subsequences and checking parental databases,
    we identify child-unique sequence content that is orthogonal to
    traditional variant calling pipelines.
  </p>
</div>

<div class="method-box">
  <h3>Algorithm Overview</h3>
  <p>
    <strong>1.</strong> Extract k-mers from child reads overlapping each
    candidate variant.<br>
    <strong>2.</strong> Build Jellyfish indexes for each parent and query
    all child k-mers.<br>
    <strong>3.</strong> K-mers absent from both parents are "child-unique"
    and provide evidence of de novo origin.<br>
    <strong>4.</strong> Count spanning fragments (DKT) and fragments with
    unique k-mers supporting the alt allele (DKA) to compute
    DKA_DKT ratio.<br>
    <strong>5.</strong> Classify the DKA-supporting reads with Kraken2 to
    measure non-human / contamination fractions (DKA_NHF, DKA_HLF, ...).<br>
    <strong>6.</strong> Apply the six-stage stratification to produce the
    final higher-confidence DNM set.
  </p>
</div>

<div class="method-box">
  <h3>Key Metrics Interpretation</h3>
  <p>
    <strong>DKU</strong> — Count of distinct child-unique k-mers at a
    variant. DKU &gt; 0 indicates potential de novo origin.<br>
    <strong>DKT</strong> — Total spanning fragments across the variant
    locus.<br>
    <strong>DKA</strong> — Alt-supporting fragments carrying at least
    one child-unique k-mer.<br>
    <strong>DKA_DKT</strong> — Ratio of DKA to DKT; the primary signal
    metric. Values &gt; 0.1 indicate higher-quality
    de novo candidates.<br>
    <strong>PKC_ALT (MAX/AVG/MIN_PKC_ALT)</strong> — ALT-allele Parental K-mer
    Count: how many times the variant's ALT-allele k-mers appear in the
    parents.  For genuine de novo variants these values should be near zero.
    MAX_PKC_ALT &lt; 1 is required for parental confirmation.<br>
    <strong>NHF (DKA_NHF)</strong> — Non-human fraction of informative
    reads. Values &ge; 0.05 trigger the contamination filter.
  </p>
</div>

<footer>
  Generated by <strong>kmer-denovo-filter</strong> v0.1.0 —
  <a href="https://github.com/jlanej/kmer_denovo_filter">GitHub</a>
</footer>

</div><!-- container -->

</body>
</html>
"""


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def generate_report(
    output_path,
    vcf_metrics_path=None,
    vcf_summary_path=None,
    vcf_path=None,
    discovery_metrics_path=None,
    discovery_summary_path=None,
):
    """Generate an interactive HTML report from pipeline output files.

    The report is fully self-contained: Plotly.js is embedded inline so it
    renders correctly in offline/HPC environments and always uses the version
    that matches the installed Python plotly package (fixes CDN version
    mismatch with plotly ≥ 6.x which requires Plotly.js 3.x).

    Parameters
    ----------
    output_path : str
        Path to write the self-contained HTML report.
    vcf_metrics_path : str, optional
        Path to VCF-mode metrics.json.
    vcf_summary_path : str, optional
        Path to VCF-mode summary.txt.
    vcf_path : str, optional
        Path to annotated VCF (for Kraken2 annotations).
    discovery_metrics_path : str, optional
        Path to discovery-mode metrics.json.
    discovery_summary_path : str, optional
        Path to discovery-mode summary.txt.

    Returns
    -------
    str
        The output_path that was written.
    """
    from jinja2 import Template

    context = {
        "mode": None,
        "vcf_metrics": None,
        "disc_metrics": None,
        "variants": [],
        "hq_variants": [],
        "variant_table_rows": [],
        "variant_table_max_rows": _VARIANT_TABLE_MAX_ROWS,
        "high_quality_threshold": _HIGH_QUALITY_DKA_DKT_THRESHOLD,
        "stratification": None,
        "disc_regions": [],
        "candidate_comparison": {},
        "dnm_evaluation": {},
        "plotly_bundle": _get_plotly_bundle(),
        "sankey_div": None,
        "funnel_div": None,
        "strat_funnel_div": None,
        "strat_sankey_div": None,
        "histogram_div": None,
        "threshold_div": None,
        "scatter_div": None,
        "heatmap_div": None,
        "pkc_box_div": None,
        "pkc_scatter_div": None,
        "contamination_div": None,
        "nhf_dist_div": None,
        "contam_funnel_div": None,
        "disc_scatter_div": None,
        "disc_size_div": None,
        "sv_evidence_div": None,
        "variant_type_div": None,
        "chrom_dist_div": None,
    }

    # ── VCF mode data ─────────────────────────────────────────────
    if vcf_metrics_path and os.path.isfile(vcf_metrics_path):
        vcf_metrics = _load_metrics(vcf_metrics_path)
        context["vcf_metrics"] = vcf_metrics
        context["mode"] = "vcf"

        context["funnel_div"] = _make_kmer_funnel_chart(vcf_metrics, mode="vcf")
        context["sankey_div"] = _make_sankey_diagram(vcf_metrics, mode="vcf")

    if vcf_summary_path and os.path.isfile(vcf_summary_path):
        variants = _load_summary_variants(vcf_summary_path)
        # Annotate each variant with its type for the table
        for v in variants:
            v["vtype"] = _classify_variant_type(v["label"])
        context["variants"] = variants

        # ── Kraken2 / contamination annotations (must be merged BEFORE
        # stratification so the NHF contamination filter can be applied) ───────
        kraken2_data = []
        if vcf_path and os.path.isfile(vcf_path):
            kraken2_data = _load_vcf_kraken2_annotations(vcf_path)
            if kraken2_data:
                _merge_kraken2_into_variants(variants, kraken2_data)

        if variants:
            # Compute stratification for both summary and downstream plots
            stratification = _compute_stratification(variants)
            context["stratification"] = stratification

            context["strat_funnel_div"] = _make_stratification_funnel(
                stratification,
            )
            context["strat_sankey_div"] = _make_stratification_sankey(
                stratification,
            )

            context["histogram_div"] = _make_dka_dkt_histogram(variants)
            context["threshold_div"] = _make_threshold_sensitivity(variants)
            context["scatter_div"] = _make_dka_vs_dkt_scatter(variants)
            context["heatmap_div"] = _make_evidence_heatmap(variants)
            context["pkc_box_div"] = _make_pkc_boxplot(variants)
            context["pkc_scatter_div"] = _make_pkc_vs_dka_dkt_scatter(variants)
            context["variant_type_div"] = _make_variant_type_breakdown(variants)
            context["chrom_dist_div"] = _make_chromosomal_distribution(variants)

            # Per-variant detail table: restrict to the final DNM set
            # (stage 5 = all filters passed).  The full per-variant list
            # is often 30 000+ rows on a real trio and is not informative
            # for a reviewer.
            hq_variants = [v for v in variants if v["stage"] >= 5]
            hq_variants.sort(
                key=lambda v: (-(v.get("stage", 0)), -v.get("dka_dkt", 0.0)),
            )
            context["hq_variants"] = hq_variants
            context["variant_table_rows"] = hq_variants[:_VARIANT_TABLE_MAX_ROWS]

            # Contamination plots (require Kraken2 data)
            if kraken2_data:
                context["contamination_div"] = _make_contamination_bar(
                    variants, kraken2_data,
                )
                context["nhf_dist_div"] = _make_nhf_distribution_plot(variants)
                context["contam_funnel_div"] = _make_contamination_funnel(
                    stratification, variants,
                )

    # ── Discovery mode data ───────────────────────────────────────
    if discovery_metrics_path and os.path.isfile(discovery_metrics_path):
        disc_metrics = _load_metrics(discovery_metrics_path)
        context["disc_metrics"] = disc_metrics
        if context["mode"] is None:
            context["mode"] = "discovery"
        else:
            context["mode"] = "combined"

        regions = disc_metrics.get("regions", [])
        context["disc_regions"] = regions
        context["candidate_comparison"] = disc_metrics.get(
            "candidate_comparison", {},
        )
        context["dnm_evaluation"] = disc_metrics.get("dnm_evaluation", {})

        if not context.get("funnel_div"):
            context["funnel_div"] = _make_kmer_funnel_chart(
                disc_metrics, mode="discovery",
            )
        if not context.get("sankey_div"):
            context["sankey_div"] = _make_sankey_diagram(
                disc_metrics, mode="discovery",
            )

        if regions:
            context["disc_scatter_div"] = _make_discovery_region_scatter(
                regions,
            )
            context["disc_size_div"] = _make_discovery_size_histogram(regions)
            context["sv_evidence_div"] = _make_sv_evidence_chart(regions)

    # ── Render ────────────────────────────────────────────────────
    template = Template(_HTML_TEMPLATE)
    html = template.render(**context)

    os.makedirs(os.path.dirname(os.path.abspath(output_path)) or ".", exist_ok=True)
    with open(output_path, "w") as fh:
        fh.write(html)

    logger.info("Report written to: %s", output_path)
    return output_path
