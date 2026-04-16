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

logger = logging.getLogger(__name__)


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
# Plot generators (return Plotly JSON strings)
# ---------------------------------------------------------------------------

def _plotly_json(fig):
    """Serialize a Plotly figure to a JSON string for embedding."""
    return fig.to_json()


def _make_kmer_funnel_chart(metrics, mode="vcf"):
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

    return _plotly_json(fig)


def _make_sankey_diagram(metrics, mode="vcf"):
    """Create a Sankey diagram showing the filtering flow."""
    import plotly.graph_objects as go

    if mode == "vcf":
        total = metrics.get("total_child_kmers", 0)
        parent_found = metrics.get("parent_found_kmers", 0)
        unique = metrics.get("child_unique_kmers", 0)
        with_reads = metrics.get("variants_with_unique_reads", 0)
        total_variants = metrics.get("total_variants", 0)
        without_reads = total_variants - with_reads

        node_labels = [
            f"Total Child K-mers ({total:,})",
            f"Found in Parents ({parent_found:,})",
            f"Child-Unique ({unique:,})",
            f"Variants with Unique Reads ({with_reads})",
            f"Variants without Unique Reads ({without_reads})",
        ]
        node_colors = ["#4C78A8", "#E45756", "#54A24B", "#72B7B2", "#BAB0AC"]

        source = [0, 0, 2, 2]
        target = [1, 2, 3, 4]
        value = [parent_found, unique,
                 max(1, with_reads), max(1, without_reads)]
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
    return _plotly_json(fig)


def _make_dka_dkt_histogram(variants):
    """Create a histogram of DKA_DKT ratios with threshold marker."""
    import plotly.graph_objects as go

    dka_dkt_values = [v["dka_dkt"] for v in variants]

    fig = go.Figure()
    fig.add_trace(go.Histogram(
        x=dka_dkt_values,
        nbinsx=30,
        marker_color="#4C78A8",
        opacity=0.85,
        hovertemplate="DKA_DKT: %{x:.3f}<br>Count: %{y}<extra></extra>",
    ))
    # Threshold line at 0.25
    fig.add_vline(
        x=0.25, line_dash="dash", line_color="#E45756", line_width=2,
        annotation_text="High-quality threshold (0.25)",
        annotation_position="top right",
        annotation_font=dict(size=11, color="#E45756"),
    )
    fig.update_layout(
        title=dict(
            text="DKA/DKT Ratio Distribution",
            font=dict(size=18),
        ),
        xaxis_title="DKA_DKT Ratio",
        yaxis_title="Number of Variants",
        template="plotly_white",
        height=400,
        margin=dict(t=60, b=40),
    )
    return _plotly_json(fig)


def _make_dka_vs_dkt_scatter(variants):
    """Create a scatter plot of DKA vs DKT colored by DKA_DKT ratio."""
    import plotly.graph_objects as go

    fig = go.Figure()
    fig.add_trace(go.Scatter(
        x=[v["dkt"] for v in variants],
        y=[v["dka"] for v in variants],
        mode="markers",
        marker=dict(
            size=[max(6, min(30, v["dku"] * 3)) for v in variants],
            color=[v["dka_dkt"] for v in variants],
            colorscale="Viridis",
            colorbar=dict(title="DKA_DKT"),
            showscale=True,
            line=dict(width=1, color="#333"),
        ),
        text=[v["label"] for v in variants],
        hovertemplate=(
            "<b>%{text}</b><br>"
            "DKT: %{x}<br>DKA: %{y}<br>"
            "DKU: %{customdata[0]}<br>"
            "DKA_DKT: %{customdata[1]:.4f}<br>"
            "Call: %{customdata[2]}"
            "<extra></extra>"
        ),
        customdata=[[v["dku"], v["dka_dkt"], v["call"]] for v in variants],
    ))

    # Add quadrant lines
    fig.add_hline(y=10, line_dash="dot", line_color="#ccc", line_width=1)
    fig.add_vline(x=0, line_dash="dot", line_color="#ccc", line_width=1)

    fig.update_layout(
        title=dict(
            text="DKA vs. DKT (size = DKU, color = DKA_DKT ratio)",
            font=dict(size=18),
        ),
        xaxis_title="DKT (Total Spanning Fragments)",
        yaxis_title="DKA (Alt-Supporting Fragments with Unique K-mers)",
        template="plotly_white",
        height=500,
        margin=dict(t=60, b=40),
    )
    return _plotly_json(fig)


def _make_evidence_heatmap(variants):
    """Create a clustered heatmap of per-variant evidence fields."""
    import plotly.graph_objects as go

    fields = [
        "dku", "dkt", "dka", "dku_dkt", "dka_dkt",
        "max_pkc", "avg_pkc", "min_pkc",
    ]
    display_fields = [
        "DKU", "DKT", "DKA", "DKU_DKT", "DKA_DKT",
        "MAX_PKC", "AVG_PKC", "MIN_PKC",
    ]
    labels = [v["label"] for v in variants]

    # Build raw data matrix
    raw = []
    for v in variants:
        raw.append([v[f] for f in fields])

    # Z-score normalization per column for visual comparability
    import statistics as stats
    n_cols = len(fields)
    n_rows = len(variants)
    z_data = []
    for r in range(n_rows):
        z_data.append([0.0] * n_cols)

    for c in range(n_cols):
        col_vals = [raw[r][c] for r in range(n_rows)]
        mean_val = stats.mean(col_vals) if col_vals else 0
        std_val = stats.pstdev(col_vals) if col_vals else 1
        if std_val == 0:
            std_val = 1
        for r in range(n_rows):
            z_data[r][c] = (raw[r][c] - mean_val) / std_val

    # Build hover text with raw values
    hover_text = []
    for r in range(n_rows):
        row_hover = []
        for c in range(n_cols):
            row_hover.append(
                f"{labels[r]}<br>{display_fields[c]}: {raw[r][c]}"
                f"<br>Z-score: {z_data[r][c]:.2f}"
            )
        hover_text.append(row_hover)

    fig = go.Figure(data=go.Heatmap(
        z=z_data,
        x=display_fields,
        y=labels,
        colorscale="RdBu_r",
        zmid=0,
        text=hover_text,
        hoverinfo="text",
        colorbar=dict(title="Z-score"),
    ))
    fig.update_layout(
        title=dict(
            text="Per-Variant Evidence Heatmap (Z-score normalized)",
            font=dict(size=18),
        ),
        template="plotly_white",
        height=max(400, 30 * len(variants) + 100),
        margin=dict(t=60, b=40, l=250),
        yaxis=dict(autorange="reversed"),
    )
    return _plotly_json(fig)


def _make_pkc_boxplot(variants):
    """Create box plots of PKC metrics by call type."""
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
            ("max_pkc", "MAX_PKC"),
            ("avg_pkc", "AVG_PKC"),
            ("min_pkc", "MIN_PKC"),
        ]:
            fig.add_trace(go.Box(
                y=[v[metric] for v in group],
                name=f"{name}<br>({label_group})",
                marker_color=color,
                boxmean=True,
            ))

    fig.update_layout(
        title=dict(
            text="Parental K-mer Count (PKC) by Call Type",
            font=dict(size=18),
        ),
        yaxis_title="K-mer Count in Parents",
        template="plotly_white",
        height=450,
        margin=dict(t=60, b=40),
        showlegend=False,
    )
    return _plotly_json(fig)


def _make_pkc_vs_dka_dkt_scatter(variants):
    """Create AVG_PKC vs DKA_DKT scatter plot."""
    import plotly.graph_objects as go

    colors = ["#54A24B" if v["call"] == "DE_NOVO" else "#E45756"
              for v in variants]

    fig = go.Figure()
    fig.add_trace(go.Scatter(
        x=[v["dka_dkt"] for v in variants],
        y=[v["avg_pkc"] for v in variants],
        mode="markers",
        marker=dict(size=10, color=colors, line=dict(width=1, color="#333")),
        text=[f"{v['label']}<br>Call: {v['call']}" for v in variants],
        hovertemplate=(
            "<b>%{text}</b><br>"
            "DKA_DKT: %{x:.4f}<br>AVG_PKC: %{y:.1f}"
            "<extra></extra>"
        ),
    ))
    fig.update_layout(
        title=dict(
            text="AVG_PKC vs. DKA_DKT Ratio",
            font=dict(size=18),
        ),
        xaxis_title="DKA_DKT Ratio",
        yaxis_title="AVG_PKC (Average Parental K-mer Count)",
        template="plotly_white",
        height=450,
        margin=dict(t=60, b=40),
    )
    return _plotly_json(fig)


def _make_contamination_bar(variants, kraken2_data):
    """Create stacked bar chart of Kraken2 read classification fractions."""
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
        k = kraken_map.get(v["label"])
        if k and any(k.get(f, 0) > 0 for f in
                      ["DKA_HLF", "DKA_NHF", "DKA_UCF", "DKA_UF"]):
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
            text="Kraken2 Read Classification per Variant",
            font=dict(size=18),
        ),
        yaxis_title="Fraction of DKA Reads",
        template="plotly_white",
        height=450,
        margin=dict(t=60, b=120),
        xaxis_tickangle=-45,
        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1),
    )
    return _plotly_json(fig)


def _make_discovery_region_scatter(regions):
    """Create scatter plot of discovery regions: reads vs k-mers."""
    import plotly.graph_objects as go

    class_colors = {"SMALL": "#4C78A8", "AMBIGUOUS": "#F58518", "SV": "#E45756"}

    fig = go.Figure()
    for cls in ["SMALL", "AMBIGUOUS", "SV"]:
        cls_regions = [r for r in regions if r.get("class") == cls]
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

    fig.update_layout(
        title=dict(
            text="Discovery Regions: Reads vs. Distinct K-mers",
            font=dict(size=18),
        ),
        xaxis_title="Supporting Reads",
        yaxis_title="Distinct Proband-Unique K-mers",
        template="plotly_white",
        height=450,
        margin=dict(t=60, b=40),
    )
    return _plotly_json(fig)


def _make_discovery_size_histogram(regions):
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
    return _plotly_json(fig)


def _make_sv_evidence_chart(regions):
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
    return _plotly_json(fig)


def _make_threshold_sensitivity(variants):
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
        x=0.25, line_dash="dash", line_color="#E45756", line_width=2,
        annotation_text="0.25",
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
    return _plotly_json(fig)


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
<script src="https://cdn.plot.ly/plotly-2.35.2.min.js" charset="utf-8"></script>
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
     Section 1: Executive Summary
     ═══════════════════════════════════════════════════════════════════ -->
<h2>1. Executive Summary</h2>
<p class="description">
  This dashboard summarizes the k-mer-based de novo variant filtering
  results. K-mers (short DNA subsequences of length <em>k</em>) from
  the child are compared against both parents.  K-mers present in the
  child but absent from both parents represent candidate de novo mutations.
</p>

<div class="metric-cards">
  {% if vcf_metrics %}
  <div class="metric-card">
    <span class="value">{{ "{:,}".format(vcf_metrics.get("total_variants", 0)) }}</span>
    <span class="label">Total Candidates</span>
  </div>
  <div class="metric-card green">
    <span class="value">{{ vcf_metrics.get("variants_with_unique_reads", 0) }}</span>
    <span class="label">Likely De Novo (DKU &gt; 0)</span>
  </div>
  <div class="metric-card">
    <span class="value">{{ "{:,}".format(vcf_metrics.get("total_child_kmers", 0)) }}</span>
    <span class="label">Total Child K-mers</span>
  </div>
  <div class="metric-card green">
    <span class="value">{{ "{:,}".format(vcf_metrics.get("child_unique_kmers", 0)) }}</span>
    <span class="label">Child-Unique K-mers</span>
  </div>
  <div class="metric-card">
    <span class="value">
      {% if vcf_metrics.get("total_child_kmers", 0) > 0 %}
        {{ "%.1f"|format(100 * vcf_metrics.get("child_unique_kmers", 0) / vcf_metrics.get("total_child_kmers", 1)) }}%
      {% else %}0%{% endif %}
    </span>
    <span class="label">Unique Rate</span>
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

{% if sankey_json %}
<div class="plot-container">
  <div id="sankey-plot" class="plot-div"></div>
</div>
{% endif %}

<!-- ═══════════════════════════════════════════════════════════════════
     Section 2: K-mer Filtering Funnel
     ═══════════════════════════════════════════════════════════════════ -->
{% if funnel_json %}
<h2>2. K-mer Filtering Funnel</h2>
<div class="section-rationale">
  <strong>Scientific rationale:</strong> A skeptic's first question is
  "how aggressive is this filter?"  The step-wise attrition plot shows
  the method is not over- or under-filtering.  This mirrors the filtering
  cascade standard in variant calling benchmarks (e.g., GIAB Mendelian
  violation rates).
</div>
<div class="plot-container">
  <div id="funnel-plot" class="plot-div"></div>
</div>
{% endif %}

<!-- ═══════════════════════════════════════════════════════════════════
     Section 3: DKA/DKT Ratio Distribution
     ═══════════════════════════════════════════════════════════════════ -->
{% if variants %}
<h2>3. DKA/DKT Ratio Distribution &amp; Threshold Sensitivity</h2>
<div class="section-rationale">
  <strong>Scientific rationale:</strong> The DKA_DKT ratio captures what
  fraction of spanning fragments carry child-unique allele-supporting
  k-mers &mdash; a direct measure of de novo evidence strength, analogous
  to variant allele fraction (VAF) in somatic callers.  Variants with
  DKA_DKT &ge; 0.25 and DKA &ge; 10 are considered high-quality
  candidates.
</div>

{% if histogram_json %}
<div class="plot-container">
  <div id="histogram-plot" class="plot-div"></div>
</div>
{% endif %}

{% if threshold_json %}
<div class="plot-container">
  <div id="threshold-plot" class="plot-div"></div>
</div>
{% endif %}

{% if scatter_json %}
<div class="plot-container">
  <div id="scatter-plot" class="plot-div"></div>
</div>
{% endif %}
{% endif %}

<!-- ═══════════════════════════════════════════════════════════════════
     Section 4: Per-Variant Evidence Heatmap
     ═══════════════════════════════════════════════════════════════════ -->
{% if heatmap_json %}
<h2>4. Per-Variant Evidence Heatmap</h2>
<div class="section-rationale">
  <strong>Scientific rationale:</strong> This is the genomics equivalent
  of a gene expression heatmap — it reveals variant groupings and outlier
  patterns that are invisible in tabular format.  Z-score normalization
  per column ensures visual comparability across fields with different
  scales.
</div>
<div class="plot-container">
  <div id="heatmap-plot" class="plot-div"></div>
</div>
{% endif %}

<!-- ═══════════════════════════════════════════════════════════════════
     Section 5: Parental K-mer Count (PKC) Analysis
     ═══════════════════════════════════════════════════════════════════ -->
{% if pkc_box_json %}
<h2>5. Parental K-mer Count (PKC) Analysis</h2>
<div class="section-rationale">
  <strong>Scientific rationale:</strong> If k-mers are absent from parents,
  is that a coverage gap or true absence?  High parental k-mer counts at
  inherited sites confirm the parents are well-sequenced; absence at
  de novo sites is therefore meaningful.  This directly addresses the
  null hypothesis that child-unique k-mers are artifacts of parental
  under-sequencing.
</div>
<div class="plot-container">
  <div id="pkc-box-plot" class="plot-div"></div>
</div>

{% if pkc_scatter_json %}
<div class="plot-container">
  <div id="pkc-scatter-plot" class="plot-div"></div>
</div>
{% endif %}
{% endif %}

<!-- ═══════════════════════════════════════════════════════════════════
     Section 6: Non-Human Contamination Profile
     ═══════════════════════════════════════════════════════════════════ -->
{% if contamination_json %}
<h2>6. Non-Human Contamination Profile (Kraken2)</h2>
<div class="section-rationale">
  <strong>Scientific rationale:</strong> Microbial contamination is a
  known source of false-positive de novo calls.  The Kraken2
  classification provides a taxonomic audit trail for every informative
  read, addressing the criticism that novel k-mers could come from
  contamination rather than true mutation.
</div>
<div class="plot-container">
  <div id="contamination-plot" class="plot-div"></div>
</div>
{% endif %}

<!-- ═══════════════════════════════════════════════════════════════════
     Section 7: Discovery Mode Results
     ═══════════════════════════════════════════════════════════════════ -->
{% if disc_regions %}
<h2>7. Discovery Mode: Region Landscape</h2>
<div class="section-rationale">
  <strong>Scientific rationale:</strong> The VCF-free discovery mode
  identifies genomic regions enriched for proband-unique k-mers without
  prior variant knowledge.  The distribution and characteristics of
  discovered regions provide independent confirmation of de novo signal.
</div>

{% if disc_scatter_json %}
<div class="plot-container">
  <div id="disc-scatter-plot" class="plot-div"></div>
</div>
{% endif %}

{% if disc_size_json %}
<div class="plot-container">
  <div id="disc-size-plot" class="plot-div"></div>
</div>
{% endif %}

{% if sv_evidence_json %}
<div class="plot-container">
  <div id="sv-evidence-plot" class="plot-div"></div>
</div>
{% endif %}
{% endif %}

<!-- ═══════════════════════════════════════════════════════════════════
     Section 8: Discovery vs VCF Concordance
     ═══════════════════════════════════════════════════════════════════ -->
{% if candidate_comparison and candidate_comparison.get("candidates") %}
<h2>8. Discovery vs. VCF Mode Concordance</h2>
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
     Section 9: Curated DNM Evaluation
     ═══════════════════════════════════════════════════════════════════ -->
{% if dnm_evaluation and dnm_evaluation.get("loci") %}
<h2>9. Curated DNM Region Evaluation (Sulovari et al. 2023)</h2>
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
     Section 10: Per-Variant Detail Table
     ═══════════════════════════════════════════════════════════════════ -->
{% if variants %}
<h2>10. Per-Variant Detail Table</h2>
<table>
  <thead>
    <tr>
      <th>Variant</th><th>DKU</th><th>DKT</th><th>DKA</th>
      <th>DKU_DKT</th><th>DKA_DKT</th><th>MAX_PKC</th>
      <th>AVG_PKC</th><th>MIN_PKC</th><th>Call</th>
    </tr>
  </thead>
  <tbody>
    {% for v in variants %}
    <tr>
      <td>{{ v.label }}</td>
      <td>{{ v.dku }}</td><td>{{ v.dkt }}</td><td>{{ v.dka }}</td>
      <td>{{ "%.4f"|format(v.dku_dkt) }}</td>
      <td>{{ "%.4f"|format(v.dka_dkt) }}</td>
      <td>{{ v.max_pkc }}</td>
      <td>{{ "%.2f"|format(v.avg_pkc) }}</td>
      <td>{{ v.min_pkc }}</td>
      <td><span class="badge {% if v.call == 'DE_NOVO' %}badge-green{% else %}badge-orange{% endif %}">
        {{ v.call }}</span></td>
    </tr>
    {% endfor %}
  </tbody>
</table>
{% endif %}

<!-- ═══════════════════════════════════════════════════════════════════
     Section 11: Method Overview
     ═══════════════════════════════════════════════════════════════════ -->
<h2>11. Method Overview &amp; Interpretation Guide</h2>

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
    DKA_DKT ratio.
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
    metric. Values &ge; 0.25 with DKA &ge; 10 indicate high-quality
    de novo candidates.<br>
    <strong>PKC (MAX/AVG/MIN)</strong> — Parental K-mer Count: how many
    times the variant's k-mers appear in parents. High PKC at inherited
    sites confirms parental sequencing depth is adequate.<br>
    <strong>NHF (DKA_NHF)</strong> — Non-human fraction of informative
    reads. Values &gt; 0.1 warrant caution as the signal may be driven
    by microbial contamination.
  </p>
</div>

<footer>
  Generated by <strong>kmer-denovo-filter</strong> v0.1.0 —
  <a href="https://github.com/jlanej/kmer_denovo_filter">GitHub</a>
</footer>

</div><!-- container -->

<!-- ═══════════════════════════════════════════════════════════════════
     Plotly rendering
     ═══════════════════════════════════════════════════════════════════ -->
<script>
{% if sankey_json %}
(function() {
  var fig = JSON.parse({{ sankey_json | tojson }});
  Plotly.newPlot("sankey-plot", fig.data, fig.layout, {responsive: true});
})();
{% endif %}

{% if funnel_json %}
(function() {
  var fig = JSON.parse({{ funnel_json | tojson }});
  Plotly.newPlot("funnel-plot", fig.data, fig.layout, {responsive: true});
})();
{% endif %}

{% if histogram_json %}
(function() {
  var fig = JSON.parse({{ histogram_json | tojson }});
  Plotly.newPlot("histogram-plot", fig.data, fig.layout, {responsive: true});
})();
{% endif %}

{% if threshold_json %}
(function() {
  var fig = JSON.parse({{ threshold_json | tojson }});
  Plotly.newPlot("threshold-plot", fig.data, fig.layout, {responsive: true});
})();
{% endif %}

{% if scatter_json %}
(function() {
  var fig = JSON.parse({{ scatter_json | tojson }});
  Plotly.newPlot("scatter-plot", fig.data, fig.layout, {responsive: true});
})();
{% endif %}

{% if heatmap_json %}
(function() {
  var fig = JSON.parse({{ heatmap_json | tojson }});
  Plotly.newPlot("heatmap-plot", fig.data, fig.layout, {responsive: true});
})();
{% endif %}

{% if pkc_box_json %}
(function() {
  var fig = JSON.parse({{ pkc_box_json | tojson }});
  Plotly.newPlot("pkc-box-plot", fig.data, fig.layout, {responsive: true});
})();
{% endif %}

{% if pkc_scatter_json %}
(function() {
  var fig = JSON.parse({{ pkc_scatter_json | tojson }});
  Plotly.newPlot("pkc-scatter-plot", fig.data, fig.layout, {responsive: true});
})();
{% endif %}

{% if contamination_json %}
(function() {
  var fig = JSON.parse({{ contamination_json | tojson }});
  Plotly.newPlot("contamination-plot", fig.data, fig.layout, {responsive: true});
})();
{% endif %}

{% if disc_scatter_json %}
(function() {
  var fig = JSON.parse({{ disc_scatter_json | tojson }});
  Plotly.newPlot("disc-scatter-plot", fig.data, fig.layout, {responsive: true});
})();
{% endif %}

{% if disc_size_json %}
(function() {
  var fig = JSON.parse({{ disc_size_json | tojson }});
  Plotly.newPlot("disc-size-plot", fig.data, fig.layout, {responsive: true});
})();
{% endif %}

{% if sv_evidence_json %}
(function() {
  var fig = JSON.parse({{ sv_evidence_json | tojson }});
  Plotly.newPlot("sv-evidence-plot", fig.data, fig.layout, {responsive: true});
})();
{% endif %}
</script>

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
        "disc_regions": [],
        "candidate_comparison": {},
        "dnm_evaluation": {},
        "sankey_json": None,
        "funnel_json": None,
        "histogram_json": None,
        "threshold_json": None,
        "scatter_json": None,
        "heatmap_json": None,
        "pkc_box_json": None,
        "pkc_scatter_json": None,
        "contamination_json": None,
        "disc_scatter_json": None,
        "disc_size_json": None,
        "sv_evidence_json": None,
    }

    # ── VCF mode data ─────────────────────────────────────────────
    if vcf_metrics_path and os.path.isfile(vcf_metrics_path):
        vcf_metrics = _load_metrics(vcf_metrics_path)
        context["vcf_metrics"] = vcf_metrics
        context["mode"] = "vcf"

        context["funnel_json"] = _make_kmer_funnel_chart(vcf_metrics, mode="vcf")
        context["sankey_json"] = _make_sankey_diagram(vcf_metrics, mode="vcf")

    if vcf_summary_path and os.path.isfile(vcf_summary_path):
        variants = _load_summary_variants(vcf_summary_path)
        context["variants"] = variants

        if variants:
            context["histogram_json"] = _make_dka_dkt_histogram(variants)
            context["threshold_json"] = _make_threshold_sensitivity(variants)
            context["scatter_json"] = _make_dka_vs_dkt_scatter(variants)
            context["heatmap_json"] = _make_evidence_heatmap(variants)
            context["pkc_box_json"] = _make_pkc_boxplot(variants)
            context["pkc_scatter_json"] = _make_pkc_vs_dka_dkt_scatter(variants)

    # Kraken2 annotations from VCF
    if vcf_path and os.path.isfile(vcf_path):
        kraken2_data = _load_vcf_kraken2_annotations(vcf_path)
        if kraken2_data and context["variants"]:
            contamination_fig = _make_contamination_bar(
                context["variants"], kraken2_data,
            )
            context["contamination_json"] = contamination_fig

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

        if not context.get("funnel_json"):
            context["funnel_json"] = _make_kmer_funnel_chart(
                disc_metrics, mode="discovery",
            )
        if not context.get("sankey_json"):
            context["sankey_json"] = _make_sankey_diagram(
                disc_metrics, mode="discovery",
            )

        if regions:
            context["disc_scatter_json"] = _make_discovery_region_scatter(
                regions,
            )
            context["disc_size_json"] = _make_discovery_size_histogram(regions)
            context["sv_evidence_json"] = _make_sv_evidence_chart(regions)

    # ── Render ────────────────────────────────────────────────────
    template = Template(_HTML_TEMPLATE)
    html = template.render(**context)

    os.makedirs(os.path.dirname(os.path.abspath(output_path)) or ".", exist_ok=True)
    with open(output_path, "w") as fh:
        fh.write(html)

    logger.info("Report written to: %s", output_path)
    return output_path
