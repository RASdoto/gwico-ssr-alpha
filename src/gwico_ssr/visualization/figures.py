"""Visualization builders for GWICO-SSR.

Generates publication-ready static figures (PNG/SVG via matplotlib/seaborn)
and interactive HTML figures (via plotly) from database-backed pipeline state.

All functions accept pre-queried data and an output directory,
keeping figure logic decoupled from DB access.
"""

from __future__ import annotations

import os
from pathlib import Path
from typing import Sequence

import matplotlib
matplotlib.use("Agg")  # Non-interactive backend for headless rendering

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------

_MOTIF_SIZE_ORDER = ["mono", "di", "tri", "tetra", "penta", "hexa"]
_MOTIF_SIZE_LABELS = {
    "mono": "Mono (1bp)", "di": "Di (2bp)", "tri": "Tri (3bp)",
    "tetra": "Tetra (4bp)", "penta": "Penta (5bp)", "hexa": "Hexa (6bp)",
}


def _ensure_dir(path: str | Path) -> Path:
    p = Path(path)
    p.mkdir(parents=True, exist_ok=True)
    return p


# ---------------------------------------------------------------------------
# 1. Motif size distribution bar chart (static)
# ---------------------------------------------------------------------------

def plot_motif_size_distribution(
    motif_size_data: list[dict],
    output_dir: str | Path,
    fmt: str = "png",
) -> Path:
    """Bar chart of SSR counts by motif size class.

    *motif_size_data*: list of dicts with keys ``motif_size`` (str) and
    ``total_ssrs`` (int), as returned by ``aggregate_by_motif_size``.
    """
    out = _ensure_dir(output_dir)
    # Preserve canonical order
    ordered = {d["motif_size"]: d["total_ssrs"] for d in motif_size_data}
    labels = [s for s in _MOTIF_SIZE_ORDER if s in ordered]
    values = [ordered[s] for s in labels]
    display_labels = [_MOTIF_SIZE_LABELS.get(s, s) for s in labels]

    fig, ax = plt.subplots(figsize=(8, 5))
    bars = ax.bar(display_labels, values, color=sns.color_palette("viridis", len(labels)))
    ax.set_xlabel("Motif Size Class")
    ax.set_ylabel("Total SSR Count")
    ax.set_title("SSR Distribution by Motif Size")
    ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
    for bar, val in zip(bars, values):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height(),
                f"{val:,}", ha="center", va="bottom", fontsize=8)
    fig.tight_layout()
    path = out / f"motif_size_distribution.{fmt}"
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path


# ---------------------------------------------------------------------------
# 2. Top motif frequency bar chart (static)
# ---------------------------------------------------------------------------

def plot_top_motifs(
    motif_freq_data: list[dict],
    output_dir: str | Path,
    top_n: int = 20,
    fmt: str = "png",
) -> Path:
    """Horizontal bar chart of the most frequent canonical motifs.

    *motif_freq_data*: list of dicts with ``motif_canonical``, ``total_count``.
    """
    out = _ensure_dir(output_dir)
    sorted_data = sorted(motif_freq_data, key=lambda d: d["total_count"], reverse=True)[:top_n]
    motifs = [d["motif_canonical"] for d in sorted_data][::-1]
    counts = [d["total_count"] for d in sorted_data][::-1]

    fig, ax = plt.subplots(figsize=(8, max(4, len(motifs) * 0.35)))
    ax.barh(motifs, counts, color=sns.color_palette("mako", len(motifs)))
    ax.set_xlabel("Total Occurrences")
    ax.set_title(f"Top {min(top_n, len(motifs))} Most Frequent Motifs")
    ax.ticklabel_format(axis="x", style="sci", scilimits=(0, 0))
    fig.tight_layout()
    path = out / f"top_motifs.{fmt}"
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path


# ---------------------------------------------------------------------------
# 3. Gene distribution bar chart (static)
# ---------------------------------------------------------------------------

def plot_gene_distribution(
    gene_data: list[dict],
    output_dir: str | Path,
    fmt: str = "png",
) -> Path:
    """Bar chart of SSR counts per gene.

    *gene_data*: list of dicts with ``gene_name`` and ``ssr_count``.
    """
    out = _ensure_dir(output_dir)
    sorted_data = sorted(gene_data, key=lambda d: d["ssr_count"], reverse=True)
    genes = [d["gene_name"] for d in sorted_data]
    counts = [d["ssr_count"] for d in sorted_data]

    fig, ax = plt.subplots(figsize=(max(6, len(genes) * 0.6), 5))
    bars = ax.bar(genes, counts, color=sns.color_palette("Set2", len(genes)))
    ax.set_xlabel("Gene")
    ax.set_ylabel("SSR Count")
    ax.set_title("SSR Distribution by Gene")
    ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
    plt.xticks(rotation=45, ha="right")
    fig.tight_layout()
    path = out / f"gene_distribution.{fmt}"
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path


# ---------------------------------------------------------------------------
# 4. Country SSR burden boxplot (static)
# ---------------------------------------------------------------------------

def plot_country_boxplot(
    country_metrics: dict[str, list[float]],
    output_dir: str | Path,
    metric_label: str = "SSR Count",
    fmt: str = "png",
) -> Path:
    """Boxplot of SSR metric distribution across top countries.

    *country_metrics*: dict mapping country name -> list of per-accession values.
    """
    out = _ensure_dir(output_dir)
    # Sort by median descending, take top 20
    sorted_countries = sorted(
        country_metrics.items(), key=lambda kv: np.median(kv[1]), reverse=True
    )[:20]

    fig, ax = plt.subplots(figsize=(max(8, len(sorted_countries) * 0.5), 6))
    data = [vals for _, vals in sorted_countries]
    labels = [name for name, _ in sorted_countries]
    bp = ax.boxplot(data, tick_labels=labels, patch_artist=True, showfliers=False)
    palette = sns.color_palette("coolwarm", len(data))
    for patch, color in zip(bp["boxes"], palette):
        patch.set_facecolor(color)
    ax.set_xlabel("Country")
    ax.set_ylabel(metric_label)
    ax.set_title(f"{metric_label} by Country")
    plt.xticks(rotation=45, ha="right")
    fig.tight_layout()
    path = out / f"country_boxplot.{fmt}"
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path


# ---------------------------------------------------------------------------
# 5. Correlation scatter (static)
# ---------------------------------------------------------------------------

def plot_correlation_scatter(
    x: list[float],
    y: list[float],
    x_label: str,
    y_label: str,
    output_dir: str | Path,
    filename: str = "correlation_scatter",
    fmt: str = "png",
) -> Path:
    """Scatter plot with regression line for two numeric arrays."""
    out = _ensure_dir(output_dir)
    fig, ax = plt.subplots(figsize=(7, 6))
    ax.scatter(x, y, alpha=0.3, s=8, color="steelblue")
    # Regression line
    if len(x) >= 2:
        z = np.polyfit(x, y, 1)
        p = np.poly1d(z)
        x_sorted = sorted(x)
        ax.plot(x_sorted, [p(xi) for xi in x_sorted], "r--", linewidth=1.5,
                label=f"y = {z[0]:.4f}x + {z[1]:.2f}")
        ax.legend(fontsize=9)
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    ax.set_title(f"{y_label} vs {x_label}")
    fig.tight_layout()
    path = out / f"{filename}.{fmt}"
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path


# ---------------------------------------------------------------------------
# 6. Motif size heatmap by country (static)
# ---------------------------------------------------------------------------

def plot_country_motif_heatmap(
    heatmap_data: dict[str, dict[str, int]],
    output_dir: str | Path,
    fmt: str = "png",
) -> Path:
    """Heatmap of motif-size counts per country.

    *heatmap_data*: dict mapping country -> {motif_size_name: count}.
    """
    out = _ensure_dir(output_dir)
    countries = sorted(heatmap_data.keys())
    sizes = _MOTIF_SIZE_ORDER

    matrix = np.zeros((len(countries), len(sizes)))
    for i, c in enumerate(countries):
        for j, s in enumerate(sizes):
            matrix[i, j] = heatmap_data[c].get(s, 0)

    fig, ax = plt.subplots(figsize=(8, max(4, len(countries) * 0.35)))
    sns.heatmap(
        matrix, ax=ax, xticklabels=[_MOTIF_SIZE_LABELS.get(s, s) for s in sizes],
        yticklabels=countries, cmap="YlOrRd", annot=len(countries) <= 20,
        fmt=".0f", linewidths=0.5,
    )
    ax.set_title("Motif Size Distribution by Country")
    fig.tight_layout()
    path = out / f"country_motif_heatmap.{fmt}"
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path


# ---------------------------------------------------------------------------
# 7. Interactive choropleth map (plotly HTML)
# ---------------------------------------------------------------------------

def plot_choropleth(
    country_data: list[dict],
    output_dir: str | Path,
    metric_key: str = "total_ssrs",
    title: str = "Total SSRs by Country",
) -> Path:
    """Interactive choropleth map using Plotly.

    *country_data*: list of dicts with ``country`` and a numeric metric key.
    """
    import plotly.express as px

    out = _ensure_dir(output_dir)
    countries = [d["country"] for d in country_data]
    values = [d.get(metric_key, 0) for d in country_data]

    fig = px.choropleth(
        locations=countries,
        locationmode="country names",
        color=values,
        color_continuous_scale="Viridis",
        labels={"color": metric_key},
        title=title,
    )
    fig.update_layout(margin=dict(l=0, r=0, t=40, b=0))
    path = out / "choropleth_map.html"
    fig.write_html(str(path))
    return path


# ---------------------------------------------------------------------------
# 8. Interactive motif sunburst (plotly HTML)
# ---------------------------------------------------------------------------

def plot_motif_sunburst(
    motif_freq_data: list[dict],
    output_dir: str | Path,
) -> Path:
    """Sunburst chart of motif frequencies grouped by size.

    *motif_freq_data*: list of dicts with ``motif_canonical``, ``motif_size``,
    ``total_count``.
    """
    import plotly.express as px

    _SIZE_NAMES = {1: "mono", 2: "di", 3: "tri", 4: "tetra", 5: "penta", 6: "hexa"}
    out = _ensure_dir(output_dir)

    parents = []
    labels = []
    values = []
    for d in motif_freq_data:
        size_name = _SIZE_NAMES.get(d["motif_size"], f"size_{d['motif_size']}")
        parents.append(size_name)
        labels.append(d["motif_canonical"])
        values.append(d["total_count"])

    # Add size-level totals
    size_totals: dict[str, int] = {}
    for p, v in zip(parents, values):
        size_totals[p] = size_totals.get(p, 0) + v
    for size_name, total in size_totals.items():
        parents.append("")
        labels.append(size_name)
        values.append(total)

    import plotly.graph_objects as go
    fig = go.Figure(go.Sunburst(
        labels=labels, parents=parents, values=values,
        branchvalues="total",
    ))
    fig.update_layout(title="Motif Frequency Sunburst", margin=dict(t=40, l=0, r=0, b=0))
    path = out / "motif_sunburst.html"
    fig.write_html(str(path))
    return path


# ---------------------------------------------------------------------------
# 9. Statistical results summary table (static)
# ---------------------------------------------------------------------------

def plot_stats_summary_table(
    stats_results: list[dict],
    output_dir: str | Path,
    fmt: str = "png",
) -> Path:
    """Render key statistical results as a table figure.

    *stats_results*: list of dicts (from AnalysisSuite.to_list()).
    """
    out = _ensure_dir(output_dir)
    # Filter to tests with p-values
    rows = []
    for r in stats_results:
        if r.get("p_value") is not None:
            rows.append([
                r.get("analysis_name", ""),
                r.get("test_name", ""),
                f"{r['statistic']:.4f}" if r.get("statistic") is not None else "—",
                f"{r['p_value']:.2e}",
                f"{r['p_value_corrected']:.2e}" if r.get("p_value_corrected") is not None else "—",
                f"{r['effect_size']:.4f}" if r.get("effect_size") is not None else "—",
                str(r.get("n", "—")),
            ])

    if not rows:
        return out / f"stats_summary.{fmt}"

    col_labels = ["Analysis", "Test", "Statistic", "p-value", "p (FDR)", "Effect Size", "N"]
    fig, ax = plt.subplots(figsize=(14, max(2, len(rows) * 0.4 + 1)))
    ax.axis("off")
    table = ax.table(
        cellText=rows, colLabels=col_labels, loc="center",
        cellLoc="center",
    )
    table.auto_set_font_size(False)
    table.set_fontsize(8)
    table.scale(1, 1.3)
    # Header styling
    for j in range(len(col_labels)):
        table[0, j].set_facecolor("#4472C4")
        table[0, j].set_text_props(color="white", fontweight="bold")
    ax.set_title("Statistical Analysis Results", fontsize=12, fontweight="bold", pad=20)
    fig.tight_layout()
    path = out / f"stats_summary.{fmt}"
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return path


# ---------------------------------------------------------------------------
# Orchestrator: generate all figures for a run
# ---------------------------------------------------------------------------

def generate_all_figures(
    motif_size_data: list[dict],
    motif_freq_data: list[dict],
    gene_data: list[dict],
    country_summary: list[dict],
    country_metrics: dict[str, list[float]] | None,
    correlation_data: dict[str, tuple[list[float], list[float]]] | None,
    heatmap_data: dict[str, dict[str, int]] | None,
    stats_results: list[dict] | None,
    output_dir: str | Path,
    fmt: str = "png",
) -> list[Path]:
    """Generate all available figures and return paths to created files."""
    out = _ensure_dir(output_dir)
    paths: list[Path] = []

    # Static plots
    if motif_size_data:
        paths.append(plot_motif_size_distribution(motif_size_data, out, fmt))

    if motif_freq_data:
        paths.append(plot_top_motifs(motif_freq_data, out, fmt=fmt))

    if gene_data:
        paths.append(plot_gene_distribution(gene_data, out, fmt))

    if country_metrics:
        paths.append(plot_country_boxplot(country_metrics, out, fmt=fmt))

    if correlation_data:
        for name, (x, y) in correlation_data.items():
            if x and y:
                x_lab, y_lab = name.split("_vs_")
                paths.append(plot_correlation_scatter(
                    x, y, x_lab, y_lab, out,
                    filename=f"scatter_{name}", fmt=fmt,
                ))

    if heatmap_data:
        paths.append(plot_country_motif_heatmap(heatmap_data, out, fmt))

    if stats_results:
        paths.append(plot_stats_summary_table(stats_results, out, fmt))

    # Interactive plots
    if country_summary:
        paths.append(plot_choropleth(country_summary, out))

    if motif_freq_data:
        paths.append(plot_motif_sunburst(motif_freq_data, out))

    return paths
