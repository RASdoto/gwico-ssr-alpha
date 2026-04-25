"""Lineage visualization functions for GWICO-SSR.

Generates publication-quality figures for lineage-aware analysis including
distributions, trends, compositions, and interactive dashboards.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Optional

import matplotlib.pyplot as plt
import numpy as np
import plotly.graph_objects as go
import seaborn as sns

from gwico_ssr.analysis.lineage import LineageSummary

logger = logging.getLogger(__name__)


def plot_lineage_distribution(
    summaries: dict[str, LineageSummary],
    title: str = "SSR Distribution by Lineage",
    save_path: Optional[str | Path] = None,
) -> None:
    """Bar chart of accession counts by lineage.

    Args:
        summaries: Dict of lineage_label → LineageSummary
        title: Figure title
        save_path: Optional file to save figure
    """
    lineages = sorted(summaries.keys())
    counts = [summaries[l].accession_count for l in lineages]

    fig, ax = plt.subplots(figsize=(12, 6))
    ax.bar(lineages, counts, color="steelblue", edgecolor="black", linewidth=0.5)
    ax.set_xlabel("Lineage", fontsize=12, fontweight="bold")
    ax.set_ylabel("Accession Count", fontsize=12, fontweight="bold")
    ax.set_title(title, fontsize=14, fontweight="bold")
    ax.tick_params(axis="x", rotation=45)
    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches="tight")
        logger.info(f"Saved lineage distribution plot to {save_path}")
    else:
        plt.show()

    plt.close()


def plot_ssr_by_lineage(
    summaries: dict[str, LineageSummary],
    title: str = "Total SSR Count by Lineage",
    save_path: Optional[str | Path] = None,
) -> None:
    """Bar chart of total SSR count by lineage.

    Args:
        summaries: Dict of lineage_label → LineageSummary
        title: Figure title
        save_path: Optional file to save figure
    """
    lineages = sorted(summaries.keys())
    ssr_counts = [summaries[l].ssr_count_total for l in lineages]

    fig, ax = plt.subplots(figsize=(12, 6))
    ax.bar(lineages, ssr_counts, color="coral", edgecolor="black", linewidth=0.5)
    ax.set_xlabel("Lineage", fontsize=12, fontweight="bold")
    ax.set_ylabel("Total SSR Count", fontsize=12, fontweight="bold")
    ax.set_title(title, fontsize=14, fontweight="bold")
    ax.tick_params(axis="x", rotation=45)
    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches="tight")
        logger.info(f"Saved SSR by lineage plot to {save_path}")
    else:
        plt.show()

    plt.close()


def plot_lineage_class_composition(
    summaries: dict[str, LineageSummary],
    title: str = "Repeat Class Composition by Lineage",
    save_path: Optional[str | Path] = None,
) -> None:
    """Stacked bar chart showing perfect/imperfect/compound composition by lineage.

    Args:
        summaries: Dict of lineage_label → LineageSummary
        title: Figure title
        save_path: Optional file to save figure
    """
    lineages = sorted(summaries.keys())
    perfect = np.array([summaries[l].perfect_count for l in lineages])
    imperfect = np.array([summaries[l].imperfect_count for l in lineages])
    compound = np.array([summaries[l].compound_component_count for l in lineages])

    fig, ax = plt.subplots(figsize=(12, 6))
    x = np.arange(len(lineages))
    width = 0.6

    ax.bar(x, perfect, width, label="Perfect", color="#2ecc71", edgecolor="black", linewidth=0.5)
    ax.bar(x, imperfect, width, bottom=perfect, label="Imperfect", color="#e74c3c", edgecolor="black", linewidth=0.5)
    ax.bar(x, compound, width, bottom=perfect + imperfect, label="Compound", color="#3498db", edgecolor="black", linewidth=0.5)

    ax.set_xlabel("Lineage", fontsize=12, fontweight="bold")
    ax.set_ylabel("SSR Count", fontsize=12, fontweight="bold")
    ax.set_title(title, fontsize=14, fontweight="bold")
    ax.set_xticks(x)
    ax.set_xticklabels(lineages, rotation=45)
    ax.legend(loc="upper right")
    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches="tight")
        logger.info(f"Saved lineage class composition plot to {save_path}")
    else:
        plt.show()

    plt.close()


def plot_lineage_metrics_heatmap(
    summaries: dict[str, LineageSummary],
    title: str = "Lineage Metrics Heatmap (Normalized)",
    save_path: Optional[str | Path] = None,
) -> None:
    """Heatmap of normalized metrics by lineage.

    Args:
        summaries: Dict of lineage_label → LineageSummary
        title: Figure title
        save_path: Optional file to save figure
    """
    lineages = sorted(summaries.keys())
    metrics = ["accession_count", "ssr_count_total", "perfect_count", "imperfect_count", "compound_component_count"]

    data = []
    for lineage in lineages:
        row = [
            summaries[lineage].accession_count,
            summaries[lineage].ssr_count_total,
            summaries[lineage].perfect_count,
            summaries[lineage].imperfect_count,
            summaries[lineage].compound_component_count,
        ]
        data.append(row)

    data = np.array(data, dtype=float)

    # Normalize by column (z-score normalization)
    data_norm = (data - data.mean(axis=0)) / (data.std(axis=0) + 1e-10)

    fig, ax = plt.subplots(figsize=(10, len(lineages) * 0.4))
    im = ax.imshow(data_norm, cmap="RdYlGn", aspect="auto")

    ax.set_xticks(np.arange(len(metrics)))
    ax.set_yticks(np.arange(len(lineages)))
    ax.set_xticklabels(metrics, rotation=45, ha="right")
    ax.set_yticklabels(lineages)
    ax.set_title(title, fontsize=14, fontweight="bold")

    # Add colorbar
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label("Normalized Value", fontsize=12)

    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches="tight")
        logger.info(f"Saved lineage metrics heatmap to {save_path}")
    else:
        plt.show()

    plt.close()


def create_interactive_lineage_plot(
    summaries: dict[str, LineageSummary],
    title: str = "Interactive Lineage Analysis Dashboard",
    save_path: Optional[str | Path] = None,
) -> None:
    """Interactive Plotly figure with lineage dropdown filtering.

    Args:
        summaries: Dict of lineage_label → LineageSummary
        title: Figure title
        save_path: Optional file to save figure
    """
    lineages = sorted(summaries.keys())

    # Create traces for each metric
    fig = go.Figure()

    # Bar chart for accession count
    fig.add_trace(go.Bar(
        x=lineages,
        y=[summaries[l].accession_count for l in lineages],
        name="Accessions",
        marker_color="steelblue",
        visible=True,
    ))

    # Bar chart for SSR count
    fig.add_trace(go.Bar(
        x=lineages,
        y=[summaries[l].ssr_count_total for l in lineages],
        name="Total SSRs",
        marker_color="coral",
        visible=False,
    ))

    # Stacked bar: repeat classes
    fig.add_trace(go.Bar(
        x=lineages,
        y=[summaries[l].perfect_count for l in lineages],
        name="Perfect",
        marker_color="#2ecc71",
        visible=False,
    ))

    fig.add_trace(go.Bar(
        x=lineages,
        y=[summaries[l].imperfect_count for l in lineages],
        name="Imperfect",
        marker_color="#e74c3c",
        visible=False,
    ))

    fig.add_trace(go.Bar(
        x=lineages,
        y=[summaries[l].compound_component_count for l in lineages],
        name="Compound",
        marker_color="#3498db",
        visible=False,
    ))

    # Create dropdown buttons
    buttons = [
        dict(
            label="Accession Count",
            method="update",
            args=[{"visible": [True, False, False, False, False]},
                  {"title": "Accession Distribution by Lineage"}],
        ),
        dict(
            label="Total SSR Count",
            method="update",
            args=[{"visible": [False, True, False, False, False]},
                  {"title": "Total SSR Count by Lineage"}],
        ),
        dict(
            label="Repeat Classes (Stacked)",
            method="update",
            args=[{"visible": [False, False, True, True, True]},
                  {"title": "Repeat Class Composition by Lineage"}],
        ),
    ]

    fig.update_layout(
        updatemenus=[
            dict(
                buttons=buttons,
                direction="down",
                pad={"r": 10, "t": 10},
                showactive=True,
                x=0.0,
                xanchor="left",
                y=1.15,
                yanchor="top",
            ),
        ],
        title=title,
        xaxis_title="Lineage",
        yaxis_title="Count",
        hovermode="x unified",
        height=600,
    )

    if save_path:
        fig.write_html(str(save_path))
        logger.info(f"Saved interactive lineage plot to {save_path}")
    else:
        fig.show()
