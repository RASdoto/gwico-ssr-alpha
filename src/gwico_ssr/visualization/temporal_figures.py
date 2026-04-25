"""Temporal visualization for GWICO-SSR.

Generates time-series figures from temporal analysis results.
Supports both static (matplotlib/seaborn) and interactive (plotly) outputs.
"""

from __future__ import annotations

import logging
from typing import Optional

import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.dates import DateFormatter, MonthLocator, YearLocator
import plotly.graph_objects as go
import plotly.express as px

from gwico_ssr.analysis.temporal import DatePrecision, TemporalSummary

logger = logging.getLogger(__name__)


def plot_temporal_distribution(
    summary: TemporalSummary,
    figsize: tuple[int, int] = (12, 6),
    save_path: Optional[str] = None,
) -> plt.Figure:
    """Plot accession count distribution over time.

    Args:
        summary: TemporalSummary object with populated bins
        figsize: Figure dimensions
        save_path: Optional path to save figure

    Returns:
        Matplotlib Figure object
    """
    if not summary.bins:
        logger.warning("No bins to plot")
        fig, ax = plt.subplots(figsize=figsize)
        ax.text(0.5, 0.5, "No temporal data available", ha="center", va="center")
        return fig

    fig, ax = plt.subplots(figsize=figsize)

    bin_ids = [b.bin_id for b in summary.bins]
    accession_counts = [b.accession_count for b in summary.bins]

    ax.plot(range(len(bin_ids)), accession_counts, marker="o", linewidth=2, markersize=6)
    ax.fill_between(range(len(bin_ids)), accession_counts, alpha=0.3)

    ax.set_xlabel("Time Bin")
    ax.set_ylabel("Accession Count")
    ax.set_title(f"Accession Distribution Over Time ({summary.date_field}, {summary.precision.value})")
    ax.set_xticks(range(len(bin_ids)))
    ax.set_xticklabels(bin_ids, rotation=45, ha="right")
    ax.grid(True, alpha=0.3)

    fig.tight_layout()

    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches="tight")
        logger.info(f"Saved temporal distribution plot to {save_path}")

    return fig


def plot_temporal_ssr_trends(
    summary: TemporalSummary,
    figsize: tuple[int, int] = (12, 6),
    save_path: Optional[str] = None,
) -> plt.Figure:
    """Plot SSR counts by repeat class over time.

    Shows trends for perfect, imperfect, and compound_component SSRs.

    Args:
        summary: TemporalSummary object with populated bins
        figsize: Figure dimensions
        save_path: Optional path to save figure

    Returns:
        Matplotlib Figure object
    """
    if not summary.bins:
        logger.warning("No bins to plot")
        fig, ax = plt.subplots(figsize=figsize)
        ax.text(0.5, 0.5, "No temporal data available", ha="center", va="center")
        return fig

    fig, ax = plt.subplots(figsize=figsize)

    bin_ids = [b.bin_id for b in summary.bins]
    x = range(len(bin_ids))

    perfect = [b.perfect_count for b in summary.bins]
    imperfect = [b.imperfect_count for b in summary.bins]
    compound = [b.compound_component_count for b in summary.bins]

    ax.plot(x, perfect, marker="o", label="Perfect", linewidth=2, markersize=6)
    ax.plot(x, imperfect, marker="s", label="Imperfect", linewidth=2, markersize=6)
    ax.plot(x, compound, marker="^", label="Compound Component", linewidth=2, markersize=6)

    ax.set_xlabel("Time Bin")
    ax.set_ylabel("SSR Count")
    ax.set_title(f"SSR Repeat Class Trends Over Time ({summary.date_field}, {summary.precision.value})")
    ax.set_xticks(x)
    ax.set_xticklabels(bin_ids, rotation=45, ha="right")
    ax.legend()
    ax.grid(True, alpha=0.3)

    fig.tight_layout()

    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches="tight")
        logger.info(f"Saved temporal SSR trends plot to {save_path}")

    return fig


def plot_temporal_bp_trends(
    summary: TemporalSummary,
    figsize: tuple[int, int] = (12, 6),
    save_path: Optional[str] = None,
) -> plt.Figure:
    """Plot total repeat length (bp) by repeat class over time.

    Args:
        summary: TemporalSummary object with populated bins
        figsize: Figure dimensions
        save_path: Optional path to save figure

    Returns:
        Matplotlib Figure object
    """
    if not summary.bins:
        logger.warning("No bins to plot")
        fig, ax = plt.subplots(figsize=figsize)
        ax.text(0.5, 0.5, "No temporal data available", ha="center", va="center")
        return fig

    fig, ax = plt.subplots(figsize=figsize)

    bin_ids = [b.bin_id for b in summary.bins]
    x = range(len(bin_ids))

    perfect_bp = [b.perfect_bp_total for b in summary.bins]
    imperfect_bp = [b.imperfect_bp_total for b in summary.bins]
    compound_bp = [b.compound_component_bp_total for b in summary.bins]

    ax.plot(x, perfect_bp, marker="o", label="Perfect", linewidth=2, markersize=6)
    ax.plot(x, imperfect_bp, marker="s", label="Imperfect", linewidth=2, markersize=6)
    ax.plot(x, compound_bp, marker="^", label="Compound Component", linewidth=2, markersize=6)

    ax.set_xlabel("Time Bin")
    ax.set_ylabel("Total Repeat Length (bp)")
    ax.set_title(f"Repeat Length Trends Over Time ({summary.date_field}, {summary.precision.value})")
    ax.set_xticks(x)
    ax.set_xticklabels(bin_ids, rotation=45, ha="right")
    ax.legend()
    ax.grid(True, alpha=0.3)

    fig.tight_layout()

    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches="tight")
        logger.info(f"Saved temporal bp trends plot to {save_path}")

    return fig


def plot_temporal_metrics_trends(
    summary: TemporalSummary,
    figsize: tuple[int, int] = (12, 6),
    save_path: Optional[str] = None,
) -> plt.Figure:
    """Plot mean RA/RD metrics over time.

    RA = Repeat Abundance (SSR count per Mb)
    RD = Repeat Density (total repeat bp per Mb)

    Args:
        summary: TemporalSummary object with populated bins
        figsize: Figure dimensions
        save_path: Optional path to save figure

    Returns:
        Matplotlib Figure object
    """
    if not summary.bins:
        logger.warning("No bins to plot")
        fig, ax = plt.subplots(figsize=figsize)
        ax.text(0.5, 0.5, "No temporal data available", ha="center", va="center")
        return fig

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=figsize)

    bin_ids = [b.bin_id for b in summary.bins]
    x = range(len(bin_ids))

    # Filter out None values
    mean_ra = [b.mean_ra if b.mean_ra is not None else 0 for b in summary.bins]
    mean_rd = [b.mean_rd if b.mean_rd is not None else 0 for b in summary.bins]

    ax1.plot(x, mean_ra, marker="o", color="steelblue", linewidth=2, markersize=6)
    ax1.set_ylabel("Mean RA")
    ax1.set_title(f"Mean Repeat Abundance Over Time ({summary.date_field})")
    ax1.set_xticks(x)
    ax1.set_xticklabels(bin_ids, rotation=45, ha="right")
    ax1.grid(True, alpha=0.3)

    ax2.plot(x, mean_rd, marker="o", color="coral", linewidth=2, markersize=6)
    ax2.set_xlabel("Time Bin")
    ax2.set_ylabel("Mean RD")
    ax2.set_title(f"Mean Repeat Density Over Time ({summary.date_field})")
    ax2.set_xticks(x)
    ax2.set_xticklabels(bin_ids, rotation=45, ha="right")
    ax2.grid(True, alpha=0.3)

    fig.tight_layout()

    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches="tight")
        logger.info(f"Saved temporal metrics trends plot to {save_path}")

    return fig


def create_interactive_temporal_plot(
    summary: TemporalSummary,
    save_path: Optional[str] = None,
) -> go.Figure:
    """Create interactive plotly figure showing accession and SSR trends.

    Includes:
    - Accession count (bar chart, primary y-axis)
    - SSR counts by class (line chart, secondary y-axis)

    Args:
        summary: TemporalSummary object with populated bins
        save_path: Optional path to save HTML figure

    Returns:
        Plotly Figure object
    """
    if not summary.bins:
        logger.warning("No bins to plot")
        fig = go.Figure()
        fig.add_annotation(text="No temporal data available", xref="paper", yref="paper",
                          x=0.5, y=0.5, showarrow=False)
        return fig

    fig = go.Figure()

    bin_ids = [b.bin_id for b in summary.bins]

    # Accession count (bar chart)
    accession_counts = [b.accession_count for b in summary.bins]
    fig.add_trace(go.Bar(
        x=bin_ids,
        y=accession_counts,
        name="Accessions",
        marker_color="lightblue",
        yaxis="y",
    ))

    # SSR counts by class (lines)
    perfect = [b.perfect_count for b in summary.bins]
    imperfect = [b.imperfect_count for b in summary.bins]
    compound = [b.compound_component_count for b in summary.bins]

    fig.add_trace(go.Scatter(
        x=bin_ids,
        y=perfect,
        name="Perfect SSRs",
        mode="lines+markers",
        yaxis="y2",
    ))

    fig.add_trace(go.Scatter(
        x=bin_ids,
        y=imperfect,
        name="Imperfect SSRs",
        mode="lines+markers",
        yaxis="y2",
    ))

    fig.add_trace(go.Scatter(
        x=bin_ids,
        y=compound,
        name="Compound SSRs",
        mode="lines+markers",
        yaxis="y2",
    ))

    # Update layout with dual y-axes
    fig.update_layout(
        title=f"Temporal SSR Analysis ({summary.date_field}, {summary.precision.value})",
        xaxis=dict(title="Time Bin"),
        yaxis=dict(title="Accession Count", side="left"),
        yaxis2=dict(title="SSR Count", overlaying="y", side="right"),
        hovermode="x unified",
        template="plotly_white",
        height=600,
        width=1000,
    )

    if save_path:
        fig.write_html(save_path)
        logger.info(f"Saved interactive temporal plot to {save_path}")

    return fig


def plot_temporal_class_stacked(
    summary: TemporalSummary,
    figsize: tuple[int, int] = (12, 6),
    save_path: Optional[str] = None,
) -> plt.Figure:
    """Plot stacked bar chart of SSR counts by repeat class over time.

    Args:
        summary: TemporalSummary object with populated bins
        figsize: Figure dimensions
        save_path: Optional path to save figure

    Returns:
        Matplotlib Figure object
    """
    if not summary.bins:
        logger.warning("No bins to plot")
        fig, ax = plt.subplots(figsize=figsize)
        ax.text(0.5, 0.5, "No temporal data available", ha="center", va="center")
        return fig

    fig, ax = plt.subplots(figsize=figsize)

    bin_ids = [b.bin_id for b in summary.bins]
    x = range(len(bin_ids))

    perfect = [b.perfect_count for b in summary.bins]
    imperfect = [b.imperfect_count for b in summary.bins]
    compound = [b.compound_component_count for b in summary.bins]

    ax.bar(x, perfect, label="Perfect", color="steelblue")
    ax.bar(x, imperfect, bottom=perfect, label="Imperfect", color="orange")
    
    compound_bottom = [p + i for p, i in zip(perfect, imperfect)]
    ax.bar(x, compound, bottom=compound_bottom, label="Compound Component", color="green")

    ax.set_xlabel("Time Bin")
    ax.set_ylabel("SSR Count")
    ax.set_title(f"SSR Class Composition Over Time ({summary.date_field}, {summary.precision.value})")
    ax.set_xticks(x)
    ax.set_xticklabels(bin_ids, rotation=45, ha="right")
    ax.legend()
    ax.grid(True, alpha=0.3, axis="y")

    fig.tight_layout()

    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches="tight")
        logger.info(f"Saved temporal class stacked plot to {save_path}")

    return fig
