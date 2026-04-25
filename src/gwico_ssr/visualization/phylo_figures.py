"""Phylogenetic tree visualization functions for GWICO-SSR.

Provides 5 visualization functions for phylo-aware SSR analysis:
- Tree structure with SSR overlay
- Clade metrics heatmap
- SSR distribution by clade
- Repeat class composition by clade
- Interactive Plotly tree visualization
"""

from __future__ import annotations

import logging
from typing import Optional

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import seaborn as sns
from plotly import graph_objects as go

from gwico_ssr.phylo.phylo_analysis import TreeCladeMetrics

logger = logging.getLogger(__name__)


def plot_tree_with_ssr_overlay(
    metrics_dict: dict[str, TreeCladeMetrics],
    figsize: tuple[int, int] = (14, 10),
) -> plt.Figure:
    """Plot tree structure with SSR metrics overlaid at tips.

    Args:
        metrics_dict: Dictionary of clade_name → TreeCladeMetrics
        figsize: Figure size (width, height)

    Returns:
        Matplotlib figure object
    """
    fig, ax = plt.subplots(figsize=figsize)

    if not metrics_dict:
        ax.text(0.5, 0.5, "No tree metrics available", ha="center", va="center")
        return fig

    # Extract data
    clades = list(metrics_dict.keys())
    ssr_counts = [metrics_dict[c].ssr_count_total for c in clades]
    accession_counts = [metrics_dict[c].accession_count for c in clades]

    # Create scatter plot with annotations
    scatter = ax.scatter(
        range(len(clades)),
        ssr_counts,
        s=[a * 10 + 100 for a in accession_counts],
        c=accession_counts,
        cmap="viridis",
        alpha=0.6,
        edgecolors="black",
        linewidth=1.5,
    )

    ax.set_xlabel("Clade (Tip Labels)", fontsize=12, fontweight="bold")
    ax.set_ylabel("Total SSR Count", fontsize=12, fontweight="bold")
    ax.set_title("Phylogenetic Tree with SSR Overlay", fontsize=14, fontweight="bold")
    ax.set_xticks(range(len(clades)))
    ax.set_xticklabels(clades, rotation=45, ha="right")
    ax.grid(True, alpha=0.3)

    # Add colorbar
    cbar = plt.colorbar(scatter, ax=ax)
    cbar.set_label("Accession Count", fontsize=10)

    # Add annotations
    for i, (clade, ssr_count, acc_count) in enumerate(zip(clades, ssr_counts, accession_counts)):
        ax.annotate(
            f"SSR: {int(ssr_count)}\nAcc: {int(acc_count)}",
            xy=(i, ssr_count),
            xytext=(0, 10),
            textcoords="offset points",
            ha="center",
            fontsize=8,
        )

    plt.tight_layout()
    return fig


def plot_clade_ssr_heatmap(
    metrics_dict: dict[str, TreeCladeMetrics],
    figsize: tuple[int, int] = (12, 8),
) -> plt.Figure:
    """Plot heatmap of SSR metrics across clades.

    Rows: clades, columns: metric types (total, perfect, imperfect, compound_component).

    Args:
        metrics_dict: Dictionary of clade_name → TreeCladeMetrics
        figsize: Figure size (width, height)

    Returns:
        Matplotlib figure object
    """
    fig, ax = plt.subplots(figsize=figsize)

    if not metrics_dict:
        ax.text(0.5, 0.5, "No tree metrics available", ha="center", va="center")
        return fig

    # Build matrix
    clades = list(metrics_dict.keys())
    metrics = ["ssr_count_total", "perfect_count", "imperfect_count", "compound_component_count"]
    data_matrix = np.zeros((len(clades), len(metrics)))

    for i, clade in enumerate(clades):
        m = metrics_dict[clade]
        data_matrix[i, 0] = m.ssr_count_total
        data_matrix[i, 1] = m.perfect_count
        data_matrix[i, 2] = m.imperfect_count
        data_matrix[i, 3] = m.compound_component_count

    # Create heatmap
    sns.heatmap(
        data_matrix,
        annot=True,
        fmt=".0f",
        cmap="YlOrRd",
        xticklabels=metrics,
        yticklabels=clades,
        ax=ax,
        cbar_kws={"label": "Count"},
    )

    ax.set_xlabel("Repeat Class", fontsize=12, fontweight="bold")
    ax.set_ylabel("Clade", fontsize=12, fontweight="bold")
    ax.set_title("SSR Metrics Heatmap by Clade", fontsize=14, fontweight="bold")
    plt.tight_layout()
    return fig


def plot_phylo_distribution(
    metrics_dict: dict[str, TreeCladeMetrics],
    figsize: tuple[int, int] = (12, 6),
) -> plt.Figure:
    """Bar chart of SSR counts by clade.

    Args:
        metrics_dict: Dictionary of clade_name → TreeCladeMetrics
        figsize: Figure size (width, height)

    Returns:
        Matplotlib figure object
    """
    fig, ax = plt.subplots(figsize=figsize)

    if not metrics_dict:
        ax.text(0.5, 0.5, "No tree metrics available", ha="center", va="center")
        return fig

    clades = list(metrics_dict.keys())
    ssr_counts = [metrics_dict[c].ssr_count_total for c in clades]
    accession_counts = [metrics_dict[c].accession_count for c in clades]

    x = np.arange(len(clades))
    width = 0.35

    bars1 = ax.bar(x - width / 2, ssr_counts, width, label="SSR Count", color="steelblue", alpha=0.8)
    bars2 = ax.bar(x + width / 2, accession_counts, width, label="Accession Count", color="coral", alpha=0.8)

    ax.set_xlabel("Clade", fontsize=12, fontweight="bold")
    ax.set_ylabel("Count", fontsize=12, fontweight="bold")
    ax.set_title("SSR Distribution by Clade", fontsize=14, fontweight="bold")
    ax.set_xticks(x)
    ax.set_xticklabels(clades, rotation=45, ha="right")
    ax.legend()
    ax.grid(True, alpha=0.3, axis="y")

    # Add value labels on bars
    for bars in [bars1, bars2]:
        for bar in bars:
            height = bar.get_height()
            ax.annotate(
                f"{int(height)}",
                xy=(bar.get_x() + bar.get_width() / 2, height),
                xytext=(0, 3),
                textcoords="offset points",
                ha="center",
                fontsize=8,
            )

    plt.tight_layout()
    return fig


def plot_phylo_class_composition(
    metrics_dict: dict[str, TreeCladeMetrics],
    figsize: tuple[int, int] = (12, 6),
) -> plt.Figure:
    """Stacked bar chart of repeat class composition by clade.

    Args:
        metrics_dict: Dictionary of clade_name → TreeCladeMetrics
        figsize: Figure size (width, height)

    Returns:
        Matplotlib figure object
    """
    fig, ax = plt.subplots(figsize=figsize)

    if not metrics_dict:
        ax.text(0.5, 0.5, "No tree metrics available", ha="center", va="center")
        return fig

    clades = list(metrics_dict.keys())
    perfect_counts = [metrics_dict[c].perfect_count for c in clades]
    imperfect_counts = [metrics_dict[c].imperfect_count for c in clades]
    compound_counts = [metrics_dict[c].compound_component_count for c in clades]

    x = np.arange(len(clades))

    p1 = ax.bar(x, perfect_counts, label="Perfect", color="#2ecc71", alpha=0.8)
    p2 = ax.bar(x, imperfect_counts, bottom=perfect_counts, label="Imperfect", color="#3498db", alpha=0.8)
    p3 = ax.bar(
        x,
        compound_counts,
        bottom=np.array(perfect_counts) + np.array(imperfect_counts),
        label="Compound Component",
        color="#e74c3c",
        alpha=0.8,
    )

    ax.set_xlabel("Clade", fontsize=12, fontweight="bold")
    ax.set_ylabel("Count", fontsize=12, fontweight="bold")
    ax.set_title("Repeat Class Composition by Clade", fontsize=14, fontweight="bold")
    ax.set_xticks(x)
    ax.set_xticklabels(clades, rotation=45, ha="right")
    ax.legend()
    ax.grid(True, alpha=0.3, axis="y")

    plt.tight_layout()
    return fig


def create_interactive_phylo_plot(
    metrics_dict: dict[str, TreeCladeMetrics],
    title: str = "Interactive Phylogenetic SSR Analysis",
) -> go.Figure:
    """Create interactive Plotly figure with tree and metrics.

    Args:
        metrics_dict: Dictionary of clade_name → TreeCladeMetrics
        title: Figure title

    Returns:
        Plotly figure object
    """
    if not metrics_dict:
        fig = go.Figure()
        fig.add_annotation(text="No tree metrics available", showarrow=False)
        return fig

    clades = list(metrics_dict.keys())
    ssr_counts = [metrics_dict[c].ssr_count_total for c in clades]
    accession_counts = [metrics_dict[c].accession_count for c in clades]
    perfect_counts = [metrics_dict[c].perfect_count for c in clades]
    imperfect_counts = [metrics_dict[c].imperfect_count for c in clades]
    compound_counts = [metrics_dict[c].compound_component_count for c in clades]

    # Create scatter plot
    fig = go.Figure()

    fig.add_trace(
        go.Scatter(
            x=clades,
            y=ssr_counts,
            mode="markers",
            marker=dict(
                size=[max(5, a / 2) for a in accession_counts],
                color=accession_counts,
                colorscale="Viridis",
                showscale=True,
                colorbar=dict(title="Accession Count"),
                line=dict(width=1, color="black"),
            ),
            text=[f"Clade: {c}<br>SSR: {s}<br>Accessions: {a}" for c, s, a in zip(clades, ssr_counts, accession_counts)],
            hovertemplate="%{text}<extra></extra>",
            name="SSR Distribution",
        )
    )

    # Add bar chart for composition
    fig.add_trace(
        go.Bar(
            x=clades,
            y=perfect_counts,
            name="Perfect",
            marker=dict(color="#2ecc71"),
            visible=False,
        )
    )

    fig.add_trace(
        go.Bar(
            x=clades,
            y=imperfect_counts,
            name="Imperfect",
            marker=dict(color="#3498db"),
            visible=False,
        )
    )

    fig.add_trace(
        go.Bar(
            x=clades,
            y=compound_counts,
            name="Compound Component",
            marker=dict(color="#e74c3c"),
            visible=False,
        )
    )

    # Add buttons for view switching
    fig.update_layout(
        updatemenus=[
            dict(
                type="buttons",
                direction="left",
                buttons=list(
                    [
                        dict(
                            args=[
                                {"visible": [True, False, False, False]},
                                {"title": title},
                            ],
                            label="Scatter View",
                            method="update",
                        ),
                        dict(
                            args=[
                                {"visible": [False, True, True, True]},
                                {"title": title + " - Repeat Composition"},
                            ],
                            label="Composition View",
                            method="update",
                        ),
                    ]
                ),
                pad={"r": 10, "t": 10},
                showactive=True,
                x=0.11,
                xanchor="left",
                y=1.15,
                yanchor="top",
            ),
        ]
    )

    fig.update_layout(
        title=title,
        xaxis_title="Clade",
        yaxis_title="Count",
        hovermode="closest",
        height=600,
        template="plotly_white",
    )

    return fig
