"""Visualization layer for GWICO-SSR."""

from gwico_ssr.visualization.figures import (
    generate_all_figures,
    plot_choropleth,
    plot_correlation_scatter,
    plot_country_boxplot,
    plot_country_motif_heatmap,
    plot_gene_distribution,
    plot_motif_size_distribution,
    plot_motif_sunburst,
    plot_stats_summary_table,
    plot_top_motifs,
)

__all__ = [
    "generate_all_figures",
    "plot_choropleth",
    "plot_correlation_scatter",
    "plot_country_boxplot",
    "plot_country_motif_heatmap",
    "plot_gene_distribution",
    "plot_motif_size_distribution",
    "plot_motif_sunburst",
    "plot_stats_summary_table",
    "plot_top_motifs",
]
