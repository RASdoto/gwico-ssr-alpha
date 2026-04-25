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
from gwico_ssr.visualization.temporal_figures import (
    plot_temporal_distribution,
    plot_temporal_ssr_trends,
    plot_temporal_bp_trends,
    plot_temporal_metrics_trends,
    create_interactive_temporal_plot,
    plot_temporal_class_stacked,
)
from gwico_ssr.visualization.lineage_figures import (
    plot_lineage_distribution,
    plot_ssr_by_lineage,
    plot_lineage_class_composition,
    plot_lineage_metrics_heatmap,
    create_interactive_lineage_plot,
)
from gwico_ssr.visualization.phylo_figures import (
    plot_tree_with_ssr_overlay,
    plot_clade_ssr_heatmap,
    plot_phylo_distribution,
    plot_phylo_class_composition,
    create_interactive_phylo_plot,
)

__all__ = [
    # Existing figures
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
    # Temporal figures
    "plot_temporal_distribution",
    "plot_temporal_ssr_trends",
    "plot_temporal_bp_trends",
    "plot_temporal_metrics_trends",
    "create_interactive_temporal_plot",
    "plot_temporal_class_stacked",
    # Lineage figures
    "plot_lineage_distribution",
    "plot_ssr_by_lineage",
    "plot_lineage_class_composition",
    "plot_lineage_metrics_heatmap",
    "create_interactive_lineage_plot",
    # Phylogenetic figures (Chunk 11)
    "plot_tree_with_ssr_overlay",
    "plot_clade_ssr_heatmap",
    "plot_phylo_distribution",
    "plot_phylo_class_composition",
    "create_interactive_phylo_plot",
]
