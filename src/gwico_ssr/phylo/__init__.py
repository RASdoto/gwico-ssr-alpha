"""Phylogenetic tree integration module for Chunk 11.

Provides tree parsing, ingestion, mapping, and phylo-aware metrics computation.
"""

from gwico_ssr.phylo.tree_parser import (
    Tree,
    TreeNode,
    parse_newick_file,
    validate_tree_structure,
)
from gwico_ssr.phylo.phylo_analysis import (
    TreeSource,
    TreeMapping,
    TreeCladeMetrics,
    TreeMappingReport,
    ingest_newick_tree,
    map_accessions_to_tips,
    compute_ssr_metrics_by_clade,
    store_tree_metrics,
)

__all__ = [
    "Tree",
    "TreeNode",
    "parse_newick_file",
    "validate_tree_structure",
    "TreeSource",
    "TreeMapping",
    "TreeCladeMetrics",
    "TreeMappingReport",
    "ingest_newick_tree",
    "map_accessions_to_tips",
    "compute_ssr_metrics_by_clade",
    "store_tree_metrics",
]
