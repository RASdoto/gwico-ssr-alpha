"""Statistical analysis layer for GWICO-SSR."""

from gwico_ssr.analysis.statistics import (
    AnalysisResult,
    AnalysisSuite,
    base_composition_test,
    chi_square_gene_country,
    chi_square_motif_country,
    correlation_gc_ssr,
    correlation_length_ssr,
    kruskal_wallis_by_country,
    run_all_analyses,
    shannon_entropy_by_country,
    shannon_entropy_motifs,
)

__all__ = [
    "AnalysisResult",
    "AnalysisSuite",
    "base_composition_test",
    "chi_square_gene_country",
    "chi_square_motif_country",
    "correlation_gc_ssr",
    "correlation_length_ssr",
    "kruskal_wallis_by_country",
    "run_all_analyses",
    "shannon_entropy_by_country",
    "shannon_entropy_motifs",
]
