"""SSR metric computation and cohort-level aggregation."""

from gwico_ssr.metrics.calculator import (
    AccessionMetricsResult,
    compute_accession_metrics,
    compute_dominant_motif,
    compute_metrics_for_accession,
    compute_motif_size_counts,
    compute_ra,
    compute_rd,
)
from gwico_ssr.metrics.aggregator import (
    CohortSummary,
    GeneSSRSummary,
    MotifFrequency,
    aggregate_by_country,
    aggregate_by_gene,
    aggregate_by_motif_size,
    aggregate_motif_frequencies,
    get_dataset_summary,
)

__all__ = [
    "AccessionMetricsResult",
    "CohortSummary",
    "GeneSSRSummary",
    "MotifFrequency",
    "aggregate_by_country",
    "aggregate_by_gene",
    "aggregate_by_motif_size",
    "aggregate_motif_frequencies",
    "compute_accession_metrics",
    "compute_dominant_motif",
    "compute_metrics_for_accession",
    "compute_motif_size_counts",
    "compute_ra",
    "compute_rd",
    "get_dataset_summary",
]
