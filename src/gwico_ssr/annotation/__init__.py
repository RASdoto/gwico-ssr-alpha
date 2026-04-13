"""SSR-to-feature annotation mapping using interval trees."""

from gwico_ssr.annotation.mapper import (
    AnnotationResult,
    FeatureInterval,
    SSRAnnotationRecord,
    annotate_accession,
    annotate_ssrs,
    annotate_ssrs_naive,
    benchmark_compare,
    build_feature_tree,
    features_from_db,
)

__all__ = [
    "AnnotationResult",
    "FeatureInterval",
    "SSRAnnotationRecord",
    "annotate_accession",
    "annotate_ssrs",
    "annotate_ssrs_naive",
    "benchmark_compare",
    "build_feature_tree",
    "features_from_db",
]
