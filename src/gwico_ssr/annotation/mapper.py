"""Interval-based SSR-to-feature annotation mapping.

Uses an interval tree (from the ``intervaltree`` package) to efficiently
map SSR coordinates to overlapping genomic features. SSRs that do not
overlap any feature are classified as **intergenic**.

Coordinate convention:
    All intervals are 0-based half-open [start, end).

Performance:
    Interval tree construction is O(n log n).  Each overlap query is
    O(log n + k) where k is the number of overlapping features.  This
    replaces the O(n * m) nested-loop approach used previously.
"""

from __future__ import annotations

import logging
import time
from dataclasses import dataclass, field
from typing import Optional, Sequence

from intervaltree import IntervalTree

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------

@dataclass
class FeatureInterval:
    """A genomic feature interval used for building the interval tree."""

    feature_id: int
    accession: str
    feature_type: str
    start: int   # 0-based inclusive
    end: int     # 0-based exclusive
    strand: str | None
    gene_name: str | None


@dataclass
class SSRAnnotationRecord:
    """An SSR-to-feature annotation result ready for persistence.
    
    Chunk 8: Now includes repeat_class for batch-aware downstream processing.
    """

    ssr_id: int
    accession: str
    feature_id: int | None
    gene_name: str | None
    region_class: str   # CDS, gene, mRNA, intergenic, other
    overlap_bp: int | None
    repeat_class: str = "perfect"  # perfect, imperfect, compound_component


@dataclass
class AnnotationResult:
    """Result of annotation mapping for a single accession."""

    accession: str
    total_ssrs: int
    annotated: int = 0
    intergenic: int = 0
    annotations: list[SSRAnnotationRecord] = field(default_factory=list)

    def to_dicts(self) -> list[dict]:
        """Convert annotations to dicts suitable for insert_ssr_annotations.
        
        Chunk 8: Includes repeat_class for batch-aware downstream operations.
        """
        return [
            {
                "ssr_id": a.ssr_id,
                "accession": a.accession,
                "feature_id": a.feature_id,
                "gene_name": a.gene_name,
                "region_class": a.region_class,
                "overlap_bp": a.overlap_bp,
                "repeat_class": a.repeat_class,
            }
            for a in self.annotations
        ]


# ---------------------------------------------------------------------------
# Interval tree builder
# ---------------------------------------------------------------------------

def build_feature_tree(features: Sequence[FeatureInterval]) -> IntervalTree:
    """Build an interval tree from feature intervals.

    Args:
        features: Sequence of FeatureInterval objects.

    Returns:
        An IntervalTree indexed on [start, end) with FeatureInterval as data.
    """
    tree = IntervalTree()
    for f in features:
        if f.end > f.start:  # skip zero-length features
            tree.addi(f.start, f.end, f)
    return tree


def features_from_db(feature_records) -> list[FeatureInterval]:
    """Convert ORM FeatureRecord objects to FeatureInterval dataclasses.

    Args:
        feature_records: Sequence of FeatureRecord ORM objects.

    Returns:
        List of FeatureInterval objects.
    """
    return [
        FeatureInterval(
            feature_id=r.feature_id,
            accession=r.accession,
            feature_type=r.feature_type,
            start=r.start,
            end=r.end,
            strand=r.strand,
            gene_name=r.gene_name,
        )
        for r in feature_records
    ]


# ---------------------------------------------------------------------------
# Overlap calculation
# ---------------------------------------------------------------------------

def _compute_overlap_bp(ssr_start: int, ssr_end: int,
                        feat_start: int, feat_end: int) -> int:
    """Compute the number of overlapping base pairs between two intervals.

    Both intervals are [start, end) half-open.
    """
    overlap_start = max(ssr_start, feat_start)
    overlap_end = min(ssr_end, feat_end)
    return max(0, overlap_end - overlap_start)


def _classify_region(feature_type: str) -> str:
    """Map a feature_type string to a region_class.

    Returns the feature_type directly for recognized types,
    or 'other' for unrecognized ones.
    """
    known = {"CDS", "gene", "mRNA", "tRNA", "rRNA", "ncRNA", "exon", "misc_feature"}
    return feature_type if feature_type in known else "other"


# ---------------------------------------------------------------------------
# Core mapping function
# ---------------------------------------------------------------------------

def annotate_ssrs(
    ssr_records,
    feature_tree: IntervalTree,
    accession: str,
) -> AnnotationResult:
    """Map SSR records to overlapping genomic features using an interval tree.

    Each SSR may overlap zero or more features (one-to-many).
    SSRs with no overlapping features are classified as intergenic.
    
    Chunk 8: Propagates repeat_class from SSRRecord to annotations.

    Args:
        ssr_records: Sequence of SSRRecord ORM objects (must have ssr_id,
            start, end, repeat_class attributes).
        feature_tree: Pre-built IntervalTree from build_feature_tree().
        accession: Accession identifier.

    Returns:
        AnnotationResult with all annotation records, including repeat_class.
    """
    result = AnnotationResult(accession=accession, total_ssrs=len(ssr_records))

    for ssr in ssr_records:
        # Get repeat_class from SSR record; default to "perfect" for backward compatibility
        repeat_class = getattr(ssr, "repeat_class", "perfect")
        
        overlaps = feature_tree.overlap(ssr.start, ssr.end)

        if not overlaps:
            # Intergenic — no overlapping features
            result.annotations.append(SSRAnnotationRecord(
                ssr_id=ssr.ssr_id,
                accession=accession,
                feature_id=None,
                gene_name=None,
                region_class="intergenic",
                overlap_bp=None,
                repeat_class=repeat_class,
            ))
            result.intergenic += 1
        else:
            # One or more overlapping features
            for iv in overlaps:
                feat: FeatureInterval = iv.data
                overlap_bp = _compute_overlap_bp(ssr.start, ssr.end,
                                                  feat.start, feat.end)
                result.annotations.append(SSRAnnotationRecord(
                    ssr_id=ssr.ssr_id,
                    accession=accession,
                    feature_id=feat.feature_id,
                    gene_name=feat.gene_name,
                    region_class=_classify_region(feat.feature_type),
                    overlap_bp=overlap_bp,
                    repeat_class=repeat_class,
                ))
            result.annotated += 1

    return result


# ---------------------------------------------------------------------------
# High-level pipeline function
# ---------------------------------------------------------------------------

def annotate_accession(session, accession: str) -> AnnotationResult:
    """Full annotation pipeline for a single accession.

    Loads features and SSRs from the database, builds an interval tree,
    maps SSRs to features, and returns the annotation result.

    Args:
        session: SQLAlchemy session.
        accession: Accession identifier.

    Returns:
        AnnotationResult with all mappings.
    """
    from gwico_ssr.db.repository import (
        get_features_for_accession,
        get_ssr_records_for_accession,
    )

    features_db = get_features_for_accession(session, accession)
    ssr_records = get_ssr_records_for_accession(session, accession)

    if not ssr_records:
        logger.debug("No SSRs for %s — skipping annotation", accession)
        return AnnotationResult(accession=accession, total_ssrs=0)

    feature_intervals = features_from_db(features_db)
    tree = build_feature_tree(feature_intervals)

    logger.debug(
        "Annotating %s: %d SSRs against %d features",
        accession, len(ssr_records), len(feature_intervals),
    )

    return annotate_ssrs(ssr_records, tree, accession)


# ---------------------------------------------------------------------------
# Naive mapping for comparison / benchmarking
# ---------------------------------------------------------------------------

def annotate_ssrs_naive(
    ssr_records,
    features: Sequence[FeatureInterval],
    accession: str,
) -> AnnotationResult:
    """Naive O(n*m) nested-loop SSR annotation for benchmarking comparison.

    This function is intentionally slow and exists only so tests can
    demonstrate the interval-tree approach is materially faster.
    """
    result = AnnotationResult(accession=accession, total_ssrs=len(ssr_records))

    for ssr in ssr_records:
        found = False
        for feat in features:
            overlap_bp = _compute_overlap_bp(ssr.start, ssr.end,
                                              feat.start, feat.end)
            if overlap_bp > 0:
                result.annotations.append(SSRAnnotationRecord(
                    ssr_id=ssr.ssr_id,
                    accession=accession,
                    feature_id=feat.feature_id,
                    gene_name=feat.gene_name,
                    region_class=_classify_region(feat.feature_type),
                    overlap_bp=overlap_bp,
                ))
                found = True

        if not found:
            result.annotations.append(SSRAnnotationRecord(
                ssr_id=ssr.ssr_id,
                accession=accession,
                feature_id=None,
                gene_name=None,
                region_class="intergenic",
                overlap_bp=None,
            ))
            result.intergenic += 1
        else:
            result.annotated += 1

    return result


def benchmark_compare(
    ssr_records,
    features: Sequence[FeatureInterval],
    accession: str,
) -> dict:
    """Run both naive and interval-tree approaches and compare timing.

    Returns a dict with timing info and correctness comparison.
    """
    # Naive
    t0 = time.perf_counter()
    naive_result = annotate_ssrs_naive(ssr_records, features, accession)
    naive_time = time.perf_counter() - t0

    # Interval tree
    t0 = time.perf_counter()
    tree = build_feature_tree(features)
    tree_result = annotate_ssrs(ssr_records, tree, accession)
    tree_time = time.perf_counter() - t0

    return {
        "naive_time_ms": naive_time * 1000,
        "tree_time_ms": tree_time * 1000,
        "speedup": naive_time / tree_time if tree_time > 0 else float("inf"),
        "naive_annotations": len(naive_result.annotations),
        "tree_annotations": len(tree_result.annotations),
        "results_match": len(naive_result.annotations) == len(tree_result.annotations),
    }
