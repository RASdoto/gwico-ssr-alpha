"""Chunk 8: Tests for batch-aware annotation, metrics, and export downstream.

Tests verify:
1. repeat_class propagates through annotation pipeline
2. Metrics calculations distinguish repeat classes
3. Aggregations group by repeat_class correctly
4. Exports include repeat_class and provenance fields
5. Backward compatibility maintained for perfect-only workflows
"""

from __future__ import annotations

import pytest
from dataclasses import dataclass

from gwico_ssr.annotation.mapper import (
    AnnotationResult,
    SSRAnnotationRecord,
    annotate_ssrs,
    build_feature_tree,
    FeatureInterval,
)
from gwico_ssr.metrics.calculator import (
    AccessionMetricsResult,
    compute_accession_metrics,
    compute_repeat_class_breakdown,
    compute_repeat_class_bp_breakdown,
)


# ---------------------------------------------------------------------------
# Mock SSR and Feature records for testing
# ---------------------------------------------------------------------------

@dataclass
class MockSSR:
    """Mock SSRRecord for testing."""
    ssr_id: int
    accession: str
    start: int
    end: int
    motif_size: int
    repeat_length_bp: int
    motif_canonical: str = "AAG"
    repeat_class: str = "perfect"  # perfect, imperfect, compound_component
    is_perfect: bool = True
    imperfection_pct: float = None


@dataclass
class MockFeature:
    """Mock FeatureRecord for testing."""
    feature_id: int
    accession: str
    feature_type: str
    start: int
    end: int
    strand: str
    gene_name: str


# ---------------------------------------------------------------------------
# Test annotation repeat_class propagation
# ---------------------------------------------------------------------------

class TestAnnotationRepeatClassPropagation:
    """Verify repeat_class propagates from SSRRecord to annotations."""

    def test_annotation_record_has_repeat_class_field(self):
        """SSRAnnotationRecord should have repeat_class field."""
        record = SSRAnnotationRecord(
            ssr_id=1,
            accession="test",
            feature_id=None,
            gene_name=None,
            region_class="intergenic",
            overlap_bp=None,
            repeat_class="imperfect",
        )
        assert record.repeat_class == "imperfect"

    def test_annotation_record_to_dicts_includes_repeat_class(self):
        """to_dicts() should include repeat_class in output."""
        result = AnnotationResult(accession="test", total_ssrs=1)
        result.annotations.append(SSRAnnotationRecord(
            ssr_id=1,
            accession="test",
            feature_id=None,
            gene_name=None,
            region_class="intergenic",
            overlap_bp=None,
            repeat_class="compound_component",
        ))
        dicts = result.to_dicts()
        assert len(dicts) == 1
        assert dicts[0]["repeat_class"] == "compound_component"

    def test_annotate_ssrs_preserves_perfect_class(self):
        """annotate_ssrs should preserve perfect SSRs as perfect."""
        ssrs = [
            MockSSR(
                ssr_id=1, accession="test", start=0, end=9,
                motif_size=3, repeat_length_bp=9,
                repeat_class="perfect"
            )
        ]
        features = []  # No features - will be intergenic
        tree = build_feature_tree([])
        
        result = annotate_ssrs(ssrs, tree, "test")
        
        assert len(result.annotations) == 1
        assert result.annotations[0].repeat_class == "perfect"
        assert result.annotations[0].region_class == "intergenic"

    def test_annotate_ssrs_preserves_imperfect_class(self):
        """annotate_ssrs should preserve imperfect SSRs as imperfect."""
        ssrs = [
            MockSSR(
                ssr_id=1, accession="test", start=0, end=9,
                motif_size=3, repeat_length_bp=9,
                repeat_class="imperfect",
                is_perfect=False,
                imperfection_pct=5.5,
            )
        ]
        tree = build_feature_tree([])
        
        result = annotate_ssrs(ssrs, tree, "test")
        
        assert len(result.annotations) == 1
        assert result.annotations[0].repeat_class == "imperfect"

    def test_annotate_ssrs_preserves_compound_component_class(self):
        """annotate_ssrs should preserve compound_component SSRs."""
        ssrs = [
            MockSSR(
                ssr_id=1, accession="test", start=0, end=9,
                motif_size=3, repeat_length_bp=9,
                repeat_class="compound_component",
            )
        ]
        tree = build_feature_tree([])
        
        result = annotate_ssrs(ssrs, tree, "test")
        
        assert len(result.annotations) == 1
        assert result.annotations[0].repeat_class == "compound_component"

    def test_annotate_ssrs_with_feature_overlap_preserves_class(self):
        """annotate_ssrs should preserve repeat_class when overlapping features."""
        ssrs = [
            MockSSR(
                ssr_id=1, accession="test", start=50, end=60,
                motif_size=3, repeat_length_bp=10,
                repeat_class="imperfect",
            )
        ]
        features = [
            FeatureInterval(
                feature_id=1, accession="test", feature_type="CDS",
                start=40, end=70, strand="+", gene_name="geneA"
            )
        ]
        tree = build_feature_tree(features)
        
        result = annotate_ssrs(ssrs, tree, "test")
        
        assert len(result.annotations) == 1
        assert result.annotations[0].repeat_class == "imperfect"
        assert result.annotations[0].region_class == "CDS"
        assert result.annotations[0].gene_name == "geneA"


# ---------------------------------------------------------------------------
# Test metrics repeat_class breakdown
# ---------------------------------------------------------------------------

class TestMetricsRepeatClassBreakdown:
    """Verify metrics calculator distinguishes repeat classes."""

    def test_accession_metrics_result_has_class_fields(self):
        """AccessionMetricsResult should have repeat_class breakdown fields."""
        result = AccessionMetricsResult(accession="test")
        assert hasattr(result, "perfect_count")
        assert hasattr(result, "imperfect_count")
        assert hasattr(result, "compound_component_count")
        assert hasattr(result, "perfect_bp_total")
        assert hasattr(result, "imperfect_bp_total")
        assert hasattr(result, "compound_component_bp_total")

    def test_compute_repeat_class_breakdown_perfect_only(self):
        """compute_repeat_class_breakdown should count perfect SSRs."""
        ssrs = [
            MockSSR(ssr_id=1, accession="test", start=0, end=9,
                   motif_size=3, repeat_length_bp=9, repeat_class="perfect"),
            MockSSR(ssr_id=2, accession="test", start=20, end=29,
                   motif_size=3, repeat_length_bp=9, repeat_class="perfect"),
        ]
        counts = compute_repeat_class_breakdown(ssrs)
        assert counts["perfect_count"] == 2
        assert counts["imperfect_count"] == 0
        assert counts["compound_component_count"] == 0

    def test_compute_repeat_class_breakdown_mixed(self):
        """compute_repeat_class_breakdown should count all classes."""
        ssrs = [
            MockSSR(ssr_id=1, accession="test", start=0, end=9,
                   motif_size=3, repeat_length_bp=9, repeat_class="perfect"),
            MockSSR(ssr_id=2, accession="test", start=20, end=29,
                   motif_size=3, repeat_length_bp=9, repeat_class="imperfect"),
            MockSSR(ssr_id=3, accession="test", start=40, end=49,
                   motif_size=3, repeat_length_bp=9, repeat_class="compound_component"),
            MockSSR(ssr_id=4, accession="test", start=60, end=69,
                   motif_size=3, repeat_length_bp=9, repeat_class="imperfect"),
        ]
        counts = compute_repeat_class_breakdown(ssrs)
        assert counts["perfect_count"] == 1
        assert counts["imperfect_count"] == 2
        assert counts["compound_component_count"] == 1

    def test_compute_repeat_class_bp_breakdown_perfect_only(self):
        """compute_repeat_class_bp_breakdown should sum perfect bp."""
        ssrs = [
            MockSSR(ssr_id=1, accession="test", start=0, end=15,
                   motif_size=3, repeat_length_bp=15, repeat_class="perfect"),
            MockSSR(ssr_id=2, accession="test", start=20, end=35,
                   motif_size=3, repeat_length_bp=15, repeat_class="perfect"),
        ]
        bp_totals = compute_repeat_class_bp_breakdown(ssrs)
        assert bp_totals["perfect_bp_total"] == 30
        assert bp_totals["imperfect_bp_total"] == 0
        assert bp_totals["compound_component_bp_total"] == 0

    def test_compute_repeat_class_bp_breakdown_mixed(self):
        """compute_repeat_class_bp_breakdown should sum all classes."""
        ssrs = [
            MockSSR(ssr_id=1, accession="test", start=0, end=10,
                   motif_size=3, repeat_length_bp=10, repeat_class="perfect"),
            MockSSR(ssr_id=2, accession="test", start=20, end=35,
                   motif_size=3, repeat_length_bp=15, repeat_class="imperfect"),
            MockSSR(ssr_id=3, accession="test", start=40, end=52,
                   motif_size=3, repeat_length_bp=12, repeat_class="compound_component"),
        ]
        bp_totals = compute_repeat_class_bp_breakdown(ssrs)
        assert bp_totals["perfect_bp_total"] == 10
        assert bp_totals["imperfect_bp_total"] == 15
        assert bp_totals["compound_component_bp_total"] == 12

    def test_compute_accession_metrics_includes_class_breakdown(self):
        """compute_accession_metrics should populate class breakdown fields."""
        ssrs = [
            MockSSR(ssr_id=1, accession="test", start=0, end=9,
                   motif_size=3, repeat_length_bp=9, repeat_class="perfect"),
            MockSSR(ssr_id=2, accession="test", start=20, end=34,
                   motif_size=3, repeat_length_bp=14, repeat_class="imperfect"),
        ]
        result = compute_accession_metrics("test", ssrs, genome_length=30000)
        
        assert result.ssr_count_total == 2
        assert result.ssr_bp_total == 23
        assert result.perfect_count == 1
        assert result.imperfect_count == 1
        assert result.compound_component_count == 0
        assert result.perfect_bp_total == 9
        assert result.imperfect_bp_total == 14

    def test_accession_metrics_to_dict_includes_class_fields(self):
        """AccessionMetricsResult.to_dict should include class fields."""
        result = AccessionMetricsResult(
            accession="test",
            ssr_count_total=3,
            ssr_bp_total=30,
            perfect_count=2,
            imperfect_count=1,
            compound_component_count=0,
            perfect_bp_total=18,
            imperfect_bp_total=12,
            compound_component_bp_total=0,
        )
        d = result.to_dict(run_id=1)
        
        assert d["perfect_count"] == 2
        assert d["imperfect_count"] == 1
        assert d["compound_component_count"] == 0
        assert d["perfect_bp_total"] == 18
        assert d["imperfect_bp_total"] == 12
        assert d["compound_component_bp_total"] == 0


# ---------------------------------------------------------------------------
# Test backward compatibility
# ---------------------------------------------------------------------------

class TestDownstreamBackwardCompatibility:
    """Verify batch-aware changes maintain perfect-only workflow compatibility."""

    def test_perfect_only_annotation_workflow(self):
        """Perfect-only workflows should work unchanged."""
        ssrs = [
            MockSSR(ssr_id=1, accession="acc1", start=0, end=9,
                   motif_size=3, repeat_length_bp=9, repeat_class="perfect"),
        ]
        tree = build_feature_tree([])
        
        result = annotate_ssrs(ssrs, tree, "acc1")
        
        assert result.total_ssrs == 1
        assert result.intergenic == 1
        assert result.annotated == 0

    def test_perfect_only_metrics_workflow(self):
        """Perfect-only metrics workflows should work unchanged."""
        ssrs = [
            MockSSR(ssr_id=1, accession="acc1", start=0, end=9,
                   motif_size=3, repeat_length_bp=9,
                   motif_canonical="AAG", repeat_class="perfect"),
        ]
        result = compute_accession_metrics("acc1", ssrs, genome_length=30000)
        
        assert result.ssr_count_total == 1
        assert result.ssr_bp_total == 9
        assert result.ra is not None
        assert result.dominant_motif == "AAG"

    def test_metrics_default_class_is_perfect(self):
        """SSR records without repeat_class should default to perfect."""
        @dataclass
        class MinimalSSR:
            ssr_id: int
            accession: str
            motif_size: int
            repeat_length_bp: int
            motif_canonical: str
            # No repeat_class field
        
        ssrs = [
            MinimalSSR(
                ssr_id=1, accession="acc1", motif_size=3,
                repeat_length_bp=9, motif_canonical="AAG"
            ),
        ]
        counts = compute_repeat_class_breakdown(ssrs)
        assert counts["perfect_count"] == 1
        assert counts["imperfect_count"] == 0

    def test_annotation_default_class_is_perfect(self):
        """SSR records without repeat_class should default to perfect in annotation."""
        @dataclass
        class MinimalSSR:
            ssr_id: int
            accession: str
            start: int
            end: int
            # No repeat_class field
        
        ssrs = [
            MinimalSSR(ssr_id=1, accession="test", start=0, end=9),
        ]
        tree = build_feature_tree([])
        
        result = annotate_ssrs(ssrs, tree, "test")
        
        assert result.annotations[0].repeat_class == "perfect"


# ---------------------------------------------------------------------------
# Test edge cases
# ---------------------------------------------------------------------------

class TestDownstreamEdgeCases:
    """Test edge cases in batch-aware downstream processing."""

    def test_empty_ssr_list_metrics(self):
        """Empty SSR list should return zero counts."""
        result = compute_accession_metrics("empty", [], genome_length=30000)
        
        assert result.ssr_count_total == 0
        assert result.ssr_bp_total == 0
        assert result.perfect_count == 0
        assert result.imperfect_count == 0
        assert result.compound_component_count == 0

    def test_multiple_overlapping_features_with_class(self):
        """Repeat_class should be consistent across multiple feature overlaps."""
        ssrs = [
            MockSSR(
                ssr_id=1, accession="test", start=50, end=60,
                motif_size=3, repeat_length_bp=10,
                repeat_class="imperfect",
            )
        ]
        features = [
            FeatureInterval(
                feature_id=1, accession="test", feature_type="CDS",
                start=40, end=70, strand="+", gene_name="geneA"
            ),
            FeatureInterval(
                feature_id=2, accession="test", feature_type="gene",
                start=30, end=80, strand="+", gene_name="geneB"
            ),
        ]
        tree = build_feature_tree(features)
        
        result = annotate_ssrs(ssrs, tree, "test")
        
        # Should create two annotation records, both with imperfect class
        assert len(result.annotations) == 2
        assert all(a.repeat_class == "imperfect" for a in result.annotations)

    def test_all_ssr_classes_present(self):
        """Metrics should handle all three repeat classes correctly."""
        ssrs = [
            MockSSR(ssr_id=1, accession="test", start=0, end=9,
                   motif_size=1, repeat_length_bp=9, repeat_class="perfect"),
            MockSSR(ssr_id=2, accession="test", start=20, end=34,
                   motif_size=2, repeat_length_bp=14, repeat_class="imperfect"),
            MockSSR(ssr_id=3, accession="test", start=40, end=52,
                   motif_size=3, repeat_length_bp=12, repeat_class="compound_component"),
            MockSSR(ssr_id=4, accession="test", start=60, end=63,
                   motif_size=1, repeat_length_bp=3, repeat_class="perfect"),
            MockSSR(ssr_id=5, accession="test", start=70, end=87,
                   motif_size=4, repeat_length_bp=17, repeat_class="imperfect"),
        ]
        result = compute_accession_metrics("test", ssrs, genome_length=100000)
        
        assert result.perfect_count == 2
        assert result.imperfect_count == 2
        assert result.compound_component_count == 1
        assert result.perfect_bp_total == 12
        assert result.imperfect_bp_total == 31
        assert result.compound_component_bp_total == 12
        # Totals should match sum
        assert result.ssr_count_total == 5
        assert result.ssr_bp_total == 55
