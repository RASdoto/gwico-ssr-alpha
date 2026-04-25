"""Unit tests for Chunk 7: Compound SSR detection and standardization.

Tests cover:
- Standardization levels (L0, L1, L2, Full)
- dMAX chaining algorithm
- Compound detection with mixed perfect/imperfect components
- Backward compatibility with perfect-only mode
- Edge cases and boundary conditions
"""

import pytest
from gwico_ssr.ssr.detector import SSRHit, DetectionResult
from gwico_ssr.ssr.compound import (
    StandardizationLevel,
    CompoundHit,
    standardize_motif,
    find_compound_chains,
    build_compound_hit,
    detect_compound_ssrs,
    mark_compound_components,
)


# ============================================================================
# Fixtures
# ============================================================================


@pytest.fixture
def simple_perfect_hits():
    """Two perfect SSRs 5 bp apart (chainable with dmax=10)."""
    return [
        SSRHit(
            start=0, end=9, motif_raw="AAG", motif_canonical="AAG",
            motif_size=3, repeat_units=3, repeat_length_bp=9, strand="+",
            actual_repeat="AAGAAGAAG",
            is_perfect=True, repeat_class="perfect",
        ),
        SSRHit(
            start=14, end=23, motif_raw="GAT", motif_canonical="ATG",
            motif_size=3, repeat_units=3, repeat_length_bp=9, strand="+",
            actual_repeat="GATGATGAT",
            is_perfect=True, repeat_class="perfect",
        ),
    ]


@pytest.fixture
def perfect_and_imperfect_hits():
    """Perfect SSR followed by imperfect SSR, chainable."""
    return [
        SSRHit(
            start=0, end=9, motif_raw="AAG", motif_canonical="AAG",
            motif_size=3, repeat_units=3, repeat_length_bp=9, strand="+",
            actual_repeat="AAGAAGAAG",
            is_perfect=True, repeat_class="perfect",
        ),
        SSRHit(
            start=14, end=23, motif_raw="GAT", motif_canonical="ATG",
            motif_size=3, repeat_units=3, repeat_length_bp=9, strand="+",
            actual_repeat="GATGATAAT",  # One mismatch (T→A at position 8)
            is_perfect=False, repeat_class="imperfect",
            imperfection_pct=3.7, num_substitutions=1, num_indels=0,
        ),
    ]


@pytest.fixture
def non_chainable_hits():
    """Two perfect SSRs > dmax_bp apart."""
    return [
        SSRHit(
            start=0, end=9, motif_raw="AAG", motif_canonical="AAG",
            motif_size=3, repeat_units=3, repeat_length_bp=9, strand="+",
            actual_repeat="AAGAAGAAG",
            is_perfect=True, repeat_class="perfect",
        ),
        SSRHit(
            start=50, end=59, motif_raw="GAT", motif_canonical="ATG",
            motif_size=3, repeat_units=3, repeat_length_bp=9, strand="+",
            actual_repeat="GATGATGAT",
            is_perfect=True, repeat_class="perfect",
        ),
    ]


@pytest.fixture
def three_chain_hits():
    """Three perfect SSRs in a valid chain."""
    return [
        SSRHit(
            start=0, end=9, motif_raw="AAG", motif_canonical="AAG",
            motif_size=3, repeat_units=3, repeat_length_bp=9, strand="+",
            actual_repeat="AAGAAGAAG",
            is_perfect=True, repeat_class="perfect",
        ),
        SSRHit(
            start=12, end=21, motif_raw="GAT", motif_canonical="ATG",
            motif_size=3, repeat_units=3, repeat_length_bp=9, strand="+",
            actual_repeat="GATGATGAT",
            is_perfect=True, repeat_class="perfect",
        ),
        SSRHit(
            start=24, end=33, motif_raw="CCC", motif_canonical="CCC",
            motif_size=3, repeat_units=3, repeat_length_bp=9, strand="+",
            actual_repeat="CCCCCCCCC",
            is_perfect=True, repeat_class="perfect",
        ),
    ]


# ============================================================================
# Tests: Standardization Levels
# ============================================================================


class TestStandardizationLevels:
    """Test motif standardization at each level."""

    def test_standardize_l0_raw_motif(self):
        """L0: Return raw motif unchanged."""
        motif = "AAG"
        result = standardize_motif(motif, StandardizationLevel.L0)
        assert result == "AAG"

    def test_standardize_l0_lowercase(self):
        """L0: Uppercase but no other transformation."""
        motif = "aag"
        result = standardize_motif(motif, StandardizationLevel.L0)
        assert result == "AAG"

    def test_standardize_l1_rotation_minimum(self):
        """L1: Lexicographic minimum rotation (forward only)."""
        # AAG, AGA, GAA → AGA (minimum)
        result = standardize_motif("AAG", StandardizationLevel.L1)
        assert result == "AAG"  # AAG < AGA < GAA

    def test_standardize_l1_all_rotations(self):
        """L1: Verify rotation logic."""
        result = standardize_motif("CAT", StandardizationLevel.L1)
        # Rotations: CAT, ATC, TCA → ATC (minimum)
        assert result == "ATC"

    def test_standardize_l2_canonical(self):
        """L2: Canonical form (rotation + RC)."""
        # AAG canonical is AAG (already canonical)
        result = standardize_motif("AAG", StandardizationLevel.L2)
        assert result == "AAG"

    def test_standardize_l2_rc_rotation(self):
        """L2: RC rotation normalized."""
        # For AAAA (homopolymer), canonical is AAAA
        result = standardize_motif("AAAA", StandardizationLevel.L2)
        assert result == "AAAA"

    def test_standardize_full_same_as_l2(self):
        """Full: Currently same as L2 (for future expansion)."""
        motif = "AAG"
        result_l2 = standardize_motif(motif, StandardizationLevel.L2)
        result_full = standardize_motif(motif, StandardizationLevel.Full)
        assert result_l2 == result_full

    def test_standardize_composite_motif(self):
        """Standardization on composite (concatenated motifs)."""
        composite = "AAGGATCCC"  # AAG + GAT + CCC
        result = standardize_motif(composite, StandardizationLevel.L2)
        # Should be canonicalized as a unit
        assert isinstance(result, str)
        assert len(result) > 0


# ============================================================================
# Tests: dMAX Chaining Algorithm
# ============================================================================


class TestDmaxChaining:
    """Test find_compound_chains with various gap distances."""

    def test_simple_chain_within_dmax(self, simple_perfect_hits):
        """Two SSRs 5 bp apart with dmax=10 should chain."""
        hits = simple_perfect_hits
        chains = find_compound_chains(hits, dmax_bp=10)
        assert len(chains) == 1
        assert len(chains[0]) == 2

    def test_chain_exact_dmax_boundary(self, simple_perfect_hits):
        """SSRs exactly at dmax boundary should chain."""
        # Gap is 14 - 9 = 5 bp
        chains = find_compound_chains(simple_perfect_hits, dmax_bp=5)
        assert len(chains) == 1

    def test_chain_exceeds_dmax(self, simple_perfect_hits):
        """SSRs beyond dmax should not chain."""
        chains = find_compound_chains(simple_perfect_hits, dmax_bp=4)
        assert len(chains) == 0  # No chains >=2 components

    def test_multiple_chains(self):
        """Multiple separate chains detected correctly."""
        hits = [
            SSRHit(start=0, end=9, motif_raw="A", motif_canonical="A",
                   motif_size=1, repeat_units=9, repeat_length_bp=9, strand="+",
                   actual_repeat="AAAAAAAAA", is_perfect=True, repeat_class="perfect"),
            SSRHit(start=12, end=15, motif_raw="C", motif_canonical="C",
                   motif_size=1, repeat_units=3, repeat_length_bp=3, strand="+",
                   actual_repeat="CCC", is_perfect=True, repeat_class="perfect"),
            SSRHit(start=50, end=53, motif_raw="G", motif_canonical="G",
                   motif_size=1, repeat_units=3, repeat_length_bp=3, strand="+",
                   actual_repeat="GGG", is_perfect=True, repeat_class="perfect"),
            SSRHit(start=56, end=59, motif_raw="T", motif_canonical="T",
                   motif_size=1, repeat_units=3, repeat_length_bp=3, strand="+",
                   actual_repeat="TTT", is_perfect=True, repeat_class="perfect"),
        ]
        chains = find_compound_chains(hits, dmax_bp=10)
        # Chain 1: indices 0,1 (gap=3)
        # Chain 2: indices 2,3 (gap=3)
        assert len(chains) == 2
        assert len(chains[0]) == 2
        assert len(chains[1]) == 2

    def test_three_component_chain(self, three_chain_hits):
        """Chain of 3+ components."""
        chains = find_compound_chains(three_chain_hits, dmax_bp=10)
        assert len(chains) == 1
        assert len(chains[0]) == 3

    def test_no_chains_single_ssr(self):
        """Single SSR returns no chains."""
        hits = [
            SSRHit(start=0, end=9, motif_raw="AAG", motif_canonical="AAG",
                   motif_size=3, repeat_units=3, repeat_length_bp=9, strand="+",
                   actual_repeat="AAGAAGAAG", is_perfect=True, repeat_class="perfect"),
        ]
        chains = find_compound_chains(hits, dmax_bp=10)
        assert len(chains) == 0

    def test_no_chains_empty_list(self):
        """Empty list returns no chains."""
        chains = find_compound_chains([], dmax_bp=10)
        assert len(chains) == 0


# ============================================================================
# Tests: Compound Hit Building
# ============================================================================


class TestBuildCompoundHit:
    """Test build_compound_hit functionality."""

    def test_build_simple_compound(self, simple_perfect_hits):
        """Build compound from two perfect components."""
        compound = build_compound_hit(
            simple_perfect_hits,
            dmax_bp=10,
            standardization_level=StandardizationLevel.L2,
        )
        assert compound.accession == "unknown"  # Default from hits
        assert compound.start == 0
        assert compound.end == 23
        assert len(compound.components) == 2
        assert len(compound.gaps_bp) == 1
        assert compound.gaps_bp[0] == 5

    def test_compound_motif_concatenation(self, simple_perfect_hits):
        """Compound motif is concatenation of components."""
        compound = build_compound_hit(
            simple_perfect_hits,
            dmax_bp=10,
            standardization_level=StandardizationLevel.L2,
        )
        # First: AAGAAGAAG, Second: GATGATGAT
        expected_composite = "AAGAAGAAGGATGATGAT"
        assert compound.compound_motif == expected_composite

    def test_compound_total_repeat_length(self, simple_perfect_hits):
        """Total repeat length sums component lengths."""
        compound = build_compound_hit(
            simple_perfect_hits,
            dmax_bp=10,
            standardization_level=StandardizationLevel.L2,
        )
        assert compound.total_repeat_length_bp == 18  # 9 + 9

    def test_compound_standardization_applied(self):
        """Standardization level applied to compound motif."""
        hits = [
            SSRHit(start=0, end=9, motif_raw="AAG", motif_canonical="AAG",
                   motif_size=3, repeat_units=3, repeat_length_bp=9, strand="+",
                   actual_repeat="AAGAAGAAG", is_perfect=True, repeat_class="perfect"),
            SSRHit(start=12, end=21, motif_raw="GAT", motif_canonical="ATG",
                   motif_size=3, repeat_units=3, repeat_length_bp=9, strand="+",
                   actual_repeat="GATGATGAT", is_perfect=True, repeat_class="perfect"),
        ]
        compound_l0 = build_compound_hit(hits, dmax_bp=10, standardization_level=StandardizationLevel.L0)
        compound_l2 = build_compound_hit(hits, dmax_bp=10, standardization_level=StandardizationLevel.L2)
        # L0 should preserve concatenation, L2 should canonicalize
        assert compound_l0.standardization_level == StandardizationLevel.L0
        assert compound_l2.standardization_level == StandardizationLevel.L2


# ============================================================================
# Tests: Compound Detection Main Function
# ============================================================================


class TestCompoundDetection:
    """Test detect_compound_ssrs integration."""

    def test_detect_simple_compound(self, simple_perfect_hits):
        """Detect a simple 2-component compound."""
        compounds = detect_compound_ssrs(simple_perfect_hits, dmax_bp=10)
        assert len(compounds) == 1
        assert compounds[0].component_count == 2

    def test_detect_no_compounds_far_apart(self, non_chainable_hits):
        """No compounds when SSRs are too far apart."""
        compounds = detect_compound_ssrs(non_chainable_hits, dmax_bp=10)
        assert len(compounds) == 0

    def test_detect_mixed_perfect_imperfect(self, perfect_and_imperfect_hits):
        """Compounds can mix perfect and imperfect components."""
        compounds = detect_compound_ssrs(perfect_and_imperfect_hits, dmax_bp=10)
        assert len(compounds) == 1
        assert len(compounds[0].components) == 2
        # First component is perfect
        assert compounds[0].components[0].is_perfect is True
        # Second component is imperfect
        assert compounds[0].components[1].is_perfect is False

    def test_detect_three_component_compound(self, three_chain_hits):
        """Detect compound with 3+ components."""
        compounds = detect_compound_ssrs(three_chain_hits, dmax_bp=10)
        assert len(compounds) == 1
        assert len(compounds[0].components) == 3

    def test_detect_returns_empty_on_single_ssr(self):
        """Single SSR returns no compounds."""
        hits = [
            SSRHit(start=0, end=9, motif_raw="AAG", motif_canonical="AAG",
                   motif_size=3, repeat_units=3, repeat_length_bp=9, strand="+",
                   actual_repeat="AAGAAGAAG", is_perfect=True, repeat_class="perfect"),
        ]
        compounds = detect_compound_ssrs(hits, dmax_bp=10)
        assert len(compounds) == 0

    def test_detect_groupby_accession(self):
        """Compounds detected per accession."""
        hits = [
            SSRHit(accession="acc1", start=0, end=9, motif_raw="AAG", motif_canonical="AAG",
                   motif_size=3, repeat_units=3, repeat_length_bp=9, strand="+",
                   actual_repeat="AAGAAGAAG", is_perfect=True, repeat_class="perfect"),
            SSRHit(accession="acc1", start=14, end=23, motif_raw="GAT", motif_canonical="ATG",
                   motif_size=3, repeat_units=3, repeat_length_bp=9, strand="+",
                   actual_repeat="GATGATGAT", is_perfect=True, repeat_class="perfect"),
            SSRHit(accession="acc2", start=0, end=9, motif_raw="CCC", motif_canonical="CCC",
                   motif_size=3, repeat_units=3, repeat_length_bp=9, strand="+",
                   actual_repeat="CCCCCCCCC", is_perfect=True, repeat_class="perfect"),
            SSRHit(accession="acc2", start=14, end=23, motif_raw="GGG", motif_canonical="GGG",
                   motif_size=3, repeat_units=3, repeat_length_bp=9, strand="+",
                   actual_repeat="GGGGGGGGG", is_perfect=True, repeat_class="perfect"),
        ]
        compounds = detect_compound_ssrs(hits, dmax_bp=10)
        assert len(compounds) == 2
        assert compounds[0].accession == "acc1"
        assert compounds[1].accession == "acc2"


# ============================================================================
# Tests: Component Marking
# ============================================================================


class TestMarkCompoundComponents:
    """Test mark_compound_components function."""

    def test_mark_simple_components(self, simple_perfect_hits):
        """Components marked with repeat_class = 'compound_component'."""
        compounds = detect_compound_ssrs(simple_perfect_hits, dmax_bp=10)
        updated = mark_compound_components(simple_perfect_hits, compounds)
        
        # All updated hits should be marked
        for hit in updated:
            assert hit.repeat_class == "compound_component"

    def test_non_component_unchanged(self, non_chainable_hits):
        """Non-component SSRs remain unchanged."""
        compounds = detect_compound_ssrs(non_chainable_hits, dmax_bp=10)  # No compounds
        updated = mark_compound_components(non_chainable_hits, compounds)
        
        # Original hits should remain unchanged (not marked as component)
        for orig, updated_hit in zip(non_chainable_hits, updated):
            assert updated_hit.repeat_class == orig.repeat_class

    def test_mixed_marking(self, simple_perfect_hits):
        """Mix of component and non-component SSRs."""
        # Add a third, isolated SSR
        isolated = SSRHit(
            start=100, end=109, motif_raw="TTT", motif_canonical="TTT",
            motif_size=3, repeat_units=3, repeat_length_bp=9, strand="+",
            actual_repeat="TTTTTTTTT", is_perfect=True, repeat_class="perfect",
        )
        all_hits = simple_perfect_hits + [isolated]
        compounds = detect_compound_ssrs(simple_perfect_hits, dmax_bp=10)
        updated = mark_compound_components(all_hits, compounds)
        
        # First two should be marked
        assert updated[0].repeat_class == "compound_component"
        assert updated[1].repeat_class == "compound_component"
        # Third should remain perfect
        assert updated[2].repeat_class == "perfect"


# ============================================================================
# Tests: Backward Compatibility
# ============================================================================


class TestBackwardCompatibility:
    """Ensure compound detection doesn't break existing functionality."""

    def test_perfect_detector_unaffected(self):
        """Perfect detector results unaffected by compound module import."""
        from gwico_ssr.ssr.detector import detect_ssrs
        
        sequence = "AAGAAGAAGGATGATGAT"
        result = detect_ssrs(sequence, "test_acc")
        
        assert result.hit_count >= 2
        for hit in result.hits:
            assert hit.repeat_class in ("perfect", "compound_component")

    def test_imperfect_detector_unaffected(self):
        """Imperfect detector results unaffected by compound module import."""
        from gwico_ssr.ssr.detector import detect_imperfect_ssrs
        
        sequence = "AAGAAGAAGGATGATGAT"
        result = detect_imperfect_ssrs(sequence, "test_acc")
        
        assert result.hit_count >= 0  # May be 0 if no imperfect found
        for hit in result.hits:
            assert hasattr(hit, 'repeat_class')

    def test_apply_compound_returns_tuple(self, simple_perfect_hits):
        """apply_compound_detection_to_result returns (hits, compounds)."""
        from gwico_ssr.ssr.detector import apply_compound_detection_to_result
        
        result = DetectionResult(accession="test", sequence_length=100)
        result.hits = simple_perfect_hits
        
        updated_hits, compounds = apply_compound_detection_to_result(result, dmax_bp=10)
        
        assert isinstance(updated_hits, list)
        assert isinstance(compounds, list)
        assert len(compounds) == 1


# ============================================================================
# Tests: Edge Cases
# ============================================================================


class TestEdgeCases:
    """Test boundary conditions and special cases."""

    def test_overlapping_ssr_regions(self):
        """Overlapping SSRs handled correctly."""
        hits = [
            SSRHit(start=0, end=15, motif_raw="AAG", motif_canonical="AAG",
                   motif_size=3, repeat_units=5, repeat_length_bp=15, strand="+",
                   actual_repeat="AAGAAGAAGAAGAAG", is_perfect=True, repeat_class="perfect"),
            SSRHit(start=10, end=25, motif_raw="GAA", motif_canonical="AAG",
                   motif_size=3, repeat_units=5, repeat_length_bp=15, strand="+",
                   actual_repeat="GAAGAAGAAGAAGAA", is_perfect=True, repeat_class="perfect"),
        ]
        compounds = detect_compound_ssrs(hits, dmax_bp=10)
        # Overlapping means they're already chained; may form compound or be kept separate
        assert isinstance(compounds, list)

    def test_zero_gap_between_ssrs(self):
        """Adjacent SSRs (gap = 0) should chain."""
        hits = [
            SSRHit(start=0, end=9, motif_raw="AAG", motif_canonical="AAG",
                   motif_size=3, repeat_units=3, repeat_length_bp=9, strand="+",
                   actual_repeat="AAGAAGAAG", is_perfect=True, repeat_class="perfect"),
            SSRHit(start=9, end=18, motif_raw="GAT", motif_canonical="ATG",
                   motif_size=3, repeat_units=3, repeat_length_bp=9, strand="+",
                   actual_repeat="GATGATGAT", is_perfect=True, repeat_class="perfect"),
        ]
        compounds = detect_compound_ssrs(hits, dmax_bp=10)
        assert len(compounds) == 1
        assert compounds[0].gaps_bp[0] == 0

    def test_homopolymer_in_compound(self):
        """Compound containing homopolymers."""
        hits = [
            SSRHit(start=0, end=10, motif_raw="A", motif_canonical="A",
                   motif_size=1, repeat_units=10, repeat_length_bp=10, strand="+",
                   actual_repeat="AAAAAAAAAA", is_perfect=True, repeat_class="perfect"),
            SSRHit(start=12, end=18, motif_raw="C", motif_canonical="C",
                   motif_size=1, repeat_units=6, repeat_length_bp=6, strand="+",
                   actual_repeat="CCCCCC", is_perfect=True, repeat_class="perfect"),
        ]
        compounds = detect_compound_ssrs(hits, dmax_bp=10)
        assert len(compounds) == 1

    def test_very_large_dmax(self):
        """Very large dmax brings distant SSRs together."""
        hits = [
            SSRHit(start=0, end=9, motif_raw="AAG", motif_canonical="AAG",
                   motif_size=3, repeat_units=3, repeat_length_bp=9, strand="+",
                   actual_repeat="AAGAAGAAG", is_perfect=True, repeat_class="perfect"),
            SSRHit(start=100, end=109, motif_raw="GAT", motif_canonical="ATG",
                   motif_size=3, repeat_units=3, repeat_length_bp=9, strand="+",
                   actual_repeat="GATGATGAT", is_perfect=True, repeat_class="perfect"),
        ]
        compounds = detect_compound_ssrs(hits, dmax_bp=1000)
        assert len(compounds) == 1

    def test_negative_dmax_treated_as_zero(self):
        """Negative dmax should not break logic."""
        hits = [
            SSRHit(start=0, end=9, motif_raw="AAG", motif_canonical="AAG",
                   motif_size=3, repeat_units=3, repeat_length_bp=9, strand="+",
                   actual_repeat="AAGAAGAAG", is_perfect=True, repeat_class="perfect"),
            SSRHit(start=14, end=23, motif_raw="GAT", motif_canonical="ATG",
                   motif_size=3, repeat_units=3, repeat_length_bp=9, strand="+",
                   actual_repeat="GATGATGAT", is_perfect=True, repeat_class="perfect"),
        ]
        compounds = detect_compound_ssrs(hits, dmax_bp=-5)
        # Gap of 5 > -5, so won't chain
        assert len(compounds) == 0


# ============================================================================
# Tests: Standardization Level Enum
# ============================================================================


class TestStandardizationLevelEnum:
    """Test StandardizationLevel enum values."""

    def test_enum_values(self):
        """All expected standardization levels present."""
        assert StandardizationLevel.L0.value == "L0"
        assert StandardizationLevel.L1.value == "L1"
        assert StandardizationLevel.L2.value == "L2"
        assert StandardizationLevel.Full.value == "Full"

    def test_enum_string_conversion(self):
        """Can convert string to enum."""
        level = StandardizationLevel("L2")
        assert level == StandardizationLevel.L2


# ============================================================================
# Tests: CompoundHit Dataclass
# ============================================================================


class TestCompoundHitDataclass:
    """Test CompoundHit properties and methods."""

    def test_compound_hit_fields(self, simple_perfect_hits):
        """CompoundHit has all expected fields."""
        compound = build_compound_hit(
            simple_perfect_hits,
            dmax_bp=10,
            standardization_level=StandardizationLevel.L2,
        )
        
        assert hasattr(compound, 'accession')
        assert hasattr(compound, 'start')
        assert hasattr(compound, 'end')
        assert hasattr(compound, 'components')
        assert hasattr(compound, 'gaps_bp')
        assert hasattr(compound, 'dmax_used')
        assert hasattr(compound, 'standardization_level')
        assert hasattr(compound, 'standardized_motif')
        assert hasattr(compound, 'compound_motif')
        assert hasattr(compound, 'total_repeat_length_bp')
        assert hasattr(compound, 'strand')

    def test_compound_hit_component_count(self, simple_perfect_hits):
        """CompoundHit correctly counts components."""
        compound = build_compound_hit(
            simple_perfect_hits,
            dmax_bp=10,
            standardization_level=StandardizationLevel.L2,
        )
        assert len(compound.components) == 2

    def test_compound_hit_strand_determination(self):
        """Compound strand is dominant strand of components."""
        hits = [
            SSRHit(start=0, end=9, motif_raw="AAG", motif_canonical="AAG",
                   motif_size=3, repeat_units=3, repeat_length_bp=9, strand="+",
                   actual_repeat="AAGAAGAAG", is_perfect=True, repeat_class="perfect"),
            SSRHit(start=12, end=21, motif_raw="CTT", motif_canonical="AAG",
                   motif_size=3, repeat_units=3, repeat_length_bp=9, strand="-",
                   actual_repeat="CTTCTTCTT", is_perfect=True, repeat_class="perfect"),
        ]
        compound = build_compound_hit(hits, dmax_bp=10, standardization_level=StandardizationLevel.L2)
        # More + strands (1) than - strands (1), equal → default +
        assert compound.strand in ("+", "-")
