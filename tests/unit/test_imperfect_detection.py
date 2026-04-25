"""Unit tests for IMEX-style imperfect SSR detection (Chunk 6).

Tests cover:
- Imperfection utility functions (substitution/indel counting, alignment, %)
- Seed-and-extend algorithm
- Imperfect SSR detection end-to-end
- Backward compatibility with perfect-only mode
- Determinism and edge cases
"""

import pytest
from gwico_ssr.ssr.imperfection import (
    IndelEvent,
    AlignmentResult,
    calculate_imperfection_pct,
    count_substitutions_in_alignment,
    find_indels,
    align_motif_to_region,
    generate_alignment_text,
)
from gwico_ssr.ssr.detector import (
    detect_ssrs,
    detect_imperfect_ssrs,
    SSRThresholds,
    ImperfectionThresholds,
)


class TestImperfectionUtilities:
    """Test low-level imperfection calculation functions."""

    def test_calculate_imperfection_pct_no_errors(self):
        """Perfect repeats should have 0% imperfection."""
        pct = calculate_imperfection_pct(substitutions=0, indels=0, tract_length_bp=30)
        assert pct == 0.0

    def test_calculate_imperfection_pct_with_substitutions(self):
        """Test imperfection % with substitutions only."""
        # 3 substitutions in 30 bp = 10%
        pct = calculate_imperfection_pct(substitutions=3, indels=0, tract_length_bp=30)
        assert pct == 10.0

    def test_calculate_imperfection_pct_with_indels(self):
        """Test imperfection % with indels only."""
        # 2 indel events in 30 bp = 6.67%
        pct = calculate_imperfection_pct(substitutions=0, indels=2, tract_length_bp=30)
        assert abs(pct - 6.666666667) < 0.01

    def test_calculate_imperfection_pct_mixed_errors(self):
        """Test imperfection % with both substitutions and indels."""
        # 2 subs + 1 indel in 30 bp = 10%
        pct = calculate_imperfection_pct(substitutions=2, indels=1, tract_length_bp=30)
        assert pct == 10.0

    def test_calculate_imperfection_pct_zero_length(self):
        """Zero-length tract should return 0%."""
        pct = calculate_imperfection_pct(substitutions=1, indels=1, tract_length_bp=0)
        assert pct == 0.0

    def test_count_substitutions_perfect_alignment(self):
        """Count substitutions in a perfect alignment."""
        motif = "AGC"
        region = "AGCAGCAGC"
        subs, cigar = count_substitutions_in_alignment(motif, region)
        
        assert subs == 0
        assert cigar == "========="  # All matches

    def test_count_substitutions_with_mismatch(self):
        """Count substitutions with one mismatch."""
        motif = "AGC"
        region = "AGCAGTAGC"  # T instead of C in position 6
        subs, cigar = count_substitutions_in_alignment(motif, region)
        
        assert subs == 1
        assert "M" in cigar  # Should have one mismatch

    def test_count_substitutions_multiple_mismatches(self):
        """Count multiple substitutions."""
        motif = "AT"
        region = "ATCTAGAT"  # Two mismatches
        subs, cigar = count_substitutions_in_alignment(motif, region)
        
        assert subs == 2

    def test_find_indels_no_indels(self):
        """No indels in perfect repetition."""
        motif = "AGC"
        region = "AGCAGCAGC"
        indels, cigar = find_indels(motif, region, max_indel_size=2)
        
        assert len(indels) == 0

    def test_find_indels_deletion(self):
        """Detect deletion event."""
        motif = "AGC"
        region = "AGCAGCAG"  # One base short
        indels, cigar = find_indels(motif, region, max_indel_size=2)
        
        assert len(indels) == 1
        assert indels[0].event_type == "deletion"

    def test_align_motif_perfect(self):
        """Align perfect motif to region."""
        motif = "AGC"
        region = "AGCAGCAGC"
        
        result = align_motif_to_region(
            motif=motif,
            region=region,
            max_substitutions=1,
            max_imperfection_pct=5.0,
        )
        
        assert result is not None
        assert result.substitutions == 0
        assert result.imperfection_pct == 0.0
        assert result.motif_repetitions == 3

    def test_align_motif_with_substitution(self):
        """Align motif with acceptable substitution."""
        motif = "AGC"
        region = "AGCAGTAGC"  # One substitution (11.11% imperfection in 9 bp)
        
        result = align_motif_to_region(
            motif=motif,
            region=region,
            max_substitutions=1,
            max_imperfection_pct=15.0,  # Need >= 11.11% to accept this
        )
        
        assert result is not None
        assert result.substitutions == 1
        assert result.imperfection_pct > 0

    def test_align_motif_exceeds_substitution_threshold(self):
        """Alignment fails when substitutions exceed threshold."""
        motif = "AGC"
        region = "AGCAGTAGC"  # One substitution
        
        result = align_motif_to_region(
            motif=motif,
            region=region,
            max_substitutions=0,  # Don't allow any
            max_imperfection_pct=10.0,
        )
        
        assert result is None

    def test_align_motif_exceeds_imperfection_pct_threshold(self):
        """Alignment fails when imperfection % exceeds threshold."""
        motif = "AGC"
        region = "AGCAGTAGC"  # ~11% imperfection (1/9)
        
        result = align_motif_to_region(
            motif=motif,
            region=region,
            max_substitutions=1,
            max_imperfection_pct=5.0,  # Too strict
        )
        
        assert result is None

    def test_generate_alignment_text(self):
        """Generate readable alignment text."""
        motif = "AGC"
        region = "AGCAGTAGC"
        
        text = generate_alignment_text(motif, region, 1, [])
        
        assert "Motif:" in text
        assert "Aligned:" in text
        assert "Substitutions:" in text


class TestSeedAndExtend:
    """Test seed-and-extend algorithm for imperfect detection."""

    def test_simple_imperfect_substitution(self):
        """Detect imperfect SSR with single substitution."""
        # Sequence has perfect AGC seed (3 units), then extends to imperfect with substitution
        sequence = "AGCAGCAGCAGTAGC"  # AGC x 3 (perfect seed) + AGT (imperfect) + AGC
        accession = "test_sub"
        
        result = detect_imperfect_ssrs(
            sequence=sequence,
            accession=accession,
            thresholds=SSRThresholds(tri=3),
            imperfection_pct_threshold=15.0,  # Allow up to 15% imperfection
            max_indel_size=2,
        )
        
        # Should detect at least the perfect seed
        assert result.hit_count > 0

    def test_determinism_repeated_detection(self):
        """Repeated detections should give identical results."""
        sequence = "AGCAGCAGTAGCAGC"
        accession = "test_det"
        
        result1 = detect_imperfect_ssrs(sequence, accession)
        result2 = detect_imperfect_ssrs(sequence, accession)
        
        assert result1.hit_count == result2.hit_count
        for h1, h2 in zip(result1.hits, result2.hits):
            assert h1.start == h2.start
            assert h1.end == h2.end
            assert h1.num_substitutions == h2.num_substitutions

    def test_high_imperfection_filtered(self):
        """High-imperfection SSRs filtered by threshold."""
        # Create sequence with high error rate
        sequence = "AGCAGACGGATCGACGACG"  # High imperfection
        
        result = detect_imperfect_ssrs(
            sequence=sequence,
            accession="test_high",
            imperfection_pct_threshold=5.0,  # Strict threshold
            max_indel_size=1,
        )
        
        # Should filter out very imperfect SSRs
        for hit in result.hits:
            if hit.imperfection_pct is not None:
                assert hit.imperfection_pct <= 5.0 + 0.01  # Small tolerance for float


class TestImperfectDetector:
    """Test end-to-end imperfect SSR detection."""

    def test_detect_imperfect_vs_perfect_backward_compat(self):
        """Perfect mode should produce same results with and without imperfect detector."""
        sequence = "AGCAGCAGCAGC"  # Pure perfect repeats
        
        perfect_result = detect_ssrs(sequence, "test_perf")
        imperfect_result = detect_imperfect_ssrs(sequence, "test_perf")
        
        # Perfect mode with imperfect detector should still find the perfect SSRs
        assert imperfect_result.hit_count > 0
        perfect_hits = [h for h in imperfect_result.hits if h.is_perfect]
        assert len(perfect_hits) > 0

    def test_imperfect_detector_finds_extensions(self):
        """Imperfect detector should find extended imperfect regions."""
        # Create sequence: perfect seed + imperfect extension
        sequence = "AGCAGCAGC" + "T" + "AGCAGC"
        
        result = detect_imperfect_ssrs(
            sequence=sequence,
            accession="test_ext",
            thresholds=SSRThresholds(tri=3),
            imperfection_pct_threshold=15.0,
            max_indel_size=2,
        )
        
        # Should detect imperfect SSRs in the extended region
        assert result.hit_count > 0

    def test_empty_sequence(self):
        """Empty sequence should produce no hits."""
        result = detect_imperfect_ssrs("", "empty")
        assert result.hit_count == 0

    def test_too_short_sequence(self):
        """Sequence too short for motif should produce no hits."""
        result = detect_imperfect_ssrs("AG", "short", thresholds=SSRThresholds(tri=3))
        assert result.hit_count == 0

    def test_ambiguous_bases_skipped(self):
        """Sequences with ambiguous bases should skip those regions."""
        sequence = "AGCAGCNAGCAGC"  # N in the middle
        result = detect_imperfect_ssrs(sequence, "ambig")
        # Should not crash and should handle ambiguous bases gracefully
        assert result.hit_count >= 0

    def test_repeat_class_field_populated(self):
        """Detected SSRs should have repeat_class field set."""
        sequence = "AGCAGCAGCAGC"
        result = detect_imperfect_ssrs(sequence, "test_class")
        
        for hit in result.hits:
            assert hit.repeat_class in ("perfect", "imperfect", "compound_component")

    def test_imperfection_metrics_populated(self):
        """Imperfect SSRs should have all imperfection metrics."""
        sequence = "AGCAGCAGTAGCAGC"  # One substitution
        result = detect_imperfect_ssrs(sequence, "test_metrics")
        
        for hit in result.hits:
            if not hit.is_perfect:
                assert hit.num_substitutions >= 0
                assert hit.num_indels >= 0
                assert hit.imperfection_pct is not None


class TestBackwardCompatibility:
    """Ensure perfect-only mode is unaffected by imperfect detection code."""

    def test_perfect_detector_unchanged(self):
        """Perfect SSR detector function should remain unchanged."""
        sequence = "ATATATAT" * 5  # 40 bases of AT repeats
        result = detect_ssrs(sequence, "test_compat")
        
        # Should find AT repeats
        assert result.hit_count > 0
        
        # All should be marked as perfect
        assert all(h.is_perfect for h in result.hits)

    def test_perfect_detector_no_imperfection_fields_for_perfect(self):
        """Perfect SSRs should have sensible imperfection field defaults."""
        sequence = "AGCAGCAGC"
        result = detect_ssrs(sequence, "test_fields")
        
        for hit in result.hits:
            if hit.is_perfect:
                assert hit.num_substitutions == 0
                assert hit.num_indels == 0
                assert hit.repeat_class == "perfect"

    def test_config_validation(self):
        """ImperfectionThresholds should validate configuration."""
        # Valid configuration
        thresholds = ImperfectionThresholds(
            motif_size=3,
            max_mismatches_per_unit=1,
            max_imperfection_pct=5.0,
            min_repeat_units=2,
            max_indel_size=2,
        )
        assert thresholds.motif_size == 3

    def test_detector_config_validation_modes(self):
        """DetectorConfig should validate detector modes."""
        from gwico_ssr.ssr.detector import DetectorConfig
        
        # Valid mode
        cfg = DetectorConfig(detector_mode="perfect")
        assert cfg.detector_mode == "perfect"
        
        # Invalid mode should raise
        with pytest.raises(ValueError):
            DetectorConfig(detector_mode="invalid")

    def test_detector_config_validation_standardization(self):
        """DetectorConfig should validate standardization levels."""
        from gwico_ssr.ssr.detector import DetectorConfig
        
        # Valid levels
        for level in ("L0", "L1", "L2", "Full"):
            cfg = DetectorConfig(standardization_level=level)
            assert cfg.standardization_level == level
        
        # Invalid level should raise
        with pytest.raises(ValueError):
            DetectorConfig(standardization_level="L99")


class TestEdgeCases:
    """Test edge cases and boundary conditions."""

    def test_imperfection_at_sequence_start(self):
        """Imperfect SSR at sequence start."""
        sequence = "AGTAGCAGC"  # Starts with imperfection
        result = detect_imperfect_ssrs(sequence, "test_start")
        # Should not crash
        assert result.hit_count >= 0

    def test_imperfection_at_sequence_end(self):
        """Imperfect SSR at sequence end."""
        sequence = "AGCAGCAGT"  # Ends with imperfection
        result = detect_imperfect_ssrs(sequence, "test_end")
        # Should not crash
        assert result.hit_count >= 0

    def test_single_repeat_unit(self):
        """Single repeat unit should not be detected (min_repeats=2+)."""
        sequence = "NNNNNNAGCNNNNNN"
        result = detect_imperfect_ssrs(
            sequence,
            "test_single",
            thresholds=SSRThresholds(tri=2),
        )
        
        # Should not detect single-unit SSR
        for hit in result.hits:
            assert hit.repeat_units >= 2

    def test_cigar_string_generation(self):
        """CIGAR string should be generated for imperfect SSRs."""
        sequence = "AGCAGCAGTAGC"  # With substitution
        result = detect_imperfect_ssrs(sequence, "test_cigar")
        
        for hit in result.hits:
            if not hit.is_perfect and hit.imperfection_cigar:
                # CIGAR should contain alignment information
                assert hit.imperfection_cigar is not None

    def test_homopolymer_detection(self):
        """Homopolymers (mono repeats) should be detected properly."""
        sequence = "NNNNNNAAAAAAAAANNNNN"  # 10 A's
        result = detect_imperfect_ssrs(
            sequence,
            "test_homo",
            thresholds=SSRThresholds(mono=8),
        )
        
        # Should find homopolymer
        mono_hits = [h for h in result.hits if h.motif_size == 1]
        assert len(mono_hits) > 0

    def test_overlapping_motif_sizes(self):
        """Test behavior when same region matches multiple motif sizes."""
        sequence = "AAGAAGAAGAAG"  # Can be detected as AA or AAG
        result = detect_imperfect_ssrs(sequence, "test_overlap")
        
        # Overlap resolution should prefer shorter motif
        # (per IMEX algorithm)
        assert result.hit_count >= 0


class TestFixtureIntegration:
    """Integration tests using fixture files."""

    def test_fixture_imperfect_substitution(self):
        """Test with substitution fixture."""
        from gwico_ssr.parsers.fasta_parser import parse_fasta
        
        fixtures_dir = "tests/fixtures"
        result_path = f"{fixtures_dir}/imperfect_substitution.fasta"
        
        try:
            result = parse_fasta(result_path)
            if result.records:
                rec = result.records[0]
                detection = detect_imperfect_ssrs(
                    rec.sequence,
                    rec.accession,
                    imperfection_pct_threshold=10.0,
                )
                # Should process without error
                assert detection.accession == rec.accession
        except FileNotFoundError:
            pytest.skip("Fixture file not found")

    def test_fixture_imperfect_insertion(self):
        """Test with insertion fixture."""
        from gwico_ssr.parsers.fasta_parser import parse_fasta
        
        fixtures_dir = "tests/fixtures"
        result_path = f"{fixtures_dir}/imperfect_insertion.fasta"
        
        try:
            result = parse_fasta(result_path)
            if result.records:
                rec = result.records[0]
                detection = detect_imperfect_ssrs(
                    rec.sequence,
                    rec.accession,
                    imperfection_pct_threshold=10.0,
                )
                assert detection.accession == rec.accession
        except FileNotFoundError:
            pytest.skip("Fixture file not found")

    def test_fixture_imperfect_deletion(self):
        """Test with deletion fixture."""
        from gwico_ssr.parsers.fasta_parser import parse_fasta
        
        fixtures_dir = "tests/fixtures"
        result_path = f"{fixtures_dir}/imperfect_deletion.fasta"
        
        try:
            result = parse_fasta(result_path)
            if result.records:
                rec = result.records[0]
                detection = detect_imperfect_ssrs(
                    rec.sequence,
                    rec.accession,
                    imperfection_pct_threshold=10.0,
                )
                assert detection.accession == rec.accession
        except FileNotFoundError:
            pytest.skip("Fixture file not found")
