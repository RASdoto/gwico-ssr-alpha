"""Validation tests for GWICO-SSR: gold-standard correctness and reproducibility.

Tests SSR detection against curated gold-standard sequences with known SSR
positions, verifying exact match of all detection attributes. Also tests
reproducibility across multiple runs.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from gwico_ssr.parsers.fasta_parser import parse_fasta
from gwico_ssr.ssr.detector import DetectionResult, SSRThresholds, detect_ssrs

GOLD_DIR = Path(__file__).parent.parent / "fixtures" / "gold_standard"
EXPECTED_PATH = GOLD_DIR / "gold_standard_expected.json"
FASTA_PATH = GOLD_DIR / "gold_standard.fasta"


@pytest.fixture
def gold_expected() -> dict:
    """Load gold-standard expected results."""
    return json.loads(EXPECTED_PATH.read_text())


@pytest.fixture
def gold_sequences() -> dict[str, str]:
    """Parse gold-standard FASTA into accession → sequence dict."""
    result = parse_fasta(FASTA_PATH)
    return {r.accession: r.sequence for r in result.records}


# ===========================================================================
# GOLD-STANDARD CORRECTNESS TESTS
# ===========================================================================


class TestGoldStandardCorrectness:
    """Verify SSR detection exactly matches gold-standard expected results."""

    def test_gold_standard_files_exist(self):
        """Ensure gold-standard fixtures are present."""
        assert FASTA_PATH.exists(), f"Missing {FASTA_PATH}"
        assert EXPECTED_PATH.exists(), f"Missing {EXPECTED_PATH}"

    def test_gold_001_hit_count(self, gold_sequences, gold_expected):
        """GOLD_001: correct number of SSR hits detected."""
        seq = gold_sequences["GOLD_001"]
        thresholds = SSRThresholds(**gold_expected["thresholds"])
        result = detect_ssrs(seq, "GOLD_001", thresholds)
        expected = gold_expected["sequences"]["GOLD_001"]
        assert result.hit_count == expected["total_hits"]

    def test_gold_001_exact_positions(self, gold_sequences, gold_expected):
        """GOLD_001: each hit has exact start/end positions."""
        seq = gold_sequences["GOLD_001"]
        thresholds = SSRThresholds(**gold_expected["thresholds"])
        result = detect_ssrs(seq, "GOLD_001", thresholds)
        expected_hits = gold_expected["sequences"]["GOLD_001"]["hits"]

        for i, (hit, exp) in enumerate(zip(result.hits, expected_hits)):
            assert hit.start == exp["start"], (
                f"Hit {i}: start {hit.start} != expected {exp['start']}"
            )
            assert hit.end == exp["end"], (
                f"Hit {i}: end {hit.end} != expected {exp['end']}"
            )

    def test_gold_001_exact_motifs(self, gold_sequences, gold_expected):
        """GOLD_001: each hit has correct motif and canonical form."""
        seq = gold_sequences["GOLD_001"]
        thresholds = SSRThresholds(**gold_expected["thresholds"])
        result = detect_ssrs(seq, "GOLD_001", thresholds)
        expected_hits = gold_expected["sequences"]["GOLD_001"]["hits"]

        for i, (hit, exp) in enumerate(zip(result.hits, expected_hits)):
            assert hit.motif_raw == exp["motif_raw"], (
                f"Hit {i}: motif_raw {hit.motif_raw} != {exp['motif_raw']}"
            )
            assert hit.motif_canonical == exp["motif_canonical"], (
                f"Hit {i}: canonical {hit.motif_canonical} != {exp['motif_canonical']}"
            )
            assert hit.motif_size == exp["motif_size"], (
                f"Hit {i}: size {hit.motif_size} != {exp['motif_size']}"
            )

    def test_gold_001_repeat_attributes(self, gold_sequences, gold_expected):
        """GOLD_001: repeat units, length, strand, actual_repeat."""
        seq = gold_sequences["GOLD_001"]
        thresholds = SSRThresholds(**gold_expected["thresholds"])
        result = detect_ssrs(seq, "GOLD_001", thresholds)
        expected_hits = gold_expected["sequences"]["GOLD_001"]["hits"]

        for i, (hit, exp) in enumerate(zip(result.hits, expected_hits)):
            assert hit.repeat_units == exp["repeat_units"], (
                f"Hit {i}: units {hit.repeat_units} != {exp['repeat_units']}"
            )
            assert hit.repeat_length_bp == exp["repeat_length_bp"], (
                f"Hit {i}: length {hit.repeat_length_bp} != {exp['repeat_length_bp']}"
            )
            assert hit.strand == exp["strand"], (
                f"Hit {i}: strand {hit.strand} != {exp['strand']}"
            )
            assert hit.actual_repeat == exp["actual_repeat"], (
                f"Hit {i}: actual {hit.actual_repeat} != {exp['actual_repeat']}"
            )

    def test_gold_001_all_motif_sizes_covered(self, gold_sequences, gold_expected):
        """GOLD_001: SSRs span all motif sizes 1-6."""
        seq = gold_sequences["GOLD_001"]
        thresholds = SSRThresholds(**gold_expected["thresholds"])
        result = detect_ssrs(seq, "GOLD_001", thresholds)
        sizes = {h.motif_size for h in result.hits}
        assert sizes == {1, 2, 3, 4, 5, 6}, f"Missing sizes: {set(range(1,7)) - sizes}"

    def test_gold_002_hit_count(self, gold_sequences, gold_expected):
        """GOLD_002: correct number of SSR hits detected."""
        seq = gold_sequences["GOLD_002"]
        thresholds = SSRThresholds(**gold_expected["thresholds"])
        result = detect_ssrs(seq, "GOLD_002", thresholds)
        expected = gold_expected["sequences"]["GOLD_002"]
        assert result.hit_count == expected["total_hits"]

    def test_gold_002_exact_positions(self, gold_sequences, gold_expected):
        """GOLD_002: each hit has exact start/end positions."""
        seq = gold_sequences["GOLD_002"]
        thresholds = SSRThresholds(**gold_expected["thresholds"])
        result = detect_ssrs(seq, "GOLD_002", thresholds)
        expected_hits = gold_expected["sequences"]["GOLD_002"]["hits"]

        for i, (hit, exp) in enumerate(zip(result.hits, expected_hits)):
            assert hit.start == exp["start"], (
                f"Hit {i}: start {hit.start} != expected {exp['start']}"
            )
            assert hit.end == exp["end"], (
                f"Hit {i}: end {hit.end} != expected {exp['end']}"
            )

    def test_gold_002_exact_motifs(self, gold_sequences, gold_expected):
        """GOLD_002: each hit has correct canonical motif and strand."""
        seq = gold_sequences["GOLD_002"]
        thresholds = SSRThresholds(**gold_expected["thresholds"])
        result = detect_ssrs(seq, "GOLD_002", thresholds)
        expected_hits = gold_expected["sequences"]["GOLD_002"]["hits"]

        for i, (hit, exp) in enumerate(zip(result.hits, expected_hits)):
            assert hit.motif_canonical == exp["motif_canonical"]
            assert hit.strand == exp["strand"]
            assert hit.actual_repeat == exp["actual_repeat"]

    def test_gold_002_repeat_attributes(self, gold_sequences, gold_expected):
        """GOLD_002: repeat units and lengths match."""
        seq = gold_sequences["GOLD_002"]
        thresholds = SSRThresholds(**gold_expected["thresholds"])
        result = detect_ssrs(seq, "GOLD_002", thresholds)
        expected_hits = gold_expected["sequences"]["GOLD_002"]["hits"]

        for i, (hit, exp) in enumerate(zip(result.hits, expected_hits)):
            assert hit.repeat_units == exp["repeat_units"]
            assert hit.repeat_length_bp == exp["repeat_length_bp"]

    def test_sequence_lengths_match(self, gold_sequences, gold_expected):
        """Parsed sequences match expected lengths."""
        for acc, spec in gold_expected["sequences"].items():
            assert len(gold_sequences[acc]) == spec["length"], (
                f"{acc}: length {len(gold_sequences[acc])} != {spec['length']}"
            )

    def test_no_hits_on_empty_sequence(self):
        """Empty sequence produces zero hits."""
        result = detect_ssrs("", "EMPTY", SSRThresholds())
        assert result.hit_count == 0

    def test_no_hits_below_threshold(self):
        """Sequence with repeats below threshold produces zero hits."""
        # 9 A's (threshold is 10 for mono)
        result = detect_ssrs("A" * 9 + "GCGCGC", "SHORT", SSRThresholds())
        # Just 9 A's + non-repeating => 0 mono hits, and the GC only 3 di-repeats (threshold 5)
        assert result.hit_count == 0


# ===========================================================================
# REPRODUCIBILITY TESTS
# ===========================================================================


class TestReproducibility:
    """Verify detection is deterministic across repeated runs."""

    def test_repeated_runs_identical_gold_001(self, gold_sequences, gold_expected):
        """Running detection 5 times on GOLD_001 yields identical results."""
        seq = gold_sequences["GOLD_001"]
        thresholds = SSRThresholds(**gold_expected["thresholds"])

        runs = [detect_ssrs(seq, "GOLD_001", thresholds) for _ in range(5)]

        first = runs[0]
        for run in runs[1:]:
            assert run.hit_count == first.hit_count
            for h1, h2 in zip(first.hits, run.hits):
                assert h1.start == h2.start
                assert h1.end == h2.end
                assert h1.motif_raw == h2.motif_raw
                assert h1.motif_canonical == h2.motif_canonical
                assert h1.repeat_units == h2.repeat_units
                assert h1.strand == h2.strand

    def test_repeated_runs_identical_gold_002(self, gold_sequences, gold_expected):
        """Running detection 5 times on GOLD_002 yields identical results."""
        seq = gold_sequences["GOLD_002"]
        thresholds = SSRThresholds(**gold_expected["thresholds"])

        runs = [detect_ssrs(seq, "GOLD_002", thresholds) for _ in range(5)]

        first = runs[0]
        for run in runs[1:]:
            assert run.hit_count == first.hit_count
            for h1, h2 in zip(first.hits, run.hits):
                assert h1.start == h2.start
                assert h1.end == h2.end
                assert h1.motif_canonical == h2.motif_canonical

    def test_different_accession_same_sequence(self, gold_sequences, gold_expected):
        """Same sequence with different accession IDs produces same hits."""
        seq = gold_sequences["GOLD_001"]
        thresholds = SSRThresholds(**gold_expected["thresholds"])

        r1 = detect_ssrs(seq, "ACC_A", thresholds)
        r2 = detect_ssrs(seq, "ACC_B", thresholds)

        assert r1.hit_count == r2.hit_count
        for h1, h2 in zip(r1.hits, r2.hits):
            assert h1.start == h2.start
            assert h1.end == h2.end
            assert h1.motif_canonical == h2.motif_canonical

    def test_synthetic_genome_reproducibility(self):
        """Synthetic genome detection is reproducible across 3 runs."""
        from tests.benchmarks.benchmark_harness import (
            check_detection_reproducibility,
        )

        result = check_detection_reproducibility(
            genome_size=10000, num_runs=3, seed=42
        )
        assert result["all_identical"], "Detection was not reproducible!"
        assert result["hits_per_run"] > 0, "Expected at least some SSRs"


# ===========================================================================
# THRESHOLD SENSITIVITY TESTS
# ===========================================================================


class TestThresholdSensitivity:
    """Verify detection respects threshold changes correctly."""

    def test_lower_thresholds_more_hits(self, gold_sequences):
        """Lower thresholds should detect equal or more SSRs."""
        seq = gold_sequences["GOLD_001"]
        default = SSRThresholds()
        lower = SSRThresholds(mono=5, di=3, tri=2, tetra=2, penta=2, hexa=2)

        r_default = detect_ssrs(seq, "GOLD_001", default)
        r_lower = detect_ssrs(seq, "GOLD_001", lower)

        assert r_lower.hit_count >= r_default.hit_count

    def test_higher_thresholds_fewer_hits(self, gold_sequences):
        """Higher thresholds should detect equal or fewer SSRs."""
        seq = gold_sequences["GOLD_001"]
        default = SSRThresholds()
        higher = SSRThresholds(mono=20, di=10, tri=8, tetra=6, penta=5, hexa=4)

        r_default = detect_ssrs(seq, "GOLD_001", default)
        r_higher = detect_ssrs(seq, "GOLD_001", higher)

        assert r_higher.hit_count <= r_default.hit_count

    def test_extreme_thresholds_zero_hits(self, gold_sequences):
        """Extremely high thresholds should yield zero hits on short sequences."""
        seq = gold_sequences["GOLD_001"]
        extreme = SSRThresholds(mono=100, di=100, tri=100, tetra=100, penta=100, hexa=100)
        result = detect_ssrs(seq, "GOLD_001", extreme)
        assert result.hit_count == 0


# ===========================================================================
# COORDINATE INVARIANT TESTS
# ===========================================================================


class TestCoordinateInvariants:
    """Verify structural invariants of SSR detection results."""

    def test_half_open_coordinates(self, gold_sequences, gold_expected):
        """All hits use 0-based half-open coordinates: end - start == repeat_length_bp."""
        seq = gold_sequences["GOLD_001"]
        thresholds = SSRThresholds(**gold_expected["thresholds"])
        result = detect_ssrs(seq, "GOLD_001", thresholds)

        for h in result.hits:
            assert h.end - h.start == h.repeat_length_bp, (
                f"[{h.start}:{h.end}) length {h.end - h.start} "
                f"!= repeat_length_bp {h.repeat_length_bp}"
            )

    def test_repeat_length_equals_motif_times_units(self, gold_sequences, gold_expected):
        """repeat_length_bp == motif_size * repeat_units for all hits."""
        for acc in gold_expected["sequences"]:
            seq = gold_sequences[acc]
            thresholds = SSRThresholds(**gold_expected["thresholds"])
            result = detect_ssrs(seq, acc, thresholds)
            for h in result.hits:
                assert h.repeat_length_bp == h.motif_size * h.repeat_units

    def test_actual_repeat_matches_sequence(self, gold_sequences, gold_expected):
        """actual_repeat field matches the sequence at [start:end)."""
        for acc in gold_expected["sequences"]:
            seq = gold_sequences[acc]
            thresholds = SSRThresholds(**gold_expected["thresholds"])
            result = detect_ssrs(seq, acc, thresholds)
            for h in result.hits:
                assert h.actual_repeat == seq[h.start : h.end].upper()

    def test_hits_within_sequence_bounds(self, gold_sequences, gold_expected):
        """All hit coordinates are within [0, sequence_length)."""
        for acc in gold_expected["sequences"]:
            seq = gold_sequences[acc]
            thresholds = SSRThresholds(**gold_expected["thresholds"])
            result = detect_ssrs(seq, acc, thresholds)
            for h in result.hits:
                assert h.start >= 0
                assert h.end <= len(seq)
                assert h.start < h.end

    def test_hits_sorted_by_position(self, gold_sequences, gold_expected):
        """Hits are returned sorted by start position."""
        for acc in gold_expected["sequences"]:
            seq = gold_sequences[acc]
            thresholds = SSRThresholds(**gold_expected["thresholds"])
            result = detect_ssrs(seq, acc, thresholds)
            starts = [h.start for h in result.hits]
            assert starts == sorted(starts)

    def test_motif_size_range(self, gold_sequences, gold_expected):
        """All motif sizes are in range 1-6."""
        for acc in gold_expected["sequences"]:
            seq = gold_sequences[acc]
            thresholds = SSRThresholds(**gold_expected["thresholds"])
            result = detect_ssrs(seq, acc, thresholds)
            for h in result.hits:
                assert 1 <= h.motif_size <= 6

    def test_strand_valid(self, gold_sequences, gold_expected):
        """All strands are '+' or '-'."""
        for acc in gold_expected["sequences"]:
            seq = gold_sequences[acc]
            thresholds = SSRThresholds(**gold_expected["thresholds"])
            result = detect_ssrs(seq, acc, thresholds)
            for h in result.hits:
                assert h.strand in ("+", "-")

    def test_repeat_units_meet_threshold(self, gold_sequences, gold_expected):
        """All detected SSRs have repeat_units >= threshold for their motif_size."""
        thresholds = SSRThresholds(**gold_expected["thresholds"])
        for acc in gold_expected["sequences"]:
            seq = gold_sequences[acc]
            result = detect_ssrs(seq, acc, thresholds)
            for h in result.hits:
                min_req = thresholds.for_size(h.motif_size)
                assert h.repeat_units >= min_req, (
                    f"Hit at [{h.start}:{h.end}) has {h.repeat_units} units "
                    f"but threshold requires {min_req}"
                )
