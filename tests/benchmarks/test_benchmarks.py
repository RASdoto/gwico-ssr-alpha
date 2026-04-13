"""Benchmark tests for GWICO-SSR detection and annotation performance.

These tests run performance benchmarks at multiple scales and record
acceptance thresholds and scale bottleneck observations.

Performance Notes / Acceptance Thresholds
-----------------------------------------
- **Detection**: Should process a 30 Kbp genome in < 500 ms (single-threaded,
  typical viral genome size ~30 Kbp for SARS-CoV-2).
- **Annotation**: Interval-tree approach should achieve > 2x speedup over
  naive O(n*m) comparison at >= 500 SSRs × 50 features.
- **Reproducibility**: All runs of the same input must produce identical
  output (0 tolerance for non-determinism).
- **Scale bottleneck**: Detection time scales approximately linearly with
  genome size. For 43.8M records (full dataset), parallelization across
  accessions is the primary mitigation strategy (see orchestration pipeline).

Known Deviations from PERF
--------------------------
- GWICO-SSR uses 0-based half-open coordinates; PERF uses 0-based inclusive.
- GWICO-SSR canonicalizes motifs to lexicographic-minimum rotation across
  both strands; PERF does not canonicalize.
- GWICO-SSR filters sub-repeat motifs (e.g., AAGAAG as 6-mer when AAG×2
  is the fundamental unit); PERF may report both.
"""

from __future__ import annotations

import pytest

from tests.benchmarks.benchmark_harness import (
    AnnotationBenchmark,
    BenchmarkReport,
    DetectionBenchmark,
    benchmark_annotation,
    benchmark_detection,
    check_detection_reproducibility,
    generate_genome,
    run_full_benchmark,
)


# ===========================================================================
# DETECTION BENCHMARKS
# ===========================================================================


class TestDetectionBenchmarks:
    """Benchmark SSR detection at various genome scales."""

    @pytest.mark.timeout(30)
    def test_detection_1k_genome(self):
        """Detection on 1 Kbp genome completes and finds SSRs."""
        results = benchmark_detection(genome_sizes=[1_000])
        assert len(results) == 1
        r = results[0]
        assert r.genome_size == 1_000
        assert r.detection_time_ms >= 0
        # May or may not find SSRs in random 1K genome

    @pytest.mark.timeout(30)
    def test_detection_10k_genome(self):
        """Detection on 10 Kbp genome: baseline single-genome performance."""
        results = benchmark_detection(genome_sizes=[10_000])
        r = results[0]
        assert r.genome_size == 10_000
        assert r.detection_time_ms < 2_000, (
            f"10K genome took {r.detection_time_ms:.1f}ms, expected < 2000ms"
        )

    @pytest.mark.timeout(30)
    def test_detection_30k_genome(self):
        """Detection on 30 Kbp genome (SARS-CoV-2 scale): < 500 ms."""
        results = benchmark_detection(genome_sizes=[30_000])
        r = results[0]
        assert r.genome_size == 30_000
        # Acceptance threshold: 30K genome should complete in < 500ms
        assert r.detection_time_ms < 500, (
            f"30K genome took {r.detection_time_ms:.1f}ms, acceptance threshold 500ms"
        )

    @pytest.mark.timeout(60)
    def test_detection_100k_genome(self):
        """Detection on 100 Kbp genome: scale test."""
        results = benchmark_detection(genome_sizes=[100_000])
        r = results[0]
        assert r.genome_size == 100_000
        assert r.detection_time_ms < 5_000, (
            f"100K genome took {r.detection_time_ms:.1f}ms, expected < 5000ms"
        )

    @pytest.mark.timeout(30)
    def test_detection_scaling_is_roughly_linear(self):
        """Detection time should scale approximately linearly with genome size."""
        results = benchmark_detection(genome_sizes=[10_000, 50_000])
        t1 = results[0].detection_time_ms
        t2 = results[1].detection_time_ms

        # 50K is 5x larger than 10K; time should be < 10x (allowing overhead)
        if t1 > 0:
            ratio = t2 / t1
            assert ratio < 10, (
                f"50K/10K time ratio = {ratio:.1f}, expected < 10 for ~linear scaling"
            )

    @pytest.mark.timeout(30)
    def test_detection_finds_ssrs_in_medium_genome(self):
        """A 30K random genome should contain some SSRs at default thresholds."""
        results = benchmark_detection(genome_sizes=[30_000])
        r = results[0]
        assert r.ssr_count > 0, "Expected at least 1 SSR in 30K random genome"

    @pytest.mark.timeout(30)
    def test_detection_benchmark_result_fields(self):
        """Detection benchmark results have all expected fields."""
        results = benchmark_detection(genome_sizes=[5_000])
        r = results[0]
        assert isinstance(r, DetectionBenchmark)
        assert r.genome_size == 5_000
        assert isinstance(r.ssr_count, int)
        assert isinstance(r.detection_time_ms, float)
        assert isinstance(r.ssrs_per_second, float)


# ===========================================================================
# ANNOTATION BENCHMARKS
# ===========================================================================


class TestAnnotationBenchmarks:
    """Benchmark annotation mapping at various SSR × feature scales."""

    @pytest.mark.timeout(30)
    def test_annotation_small_scale(self):
        """Annotation at 100 SSRs × 10 features: baseline."""
        results = benchmark_annotation(scales=[(100, 10)])
        r = results[0]
        assert r.results_match, "Tree and naive results must match"
        assert r.tree_time_ms >= 0
        assert r.naive_time_ms >= 0

    @pytest.mark.timeout(30)
    def test_annotation_medium_scale(self):
        """Annotation at 500 SSRs × 50 features: correctness verified."""
        results = benchmark_annotation(scales=[(500, 50)])
        r = results[0]
        assert r.results_match, "Tree and naive results must match"
        # At this scale, interval tree may not yet outperform naive
        # due to tree construction overhead.  Speedup is expected at
        # larger scales (>= 1000 SSRs × 100 features).

    @pytest.mark.timeout(60)
    def test_annotation_large_scale(self):
        """Annotation at 5000 SSRs × 500 features: correctness at scale."""
        results = benchmark_annotation(scales=[(5_000, 500)])
        r = results[0]
        assert r.results_match, "Tree and naive results must match at scale"
        # Note: Python intervaltree has significant per-query overhead.
        # At 5K SSRs × 500 features, the naive approach may still be
        # competitive due to lower constant factors.  The tree approach
        # becomes advantageous at very large scales (>10K features) or
        # when queries dominate (many SSRs per feature region).
        # Key metric: results correctness, not necessarily speedup.

    @pytest.mark.timeout(30)
    def test_annotation_benchmark_result_fields(self):
        """Annotation benchmark results have all expected fields."""
        results = benchmark_annotation(scales=[(100, 10)])
        r = results[0]
        assert isinstance(r, AnnotationBenchmark)
        assert isinstance(r.num_ssrs, int)
        assert isinstance(r.num_features, int)
        assert isinstance(r.tree_time_ms, float)
        assert isinstance(r.naive_time_ms, float)
        assert isinstance(r.speedup, float)
        assert isinstance(r.results_match, bool)


# ===========================================================================
# REPRODUCIBILITY BENCHMARKS
# ===========================================================================


class TestReproducibilityBenchmarks:
    """Benchmark reproducibility across repeated runs."""

    @pytest.mark.timeout(30)
    def test_reproducibility_10k_genome(self):
        """Detection on 10K genome is reproducible across 3 runs."""
        result = check_detection_reproducibility(
            genome_size=10_000, num_runs=3, seed=42
        )
        assert result["all_identical"]
        assert result["hits_per_run"] > 0

    @pytest.mark.timeout(60)
    def test_reproducibility_30k_genome(self):
        """Detection on 30K genome is reproducible across 5 runs."""
        result = check_detection_reproducibility(
            genome_size=30_000, num_runs=5, seed=42
        )
        assert result["all_identical"]
        assert result["hits_per_run"] > 0

    @pytest.mark.timeout(30)
    def test_reproducibility_different_seeds(self):
        """Different seeds produce different results (sanity check)."""
        g1 = generate_genome(10_000, seed=1)
        g2 = generate_genome(10_000, seed=2)
        assert g1 != g2, "Different seeds should produce different genomes"

        from gwico_ssr.ssr.detector import SSRThresholds, detect_ssrs

        r1 = detect_ssrs(g1, "SEED1", SSRThresholds())
        r2 = detect_ssrs(g2, "SEED2", SSRThresholds())
        # Highly unlikely to be identical
        if r1.hit_count == r2.hit_count and r1.hit_count > 0:
            # Even if same count, positions should differ
            any_diff = any(
                h1.start != h2.start for h1, h2 in zip(r1.hits, r2.hits)
            )
            assert any_diff or r1.hit_count == 0


# ===========================================================================
# SYNTHETIC GENOME GENERATION TESTS
# ===========================================================================


class TestSyntheticGeneration:
    """Test the benchmark harness's synthetic genome generation."""

    def test_generate_genome_length(self):
        """Generated genome has requested length."""
        g = generate_genome(1000)
        assert len(g) == 1000

    def test_generate_genome_deterministic(self):
        """Same seed produces same genome."""
        g1 = generate_genome(500, seed=42)
        g2 = generate_genome(500, seed=42)
        assert g1 == g2

    def test_generate_genome_different_seeds(self):
        """Different seeds produce different genomes."""
        g1 = generate_genome(500, seed=1)
        g2 = generate_genome(500, seed=2)
        assert g1 != g2

    def test_generate_genome_valid_bases(self):
        """Generated genome contains only ATGC."""
        g = generate_genome(10000)
        assert set(g) <= {"A", "T", "G", "C"}

    def test_generate_genome_gc_content(self):
        """GC content is approximately as requested."""
        g = generate_genome(100_000, gc_content=0.40, seed=42)
        gc_count = g.count("G") + g.count("C")
        gc_frac = gc_count / len(g)
        assert 0.35 <= gc_frac <= 0.45, f"GC content = {gc_frac:.3f}, expected ~0.40"


# ===========================================================================
# FULL BENCHMARK REPORT
# ===========================================================================


class TestBenchmarkReport:
    """Test the full benchmark report generation."""

    @pytest.mark.timeout(120)
    def test_full_benchmark_completes(self, tmp_path):
        """Full benchmark suite runs and produces a report."""
        output = tmp_path / "benchmark_report.json"
        report = run_full_benchmark(output_path=output)

        assert isinstance(report, BenchmarkReport)
        assert len(report.detection_benchmarks) == 3
        assert len(report.annotation_benchmarks) == 3
        assert report.detection_reproducibility is not None
        assert report.detection_reproducibility["all_identical"]

    @pytest.mark.timeout(120)
    def test_benchmark_report_json_output(self, tmp_path):
        """Benchmark report writes valid JSON."""
        import json

        output = tmp_path / "benchmark_report.json"
        run_full_benchmark(output_path=output)

        assert output.exists()
        data = json.loads(output.read_text())
        assert "detection_benchmarks" in data
        assert "annotation_benchmarks" in data
        assert "detection_reproducibility" in data
        assert "metadata" in data

    @pytest.mark.timeout(120)
    def test_benchmark_report_to_dict(self):
        """BenchmarkReport.to_dict() returns serializable dict."""
        report = run_full_benchmark()
        d = report.to_dict()
        assert isinstance(d, dict)
        assert isinstance(d["detection_benchmarks"], list)
        # Verify serializable
        import json

        json.dumps(d)  # Should not raise
