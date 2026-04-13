"""Benchmark harness for GWICO-SSR pipeline stages.

Provides utilities for generating synthetic genomes at various scales,
running detection + annotation benchmarks, and producing structured
benchmark reports.
"""

from __future__ import annotations

import json
import random
import time
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Optional

from gwico_ssr.annotation.mapper import (
    FeatureInterval,
    annotate_ssrs,
    benchmark_compare,
    build_feature_tree,
)
from gwico_ssr.ssr.detector import DetectionResult, SSRThresholds, detect_ssrs


# ---------------------------------------------------------------------------
# Synthetic genome generation
# ---------------------------------------------------------------------------


def generate_genome(length: int, gc_content: float = 0.38, seed: int = 42) -> str:
    """Generate a synthetic genome of specified length.

    Parameters
    ----------
    length
        Total bases in the generated sequence.
    gc_content
        Fraction of G+C bases (0.0–1.0).
    seed
        Random seed for reproducibility.
    """
    rng = random.Random(seed)
    at_weight = (1.0 - gc_content) / 2.0
    gc_weight = gc_content / 2.0
    bases = rng.choices(
        "ATGC", weights=[at_weight, at_weight, gc_weight, gc_weight], k=length
    )
    return "".join(bases)


def inject_ssrs(
    sequence: str,
    ssrs: list[dict],
    seed: int = 42,
) -> str:
    """Inject known SSR patterns into a sequence at specified positions.

    Each SSR dict should have: motif, repeat_units, position.
    """
    seq_list = list(sequence)
    for ssr in ssrs:
        motif = ssr["motif"]
        repeat_units = ssr["repeat_units"]
        pos = ssr["position"]
        repeat_seq = motif * repeat_units
        for i, base in enumerate(repeat_seq):
            if pos + i < len(seq_list):
                seq_list[pos + i] = base
    return "".join(seq_list)


def generate_features(
    genome_length: int, num_features: int, seed: int = 42
) -> list[FeatureInterval]:
    """Generate synthetic gene features distributed across a genome."""
    rng = random.Random(seed)
    features = []
    gene_names = [
        "ORF1ab", "S", "ORF3a", "E", "M", "ORF6", "ORF7a", "ORF8", "N", "ORF10",
    ]
    spacing = genome_length // (num_features + 1)
    for i in range(num_features):
        start = spacing * (i + 1) - spacing // 4
        end = start + rng.randint(200, spacing // 2)
        if end >= genome_length:
            end = genome_length - 1
        features.append(
            FeatureInterval(
                feature_id=i + 1,
                accession="BENCH_ACC",
                feature_type="CDS",
                start=max(0, start),
                end=end,
                strand="+",
                gene_name=gene_names[i % len(gene_names)],
            )
        )
    return features


# ---------------------------------------------------------------------------
# Benchmark result dataclasses
# ---------------------------------------------------------------------------


@dataclass
class DetectionBenchmark:
    """Result of a single SSR detection benchmark run."""

    genome_size: int
    ssr_count: int
    detection_time_ms: float
    ssrs_per_second: float


@dataclass
class AnnotationBenchmark:
    """Result of a single annotation benchmark run."""

    num_ssrs: int
    num_features: int
    tree_time_ms: float
    naive_time_ms: float
    speedup: float
    results_match: bool


@dataclass
class BenchmarkReport:
    """Complete benchmark report for a suite of tests."""

    detection_benchmarks: list[DetectionBenchmark] = field(default_factory=list)
    annotation_benchmarks: list[AnnotationBenchmark] = field(default_factory=list)
    detection_reproducibility: Optional[dict] = None
    metadata: dict = field(default_factory=dict)

    def to_dict(self) -> dict:
        return {
            "detection_benchmarks": [asdict(b) for b in self.detection_benchmarks],
            "annotation_benchmarks": [asdict(b) for b in self.annotation_benchmarks],
            "detection_reproducibility": self.detection_reproducibility,
            "metadata": self.metadata,
        }

    def to_json(self, path: Path) -> None:
        path.write_text(json.dumps(self.to_dict(), indent=2))


# ---------------------------------------------------------------------------
# Benchmark runners
# ---------------------------------------------------------------------------


def benchmark_detection(
    genome_sizes: list[int],
    thresholds: Optional[SSRThresholds] = None,
    gc_content: float = 0.38,
    seed: int = 42,
) -> list[DetectionBenchmark]:
    """Benchmark SSR detection at various genome sizes.

    Parameters
    ----------
    genome_sizes
        List of genome sizes to benchmark (in base pairs).
    thresholds
        Detection thresholds. Uses defaults if None.
    gc_content
        GC content for synthetic genomes.
    seed
        Random seed for reproducibility.

    Returns
    -------
    List of DetectionBenchmark results, one per genome size.
    """
    if thresholds is None:
        thresholds = SSRThresholds()

    results = []
    for size in genome_sizes:
        genome = generate_genome(size, gc_content=gc_content, seed=seed)

        t0 = time.perf_counter()
        detection = detect_ssrs(genome, f"BENCH_{size}", thresholds)
        elapsed = time.perf_counter() - t0

        elapsed_ms = elapsed * 1000
        ssrs_sec = detection.hit_count / elapsed if elapsed > 0 else 0

        results.append(
            DetectionBenchmark(
                genome_size=size,
                ssr_count=detection.hit_count,
                detection_time_ms=round(elapsed_ms, 2),
                ssrs_per_second=round(ssrs_sec, 1),
            )
        )

    return results


def benchmark_annotation(
    scales: list[tuple[int, int]],
    seed: int = 42,
) -> list[AnnotationBenchmark]:
    """Benchmark annotation mapping at various SSR × feature scales.

    Parameters
    ----------
    scales
        List of (num_ssrs, num_features) tuples to benchmark.
    seed
        Random seed for reproducibility.

    Returns
    -------
    List of AnnotationBenchmark results.
    """
    rng = random.Random(seed)
    results = []

    for num_ssrs, num_features in scales:
        genome_length = max(num_features * 200, 50000)

        features = [
            FeatureInterval(
                feature_id=i,
                accession="BENCH",
                feature_type="CDS",
                start=i * (genome_length // num_features),
                end=i * (genome_length // num_features) + 80,
                strand="+",
                gene_name=f"gene_{i}",
            )
            for i in range(num_features)
        ]

        class _FakeSSR:
            def __init__(self, ssr_id, start, end):
                self.ssr_id = ssr_id
                self.start = start
                self.end = end

        ssrs = [
            _FakeSSR(
                ssr_id=i,
                start=rng.randint(0, genome_length - 20),
                end=rng.randint(0, genome_length - 20) + 10,
            )
            for i in range(num_ssrs)
        ]

        comparison = benchmark_compare(ssrs, features, "BENCH")

        results.append(
            AnnotationBenchmark(
                num_ssrs=num_ssrs,
                num_features=num_features,
                tree_time_ms=round(comparison["tree_time_ms"], 2),
                naive_time_ms=round(comparison["naive_time_ms"], 2),
                speedup=round(comparison["speedup"], 2),
                results_match=comparison["results_match"],
            )
        )

    return results


def check_detection_reproducibility(
    genome_size: int = 30000,
    num_runs: int = 3,
    seed: int = 42,
) -> dict:
    """Verify that repeated detection runs produce identical results.

    Returns a dict with run details and whether all runs matched.
    """
    genome = generate_genome(genome_size, seed=seed)
    thresholds = SSRThresholds()

    results: list[DetectionResult] = []
    for _ in range(num_runs):
        result = detect_ssrs(genome, "REPRO_TEST", thresholds)
        results.append(result)

    # Compare all runs to the first
    first = results[0]
    all_match = True
    for i, r in enumerate(results[1:], 1):
        if r.hit_count != first.hit_count:
            all_match = False
            break
        for h1, h2 in zip(first.hits, r.hits):
            if (
                h1.start != h2.start
                or h1.end != h2.end
                or h1.motif_canonical != h2.motif_canonical
                or h1.repeat_units != h2.repeat_units
            ):
                all_match = False
                break

    return {
        "genome_size": genome_size,
        "num_runs": num_runs,
        "hits_per_run": first.hit_count,
        "all_identical": all_match,
    }


def run_full_benchmark(
    output_path: Optional[Path] = None,
) -> BenchmarkReport:
    """Run the complete benchmark suite and return a report.

    Benchmarks SSR detection at 1K, 10K, 30K, 100K genome sizes,
    annotation mapping at multiple SSR × feature scales,
    and reproducibility verification.
    """
    import platform

    report = BenchmarkReport(
        metadata={
            "platform": platform.platform(),
            "python": platform.python_version(),
            "processor": platform.processor(),
        }
    )

    # Detection benchmarks
    report.detection_benchmarks = benchmark_detection(
        genome_sizes=[1_000, 10_000, 30_000]
    )

    # Annotation benchmarks
    report.annotation_benchmarks = benchmark_annotation(
        scales=[
            (100, 10),
            (500, 50),
            (1_000, 100),
        ]
    )

    # Reproducibility
    report.detection_reproducibility = check_detection_reproducibility()

    if output_path:
        report.to_json(output_path)

    return report
