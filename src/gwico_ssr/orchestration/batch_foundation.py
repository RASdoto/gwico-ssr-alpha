"""Batch-first foundation helpers for planning artifact-oriented workflows."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from gwico_ssr.config import Settings


@dataclass(frozen=True)
class BatchFoundationPlan:
    request_batch_size: int
    artifact_batch_size: int
    parse_chunk_size: int
    retain_raw_artifacts: bool
    manifest_policy: str
    duplicate_handling: str


def build_batch_foundation_plan(settings: Settings) -> BatchFoundationPlan:
    """Build a resolved batch foundation plan from loaded settings."""
    request_batch_size = settings.ncbi.request_batch_size or settings.ncbi.batch_size
    return BatchFoundationPlan(
        request_batch_size=request_batch_size,
        artifact_batch_size=settings.batch.artifact_batch_size,
        parse_chunk_size=settings.batch.parse_chunk_size,
        retain_raw_artifacts=settings.batch.retain_raw_artifacts,
        manifest_policy=settings.batch.manifest_policy,
        duplicate_handling=settings.batch.duplicate_handling,
    )


def make_artifact_base_dir(output_dir: str | Path) -> Path:
    """Return deterministic base directory for raw batch artifacts."""
    return Path(output_dir) / "batches" / "raw"
