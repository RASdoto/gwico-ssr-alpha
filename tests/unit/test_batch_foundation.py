"""Tests for Chunk 1 batch-foundation orchestration helpers."""

from __future__ import annotations

from gwico_ssr.config import load_settings
from gwico_ssr.orchestration.batch_foundation import (
    build_batch_foundation_plan,
    make_artifact_base_dir,
)


def test_batch_foundation_plan_defaults():
    settings = load_settings("/nonexistent/path.toml")
    plan = build_batch_foundation_plan(settings)

    assert plan.request_batch_size == 500
    assert plan.artifact_batch_size == 10000
    assert plan.parse_chunk_size == 1000
    assert plan.retain_raw_artifacts is True
    assert plan.manifest_policy == "required"
    assert plan.duplicate_handling == "error"


def test_batch_foundation_plan_from_toml(sample_toml):
    settings = load_settings(str(sample_toml))
    plan = build_batch_foundation_plan(settings)

    assert plan.request_batch_size == 700
    assert plan.artifact_batch_size == 12000
    assert plan.parse_chunk_size == 1500


def test_artifact_base_dir_is_deterministic():
    base = make_artifact_base_dir("outputs")
    assert str(base).replace("\\", "/") == "outputs/batches/raw"
