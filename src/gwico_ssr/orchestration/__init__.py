"""Orchestration, checkpointing, and recovery for GWICO-SSR pipeline."""

from gwico_ssr.orchestration.batch_foundation import (
    BatchFoundationPlan,
    build_batch_foundation_plan,
    make_artifact_base_dir,
)

from gwico_ssr.orchestration.pipeline import (
    CORE_STAGES,
    STAGE_ORDER,
    PipelineOrchestrator,
    PipelineResult,
    PipelineStage,
    StageResult,
    execute_analyze,
    execute_annotate,
    execute_detect,
    execute_metrics,
    get_checkpoint,
    get_completed_stages,
    get_failed_accessions,
    get_last_completed_stage,
    record_checkpoint,
    record_failed_accession,
    resolve_failed_accession,
)

__all__ = [
    "BatchFoundationPlan",
    "CORE_STAGES",
    "STAGE_ORDER",
    "build_batch_foundation_plan",
    "make_artifact_base_dir",
    "PipelineOrchestrator",
    "PipelineResult",
    "PipelineStage",
    "StageResult",
    "execute_analyze",
    "execute_annotate",
    "execute_detect",
    "execute_metrics",
    "get_checkpoint",
    "get_completed_stages",
    "get_failed_accessions",
    "get_last_completed_stage",
    "record_checkpoint",
    "record_failed_accession",
    "resolve_failed_accession",
]
