"""Pipeline orchestration, checkpointing, and recovery.

Defines the canonical pipeline stage ordering, provides checkpoint tracking
per run, failed-accession retry support, dry-run mode, and an orchestrator
that can execute or resume a multi-stage pipeline run.
"""

from __future__ import annotations

import enum
import logging
import time
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional, Sequence

from sqlalchemy import delete, func, select
from sqlalchemy.orm import Session

from gwico_ssr.models.schema import (
    Accession,
    AccessionMetrics,
    FailedAccession as FailedAccessionModel,
    FeatureRecord,
    SSRAnnotation as SSRAnnotationModel,
    SSRRecord as SSRRecordModel,
    StageCheckpoint as StageCheckpointModel,
    StatisticalResult,
)
from gwico_ssr.db.repository import (
    create_run,
    get_run,
    insert_ssr_annotations,
    insert_ssr_records,
    update_run_status,
    upsert_accession_metrics,
)

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Pipeline stage definitions
# ---------------------------------------------------------------------------


class PipelineStage(enum.Enum):
    """Ordered pipeline stages."""

    INGEST = "ingest"
    DOWNLOAD = "download"
    PARSE = "parse"
    DETECT = "detect"
    ANNOTATE = "annotate"
    METRICS = "metrics"
    ANALYZE = "analyze"
    VISUALIZE = "visualize"
    EXPORT = "export"


STAGE_ORDER: list[PipelineStage] = list(PipelineStage)

# Core stages that the orchestrator manages by default
CORE_STAGES: list[PipelineStage] = [
    PipelineStage.DETECT,
    PipelineStage.ANNOTATE,
    PipelineStage.METRICS,
    PipelineStage.ANALYZE,
]

# Stages that process accessions one at a time
PER_ACCESSION_STAGES = {
    PipelineStage.DETECT,
    PipelineStage.ANNOTATE,
    PipelineStage.METRICS,
}


def get_stage_index(stage: PipelineStage) -> int:
    """Return the ordinal index of a stage in the canonical ordering."""
    return STAGE_ORDER.index(stage)


# ---------------------------------------------------------------------------
# Result dataclasses
# ---------------------------------------------------------------------------


@dataclass
class StageResult:
    """Outcome of executing a single pipeline stage."""

    stage: PipelineStage
    status: str  # completed, failed, skipped, dry-run
    processed: int = 0
    skipped: int = 0
    failed: int = 0
    errors: list[str] = field(default_factory=list)
    duration_seconds: float = 0.0


@dataclass
class PipelineResult:
    """Outcome of a full or partial pipeline run."""

    run_id: int
    dataset_id: int
    stages: list[StageResult] = field(default_factory=list)
    status: str = "completed"  # completed, failed, dry-run


# ---------------------------------------------------------------------------
# Checkpoint helpers
# ---------------------------------------------------------------------------


def record_checkpoint(
    session: Session,
    run_id: int,
    stage: PipelineStage,
    status: str,
    processed: int = 0,
    skipped: int = 0,
    failed: int = 0,
) -> None:
    """Record or update a stage checkpoint for a run."""
    now = datetime.now(timezone.utc)
    existing = session.execute(
        select(StageCheckpointModel).where(
            StageCheckpointModel.run_id == run_id,
            StageCheckpointModel.stage == stage.value,
        )
    ).scalar_one_or_none()

    if existing:
        existing.status = status
        existing.finished_at = now if status in ("completed", "failed") else None
        existing.processed_count = processed
        existing.skipped_count = skipped
        existing.failed_count = failed
    else:
        cp = StageCheckpointModel(
            run_id=run_id,
            stage=stage.value,
            status=status,
            started_at=now,
            finished_at=now if status in ("completed", "failed") else None,
            processed_count=processed,
            skipped_count=skipped,
            failed_count=failed,
        )
        session.add(cp)
    session.flush()


def get_checkpoint(
    session: Session, run_id: int, stage: PipelineStage
) -> Optional[StageCheckpointModel]:
    """Get checkpoint for a specific stage of a run."""
    return session.execute(
        select(StageCheckpointModel).where(
            StageCheckpointModel.run_id == run_id,
            StageCheckpointModel.stage == stage.value,
        )
    ).scalar_one_or_none()


def get_completed_stages(session: Session, run_id: int) -> list[PipelineStage]:
    """Get all completed stages for a run, in canonical order."""
    rows = (
        session.execute(
            select(StageCheckpointModel.stage).where(
                StageCheckpointModel.run_id == run_id,
                StageCheckpointModel.status == "completed",
            )
        )
        .scalars()
        .all()
    )
    completed = set(rows)
    return [s for s in STAGE_ORDER if s.value in completed]


def get_last_completed_stage(
    session: Session, run_id: int
) -> Optional[PipelineStage]:
    """Get the last completed stage for a run."""
    completed = get_completed_stages(session, run_id)
    return completed[-1] if completed else None


# ---------------------------------------------------------------------------
# Failed accession helpers
# ---------------------------------------------------------------------------


def record_failed_accession(
    session: Session,
    run_id: int,
    accession: str,
    stage: PipelineStage,
    error_message: str,
) -> None:
    """Record or update a failed accession for later retry."""
    existing = session.execute(
        select(FailedAccessionModel).where(
            FailedAccessionModel.run_id == run_id,
            FailedAccessionModel.accession == accession,
            FailedAccessionModel.stage == stage.value,
        )
    ).scalar_one_or_none()

    if existing:
        existing.retry_count += 1
        existing.error_message = error_message
        existing.resolved = False
    else:
        fa = FailedAccessionModel(
            run_id=run_id,
            accession=accession,
            stage=stage.value,
            error_message=error_message,
            retry_count=0,
            resolved=False,
        )
        session.add(fa)
    session.flush()


def get_failed_accessions(
    session: Session,
    run_id: int,
    stage: Optional[PipelineStage] = None,
    unresolved_only: bool = True,
) -> list[FailedAccessionModel]:
    """Get failed accessions, optionally filtered by stage."""
    stmt = select(FailedAccessionModel).where(
        FailedAccessionModel.run_id == run_id
    )
    if stage:
        stmt = stmt.where(FailedAccessionModel.stage == stage.value)
    if unresolved_only:
        stmt = stmt.where(FailedAccessionModel.resolved == False)  # noqa: E712
    return list(session.execute(stmt).scalars().all())


def resolve_failed_accession(
    session: Session,
    run_id: int,
    accession: str,
    stage: PipelineStage,
) -> None:
    """Mark a failed accession as resolved after successful retry."""
    existing = session.execute(
        select(FailedAccessionModel).where(
            FailedAccessionModel.run_id == run_id,
            FailedAccessionModel.accession == accession,
            FailedAccessionModel.stage == stage.value,
        )
    ).scalar_one_or_none()
    if existing:
        existing.resolved = True
        session.flush()


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


def _get_dataset_accessions(session: Session, dataset_id: int) -> list[str]:
    """Get all accession IDs for a dataset."""
    return list(
        session.execute(
            select(Accession.accession).where(
                Accession.dataset_id == dataset_id
            )
        )
        .scalars()
        .all()
    )


# ---------------------------------------------------------------------------
# Stage executors
# ---------------------------------------------------------------------------


def execute_detect(
    session: Session,
    run_id: int,
    dataset_id: int,
    data_dir: str,
    thresholds: Optional[object] = None,
    force: bool = False,
    dry_run: bool = False,
) -> StageResult:
    """Execute SSR detection for all accessions in a dataset."""
    from gwico_ssr.parsers.fasta_parser import parse_fasta
    from gwico_ssr.ssr.detector import SSRThresholds, detect_ssrs

    if thresholds is None:
        thresholds = SSRThresholds()

    t0 = time.monotonic()
    accessions = _get_dataset_accessions(session, dataset_id)
    processed = 0
    skipped_count = 0
    failed_count = 0
    errors: list[str] = []

    for acc_id in accessions:
        if not force:
            existing = session.execute(
                select(func.count()).where(SSRRecordModel.accession == acc_id)
            ).scalar()
            if existing and existing > 0:
                skipped_count += 1
                continue

        if dry_run:
            processed += 1
            continue

        try:
            fasta_path = Path(data_dir) / "sequences" / "fasta" / f"{acc_id}.fasta"
            if not fasta_path.exists():
                skipped_count += 1
                continue

            fasta_result = parse_fasta(str(fasta_path))
            if not fasta_result.records:
                skipped_count += 1
                continue

            seq_rec = fasta_result.records[0]
            for r in fasta_result.records:
                if r.accession == acc_id:
                    seq_rec = r
                    break

            if force:
                session.execute(
                    delete(SSRRecordModel).where(
                        SSRRecordModel.accession == acc_id
                    )
                )
                session.flush()

            detection = detect_ssrs(seq_rec.sequence, acc_id, thresholds)
            if detection.hits:
                dicts = detection.to_dicts()
                for d in dicts:
                    d["run_id"] = run_id
                insert_ssr_records(session, dicts)

            resolve_failed_accession(session, run_id, acc_id, PipelineStage.DETECT)
            processed += 1

        except Exception as exc:
            failed_count += 1
            msg = f"{acc_id}: {exc}"
            errors.append(msg)
            record_failed_accession(
                session, run_id, acc_id, PipelineStage.DETECT, str(exc)
            )

    status = "dry-run" if dry_run else ("completed" if failed_count == 0 else "failed")
    return StageResult(
        stage=PipelineStage.DETECT,
        status=status,
        processed=processed,
        skipped=skipped_count,
        failed=failed_count,
        errors=errors,
        duration_seconds=time.monotonic() - t0,
    )


def execute_annotate(
    session: Session,
    run_id: int,
    dataset_id: int,
    force: bool = False,
    dry_run: bool = False,
) -> StageResult:
    """Execute annotation mapping for all accessions."""
    from gwico_ssr.annotation.mapper import annotate_accession

    t0 = time.monotonic()
    accessions = _get_dataset_accessions(session, dataset_id)
    processed = 0
    skipped_count = 0
    failed_count = 0
    errors: list[str] = []

    for acc_id in accessions:
        # Skip if no SSRs
        ssr_count = session.execute(
            select(func.count()).where(SSRRecordModel.accession == acc_id)
        ).scalar()
        if not ssr_count or ssr_count == 0:
            skipped_count += 1
            continue

        if not force:
            existing = session.execute(
                select(func.count()).where(
                    SSRAnnotationModel.accession == acc_id
                )
            ).scalar()
            if existing and existing > 0:
                skipped_count += 1
                continue

        if dry_run:
            processed += 1
            continue

        try:
            if force:
                session.execute(
                    delete(SSRAnnotationModel).where(
                        SSRAnnotationModel.accession == acc_id
                    )
                )
                session.flush()

            result = annotate_accession(session, acc_id)
            if result.annotations:
                insert_ssr_annotations(session, result.to_dicts())

            resolve_failed_accession(
                session, run_id, acc_id, PipelineStage.ANNOTATE
            )
            processed += 1

        except Exception as exc:
            failed_count += 1
            msg = f"{acc_id}: {exc}"
            errors.append(msg)
            record_failed_accession(
                session, run_id, acc_id, PipelineStage.ANNOTATE, str(exc)
            )

    status = "dry-run" if dry_run else ("completed" if failed_count == 0 else "failed")
    return StageResult(
        stage=PipelineStage.ANNOTATE,
        status=status,
        processed=processed,
        skipped=skipped_count,
        failed=failed_count,
        errors=errors,
        duration_seconds=time.monotonic() - t0,
    )


def execute_metrics(
    session: Session,
    run_id: int,
    dataset_id: int,
    force: bool = False,
    dry_run: bool = False,
) -> StageResult:
    """Compute per-accession metrics."""
    from gwico_ssr.db.repository import has_accession_metrics
    from gwico_ssr.metrics.calculator import compute_metrics_for_accession

    t0 = time.monotonic()
    # Only accessions that have SSRs
    acc_ids = list(
        session.execute(
            select(SSRRecordModel.accession.distinct())
        )
        .scalars()
        .all()
    )
    processed = 0
    skipped_count = 0
    failed_count = 0
    errors: list[str] = []

    for acc_id in acc_ids:
        if not force and has_accession_metrics(session, acc_id, run_id):
            skipped_count += 1
            continue

        if dry_run:
            processed += 1
            continue

        try:
            result = compute_metrics_for_accession(session, acc_id)
            upsert_accession_metrics(session, **result.to_dict(run_id))
            resolve_failed_accession(
                session, run_id, acc_id, PipelineStage.METRICS
            )
            processed += 1

        except Exception as exc:
            failed_count += 1
            msg = f"{acc_id}: {exc}"
            errors.append(msg)
            record_failed_accession(
                session, run_id, acc_id, PipelineStage.METRICS, str(exc)
            )

    status = "dry-run" if dry_run else ("completed" if failed_count == 0 else "failed")
    return StageResult(
        stage=PipelineStage.METRICS,
        status=status,
        processed=processed,
        skipped=skipped_count,
        failed=failed_count,
        errors=errors,
        duration_seconds=time.monotonic() - t0,
    )


def execute_analyze(
    session: Session,
    run_id: int,
    dataset_id: int,
    force: bool = False,
    dry_run: bool = False,
) -> StageResult:
    """Run statistical analyses on computed metrics."""
    from gwico_ssr.analysis.statistics import run_all_analyses

    t0 = time.monotonic()

    if dry_run:
        return StageResult(
            stage=PipelineStage.ANALYZE,
            status="dry-run",
            processed=1,
            duration_seconds=time.monotonic() - t0,
        )

    try:
        # Check for existing results
        existing = session.execute(
            select(func.count(StatisticalResult.result_id)).where(
                StatisticalResult.run_id == run_id
            )
        ).scalar() or 0

        if existing > 0 and not force:
            return StageResult(
                stage=PipelineStage.ANALYZE,
                status="skipped",
                skipped=1,
                duration_seconds=time.monotonic() - t0,
            )

        if force and existing > 0:
            session.execute(
                delete(StatisticalResult).where(
                    StatisticalResult.run_id == run_id
                )
            )
            session.flush()

        suite = run_all_analyses(
            session, run_id=run_id, dataset_id=dataset_id
        )
        suite.apply_fdr()
        suite.persist(session, run_id)

        return StageResult(
            stage=PipelineStage.ANALYZE,
            status="completed",
            processed=len(suite.results),
            duration_seconds=time.monotonic() - t0,
        )
    except Exception as exc:
        return StageResult(
            stage=PipelineStage.ANALYZE,
            status="failed",
            failed=1,
            errors=[str(exc)],
            duration_seconds=time.monotonic() - t0,
        )


# ---------------------------------------------------------------------------
# Pipeline orchestrator
# ---------------------------------------------------------------------------

# Map stages to executor functions
_STAGE_EXECUTORS = {
    PipelineStage.DETECT: execute_detect,
    PipelineStage.ANNOTATE: execute_annotate,
    PipelineStage.METRICS: execute_metrics,
    PipelineStage.ANALYZE: execute_analyze,
}


class PipelineOrchestrator:
    """Orchestrates multi-stage pipeline execution with checkpointing.

    Usage::

        orch = PipelineOrchestrator(session, dataset_id=1, data_dir="outputs")
        result = orch.run()           # full run
        result = orch.resume(run_id)  # resume interrupted run
        result = orch.run(dry_run=True)  # plan without executing
    """

    def __init__(
        self,
        session: Session,
        dataset_id: int,
        data_dir: str = "outputs",
        thresholds: Optional[object] = None,
    ) -> None:
        self.session = session
        self.dataset_id = dataset_id
        self.data_dir = data_dir
        self.thresholds = thresholds

    def run(
        self,
        stages: Optional[Sequence[PipelineStage]] = None,
        force: bool = False,
        dry_run: bool = False,
        retry_failed: bool = False,
    ) -> PipelineResult:
        """Execute a new pipeline run through the specified stages.

        Parameters
        ----------
        stages
            Which stages to execute, in order. Defaults to CORE_STAGES.
        force
            Re-process even if results already exist.
        dry_run
            Plan the run without making changes.
        retry_failed
            If True, retry previously failed accessions.
        """
        target_stages = list(stages or CORE_STAGES)

        run = create_run(
            self.session, dataset_id=self.dataset_id, stage="pipeline"
        )
        update_run_status(
            self.session, run.run_id, "running", stage="pipeline"
        )

        result = PipelineResult(
            run_id=run.run_id,
            dataset_id=self.dataset_id,
        )

        if dry_run:
            result.status = "dry-run"

        for stage in target_stages:
            stage_result = self._execute_stage(
                run.run_id, stage, force=force, dry_run=dry_run
            )
            result.stages.append(stage_result)

            if stage_result.status == "failed" and not dry_run:
                result.status = "failed"
                update_run_status(
                    self.session,
                    run.run_id,
                    "failed",
                    stage=stage.value,
                    notes="; ".join(stage_result.errors[:5]),
                )
                break

        if result.status not in ("failed", "dry-run"):
            result.status = "completed"
            update_run_status(
                self.session, run.run_id, "completed", stage="pipeline"
            )

        return result

    def resume(
        self,
        run_id: int,
        force: bool = False,
        dry_run: bool = False,
        retry_failed: bool = False,
    ) -> PipelineResult:
        """Resume an interrupted or failed pipeline run.

        Skips stages that already completed unless *force* is True.
        """
        run = get_run(self.session, run_id)
        if run is None:
            raise ValueError(f"Run {run_id} not found")

        completed = set(get_completed_stages(self.session, run_id))

        update_run_status(
            self.session, run_id, "running", stage="pipeline"
        )

        result = PipelineResult(
            run_id=run_id,
            dataset_id=run.dataset_id,
        )

        if dry_run:
            result.status = "dry-run"

        for stage in CORE_STAGES:
            if stage in completed and not force:
                result.stages.append(
                    StageResult(stage=stage, status="skipped", skipped=1)
                )
                continue

            stage_result = self._execute_stage(
                run_id,
                stage,
                force=force,
                dry_run=dry_run,
            )
            result.stages.append(stage_result)

            if stage_result.status == "failed" and not dry_run:
                result.status = "failed"
                update_run_status(
                    self.session,
                    run_id,
                    "failed",
                    stage=stage.value,
                    notes="; ".join(stage_result.errors[:5]),
                )
                break

        if result.status not in ("failed", "dry-run"):
            result.status = "completed"
            update_run_status(
                self.session, run_id, "completed", stage="pipeline"
            )

        return result

    def get_status(self, run_id: int) -> dict:
        """Get the current status of a pipeline run."""
        run = get_run(self.session, run_id)
        if run is None:
            return {"error": f"Run {run_id} not found"}

        completed = get_completed_stages(self.session, run_id)
        failed = get_failed_accessions(self.session, run_id)

        return {
            "run_id": run_id,
            "status": run.status,
            "stage": run.stage,
            "completed_stages": [s.value for s in completed],
            "failed_accessions": len(failed),
            "started_at": str(run.started_at) if run.started_at else None,
            "finished_at": str(run.finished_at) if run.finished_at else None,
        }

    def _execute_stage(
        self,
        run_id: int,
        stage: PipelineStage,
        force: bool = False,
        dry_run: bool = False,
    ) -> StageResult:
        """Execute a single pipeline stage with checkpoint tracking."""
        executor = _STAGE_EXECUTORS.get(stage)
        if executor is None:
            return StageResult(
                stage=stage,
                status="skipped",
                skipped=1,
            )

        logger.info("Starting stage %s for run %d", stage.value, run_id)
        record_checkpoint(self.session, run_id, stage, "running")
        update_run_status(self.session, run_id, "running", stage=stage.value)

        # Build kwargs for the executor
        kwargs: dict = {
            "session": self.session,
            "run_id": run_id,
            "dataset_id": self.dataset_id,
            "force": force,
            "dry_run": dry_run,
        }

        # Stage-specific kwargs
        if stage == PipelineStage.DETECT:
            kwargs["data_dir"] = self.data_dir
            kwargs["thresholds"] = self.thresholds

        stage_result = executor(**kwargs)

        record_checkpoint(
            self.session,
            run_id,
            stage,
            stage_result.status,
            processed=stage_result.processed,
            skipped=stage_result.skipped,
            failed=stage_result.failed,
        )

        logger.info(
            "Stage %s: status=%s processed=%d skipped=%d failed=%d (%.1fs)",
            stage.value,
            stage_result.status,
            stage_result.processed,
            stage_result.skipped,
            stage_result.failed,
            stage_result.duration_seconds,
        )

        return stage_result
