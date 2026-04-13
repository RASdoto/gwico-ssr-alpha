"""Tests for GWICO-SSR orchestration, checkpointing, and recovery."""

from __future__ import annotations

import json
from pathlib import Path

import pytest
from sqlalchemy import create_engine, select, func
from sqlalchemy.orm import sessionmaker

from gwico_ssr.db import create_tables
from gwico_ssr.db.repository import (
    create_run,
    get_or_create_dataset,
    insert_features,
    insert_ssr_records,
    update_run_status,
    upsert_accession,
    upsert_accession_metrics,
)
from gwico_ssr.models.schema import (
    AccessionMetrics,
    FailedAccession as FailedAccessionModel,
    SSRAnnotation as SSRAnnotationModel,
    SSRRecord as SSRRecordModel,
    StageCheckpoint as StageCheckpointModel,
    StatisticalResult,
)
from gwico_ssr.orchestration.pipeline import (
    CORE_STAGES,
    STAGE_ORDER,
    PipelineOrchestrator,
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


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def db_session():
    engine = create_engine("sqlite:///:memory:")
    create_tables(engine)
    factory = sessionmaker(bind=engine)
    session = factory()
    yield session
    session.close()


@pytest.fixture
def dataset_with_accessions(db_session):
    """Dataset with 3 accessions, features, and FASTA files (no SSRs yet)."""
    session = db_session
    ds = get_or_create_dataset(session, "orch_test", "csv")

    for i in range(3):
        acc_id = f"ACC_{i:03d}"
        upsert_accession(
            session,
            accession=acc_id,
            dataset_id=ds.dataset_id,
            country=["USA", "China", "USA"][i],
            genome_length=29000 + i * 100,
            gc_content=0.38,
        )
        # Add CDS feature for annotation
        insert_features(
            session,
            [
                {
                    "accession": acc_id,
                    "feature_type": "CDS",
                    "start": 0,
                    "end": 15000,
                    "strand": "+",
                    "gene_name": "ORF1ab",
                    "product": "polyprotein",
                    "locus_tag": None,
                    "annotation_source": "genbank",
                },
                {
                    "accession": acc_id,
                    "feature_type": "CDS",
                    "start": 20000,
                    "end": 25000,
                    "strand": "+",
                    "gene_name": "S",
                    "product": "spike protein",
                    "locus_tag": None,
                    "annotation_source": "genbank",
                },
            ],
        )
    session.flush()
    return session, ds


@pytest.fixture
def dataset_with_ssrs(dataset_with_accessions):
    """Dataset with accessions that already have SSRs and features."""
    session, ds = dataset_with_accessions
    run = create_run(session, dataset_id=ds.dataset_id, stage="pipeline")

    ssrs = []
    for i in range(3):
        acc_id = f"ACC_{i:03d}"
        # SSR inside ORF1ab
        ssrs.append(
            {
                "accession": acc_id,
                "start": 100,
                "end": 109,
                "motif_raw": "AAG",
                "motif_canonical": "AAG",
                "motif_size": 3,
                "repeat_units": 3,
                "repeat_length_bp": 9,
                "strand": "+",
                "actual_repeat": "AAGAAGAAG",
                "detector_version": "gwico-ssr-1.0",
                "run_id": run.run_id,
            }
        )
        # SSR in intergenic region
        ssrs.append(
            {
                "accession": acc_id,
                "start": 16000,
                "end": 16010,
                "motif_raw": "AT",
                "motif_canonical": "AT",
                "motif_size": 2,
                "repeat_units": 5,
                "repeat_length_bp": 10,
                "strand": "+",
                "actual_repeat": "ATATATATAT",
                "detector_version": "gwico-ssr-1.0",
                "run_id": run.run_id,
            }
        )
    insert_ssr_records(session, ssrs)
    session.flush()
    return session, ds, run


@pytest.fixture
def fasta_dir(tmp_path):
    """Create minimal FASTA files for detection tests."""
    fasta_base = tmp_path / "sequences" / "fasta"
    fasta_base.mkdir(parents=True)

    # A short sequence with a known SSR
    seq = "ATCG" * 50 + "AAGAAGAAGAAGAAG" + "GCTA" * 50  # AAG x 5

    for i in range(3):
        acc_id = f"ACC_{i:03d}"
        fasta_file = fasta_base / f"{acc_id}.fasta"
        fasta_file.write_text(f">{acc_id} test sequence\n{seq}\n")

    return str(tmp_path)


# ---------------------------------------------------------------------------
# Test: Pipeline Stage Ordering
# ---------------------------------------------------------------------------


class TestPipelineStage:
    def test_stage_order(self):
        assert STAGE_ORDER[0] == PipelineStage.INGEST
        assert STAGE_ORDER[-1] == PipelineStage.EXPORT
        assert len(STAGE_ORDER) == 9

    def test_core_stages(self):
        assert CORE_STAGES == [
            PipelineStage.DETECT,
            PipelineStage.ANNOTATE,
            PipelineStage.METRICS,
            PipelineStage.ANALYZE,
        ]

    def test_stage_values(self):
        assert PipelineStage.DETECT.value == "detect"
        assert PipelineStage.ANNOTATE.value == "annotate"
        assert PipelineStage.METRICS.value == "metrics"
        assert PipelineStage.ANALYZE.value == "analyze"


# ---------------------------------------------------------------------------
# Test: Checkpoint Tracking
# ---------------------------------------------------------------------------


class TestCheckpoint:
    def test_record_and_get_checkpoint(self, db_session):
        ds = get_or_create_dataset(db_session, "cp_test", "csv")
        run = create_run(db_session, dataset_id=ds.dataset_id, stage="test")

        record_checkpoint(
            db_session,
            run.run_id,
            PipelineStage.DETECT,
            "completed",
            processed=10,
            skipped=2,
            failed=1,
        )

        cp = get_checkpoint(db_session, run.run_id, PipelineStage.DETECT)
        assert cp is not None
        assert cp.stage == "detect"
        assert cp.status == "completed"
        assert cp.processed_count == 10
        assert cp.skipped_count == 2
        assert cp.failed_count == 1
        assert cp.finished_at is not None

    def test_checkpoint_update(self, db_session):
        ds = get_or_create_dataset(db_session, "cp_update", "csv")
        run = create_run(db_session, dataset_id=ds.dataset_id)

        record_checkpoint(
            db_session, run.run_id, PipelineStage.DETECT, "running"
        )
        cp = get_checkpoint(db_session, run.run_id, PipelineStage.DETECT)
        assert cp.status == "running"
        assert cp.finished_at is None

        record_checkpoint(
            db_session,
            run.run_id,
            PipelineStage.DETECT,
            "completed",
            processed=5,
        )
        cp = get_checkpoint(db_session, run.run_id, PipelineStage.DETECT)
        assert cp.status == "completed"
        assert cp.processed_count == 5
        assert cp.finished_at is not None

    def test_completed_stages(self, db_session):
        ds = get_or_create_dataset(db_session, "cp_completed", "csv")
        run = create_run(db_session, dataset_id=ds.dataset_id)

        record_checkpoint(
            db_session, run.run_id, PipelineStage.DETECT, "completed"
        )
        record_checkpoint(
            db_session, run.run_id, PipelineStage.ANNOTATE, "completed"
        )
        record_checkpoint(
            db_session, run.run_id, PipelineStage.METRICS, "running"
        )

        completed = get_completed_stages(db_session, run.run_id)
        assert len(completed) == 2
        assert PipelineStage.DETECT in completed
        assert PipelineStage.ANNOTATE in completed
        assert PipelineStage.METRICS not in completed

    def test_last_completed_stage(self, db_session):
        ds = get_or_create_dataset(db_session, "cp_last", "csv")
        run = create_run(db_session, dataset_id=ds.dataset_id)

        assert get_last_completed_stage(db_session, run.run_id) is None

        record_checkpoint(
            db_session, run.run_id, PipelineStage.DETECT, "completed"
        )
        assert (
            get_last_completed_stage(db_session, run.run_id)
            == PipelineStage.DETECT
        )

        record_checkpoint(
            db_session, run.run_id, PipelineStage.METRICS, "completed"
        )
        assert (
            get_last_completed_stage(db_session, run.run_id)
            == PipelineStage.METRICS
        )

    def test_completed_stages_order(self, db_session):
        """Completed stages are returned in canonical order."""
        ds = get_or_create_dataset(db_session, "cp_order", "csv")
        run = create_run(db_session, dataset_id=ds.dataset_id)

        # Record in reverse order
        record_checkpoint(
            db_session, run.run_id, PipelineStage.ANALYZE, "completed"
        )
        record_checkpoint(
            db_session, run.run_id, PipelineStage.DETECT, "completed"
        )

        completed = get_completed_stages(db_session, run.run_id)
        assert completed == [PipelineStage.DETECT, PipelineStage.ANALYZE]


# ---------------------------------------------------------------------------
# Test: Failed Accession Tracking
# ---------------------------------------------------------------------------


class TestFailedAccession:
    def test_record_and_query(self, db_session):
        ds = get_or_create_dataset(db_session, "fail_test", "csv")
        run = create_run(db_session, dataset_id=ds.dataset_id)

        record_failed_accession(
            db_session,
            run.run_id,
            "ACC_001",
            PipelineStage.DETECT,
            "FASTA not found",
        )

        failed = get_failed_accessions(db_session, run.run_id)
        assert len(failed) == 1
        assert failed[0].accession == "ACC_001"
        assert failed[0].stage == "detect"
        assert failed[0].error_message == "FASTA not found"
        assert failed[0].retry_count == 0
        assert failed[0].resolved is False

    def test_retry_increments_count(self, db_session):
        ds = get_or_create_dataset(db_session, "retry_test", "csv")
        run = create_run(db_session, dataset_id=ds.dataset_id)

        record_failed_accession(
            db_session,
            run.run_id,
            "ACC_001",
            PipelineStage.DETECT,
            "error 1",
        )
        record_failed_accession(
            db_session,
            run.run_id,
            "ACC_001",
            PipelineStage.DETECT,
            "error 2",
        )

        failed = get_failed_accessions(db_session, run.run_id)
        assert len(failed) == 1
        assert failed[0].retry_count == 1
        assert failed[0].error_message == "error 2"

    def test_resolve(self, db_session):
        ds = get_or_create_dataset(db_session, "resolve_test", "csv")
        run = create_run(db_session, dataset_id=ds.dataset_id)

        record_failed_accession(
            db_session,
            run.run_id,
            "ACC_001",
            PipelineStage.DETECT,
            "error",
        )
        resolve_failed_accession(
            db_session, run.run_id, "ACC_001", PipelineStage.DETECT
        )

        # Should not appear in unresolved queries
        failed = get_failed_accessions(
            db_session, run.run_id, unresolved_only=True
        )
        assert len(failed) == 0

        # Should appear when including resolved
        all_failed = get_failed_accessions(
            db_session, run.run_id, unresolved_only=False
        )
        assert len(all_failed) == 1
        assert all_failed[0].resolved is True

    def test_filter_by_stage(self, db_session):
        ds = get_or_create_dataset(db_session, "stage_filter", "csv")
        run = create_run(db_session, dataset_id=ds.dataset_id)

        record_failed_accession(
            db_session,
            run.run_id,
            "ACC_001",
            PipelineStage.DETECT,
            "err1",
        )
        record_failed_accession(
            db_session,
            run.run_id,
            "ACC_002",
            PipelineStage.ANNOTATE,
            "err2",
        )

        detect_fails = get_failed_accessions(
            db_session, run.run_id, stage=PipelineStage.DETECT
        )
        assert len(detect_fails) == 1
        assert detect_fails[0].accession == "ACC_001"

    def test_multiple_accessions(self, db_session):
        ds = get_or_create_dataset(db_session, "multi_fail", "csv")
        run = create_run(db_session, dataset_id=ds.dataset_id)

        for i in range(5):
            record_failed_accession(
                db_session,
                run.run_id,
                f"ACC_{i:03d}",
                PipelineStage.DETECT,
                f"error {i}",
            )

        failed = get_failed_accessions(db_session, run.run_id)
        assert len(failed) == 5


# ---------------------------------------------------------------------------
# Test: Stage Executors
# ---------------------------------------------------------------------------


class TestExecuteAnnotate:
    def test_annotate_creates_annotations(self, dataset_with_ssrs):
        session, ds, run = dataset_with_ssrs

        result = execute_annotate(
            session, run.run_id, ds.dataset_id, force=False, dry_run=False
        )

        assert result.status == "completed"
        assert result.processed == 3
        assert result.failed == 0

        # Verify annotations were created
        count = session.execute(
            select(func.count()).select_from(SSRAnnotationModel)
        ).scalar()
        assert count > 0

    def test_annotate_skips_existing(self, dataset_with_ssrs):
        session, ds, run = dataset_with_ssrs

        # First run
        execute_annotate(session, run.run_id, ds.dataset_id)

        # Second run — should skip
        result = execute_annotate(
            session, run.run_id, ds.dataset_id, force=False
        )
        assert result.skipped == 3
        assert result.processed == 0

    def test_annotate_force_rerun(self, dataset_with_ssrs):
        session, ds, run = dataset_with_ssrs

        execute_annotate(session, run.run_id, ds.dataset_id)

        result = execute_annotate(
            session, run.run_id, ds.dataset_id, force=True
        )
        assert result.processed == 3

    def test_annotate_dry_run(self, dataset_with_ssrs):
        session, ds, run = dataset_with_ssrs

        result = execute_annotate(
            session, run.run_id, ds.dataset_id, dry_run=True
        )
        assert result.status == "dry-run"
        assert result.processed == 3

        # Verify no annotations were actually created
        count = session.execute(
            select(func.count()).select_from(SSRAnnotationModel)
        ).scalar()
        assert count == 0


class TestExecuteMetrics:
    def test_metrics_computes(self, dataset_with_ssrs):
        session, ds, run = dataset_with_ssrs

        result = execute_metrics(
            session, run.run_id, ds.dataset_id, force=False
        )

        assert result.status == "completed"
        assert result.processed == 3
        assert result.failed == 0

        # Verify metrics were created
        count = session.execute(
            select(func.count()).select_from(AccessionMetrics)
        ).scalar()
        assert count == 3

    def test_metrics_skips_existing(self, dataset_with_ssrs):
        session, ds, run = dataset_with_ssrs

        execute_metrics(session, run.run_id, ds.dataset_id)
        result = execute_metrics(
            session, run.run_id, ds.dataset_id, force=False
        )
        assert result.skipped == 3
        assert result.processed == 0

    def test_metrics_dry_run(self, dataset_with_ssrs):
        session, ds, run = dataset_with_ssrs

        result = execute_metrics(
            session, run.run_id, ds.dataset_id, dry_run=True
        )
        assert result.status == "dry-run"

        count = session.execute(
            select(func.count()).select_from(AccessionMetrics)
        ).scalar()
        assert count == 0


class TestExecuteDetect:
    def test_detect_with_fasta(self, dataset_with_accessions, fasta_dir):
        session, ds = dataset_with_accessions
        run = create_run(session, dataset_id=ds.dataset_id, stage="pipeline")

        result = execute_detect(
            session, run.run_id, ds.dataset_id, data_dir=fasta_dir
        )

        assert result.status == "completed"
        assert result.processed > 0

        ssr_count = session.execute(
            select(func.count()).select_from(SSRRecordModel)
        ).scalar()
        assert ssr_count > 0

    def test_detect_skips_existing(self, dataset_with_ssrs, fasta_dir):
        session, ds, run = dataset_with_ssrs

        result = execute_detect(
            session, run.run_id, ds.dataset_id, data_dir=fasta_dir, force=False
        )
        assert result.skipped == 3
        assert result.processed == 0

    def test_detect_dry_run(self, dataset_with_accessions, fasta_dir):
        session, ds = dataset_with_accessions
        run = create_run(session, dataset_id=ds.dataset_id, stage="pipeline")

        result = execute_detect(
            session, run.run_id, ds.dataset_id, data_dir=fasta_dir, dry_run=True
        )
        assert result.status == "dry-run"

        ssr_count = session.execute(
            select(func.count()).select_from(SSRRecordModel)
        ).scalar()
        assert ssr_count == 0

    def test_detect_missing_fasta_skips(self, dataset_with_accessions):
        """Accessions with no FASTA file are gracefully skipped."""
        session, ds = dataset_with_accessions
        run = create_run(session, dataset_id=ds.dataset_id, stage="pipeline")

        result = execute_detect(
            session, run.run_id, ds.dataset_id, data_dir="/nonexistent"
        )
        assert result.skipped == 3
        assert result.processed == 0


class TestExecuteAnalyze:
    def test_analyze_runs(self, dataset_with_ssrs):
        session, ds, run = dataset_with_ssrs

        # Need metrics first
        execute_metrics(session, run.run_id, ds.dataset_id)

        result = execute_analyze(
            session, run.run_id, ds.dataset_id, force=False
        )
        assert result.status == "completed"
        assert result.processed > 0

    def test_analyze_skips_existing(self, dataset_with_ssrs):
        session, ds, run = dataset_with_ssrs

        execute_metrics(session, run.run_id, ds.dataset_id)
        execute_analyze(session, run.run_id, ds.dataset_id)

        result = execute_analyze(
            session, run.run_id, ds.dataset_id, force=False
        )
        assert result.status == "skipped"

    def test_analyze_dry_run(self, dataset_with_ssrs):
        session, ds, run = dataset_with_ssrs

        execute_metrics(session, run.run_id, ds.dataset_id)

        result = execute_analyze(
            session, run.run_id, ds.dataset_id, dry_run=True
        )
        assert result.status == "dry-run"

        count = session.execute(
            select(func.count()).select_from(StatisticalResult)
        ).scalar()
        assert count == 0


# ---------------------------------------------------------------------------
# Test: Pipeline Orchestrator
# ---------------------------------------------------------------------------


class TestOrchestrator:
    def test_full_run_annotate_through_analyze(self, dataset_with_ssrs):
        """Run annotate → metrics → analyze as a pipeline."""
        session, ds, run_obj = dataset_with_ssrs

        orch = PipelineOrchestrator(session, ds.dataset_id)
        result = orch.run(
            stages=[
                PipelineStage.ANNOTATE,
                PipelineStage.METRICS,
                PipelineStage.ANALYZE,
            ]
        )

        assert result.status == "completed"
        assert len(result.stages) == 3
        assert all(sr.status == "completed" for sr in result.stages)

        # Verify checkpoints were recorded
        completed = get_completed_stages(session, result.run_id)
        assert PipelineStage.ANNOTATE in completed
        assert PipelineStage.METRICS in completed
        assert PipelineStage.ANALYZE in completed

    def test_full_run_with_detect(
        self, dataset_with_accessions, fasta_dir
    ):
        """Run full pipeline from detect through analyze."""
        session, ds = dataset_with_accessions

        orch = PipelineOrchestrator(
            session, ds.dataset_id, data_dir=fasta_dir
        )
        result = orch.run()

        assert result.status == "completed"
        assert len(result.stages) == 4
        for sr in result.stages:
            assert sr.status == "completed"

    def test_dry_run(self, dataset_with_ssrs):
        session, ds, _ = dataset_with_ssrs

        orch = PipelineOrchestrator(session, ds.dataset_id)
        result = orch.run(
            stages=[PipelineStage.ANNOTATE, PipelineStage.METRICS],
            dry_run=True,
        )

        assert result.status == "dry-run"
        for sr in result.stages:
            assert sr.status == "dry-run"

    def test_force_rerun(self, dataset_with_ssrs):
        session, ds, _ = dataset_with_ssrs

        orch = PipelineOrchestrator(session, ds.dataset_id)

        # First run
        r1 = orch.run(stages=[PipelineStage.ANNOTATE, PipelineStage.METRICS])
        assert r1.status == "completed"

        # Second run with force
        r2 = orch.run(
            stages=[PipelineStage.ANNOTATE, PipelineStage.METRICS],
            force=True,
        )
        assert r2.status == "completed"
        # Annotations should have been re-processed
        annot_stage = r2.stages[0]
        assert annot_stage.processed == 3


# ---------------------------------------------------------------------------
# Test: Resume and Recovery
# ---------------------------------------------------------------------------


class TestResume:
    def test_resume_skips_completed(self, dataset_with_ssrs):
        """Resume skips stages that already completed."""
        session, ds, _ = dataset_with_ssrs

        orch = PipelineOrchestrator(session, ds.dataset_id)

        # Run only annotate
        r1 = orch.run(stages=[PipelineStage.ANNOTATE])
        assert r1.status == "completed"

        # Manually mark the run as if it was interrupted
        # Record annotate as completed for this run
        # (already done by orch.run)

        # Resume: annotate should be skipped, rest should run
        r2 = orch.resume(r1.run_id)
        assert r2.status == "completed"

        # Find the annotate stage result — it should be skipped
        annotate_result = next(
            sr for sr in r2.stages if sr.stage == PipelineStage.ANNOTATE
        )
        assert annotate_result.status == "skipped"

        # Metrics should have been processed
        metrics_result = next(
            sr for sr in r2.stages if sr.stage == PipelineStage.METRICS
        )
        assert metrics_result.status == "completed"

    def test_resume_force_reprocesses(self, dataset_with_ssrs):
        """Resume with force re-runs completed stages."""
        session, ds, _ = dataset_with_ssrs

        orch = PipelineOrchestrator(session, ds.dataset_id)

        r1 = orch.run(stages=[PipelineStage.ANNOTATE])
        assert r1.status == "completed"

        r2 = orch.resume(r1.run_id, force=True)
        annotate_result = next(
            sr for sr in r2.stages if sr.stage == PipelineStage.ANNOTATE
        )
        # With force, annotate should be re-processed, not skipped
        assert annotate_result.status == "completed"
        assert annotate_result.processed == 3

    def test_resume_nonexistent_run(self, db_session):
        ds = get_or_create_dataset(db_session, "none", "csv")
        orch = PipelineOrchestrator(db_session, ds.dataset_id)

        with pytest.raises(ValueError, match="Run 9999 not found"):
            orch.resume(9999)

    def test_simulate_interruption_and_resume(self, dataset_with_ssrs):
        """Simulate: annotate completes, metrics 'interrupts', then resume."""
        session, ds, _ = dataset_with_ssrs

        # Step 1: Run only annotate
        orch = PipelineOrchestrator(session, ds.dataset_id)
        r1 = orch.run(stages=[PipelineStage.ANNOTATE])

        # Mark run as failed (simulating interruption during metrics)
        update_run_status(session, r1.run_id, "failed", stage="metrics")

        # Step 2: Resume — annotate should be skipped, rest should run
        r2 = orch.resume(r1.run_id)
        assert r2.status == "completed"

        annotate_sr = next(
            sr for sr in r2.stages if sr.stage == PipelineStage.ANNOTATE
        )
        assert annotate_sr.status == "skipped"

        metrics_sr = next(
            sr for sr in r2.stages if sr.stage == PipelineStage.METRICS
        )
        assert metrics_sr.status == "completed"

        analyze_sr = next(
            sr for sr in r2.stages if sr.stage == PipelineStage.ANALYZE
        )
        assert analyze_sr.status == "completed"


# ---------------------------------------------------------------------------
# Test: Idempotent Execution
# ---------------------------------------------------------------------------


class TestIdempotent:
    def test_no_duplicate_annotations(self, dataset_with_ssrs):
        """Running annotate twice without force does not duplicate data."""
        session, ds, run = dataset_with_ssrs

        execute_annotate(session, run.run_id, ds.dataset_id)
        count1 = session.execute(
            select(func.count()).select_from(SSRAnnotationModel)
        ).scalar()

        execute_annotate(session, run.run_id, ds.dataset_id, force=False)
        count2 = session.execute(
            select(func.count()).select_from(SSRAnnotationModel)
        ).scalar()

        assert count1 == count2

    def test_no_duplicate_metrics(self, dataset_with_ssrs):
        """Running metrics twice without force does not duplicate data."""
        session, ds, run = dataset_with_ssrs

        execute_metrics(session, run.run_id, ds.dataset_id)
        count1 = session.execute(
            select(func.count()).select_from(AccessionMetrics)
        ).scalar()

        execute_metrics(session, run.run_id, ds.dataset_id, force=False)
        count2 = session.execute(
            select(func.count()).select_from(AccessionMetrics)
        ).scalar()

        assert count1 == count2

    def test_no_duplicate_stats(self, dataset_with_ssrs):
        """Running analyze twice without force does not duplicate data."""
        session, ds, run = dataset_with_ssrs

        execute_metrics(session, run.run_id, ds.dataset_id)
        execute_analyze(session, run.run_id, ds.dataset_id)
        count1 = session.execute(
            select(func.count()).select_from(StatisticalResult)
        ).scalar()

        execute_analyze(session, run.run_id, ds.dataset_id, force=False)
        count2 = session.execute(
            select(func.count()).select_from(StatisticalResult)
        ).scalar()

        assert count1 == count2


# ---------------------------------------------------------------------------
# Test: Orchestrator Status
# ---------------------------------------------------------------------------


class TestOrchestratorStatus:
    def test_get_status(self, dataset_with_ssrs):
        session, ds, _ = dataset_with_ssrs

        orch = PipelineOrchestrator(session, ds.dataset_id)
        r = orch.run(stages=[PipelineStage.ANNOTATE])

        status = orch.get_status(r.run_id)
        assert status["run_id"] == r.run_id
        assert status["status"] == "completed"
        assert "annotate" in status["completed_stages"]

    def test_get_status_nonexistent(self, db_session):
        ds = get_or_create_dataset(db_session, "none2", "csv")
        orch = PipelineOrchestrator(db_session, ds.dataset_id)
        status = orch.get_status(99999)
        assert "error" in status


# ---------------------------------------------------------------------------
# Test: CLI run command
# ---------------------------------------------------------------------------


class TestCLIRun:
    def test_run_help(self):
        from click.testing import CliRunner
        from gwico_ssr.cli import cli

        runner = CliRunner()
        result = runner.invoke(cli, ["run", "--help"])
        assert result.exit_code == 0
        assert "--dry-run" in result.output
        assert "--resume" in result.output
        assert "--force" in result.output
        assert "--retry-failed" in result.output

    def test_run_dry_run(self, tmp_path):
        from click.testing import CliRunner
        from gwico_ssr.cli import cli

        db_path = str(tmp_path / "test.db")
        db_url = f"sqlite:///{db_path}"

        runner = CliRunner()

        # Init DB
        runner.invoke(
            cli,
            ["--config", "/nonexistent.toml", "init-db"],
            env={"GWICO_SSR_DB_URL": db_url},
        )

        # Ingest a simple accession
        csv_file = tmp_path / "test.csv"
        csv_file.write_text(
            "Accession,Release_Date,Species,Length,Nuc_Completeness,"
            "Geo_Location,USA,Host,Isolation_Source,Collection_Date\n"
            "NC_045512.2,2020-01-13,SARS-CoV-2,29903,complete,China,,Homo sapiens,,2019-12\n"
        )
        runner.invoke(
            cli,
            [
                "--config", "/nonexistent.toml",
                "ingest", str(csv_file), "--dataset-name", "dryrun_ds",
            ],
            env={"GWICO_SSR_DB_URL": db_url},
        )

        # Run pipeline in dry-run mode
        result = runner.invoke(
            cli,
            [
                "--config", "/nonexistent.toml", "--log-level", "ERROR",
                "run", "dryrun_ds", "--dry-run", "--json-summary",
            ],
            env={"GWICO_SSR_DB_URL": db_url},
        )
        assert result.exit_code == 0
        data = json.loads(result.output)
        assert data["status"] == "dry-run"
        assert len(data["stages"]) == 4

    def test_run_unknown_dataset(self, tmp_path):
        from click.testing import CliRunner
        from gwico_ssr.cli import cli

        db_path = str(tmp_path / "test2.db")
        db_url = f"sqlite:///{db_path}"

        runner = CliRunner()
        runner.invoke(
            cli,
            ["--config", "/nonexistent.toml", "init-db"],
            env={"GWICO_SSR_DB_URL": db_url},
        )

        result = runner.invoke(
            cli,
            [
                "--config", "/nonexistent.toml",
                "run", "nonexistent_dataset",
            ],
            env={"GWICO_SSR_DB_URL": db_url},
        )
        assert "Dataset not found" in result.output
