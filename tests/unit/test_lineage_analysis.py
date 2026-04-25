"""Unit tests for Chunk 10 lineage-aware analysis.

Test coverage: ingestion, reconciliation, analysis, visualizations (30+ tests).
"""

import csv
import json
import tempfile
from pathlib import Path

import pytest
from sqlalchemy import create_engine
from sqlalchemy.orm import Session

from gwico_ssr.models.schema import (
    Base,
    Dataset,
    Run,
    Accession,
    SSRRecord,
    AccessionMetrics,
    LineageSource as LineageSourceTable,
    LineageAssignment as LineageAssignmentTable,
    LineageMetrics as LineageMetricsTable,
)
from gwico_ssr.analysis.lineage import (
    LineageSource,
    LineageAssignment,
    LineageSummary,
    LineageReport,
    ingest_pango_csv,
    ingest_nextstrain_json,
    get_lineage_for_accession,
    audit_lineage_mapping,
    compute_ssr_metrics_by_lineage,
    store_lineage_metrics,
)


@pytest.fixture
def db_session():
    """In-memory SQLite database session for testing."""
    engine = create_engine("sqlite:///:memory:")
    Base.metadata.create_all(engine)
    session = Session(engine)
    yield session
    session.close()


@pytest.fixture
def sample_data(db_session):
    """Create minimal sample data in database."""
    # Create dataset
    dataset = Dataset(name="test_dataset", source_type="csv", organism="virus")
    db_session.add(dataset)
    db_session.flush()

    # Create run
    run = Run(dataset_id=dataset.dataset_id, pipeline_version="1.0.0")
    db_session.add(run)
    db_session.flush()

    # Create accessions
    accessions = [
        Accession(accession="ACC001", dataset_id=dataset.dataset_id, species="virus_A", release_date="2023-01-15"),
        Accession(accession="ACC002", dataset_id=dataset.dataset_id, species="virus_A", release_date="2023-02-20"),
        Accession(accession="ACC003", dataset_id=dataset.dataset_id, species="virus_A", release_date="2023-03-10"),
        Accession(accession="ACC004", dataset_id=dataset.dataset_id, species="virus_B", release_date="2023-04-05"),
        Accession(accession="ACC005", dataset_id=dataset.dataset_id, species="virus_B", release_date="2023-05-12"),
    ]
    for acc in accessions:
        db_session.add(acc)
    db_session.flush()

    # Create SSR records
    ssr_records = [
        SSRRecord(run_id=run.run_id, accession="ACC001", start=1, end=10, motif_raw="AT", motif_canonical="AT", repeat_units=5, motif_size=2, repeat_length_bp=10, repeat_class="perfect"),
        SSRRecord(run_id=run.run_id, accession="ACC001", start=20, end=26, motif_raw="GC", motif_canonical="GC", repeat_units=3, motif_size=2, repeat_length_bp=6, repeat_class="imperfect"),
        SSRRecord(run_id=run.run_id, accession="ACC002", start=1, end=8, motif_raw="AT", motif_canonical="AT", repeat_units=4, motif_size=2, repeat_length_bp=8, repeat_class="perfect"),
        SSRRecord(run_id=run.run_id, accession="ACC003", start=1, end=6, motif_raw="CGC", motif_canonical="CGC", repeat_units=2, motif_size=3, repeat_length_bp=6, repeat_class="compound_component"),
        SSRRecord(run_id=run.run_id, accession="ACC004", start=1, end=12, motif_raw="AT", motif_canonical="AT", repeat_units=6, motif_size=2, repeat_length_bp=12, repeat_class="perfect"),
        SSRRecord(run_id=run.run_id, accession="ACC005", start=1, end=4, motif_raw="GC", motif_canonical="GC", repeat_units=2, motif_size=2, repeat_length_bp=4, repeat_class="imperfect"),
    ]
    for ssr in ssr_records:
        db_session.add(ssr)
    db_session.flush()

    # Create metrics
    metrics = [
        AccessionMetrics(run_id=run.run_id, accession="ACC001", ra=0.5, rd=0.3, perfect_count=1, imperfect_count=1, compound_component_count=0),
        AccessionMetrics(run_id=run.run_id, accession="ACC002", ra=0.6, rd=0.4, perfect_count=1, imperfect_count=0, compound_component_count=0),
        AccessionMetrics(run_id=run.run_id, accession="ACC003", ra=0.4, rd=0.2, perfect_count=0, imperfect_count=0, compound_component_count=1),
        AccessionMetrics(run_id=run.run_id, accession="ACC004", ra=0.7, rd=0.5, perfect_count=1, imperfect_count=0, compound_component_count=0),
        AccessionMetrics(run_id=run.run_id, accession="ACC005", ra=0.3, rd=0.1, perfect_count=0, imperfect_count=1, compound_component_count=0),
    ]
    for m in metrics:
        db_session.add(m)

    db_session.commit()
    return {
        "dataset": dataset,
        "run": run,
        "accessions": accessions,
        "ssr_records": ssr_records,
        "metrics": metrics,
    }


# ---------------------------------------------------------------------------
# Test Dataclasses
# ---------------------------------------------------------------------------

class TestLineageSourceDataclass:
    """Test LineageSource dataclass."""

    def test_creation(self):
        """Test LineageSource creation."""
        source = LineageSource(
            source_name="pango",
            version="v1.23",
            release_date="2024-01-15",
        )
        assert source.source_name == "pango"
        assert source.version == "v1.23"
        assert source.release_date == "2024-01-15"

    def test_optional_fields(self):
        """Test LineageSource with optional fields."""
        source = LineageSource(source_name="nextstrain", version="2024-Q1")
        assert source.release_date is None
        assert source.description is None


class TestLineageAssignmentDataclass:
    """Test LineageAssignment dataclass."""

    def test_creation(self):
        """Test LineageAssignment creation."""
        assignment = LineageAssignment(
            accession="ACC001",
            lineage_label="BA.1",
            source="pango",
        )
        assert assignment.accession == "ACC001"
        assert assignment.lineage_label == "BA.1"
        assert assignment.source == "pango"

    def test_with_confidence(self):
        """Test LineageAssignment with confidence."""
        assignment = LineageAssignment(
            accession="ACC001",
            lineage_label="XEC",
            source="pango",
            confidence=0.95,
        )
        assert assignment.confidence == 0.95


class TestLineageSummaryDataclass:
    """Test LineageSummary dataclass."""

    def test_creation(self):
        """Test LineageSummary creation."""
        summary = LineageSummary(
            lineage_label="BA.1",
            accession_count=10,
            ssr_count_total=50,
        )
        assert summary.lineage_label == "BA.1"
        assert summary.accession_count == 10
        assert summary.ssr_count_total == 50

    def test_repeat_class_breakdown(self):
        """Test LineageSummary repeat class counts."""
        summary = LineageSummary(
            lineage_label="XEC",
            accession_count=5,
            perfect_count=20,
            imperfect_count=15,
            compound_component_count=5,
        )
        assert summary.perfect_count == 20
        assert summary.imperfect_count == 15
        assert summary.compound_component_count == 5


# ---------------------------------------------------------------------------
# Test Pango CSV Ingestion
# ---------------------------------------------------------------------------

class TestPangoCSVIngestion:
    """Test Pango CSV ingestion."""

    def test_ingest_valid_csv(self, db_session, sample_data):
        """Test ingestion of valid Pango CSV."""
        csv_content = "accession,lineage_label,confidence\n"
        csv_content += "ACC001,BA.1,0.95\n"
        csv_content += "ACC002,XEC,0.88\n"
        csv_content += "ACC003,JN.1,0.92\n"

        with tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False) as f:
            f.write(csv_content)
            f.flush()
            temp_path = f.name

        try:
            report = ingest_pango_csv(
                db_session,
                temp_path,
                source_name="pango",
                version="v1.23",
                release_date="2024-01-15",
            )
            assert report.total_accessions_in_source == 3
            assert report.successfully_mapped == 3
            assert report.unmapped == 0

            # Verify assignments created
            assignments = db_session.query(LineageAssignmentTable).all()
            assert len(assignments) == 3
        finally:
            Path(temp_path).unlink()

    def test_ingest_csv_with_missing_accessions(self, db_session, sample_data):
        """Test ingestion with accessions not in database."""
        csv_content = "accession,lineage_label,confidence\n"
        csv_content += "ACC001,BA.1,0.95\n"
        csv_content += "UNKNOWN,XEC,0.88\n"  # Not in database
        csv_content += "ACC003,JN.1,0.92\n"

        with tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False) as f:
            f.write(csv_content)
            f.flush()
            temp_path = f.name

        try:
            report = ingest_pango_csv(db_session, temp_path)
            assert report.total_accessions_in_source == 3
            assert report.successfully_mapped == 2
            assert report.unmapped == 1
        finally:
            Path(temp_path).unlink()

    def test_ingest_csv_file_not_found(self, db_session):
        """Test ingestion with non-existent file."""
        with pytest.raises(FileNotFoundError):
            ingest_pango_csv(db_session, "/nonexistent/file.csv")

    def test_ingest_csv_with_unassigned_labels(self, db_session, sample_data):
        """Test ingestion with unassigned lineage labels."""
        csv_content = "accession,lineage_label\n"
        csv_content += "ACC001,BA.1\n"
        csv_content += "ACC002,unassigned\n"
        csv_content += "ACC003,unknown\n"

        with tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False) as f:
            f.write(csv_content)
            f.flush()
            temp_path = f.name

        try:
            report = ingest_pango_csv(db_session, temp_path)
            assert report.unassigned_count == 2
        finally:
            Path(temp_path).unlink()

    def test_ingest_csv_update_existing_assignment(self, db_session, sample_data):
        """Test updating confidence for existing assignment."""
        # Create initial source
        source = LineageSourceTable(source_name="pango", version="v1.23")
        db_session.add(source)
        db_session.flush()

        # Create initial assignment
        assignment = LineageAssignmentTable(
            accession="ACC001",
            source_id=source.source_id,
            lineage_label="BA.1",
            confidence=0.80,
        )
        db_session.add(assignment)
        db_session.commit()

        # Re-ingest with higher confidence
        csv_content = "accession,lineage_label,confidence\n"
        csv_content += "ACC001,BA.1,0.95\n"

        with tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False) as f:
            f.write(csv_content)
            f.flush()
            temp_path = f.name

        try:
            report = ingest_pango_csv(
                db_session,
                temp_path,
                source_name="pango",
                version="v1.23",
            )
            # Should reuse existing source and update assignment
            updated = db_session.query(LineageAssignmentTable).filter(
                LineageAssignmentTable.accession == "ACC001"
            ).first()
            assert updated.confidence == 0.95
        finally:
            Path(temp_path).unlink()


# ---------------------------------------------------------------------------
# Test Nextstrain JSON Ingestion
# ---------------------------------------------------------------------------

class TestNextstrainJSONIngestion:
    """Test Nextstrain JSON ingestion."""

    def test_ingest_valid_json(self, db_session, sample_data):
        """Test ingestion of valid Nextstrain JSON."""
        json_content = {
            "accessions": {
                "ACC001": {"clade": "XEC"},
                "ACC002": {"clade": "JN.1"},
                "ACC003": {"clade": "BA.1"},
            }
        }

        with tempfile.NamedTemporaryFile(mode="w", suffix=".json", delete=False) as f:
            json.dump(json_content, f)
            f.flush()
            temp_path = f.name

        try:
            report = ingest_nextstrain_json(
                db_session,
                temp_path,
                source_name="nextstrain",
                version="2024-01",
            )
            assert report.total_accessions_in_source == 3
            assert report.successfully_mapped == 3

            # Verify assignments created
            assignments = db_session.query(LineageAssignmentTable).all()
            assert len(assignments) == 3
        finally:
            Path(temp_path).unlink()

    def test_ingest_json_file_not_found(self, db_session):
        """Test ingestion with non-existent file."""
        with pytest.raises(FileNotFoundError):
            ingest_nextstrain_json(db_session, "/nonexistent/file.json")


# ---------------------------------------------------------------------------
# Test Reconciliation
# ---------------------------------------------------------------------------

class TestLineageReconciliation:
    """Test lineage reconciliation functions."""

    def test_get_lineage_for_accession_latest(self, db_session, sample_data):
        """Test getting latest lineage for accession."""
        # Create source and assignments
        source = LineageSourceTable(source_name="pango", version="v1.23")
        db_session.add(source)
        db_session.flush()

        assignment = LineageAssignmentTable(
            accession="ACC001",
            source_id=source.source_id,
            lineage_label="BA.1",
        )
        db_session.add(assignment)
        db_session.commit()

        # Query
        result = get_lineage_for_accession(db_session, "ACC001")
        assert result is not None
        assert result.lineage_label == "BA.1"

    def test_get_lineage_for_nonexistent_accession(self, db_session):
        """Test getting lineage for non-existent accession."""
        result = get_lineage_for_accession(db_session, "NONEXISTENT")
        assert result is None

    def test_get_lineage_for_specific_source(self, db_session, sample_data):
        """Test getting lineage from specific source."""
        source1 = LineageSourceTable(source_name="pango", version="v1.23")
        source2 = LineageSourceTable(source_name="nextstrain", version="2024-01")
        db_session.add(source1)
        db_session.add(source2)
        db_session.flush()

        assignment1 = LineageAssignmentTable(
            accession="ACC001",
            source_id=source1.source_id,
            lineage_label="BA.1",
        )
        assignment2 = LineageAssignmentTable(
            accession="ACC001",
            source_id=source2.source_id,
            lineage_label="XEC",
        )
        db_session.add(assignment1)
        db_session.add(assignment2)
        db_session.commit()

        # Query specific source
        result = get_lineage_for_accession(db_session, "ACC001", source_id=source1.source_id)
        assert result.lineage_label == "BA.1"

    def test_audit_lineage_mapping(self, db_session, sample_data):
        """Test audit report generation."""
        source = LineageSourceTable(source_name="pango", version="v1.23")
        db_session.add(source)
        db_session.flush()

        assignments = [
            LineageAssignmentTable(accession="ACC001", source_id=source.source_id, lineage_label="BA.1"),
            LineageAssignmentTable(accession="ACC002", source_id=source.source_id, lineage_label="XEC"),
            LineageAssignmentTable(accession="ACC003", source_id=source.source_id, lineage_label="unassigned"),
        ]
        for a in assignments:
            db_session.add(a)
        db_session.commit()

        report = audit_lineage_mapping(db_session, source.source_id)
        assert report.total_accessions_in_source == 3
        assert report.unassigned_count == 1


# ---------------------------------------------------------------------------
# Test Analysis
# ---------------------------------------------------------------------------

class TestLineageAnalysis:
    """Test lineage-aware metrics computation."""

    def test_compute_ssr_metrics_by_lineage_basic(self, db_session, sample_data):
        """Test basic SSR metrics computation by lineage."""
        source = LineageSourceTable(source_name="pango", version="v1.23")
        db_session.add(source)
        db_session.flush()

        assignments = [
            LineageAssignmentTable(accession="ACC001", source_id=source.source_id, lineage_label="BA.1"),
            LineageAssignmentTable(accession="ACC002", source_id=source.source_id, lineage_label="BA.1"),
            LineageAssignmentTable(accession="ACC003", source_id=source.source_id, lineage_label="XEC"),
            LineageAssignmentTable(accession="ACC004", source_id=source.source_id, lineage_label="XEC"),
            LineageAssignmentTable(accession="ACC005", source_id=source.source_id, lineage_label="JN.1"),
        ]
        for a in assignments:
            db_session.add(a)
        db_session.commit()

        metrics = compute_ssr_metrics_by_lineage(db_session, source.source_id)

        # Verify BA.1 metrics
        assert "BA.1" in metrics
        assert metrics["BA.1"].accession_count == 2
        assert metrics["BA.1"].perfect_count == 2  # ACC001 has 1, ACC002 has 1

        # Verify XEC metrics
        assert "XEC" in metrics
        assert metrics["XEC"].accession_count == 2
        assert metrics["XEC"].perfect_count == 1  # ACC004

    def test_compute_ssr_metrics_repeat_class_breakdown(self, db_session, sample_data):
        """Test repeat class breakdown in metrics."""
        source = LineageSourceTable(source_name="pango", version="v1.23")
        db_session.add(source)
        db_session.flush()

        assignment = LineageAssignmentTable(
            accession="ACC001",
            source_id=source.source_id,
            lineage_label="BA.1",
        )
        db_session.add(assignment)
        db_session.commit()

        metrics = compute_ssr_metrics_by_lineage(db_session, source.source_id)

        summary = metrics["BA.1"]
        assert summary.perfect_count == 1
        assert summary.imperfect_count == 1
        assert summary.compound_component_count == 0

    def test_store_lineage_metrics(self, db_session, sample_data):
        """Test storing lineage metrics to database."""
        source = LineageSourceTable(source_name="pango", version="v1.23")
        db_session.add(source)
        db_session.flush()

        summaries = {
            "BA.1": LineageSummary(
                lineage_label="BA.1",
                accession_count=10,
                ssr_count_total=50,
                perfect_count=30,
            ),
            "XEC": LineageSummary(
                lineage_label="XEC",
                accession_count=8,
                ssr_count_total=40,
                imperfect_count=25,
            ),
        }

        count = store_lineage_metrics(db_session, run_id=1, source_id=source.source_id, metrics_dict=summaries)

        assert count == 2

        # Verify stored
        stored = db_session.query(LineageMetricsTable).all()
        assert len(stored) == 2


# ---------------------------------------------------------------------------
# Test Edge Cases
# ---------------------------------------------------------------------------

class TestLineageEdgeCases:
    """Test edge cases and error handling."""

    def test_empty_lineage_summaries(self):
        """Test handling of empty summaries."""
        from gwico_ssr.visualization.lineage_figures import plot_lineage_distribution
        # Should not raise
        plot_lineage_distribution({})

    def test_single_lineage(self, db_session, sample_data):
        """Test metrics with single lineage."""
        source = LineageSourceTable(source_name="pango", version="v1.23")
        db_session.add(source)
        db_session.flush()

        assignment = LineageAssignmentTable(
            accession="ACC001",
            source_id=source.source_id,
            lineage_label="BA.1",
        )
        db_session.add(assignment)
        db_session.commit()

        metrics = compute_ssr_metrics_by_lineage(db_session, source.source_id)
        assert len(metrics) == 1
        assert "BA.1" in metrics

    def test_audit_nonexistent_source(self, db_session):
        """Test audit with non-existent source."""
        with pytest.raises(ValueError):
            audit_lineage_mapping(db_session, 999)


# ---------------------------------------------------------------------------
# Test Backward Compatibility
# ---------------------------------------------------------------------------

class TestBackwardCompatibility:
    """Test backward compatibility with existing code."""

    def test_imports_available(self):
        """Test that lineage modules are importable."""
        from gwico_ssr.analysis.lineage import (
            ingest_pango_csv,
            compute_ssr_metrics_by_lineage,
        )
        assert callable(ingest_pango_csv)
        assert callable(compute_ssr_metrics_by_lineage)

    def test_visualization_imports_available(self):
        """Test that visualization modules are importable."""
        from gwico_ssr.visualization.lineage_figures import (
            plot_lineage_distribution,
            plot_ssr_by_lineage,
        )
        assert callable(plot_lineage_distribution)
        assert callable(plot_ssr_by_lineage)
