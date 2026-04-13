"""Tests for ORM models and schema creation."""

from __future__ import annotations

from sqlalchemy import inspect, text

from gwico_ssr.db import create_tables, get_engine
from gwico_ssr.models.schema import Base


def test_create_all_tables(tmp_path):
    """All core tables are created with correct names."""
    db_path = tmp_path / "test.db"
    engine = get_engine(f"sqlite:///{db_path}", echo=False)
    create_tables(engine)

    inspector = inspect(engine)
    table_names = set(inspector.get_table_names())

    expected = {
        "datasets",
        "runs",
        "accessions",
        "sequence_records",
        "feature_records",
        "ssr_records",
        "ssr_annotations",
        "accession_metrics",
        "statistical_results",
        "stage_checkpoints",
        "failed_accessions",
    }
    assert expected == table_names


def test_indexes_exist(tmp_path):
    """Key indexes are created on the tables."""
    db_path = tmp_path / "test.db"
    engine = get_engine(f"sqlite:///{db_path}")
    create_tables(engine)

    inspector = inspect(engine)

    # Accession indexes
    acc_indexes = {idx["name"] for idx in inspector.get_indexes("accessions")}
    assert "ix_accessions_country" in acc_indexes
    assert "ix_accessions_collection_date" in acc_indexes

    # SSR indexes
    ssr_indexes = {idx["name"] for idx in inspector.get_indexes("ssr_records")}
    assert "ix_ssr_motif_canonical" in ssr_indexes
    assert "ix_ssr_accession_start_end" in ssr_indexes
    assert "ix_ssr_run_id" in ssr_indexes

    # Feature indexes
    feat_indexes = {idx["name"] for idx in inspector.get_indexes("feature_records")}
    assert "ix_features_accession_start_end" in feat_indexes

    # Metrics unique constraint
    unique_constraints = inspector.get_unique_constraints("accession_metrics")
    uq_names = {uc["name"] for uc in unique_constraints}
    assert "uq_accession_metrics_acc_run" in uq_names


def test_foreign_keys_enabled(tmp_path):
    """SQLite foreign keys are enabled via engine pragmas."""
    db_path = tmp_path / "test.db"
    engine = get_engine(f"sqlite:///{db_path}")
    create_tables(engine)

    with engine.connect() as conn:
        result = conn.execute(text("PRAGMA foreign_keys")).scalar()
        assert result == 1


def test_wal_mode(tmp_path):
    """SQLite WAL journal mode is set."""
    db_path = tmp_path / "test.db"
    engine = get_engine(f"sqlite:///{db_path}")
    create_tables(engine)

    with engine.connect() as conn:
        result = conn.execute(text("PRAGMA journal_mode")).scalar()
        assert result == "wal"


def test_idempotent_create(tmp_path):
    """Calling create_tables twice does not error."""
    db_path = tmp_path / "test.db"
    engine = get_engine(f"sqlite:///{db_path}")
    create_tables(engine)
    create_tables(engine)  # Should not raise

    inspector = inspect(engine)
    assert "datasets" in inspector.get_table_names()


def test_column_counts(tmp_path):
    """Tables have the expected number of columns."""
    db_path = tmp_path / "test.db"
    engine = get_engine(f"sqlite:///{db_path}")
    create_tables(engine)

    inspector = inspect(engine)
    columns = {
        "datasets": 7,
        "runs": 10,
        "accessions": 13,
        "sequence_records": 9,
        "feature_records": 10,
        "ssr_records": 13,
        "ssr_annotations": 7,
        "accession_metrics": 14,  # id + accession + run_id + 11 metric fields
        "statistical_results": 14,
    }
    for table, expected_count in columns.items():
        actual = len(inspector.get_columns(table))
        assert actual == expected_count, f"{table} has {actual} columns, expected {expected_count}"
