"""Tests for repository helpers — CRUD, upserts, and deduplication."""

from __future__ import annotations

import pytest
from sqlalchemy import select

from gwico_ssr.db import create_tables, get_engine, get_session_factory, session_scope
from gwico_ssr.db.repository import (
    bulk_upsert_accessions,
    create_run,
    get_accession,
    get_dataset_by_name,
    get_features_for_accession,
    get_or_create_dataset,
    get_run,
    get_ssr_records_for_accession,
    insert_features,
    insert_ssr_annotations,
    insert_ssr_records,
    insert_statistical_result,
    list_accessions,
    update_run_status,
    upsert_accession,
    upsert_accession_metrics,
    upsert_sequence_record,
)
from gwico_ssr.models.schema import Accession, Dataset, Run


@pytest.fixture
def db_session(tmp_path):
    """Provide a database session with all tables created."""
    db_path = tmp_path / "test.db"
    engine = get_engine(f"sqlite:///{db_path}")
    create_tables(engine)
    factory = get_session_factory(engine)
    with session_scope(factory) as session:
        yield session


# ---------------------------------------------------------------------------
# Dataset tests
# ---------------------------------------------------------------------------

def test_get_or_create_dataset(db_session):
    """Create a dataset and retrieve it."""
    ds = get_or_create_dataset(
        db_session, name="test_ds", source_type="csv", organism="SARS-CoV-2"
    )
    assert ds.dataset_id is not None
    assert ds.name == "test_ds"
    assert ds.source_type == "csv"


def test_get_or_create_dataset_idempotent(db_session):
    """Calling get_or_create twice returns the same record."""
    ds1 = get_or_create_dataset(db_session, name="dup", source_type="csv")
    ds2 = get_or_create_dataset(db_session, name="dup", source_type="fasta")
    assert ds1.dataset_id == ds2.dataset_id
    # source_type should NOT be updated since it was found by name
    assert ds2.source_type == "csv"


def test_get_dataset_by_name(db_session):
    """Lookup dataset by name."""
    get_or_create_dataset(db_session, name="lookup_test", source_type="ncbi")
    found = get_dataset_by_name(db_session, "lookup_test")
    assert found is not None
    assert found.name == "lookup_test"


def test_get_dataset_by_name_not_found(db_session):
    """Non-existent dataset returns None."""
    assert get_dataset_by_name(db_session, "nonexistent") is None


# ---------------------------------------------------------------------------
# Run tests
# ---------------------------------------------------------------------------

def test_create_run(db_session):
    """Create a run with correct defaults."""
    ds = get_or_create_dataset(db_session, name="run_ds", source_type="csv")
    run = create_run(db_session, dataset_id=ds.dataset_id, config_hash="abc123")
    assert run.run_id is not None
    assert run.status == "started"
    assert run.pipeline_version is not None
    assert run.hostname is not None


def test_update_run_status(db_session):
    """Run status transitions work correctly."""
    ds = get_or_create_dataset(db_session, name="status_ds", source_type="csv")
    run = create_run(db_session, dataset_id=ds.dataset_id)

    updated = update_run_status(db_session, run.run_id, "completed", stage="export")
    assert updated is not None
    assert updated.status == "completed"
    assert updated.stage == "export"
    assert updated.finished_at is not None


def test_update_run_status_not_found(db_session):
    """Updating a non-existent run returns None."""
    assert update_run_status(db_session, 99999, "failed") is None


def test_get_run(db_session):
    """Get run by ID."""
    ds = get_or_create_dataset(db_session, name="get_run_ds", source_type="csv")
    run = create_run(db_session, dataset_id=ds.dataset_id)
    found = get_run(db_session, run.run_id)
    assert found is not None
    assert found.run_id == run.run_id


# ---------------------------------------------------------------------------
# Accession tests
# ---------------------------------------------------------------------------

def test_upsert_accession_insert(db_session):
    """Insert a new accession."""
    ds = get_or_create_dataset(db_session, name="acc_ds", source_type="csv")
    acc = upsert_accession(
        db_session,
        accession="NC_045512.2",
        dataset_id=ds.dataset_id,
        species="SARS-CoV-2",
        country="China",
        genome_length=29903,
    )
    assert acc.accession == "NC_045512.2"
    assert acc.country == "China"
    assert acc.genome_length == 29903


def test_upsert_accession_update(db_session):
    """Upserting an existing accession updates fields."""
    ds = get_or_create_dataset(db_session, name="upsert_ds", source_type="csv")
    upsert_accession(
        db_session,
        accession="TEST001",
        dataset_id=ds.dataset_id,
        country="USA",
        source_priority=1,
    )
    upsert_accession(
        db_session,
        accession="TEST001",
        dataset_id=ds.dataset_id,
        country="Canada",
        source_priority=2,
    )
    acc = get_accession(db_session, "TEST001")
    assert acc is not None
    # Higher priority overwrites
    assert acc.country == "Canada"


def test_upsert_accession_lower_priority_ignored(db_session):
    """Lower priority does not overwrite existing values."""
    ds = get_or_create_dataset(db_session, name="prio_ds", source_type="csv")
    upsert_accession(
        db_session,
        accession="PRIO001",
        dataset_id=ds.dataset_id,
        country="Japan",
        source_priority=5,
    )
    upsert_accession(
        db_session,
        accession="PRIO001",
        dataset_id=ds.dataset_id,
        country="Korea",
        source_priority=1,
    )
    acc = get_accession(db_session, "PRIO001")
    assert acc is not None
    assert acc.country == "Japan"  # Higher priority preserved


def test_bulk_upsert_accessions(db_session):
    """Bulk upsert processes multiple records."""
    ds = get_or_create_dataset(db_session, name="bulk_ds", source_type="csv")
    records = [
        {"accession": f"BULK{i:03d}", "dataset_id": ds.dataset_id, "country": "US"}
        for i in range(10)
    ]
    count = bulk_upsert_accessions(db_session, records)
    assert count == 10

    found = list_accessions(db_session, dataset_id=ds.dataset_id)
    assert len(found) == 10


def test_list_accessions_with_filters(db_session):
    """List accessions with country filter."""
    ds = get_or_create_dataset(db_session, name="filter_ds", source_type="csv")
    upsert_accession(db_session, accession="F001", dataset_id=ds.dataset_id, country="USA")
    upsert_accession(db_session, accession="F002", dataset_id=ds.dataset_id, country="China")
    upsert_accession(db_session, accession="F003", dataset_id=ds.dataset_id, country="USA")

    usa = list_accessions(db_session, country="USA")
    assert len(usa) == 2

    limited = list_accessions(db_session, limit=1)
    assert len(limited) == 1


# ---------------------------------------------------------------------------
# SequenceRecord tests
# ---------------------------------------------------------------------------

def test_upsert_sequence_record(db_session):
    """Insert and update a sequence record."""
    ds = get_or_create_dataset(db_session, name="seq_ds", source_type="csv")
    upsert_accession(db_session, accession="SEQ001", dataset_id=ds.dataset_id)

    sr = upsert_sequence_record(
        db_session,
        accession="SEQ001",
        download_status="downloaded",
        fasta_path="/data/SEQ001.fasta",
    )
    assert sr.download_status == "downloaded"

    # Update
    sr2 = upsert_sequence_record(
        db_session,
        accession="SEQ001",
        parse_status="parsed",
        sequence_length=29903,
    )
    assert sr2.parse_status == "parsed"
    assert sr2.sequence_length == 29903
    assert sr2.fasta_path == "/data/SEQ001.fasta"  # Preserved


# ---------------------------------------------------------------------------
# FeatureRecord tests
# ---------------------------------------------------------------------------

def test_insert_features(db_session):
    """Bulk insert features."""
    ds = get_or_create_dataset(db_session, name="feat_ds", source_type="csv")
    upsert_accession(db_session, accession="FEAT001", dataset_id=ds.dataset_id)

    features = [
        {"accession": "FEAT001", "feature_type": "CDS", "start": 266, "end": 21555, "gene_name": "ORF1ab", "strand": "+"},
        {"accession": "FEAT001", "feature_type": "CDS", "start": 21563, "end": 25384, "gene_name": "S", "strand": "+"},
    ]
    count = insert_features(db_session, features)
    assert count == 2

    found = get_features_for_accession(db_session, "FEAT001")
    assert len(found) == 2
    gene_names = {f.gene_name for f in found}
    assert "ORF1ab" in gene_names
    assert "S" in gene_names


# ---------------------------------------------------------------------------
# SSRRecord tests
# ---------------------------------------------------------------------------

def test_insert_ssr_records(db_session):
    """Bulk insert SSR records."""
    ds = get_or_create_dataset(db_session, name="ssr_ds", source_type="csv")
    upsert_accession(db_session, accession="SSR001", dataset_id=ds.dataset_id)
    run = create_run(db_session, dataset_id=ds.dataset_id)

    ssrs = [
        {
            "accession": "SSR001",
            "start": 626,
            "end": 635,
            "motif_raw": "AAG",
            "motif_canonical": "AAG",
            "repeat_units": 3,
            "motif_size": 3,
            "repeat_length_bp": 9,
            "strand": "+",
            "actual_repeat": "AAGAAGAAG",
            "run_id": run.run_id,
        },
    ]
    count = insert_ssr_records(db_session, ssrs)
    assert count == 1

    found = get_ssr_records_for_accession(db_session, "SSR001")
    assert len(found) == 1
    assert found[0].motif_raw == "AAG"
    assert found[0].repeat_length_bp == 9


# ---------------------------------------------------------------------------
# SSRAnnotation tests
# ---------------------------------------------------------------------------

def test_insert_ssr_annotations(db_session):
    """Insert SSR annotations linking SSRs to features."""
    ds = get_or_create_dataset(db_session, name="annot_ds", source_type="csv")
    upsert_accession(db_session, accession="ANN001", dataset_id=ds.dataset_id)
    run = create_run(db_session, dataset_id=ds.dataset_id)

    insert_features(db_session, [
        {"accession": "ANN001", "feature_type": "CDS", "start": 100, "end": 1000, "gene_name": "gene1"},
    ])
    features = get_features_for_accession(db_session, "ANN001")

    insert_ssr_records(db_session, [{
        "accession": "ANN001", "start": 200, "end": 209, "motif_raw": "ATG",
        "motif_canonical": "ATG", "repeat_units": 3, "motif_size": 3,
        "repeat_length_bp": 9, "strand": "+", "run_id": run.run_id,
    }])
    ssrs = get_ssr_records_for_accession(db_session, "ANN001")

    count = insert_ssr_annotations(db_session, [{
        "ssr_id": ssrs[0].ssr_id,
        "accession": "ANN001",
        "feature_id": features[0].feature_id,
        "gene_name": "gene1",
        "region_class": "CDS",
        "overlap_bp": 9,
    }])
    assert count == 1


# ---------------------------------------------------------------------------
# AccessionMetrics tests
# ---------------------------------------------------------------------------

def test_upsert_accession_metrics(db_session):
    """Insert and update accession metrics."""
    ds = get_or_create_dataset(db_session, name="met_ds", source_type="csv")
    upsert_accession(db_session, accession="MET001", dataset_id=ds.dataset_id)
    run = create_run(db_session, dataset_id=ds.dataset_id)

    m = upsert_accession_metrics(
        db_session,
        accession="MET001",
        run_id=run.run_id,
        ssr_count_total=57,
        ssr_bp_total=400,
        tri_count=30,
    )
    assert m.ssr_count_total == 57

    # Upsert same (accession, run_id) updates
    m2 = upsert_accession_metrics(
        db_session,
        accession="MET001",
        run_id=run.run_id,
        ssr_count_total=60,
    )
    assert m2.ssr_count_total == 60
    assert m2.tri_count == 30  # Preserved


# ---------------------------------------------------------------------------
# StatisticalResult tests
# ---------------------------------------------------------------------------

def test_insert_statistical_result(db_session):
    """Insert a statistical result."""
    ds = get_or_create_dataset(db_session, name="stats_ds", source_type="csv")
    run = create_run(db_session, dataset_id=ds.dataset_id)

    result = insert_statistical_result(
        db_session,
        run_id=run.run_id,
        analysis_name="chi_square_motif_country",
        test_name="chi_square",
        statistic=784.8,
        p_value=2.9e-114,
        p_value_corrected=3.5e-113,
        n=82,
        metadata_json='{"df": 81}',
    )
    assert result.result_id is not None
    assert result.statistic == 784.8


# ---------------------------------------------------------------------------
# Session scope tests
# ---------------------------------------------------------------------------

def test_session_scope_commit(tmp_path):
    """session_scope commits on clean exit."""
    db_path = tmp_path / "commit.db"
    engine = get_engine(f"sqlite:///{db_path}")
    create_tables(engine)
    factory = get_session_factory(engine)

    with session_scope(factory) as session:
        get_or_create_dataset(session, name="committed", source_type="csv")

    # Verify data persisted in a new session
    with session_scope(factory) as session:
        ds = get_dataset_by_name(session, "committed")
        assert ds is not None


def test_session_scope_rollback(tmp_path):
    """session_scope rolls back on exception."""
    db_path = tmp_path / "rollback.db"
    engine = get_engine(f"sqlite:///{db_path}")
    create_tables(engine)
    factory = get_session_factory(engine)

    with pytest.raises(ValueError):
        with session_scope(factory) as session:
            get_or_create_dataset(session, name="rolled_back", source_type="csv")
            raise ValueError("deliberate error")

    # Verify data was NOT persisted
    with session_scope(factory) as session:
        ds = get_dataset_by_name(session, "rolled_back")
        assert ds is None


# ---------------------------------------------------------------------------
# CLI init-db test
# ---------------------------------------------------------------------------

def test_init_db_cli(tmp_path):
    """init-db CLI command creates the database."""
    from click.testing import CliRunner
    from gwico_ssr.cli import cli

    db_path = tmp_path / "cli_test.db"
    runner = CliRunner()
    result = runner.invoke(cli, [
        "--config", "/nonexistent.toml",
        "init-db",
    ], env={"GWICO_SSR_DB_URL": f"sqlite:///{db_path}"})

    assert result.exit_code == 0, result.output
    assert "Database initialized" in result.output
    assert db_path.exists()
