"""Tests for the GWICO-SSR ingestion module (Chunk 2).

Covers normalizers, CSV parsing, validation, duplicate handling,
local file registration, and CLI ingest command.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest
from click.testing import CliRunner
from sqlalchemy import create_engine, event
from sqlalchemy.orm import Session, sessionmaker

from gwico_ssr.cli import cli
from gwico_ssr.db import create_tables, get_session_factory, session_scope
from gwico_ssr.db.repository import get_accession, get_or_create_dataset, list_accessions
from gwico_ssr.ingest.csv_loader import (
    IngestSummary,
    RowError,
    load_csv,
    parse_csv_row,
    register_local_files,
    validate_header,
)
from gwico_ssr.ingest.normalizers import (
    normalize_accession,
    normalize_completeness,
    normalize_date,
    normalize_geo_location,
    normalize_length,
)
from gwico_ssr.models.schema import Base


# ---- Fixtures ----

FIXTURES_DIR = Path(__file__).resolve().parent.parent / "fixtures"


@pytest.fixture
def db_session(tmp_path):
    """Create an in-memory SQLite DB with all tables and yield a session."""
    url = "sqlite:///:memory:"
    engine = create_engine(url)

    @event.listens_for(engine, "connect")
    def _pragmas(dbapi_conn, _rec):
        cursor = dbapi_conn.cursor()
        cursor.execute("PRAGMA foreign_keys=ON")
        cursor.close()

    Base.metadata.create_all(engine)
    factory = sessionmaker(bind=engine, expire_on_commit=False)
    session = factory()
    try:
        yield session
        session.commit()
    except Exception:
        session.rollback()
        raise
    finally:
        session.close()


@pytest.fixture
def db_engine_and_factory(tmp_path):
    """Create a file-backed SQLite DB and return (engine, factory, db_path)."""
    db_path = tmp_path / "test.db"
    url = f"sqlite:///{db_path}"
    engine = create_engine(url)

    @event.listens_for(engine, "connect")
    def _pragmas(dbapi_conn, _rec):
        cursor = dbapi_conn.cursor()
        cursor.execute("PRAGMA foreign_keys=ON")
        cursor.close()

    Base.metadata.create_all(engine)
    factory = sessionmaker(bind=engine, expire_on_commit=False)
    return engine, factory, db_path


# ==========================================================================
# Normalizer tests
# ==========================================================================


class TestNormalizeDate:
    def test_iso_full(self):
        assert normalize_date("2020-01-13T00:00:00Z") == "2020-01-13"

    def test_iso_no_time(self):
        assert normalize_date("2020-01-13") == "2020-01-13"

    def test_year_month(self):
        assert normalize_date("2019-12") == "2019-12"

    def test_year_only(self):
        assert normalize_date("2020") == "2020"

    def test_empty(self):
        assert normalize_date("") is None
        assert normalize_date(None) is None

    def test_garbage(self):
        assert normalize_date("not-a-date") is None

    def test_whitespace(self):
        assert normalize_date("  2020-01-13  ") == "2020-01-13"


class TestNormalizeCompleteness:
    def test_complete(self):
        assert normalize_completeness("complete") is True

    def test_complete_genome(self):
        assert normalize_completeness("complete genome") is True

    def test_partial(self):
        assert normalize_completeness("partial") is False

    def test_empty(self):
        assert normalize_completeness("") is None
        assert normalize_completeness(None) is None

    def test_case_insensitive(self):
        assert normalize_completeness("COMPLETE") is True
        assert normalize_completeness("Complete Genome") is True


class TestNormalizeGeoLocation:
    def test_country_only(self):
        assert normalize_geo_location("China") == ("China", None)

    def test_country_and_region(self):
        assert normalize_geo_location("USA: California") == ("USA", "California")

    def test_country_with_spaces(self):
        assert normalize_geo_location("South Korea: Seoul") == ("South Korea", "Seoul")

    def test_empty(self):
        assert normalize_geo_location("") == (None, None)
        assert normalize_geo_location(None) == (None, None)

    def test_colon_only(self):
        c, r = normalize_geo_location(":")
        assert c is None
        assert r is None


class TestNormalizeLength:
    def test_valid_int(self):
        assert normalize_length("29903") == 29903

    def test_negative(self):
        assert normalize_length("-5") is None

    def test_zero(self):
        assert normalize_length("0") is None

    def test_non_numeric(self):
        assert normalize_length("abc") is None

    def test_empty(self):
        assert normalize_length("") is None
        assert normalize_length(None) is None


class TestNormalizeAccession:
    def test_valid(self):
        assert normalize_accession("NC_045512.2") == "NC_045512.2"

    def test_stripped(self):
        assert normalize_accession("  OQ123456.1  ") == "OQ123456.1"

    def test_invalid_chars(self):
        assert normalize_accession("@INVALID!") is None

    def test_empty(self):
        assert normalize_accession("") is None
        assert normalize_accession(None) is None


# ==========================================================================
# Header validation tests
# ==========================================================================


class TestValidateHeader:
    def test_full_header(self):
        header = [
            "Accession", "Release_Date", "Species", "Length",
            "Nuc_Completeness", "Geo_Location", "USA", "Host",
            "Isolation_Source", "Collection_Date",
        ]
        col_map, warnings = validate_header(header)
        assert "accession" in col_map
        assert col_map["accession"] == 0
        assert len(warnings) == 0

    def test_missing_required(self):
        header = ["Release_Date", "Species"]
        with pytest.raises(ValueError, match="Missing required columns.*accession"):
            validate_header(header)

    def test_extra_columns(self):
        header = ["Accession", "Extra_Column"]
        col_map, warnings = validate_header(header)
        assert "accession" in col_map
        assert any("Unexpected" in w for w in warnings)

    def test_missing_optional(self):
        header = ["Accession"]
        col_map, warnings = validate_header(header)
        assert "accession" in col_map
        assert any("Optional" in w for w in warnings)


# ==========================================================================
# Row parsing tests
# ==========================================================================


class TestParseRow:
    def test_valid_full_row(self):
        header = [
            "Accession", "Release_Date", "Species", "Length",
            "Nuc_Completeness", "Geo_Location", "USA", "Host",
            "Isolation_Source", "Collection_Date",
        ]
        col_map, _ = validate_header(header)
        row = [
            "NC_045512.2", "2020-01-13T00:00:00Z",
            "SARS-CoV-2", "29903", "complete",
            "China", "", "Homo sapiens", "", "2019-12",
        ]
        record, errors = parse_csv_row(row, col_map, 2)
        assert record is not None
        assert record["accession"] == "NC_045512.2"
        assert record["release_date"] == "2020-01-13"
        assert record["genome_length"] == 29903
        assert record["is_complete"] is True
        assert record["country"] == "China"
        assert record["collection_date"] == "2019-12"
        # Non-critical errors (warnings) for bad dates are in errors list,
        # but row is still valid
        critical = [e for e in errors if e.field == "accession"]
        assert len(critical) == 0

    def test_missing_accession(self):
        col_map = {"accession": 0, "species": 1}
        row = ["", "SARS-CoV-2"]
        record, errors = parse_csv_row(row, col_map, 3)
        assert record is None
        assert any(e.field == "accession" for e in errors)

    def test_invalid_accession(self):
        col_map = {"accession": 0}
        row = ["@BAD!"]
        record, errors = parse_csv_row(row, col_map, 4)
        assert record is None

    def test_bad_date_still_valid_row(self):
        col_map = {"accession": 0, "release_date": 1}
        row = ["OQ111111.1", "not-a-date"]
        record, errors = parse_csv_row(row, col_map, 5)
        assert record is not None
        assert record["release_date"] is None
        assert any(e.field == "release_date" for e in errors)

    def test_usa_column_fills_region(self):
        col_map = {"accession": 0, "geo_location": 1, "usa": 2}
        row = ["OQ111111.1", "USA", "Texas"]
        record, errors = parse_csv_row(row, col_map, 6)
        assert record is not None
        assert record["country"] == "USA"
        assert record["region"] == "Texas"


# ==========================================================================
# CSV loading tests (file-level)
# ==========================================================================


class TestLoadCSV:
    def test_load_valid_csv(self, db_session):
        csv_path = FIXTURES_DIR / "sample_metadata.csv"
        summary = load_csv(csv_path, db_session)
        assert summary.total_rows == 5
        assert summary.valid_rows == 5
        assert summary.invalid_rows == 0
        assert summary.upserted_rows == 5

    def test_accessions_persisted(self, db_session):
        csv_path = FIXTURES_DIR / "sample_metadata.csv"
        load_csv(csv_path, db_session)
        acc = get_accession(db_session, "NC_045512.2")
        assert acc is not None
        assert acc.country == "China"
        assert acc.genome_length == 29903
        assert acc.is_complete is True

    def test_usa_region_populated(self, db_session):
        csv_path = FIXTURES_DIR / "sample_metadata.csv"
        load_csv(csv_path, db_session)
        acc = get_accession(db_session, "OQ123456.1")
        assert acc is not None
        assert acc.country == "USA"
        assert acc.region == "California"

    def test_malformed_csv(self, db_session):
        csv_path = FIXTURES_DIR / "malformed_metadata.csv"
        summary = load_csv(csv_path, db_session)
        # Row 2: missing accession → invalid
        # Row 3: bad date + bad length → valid (non-critical errors)
        # Row 4: bad date → valid
        # Row 5: @INVALID! accession → invalid
        # Row 6: valid
        # Row 7: duplicate of row 6 → skipped
        # Rows 8-9: empty → skipped
        assert summary.invalid_rows == 2  # missing accession + invalid accession
        assert summary.duplicate_rows == 1
        assert len(summary.errors) > 0

    def test_duplicate_within_file(self, db_session):
        csv_path = FIXTURES_DIR / "malformed_metadata.csv"
        summary = load_csv(csv_path, db_session)
        assert summary.duplicate_rows == 1
        assert any("duplicate" in w.lower() for w in summary.warnings)

    def test_file_not_found(self, db_session):
        summary = load_csv(Path("/nonexistent/file.csv"), db_session)
        assert len(summary.errors) > 0
        assert summary.errors[0].field == "file"

    def test_empty_csv(self, db_session, tmp_path):
        empty_csv = tmp_path / "empty.csv"
        empty_csv.write_text("", encoding="utf-8")
        summary = load_csv(empty_csv, db_session)
        assert len(summary.errors) > 0
        assert summary.errors[0].field == "file"

    def test_header_only_csv(self, db_session, tmp_path):
        csv_file = tmp_path / "header_only.csv"
        csv_file.write_text(
            "Accession,Release_Date,Species,Length,Nuc_Completeness,"
            "Geo_Location,USA,Host,Isolation_Source,Collection_Date\n",
            encoding="utf-8",
        )
        summary = load_csv(csv_file, db_session)
        assert summary.total_rows == 0
        assert summary.valid_rows == 0

    def test_missing_required_header(self, db_session, tmp_path):
        csv_file = tmp_path / "bad_header.csv"
        csv_file.write_text("Species,Length\nSARS,29903\n", encoding="utf-8")
        summary = load_csv(csv_file, db_session)
        assert len(summary.errors) > 0
        assert summary.errors[0].field == "header"

    def test_dataset_created(self, db_session):
        csv_path = FIXTURES_DIR / "sample_metadata.csv"
        load_csv(csv_path, db_session, dataset_name="test_dataset")
        from gwico_ssr.db.repository import get_dataset_by_name
        ds = get_dataset_by_name(db_session, "test_dataset")
        assert ds is not None
        assert ds.source_type == "csv"
        assert ds.input_manifest_hash is not None

    def test_idempotent_reimport(self, db_session):
        csv_path = FIXTURES_DIR / "sample_metadata.csv"
        s1 = load_csv(csv_path, db_session, dataset_name="ds1")
        s2 = load_csv(csv_path, db_session, dataset_name="ds1")
        # Second import should upsert same accessions
        assert s1.upserted_rows == 5
        assert s2.upserted_rows == 5
        # Total in DB should still be 5 (upserts, not inserts)
        accs = list_accessions(db_session, dataset_id=1)
        assert len(accs) == 5

    def test_source_priority_upsert(self, db_session):
        csv_path = FIXTURES_DIR / "sample_metadata.csv"
        load_csv(csv_path, db_session, dataset_name="ds1", source_priority=1)
        acc = get_accession(db_session, "NC_045512.2")
        assert acc.country == "China"

        # Create a CSV with higher priority that changes the country
        import tempfile
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".csv", delete=False, encoding="utf-8"
        ) as f:
            f.write("Accession,Geo_Location\nNC_045512.2,Taiwan\n")
            override_path = Path(f.name)

        try:
            load_csv(override_path, db_session, dataset_name="ds2", source_priority=5)
            db_session.expire_all()
            acc = get_accession(db_session, "NC_045512.2")
            assert acc.country == "Taiwan"
        finally:
            override_path.unlink(missing_ok=True)

    def test_summary_to_dict(self, db_session):
        csv_path = FIXTURES_DIR / "sample_metadata.csv"
        summary = load_csv(csv_path, db_session)
        d = summary.to_dict()
        assert isinstance(d, dict)
        assert d["total_rows"] == 5
        assert d["valid_rows"] == 5
        assert "errors" in d
        assert "warnings" in d


# ==========================================================================
# Local file registration tests
# ==========================================================================


class TestRegisterLocalFiles:
    def test_register_fasta(self, db_session, tmp_path):
        # Create an accession first
        ds = get_or_create_dataset(db_session, name="test", source_type="csv")
        from gwico_ssr.db.repository import upsert_accession
        upsert_accession(
            db_session,
            accession="NC_045512.2",
            dataset_id=ds.dataset_id,
        )

        fasta_file = tmp_path / "NC_045512.2.fasta"
        fasta_file.write_text(">NC_045512.2\nATGC\n", encoding="utf-8")

        success, errors = register_local_files(
            db_session,
            accession="NC_045512.2",
            fasta_path=fasta_file,
        )
        assert success is True
        assert len(errors) == 0

        from gwico_ssr.db.repository import upsert_sequence_record
        from gwico_ssr.models.schema import SequenceRecord
        sr = db_session.get(SequenceRecord, "NC_045512.2")
        assert sr is not None
        assert sr.sequence_source == "local"
        assert sr.fasta_path is not None

    def test_register_missing_file(self, db_session, tmp_path):
        ds = get_or_create_dataset(db_session, name="test", source_type="csv")
        from gwico_ssr.db.repository import upsert_accession
        upsert_accession(
            db_session,
            accession="NC_045512.2",
            dataset_id=ds.dataset_id,
        )

        success, errors = register_local_files(
            db_session,
            accession="NC_045512.2",
            fasta_path=tmp_path / "nonexistent.fasta",
        )
        assert success is False
        assert len(errors) > 0

    def test_register_genbank(self, db_session, tmp_path):
        ds = get_or_create_dataset(db_session, name="test", source_type="csv")
        from gwico_ssr.db.repository import upsert_accession
        upsert_accession(
            db_session,
            accession="NC_045512.2",
            dataset_id=ds.dataset_id,
        )

        gb_file = tmp_path / "NC_045512.2.gb"
        gb_file.write_text("LOCUS NC_045512\n", encoding="utf-8")

        success, errors = register_local_files(
            db_session,
            accession="NC_045512.2",
            genbank_path=gb_file,
        )
        assert success is True
        from gwico_ssr.models.schema import SequenceRecord
        sr = db_session.get(SequenceRecord, "NC_045512.2")
        assert sr.genbank_path is not None


# ==========================================================================
# CLI ingest command tests
# ==========================================================================


class TestIngestCLI:
    def test_ingest_command_help(self):
        runner = CliRunner()
        result = runner.invoke(cli, ["ingest", "--help"])
        assert result.exit_code == 0
        assert "Ingest accession metadata" in result.output

    def test_ingest_valid_csv(self, tmp_path):
        runner = CliRunner()
        csv_path = FIXTURES_DIR / "sample_metadata.csv"
        db_path = tmp_path / "test.db"
        config = tmp_path / "config.toml"
        config.write_text(
            f'[database]\nurl = "sqlite:///{db_path.as_posix()}"\n'
            "[logging]\nlevel = \"WARNING\"\nformat = \"text\"\nfile = \"\"\n"
            "[ncbi]\nemail = \"test@test.com\"\n"
            "[ssr]\n[output]\ndir = \"outputs\"\n",
            encoding="utf-8",
        )
        result = runner.invoke(cli, ["--config", str(config), "ingest", str(csv_path)])
        assert result.exit_code == 0
        assert "Total rows: 5" in result.output
        assert "Valid: 5" in result.output
        assert "Upserted to DB: 5" in result.output

    def test_ingest_json_summary(self, tmp_path):
        runner = CliRunner()
        csv_path = FIXTURES_DIR / "sample_metadata.csv"
        db_path = tmp_path / "test.db"
        config = tmp_path / "config.toml"
        config.write_text(
            f'[database]\nurl = "sqlite:///{db_path.as_posix()}"\n'
            "[logging]\nlevel = \"WARNING\"\nformat = \"text\"\nfile = \"\"\n"
            "[ncbi]\nemail = \"test@test.com\"\n"
            "[ssr]\n[output]\ndir = \"outputs\"\n",
            encoding="utf-8",
        )
        result = runner.invoke(
            cli, ["--config", str(config), "ingest", "--json-summary", str(csv_path)]
        )
        assert result.exit_code == 0
        data = json.loads(result.output)
        assert data["total_rows"] == 5
        assert data["valid_rows"] == 5

    def test_ingest_malformed_csv(self, tmp_path):
        runner = CliRunner()
        csv_path = FIXTURES_DIR / "malformed_metadata.csv"
        db_path = tmp_path / "test.db"
        config = tmp_path / "config.toml"
        config.write_text(
            f'[database]\nurl = "sqlite:///{db_path.as_posix()}"\n'
            "[logging]\nlevel = \"WARNING\"\nformat = \"text\"\nfile = \"\"\n"
            "[ncbi]\nemail = \"test@test.com\"\n"
            "[ssr]\n[output]\ndir = \"outputs\"\n",
            encoding="utf-8",
        )
        result = runner.invoke(cli, ["--config", str(config), "ingest", str(csv_path)])
        assert result.exit_code == 0
        assert "Invalid:" in result.output
        assert "Errors:" in result.output
