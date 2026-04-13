"""Tests for the GWICO-SSR NCBI acquisition pipeline (Chunk 3).

Covers EntrezClient, download orchestration, retry manifests.
All NCBI network I/O is mocked — no real Entrez calls.
"""

from __future__ import annotations

import json
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest
from click.testing import CliRunner
from sqlalchemy import create_engine, event
from sqlalchemy.orm import sessionmaker

from gwico_ssr.cli import cli
from gwico_ssr.db.repository import (
    get_accession,
    get_or_create_dataset,
    upsert_accession,
    upsert_sequence_record,
)
from gwico_ssr.ingest.downloader import (
    DownloadSummary,
    download_accessions,
    load_retry_manifest,
    write_retry_manifest,
    _fasta_path,
    _genbank_path,
    _validate_fasta,
    _validate_genbank,
)
from gwico_ssr.ingest.entrez_client import (
    EntrezClient,
    EntrezConfig,
    EntrezResult,
    _extract_accession_from_header,
    _split_fasta,
)
from gwico_ssr.models.schema import Base, SequenceRecord


# ---- Fixtures ----

@pytest.fixture
def db_session():
    """In-memory SQLite DB with all tables."""
    engine = create_engine("sqlite:///:memory:")

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
def sample_dataset(db_session):
    """Create a dataset with 3 accessions."""
    ds = get_or_create_dataset(db_session, name="test_ds", source_type="csv")
    for acc_id in ["NC_045512.2", "OQ111111.1", "OQ222222.1"]:
        upsert_accession(db_session, accession=acc_id, dataset_id=ds.dataset_id)
    return ds


SAMPLE_FASTA = """>NC_045512.2 Severe acute respiratory syndrome coronavirus 2
ATTAAAGGTTTATACCTTCCCAGGTAACAAACCAACCAACTTTCGATCTCTTGTAGATCTG
TTCTCTAAACGAACTTTAAAATCTGTGTGGCTGTCACTCGGCTGCATGCTTAGTGCACTCA
"""

SAMPLE_GENBANK = """LOCUS       NC_045512              29903 bp    ss-RNA  linear   VRL 18-JUL-2020
DEFINITION  Severe acute respiratory syndrome coronavirus 2 isolate Wuhan-Hu-1,
            complete genome.
ACCESSION   NC_045512
VERSION     NC_045512.2
FEATURES             Location/Qualifiers
     source          1..29903
                     /organism="Severe acute respiratory syndrome coronavirus 2"
ORIGIN
        1 attaaaggtt tataccttcc caggtaacaa accaaccaac tttcgatctc ttgtagatct
//
"""


# ==========================================================================
# EntrezClient unit tests
# ==========================================================================


class TestEntrezConfig:
    def test_defaults(self):
        config = EntrezConfig()
        assert config.max_retries == 3
        assert config.rate_limit == 3.0
        assert config.batch_size == 500


class TestEntrezClient:
    def test_min_interval(self):
        config = EntrezConfig(rate_limit=10.0)
        client = EntrezClient(config)
        assert client.min_interval == pytest.approx(0.1, abs=0.01)

    def test_min_interval_zero_rate(self):
        config = EntrezConfig(rate_limit=0.0)
        client = EntrezClient(config)
        assert client.min_interval == 0.0

    @patch("gwico_ssr.ingest.entrez_client.Entrez.efetch")
    def test_fetch_fasta_success(self, mock_efetch):
        mock_handle = MagicMock()
        mock_handle.read.return_value = SAMPLE_FASTA
        mock_efetch.return_value = mock_handle

        config = EntrezConfig(rate_limit=0.0, max_retries=1)
        client = EntrezClient(config)
        result = client.fetch_fasta("NC_045512.2")

        assert result.success is True
        assert result.file_type == "fasta"
        assert result.data is not None
        assert result.data.startswith(">")

    @patch("gwico_ssr.ingest.entrez_client.Entrez.efetch")
    def test_fetch_genbank_success(self, mock_efetch):
        mock_handle = MagicMock()
        mock_handle.read.return_value = SAMPLE_GENBANK
        mock_efetch.return_value = mock_handle

        config = EntrezConfig(rate_limit=0.0, max_retries=1)
        client = EntrezClient(config)
        result = client.fetch_genbank("NC_045512.2")

        assert result.success is True
        assert result.file_type == "genbank"
        assert "LOCUS" in result.data

    @patch("gwico_ssr.ingest.entrez_client.Entrez.efetch")
    def test_fetch_fasta_failure_retries(self, mock_efetch):
        mock_efetch.side_effect = Exception("Network error")

        config = EntrezConfig(rate_limit=0.0, max_retries=2)
        client = EntrezClient(config)

        # Patch time.sleep to not actually wait
        with patch("gwico_ssr.ingest.entrez_client.time.sleep"):
            result = client.fetch_fasta("NC_045512.2")

        assert result.success is False
        assert "2 attempts" in result.error
        assert mock_efetch.call_count == 2

    @patch("gwico_ssr.ingest.entrez_client.Entrez.efetch")
    def test_fetch_empty_response(self, mock_efetch):
        mock_handle = MagicMock()
        mock_handle.read.return_value = ""
        mock_efetch.return_value = mock_handle

        config = EntrezConfig(rate_limit=0.0, max_retries=1)
        client = EntrezClient(config)

        with patch("gwico_ssr.ingest.entrez_client.time.sleep"):
            result = client.fetch_fasta("NC_045512.2")

        assert result.success is False

    @patch("gwico_ssr.ingest.entrez_client.Entrez.efetch")
    def test_fetch_bytes_response(self, mock_efetch):
        mock_handle = MagicMock()
        mock_handle.read.return_value = SAMPLE_FASTA.encode("utf-8")
        mock_efetch.return_value = mock_handle

        config = EntrezConfig(rate_limit=0.0, max_retries=1)
        client = EntrezClient(config)
        result = client.fetch_fasta("NC_045512.2")

        assert result.success is True
        assert isinstance(result.data, str)


class TestFastaHelpers:
    def test_split_fasta_single(self):
        records = _split_fasta(SAMPLE_FASTA)
        assert len(records) == 1
        assert records[0].startswith(">NC_045512.2")

    def test_split_fasta_multi(self):
        multi = SAMPLE_FASTA + "\n>OQ111111.1 Another seq\nATGCATGC\n"
        records = _split_fasta(multi)
        assert len(records) == 2

    def test_split_fasta_empty(self):
        assert _split_fasta("") == []

    def test_extract_accession_from_header(self):
        acc = _extract_accession_from_header(
            ">NC_045512.2 Severe acute...",
            ["NC_045512.2", "OQ111111.1"],
        )
        assert acc == "NC_045512.2"

    def test_extract_accession_unknown(self):
        acc = _extract_accession_from_header(
            ">UNKNOWN_ACC Some desc",
            ["NC_045512.2"],
        )
        assert acc == "UNKNOWN_ACC"

    def test_extract_accession_no_header(self):
        assert _extract_accession_from_header("ATGC", []) is None


# ==========================================================================
# Validation helpers
# ==========================================================================


class TestValidation:
    def test_validate_fasta_ok(self):
        assert _validate_fasta(SAMPLE_FASTA) is True

    def test_validate_fasta_bad(self):
        assert _validate_fasta("not fasta") is False
        assert _validate_fasta(">short") is False

    def test_validate_genbank_ok(self):
        assert _validate_genbank(SAMPLE_GENBANK) is True

    def test_validate_genbank_bad(self):
        assert _validate_genbank("not genbank") is False
        assert _validate_genbank("LOCUS short") is False


# ==========================================================================
# Download orchestrator tests
# ==========================================================================


class TestDownloadAccessions:
    def _make_client(self) -> EntrezClient:
        config = EntrezConfig(rate_limit=0.0, max_retries=1)
        return EntrezClient(config)

    @patch.object(EntrezClient, "fetch_fasta")
    @patch.object(EntrezClient, "fetch_genbank")
    def test_download_success(self, mock_gb, mock_fa, db_session, sample_dataset, tmp_path):
        mock_fa.return_value = EntrezResult(
            accession="NC_045512.2", file_type="fasta", success=True, data=SAMPLE_FASTA,
        )
        mock_gb.return_value = EntrezResult(
            accession="NC_045512.2", file_type="genbank", success=True, data=SAMPLE_GENBANK,
        )

        client = self._make_client()
        summary = download_accessions(
            session=db_session,
            client=client,
            accession_ids=["NC_045512.2"],
            output_dir=str(tmp_path),
        )

        assert summary.downloaded_ok == 1
        assert summary.failed == 0

        # Files created
        fasta_file = _fasta_path(tmp_path / "sequences", "NC_045512.2")
        genbank_file = _genbank_path(tmp_path / "sequences", "NC_045512.2")
        assert fasta_file.exists()
        assert genbank_file.exists()

        # DB updated
        sr = db_session.get(SequenceRecord, "NC_045512.2")
        assert sr is not None
        assert sr.download_status == "downloaded"
        assert sr.sequence_source == "ncbi"

    @patch.object(EntrezClient, "fetch_fasta")
    @patch.object(EntrezClient, "fetch_genbank")
    def test_download_failure(self, mock_gb, mock_fa, db_session, sample_dataset, tmp_path):
        mock_fa.return_value = EntrezResult(
            accession="NC_045512.2", file_type="fasta", success=False, error="Network error",
        )
        mock_gb.return_value = EntrezResult(
            accession="NC_045512.2", file_type="genbank", success=True, data=SAMPLE_GENBANK,
        )

        client = self._make_client()
        summary = download_accessions(
            session=db_session,
            client=client,
            accession_ids=["NC_045512.2"],
            output_dir=str(tmp_path),
        )

        assert summary.failed == 1
        assert len(summary.failures) == 1

        sr = db_session.get(SequenceRecord, "NC_045512.2")
        assert sr is not None
        assert sr.download_status == "failed"

    @patch.object(EntrezClient, "fetch_fasta")
    @patch.object(EntrezClient, "fetch_genbank")
    def test_skip_already_downloaded(self, mock_gb, mock_fa, db_session, sample_dataset, tmp_path):
        # Mark as already downloaded
        upsert_sequence_record(
            db_session,
            accession="NC_045512.2",
            sequence_source="ncbi",
            download_status="downloaded",
            fasta_path=str(tmp_path / "sequences" / "fasta" / "NC_045512.2.fasta"),
            genbank_path=str(tmp_path / "sequences" / "genbank" / "NC_045512.2.gb"),
        )
        # Create the files on disk
        fasta_p = tmp_path / "sequences" / "fasta" / "NC_045512.2.fasta"
        genbank_p = tmp_path / "sequences" / "genbank" / "NC_045512.2.gb"
        fasta_p.parent.mkdir(parents=True, exist_ok=True)
        genbank_p.parent.mkdir(parents=True, exist_ok=True)
        fasta_p.write_text(SAMPLE_FASTA, encoding="utf-8")
        genbank_p.write_text(SAMPLE_GENBANK, encoding="utf-8")

        client = self._make_client()
        summary = download_accessions(
            session=db_session,
            client=client,
            accession_ids=["NC_045512.2"],
            output_dir=str(tmp_path),
        )

        assert summary.already_downloaded == 1
        assert summary.downloaded_ok == 0
        mock_fa.assert_not_called()
        mock_gb.assert_not_called()

    @patch.object(EntrezClient, "fetch_fasta")
    @patch.object(EntrezClient, "fetch_genbank")
    def test_force_redownload(self, mock_gb, mock_fa, db_session, sample_dataset, tmp_path):
        # Mark as already downloaded
        upsert_sequence_record(
            db_session,
            accession="NC_045512.2",
            sequence_source="ncbi",
            download_status="downloaded",
            fasta_path=str(tmp_path / "fake.fasta"),
        )

        mock_fa.return_value = EntrezResult(
            accession="NC_045512.2", file_type="fasta", success=True, data=SAMPLE_FASTA,
        )
        mock_gb.return_value = EntrezResult(
            accession="NC_045512.2", file_type="genbank", success=True, data=SAMPLE_GENBANK,
        )

        client = self._make_client()
        summary = download_accessions(
            session=db_session,
            client=client,
            accession_ids=["NC_045512.2"],
            output_dir=str(tmp_path),
            force=True,
        )

        assert summary.downloaded_ok == 1
        assert summary.already_downloaded == 0
        mock_fa.assert_called_once()

    def test_accession_not_in_db(self, db_session, sample_dataset, tmp_path):
        client = self._make_client()
        summary = download_accessions(
            session=db_session,
            client=client,
            accession_ids=["NONEXISTENT.1"],
            output_dir=str(tmp_path),
        )

        assert summary.failed == 1
        assert "not found" in summary.failures[0]["error"].lower()

    @patch.object(EntrezClient, "fetch_fasta")
    def test_fasta_only(self, mock_fa, db_session, sample_dataset, tmp_path):
        mock_fa.return_value = EntrezResult(
            accession="NC_045512.2", file_type="fasta", success=True, data=SAMPLE_FASTA,
        )

        client = self._make_client()
        summary = download_accessions(
            session=db_session,
            client=client,
            accession_ids=["NC_045512.2"],
            output_dir=str(tmp_path),
            file_types=["fasta"],
        )

        assert summary.downloaded_ok == 1
        fasta_file = _fasta_path(tmp_path / "sequences", "NC_045512.2")
        assert fasta_file.exists()

    @patch.object(EntrezClient, "fetch_fasta")
    @patch.object(EntrezClient, "fetch_genbank")
    def test_multiple_accessions(self, mock_gb, mock_fa, db_session, sample_dataset, tmp_path):
        def fa_side_effect(acc):
            return EntrezResult(acc, "fasta", True, f">{acc}\nATGC\nATGCATGC\n")

        def gb_side_effect(acc):
            return EntrezResult(
                acc, "genbank", True,
                f"LOCUS       {acc}              29903 bp\nDEFINITION  test\n//\n",
            )

        mock_fa.side_effect = fa_side_effect
        mock_gb.side_effect = gb_side_effect

        client = self._make_client()
        summary = download_accessions(
            session=db_session,
            client=client,
            accession_ids=["NC_045512.2", "OQ111111.1", "OQ222222.1"],
            output_dir=str(tmp_path),
        )

        assert summary.downloaded_ok == 3
        assert summary.failed == 0


# ==========================================================================
# Retry manifest tests
# ==========================================================================


class TestRetryManifest:
    def test_write_manifest(self, tmp_path):
        summary = DownloadSummary(
            total_requested=3,
            failed=1,
            failures=[{"accession": "OQ111111.1", "error": "Timeout"}],
        )
        path = write_retry_manifest(summary, tmp_path)
        assert path is not None
        assert path.exists()

        data = json.loads(path.read_text(encoding="utf-8"))
        assert data["failed_count"] == 1
        assert "OQ111111.1" in data["accessions"]

    def test_no_manifest_on_success(self, tmp_path):
        summary = DownloadSummary(total_requested=3, downloaded_ok=3)
        path = write_retry_manifest(summary, tmp_path)
        assert path is None

    def test_load_manifest(self, tmp_path):
        manifest = {
            "failed_count": 2,
            "accessions": ["OQ111111.1", "OQ222222.1"],
            "details": [],
        }
        manifest_path = tmp_path / "retry_manifest.json"
        manifest_path.write_text(json.dumps(manifest), encoding="utf-8")

        acc_ids = load_retry_manifest(manifest_path)
        assert acc_ids == ["OQ111111.1", "OQ222222.1"]

    def test_load_missing_manifest(self, tmp_path):
        acc_ids = load_retry_manifest(tmp_path / "nonexistent.json")
        assert acc_ids == []


# ==========================================================================
# Deterministic path tests
# ==========================================================================


class TestDeterministicPaths:
    def test_fasta_path(self, tmp_path):
        p = _fasta_path(tmp_path, "NC_045512.2")
        assert p == tmp_path / "fasta" / "NC_045512.2.fasta"

    def test_genbank_path(self, tmp_path):
        p = _genbank_path(tmp_path, "NC_045512.2")
        assert p == tmp_path / "genbank" / "NC_045512.2.gb"


# ==========================================================================
# DownloadSummary tests
# ==========================================================================


class TestDownloadSummary:
    def test_to_dict(self):
        summary = DownloadSummary(
            total_requested=10,
            already_downloaded=3,
            downloaded_ok=5,
            failed=2,
            failures=[{"accession": "X", "error": "E"}],
        )
        d = summary.to_dict()
        assert d["total_requested"] == 10
        assert d["already_downloaded"] == 3
        assert d["downloaded_ok"] == 5
        assert d["failed"] == 2
        assert len(d["failures"]) == 1


# ==========================================================================
# CLI download command tests
# ==========================================================================


class TestDownloadCLI:
    def test_download_help(self):
        runner = CliRunner()
        result = runner.invoke(cli, ["download", "--help"])
        assert result.exit_code == 0
        assert "Download FASTA/GenBank" in result.output

    @patch.object(EntrezClient, "fetch_fasta")
    @patch.object(EntrezClient, "fetch_genbank")
    def test_download_with_accessions(self, mock_gb, mock_fa, tmp_path):
        mock_fa.return_value = EntrezResult(
            accession="NC_045512.2", file_type="fasta", success=True, data=SAMPLE_FASTA,
        )
        mock_gb.return_value = EntrezResult(
            accession="NC_045512.2", file_type="genbank", success=True, data=SAMPLE_GENBANK,
        )

        runner = CliRunner()
        db_path = tmp_path / "test.db"
        config = tmp_path / "config.toml"
        config.write_text(
            f'[database]\nurl = "sqlite:///{db_path.as_posix()}"\n'
            "[logging]\nlevel = \"WARNING\"\nformat = \"text\"\nfile = \"\"\n"
            "[ncbi]\nemail = \"test@test.com\"\nrate_limit = 0.0\nmax_retries = 1\n"
            "[ssr]\n[output]\n"
            f'dir = "{tmp_path.as_posix()}"\n',
            encoding="utf-8",
        )

        # First ingest some accessions
        fixtures_dir = Path(__file__).resolve().parent.parent / "fixtures"
        csv_path = fixtures_dir / "sample_metadata.csv"
        result = runner.invoke(cli, [
            "--config", str(config),
            "ingest", str(csv_path), "--dataset-name", "test_ds",
        ])
        assert result.exit_code == 0

        # Now download
        result = runner.invoke(cli, [
            "--config", str(config),
            "download", "test_ds", "--accessions", "NC_045512.2",
        ])
        assert result.exit_code == 0
        assert "Downloaded OK: 1" in result.output

    def test_download_dataset_not_found(self, tmp_path):
        runner = CliRunner()
        db_path = tmp_path / "test.db"
        config = tmp_path / "config.toml"
        config.write_text(
            f'[database]\nurl = "sqlite:///{db_path.as_posix()}"\n'
            "[logging]\nlevel = \"WARNING\"\nformat = \"text\"\nfile = \"\"\n"
            "[ncbi]\nemail = \"test@test.com\"\n"
            "[ssr]\n[output]\ndir = \"outputs\"\n",
            encoding="utf-8",
        )
        result = runner.invoke(cli, [
            "--config", str(config),
            "download", "nonexistent_ds",
        ])
        assert "not found" in result.output.lower() or result.exit_code != 0

    @patch.object(EntrezClient, "fetch_fasta")
    @patch.object(EntrezClient, "fetch_genbank")
    def test_download_json_summary(self, mock_gb, mock_fa, tmp_path):
        mock_fa.return_value = EntrezResult(
            accession="NC_045512.2", file_type="fasta", success=True, data=SAMPLE_FASTA,
        )
        mock_gb.return_value = EntrezResult(
            accession="NC_045512.2", file_type="genbank", success=True, data=SAMPLE_GENBANK,
        )

        runner = CliRunner()
        db_path = tmp_path / "test.db"
        config = tmp_path / "config.toml"
        config.write_text(
            f'[database]\nurl = "sqlite:///{db_path.as_posix()}"\n'
            "[logging]\nlevel = \"WARNING\"\nformat = \"text\"\nfile = \"\"\n"
            "[ncbi]\nemail = \"test@test.com\"\nrate_limit = 0.0\nmax_retries = 1\n"
            "[ssr]\n[output]\n"
            f'dir = "{tmp_path.as_posix()}"\n',
            encoding="utf-8",
        )

        fixtures_dir = Path(__file__).resolve().parent.parent / "fixtures"
        csv_path = fixtures_dir / "sample_metadata.csv"
        runner.invoke(cli, [
            "--config", str(config),
            "ingest", str(csv_path), "--dataset-name", "test_ds",
        ])

        result = runner.invoke(cli, [
            "--config", str(config),
            "download", "test_ds",
            "--accessions", "NC_045512.2",
            "--json-summary",
        ])
        assert result.exit_code == 0
        data = json.loads(result.output)
        assert data["downloaded_ok"] == 1
