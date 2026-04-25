"""Tests for batch-aware parsing (Chunk 4).

Tests the new batch parsing functionality that integrates with
composite handling from Chunk 3.
"""

import tempfile
from pathlib import Path

import pytest
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

from gwico_ssr.db import create_tables
from gwico_ssr.db.repository import (
    get_accession,
    get_or_create_dataset,
)
from gwico_ssr.ingest import split_composite_fasta, split_composite_genbank
from gwico_ssr.parsers.persist import (
    parse_batch_directory,
    parse_composite_artifact,
    parse_dataset_accessions,
)


@pytest.fixture
def fixtures_dir() -> Path:
    """Return path to test fixtures directory."""
    return Path(__file__).parent.parent / "fixtures"


@pytest.fixture
def composite_fasta_path(fixtures_dir: Path) -> Path:
    """Return path to composite FASTA fixture."""
    return fixtures_dir / "multi_record.fasta"


@pytest.fixture
def composite_genbank_path(fixtures_dir: Path) -> Path:
    """Return path to composite GenBank fixture."""
    return fixtures_dir / "multi_record.gb"


@pytest.fixture
def test_db():
    """Create an in-memory test database."""
    engine = create_engine("sqlite:///:memory:")
    create_tables(engine)
    factory = sessionmaker(bind=engine)
    return factory


@pytest.fixture
def test_dataset(test_db):
    """Create a test dataset with accessions."""
    Session = test_db
    with Session() as session:
        ds = get_or_create_dataset(session, "batch_test_dataset", "csv")

        # Create accessions for the composite fixtures
        accessions = ["NC_045512.2", "OL672836.1"]
        from gwico_ssr.db.repository import upsert_accession
        for acc in accessions:
            upsert_accession(
                session,
                accession=acc,
                dataset_id=ds.dataset_id,
                species="SARS-CoV-2",
            )
        session.commit()
        return ds.dataset_id


class TestParseCompositeArtifact:
    """Tests for parse_composite_artifact function."""

    def test_parse_composite_fasta(self, test_db, test_dataset, composite_fasta_path):
        """Test parsing a composite FASTA file."""
        Session = test_db
        with tempfile.TemporaryDirectory() as tmpdir:
            with Session() as session:
                result = parse_composite_artifact(
                    session,
                    composite_fasta_path=composite_fasta_path,
                    normalization_dir=tmpdir,
                )

                assert result.total_accessions == 2
                assert result.fasta_parsed == 2
                assert result.fasta_failed == 0
                assert len(result.errors) == 0

    def test_parse_composite_fasta_creates_normalized(
        self, test_db, test_dataset, composite_fasta_path
    ):
        """Test that parsing composite FASTA creates normalized files."""
        Session = test_db
        with tempfile.TemporaryDirectory() as tmpdir:
            with Session() as session:
                result = parse_composite_artifact(
                    session,
                    composite_fasta_path=composite_fasta_path,
                    normalization_dir=tmpdir,
                )

                # Check that normalized files were created
                normalized_dir = Path(tmpdir)
                assert (normalized_dir / "NC_045512.2.fasta").exists()
                assert (normalized_dir / "OL672836.1.fasta").exists()

    def test_parse_composite_genbank(self, test_db, test_dataset, composite_genbank_path):
        """Test parsing a composite GenBank file."""
        Session = test_db
        with tempfile.TemporaryDirectory() as tmpdir:
            with Session() as session:
                result = parse_composite_artifact(
                    session,
                    composite_genbank_path=composite_genbank_path,
                    normalization_dir=tmpdir,
                )

                assert result.total_accessions == 2
                assert result.genbank_parsed == 2
                assert result.genbank_failed == 0

    def test_parse_composite_both_fasta_and_genbank(
        self,
        test_db,
        test_dataset,
        composite_fasta_path,
        composite_genbank_path,
    ):
        """Test parsing composite FASTA and GenBank together."""
        Session = test_db
        with tempfile.TemporaryDirectory() as tmpdir:
            with Session() as session:
                result = parse_composite_artifact(
                    session,
                    composite_fasta_path=composite_fasta_path,
                    composite_genbank_path=composite_genbank_path,
                    normalization_dir=tmpdir,
                )

                assert result.total_accessions == 2
                assert result.fasta_parsed == 2
                assert result.genbank_parsed == 2
                assert result.features_inserted > 0

    def test_parse_composite_persists_to_db(
        self, test_db, test_dataset, composite_fasta_path
    ):
        """Test that parsing composite artifact updates the database."""
        Session = test_db
        with tempfile.TemporaryDirectory() as tmpdir:
            with Session() as session:
                # Parse the composite file
                result = parse_composite_artifact(
                    session,
                    composite_fasta_path=composite_fasta_path,
                    normalization_dir=tmpdir,
                )

                # Verify database was updated
                acc1 = get_accession(session, "NC_045512.2")
                acc2 = get_accession(session, "OL672836.1")

                assert acc1 is not None
                assert acc1.gc_content is not None
                assert acc1.genome_length is not None

                assert acc2 is not None
                assert acc2.gc_content is not None
                assert acc2.genome_length is not None


class TestParseBatchDirectory:
    """Tests for parse_batch_directory function."""

    def test_parse_fasta_batch(self, test_db, test_dataset, composite_fasta_path):
        """Test parsing batch directory with FASTA files."""
        Session = test_db
        with tempfile.TemporaryDirectory() as tmpdir:
            # Set up normalized batch directory
            batch_dir = Path(tmpdir) / "batch"
            (batch_dir / "fasta").mkdir(parents=True, exist_ok=True)

            # Split composite into normalized batch
            split_result = split_composite_fasta(
                composite_fasta_path,
                batch_dir / "fasta",
            )

            # Parse the batch directory
            with Session() as session:
                result = parse_batch_directory(
                    session,
                    batch_dir=batch_dir,
                    file_type="fasta",
                )

                assert result.total_accessions == 2
                assert result.fasta_parsed == 2

    def test_parse_genbank_batch(self, test_db, test_dataset, composite_genbank_path):
        """Test parsing batch directory with GenBank files."""
        Session = test_db
        with tempfile.TemporaryDirectory() as tmpdir:
            # Set up normalized batch directory
            batch_dir = Path(tmpdir) / "batch"
            (batch_dir / "genbank").mkdir(parents=True, exist_ok=True)

            # Split composite into normalized batch
            split_result = split_composite_genbank(
                composite_genbank_path,
                batch_dir / "genbank",
            )

            # Parse the batch directory
            with Session() as session:
                result = parse_batch_directory(
                    session,
                    batch_dir=batch_dir,
                    file_type="genbank",
                )

                assert result.total_accessions == 2
                assert result.genbank_parsed == 2

    def test_parse_mixed_batch(
        self,
        test_db,
        test_dataset,
        composite_fasta_path,
        composite_genbank_path,
    ):
        """Test parsing batch directory with both FASTA and GenBank."""
        Session = test_db
        with tempfile.TemporaryDirectory() as tmpdir:
            # Set up normalized batch directory
            batch_dir = Path(tmpdir) / "batch"
            (batch_dir / "fasta").mkdir(parents=True, exist_ok=True)
            (batch_dir / "genbank").mkdir(parents=True, exist_ok=True)

            # Split composites
            split_composite_fasta(composite_fasta_path, batch_dir / "fasta")
            split_composite_genbank(composite_genbank_path, batch_dir / "genbank")

            # Parse the batch directory
            with Session() as session:
                result = parse_batch_directory(
                    session,
                    batch_dir=batch_dir,
                    file_type="both",
                )

                assert result.total_accessions == 2
                assert result.fasta_parsed == 2
                assert result.genbank_parsed == 2


class TestBackwardCompatibility:
    """Tests ensuring single-accession parsing still works."""

    def test_parse_dataset_accessions_still_works(
        self, test_db, test_dataset, fixtures_dir
    ):
        """Test that original parse_dataset_accessions function still works."""
        Session = test_db
        with tempfile.TemporaryDirectory() as tmpdir:
            # Set up deterministic file layout
            output_dir = Path(tmpdir)
            (output_dir / "sequences" / "fasta").mkdir(parents=True, exist_ok=True)
            (output_dir / "sequences" / "genbank").mkdir(parents=True, exist_ok=True)

            # Copy fixture files
            sample_fasta = fixtures_dir / "sample.fasta"
            sample_gb = fixtures_dir / "sample.gb"

            import shutil
            shutil.copy(sample_fasta, output_dir / "sequences" / "fasta" / "sample.fasta")
            shutil.copy(sample_gb, output_dir / "sequences" / "genbank" / "sample.gb")

            # Create accession in dataset
            with Session() as session:
                upsert_accession = __import__(
                    "gwico_ssr.db.repository", fromlist=["upsert_accession"]
                ).upsert_accession
                upsert_accession(
                    session,
                    accession="sample",
                    dataset_id=test_dataset,
                    species="test",
                )
                session.commit()

            # Parse using original function
            with Session() as session:
                result = parse_dataset_accessions(
                    session,
                    accession_ids=["sample"],
                    data_dir=output_dir,
                )

                # Verify it still works
                assert result.total_accessions == 1
                assert result.fasta_parsed + result.genbank_parsed > 0
