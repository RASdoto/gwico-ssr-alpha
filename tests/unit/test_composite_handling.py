"""Tests for composite FASTA and GenBank file handling.

Tests the new Chunk 3 functionality for splitting and normalizing
multi-record composite files into individual accession-level records.
"""

import tempfile
from pathlib import Path

import pytest

from gwico_ssr.ingest.normalizers import (
    split_composite_fasta,
    split_composite_genbank,
)
from gwico_ssr.parsers.fasta_parser import parse_fasta
from gwico_ssr.parsers.genbank_parser import parse_genbank_composite


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


class TestSplitCompositeFasta:
    """Tests for split_composite_fasta function."""

    def test_split_composite_fasta_basic(self, composite_fasta_path: Path):
        """Test basic composite FASTA splitting."""
        with tempfile.TemporaryDirectory() as tmpdir:
            result = split_composite_fasta(composite_fasta_path, tmpdir)

            assert result.file_type == "fasta"
            assert result.record_count == 2
            assert result.accession_count == 2
            assert result.checksum_sha256  # Should be non-empty
            assert len(result.normalized_records) == 2
            assert len(result.errors) == 0

            # Check normalized files exist
            for rec in result.normalized_records:
                assert Path(rec["output_path"]).exists()
                assert Path(rec["output_path"]).stat().st_size > 0

    def test_split_composite_fasta_accession_extraction(self, composite_fasta_path: Path):
        """Test that accessions are correctly extracted."""
        with tempfile.TemporaryDirectory() as tmpdir:
            result = split_composite_fasta(composite_fasta_path, tmpdir)

            accessions = [rec["accession"] for rec in result.normalized_records]
            assert "NC_045512.2" in accessions
            assert "OL672836.1" in accessions

    def test_split_composite_fasta_output_files(self, composite_fasta_path: Path):
        """Test that output files can be parsed independently."""
        with tempfile.TemporaryDirectory() as tmpdir:
            result = split_composite_fasta(composite_fasta_path, tmpdir)

            for rec in result.normalized_records:
                output_path = Path(rec["output_path"])
                # Parse the output file
                parse_result = parse_fasta(output_path)
                assert parse_result.success
                assert parse_result.record_count == 1
                assert parse_result.records[0].accession == rec["accession"]

    def test_split_composite_fasta_file_not_found(self):
        """Test error handling when file doesn't exist."""
        result = split_composite_fasta("/nonexistent/file.fasta", "/tmp")
        assert len(result.errors) > 0
        assert "not found" in result.errors[0]["message"].lower()

    def test_split_composite_fasta_empty_file(self):
        """Test error handling for empty file."""
        with tempfile.TemporaryDirectory() as tmpdir:
            empty_file = Path(tmpdir) / "empty.fasta"
            empty_file.write_text("")
            
            result = split_composite_fasta(empty_file, tmpdir)
            assert len(result.errors) > 0
            assert "empty" in result.errors[0]["message"].lower()

    def test_split_composite_fasta_duplicate_handling_error(self, fixtures_dir: Path):
        """Test duplicate handling with error policy."""
        # Create a FASTA file with duplicate accessions
        dup_fasta = fixtures_dir / "duplicate_accessions.fasta"
        dup_fasta.write_text(">NC_045512.2\nACGT\n>NC_045512.2\nGGGG\n")

        try:
            with tempfile.TemporaryDirectory() as tmpdir:
                result = split_composite_fasta(
                    dup_fasta,
                    tmpdir,
                    duplicate_policy="error",
                )
                # Should have an error for duplicate
                assert len(result.errors) > 0
                assert any("duplicate" in e["message"].lower() for e in result.errors)
        finally:
            dup_fasta.unlink()

    def test_split_composite_fasta_duplicate_handling_skip(self, fixtures_dir: Path):
        """Test duplicate handling with skip policy."""
        dup_fasta = fixtures_dir / "duplicate_accessions_skip.fasta"
        dup_fasta.write_text(">NC_045512.2\nACGT\n>NC_045512.2\nGGGG\n")

        try:
            with tempfile.TemporaryDirectory() as tmpdir:
                result = split_composite_fasta(
                    dup_fasta,
                    tmpdir,
                    duplicate_policy="skip",
                )
                # Should have duplicates recorded but only 1 record processed
                assert result.record_count == 2
                assert result.accession_count == 1  # Only one unique
                assert len(result.normalized_records) == 1
                assert len(result.duplicates) == 1
        finally:
            dup_fasta.unlink()

    def test_split_composite_fasta_duplicate_handling_first(self, fixtures_dir: Path):
        """Test duplicate handling with first policy (default)."""
        dup_fasta = fixtures_dir / "duplicate_accessions_first.fasta"
        dup_fasta.write_text(">NC_045512.2\nACGT\n>NC_045512.2\nGGGG\n")

        try:
            with tempfile.TemporaryDirectory() as tmpdir:
                result = split_composite_fasta(
                    dup_fasta,
                    tmpdir,
                    duplicate_policy="first",
                )
                # Should silently keep first, skip second
                assert result.record_count == 2
                assert result.accession_count == 1
                assert len(result.normalized_records) == 1
        finally:
            dup_fasta.unlink()


class TestSplitCompositeGenBank:
    """Tests for split_composite_genbank function."""

    def test_split_composite_genbank_basic(self, composite_genbank_path: Path):
        """Test basic composite GenBank splitting."""
        with tempfile.TemporaryDirectory() as tmpdir:
            result = split_composite_genbank(composite_genbank_path, tmpdir)

            assert result.file_type == "genbank"
            assert result.record_count == 2
            assert result.accession_count == 2
            assert result.checksum_sha256  # Should be non-empty
            assert len(result.normalized_records) == 2
            assert len(result.errors) == 0

            # Check normalized files exist
            for rec in result.normalized_records:
                assert Path(rec["output_path"]).exists()
                assert Path(rec["output_path"]).stat().st_size > 0

    def test_split_composite_genbank_accession_extraction(self, composite_genbank_path: Path):
        """Test that accessions are correctly extracted from GenBank."""
        with tempfile.TemporaryDirectory() as tmpdir:
            result = split_composite_genbank(composite_genbank_path, tmpdir)

            accessions = [rec["accession"] for rec in result.normalized_records]
            # GenBank records include version numbers
            assert any("NC_045512" in acc for acc in accessions)
            assert any("OL672836" in acc for acc in accessions)

    def test_split_composite_genbank_output_files(self, composite_genbank_path: Path):
        """Test that output files can be parsed independently."""
        with tempfile.TemporaryDirectory() as tmpdir:
            result = split_composite_genbank(composite_genbank_path, tmpdir)

            for rec in result.normalized_records:
                output_path = Path(rec["output_path"])
                # Parse the output file
                parse_results = parse_genbank_composite(output_path)
                assert len(parse_results) >= 1
                assert parse_results[0].success
                assert parse_results[0].sequence_info is not None
                assert parse_results[0].sequence_info.accession == rec["accession"]

    def test_split_composite_genbank_file_not_found(self):
        """Test error handling when file doesn't exist."""
        result = split_composite_genbank("/nonexistent/file.gb", "/tmp")
        assert len(result.errors) > 0
        assert "not found" in result.errors[0]["message"].lower()

    def test_split_composite_genbank_empty_file(self):
        """Test error handling for empty file."""
        with tempfile.TemporaryDirectory() as tmpdir:
            empty_file = Path(tmpdir) / "empty.gb"
            empty_file.write_text("")
            
            result = split_composite_genbank(empty_file, tmpdir)
            assert len(result.errors) > 0
            assert "empty" in result.errors[0]["message"].lower()


class TestParseGenbankComposite:
    """Tests for the new parse_genbank_composite function."""

    def test_parse_genbank_composite_basic(self, composite_genbank_path: Path):
        """Test basic multi-record GenBank parsing."""
        results = parse_genbank_composite(composite_genbank_path)

        assert len(results) == 2
        assert results[0].success
        assert results[1].success
        assert results[0].sequence_info is not None
        assert results[1].sequence_info is not None

    def test_parse_genbank_composite_accessions(self, composite_genbank_path: Path):
        """Test that accessions are extracted correctly."""
        results = parse_genbank_composite(composite_genbank_path)

        accessions = [r.sequence_info.accession for r in results if r.sequence_info]
        # GenBank records include version numbers
        assert any("NC_045512" in acc for acc in accessions)
        assert any("OL672836" in acc for acc in accessions)

    def test_parse_genbank_composite_features(self, composite_genbank_path: Path):
        """Test that features are extracted for each record."""
        results = parse_genbank_composite(composite_genbank_path)

        for result in results:
            assert result.success
            # Should have at least one feature (CDS)
            assert result.feature_count > 0
            assert result.cds_count > 0


class TestCompositeIntegration:
    """Integration tests combining normalization and parsing."""

    def test_roundtrip_split_and_reparse_fasta(self, composite_fasta_path: Path):
        """Test that splitting and re-parsing FASTA preserves data."""
        with tempfile.TemporaryDirectory() as tmpdir:
            # Split the composite file
            split_result = split_composite_fasta(composite_fasta_path, tmpdir)
            assert split_result.accession_count == 2

            # Re-parse each normalized file
            for rec in split_result.normalized_records:
                parse_result = parse_fasta(rec["output_path"])
                assert parse_result.success
                assert len(parse_result.records) == 1
                assert parse_result.records[0].accession == rec["accession"]

    def test_roundtrip_split_and_reparse_genbank(self, composite_genbank_path: Path):
        """Test that splitting and re-parsing GenBank preserves data."""
        with tempfile.TemporaryDirectory() as tmpdir:
            # Split the composite file
            split_result = split_composite_genbank(composite_genbank_path, tmpdir)
            assert split_result.accession_count == 2

            # Re-parse each normalized file
            for rec in split_result.normalized_records:
                parse_results = parse_genbank_composite(rec["output_path"])
                assert len(parse_results) >= 1
                assert parse_results[0].success
                assert parse_results[0].sequence_info.accession == rec["accession"]
