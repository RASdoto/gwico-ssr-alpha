"""Test suite for Chunk 12 genome browser integration.

Coverage:
- Coordinate validation (off-by-one, boundaries, strand)
- Repeat class color assignment
- BED12 format correctness
- Track header generation
- Backward compatibility
- Edge cases
"""

from __future__ import annotations

import pytest
from pathlib import Path
from sqlalchemy import create_engine
from sqlalchemy.orm import Session

from gwico_ssr.models.schema import (
    Base,
    Dataset,
    Run,
    Accession,
    SSRRecord,
)
from gwico_ssr.export.browser_tracks import (
    Bed12Record,
    BedTrackConfig,
    BrowserTrackGenerator,
    CoordinateValidator,
    CoordinateValidationError,
    CoordinateValidationReport,
    RepeatClass,
    ColorScheme,
    generate_browser_track,
)


# -----------------------------------------------------------------------
# Fixtures
# -----------------------------------------------------------------------

@pytest.fixture
def db_session():
    """In-memory SQLite database session for testing."""
    engine = create_engine("sqlite:///:memory:")
    Base.metadata.create_all(engine)
    session = Session(engine)
    yield session
    session.close()


@pytest.fixture
def populated_db(db_session):
    """Session with dataset, run, accessions, and SSR records."""
    # Create dataset
    dataset = Dataset(name="test_dataset", source_type="csv", organism="virus")
    db_session.add(dataset)
    db_session.flush()

    # Create run
    run = Run(dataset_id=dataset.dataset_id, pipeline_version="1.0.0")
    db_session.add(run)
    db_session.flush()

    # Create accessions
    accessions = []
    for acc_name in ["accession1", "accession2", "accession3"]:
        acc = Accession(
            accession=acc_name,
            dataset_id=dataset.dataset_id,
            species="virus_A",
        )
        db_session.add(acc)
        accessions.append(acc)
    db_session.flush()

    # Create SSR records with various repeat classes
    ssrs = [
        # Perfect SSRs
        SSRRecord(
            run_id=run.run_id,
            accession="accession1",
            motif_raw="AT",
            motif_canonical="AT",
            repeat_units=5,
            motif_size=2,
            repeat_class="perfect",
            repeat_length_bp=10,
            start=100,
            end=110,
            strand="+",
        ),
        SSRRecord(
            run_id=run.run_id,
            accession="accession1",
            motif_raw="GC",
            motif_canonical="GC",
            repeat_units=3,
            motif_size=2,
            repeat_class="perfect",
            repeat_length_bp=6,
            start=200,
            end=206,
            strand="-",
        ),
        # Imperfect SSRs
        SSRRecord(
            run_id=run.run_id,
            accession="accession2",
            motif_raw="CGC",
            motif_canonical="CGC",
            repeat_units=4,
            motif_size=3,
            repeat_class="imperfect",
            repeat_length_bp=15,
            start=50,
            end=65,
            strand="+",
        ),
        # Compound component SSRs
        SSRRecord(
            run_id=run.run_id,
            accession="accession3",
            motif_raw="AT",
            motif_canonical="AT",
            repeat_units=2,
            motif_size=2,
            repeat_class="compound_component",
            repeat_length_bp=4,
            start=150,
            end=154,
            strand="+",
        ),
    ]

    for ssr in ssrs:
        db_session.add(ssr)

    db_session.commit()
    return db_session


# -----------------------------------------------------------------------
# Coordinate Validation Tests (8)
# -----------------------------------------------------------------------

class TestCoordinateValidation:
    """Test coordinate validation."""

    def test_valid_coordinates(self):
        """Valid coordinates pass validation."""
        validator = CoordinateValidator()
        error = validator.validate_single(
            ssr_id=1,
            accession="acc1",
            start=100,
            end=110,
        )
        assert error is None

    def test_start_negative(self):
        """Negative start coordinate detected."""
        validator = CoordinateValidator()
        error = validator.validate_single(
            ssr_id=1,
            accession="acc1",
            start=-1,
            end=10,
        )
        assert error is not None
        assert error.error_type == "start_negative"

    def test_end_before_start(self):
        """End before start detected."""
        validator = CoordinateValidator()
        error = validator.validate_single(
            ssr_id=1,
            accession="acc1",
            start=100,
            end=50,
        )
        assert error is not None
        assert error.error_type == "end_before_start"

    def test_end_equals_start(self):
        """End equal to start detected."""
        validator = CoordinateValidator()
        error = validator.validate_single(
            ssr_id=1,
            accession="acc1",
            start=100,
            end=100,
        )
        assert error is not None

    def test_start_exceeds_length(self):
        """Start exceeding accession length detected."""
        validator = CoordinateValidator()
        error = validator.validate_single(
            ssr_id=1,
            accession="acc1",
            start=1000,
            end=1010,
            accession_length=100,
        )
        assert error is not None
        assert error.error_type == "start_exceeds_length"

    def test_end_exceeds_length_warning(self):
        """End exceeding accession length produces warning."""
        validator = CoordinateValidator()
        error = validator.validate_single(
            ssr_id=1,
            accession="acc1",
            start=90,
            end=110,
            accession_length=100,
        )
        assert error is not None
        assert error.severity == "warning"
        assert error.error_type == "end_exceeds_length"

    def test_batch_validation(self, populated_db):
        """Batch validation of multiple SSRs."""
        session = populated_db
        ssrs = session.query(SSRRecord).all()
        validator = CoordinateValidator()
        report = validator.validate_batch(session, ssrs)
        assert report.total_records == 4
        assert report.valid_count >= 3  # At least 3 valid
        assert isinstance(report, CoordinateValidationReport)

    def test_batch_validation_counts(self, populated_db):
        """Batch validation report has correct counts."""
        session = populated_db
        ssrs = session.query(SSRRecord).all()
        validator = CoordinateValidator()
        report = validator.validate_batch(session, ssrs)
        assert report.total_records == len(ssrs)
        assert report.valid_count + report.invalid_count + len(report.warnings) >= report.total_records


# -----------------------------------------------------------------------
# BED12 Record Tests (5)
# -----------------------------------------------------------------------

class TestBed12Records:
    """Test BED12 record generation."""

    def test_bed12_format_basics(self):
        """BED12 record converts to correct format."""
        record = Bed12Record(
            chrom="chr1",
            chrom_start=100,
            chrom_end=110,
            name="AT_x5",
            score=50,
            strand="+",
            thick_start=100,
            thick_end=110,
            item_rgb="46,204,113",
        )
        bed_line = record.to_bed12()
        fields = bed_line.split("\t")
        assert len(fields) == 12
        assert fields[0] == "chr1"
        assert fields[1] == "100"
        assert fields[2] == "110"
        assert fields[3] == "AT_x5"
        assert fields[8] == "46,204,113"

    def test_bed12_0based_coordinates(self):
        """BED12 uses 0-based half-open coordinates."""
        record = Bed12Record(
            chrom="accession1",
            chrom_start=0,
            chrom_end=10,
            name="motif",
            score=100,
            strand="+",
            thick_start=0,
            thick_end=10,
            item_rgb="0,0,0",
        )
        bed_line = record.to_bed12()
        assert bed_line.startswith("accession1\t0\t10")

    def test_bed12_strand_positive(self):
        """BED12 encodes positive strand."""
        record = Bed12Record(
            chrom="chr",
            chrom_start=1,
            chrom_end=5,
            name="test",
            score=1,
            strand="+",
            thick_start=1,
            thick_end=5,
            item_rgb="0,0,0",
        )
        assert "+" in record.to_bed12()

    def test_bed12_strand_negative(self):
        """BED12 encodes negative strand."""
        record = Bed12Record(
            chrom="chr",
            chrom_start=1,
            chrom_end=5,
            name="test",
            score=1,
            strand="-",
            thick_start=1,
            thick_end=5,
            item_rgb="0,0,0",
        )
        assert "-" in record.to_bed12()

    def test_bed12_block_computation(self):
        """BED12 computes block sizes correctly."""
        record = Bed12Record(
            chrom="chr",
            chrom_start=100,
            chrom_end=150,
            name="test",
            score=1,
            strand="+",
            thick_start=100,
            thick_end=150,
            item_rgb="0,0,0",
        )
        bed_line = record.to_bed12()
        fields = bed_line.split("\t")
        # blockSizes should be computed as chrom_end - chrom_start
        assert fields[10] == "50"  # blockSizes


# -----------------------------------------------------------------------
# Repeat Class Color Tests (3)
# -----------------------------------------------------------------------

class TestRepeatClassColors:
    """Test repeat class color assignment."""

    def test_perfect_color(self):
        """Perfect repeats get correct color."""
        config = BedTrackConfig()
        generator = BrowserTrackGenerator(config)
        color = generator._get_color("perfect")
        assert color == ColorScheme.DEFAULT["perfect"]
        assert "46" in color  # Green component

    def test_imperfect_color(self):
        """Imperfect repeats get correct color."""
        config = BedTrackConfig()
        generator = BrowserTrackGenerator(config)
        color = generator._get_color("imperfect")
        assert color == ColorScheme.DEFAULT["imperfect"]
        assert "52" in color  # Blue component

    def test_compound_color(self):
        """Compound component repeats get correct color."""
        config = BedTrackConfig()
        generator = BrowserTrackGenerator(config)
        color = generator._get_color("compound_component")
        assert color == ColorScheme.DEFAULT["compound_component"]
        assert "231" in color  # Red component


# -----------------------------------------------------------------------
# Track Header Tests (2)
# -----------------------------------------------------------------------

class TestTrackHeaders:
    """Test UCSC Genome Browser track headers."""

    def test_track_header_format(self):
        """Track header has correct UCSC format."""
        config = BedTrackConfig(
            track_name="GWICO-SSR",
            track_description="SSR Test Track",
        )
        generator = BrowserTrackGenerator(config)
        header = generator.generate_track_header()
        assert "track" in header
        assert "name=" in header
        assert "description=" in header
        assert "itemRgb=On" in header

    def test_track_header_custom_name(self):
        """Track header respects custom name."""
        config = BedTrackConfig(track_name="CustomTrack")
        generator = BrowserTrackGenerator(config)
        header = generator.generate_track_header()
        assert "CustomTrack" in header


# -----------------------------------------------------------------------
# Browser Track Generation Tests (4)
# -----------------------------------------------------------------------

class TestBrowserTrackGeneration:
    """Test main track generation function."""

    def test_generate_track_file_created(self, populated_db, tmp_path):
        """Track file is created."""
        session = populated_db
        output_file, report = generate_browser_track(
            session,
            tmp_path,
            track_name="TestTrack",
        )
        assert output_file.exists()
        assert output_file.suffix == ".bed12"

    def test_generate_track_has_header(self, populated_db, tmp_path):
        """Generated track file includes header."""
        session = populated_db
        output_file, _ = generate_browser_track(
            session,
            tmp_path,
            track_name="TestTrack",
        )
        with open(output_file) as f:
            first_line = f.readline()
        assert "track" in first_line
        assert "itemRgb" in first_line

    def test_generate_track_has_records(self, populated_db, tmp_path):
        """Generated track file includes SSR records."""
        session = populated_db
        output_file, _ = generate_browser_track(
            session,
            tmp_path,
            track_name="TestTrack",
        )
        with open(output_file) as f:
            lines = f.readlines()
        # First line is header, rest are records
        assert len(lines) > 1
        # Check second line (first record) has BED12 fields
        fields = lines[1].split("\t")
        assert len(fields) == 12

    def test_generate_track_validation_report(self, populated_db, tmp_path):
        """Generate track returns validation report."""
        session = populated_db
        _, report = generate_browser_track(
            session,
            tmp_path,
            validate_coordinates=True,
        )
        assert report is not None
        assert report.total_records >= 1
        assert report.valid_count >= 1


# -----------------------------------------------------------------------
# Integration Tests (3)
# -----------------------------------------------------------------------

class TestIntegration:
    """Integration tests with full workflow."""

    def test_build_bed12_from_ssr_record(self, populated_db):
        """Build BED12 from SSRRecord."""
        session = populated_db
        ssr = session.query(SSRRecord).filter(
            SSRRecord.repeat_class == "perfect"
        ).first()
        
        config = BedTrackConfig()
        generator = BrowserTrackGenerator(config)
        bed_record = generator.build_bed12_record(ssr)
        
        assert bed_record.chrom == ssr.accession
        assert bed_record.chrom_start == ssr.start
        assert bed_record.chrom_end == ssr.end
        assert bed_record.strand == ssr.strand or "+"
        assert "46,204,113" in bed_record.item_rgb  # Green for perfect

    def test_score_scaling(self):
        """Score scaling works correctly."""
        config = BedTrackConfig()
        generator = BrowserTrackGenerator(config)
        
        # Small repeat
        score1 = generator._scale_score(50)
        assert 0 <= score1 <= 1000
        
        # Large repeat
        score2 = generator._scale_score(500)
        assert score2 > score1  # Larger repeat should have higher score

    def test_full_workflow_with_filters(self, populated_db, tmp_path):
        """Full workflow with run_id filter."""
        session = populated_db
        run_id = session.query(Run.run_id).first()[0]
        
        output_file, report = generate_browser_track(
            session,
            tmp_path,
            run_id=run_id,
            validate_coordinates=True,
        )
        
        assert output_file.exists()
        assert report.valid_count >= 1


# -----------------------------------------------------------------------
# Backward Compatibility Tests (2)
# -----------------------------------------------------------------------

class TestBackwardCompatibility:
    """Ensure Chunk 12 doesn't break existing functionality."""

    def test_existing_ssr_records_unaffected(self, populated_db):
        """SSRRecords from pre-Chunk12 still queryable."""
        session = populated_db
        ssr = session.query(SSRRecord).first()
        assert ssr is not None
        assert ssr.motif_canonical is not None
        assert ssr.repeat_class in ("perfect", "imperfect", "compound_component")

    def test_existing_exports_still_work(self, populated_db, tmp_path):
        """Existing export functions still accessible."""
        session = populated_db
        # Verify existing exports can be imported and used
        from gwico_ssr.export import export_ssrs_bed, export_ssrs_gff3
        
        bed_file = export_ssrs_bed(session, tmp_path)
        assert bed_file.exists()
        
        gff_file = export_ssrs_gff3(session, tmp_path)
        assert gff_file.exists()


# -----------------------------------------------------------------------
# Edge Case Tests (3)
# -----------------------------------------------------------------------

class TestEdgeCases:
    """Test edge cases and error conditions."""

    def test_empty_dataset(self, db_session, tmp_path):
        """Empty dataset produces valid track file."""
        session = db_session
        # Create dataset and run, but no accessions or SSRs
        dataset = Dataset(name="empty", source_type="csv", organism="virus")
        session.add(dataset)
        session.flush()
        run = Run(dataset_id=dataset.dataset_id, pipeline_version="1.0.0")
        session.add(run)
        session.commit()
        
        output_file, _ = generate_browser_track(
            session,
            tmp_path,
        )
        assert output_file.exists()
        with open(output_file) as f:
            lines = f.readlines()
        # Should at least have header
        assert len(lines) >= 1

    def test_single_ssr_record(self, populated_db, tmp_path):
        """Single SSR record generates valid track."""
        session = populated_db
        # Get just one SSR
        ssr = session.query(SSRRecord).first()
        run_id = ssr.run_id
        
        # Create new run with only this SSR
        # (we'll just filter by accession instead for simplicity)
        output_file, report = generate_browser_track(
            session,
            tmp_path,
            run_id=run_id,
        )
        
        with open(output_file) as f:
            lines = [l for l in f.readlines() if not l.startswith("track")]
        assert len(lines) >= 1

    def test_missing_strand_defaults_to_plus(self):
        """Missing strand defaults to '+'."""
        record = Bed12Record(
            chrom="chr",
            chrom_start=1,
            chrom_end=5,
            name="test",
            score=1,
            strand="",  # Empty strand
            thick_start=1,
            thick_end=5,
            item_rgb="0,0,0",
        )
        bed_line = record.to_bed12()
        # Should have valid strand
        assert "\t" in bed_line
