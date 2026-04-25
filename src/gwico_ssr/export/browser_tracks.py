"""Genome browser track generation for GWICO-SSR.

Generates browser-ready BED12 tracks with color coding by repeat class,
coordinate validation, and support for UCSC, IGV, and JBrowse formats.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path
from typing import Optional

from sqlalchemy import select
from sqlalchemy.orm import Session

from gwico_ssr.models.schema import Accession, SSRRecord

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Enums and Constants
# ---------------------------------------------------------------------------

class RepeatClass(str, Enum):
    """SSR repeat class enumeration."""

    PERFECT = "perfect"
    IMPERFECT = "imperfect"
    COMPOUND_COMPONENT = "compound_component"


class ColorScheme(dict):
    """Color scheme for repeat classes in RGB hex format."""

    # Based on visualization palette from Chunks 9-11
    DEFAULT = {
        RepeatClass.PERFECT: "46,204,113",  # Green (#2ecc71)
        RepeatClass.IMPERFECT: "52,152,219",  # Blue (#3498db)
        RepeatClass.COMPOUND_COMPONENT: "231,76,60",  # Red (#e74c3c)
    }


# ---------------------------------------------------------------------------
# Dataclasses
# ---------------------------------------------------------------------------

@dataclass
class CoordinateValidationError:
    """Single coordinate validation error."""

    ssr_id: int
    accession: str
    error_type: str  # 'start_negative', 'end_before_start', 'bounds_exceeded', etc.
    detail: str
    severity: str = "error"  # 'error' or 'warning'


@dataclass
class CoordinateValidationReport:
    """Report from bulk coordinate validation."""

    total_records: int
    valid_count: int
    invalid_count: int
    errors: list[CoordinateValidationError] = field(default_factory=list)
    warnings: list[CoordinateValidationError] = field(default_factory=list)

    @property
    def pass_rate(self) -> float:
        """Percentage of valid records."""
        return (self.valid_count / self.total_records * 100) if self.total_records > 0 else 0.0

    @property
    def total_issues(self) -> int:
        """Total errors + warnings."""
        return len(self.errors) + len(self.warnings)


@dataclass
class BedTrackConfig:
    """Configuration for BED track generation."""

    track_name: str = "GWICO-SSR"
    track_description: str = "SSR records with repeat class coloring"
    color_scheme: dict = field(default_factory=lambda: ColorScheme.DEFAULT)
    include_header: bool = True
    use_thick_region: bool = True  # Use thickStart=start, thickEnd=end for full SSR region
    score_field: str = "repeat_length_bp"  # Field to use for score (0-1000 scaled)
    strict_validation: bool = False  # If True, skip invalid coordinates; if False, warn only


@dataclass
class Bed12Record:
    """Single BED12 formatted record."""

    chrom: str
    chrom_start: int  # 0-based
    chrom_end: int  # 0-based, half-open
    name: str
    score: int  # 0-1000
    strand: str  # +/-/.
    thick_start: int  # 0-based
    thick_end: int  # 0-based, half-open
    item_rgb: str  # R,G,B
    block_count: int = 1
    block_sizes: str = "1"  # Will be computed as "chrom_end - chrom_start"
    block_starts: str = "0"

    def to_bed12(self) -> str:
        """Convert to BED12 format line."""
        self.block_sizes = str(self.chrom_end - self.chrom_start)
        return (
            f"{self.chrom}\t{self.chrom_start}\t{self.chrom_end}\t{self.name}\t{self.score}\t"
            f"{self.strand}\t{self.thick_start}\t{self.thick_end}\t{self.item_rgb}\t"
            f"{self.block_count}\t{self.block_sizes}\t{self.block_starts}"
        )


# ---------------------------------------------------------------------------
# Coordinate Validation
# ---------------------------------------------------------------------------

class CoordinateValidator:
    """Validates SSR coordinates against stored ranges and constraints."""

    def __init__(self, strict: bool = False):
        """Initialize validator.

        Args:
            strict: If True, raise on invalid coordinates. If False, collect warnings.
        """
        self.strict = strict

    def validate_single(
        self,
        ssr_id: int,
        accession: str,
        start: int,
        end: int,
        accession_length: Optional[int] = None,
    ) -> Optional[CoordinateValidationError]:
        """Validate single SSR record coordinates.

        Returns:
            CoordinateValidationError if invalid, None if valid.
        """
        errors = []

        # Check basic constraints
        if start < 0:
            return CoordinateValidationError(
                ssr_id=ssr_id,
                accession=accession,
                error_type="start_negative",
                detail=f"start={start} is negative",
            )

        if end <= start:
            return CoordinateValidationError(
                ssr_id=ssr_id,
                accession=accession,
                error_type="end_before_start",
                detail=f"end={end} <= start={start}",
            )

        # Check bounds if accession length provided
        if accession_length is not None:
            if start >= accession_length:
                return CoordinateValidationError(
                    ssr_id=ssr_id,
                    accession=accession,
                    error_type="start_exceeds_length",
                    detail=f"start={start} >= accession_length={accession_length}",
                )

            if end > accession_length:
                return CoordinateValidationError(
                    ssr_id=ssr_id,
                    accession=accession,
                    error_type="end_exceeds_length",
                    detail=f"end={end} > accession_length={accession_length}",
                    severity="warning",
                )

        return None

    def validate_batch(
        self,
        session: Session,
        ssrs: list[SSRRecord],
    ) -> CoordinateValidationReport:
        """Validate a batch of SSR records.

        Args:
            session: SQLAlchemy session
            ssrs: List of SSRRecord objects

        Returns:
            CoordinateValidationReport with results
        """
        report = CoordinateValidationReport(
            total_records=len(ssrs),
            valid_count=0,
            invalid_count=0,
        )

        # Pre-fetch accession lengths
        acc_lengths = {}
        for ssr in ssrs:
            if ssr.accession not in acc_lengths:
                acc_obj = session.query(Accession).filter(
                    Accession.accession == ssr.accession
                ).first()
                if acc_obj and hasattr(acc_obj, "sequence_length"):
                    acc_lengths[ssr.accession] = acc_obj.sequence_length
                else:
                    acc_lengths[ssr.accession] = None

        # Validate each SSR
        for ssr in ssrs:
            error = self.validate_single(
                ssr_id=ssr.ssr_id,
                accession=ssr.accession,
                start=ssr.start,
                end=ssr.end,
                accession_length=acc_lengths.get(ssr.accession),
            )

            if error is None:
                report.valid_count += 1
            else:
                if error.severity == "error":
                    report.invalid_count += 1
                    report.errors.append(error)
                else:
                    report.warnings.append(error)

        return report


# ---------------------------------------------------------------------------
# Browser Track Generation
# ---------------------------------------------------------------------------

class BrowserTrackGenerator:
    """Generates genome browser tracks in various formats."""

    def __init__(self, config: BedTrackConfig = None):
        """Initialize generator.

        Args:
            config: BedTrackConfig with generation options.
        """
        self.config = config or BedTrackConfig()
        self.validator = CoordinateValidator(strict=self.config.strict_validation)

    def _scale_score(self, repeat_length_bp: int, max_length: int = 500) -> int:
        """Scale repeat length to 0-1000 for BED score field.

        Args:
            repeat_length_bp: Repeat length in base pairs
            max_length: Maximum repeat length to consider (default 500bp = score 1000)

        Returns:
            Scaled score (0-1000)
        """
        if max_length <= 0:
            return 0
        scaled = min(1000, int((repeat_length_bp / max_length) * 1000))
        return max(0, scaled)

    def _get_color(self, repeat_class: str) -> str:
        """Get RGB color for repeat class.

        Args:
            repeat_class: repeat_class value (perfect/imperfect/compound_component)

        Returns:
            RGB hex string (R,G,B)
        """
        return self.config.color_scheme.get(
            repeat_class,
            self.config.color_scheme.get(RepeatClass.PERFECT, "0,0,0"),
        )

    def build_bed12_record(self, ssr: SSRRecord) -> Bed12Record:
        """Build BED12 record from SSRRecord.

        Args:
            ssr: SSRRecord object from database

        Returns:
            Bed12Record formatted for output
        """
        # Determine score
        if self.config.score_field == "repeat_length_bp":
            score = self._scale_score(ssr.repeat_length_bp)
        else:
            score = 500  # Default middle value

        # Determine strand
        strand = ssr.strand or "+"
        if strand not in ("+", "-", "."):
            strand = "."

        # Get color
        repeat_class = getattr(ssr, "repeat_class", "perfect")
        color = self._get_color(repeat_class)

        # Build name
        name = f"{ssr.motif_canonical}x{ssr.repeat_units}"

        # Thick region (full SSR if using thick region)
        if self.config.use_thick_region:
            thick_start = ssr.start
            thick_end = ssr.end
        else:
            thick_start = ssr.start
            thick_end = ssr.start + 1  # Minimal thick region

        return Bed12Record(
            chrom=ssr.accession,
            chrom_start=ssr.start,
            chrom_end=ssr.end,
            name=name,
            score=score,
            strand=strand,
            thick_start=thick_start,
            thick_end=thick_end,
            item_rgb=color,
        )

    def generate_track_header(self) -> str:
        """Generate UCSC Genome Browser track header line.

        Returns:
            Track header string
        """
        header = (
            f'track name="{self.config.track_name}" '
            f'description="{self.config.track_description}" '
            f'itemRgb=On\n'
        )
        return header


# ---------------------------------------------------------------------------
# Main Export Function
# ---------------------------------------------------------------------------

def generate_browser_track(
    session: Session,
    output_dir: str | Path,
    *,
    track_name: str = "GWICO-SSR",
    run_id: int | None = None,
    dataset_id: int | None = None,
    validate_coordinates: bool = True,
    strict_validation: bool = False,
) -> tuple[Path, CoordinateValidationReport]:
    """Generate browser-ready BED12 track from SSRRecords.

    Args:
        session: SQLAlchemy session
        output_dir: Output directory for track file
        track_name: Name for track
        run_id: Optional run_id filter
        dataset_id: Optional dataset_id filter
        validate_coordinates: Whether to validate coordinates
        strict_validation: If True, skip invalid records; if False, warn only

    Returns:
        Tuple of (output_file_path, validation_report)
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    output_file = output_dir / f"{track_name.lower().replace(' ', '_')}.bed12"

    # Build query
    stmt = select(SSRRecord).order_by(SSRRecord.accession, SSRRecord.start)

    if run_id is not None:
        stmt = stmt.where(SSRRecord.run_id == run_id)

    if dataset_id is not None:
        stmt = stmt.join(Accession).where(Accession.dataset_id == dataset_id)

    ssrs = session.execute(stmt).scalars().all()

    # Validate if requested
    validation_report = None
    if validate_coordinates:
        validator = CoordinateValidator(strict=strict_validation)
        validation_report = validator.validate_batch(session, ssrs)
        logger.info(
            f"Coordinate validation: {validation_report.valid_count}/{validation_report.total_records} "
            f"valid ({validation_report.pass_rate:.1f}%)"
        )

    # Generate track
    config = BedTrackConfig(
        track_name=track_name,
        strict_validation=strict_validation,
    )
    generator = BrowserTrackGenerator(config)

    with open(output_file, "w", encoding="utf-8") as f:
        # Write header
        if config.include_header:
            f.write(generator.generate_track_header())

        # Write records
        for ssr in ssrs:
            # Skip invalid if strict
            if strict_validation and validation_report:
                if any(e.ssr_id == ssr.ssr_id for e in validation_report.errors):
                    continue

            bed_record = generator.build_bed12_record(ssr)
            f.write(bed_record.to_bed12() + "\n")

    logger.info(f"Generated browser track: {output_file}")
    return output_file, validation_report
