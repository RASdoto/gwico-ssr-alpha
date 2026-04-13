"""CSV ingestion and validation for GWICO-SSR accession metadata.

Parses accession metadata CSVs, validates fields, normalizes values,
detects duplicates, and persists records to the database.
"""

from __future__ import annotations

import csv
import hashlib
import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional, Sequence

from sqlalchemy.orm import Session

from gwico_ssr.db.repository import (
    bulk_upsert_accessions,
    get_or_create_dataset,
    upsert_sequence_record,
)
from gwico_ssr.ingest.normalizers import (
    normalize_accession,
    normalize_completeness,
    normalize_date,
    normalize_geo_location,
    normalize_length,
)

logger = logging.getLogger(__name__)

# Expected CSV header columns (case-insensitive matching)
EXPECTED_COLUMNS = {
    "accession",
    "release_date",
    "species",
    "length",
    "nuc_completeness",
    "geo_location",
    "usa",
    "host",
    "isolation_source",
    "collection_date",
}

# Minimum required columns for a valid CSV
REQUIRED_COLUMNS = {"accession"}


@dataclass
class RowError:
    """A validation error for a single CSV row."""

    row_number: int
    field: str
    message: str
    raw_value: Optional[str] = None


@dataclass
class IngestSummary:
    """Summary of an ingestion operation."""

    source_file: str
    total_rows: int = 0
    valid_rows: int = 0
    invalid_rows: int = 0
    duplicate_rows: int = 0
    upserted_rows: int = 0
    skipped_rows: int = 0
    errors: list[RowError] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)

    def to_dict(self) -> dict:
        """Convert summary to a dictionary for JSON serialization."""
        return {
            "source_file": self.source_file,
            "total_rows": self.total_rows,
            "valid_rows": self.valid_rows,
            "invalid_rows": self.invalid_rows,
            "duplicate_rows": self.duplicate_rows,
            "upserted_rows": self.upserted_rows,
            "skipped_rows": self.skipped_rows,
            "error_count": len(self.errors),
            "errors": [
                {
                    "row": e.row_number,
                    "field": e.field,
                    "message": e.message,
                    "raw_value": e.raw_value,
                }
                for e in self.errors
            ],
            "warnings": self.warnings,
        }


def _normalize_header(raw_header: str) -> str:
    """Normalize a CSV header name to lowercase with underscores."""
    return raw_header.strip().lower().replace(" ", "_")


def validate_header(header: list[str]) -> tuple[dict[str, int], list[str]]:
    """Validate and map CSV header columns.

    Returns a mapping of normalized column name → index, and a list of warnings.
    Raises ``ValueError`` if required columns are missing.
    """
    normalized = [_normalize_header(h) for h in header]
    col_map = {name: idx for idx, name in enumerate(normalized)}

    missing_required = REQUIRED_COLUMNS - set(normalized)
    if missing_required:
        raise ValueError(
            f"Missing required columns: {', '.join(sorted(missing_required))}"
        )

    warnings = []
    extra = set(normalized) - EXPECTED_COLUMNS
    if extra:
        warnings.append(f"Unexpected columns (ignored): {', '.join(sorted(extra))}")

    missing_optional = EXPECTED_COLUMNS - set(normalized)
    if missing_optional:
        warnings.append(
            f"Optional columns not present: {', '.join(sorted(missing_optional))}"
        )

    return col_map, warnings


def _get_field(row: list[str], col_map: dict[str, int], field_name: str) -> Optional[str]:
    """Safely get a field value from a row by column name."""
    idx = col_map.get(field_name)
    if idx is None or idx >= len(row):
        return None
    val = row[idx].strip()
    return val if val else None


def parse_csv_row(
    row: list[str],
    col_map: dict[str, int],
    row_number: int,
) -> tuple[Optional[dict], list[RowError]]:
    """Parse and normalize a single CSV row into an accession dict.

    Returns ``(record_dict, errors)``. If ``record_dict`` is ``None``,
    the row is invalid due to critical errors.
    """
    errors: list[RowError] = []

    # --- Accession (required) ---
    raw_accession = _get_field(row, col_map, "accession")
    accession = normalize_accession(raw_accession)
    if accession is None:
        errors.append(
            RowError(row_number, "accession", "Missing or invalid accession", raw_accession)
        )
        return None, errors

    # --- Release date ---
    raw_release = _get_field(row, col_map, "release_date")
    release_date = normalize_date(raw_release)
    if raw_release and release_date is None:
        errors.append(
            RowError(row_number, "release_date", "Unparseable date", raw_release)
        )

    # --- Species ---
    species = _get_field(row, col_map, "species")

    # --- Length ---
    raw_length = _get_field(row, col_map, "length")
    genome_length = normalize_length(raw_length)
    if raw_length and genome_length is None:
        errors.append(
            RowError(row_number, "length", "Invalid length value", raw_length)
        )

    # --- Completeness ---
    raw_complete = _get_field(row, col_map, "nuc_completeness")
    is_complete = normalize_completeness(raw_complete)

    # --- Geo location ---
    raw_geo = _get_field(row, col_map, "geo_location")
    country, region = normalize_geo_location(raw_geo)

    # --- USA column (supplementary geo for US states) ---
    raw_usa = _get_field(row, col_map, "usa")
    if raw_usa and country == "USA" and region is None:
        region = raw_usa

    # --- Host ---
    host = _get_field(row, col_map, "host")

    # --- Collection date ---
    raw_collection = _get_field(row, col_map, "collection_date")
    collection_date = normalize_date(raw_collection)
    if raw_collection and collection_date is None:
        errors.append(
            RowError(row_number, "collection_date", "Unparseable date", raw_collection)
        )

    record = {
        "accession": accession,
        "species": species,
        "release_date": release_date,
        "collection_date": collection_date,
        "geo_location_raw": raw_geo,
        "country": country,
        "region": region,
        "host": host,
        "is_complete": is_complete,
        "genome_length": genome_length,
        "source_priority": 1,  # CSV metadata source priority
    }

    return record, errors


def compute_file_hash(file_path: Path) -> str:
    """Compute SHA-256 hash of a file's contents."""
    h = hashlib.sha256()
    with open(file_path, "rb") as f:
        for chunk in iter(lambda: f.read(8192), b""):
            h.update(chunk)
    return h.hexdigest()


def load_csv(
    file_path: Path,
    session: Session,
    dataset_name: Optional[str] = None,
    source_priority: int = 1,
) -> IngestSummary:
    """Load an accession metadata CSV into the database.

    Args:
        file_path: Path to the CSV file.
        session: Active database session.
        dataset_name: Name for the dataset. Defaults to the filename stem.
        source_priority: Priority for dedup resolution (higher wins).

    Returns:
        An ``IngestSummary`` with counts and any errors.
    """
    file_path = Path(file_path)
    summary = IngestSummary(source_file=str(file_path))

    if not file_path.exists():
        summary.errors.append(
            RowError(0, "file", f"File not found: {file_path}", None)
        )
        return summary

    if not file_path.suffix.lower() == ".csv":
        summary.warnings.append(f"File does not have .csv extension: {file_path.name}")

    # Compute manifest hash
    file_hash = compute_file_hash(file_path)

    # Create or get dataset
    ds_name = dataset_name or file_path.stem
    dataset = get_or_create_dataset(
        session,
        name=ds_name,
        source_type="csv",
        input_manifest_hash=file_hash,
    )

    # Read and parse CSV
    seen_accessions: set[str] = set()
    records_to_upsert: list[dict] = []

    with open(file_path, "r", encoding="utf-8", errors="replace") as f:
        reader = csv.reader(f)

        # Read header
        try:
            raw_header = next(reader)
        except StopIteration:
            summary.errors.append(
                RowError(0, "file", "Empty CSV file", None)
            )
            return summary

        try:
            col_map, header_warnings = validate_header(raw_header)
        except ValueError as exc:
            summary.errors.append(
                RowError(0, "header", str(exc), None)
            )
            return summary
        summary.warnings.extend(header_warnings)

        # Parse rows
        for row_idx, row in enumerate(reader, start=2):  # row 1 = header
            summary.total_rows += 1

            if not any(cell.strip() for cell in row):
                summary.skipped_rows += 1
                continue

            record, errors = parse_csv_row(row, col_map, row_idx)

            if errors:
                summary.errors.extend(errors)

            if record is None:
                summary.invalid_rows += 1
                continue

            # Duplicate detection within this import
            acc_id = record["accession"]
            if acc_id in seen_accessions:
                summary.duplicate_rows += 1
                summary.warnings.append(
                    f"Row {row_idx}: duplicate accession '{acc_id}' within same file (skipped)"
                )
                continue
            seen_accessions.add(acc_id)

            # Apply dataset_id and source_priority
            record["dataset_id"] = dataset.dataset_id
            record["source_priority"] = source_priority

            records_to_upsert.append(record)
            summary.valid_rows += 1

    # Bulk upsert to database
    if records_to_upsert:
        summary.upserted_rows = bulk_upsert_accessions(session, records_to_upsert)

    logger.info(
        "CSV ingestion complete",
        extra={
            "file": str(file_path),
            "total": summary.total_rows,
            "valid": summary.valid_rows,
            "invalid": summary.invalid_rows,
            "duplicates": summary.duplicate_rows,
            "upserted": summary.upserted_rows,
        },
    )

    return summary


# ---------------------------------------------------------------------------
# Local file registration
# ---------------------------------------------------------------------------

_SEQUENCE_EXTENSIONS = {
    ".fasta": "fasta",
    ".fa": "fasta",
    ".fna": "fasta",
    ".gb": "genbank",
    ".gbk": "genbank",
    ".genbank": "genbank",
    ".gff": "gff3",
    ".gff3": "gff3",
}


def register_local_files(
    session: Session,
    accession: str,
    fasta_path: Optional[Path] = None,
    genbank_path: Optional[Path] = None,
    gff3_path: Optional[Path] = None,
) -> tuple[bool, list[str]]:
    """Register local sequence/annotation files for an accession.

    Validates that paths exist but does not parse file contents (deferred
    to the parser chunk). Creates or updates the SequenceRecord entry.

    Returns ``(success, errors)``.
    """
    errors: list[str] = []

    kwargs: dict = {
        "accession": accession,
        "sequence_source": "local",
        "download_status": "skipped",
    }

    if fasta_path is not None:
        fasta_path = Path(fasta_path)
        if not fasta_path.exists():
            errors.append(f"FASTA file not found: {fasta_path}")
        else:
            kwargs["fasta_path"] = str(fasta_path.resolve())

    if genbank_path is not None:
        genbank_path = Path(genbank_path)
        if not genbank_path.exists():
            errors.append(f"GenBank file not found: {genbank_path}")
        else:
            kwargs["genbank_path"] = str(genbank_path.resolve())

    if gff3_path is not None:
        # Store GFF3 path — the SequenceRecord model doesn't have a gff3_path
        # column, but we validate the path and log it for now.
        gff3_path = Path(gff3_path)
        if not gff3_path.exists():
            errors.append(f"GFF3 file not found: {gff3_path}")

    if errors:
        return False, errors

    upsert_sequence_record(session, **kwargs)

    logger.info(
        "Registered local files for accession",
        extra={
            "accession": accession,
            "fasta": str(fasta_path) if fasta_path else None,
            "genbank": str(genbank_path) if genbank_path else None,
            "gff3": str(gff3_path) if gff3_path else None,
        },
    )

    return True, []
