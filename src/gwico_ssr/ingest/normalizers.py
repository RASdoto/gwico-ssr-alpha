"""Metadata normalization routines for GWICO-SSR ingestion.

Handles date parsing, completeness flag extraction, country/region splitting,
and other field-level normalization for accession metadata.
Also provides composite FASTA and GenBank file splitting utilities.
"""

from __future__ import annotations

import hashlib
import logging
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional, Sequence

from Bio import SeqIO

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Date normalization
# ---------------------------------------------------------------------------

_ISO_FULL_RE = re.compile(
    r"^(\d{4})-(\d{2})-(\d{2})"  # YYYY-MM-DD
    r"(?:T\d{2}:\d{2}:\d{2}Z?)?$"  # optional T00:00:00Z
)
_YEAR_MONTH_RE = re.compile(r"^(\d{4})-(\d{2})$")
_YEAR_ONLY_RE = re.compile(r"^(\d{4})$")


def normalize_date(raw: Optional[str]) -> Optional[str]:
    """Normalize a date string to a consistent format.

    Accepts: ``2020-01-13T00:00:00Z``, ``2020-01-13``, ``2020-01``, ``2020``.
    Returns:  ``YYYY-MM-DD``, ``YYYY-MM``, or ``YYYY`` depending on precision.
    Returns ``None`` for empty or unparseable values.
    """
    if not raw or not raw.strip():
        return None
    raw = raw.strip()

    m = _ISO_FULL_RE.match(raw)
    if m:
        return f"{m.group(1)}-{m.group(2)}-{m.group(3)}"

    m = _YEAR_MONTH_RE.match(raw)
    if m:
        return f"{m.group(1)}-{m.group(2)}"

    m = _YEAR_ONLY_RE.match(raw)
    if m:
        return m.group(1)

    return None


# ---------------------------------------------------------------------------
# Completeness normalization
# ---------------------------------------------------------------------------

_COMPLETE_VALUES = frozenset({"complete", "complete genome"})


def normalize_completeness(raw: Optional[str]) -> Optional[bool]:
    """Convert nucleotide completeness string to a boolean.

    ``"complete"`` or ``"complete genome"`` → ``True``.
    Any other non-empty value → ``False``.
    Empty / missing → ``None``.
    """
    if not raw or not raw.strip():
        return None
    val = raw.strip().lower()
    return val in _COMPLETE_VALUES


# ---------------------------------------------------------------------------
# Geo-location normalization
# ---------------------------------------------------------------------------

def normalize_geo_location(raw: Optional[str]) -> tuple[Optional[str], Optional[str]]:
    """Parse geo_location into (country, region).

    Expected formats:
    - ``"China"`` → ``("China", None)``
    - ``"USA: California"`` → ``("USA", "California")``
    - ``"South Korea: Seoul"`` → ``("South Korea", "Seoul")``

    Returns ``(None, None)`` for empty / missing values.
    """
    if not raw or not raw.strip():
        return None, None
    raw = raw.strip()

    if ":" in raw:
        parts = raw.split(":", 1)
        country = parts[0].strip() or None
        region = parts[1].strip() or None
        return country, region

    return raw, None


# ---------------------------------------------------------------------------
# Genome length normalization
# ---------------------------------------------------------------------------

def normalize_length(raw: Optional[str]) -> Optional[int]:
    """Convert length field to an integer.

    Returns ``None`` for empty, non-numeric, or negative values.
    """
    if not raw or not str(raw).strip():
        return None
    try:
        val = int(str(raw).strip())
        return val if val > 0 else None
    except (ValueError, TypeError):
        return None


# ---------------------------------------------------------------------------
# Accession normalization
# ---------------------------------------------------------------------------

_ACCESSION_RE = re.compile(r"^[A-Za-z0-9_.]+$")


def normalize_accession(raw: Optional[str]) -> Optional[str]:
    """Validate and normalize an accession identifier.

    Strips whitespace and verifies format.
    Returns ``None`` for empty or invalid accession strings.
    """
    if not raw or not raw.strip():
        return None
    val = raw.strip()
    if not _ACCESSION_RE.match(val):
        return None
    return val


# ---------------------------------------------------------------------------
# Composite File Normalization
# ---------------------------------------------------------------------------

@dataclass
class CompositeNormalizationResult:
    """Result of normalizing a composite FASTA or GenBank file."""

    source_path: str
    file_type: str  # "fasta" or "genbank"
    checksum_sha256: str
    record_count: int
    accession_count: int
    normalized_records: list[dict] = field(default_factory=list)  # List of {accession, output_path}
    errors: list[dict] = field(default_factory=list)
    duplicates: list[dict] = field(default_factory=list)

    def to_dict(self) -> dict:
        return {
            "source_path": self.source_path,
            "file_type": self.file_type,
            "checksum_sha256": self.checksum_sha256,
            "record_count": self.record_count,
            "accession_count": self.accession_count,
            "normalized_records": self.normalized_records,
            "errors": self.errors,
            "duplicates": self.duplicates,
        }


def _extract_accession_from_fasta_id(record_id: str) -> str:
    """Extract accession from FASTA record ID (first whitespace-delimited token)."""
    return record_id.split()[0].strip()


def _compute_sha256(content: str) -> str:
    """Compute SHA-256 hash of content."""
    return hashlib.sha256(content.encode("utf-8", errors="ignore")).hexdigest()


def split_composite_fasta(
    composite_path: str | Path,
    output_dir: str | Path,
    overwrite: bool = False,
    duplicate_policy: str = "first",
) -> CompositeNormalizationResult:
    """Split a composite (multi-record) FASTA file into individual normalized files.

    Args:
        composite_path: Path to the composite FASTA file.
        output_dir: Directory where normalized single-record files will be written.
        overwrite: If True, overwrite existing normalized files; if False, skip duplicates.
        duplicate_policy: "first" (keep first), "error" (fail on duplicates), or "skip" (ignore extras).

    Returns:
        CompositeNormalizationResult with per-record output paths and any issues.
    """
    composite_path = Path(composite_path)
    output_dir = Path(output_dir)
    result = CompositeNormalizationResult(
        source_path=str(composite_path.resolve()),
        file_type="fasta",
        checksum_sha256="",
        record_count=0,
        accession_count=0,
    )

    if not composite_path.exists():
        result.errors.append({
            "accession": None,
            "message": f"File not found: {composite_path}",
        })
        return result

    if composite_path.stat().st_size == 0:
        result.errors.append({
            "accession": None,
            "message": "File is empty",
        })
        return result

    # Read entire file for checksum
    try:
        content = composite_path.read_text(encoding="utf-8", errors="ignore")
        result.checksum_sha256 = _compute_sha256(content)
    except Exception as e:
        result.errors.append({
            "accession": None,
            "message": f"Error reading file for checksum: {e}",
        })
        return result

    # Parse and split records
    seen_accessions = set()
    try:
        output_dir.mkdir(parents=True, exist_ok=True)
        for record in SeqIO.parse(str(composite_path), "fasta"):
            result.record_count += 1
            accession = _extract_accession_from_fasta_id(record.id)

            if accession in seen_accessions:
                result.duplicates.append({
                    "accession": accession,
                    "record_index": result.record_count - 1,
                    "policy": duplicate_policy,
                })
                if duplicate_policy == "error":
                    result.errors.append({
                        "accession": accession,
                        "message": f"Duplicate accession in composite file (duplicate_policy={duplicate_policy})",
                    })
                    continue
                elif duplicate_policy == "skip":
                    continue
                # "first" policy: silently skip additional occurrences

            seen_accessions.add(accession)

            # Write normalized file
            try:
                output_path = output_dir / f"{accession}.fasta"
                if output_path.exists() and not overwrite:
                    result.errors.append({
                        "accession": accession,
                        "message": f"Output file already exists (set overwrite=True to replace): {output_path}",
                    })
                    continue

                # Write single record to file
                SeqIO.write([record], str(output_path), "fasta")
                result.normalized_records.append({
                    "accession": accession,
                    "output_path": str(output_path.resolve()),
                })
                logger.debug("Split composite FASTA: %s -> %s", accession, output_path)

            except Exception as e:
                result.errors.append({
                    "accession": accession,
                    "message": f"Error writing normalized FASTA: {e}",
                })

    except Exception as e:
        result.errors.append({
            "accession": None,
            "message": f"Error parsing composite FASTA: {e}",
        })

    result.accession_count = len(seen_accessions)
    return result


def split_composite_genbank(
    composite_path: str | Path,
    output_dir: str | Path,
    overwrite: bool = False,
    duplicate_policy: str = "first",
) -> CompositeNormalizationResult:
    """Split a composite (multi-record) GenBank file into individual normalized files.

    Args:
        composite_path: Path to the composite GenBank file.
        output_dir: Directory where normalized single-record files will be written.
        overwrite: If True, overwrite existing normalized files; if False, skip duplicates.
        duplicate_policy: "first" (keep first), "error" (fail on duplicates), or "skip" (ignore extras).

    Returns:
        CompositeNormalizationResult with per-record output paths and any issues.
    """
    composite_path = Path(composite_path)
    output_dir = Path(output_dir)
    result = CompositeNormalizationResult(
        source_path=str(composite_path.resolve()),
        file_type="genbank",
        checksum_sha256="",
        record_count=0,
        accession_count=0,
    )

    if not composite_path.exists():
        result.errors.append({
            "accession": None,
            "message": f"File not found: {composite_path}",
        })
        return result

    if composite_path.stat().st_size == 0:
        result.errors.append({
            "accession": None,
            "message": "File is empty",
        })
        return result

    # Read entire file for checksum
    try:
        content = composite_path.read_text(encoding="utf-8", errors="ignore")
        result.checksum_sha256 = _compute_sha256(content)
    except Exception as e:
        result.errors.append({
            "accession": None,
            "message": f"Error reading file for checksum: {e}",
        })
        return result

    # Parse and split records
    seen_accessions = set()
    try:
        output_dir.mkdir(parents=True, exist_ok=True)
        for record in SeqIO.parse(str(composite_path), "genbank"):
            result.record_count += 1
            accession = record.id if record.id and record.id != "<unknown id>" else record.name

            if not accession:
                result.errors.append({
                    "accession": None,
                    "message": f"Record {result.record_count}: Cannot extract accession from GenBank record",
                })
                continue

            if accession in seen_accessions:
                result.duplicates.append({
                    "accession": accession,
                    "record_index": result.record_count - 1,
                    "policy": duplicate_policy,
                })
                if duplicate_policy == "error":
                    result.errors.append({
                        "accession": accession,
                        "message": f"Duplicate accession in composite file (duplicate_policy={duplicate_policy})",
                    })
                    continue
                elif duplicate_policy == "skip":
                    continue
                # "first" policy: silently skip additional occurrences

            seen_accessions.add(accession)

            # Write normalized file
            try:
                output_path = output_dir / f"{accession}.gb"
                if output_path.exists() and not overwrite:
                    result.errors.append({
                        "accession": accession,
                        "message": f"Output file already exists (set overwrite=True to replace): {output_path}",
                    })
                    continue

                # Write single record to file
                SeqIO.write([record], str(output_path), "genbank")
                result.normalized_records.append({
                    "accession": accession,
                    "output_path": str(output_path.resolve()),
                })
                logger.debug("Split composite GenBank: %s -> %s", accession, output_path)

            except Exception as e:
                result.errors.append({
                    "accession": accession,
                    "message": f"Error writing normalized GenBank: {e}",
                })

    except Exception as e:
        result.errors.append({
            "accession": None,
            "message": f"Error parsing composite GenBank: {e}",
        })

    result.accession_count = len(seen_accessions)
    return result

