"""FASTA file parser for GWICO-SSR.

Parses FASTA files to extract sequence data, compute sequence hash,
length, and GC content. Results are returned as structured objects
suitable for database persistence.
"""

from __future__ import annotations

import hashlib
import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional, Sequence

from Bio import SeqIO

logger = logging.getLogger(__name__)


@dataclass
class ParsedSequence:
    """Parsed data from a single FASTA record."""

    accession: str
    sequence: str
    sequence_length: int
    sequence_hash: str
    gc_content: float
    description: str = ""


@dataclass
class ParseError:
    """Structured error from parsing."""

    file_path: str
    accession: Optional[str]
    message: str


@dataclass
class FastaParseResult:
    """Result of parsing a FASTA file (may contain multiple records)."""

    file_path: str
    records: list[ParsedSequence] = field(default_factory=list)
    errors: list[ParseError] = field(default_factory=list)

    @property
    def success(self) -> bool:
        return len(self.records) > 0 and len(self.errors) == 0

    @property
    def record_count(self) -> int:
        return len(self.records)

    def to_dict(self) -> dict:
        return {
            "file_path": self.file_path,
            "record_count": self.record_count,
            "success": self.success,
            "errors": [
                {"file_path": e.file_path, "accession": e.accession, "message": e.message}
                for e in self.errors
            ],
        }


def _compute_gc_content(sequence: str) -> float:
    """Compute GC content as a fraction (0.0 to 1.0)."""
    if not sequence:
        return 0.0
    upper = sequence.upper()
    gc = upper.count("G") + upper.count("C")
    # Only count definite bases (A, T, G, C)
    total = upper.count("A") + upper.count("T") + upper.count("G") + upper.count("C")
    if total == 0:
        return 0.0
    return gc / total


def _compute_hash(sequence: str) -> str:
    """Compute SHA-256 hash of the uppercase sequence."""
    return hashlib.sha256(sequence.upper().encode("ascii", errors="ignore")).hexdigest()


def _extract_accession(record_id: str) -> str:
    """Extract accession from a FASTA record ID.

    Handles formats like 'NC_045512.2' or 'NC_045512.2 some description'.
    """
    # BioPython record.id is the first whitespace-delimited token
    return record_id.split()[0].strip()


def parse_fasta(file_path: str | Path) -> FastaParseResult:
    """Parse a FASTA file and extract sequence records.

    Args:
        file_path: Path to the FASTA file.

    Returns:
        FastaParseResult with parsed records and any errors.
    """
    file_path = Path(file_path)
    result = FastaParseResult(file_path=str(file_path))

    if not file_path.exists():
        result.errors.append(ParseError(
            file_path=str(file_path),
            accession=None,
            message=f"File not found: {file_path}",
        ))
        return result

    if file_path.stat().st_size == 0:
        result.errors.append(ParseError(
            file_path=str(file_path),
            accession=None,
            message="File is empty",
        ))
        return result

    try:
        for record in SeqIO.parse(str(file_path), "fasta"):
            try:
                accession = _extract_accession(record.id)
                seq_str = str(record.seq)

                if not seq_str:
                    result.errors.append(ParseError(
                        file_path=str(file_path),
                        accession=accession,
                        message="Empty sequence",
                    ))
                    continue

                parsed = ParsedSequence(
                    accession=accession,
                    sequence=seq_str,
                    sequence_length=len(seq_str),
                    sequence_hash=_compute_hash(seq_str),
                    gc_content=_compute_gc_content(seq_str),
                    description=record.description,
                )
                result.records.append(parsed)
                logger.debug("Parsed FASTA record: %s (%d bp)", accession, len(seq_str))

            except Exception as e:
                result.errors.append(ParseError(
                    file_path=str(file_path),
                    accession=getattr(record, "id", None),
                    message=f"Error parsing record: {e}",
                ))
    except Exception as e:
        result.errors.append(ParseError(
            file_path=str(file_path),
            accession=None,
            message=f"Error reading FASTA file: {e}",
        ))

    return result


def parse_fasta_multi(file_paths: Sequence[str | Path]) -> list[FastaParseResult]:
    """Parse multiple FASTA files.

    Args:
        file_paths: Paths to the FASTA files.

    Returns:
        List of FastaParseResult, one per file.
    """
    return [parse_fasta(fp) for fp in file_paths]
