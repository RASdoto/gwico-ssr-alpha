"""GenBank file parser for GWICO-SSR.

Parses GenBank files to extract CDS and other genomic features,
normalizing gene names, coordinates, and strand information.
Also extracts sequence data when present.
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
class ParsedFeature:
    """A normalized genomic feature extracted from GenBank."""

    accession: str
    feature_type: str
    start: int
    end: int
    strand: str  # "+", "-", or "."
    gene_name: Optional[str] = None
    product: Optional[str] = None
    locus_tag: Optional[str] = None
    annotation_source: str = "genbank"


@dataclass
class GenBankSequenceInfo:
    """Sequence-level data extracted from GenBank."""

    accession: str
    sequence_length: int
    sequence_hash: Optional[str] = None
    gc_content: Optional[float] = None
    is_circular: bool = False
    organism: Optional[str] = None
    description: str = ""


@dataclass
class GenBankParseResult:
    """Result of parsing a GenBank file."""

    file_path: str
    sequence_info: Optional[GenBankSequenceInfo] = None
    features: list[ParsedFeature] = field(default_factory=list)
    errors: list[dict] = field(default_factory=list)

    @property
    def success(self) -> bool:
        return self.sequence_info is not None and len(self.errors) == 0

    @property
    def feature_count(self) -> int:
        return len(self.features)

    @property
    def cds_count(self) -> int:
        return sum(1 for f in self.features if f.feature_type == "CDS")

    def to_dict(self) -> dict:
        return {
            "file_path": self.file_path,
            "accession": self.sequence_info.accession if self.sequence_info else None,
            "sequence_length": self.sequence_info.sequence_length if self.sequence_info else None,
            "feature_count": self.feature_count,
            "cds_count": self.cds_count,
            "success": self.success,
            "errors": self.errors,
        }


# Feature types we extract from GenBank files
FEATURE_TYPES_OF_INTEREST = {"CDS", "gene", "mRNA", "tRNA", "rRNA", "ncRNA", "misc_feature"}


def _normalize_strand(strand_int: Optional[int]) -> str:
    """Convert BioPython strand integer to string representation."""
    if strand_int == 1:
        return "+"
    elif strand_int == -1:
        return "-"
    return "."


def _extract_qualifier(feature, *keys: str) -> Optional[str]:
    """Extract the first matching qualifier value from a BioPython feature."""
    for key in keys:
        vals = feature.qualifiers.get(key, [])
        if vals:
            return vals[0]
    return None


def _compute_gc_content(sequence: str) -> float:
    """Compute GC content as a fraction (0.0 to 1.0)."""
    if not sequence:
        return 0.0
    upper = sequence.upper()
    gc = upper.count("G") + upper.count("C")
    total = upper.count("A") + upper.count("T") + upper.count("G") + upper.count("C")
    if total == 0:
        return 0.0
    return gc / total


def parse_genbank(file_path: str | Path) -> GenBankParseResult:
    """Parse a GenBank file and extract sequence info and features.

    Extracts CDS, gene, mRNA, and other annotation features, normalizing
    coordinates and strand. For compound/join locations, uses the overall
    span (min start, max end).

    Args:
        file_path: Path to the GenBank file.

    Returns:
        GenBankParseResult with parsed sequence info, features, and errors.
    """
    file_path = Path(file_path)
    result = GenBankParseResult(file_path=str(file_path))

    if not file_path.exists():
        result.errors.append({
            "file_path": str(file_path),
            "accession": None,
            "message": f"File not found: {file_path}",
        })
        return result

    if file_path.stat().st_size == 0:
        result.errors.append({
            "file_path": str(file_path),
            "accession": None,
            "message": "File is empty",
        })
        return result

    try:
        # GenBank files for single accessions typically have one record
        records = list(SeqIO.parse(str(file_path), "genbank"))
        if not records:
            result.errors.append({
                "file_path": str(file_path),
                "accession": None,
                "message": "No records found in GenBank file",
            })
            return result

        record = records[0]
        accession = record.id if record.id and record.id != "<unknown id>" else record.name

        # Sequence info
        seq_str = str(record.seq) if record.seq else ""
        seq_hash = hashlib.sha256(seq_str.upper().encode("ascii", errors="ignore")).hexdigest() if seq_str else None
        gc = _compute_gc_content(seq_str) if seq_str else None

        # Check topology for circularity
        is_circular = False
        annotations = getattr(record, "annotations", {})
        topology = annotations.get("topology", "")
        if topology == "circular":
            is_circular = True

        organism = annotations.get("organism", None)

        result.sequence_info = GenBankSequenceInfo(
            accession=accession,
            sequence_length=len(seq_str) if seq_str else 0,
            sequence_hash=seq_hash,
            gc_content=gc,
            is_circular=is_circular,
            organism=organism,
            description=record.description or "",
        )

        # Extract features
        for feature in record.features:
            if feature.type not in FEATURE_TYPES_OF_INTEREST:
                continue

            try:
                # Use int() on location positions for BioPython FeatureLocation
                start = int(feature.location.start)
                end = int(feature.location.end)
                strand = _normalize_strand(feature.location.strand)

                gene_name = _extract_qualifier(feature, "gene", "gene_synonym")
                product = _extract_qualifier(feature, "product")
                locus_tag = _extract_qualifier(feature, "locus_tag")

                parsed_feat = ParsedFeature(
                    accession=accession,
                    feature_type=feature.type,
                    start=start,
                    end=end,
                    strand=strand,
                    gene_name=gene_name,
                    product=product,
                    locus_tag=locus_tag,
                    annotation_source="genbank",
                )
                result.features.append(parsed_feat)

            except Exception as e:
                result.errors.append({
                    "file_path": str(file_path),
                    "accession": accession,
                    "message": f"Error parsing feature {feature.type}: {e}",
                })

        logger.debug(
            "Parsed GenBank: %s (%d bp, %d features, %d CDS)",
            accession,
            result.sequence_info.sequence_length,
            result.feature_count,
            result.cds_count,
        )

    except Exception as e:
        result.errors.append({
            "file_path": str(file_path),
            "accession": None,
            "message": f"Error reading GenBank file: {e}",
        })

    return result


def parse_genbank_multi(file_paths: Sequence[str | Path]) -> list[GenBankParseResult]:
    """Parse multiple GenBank files.

    Args:
        file_paths: Paths to GenBank files.

    Returns:
        List of GenBankParseResult, one per file.
    """
    return [parse_genbank(fp) for fp in file_paths]
