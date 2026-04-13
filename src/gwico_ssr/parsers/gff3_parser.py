"""GFF3 annotation file parser for GWICO-SSR.

Parses GFF3 files to extract genomic features (CDS, gene, mRNA, etc.)
into the same normalized ParsedFeature structure used by the GenBank parser.
Uses a lightweight custom parser — no heavy GFF library dependency.

GFF3 spec reference: https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md
Columns: seqid, source, type, start, end, score, strand, phase, attributes
"""

from __future__ import annotations

import logging
import urllib.parse
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional, Sequence

logger = logging.getLogger(__name__)

# Reuse ParsedFeature from genbank_parser for normalized output
from gwico_ssr.parsers.genbank_parser import ParsedFeature

# GFF3 feature types we extract (same as GenBank)
GFF3_FEATURE_TYPES = {"CDS", "gene", "mRNA", "tRNA", "rRNA", "ncRNA", "exon", "misc_feature"}


@dataclass
class GFF3ParseResult:
    """Result of parsing a GFF3 file."""

    file_path: str
    features: list[ParsedFeature] = field(default_factory=list)
    errors: list[dict] = field(default_factory=list)
    directives: list[str] = field(default_factory=list)  # ##directives found

    @property
    def success(self) -> bool:
        return len(self.features) > 0 and len(self.errors) == 0

    @property
    def feature_count(self) -> int:
        return len(self.features)

    @property
    def cds_count(self) -> int:
        return sum(1 for f in self.features if f.feature_type == "CDS")

    def to_dict(self) -> dict:
        return {
            "file_path": self.file_path,
            "feature_count": self.feature_count,
            "cds_count": self.cds_count,
            "success": self.success,
            "errors": self.errors,
        }


def _parse_attributes(attr_str: str) -> dict[str, str]:
    """Parse GFF3 column 9 attributes.

    Format: key1=value1;key2=value2
    Values are URL-encoded per the GFF3 spec.
    """
    attrs: dict[str, str] = {}
    if not attr_str or attr_str == ".":
        return attrs
    for pair in attr_str.split(";"):
        pair = pair.strip()
        if not pair:
            continue
        if "=" in pair:
            key, _, value = pair.partition("=")
            attrs[key.strip()] = urllib.parse.unquote(value.strip())
        else:
            # Bare attribute without value
            attrs[pair] = ""
    return attrs


def _normalize_strand(strand_char: str) -> str:
    """Normalize GFF3 strand character."""
    if strand_char == "+":
        return "+"
    elif strand_char == "-":
        return "-"
    return "."


def parse_gff3(
    file_path: str | Path,
    default_accession: Optional[str] = None,
) -> GFF3ParseResult:
    """Parse a GFF3 file and extract features.

    GFF3 coordinates are 1-based inclusive; we convert to 0-based half-open
    (matching BioPython/GenBank convention) for consistency.

    Args:
        file_path: Path to the GFF3 file.
        default_accession: Accession to use if seqid column doesn't match
            an expected accession format. If None, uses seqid as-is.

    Returns:
        GFF3ParseResult with parsed features and errors.
    """
    file_path = Path(file_path)
    result = GFF3ParseResult(file_path=str(file_path))

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
        with open(file_path, "r", encoding="utf-8", errors="replace") as fh:
            for line_num, line in enumerate(fh, start=1):
                line = line.rstrip("\n\r")

                # Skip empty lines
                if not line:
                    continue

                # Directives (## lines)
                if line.startswith("##"):
                    result.directives.append(line)
                    if line == "##FASTA":
                        # FASTA section follows; stop parsing features
                        break
                    continue

                # Comment lines
                if line.startswith("#"):
                    continue

                # Parse feature line
                parts = line.split("\t")
                if len(parts) != 9:
                    result.errors.append({
                        "file_path": str(file_path),
                        "accession": None,
                        "message": f"Line {line_num}: expected 9 tab-separated columns, got {len(parts)}",
                    })
                    continue

                seqid, source, feat_type, start_str, end_str, score, strand_char, phase, attr_str = parts

                # Filter to feature types of interest
                if feat_type not in GFF3_FEATURE_TYPES:
                    continue

                try:
                    # GFF3 uses 1-based inclusive coordinates
                    # Convert to 0-based half-open for internal consistency
                    start = int(start_str) - 1
                    end = int(end_str)
                except ValueError:
                    result.errors.append({
                        "file_path": str(file_path),
                        "accession": seqid,
                        "message": f"Line {line_num}: invalid coordinates start={start_str}, end={end_str}",
                    })
                    continue

                attrs = _parse_attributes(attr_str)
                accession = default_accession if default_accession else seqid

                gene_name = attrs.get("gene") or attrs.get("Name") or attrs.get("gene_name")
                product = attrs.get("product")
                locus_tag = attrs.get("locus_tag")

                feature = ParsedFeature(
                    accession=accession,
                    feature_type=feat_type,
                    start=start,
                    end=end,
                    strand=_normalize_strand(strand_char),
                    gene_name=gene_name,
                    product=product,
                    locus_tag=locus_tag,
                    annotation_source="gff3",
                )
                result.features.append(feature)

    except Exception as e:
        result.errors.append({
            "file_path": str(file_path),
            "accession": None,
            "message": f"Error reading GFF3 file: {e}",
        })

    if result.features:
        logger.debug(
            "Parsed GFF3: %s (%d features, %d CDS)",
            file_path.name,
            result.feature_count,
            result.cds_count,
        )

    return result


def parse_gff3_multi(
    file_paths: Sequence[str | Path],
    default_accession: Optional[str] = None,
) -> list[GFF3ParseResult]:
    """Parse multiple GFF3 files.

    Args:
        file_paths: Paths to GFF3 files.
        default_accession: Default accession if seqid doesn't match.

    Returns:
        List of GFF3ParseResult, one per file.
    """
    return [parse_gff3(fp, default_accession=default_accession) for fp in file_paths]
