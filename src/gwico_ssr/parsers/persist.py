"""Persistence orchestrator for parsed sequence and feature data.

Coordinates parsing of FASTA/GenBank/GFF3 files and persisting
results to the database via the repository layer.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

from sqlalchemy.orm import Session

from gwico_ssr.db.repository import (
    get_accession,
    get_features_for_accession,
    insert_features,
    upsert_sequence_record,
)
from gwico_ssr.parsers.fasta_parser import ParsedSequence, parse_fasta
from gwico_ssr.parsers.genbank_parser import (
    GenBankParseResult,
    ParsedFeature,
    parse_genbank,
    parse_genbank_composite,
)
from gwico_ssr.parsers.gff3_parser import GFF3ParseResult, parse_gff3

logger = logging.getLogger(__name__)


@dataclass
class ParseSummary:
    """Summary of a parse batch operation."""

    total_accessions: int = 0
    fasta_parsed: int = 0
    fasta_failed: int = 0
    genbank_parsed: int = 0
    genbank_failed: int = 0
    gff3_parsed: int = 0
    gff3_failed: int = 0
    features_inserted: int = 0
    errors: list[dict] = field(default_factory=list)

    def to_dict(self) -> dict:
        return {
            "total_accessions": self.total_accessions,
            "fasta_parsed": self.fasta_parsed,
            "fasta_failed": self.fasta_failed,
            "genbank_parsed": self.genbank_parsed,
            "genbank_failed": self.genbank_failed,
            "gff3_parsed": self.gff3_parsed,
            "gff3_failed": self.gff3_failed,
            "features_inserted": self.features_inserted,
            "errors": self.errors[:50],  # cap for serialization
        }


def persist_parsed_fasta(
    session: Session,
    accession: str,
    parsed: ParsedSequence,
) -> None:
    """Persist parsed FASTA data to SequenceRecord and Accession."""
    upsert_sequence_record(
        session,
        accession=accession,
        sequence_hash=parsed.sequence_hash,
        sequence_length=parsed.sequence_length,
        parse_status="parsed",
    )

    # Update GC content and genome length on the Accession
    acc_obj = get_accession(session, accession)
    if acc_obj is not None:
        acc_obj.gc_content = parsed.gc_content
        acc_obj.genome_length = parsed.sequence_length
        session.flush()


def persist_parsed_features(
    session: Session,
    features: list[ParsedFeature],
    accession: str,
    replace_existing: bool = False,
) -> int:
    """Persist parsed features to the database.

    Args:
        session: DB session.
        features: List of parsed features to insert.
        accession: Accession these features belong to.
        replace_existing: If True, delete existing features first.

    Returns:
        Number of features inserted.
    """
    if replace_existing:
        existing = get_features_for_accession(session, accession)
        for feat in existing:
            session.delete(feat)
        session.flush()

    feature_dicts = [
        {
            "accession": f.accession,
            "feature_type": f.feature_type,
            "start": f.start,
            "end": f.end,
            "strand": f.strand,
            "gene_name": f.gene_name,
            "product": f.product,
            "locus_tag": f.locus_tag,
            "annotation_source": f.annotation_source,
        }
        for f in features
    ]
    return insert_features(session, feature_dicts)


def parse_and_persist(
    session: Session,
    accession: str,
    fasta_path: Optional[str] = None,
    genbank_path: Optional[str] = None,
    gff3_path: Optional[str] = None,
    force: bool = False,
) -> dict:
    """Parse available files for a single accession and persist results.

    Args:
        session: DB session.
        accession: Accession ID.
        fasta_path: Path to FASTA file.
        genbank_path: Path to GenBank file.
        gff3_path: Path to GFF3 file.
        force: Re-parse even if already parsed.

    Returns:
        Dict with parse status info.
    """
    result = {
        "accession": accession,
        "fasta_parsed": False,
        "genbank_parsed": False,
        "gff3_parsed": False,
        "features_inserted": 0,
        "errors": [],
    }

    # Check current parse status
    from gwico_ssr.models.schema import SequenceRecord
    seq_rec = session.get(SequenceRecord, accession)
    if seq_rec and seq_rec.parse_status == "parsed" and not force:
        result["skipped"] = True
        return result

    # Parse FASTA
    if fasta_path and Path(fasta_path).exists():
        fasta_result = parse_fasta(fasta_path)
        if fasta_result.records:
            # Find matching record
            matching = [r for r in fasta_result.records if r.accession == accession]
            if not matching:
                # If single record, use it regardless of ID
                matching = fasta_result.records[:1]
            if matching:
                persist_parsed_fasta(session, accession, matching[0])
                result["fasta_parsed"] = True
        for err in fasta_result.errors:
            result["errors"].append({"source": "fasta", "message": err.message})

    # Parse GenBank (features + optional sequence info)
    features_to_insert: list[ParsedFeature] = []
    if genbank_path and Path(genbank_path).exists():
        gb_result = parse_genbank(genbank_path)
        if gb_result.sequence_info:
            # Update sequence info if FASTA wasn't parsed
            if not result["fasta_parsed"]:
                upsert_sequence_record(
                    session,
                    accession=accession,
                    sequence_hash=gb_result.sequence_info.sequence_hash,
                    sequence_length=gb_result.sequence_info.sequence_length,
                    parse_status="parsed",
                )
                acc_obj = get_accession(session, accession)
                if acc_obj and gb_result.sequence_info.gc_content is not None:
                    acc_obj.gc_content = gb_result.sequence_info.gc_content
                    acc_obj.genome_length = gb_result.sequence_info.sequence_length
                    session.flush()

            # Update circularity info
            upsert_sequence_record(
                session,
                accession=accession,
                is_circular=gb_result.sequence_info.is_circular,
            )

            result["genbank_parsed"] = True
            features_to_insert.extend(gb_result.features)

        for err in gb_result.errors:
            result["errors"].append({"source": "genbank", "message": err.get("message", str(err))})

    # Parse GFF3 (features only)
    if gff3_path and Path(gff3_path).exists():
        gff3_result = parse_gff3(gff3_path, default_accession=accession)
        if gff3_result.features:
            result["gff3_parsed"] = True
            features_to_insert.extend(gff3_result.features)

        for err in gff3_result.errors:
            result["errors"].append({"source": "gff3", "message": err.get("message", str(err))})

    # Persist features (replace if force)
    if features_to_insert:
        count = persist_parsed_features(
            session, features_to_insert, accession, replace_existing=force
        )
        result["features_inserted"] = count

    # Update parse_status
    if result["fasta_parsed"] or result["genbank_parsed"]:
        upsert_sequence_record(session, accession=accession, parse_status="parsed")
    elif result["errors"]:
        upsert_sequence_record(session, accession=accession, parse_status="failed")

    return result


def parse_dataset_accessions(
    session: Session,
    accession_ids: list[str],
    data_dir: str | Path,
    gff3_dir: Optional[str | Path] = None,
    force: bool = False,
) -> ParseSummary:
    """Parse downloaded files for a list of accessions.

    Looks for FASTA and GenBank files in the deterministic layout
    created by the downloader: data_dir/sequences/fasta/{acc}.fasta
    and data_dir/sequences/genbank/{acc}.gb.

    GFF3 files are looked up in gff3_dir/{acc}.gff3 if provided.

    Args:
        session: DB session.
        accession_ids: List of accession IDs to parse.
        data_dir: Base data directory.
        gff3_dir: Optional directory containing GFF3 files.
        force: Re-parse even if already parsed.

    Returns:
        ParseSummary with counts and errors.
    """
    data_dir = Path(data_dir)
    summary = ParseSummary(total_accessions=len(accession_ids))

    for acc in accession_ids:
        fasta_path = data_dir / "sequences" / "fasta" / f"{acc}.fasta"
        genbank_path = data_dir / "sequences" / "genbank" / f"{acc}.gb"
        gff3_path = None
        if gff3_dir:
            gff3_path = str(Path(gff3_dir) / f"{acc}.gff3")

        result = parse_and_persist(
            session=session,
            accession=acc,
            fasta_path=str(fasta_path) if fasta_path.exists() else None,
            genbank_path=str(genbank_path) if genbank_path.exists() else None,
            gff3_path=gff3_path,
            force=force,
        )

        if result.get("skipped"):
            continue

        if result["fasta_parsed"]:
            summary.fasta_parsed += 1
        elif str(fasta_path) and fasta_path.exists():
            summary.fasta_failed += 1

        if result["genbank_parsed"]:
            summary.genbank_parsed += 1
        elif str(genbank_path) and genbank_path.exists():
            summary.genbank_failed += 1

        if result["gff3_parsed"]:
            summary.gff3_parsed += 1
        elif gff3_path and Path(gff3_path).exists():
            summary.gff3_failed += 1

        summary.features_inserted += result["features_inserted"]
        summary.errors.extend(result["errors"])

    return summary


# ---------------------------------------------------------------------------
# Batch-Aware Parsing (Chunk 4)
# ---------------------------------------------------------------------------


def parse_composite_artifact(
    session: Session,
    composite_fasta_path: Optional[str | Path] = None,
    composite_genbank_path: Optional[str | Path] = None,
    normalization_dir: Optional[str | Path] = None,
    force: bool = False,
) -> ParseSummary:
    """Parse composite FASTA/GenBank artifact and persist all records.

    This is the batch-aware entry point for parsing. It can handle:
    1. Composite files directly (auto-split via normalization)
    2. Pre-normalized directories (split files from Chunk 3)

    Args:
        session: DB session.
        composite_fasta_path: Path to composite FASTA file (will be split).
        composite_genbank_path: Path to composite GenBank file (will be split).
        normalization_dir: If provided and composite files given, split into this dir.
        force: Re-parse even if already parsed.

    Returns:
        ParseSummary with counts and errors.
    """
    from gwico_ssr.ingest.normalizers import (
        split_composite_fasta,
        split_composite_genbank,
    )

    summary = ParseSummary()
    normalized_fasta_files = {}  # accession -> path
    normalized_genbank_files = {}  # accession -> path

    # Handle composite FASTA
    if composite_fasta_path:
        composite_fasta_path = Path(composite_fasta_path)
        if not normalization_dir:
            normalization_dir = composite_fasta_path.parent / "normalized"

        split_result = split_composite_fasta(
            composite_fasta_path,
            normalization_dir,
            overwrite=force,
        )

        for rec in split_result.normalized_records:
            normalized_fasta_files[rec["accession"]] = rec["output_path"]

        for err in split_result.errors:
            summary.errors.append({
                "source": "composite_fasta",
                "accession": err.get("accession"),
                "message": err.get("message", str(err)),
            })

    # Handle composite GenBank
    if composite_genbank_path:
        composite_genbank_path = Path(composite_genbank_path)
        if not normalization_dir:
            normalization_dir = composite_genbank_path.parent / "normalized"

        split_result = split_composite_genbank(
            composite_genbank_path,
            normalization_dir,
            overwrite=force,
        )

        for rec in split_result.normalized_records:
            normalized_genbank_files[rec["accession"]] = rec["output_path"]

        for err in split_result.errors:
            summary.errors.append({
                "source": "composite_genbank",
                "accession": err.get("accession"),
                "message": err.get("message", str(err)),
            })

    # Parse all normalized files
    all_accessions = set(normalized_fasta_files.keys()) | set(
        normalized_genbank_files.keys()
    )
    summary.total_accessions = len(all_accessions)

    for acc in all_accessions:
        fasta_path = normalized_fasta_files.get(acc)
        genbank_path = normalized_genbank_files.get(acc)

        result = parse_and_persist(
            session=session,
            accession=acc,
            fasta_path=fasta_path,
            genbank_path=genbank_path,
            gff3_path=None,
            force=force,
        )

        if result.get("skipped"):
            continue

        if result["fasta_parsed"]:
            summary.fasta_parsed += 1
        elif fasta_path:
            summary.fasta_failed += 1

        if result["genbank_parsed"]:
            summary.genbank_parsed += 1
        elif genbank_path:
            summary.genbank_failed += 1

        summary.features_inserted += result["features_inserted"]
        summary.errors.extend(result["errors"])

    return summary


def parse_batch_directory(
    session: Session,
    batch_dir: str | Path,
    file_type: str = "fasta",
    gff3_dir: Optional[str | Path] = None,
    force: bool = False,
) -> ParseSummary:
    """Parse all normalized files in a batch directory.

    Looks for normalized single-record files in batch_dir/{file_type}/*.{fasta,gb}.

    Args:
        session: DB session.
        batch_dir: Directory containing normalized files.
        file_type: "fasta" or "genbank" (or "both").
        gff3_dir: Optional directory containing GFF3 files.
        force: Re-parse even if already parsed.

    Returns:
        ParseSummary with counts and errors.
    """
    batch_dir = Path(batch_dir)
    summary = ParseSummary()

    # Collect all accessions from available file types
    accessions_to_parse = set()

    if file_type in ("fasta", "both"):
        fasta_dir = batch_dir / "fasta"
        if fasta_dir.exists():
            for fasta_file in fasta_dir.glob("*.fasta"):
                accession = fasta_file.stem
                accessions_to_parse.add(accession)

    if file_type in ("genbank", "both"):
        genbank_dir = batch_dir / "genbank"
        if genbank_dir.exists():
            for gb_file in genbank_dir.glob("*.gb"):
                accession = gb_file.stem
                accessions_to_parse.add(accession)

    summary.total_accessions = len(accessions_to_parse)

    # Parse each accession
    for acc in sorted(accessions_to_parse):
        fasta_path = batch_dir / "fasta" / f"{acc}.fasta"
        genbank_path = batch_dir / "genbank" / f"{acc}.gb"
        gff3_path = None
        if gff3_dir:
            gff3_path = str(Path(gff3_dir) / f"{acc}.gff3")

        result = parse_and_persist(
            session=session,
            accession=acc,
            fasta_path=str(fasta_path) if fasta_path.exists() else None,
            genbank_path=str(genbank_path) if genbank_path.exists() else None,
            gff3_path=gff3_path,
            force=force,
        )

        if result.get("skipped"):
            continue

        if result["fasta_parsed"]:
            summary.fasta_parsed += 1
        elif fasta_path.exists():
            summary.fasta_failed += 1

        if result["genbank_parsed"]:
            summary.genbank_parsed += 1
        elif genbank_path.exists():
            summary.genbank_failed += 1

        if result["gff3_parsed"]:
            summary.gff3_parsed += 1
        elif gff3_path and Path(gff3_path).exists():
            summary.gff3_failed += 1

        summary.features_inserted += result["features_inserted"]
        summary.errors.extend(result["errors"])

    return summary
