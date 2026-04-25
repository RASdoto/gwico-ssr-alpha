"""Lineage-aware analysis for GWICO-SSR.

Implements lineage ingestion, reconciliation, and lineage-aware metrics with
full audit trail and reproducibility guarantees.
"""

from __future__ import annotations

import csv
import json
import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

from sqlalchemy import func, select, case
from sqlalchemy.orm import Session

from gwico_ssr.models.schema import (
    Accession,
    AccessionMetrics,
    LineageAssignment as LineageAssignmentTable,
    LineageMetrics as LineageMetricsTable,
    LineageSource as LineageSourceTable,
    SSRRecord,
)

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Dataclasses
# ---------------------------------------------------------------------------

@dataclass
class LineageSource:
    """Metadata about a lineage reference source."""
    source_name: str  # pango, nextstrain, user-defined
    version: str  # e.g., "v1.23", "2024-Q1"
    release_date: Optional[str] = None  # ISO 8601 date
    description: Optional[str] = None


@dataclass
class LineageAssignment:
    """Single accession-to-lineage mapping."""
    accession: str
    lineage_label: str  # e.g., "BA.1", "XEC", "unassigned"
    source: str  # source_name from LineageSource
    confidence: Optional[float] = None  # 0.0-1.0 if available
    assignment_date: Optional[str] = None
    is_deprecated: bool = False
    notes: Optional[str] = None


@dataclass
class LineageSummary:
    """Aggregated metrics for a single lineage."""
    lineage_label: str
    accession_count: int = 0
    ssr_count_total: int = 0
    ssr_bp_total: int = 0
    perfect_count: int = 0
    imperfect_count: int = 0
    compound_component_count: int = 0
    perfect_bp_total: int = 0
    imperfect_bp_total: int = 0
    compound_component_bp_total: int = 0
    mean_ra: Optional[float] = None
    mean_rd: Optional[float] = None


@dataclass
class LineageReport:
    """Complete lineage assignment report with reconciliation info."""
    source: LineageSource
    total_accessions_in_source: int
    total_accessions_in_database: int
    successfully_mapped: int
    unmapped: int
    deprecated_labels: list[str] = field(default_factory=list)
    unassigned_count: int = 0
    conflicts: list[str] = field(default_factory=list)
    analysis_notes: list[str] = field(default_factory=list)


# ---------------------------------------------------------------------------
# Lineage ingestion
# ---------------------------------------------------------------------------

def ingest_pango_csv(
    session: Session,
    csv_path: str | Path,
    source_name: str = "pango",
    version: str = "user-provided",
    release_date: Optional[str] = None,
) -> LineageReport:
    """Ingest Pango lineage assignments from CSV.

    Expected CSV columns: accession, lineage_label, confidence (optional)

    Args:
        session: SQLAlchemy session
        csv_path: Path to CSV file
        source_name: Source identifier (e.g., "pango", "pango-designated")
        version: Version string
        release_date: Optional ISO 8601 release date

    Returns:
        LineageReport with ingestion summary
    """
    csv_path = Path(csv_path)
    if not csv_path.exists():
        raise FileNotFoundError(f"CSV file not found: {csv_path}")

    # Create or get lineage source
    source_obj = session.query(LineageSourceTable).filter(
        LineageSourceTable.source_name == source_name,
        LineageSourceTable.version == version,
    ).first()

    if not source_obj:
        source_obj = LineageSourceTable(
            source_name=source_name,
            version=version,
            release_date=release_date,
        )
        session.add(source_obj)
        session.flush()

    # Read CSV and ingest assignments
    report = LineageReport(
        source=LineageSource(
            source_name=source_name,
            version=version,
            release_date=release_date,
        ),
        total_accessions_in_source=0,
        total_accessions_in_database=0,
        successfully_mapped=0,
        unmapped=0,
    )

    with open(csv_path, "r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            accession = row.get("accession")
            lineage_label = row.get("lineage_label") or row.get("pango_lineage")
            confidence = float(row.get("confidence", 0.0)) if row.get("confidence") else None

            if not accession or not lineage_label:
                continue

            report.total_accessions_in_source += 1

            # Check if accession exists in database
            acc_obj = session.query(Accession).filter(Accession.accession == accession).first()
            if not acc_obj:
                report.unmapped += 1
                report.analysis_notes.append(f"Accession {accession} not in database")
                continue

            report.total_accessions_in_database += 1

            # Create or update assignment
            existing = session.query(LineageAssignmentTable).filter(
                LineageAssignmentTable.accession == accession,
                LineageAssignmentTable.source_id == source_obj.source_id,
            ).first()

            if existing:
                # Update confidence if provided
                if confidence is not None:
                    existing.confidence = confidence
                logger.info(f"Updated assignment for {accession}")
            else:
                assignment = LineageAssignmentTable(
                    accession=accession,
                    source_id=source_obj.source_id,
                    lineage_label=lineage_label,
                    confidence=confidence,
                )
                session.add(assignment)
                report.successfully_mapped += 1

            # Track unassigned
            if lineage_label.lower() in ("unassigned", "unknown", "none"):
                report.unassigned_count += 1

    session.commit()
    logger.info(f"Ingested {report.successfully_mapped} lineage assignments from {csv_path}")

    return report


def ingest_nextstrain_json(
    session: Session,
    json_path: str | Path,
    source_name: str = "nextstrain",
    version: str = "user-provided",
    release_date: Optional[str] = None,
) -> LineageReport:
    """Ingest lineage assignments from Nextstrain JSON export.

    Expected structure: {accessions: {accession: {clade: "XYZ", ...}, ...}}

    Args:
        session: SQLAlchemy session
        json_path: Path to JSON file
        source_name: Source identifier
        version: Version string
        release_date: Optional ISO 8601 release date

    Returns:
        LineageReport with ingestion summary
    """
    json_path = Path(json_path)
    if not json_path.exists():
        raise FileNotFoundError(f"JSON file not found: {json_path}")

    # Create or get lineage source
    source_obj = session.query(LineageSourceTable).filter(
        LineageSourceTable.source_name == source_name,
        LineageSourceTable.version == version,
    ).first()

    if not source_obj:
        source_obj = LineageSourceTable(
            source_name=source_name,
            version=version,
            release_date=release_date,
        )
        session.add(source_obj)
        session.flush()

    # Load and parse JSON
    report = LineageReport(
        source=LineageSource(
            source_name=source_name,
            version=version,
            release_date=release_date,
        ),
        total_accessions_in_source=0,
        total_accessions_in_database=0,
        successfully_mapped=0,
        unmapped=0,
    )

    with open(json_path, "r") as f:
        data = json.load(f)

    accessions_dict = data.get("accessions", {})
    for accession, metadata in accessions_dict.items():
        clade = metadata.get("clade")
        if not clade:
            continue

        report.total_accessions_in_source += 1

        # Check if accession exists
        acc_obj = session.query(Accession).filter(Accession.accession == accession).first()
        if not acc_obj:
            report.unmapped += 1
            report.analysis_notes.append(f"Accession {accession} not in database")
            continue

        report.total_accessions_in_database += 1

        # Create or update assignment
        existing = session.query(LineageAssignmentTable).filter(
            LineageAssignmentTable.accession == accession,
            LineageAssignmentTable.source_id == source_obj.source_id,
        ).first()

        if existing:
            logger.info(f"Updated assignment for {accession}")
        else:
            assignment = LineageAssignmentTable(
                accession=accession,
                source_id=source_obj.source_id,
                lineage_label=clade,
                confidence=metadata.get("confidence"),
            )
            session.add(assignment)
            report.successfully_mapped += 1

    session.commit()
    logger.info(f"Ingested {report.successfully_mapped} lineage assignments from {json_path}")

    return report


# ---------------------------------------------------------------------------
# Lineage reconciliation
# ---------------------------------------------------------------------------

def get_lineage_for_accession(
    session: Session,
    accession: str,
    source_id: Optional[int] = None,
) -> Optional[LineageAssignmentTable]:
    """Get the most recent lineage assignment for an accession.

    Args:
        session: SQLAlchemy session
        accession: Accession identifier
        source_id: Optional filter for specific source

    Returns:
        Most recent LineageAssignment or None
    """
    query = session.query(LineageAssignmentTable).filter(
        LineageAssignmentTable.accession == accession
    )

    if source_id is not None:
        query = query.filter(LineageAssignmentTable.source_id == source_id)

    # Get most recent (latest created_at)
    result = query.order_by(LineageAssignmentTable.created_at.desc()).first()
    return result


def audit_lineage_mapping(
    session: Session,
    source_id: int,
) -> LineageReport:
    """Generate audit report for a lineage source.

    Args:
        session: SQLAlchemy session
        source_id: LineageSource ID

    Returns:
        LineageReport with audit details
    """
    source_obj = session.query(LineageSourceTable).filter(
        LineageSourceTable.source_id == source_id
    ).first()

    if not source_obj:
        raise ValueError(f"LineageSource {source_id} not found")

    # Count assignments
    total_assigned = session.query(func.count(LineageAssignmentTable.assignment_id)).filter(
        LineageAssignmentTable.source_id == source_id
    ).scalar()

    deprecated_count = session.query(func.count(LineageAssignmentTable.assignment_id)).filter(
        LineageAssignmentTable.source_id == source_id,
        LineageAssignmentTable.is_deprecated == True,
    ).scalar()

    unassigned_count = session.query(func.count(LineageAssignmentTable.assignment_id)).filter(
        LineageAssignmentTable.source_id == source_id,
        LineageAssignmentTable.lineage_label.in_(["unassigned", "unknown", "none"]),
    ).scalar()

    # Unique lineage labels
    unique_labels = session.query(
        func.distinct(LineageAssignmentTable.lineage_label)
    ).filter(
        LineageAssignmentTable.source_id == source_id
    ).all()

    report = LineageReport(
        source=LineageSource(
            source_name=source_obj.source_name,
            version=source_obj.version,
            release_date=source_obj.release_date,
        ),
        total_accessions_in_source=total_assigned,
        total_accessions_in_database=total_assigned,
        successfully_mapped=total_assigned - deprecated_count,
        unmapped=deprecated_count,
        unassigned_count=unassigned_count,
        deprecated_labels=[label[0] for label in session.query(
            LineageAssignmentTable.lineage_label
        ).filter(
            LineageAssignmentTable.source_id == source_id,
            LineageAssignmentTable.is_deprecated == True,
        ).all()],
    )

    return report


# ---------------------------------------------------------------------------
# Lineage-aware metrics
# ---------------------------------------------------------------------------

def compute_ssr_metrics_by_lineage(
    session: Session,
    source_id: int,
    run_id: Optional[int] = None,
) -> dict[str, LineageSummary]:
    """Compute SSR metrics aggregated by lineage with repeat-class breakdown.

    Args:
        session: SQLAlchemy session
        source_id: LineageSource ID
        run_id: Optional run_id filter

    Returns:
        Dict mapping lineage_label → LineageSummary
    """
    # Get all lineage labels for source
    labels = session.query(
        func.distinct(LineageAssignmentTable.lineage_label)
    ).filter(
        LineageAssignmentTable.source_id == source_id
    ).all()

    results = {}

    for label_tuple in labels:
        label = label_tuple[0]

        # Get accessions for this lineage
        accessions_for_lineage = session.query(LineageAssignmentTable.accession).filter(
            LineageAssignmentTable.source_id == source_id,
            LineageAssignmentTable.lineage_label == label,
        ).all()

        accession_list = [a[0] for a in accessions_for_lineage]
        accession_count = len(accession_list)

        if accession_count == 0:
            continue

        # Query SSR metrics for these accessions
        stmt = select(
            func.count(SSRRecord.ssr_id).label("ssr_count"),
            func.sum(SSRRecord.repeat_length_bp).label("ssr_bp_total"),
            func.sum(case((SSRRecord.repeat_class == "perfect", 1), else_=0)).label("perfect_count"),
            func.sum(case((SSRRecord.repeat_class == "imperfect", 1), else_=0)).label("imperfect_count"),
            func.sum(case((SSRRecord.repeat_class == "compound_component", 1), else_=0)).label("compound_count"),
            func.sum(case((SSRRecord.repeat_class == "perfect", SSRRecord.repeat_length_bp), else_=0)).label("perfect_bp"),
            func.sum(case((SSRRecord.repeat_class == "imperfect", SSRRecord.repeat_length_bp), else_=0)).label("imperfect_bp"),
            func.sum(case((SSRRecord.repeat_class == "compound_component", SSRRecord.repeat_length_bp), else_=0)).label("compound_bp"),
        ).where(
            SSRRecord.accession.in_(accession_list)
        )

        if run_id is not None:
            stmt = stmt.where(SSRRecord.run_id == run_id)

        ssr_result = session.execute(stmt).one()

        # Query metrics for mean RA/RD
        stmt_metrics = select(
            func.avg(AccessionMetrics.ra).label("mean_ra"),
            func.avg(AccessionMetrics.rd).label("mean_rd"),
        ).where(
            AccessionMetrics.accession.in_(accession_list)
        )

        if run_id is not None:
            stmt_metrics = stmt_metrics.where(AccessionMetrics.run_id == run_id)

        metrics_result = session.execute(stmt_metrics).one()

        summary = LineageSummary(
            lineage_label=label,
            accession_count=accession_count,
            ssr_count_total=ssr_result.ssr_count or 0,
            ssr_bp_total=ssr_result.ssr_bp_total or 0,
            perfect_count=ssr_result.perfect_count or 0,
            imperfect_count=ssr_result.imperfect_count or 0,
            compound_component_count=ssr_result.compound_count or 0,
            perfect_bp_total=int(ssr_result.perfect_bp) if ssr_result.perfect_bp else 0,
            imperfect_bp_total=int(ssr_result.imperfect_bp) if ssr_result.imperfect_bp else 0,
            compound_component_bp_total=int(ssr_result.compound_bp) if ssr_result.compound_bp else 0,
            mean_ra=metrics_result.mean_ra,
            mean_rd=metrics_result.mean_rd,
        )

        results[label] = summary

    return results


def store_lineage_metrics(
    session: Session,
    run_id: int,
    source_id: int,
    metrics_dict: dict[str, LineageSummary],
) -> int:
    """Store lineage metrics to database for persistence and caching.

    Args:
        session: SQLAlchemy session
        run_id: Run ID
        source_id: LineageSource ID
        metrics_dict: Dictionary of lineage_label → LineageSummary

    Returns:
        Number of metric records created
    """
    count = 0
    for lineage_label, summary in metrics_dict.items():
        metric = LineageMetricsTable(
            run_id=run_id,
            source_id=source_id,
            lineage_label=lineage_label,
            accession_count=summary.accession_count,
            ssr_count_total=summary.ssr_count_total,
            ssr_bp_total=summary.ssr_bp_total,
            perfect_count=summary.perfect_count,
            imperfect_count=summary.imperfect_count,
            compound_component_count=summary.compound_component_count,
            perfect_bp_total=summary.perfect_bp_total,
            imperfect_bp_total=summary.imperfect_bp_total,
            compound_component_bp_total=summary.compound_component_bp_total,
            mean_ra=summary.mean_ra,
            mean_rd=summary.mean_rd,
        )
        session.add(metric)
        count += 1

    session.commit()
    logger.info(f"Stored {count} lineage metrics records for run {run_id}")
    return count
