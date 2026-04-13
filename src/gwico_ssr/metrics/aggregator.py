"""Cohort-level SSR metric aggregations.

Aggregates per-accession metrics into cohort summaries grouped by
country, gene, motif, motif_size, and dataset. All queries run against
the database to avoid loading full datasets into memory.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import Optional

from sqlalchemy import func, select
from sqlalchemy.orm import Session

from gwico_ssr.models.schema import (
    Accession,
    AccessionMetrics,
    SSRAnnotation,
    SSRRecord,
)

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Cohort summary data structures
# ---------------------------------------------------------------------------

@dataclass
class CohortSummary:
    """Generic grouping summary."""

    group_key: str
    group_value: str
    count: int
    total_ssrs: int = 0
    total_ssr_bp: int = 0
    mean_ra: float | None = None
    mean_rd: float | None = None


@dataclass
class MotifFrequency:
    """Motif frequency across the cohort."""

    motif_canonical: str
    motif_size: int
    total_count: int
    total_bp: int


@dataclass
class GeneSSRSummary:
    """SSR count per gene across the cohort."""

    gene_name: str
    ssr_count: int
    accession_count: int


# ---------------------------------------------------------------------------
# Database-side aggregation queries
# ---------------------------------------------------------------------------

def aggregate_by_country(session: Session, run_id: int) -> list[CohortSummary]:
    """Aggregate accession metrics by country for a given run.

    Returns one CohortSummary per country with accession count,
    total SSRs, total SSR bp, and mean RA/RD.
    """
    stmt = (
        select(
            Accession.country,
            func.count(AccessionMetrics.id).label("count"),
            func.sum(AccessionMetrics.ssr_count_total).label("total_ssrs"),
            func.sum(AccessionMetrics.ssr_bp_total).label("total_ssr_bp"),
            func.avg(AccessionMetrics.ra).label("mean_ra"),
            func.avg(AccessionMetrics.rd).label("mean_rd"),
        )
        .join(AccessionMetrics, Accession.accession == AccessionMetrics.accession)
        .where(AccessionMetrics.run_id == run_id)
        .where(Accession.country.isnot(None))
        .group_by(Accession.country)
        .order_by(func.count(AccessionMetrics.id).desc())
    )
    rows = session.execute(stmt).all()
    return [
        CohortSummary(
            group_key="country",
            group_value=r.country,
            count=r.count,
            total_ssrs=r.total_ssrs or 0,
            total_ssr_bp=r.total_ssr_bp or 0,
            mean_ra=r.mean_ra,
            mean_rd=r.mean_rd,
        )
        for r in rows
    ]


def aggregate_by_motif_size(session: Session, run_id: int) -> list[CohortSummary]:
    """Aggregate SSR counts by motif size across all accessions in a run."""
    # Use the per-accession motif-size counts from AccessionMetrics
    size_labels = {
        1: "mono",
        2: "di",
        3: "tri",
        4: "tetra",
        5: "penta",
        6: "hexa",
    }
    size_columns = {
        1: AccessionMetrics.mono_count,
        2: AccessionMetrics.di_count,
        3: AccessionMetrics.tri_count,
        4: AccessionMetrics.tetra_count,
        5: AccessionMetrics.penta_count,
        6: AccessionMetrics.hexa_count,
    }
    results = []
    for size, col in size_columns.items():
        stmt = (
            select(
                func.count(AccessionMetrics.id).label("count"),
                func.sum(col).label("total_ssrs"),
            )
            .where(AccessionMetrics.run_id == run_id)
        )
        row = session.execute(stmt).one()
        results.append(CohortSummary(
            group_key="motif_size",
            group_value=size_labels[size],
            count=row.count,
            total_ssrs=row.total_ssrs or 0,
        ))
    return results


def aggregate_motif_frequencies(session: Session, dataset_id: int | None = None) -> list[MotifFrequency]:
    """Count SSR occurrences by canonical motif across the dataset.

    Returns motifs sorted by total count descending.
    """
    stmt = (
        select(
            SSRRecord.motif_canonical,
            SSRRecord.motif_size,
            func.count(SSRRecord.ssr_id).label("total_count"),
            func.sum(SSRRecord.repeat_length_bp).label("total_bp"),
        )
        .group_by(SSRRecord.motif_canonical, SSRRecord.motif_size)
        .order_by(func.count(SSRRecord.ssr_id).desc())
    )
    if dataset_id is not None:
        stmt = stmt.join(Accession, SSRRecord.accession == Accession.accession).where(
            Accession.dataset_id == dataset_id
        )
    rows = session.execute(stmt).all()
    return [
        MotifFrequency(
            motif_canonical=r.motif_canonical,
            motif_size=r.motif_size,
            total_count=r.total_count,
            total_bp=r.total_bp or 0,
        )
        for r in rows
    ]


def aggregate_by_gene(session: Session, dataset_id: int | None = None) -> list[GeneSSRSummary]:
    """Count SSRs per gene across the dataset, using SSRAnnotation records.

    Returns genes sorted by SSR count descending.
    """
    stmt = (
        select(
            SSRAnnotation.gene_name,
            func.count(SSRAnnotation.annotation_id).label("ssr_count"),
            func.count(func.distinct(SSRAnnotation.accession)).label("accession_count"),
        )
        .where(SSRAnnotation.gene_name.isnot(None))
        .group_by(SSRAnnotation.gene_name)
        .order_by(func.count(SSRAnnotation.annotation_id).desc())
    )
    if dataset_id is not None:
        stmt = stmt.join(Accession, SSRAnnotation.accession == Accession.accession).where(
            Accession.dataset_id == dataset_id
        )
    rows = session.execute(stmt).all()
    return [
        GeneSSRSummary(
            gene_name=r.gene_name,
            ssr_count=r.ssr_count,
            accession_count=r.accession_count,
        )
        for r in rows
    ]


def get_dataset_summary(session: Session, run_id: int) -> dict:
    """Build a comprehensive dataset summary for a run.

    Returns a dict with total accessions, SSR counts, motif-size breakdown,
    country breakdown, gene breakdown, and top motifs.
    """
    # Total metrics
    totals_stmt = (
        select(
            func.count(AccessionMetrics.id).label("total_accessions"),
            func.sum(AccessionMetrics.ssr_count_total).label("total_ssrs"),
            func.sum(AccessionMetrics.ssr_bp_total).label("total_ssr_bp"),
            func.avg(AccessionMetrics.ra).label("mean_ra"),
            func.avg(AccessionMetrics.rd).label("mean_rd"),
            func.avg(AccessionMetrics.ssr_count_total).label("mean_ssrs_per_accession"),
        )
        .where(AccessionMetrics.run_id == run_id)
    )
    totals = session.execute(totals_stmt).one()

    by_size = aggregate_by_motif_size(session, run_id)
    by_country = aggregate_by_country(session, run_id)

    return {
        "run_id": run_id,
        "total_accessions": totals.total_accessions,
        "total_ssrs": totals.total_ssrs or 0,
        "total_ssr_bp": totals.total_ssr_bp or 0,
        "mean_ra": totals.mean_ra,
        "mean_rd": totals.mean_rd,
        "mean_ssrs_per_accession": totals.mean_ssrs_per_accession,
        "by_motif_size": [
            {"motif_size": s.group_value, "total_ssrs": s.total_ssrs}
            for s in by_size
        ],
        "by_country_count": len(by_country),
        "top_countries": [
            {"country": c.group_value, "accessions": c.count, "total_ssrs": c.total_ssrs}
            for c in by_country[:10]
        ],
    }
