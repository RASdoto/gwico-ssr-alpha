"""Temporal analysis for GWICO-SSR.

Implements time-binned SSR analysis on persisted run state, including:
- Date parsing and normalization with precision levels
- Configurable time binning (year, month, quarter, week, day)
- Temporal distribution queries
- Repeat-class-aware temporal metrics
- Motif temporal trends

All queries run on database side (no full-dataset loads).
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from datetime import datetime, timedelta
from enum import Enum
from typing import Optional

from sqlalchemy import func, select
from sqlalchemy.orm import Session

from gwico_ssr.models.schema import (
    Accession,
    AccessionMetrics,
    SSRRecord,
)

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Enums and dataclasses
# ---------------------------------------------------------------------------

class DatePrecision(Enum):
    """Temporal precision levels for binning."""
    YEAR = "year"
    QUARTER = "quarter"
    MONTH = "month"
    WEEK = "week"
    DAY = "day"


@dataclass
class TimeBin:
    """A time bin with aggregated SSR metrics."""

    bin_id: str  # e.g., "2024-Q1", "2024-01", "2024"
    start_date: datetime
    end_date: datetime
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
    excluded_accessions: int = 0


@dataclass
class TemporalSummary:
    """Collection of temporal bins with metadata."""

    date_field: str  # "release_date" or "collection_date"
    precision: DatePrecision
    bins: list[TimeBin] = field(default_factory=list)
    total_accessions: int = 0
    excluded_accessions: int = 0
    date_range_start: Optional[datetime] = None
    date_range_end: Optional[datetime] = None
    analysis_notes: list[str] = field(default_factory=list)

    def total_ssrs(self) -> int:
        """Sum total SSRs across all bins."""
        return sum(b.ssr_count_total for b in self.bins)

    def total_ssrs_bp(self) -> int:
        """Sum total SSR bp across all bins."""
        return sum(b.ssr_bp_total for b in self.bins)


# ---------------------------------------------------------------------------
# Date parsing and normalization
# ---------------------------------------------------------------------------

def parse_date_string(date_str: str) -> Optional[datetime]:
    """Parse ISO 8601 or partial date strings to datetime.

    Supports:
    - Full: 2024-04-25, 2024-04-25T10:30:00
    - Partial: 2024-04, 2024-Q1, 2024

    Returns None if parsing fails.
    """
    if not date_str or not isinstance(date_str, str):
        return None

    date_str = date_str.strip()

    # Full date with time
    for fmt in ["%Y-%m-%dT%H:%M:%S", "%Y-%m-%d %H:%M:%S", "%Y-%m-%d"]:
        try:
            return datetime.strptime(date_str, fmt)
        except ValueError:
            continue

    # Partial dates
    try:
        # YYYY-MM
        if len(date_str) == 7 and date_str[4] == "-":
            # Check if it's YYYY-MM or YYYY-Qq
            if date_str[5] == "Q":
                # Quarter format
                year = int(date_str[:4])
                quarter = int(date_str[6])
                if quarter < 1 or quarter > 4:
                    return None
                month = (quarter - 1) * 3 + 1
                return datetime(year, month, 1)
            else:
                # Month format
                return datetime.strptime(date_str, "%Y-%m")
        # YYYY only
        if len(date_str) == 4 and date_str.isdigit():
            return datetime.strptime(date_str, "%Y")
    except (ValueError, IndexError):
        pass

    return None


def normalize_date_to_precision(
    date: datetime,
    precision: DatePrecision,
) -> datetime:
    """Normalize a datetime to a precision level.

    YEAR: Jan 1 of year
    QUARTER: First day of quarter (Jan 1, Apr 1, Jul 1, Oct 1)
    MONTH: First day of month
    WEEK: Monday of that week
    DAY: Midnight of that day
    """
    if precision == DatePrecision.YEAR:
        return date.replace(month=1, day=1, hour=0, minute=0, second=0, microsecond=0)
    elif precision == DatePrecision.QUARTER:
        quarter = (date.month - 1) // 3
        month = quarter * 3 + 1
        return date.replace(month=month, day=1, hour=0, minute=0, second=0, microsecond=0)
    elif precision == DatePrecision.MONTH:
        return date.replace(day=1, hour=0, minute=0, second=0, microsecond=0)
    elif precision == DatePrecision.WEEK:
        # Monday = weekday 0
        days_since_monday = date.weekday()
        week_start = date - timedelta(days=days_since_monday)
        return week_start.replace(hour=0, minute=0, second=0, microsecond=0)
    elif precision == DatePrecision.DAY:
        return date.replace(hour=0, minute=0, second=0, microsecond=0)
    else:
        return date


def bin_id_from_date(date: datetime, precision: DatePrecision) -> str:
    """Generate a bin ID string from a date and precision.

    Examples:
    - YEAR: "2024"
    - QUARTER: "2024-Q1"
    - MONTH: "2024-04"
    - WEEK: "2024-W17"
    - DAY: "2024-04-25"
    """
    if precision == DatePrecision.YEAR:
        return date.strftime("%Y")
    elif precision == DatePrecision.QUARTER:
        quarter = (date.month - 1) // 3 + 1
        return f"{date.year}-Q{quarter}"
    elif precision == DatePrecision.MONTH:
        return date.strftime("%Y-%m")
    elif precision == DatePrecision.WEEK:
        return date.strftime("%Y-W%V")
    elif precision == DatePrecision.DAY:
        return date.strftime("%Y-%m-%d")
    else:
        return str(date)


# ---------------------------------------------------------------------------
# Time bin generation
# ---------------------------------------------------------------------------

def create_time_bins(
    start_date: datetime,
    end_date: datetime,
    precision: DatePrecision,
    min_accessions: int = 0,
) -> list[tuple[datetime, datetime]]:
    """Generate time bin boundaries from start to end date.

    Returns list of (bin_start, bin_end) tuples covering the date range.
    Each bin spans to the first moment of the next bin (half-open intervals).

    Args:
        start_date: Bin generation starts here
        end_date: Bin generation ends here
        precision: Granularity (YEAR, MONTH, WEEK, DAY, etc.)
        min_accessions: Not used for generation; for filtering bins later

    Returns:
        List of (start, end) datetime tuples representing bin boundaries
    """
    bins = []
    current = normalize_date_to_precision(start_date, precision)

    while current <= end_date:
        if precision == DatePrecision.YEAR:
            next_date = current.replace(year=current.year + 1)
        elif precision == DatePrecision.QUARTER:
            quarter = (current.month - 1) // 3
            if quarter == 3:
                next_date = current.replace(year=current.year + 1, month=1)
            else:
                next_date = current.replace(month=current.month + 3)
        elif precision == DatePrecision.MONTH:
            if current.month == 12:
                next_date = current.replace(year=current.year + 1, month=1)
            else:
                next_date = current.replace(month=current.month + 1)
        elif precision == DatePrecision.WEEK:
            next_date = current + timedelta(days=7)
        elif precision == DatePrecision.DAY:
            next_date = current + timedelta(days=1)
        else:
            break

        bins.append((current, next_date))
        current = next_date

    return bins


# ---------------------------------------------------------------------------
# Temporal queries
# ---------------------------------------------------------------------------

def get_date_range_for_field(
    session: Session,
    date_field: str,
    run_id: Optional[int] = None,
) -> tuple[Optional[datetime], Optional[datetime]]:
    """Find min and max dates in a date field.

    Args:
        session: SQLAlchemy session
        date_field: "release_date" or "collection_date"
        run_id: Optional run_id to filter SSRs

    Returns:
        (min_date, max_date) or (None, None) if no data
    """
    col = getattr(Accession, date_field)
    stmt = select(
        func.min(col).label("min_date"),
        func.max(col).label("max_date"),
    )

    if run_id is not None:
        stmt = stmt.join(SSRRecord, Accession.accession == SSRRecord.accession).where(
            SSRRecord.run_id == run_id
        )

    result = session.execute(stmt).one()
    min_date = result.min_date
    max_date = result.max_date

    if min_date and isinstance(min_date, str):
        min_date = parse_date_string(min_date)
    if max_date and isinstance(max_date, str):
        max_date = parse_date_string(max_date)

    return min_date, max_date


def query_accessions_in_time_bin(
    session: Session,
    date_field: str,
    bin_start: datetime,
    bin_end: datetime,
    run_id: Optional[int] = None,
) -> int:
    """Count unique accessions in a time bin.

    Args:
        session: SQLAlchemy session
        date_field: "release_date" or "collection_date"
        bin_start: Bin start (inclusive)
        bin_end: Bin end (exclusive)
        run_id: Optional run_id filter

    Returns:
        Count of accessions
    """
    col = getattr(Accession, date_field)
    stmt = select(func.count(func.distinct(Accession.accession))).where(
        (col >= bin_start.isoformat()) & (col < bin_end.isoformat())
    )

    if run_id is not None:
        stmt = stmt.join(SSRRecord, Accession.accession == SSRRecord.accession).where(
            SSRRecord.run_id == run_id
        )

    count = session.execute(stmt).scalar()
    return count or 0


def query_temporal_ssr_metrics(
    session: Session,
    date_field: str,
    bin_start: datetime,
    bin_end: datetime,
    run_id: Optional[int] = None,
) -> dict:
    """Query SSR metrics for accessions in a time bin.

    Returns dict with: ssr_count, ssr_bp_total, perfect_count, imperfect_count,
    compound_component_count, perfect_bp_total, imperfect_bp_total,
    compound_component_bp_total, mean_ra, mean_rd
    """
    from sqlalchemy import case
    
    col = getattr(Accession, date_field)
    
    # Query all SSR metrics for the bin
    stmt = select(
        func.count(SSRRecord.ssr_id).label("ssr_count"),
        func.sum(SSRRecord.repeat_length_bp).label("ssr_bp_total"),
        func.sum(case((SSRRecord.repeat_class == "perfect", 1), else_=0)).label("perfect_count"),
        func.sum(case((SSRRecord.repeat_class == "imperfect", 1), else_=0)).label("imperfect_count"),
        func.sum(case((SSRRecord.repeat_class == "compound_component", 1), else_=0)).label("compound_count"),
        func.sum(case((SSRRecord.repeat_class == "perfect", SSRRecord.repeat_length_bp), else_=0)).label("perfect_bp"),
        func.sum(case((SSRRecord.repeat_class == "imperfect", SSRRecord.repeat_length_bp), else_=0)).label("imperfect_bp"),
        func.sum(case((SSRRecord.repeat_class == "compound_component", SSRRecord.repeat_length_bp), else_=0)).label("compound_bp"),
    ).join(
        Accession, SSRRecord.accession == Accession.accession
    ).where(
        (col >= bin_start.isoformat()) & (col < bin_end.isoformat())
    )

    if run_id is not None:
        stmt = stmt.where(SSRRecord.run_id == run_id)

    result = session.execute(stmt).one()

    # Mean RA/RD from metrics table
    stmt_metrics = select(
        func.avg(AccessionMetrics.ra).label("mean_ra"),
        func.avg(AccessionMetrics.rd).label("mean_rd"),
    ).join(
        Accession, AccessionMetrics.accession == Accession.accession
    ).where(
        (col >= bin_start.isoformat()) & (col < bin_end.isoformat())
    )
    if run_id is not None:
        stmt_metrics = stmt_metrics.where(AccessionMetrics.run_id == run_id)

    result_metrics = session.execute(stmt_metrics).one()

    return {
        "ssr_count_total": result.ssr_count or 0,
        "ssr_bp_total": result.ssr_bp_total or 0,
        "perfect_count": result.perfect_count or 0,
        "imperfect_count": result.imperfect_count or 0,
        "compound_component_count": result.compound_count or 0,
        "perfect_bp_total": int(result.perfect_bp) if result.perfect_bp else 0,
        "imperfect_bp_total": int(result.imperfect_bp) if result.imperfect_bp else 0,
        "compound_component_bp_total": int(result.compound_bp) if result.compound_bp else 0,
        "mean_ra": result_metrics.mean_ra,
        "mean_rd": result_metrics.mean_rd,
    }


def compute_temporal_summary(
    session: Session,
    date_field: str,
    precision: DatePrecision,
    run_id: Optional[int] = None,
) -> TemporalSummary:
    """Compute temporal summary with SSR metrics across time bins.

    Args:
        session: SQLAlchemy session
        date_field: "release_date" or "collection_date"
        precision: Temporal precision (YEAR, MONTH, etc.)
        run_id: Optional run_id filter

    Returns:
        TemporalSummary with populated bins
    """
    # Get date range
    min_date, max_date = get_date_range_for_field(session, date_field, run_id)

    if not min_date or not max_date:
        logger.warning(f"No dates found in {date_field}")
        return TemporalSummary(
            date_field=date_field,
            precision=precision,
            analysis_notes=["No data available for temporal analysis"],
        )

    # Generate bins
    bin_boundaries = create_time_bins(min_date, max_date, precision)

    summary = TemporalSummary(
        date_field=date_field,
        precision=precision,
        date_range_start=min_date,
        date_range_end=max_date,
    )

    for bin_start, bin_end in bin_boundaries:
        bin_id = bin_id_from_date(bin_start, precision)
        accession_count = query_accessions_in_time_bin(
            session, date_field, bin_start, bin_end, run_id
        )

        if accession_count == 0:
            summary.analysis_notes.append(f"Bin {bin_id}: 0 accessions, skipped")
            continue

        metrics = query_temporal_ssr_metrics(
            session, date_field, bin_start, bin_end, run_id
        )

        bin_obj = TimeBin(
            bin_id=bin_id,
            start_date=bin_start,
            end_date=bin_end,
            accession_count=accession_count,
            ssr_count_total=metrics["ssr_count_total"],
            ssr_bp_total=metrics["ssr_bp_total"],
            perfect_count=metrics["perfect_count"],
            imperfect_count=metrics["imperfect_count"],
            compound_component_count=metrics["compound_component_count"],
            perfect_bp_total=metrics["perfect_bp_total"],
            imperfect_bp_total=metrics["imperfect_bp_total"],
            compound_component_bp_total=metrics["compound_component_bp_total"],
            mean_ra=metrics["mean_ra"],
            mean_rd=metrics["mean_rd"],
        )
        summary.bins.append(bin_obj)
        summary.total_accessions += accession_count

    if summary.analysis_notes:
        logger.info(f"Temporal analysis notes: {'; '.join(summary.analysis_notes)}")

    return summary
