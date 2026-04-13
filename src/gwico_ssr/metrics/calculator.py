"""Per-accession SSR metric computation.

Metric definitions:
    ssr_count_total: Total number of SSR loci detected.
    ssr_bp_total: Total base pairs covered by SSR loci.
    mono_count .. hexa_count: SSR count by motif size (1-6).
    ra: Relative Abundance = ssr_count_total / (genome_length / 1000).
        Units: SSRs per kilobase.
    rd: Relative Density = ssr_bp_total / (genome_length / 1_000_000).
        Units: SSR bp per megabase.
    dominant_motif: The canonical motif with the highest count for this accession.

All metrics are derived from database state (SSRRecord table) for a given accession.
"""

from __future__ import annotations

import logging
from collections import Counter
from dataclasses import dataclass, field
from typing import Optional, Sequence

logger = logging.getLogger(__name__)


@dataclass
class AccessionMetricsResult:
    """Computed metrics for a single accession."""

    accession: str
    ssr_count_total: int = 0
    ssr_bp_total: int = 0
    mono_count: int = 0
    di_count: int = 0
    tri_count: int = 0
    tetra_count: int = 0
    penta_count: int = 0
    hexa_count: int = 0
    ra: float | None = None
    rd: float | None = None
    dominant_motif: str | None = None
    genome_length: int | None = None
    gc_content: float | None = None

    def to_dict(self, run_id: int) -> dict:
        """Convert to dict suitable for upsert_accession_metrics."""
        return {
            "accession": self.accession,
            "run_id": run_id,
            "ssr_count_total": self.ssr_count_total,
            "ssr_bp_total": self.ssr_bp_total,
            "ra": self.ra,
            "rd": self.rd,
            "mono_count": self.mono_count,
            "di_count": self.di_count,
            "tri_count": self.tri_count,
            "tetra_count": self.tetra_count,
            "penta_count": self.penta_count,
            "hexa_count": self.hexa_count,
            "dominant_motif": self.dominant_motif,
        }


# ---------------------------------------------------------------------------
# Motif-size count helpers
# ---------------------------------------------------------------------------

_SIZE_ATTR = {
    1: "mono_count",
    2: "di_count",
    3: "tri_count",
    4: "tetra_count",
    5: "penta_count",
    6: "hexa_count",
}


def compute_motif_size_counts(ssr_records) -> dict[str, int]:
    """Count SSRs by motif size from a sequence of SSRRecord objects.

    Returns a dict with keys mono_count..hexa_count.
    """
    counts = {attr: 0 for attr in _SIZE_ATTR.values()}
    for ssr in ssr_records:
        attr = _SIZE_ATTR.get(ssr.motif_size)
        if attr:
            counts[attr] += 1
    return counts


def compute_dominant_motif(ssr_records) -> str | None:
    """Find the canonical motif with the highest count.

    Returns None if no SSRs exist.
    """
    counter: Counter[str] = Counter()
    for ssr in ssr_records:
        counter[ssr.motif_canonical] += 1
    if not counter:
        return None
    return counter.most_common(1)[0][0]


# ---------------------------------------------------------------------------
# RA / RD formulas
# ---------------------------------------------------------------------------

def compute_ra(ssr_count: int, genome_length: int | None) -> float | None:
    """Relative Abundance = SSR count / genome size in kb.

    Returns None if genome_length is missing or zero.
    """
    if not genome_length or genome_length <= 0:
        return None
    return ssr_count / (genome_length / 1000.0)


def compute_rd(ssr_bp_total: int, genome_length: int | None) -> float | None:
    """Relative Density = total SSR bp / genome size in Mb.

    Returns None if genome_length is missing or zero.
    """
    if not genome_length or genome_length <= 0:
        return None
    return ssr_bp_total / (genome_length / 1_000_000.0)


# ---------------------------------------------------------------------------
# Per-accession metric computation
# ---------------------------------------------------------------------------

def compute_accession_metrics(
    accession: str,
    ssr_records,
    genome_length: int | None = None,
    gc_content: float | None = None,
) -> AccessionMetricsResult:
    """Compute all metrics for a single accession from its SSR records.

    Args:
        accession: Accession identifier.
        ssr_records: Sequence of SSRRecord ORM objects.
        genome_length: Genome length from Accession or SequenceRecord.
        gc_content: GC content from Accession.

    Returns:
        AccessionMetricsResult with all computed metrics.
    """
    ssr_list = list(ssr_records)

    ssr_count_total = len(ssr_list)
    ssr_bp_total = sum(ssr.repeat_length_bp for ssr in ssr_list)

    size_counts = compute_motif_size_counts(ssr_list)
    dominant = compute_dominant_motif(ssr_list)
    ra = compute_ra(ssr_count_total, genome_length)
    rd = compute_rd(ssr_bp_total, genome_length)

    return AccessionMetricsResult(
        accession=accession,
        ssr_count_total=ssr_count_total,
        ssr_bp_total=ssr_bp_total,
        ra=ra,
        rd=rd,
        dominant_motif=dominant,
        genome_length=genome_length,
        gc_content=gc_content,
        **size_counts,
    )


def compute_metrics_for_accession(session, accession: str) -> AccessionMetricsResult:
    """Load SSR records and accession data from DB and compute metrics.

    Args:
        session: SQLAlchemy session.
        accession: Accession identifier.

    Returns:
        AccessionMetricsResult.
    """
    from gwico_ssr.db.repository import get_accession, get_ssr_records_for_accession

    acc_obj = get_accession(session, accession)
    ssr_records = get_ssr_records_for_accession(session, accession)

    genome_length = None
    gc_content = None
    if acc_obj:
        genome_length = acc_obj.genome_length
        gc_content = acc_obj.gc_content
        # Fall back to sequence_record.sequence_length if genome_length not set
        if genome_length is None and acc_obj.sequence_record:
            genome_length = acc_obj.sequence_record.sequence_length

    return compute_accession_metrics(accession, ssr_records, genome_length, gc_content)
