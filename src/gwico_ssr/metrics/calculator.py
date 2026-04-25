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
    """Computed metrics for a single accession.
    
    Chunk 8: Includes repeat_class breakdowns for batch-aware downstream analysis.
    """

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
    # Chunk 8: Repeat class breakdowns
    perfect_count: int = 0
    imperfect_count: int = 0
    compound_component_count: int = 0
    perfect_bp_total: int = 0
    imperfect_bp_total: int = 0
    compound_component_bp_total: int = 0

    def to_dict(self, run_id: int) -> dict:
        """Convert to dict suitable for upsert_accession_metrics.
        
        Chunk 8: Includes repeat_class breakdown fields.
        """
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
            "perfect_count": self.perfect_count,
            "imperfect_count": self.imperfect_count,
            "compound_component_count": self.compound_component_count,
            "perfect_bp_total": self.perfect_bp_total,
            "imperfect_bp_total": self.imperfect_bp_total,
            "compound_component_bp_total": self.compound_component_bp_total,
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
# Chunk 8: Repeat class breakdown helpers
# ---------------------------------------------------------------------------

def compute_repeat_class_breakdown(ssr_records) -> dict[str, int]:
    """Count SSRs by repeat_class.
    
    Returns dict with keys: perfect_count, imperfect_count, compound_component_count.
    """
    counts = {
        "perfect_count": 0,
        "imperfect_count": 0,
        "compound_component_count": 0,
    }
    for ssr in ssr_records:
        repeat_class = getattr(ssr, "repeat_class", "perfect")  # backward compat
        if repeat_class == "perfect":
            counts["perfect_count"] += 1
        elif repeat_class == "imperfect":
            counts["imperfect_count"] += 1
        elif repeat_class == "compound_component":
            counts["compound_component_count"] += 1
    return counts


def compute_repeat_class_bp_breakdown(ssr_records) -> dict[str, int]:
    """Sum repeat_length_bp by repeat_class.
    
    Returns dict with keys: perfect_bp_total, imperfect_bp_total, compound_component_bp_total.
    """
    bp_totals = {
        "perfect_bp_total": 0,
        "imperfect_bp_total": 0,
        "compound_component_bp_total": 0,
    }
    for ssr in ssr_records:
        repeat_class = getattr(ssr, "repeat_class", "perfect")  # backward compat
        bp = ssr.repeat_length_bp
        if repeat_class == "perfect":
            bp_totals["perfect_bp_total"] += bp
        elif repeat_class == "imperfect":
            bp_totals["imperfect_bp_total"] += bp
        elif repeat_class == "compound_component":
            bp_totals["compound_component_bp_total"] += bp
    return bp_totals


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
    
    Chunk 8: Includes repeat_class breakdown computation.

    Args:
        accession: Accession identifier.
        ssr_records: Sequence of SSRRecord ORM objects.
        genome_length: Genome length from Accession or SequenceRecord.
        gc_content: GC content from Accession.

    Returns:
        AccessionMetricsResult with all computed metrics, including class breakdowns.
    """
    ssr_list = list(ssr_records)

    ssr_count_total = len(ssr_list)
    ssr_bp_total = sum(ssr.repeat_length_bp for ssr in ssr_list)

    size_counts = compute_motif_size_counts(ssr_list)
    dominant = compute_dominant_motif(ssr_list)
    ra = compute_ra(ssr_count_total, genome_length)
    rd = compute_rd(ssr_bp_total, genome_length)
    
    # Chunk 8: Compute repeat class breakdowns
    class_counts = compute_repeat_class_breakdown(ssr_list)
    class_bp_totals = compute_repeat_class_bp_breakdown(ssr_list)

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
        **class_counts,
        **class_bp_totals,
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
