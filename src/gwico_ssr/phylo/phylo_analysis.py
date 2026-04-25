"""Phylogenetic tree-aware analysis for GWICO-SSR.

Implements tree ingestion, accession-to-tip mapping, and phylo-aware metrics
with full audit trail and reproducibility guarantees.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

from sqlalchemy import func, select, case
from sqlalchemy.orm import Session

from gwico_ssr.models.schema import (
    Accession,
    AccessionMetrics,
    SSRRecord,
    TreeMapping as TreeMappingTable,
    TreeMetrics as TreeMetricsTable,
    TreeSource as TreeSourceTable,
)
from gwico_ssr.phylo.tree_parser import Tree, parse_newick_file, validate_tree_structure

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Dataclasses
# ---------------------------------------------------------------------------

@dataclass
class TreeSource:
    """Metadata about a phylogenetic tree source."""

    source_name: str  # user-supplied, nextstrain, etc.
    version: str  # v1, 2024-01, etc.
    file_name: Optional[str] = None
    file_hash: Optional[str] = None
    num_tips: int = 0
    release_date: Optional[str] = None
    description: Optional[str] = None


@dataclass
class TreeMapping:
    """Single accession-to-tree-tip mapping."""

    accession: str
    tip_label: str
    source: str  # source_name from TreeSource
    is_confident: bool = True
    conflict_notes: Optional[str] = None


@dataclass
class TreeCladeMetrics:
    """Aggregated metrics for a single clade/tip in tree."""

    clade_name: str
    tip_count: int = 0
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
class TreeMappingReport:
    """Complete tree mapping report with audit info."""

    source: TreeSource
    total_tips_in_tree: int
    total_accessions_in_database: int
    successfully_mapped: int
    unmapped: int
    conflicted: int = 0
    conflicts: list[str] = field(default_factory=list)
    analysis_notes: list[str] = field(default_factory=list)


# ---------------------------------------------------------------------------
# Tree ingestion
# ---------------------------------------------------------------------------

def ingest_newick_tree(
    session: Session,
    tree_file: str | Path,
    source_name: str = "user-supplied",
    version: str = "v1",
    release_date: Optional[str] = None,
) -> tuple[Tree, TreeMappingReport]:
    """Ingest phylogenetic tree from Newick file.

    Args:
        session: SQLAlchemy session
        tree_file: Path to Newick format file
        source_name: Source identifier
        version: Version string
        release_date: Optional ISO 8601 release date

    Returns:
        Tuple of (parsed Tree object, TreeMappingReport)

    Raises:
        FileNotFoundError: If tree file not found
        ValueError: If Newick format invalid
    """
    tree_file = Path(tree_file)

    # Parse tree file
    tree, file_hash = parse_newick_file(tree_file)

    # Validate tree structure
    validation_msgs = validate_tree_structure(tree)
    if validation_msgs:
        logger.warning(f"Tree validation: {validation_msgs}")

    # Create or get tree source
    source_obj = session.query(TreeSourceTable).filter(
        TreeSourceTable.source_name == source_name,
        TreeSourceTable.version == version,
    ).first()

    if not source_obj:
        source_obj = TreeSourceTable(
            source_name=source_name,
            version=version,
            file_name=tree_file.name,
            file_hash=file_hash,
            num_tips=tree.num_tips,
            release_date=release_date,
        )
        session.add(source_obj)
        session.flush()

    report = TreeMappingReport(
        source=TreeSource(
            source_name=source_name,
            version=version,
            file_name=tree_file.name,
            file_hash=file_hash,
            num_tips=tree.num_tips,
            release_date=release_date,
        ),
        total_tips_in_tree=tree.num_tips,
        total_accessions_in_database=0,
        successfully_mapped=0,
        unmapped=0,
    )

    logger.info(f"Ingested tree: {tree.num_tips} tips from {tree_file.name}")

    return tree, report


def map_accessions_to_tips(
    session: Session,
    tree: Tree,
    source_id: int,
) -> TreeMappingReport:
    """Map accessions in database to tree tips.

    Args:
        session: SQLAlchemy session
        tree: Parsed tree object
        source_id: TreeSource ID

    Returns:
        TreeMappingReport with mapping results
    """
    source_obj = session.query(TreeSourceTable).filter(
        TreeSourceTable.source_id == source_id
    ).first()

    if not source_obj:
        raise ValueError(f"TreeSource {source_id} not found")

    # Get all accessions
    accessions = session.query(Accession.accession).all()
    accession_list = [a[0] for a in accessions]

    # Get tree tips
    tree_tips = set(tree.get_all_tips())

    report = TreeMappingReport(
        source=TreeSource(
            source_name=source_obj.source_name,
            version=source_obj.version,
            file_name=source_obj.file_name,
            file_hash=source_obj.file_hash,
            num_tips=source_obj.num_tips,
            release_date=source_obj.release_date,
        ),
        total_tips_in_tree=tree.num_tips,
        total_accessions_in_database=len(accession_list),
        successfully_mapped=0,
        unmapped=0,
    )

    # Create mappings where accession matches tip label
    for accession in accession_list:
        if accession in tree_tips:
            # Check for existing mapping
            existing = session.query(TreeMappingTable).filter(
                TreeMappingTable.accession == accession,
                TreeMappingTable.source_id == source_id,
            ).first()

            if not existing:
                mapping = TreeMappingTable(
                    accession=accession,
                    source_id=source_id,
                    tip_label=accession,
                    is_confident=True,
                )
                session.add(mapping)
                report.successfully_mapped += 1
        else:
            report.unmapped += 1
            report.analysis_notes.append(f"Accession {accession} not found in tree tips")

    session.commit()
    logger.info(f"Mapped {report.successfully_mapped} accessions to tree tips")

    return report


def audit_tree_mapping(
    session: Session,
    source_id: int,
) -> TreeMappingReport:
    """Generate audit report for tree mapping.

    Args:
        session: SQLAlchemy session
        source_id: TreeSource ID

    Returns:
        TreeMappingReport with audit details
    """
    source_obj = session.query(TreeSourceTable).filter(
        TreeSourceTable.source_id == source_id
    ).first()

    if not source_obj:
        raise ValueError(f"TreeSource {source_id} not found")

    total_mapped = session.query(func.count(TreeMappingTable.mapping_id)).filter(
        TreeMappingTable.source_id == source_id
    ).scalar()

    confident_count = session.query(func.count(TreeMappingTable.mapping_id)).filter(
        TreeMappingTable.source_id == source_id,
        TreeMappingTable.is_confident == True,
    ).scalar()

    conflicted_count = total_mapped - confident_count if total_mapped else 0

    report = TreeMappingReport(
        source=TreeSource(
            source_name=source_obj.source_name,
            version=source_obj.version,
            file_name=source_obj.file_name,
            file_hash=source_obj.file_hash,
            num_tips=source_obj.num_tips,
            release_date=source_obj.release_date,
        ),
        total_tips_in_tree=source_obj.num_tips,
        total_accessions_in_database=0,
        successfully_mapped=total_mapped or 0,
        unmapped=0,
        conflicted=conflicted_count,
    )

    return report


# ---------------------------------------------------------------------------
# Phylo-aware metrics
# ---------------------------------------------------------------------------

def compute_ssr_metrics_by_clade(
    session: Session,
    source_id: int,
    run_id: Optional[int] = None,
) -> dict[str, TreeCladeMetrics]:
    """Compute SSR metrics aggregated by clade/tip with repeat-class breakdown.

    Args:
        session: SQLAlchemy session
        source_id: TreeSource ID
        run_id: Optional run_id filter

    Returns:
        Dict mapping clade_name → TreeCladeMetrics
    """
    # Get all unique tips in mappings
    tips = session.query(
        func.distinct(TreeMappingTable.tip_label)
    ).filter(
        TreeMappingTable.source_id == source_id
    ).all()

    results = {}

    for tip_tuple in tips:
        tip_label = tip_tuple[0]

        # Get accessions mapped to this tip
        accessions_for_tip = session.query(TreeMappingTable.accession).filter(
            TreeMappingTable.source_id == source_id,
            TreeMappingTable.tip_label == tip_label,
        ).all()

        accession_list = [a[0] for a in accessions_for_tip]
        accession_count = len(accession_list)

        if accession_count == 0:
            continue

        # Query SSR metrics
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

        summary = TreeCladeMetrics(
            clade_name=tip_label,
            tip_count=1,  # Single tip in this case (leaf nodes)
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

        results[tip_label] = summary

    return results


def store_tree_metrics(
    session: Session,
    run_id: int,
    source_id: int,
    metrics_dict: dict[str, TreeCladeMetrics],
) -> int:
    """Store tree metrics to database for persistence and caching.

    Args:
        session: SQLAlchemy session
        run_id: Run ID
        source_id: TreeSource ID
        metrics_dict: Dictionary of clade_name → TreeCladeMetrics

    Returns:
        Number of metric records created
    """
    count = 0
    for clade_name, summary in metrics_dict.items():
        metric = TreeMetricsTable(
            run_id=run_id,
            source_id=source_id,
            clade_name=clade_name,
            tip_count=summary.tip_count,
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
    logger.info(f"Stored {count} tree metrics records for run {run_id}")
    return count
