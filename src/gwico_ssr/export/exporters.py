"""Export builders for GWICO-SSR.

Generates CSV, JSON, BED, and GFF3 exports from database-backed pipeline
state, plus run manifests for provenance tracking.
"""

from __future__ import annotations

import csv
import hashlib
import json
import platform
from datetime import datetime, timezone
from io import StringIO
from pathlib import Path
from typing import Any, Sequence

from sqlalchemy import select
from sqlalchemy.orm import Session

from gwico_ssr import __version__
from gwico_ssr.models.schema import (
    Accession,
    AccessionMetrics,
    Dataset,
    Run,
    SSRAnnotation,
    SSRRecord,
    StatisticalResult,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _ensure_dir(path: str | Path) -> Path:
    p = Path(path)
    p.mkdir(parents=True, exist_ok=True)
    return p


def _file_sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(8192), b""):
            h.update(chunk)
    return h.hexdigest()


# ---------------------------------------------------------------------------
# Filters
# ---------------------------------------------------------------------------

def _apply_ssr_filters(
    stmt,
    *,
    dataset_id: int | None = None,
    country: str | None = None,
    gene: str | None = None,
    motif: str | None = None,
    run_id: int | None = None,
):
    """Apply optional filters to an SSR query. Returns modified statement."""
    _joined_accession = False
    if run_id is not None:
        stmt = stmt.where(SSRRecord.run_id == run_id)
    if dataset_id is not None:
        stmt = stmt.join(Accession, SSRRecord.accession == Accession.accession).where(
            Accession.dataset_id == dataset_id
        )
        _joined_accession = True
    if country is not None:
        if not _joined_accession:
            stmt = stmt.join(Accession, SSRRecord.accession == Accession.accession)
            _joined_accession = True
        stmt = stmt.where(Accession.country == country)
    if motif is not None:
        stmt = stmt.where(SSRRecord.motif_canonical == motif)
    return stmt


# ---------------------------------------------------------------------------
# CSV Export
# ---------------------------------------------------------------------------

def export_ssrs_csv(
    session: Session,
    output_dir: str | Path,
    *,
    dataset_id: int | None = None,
    run_id: int | None = None,
    country: str | None = None,
    motif: str | None = None,
) -> Path:
    """Export SSR records to CSV.
    
    Chunk 8: Includes repeat_class for downstream batch-aware processing.
    """
    out = _ensure_dir(output_dir)
    path = out / "ssr_records.csv"

    stmt = select(SSRRecord)
    stmt = _apply_ssr_filters(stmt, dataset_id=dataset_id, run_id=run_id,
                               country=country, motif=motif)

    rows = session.execute(stmt).scalars().all()

    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow([
            "ssr_id", "accession", "start", "end", "motif_raw",
            "motif_canonical", "repeat_units", "motif_size",
            "repeat_length_bp", "strand", "actual_repeat", "detector_version",
            "repeat_class", "is_perfect", "imperfection_pct",
        ])
        for r in rows:
            writer.writerow([
                r.ssr_id, r.accession, r.start, r.end, r.motif_raw,
                r.motif_canonical, r.repeat_units, r.motif_size,
                r.repeat_length_bp, r.strand, r.actual_repeat,
                r.detector_version, getattr(r, "repeat_class", "perfect"),
                getattr(r, "is_perfect", True),
                getattr(r, "imperfection_pct", None),
            ])
    return path


def export_metrics_csv(
    session: Session,
    output_dir: str | Path,
    *,
    run_id: int | None = None,
) -> Path:
    """Export per-accession metrics to CSV.
    
    Chunk 8: Includes repeat_class breakdown fields for batch-aware analysis.
    """
    out = _ensure_dir(output_dir)
    path = out / "accession_metrics.csv"

    stmt = select(AccessionMetrics)
    if run_id is not None:
        stmt = stmt.where(AccessionMetrics.run_id == run_id)

    rows = session.execute(stmt).scalars().all()

    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow([
            "accession", "run_id", "ssr_count_total", "ssr_bp_total",
            "ra", "rd", "mono_count", "di_count", "tri_count",
            "tetra_count", "penta_count", "hexa_count", "dominant_motif",
            "perfect_count", "imperfect_count", "compound_component_count",
            "perfect_bp_total", "imperfect_bp_total", "compound_component_bp_total",
        ])
        for m in rows:
            writer.writerow([
                m.accession, m.run_id, m.ssr_count_total, m.ssr_bp_total,
                m.ra, m.rd, m.mono_count, m.di_count, m.tri_count,
                m.tetra_count, m.penta_count, m.hexa_count, m.dominant_motif,
                getattr(m, "perfect_count", 0),
                getattr(m, "imperfect_count", 0),
                getattr(m, "compound_component_count", 0),
                getattr(m, "perfect_bp_total", 0),
                getattr(m, "imperfect_bp_total", 0),
                getattr(m, "compound_component_bp_total", 0),
            ])
    return path


def export_stats_csv(
    session: Session,
    output_dir: str | Path,
    *,
    run_id: int | None = None,
) -> Path:
    """Export statistical results to CSV."""
    out = _ensure_dir(output_dir)
    path = out / "statistical_results.csv"

    stmt = select(StatisticalResult)
    if run_id is not None:
        stmt = stmt.where(StatisticalResult.run_id == run_id)

    rows = session.execute(stmt).scalars().all()

    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow([
            "result_id", "run_id", "analysis_name", "grouping", "metric",
            "test_name", "statistic", "p_value", "p_value_corrected",
            "effect_size", "ci_lower", "ci_upper", "n", "metadata_json",
        ])
        for s in rows:
            writer.writerow([
                s.result_id, s.run_id, s.analysis_name, s.grouping,
                s.metric, s.test_name, s.statistic, s.p_value,
                s.p_value_corrected, s.effect_size, s.ci_lower, s.ci_upper,
                s.n, s.metadata_json,
            ])
    return path


# ---------------------------------------------------------------------------
# JSON Export
# ---------------------------------------------------------------------------

def export_ssrs_json(
    session: Session,
    output_dir: str | Path,
    *,
    dataset_id: int | None = None,
    run_id: int | None = None,
    country: str | None = None,
    motif: str | None = None,
) -> Path:
    """Export SSR records to JSON.
    
    Chunk 8: Includes repeat_class and imperfection fields.
    """
    out = _ensure_dir(output_dir)
    path = out / "ssr_records.json"

    stmt = select(SSRRecord)
    stmt = _apply_ssr_filters(stmt, dataset_id=dataset_id, run_id=run_id,
                               country=country, motif=motif)

    rows = session.execute(stmt).scalars().all()

    records = [
        {
            "ssr_id": r.ssr_id,
            "accession": r.accession,
            "start": r.start,
            "end": r.end,
            "motif_raw": r.motif_raw,
            "motif_canonical": r.motif_canonical,
            "repeat_units": r.repeat_units,
            "motif_size": r.motif_size,
            "repeat_length_bp": r.repeat_length_bp,
            "strand": r.strand,
            "actual_repeat": r.actual_repeat,
            "detector_version": r.detector_version,
            "repeat_class": getattr(r, "repeat_class", "perfect"),
            "is_perfect": getattr(r, "is_perfect", True),
            "imperfection_pct": getattr(r, "imperfection_pct", None),
        }
        for r in rows
    ]

    with open(path, "w", encoding="utf-8") as f:
        json.dump(records, f, indent=2, default=str)
    return path


def export_stats_json(
    session: Session,
    output_dir: str | Path,
    *,
    run_id: int | None = None,
) -> Path:
    """Export statistical results to JSON."""
    out = _ensure_dir(output_dir)
    path = out / "statistical_results.json"

    stmt = select(StatisticalResult)
    if run_id is not None:
        stmt = stmt.where(StatisticalResult.run_id == run_id)

    rows = session.execute(stmt).scalars().all()

    records = [
        {
            "result_id": s.result_id,
            "run_id": s.run_id,
            "analysis_name": s.analysis_name,
            "grouping": s.grouping,
            "metric": s.metric,
            "test_name": s.test_name,
            "statistic": s.statistic,
            "p_value": s.p_value,
            "p_value_corrected": s.p_value_corrected,
            "effect_size": s.effect_size,
            "ci_lower": s.ci_lower,
            "ci_upper": s.ci_upper,
            "n": s.n,
            "metadata": json.loads(s.metadata_json) if s.metadata_json else None,
        }
        for s in rows
    ]

    with open(path, "w", encoding="utf-8") as f:
        json.dump(records, f, indent=2, default=str)
    return path


# ---------------------------------------------------------------------------
# BED Export (0-based, half-open)
# ---------------------------------------------------------------------------

def export_ssrs_bed(
    session: Session,
    output_dir: str | Path,
    *,
    dataset_id: int | None = None,
    run_id: int | None = None,
    country: str | None = None,
    motif: str | None = None,
) -> Path:
    """Export SSR records in BED6 format.

    BED columns: chrom, chromStart, chromEnd, name, score, strand.
    Coordinates are 0-based, half-open (matching SSRRecord storage).
    """
    out = _ensure_dir(output_dir)
    path = out / "ssr_records.bed"

    stmt = select(SSRRecord).order_by(SSRRecord.accession, SSRRecord.start)
    stmt = _apply_ssr_filters(stmt, dataset_id=dataset_id, run_id=run_id,
                               country=country, motif=motif)

    rows = session.execute(stmt).scalars().all()

    with open(path, "w", encoding="utf-8") as f:
        for r in rows:
            # name = canonical_motif x repeat_units
            name = f"{r.motif_canonical}x{r.repeat_units}"
            score = r.repeat_length_bp  # Use repeat length as score
            strand = r.strand or "+"
            f.write(f"{r.accession}\t{r.start}\t{r.end}\t{name}\t{score}\t{strand}\n")
    return path


# ---------------------------------------------------------------------------
# GFF3 Export
# ---------------------------------------------------------------------------

def export_ssrs_gff3(
    session: Session,
    output_dir: str | Path,
    *,
    dataset_id: int | None = None,
    run_id: int | None = None,
    country: str | None = None,
    motif: str | None = None,
) -> Path:
    """Export SSR records in GFF3 format.

    GFF3 uses 1-based, fully-closed coordinates.
    SSRRecord stores 0-based, half-open → convert: start+1, end unchanged.
    
    Chunk 8: Includes repeat_class in attributes.
    """
    out = _ensure_dir(output_dir)
    path = out / "ssr_records.gff3"

    stmt = select(SSRRecord).order_by(SSRRecord.accession, SSRRecord.start)
    stmt = _apply_ssr_filters(stmt, dataset_id=dataset_id, run_id=run_id,
                               country=country, motif=motif)

    rows = session.execute(stmt).scalars().all()

    with open(path, "w", encoding="utf-8") as f:
        f.write("##gff-version 3\n")
        for r in rows:
            # 0-based half-open → 1-based closed
            gff_start = r.start + 1
            gff_end = r.end  # half-open end == closed end in 1-based
            strand = r.strand or "+"
            attrs = (
                f"ID=ssr_{r.ssr_id};"
                f"Name={r.motif_canonical}x{r.repeat_units};"
                f"motif_raw={r.motif_raw};"
                f"motif_canonical={r.motif_canonical};"
                f"motif_size={r.motif_size};"
                f"repeat_units={r.repeat_units};"
                f"repeat_length_bp={r.repeat_length_bp};"
                f"repeat_class={getattr(r, 'repeat_class', 'perfect')}"
            )
            if r.actual_repeat:
                attrs += f";actual_repeat={r.actual_repeat}"
            if hasattr(r, "imperfection_pct") and r.imperfection_pct is not None:
                attrs += f";imperfection_pct={r.imperfection_pct}"
            f.write(
                f"{r.accession}\tgwico-ssr\tmicrosatellite\t{gff_start}\t{gff_end}"
                f"\t{r.repeat_length_bp}\t{strand}\t.\t{attrs}\n"
            )
    return path


# ---------------------------------------------------------------------------
# Publication Tables
# ---------------------------------------------------------------------------

def export_publication_tables(
    session: Session,
    output_dir: str | Path,
    *,
    run_id: int,
    dataset_id: int,
) -> list[Path]:
    """Export publication-ready summary tables as CSV files.

    Returns list of created file paths.
    """
    from gwico_ssr.metrics.aggregator import (
        aggregate_by_country,
        aggregate_by_gene,
        aggregate_motif_frequencies,
        get_dataset_summary,
    )

    out = _ensure_dir(output_dir)
    paths: list[Path] = []

    # 1. Dataset summary
    summary = get_dataset_summary(session, run_id)
    summary_path = out / "dataset_summary.json"
    with open(summary_path, "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2, default=str)
    paths.append(summary_path)

    # 2. Country summary table
    country_data = aggregate_by_country(session, run_id)
    country_path = out / "country_summary.csv"
    with open(country_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["country", "accession_count", "total_ssrs",
                         "total_ssr_bp", "mean_ra", "mean_rd"])
        for c in country_data:
            writer.writerow([c.group_value, c.count, c.total_ssrs,
                             c.total_ssr_bp, c.mean_ra, c.mean_rd])
    paths.append(country_path)

    # 3. Motif frequency table
    motif_data = aggregate_motif_frequencies(session, dataset_id)
    motif_path = out / "motif_frequencies.csv"
    with open(motif_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["motif_canonical", "motif_size", "total_count", "total_bp"])
        for m in motif_data:
            writer.writerow([m.motif_canonical, m.motif_size,
                             m.total_count, m.total_bp])
    paths.append(motif_path)

    # 4. Gene distribution table
    gene_data = aggregate_by_gene(session, dataset_id)
    gene_path = out / "gene_distribution.csv"
    with open(gene_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["gene_name", "ssr_count", "accession_count"])
        for g in gene_data:
            writer.writerow([g.gene_name, g.ssr_count, g.accession_count])
    paths.append(gene_path)

    return paths


# ---------------------------------------------------------------------------
# Run Manifest
# ---------------------------------------------------------------------------

def generate_run_manifest(
    session: Session,
    output_dir: str | Path,
    *,
    run_id: int,
    config_snapshot: dict[str, Any] | None = None,
    output_files: list[Path] | None = None,
) -> Path:
    """Generate a run manifest JSON with provenance metadata.

    Ties every output artifact to the run that produced it.
    """
    out = _ensure_dir(output_dir)
    path = out / "run_manifest.json"

    run = session.execute(
        select(Run).where(Run.run_id == run_id)
    ).scalar_one_or_none()

    if run is None:
        raise ValueError(f"Run {run_id} not found")

    dataset = session.execute(
        select(Dataset).where(Dataset.dataset_id == run.dataset_id)
    ).scalar_one_or_none()

    # Count records
    from sqlalchemy import func as sa_func

    ssr_count = session.execute(
        select(sa_func.count(SSRRecord.ssr_id))
        .where(SSRRecord.run_id == run_id)
    ).scalar() or 0

    annotation_count = session.execute(
        select(sa_func.count(SSRAnnotation.annotation_id))
        .join(SSRRecord, SSRAnnotation.ssr_id == SSRRecord.ssr_id)
        .where(SSRRecord.run_id == run_id)
    ).scalar() or 0

    metrics_count = session.execute(
        select(sa_func.count(AccessionMetrics.id))
        .where(AccessionMetrics.run_id == run_id)
    ).scalar() or 0

    # Build file hashes for output files
    file_manifest = []
    if output_files:
        for fp in output_files:
            if fp.exists():
                file_manifest.append({
                    "path": str(fp),
                    "size_bytes": fp.stat().st_size,
                    "sha256": _file_sha256(fp),
                })

    manifest = {
        "gwico_ssr_version": __version__,
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "hostname": platform.node(),
        "python_version": platform.python_version(),
        "run": {
            "run_id": run.run_id,
            "dataset_id": run.dataset_id,
            "dataset_name": dataset.name if dataset else None,
            "stage": run.stage,
            "status": run.status,
            "pipeline_version": run.pipeline_version,
            "config_hash": run.config_hash,
            "started_at": run.started_at.isoformat() if run.started_at else None,
            "finished_at": run.finished_at.isoformat() if run.finished_at else None,
        },
        "counts": {
            "ssr_records": ssr_count,
            "annotations": annotation_count,
            "accession_metrics": metrics_count,
        },
        "config_snapshot": config_snapshot,
        "output_files": file_manifest,
    }

    with open(path, "w", encoding="utf-8") as f:
        json.dump(manifest, f, indent=2, default=str)
    return path


# ---------------------------------------------------------------------------
# Orchestrator: export all formats for a run
# ---------------------------------------------------------------------------

def export_all(
    session: Session,
    output_dir: str | Path,
    *,
    run_id: int,
    dataset_id: int,
    formats: set[str] | None = None,
    country: str | None = None,
    motif: str | None = None,
    config_snapshot: dict[str, Any] | None = None,
) -> list[Path]:
    """Export all requested formats and generate a run manifest.

    *formats*: set of format names. Default: all.
    Supported: ``csv``, ``json``, ``bed``, ``gff3``, ``tables``.
    """
    all_formats = {"csv", "json", "bed", "gff3", "tables"}
    if formats is None:
        formats = all_formats
    else:
        formats = formats & all_formats

    out = _ensure_dir(output_dir)
    paths: list[Path] = []

    filter_kwargs = dict(
        dataset_id=dataset_id, run_id=run_id,
        country=country, motif=motif,
    )

    if "csv" in formats:
        paths.append(export_ssrs_csv(session, out, **filter_kwargs))
        paths.append(export_metrics_csv(session, out, run_id=run_id))
        paths.append(export_stats_csv(session, out, run_id=run_id))

    if "json" in formats:
        paths.append(export_ssrs_json(session, out, **filter_kwargs))
        paths.append(export_stats_json(session, out, run_id=run_id))

    if "bed" in formats:
        paths.append(export_ssrs_bed(session, out, **filter_kwargs))

    if "gff3" in formats:
        paths.append(export_ssrs_gff3(session, out, **filter_kwargs))

    if "tables" in formats:
        paths.extend(export_publication_tables(
            session, out, run_id=run_id, dataset_id=dataset_id,
        ))

    # Generate manifest last, so it can include file hashes
    manifest_path = generate_run_manifest(
        session, out, run_id=run_id,
        config_snapshot=config_snapshot,
        output_files=paths,
    )
    paths.append(manifest_path)

    return paths
