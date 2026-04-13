"""Repository helpers for GWICO-SSR database operations.

Provides insert, upsert, lookup, and bulk operations for core entities.
"""

from __future__ import annotations

import datetime
import platform
from typing import Optional, Sequence

from sqlalchemy import select
from sqlalchemy.orm import Session

from gwico_ssr import __version__
from gwico_ssr.models.schema import (
    Accession,
    AccessionMetrics,
    Dataset,
    FeatureRecord,
    Run,
    SequenceRecord,
    SSRAnnotation,
    SSRRecord,
    StatisticalResult,
)


# ---------------------------------------------------------------------------
# Dataset
# ---------------------------------------------------------------------------

def get_or_create_dataset(
    session: Session,
    name: str,
    source_type: str,
    organism: Optional[str] = None,
    description: Optional[str] = None,
    input_manifest_hash: Optional[str] = None,
) -> Dataset:
    """Return existing dataset by name, or create a new one."""
    stmt = select(Dataset).where(Dataset.name == name)
    dataset = session.execute(stmt).scalar_one_or_none()
    if dataset is not None:
        return dataset
    dataset = Dataset(
        name=name,
        source_type=source_type,
        organism=organism,
        description=description,
        input_manifest_hash=input_manifest_hash,
    )
    session.add(dataset)
    session.flush()
    return dataset


def get_dataset_by_name(session: Session, name: str) -> Optional[Dataset]:
    """Lookup a dataset by name."""
    stmt = select(Dataset).where(Dataset.name == name)
    return session.execute(stmt).scalar_one_or_none()


# ---------------------------------------------------------------------------
# Run
# ---------------------------------------------------------------------------

def create_run(
    session: Session,
    dataset_id: int,
    config_hash: Optional[str] = None,
    stage: Optional[str] = None,
    notes: Optional[str] = None,
) -> Run:
    """Create a new pipeline run record."""
    run = Run(
        dataset_id=dataset_id,
        pipeline_version=__version__,
        config_hash=config_hash,
        started_at=datetime.datetime.now(datetime.UTC),
        status="started",
        stage=stage,
        hostname=platform.node(),
        notes=notes,
    )
    session.add(run)
    session.flush()
    return run


def update_run_status(
    session: Session,
    run_id: int,
    status: str,
    stage: Optional[str] = None,
    notes: Optional[str] = None,
) -> Optional[Run]:
    """Update run status and optionally stage/notes."""
    run = session.get(Run, run_id)
    if run is None:
        return None
    run.status = status
    if stage is not None:
        run.stage = stage
    if notes is not None:
        run.notes = notes
    if status in ("completed", "failed"):
        run.finished_at = datetime.datetime.now(datetime.UTC)
    session.flush()
    return run


def get_run(session: Session, run_id: int) -> Optional[Run]:
    """Lookup a run by ID."""
    return session.get(Run, run_id)


# ---------------------------------------------------------------------------
# Accession
# ---------------------------------------------------------------------------

def upsert_accession(session: Session, **kwargs: object) -> Accession:
    """Insert or update an accession record.

    If the accession already exists, fields with non-None values in kwargs
    are updated. source_priority is compared: higher priority wins.
    """
    acc_id = kwargs.get("accession")
    existing = session.get(Accession, acc_id)
    if existing is not None:
        new_priority = kwargs.get("source_priority", 0) or 0
        old_priority = existing.source_priority or 0
        if new_priority >= old_priority:
            for key, value in kwargs.items():
                if key != "accession" and value is not None:
                    setattr(existing, key, value)
        session.flush()
        return existing
    acc = Accession(**kwargs)  # type: ignore[arg-type]
    session.add(acc)
    session.flush()
    return acc


def bulk_upsert_accessions(session: Session, records: Sequence[dict]) -> int:
    """Upsert multiple accession records. Returns count of processed records."""
    count = 0
    for rec in records:
        upsert_accession(session, **rec)
        count += 1
    return count


def get_accession(session: Session, accession: str) -> Optional[Accession]:
    """Lookup an accession by ID."""
    return session.get(Accession, accession)


def list_accessions(
    session: Session,
    dataset_id: Optional[int] = None,
    country: Optional[str] = None,
    limit: Optional[int] = None,
) -> list[Accession]:
    """List accessions with optional filters."""
    stmt = select(Accession)
    if dataset_id is not None:
        stmt = stmt.where(Accession.dataset_id == dataset_id)
    if country is not None:
        stmt = stmt.where(Accession.country == country)
    if limit is not None:
        stmt = stmt.limit(limit)
    return list(session.execute(stmt).scalars().all())


# ---------------------------------------------------------------------------
# SequenceRecord
# ---------------------------------------------------------------------------

def upsert_sequence_record(session: Session, **kwargs: object) -> SequenceRecord:
    """Insert or update a sequence record keyed by accession."""
    acc_id = kwargs.get("accession")
    existing = session.get(SequenceRecord, acc_id)
    if existing is not None:
        for key, value in kwargs.items():
            if key != "accession" and value is not None:
                setattr(existing, key, value)
        session.flush()
        return existing
    rec = SequenceRecord(**kwargs)  # type: ignore[arg-type]
    session.add(rec)
    session.flush()
    return rec


# ---------------------------------------------------------------------------
# FeatureRecord
# ---------------------------------------------------------------------------

def insert_features(session: Session, features: Sequence[dict]) -> int:
    """Bulk insert feature records. Returns count."""
    objects = [FeatureRecord(**f) for f in features]
    session.add_all(objects)
    session.flush()
    return len(objects)


def get_features_for_accession(session: Session, accession: str) -> list[FeatureRecord]:
    """Get all feature records for an accession."""
    stmt = select(FeatureRecord).where(FeatureRecord.accession == accession)
    return list(session.execute(stmt).scalars().all())


# ---------------------------------------------------------------------------
# SSRRecord
# ---------------------------------------------------------------------------

def insert_ssr_records(session: Session, records: Sequence[dict]) -> int:
    """Bulk insert SSR records. Returns count."""
    objects = [SSRRecord(**r) for r in records]
    session.add_all(objects)
    session.flush()
    return len(objects)


def get_ssr_records_for_accession(session: Session, accession: str) -> list[SSRRecord]:
    """Get all SSR records for an accession."""
    stmt = select(SSRRecord).where(SSRRecord.accession == accession)
    return list(session.execute(stmt).scalars().all())


# ---------------------------------------------------------------------------
# SSRAnnotation
# ---------------------------------------------------------------------------

def insert_ssr_annotations(session: Session, annotations: Sequence[dict]) -> int:
    """Bulk insert SSR annotations. Returns count."""
    objects = [SSRAnnotation(**a) for a in annotations]
    session.add_all(objects)
    session.flush()
    return len(objects)


# ---------------------------------------------------------------------------
# AccessionMetrics
# ---------------------------------------------------------------------------

def upsert_accession_metrics(session: Session, **kwargs: object) -> AccessionMetrics:
    """Insert or replace accession metrics for a given (accession, run_id) pair."""
    acc = kwargs.get("accession")
    rid = kwargs.get("run_id")
    stmt = select(AccessionMetrics).where(
        AccessionMetrics.accession == acc,
        AccessionMetrics.run_id == rid,
    )
    existing = session.execute(stmt).scalar_one_or_none()
    if existing is not None:
        for key, value in kwargs.items():
            if key not in ("accession", "run_id") and value is not None:
                setattr(existing, key, value)
        session.flush()
        return existing
    m = AccessionMetrics(**kwargs)  # type: ignore[arg-type]
    session.add(m)
    session.flush()
    return m


def get_accession_metrics(
    session: Session, accession: str, run_id: int | None = None
) -> list[AccessionMetrics]:
    """Get metrics for an accession, optionally filtered by run_id."""
    stmt = select(AccessionMetrics).where(AccessionMetrics.accession == accession)
    if run_id is not None:
        stmt = stmt.where(AccessionMetrics.run_id == run_id)
    return list(session.execute(stmt).scalars().all())


def has_accession_metrics(session: Session, accession: str, run_id: int) -> bool:
    """Check if metrics already exist for an accession + run_id pair."""
    from sqlalchemy import func as sqla_func
    count = session.execute(
        select(sqla_func.count()).where(
            AccessionMetrics.accession == accession,
            AccessionMetrics.run_id == run_id,
        )
    ).scalar()
    return bool(count and count > 0)


# ---------------------------------------------------------------------------
# StatisticalResult
# ---------------------------------------------------------------------------

def insert_statistical_result(session: Session, **kwargs: object) -> StatisticalResult:
    """Insert a statistical result record."""
    result = StatisticalResult(**kwargs)  # type: ignore[arg-type]
    session.add(result)
    session.flush()
    return result
