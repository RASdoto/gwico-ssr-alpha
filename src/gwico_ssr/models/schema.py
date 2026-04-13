"""SQLAlchemy ORM models for GWICO-SSR canonical data model.

Nine core entities: Dataset, Run, Accession, SequenceRecord, FeatureRecord,
SSRRecord, SSRAnnotation, AccessionMetrics, StatisticalResult.
"""

from __future__ import annotations

import datetime
from typing import Optional

from sqlalchemy import (
    Boolean,
    DateTime,
    Float,
    ForeignKey,
    Index,
    Integer,
    String,
    Text,
    UniqueConstraint,
)
from sqlalchemy.orm import (
    DeclarativeBase,
    Mapped,
    mapped_column,
    relationship,
)


class Base(DeclarativeBase):
    """Declarative base for all GWICO-SSR ORM models."""
    pass


class Dataset(Base):
    __tablename__ = "datasets"

    dataset_id: Mapped[int] = mapped_column(Integer, primary_key=True, autoincrement=True)
    name: Mapped[str] = mapped_column(String(255), nullable=False, unique=True)
    source_type: Mapped[str] = mapped_column(String(50), nullable=False)  # csv, fasta, ncbi, mixed
    organism: Mapped[Optional[str]] = mapped_column(String(255), nullable=True)
    created_at: Mapped[datetime.datetime] = mapped_column(
        DateTime, nullable=False, default=lambda: datetime.datetime.now(datetime.UTC)
    )
    description: Mapped[Optional[str]] = mapped_column(Text, nullable=True)
    input_manifest_hash: Mapped[Optional[str]] = mapped_column(String(128), nullable=True)

    runs: Mapped[list[Run]] = relationship("Run", back_populates="dataset")
    accessions: Mapped[list[Accession]] = relationship("Accession", back_populates="dataset")


class Run(Base):
    __tablename__ = "runs"

    run_id: Mapped[int] = mapped_column(Integer, primary_key=True, autoincrement=True)
    dataset_id: Mapped[int] = mapped_column(Integer, ForeignKey("datasets.dataset_id"), nullable=False)
    pipeline_version: Mapped[str] = mapped_column(String(50), nullable=False)
    config_hash: Mapped[Optional[str]] = mapped_column(String(128), nullable=True)
    started_at: Mapped[datetime.datetime] = mapped_column(
        DateTime, nullable=False, default=lambda: datetime.datetime.now(datetime.UTC)
    )
    finished_at: Mapped[Optional[datetime.datetime]] = mapped_column(DateTime, nullable=True)
    status: Mapped[str] = mapped_column(String(50), nullable=False, default="started")  # started, running, completed, failed
    stage: Mapped[Optional[str]] = mapped_column(String(100), nullable=True)
    hostname: Mapped[Optional[str]] = mapped_column(String(255), nullable=True)
    notes: Mapped[Optional[str]] = mapped_column(Text, nullable=True)

    dataset: Mapped[Dataset] = relationship("Dataset", back_populates="runs")
    ssr_records: Mapped[list[SSRRecord]] = relationship("SSRRecord", back_populates="run")
    accession_metrics: Mapped[list[AccessionMetrics]] = relationship("AccessionMetrics", back_populates="run")
    statistical_results: Mapped[list[StatisticalResult]] = relationship("StatisticalResult", back_populates="run")

    __table_args__ = (
        Index("ix_runs_dataset_id", "dataset_id"),
        Index("ix_runs_status", "status"),
    )


class Accession(Base):
    __tablename__ = "accessions"

    accession: Mapped[str] = mapped_column(String(50), primary_key=True)
    dataset_id: Mapped[int] = mapped_column(Integer, ForeignKey("datasets.dataset_id"), nullable=False)
    species: Mapped[Optional[str]] = mapped_column(String(255), nullable=True)
    release_date: Mapped[Optional[str]] = mapped_column(String(50), nullable=True)
    collection_date: Mapped[Optional[str]] = mapped_column(String(50), nullable=True)
    geo_location_raw: Mapped[Optional[str]] = mapped_column(String(500), nullable=True)
    country: Mapped[Optional[str]] = mapped_column(String(100), nullable=True)
    region: Mapped[Optional[str]] = mapped_column(String(255), nullable=True)
    host: Mapped[Optional[str]] = mapped_column(String(255), nullable=True)
    is_complete: Mapped[Optional[bool]] = mapped_column(Boolean, nullable=True)
    genome_length: Mapped[Optional[int]] = mapped_column(Integer, nullable=True)
    gc_content: Mapped[Optional[float]] = mapped_column(Float, nullable=True)
    source_priority: Mapped[Optional[int]] = mapped_column(Integer, nullable=True, default=0)

    dataset: Mapped[Dataset] = relationship("Dataset", back_populates="accessions")
    sequence_record: Mapped[Optional[SequenceRecord]] = relationship(
        "SequenceRecord", back_populates="accession_rel", uselist=False
    )
    features: Mapped[list[FeatureRecord]] = relationship("FeatureRecord", back_populates="accession_rel")
    ssr_records: Mapped[list[SSRRecord]] = relationship("SSRRecord", back_populates="accession_rel")
    ssr_annotations: Mapped[list[SSRAnnotation]] = relationship("SSRAnnotation", back_populates="accession_rel")
    metrics: Mapped[list[AccessionMetrics]] = relationship("AccessionMetrics", back_populates="accession_rel")

    __table_args__ = (
        Index("ix_accessions_dataset_id", "dataset_id"),
        Index("ix_accessions_country", "country"),
        Index("ix_accessions_collection_date", "collection_date"),
        Index("ix_accessions_species", "species"),
    )


class SequenceRecord(Base):
    __tablename__ = "sequence_records"

    accession: Mapped[str] = mapped_column(
        String(50), ForeignKey("accessions.accession"), primary_key=True
    )
    sequence_hash: Mapped[Optional[str]] = mapped_column(String(128), nullable=True)
    sequence_length: Mapped[Optional[int]] = mapped_column(Integer, nullable=True)
    sequence_source: Mapped[Optional[str]] = mapped_column(String(50), nullable=True)  # ncbi, local
    fasta_path: Mapped[Optional[str]] = mapped_column(String(1000), nullable=True)
    genbank_path: Mapped[Optional[str]] = mapped_column(String(1000), nullable=True)
    is_circular: Mapped[Optional[bool]] = mapped_column(Boolean, nullable=True, default=False)
    download_status: Mapped[Optional[str]] = mapped_column(
        String(50), nullable=True, default="pending"
    )  # pending, downloaded, failed, skipped
    parse_status: Mapped[Optional[str]] = mapped_column(
        String(50), nullable=True, default="pending"
    )  # pending, parsed, failed

    accession_rel: Mapped[Accession] = relationship("Accession", back_populates="sequence_record")


class FeatureRecord(Base):
    __tablename__ = "feature_records"

    feature_id: Mapped[int] = mapped_column(Integer, primary_key=True, autoincrement=True)
    accession: Mapped[str] = mapped_column(
        String(50), ForeignKey("accessions.accession"), nullable=False
    )
    feature_type: Mapped[str] = mapped_column(String(50), nullable=False)  # CDS, gene, mRNA, etc.
    start: Mapped[int] = mapped_column(Integer, nullable=False)
    end: Mapped[int] = mapped_column(Integer, nullable=False)
    strand: Mapped[Optional[str]] = mapped_column(String(1), nullable=True)  # +, -, .
    gene_name: Mapped[Optional[str]] = mapped_column(String(255), nullable=True)
    product: Mapped[Optional[str]] = mapped_column(String(500), nullable=True)
    locus_tag: Mapped[Optional[str]] = mapped_column(String(255), nullable=True)
    annotation_source: Mapped[Optional[str]] = mapped_column(String(50), nullable=True)  # genbank, gff3

    accession_rel: Mapped[Accession] = relationship("Accession", back_populates="features")
    ssr_annotations: Mapped[list[SSRAnnotation]] = relationship("SSRAnnotation", back_populates="feature")

    __table_args__ = (
        Index("ix_features_accession", "accession"),
        Index("ix_features_accession_start_end", "accession", "start", "end"),
    )


class SSRRecord(Base):
    __tablename__ = "ssr_records"

    ssr_id: Mapped[int] = mapped_column(Integer, primary_key=True, autoincrement=True)
    accession: Mapped[str] = mapped_column(
        String(50), ForeignKey("accessions.accession"), nullable=False
    )
    start: Mapped[int] = mapped_column(Integer, nullable=False)
    end: Mapped[int] = mapped_column(Integer, nullable=False)
    motif_raw: Mapped[str] = mapped_column(String(6), nullable=False)
    motif_canonical: Mapped[str] = mapped_column(String(6), nullable=False)
    repeat_units: Mapped[int] = mapped_column(Integer, nullable=False)
    motif_size: Mapped[int] = mapped_column(Integer, nullable=False)  # 1-6
    repeat_length_bp: Mapped[int] = mapped_column(Integer, nullable=False)
    strand: Mapped[str] = mapped_column(String(1), nullable=False, default="+")
    actual_repeat: Mapped[Optional[str]] = mapped_column(Text, nullable=True)
    detector_version: Mapped[Optional[str]] = mapped_column(String(50), nullable=True)
    run_id: Mapped[Optional[int]] = mapped_column(
        Integer, ForeignKey("runs.run_id"), nullable=True
    )

    accession_rel: Mapped[Accession] = relationship("Accession", back_populates="ssr_records")
    run: Mapped[Optional[Run]] = relationship("Run", back_populates="ssr_records")
    annotations: Mapped[list[SSRAnnotation]] = relationship("SSRAnnotation", back_populates="ssr")

    __table_args__ = (
        Index("ix_ssr_accession", "accession"),
        Index("ix_ssr_accession_start_end", "accession", "start", "end"),
        Index("ix_ssr_motif_canonical", "motif_canonical"),
        Index("ix_ssr_motif_size", "motif_size"),
        Index("ix_ssr_run_id", "run_id"),
    )


class SSRAnnotation(Base):
    __tablename__ = "ssr_annotations"

    annotation_id: Mapped[int] = mapped_column(Integer, primary_key=True, autoincrement=True)
    ssr_id: Mapped[int] = mapped_column(
        Integer, ForeignKey("ssr_records.ssr_id"), nullable=False
    )
    accession: Mapped[str] = mapped_column(
        String(50), ForeignKey("accessions.accession"), nullable=False
    )
    feature_id: Mapped[Optional[int]] = mapped_column(
        Integer, ForeignKey("feature_records.feature_id"), nullable=True
    )
    gene_name: Mapped[Optional[str]] = mapped_column(String(255), nullable=True)
    region_class: Mapped[str] = mapped_column(
        String(50), nullable=False, default="intergenic"
    )  # CDS, UTR5, UTR3, intergenic, other
    overlap_bp: Mapped[Optional[int]] = mapped_column(Integer, nullable=True)

    ssr: Mapped[SSRRecord] = relationship("SSRRecord", back_populates="annotations")
    accession_rel: Mapped[Accession] = relationship("Accession", back_populates="ssr_annotations")
    feature: Mapped[Optional[FeatureRecord]] = relationship("FeatureRecord", back_populates="ssr_annotations")

    __table_args__ = (
        Index("ix_ssr_annot_ssr_id", "ssr_id"),
        Index("ix_ssr_annot_accession", "accession"),
        Index("ix_ssr_annot_gene_name", "gene_name"),
    )


class AccessionMetrics(Base):
    __tablename__ = "accession_metrics"

    id: Mapped[int] = mapped_column(Integer, primary_key=True, autoincrement=True)
    accession: Mapped[str] = mapped_column(
        String(50), ForeignKey("accessions.accession"), nullable=False
    )
    run_id: Mapped[int] = mapped_column(
        Integer, ForeignKey("runs.run_id"), nullable=False
    )
    ssr_count_total: Mapped[int] = mapped_column(Integer, nullable=False, default=0)
    ssr_bp_total: Mapped[int] = mapped_column(Integer, nullable=False, default=0)
    ra: Mapped[Optional[float]] = mapped_column(Float, nullable=True)  # SSR count / genome size (kb)
    rd: Mapped[Optional[float]] = mapped_column(Float, nullable=True)  # SSR bp / genome size (Mb)
    mono_count: Mapped[int] = mapped_column(Integer, nullable=False, default=0)
    di_count: Mapped[int] = mapped_column(Integer, nullable=False, default=0)
    tri_count: Mapped[int] = mapped_column(Integer, nullable=False, default=0)
    tetra_count: Mapped[int] = mapped_column(Integer, nullable=False, default=0)
    penta_count: Mapped[int] = mapped_column(Integer, nullable=False, default=0)
    hexa_count: Mapped[int] = mapped_column(Integer, nullable=False, default=0)
    dominant_motif: Mapped[Optional[str]] = mapped_column(String(6), nullable=True)

    accession_rel: Mapped[Accession] = relationship("Accession", back_populates="metrics")
    run: Mapped[Run] = relationship("Run", back_populates="accession_metrics")

    __table_args__ = (
        UniqueConstraint("accession", "run_id", name="uq_accession_metrics_acc_run"),
        Index("ix_metrics_accession", "accession"),
        Index("ix_metrics_run_id", "run_id"),
    )


class StatisticalResult(Base):
    __tablename__ = "statistical_results"

    result_id: Mapped[int] = mapped_column(Integer, primary_key=True, autoincrement=True)
    run_id: Mapped[int] = mapped_column(
        Integer, ForeignKey("runs.run_id"), nullable=False
    )
    analysis_name: Mapped[str] = mapped_column(String(255), nullable=False)
    grouping: Mapped[Optional[str]] = mapped_column(String(255), nullable=True)
    metric: Mapped[Optional[str]] = mapped_column(String(255), nullable=True)
    test_name: Mapped[str] = mapped_column(String(100), nullable=False)
    statistic: Mapped[Optional[float]] = mapped_column(Float, nullable=True)
    p_value: Mapped[Optional[float]] = mapped_column(Float, nullable=True)
    p_value_corrected: Mapped[Optional[float]] = mapped_column(Float, nullable=True)
    effect_size: Mapped[Optional[float]] = mapped_column(Float, nullable=True)
    ci_lower: Mapped[Optional[float]] = mapped_column(Float, nullable=True)
    ci_upper: Mapped[Optional[float]] = mapped_column(Float, nullable=True)
    n: Mapped[Optional[int]] = mapped_column(Integer, nullable=True)
    metadata_json: Mapped[Optional[str]] = mapped_column(Text, nullable=True)

    run: Mapped[Run] = relationship("Run", back_populates="statistical_results")

    __table_args__ = (
        Index("ix_stats_run_id", "run_id"),
        Index("ix_stats_analysis_name", "analysis_name"),
    )


class StageCheckpoint(Base):
    """Tracks completion status of each pipeline stage within a run."""

    __tablename__ = "stage_checkpoints"

    id: Mapped[int] = mapped_column(Integer, primary_key=True, autoincrement=True)
    run_id: Mapped[int] = mapped_column(
        Integer, ForeignKey("runs.run_id"), nullable=False
    )
    stage: Mapped[str] = mapped_column(String(50), nullable=False)
    status: Mapped[str] = mapped_column(
        String(50), nullable=False, default="pending"
    )  # pending, running, completed, failed, skipped
    started_at: Mapped[Optional[datetime.datetime]] = mapped_column(DateTime, nullable=True)
    finished_at: Mapped[Optional[datetime.datetime]] = mapped_column(DateTime, nullable=True)
    processed_count: Mapped[int] = mapped_column(Integer, nullable=False, default=0)
    skipped_count: Mapped[int] = mapped_column(Integer, nullable=False, default=0)
    failed_count: Mapped[int] = mapped_column(Integer, nullable=False, default=0)

    __table_args__ = (
        UniqueConstraint("run_id", "stage", name="uq_run_stage"),
        Index("ix_stage_checkpoints_run_id", "run_id"),
    )


class FailedAccession(Base):
    """Tracks accessions that failed at a specific pipeline stage for retry."""

    __tablename__ = "failed_accessions"

    id: Mapped[int] = mapped_column(Integer, primary_key=True, autoincrement=True)
    run_id: Mapped[int] = mapped_column(
        Integer, ForeignKey("runs.run_id"), nullable=False
    )
    accession: Mapped[str] = mapped_column(String(50), nullable=False)
    stage: Mapped[str] = mapped_column(String(50), nullable=False)
    error_message: Mapped[Optional[str]] = mapped_column(Text, nullable=True)
    retry_count: Mapped[int] = mapped_column(Integer, nullable=False, default=0)
    resolved: Mapped[bool] = mapped_column(Boolean, nullable=False, default=False)
    created_at: Mapped[Optional[datetime.datetime]] = mapped_column(
        DateTime, nullable=True,
        default=lambda: datetime.datetime.now(datetime.UTC),
    )

    __table_args__ = (
        Index("ix_failed_accessions_run_id", "run_id"),
        Index("ix_failed_accessions_stage", "stage"),
    )
