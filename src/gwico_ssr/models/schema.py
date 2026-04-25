"""SQLAlchemy ORM models for GWICO-SSR canonical data model.

Core entities: Dataset, Run, Accession, SequenceRecord, FeatureRecord,
SSRRecord, SSRAnnotation, AccessionMetrics, StatisticalResult, BatchArtifact,
BatchNormalizedRecord, StageCheckpoint, FailedAccession, LineageSource,
LineageAssignment, LineageMetrics (Chunk 10), TreeSource, TreeMapping,
TreeMetrics (Chunk 11).
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
    lineage_assignments: Mapped[list[LineageAssignment]] = relationship("LineageAssignment", back_populates="accession_rel")
    tree_mappings: Mapped[list[TreeMapping]] = relationship("TreeMapping", back_populates="accession_rel")

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
    # Chunk 5: IMEX detector fields (backward compatible, nullable)
    is_perfect: Mapped[Optional[bool]] = mapped_column(Boolean, nullable=True, default=True)
    repeat_class: Mapped[str] = mapped_column(String(20), nullable=False, default="perfect")  # perfect, imperfect, compound_component
    imperfection_pct: Mapped[Optional[float]] = mapped_column(Float, nullable=True)  # % mismatch
    num_substitutions: Mapped[Optional[int]] = mapped_column(Integer, nullable=True, default=0)
    num_indels: Mapped[Optional[int]] = mapped_column(Integer, nullable=True, default=0)
    imperfection_cigar: Mapped[Optional[str]] = mapped_column(Text, nullable=True)

    accession_rel: Mapped[Accession] = relationship("Accession", back_populates="ssr_records")
    run: Mapped[Optional[Run]] = relationship("Run", back_populates="ssr_records")
    annotations: Mapped[list[SSRAnnotation]] = relationship("SSRAnnotation", back_populates="ssr")
    compound_components: Mapped[list[CompoundSSRComponent]] = relationship(
        "CompoundSSRComponent", back_populates="ssr_record", foreign_keys="CompoundSSRComponent.ssr_id"
    )

    __table_args__ = (
        Index("ix_ssr_accession", "accession"),
        Index("ix_ssr_accession_start_end", "accession", "start", "end"),
        Index("ix_ssr_motif_canonical", "motif_canonical"),
        Index("ix_ssr_motif_size", "motif_size"),
        Index("ix_ssr_run_id", "run_id"),
        Index("ix_ssr_repeat_class", "repeat_class"),
        Index("ix_ssr_is_perfect", "is_perfect"),
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
    # Chunk 8: Batch-aware downstream
    repeat_class: Mapped[str] = mapped_column(String(20), nullable=False, default="perfect")  # perfect, imperfect, compound_component

    ssr: Mapped[SSRRecord] = relationship("SSRRecord", back_populates="annotations")
    accession_rel: Mapped[Accession] = relationship("Accession", back_populates="ssr_annotations")
    feature: Mapped[Optional[FeatureRecord]] = relationship("FeatureRecord", back_populates="ssr_annotations")

    __table_args__ = (
        Index("ix_ssr_annot_ssr_id", "ssr_id"),
        Index("ix_ssr_annot_accession", "accession"),
        Index("ix_ssr_annot_gene_name", "gene_name"),
    )


class CompoundSSR(Base):
    """Chunk 5: Compound SSR records for chained imperfect/perfect repeats.
    
    Represents a sequence of multiple SSRs chained together with gaps <= dmax.
    """
    __tablename__ = "compound_ssrs"

    compound_id: Mapped[int] = mapped_column(Integer, primary_key=True, autoincrement=True)
    run_id: Mapped[Optional[int]] = mapped_column(
        Integer, ForeignKey("runs.run_id"), nullable=True
    )
    accession: Mapped[str] = mapped_column(
        String(50), ForeignKey("accessions.accession"), nullable=False
    )
    start: Mapped[int] = mapped_column(Integer, nullable=False)  # start of first component
    end: Mapped[int] = mapped_column(Integer, nullable=False)  # end of last component
    component_count: Mapped[int] = mapped_column(Integer, nullable=False)
    total_repeat_length_bp: Mapped[int] = mapped_column(Integer, nullable=False)
    dmax_used: Mapped[int] = mapped_column(Integer, nullable=False)  # gap threshold used
    standardization_level: Mapped[str] = mapped_column(String(10), nullable=False, default="L2")  # L0, L1, L2, Full
    strand: Mapped[str] = mapped_column(String(1), nullable=False, default="+")

    components: Mapped[list[CompoundSSRComponent]] = relationship(
        "CompoundSSRComponent", back_populates="compound", cascade="all, delete-orphan"
    )

    __table_args__ = (
        Index("ix_compound_accession", "accession"),
        Index("ix_compound_accession_start_end", "accession", "start", "end"),
        Index("ix_compound_run_id", "run_id"),
        Index("ix_compound_standardization_level", "standardization_level"),
    )


class CompoundSSRComponent(Base):
    """Chunk 5: Individual components of a compound SSR.
    
    Each row represents one SSR within a compound SSR chain.
    """
    __tablename__ = "compound_ssr_components"

    component_id: Mapped[int] = mapped_column(Integer, primary_key=True, autoincrement=True)
    compound_id: Mapped[int] = mapped_column(
        Integer, ForeignKey("compound_ssrs.compound_id"), nullable=False
    )
    ssr_id: Mapped[int] = mapped_column(
        Integer, ForeignKey("ssr_records.ssr_id"), nullable=False
    )
    component_order: Mapped[int] = mapped_column(Integer, nullable=False)  # 0, 1, 2, ...
    gap_to_next_bp: Mapped[Optional[int]] = mapped_column(Integer, nullable=True)  # None if last component

    compound: Mapped[CompoundSSR] = relationship("CompoundSSR", back_populates="components")
    ssr_record: Mapped[SSRRecord] = relationship(
        "SSRRecord", back_populates="compound_components", foreign_keys=[ssr_id]
    )

    __table_args__ = (
        Index("ix_comp_compound_id", "compound_id"),
        Index("ix_comp_ssr_id", "ssr_id"),
        Index("ix_comp_order", "compound_id", "component_order"),
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
    # Chunk 8: Batch-aware downstream repeat_class breakdown
    perfect_count: Mapped[int] = mapped_column(Integer, nullable=False, default=0)
    imperfect_count: Mapped[int] = mapped_column(Integer, nullable=False, default=0)
    compound_component_count: Mapped[int] = mapped_column(Integer, nullable=False, default=0)
    perfect_bp_total: Mapped[int] = mapped_column(Integer, nullable=False, default=0)
    imperfect_bp_total: Mapped[int] = mapped_column(Integer, nullable=False, default=0)
    compound_component_bp_total: Mapped[int] = mapped_column(Integer, nullable=False, default=0)

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


class BatchArtifact(Base):
    """Tracks persisted raw composite artifacts used in batch-first workflows."""

    __tablename__ = "batch_artifacts"

    artifact_id: Mapped[int] = mapped_column(Integer, primary_key=True, autoincrement=True)
    dataset_id: Mapped[int] = mapped_column(Integer, ForeignKey("datasets.dataset_id"), nullable=False)
    run_id: Mapped[Optional[int]] = mapped_column(Integer, ForeignKey("runs.run_id"), nullable=True)
    file_type: Mapped[str] = mapped_column(String(20), nullable=False)  # fasta, genbank
    source_type: Mapped[str] = mapped_column(String(30), nullable=False)  # ncbi, local
    artifact_path: Mapped[str] = mapped_column(String(1000), nullable=False)
    checksum_sha256: Mapped[Optional[str]] = mapped_column(String(128), nullable=True)
    accession_count: Mapped[Optional[int]] = mapped_column(Integer, nullable=True)
    request_batch_size: Mapped[Optional[int]] = mapped_column(Integer, nullable=True)
    artifact_batch_size: Mapped[Optional[int]] = mapped_column(Integer, nullable=True)
    manifest_path: Mapped[Optional[str]] = mapped_column(String(1000), nullable=True)
    status: Mapped[str] = mapped_column(String(30), nullable=False, default="created")
    created_at: Mapped[datetime.datetime] = mapped_column(
        DateTime, nullable=False, default=lambda: datetime.datetime.now(datetime.UTC)
    )

    __table_args__ = (
        Index("ix_batch_artifacts_dataset_id", "dataset_id"),
        Index("ix_batch_artifacts_run_id", "run_id"),
        Index("ix_batch_artifacts_file_type", "file_type"),
        Index("ix_batch_artifacts_status", "status"),
        UniqueConstraint("artifact_path", name="uq_batch_artifacts_artifact_path"),
    )


class BatchNormalizedRecord(Base):
    """Links normalized accession outputs to their source batch artifact."""

    __tablename__ = "batch_normalized_records"

    normalized_id: Mapped[int] = mapped_column(Integer, primary_key=True, autoincrement=True)
    artifact_id: Mapped[int] = mapped_column(
        Integer, ForeignKey("batch_artifacts.artifact_id"), nullable=False
    )
    accession: Mapped[str] = mapped_column(String(50), nullable=False)
    record_index: Mapped[int] = mapped_column(Integer, nullable=False)
    normalized_fasta_path: Mapped[Optional[str]] = mapped_column(String(1000), nullable=True)
    normalized_genbank_path: Mapped[Optional[str]] = mapped_column(String(1000), nullable=True)
    checksum_sha256: Mapped[Optional[str]] = mapped_column(String(128), nullable=True)
    duplicate_policy: Mapped[Optional[str]] = mapped_column(String(30), nullable=True)
    parse_status: Mapped[str] = mapped_column(String(30), nullable=False, default="pending")
    created_at: Mapped[datetime.datetime] = mapped_column(
        DateTime, nullable=False, default=lambda: datetime.datetime.now(datetime.UTC)
    )

    __table_args__ = (
        Index("ix_batch_norm_records_artifact_id", "artifact_id"),
        Index("ix_batch_norm_records_accession", "accession"),
        Index("ix_batch_norm_records_parse_status", "parse_status"),
        UniqueConstraint(
            "artifact_id",
            "accession",
            "record_index",
            name="uq_batch_norm_artifact_acc_record",
        ),
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


class LineageSource(Base):
    """Chunk 10: Metadata about lineage assignment source and version."""

    __tablename__ = "lineage_sources"

    source_id: Mapped[int] = mapped_column(Integer, primary_key=True, autoincrement=True)
    source_name: Mapped[str] = mapped_column(String(50), nullable=False)  # pango, nextstrain, user-defined
    version: Mapped[str] = mapped_column(String(50), nullable=False)  # e.g., "v1.23", "2024-Q1"
    release_date: Mapped[Optional[str]] = mapped_column(String(50), nullable=True)  # ISO 8601 date
    description: Mapped[Optional[str]] = mapped_column(Text, nullable=True)
    created_at: Mapped[datetime.datetime] = mapped_column(
        DateTime, nullable=False, default=lambda: datetime.datetime.now(datetime.UTC)
    )

    assignments: Mapped[list[LineageAssignment]] = relationship(
        "LineageAssignment", back_populates="source"
    )

    __table_args__ = (
        Index("ix_lineage_source_name_version", "source_name", "version"),
    )


class LineageAssignment(Base):
    """Chunk 10: Accession-to-lineage mapping with provenance and audit trail."""

    __tablename__ = "lineage_assignments"

    assignment_id: Mapped[int] = mapped_column(Integer, primary_key=True, autoincrement=True)
    accession: Mapped[str] = mapped_column(
        String(50), ForeignKey("accessions.accession"), nullable=False
    )
    source_id: Mapped[int] = mapped_column(
        Integer, ForeignKey("lineage_sources.source_id"), nullable=False
    )
    lineage_label: Mapped[str] = mapped_column(String(100), nullable=False)  # e.g., "BA.1", "XEC", "unassigned"
    confidence: Mapped[Optional[float]] = mapped_column(Float, nullable=True)  # 0.0-1.0 if available
    assignment_date: Mapped[Optional[str]] = mapped_column(String(50), nullable=True)  # When assignment was made
    is_deprecated: Mapped[bool] = mapped_column(Boolean, nullable=False, default=False)  # Mark if lineage deprecated
    notes: Mapped[Optional[str]] = mapped_column(Text, nullable=True)  # User notes, conflict info, etc.
    created_at: Mapped[datetime.datetime] = mapped_column(
        DateTime, nullable=False, default=lambda: datetime.datetime.now(datetime.UTC)
    )

    accession_rel: Mapped[Accession] = relationship("Accession", back_populates="lineage_assignments")
    source: Mapped[LineageSource] = relationship("LineageSource", back_populates="assignments")

    __table_args__ = (
        Index("ix_lineage_assign_accession", "accession"),
        Index("ix_lineage_assign_source_id", "source_id"),
        Index("ix_lineage_assign_label", "lineage_label"),
        Index("ix_lineage_assign_deprecated", "is_deprecated"),
    )


class LineageMetrics(Base):
    """Chunk 10: SSR metrics aggregated by lineage and repeat class."""

    __tablename__ = "lineage_metrics"

    metric_id: Mapped[int] = mapped_column(Integer, primary_key=True, autoincrement=True)
    run_id: Mapped[int] = mapped_column(
        Integer, ForeignKey("runs.run_id"), nullable=False
    )
    source_id: Mapped[int] = mapped_column(
        Integer, ForeignKey("lineage_sources.source_id"), nullable=False
    )
    lineage_label: Mapped[str] = mapped_column(String(100), nullable=False)
    accession_count: Mapped[int] = mapped_column(Integer, nullable=False, default=0)
    ssr_count_total: Mapped[int] = mapped_column(Integer, nullable=False, default=0)
    ssr_bp_total: Mapped[int] = mapped_column(Integer, nullable=False, default=0)
    perfect_count: Mapped[int] = mapped_column(Integer, nullable=False, default=0)
    imperfect_count: Mapped[int] = mapped_column(Integer, nullable=False, default=0)
    compound_component_count: Mapped[int] = mapped_column(Integer, nullable=False, default=0)
    perfect_bp_total: Mapped[int] = mapped_column(Integer, nullable=False, default=0)
    imperfect_bp_total: Mapped[int] = mapped_column(Integer, nullable=False, default=0)
    compound_component_bp_total: Mapped[int] = mapped_column(Integer, nullable=False, default=0)
    mean_ra: Mapped[Optional[float]] = mapped_column(Float, nullable=True)  # Mean repeat abundance
    mean_rd: Mapped[Optional[float]] = mapped_column(Float, nullable=True)  # Mean repeat density

    run: Mapped[Run] = relationship("Run")
    source: Mapped[LineageSource] = relationship("LineageSource")

    __table_args__ = (
        Index("ix_lineage_metrics_run_id", "run_id"),
        Index("ix_lineage_metrics_source_id", "source_id"),
        Index("ix_lineage_metrics_label", "lineage_label"),
    )


class TreeSource(Base):
    """Chunk 11: Metadata about phylogenetic tree source and version."""

    __tablename__ = "tree_sources"

    source_id: Mapped[int] = mapped_column(Integer, primary_key=True, autoincrement=True)
    source_name: Mapped[str] = mapped_column(String(100), nullable=False)
    version: Mapped[str] = mapped_column(String(50), nullable=False)
    file_name: Mapped[Optional[str]] = mapped_column(String(255), nullable=True)
    file_hash: Mapped[Optional[str]] = mapped_column(String(128), nullable=True)
    num_tips: Mapped[int] = mapped_column(Integer, nullable=False, default=0)
    release_date: Mapped[Optional[str]] = mapped_column(String(50), nullable=True)
    description: Mapped[Optional[str]] = mapped_column(Text, nullable=True)
    created_at: Mapped[datetime.datetime] = mapped_column(
        DateTime, nullable=False, default=lambda: datetime.datetime.now(datetime.UTC)
    )

    mappings: Mapped[list[TreeMapping]] = relationship("TreeMapping", back_populates="source")
    metrics: Mapped[list[TreeMetrics]] = relationship("TreeMetrics", back_populates="source")

    __table_args__ = (
        Index("ix_tree_source_name_version", "source_name", "version"),
    )


class TreeMapping(Base):
    """Chunk 11: Accession-to-tree-tip mapping with provenance and conflict tracking."""

    __tablename__ = "tree_mappings"

    mapping_id: Mapped[int] = mapped_column(Integer, primary_key=True, autoincrement=True)
    accession: Mapped[str] = mapped_column(
        String(50), ForeignKey("accessions.accession"), nullable=False
    )
    source_id: Mapped[int] = mapped_column(
        Integer, ForeignKey("tree_sources.source_id"), nullable=False
    )
    tip_label: Mapped[str] = mapped_column(String(255), nullable=False)
    is_confident: Mapped[bool] = mapped_column(Boolean, nullable=False, default=True)
    conflict_notes: Mapped[Optional[str]] = mapped_column(Text, nullable=True)
    created_at: Mapped[datetime.datetime] = mapped_column(
        DateTime, nullable=False, default=lambda: datetime.datetime.now(datetime.UTC)
    )

    accession_rel: Mapped[Accession] = relationship("Accession", back_populates="tree_mappings")
    source: Mapped[TreeSource] = relationship("TreeSource", back_populates="mappings")

    __table_args__ = (
        Index("ix_tree_mapping_accession", "accession"),
        Index("ix_tree_mapping_source_id", "source_id"),
        Index("ix_tree_mapping_tip_label", "tip_label"),
    )


class TreeMetrics(Base):
    """Chunk 11: SSR metrics aggregated by clade in phylogenetic tree."""

    __tablename__ = "tree_metrics"

    metric_id: Mapped[int] = mapped_column(Integer, primary_key=True, autoincrement=True)
    run_id: Mapped[int] = mapped_column(
        Integer, ForeignKey("runs.run_id"), nullable=False
    )
    source_id: Mapped[int] = mapped_column(
        Integer, ForeignKey("tree_sources.source_id"), nullable=False
    )
    clade_name: Mapped[str] = mapped_column(String(255), nullable=False)
    tip_count: Mapped[int] = mapped_column(Integer, nullable=False, default=0)
    accession_count: Mapped[int] = mapped_column(Integer, nullable=False, default=0)
    ssr_count_total: Mapped[int] = mapped_column(Integer, nullable=False, default=0)
    ssr_bp_total: Mapped[int] = mapped_column(Integer, nullable=False, default=0)
    perfect_count: Mapped[int] = mapped_column(Integer, nullable=False, default=0)
    imperfect_count: Mapped[int] = mapped_column(Integer, nullable=False, default=0)
    compound_component_count: Mapped[int] = mapped_column(Integer, nullable=False, default=0)
    perfect_bp_total: Mapped[int] = mapped_column(Integer, nullable=False, default=0)
    imperfect_bp_total: Mapped[int] = mapped_column(Integer, nullable=False, default=0)
    compound_component_bp_total: Mapped[int] = mapped_column(Integer, nullable=False, default=0)
    mean_ra: Mapped[Optional[float]] = mapped_column(Float, nullable=True)
    mean_rd: Mapped[Optional[float]] = mapped_column(Float, nullable=True)

    run: Mapped[Run] = relationship("Run")
    source: Mapped[TreeSource] = relationship("TreeSource", back_populates="metrics")

    __table_args__ = (
        Index("ix_tree_metrics_run_id", "run_id"),
        Index("ix_tree_metrics_source_id", "source_id"),
        Index("ix_tree_metrics_clade_name", "clade_name"),
    )
