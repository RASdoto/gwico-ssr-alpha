"""Tests for GWICO-SSR parsers: FASTA, GenBank, GFF3, and persistence."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pytest
from sqlalchemy import create_engine
from sqlalchemy.orm import Session, sessionmaker

from gwico_ssr.db import create_tables
from gwico_ssr.db.repository import (
    get_accession,
    get_features_for_accession,
    get_or_create_dataset,
    upsert_accession,
    upsert_sequence_record,
)
from gwico_ssr.models.schema import Base, SequenceRecord

FIXTURES_DIR = Path(__file__).parent.parent / "fixtures"


# ---------------------------------------------------------------------------
# Shared DB fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def db_session(tmp_path):
    """Create an in-memory SQLite session with schema."""
    engine = create_engine("sqlite:///:memory:")
    create_tables(engine)
    factory = sessionmaker(bind=engine)
    session = factory()
    yield session
    session.close()


@pytest.fixture
def db_session_with_accession(db_session):
    """DB session with a test dataset and accession pre-created."""
    ds = get_or_create_dataset(db_session, name="test_ds", source_type="csv")
    upsert_accession(
        db_session,
        accession="NC_045512.2",
        dataset_id=ds.dataset_id,
        species="SARS-CoV-2",
    )
    upsert_sequence_record(
        db_session,
        accession="NC_045512.2",
        download_status="downloaded",
    )
    return db_session


# ===========================================================================
# FASTA PARSER TESTS
# ===========================================================================

class TestFastaParser:

    def test_parse_valid_fasta(self):
        from gwico_ssr.parsers.fasta_parser import parse_fasta

        result = parse_fasta(FIXTURES_DIR / "sample.fasta")
        assert result.success
        assert result.record_count == 1
        rec = result.records[0]
        assert rec.accession == "NC_045512.2"
        assert rec.sequence_length > 0
        assert 0.0 < rec.gc_content < 1.0
        assert len(rec.sequence_hash) == 64  # SHA-256 hex

    def test_parse_multi_record_fasta(self):
        from gwico_ssr.parsers.fasta_parser import parse_fasta

        result = parse_fasta(FIXTURES_DIR / "multi_record.fasta")
        assert result.success
        assert result.record_count == 2
        accessions = [r.accession for r in result.records]
        assert "NC_045512.2" in accessions
        assert "OL672836.1" in accessions

    def test_parse_empty_sequence(self):
        from gwico_ssr.parsers.fasta_parser import parse_fasta

        result = parse_fasta(FIXTURES_DIR / "empty_sequence.fasta")
        # Empty sequence record should produce an error
        assert not result.success or result.record_count == 0

    def test_parse_malformed_fasta(self):
        from gwico_ssr.parsers.fasta_parser import parse_fasta

        result = parse_fasta(FIXTURES_DIR / "malformed.fasta")
        # Malformed file may parse with no valid records
        assert result.record_count == 0

    def test_parse_missing_file(self):
        from gwico_ssr.parsers.fasta_parser import parse_fasta

        result = parse_fasta("/nonexistent/file.fasta")
        assert not result.success
        assert len(result.errors) == 1
        assert "not found" in result.errors[0].message.lower()

    def test_parse_empty_file(self, tmp_path):
        from gwico_ssr.parsers.fasta_parser import parse_fasta

        empty = tmp_path / "empty.fasta"
        empty.write_text("")
        result = parse_fasta(empty)
        assert not result.success
        assert len(result.errors) == 1

    def test_gc_content_calculation(self):
        from gwico_ssr.parsers.fasta_parser import _compute_gc_content

        assert _compute_gc_content("GGCC") == 1.0
        assert _compute_gc_content("AATT") == 0.0
        assert _compute_gc_content("ATGC") == 0.5
        assert _compute_gc_content("") == 0.0
        assert _compute_gc_content("NNNN") == 0.0  # No definite bases

    def test_hash_deterministic(self):
        from gwico_ssr.parsers.fasta_parser import _compute_hash

        h1 = _compute_hash("ATCG")
        h2 = _compute_hash("ATCG")
        h3 = _compute_hash("atcg")  # case insensitive
        assert h1 == h2
        assert h1 == h3

    def test_extract_accession(self):
        from gwico_ssr.parsers.fasta_parser import _extract_accession

        assert _extract_accession("NC_045512.2") == "NC_045512.2"
        assert _extract_accession("NC_045512.2 some desc") == "NC_045512.2"

    def test_parse_fasta_multi(self):
        from gwico_ssr.parsers.fasta_parser import parse_fasta_multi

        results = parse_fasta_multi([
            FIXTURES_DIR / "sample.fasta",
            FIXTURES_DIR / "multi_record.fasta",
        ])
        assert len(results) == 2
        assert results[0].record_count == 1
        assert results[1].record_count == 2

    def test_to_dict(self):
        from gwico_ssr.parsers.fasta_parser import parse_fasta

        result = parse_fasta(FIXTURES_DIR / "sample.fasta")
        d = result.to_dict()
        assert "record_count" in d
        assert d["success"] is True


# ===========================================================================
# GENBANK PARSER TESTS
# ===========================================================================

class TestGenBankParser:

    def test_parse_valid_genbank(self):
        from gwico_ssr.parsers.genbank_parser import parse_genbank

        result = parse_genbank(FIXTURES_DIR / "sample.gb")
        assert result.success
        assert result.sequence_info is not None
        assert result.sequence_info.accession == "NC_045512.2"
        assert result.sequence_info.sequence_length == 444
        assert result.sequence_info.organism == "Severe acute respiratory syndrome coronavirus 2"

    def test_genbank_cds_extraction(self):
        from gwico_ssr.parsers.genbank_parser import parse_genbank

        result = parse_genbank(FIXTURES_DIR / "sample.gb")
        assert result.cds_count == 3  # ORF1ab, S, E
        cds_genes = [f.gene_name for f in result.features if f.feature_type == "CDS"]
        assert "ORF1ab" in cds_genes
        assert "S" in cds_genes
        assert "E" in cds_genes

    def test_genbank_feature_coordinates(self):
        from gwico_ssr.parsers.genbank_parser import parse_genbank

        result = parse_genbank(FIXTURES_DIR / "sample.gb")
        # Find ORF1ab CDS — GenBank 10..200 becomes 0-based: start=9, end=200
        orf1ab = [f for f in result.features if f.gene_name == "ORF1ab" and f.feature_type == "CDS"]
        assert len(orf1ab) == 1
        assert orf1ab[0].start == 9  # BioPython converts to 0-based
        assert orf1ab[0].end == 200
        assert orf1ab[0].strand == "+"

    def test_genbank_all_features(self):
        from gwico_ssr.parsers.genbank_parser import parse_genbank

        result = parse_genbank(FIXTURES_DIR / "sample.gb")
        # 3 CDS + 3 gene = 6 features of interest
        assert result.feature_count == 6

    def test_genbank_locus_tags(self):
        from gwico_ssr.parsers.genbank_parser import parse_genbank

        result = parse_genbank(FIXTURES_DIR / "sample.gb")
        cds_feats = [f for f in result.features if f.feature_type == "CDS"]
        for f in cds_feats:
            assert f.locus_tag is not None
            assert f.locus_tag.startswith("GU280_gp")

    def test_genbank_annotation_source(self):
        from gwico_ssr.parsers.genbank_parser import parse_genbank

        result = parse_genbank(FIXTURES_DIR / "sample.gb")
        for f in result.features:
            assert f.annotation_source == "genbank"

    def test_genbank_gc_content(self):
        from gwico_ssr.parsers.genbank_parser import parse_genbank

        result = parse_genbank(FIXTURES_DIR / "sample.gb")
        assert result.sequence_info.gc_content is not None
        assert 0.0 < result.sequence_info.gc_content < 1.0

    def test_genbank_hash(self):
        from gwico_ssr.parsers.genbank_parser import parse_genbank

        result = parse_genbank(FIXTURES_DIR / "sample.gb")
        assert result.sequence_info.sequence_hash is not None
        assert len(result.sequence_info.sequence_hash) == 64

    def test_genbank_missing_file(self):
        from gwico_ssr.parsers.genbank_parser import parse_genbank

        result = parse_genbank("/nonexistent/file.gb")
        assert not result.success
        assert len(result.errors) == 1

    def test_genbank_empty_file(self, tmp_path):
        from gwico_ssr.parsers.genbank_parser import parse_genbank

        empty = tmp_path / "empty.gb"
        empty.write_text("")
        result = parse_genbank(empty)
        assert not result.success

    def test_genbank_to_dict(self):
        from gwico_ssr.parsers.genbank_parser import parse_genbank

        result = parse_genbank(FIXTURES_DIR / "sample.gb")
        d = result.to_dict()
        assert d["success"] is True
        assert d["cds_count"] == 3
        assert d["accession"] == "NC_045512.2"

    def test_normalize_strand(self):
        from gwico_ssr.parsers.genbank_parser import _normalize_strand

        assert _normalize_strand(1) == "+"
        assert _normalize_strand(-1) == "-"
        assert _normalize_strand(0) == "."
        assert _normalize_strand(None) == "."

    def test_parse_genbank_multi(self):
        from gwico_ssr.parsers.genbank_parser import parse_genbank_multi

        results = parse_genbank_multi([FIXTURES_DIR / "sample.gb"])
        assert len(results) == 1
        assert results[0].success


# ===========================================================================
# GFF3 PARSER TESTS
# ===========================================================================

class TestGFF3Parser:

    def test_parse_valid_gff3(self):
        from gwico_ssr.parsers.gff3_parser import parse_gff3

        result = parse_gff3(FIXTURES_DIR / "sample.gff3")
        assert result.success
        assert result.feature_count > 0

    def test_gff3_cds_extraction(self):
        from gwico_ssr.parsers.gff3_parser import parse_gff3

        result = parse_gff3(FIXTURES_DIR / "sample.gff3")
        assert result.cds_count == 4  # ORF1ab, S, E, N
        cds_genes = [f.gene_name for f in result.features if f.feature_type == "CDS"]
        assert "ORF1ab" in cds_genes
        assert "S" in cds_genes
        assert "E" in cds_genes
        assert "N" in cds_genes

    def test_gff3_coordinate_conversion(self):
        from gwico_ssr.parsers.gff3_parser import parse_gff3

        result = parse_gff3(FIXTURES_DIR / "sample.gff3")
        # GFF3 line: NC_045512.2 RefSeq gene 266 21555 . + . ...
        # Converted to 0-based: start=265, end=21555
        orf1ab_genes = [f for f in result.features if f.gene_name == "ORF1ab" and f.feature_type == "gene"]
        assert len(orf1ab_genes) == 1
        assert orf1ab_genes[0].start == 265
        assert orf1ab_genes[0].end == 21555

    def test_gff3_strand_normalization(self):
        from gwico_ssr.parsers.gff3_parser import parse_gff3

        result = parse_gff3(FIXTURES_DIR / "sample.gff3")
        # N gene is on minus strand
        n_feats = [f for f in result.features if f.gene_name == "N"]
        for f in n_feats:
            assert f.strand == "-"
        # Other genes are on plus strand
        s_feats = [f for f in result.features if f.gene_name == "S"]
        for f in s_feats:
            assert f.strand == "+"

    def test_gff3_annotation_source(self):
        from gwico_ssr.parsers.gff3_parser import parse_gff3

        result = parse_gff3(FIXTURES_DIR / "sample.gff3")
        for f in result.features:
            assert f.annotation_source == "gff3"

    def test_gff3_default_accession(self):
        from gwico_ssr.parsers.gff3_parser import parse_gff3

        result = parse_gff3(FIXTURES_DIR / "sample.gff3", default_accession="CUSTOM_ACC")
        for f in result.features:
            assert f.accession == "CUSTOM_ACC"

    def test_gff3_without_default_accession(self):
        from gwico_ssr.parsers.gff3_parser import parse_gff3

        result = parse_gff3(FIXTURES_DIR / "sample.gff3")
        for f in result.features:
            assert f.accession == "NC_045512.2"  # seqid column

    def test_gff3_locus_tags(self):
        from gwico_ssr.parsers.gff3_parser import parse_gff3

        result = parse_gff3(FIXTURES_DIR / "sample.gff3")
        cds_feats = [f for f in result.features if f.feature_type == "CDS"]
        for f in cds_feats:
            assert f.locus_tag is not None

    def test_gff3_products(self):
        from gwico_ssr.parsers.gff3_parser import parse_gff3

        result = parse_gff3(FIXTURES_DIR / "sample.gff3")
        cds_feats = [f for f in result.features if f.feature_type == "CDS"]
        products = [f.product for f in cds_feats]
        assert "surface glycoprotein" in products
        assert "envelope protein" in products

    def test_gff3_missing_file(self):
        from gwico_ssr.parsers.gff3_parser import parse_gff3

        result = parse_gff3("/nonexistent/file.gff3")
        assert not result.success
        assert len(result.errors) == 1

    def test_gff3_empty_file(self, tmp_path):
        from gwico_ssr.parsers.gff3_parser import parse_gff3

        empty = tmp_path / "empty.gff3"
        empty.write_text("")
        result = parse_gff3(empty)
        assert not result.success

    def test_gff3_malformed(self):
        from gwico_ssr.parsers.gff3_parser import parse_gff3

        result = parse_gff3(FIXTURES_DIR / "malformed.gff3")
        # Should have errors for bad column counts
        assert len(result.errors) > 0

    def test_gff3_to_dict(self):
        from gwico_ssr.parsers.gff3_parser import parse_gff3

        result = parse_gff3(FIXTURES_DIR / "sample.gff3")
        d = result.to_dict()
        assert d["success"] is True
        assert d["cds_count"] == 4

    def test_gff3_directives(self):
        from gwico_ssr.parsers.gff3_parser import parse_gff3

        result = parse_gff3(FIXTURES_DIR / "sample.gff3")
        assert any("gff-version" in d for d in result.directives)
        assert any("sequence-region" in d for d in result.directives)

    def test_parse_attributes(self):
        from gwico_ssr.parsers.gff3_parser import _parse_attributes

        attrs = _parse_attributes("gene=ORF1ab;product=ORF1ab%20polyprotein;locus_tag=GU280_gp01")
        assert attrs["gene"] == "ORF1ab"
        assert attrs["product"] == "ORF1ab polyprotein"  # URL-decoded
        assert attrs["locus_tag"] == "GU280_gp01"

    def test_parse_attributes_empty(self):
        from gwico_ssr.parsers.gff3_parser import _parse_attributes

        assert _parse_attributes("") == {}
        assert _parse_attributes(".") == {}

    def test_parse_gff3_multi(self):
        from gwico_ssr.parsers.gff3_parser import parse_gff3_multi

        results = parse_gff3_multi([FIXTURES_DIR / "sample.gff3"])
        assert len(results) == 1
        assert results[0].success


# ===========================================================================
# PERSISTENCE TESTS
# ===========================================================================

class TestParsePersistence:

    def test_persist_fasta(self, db_session_with_accession):
        from gwico_ssr.parsers.fasta_parser import ParsedSequence
        from gwico_ssr.parsers.persist import persist_parsed_fasta

        parsed = ParsedSequence(
            accession="NC_045512.2",
            sequence="ATCGATCG",
            sequence_length=8,
            sequence_hash="abc123",
            gc_content=0.5,
        )
        persist_parsed_fasta(db_session_with_accession, "NC_045512.2", parsed)

        seq_rec = db_session_with_accession.get(SequenceRecord, "NC_045512.2")
        assert seq_rec.sequence_hash == "abc123"
        assert seq_rec.sequence_length == 8
        assert seq_rec.parse_status == "parsed"

        acc = get_accession(db_session_with_accession, "NC_045512.2")
        assert acc.gc_content == 0.5
        assert acc.genome_length == 8

    def test_persist_features(self, db_session_with_accession):
        from gwico_ssr.parsers.genbank_parser import ParsedFeature
        from gwico_ssr.parsers.persist import persist_parsed_features

        features = [
            ParsedFeature(
                accession="NC_045512.2",
                feature_type="CDS",
                start=10,
                end=200,
                strand="+",
                gene_name="ORF1ab",
                product="polyprotein",
                locus_tag="gp01",
            ),
            ParsedFeature(
                accession="NC_045512.2",
                feature_type="CDS",
                start=220,
                end=400,
                strand="+",
                gene_name="S",
                product="spike",
                locus_tag="gp02",
            ),
        ]
        count = persist_parsed_features(db_session_with_accession, features, "NC_045512.2")
        assert count == 2

        db_features = get_features_for_accession(db_session_with_accession, "NC_045512.2")
        assert len(db_features) == 2
        gene_names = [f.gene_name for f in db_features]
        assert "ORF1ab" in gene_names
        assert "S" in gene_names

    def test_persist_features_replace(self, db_session_with_accession):
        from gwico_ssr.parsers.genbank_parser import ParsedFeature
        from gwico_ssr.parsers.persist import persist_parsed_features

        feat1 = [ParsedFeature(
            accession="NC_045512.2", feature_type="CDS",
            start=10, end=200, strand="+", gene_name="OLD",
        )]
        persist_parsed_features(db_session_with_accession, feat1, "NC_045512.2")

        feat2 = [ParsedFeature(
            accession="NC_045512.2", feature_type="CDS",
            start=10, end=200, strand="+", gene_name="NEW",
        )]
        persist_parsed_features(
            db_session_with_accession, feat2, "NC_045512.2", replace_existing=True
        )

        db_features = get_features_for_accession(db_session_with_accession, "NC_045512.2")
        assert len(db_features) == 1
        assert db_features[0].gene_name == "NEW"

    def test_parse_and_persist_fasta_genbank(self, db_session_with_accession, tmp_path):
        """End-to-end parse and persist from fixture files."""
        import shutil
        from gwico_ssr.parsers.persist import parse_and_persist

        # Copy fixtures to a deterministic layout
        fasta_dir = tmp_path / "sequences" / "fasta"
        gb_dir = tmp_path / "sequences" / "genbank"
        fasta_dir.mkdir(parents=True)
        gb_dir.mkdir(parents=True)

        # We need to use the accession from the genbank file (NC_045512)
        # Create a matching accession record
        ds = get_or_create_dataset(db_session_with_accession, name="test_ds2", source_type="csv")
        upsert_accession(
            db_session_with_accession,
            accession="NC_045512",
            dataset_id=ds.dataset_id,
        )
        upsert_sequence_record(
            db_session_with_accession,
            accession="NC_045512",
            download_status="downloaded",
        )

        shutil.copy(FIXTURES_DIR / "sample.fasta", fasta_dir / "NC_045512.fasta")
        shutil.copy(FIXTURES_DIR / "sample.gb", gb_dir / "NC_045512.gb")

        result = parse_and_persist(
            session=db_session_with_accession,
            accession="NC_045512",
            fasta_path=str(fasta_dir / "NC_045512.fasta"),
            genbank_path=str(gb_dir / "NC_045512.gb"),
        )

        # FASTA record ID is NC_045512.2, so it won't match "NC_045512" exactly
        # but single-record fallback should handle it
        assert result["fasta_parsed"] or result["genbank_parsed"]
        assert result["features_inserted"] > 0

    def test_parse_and_persist_skip_already_parsed(self, db_session_with_accession):
        from gwico_ssr.parsers.persist import parse_and_persist

        # Mark as already parsed
        upsert_sequence_record(
            db_session_with_accession,
            accession="NC_045512.2",
            parse_status="parsed",
        )

        result = parse_and_persist(
            session=db_session_with_accession,
            accession="NC_045512.2",
        )
        assert result.get("skipped") is True

    def test_parse_and_persist_force(self, db_session_with_accession, tmp_path):
        from gwico_ssr.parsers.persist import parse_and_persist
        import shutil

        upsert_sequence_record(
            db_session_with_accession,
            accession="NC_045512.2",
            parse_status="parsed",
        )

        fasta_dir = tmp_path / "sequences" / "fasta"
        fasta_dir.mkdir(parents=True)
        shutil.copy(FIXTURES_DIR / "sample.fasta", fasta_dir / "NC_045512.2.fasta")

        result = parse_and_persist(
            session=db_session_with_accession,
            accession="NC_045512.2",
            fasta_path=str(fasta_dir / "NC_045512.2.fasta"),
            force=True,
        )
        assert result.get("skipped") is not True
        assert result["fasta_parsed"] is True

    def test_parse_dataset_accessions(self, db_session_with_accession, tmp_path):
        import shutil
        from gwico_ssr.parsers.persist import parse_dataset_accessions

        fasta_dir = tmp_path / "sequences" / "fasta"
        gb_dir = tmp_path / "sequences" / "genbank"
        fasta_dir.mkdir(parents=True)
        gb_dir.mkdir(parents=True)

        shutil.copy(FIXTURES_DIR / "sample.fasta", fasta_dir / "NC_045512.2.fasta")
        shutil.copy(FIXTURES_DIR / "sample.gb", gb_dir / "NC_045512.2.gb")

        summary = parse_dataset_accessions(
            session=db_session_with_accession,
            accession_ids=["NC_045512.2"],
            data_dir=str(tmp_path),
        )

        assert summary.total_accessions == 1
        assert summary.fasta_parsed == 1

    def test_parse_summary_to_dict(self):
        from gwico_ssr.parsers.persist import ParseSummary

        s = ParseSummary(
            total_accessions=5,
            fasta_parsed=3,
            genbank_parsed=2,
            features_inserted=10,
        )
        d = s.to_dict()
        assert d["total_accessions"] == 5
        assert d["fasta_parsed"] == 3
        assert d["features_inserted"] == 10


# ===========================================================================
# CLI PARSE COMMAND TESTS
# ===========================================================================

class TestParseCLI:

    def test_parse_help(self, tmp_path, sample_toml):
        from click.testing import CliRunner
        from gwico_ssr.cli import cli

        runner = CliRunner()
        result = runner.invoke(cli, ["--config", str(sample_toml), "parse", "--help"])
        assert result.exit_code == 0
        assert "Parse downloaded" in result.output

    def test_parse_no_dataset(self, tmp_path, sample_toml):
        from click.testing import CliRunner
        from gwico_ssr.cli import cli

        # Point DB to tmp_path so it creates fresh
        config = tmp_path / "cfg.toml"
        db_url = str(tmp_path / 'test.db').replace('\\', '/')
        out_dir = str(tmp_path).replace('\\', '/')
        config.write_text(f"""
[database]
url = "sqlite:///{db_url}"
[ncbi]
email = "test@example.com"
[ssr]
[logging]
level = "WARNING"
format = "text"
file = ""
[output]
dir = "{out_dir}"
""", encoding="utf-8")

        runner = CliRunner()
        result = runner.invoke(cli, ["--config", str(config), "parse", "nonexistent"])
        assert "Dataset not found" in result.output or "No accessions" in result.output

    def test_parse_with_accessions(self, tmp_path, sample_toml):
        """Test parse command with explicit accession list."""
        import shutil
        from click.testing import CliRunner
        from gwico_ssr.cli import cli

        db_path = tmp_path / "test.db"
        config = tmp_path / "cfg.toml"
        db_url = str(db_path).replace('\\', '/')
        out_dir = str(tmp_path).replace('\\', '/')
        config.write_text(f"""
[database]
url = "sqlite:///{db_url}"
[ncbi]
email = "test@example.com"
[ssr]
[logging]
level = "WARNING"
format = "text"
file = ""
[output]
dir = "{out_dir}"
""", encoding="utf-8")

        # Set up a dataset + accession + downloaded files
        from sqlalchemy import create_engine
        from sqlalchemy.orm import sessionmaker
        from gwico_ssr.db import create_tables
        from gwico_ssr.db.repository import get_or_create_dataset, upsert_accession, upsert_sequence_record

        engine = create_engine(f"sqlite:///{db_path}")
        create_tables(engine)
        factory = sessionmaker(bind=engine)
        sess = factory()
        ds = get_or_create_dataset(sess, "myds", "csv")
        upsert_accession(sess, accession="NC_045512.2", dataset_id=ds.dataset_id)
        upsert_sequence_record(sess, accession="NC_045512.2", download_status="downloaded")
        sess.commit()
        sess.close()
        engine.dispose()

        fasta_dir = tmp_path / "sequences" / "fasta"
        fasta_dir.mkdir(parents=True)
        shutil.copy(FIXTURES_DIR / "sample.fasta", fasta_dir / "NC_045512.2.fasta")

        runner = CliRunner()
        result = runner.invoke(cli, [
            "--config", str(config),
            "parse", "myds",
            "--accessions", "NC_045512.2",
        ])
        assert result.exit_code == 0
        assert "FASTA parsed" in result.output

    def test_parse_json_summary(self, tmp_path, sample_toml):
        """Test parse command with JSON output."""
        import shutil
        from click.testing import CliRunner
        from gwico_ssr.cli import cli

        db_path = tmp_path / "test.db"
        config = tmp_path / "cfg.toml"
        db_url = str(db_path).replace('\\', '/')
        out_dir = str(tmp_path).replace('\\', '/')
        config.write_text(f"""
[database]
url = "sqlite:///{db_url}"
[ncbi]
email = "test@example.com"
[ssr]
[logging]
level = "WARNING"
format = "text"
file = ""
[output]
dir = "{out_dir}"
""", encoding="utf-8")

        from sqlalchemy import create_engine
        from sqlalchemy.orm import sessionmaker
        from gwico_ssr.db import create_tables
        from gwico_ssr.db.repository import get_or_create_dataset, upsert_accession, upsert_sequence_record

        engine = create_engine(f"sqlite:///{db_path}")
        create_tables(engine)
        factory = sessionmaker(bind=engine)
        sess = factory()
        ds = get_or_create_dataset(sess, "myds", "csv")
        upsert_accession(sess, accession="NC_045512.2", dataset_id=ds.dataset_id)
        upsert_sequence_record(sess, accession="NC_045512.2", download_status="downloaded")
        sess.commit()
        sess.close()
        engine.dispose()

        fasta_dir = tmp_path / "sequences" / "fasta"
        fasta_dir.mkdir(parents=True)
        shutil.copy(FIXTURES_DIR / "sample.fasta", fasta_dir / "NC_045512.2.fasta")

        runner = CliRunner()
        result = runner.invoke(cli, [
            "--config", str(config),
            "parse", "myds",
            "--accessions", "NC_045512.2",
            "--json-summary",
        ])
        assert result.exit_code == 0
        data = json.loads(result.output)
        assert data["total_accessions"] == 1
        assert data["fasta_parsed"] == 1


# ===========================================================================
# PACKAGE IMPORT TESTS
# ===========================================================================

class TestParserImports:

    def test_import_parsers_package(self):
        from gwico_ssr.parsers import (
            FastaParseResult,
            GFF3ParseResult,
            GenBankParseResult,
            ParseError,
            ParseSummary,
            ParsedFeature,
            ParsedSequence,
            parse_and_persist,
            parse_dataset_accessions,
            parse_fasta,
            parse_genbank,
            parse_gff3,
        )
        # All imports should succeed
        assert parse_fasta is not None
        assert parse_genbank is not None
        assert parse_gff3 is not None
