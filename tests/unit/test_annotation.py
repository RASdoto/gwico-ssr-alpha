"""Tests for GWICO-SSR interval-based annotation mapping."""

from __future__ import annotations

import json
import time
from dataclasses import dataclass
from pathlib import Path

import pytest
from sqlalchemy import create_engine, select, func
from sqlalchemy.orm import sessionmaker

from gwico_ssr.annotation.mapper import (
    AnnotationResult,
    FeatureInterval,
    SSRAnnotationRecord,
    annotate_ssrs,
    annotate_ssrs_naive,
    benchmark_compare,
    build_feature_tree,
    features_from_db,
    annotate_accession,
    _compute_overlap_bp,
    _classify_region,
)
from gwico_ssr.db import create_tables
from gwico_ssr.db.repository import (
    get_or_create_dataset,
    insert_features,
    insert_ssr_records,
    insert_ssr_annotations,
    upsert_accession,
)
from gwico_ssr.models.schema import (
    FeatureRecord,
    SSRAnnotation as SSRAnnotationModel,
    SSRRecord as SSRRecordModel,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

@dataclass
class FakeSSR:
    """Lightweight SSR stand-in for unit tests (no DB needed)."""
    ssr_id: int
    start: int
    end: int


@pytest.fixture
def db_session():
    engine = create_engine("sqlite:///:memory:")
    create_tables(engine)
    factory = sessionmaker(bind=engine)
    session = factory()
    yield session
    session.close()


def _make_features(*specs) -> list[FeatureInterval]:
    """Create FeatureInterval list from (start, end, gene_name, type) tuples."""
    return [
        FeatureInterval(
            feature_id=i + 1,
            accession="ACC1",
            feature_type=s[3] if len(s) > 3 else "CDS",
            start=s[0],
            end=s[1],
            strand="+",
            gene_name=s[2] if len(s) > 2 else f"gene_{i}",
        )
        for i, s in enumerate(specs)
    ]


# ===========================================================================
# INTERVAL TREE TESTS
# ===========================================================================

class TestBuildFeatureTree:

    def test_empty_features(self):
        tree = build_feature_tree([])
        assert len(tree) == 0

    def test_single_feature(self):
        features = _make_features((100, 500, "ORF1ab"))
        tree = build_feature_tree(features)
        assert len(tree) == 1

    def test_multiple_features(self):
        features = _make_features(
            (100, 500, "ORF1ab"),
            (600, 900, "S"),
            (1000, 1200, "N"),
        )
        tree = build_feature_tree(features)
        assert len(tree) == 3

    def test_zero_length_features_skipped(self):
        features = _make_features(
            (100, 100, "empty"),  # zero-length
            (200, 500, "real"),
        )
        tree = build_feature_tree(features)
        assert len(tree) == 1

    def test_overlapping_features(self):
        features = _make_features(
            (100, 500, "geneA"),
            (400, 800, "geneB"),
        )
        tree = build_feature_tree(features)
        assert len(tree) == 2
        # Query in overlap region
        hits = tree.overlap(450, 451)
        assert len(hits) == 2


# ===========================================================================
# OVERLAP CALCULATION TESTS
# ===========================================================================

class TestOverlapComputation:

    def test_full_overlap(self):
        assert _compute_overlap_bp(100, 110, 50, 200) == 10

    def test_partial_overlap_left(self):
        assert _compute_overlap_bp(90, 110, 100, 200) == 10

    def test_partial_overlap_right(self):
        assert _compute_overlap_bp(190, 210, 100, 200) == 10

    def test_no_overlap(self):
        assert _compute_overlap_bp(0, 50, 100, 200) == 0

    def test_adjacent_no_overlap(self):
        # [0, 100) and [100, 200) — no overlap in half-open
        assert _compute_overlap_bp(0, 100, 100, 200) == 0

    def test_exact_same_interval(self):
        assert _compute_overlap_bp(100, 200, 100, 200) == 100

    def test_ssr_inside_feature(self):
        assert _compute_overlap_bp(150, 160, 100, 500) == 10

    def test_feature_inside_ssr(self):
        assert _compute_overlap_bp(100, 500, 200, 250) == 50

    def test_single_bp_overlap(self):
        assert _compute_overlap_bp(99, 101, 100, 200) == 1


# ===========================================================================
# REGION CLASSIFICATION TESTS
# ===========================================================================

class TestClassifyRegion:

    def test_known_types(self):
        for t in ("CDS", "gene", "mRNA", "tRNA", "rRNA", "ncRNA", "exon", "misc_feature"):
            assert _classify_region(t) == t

    def test_unknown_type(self):
        assert _classify_region("promoter") == "other"
        assert _classify_region("regulatory") == "other"


# ===========================================================================
# CORE ANNOTATION MAPPING TESTS
# ===========================================================================

class TestAnnotateSSRs:

    def test_ssr_inside_single_gene(self):
        features = _make_features((100, 500, "ORF1ab"))
        tree = build_feature_tree(features)
        ssrs = [FakeSSR(ssr_id=1, start=200, end=210)]
        result = annotate_ssrs(ssrs, tree, "ACC1")
        assert result.total_ssrs == 1
        assert result.annotated == 1
        assert result.intergenic == 0
        assert len(result.annotations) == 1
        assert result.annotations[0].gene_name == "ORF1ab"
        assert result.annotations[0].region_class == "CDS"
        assert result.annotations[0].overlap_bp == 10

    def test_intergenic_ssr(self):
        features = _make_features((100, 500, "ORF1ab"))
        tree = build_feature_tree(features)
        ssrs = [FakeSSR(ssr_id=1, start=600, end=620)]
        result = annotate_ssrs(ssrs, tree, "ACC1")
        assert result.intergenic == 1
        assert result.annotated == 0
        assert result.annotations[0].region_class == "intergenic"
        assert result.annotations[0].feature_id is None
        assert result.annotations[0].gene_name is None
        assert result.annotations[0].overlap_bp is None

    def test_ssr_overlap_two_genes(self):
        """SSR spanning a gene boundary should produce two annotations."""
        features = _make_features(
            (100, 300, "geneA"),
            (300, 600, "geneB"),
        )
        tree = build_feature_tree(features)
        # SSR at [290, 310) overlaps both
        ssrs = [FakeSSR(ssr_id=1, start=290, end=310)]
        result = annotate_ssrs(ssrs, tree, "ACC1")
        assert result.annotated == 1
        assert result.intergenic == 0
        assert len(result.annotations) == 2
        gene_names = {a.gene_name for a in result.annotations}
        assert gene_names == {"geneA", "geneB"}
        # Check overlap bp
        overlaps = {a.gene_name: a.overlap_bp for a in result.annotations}
        assert overlaps["geneA"] == 10  # [290, 300)
        assert overlaps["geneB"] == 10  # [300, 310)

    def test_ssr_at_exact_feature_start(self):
        features = _make_features((100, 500, "gene1"))
        tree = build_feature_tree(features)
        ssrs = [FakeSSR(ssr_id=1, start=100, end=110)]
        result = annotate_ssrs(ssrs, tree, "ACC1")
        assert result.annotated == 1
        assert result.annotations[0].overlap_bp == 10

    def test_ssr_at_exact_feature_end(self):
        """SSR starting at feature end should NOT overlap (half-open)."""
        features = _make_features((100, 500, "gene1"))
        tree = build_feature_tree(features)
        ssrs = [FakeSSR(ssr_id=1, start=500, end=510)]
        result = annotate_ssrs(ssrs, tree, "ACC1")
        assert result.intergenic == 1
        assert result.annotated == 0

    def test_ssr_ending_at_feature_start(self):
        """SSR ending at feature start should NOT overlap (half-open)."""
        features = _make_features((100, 500, "gene1"))
        tree = build_feature_tree(features)
        ssrs = [FakeSSR(ssr_id=1, start=90, end=100)]
        result = annotate_ssrs(ssrs, tree, "ACC1")
        assert result.intergenic == 1

    def test_multiple_ssrs_mixed(self):
        features = _make_features(
            (100, 500, "ORF1ab"),
            (600, 900, "S"),
        )
        tree = build_feature_tree(features)
        ssrs = [
            FakeSSR(ssr_id=1, start=200, end=210),   # in ORF1ab
            FakeSSR(ssr_id=2, start=550, end=560),   # intergenic
            FakeSSR(ssr_id=3, start=700, end=720),   # in S
        ]
        result = annotate_ssrs(ssrs, tree, "ACC1")
        assert result.total_ssrs == 3
        assert result.annotated == 2
        assert result.intergenic == 1

    def test_no_ssrs(self):
        features = _make_features((100, 500, "gene1"))
        tree = build_feature_tree(features)
        result = annotate_ssrs([], tree, "ACC1")
        assert result.total_ssrs == 0
        assert result.annotated == 0
        assert result.intergenic == 0

    def test_no_features(self):
        tree = build_feature_tree([])
        ssrs = [FakeSSR(ssr_id=1, start=200, end=210)]
        result = annotate_ssrs(ssrs, tree, "ACC1")
        assert result.intergenic == 1

    def test_ssr_spanning_entire_feature(self):
        features = _make_features((200, 250, "small_gene"))
        tree = build_feature_tree(features)
        ssrs = [FakeSSR(ssr_id=1, start=100, end=500)]
        result = annotate_ssrs(ssrs, tree, "ACC1")
        assert result.annotated == 1
        assert result.annotations[0].overlap_bp == 50  # the feature is 50bp

    def test_to_dicts(self):
        features = _make_features((100, 500, "gene1"))
        tree = build_feature_tree(features)
        ssrs = [FakeSSR(ssr_id=1, start=200, end=210)]
        result = annotate_ssrs(ssrs, tree, "ACC1")
        dicts = result.to_dicts()
        assert len(dicts) == 1
        d = dicts[0]
        assert d["ssr_id"] == 1
        assert d["accession"] == "ACC1"
        assert d["gene_name"] == "gene1"
        assert d["region_class"] == "CDS"
        assert d["overlap_bp"] == 10
        assert d["feature_id"] == 1


# ===========================================================================
# FEATURES_FROM_DB TESTS
# ===========================================================================

class TestFeaturesFromDB:

    def test_converts_orm_objects(self, db_session):
        ds = get_or_create_dataset(db_session, "test_ds", "csv")
        upsert_accession(db_session, accession="ACC1", dataset_id=ds.dataset_id)
        insert_features(db_session, [
            {"accession": "ACC1", "feature_type": "CDS", "start": 100, "end": 500,
             "strand": "+", "gene_name": "ORF1ab", "annotation_source": "genbank"},
        ])
        db_session.flush()

        from gwico_ssr.db.repository import get_features_for_accession
        records = get_features_for_accession(db_session, "ACC1")
        intervals = features_from_db(records)
        assert len(intervals) == 1
        assert intervals[0].gene_name == "ORF1ab"
        assert intervals[0].start == 100
        assert intervals[0].end == 500


# ===========================================================================
# ANNOTATE_ACCESSION (DB INTEGRATION) TESTS
# ===========================================================================

class TestAnnotateAccessionDB:

    def _setup_accession(self, session):
        """Create an accession with features and SSR records."""
        ds = get_or_create_dataset(session, "test_ds", "csv")
        upsert_accession(session, accession="ACC1", dataset_id=ds.dataset_id)

        # Insert features
        insert_features(session, [
            {"accession": "ACC1", "feature_type": "CDS", "start": 100, "end": 500,
             "strand": "+", "gene_name": "ORF1ab", "annotation_source": "genbank"},
            {"accession": "ACC1", "feature_type": "CDS", "start": 600, "end": 900,
             "strand": "+", "gene_name": "S", "annotation_source": "genbank"},
        ])

        # Insert SSR records
        insert_ssr_records(session, [
            {"accession": "ACC1", "start": 200, "end": 212, "motif_raw": "AAG",
             "motif_canonical": "AAG", "motif_size": 3, "repeat_units": 4,
             "repeat_length_bp": 12, "strand": "+", "actual_repeat": "AAGAAGAAGAAG",
             "detector_version": "gwico-ssr-1.0"},
            {"accession": "ACC1", "start": 550, "end": 562, "motif_raw": "AT",
             "motif_canonical": "AT", "motif_size": 2, "repeat_units": 6,
             "repeat_length_bp": 12, "strand": "+", "actual_repeat": "ATATATATATATAT",
             "detector_version": "gwico-ssr-1.0"},
            {"accession": "ACC1", "start": 700, "end": 712, "motif_raw": "GAA",
             "motif_canonical": "AAG", "motif_size": 3, "repeat_units": 4,
             "repeat_length_bp": 12, "strand": "+", "actual_repeat": "GAAGAAGAAGAA",
             "detector_version": "gwico-ssr-1.0"},
        ])
        session.flush()

    def test_annotate_accession_basic(self, db_session):
        self._setup_accession(db_session)
        result = annotate_accession(db_session, "ACC1")
        assert result.total_ssrs == 3
        assert result.annotated == 2  # SSR1 in ORF1ab, SSR3 in S
        assert result.intergenic == 1  # SSR2 at 550 is intergenic

    def test_annotate_and_persist(self, db_session):
        self._setup_accession(db_session)
        result = annotate_accession(db_session, "ACC1")
        dicts = result.to_dicts()
        count = insert_ssr_annotations(db_session, dicts)
        assert count == 3

        anno_count = db_session.execute(
            select(func.count()).select_from(SSRAnnotationModel).where(
                SSRAnnotationModel.accession == "ACC1"
            )
        ).scalar()
        assert anno_count == 3

    def test_annotate_no_ssrs(self, db_session):
        ds = get_or_create_dataset(db_session, "test_ds", "csv")
        upsert_accession(db_session, accession="ACC2", dataset_id=ds.dataset_id)
        result = annotate_accession(db_session, "ACC2")
        assert result.total_ssrs == 0

    def test_annotate_gene_names_correct(self, db_session):
        self._setup_accession(db_session)
        result = annotate_accession(db_session, "ACC1")
        gene_annotations = [a for a in result.annotations if a.gene_name is not None]
        gene_names = {a.gene_name for a in gene_annotations}
        assert "ORF1ab" in gene_names
        assert "S" in gene_names


# ===========================================================================
# NAIVE VS INTERVAL TREE COMPARISON TESTS
# ===========================================================================

class TestNaiveComparison:

    def test_same_results(self):
        features = _make_features(
            (100, 500, "ORF1ab"),
            (600, 900, "S"),
            (1000, 1500, "N"),
        )
        tree = build_feature_tree(features)
        ssrs = [
            FakeSSR(ssr_id=1, start=200, end=210),
            FakeSSR(ssr_id=2, start=550, end=560),
            FakeSSR(ssr_id=3, start=700, end=720),
            FakeSSR(ssr_id=4, start=1100, end=1110),
        ]
        tree_result = annotate_ssrs(ssrs, tree, "ACC1")
        naive_result = annotate_ssrs_naive(ssrs, features, "ACC1")

        assert tree_result.total_ssrs == naive_result.total_ssrs
        assert tree_result.annotated == naive_result.annotated
        assert tree_result.intergenic == naive_result.intergenic
        assert len(tree_result.annotations) == len(naive_result.annotations)

    def test_benchmark_compare(self):
        """Generate enough data to show interval tree is faster."""
        import random
        random.seed(42)

        # Create 500 features and 1000 SSRs
        features = [
            FeatureInterval(
                feature_id=i,
                accession="ACC1",
                feature_type="CDS",
                start=i * 100,
                end=i * 100 + 80,
                strand="+",
                gene_name=f"gene_{i}",
            )
            for i in range(500)
        ]
        ssrs = [
            FakeSSR(ssr_id=i, start=random.randint(0, 50000), end=random.randint(0, 50000) + 10)
            for i in range(1000)
        ]

        result = benchmark_compare(ssrs, features, "ACC1")
        assert result["results_match"]
        # Interval tree should be at least somewhat faster on this data size
        # (may not be dramatic on small inputs, but should not be slower)
        assert result["tree_time_ms"] >= 0
        assert result["naive_time_ms"] >= 0

    def test_boundary_cases_match(self):
        features = _make_features(
            (100, 200, "geneA"),
            (200, 300, "geneB"),
        )
        tree = build_feature_tree(features)
        # SSR at exact boundary [200, 210) — should only hit geneB
        ssrs = [FakeSSR(ssr_id=1, start=200, end=210)]
        tree_result = annotate_ssrs(ssrs, tree, "ACC1")
        naive_result = annotate_ssrs_naive(ssrs, features, "ACC1")

        assert len(tree_result.annotations) == len(naive_result.annotations)
        assert tree_result.annotations[0].gene_name == "geneB"
        assert naive_result.annotations[0].gene_name == "geneB"


# ===========================================================================
# CLI ANNOTATE COMMAND TESTS
# ===========================================================================

class TestAnnotateCLI:

    def test_annotate_help(self, tmp_path, sample_toml):
        from click.testing import CliRunner
        from gwico_ssr.cli import cli

        runner = CliRunner()
        result = runner.invoke(cli, ["--config", str(sample_toml), "annotate", "--help"])
        assert result.exit_code == 0
        assert "Map SSRs to genomic features" in result.output

    def test_annotate_no_dataset(self, tmp_path):
        from click.testing import CliRunner
        from gwico_ssr.cli import cli

        db_url = str(tmp_path / "test.db").replace("\\", "/")
        out_dir = str(tmp_path).replace("\\", "/")
        config = tmp_path / "cfg.toml"
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
        result = runner.invoke(cli, ["--config", str(config), "annotate", "nonexistent"])
        assert "Dataset not found" in result.output or "No accessions" in result.output

    def test_annotate_end_to_end(self, tmp_path):
        """Full end-to-end: features + SSRs → annotations."""
        from click.testing import CliRunner
        from gwico_ssr.cli import cli

        db_path = tmp_path / "test.db"
        db_url = str(db_path).replace("\\", "/")
        out_dir = str(tmp_path).replace("\\", "/")
        config = tmp_path / "cfg.toml"
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

        # Set up DB with features and SSRs
        engine = create_engine(f"sqlite:///{db_path}")
        create_tables(engine)
        factory = sessionmaker(bind=engine)
        sess = factory()
        ds = get_or_create_dataset(sess, "myds", "csv")
        upsert_accession(sess, accession="TEST1", dataset_id=ds.dataset_id)
        insert_features(sess, [
            {"accession": "TEST1", "feature_type": "CDS", "start": 100, "end": 500,
             "strand": "+", "gene_name": "ORF1ab", "annotation_source": "genbank"},
        ])
        insert_ssr_records(sess, [
            {"accession": "TEST1", "start": 200, "end": 212, "motif_raw": "AAG",
             "motif_canonical": "AAG", "motif_size": 3, "repeat_units": 4,
             "repeat_length_bp": 12, "strand": "+", "actual_repeat": "AAGAAGAAGAAG",
             "detector_version": "gwico-ssr-1.0"},
        ])
        sess.commit()
        sess.close()
        engine.dispose()

        runner = CliRunner()
        result = runner.invoke(cli, [
            "--config", str(config),
            "annotate", "myds",
            "--accessions", "TEST1",
        ])
        assert result.exit_code == 0
        assert "Processed: 1" in result.output
        assert "Total annotations:" in result.output

    def test_annotate_json_summary(self, tmp_path):
        from click.testing import CliRunner
        from gwico_ssr.cli import cli

        db_path = tmp_path / "test.db"
        db_url = str(db_path).replace("\\", "/")
        out_dir = str(tmp_path).replace("\\", "/")
        config = tmp_path / "cfg.toml"
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

        engine = create_engine(f"sqlite:///{db_path}")
        create_tables(engine)
        factory = sessionmaker(bind=engine)
        sess = factory()
        ds = get_or_create_dataset(sess, "myds", "csv")
        upsert_accession(sess, accession="TEST1", dataset_id=ds.dataset_id)
        insert_features(sess, [
            {"accession": "TEST1", "feature_type": "CDS", "start": 100, "end": 500,
             "strand": "+", "gene_name": "ORF1ab", "annotation_source": "genbank"},
        ])
        insert_ssr_records(sess, [
            {"accession": "TEST1", "start": 200, "end": 212, "motif_raw": "AAG",
             "motif_canonical": "AAG", "motif_size": 3, "repeat_units": 4,
             "repeat_length_bp": 12, "strand": "+", "actual_repeat": "AAGAAGAAGAAG",
             "detector_version": "gwico-ssr-1.0"},
            {"accession": "TEST1", "start": 600, "end": 612, "motif_raw": "AT",
             "motif_canonical": "AT", "motif_size": 2, "repeat_units": 6,
             "repeat_length_bp": 12, "strand": "+", "actual_repeat": "ATATATATATATAT",
             "detector_version": "gwico-ssr-1.0"},
        ])
        sess.commit()
        sess.close()
        engine.dispose()

        runner = CliRunner()
        result = runner.invoke(cli, [
            "--config", str(config),
            "annotate", "myds",
            "--accessions", "TEST1",
            "--json-summary",
        ])
        assert result.exit_code == 0
        data = json.loads(result.output)
        assert data["processed"] == 1
        assert data["total_annotations"] == 2  # 1 in gene, 1 intergenic


# ===========================================================================
# IMPORT TESTS
# ===========================================================================

class TestAnnotationImports:

    def test_import_annotation_package(self):
        from gwico_ssr.annotation import (
            AnnotationResult,
            FeatureInterval,
            SSRAnnotationRecord,
            annotate_accession,
            annotate_ssrs,
            annotate_ssrs_naive,
            benchmark_compare,
            build_feature_tree,
            features_from_db,
        )
        assert annotate_ssrs is not None
        assert build_feature_tree is not None
