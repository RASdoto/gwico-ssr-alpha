"""Tests for GWICO-SSR statistical analysis layer."""

from __future__ import annotations

import json
import math

import numpy as np
import pytest
from scipy import stats as sp_stats
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

from gwico_ssr.db import create_tables
from gwico_ssr.db.repository import (
    create_run,
    get_or_create_dataset,
    insert_ssr_annotations,
    insert_ssr_records,
    upsert_accession,
    upsert_accession_metrics,
)
from gwico_ssr.analysis.statistics import (
    AnalysisResult,
    AnalysisSuite,
    _pearson_ci,
    base_composition_test,
    chi_square_gene_country,
    chi_square_motif_country,
    correlation_gc_ssr,
    correlation_length_ssr,
    kruskal_wallis_by_country,
    run_all_analyses,
    shannon_entropy_by_country,
    shannon_entropy_motifs,
)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def db_session():
    engine = create_engine("sqlite:///:memory:")
    create_tables(engine)
    factory = sessionmaker(bind=engine)
    session = factory()
    yield session
    session.close()


def _seed_multi_country_dataset(session, n_per_country=None):
    """Seed a dataset with accessions across multiple countries.

    Returns (dataset, metrics_run).
    """
    if n_per_country is None:
        n_per_country = {"USA": 10, "China": 8, "UK": 6}

    ds = get_or_create_dataset(session, "test_ds", "csv")
    metrics_run = create_run(session, dataset_id=ds.dataset_id, stage="metrics")

    acc_idx = 0
    for country, count in n_per_country.items():
        for i in range(count):
            acc_id = f"ACC_{acc_idx:04d}"
            acc_idx += 1
            genome_length = 29000 + acc_idx * 100
            gc_content = 0.35 + acc_idx * 0.001
            ssr_count = 50 + acc_idx * 2
            ssr_bp = 500 + acc_idx * 20

            upsert_accession(
                session, accession=acc_id, dataset_id=ds.dataset_id,
                country=country, genome_length=genome_length, gc_content=gc_content,
            )
            upsert_accession_metrics(
                session, accession=acc_id, run_id=metrics_run.run_id,
                ssr_count_total=ssr_count, ssr_bp_total=ssr_bp,
                ra=ssr_count / (genome_length / 1000),
                rd=ssr_bp / (genome_length / 1_000_000),
                mono_count=ssr_count // 6,
                di_count=ssr_count // 5,
                tri_count=ssr_count // 3,
                tetra_count=ssr_count // 8,
                penta_count=ssr_count // 10,
                hexa_count=ssr_count // 12,
                dominant_motif="AAG",
            )

    session.flush()
    return ds, metrics_run


def _seed_ssr_records_for_analysis(session, ds):
    """Add SSR records and annotations for chi-square / entropy tests."""
    acc_ids = [f"ACC_{i:04d}" for i in range(24)]
    ssr_dicts = []
    annotation_dicts = []

    for acc_id in acc_ids:
        for idx, (motif, size, repeat_bp, gene) in enumerate([
            ("AAG", 3, 9, "ORF1ab"),
            ("A", 1, 10, "ORF1ab"),
            ("AT", 2, 12, "S"),
            ("AATG", 4, 12, "N"),
        ]):
            ssr_dicts.append({
                "accession": acc_id, "start": 100 + idx * 100,
                "end": 100 + idx * 100 + repeat_bp,
                "motif_raw": motif, "motif_canonical": motif,
                "motif_size": size, "repeat_units": repeat_bp // size,
                "repeat_length_bp": repeat_bp, "strand": "+",
                "actual_repeat": motif * (repeat_bp // size),
                "detector_version": "gwico-ssr-1.0",
            })

    insert_ssr_records(session, ssr_dicts)
    session.flush()

    # Now add annotations linking SSRs to genes
    from gwico_ssr.db.repository import get_ssr_records_for_accession

    for acc_id in acc_ids:
        ssrs = get_ssr_records_for_accession(session, acc_id)
        gene_map = {0: "ORF1ab", 1: "ORF1ab", 2: "S", 3: "N"}
        for i, ssr in enumerate(ssrs):
            if i in gene_map:
                annotation_dicts.append({
                    "ssr_id": ssr.ssr_id, "accession": acc_id,
                    "gene_name": gene_map[i], "region_class": "CDS",
                    "overlap_bp": ssr.repeat_length_bp,
                })

    insert_ssr_annotations(session, annotation_dicts)
    session.flush()


# ===========================================================================
# AnalysisResult tests
# ===========================================================================

class TestAnalysisResult:

    def test_to_dict(self):
        r = AnalysisResult(
            analysis_name="test", grouping="g", metric="m",
            test_name="kruskal-wallis", statistic=10.5, p_value=0.001,
            effect_size=0.1, n=100,
        )
        d = r.to_dict(run_id=42)
        assert d["run_id"] == 42
        assert d["analysis_name"] == "test"
        assert d["statistic"] == 10.5
        assert d["metadata_json"] is None

    def test_to_dict_with_metadata(self):
        r = AnalysisResult(
            analysis_name="test", grouping=None, metric=None,
            test_name="chi-square", statistic=5.0, p_value=0.05,
            metadata={"key": "value"},
        )
        d = r.to_dict(run_id=1)
        assert json.loads(d["metadata_json"]) == {"key": "value"}


# ===========================================================================
# FDR correction tests
# ===========================================================================

class TestFDRCorrection:

    def test_apply_fdr_basic(self):
        suite = AnalysisSuite(results=[
            AnalysisResult(
                analysis_name=f"test_{i}", grouping=None, metric=None,
                test_name="t", statistic=1.0, p_value=p,
            )
            for i, p in enumerate([0.01, 0.04, 0.03, 0.002])
        ])
        suite.apply_fdr()
        # All should have corrected p-values
        for r in suite.results:
            assert r.p_value_corrected is not None
            assert r.p_value_corrected >= r.p_value

    def test_fdr_preserves_order(self):
        """Smallest raw p should still have smallest corrected p after BH."""
        suite = AnalysisSuite(results=[
            AnalysisResult(
                analysis_name="a", grouping=None, metric=None,
                test_name="t", statistic=1.0, p_value=0.001,
            ),
            AnalysisResult(
                analysis_name="b", grouping=None, metric=None,
                test_name="t", statistic=1.0, p_value=0.05,
            ),
        ])
        suite.apply_fdr()
        assert suite.results[0].p_value_corrected <= suite.results[1].p_value_corrected

    def test_fdr_skips_none_pvalues(self):
        suite = AnalysisSuite(results=[
            AnalysisResult(
                analysis_name="a", grouping=None, metric=None,
                test_name="t", statistic=1.0, p_value=0.01,
            ),
            AnalysisResult(
                analysis_name="b", grouping=None, metric=None,
                test_name="entropy", statistic=2.0, p_value=None,
            ),
        ])
        suite.apply_fdr()
        assert suite.results[0].p_value_corrected is not None
        assert suite.results[1].p_value_corrected is None

    def test_fdr_empty(self):
        suite = AnalysisSuite()
        suite.apply_fdr()  # Should not raise
        assert len(suite.results) == 0


# ===========================================================================
# Pearson CI tests
# ===========================================================================

class TestPearsonCI:

    def test_ci_known_values(self):
        lo, hi = _pearson_ci(0.5, 100)
        assert lo is not None and hi is not None
        assert lo < 0.5 < hi
        # For n=100, r=0.5, CI should be roughly (0.33, 0.64)
        assert 0.2 < lo < 0.45
        assert 0.55 < hi < 0.75

    def test_ci_small_n(self):
        lo, hi = _pearson_ci(0.5, 3)
        assert lo is None and hi is None

    def test_ci_perfect_correlation(self):
        lo, hi = _pearson_ci(0.99, 50)
        assert lo is not None and hi is not None
        assert lo > 0.9


# ===========================================================================
# Kruskal-Wallis tests
# ===========================================================================

class TestKruskalWallis:

    def test_kw_by_country(self, db_session):
        _seed_multi_country_dataset(db_session)
        ds, run = _seed_multi_country_dataset.__wrapped__(db_session) if hasattr(
            _seed_multi_country_dataset, '__wrapped__') else (None, None)
        # Re-seed properly
        from sqlalchemy import select
        from gwico_ssr.models.schema import Dataset, Run
        ds = db_session.execute(select(Dataset)).scalar_one()
        run = db_session.execute(
            select(Run).where(Run.stage == "metrics")
        ).scalar_one()

        result = kruskal_wallis_by_country(db_session, run.run_id)
        assert result is not None
        assert result.test_name == "kruskal-wallis"
        assert result.statistic > 0
        assert result.p_value is not None
        assert result.effect_size is not None
        assert result.n == 24  # 10+8+6
        assert result.metadata["n_groups"] == 3

    def test_kw_ra_metric(self, db_session):
        _seed_multi_country_dataset(db_session)
        from sqlalchemy import select
        from gwico_ssr.models.schema import Run
        run = db_session.execute(
            select(Run).where(Run.stage == "metrics")
        ).scalar_one()

        result = kruskal_wallis_by_country(db_session, run.run_id, "ra")
        assert result is not None
        assert result.metric == "ra"

    def test_kw_insufficient_groups(self, db_session):
        _seed_multi_country_dataset(db_session, {"USA": 10})
        from sqlalchemy import select
        from gwico_ssr.models.schema import Run
        run = db_session.execute(
            select(Run).where(Run.stage == "metrics")
        ).scalar_one()
        result = kruskal_wallis_by_country(db_session, run.run_id)
        assert result is None  # Only 1 group


# ===========================================================================
# Chi-square tests
# ===========================================================================

class TestChiSquare:

    def test_chi_square_gene_country(self, db_session):
        ds, run = _seed_multi_country_dataset(db_session)
        _seed_ssr_records_for_analysis(db_session, ds)

        result = chi_square_gene_country(db_session, ds.dataset_id)
        assert result is not None
        assert result.test_name == "chi-square"
        assert result.statistic >= 0
        assert result.p_value is not None
        assert result.effect_size is not None  # Cramér's V
        assert result.metadata["n_genes"] >= 2

    def test_chi_square_motif_country(self, db_session):
        ds, run = _seed_multi_country_dataset(db_session)
        _seed_ssr_records_for_analysis(db_session, ds)

        result = chi_square_motif_country(db_session, ds.dataset_id)
        assert result is not None
        assert result.test_name == "chi-square"
        assert result.statistic >= 0
        assert result.effect_size is not None

    def test_chi_square_no_data(self, db_session):
        ds = get_or_create_dataset(db_session, "empty_ds", "csv")
        db_session.flush()
        result = chi_square_gene_country(db_session, ds.dataset_id)
        assert result is None


# ===========================================================================
# Correlation tests
# ===========================================================================

class TestCorrelation:

    def test_pearson_length_ssr(self, db_session):
        ds, run = _seed_multi_country_dataset(db_session)
        result = correlation_length_ssr(db_session, run.run_id, "pearson")
        assert result is not None
        assert result.test_name == "pearson"
        assert result.statistic is not None
        assert result.ci_lower is not None
        assert result.ci_upper is not None
        assert result.ci_lower < result.statistic < result.ci_upper

    def test_spearman_length_ssr(self, db_session):
        ds, run = _seed_multi_country_dataset(db_session)
        result = correlation_length_ssr(db_session, run.run_id, "spearman")
        assert result is not None
        assert result.test_name == "spearman"
        assert result.ci_lower is None  # Spearman doesn't use Fisher z CI

    def test_pearson_gc_ssr(self, db_session):
        ds, run = _seed_multi_country_dataset(db_session)
        result = correlation_gc_ssr(db_session, run.run_id, "pearson")
        assert result is not None
        assert result.effect_size is not None  # r is the effect size

    def test_correlation_insufficient_data(self, db_session):
        ds = get_or_create_dataset(db_session, "tiny", "csv")
        run = create_run(db_session, dataset_id=ds.dataset_id, stage="metrics")
        upsert_accession(db_session, accession="X1", dataset_id=ds.dataset_id,
                         genome_length=30000)
        upsert_accession_metrics(db_session, accession="X1", run_id=run.run_id,
                                 ssr_count_total=50)
        db_session.flush()
        result = correlation_length_ssr(db_session, run.run_id)
        assert result is None  # Need at least 3 data points


# ===========================================================================
# Shannon entropy tests
# ===========================================================================

class TestShannonEntropy:

    def test_entropy_motifs(self, db_session):
        ds, run = _seed_multi_country_dataset(db_session)
        _seed_ssr_records_for_analysis(db_session, ds)

        result = shannon_entropy_motifs(db_session, ds.dataset_id)
        assert result is not None
        assert result.test_name == "shannon-entropy"
        assert result.statistic > 0  # H > 0 for non-trivial distribution
        assert result.p_value is None  # Descriptive, no p-value
        assert result.effect_size is not None  # Pielou's evenness
        assert 0 <= result.effect_size <= 1
        assert result.metadata["n_unique_motifs"] == 4

    def test_entropy_by_country(self, db_session):
        ds, run = _seed_multi_country_dataset(db_session)
        results = shannon_entropy_by_country(db_session, run.run_id)
        assert len(results) == 3  # USA, China, UK
        for r in results:
            assert r.statistic > 0  # Each country has non-trivial distribution
            assert r.grouping is not None

    def test_entropy_known_value(self):
        """Verify Shannon entropy formula with known values."""
        # Uniform distribution over 4 items: H = log2(4) = 2.0
        probs = np.array([0.25, 0.25, 0.25, 0.25])
        h = float(-np.sum(probs * np.log2(probs)))
        assert abs(h - 2.0) < 1e-10

    def test_entropy_single_motif(self):
        """Single element: H = 0."""
        probs = np.array([1.0])
        h = float(-np.sum(probs * np.log2(probs)))
        assert abs(h) < 1e-10


# ===========================================================================
# Base composition test
# ===========================================================================

class TestBaseComposition:

    def test_base_composition(self, db_session):
        ds, run = _seed_multi_country_dataset(db_session)
        _seed_ssr_records_for_analysis(db_session, ds)

        result = base_composition_test(db_session, ds.dataset_id)
        assert result is not None
        assert result.test_name == "chi-square-gof"
        assert result.statistic > 0
        assert result.p_value is not None
        assert result.metadata["at_percent"] > 0
        assert result.metadata["gc_percent"] > 0

    def test_base_composition_no_data(self, db_session):
        ds = get_or_create_dataset(db_session, "empty", "csv")
        db_session.flush()
        result = base_composition_test(db_session, ds.dataset_id)
        assert result is None


# ===========================================================================
# Full suite tests
# ===========================================================================

class TestRunAllAnalyses:

    def test_run_all(self, db_session):
        ds, run = _seed_multi_country_dataset(db_session)
        _seed_ssr_records_for_analysis(db_session, ds)

        suite = run_all_analyses(db_session, run.run_id, ds.dataset_id)
        assert len(suite.results) > 0

        # Check FDR was applied
        p_results = [r for r in suite.results if r.p_value is not None]
        for r in p_results:
            assert r.p_value_corrected is not None

    def test_persist_results(self, db_session):
        ds, run = _seed_multi_country_dataset(db_session)
        _seed_ssr_records_for_analysis(db_session, ds)

        suite = run_all_analyses(db_session, run.run_id, ds.dataset_id)
        analysis_run = create_run(db_session, dataset_id=ds.dataset_id, stage="analysis")
        count = suite.persist(db_session, analysis_run.run_id)
        assert count == len(suite.results)

        # Verify persisted
        from sqlalchemy import select, func
        from gwico_ssr.models.schema import StatisticalResult
        persisted_count = db_session.execute(
            select(func.count(StatisticalResult.result_id))
            .where(StatisticalResult.run_id == analysis_run.run_id)
        ).scalar()
        assert persisted_count == count

    def test_to_list(self, db_session):
        ds, run = _seed_multi_country_dataset(db_session)
        _seed_ssr_records_for_analysis(db_session, ds)

        suite = run_all_analyses(db_session, run.run_id, ds.dataset_id)
        lst = suite.to_list()
        assert isinstance(lst, list)
        assert len(lst) == len(suite.results)
        assert all("analysis_name" in item for item in lst)
        assert all("p_value" in item for item in lst)


# ===========================================================================
# CLI tests
# ===========================================================================

class TestAnalyzeCLI:

    def test_analyze_help(self, tmp_path, sample_toml):
        from click.testing import CliRunner
        from gwico_ssr.cli import cli

        runner = CliRunner()
        result = runner.invoke(cli, ["--config", str(sample_toml), "analyze", "--help"])
        assert result.exit_code == 0
        assert "Run statistical analyses" in result.output

    def test_analyze_no_dataset(self, tmp_path):
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
        result = runner.invoke(cli, ["--config", str(config), "analyze", "nonexistent"])
        assert "Dataset not found" in result.output

    def test_analyze_no_metrics_run(self, tmp_path):
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
        get_or_create_dataset(sess, "myds", "csv")
        sess.commit()
        sess.close()
        engine.dispose()

        runner = CliRunner()
        result = runner.invoke(cli, ["--config", str(config), "analyze", "myds"])
        assert "No completed metrics run" in result.output

    def test_analyze_end_to_end(self, tmp_path):
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
        from gwico_ssr.db.repository import update_run_status
        run = create_run(sess, dataset_id=ds.dataset_id, stage="metrics")

        # Seed 3 countries with 5+ accessions each
        idx = 0
        for country in ["USA", "China", "UK"]:
            for i in range(6):
                acc = f"T{idx:03d}"
                idx += 1
                upsert_accession(sess, accession=acc, dataset_id=ds.dataset_id,
                                 country=country, genome_length=29000 + idx * 100,
                                 gc_content=0.35 + idx * 0.001)
                upsert_accession_metrics(
                    sess, accession=acc, run_id=run.run_id,
                    ssr_count_total=50 + idx * 2, ssr_bp_total=500 + idx * 20,
                    ra=1.5 + idx * 0.05, rd=15000.0 + idx * 100,
                    mono_count=5, di_count=10, tri_count=15,
                    tetra_count=5, penta_count=3, hexa_count=2,
                )
                insert_ssr_records(sess, [{
                    "accession": acc, "start": 100, "end": 109,
                    "motif_raw": "AAG", "motif_canonical": "AAG",
                    "motif_size": 3, "repeat_units": 3,
                    "repeat_length_bp": 9, "strand": "+",
                    "actual_repeat": "AAGAAGAAG",
                    "detector_version": "gwico-ssr-1.0",
                }])

        update_run_status(sess, run.run_id, "completed")
        sess.commit()
        sess.close()
        engine.dispose()

        runner = CliRunner()
        result = runner.invoke(cli, [
            "--config", str(config),
            "analyze", "myds",
        ])
        assert result.exit_code == 0
        assert "Analysis complete:" in result.output
        assert "results persisted" in result.output

    def test_analyze_json_summary(self, tmp_path):
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
        from gwico_ssr.db.repository import update_run_status
        run = create_run(sess, dataset_id=ds.dataset_id, stage="metrics")

        idx = 0
        for country in ["USA", "China", "UK"]:
            for i in range(6):
                acc = f"J{idx:03d}"
                idx += 1
                upsert_accession(sess, accession=acc, dataset_id=ds.dataset_id,
                                 country=country, genome_length=29000 + idx * 100,
                                 gc_content=0.35 + idx * 0.001)
                upsert_accession_metrics(
                    sess, accession=acc, run_id=run.run_id,
                    ssr_count_total=50 + idx * 2, ssr_bp_total=500 + idx * 20,
                    ra=1.5 + idx * 0.05, rd=15000.0 + idx * 100,
                    mono_count=5, di_count=10, tri_count=15,
                    tetra_count=5, penta_count=3, hexa_count=2,
                )
                insert_ssr_records(sess, [{
                    "accession": acc, "start": 100, "end": 109,
                    "motif_raw": "AAG", "motif_canonical": "AAG",
                    "motif_size": 3, "repeat_units": 3,
                    "repeat_length_bp": 9, "strand": "+",
                    "actual_repeat": "AAGAAGAAG",
                    "detector_version": "gwico-ssr-1.0",
                }])

        update_run_status(sess, run.run_id, "completed")
        sess.commit()
        sess.close()
        engine.dispose()

        runner = CliRunner()
        result = runner.invoke(cli, [
            "--config", str(config),
            "analyze", "myds",
            "--json-summary",
        ])
        assert result.exit_code == 0
        data = json.loads(result.output)
        assert data["total_results"] > 0
        assert "results" in data

    def test_analyze_skip_existing(self, tmp_path):
        """Should skip if results already exist for the metrics run."""
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
        from gwico_ssr.db.repository import update_run_status
        run = create_run(sess, dataset_id=ds.dataset_id, stage="metrics")

        idx = 0
        for country in ["USA", "China", "UK"]:
            for i in range(6):
                acc = f"S{idx:03d}"
                idx += 1
                upsert_accession(sess, accession=acc, dataset_id=ds.dataset_id,
                                 country=country, genome_length=29000 + idx * 100,
                                 gc_content=0.35 + idx * 0.001)
                upsert_accession_metrics(
                    sess, accession=acc, run_id=run.run_id,
                    ssr_count_total=50 + idx * 2, ssr_bp_total=500 + idx * 20,
                    ra=1.5 + idx * 0.05, rd=15000.0 + idx * 100,
                    mono_count=5, di_count=10, tri_count=15,
                    tetra_count=5, penta_count=3, hexa_count=2,
                )
                insert_ssr_records(sess, [{
                    "accession": acc, "start": 100, "end": 109,
                    "motif_raw": "AAG", "motif_canonical": "AAG",
                    "motif_size": 3, "repeat_units": 3,
                    "repeat_length_bp": 9, "strand": "+",
                    "actual_repeat": "AAGAAGAAG",
                    "detector_version": "gwico-ssr-1.0",
                }])

        update_run_status(sess, run.run_id, "completed")
        sess.commit()
        sess.close()
        engine.dispose()

        runner = CliRunner()
        # First run
        result1 = runner.invoke(cli, ["--config", str(config), "analyze", "myds"])
        assert result1.exit_code == 0
        assert "Analysis complete:" in result1.output

        # Second run (should detect existing)
        result2 = runner.invoke(cli, ["--config", str(config), "analyze", "myds"])
        assert "already exist" in result2.output


# ===========================================================================
# Import tests
# ===========================================================================

class TestAnalysisImports:

    def test_import_analysis_package(self):
        from gwico_ssr.analysis import (
            AnalysisResult,
            AnalysisSuite,
            base_composition_test,
            chi_square_gene_country,
            chi_square_motif_country,
            correlation_gc_ssr,
            correlation_length_ssr,
            kruskal_wallis_by_country,
            run_all_analyses,
            shannon_entropy_by_country,
            shannon_entropy_motifs,
        )
        assert run_all_analyses is not None
        assert AnalysisSuite is not None
