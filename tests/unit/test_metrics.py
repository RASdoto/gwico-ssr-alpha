"""Tests for GWICO-SSR metric computation and cohort aggregation."""

from __future__ import annotations

import json
from dataclasses import dataclass

import pytest
from sqlalchemy import create_engine, select, func
from sqlalchemy.orm import sessionmaker

from gwico_ssr.db import create_tables
from gwico_ssr.db.repository import (
    create_run,
    get_or_create_dataset,
    has_accession_metrics,
    insert_features,
    insert_ssr_annotations,
    insert_ssr_records,
    upsert_accession,
    upsert_accession_metrics,
    get_accession_metrics,
)
from gwico_ssr.metrics.calculator import (
    AccessionMetricsResult,
    compute_accession_metrics,
    compute_dominant_motif,
    compute_metrics_for_accession,
    compute_motif_size_counts,
    compute_ra,
    compute_rd,
)
from gwico_ssr.metrics.aggregator import (
    aggregate_by_country,
    aggregate_by_gene,
    aggregate_by_motif_size,
    aggregate_motif_frequencies,
    get_dataset_summary,
)
from gwico_ssr.models.schema import AccessionMetrics


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

@dataclass
class FakeSSR:
    """Lightweight SSR stand-in for formula tests (no DB needed)."""
    motif_size: int
    motif_canonical: str
    repeat_length_bp: int


@pytest.fixture
def db_session():
    engine = create_engine("sqlite:///:memory:")
    create_tables(engine)
    factory = sessionmaker(bind=engine)
    session = factory()
    yield session
    session.close()


def _seed_accession_with_ssrs(session, accession="ACC1", country="USA",
                                genome_length=30000, gc_content=0.38):
    """Create an accession with SSR records for testing."""
    ds = get_or_create_dataset(session, "test_ds", "csv")
    upsert_accession(
        session, accession=accession, dataset_id=ds.dataset_id,
        country=country, genome_length=genome_length, gc_content=gc_content,
    )
    # Insert diverse SSR records
    ssr_dicts = [
        {"accession": accession, "start": 100, "end": 112, "motif_raw": "AAG",
         "motif_canonical": "AAG", "motif_size": 3, "repeat_units": 4,
         "repeat_length_bp": 12, "strand": "+", "actual_repeat": "AAGAAGAAGAAG",
         "detector_version": "gwico-ssr-1.0"},
        {"accession": accession, "start": 200, "end": 210, "motif_raw": "A",
         "motif_canonical": "A", "motif_size": 1, "repeat_units": 10,
         "repeat_length_bp": 10, "strand": "+", "actual_repeat": "A" * 10,
         "detector_version": "gwico-ssr-1.0"},
        {"accession": accession, "start": 300, "end": 312, "motif_raw": "AT",
         "motif_canonical": "AT", "motif_size": 2, "repeat_units": 6,
         "repeat_length_bp": 12, "strand": "+", "actual_repeat": "ATATATATATAT",
         "detector_version": "gwico-ssr-1.0"},
        {"accession": accession, "start": 500, "end": 512, "motif_raw": "AAG",
         "motif_canonical": "AAG", "motif_size": 3, "repeat_units": 4,
         "repeat_length_bp": 12, "strand": "+", "actual_repeat": "AAGAAGAAGAAG",
         "detector_version": "gwico-ssr-1.0"},
        {"accession": accession, "start": 600, "end": 612, "motif_raw": "AATG",
         "motif_canonical": "AATG", "motif_size": 4, "repeat_units": 3,
         "repeat_length_bp": 12, "strand": "+", "actual_repeat": "AATGAATGAATG",
         "detector_version": "gwico-ssr-1.0"},
    ]
    insert_ssr_records(session, ssr_dicts)
    session.flush()
    return ds


# ===========================================================================
# FORMULA UNIT TESTS
# ===========================================================================

class TestRAFormula:

    def test_ra_basic(self):
        # 57 SSRs / (29903 / 1000) = 57 / 29.903 ≈ 1.9063
        result = compute_ra(57, 29903)
        assert result is not None
        assert abs(result - 57 / 29.903) < 0.001

    def test_ra_zero_genome(self):
        assert compute_ra(10, 0) is None

    def test_ra_none_genome(self):
        assert compute_ra(10, None) is None

    def test_ra_zero_ssrs(self):
        result = compute_ra(0, 30000)
        assert result == 0.0

    def test_ra_known_value(self):
        # 5 SSRs / (30000 / 1000) = 5 / 30 = 0.1667
        result = compute_ra(5, 30000)
        assert result is not None
        assert abs(result - 5 / 30.0) < 0.0001


class TestRDFormula:

    def test_rd_basic(self):
        # 500 bp / (29903 / 1_000_000) = 500 / 0.029903 ≈ 16720.7
        result = compute_rd(500, 29903)
        assert result is not None
        assert abs(result - 500 / 0.029903) < 1.0

    def test_rd_zero_genome(self):
        assert compute_rd(100, 0) is None

    def test_rd_none_genome(self):
        assert compute_rd(100, None) is None

    def test_rd_zero_bp(self):
        result = compute_rd(0, 30000)
        assert result == 0.0

    def test_rd_known_value(self):
        # 58 bp / (30000 / 1_000_000) = 58 / 0.03 = 1933.33
        result = compute_rd(58, 30000)
        assert result is not None
        assert abs(result - 58 / 0.03) < 0.1


class TestMotifSizeCounts:

    def test_counts(self):
        ssrs = [
            FakeSSR(motif_size=1, motif_canonical="A", repeat_length_bp=10),
            FakeSSR(motif_size=1, motif_canonical="A", repeat_length_bp=12),
            FakeSSR(motif_size=3, motif_canonical="AAG", repeat_length_bp=9),
            FakeSSR(motif_size=6, motif_canonical="AAGCTG", repeat_length_bp=18),
        ]
        counts = compute_motif_size_counts(ssrs)
        assert counts["mono_count"] == 2
        assert counts["tri_count"] == 1
        assert counts["hexa_count"] == 1
        assert counts["di_count"] == 0

    def test_empty(self):
        counts = compute_motif_size_counts([])
        assert all(v == 0 for v in counts.values())


class TestDominantMotif:

    def test_single_motif(self):
        ssrs = [
            FakeSSR(motif_size=3, motif_canonical="AAG", repeat_length_bp=9),
            FakeSSR(motif_size=3, motif_canonical="AAG", repeat_length_bp=12),
        ]
        assert compute_dominant_motif(ssrs) == "AAG"

    def test_multiple_motifs(self):
        ssrs = [
            FakeSSR(motif_size=3, motif_canonical="AAG", repeat_length_bp=9),
            FakeSSR(motif_size=3, motif_canonical="AAG", repeat_length_bp=12),
            FakeSSR(motif_size=1, motif_canonical="A", repeat_length_bp=10),
        ]
        assert compute_dominant_motif(ssrs) == "AAG"  # 2 vs 1

    def test_empty(self):
        assert compute_dominant_motif([]) is None


# ===========================================================================
# COMPUTE_ACCESSION_METRICS TESTS
# ===========================================================================

class TestComputeAccessionMetrics:

    def test_basic_metrics(self):
        ssrs = [
            FakeSSR(motif_size=1, motif_canonical="A", repeat_length_bp=10),
            FakeSSR(motif_size=3, motif_canonical="AAG", repeat_length_bp=12),
            FakeSSR(motif_size=3, motif_canonical="AAG", repeat_length_bp=9),
            FakeSSR(motif_size=2, motif_canonical="AT", repeat_length_bp=12),
            FakeSSR(motif_size=4, motif_canonical="AATG", repeat_length_bp=12),
        ]
        result = compute_accession_metrics("ACC1", ssrs, genome_length=30000, gc_content=0.38)
        assert result.accession == "ACC1"
        assert result.ssr_count_total == 5
        assert result.ssr_bp_total == 55
        assert result.mono_count == 1
        assert result.di_count == 1
        assert result.tri_count == 2
        assert result.tetra_count == 1
        assert result.penta_count == 0
        assert result.hexa_count == 0
        assert result.dominant_motif == "AAG"
        assert result.ra is not None
        assert abs(result.ra - 5 / 30.0) < 0.001
        assert result.rd is not None
        assert abs(result.rd - 55 / 0.03) < 0.1
        assert result.genome_length == 30000
        assert result.gc_content == 0.38

    def test_no_ssrs(self):
        result = compute_accession_metrics("ACC1", [], genome_length=30000)
        assert result.ssr_count_total == 0
        assert result.ssr_bp_total == 0
        assert result.ra == 0.0
        assert result.rd == 0.0
        assert result.dominant_motif is None

    def test_no_genome_length(self):
        ssrs = [FakeSSR(motif_size=1, motif_canonical="A", repeat_length_bp=10)]
        result = compute_accession_metrics("ACC1", ssrs, genome_length=None)
        assert result.ra is None
        assert result.rd is None

    def test_to_dict(self):
        ssrs = [FakeSSR(motif_size=3, motif_canonical="AAG", repeat_length_bp=9)]
        result = compute_accession_metrics("ACC1", ssrs, genome_length=30000)
        d = result.to_dict(run_id=1)
        assert d["accession"] == "ACC1"
        assert d["run_id"] == 1
        assert "ssr_count_total" in d
        assert "ra" in d
        assert "dominant_motif" in d


# ===========================================================================
# DB INTEGRATION TESTS
# ===========================================================================

class TestComputeMetricsFromDB:

    def test_compute_metrics_for_accession(self, db_session):
        _seed_accession_with_ssrs(db_session)
        result = compute_metrics_for_accession(db_session, "ACC1")
        assert result.accession == "ACC1"
        assert result.ssr_count_total == 5
        assert result.ssr_bp_total == 58  # 12+10+12+12+12
        assert result.mono_count == 1
        assert result.di_count == 1
        assert result.tri_count == 2
        assert result.tetra_count == 1
        assert result.dominant_motif == "AAG"
        assert result.ra is not None
        assert result.rd is not None

    def test_persist_and_retrieve_metrics(self, db_session):
        ds = _seed_accession_with_ssrs(db_session)
        run = create_run(db_session, dataset_id=ds.dataset_id, stage="metrics")
        result = compute_metrics_for_accession(db_session, "ACC1")
        upsert_accession_metrics(db_session, **result.to_dict(run.run_id))

        metrics = get_accession_metrics(db_session, "ACC1", run.run_id)
        assert len(metrics) == 1
        m = metrics[0]
        assert m.ssr_count_total == 5
        assert m.dominant_motif == "AAG"
        assert m.ra is not None

    def test_has_accession_metrics(self, db_session):
        ds = _seed_accession_with_ssrs(db_session)
        run = create_run(db_session, dataset_id=ds.dataset_id, stage="metrics")
        assert has_accession_metrics(db_session, "ACC1", run.run_id) is False
        result = compute_metrics_for_accession(db_session, "ACC1")
        upsert_accession_metrics(db_session, **result.to_dict(run.run_id))
        assert has_accession_metrics(db_session, "ACC1", run.run_id) is True

    def test_upsert_replaces_metrics(self, db_session):
        ds = _seed_accession_with_ssrs(db_session)
        run = create_run(db_session, dataset_id=ds.dataset_id, stage="metrics")
        upsert_accession_metrics(
            db_session, accession="ACC1", run_id=run.run_id,
            ssr_count_total=10, ssr_bp_total=100,
        )
        upsert_accession_metrics(
            db_session, accession="ACC1", run_id=run.run_id,
            ssr_count_total=20,
        )
        metrics = get_accession_metrics(db_session, "ACC1", run.run_id)
        assert len(metrics) == 1
        assert metrics[0].ssr_count_total == 20

    def test_no_accession(self, db_session):
        result = compute_metrics_for_accession(db_session, "NONEXISTENT")
        assert result.ssr_count_total == 0
        assert result.ra is None


# ===========================================================================
# COHORT AGGREGATION TESTS
# ===========================================================================

class TestCohortAggregations:

    def _seed_multi_country(self, session):
        ds = get_or_create_dataset(session, "test_ds", "csv")
        run = create_run(session, dataset_id=ds.dataset_id, stage="metrics")

        for acc_id, country, ssr_count, ssr_bp in [
            ("A1", "USA", 50, 500),
            ("A2", "USA", 60, 600),
            ("A3", "China", 40, 400),
        ]:
            upsert_accession(session, accession=acc_id, dataset_id=ds.dataset_id,
                             country=country, genome_length=30000)
            upsert_accession_metrics(
                session, accession=acc_id, run_id=run.run_id,
                ssr_count_total=ssr_count, ssr_bp_total=ssr_bp,
                ra=ssr_count / 30.0, rd=ssr_bp / 0.03,
                mono_count=5, di_count=10, tri_count=20,
                tetra_count=5, penta_count=5, hexa_count=5,
            )
        session.flush()
        return ds, run

    def test_aggregate_by_country(self, db_session):
        _, run = self._seed_multi_country(db_session)
        results = aggregate_by_country(db_session, run.run_id)
        assert len(results) == 2
        usa = [r for r in results if r.group_value == "USA"]
        assert len(usa) == 1
        assert usa[0].count == 2
        assert usa[0].total_ssrs == 110

    def test_aggregate_by_motif_size(self, db_session):
        _, run = self._seed_multi_country(db_session)
        results = aggregate_by_motif_size(db_session, run.run_id)
        assert len(results) == 6
        mono = [r for r in results if r.group_value == "mono"]
        assert len(mono) == 1
        assert mono[0].total_ssrs == 15  # 5*3

    def test_aggregate_motif_frequencies(self, db_session):
        ds = _seed_accession_with_ssrs(db_session)
        results = aggregate_motif_frequencies(db_session, dataset_id=ds.dataset_id)
        assert len(results) > 0
        # AAG should be the top motif (2 occurrences)
        assert results[0].motif_canonical == "AAG"
        assert results[0].total_count == 2

    def test_aggregate_by_gene(self, db_session):
        ds = _seed_accession_with_ssrs(db_session)
        # Need to add annotations for gene aggregation
        from gwico_ssr.db.repository import get_ssr_records_for_accession
        ssrs = get_ssr_records_for_accession(db_session, "ACC1")
        insert_ssr_annotations(db_session, [
            {"ssr_id": ssrs[0].ssr_id, "accession": "ACC1",
             "gene_name": "ORF1ab", "region_class": "CDS", "overlap_bp": 12},
            {"ssr_id": ssrs[1].ssr_id, "accession": "ACC1",
             "gene_name": "ORF1ab", "region_class": "CDS", "overlap_bp": 10},
            {"ssr_id": ssrs[2].ssr_id, "accession": "ACC1",
             "gene_name": "S", "region_class": "CDS", "overlap_bp": 12},
        ])
        db_session.flush()

        results = aggregate_by_gene(db_session, dataset_id=ds.dataset_id)
        assert len(results) == 2
        assert results[0].gene_name == "ORF1ab"
        assert results[0].ssr_count == 2

    def test_get_dataset_summary(self, db_session):
        _, run = self._seed_multi_country(db_session)
        summary = get_dataset_summary(db_session, run.run_id)
        assert summary["total_accessions"] == 3
        assert summary["total_ssrs"] == 150  # 50+60+40
        assert summary["mean_ra"] is not None
        assert len(summary["by_motif_size"]) == 6
        assert len(summary["top_countries"]) == 2


# ===========================================================================
# CLI METRICS COMMAND TESTS
# ===========================================================================

class TestMetricsCLI:

    def test_metrics_help(self, tmp_path, sample_toml):
        from click.testing import CliRunner
        from gwico_ssr.cli import cli

        runner = CliRunner()
        result = runner.invoke(cli, ["--config", str(sample_toml), "metrics", "--help"])
        assert result.exit_code == 0
        assert "Compute per-accession SSR metrics" in result.output

    def test_metrics_no_dataset(self, tmp_path):
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
        result = runner.invoke(cli, ["--config", str(config), "metrics", "nonexistent"])
        assert "Dataset not found" in result.output

    def test_metrics_end_to_end(self, tmp_path):
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
        upsert_accession(sess, accession="TEST1", dataset_id=ds.dataset_id,
                         genome_length=30000, gc_content=0.38)
        insert_ssr_records(sess, [
            {"accession": "TEST1", "start": 100, "end": 112, "motif_raw": "AAG",
             "motif_canonical": "AAG", "motif_size": 3, "repeat_units": 4,
             "repeat_length_bp": 12, "strand": "+", "actual_repeat": "AAGAAGAAGAAG",
             "detector_version": "gwico-ssr-1.0"},
            {"accession": "TEST1", "start": 200, "end": 210, "motif_raw": "A",
             "motif_canonical": "A", "motif_size": 1, "repeat_units": 10,
             "repeat_length_bp": 10, "strand": "+", "actual_repeat": "A" * 10,
             "detector_version": "gwico-ssr-1.0"},
        ])
        sess.commit()
        sess.close()
        engine.dispose()

        runner = CliRunner()
        result = runner.invoke(cli, [
            "--config", str(config),
            "metrics", "myds",
            "--accessions", "TEST1",
        ])
        assert result.exit_code == 0
        assert "Computed: 1" in result.output
        assert "Total SSRs across dataset:" in result.output

    def test_metrics_json_summary(self, tmp_path):
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
        upsert_accession(sess, accession="TEST1", dataset_id=ds.dataset_id,
                         genome_length=30000)
        insert_ssr_records(sess, [
            {"accession": "TEST1", "start": 100, "end": 112, "motif_raw": "AAG",
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
            "metrics", "myds",
            "--accessions", "TEST1",
            "--json-summary",
        ])
        assert result.exit_code == 0
        data = json.loads(result.output)
        assert data["computed"] == 1
        assert data["total_ssrs"] == 1


# ===========================================================================
# IMPORT TESTS
# ===========================================================================

class TestMetricsImports:

    def test_import_metrics_package(self):
        from gwico_ssr.metrics import (
            AccessionMetricsResult,
            CohortSummary,
            GeneSSRSummary,
            MotifFrequency,
            aggregate_by_country,
            aggregate_by_gene,
            aggregate_by_motif_size,
            aggregate_motif_frequencies,
            compute_accession_metrics,
            compute_dominant_motif,
            compute_metrics_for_accession,
            compute_motif_size_counts,
            compute_ra,
            compute_rd,
            get_dataset_summary,
        )
        assert compute_ra is not None
        assert aggregate_by_country is not None
