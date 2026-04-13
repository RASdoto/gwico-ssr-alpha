"""Tests for GWICO-SSR visualization layer."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from gwico_ssr.visualization.figures import (
    generate_all_figures,
    plot_choropleth,
    plot_correlation_scatter,
    plot_country_boxplot,
    plot_country_motif_heatmap,
    plot_gene_distribution,
    plot_motif_size_distribution,
    plot_motif_sunburst,
    plot_stats_summary_table,
    plot_top_motifs,
)


# ---------------------------------------------------------------------------
# Sample data factories
# ---------------------------------------------------------------------------

def _motif_size_data():
    return [
        {"motif_size": "mono", "total_ssrs": 1200},
        {"motif_size": "di", "total_ssrs": 3500},
        {"motif_size": "tri", "total_ssrs": 5000},
        {"motif_size": "tetra", "total_ssrs": 800},
        {"motif_size": "penta", "total_ssrs": 400},
        {"motif_size": "hexa", "total_ssrs": 200},
    ]


def _motif_freq_data():
    return [
        {"motif_canonical": "A", "motif_size": 1, "total_count": 500, "total_bp": 5000},
        {"motif_canonical": "T", "motif_size": 1, "total_count": 480, "total_bp": 4800},
        {"motif_canonical": "AT", "motif_size": 2, "total_count": 1200, "total_bp": 12000},
        {"motif_canonical": "AG", "motif_size": 2, "total_count": 800, "total_bp": 8000},
        {"motif_canonical": "AAG", "motif_size": 3, "total_count": 2000, "total_bp": 18000},
        {"motif_canonical": "AAT", "motif_size": 3, "total_count": 1500, "total_bp": 13500},
        {"motif_canonical": "AATG", "motif_size": 4, "total_count": 300, "total_bp": 3600},
        {"motif_canonical": "AAAG", "motif_size": 4, "total_count": 250, "total_bp": 3000},
    ]


def _gene_data():
    return [
        {"gene_name": "ORF1ab", "ssr_count": 3000, "accession_count": 100},
        {"gene_name": "S", "ssr_count": 1500, "accession_count": 100},
        {"gene_name": "N", "ssr_count": 800, "accession_count": 80},
        {"gene_name": "M", "ssr_count": 400, "accession_count": 60},
        {"gene_name": "E", "ssr_count": 200, "accession_count": 50},
    ]


def _country_summary():
    return [
        {"country": "USA", "total_ssrs": 10000, "count": 50, "mean_ra": 12.5, "mean_rd": 150.0},
        {"country": "China", "total_ssrs": 8000, "count": 40, "mean_ra": 11.0, "mean_rd": 140.0},
        {"country": "UK", "total_ssrs": 6000, "count": 30, "mean_ra": 10.5, "mean_rd": 130.0},
        {"country": "India", "total_ssrs": 5000, "count": 25, "mean_ra": 10.0, "mean_rd": 125.0},
    ]


def _country_metrics():
    np.random.seed(42)
    return {
        "USA": list(np.random.normal(100, 20, 50).astype(float)),
        "China": list(np.random.normal(90, 15, 40).astype(float)),
        "UK": list(np.random.normal(85, 18, 30).astype(float)),
    }


def _heatmap_data():
    return {
        "USA": {"mono": 200, "di": 500, "tri": 800, "tetra": 100, "penta": 50, "hexa": 30},
        "China": {"mono": 180, "di": 450, "tri": 700, "tetra": 90, "penta": 40, "hexa": 25},
        "UK": {"mono": 150, "di": 400, "tri": 600, "tetra": 80, "penta": 35, "hexa": 20},
    }


def _stats_results():
    return [
        {
            "analysis_name": "kruskal_wallis_ssr_count_total",
            "test_name": "Kruskal-Wallis H",
            "statistic": 15.234,
            "p_value": 0.0005,
            "p_value_corrected": 0.003,
            "effect_size": 0.12,
            "n": 120,
        },
        {
            "analysis_name": "correlation_length_ssr",
            "test_name": "Pearson r",
            "statistic": 0.85,
            "p_value": 1.2e-10,
            "p_value_corrected": 8.4e-10,
            "effect_size": 0.85,
            "n": 100,
        },
        {
            "analysis_name": "shannon_entropy_motifs",
            "test_name": "Shannon entropy",
            "statistic": 2.345,
            "p_value": None,
            "p_value_corrected": None,
            "effect_size": 0.78,
            "n": 50,
        },
    ]


# ---------------------------------------------------------------------------
# Tests: Individual figure functions
# ---------------------------------------------------------------------------

class TestMotifSizeDistribution:
    def test_creates_file(self, tmp_path):
        path = plot_motif_size_distribution(_motif_size_data(), tmp_path)
        assert path.exists()
        assert path.suffix == ".png"
        assert path.stat().st_size > 0

    def test_svg_format(self, tmp_path):
        path = plot_motif_size_distribution(_motif_size_data(), tmp_path, fmt="svg")
        assert path.suffix == ".svg"
        assert path.exists()

    def test_partial_sizes(self, tmp_path):
        data = [{"motif_size": "mono", "total_ssrs": 100},
                {"motif_size": "tri", "total_ssrs": 500}]
        path = plot_motif_size_distribution(data, tmp_path)
        assert path.exists()

    def test_creates_output_dir(self, tmp_path):
        nested = tmp_path / "sub" / "dir"
        path = plot_motif_size_distribution(_motif_size_data(), nested)
        assert nested.exists()
        assert path.exists()


class TestTopMotifs:
    def test_creates_file(self, tmp_path):
        path = plot_top_motifs(_motif_freq_data(), tmp_path)
        assert path.exists()
        assert path.stat().st_size > 0

    def test_top_n_limit(self, tmp_path):
        path = plot_top_motifs(_motif_freq_data(), tmp_path, top_n=3)
        assert path.exists()

    def test_svg_format(self, tmp_path):
        path = plot_top_motifs(_motif_freq_data(), tmp_path, fmt="svg")
        assert path.suffix == ".svg"


class TestGeneDistribution:
    def test_creates_file(self, tmp_path):
        path = plot_gene_distribution(_gene_data(), tmp_path)
        assert path.exists()
        assert path.stat().st_size > 0

    def test_single_gene(self, tmp_path):
        data = [{"gene_name": "ORF1ab", "ssr_count": 3000}]
        path = plot_gene_distribution(data, tmp_path)
        assert path.exists()


class TestCountryBoxplot:
    def test_creates_file(self, tmp_path):
        path = plot_country_boxplot(_country_metrics(), tmp_path)
        assert path.exists()
        assert path.stat().st_size > 0

    def test_custom_metric_label(self, tmp_path):
        path = plot_country_boxplot(_country_metrics(), tmp_path,
                                    metric_label="Relative Abundance")
        assert path.exists()

    def test_svg_format(self, tmp_path):
        path = plot_country_boxplot(_country_metrics(), tmp_path, fmt="svg")
        assert path.suffix == ".svg"


class TestCorrelationScatter:
    def test_creates_file(self, tmp_path):
        x = list(range(100))
        y = [v * 2 + np.random.normal(0, 5) for v in x]
        path = plot_correlation_scatter(x, y, "X", "Y", tmp_path)
        assert path.exists()
        assert path.stat().st_size > 0

    def test_custom_filename(self, tmp_path):
        path = plot_correlation_scatter(
            [1, 2, 3], [4, 5, 6], "A", "B", tmp_path,
            filename="custom_scatter",
        )
        assert path.name.startswith("custom_scatter")

    def test_two_points_regression(self, tmp_path):
        path = plot_correlation_scatter([1, 2], [3, 4], "X", "Y", tmp_path)
        assert path.exists()


class TestCountryMotifHeatmap:
    def test_creates_file(self, tmp_path):
        path = plot_country_motif_heatmap(_heatmap_data(), tmp_path)
        assert path.exists()
        assert path.stat().st_size > 0

    def test_svg_format(self, tmp_path):
        path = plot_country_motif_heatmap(_heatmap_data(), tmp_path, fmt="svg")
        assert path.suffix == ".svg"

    def test_single_country(self, tmp_path):
        data = {"USA": {"mono": 200, "tri": 800}}
        path = plot_country_motif_heatmap(data, tmp_path)
        assert path.exists()


class TestChoropleth:
    def test_creates_html(self, tmp_path):
        path = plot_choropleth(_country_summary(), tmp_path)
        assert path.exists()
        assert path.suffix == ".html"
        assert path.stat().st_size > 0

    def test_custom_metric(self, tmp_path):
        path = plot_choropleth(_country_summary(), tmp_path,
                               metric_key="count", title="Accession Counts")
        assert path.exists()


class TestMotifSunburst:
    def test_creates_html(self, tmp_path):
        path = plot_motif_sunburst(_motif_freq_data(), tmp_path)
        assert path.exists()
        assert path.suffix == ".html"
        assert path.stat().st_size > 0


class TestStatsSummaryTable:
    def test_creates_file(self, tmp_path):
        path = plot_stats_summary_table(_stats_results(), tmp_path)
        assert path.exists()
        assert path.stat().st_size > 0

    def test_empty_results(self, tmp_path):
        path = plot_stats_summary_table([], tmp_path)
        # Returns path but file may not exist (no rows to render)
        assert isinstance(path, Path)

    def test_no_p_value_rows(self, tmp_path):
        data = [{"analysis_name": "test", "test_name": "t", "statistic": 1.0,
                 "p_value": None, "p_value_corrected": None, "effect_size": None, "n": 10}]
        path = plot_stats_summary_table(data, tmp_path)
        assert isinstance(path, Path)


# ---------------------------------------------------------------------------
# Tests: Orchestrator
# ---------------------------------------------------------------------------

class TestGenerateAllFigures:
    def test_generates_all_available(self, tmp_path):
        paths = generate_all_figures(
            motif_size_data=_motif_size_data(),
            motif_freq_data=_motif_freq_data(),
            gene_data=_gene_data(),
            country_summary=_country_summary(),
            country_metrics=_country_metrics(),
            correlation_data={"X_vs_Y": ([1, 2, 3, 4, 5], [2, 4, 6, 8, 10])},
            heatmap_data=_heatmap_data(),
            stats_results=_stats_results(),
            output_dir=tmp_path,
        )
        # motif_size + top_motifs + gene_dist + country_boxplot + 1 scatter
        # + heatmap + stats_table + choropleth + sunburst = 9
        assert len(paths) == 9
        for p in paths:
            assert p.exists()

    def test_minimal_data(self, tmp_path):
        paths = generate_all_figures(
            motif_size_data=[],
            motif_freq_data=[],
            gene_data=[],
            country_summary=[],
            country_metrics=None,
            correlation_data=None,
            heatmap_data=None,
            stats_results=None,
            output_dir=tmp_path,
        )
        assert paths == []

    def test_partial_data(self, tmp_path):
        paths = generate_all_figures(
            motif_size_data=_motif_size_data(),
            motif_freq_data=[],
            gene_data=[],
            country_summary=[],
            country_metrics=None,
            correlation_data=None,
            heatmap_data=None,
            stats_results=None,
            output_dir=tmp_path,
        )
        assert len(paths) == 1
        assert "motif_size_distribution" in paths[0].name

    def test_svg_format(self, tmp_path):
        paths = generate_all_figures(
            motif_size_data=_motif_size_data(),
            motif_freq_data=[],
            gene_data=_gene_data(),
            country_summary=[],
            country_metrics=None,
            correlation_data=None,
            heatmap_data=None,
            stats_results=None,
            output_dir=tmp_path,
            fmt="svg",
        )
        for p in paths:
            if p.suffix != ".html":
                assert p.suffix == ".svg"

    def test_multiple_correlations(self, tmp_path):
        paths = generate_all_figures(
            motif_size_data=[],
            motif_freq_data=[],
            gene_data=[],
            country_summary=[],
            country_metrics=None,
            correlation_data={
                "Length_vs_SSR": ([1, 2, 3], [4, 5, 6]),
                "GC_vs_SSR": ([0.3, 0.4, 0.5], [10, 20, 30]),
            },
            heatmap_data=None,
            stats_results=None,
            output_dir=tmp_path,
        )
        assert len(paths) == 2


# ---------------------------------------------------------------------------
# Tests: CLI integration
# ---------------------------------------------------------------------------

class TestVisualizeCLI:
    def test_help_text(self):
        from click.testing import CliRunner
        from gwico_ssr.cli import cli

        runner = CliRunner()
        result = runner.invoke(cli, ["visualize", "--help"])
        assert result.exit_code == 0
        assert "Generate figures" in result.output

    def test_dataset_not_found(self, tmp_path):
        from click.testing import CliRunner
        from gwico_ssr.cli import cli

        db_url = f"sqlite:///{tmp_path / 'test.db'}"
        runner = CliRunner()
        result = runner.invoke(
            cli, ["--config", "/nonexistent.toml", "visualize", "nonexistent"],
            env={"GWICO_SSR_DB_URL": db_url},
        )
        assert "Dataset not found" in result.output

    def test_no_metrics_run(self, tmp_path):
        from click.testing import CliRunner
        from sqlalchemy import create_engine
        from sqlalchemy.orm import sessionmaker

        from gwico_ssr.cli import cli
        from gwico_ssr.db import create_tables
        from gwico_ssr.db.repository import get_or_create_dataset

        db_path = tmp_path / "test.db"
        db_url = f"sqlite:///{db_path}"
        engine = create_engine(db_url)
        create_tables(engine)
        factory = sessionmaker(bind=engine)
        with factory() as session:
            get_or_create_dataset(session, "test_ds", "csv")
            session.commit()

        runner = CliRunner()
        result = runner.invoke(
            cli, ["--config", "/nonexistent.toml", "visualize", "test_ds"],
            env={"GWICO_SSR_DB_URL": db_url},
        )
        assert "No completed metrics run" in result.output

    def test_end_to_end(self, tmp_path):
        """Full end-to-end: seed data and generate all figures."""
        from click.testing import CliRunner
        from sqlalchemy import create_engine
        from sqlalchemy.orm import sessionmaker

        from gwico_ssr.cli import cli
        from gwico_ssr.db import create_tables
        from gwico_ssr.db.repository import (
            create_run,
            get_or_create_dataset,
            insert_ssr_records,
            update_run_status,
            upsert_accession,
            upsert_accession_metrics,
        )

        db_path = tmp_path / "test.db"
        db_url = f"sqlite:///{db_path}"
        engine = create_engine(db_url)
        create_tables(engine)
        factory = sessionmaker(bind=engine)

        with factory() as session:
            ds = get_or_create_dataset(session, "viz_ds", "csv")
            metrics_run = create_run(session, dataset_id=ds.dataset_id, stage="metrics")

            countries = {"USA": 5, "China": 4, "UK": 3}
            acc_idx = 0
            for country, count in countries.items():
                for i in range(count):
                    acc_id = f"ACC_{acc_idx:04d}"
                    acc_idx += 1
                    upsert_accession(
                        session, accession=acc_id, dataset_id=ds.dataset_id,
                        country=country, genome_length=29000 + acc_idx * 100,
                        gc_content=0.35 + acc_idx * 0.001,
                    )
                    upsert_accession_metrics(
                        session, accession=acc_id, run_id=metrics_run.run_id,
                        ssr_count_total=50 + acc_idx * 2, ssr_bp_total=500 + acc_idx * 20,
                        ra=2.0, rd=20.0,
                        mono_count=5, di_count=10, tri_count=15,
                        tetra_count=3, penta_count=2, hexa_count=1,
                        dominant_motif="AAG",
                    )

            # Add SSR records for gene aggregation
            ssrs = []
            for j in range(acc_idx):
                acc_id = f"ACC_{j:04d}"
                for gene_idx, (motif, size) in enumerate([("A", 1), ("AT", 2), ("AAG", 3)]):
                    ssrs.append({
                        "accession": acc_id, "start": 100 + gene_idx * 100,
                        "end": 100 + gene_idx * 100 + size * 3,
                        "motif_raw": motif, "motif_canonical": motif,
                        "motif_size": size, "repeat_units": 3,
                        "repeat_length_bp": size * 3, "strand": "+",
                        "actual_repeat": motif * 3,
                        "detector_version": "gwico-ssr-1.0",
                    })
            insert_ssr_records(session, ssrs)

            update_run_status(session, metrics_run.run_id, "completed")
            session.commit()

        out_dir = tmp_path / "figures"
        runner = CliRunner()
        result = runner.invoke(
            cli,
            ["--config", "/nonexistent.toml", "visualize", "viz_ds", "--output-dir", str(out_dir), "--format", "png"],
            env={"GWICO_SSR_DB_URL": db_url},
        )
        assert result.exit_code == 0, result.output
        assert "Generated" in result.output
        assert out_dir.exists()
        # At least some figure files should exist
        files = list(out_dir.iterdir())
        assert len(files) >= 1
