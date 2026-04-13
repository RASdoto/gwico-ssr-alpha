"""Tests for GWICO-SSR export and provenance layer."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import pytest
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

from gwico_ssr.db import create_tables
from gwico_ssr.db.repository import (
    create_run,
    get_or_create_dataset,
    insert_ssr_records,
    insert_statistical_result,
    update_run_status,
    upsert_accession,
    upsert_accession_metrics,
)
from gwico_ssr.export.exporters import (
    export_all,
    export_metrics_csv,
    export_ssrs_bed,
    export_ssrs_csv,
    export_ssrs_gff3,
    export_ssrs_json,
    export_stats_csv,
    export_stats_json,
    export_publication_tables,
    generate_run_manifest,
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


@pytest.fixture
def seeded_session(db_session):
    """Session with dataset, accessions, SSRs, metrics, and stats."""
    session = db_session
    ds = get_or_create_dataset(session, "export_ds", "csv")
    run = create_run(session, dataset_id=ds.dataset_id, stage="metrics")

    countries = {"USA": 3, "China": 2}
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
                session, accession=acc_id, run_id=run.run_id,
                ssr_count_total=50 + acc_idx * 2, ssr_bp_total=500 + acc_idx * 20,
                ra=2.0, rd=20.0,
                mono_count=5, di_count=10, tri_count=15,
                tetra_count=3, penta_count=2, hexa_count=1,
                dominant_motif="AAG",
            )

    # SSR records
    ssrs = []
    for j in range(acc_idx):
        acc_id = f"ACC_{j:04d}"
        for idx, (motif, size, bp) in enumerate([("A", 1, 10), ("AT", 2, 12), ("AAG", 3, 9)]):
            ssrs.append({
                "accession": acc_id, "start": 100 + idx * 200,
                "end": 100 + idx * 200 + bp,
                "motif_raw": motif, "motif_canonical": motif,
                "motif_size": size, "repeat_units": bp // size,
                "repeat_length_bp": bp, "strand": "+",
                "actual_repeat": motif * (bp // size),
                "detector_version": "gwico-ssr-1.0",
                "run_id": run.run_id,
            })
    insert_ssr_records(session, ssrs)

    # Statistical results
    insert_statistical_result(
        session, run_id=run.run_id,
        analysis_name="test_analysis", grouping="country",
        metric="ssr_count_total", test_name="Kruskal-Wallis H",
        statistic=15.0, p_value=0.001, p_value_corrected=0.005,
        effect_size=0.12, ci_lower=None, ci_upper=None, n=5,
        metadata_json=json.dumps({"note": "test"}),
    )

    update_run_status(session, run.run_id, "completed")
    session.flush()

    return session, ds, run


# ---------------------------------------------------------------------------
# CSV Export Tests
# ---------------------------------------------------------------------------

class TestCSVExport:
    def test_ssr_csv(self, seeded_session, tmp_path):
        session, ds, run = seeded_session
        path = export_ssrs_csv(session, tmp_path, run_id=run.run_id)
        assert path.exists()
        assert path.name == "ssr_records.csv"
        with open(path) as f:
            reader = csv.DictReader(f)
            rows = list(reader)
        assert len(rows) == 15  # 5 accessions * 3 SSRs
        assert "motif_canonical" in rows[0]
        assert "accession" in rows[0]

    def test_ssr_csv_filter_motif(self, seeded_session, tmp_path):
        session, ds, run = seeded_session
        path = export_ssrs_csv(session, tmp_path, run_id=run.run_id, motif="A")
        with open(path) as f:
            rows = list(csv.DictReader(f))
        assert len(rows) == 5  # Only mono "A"
        assert all(r["motif_canonical"] == "A" for r in rows)

    def test_metrics_csv(self, seeded_session, tmp_path):
        session, ds, run = seeded_session
        path = export_metrics_csv(session, tmp_path, run_id=run.run_id)
        assert path.exists()
        assert path.name == "accession_metrics.csv"
        with open(path) as f:
            rows = list(csv.DictReader(f))
        assert len(rows) == 5

    def test_stats_csv(self, seeded_session, tmp_path):
        session, ds, run = seeded_session
        path = export_stats_csv(session, tmp_path, run_id=run.run_id)
        assert path.exists()
        with open(path) as f:
            rows = list(csv.DictReader(f))
        assert len(rows) == 1
        assert rows[0]["analysis_name"] == "test_analysis"


# ---------------------------------------------------------------------------
# JSON Export Tests
# ---------------------------------------------------------------------------

class TestJSONExport:
    def test_ssr_json(self, seeded_session, tmp_path):
        session, ds, run = seeded_session
        path = export_ssrs_json(session, tmp_path, run_id=run.run_id)
        assert path.exists()
        data = json.loads(path.read_text())
        assert len(data) == 15
        assert "motif_canonical" in data[0]

    def test_stats_json(self, seeded_session, tmp_path):
        session, ds, run = seeded_session
        path = export_stats_json(session, tmp_path, run_id=run.run_id)
        assert path.exists()
        data = json.loads(path.read_text())
        assert len(data) == 1
        assert data[0]["analysis_name"] == "test_analysis"
        assert data[0]["metadata"] == {"note": "test"}


# ---------------------------------------------------------------------------
# BED Export Tests
# ---------------------------------------------------------------------------

class TestBEDExport:
    def test_bed_format(self, seeded_session, tmp_path):
        session, ds, run = seeded_session
        path = export_ssrs_bed(session, tmp_path, run_id=run.run_id)
        assert path.exists()
        assert path.suffix == ".bed"
        lines = path.read_text().strip().split("\n")
        assert len(lines) == 15
        # Verify BED6 columns
        fields = lines[0].split("\t")
        assert len(fields) == 6
        # chrom, start, end, name, score, strand
        assert fields[0].startswith("ACC_")
        assert int(fields[1]) >= 0  # 0-based start
        assert int(fields[2]) > int(fields[1])  # end > start
        assert "x" in fields[3]  # name like "Ax10"
        assert fields[5] in ("+", "-")

    def test_bed_sorted(self, seeded_session, tmp_path):
        session, ds, run = seeded_session
        path = export_ssrs_bed(session, tmp_path, run_id=run.run_id)
        lines = path.read_text().strip().split("\n")
        prev_chrom, prev_start = "", -1
        for line in lines:
            fields = line.split("\t")
            chrom = fields[0]
            start = int(fields[1])
            if chrom == prev_chrom:
                assert start >= prev_start
            prev_chrom, prev_start = chrom, start


# ---------------------------------------------------------------------------
# GFF3 Export Tests
# ---------------------------------------------------------------------------

class TestGFF3Export:
    def test_gff3_format(self, seeded_session, tmp_path):
        session, ds, run = seeded_session
        path = export_ssrs_gff3(session, tmp_path, run_id=run.run_id)
        assert path.exists()
        assert path.suffix == ".gff3"
        lines = path.read_text().strip().split("\n")
        assert lines[0] == "##gff-version 3"
        data_lines = [l for l in lines if not l.startswith("#")]
        assert len(data_lines) == 15

    def test_gff3_coordinates(self, seeded_session, tmp_path):
        session, ds, run = seeded_session
        path = export_ssrs_gff3(session, tmp_path, run_id=run.run_id)
        lines = path.read_text().strip().split("\n")
        data_lines = [l for l in lines if not l.startswith("#")]
        # GFF3 is 1-based closed
        fields = data_lines[0].split("\t")
        assert len(fields) == 9
        gff_start = int(fields[3])
        assert gff_start >= 1  # 1-based

    def test_gff3_attributes(self, seeded_session, tmp_path):
        session, ds, run = seeded_session
        path = export_ssrs_gff3(session, tmp_path, run_id=run.run_id)
        lines = path.read_text().strip().split("\n")
        data_lines = [l for l in lines if not l.startswith("#")]
        attrs = data_lines[0].split("\t")[8]
        assert "ID=ssr_" in attrs
        assert "motif_canonical=" in attrs
        assert "repeat_units=" in attrs

    def test_gff3_source_column(self, seeded_session, tmp_path):
        session, ds, run = seeded_session
        path = export_ssrs_gff3(session, tmp_path, run_id=run.run_id)
        lines = path.read_text().strip().split("\n")
        data_lines = [l for l in lines if not l.startswith("#")]
        fields = data_lines[0].split("\t")
        assert fields[1] == "gwico-ssr"
        assert fields[2] == "microsatellite"


# ---------------------------------------------------------------------------
# Publication Tables Tests
# ---------------------------------------------------------------------------

class TestPublicationTables:
    def test_generates_tables(self, seeded_session, tmp_path):
        session, ds, run = seeded_session
        paths = export_publication_tables(
            session, tmp_path, run_id=run.run_id, dataset_id=ds.dataset_id,
        )
        assert len(paths) == 4
        names = {p.name for p in paths}
        assert "dataset_summary.json" in names
        assert "country_summary.csv" in names
        assert "motif_frequencies.csv" in names
        assert "gene_distribution.csv" in names

    def test_country_summary_content(self, seeded_session, tmp_path):
        session, ds, run = seeded_session
        paths = export_publication_tables(
            session, tmp_path, run_id=run.run_id, dataset_id=ds.dataset_id,
        )
        country_path = [p for p in paths if p.name == "country_summary.csv"][0]
        with open(country_path) as f:
            rows = list(csv.DictReader(f))
        countries = {r["country"] for r in rows}
        assert "USA" in countries
        assert "China" in countries


# ---------------------------------------------------------------------------
# Run Manifest Tests
# ---------------------------------------------------------------------------

class TestRunManifest:
    def test_manifest_structure(self, seeded_session, tmp_path):
        session, ds, run = seeded_session
        path = generate_run_manifest(session, tmp_path, run_id=run.run_id)
        assert path.exists()
        assert path.name == "run_manifest.json"
        data = json.loads(path.read_text())
        assert "gwico_ssr_version" in data
        assert "generated_at" in data
        assert "run" in data
        assert data["run"]["run_id"] == run.run_id
        assert data["run"]["dataset_name"] == "export_ds"
        assert "counts" in data

    def test_manifest_with_files(self, seeded_session, tmp_path):
        session, ds, run = seeded_session
        # Create a dummy file
        dummy = tmp_path / "dummy.csv"
        dummy.write_text("a,b\n1,2\n")
        path = generate_run_manifest(
            session, tmp_path, run_id=run.run_id,
            output_files=[dummy],
        )
        data = json.loads(path.read_text())
        assert len(data["output_files"]) == 1
        assert data["output_files"][0]["sha256"]
        assert data["output_files"][0]["size_bytes"] > 0

    def test_manifest_with_config(self, seeded_session, tmp_path):
        session, ds, run = seeded_session
        path = generate_run_manifest(
            session, tmp_path, run_id=run.run_id,
            config_snapshot={"db": "sqlite:///test.db"},
        )
        data = json.loads(path.read_text())
        assert data["config_snapshot"]["db"] == "sqlite:///test.db"

    def test_manifest_invalid_run(self, db_session, tmp_path):
        with pytest.raises(ValueError, match="Run 9999 not found"):
            generate_run_manifest(db_session, tmp_path, run_id=9999)


# ---------------------------------------------------------------------------
# Export All Orchestrator Tests
# ---------------------------------------------------------------------------

class TestExportAll:
    def test_all_formats(self, seeded_session, tmp_path):
        session, ds, run = seeded_session
        paths = export_all(
            session, tmp_path,
            run_id=run.run_id, dataset_id=ds.dataset_id,
        )
        names = {p.name for p in paths}
        assert "ssr_records.csv" in names
        assert "accession_metrics.csv" in names
        assert "statistical_results.csv" in names
        assert "ssr_records.json" in names
        assert "statistical_results.json" in names
        assert "ssr_records.bed" in names
        assert "ssr_records.gff3" in names
        assert "dataset_summary.json" in names
        assert "run_manifest.json" in names

    def test_selective_formats(self, seeded_session, tmp_path):
        session, ds, run = seeded_session
        paths = export_all(
            session, tmp_path,
            run_id=run.run_id, dataset_id=ds.dataset_id,
            formats={"bed", "gff3"},
        )
        names = {p.name for p in paths}
        assert "ssr_records.bed" in names
        assert "ssr_records.gff3" in names
        assert "run_manifest.json" in names  # Always included
        assert "ssr_records.csv" not in names

    def test_filter_by_country(self, seeded_session, tmp_path):
        session, ds, run = seeded_session
        paths = export_all(
            session, tmp_path,
            run_id=run.run_id, dataset_id=ds.dataset_id,
            formats={"csv"},
            country="USA",
        )
        csv_path = [p for p in paths if p.name == "ssr_records.csv"][0]
        with open(csv_path) as f:
            rows = list(csv.DictReader(f))
        assert len(rows) == 9  # 3 USA accessions * 3 SSRs


# ---------------------------------------------------------------------------
# CLI Tests
# ---------------------------------------------------------------------------

class TestExportCLI:
    def test_help_text(self):
        from click.testing import CliRunner
        from gwico_ssr.cli import cli

        runner = CliRunner()
        result = runner.invoke(cli, ["export", "--help"])
        assert result.exit_code == 0
        assert "Export analysis outputs" in result.output

    def test_dataset_not_found(self, tmp_path):
        from click.testing import CliRunner
        from gwico_ssr.cli import cli

        db_url = f"sqlite:///{tmp_path / 'test.db'}"
        runner = CliRunner()
        result = runner.invoke(
            cli, ["--config", "/nonexistent.toml", "export", "nonexistent"],
            env={"GWICO_SSR_DB_URL": db_url},
        )
        assert "Dataset not found" in result.output

    def test_end_to_end(self, tmp_path):
        from click.testing import CliRunner
        from sqlalchemy import create_engine
        from sqlalchemy.orm import sessionmaker

        from gwico_ssr.cli import cli
        from gwico_ssr.db import create_tables as ct
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
        ct(engine)
        factory = sessionmaker(bind=engine)

        with factory() as session:
            ds = get_or_create_dataset(session, "exp_ds", "csv")
            run = create_run(session, dataset_id=ds.dataset_id, stage="metrics")
            upsert_accession(
                session, accession="ACC_0001", dataset_id=ds.dataset_id,
                country="USA", genome_length=29000, gc_content=0.38,
            )
            upsert_accession_metrics(
                session, accession="ACC_0001", run_id=run.run_id,
                ssr_count_total=50, ssr_bp_total=500, ra=2.0, rd=20.0,
                mono_count=5, di_count=10, tri_count=15,
                tetra_count=3, penta_count=2, hexa_count=1,
                dominant_motif="AAG",
            )
            insert_ssr_records(session, [{
                "accession": "ACC_0001", "start": 100, "end": 109,
                "motif_raw": "AAG", "motif_canonical": "AAG",
                "motif_size": 3, "repeat_units": 3,
                "repeat_length_bp": 9, "strand": "+",
                "actual_repeat": "AAGAAGAAG",
                "detector_version": "gwico-ssr-1.0",
                "run_id": run.run_id,
            }])
            update_run_status(session, run.run_id, "completed")
            session.commit()

        out_dir = tmp_path / "exports"
        runner = CliRunner()
        result = runner.invoke(
            cli,
            ["--config", "/nonexistent.toml", "export", "exp_ds",
             "--output-dir", str(out_dir), "--formats", "csv,bed,gff3"],
            env={"GWICO_SSR_DB_URL": db_url},
        )
        assert result.exit_code == 0, result.output
        assert "Exported" in result.output
        assert out_dir.exists()
        files = {f.name for f in out_dir.iterdir()}
        assert "ssr_records.csv" in files
        assert "ssr_records.bed" in files
        assert "ssr_records.gff3" in files
        assert "run_manifest.json" in files
