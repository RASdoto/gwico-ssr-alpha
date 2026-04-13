"""Tests for GWICO-SSR SSR detection engine: motif utils, detector, and CLI."""

from __future__ import annotations

import json
import shutil
from pathlib import Path

import pytest
from sqlalchemy import create_engine, select, func
from sqlalchemy.orm import sessionmaker

from gwico_ssr.db import create_tables
from gwico_ssr.db.repository import (
    get_or_create_dataset,
    insert_ssr_records,
    upsert_accession,
    upsert_sequence_record,
)
from gwico_ssr.models.schema import SSRRecord as SSRRecordModel

FIXTURES_DIR = Path(__file__).parent.parent / "fixtures"


@pytest.fixture
def db_session():
    engine = create_engine("sqlite:///:memory:")
    create_tables(engine)
    factory = sessionmaker(bind=engine)
    session = factory()
    yield session
    session.close()


# ===========================================================================
# MOTIF CANONICALIZATION TESTS
# ===========================================================================

class TestMotifUtils:

    def test_reverse_complement(self):
        from gwico_ssr.ssr.motif import reverse_complement
        assert reverse_complement("ATCG") == "CGAT"
        assert reverse_complement("AAG") == "CTT"
        assert reverse_complement("A") == "T"
        assert reverse_complement("AT") == "AT"  # palindrome
        assert reverse_complement("GC") == "GC"  # palindrome

    def test_all_rotations(self):
        from gwico_ssr.ssr.motif import all_rotations
        rots = all_rotations("AAG")
        assert set(rots) == {"AAG", "AGA", "GAA"}
        rots2 = all_rotations("AT")
        assert set(rots2) == {"AT", "TA"}
        rots1 = all_rotations("A")
        assert rots1 == ["A"]

    def test_canonicalize_motif_rotations(self):
        from gwico_ssr.ssr.motif import canonicalize_motif
        # All rotations of AAG should yield the same canonical
        assert canonicalize_motif("AAG") == "AAG"
        assert canonicalize_motif("AGA") == "AAG"
        assert canonicalize_motif("GAA") == "AAG"

    def test_canonicalize_motif_reverse_complement(self):
        from gwico_ssr.ssr.motif import canonicalize_motif
        # RC of AAG = CTT; rotations CTT, TTC, TCT
        assert canonicalize_motif("CTT") == "AAG"
        assert canonicalize_motif("TTC") == "AAG"
        assert canonicalize_motif("TCT") == "AAG"

    def test_canonicalize_motif_mono(self):
        from gwico_ssr.ssr.motif import canonicalize_motif
        assert canonicalize_motif("A") == "A"
        assert canonicalize_motif("T") == "A"  # RC of A = T, min("A","T") = "A"
        assert canonicalize_motif("C") == "C"
        assert canonicalize_motif("G") == "C"  # RC of C = G, min("C","G") = "C"

    def test_canonicalize_motif_di(self):
        from gwico_ssr.ssr.motif import canonicalize_motif
        assert canonicalize_motif("AT") == "AT"
        assert canonicalize_motif("TA") == "AT"  # rotation of AT
        assert canonicalize_motif("CG") == "CG"
        assert canonicalize_motif("GC") == "CG"  # rotation of CG

    def test_canonicalize_motif_hexa(self):
        from gwico_ssr.ssr.motif import canonicalize_motif
        # Check that a hexamer canonicalizes consistently
        c1 = canonicalize_motif("AAGCTG")
        c2 = canonicalize_motif("AGCTGA")  # rotation
        c3 = canonicalize_motif("CAGCTT")  # RC of AAGCTG
        assert c1 == c2
        assert c1 == c3

    def test_determine_strand_forward(self):
        from gwico_ssr.ssr.motif import determine_strand
        assert determine_strand("AAG", "AAG") == "+"
        assert determine_strand("AGA", "AAG") == "+"  # rotation of canonical

    def test_determine_strand_reverse(self):
        from gwico_ssr.ssr.motif import determine_strand
        assert determine_strand("CTT", "AAG") == "-"  # RC
        assert determine_strand("TTC", "AAG") == "-"  # rotation of RC

    def test_is_sub_repeat(self):
        from gwico_ssr.ssr.motif import is_sub_repeat
        assert is_sub_repeat("AAGAAG") is True   # 2× AAG
        assert is_sub_repeat("AAAA") is True      # 4× A
        assert is_sub_repeat("ATATAT") is True    # 3× AT
        assert is_sub_repeat("AATAAT") is True    # 2× AAT
        assert is_sub_repeat("AAG") is False       # not a sub-repeat
        assert is_sub_repeat("AAGCTG") is False    # not a sub-repeat
        assert is_sub_repeat("A") is False         # single base
        assert is_sub_repeat("AT") is False        # dimer, not a sub-repeat

    def test_canonicalize_case_insensitive(self):
        from gwico_ssr.ssr.motif import canonicalize_motif
        assert canonicalize_motif("aag") == canonicalize_motif("AAG")
        assert canonicalize_motif("Ctt") == canonicalize_motif("CTT")


# ===========================================================================
# SSR DETECTOR TESTS
# ===========================================================================

class TestSSRDetector:

    def test_detect_mononucleotide(self):
        from gwico_ssr.ssr.detector import detect_ssrs, SSRThresholds
        # 12 consecutive A's → mono SSR
        seq = "NNNN" + "A" * 12 + "NNNN"
        result = detect_ssrs(seq, "test", SSRThresholds(mono=10))
        assert result.hit_count >= 1
        mono_hits = [h for h in result.hits if h.motif_size == 1]
        assert len(mono_hits) == 1
        assert mono_hits[0].motif_raw == "A"
        assert mono_hits[0].repeat_units == 12
        assert mono_hits[0].start == 4
        assert mono_hits[0].end == 16

    def test_detect_dinucleotide(self):
        from gwico_ssr.ssr.detector import detect_ssrs, SSRThresholds
        seq = "GCGC" + "AT" * 6 + "GCGC"
        result = detect_ssrs(seq, "test", SSRThresholds(di=5))
        di_hits = [h for h in result.hits if h.motif_size == 2]
        assert len(di_hits) >= 1
        assert di_hits[0].repeat_units == 6
        assert di_hits[0].motif_raw == "AT"
        assert di_hits[0].actual_repeat == "ATATATATATAT"

    def test_detect_trinucleotide(self):
        from gwico_ssr.ssr.detector import detect_ssrs, SSRThresholds
        seq = "TTTT" + "AAG" * 4 + "TTTT"
        result = detect_ssrs(seq, "test", SSRThresholds(tri=3))
        tri_hits = [h for h in result.hits if h.motif_size == 3]
        assert len(tri_hits) >= 1
        h = tri_hits[0]
        assert h.repeat_units == 4
        assert h.motif_raw == "AAG"
        assert h.motif_canonical == "AAG"
        assert h.repeat_length_bp == 12

    def test_detect_tetranucleotide(self):
        from gwico_ssr.ssr.detector import detect_ssrs, SSRThresholds
        seq = "CC" + "AATG" * 4 + "CC"
        result = detect_ssrs(seq, "test", SSRThresholds(tetra=3))
        tetra_hits = [h for h in result.hits if h.motif_size == 4]
        assert len(tetra_hits) >= 1
        assert tetra_hits[0].repeat_units == 4

    def test_detect_pentanucleotide(self):
        from gwico_ssr.ssr.detector import detect_ssrs, SSRThresholds
        seq = "GG" + "AATGC" * 4 + "GG"
        result = detect_ssrs(seq, "test", SSRThresholds(penta=3))
        penta_hits = [h for h in result.hits if h.motif_size == 5]
        assert len(penta_hits) >= 1
        assert penta_hits[0].repeat_units == 4

    def test_detect_hexanucleotide(self):
        from gwico_ssr.ssr.detector import detect_ssrs, SSRThresholds
        seq = "TT" + "AATGCG" * 3 + "TT"
        result = detect_ssrs(seq, "test", SSRThresholds(hexa=2))
        hexa_hits = [h for h in result.hits if h.motif_size == 6]
        assert len(hexa_hits) >= 1
        assert hexa_hits[0].repeat_units == 3

    def test_no_detection_below_threshold(self):
        from gwico_ssr.ssr.detector import detect_ssrs, SSRThresholds
        # Only 9 A's, threshold is 10
        seq = "A" * 9
        result = detect_ssrs(seq, "test", SSRThresholds(mono=10))
        mono_hits = [h for h in result.hits if h.motif_size == 1]
        assert len(mono_hits) == 0

    def test_exact_threshold(self):
        from gwico_ssr.ssr.detector import detect_ssrs, SSRThresholds
        # Exactly 10 A's, threshold is 10
        seq = "A" * 10
        result = detect_ssrs(seq, "test", SSRThresholds(mono=10))
        mono_hits = [h for h in result.hits if h.motif_size == 1]
        assert len(mono_hits) == 1
        assert mono_hits[0].repeat_units == 10

    def test_multiple_ssrs_in_sequence(self):
        from gwico_ssr.ssr.detector import detect_ssrs, SSRThresholds
        # Two separate mono-SSRs
        seq = "A" * 12 + "GCGCGCGC" + "T" * 12
        result = detect_ssrs(seq, "test", SSRThresholds(mono=10))
        mono_hits = [h for h in result.hits if h.motif_size == 1]
        assert len(mono_hits) == 2

    def test_ambiguous_bases_skipped(self):
        from gwico_ssr.ssr.detector import detect_ssrs, SSRThresholds
        # N's in the middle should break the SSR
        seq = "A" * 5 + "N" + "A" * 5
        result = detect_ssrs(seq, "test", SSRThresholds(mono=10))
        mono_hits = [h for h in result.hits if h.motif_size == 1]
        assert len(mono_hits) == 0

    def test_coordinates_zero_based(self):
        from gwico_ssr.ssr.detector import detect_ssrs, SSRThresholds
        seq = "GC" + "AAG" * 4 + "GC"
        result = detect_ssrs(seq, "test", SSRThresholds(tri=3))
        tri_hits = [h for h in result.hits if h.motif_size == 3]
        assert len(tri_hits) == 1
        assert tri_hits[0].start == 2
        assert tri_hits[0].end == 14  # 2 + 3*4

    def test_actual_repeat_content(self):
        from gwico_ssr.ssr.detector import detect_ssrs, SSRThresholds
        seq = "TT" + "AAG" * 5 + "TT"
        result = detect_ssrs(seq, "test", SSRThresholds(tri=3))
        tri_hits = [h for h in result.hits if h.motif_size == 3]
        assert len(tri_hits) == 1
        assert tri_hits[0].actual_repeat == "AAG" * 5

    def test_sub_repeat_filtering(self):
        from gwico_ssr.ssr.detector import detect_ssrs, SSRThresholds
        # ATATAT is 3×AT (dimer) or 1×ATATAT (hexamer)
        # The sub-repeat filter should ensure only the dimer is detected
        seq = "GG" + "AT" * 6 + "GG"
        result = detect_ssrs(seq, "test", SSRThresholds(di=5, hexa=2))
        di_hits = [h for h in result.hits if h.motif_size == 2]
        hexa_hits = [h for h in result.hits if h.motif_size == 6]
        assert len(di_hits) >= 1
        # Hexamer ATATAT is a sub-repeat of AT, so should not appear
        for hh in hexa_hits:
            from gwico_ssr.ssr.motif import is_sub_repeat
            assert not is_sub_repeat(hh.motif_raw)

    def test_canonical_motif_in_results(self):
        from gwico_ssr.ssr.detector import detect_ssrs, SSRThresholds
        # CTT repeats → canonical should be AAG
        seq = "GG" + "CTT" * 4 + "GG"
        result = detect_ssrs(seq, "test", SSRThresholds(tri=3))
        tri_hits = [h for h in result.hits if h.motif_size == 3]
        assert len(tri_hits) == 1
        assert tri_hits[0].motif_raw == "CTT"
        assert tri_hits[0].motif_canonical == "AAG"
        assert tri_hits[0].strand == "-"

    def test_strand_assignment(self):
        from gwico_ssr.ssr.detector import detect_ssrs, SSRThresholds
        # AAG on + strand
        seq1 = "GG" + "AAG" * 4 + "GG"
        r1 = detect_ssrs(seq1, "test", SSRThresholds(tri=3))
        tri1 = [h for h in r1.hits if h.motif_size == 3]
        assert tri1[0].strand == "+"

        # CTT on - strand (RC of AAG)
        seq2 = "GG" + "CTT" * 4 + "GG"
        r2 = detect_ssrs(seq2, "test", SSRThresholds(tri=3))
        tri2 = [h for h in r2.hits if h.motif_size == 3]
        assert tri2[0].strand == "-"

    def test_empty_sequence(self):
        from gwico_ssr.ssr.detector import detect_ssrs
        result = detect_ssrs("", "test")
        assert result.hit_count == 0

    def test_short_sequence(self):
        from gwico_ssr.ssr.detector import detect_ssrs
        result = detect_ssrs("ATCG", "test")
        assert result.hit_count == 0

    def test_thresholds_from_config(self):
        from gwico_ssr.ssr.detector import SSRThresholds

        class FakeSSRSettings:
            min_repeats_mono = 8
            min_repeats_di = 4
            min_repeats_tri = 3
            min_repeats_tetra = 3
            min_repeats_penta = 3
            min_repeats_hexa = 2

        t = SSRThresholds.from_config(FakeSSRSettings())
        assert t.mono == 8
        assert t.di == 4
        assert t.for_size(1) == 8
        assert t.for_size(6) == 2

    def test_count_by_motif_size(self):
        from gwico_ssr.ssr.detector import detect_ssrs, SSRThresholds
        seq = "A" * 12 + "GCGC" + "AAG" * 4
        result = detect_ssrs(seq, "test", SSRThresholds(mono=10, tri=3))
        counts = result.count_by_motif_size()
        assert 1 in counts
        assert 3 in counts

    def test_to_dicts(self):
        from gwico_ssr.ssr.detector import detect_ssrs, SSRThresholds
        seq = "A" * 12
        result = detect_ssrs(seq, "test_acc", SSRThresholds(mono=10))
        dicts = result.to_dicts()
        assert len(dicts) >= 1
        d = dicts[0]
        assert d["accession"] == "test_acc"
        assert "start" in d
        assert "end" in d
        assert "motif_raw" in d
        assert "motif_canonical" in d
        assert "detector_version" in d

    def test_detect_ssrs_default_thresholds(self):
        from gwico_ssr.ssr.detector import detect_ssrs
        # 10 A's should be detected with default thresholds (mono=10)
        seq = "A" * 10
        result = detect_ssrs(seq, "test")
        assert result.hit_count >= 1

    def test_detect_ssrs_from_fasta(self):
        from gwico_ssr.ssr.detector import detect_ssrs_from_fasta, SSRThresholds

        results = detect_ssrs_from_fasta(
            str(FIXTURES_DIR / "sample.fasta"),
            SSRThresholds(mono=5, di=3, tri=3, tetra=2, penta=2, hexa=2),
        )
        assert len(results) == 1
        # The sample FASTA should have some SSRs at low thresholds
        assert results[0].sequence_length > 0

    def test_perf_style_validation(self):
        """Validate against PERF-style expected output for NC_045512.2 snippet.

        PERF output (from context.md): NC_045512.2  626  635  AAG  9  -  3  TTC
        This means an AAG motif at positions 626-635 (1-based), 9 bp, strand -,
        motif size 3, actual repeat TTC.

        We test with a synthetic sequence that should produce a known SSR.
        """
        from gwico_ssr.ssr.detector import detect_ssrs, SSRThresholds

        # Create a sequence with TTC repeated 3 times at a known position
        # TTC×3 = TTCTTCTTC (9 bp). canonical(TTC) = AAG, strand = -
        prefix = "A" * 10
        repeat = "TTC" * 3
        suffix = "G" * 10
        seq = prefix + repeat + suffix
        result = detect_ssrs(seq, "NC_045512.2", SSRThresholds(tri=3))
        tri_hits = [h for h in result.hits if h.motif_size == 3]
        assert len(tri_hits) == 1
        h = tri_hits[0]
        assert h.start == 10
        assert h.end == 19
        assert h.motif_raw == "TTC"
        assert h.motif_canonical == "AAG"
        assert h.strand == "-"
        assert h.repeat_length_bp == 9
        assert h.actual_repeat == "TTCTTCTTC"


# ===========================================================================
# DB PERSISTENCE TESTS
# ===========================================================================

class TestSSRPersistence:

    def test_insert_ssr_records(self, db_session):
        ds = get_or_create_dataset(db_session, "test_ds", "csv")
        upsert_accession(db_session, accession="ACC1", dataset_id=ds.dataset_id)

        from gwico_ssr.ssr.detector import detect_ssrs, SSRThresholds
        seq = "GG" + "AAG" * 5 + "GG"
        result = detect_ssrs(seq, "ACC1", SSRThresholds(tri=3))
        dicts = result.to_dicts()
        count = insert_ssr_records(db_session, dicts)
        assert count > 0

        ssr_count = db_session.execute(
            select(func.count()).select_from(SSRRecordModel).where(SSRRecordModel.accession == "ACC1")
        ).scalar()
        assert ssr_count > 0

    def test_detector_version_in_records(self, db_session):
        ds = get_or_create_dataset(db_session, "test_ds", "csv")
        upsert_accession(db_session, accession="ACC1", dataset_id=ds.dataset_id)

        from gwico_ssr.ssr.detector import detect_ssrs, DETECTOR_VERSION, SSRThresholds
        seq = "A" * 12
        result = detect_ssrs(seq, "ACC1", SSRThresholds(mono=10))
        dicts = result.to_dicts()
        insert_ssr_records(db_session, dicts)

        records = db_session.execute(
            select(SSRRecordModel).where(SSRRecordModel.accession == "ACC1")
        ).scalars().all()
        for r in records:
            assert r.detector_version == DETECTOR_VERSION


# ===========================================================================
# CLI DETECT COMMAND TESTS
# ===========================================================================

class TestDetectCLI:

    def test_detect_help(self, tmp_path, sample_toml):
        from click.testing import CliRunner
        from gwico_ssr.cli import cli

        runner = CliRunner()
        result = runner.invoke(cli, ["--config", str(sample_toml), "detect", "--help"])
        assert result.exit_code == 0
        assert "Detect perfect SSRs" in result.output

    def test_detect_no_dataset(self, tmp_path):
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
        result = runner.invoke(cli, ["--config", str(config), "detect", "nonexistent"])
        assert "Dataset not found" in result.output or "No accessions" in result.output

    def test_detect_with_fasta(self, tmp_path):
        """End-to-end: create accession, place FASTA, detect SSRs."""
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
min_repeats_mono = 5
min_repeats_di = 3
min_repeats_tri = 3
min_repeats_tetra = 2
min_repeats_penta = 2
min_repeats_hexa = 2
[logging]
level = "WARNING"
format = "text"
file = ""
[output]
dir = "{out_dir}"
""", encoding="utf-8")

        # Set up DB
        engine = create_engine(f"sqlite:///{db_path}")
        create_tables(engine)
        factory = sessionmaker(bind=engine)
        sess = factory()
        ds = get_or_create_dataset(sess, "myds", "csv")
        upsert_accession(sess, accession="TEST1", dataset_id=ds.dataset_id)
        sess.commit()
        sess.close()
        engine.dispose()

        # Write a FASTA with known SSRs
        fasta_dir = tmp_path / "sequences" / "fasta"
        fasta_dir.mkdir(parents=True)
        fasta_content = ">TEST1 Test sequence\n" + "AAAAAAAAAA" + "GC" * 10 + "AAGAAGAAGAAG" + "\n"
        (fasta_dir / "TEST1.fasta").write_text(fasta_content)

        runner = CliRunner()
        result = runner.invoke(cli, [
            "--config", str(config),
            "detect", "myds",
            "--accessions", "TEST1",
        ])
        assert result.exit_code == 0
        assert "Processed: 1" in result.output
        assert "Total SSRs found:" in result.output

    def test_detect_json_summary(self, tmp_path):
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
min_repeats_mono = 5
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
        sess.commit()
        sess.close()
        engine.dispose()

        fasta_dir = tmp_path / "sequences" / "fasta"
        fasta_dir.mkdir(parents=True)
        (fasta_dir / "TEST1.fasta").write_text(">TEST1\n" + "A" * 20 + "\n")

        runner = CliRunner()
        result = runner.invoke(cli, [
            "--config", str(config),
            "detect", "myds",
            "--accessions", "TEST1",
            "--json-summary",
        ])
        assert result.exit_code == 0
        data = json.loads(result.output)
        assert data["processed"] == 1
        assert data["total_ssrs_found"] > 0


# ===========================================================================
# SSR PACKAGE IMPORT TESTS
# ===========================================================================

class TestSSRImports:

    def test_import_ssr_package(self):
        from gwico_ssr.ssr import (
            DETECTOR_VERSION,
            DetectionResult,
            SSRHit,
            SSRThresholds,
            all_rotations,
            canonicalize_motif,
            detect_ssrs,
            detect_ssrs_from_fasta,
            determine_strand,
            is_sub_repeat,
            reverse_complement,
        )
        assert detect_ssrs is not None
        assert canonicalize_motif is not None
        assert DETECTOR_VERSION is not None
