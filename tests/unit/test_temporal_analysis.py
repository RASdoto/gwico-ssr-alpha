"""Tests for temporal analysis module."""

import pytest
from datetime import datetime, timedelta
from sqlalchemy import create_engine
from sqlalchemy.orm import Session

from gwico_ssr.analysis.temporal import (
    DatePrecision,
    TimeBin,
    TemporalSummary,
    parse_date_string,
    normalize_date_to_precision,
    bin_id_from_date,
    create_time_bins,
    get_date_range_for_field,
    query_accessions_in_time_bin,
    query_temporal_ssr_metrics,
    compute_temporal_summary,
)
from gwico_ssr.models.schema import Base, Accession, SSRRecord, Run, AccessionMetrics
from gwico_ssr.db.repository import create_run, get_or_create_dataset, insert_ssr_records, upsert_accession, upsert_accession_metrics


@pytest.fixture
def in_memory_db():
    """Create an in-memory SQLite database for testing."""
    engine = create_engine("sqlite:///:memory:")
    Base.metadata.create_all(engine)
    return engine


@pytest.fixture
def session(in_memory_db):
    """Create a test session."""
    with Session(in_memory_db) as sess:
        yield sess


# ---------------------------------------------------------------------------
# Date parsing tests
# ---------------------------------------------------------------------------

class TestParseDateString:
    """Test parse_date_string function."""

    def test_parse_full_iso_date(self):
        """Parse full ISO 8601 date."""
        dt = parse_date_string("2024-04-25")
        assert dt is not None
        assert dt.year == 2024
        assert dt.month == 4
        assert dt.day == 25

    def test_parse_iso_datetime(self):
        """Parse ISO 8601 datetime."""
        dt = parse_date_string("2024-04-25T10:30:00")
        assert dt is not None
        assert dt.hour == 10
        assert dt.minute == 30

    def test_parse_partial_month(self):
        """Parse partial month date."""
        dt = parse_date_string("2024-04")
        assert dt is not None
        assert dt.year == 2024
        assert dt.month == 4
        assert dt.day == 1

    def test_parse_quarter(self):
        """Parse quarter notation."""
        dt = parse_date_string("2024-Q1")
        assert dt is not None
        assert dt.year == 2024
        assert dt.month == 1

        dt = parse_date_string("2024-Q3")
        assert dt is not None
        assert dt.month == 7

    def test_parse_year_only(self):
        """Parse year-only date."""
        dt = parse_date_string("2024")
        assert dt is not None
        assert dt.year == 2024
        assert dt.month == 1
        assert dt.day == 1

    def test_parse_invalid_date(self):
        """Invalid date strings return None."""
        assert parse_date_string("not-a-date") is None
        assert parse_date_string("2024-13") is None
        assert parse_date_string("") is None
        assert parse_date_string(None) is None

    def test_parse_invalid_quarter(self):
        """Invalid quarter returns None."""
        assert parse_date_string("2024-Q5") is None
        assert parse_date_string("2024-Q0") is None


# ---------------------------------------------------------------------------
# Date precision tests
# ---------------------------------------------------------------------------

class TestNormalizeDateToPrecision:
    """Test normalize_date_to_precision function."""

    def test_normalize_to_year(self):
        """Normalize to YEAR precision."""
        dt = datetime(2024, 4, 25, 10, 30, 45)
        normalized = normalize_date_to_precision(dt, DatePrecision.YEAR)
        assert normalized.year == 2024
        assert normalized.month == 1
        assert normalized.day == 1
        assert normalized.hour == 0

    def test_normalize_to_quarter(self):
        """Normalize to QUARTER precision."""
        # Q1
        dt = datetime(2024, 2, 15)
        normalized = normalize_date_to_precision(dt, DatePrecision.QUARTER)
        assert normalized.month == 1

        # Q3
        dt = datetime(2024, 8, 15)
        normalized = normalize_date_to_precision(dt, DatePrecision.QUARTER)
        assert normalized.month == 7

    def test_normalize_to_month(self):
        """Normalize to MONTH precision."""
        dt = datetime(2024, 4, 25, 10, 30)
        normalized = normalize_date_to_precision(dt, DatePrecision.MONTH)
        assert normalized.year == 2024
        assert normalized.month == 4
        assert normalized.day == 1
        assert normalized.hour == 0

    def test_normalize_to_week(self):
        """Normalize to WEEK precision (Monday)."""
        # 2024-04-25 is a Thursday
        dt = datetime(2024, 4, 25)
        normalized = normalize_date_to_precision(dt, DatePrecision.WEEK)
        # Monday of that week should be 2024-04-22
        assert normalized.weekday() == 0  # Monday
        assert normalized.day == 22

    def test_normalize_to_day(self):
        """Normalize to DAY precision."""
        dt = datetime(2024, 4, 25, 10, 30, 45)
        normalized = normalize_date_to_precision(dt, DatePrecision.DAY)
        assert normalized.year == 2024
        assert normalized.month == 4
        assert normalized.day == 25
        assert normalized.hour == 0
        assert normalized.minute == 0


# ---------------------------------------------------------------------------
# Bin ID tests
# ---------------------------------------------------------------------------

class TestBinIdFromDate:
    """Test bin_id_from_date function."""

    def test_bin_id_year(self):
        """Generate YEAR bin ID."""
        dt = datetime(2024, 4, 25)
        bin_id = bin_id_from_date(dt, DatePrecision.YEAR)
        assert bin_id == "2024"

    def test_bin_id_quarter(self):
        """Generate QUARTER bin ID."""
        dt = datetime(2024, 4, 25)
        bin_id = bin_id_from_date(dt, DatePrecision.QUARTER)
        assert bin_id == "2024-Q2"

        dt = datetime(2024, 1, 1)
        bin_id = bin_id_from_date(dt, DatePrecision.QUARTER)
        assert bin_id == "2024-Q1"

    def test_bin_id_month(self):
        """Generate MONTH bin ID."""
        dt = datetime(2024, 4, 25)
        bin_id = bin_id_from_date(dt, DatePrecision.MONTH)
        assert bin_id == "2024-04"

    def test_bin_id_week(self):
        """Generate WEEK bin ID."""
        dt = datetime(2024, 4, 25)
        bin_id = bin_id_from_date(dt, DatePrecision.WEEK)
        # Should have format YYYY-Www
        assert "2024" in bin_id
        assert "W" in bin_id

    def test_bin_id_day(self):
        """Generate DAY bin ID."""
        dt = datetime(2024, 4, 25)
        bin_id = bin_id_from_date(dt, DatePrecision.DAY)
        assert bin_id == "2024-04-25"


# ---------------------------------------------------------------------------
# Time bin generation tests
# ---------------------------------------------------------------------------

class TestCreateTimeBins:
    """Test create_time_bins function."""

    def test_create_yearly_bins(self):
        """Generate yearly bins."""
        start = datetime(2022, 1, 1)
        end = datetime(2024, 12, 31)
        bins = create_time_bins(start, end, DatePrecision.YEAR)

        assert len(bins) == 3  # 2022, 2023, 2024
        assert bins[0][0].year == 2022
        assert bins[1][0].year == 2023
        assert bins[2][0].year == 2024

    def test_create_monthly_bins(self):
        """Generate monthly bins."""
        start = datetime(2024, 1, 1)
        end = datetime(2024, 3, 31)
        bins = create_time_bins(start, end, DatePrecision.MONTH)

        assert len(bins) == 3  # Jan, Feb, Mar
        assert bins[0][0].month == 1
        assert bins[1][0].month == 2
        assert bins[2][0].month == 3

    def test_create_quarterly_bins(self):
        """Generate quarterly bins."""
        start = datetime(2024, 1, 1)
        end = datetime(2024, 12, 31)
        bins = create_time_bins(start, end, DatePrecision.QUARTER)

        assert len(bins) == 4  # Q1, Q2, Q3, Q4

    def test_create_daily_bins(self):
        """Generate daily bins."""
        start = datetime(2024, 4, 1)
        end = datetime(2024, 4, 5)
        bins = create_time_bins(start, end, DatePrecision.DAY)

        assert len(bins) == 5  # 4 days inclusive
        assert bins[0][0].day == 1
        assert bins[1][0].day == 2

    def test_bin_boundaries_are_half_open(self):
        """Verify bins use half-open intervals [start, end)."""
        start = datetime(2024, 1, 1)
        end = datetime(2024, 1, 5)
        bins = create_time_bins(start, end, DatePrecision.DAY)

        # Each bin's end should be the next bin's start
        for i in range(len(bins) - 1):
            assert bins[i][1] == bins[i + 1][0]


# ---------------------------------------------------------------------------
# TimeBin and TemporalSummary tests
# ---------------------------------------------------------------------------

class TestTimeBinDataclass:
    """Test TimeBin dataclass."""

    def test_timebin_creation(self):
        """Create TimeBin with data."""
        bin_obj = TimeBin(
            bin_id="2024-01",
            start_date=datetime(2024, 1, 1),
            end_date=datetime(2024, 2, 1),
            accession_count=100,
            ssr_count_total=500,
        )
        assert bin_obj.bin_id == "2024-01"
        assert bin_obj.accession_count == 100
        assert bin_obj.ssr_count_total == 500

    def test_timebin_default_values(self):
        """TimeBin fields default to 0."""
        bin_obj = TimeBin(
            bin_id="2024-01",
            start_date=datetime(2024, 1, 1),
            end_date=datetime(2024, 2, 1),
        )
        assert bin_obj.accession_count == 0
        assert bin_obj.ssr_count_total == 0
        assert bin_obj.perfect_count == 0


class TestTemporalSummary:
    """Test TemporalSummary dataclass."""

    def test_temporal_summary_creation(self):
        """Create TemporalSummary."""
        summary = TemporalSummary(
            date_field="release_date",
            precision=DatePrecision.MONTH,
        )
        assert summary.date_field == "release_date"
        assert summary.precision == DatePrecision.MONTH
        assert len(summary.bins) == 0

    def test_total_ssrs(self):
        """Sum SSRs across bins."""
        summary = TemporalSummary(date_field="release_date", precision=DatePrecision.MONTH)
        summary.bins = [
            TimeBin(
                bin_id="2024-01",
                start_date=datetime(2024, 1, 1),
                end_date=datetime(2024, 2, 1),
                ssr_count_total=100,
            ),
            TimeBin(
                bin_id="2024-02",
                start_date=datetime(2024, 2, 1),
                end_date=datetime(2024, 3, 1),
                ssr_count_total=150,
            ),
        ]
        assert summary.total_ssrs() == 250

    def test_total_ssrs_bp(self):
        """Sum SSR bp across bins."""
        summary = TemporalSummary(date_field="release_date", precision=DatePrecision.MONTH)
        summary.bins = [
            TimeBin(
                bin_id="2024-01",
                start_date=datetime(2024, 1, 1),
                end_date=datetime(2024, 2, 1),
                ssr_bp_total=1000,
            ),
            TimeBin(
                bin_id="2024-02",
                start_date=datetime(2024, 2, 1),
                end_date=datetime(2024, 3, 1),
                ssr_bp_total=1500,
            ),
        ]
        assert summary.total_ssrs_bp() == 2500


# ---------------------------------------------------------------------------
# Database query tests
# ---------------------------------------------------------------------------

class TestDatabaseQueries:
    """Test database query functions."""

    @pytest.fixture
    def populated_db(self, session):
        """Populate database with test data."""
        # Create dataset
        dataset = get_or_create_dataset(session, "test_dataset", "test_description")
        session.flush()

        # Create run
        run = create_run(session, dataset_id=dataset.dataset_id, stage="detect")
        session.flush()

        # Create accessions with dates
        accessions_data = [
            ("ACC001", "2024-01-15", "2024-01-10"),
            ("ACC002", "2024-01-20", "2024-01-05"),
            ("ACC003", "2024-02-10", "2024-02-01"),
            ("ACC004", "2024-02-15", "2024-02-05"),
            ("ACC005", "2024-03-01", "2024-03-01"),
        ]

        for acc, release, collection in accessions_data:
            upsert_accession(
                session,
                accession=acc,
                dataset_id=dataset.dataset_id,
                release_date=release,
                collection_date=collection,
            )
        session.flush()

        # Create SSRs as dictionaries
        ssr_data = [
            {"accession": "ACC001", "repeat_class": "perfect", "repeat_units": 20, "repeat_length_bp": 100},
            {"accession": "ACC001", "repeat_class": "perfect", "repeat_units": 18, "repeat_length_bp": 90},
            {"accession": "ACC002", "repeat_class": "imperfect", "repeat_units": 22, "repeat_length_bp": 110},
            {"accession": "ACC003", "repeat_class": "perfect", "repeat_units": 20, "repeat_length_bp": 100},
            {"accession": "ACC003", "repeat_class": "perfect", "repeat_units": 25, "repeat_length_bp": 125},
            {"accession": "ACC004", "repeat_class": "compound_component", "repeat_units": 15, "repeat_length_bp": 75},
            {"accession": "ACC005", "repeat_class": "perfect", "repeat_units": 20, "repeat_length_bp": 100},
        ]

        ssr_records = []
        for i, ssr_data_item in enumerate(ssr_data):
            ssr = {
                "accession": ssr_data_item["accession"],
                "start": 100 + i * 200,
                "end": 100 + i * 200 + ssr_data_item["repeat_length_bp"],
                "motif_raw": "AGC",
                "motif_canonical": "AGC",
                "repeat_units": ssr_data_item["repeat_units"],
                "motif_size": 3,
                "repeat_length_bp": ssr_data_item["repeat_length_bp"],
                "strand": "+",
                "run_id": run.run_id,
                "repeat_class": ssr_data_item["repeat_class"],
            }
            ssr_records.append(ssr)
        
        insert_ssr_records(session, ssr_records)
        session.flush()

        # Create metrics
        for acc in ["ACC001", "ACC002", "ACC003", "ACC004", "ACC005"]:
            upsert_accession_metrics(
                session,
                run_id=run.run_id,
                accession=acc,
                ra=10.5,
                rd=1000.0,
            )
        session.commit()

        return session

    def test_get_date_range_for_field(self, populated_db):
        """Query date range."""
        # Get run_id from a test SSR
        run_id = populated_db.query(SSRRecord).first().run_id
        
        min_date, max_date = get_date_range_for_field(
            populated_db, "release_date", run_id=run_id
        )
        assert min_date is not None
        assert max_date is not None
        assert min_date < max_date

    def test_query_accessions_in_time_bin(self, populated_db):
        """Query accession count in time bin."""
        # Get run_id from a test SSR
        run_id = populated_db.query(SSRRecord).first().run_id
        
        bin_start = datetime(2024, 1, 1)
        bin_end = datetime(2024, 2, 1)
        count = query_accessions_in_time_bin(
            populated_db, "release_date", bin_start, bin_end, run_id=run_id
        )
        assert count == 2  # ACC001, ACC002

    def test_query_temporal_ssr_metrics(self, populated_db):
        """Query SSR metrics for time bin."""
        # Get run_id from a test SSR
        run_id = populated_db.query(SSRRecord).first().run_id
        
        bin_start = datetime(2024, 1, 1)
        bin_end = datetime(2024, 2, 1)
        metrics = query_temporal_ssr_metrics(
            populated_db, "release_date", bin_start, bin_end, run_id=run_id
        )
        assert metrics["ssr_count_total"] == 3  # 2 from ACC001 + 1 from ACC002
        assert metrics["perfect_count"] == 2
        assert metrics["imperfect_count"] == 1

    def test_compute_temporal_summary(self, populated_db):
        """Compute full temporal summary."""
        # Get run_id from a test SSR
        run_id = populated_db.query(SSRRecord).first().run_id
        
        summary = compute_temporal_summary(
            populated_db, "release_date", DatePrecision.MONTH, run_id=run_id
        )
        assert summary.date_field == "release_date"
        assert summary.precision == DatePrecision.MONTH
        assert len(summary.bins) > 0
        assert summary.total_accessions > 0

    def test_temporal_summary_with_collection_date(self, populated_db):
        """Temporal summary using collection_date field."""
        # Get run_id from a test SSR
        run_id = populated_db.query(SSRRecord).first().run_id
        
        summary = compute_temporal_summary(
            populated_db, "collection_date", DatePrecision.MONTH, run_id=run_id
        )
        assert summary.date_field == "collection_date"
        assert len(summary.bins) > 0


# ---------------------------------------------------------------------------
# Edge case tests
# ---------------------------------------------------------------------------

class TestEdgeCases:
    """Test edge cases and error handling."""

    def test_parse_date_whitespace(self):
        """Handle whitespace in date strings."""
        dt = parse_date_string("  2024-04-25  ")
        assert dt is not None
        assert dt.year == 2024

    def test_normalize_leap_year(self):
        """Handle leap year dates correctly."""
        # Feb 29, 2024 (leap year)
        dt = datetime(2024, 2, 29)
        normalized = normalize_date_to_precision(dt, DatePrecision.MONTH)
        assert normalized.month == 2
        assert normalized.day == 1

    def test_create_bins_single_day(self):
        """Create bins for single-day range."""
        start = datetime(2024, 4, 25)
        end = datetime(2024, 4, 25)
        bins = create_time_bins(start, end, DatePrecision.DAY)
        assert len(bins) == 1

    def test_empty_temporal_summary(self, session):
        """Temporal summary with no data."""
        summary = compute_temporal_summary(session, "release_date", DatePrecision.MONTH)
        assert len(summary.bins) == 0
        assert len(summary.analysis_notes) > 0

    def test_bin_id_edge_dates(self):
        """Generate bin IDs for edge dates."""
        # Year boundary
        dt = datetime(2024, 12, 31)
        bin_id = bin_id_from_date(dt, DatePrecision.YEAR)
        assert bin_id == "2024"

        # Month boundary
        dt = datetime(2024, 1, 1)
        bin_id = bin_id_from_date(dt, DatePrecision.MONTH)
        assert bin_id == "2024-01"


# ---------------------------------------------------------------------------
# Backward compatibility tests
# ---------------------------------------------------------------------------

class TestBackwardCompatibility:
    """Test backward compatibility with existing code."""

    def test_temporal_module_imports(self):
        """Verify all public functions can be imported."""
        from gwico_ssr.analysis.temporal import (
            DatePrecision,
            TimeBin,
            TemporalSummary,
            parse_date_string,
            normalize_date_to_precision,
            bin_id_from_date,
            create_time_bins,
            compute_temporal_summary,
        )
        assert all([
            DatePrecision,
            TimeBin,
            TemporalSummary,
            parse_date_string,
            normalize_date_to_precision,
            bin_id_from_date,
            create_time_bins,
            compute_temporal_summary,
        ])

    def test_visualization_imports(self):
        """Verify visualization functions can be imported."""
        from gwico_ssr.visualization.temporal_figures import (
            plot_temporal_distribution,
            plot_temporal_ssr_trends,
            create_interactive_temporal_plot,
        )
        assert all([
            plot_temporal_distribution,
            plot_temporal_ssr_trends,
            create_interactive_temporal_plot,
        ])


# ---------------------------------------------------------------------------
# Repeat class breakdown tests
# ---------------------------------------------------------------------------

class TestRepeatClassBreakdown:
    """Test repeat-class-aware temporal analysis."""

    @pytest.fixture
    def mixed_class_db(self, session):
        """Database with perfect, imperfect, and compound SSRs."""
        # Create dataset
        dataset = get_or_create_dataset(session, "mixed_class_dataset", "mixed class test")
        session.flush()

        # Create run
        run = create_run(session, dataset_id=dataset.dataset_id, stage="detect")
        session.flush()

        # Create accessions
        for i in range(1, 4):
            upsert_accession(
                session,
                accession=f"ACC{i:03d}",
                dataset_id=dataset.dataset_id,
                release_date=f"2024-{i:02d}-15",
                collection_date=f"2024-{i:02d}-01",
            )
        session.flush()

        # SSRs with mixed classes
        ssr_data = [
            ("ACC001", "perfect", 20, 100),
            ("ACC001", "imperfect", 22, 110),
            ("ACC002", "perfect", 20, 100),
            ("ACC002", "perfect", 18, 90),
            ("ACC002", "compound_component", 15, 75),
            ("ACC003", "imperfect", 25, 125),
        ]

        ssr_records = []
        for i, (accession, repeat_class, repeat_units, repeat_bp) in enumerate(ssr_data):
            ssr = {
                "accession": accession,
                "start": 100 + i * 200,
                "end": 100 + i * 200 + repeat_bp,
                "motif_raw": "AGC",
                "motif_canonical": "AGC",
                "repeat_units": repeat_units,
                "motif_size": 3,
                "repeat_length_bp": repeat_bp,
                "strand": "+",
                "run_id": run.run_id,
                "repeat_class": repeat_class,
            }
            ssr_records.append(ssr)
        
        insert_ssr_records(session, ssr_records)
        session.flush()

        # Metrics
        for i in range(1, 4):
            upsert_accession_metrics(
                session,
                run_id=run.run_id,
                accession=f"ACC{i:03d}",
                ra=10.5,
                rd=1000.0,
            )
        session.commit()

        return session

    def test_perfect_count_breakdown(self, mixed_class_db):
        """Verify perfect_count in temporal summary."""
        # Get run_id from a test SSR
        run_id = mixed_class_db.query(SSRRecord).first().run_id
        
        summary = compute_temporal_summary(
            mixed_class_db, "release_date", DatePrecision.MONTH, run_id=run_id
        )
        total_perfect = sum(b.perfect_count for b in summary.bins)
        assert total_perfect == 3

    def test_imperfect_count_breakdown(self, mixed_class_db):
        """Verify imperfect_count in temporal summary."""
        # Get run_id from a test SSR
        run_id = mixed_class_db.query(SSRRecord).first().run_id
        
        summary = compute_temporal_summary(
            mixed_class_db, "release_date", DatePrecision.MONTH, run_id=run_id
        )
        total_imperfect = sum(b.imperfect_count for b in summary.bins)
        assert total_imperfect == 2

    def test_compound_count_breakdown(self, mixed_class_db):
        """Verify compound_component_count in temporal summary."""
        # Get run_id from a test SSR
        run_id = mixed_class_db.query(SSRRecord).first().run_id
        
        summary = compute_temporal_summary(
            mixed_class_db, "release_date", DatePrecision.MONTH, run_id=run_id
        )
        total_compound = sum(b.compound_component_count for b in summary.bins)
        assert total_compound == 1
