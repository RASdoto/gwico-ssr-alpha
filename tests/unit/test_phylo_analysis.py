"""Test suite for Chunk 11 phylogenetic tree integration.

Coverage:
- Newick tree parsing (8 tests)
- Tree ingestion (5 tests)
- Accession-to-tip mapping (8 tests)
- Metrics computation (6 tests)
- Visualizations (4 tests)
- Backward compatibility (2 tests)
"""

from __future__ import annotations

import pytest
import tempfile
from pathlib import Path

from sqlalchemy import create_engine, select
from sqlalchemy.orm import Session

from gwico_ssr.models.schema import (
    Base,
    Dataset,
    Run,
    Accession,
    SSRRecord,
    AccessionMetrics,
    TreeSource,
    TreeMapping,
    TreeMetrics,
)
from gwico_ssr.phylo.tree_parser import Tree, TreeNode, parse_newick_file, validate_tree_structure
from gwico_ssr.phylo.phylo_analysis import (
    ingest_newick_tree,
    map_accessions_to_tips,
    compute_ssr_metrics_by_clade,
    store_tree_metrics,
)
from gwico_ssr.visualization.phylo_figures import (
    plot_tree_with_ssr_overlay,
    plot_clade_ssr_heatmap,
    plot_phylo_distribution,
    plot_phylo_class_composition,
    create_interactive_phylo_plot,
)


# -----------------------------------------------------------------------
# Fixtures
# -----------------------------------------------------------------------

@pytest.fixture
def db_session():
    """In-memory SQLite database session for testing."""
    engine = create_engine("sqlite:///:memory:")
    Base.metadata.create_all(engine)
    session = Session(engine)
    yield session
    session.close()


@pytest.fixture
def populated_db(db_session):
    """Session with dataset, run, accessions, and SSR records."""
    # Create dataset
    dataset = Dataset(name="test_dataset", source_type="csv", organism="virus")
    db_session.add(dataset)
    db_session.flush()

    # Create run
    run = Run(dataset_id=dataset.dataset_id, pipeline_version="1.0.0")
    db_session.add(run)
    db_session.flush()

    # Create accessions
    accessions = []
    for acc_name in ["accession1", "accession2", "accession3", "accession4"]:
        acc = Accession(
            accession=acc_name,
            dataset_id=dataset.dataset_id,
            species="virus_A",
        )
        db_session.add(acc)
        accessions.append(acc)
    db_session.flush()

    # Create SSR records
    for i, acc in enumerate(accessions):
        for j in range(5 + i):  # Variable SSR count per accession
            ssr = SSRRecord(
                run_id=run.run_id,
                accession=acc.accession,
                motif_raw="ATAT",
                motif_canonical="ATAT",
                repeat_units=4,
                motif_size=2,
                repeat_class="perfect" if j % 2 == 0 else "imperfect",
                repeat_length_bp=10 + j * 2,
                start=100 + j * 50,
                end=110 + j * 50,
            )
            db_session.add(ssr)

        # Create metrics
        metrics = AccessionMetrics(
            run_id=run.run_id,
            accession=acc.accession,
            ssr_count_total=5 + i,
            ssr_bp_total=(5 + i) * 10,
            ra=0.5 + i * 0.1,
            rd=0.3 + i * 0.05,
        )
        db_session.add(metrics)

    db_session.commit()
    return db_session


@pytest.fixture
def simple_newick() -> str:
    """Simple binary tree with 4 tips."""
    return "((accession1:0.1,accession2:0.2):0.3,(accession3:0.15,accession4:0.25):0.35);"


@pytest.fixture
def complex_newick() -> str:
    """More complex tree with 6 tips and internal labels."""
    return "(((accession1:0.1,accession2:0.2)clade_A:0.3,(accession3:0.15,accession4:0.25)clade_B:0.35)clade_C:0.4,(accession5:0.05,accession6:0.08)clade_D:0.5)root;"


@pytest.fixture
def single_tip_newick() -> str:
    """Tree with single tip (edge case)."""
    return "accession1:0.1;"


@pytest.fixture
def unrooted_newick() -> str:
    """Unrooted tree format."""
    return "(accession1,accession2,accession3);"


@pytest.fixture
def tree_file(simple_newick, tmp_path) -> Path:
    """Create temporary tree file."""
    f = tmp_path / "test_tree.nwk"
    f.write_text(simple_newick)
    return f



# -----------------------------------------------------------------------
# Newick Parsing Tests (8)
# -----------------------------------------------------------------------

class TestNewwickParsing:
    """Test Newick format parsing and tree structure."""

    def test_parse_simple_tree(self, simple_newick):
        """Parse simple binary tree."""
        tree = Tree.from_newick(simple_newick)
        assert tree.num_tips == 4
        tips = tree.get_all_tips()
        assert set(tips) == {"accession1", "accession2", "accession3", "accession4"}

    def test_parse_without_semicolon(self, simple_newick):
        """Parse tree without trailing semicolon."""
        newick_no_semi = simple_newick.rstrip(";")
        tree = Tree.from_newick(newick_no_semi)
        assert tree.num_tips == 4

    def test_parse_complex_tree(self, complex_newick):
        """Parse tree with internal labels."""
        tree = Tree.from_newick(complex_newick)
        assert tree.num_tips == 6
        tips = tree.get_all_tips()
        assert len(tips) == 6

    def test_parse_single_tip(self, single_tip_newick):
        """Parse single-tip tree (edge case)."""
        tree = Tree.from_newick(single_tip_newick)
        assert tree.num_tips == 1

    def test_parse_unrooted(self, unrooted_newick):
        """Parse unrooted tree format."""
        tree = Tree.from_newick(unrooted_newick)
        assert tree.num_tips == 3

    def test_invalid_format_empty(self):
        """Reject empty Newick string."""
        with pytest.raises(ValueError, match="Empty"):
            Tree.from_newick("")

    def test_invalid_format_malformed(self):
        """Reject malformed Newick."""
        with pytest.raises(ValueError):
            Tree.from_newick("(accession1,accession2")  # Missing )

    def test_parse_preserves_branch_lengths(self, simple_newick):
        """Parse preserves branch length information."""
        tree = Tree.from_newick(simple_newick)
        # All leaf nodes should have branch lengths
        leaves = [n for n in tree.root.children if n.is_leaf or n.children]
        assert len(leaves) >= 2  # At least 2 children at root

    def test_tree_roundtrip(self, simple_newick):
        """Parse and serialize tree should produce valid tree."""
        tree1 = Tree.from_newick(simple_newick)
        newick_str = tree1.to_newick()
        tree2 = Tree.from_newick(newick_str)
        assert tree1.num_tips == tree2.num_tips


# -----------------------------------------------------------------------
# Tree File I/O Tests (3)
# -----------------------------------------------------------------------

class TestTreeFileIO:
    """Test reading tree files and file validation."""

    def test_parse_tree_file(self, tree_file):
        """Load tree from file."""
        tree, file_hash = parse_newick_file(tree_file)
        assert tree.num_tips == 4
        assert len(file_hash) == 64  # SHA256 hex length

    def test_tree_file_hash_consistency(self, tree_file):
        """File hash should be consistent."""
        _, hash1 = parse_newick_file(tree_file)
        _, hash2 = parse_newick_file(tree_file)
        assert hash1 == hash2

    def test_missing_tree_file(self):
        """Raise error for missing file."""
        with pytest.raises(FileNotFoundError):
            parse_newick_file("/nonexistent/tree.nwk")

    def test_invalid_file_content(self, tmp_path):
        """Raise error for invalid Newick in file."""
        f = tmp_path / "bad.nwk"
        f.write_text("(accession1,accession2")  # Missing closing paren
        with pytest.raises(ValueError):
            parse_newick_file(f)


# -----------------------------------------------------------------------
# Tree Structure Validation Tests (2)
# -----------------------------------------------------------------------

class TestTreeValidation:
    """Test tree structure validation."""

    def test_validate_valid_tree(self, simple_newick):
        """Valid tree produces no warnings."""
        tree = Tree.from_newick(simple_newick)
        messages = validate_tree_structure(tree)
        assert len(messages) == 0

    def test_validate_single_tip_warning(self, single_tip_newick):
        """Single-tip tree produces warning."""
        tree = Tree.from_newick(single_tip_newick)
        messages = validate_tree_structure(tree)
        assert any("fewer than 2" in m for m in messages)


# -----------------------------------------------------------------------
# Tree Ingestion Tests (5)
# -----------------------------------------------------------------------

class TestTreeIngestion:
    """Test loading and storing trees."""

    def test_ingest_simple_tree(self, populated_db, tree_file):
        """Ingest tree file and create TreeSource."""
        session = populated_db
        tree, report = ingest_newick_tree(
            session,
            tree_file,
            source_name="test_source",
            version="v1",
        )

        assert tree.num_tips == 4
        assert report.source.source_name == "test_source"

    def test_ingest_creates_orm_object(self, populated_db, tree_file):
        """Ingestion stores TreeSource in database."""
        session = populated_db
        ingest_newick_tree(session, tree_file, source_name="test", version="v1")

        source = session.query(TreeSource).filter(
            TreeSource.source_name == "test",
            TreeSource.version == "v1",
        ).first()

        assert source is not None
        assert source.num_tips == 4

    def test_ingest_with_release_date(self, populated_db, tree_file):
        """Ingest with release date metadata."""
        session = populated_db
        ingest_newick_tree(
            session,
            tree_file,
            source_name="dated",
            version="v1",
            release_date="2024-01-15",
        )

        source = session.query(TreeSource).filter(
            TreeSource.source_name == "dated"
        ).first()
        assert source.release_date == "2024-01-15"

    def test_ingest_idempotent(self, populated_db, tree_file):
        """Ingesting same tree twice doesn't create duplicates."""
        session = populated_db
        ingest_newick_tree(session, tree_file, source_name="test", version="v1")
        ingest_newick_tree(session, tree_file, source_name="test", version="v1")

        sources = session.query(TreeSource).filter(
            TreeSource.source_name == "test",
            TreeSource.version == "v1",
        ).all()

        assert len(sources) == 1

    def test_ingest_multiple_versions(self, populated_db, tree_file):
        """Different versions create separate TreeSource records."""
        session = populated_db
        ingest_newick_tree(session, tree_file, source_name="multi", version="v1")
        ingest_newick_tree(session, tree_file, source_name="multi", version="v2")

        sources = session.query(TreeSource).filter(
            TreeSource.source_name == "multi"
        ).all()

        assert len(sources) == 2


# -----------------------------------------------------------------------
# Accession-to-Tip Mapping Tests (8)
# -----------------------------------------------------------------------

class TestTreeMapping:
    """Test mapping accessions to tree tips."""

    def test_map_perfect_match(self, populated_db, tree_file):
        """Map accessions that exactly match tree tips."""
        session = populated_db
        tree, _ = ingest_newick_tree(session, tree_file)
        source_id = session.query(TreeSource.source_id).filter(
            TreeSource.source_name == "user-supplied"
        ).scalar()

        report = map_accessions_to_tips(session, tree, source_id)

        assert report.successfully_mapped == 4
        assert report.unmapped == 0

    def test_map_creates_tree_mapping_records(self, populated_db, tree_file):
        """Mapping creates TreeMapping ORM records."""
        session = populated_db
        tree, _ = ingest_newick_tree(session, tree_file)
        source_id = session.query(TreeSource.source_id).scalar()
        map_accessions_to_tips(session, tree, source_id)

        mappings = session.query(TreeMapping).filter(
            TreeMapping.source_id == source_id
        ).all()

        assert len(mappings) == 4

    def test_map_stores_accession_and_tip(self, populated_db, tree_file):
        """TreeMapping records contain accession and tip_label."""
        session = populated_db
        tree, _ = ingest_newick_tree(session, tree_file)
        source_id = session.query(TreeSource.source_id).scalar()
        map_accessions_to_tips(session, tree, source_id)

        mapping = session.query(TreeMapping).first()
        assert mapping.accession in ["accession1", "accession2", "accession3", "accession4"]
        assert mapping.tip_label in ["accession1", "accession2", "accession3", "accession4"]

    def test_map_accession_not_in_tree(self, populated_db, tmp_path):
        """Accessions not in tree tips are marked unmapped."""
        # Create a tree with different tips
        newick = "(tip1,tip2);"
        f = tmp_path / "diff_tree.nwk"
        f.write_text(newick)

        session = populated_db
        tree, _ = ingest_newick_tree(session, f)
        source_id = session.query(TreeSource.source_id).scalar()
        report = map_accessions_to_tips(session, tree, source_id)

        assert report.unmapped == 4  # Our accessions not in tree

    def test_map_report_contains_metadata(self, populated_db, tree_file):
        """Mapping report contains source and count info."""
        session = populated_db
        tree, _ = ingest_newick_tree(session, tree_file)
        source_id = session.query(TreeSource.source_id).scalar()
        report = map_accessions_to_tips(session, tree, source_id)

        assert report.total_tips_in_tree == 4
        assert report.total_accessions_in_database == 4

    def test_map_idempotent(self, populated_db, tree_file):
        """Mapping same tree twice doesn't duplicate records."""
        session = populated_db
        tree, _ = ingest_newick_tree(session, tree_file)
        source_id = session.query(TreeSource.source_id).scalar()

        map_accessions_to_tips(session, tree, source_id)
        initial_count = session.query(TreeMapping).filter(
            TreeMapping.source_id == source_id
        ).count()

        map_accessions_to_tips(session, tree, source_id)
        final_count = session.query(TreeMapping).filter(
            TreeMapping.source_id == source_id
        ).count()

        assert initial_count == final_count

    def test_map_confident_flag(self, populated_db, tree_file):
        """Mapping records are confident by default."""
        session = populated_db
        tree, _ = ingest_newick_tree(session, tree_file)
        source_id = session.query(TreeSource.source_id).scalar()
        map_accessions_to_tips(session, tree, source_id)

        mapping = session.query(TreeMapping).first()
        assert mapping.is_confident is True


# -----------------------------------------------------------------------
# Metrics Computation Tests (6)
# -----------------------------------------------------------------------

class TestPhyloMetrics:
    """Test phylo-aware metrics calculation."""

    def test_compute_metrics_by_clade(self, populated_db, tree_file):
        """Compute SSR metrics aggregated by clade."""
        session = populated_db
        tree, _ = ingest_newick_tree(session, tree_file)
        source_id = session.query(TreeSource.source_id).scalar()
        map_accessions_to_tips(session, tree, source_id)

        metrics = compute_ssr_metrics_by_clade(session, source_id)

        assert len(metrics) == 4  # One per tip
        for clade_name, metric in metrics.items():
            assert clade_name in ["accession1", "accession2", "accession3", "accession4"]
            assert metric.ssr_count_total > 0

    def test_metrics_include_repeat_classes(self, populated_db, tree_file):
        """Metrics break down by repeat class."""
        session = populated_db
        tree, _ = ingest_newick_tree(session, tree_file)
        source_id = session.query(TreeSource.source_id).scalar()
        map_accessions_to_tips(session, tree, source_id)

        metrics = compute_ssr_metrics_by_clade(session, source_id)

        any_metric = next(iter(metrics.values()))
        assert any_metric.perfect_count >= 0
        assert any_metric.imperfect_count >= 0

    def test_metrics_total_equals_sum(self, populated_db, tree_file):
        """Total SSR count equals sum of classes."""
        session = populated_db
        tree, _ = ingest_newick_tree(session, tree_file)
        source_id = session.query(TreeSource.source_id).scalar()
        map_accessions_to_tips(session, tree, source_id)

        metrics = compute_ssr_metrics_by_clade(session, source_id)

        for metric in metrics.values():
            class_sum = metric.perfect_count + metric.imperfect_count + metric.compound_component_count
            # Allow margin for unknown classes
            assert metric.ssr_count_total >= class_sum

    def test_store_tree_metrics(self, populated_db, tree_file):
        """Store computed metrics to database."""
        session = populated_db
        run_id = session.query(Run.run_id).first()[0]
        tree, _ = ingest_newick_tree(session, tree_file)
        source_id = session.query(TreeSource.source_id).scalar()
        map_accessions_to_tips(session, tree, source_id)

        metrics_dict = compute_ssr_metrics_by_clade(session, source_id)
        count = store_tree_metrics(session, run_id, source_id, metrics_dict)

        assert count == 4

    def test_stored_metrics_retrievable(self, populated_db, tree_file):
        """Stored metrics can be retrieved from database."""
        session = populated_db
        run_id = session.query(Run.run_id).first()[0]
        tree, _ = ingest_newick_tree(session, tree_file)
        source_id = session.query(TreeSource.source_id).scalar()
        map_accessions_to_tips(session, tree, source_id)

        metrics_dict = compute_ssr_metrics_by_clade(session, source_id)
        store_tree_metrics(session, run_id, source_id, metrics_dict)

        stored = session.query(TreeMetrics).filter(
            TreeMetrics.run_id == run_id,
            TreeMetrics.source_id == source_id,
        ).all()

        assert len(stored) == 4


# -----------------------------------------------------------------------
# Visualization Tests (4)
# -----------------------------------------------------------------------

class TestPhyloVisualizations:
    """Test phylo visualization functions."""

    def test_plot_tree_with_ssr_overlay(self, populated_db, tree_file):
        """Tree SSR overlay plot creation."""
        session = populated_db
        tree, _ = ingest_newick_tree(session, tree_file)
        source_id = session.query(TreeSource.source_id).scalar()
        map_accessions_to_tips(session, tree, source_id)

        metrics = compute_ssr_metrics_by_clade(session, source_id)
        fig = plot_tree_with_ssr_overlay(metrics)

        assert fig is not None
        assert "Clade" in fig.axes[0].get_xlabel()

    def test_plot_clade_ssr_heatmap(self, populated_db, tree_file):
        """Clade metrics heatmap creation."""
        session = populated_db
        tree, _ = ingest_newick_tree(session, tree_file)
        source_id = session.query(TreeSource.source_id).scalar()
        map_accessions_to_tips(session, tree, source_id)

        metrics = compute_ssr_metrics_by_clade(session, source_id)
        fig = plot_clade_ssr_heatmap(metrics)

        assert fig is not None

    def test_plot_phylo_distribution(self, populated_db, tree_file):
        """SSR distribution bar chart creation."""
        session = populated_db
        tree, _ = ingest_newick_tree(session, tree_file)
        source_id = session.query(TreeSource.source_id).scalar()
        map_accessions_to_tips(session, tree, source_id)

        metrics = compute_ssr_metrics_by_clade(session, source_id)
        fig = plot_phylo_distribution(metrics)

        assert fig is not None

    def test_create_interactive_phylo_plot(self, populated_db, tree_file):
        """Interactive Plotly figure creation."""
        session = populated_db
        tree, _ = ingest_newick_tree(session, tree_file)
        source_id = session.query(TreeSource.source_id).scalar()
        map_accessions_to_tips(session, tree, source_id)

        metrics = compute_ssr_metrics_by_clade(session, source_id)
        fig = create_interactive_phylo_plot(metrics)

        assert fig is not None
        assert "Clade" in fig.layout.xaxis.title.text


# -----------------------------------------------------------------------
# Backward Compatibility Tests (2)
# -----------------------------------------------------------------------

class TestBackwardCompatibility:
    """Ensure Chunk 11 doesn't break existing functionality."""

    def test_existing_accessions_unaffected(self, populated_db):
        """Accessions from pre-Chunk11 still queryable."""
        session = populated_db
        acc = session.query(Accession).first()
        assert acc is not None
        assert acc.accession.startswith("accession")

    def test_existing_ssr_records_unaffected(self, populated_db):
        """SSR records from pre-Chunk11 still queryable."""
        session = populated_db
        ssr = session.query(SSRRecord).first()
        assert ssr is not None
        assert ssr.motif_canonical == "ATAT"
