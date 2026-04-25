#!/usr/bin/env python
"""GWICO-SSR Large-Cohort Example Workflow.

Demonstrates production-scale batch-first workflow for 1M-accession datasets.
Includes batch size configuration, failure recovery, and performance monitoring.

Usage:
    # Run with default settings (10K batch size)
    python examples/large_cohort_example.py
    
    # Custom batch size
    python examples/large_cohort_example.py --batch-size 5000
    
    # Resume after interruption
    python examples/large_cohort_example.py --resume

Prerequisites:
    pip install -e ".[dev]"
    
Expected Output:
    - SQLite database: outputs/large_cohort/gwico.db (~150 GB for 1M accessions)
    - Browser tracks: outputs/large_cohort/browser_tracks.bed
    - Analysis results: outputs/large_cohort/analysis_results.csv
    - Manifest: outputs/large_cohort/RUN_MANIFEST.json
    
Estimated Runtime:
    - 100K accessions: 6-8 hours (8-core system)
    - 1M accessions: 40-60 hours (8-core system) or 270 hours (single-threaded)
"""

from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
import tempfile
from datetime import datetime
from pathlib import Path
from typing import Optional

# Resolve paths
REPO_ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = REPO_ROOT / "data" / "demo"
OUTPUT_DIR = REPO_ROOT / "outputs" / "large_cohort"
DB_PATH = OUTPUT_DIR / "gwico.db"
CONFIG_PATH = OUTPUT_DIR / "large_cohort.toml"
MANIFEST_PATH = OUTPUT_DIR / "RUN_MANIFEST.json"

CLI = [sys.executable, "-m", "gwico_ssr"]


class LargeCohortWorkflow:
    """Orchestrate production-scale GWICO-SSR analysis."""

    def __init__(
        self,
        batch_size: int = 10000,
        dataset_name: str = "large_cohort_demo",
        resume: bool = False,
        dry_run: bool = False,
    ):
        """Initialize workflow.

        Args:
            batch_size: Accessions per batch (1K, 5K, 10K, 20K)
            dataset_name: Name for database dataset
            resume: Resume from last completed phase
            dry_run: Print commands without executing
        """
        self.batch_size = batch_size
        self.dataset_name = dataset_name
        self.resume = resume
        self.dry_run = dry_run
        self.start_time = datetime.now()
        self.phase_times = {}

        # Validate batch size
        assert batch_size in [1000, 5000, 10000, 20000], \
            f"Batch size must be 1K, 5K, 10K, or 20K; got {batch_size}"

    def setup(self) -> None:
        """Create output directory and configuration."""
        OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
        
        # Create configuration file
        config_content = self._generate_config()
        CONFIG_PATH.write_text(config_content)
        print(f"✓ Configuration: {CONFIG_PATH}")

    def _generate_config(self) -> str:
        """Generate TOML configuration for batch-first workflow."""
        return f"""# GWICO-SSR Large-Cohort Configuration
# Generated: {datetime.now().isoformat()}
# Batch Size: {self.batch_size:,} accessions

[database]
url = "sqlite:///{DB_PATH}"
journal_mode = "WAL"
synchronous = "NORMAL"
cache_size = 2000

[download]
request_batch_size = 500
artifact_batch_size = {self.batch_size}
retry_max_attempts = 5
retry_backoff_seconds = 60

[parsing]
chunk_size = {min(self.batch_size, 1000)}
format_detection = "auto"
malformed_strategy = "skip"
feature_extraction = true

[detection]
motif_sizes = [1, 2, 3, 4, 5, 6]
perfect = true
imperfect = true
imperfect_max_mismatches = 2
compound = true
compound_dmax = 1

[logging]
level = "INFO"
format = "json"

[performance]
# Performance tuning
use_wal_mode = true
connection_pool_size = 20
batch_insert_size = 1000
"""

    def run_cmd(
        self,
        args: list[str],
        description: str,
        phase_name: Optional[str] = None,
    ) -> None:
        """Run CLI command with status reporting.

        Args:
            args: Command arguments after 'python -m gwico_ssr'
            description: Human-readable description
            phase_name: Phase name for timing tracking
        """
        cmd = CLI + args
        print(f"\n{'='*70}")
        print(f"  {description}")
        print(f"  $ {' '.join(cmd)}")
        print(f"{'='*70}")

        if self.dry_run:
            print("  [DRY RUN - Command not executed]")
            return

        # Set environment
        env = dict(os.environ)
        env["GWICO_SSR_DB_URL"] = f"sqlite:///{DB_PATH}"
        env["GWICO_SSR_CONFIG"] = str(CONFIG_PATH)

        # Run command
        start = datetime.now()
        result = subprocess.run(cmd, env=env, capture_output=True, text=True)
        elapsed = datetime.now() - start

        # Log output (truncate for large outputs)
        if result.stdout:
            lines = result.stdout.split('\n')
            if len(lines) > 20:
                print('\n'.join(lines[:10]))
                print(f"... [{len(lines)-20} lines omitted] ...")
                print('\n'.join(lines[-10:]))
            else:
                print(result.stdout[:2000])

        if result.returncode != 0:
            print(f"\n❌ FAILED with exit code {result.returncode}")
            print(f"STDERR: {result.stderr[:1000]}")
            raise RuntimeError(f"Command failed: {' '.join(cmd)}")

        # Track phase time
        if phase_name:
            self.phase_times[phase_name] = elapsed
            print(f"\n  ⏱ Phase time: {elapsed}")

        print(f"  ✅ Success")

    def phase_1_ingest(self) -> None:
        """Phase 1: Ingest metadata (0.5 hours for 1M accessions)."""
        print("\n" + "="*70)
        print("PHASE 1: INGEST METADATA")
        print("="*70)

        # For demo, use existing demo metadata
        # In production, split large metadata files by batch size
        self.run_cmd(
            [
                "ingest",
                str(DATA_DIR / "demo_metadata.csv"),
                "--dataset-name", self.dataset_name,
                "--organism", "SARS-CoV-2-demo",
            ],
            "Ingesting metadata from CSV",
            phase_name="ingest",
        )

    def phase_2_download(self) -> None:
        """Phase 2: Download sequences from NCBI (18-24 hours for 1M)."""
        print("\n" + "="*70)
        print("PHASE 2: DOWNLOAD FROM NCBI")
        print("="*70)

        args = [
            "download",
            "--dataset-name", self.dataset_name,
            "--output-dir", str(OUTPUT_DIR / "downloads"),
            "--show-progress",
        ]

        if self.resume:
            args.append("--resume")

        self.run_cmd(
            args,
            f"Downloading sequences (batch size: {self.batch_size:,})",
            phase_name="download",
        )

    def phase_3_parse(self) -> None:
        """Phase 3: Parse sequences (3-6 hours for 1M)."""
        print("\n" + "="*70)
        print("PHASE 3: PARSE SEQUENCES")
        print("="*70)

        args = [
            "parse",
            "--dataset-name", self.dataset_name,
            "--data-dir", str(OUTPUT_DIR / "downloads"),
            "--show-progress",
        ]

        if self.resume:
            args.append("--resume")

        self.run_cmd(
            args,
            f"Parsing sequences (chunk size: {self.batch_size//10})",
            phase_name="parse",
        )

    def phase_4_detect(self) -> None:
        """Phase 4: Detect SSRs (4-8 hours for 1M single-threaded)."""
        print("\n" + "="*70)
        print("PHASE 4: DETECT SSRs")
        print("="*70)

        print("  Note: Detection is CPU-intensive")
        print("  For 8-core parallelization, use GNU parallel:")
        print(f"    parallel -j 8 'python -m gwico_ssr detect ...")
        print(f"      --dataset-name {self.dataset_name} --accession-batch {{}}' \\\n"
              f"      ::: $(seq 0 100000 1000000)")

        args = [
            "detect",
            "--dataset-name", self.dataset_name,
            "--motif-sizes", "1,2,3,4,5,6",
            "--show-progress",
        ]

        if self.resume:
            args.append("--resume")

        self.run_cmd(
            args,
            "Detecting SSRs (perfect, imperfect, compound)",
            phase_name="detect",
        )

    def phase_5_annotate(self) -> None:
        """Phase 5: Map SSRs to genes (2-4 hours for 1M)."""
        print("\n" + "="*70)
        print("PHASE 5: ANNOTATE SSRs TO GENES")
        print("="*70)

        args = [
            "annotate",
            "--dataset-name", self.dataset_name,
            "--show-progress",
        ]

        if self.resume:
            args.append("--resume")

        self.run_cmd(
            args,
            "Mapping SSRs to gene features",
            phase_name="annotate",
        )

    def phase_6_analyze(self) -> None:
        """Phase 6: Statistical analysis (<1 hour for 1M)."""
        print("\n" + "="*70)
        print("PHASE 6: STATISTICAL ANALYSIS")
        print("="*70)

        self.run_cmd(
            [
                "analyze",
                "--dataset-name", self.dataset_name,
                "--analyses", "all",
                "--fdr-correction", "benjamini-hochberg",
                "--output", str(OUTPUT_DIR / "analysis"),
                "--show-progress",
            ],
            "Running comprehensive statistical analysis",
            phase_name="analyze",
        )

    def phase_7_visualize(self) -> None:
        """Phase 7: Generate publication figures (30 min for 1M)."""
        print("\n" + "="*70)
        print("PHASE 7: VISUALIZATION")
        print("="*70)

        self.run_cmd(
            [
                "visualize",
                "--dataset-name", self.dataset_name,
                "--figure-types", "all",
                "--output", str(OUTPUT_DIR / "figures"),
                "--dpi", "300",
            ],
            "Generating publication-quality figures",
            phase_name="visualize",
        )

    def phase_8_export(self) -> None:
        """Phase 8: Export formats and create manifest (30 min for 1M)."""
        print("\n" + "="*70)
        print("PHASE 8: EXPORT AND MANIFEST")
        print("="*70)

        # Export all formats
        self.run_cmd(
            [
                "export",
                "--dataset-name", self.dataset_name,
                "--export-type", "all",
                "--output", str(OUTPUT_DIR / "exports"),
            ],
            "Exporting analysis in all formats (CSV, BED, GFF3, JSON)",
            phase_name="export",
        )

        # Create manifest
        self.run_cmd(
            [
                "export",
                "--dataset-name", self.dataset_name,
                "--export-type", "manifest",
                "--output", str(MANIFEST_PATH),
            ],
            "Creating reproducibility manifest",
            phase_name="manifest",
        )

    def create_browser_tracks(self) -> None:
        """Generate UCSC Genome Browser tracks."""
        print("\n" + "="*70)
        print("GENERATING BROWSER TRACKS")
        print("="*70)

        tracks_dir = OUTPUT_DIR / "browser_tracks"
        tracks_dir.mkdir(exist_ok=True)

        self.run_cmd(
            [
                "export",
                "--dataset-name", self.dataset_name,
                "--export-type", "browser_tracks",
                "--output", str(tracks_dir),
            ],
            "Creating UCSC Genome Browser BED12 tracks",
            phase_name="browser_tracks",
        )

    def validate(self) -> None:
        """Validate outputs and create checksum file."""
        print("\n" + "="*70)
        print("VALIDATION")
        print("="*70)

        # Check key outputs exist
        required_files = [
            OUTPUT_DIR / "gwico.db",
            OUTPUT_DIR / "RUN_MANIFEST.json",
            OUTPUT_DIR / "exports" / "ssr_coordinates.bed",
        ]

        for fpath in required_files:
            if not fpath.exists():
                print(f"  ⚠ Warning: {fpath} not found")
            else:
                size = fpath.stat().st_size
                print(f"  ✓ {fpath.name}: {size / 1e9:.2f} GB")

        # Create SHA256 checksums
        print("\n  Creating SHA256 checksums...")
        import subprocess
        result = subprocess.run(
            f"cd {OUTPUT_DIR} && find . -type f -exec sha256sum {{}} > SHA256SUMS.txt \\;",
            shell=True,
            capture_output=True,
            text=True,
        )
        if result.returncode == 0:
            print("  ✓ SHA256SUMS.txt created")

    def generate_report(self) -> None:
        """Generate execution report."""
        print("\n" + "="*70)
        print("EXECUTION REPORT")
        print("="*70)

        total_time = datetime.now() - self.start_time

        report = {
            "timestamp": self.start_time.isoformat(),
            "dataset_name": self.dataset_name,
            "batch_size": self.batch_size,
            "total_runtime": str(total_time),
            "phase_times": {
                name: str(elapsed) for name, elapsed in self.phase_times.items()
            },
            "output_directory": str(OUTPUT_DIR),
            "database": str(DB_PATH),
            "config": str(CONFIG_PATH),
        }

        print(json.dumps(report, indent=2))

        # Save report
        report_path = OUTPUT_DIR / "EXECUTION_REPORT.json"
        report_path.write_text(json.dumps(report, indent=2))
        print(f"\n  ✓ Report saved to {report_path}")

    def run(self) -> None:
        """Execute full workflow."""
        try:
            self.setup()
            
            # Determine starting phase
            if self.resume:
                print(f"\n  Resuming from last checkpoint...")
            else:
                print(f"\n  Starting large-cohort workflow")
                print(f"  Batch size: {self.batch_size:,} accessions")
                print(f"  Dataset: {self.dataset_name}")
                print(f"  Output: {OUTPUT_DIR}")

            # Run phases
            self.phase_1_ingest()
            self.phase_2_download()
            self.phase_3_parse()
            self.phase_4_detect()
            self.phase_5_annotate()
            self.phase_6_analyze()
            self.phase_7_visualize()
            self.phase_8_export()
            self.create_browser_tracks()

            # Finalize
            self.validate()
            self.generate_report()

            print("\n" + "="*70)
            print("✅ WORKFLOW COMPLETE")
            print("="*70)
            print(f"\nResults available at: {OUTPUT_DIR}")
            print(f"Database: {DB_PATH}")
            print(f"Manifest: {MANIFEST_PATH}")

        except KeyboardInterrupt:
            print("\n\n⚠ Workflow interrupted (Ctrl+C)")
            print("Resume with: python examples/large_cohort_example.py --resume")
            sys.exit(1)
        except Exception as e:
            print(f"\n❌ Workflow failed: {e}")
            sys.exit(1)


def main():
    """Command-line interface."""
    parser = argparse.ArgumentParser(
        description="Large-cohort GWICO-SSR workflow",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Run with default 10K batch size
  python examples/large_cohort_example.py
  
  # Run with 5K batch size (memory-constrained)
  python examples/large_cohort_example.py --batch-size 5000
  
  # Run with 20K batch size (high-performance)
  python examples/large_cohort_example.py --batch-size 20000
  
  # Resume after interruption
  python examples/large_cohort_example.py --resume
  
  # Dry run (print commands without executing)
  python examples/large_cohort_example.py --dry-run
        """,
    )

    parser.add_argument(
        "--batch-size",
        type=int,
        default=10000,
        choices=[1000, 5000, 10000, 20000],
        help="Accessions per batch (default: 10000)",
    )

    parser.add_argument(
        "--dataset-name",
        type=str,
        default="large_cohort_demo",
        help="Dataset name (default: large_cohort_demo)",
    )

    parser.add_argument(
        "--resume",
        action="store_true",
        help="Resume from last checkpoint",
    )

    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print commands without executing",
    )

    args = parser.parse_args()

    workflow = LargeCohortWorkflow(
        batch_size=args.batch_size,
        dataset_name=args.dataset_name,
        resume=args.resume,
        dry_run=args.dry_run,
    )

    workflow.run()


if __name__ == "__main__":
    main()
