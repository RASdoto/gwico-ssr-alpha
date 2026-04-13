#!/usr/bin/env python
"""GWICO-SSR End-to-End Demo Run.

Demonstrates a complete pipeline execution using the bundled demo dataset.
This script is the reference for reproducing manuscript-support artifacts.

Usage:
    python examples/demo_run.py

Prerequisites:
    pip install -e ".[dev]"
"""

from __future__ import annotations

import json
import shutil
import subprocess
import sys
from pathlib import Path

# Resolve paths
REPO_ROOT = Path(__file__).resolve().parent.parent
DEMO_DIR = REPO_ROOT / "data" / "demo"
OUTPUT_DIR = REPO_ROOT / "outputs" / "demo_run"
DB_PATH = OUTPUT_DIR / "demo.db"
CONFIG_PATH = REPO_ROOT / "config" / "default.toml"

CLI = [sys.executable, "-m", "gwico_ssr"]


def run_cmd(args: list[str], description: str) -> None:
    """Run a CLI command and print status."""
    cmd = CLI + args
    print(f"\n{'='*60}")
    print(f"  {description}")
    print(f"  $ {' '.join(cmd)}")
    print(f"{'='*60}")
    env = {
        **dict(__import__("os").environ),
        "GWICO_SSR_DB_URL": f"sqlite:///{DB_PATH}",
    }
    result = subprocess.run(cmd, env=env, capture_output=True, text=True)
    if result.stdout:
        print(result.stdout[:2000])
    if result.returncode != 0:
        print(f"STDERR: {result.stderr[:1000]}")
        print(f"Command failed with exit code {result.returncode}")
        sys.exit(1)


def main() -> None:
    print("GWICO-SSR End-to-End Demo Run")
    print(f"Output directory: {OUTPUT_DIR}")

    # Clean previous run
    if OUTPUT_DIR.exists():
        shutil.rmtree(OUTPUT_DIR)
    OUTPUT_DIR.mkdir(parents=True)

    # Step 1: Initialize database
    run_cmd(
        ["init-db"],
        "Step 1: Initialize the database",
    )

    # Step 2: Ingest metadata
    run_cmd(
        ["ingest", str(DEMO_DIR / "demo_metadata.csv"),
         "--dataset-name", "demo", "--json-summary"],
        "Step 2: Ingest accession metadata from CSV",
    )

    # Step 3: Parse sequences and annotations
    run_cmd(
        ["parse", "demo",
         "--data-dir", str(DEMO_DIR),
         "--json-summary"],
        "Step 3: Parse FASTA sequences and GenBank annotations",
    )

    # Step 4: Detect SSRs
    run_cmd(
        ["detect", "demo",
         "--data-dir", str(DEMO_DIR),
         "--json-summary"],
        "Step 4: Detect SSRs in parsed sequences",
    )

    # Step 5: Annotate SSRs against gene features
    run_cmd(
        ["annotate", "demo", "--json-summary"],
        "Step 5: Map SSRs to genomic features",
    )

    # Step 6: Compute metrics
    run_cmd(
        ["metrics", "demo", "--json-summary"],
        "Step 6: Compute per-accession SSR metrics",
    )

    # Step 7: Run statistical analyses
    run_cmd(
        ["analyze", "demo", "--json-summary"],
        "Step 7: Run statistical analyses with FDR correction",
    )

    # Step 8: Generate visualizations
    run_cmd(
        ["visualize", "demo",
         "--output-dir", str(OUTPUT_DIR / "figures"),
         "--format", "png"],
        "Step 8: Generate publication-ready figures",
    )

    # Step 9: Export results
    run_cmd(
        ["export", "demo",
         "--output-dir", str(OUTPUT_DIR / "exports"),
         "--formats", "csv,json,bed,gff3,tables"],
        "Step 9: Export results in all formats",
    )

    # Summary
    print(f"\n{'='*60}")
    print("  Demo run complete!")
    print(f"{'='*60}")
    print(f"Database: {DB_PATH}")
    print(f"Figures:  {OUTPUT_DIR / 'figures'}")
    print(f"Exports:  {OUTPUT_DIR / 'exports'}")
    print()
    print("Next steps:")
    print("  - Inspect the database: python -m gwico_ssr info")
    print("  - View figures in outputs/demo_run/figures/")
    print("  - Check exports in outputs/demo_run/exports/")


if __name__ == "__main__":
    main()
