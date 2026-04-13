"""Download orchestration for GWICO-SSR NCBI acquisition.

Manages the download of FASTA and GenBank files for accessions in the database,
with deterministic file layout, skip-if-exists logic, and failure tracking.
"""

from __future__ import annotations

import json
import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional, Sequence

from sqlalchemy.orm import Session

from gwico_ssr.db.repository import (
    get_accession,
    list_accessions,
    upsert_sequence_record,
)
from gwico_ssr.ingest.entrez_client import EntrezClient, EntrezConfig, EntrezResult

logger = logging.getLogger(__name__)


@dataclass
class DownloadSummary:
    """Summary of a download batch operation."""

    total_requested: int = 0
    already_downloaded: int = 0
    downloaded_ok: int = 0
    failed: int = 0
    failures: list[dict] = field(default_factory=list)

    def to_dict(self) -> dict:
        return {
            "total_requested": self.total_requested,
            "already_downloaded": self.already_downloaded,
            "downloaded_ok": self.downloaded_ok,
            "failed": self.failed,
            "failures": self.failures,
        }


def _data_dir(output_base: str | Path) -> Path:
    """Return the deterministic data directory for downloaded files."""
    return Path(output_base) / "sequences"


def _fasta_path(data_dir: Path, accession: str) -> Path:
    """Deterministic path for a FASTA file."""
    return data_dir / "fasta" / f"{accession}.fasta"


def _genbank_path(data_dir: Path, accession: str) -> Path:
    """Deterministic path for a GenBank file."""
    return data_dir / "genbank" / f"{accession}.gb"


def _write_file(path: Path, data: str) -> None:
    """Write data to a file, creating parent directories as needed."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(data, encoding="utf-8")


def _validate_fasta(data: str) -> bool:
    """Basic sanity check for FASTA data."""
    stripped = data.strip()
    return stripped.startswith(">") and len(stripped) > 10


def _validate_genbank(data: str) -> bool:
    """Basic sanity check for GenBank data."""
    stripped = data.strip()
    return stripped.startswith("LOCUS") and len(stripped) > 50


def download_accessions(
    session: Session,
    client: EntrezClient,
    accession_ids: Sequence[str],
    output_dir: str | Path,
    file_types: Sequence[str] = ("fasta", "genbank"),
    force: bool = False,
) -> DownloadSummary:
    """Download FASTA and/or GenBank files for a list of accessions.

    Args:
        session: Active database session.
        client: Configured EntrezClient.
        accession_ids: Accession IDs to download.
        output_dir: Base output directory.
        file_types: Which file types to download ("fasta", "genbank", or both).
        force: If True, re-download even if files already exist.

    Returns:
        DownloadSummary with counts and failure details.
    """
    summary = DownloadSummary(total_requested=len(accession_ids))
    data_dir = _data_dir(output_dir)

    for acc_id in accession_ids:
        # Verify accession exists in DB
        acc = get_accession(session, acc_id)
        if acc is None:
            summary.failed += 1
            summary.failures.append({
                "accession": acc_id,
                "error": "Accession not found in database",
            })
            continue

        # Check if already downloaded (unless force)
        if not force:
            seq_rec = acc.sequence_record
            if seq_rec and seq_rec.download_status == "downloaded":
                # Verify files still exist
                files_ok = True
                if "fasta" in file_types and seq_rec.fasta_path:
                    if not Path(seq_rec.fasta_path).exists():
                        files_ok = False
                if "genbank" in file_types and seq_rec.genbank_path:
                    if not Path(seq_rec.genbank_path).exists():
                        files_ok = False
                if files_ok:
                    summary.already_downloaded += 1
                    logger.debug("Skipping %s (already downloaded)", acc_id)
                    continue

        # Download requested file types
        fasta_ok = True
        genbank_ok = True
        fasta_file: Optional[Path] = None
        genbank_file: Optional[Path] = None
        error_msg: Optional[str] = None

        if "fasta" in file_types:
            fasta_file = _fasta_path(data_dir, acc_id)

            # Skip download if file exists on disk and not forcing
            if not force and fasta_file.exists() and fasta_file.stat().st_size > 0:
                logger.debug("FASTA file exists on disk: %s", fasta_file)
            else:
                result = client.fetch_fasta(acc_id)
                if result.success and result.data and _validate_fasta(result.data):
                    _write_file(fasta_file, result.data)
                else:
                    fasta_ok = False
                    error_msg = result.error or "Invalid FASTA data"

        if "genbank" in file_types:
            genbank_file = _genbank_path(data_dir, acc_id)

            if not force and genbank_file.exists() and genbank_file.stat().st_size > 0:
                logger.debug("GenBank file exists on disk: %s", genbank_file)
            else:
                result = client.fetch_genbank(acc_id)
                if result.success and result.data and _validate_genbank(result.data):
                    _write_file(genbank_file, result.data)
                else:
                    genbank_ok = False
                    error_msg = result.error or "Invalid GenBank data"

        # Update sequence record
        if fasta_ok and genbank_ok:
            sr_kwargs: dict = {
                "accession": acc_id,
                "sequence_source": "ncbi",
                "download_status": "downloaded",
            }
            if fasta_file:
                sr_kwargs["fasta_path"] = str(fasta_file.resolve())
            if genbank_file:
                sr_kwargs["genbank_path"] = str(genbank_file.resolve())
            upsert_sequence_record(session, **sr_kwargs)
            summary.downloaded_ok += 1
            logger.info("Downloaded %s", acc_id)
        else:
            upsert_sequence_record(
                session,
                accession=acc_id,
                sequence_source="ncbi",
                download_status="failed",
            )
            summary.failed += 1
            summary.failures.append({
                "accession": acc_id,
                "error": error_msg or "Unknown download failure",
            })
            logger.warning("Download failed for %s: %s", acc_id, error_msg)

    return summary


def write_retry_manifest(
    summary: DownloadSummary,
    output_dir: str | Path,
) -> Optional[Path]:
    """Write a JSON manifest of failed accessions for retry.

    Returns the manifest path, or None if there are no failures.
    """
    if not summary.failures:
        return None

    manifest_path = Path(output_dir) / "retry_manifest.json"
    manifest_path.parent.mkdir(parents=True, exist_ok=True)

    manifest = {
        "failed_count": summary.failed,
        "accessions": [f["accession"] for f in summary.failures],
        "details": summary.failures,
    }
    manifest_path.write_text(
        json.dumps(manifest, indent=2), encoding="utf-8"
    )
    logger.info("Retry manifest written to %s", manifest_path)
    return manifest_path


def load_retry_manifest(manifest_path: str | Path) -> list[str]:
    """Load accession IDs from a retry manifest file."""
    path = Path(manifest_path)
    if not path.exists():
        return []
    data = json.loads(path.read_text(encoding="utf-8"))
    return data.get("accessions", [])
