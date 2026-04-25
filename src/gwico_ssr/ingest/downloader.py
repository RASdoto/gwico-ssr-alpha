"""Download orchestration for GWICO-SSR NCBI acquisition.

Manages the download of FASTA and GenBank files for accessions in the database,
with deterministic file layout, skip-if-exists logic, and failure tracking.
"""

from __future__ import annotations

import hashlib
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
from gwico_ssr.models.schema import BatchArtifact, BatchNormalizedRecord

logger = logging.getLogger(__name__)


@dataclass
class DownloadSummary:
    """Summary of a download batch operation."""

    total_requested: int = 0
    already_downloaded: int = 0
    downloaded_ok: int = 0
    failed: int = 0
    request_batch_size: int = 1
    artifact_batch_size: int = 1
    raw_artifacts_written: int = 0
    failures: list[dict] = field(default_factory=list)

    def to_dict(self) -> dict:
        return {
            "total_requested": self.total_requested,
            "already_downloaded": self.already_downloaded,
            "downloaded_ok": self.downloaded_ok,
            "failed": self.failed,
            "request_batch_size": self.request_batch_size,
            "artifact_batch_size": self.artifact_batch_size,
            "raw_artifacts_written": self.raw_artifacts_written,
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


def _chunked(items: Sequence[str], size: int) -> list[list[str]]:
    if size <= 0:
        size = 1
    return [list(items[i:i + size]) for i in range(0, len(items), size)]


def _artifact_path(output_dir: str | Path, file_type: str, index: int) -> Path:
    suffix = "fasta" if file_type == "fasta" else "gb"
    return Path(output_dir) / "batches" / "raw" / file_type / f"batch_{index:06d}.{suffix}"


def _checksum_sha256(content: str) -> str:
    return hashlib.sha256(content.encode("utf-8", errors="ignore")).hexdigest()


def _persist_artifact_and_mappings(
    session: Session,
    dataset_id: int,
    run_id: Optional[int],
    file_type: str,
    source_type: str,
    artifact_path: Path,
    checksum: str,
    accession_order: list[str],
    normalized_paths: dict[str, Path],
    request_batch_size: int,
    artifact_batch_size: int,
    manifest_path: Optional[Path],
    duplicate_handling: str,
) -> None:
    artifact = BatchArtifact(
        dataset_id=dataset_id,
        run_id=run_id,
        file_type=file_type,
        source_type=source_type,
        artifact_path=str(artifact_path.resolve()),
        checksum_sha256=checksum,
        accession_count=len(accession_order),
        request_batch_size=request_batch_size,
        artifact_batch_size=artifact_batch_size,
        manifest_path=str(manifest_path.resolve()) if manifest_path else None,
        status="created",
    )
    session.add(artifact)
    session.flush()

    for idx, acc in enumerate(accession_order):
        norm = BatchNormalizedRecord(
            artifact_id=artifact.artifact_id,
            accession=acc,
            record_index=idx,
            normalized_fasta_path=str(normalized_paths[acc].resolve()) if file_type == "fasta" else None,
            normalized_genbank_path=str(normalized_paths[acc].resolve()) if file_type == "genbank" else None,
            checksum_sha256=None,
            duplicate_policy=duplicate_handling,
            parse_status="pending",
        )
        session.add(norm)
    session.flush()


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
    request_batch_size: int = 1,
    artifact_batch_size: int = 1,
    persist_composite_artifacts: bool = True,
    manifest_policy: str = "required",
    duplicate_handling: str = "error",
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
    summary = DownloadSummary(
        total_requested=len(accession_ids),
        request_batch_size=max(1, request_batch_size),
        artifact_batch_size=max(1, artifact_batch_size),
    )
    data_dir = _data_dir(output_dir)
    artifact_index = {"fasta": 0, "genbank": 0}

    pending: list[str] = []
    acc_to_dataset: dict[str, int] = {}

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
        acc_to_dataset[acc_id] = acc.dataset_id

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
        pending.append(acc_id)

    for artifact_chunk in _chunked(pending, summary.artifact_batch_size):
        chunk_results: dict[str, dict[str, EntrezResult]] = {
            "fasta": {},
            "genbank": {},
        }

        if "fasta" in file_types:
            for request_chunk in _chunked(artifact_chunk, summary.request_batch_size):
                if len(request_chunk) == 1:
                    res = client.fetch_fasta(request_chunk[0])
                    chunk_results["fasta"][request_chunk[0]] = res
                else:
                    for res in client.fetch_batch_fasta(request_chunk):
                        chunk_results["fasta"][res.accession] = res

        if "genbank" in file_types:
            for request_chunk in _chunked(artifact_chunk, summary.request_batch_size):
                if len(request_chunk) == 1:
                    res = client.fetch_genbank(request_chunk[0])
                    chunk_results["genbank"][request_chunk[0]] = res
                else:
                    for res in client.fetch_batch_genbank(request_chunk):
                        chunk_results["genbank"][res.accession] = res

        artifact_members: dict[str, list[str]] = {"fasta": [], "genbank": []}
        artifact_payloads: dict[str, list[str]] = {"fasta": [], "genbank": []}
        normalized_paths: dict[str, dict[str, Path]] = {"fasta": {}, "genbank": {}}

        for acc_id in artifact_chunk:
            fasta_ok = True
            genbank_ok = True
            fasta_file: Optional[Path] = None
            genbank_file: Optional[Path] = None
            error_msg: Optional[str] = None

            if "fasta" in file_types:
                fasta_file = _fasta_path(data_dir, acc_id)
                if not force and fasta_file.exists() and fasta_file.stat().st_size > 0:
                    logger.debug("FASTA file exists on disk: %s", fasta_file)
                else:
                    result = chunk_results["fasta"].get(acc_id)
                    if result and result.success and result.data and _validate_fasta(result.data):
                        _write_file(fasta_file, result.data)
                        artifact_members["fasta"].append(acc_id)
                        artifact_payloads["fasta"].append(result.data.strip())
                        normalized_paths["fasta"][acc_id] = fasta_file
                    else:
                        fasta_ok = False
                        error_msg = (result.error if result else "Missing FASTA batch result") or "Invalid FASTA data"

            if "genbank" in file_types:
                genbank_file = _genbank_path(data_dir, acc_id)
                if not force and genbank_file.exists() and genbank_file.stat().st_size > 0:
                    logger.debug("GenBank file exists on disk: %s", genbank_file)
                else:
                    result = chunk_results["genbank"].get(acc_id)
                    if result and result.success and result.data and _validate_genbank(result.data):
                        _write_file(genbank_file, result.data)
                        artifact_members["genbank"].append(acc_id)
                        artifact_payloads["genbank"].append(result.data.strip())
                        normalized_paths["genbank"][acc_id] = genbank_file
                    else:
                        genbank_ok = False
                        error_msg = (result.error if result else "Missing GenBank batch result") or "Invalid GenBank data"

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

        if persist_composite_artifacts:
            for ftype in ("fasta", "genbank"):
                if ftype not in file_types:
                    continue
                members = artifact_members[ftype]
                if not members:
                    continue

                artifact_index[ftype] += 1
                artifact_file = _artifact_path(output_dir, ftype, artifact_index[ftype])

                if ftype == "fasta":
                    content = "\n".join(artifact_payloads[ftype]).strip() + "\n"
                else:
                    content = "\n\n".join(artifact_payloads[ftype]).strip() + "\n"

                _write_file(artifact_file, content)
                checksum = _checksum_sha256(content)

                manifest_path: Optional[Path] = None
                if manifest_policy.lower() == "required":
                    manifest_path = artifact_file.with_suffix(artifact_file.suffix + ".manifest.json")
                    manifest_payload = {
                        "file_type": ftype,
                        "artifact_path": str(artifact_file),
                        "checksum_sha256": checksum,
                        "accessions": members,
                        "request_batch_size": summary.request_batch_size,
                        "artifact_batch_size": summary.artifact_batch_size,
                    }
                    _write_file(manifest_path, json.dumps(manifest_payload, indent=2))

                dataset_id = acc_to_dataset[members[0]]
                _persist_artifact_and_mappings(
                    session=session,
                    dataset_id=dataset_id,
                    run_id=None,
                    file_type=ftype,
                    source_type="ncbi",
                    artifact_path=artifact_file,
                    checksum=checksum,
                    accession_order=members,
                    normalized_paths=normalized_paths[ftype],
                    request_batch_size=summary.request_batch_size,
                    artifact_batch_size=summary.artifact_batch_size,
                    manifest_path=manifest_path,
                    duplicate_handling=duplicate_handling,
                )
                summary.raw_artifacts_written += 1

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
