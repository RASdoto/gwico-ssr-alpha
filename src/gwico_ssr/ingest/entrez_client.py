"""NCBI Entrez client wrapper for GWICO-SSR.

Provides rate-limited, retry-capable access to NCBI Entrez for downloading
FASTA and GenBank files. Designed for testability — all network calls go
through a single ``EntrezClient`` class that can be mocked.
"""

from __future__ import annotations

import logging
import time
from dataclasses import dataclass, field
from typing import Optional

from Bio import Entrez

logger = logging.getLogger(__name__)


@dataclass
class EntrezResult:
    """Result of a single Entrez fetch operation."""

    accession: str
    file_type: str  # "fasta" or "genbank"
    success: bool
    data: Optional[str] = None
    error: Optional[str] = None


@dataclass
class EntrezConfig:
    """Configuration for the Entrez client."""

    email: str = ""
    api_key: str = ""
    max_retries: int = 3
    rate_limit: float = 3.0  # requests per second (up to 10 with API key)
    batch_size: int = 500
    tool: str = "gwico-ssr"


class EntrezClient:
    """Rate-limited wrapper around BioPython's Entrez module.

    Handles retries with exponential backoff and enforces rate limits.
    All NCBI network I/O is centralized here for mockability.
    """

    def __init__(self, config: EntrezConfig) -> None:
        self._config = config
        self._last_request_time: float = 0.0

        # Configure BioPython Entrez
        Entrez.email = config.email or "gwico-ssr@example.com"
        if config.api_key:
            Entrez.api_key = config.api_key
        Entrez.tool = config.tool

    @property
    def min_interval(self) -> float:
        """Minimum seconds between requests."""
        if self._config.rate_limit <= 0:
            return 0.0
        return 1.0 / self._config.rate_limit

    def _throttle(self) -> None:
        """Enforce rate limiting between requests."""
        now = time.monotonic()
        elapsed = now - self._last_request_time
        wait = self.min_interval - elapsed
        if wait > 0:
            time.sleep(wait)
        self._last_request_time = time.monotonic()

    def _fetch_with_retry(
        self,
        accession: str,
        db: str,
        rettype: str,
        retmode: str,
    ) -> tuple[bool, Optional[str], Optional[str]]:
        """Fetch a single record with retries and exponential backoff.

        Returns ``(success, data, error_message)``.
        """
        max_retries = self._config.max_retries
        for attempt in range(1, max_retries + 1):
            try:
                self._throttle()
                handle = Entrez.efetch(
                    db=db,
                    id=accession,
                    rettype=rettype,
                    retmode=retmode,
                )
                data = handle.read()
                handle.close()

                if isinstance(data, bytes):
                    data = data.decode("utf-8", errors="replace")

                # Basic validity check
                if not data or len(data.strip()) == 0:
                    raise ValueError("Empty response from NCBI")

                return True, data, None

            except Exception as exc:
                error_msg = f"{type(exc).__name__}: {exc}"
                logger.warning(
                    "Entrez fetch failed (attempt %d/%d)",
                    attempt,
                    max_retries,
                    extra={
                        "accession": accession,
                        "rettype": rettype,
                        "attempt": attempt,
                        "error": error_msg,
                    },
                )
                if attempt < max_retries:
                    backoff = 2 ** (attempt - 1)  # 1s, 2s, 4s, ...
                    time.sleep(backoff)

        return False, None, f"Failed after {max_retries} attempts"

    def fetch_fasta(self, accession: str) -> EntrezResult:
        """Download a FASTA record for a single accession."""
        success, data, error = self._fetch_with_retry(
            accession=accession,
            db="nucleotide",
            rettype="fasta",
            retmode="text",
        )
        return EntrezResult(
            accession=accession,
            file_type="fasta",
            success=success,
            data=data,
            error=error,
        )

    def fetch_genbank(self, accession: str) -> EntrezResult:
        """Download a GenBank record for a single accession."""
        success, data, error = self._fetch_with_retry(
            accession=accession,
            db="nucleotide",
            rettype="gb",
            retmode="text",
        )
        return EntrezResult(
            accession=accession,
            file_type="genbank",
            success=success,
            data=data,
            error=error,
        )

    def fetch_batch_fasta(self, accessions: list[str]) -> list[EntrezResult]:
        """Download FASTA records for a batch of accessions.

        Uses a single efetch call with comma-joined IDs for efficiency.
        Falls back to individual fetches on batch failure.
        """
        if not accessions:
            return []

        results: list[EntrezResult] = []
        try:
            self._throttle()
            handle = Entrez.efetch(
                db="nucleotide",
                id=",".join(accessions),
                rettype="fasta",
                retmode="text",
            )
            data = handle.read()
            handle.close()

            if isinstance(data, bytes):
                data = data.decode("utf-8", errors="replace")

            # Split batch response into individual records
            records = _split_fasta(data)

            # Map records back to accessions
            found_ids = set()
            for record_text in records:
                # Extract accession from header line
                header_line = record_text.split("\n", 1)[0]
                rec_id = _extract_accession_from_header(header_line, accessions)
                if rec_id:
                    found_ids.add(rec_id)
                    results.append(EntrezResult(
                        accession=rec_id,
                        file_type="fasta",
                        success=True,
                        data=record_text,
                    ))

            # Mark missing accessions
            for acc in accessions:
                if acc not in found_ids:
                    results.append(EntrezResult(
                        accession=acc,
                        file_type="fasta",
                        success=False,
                        error="Not found in batch response",
                    ))

        except Exception as exc:
            logger.warning(
                "Batch FASTA fetch failed, falling back to individual fetches",
                extra={"batch_size": len(accessions), "error": str(exc)},
            )
            for acc in accessions:
                results.append(self.fetch_fasta(acc))

        return results


def _split_fasta(data: str) -> list[str]:
    """Split multi-FASTA text into individual records."""
    records: list[str] = []
    current: list[str] = []
    for line in data.split("\n"):
        if line.startswith(">"):
            if current:
                records.append("\n".join(current))
            current = [line]
        elif current:
            current.append(line)
    if current:
        records.append("\n".join(current))
    return records


def _extract_accession_from_header(
    header: str, known_accessions: list[str]
) -> Optional[str]:
    """Extract an accession ID from a FASTA header line.

    Tries to match against the known accession list first (exact or prefix),
    then falls back to extracting the first whitespace-delimited token.
    """
    if not header.startswith(">"):
        return None
    # Remove '>' prefix
    header_body = header[1:].strip()
    first_token = header_body.split()[0] if header_body else ""

    # Try exact match against known accessions
    for acc in known_accessions:
        if first_token == acc or first_token.startswith(acc):
            return acc

    return first_token if first_token else None
