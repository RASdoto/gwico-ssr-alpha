"""Metadata normalization routines for GWICO-SSR ingestion.

Handles date parsing, completeness flag extraction, country/region splitting,
and other field-level normalization for accession metadata.
"""

from __future__ import annotations

import re
from typing import Optional


# ---------------------------------------------------------------------------
# Date normalization
# ---------------------------------------------------------------------------

_ISO_FULL_RE = re.compile(
    r"^(\d{4})-(\d{2})-(\d{2})"  # YYYY-MM-DD
    r"(?:T\d{2}:\d{2}:\d{2}Z?)?$"  # optional T00:00:00Z
)
_YEAR_MONTH_RE = re.compile(r"^(\d{4})-(\d{2})$")
_YEAR_ONLY_RE = re.compile(r"^(\d{4})$")


def normalize_date(raw: Optional[str]) -> Optional[str]:
    """Normalize a date string to a consistent format.

    Accepts: ``2020-01-13T00:00:00Z``, ``2020-01-13``, ``2020-01``, ``2020``.
    Returns:  ``YYYY-MM-DD``, ``YYYY-MM``, or ``YYYY`` depending on precision.
    Returns ``None`` for empty or unparseable values.
    """
    if not raw or not raw.strip():
        return None
    raw = raw.strip()

    m = _ISO_FULL_RE.match(raw)
    if m:
        return f"{m.group(1)}-{m.group(2)}-{m.group(3)}"

    m = _YEAR_MONTH_RE.match(raw)
    if m:
        return f"{m.group(1)}-{m.group(2)}"

    m = _YEAR_ONLY_RE.match(raw)
    if m:
        return m.group(1)

    return None


# ---------------------------------------------------------------------------
# Completeness normalization
# ---------------------------------------------------------------------------

_COMPLETE_VALUES = frozenset({"complete", "complete genome"})


def normalize_completeness(raw: Optional[str]) -> Optional[bool]:
    """Convert nucleotide completeness string to a boolean.

    ``"complete"`` or ``"complete genome"`` → ``True``.
    Any other non-empty value → ``False``.
    Empty / missing → ``None``.
    """
    if not raw or not raw.strip():
        return None
    val = raw.strip().lower()
    return val in _COMPLETE_VALUES


# ---------------------------------------------------------------------------
# Geo-location normalization
# ---------------------------------------------------------------------------

def normalize_geo_location(raw: Optional[str]) -> tuple[Optional[str], Optional[str]]:
    """Parse geo_location into (country, region).

    Expected formats:
    - ``"China"`` → ``("China", None)``
    - ``"USA: California"`` → ``("USA", "California")``
    - ``"South Korea: Seoul"`` → ``("South Korea", "Seoul")``

    Returns ``(None, None)`` for empty / missing values.
    """
    if not raw or not raw.strip():
        return None, None
    raw = raw.strip()

    if ":" in raw:
        parts = raw.split(":", 1)
        country = parts[0].strip() or None
        region = parts[1].strip() or None
        return country, region

    return raw, None


# ---------------------------------------------------------------------------
# Genome length normalization
# ---------------------------------------------------------------------------

def normalize_length(raw: Optional[str]) -> Optional[int]:
    """Convert length field to an integer.

    Returns ``None`` for empty, non-numeric, or negative values.
    """
    if not raw or not str(raw).strip():
        return None
    try:
        val = int(str(raw).strip())
        return val if val > 0 else None
    except (ValueError, TypeError):
        return None


# ---------------------------------------------------------------------------
# Accession normalization
# ---------------------------------------------------------------------------

_ACCESSION_RE = re.compile(r"^[A-Za-z0-9_.]+$")


def normalize_accession(raw: Optional[str]) -> Optional[str]:
    """Validate and normalize an accession identifier.

    Strips whitespace and verifies format.
    Returns ``None`` for empty or invalid accession strings.
    """
    if not raw or not raw.strip():
        return None
    val = raw.strip()
    if not _ACCESSION_RE.match(val):
        return None
    return val
