"""Perfect SSR detection engine for GWICO-SSR.

Scans DNA sequences for perfect Simple Sequence Repeats (microsatellites)
with motif sizes 1-6 bp. Uses configurable minimum repeat thresholds.

Algorithm:
    For each motif_size (1 to 6):
        Slide along the sequence. At each position, check if the next
        `motif_size` bases repeat contiguously. Count consecutive repeat
        units. If the count meets the minimum threshold, emit an SSR record.

Overlap precedence:
    When a region is detected as an SSR at multiple motif sizes, longer
    motifs take priority. For example, AAGAAG detected as both a 3-mer
    repeat (AAG×2) and a 6-mer (AAGAAG×1): the 3-mer representation is
    preferred because the 6-mer is a sub-repeat. We filter out SSRs whose
    motif is a sub-repeat of a shorter motif.

Coordinates:
    0-based half-open [start, end), consistent with internal conventions.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import Optional, Sequence

from gwico_ssr.ssr.motif import canonicalize_motif, determine_strand, is_sub_repeat

logger = logging.getLogger(__name__)

# Detector version string embedded in SSR records
DETECTOR_VERSION = "gwico-ssr-1.0"


@dataclass
class SSRHit:
    """A single SSR detection result."""

    start: int          # 0-based inclusive
    end: int            # 0-based exclusive
    motif_raw: str      # The actual motif found in the sequence
    motif_canonical: str
    motif_size: int     # 1-6
    repeat_units: int
    repeat_length_bp: int
    strand: str         # "+" or "-"
    actual_repeat: str  # The actual repeated sequence from the genome


@dataclass
class DetectionResult:
    """Result of SSR detection on a single sequence."""

    accession: str
    sequence_length: int
    hits: list[SSRHit] = field(default_factory=list)
    detector_version: str = DETECTOR_VERSION

    @property
    def hit_count(self) -> int:
        return len(self.hits)

    def count_by_motif_size(self) -> dict[int, int]:
        """Count hits grouped by motif size."""
        counts: dict[int, int] = {}
        for h in self.hits:
            counts[h.motif_size] = counts.get(h.motif_size, 0) + 1
        return counts

    def to_dicts(self) -> list[dict]:
        """Convert hits to dicts suitable for insert_ssr_records."""
        return [
            {
                "accession": self.accession,
                "start": h.start,
                "end": h.end,
                "motif_raw": h.motif_raw,
                "motif_canonical": h.motif_canonical,
                "motif_size": h.motif_size,
                "repeat_units": h.repeat_units,
                "repeat_length_bp": h.repeat_length_bp,
                "strand": h.strand,
                "actual_repeat": h.actual_repeat,
                "detector_version": self.detector_version,
            }
            for h in self.hits
        ]


@dataclass
class SSRThresholds:
    """Minimum repeat counts per motif size."""

    mono: int = 10
    di: int = 5
    tri: int = 3
    tetra: int = 3
    penta: int = 3
    hexa: int = 2

    def for_size(self, motif_size: int) -> int:
        return {
            1: self.mono,
            2: self.di,
            3: self.tri,
            4: self.tetra,
            5: self.penta,
            6: self.hexa,
        }[motif_size]

    @classmethod
    def from_config(cls, ssr_settings) -> SSRThresholds:
        """Create thresholds from SSRSettings config object."""
        return cls(
            mono=ssr_settings.min_repeats_mono,
            di=ssr_settings.min_repeats_di,
            tri=ssr_settings.min_repeats_tri,
            tetra=ssr_settings.min_repeats_tetra,
            penta=ssr_settings.min_repeats_penta,
            hexa=ssr_settings.min_repeats_hexa,
        )


def _scan_for_motif_size(
    sequence: str,
    motif_size: int,
    min_repeats: int,
) -> list[SSRHit]:
    """Scan a sequence for perfect SSRs of a given motif size.

    Args:
        sequence: Uppercase DNA sequence.
        motif_size: Size of the motif (1-6).
        min_repeats: Minimum number of consecutive repeats.

    Returns:
        List of SSRHit objects found.
    """
    hits: list[SSRHit] = []
    seq_len = len(sequence)
    if seq_len < motif_size * min_repeats:
        return hits

    i = 0
    while i <= seq_len - motif_size:
        motif = sequence[i : i + motif_size]

        # Skip motifs containing ambiguous bases
        if not all(c in "ACGT" for c in motif):
            i += 1
            continue

        # Count consecutive repeats of this motif
        repeat_count = 1
        j = i + motif_size
        while j + motif_size <= seq_len and sequence[j : j + motif_size] == motif:
            repeat_count += 1
            j += motif_size

        if repeat_count >= min_repeats:
            # Check if this motif is a sub-repeat of a shorter motif
            # If so, skip it — the shorter motif scan will catch it
            if not is_sub_repeat(motif):
                canonical = canonicalize_motif(motif)
                strand = determine_strand(motif, canonical)
                end = i + motif_size * repeat_count
                actual = sequence[i:end]

                hits.append(SSRHit(
                    start=i,
                    end=end,
                    motif_raw=motif,
                    motif_canonical=canonical,
                    motif_size=motif_size,
                    repeat_units=repeat_count,
                    repeat_length_bp=end - i,
                    strand=strand,
                    actual_repeat=actual,
                ))

            # Advance past this repeat region
            i = j
        else:
            i += 1

    return hits


def _resolve_overlaps(hits: list[SSRHit]) -> list[SSRHit]:
    """Resolve overlapping SSR detections.

    When the same genomic region is detected as SSRs at different motif sizes,
    we prefer the shorter motif size (the more fundamental repeat unit).
    For example, AAGAAG as a trimer (AAG×2) is preferred over as a hexamer.

    The sub-repeat filter in _scan_for_motif_size already handles most cases.
    This function handles remaining overlaps by keeping the hit with the
    smaller motif_size, or if equal, the longer repeat.
    """
    if not hits:
        return hits

    # Sort by start position, then by motif_size (prefer shorter motifs)
    hits.sort(key=lambda h: (h.start, h.motif_size, -h.repeat_length_bp))

    resolved: list[SSRHit] = []
    for hit in hits:
        # Check if this hit overlaps with the last resolved hit
        if resolved and hit.start < resolved[-1].end:
            prev = resolved[-1]
            # If the new hit starts at the same position, prefer shorter motif
            if hit.start == prev.start:
                if hit.motif_size < prev.motif_size:
                    resolved[-1] = hit
                elif hit.motif_size == prev.motif_size and hit.repeat_length_bp > prev.repeat_length_bp:
                    resolved[-1] = hit
            # Otherwise, if it partially overlaps, keep both
            # (different SSRs can overlap at boundaries)
            else:
                resolved.append(hit)
        else:
            resolved.append(hit)

    return resolved


def detect_ssrs(
    sequence: str,
    accession: str,
    thresholds: Optional[SSRThresholds] = None,
) -> DetectionResult:
    """Detect perfect SSRs in a DNA sequence.

    Scans for all motif sizes 1-6 and resolves overlaps.

    Args:
        sequence: DNA sequence string.
        accession: Accession identifier.
        thresholds: Minimum repeat thresholds. Uses defaults if None.

    Returns:
        DetectionResult with all SSR hits.
    """
    if thresholds is None:
        thresholds = SSRThresholds()

    seq_upper = sequence.upper()
    result = DetectionResult(
        accession=accession,
        sequence_length=len(seq_upper),
    )

    all_hits: list[SSRHit] = []

    for motif_size in range(1, 7):
        min_reps = thresholds.for_size(motif_size)
        hits = _scan_for_motif_size(seq_upper, motif_size, min_reps)
        all_hits.extend(hits)

    result.hits = _resolve_overlaps(all_hits)

    logger.debug(
        "SSR detection: %s (%d bp) → %d SSRs",
        accession, len(seq_upper), result.hit_count,
    )

    return result


def detect_ssrs_from_fasta(
    fasta_path: str,
    thresholds: Optional[SSRThresholds] = None,
) -> list[DetectionResult]:
    """Detect SSRs in all records of a FASTA file.

    Args:
        fasta_path: Path to FASTA file.
        thresholds: Minimum repeat thresholds.

    Returns:
        List of DetectionResult, one per sequence record.
    """
    from gwico_ssr.parsers.fasta_parser import parse_fasta

    parse_result = parse_fasta(fasta_path)
    results: list[DetectionResult] = []

    for rec in parse_result.records:
        dr = detect_ssrs(rec.sequence, rec.accession, thresholds)
        results.append(dr)

    return results
