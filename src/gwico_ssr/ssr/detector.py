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
from gwico_ssr.ssr.imperfection import (
    align_motif_to_region,
    calculate_imperfection_pct,
)

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
    # Chunk 5: IMEX detector fields (optional, scaffold only)
    is_perfect: bool = True
    repeat_class: str = "perfect"  # perfect, imperfect, compound_component
    imperfection_pct: Optional[float] = None
    num_substitutions: int = 0
    num_indels: int = 0
    imperfection_cigar: Optional[str] = None
    # Chunk 7: Accession for compound detection
    accession: str = "unknown"  # Default for backward compatibility


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
                # Chunk 5: Include imperfection metadata if present
                "is_perfect": h.is_perfect,
                "repeat_class": h.repeat_class,
                "imperfection_pct": h.imperfection_pct,
                "num_substitutions": h.num_substitutions,
                "num_indels": h.num_indels,
                "imperfection_cigar": h.imperfection_cigar,
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


@dataclass
class DetectorConfig:
    """Chunk 5: Configuration for IMEX-compatible SSR detection modes.
    
    Supports perfect-only (alpha), imperfect (Chunk 6), and compound (Chunk 7) modes.
    """

    detector_mode: str = "perfect"  # perfect, imperfect, compound
    imperfection_threshold_pct: float = 5.0
    indel_max_size: int = 1
    compound_dmax_bp: int = 10
    standardization_level: str = "L2"  # L0, L1, L2, Full
    thresholds: Optional[SSRThresholds] = None

    def __post_init__(self):
        """Validate configuration."""
        valid_modes = ("perfect", "imperfect", "compound")
        if self.detector_mode not in valid_modes:
            raise ValueError(
                f"detector_mode must be one of {valid_modes}, got '{self.detector_mode}'"
            )

        valid_levels = ("L0", "L1", "L2", "Full")
        if self.standardization_level not in valid_levels:
            raise ValueError(
                f"standardization_level must be one of {valid_levels}, got '{self.standardization_level}'"
            )

        if self.thresholds is None:
            self.thresholds = SSRThresholds()

    @classmethod
    def from_settings(cls, detector_settings, ssr_settings) -> "DetectorConfig":
        """Create config from DetectorSettings and SSRSettings."""
        return cls(
            detector_mode=detector_settings.detector_mode,
            imperfection_threshold_pct=detector_settings.imperfection_threshold_pct,
            indel_max_size=detector_settings.indel_max_size,
            compound_dmax_bp=detector_settings.compound_dmax_bp,
            standardization_level=detector_settings.standardization_level,
            thresholds=SSRThresholds.from_config(ssr_settings),
        )


@dataclass
class ImperfectionMetadata:
    """Chunk 5: Metadata for imperfect SSR detections (scaffolding for Chunk 6).
    
    Not used in perfect-only mode but available for imperfect detection.
    """

    is_perfect: bool = True
    imperfection_pct: float = 0.0
    num_substitutions: int = 0
    num_indels: int = 0
    imperfection_cigar: Optional[str] = None  # CIGAR-like alignment string


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


def _set_accession_on_hits(hits: list[SSRHit], accession: str) -> list[SSRHit]:
    """Populate accession field on all hits."""
    for hit in hits:
        hit.accession = accession
    return hits


# ============================================================================
# Chunk 6: Imperfect SSR Detection (IMEX-style seed-and-extend strategy)
# ============================================================================


@dataclass
class ImperfectionThresholds:
    """Chunk 6: Per-motif-size thresholds for imperfect SSR detection."""

    motif_size: int
    max_mismatches_per_unit: int = 1  # Max substitutions in any single unit
    max_imperfection_pct: float = 5.0  # Max aggregate imperfection %
    min_repeat_units: int = 2         # Minimum repeat units to qualify
    max_indel_size: int = 2           # Max indel size to tolerate


def _seed_and_extend(
    sequence: str,
    motif: str,
    seed_start: int,
    seed_end: int,
    imperfection_thresholds: ImperfectionThresholds,
) -> Optional[SSRHit]:
    """Seed-and-extend algorithm: extend from a perfect seed while tolerating imperfections.

    Given a perfect motif seed within the sequence, this function attempts to extend
    left and right, accepting substitutions and small indels as long as they remain
    within configured thresholds.

    Args:
        sequence: Full DNA sequence
        motif: The repeat motif (should be perfect in the seed region)
        seed_start: Start position of the perfect seed
        seed_end: End position of the perfect seed (exclusive)
        imperfection_thresholds: Parameters controlling extension tolerance

    Returns:
        SSRHit if successful extension found, None otherwise
    """
    motif_len = len(motif)
    extended_start = seed_start
    extended_end = seed_end

    # Extend leftward
    while extended_start > 0:
        # Try to match the previous motif_len bases
        candidate_start = max(0, extended_start - motif_len)
        candidate = sequence[candidate_start:extended_start]

        if len(candidate) < motif_len:
            # Reached sequence boundary
            break

        # Try to align with tolerance
        alignment = align_motif_to_region(
            motif=motif,
            region=candidate,
            max_substitutions=imperfection_thresholds.max_mismatches_per_unit,
            max_imperfection_pct=imperfection_thresholds.max_imperfection_pct,
            start_pos=candidate_start,
            max_indel_size=imperfection_thresholds.max_indel_size,
        )

        if alignment:
            extended_start = candidate_start
        else:
            break

    # Extend rightward
    seq_len = len(sequence)
    while extended_end < seq_len:
        # Try to match the next motif_len bases
        candidate_end = min(seq_len, extended_end + motif_len)
        candidate = sequence[extended_end:candidate_end]

        if len(candidate) < motif_len:
            # Reached sequence boundary or insufficient bases
            break

        # Try to align with tolerance
        alignment = align_motif_to_region(
            motif=motif,
            region=candidate,
            max_substitutions=imperfection_thresholds.max_mismatches_per_unit,
            max_imperfection_pct=imperfection_thresholds.max_imperfection_pct,
            start_pos=extended_end,
            max_indel_size=imperfection_thresholds.max_indel_size,
        )

        if alignment:
            extended_end = candidate_end
        else:
            break

    # Check if extension resulted in enough repeats
    extended_length = extended_end - extended_start
    num_units = extended_length / motif_len
    if num_units < imperfection_thresholds.min_repeat_units:
        return None

    # Calculate final imperfection metrics
    extended_region = sequence[extended_start:extended_end]
    total_subs = 0
    total_indels = 0

    # Count all imperfections in the extended region
    full_alignment = align_motif_to_region(
        motif=motif,
        region=extended_region,
        max_substitutions=imperfection_thresholds.max_mismatches_per_unit * int(num_units),
        max_imperfection_pct=100.0,  # Don't filter; we'll check later
        start_pos=extended_start,
        max_indel_size=imperfection_thresholds.max_indel_size,
    )

    if full_alignment:
        total_subs = full_alignment.substitutions
        total_indels = full_alignment.total_indel_count
        imperfection_pct = full_alignment.imperfection_pct
        cigar = full_alignment.cigar
    else:
        imperfection_pct = 0.0
        cigar = None

    # Check against final imperfection % threshold
    if imperfection_pct > imperfection_thresholds.max_imperfection_pct:
        return None

    # Create SSRHit
    canonical = canonicalize_motif(motif)
    strand = determine_strand(motif, canonical)

    return SSRHit(
        start=extended_start,
        end=extended_end,
        motif_raw=motif,
        motif_canonical=canonical,
        motif_size=len(motif),
        repeat_units=int(num_units),
        repeat_length_bp=extended_length,
        strand=strand,
        actual_repeat=extended_region,
        is_perfect=(total_subs == 0 and total_indels == 0),
        repeat_class="imperfect" if (total_subs > 0 or total_indels > 0) else "perfect",
        imperfection_pct=imperfection_pct if imperfection_pct > 0 else None,
        num_substitutions=total_subs,
        num_indels=total_indels,
        imperfection_cigar=cigar,
    )


def _detect_imperfect_for_motif_size(
    sequence: str,
    motif_size: int,
    min_repeats: int,
    imperfection_pct_threshold: float,
    max_indel_size: int,
) -> list[SSRHit]:
    """Scan for imperfect SSRs of a given motif size using seed-and-extend.

    Strategy:
    1. Find all perfect seeds of the motif size with at least min_repeats
    2. For each seed, attempt seed-and-extend with imperfection tolerance
    3. Report results that exceed min_repeats and imperfection % threshold

    Args:
        sequence: Uppercase DNA sequence
        motif_size: Size of motif (1-6)
        min_repeats: Minimum repeat units
        imperfection_pct_threshold: Maximum imperfection % allowed
        max_indel_size: Maximum indel size

    Returns:
        List of imperfect SSRHit objects found
    """
    hits: list[SSRHit] = []
    seq_len = len(sequence)

    if seq_len < motif_size * min_repeats:
        return hits

    # First, find perfect seeds
    perfect_seeds = _scan_for_motif_size(sequence, motif_size, min_repeats)

    # For each perfect seed, try to extend with imperfection tolerance
    for seed in perfect_seeds:
        motif = seed.motif_raw
        thresholds = ImperfectionThresholds(
            motif_size=motif_size,
            max_mismatches_per_unit=1,  # Allow 1 sub per unit
            max_imperfection_pct=imperfection_pct_threshold,
            min_repeat_units=min_repeats,
            max_indel_size=max_indel_size,
        )

        extended = _seed_and_extend(
            sequence=sequence,
            motif=motif,
            seed_start=seed.start,
            seed_end=seed.end,
            imperfection_thresholds=thresholds,
        )

        if extended and extended.imperfection_pct is not None and extended.imperfection_pct > 0:
            # This is an imperfect SSR (extended beyond perfect seed)
            hits.append(extended)
        elif extended and (extended.imperfection_pct is None or extended.imperfection_pct == 0):
            # Perfect SSR found; will be deduplicated with perfect detector results
            pass

    return hits


def detect_imperfect_ssrs(
    sequence: str,
    accession: str,
    thresholds: Optional[SSRThresholds] = None,
    imperfection_pct_threshold: float = 5.0,
    max_indel_size: int = 2,
) -> DetectionResult:
    """Detect both perfect and imperfect SSRs in a DNA sequence.

    This function scans for imperfect SSRs in addition to perfect ones.
    Uses seed-and-extend strategy from IMEX algorithm.

    Args:
        sequence: DNA sequence string
        accession: Accession identifier
        thresholds: Minimum repeat thresholds (perfect mode)
        imperfection_pct_threshold: Maximum imperfection % for imperfect mode
        max_indel_size: Maximum indel size to tolerate

    Returns:
        DetectionResult with both perfect and imperfect SSR hits
    """
    if thresholds is None:
        thresholds = SSRThresholds()

    seq_upper = sequence.upper()
    result = DetectionResult(
        accession=accession,
        sequence_length=len(seq_upper),
    )

    all_hits: list[SSRHit] = []

    # Detect perfect SSRs first (serves as seeds)
    for motif_size in range(1, 7):
        min_reps = thresholds.for_size(motif_size)
        perfect_hits = _scan_for_motif_size(seq_upper, motif_size, min_reps)
        all_hits.extend(perfect_hits)

    # Detect imperfect SSRs (extend from perfect seeds)
    for motif_size in range(1, 7):
        min_reps = thresholds.for_size(motif_size)
        imperfect_hits = _detect_imperfect_for_motif_size(
            seq_upper,
            motif_size,
            min_reps,
            imperfection_pct_threshold,
            max_indel_size,
        )
        all_hits.extend(imperfect_hits)

    # Resolve overlaps, preferring perfect over imperfect at same position
    result.hits = _resolve_overlaps(all_hits)
    result.hits = _set_accession_on_hits(result.hits, accession)

    logger.debug(
        "Imperfect SSR detection: %s (%d bp) → %d SSRs (perfect + imperfect)",
        accession, len(seq_upper), result.hit_count,
    )

    return result


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
    result.hits = _set_accession_on_hits(result.hits, accession)

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


# ============================================================================
# Chunk 7: Compound SSR Detection (IMEX-style dMAX chaining)
# ============================================================================


def apply_compound_detection_to_result(
    result: DetectionResult,
    dmax_bp: int = 10,
    standardization_level: str = "L2",
) -> tuple[list[SSRHit], list]:
    """Apply compound SSR detection to a DetectionResult and update hit classifications.

    Chains nearby SSRs together using dMAX gap distance, then marks component SSRs
    as "compound_component" and returns the compound records.

    Args:
        result: DetectionResult with individual SSR hits.
        dmax_bp: Maximum gap distance for chaining (default 10 bp).
        standardization_level: Motif standardization level (L0, L1, L2, Full).

    Returns:
        Tuple of:
        - Updated SSRHit list with repeat_class = "compound_component" for components
        - List of CompoundHit objects representing detected compounds
    """
    from gwico_ssr.ssr.compound import (
        detect_compound_ssrs,
        mark_compound_components,
        StandardizationLevel,
    )

    if len(result.hits) < 2:
        return result.hits, []

    # Detect compounds
    std_level = StandardizationLevel(standardization_level)
    compounds = detect_compound_ssrs(result.hits, dmax_bp, std_level)

    if not compounds:
        return result.hits, []

    # Mark component SSRs
    updated_hits = mark_compound_components(result.hits, compounds)

    return updated_hits, compounds
