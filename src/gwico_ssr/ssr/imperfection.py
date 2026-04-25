"""Imperfection utilities for IMEX-style imperfect SSR detection.

This module implements the core algorithms for detecting imperfect SSRs:
- Substitutions (point mutations within repeat units)
- Indels (insertions/deletions at repeat unit boundaries)
- Imperfection percentage calculation
- Seed-and-extend detection strategy

References:
    Mudunuri & Nagarajaram (2007). IMEx: Search for Microsatellites in Genomes.
    Bioinformatics 23(10):1181–1187.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import Optional, Tuple, List

logger = logging.getLogger(__name__)


@dataclass
class IndelEvent:
    """A single insertion or deletion event within a repeat tract."""
    position: int           # Position in the sequence where indel occurs
    event_type: str         # "insertion" or "deletion"
    length_bp: int          # Number of base pairs involved
    sequence: Optional[str] = None  # Actual inserted/deleted bases


@dataclass
class AlignmentResult:
    """Result of aligning a motif to a sequence region with tolerance for imperfections."""
    
    motif: str                              # The expected repeat motif
    aligned_region: str                     # The actual sequence region that aligned
    start: int                              # Start position in sequence
    end: int                                # End position in sequence (exclusive)
    motif_repetitions: int                  # Number of times motif repeats
    substitutions: int                      # Count of point mutations
    indels: List[IndelEvent]                # List of indel events
    imperfection_pct: float                 # Calculated imperfection percentage
    cigar: Optional[str] = None             # CIGAR string representation
    
    @property
    def total_indel_count(self) -> int:
        """Total number of indel events."""
        return len(self.indels)
    
    @property
    def tract_length_bp(self) -> int:
        """Length of the aligned tract in base pairs."""
        return self.end - self.start


def calculate_imperfection_pct(
    substitutions: int,
    indels: int,
    tract_length_bp: int
) -> float:
    """Calculate imperfection percentage using IMEX formula.
    
    Formula: p = (S + I) / L × 100
    Where:
        S = total substitution count
        I = total indel event count
        L = tract length in bp
    
    Args:
        substitutions: Number of point mutations detected
        indels: Number of indel events detected
        tract_length_bp: Total length of the tract in base pairs
    
    Returns:
        Imperfection percentage (0-100), or 0.0 if tract_length_bp is 0
    """
    if tract_length_bp == 0:
        return 0.0
    
    total_errors = substitutions + indels
    return (total_errors / tract_length_bp) * 100.0


def count_substitutions_in_alignment(
    motif: str,
    aligned_region: str,
    start_pos: int = 0
) -> Tuple[int, str]:
    """Count point mutations when motif is aligned to a region.
    
    Assumes aligned_region is a multiple of motif length.
    
    Args:
        motif: The repeat unit to check against
        aligned_region: The actual sequence region
        start_pos: Starting position for CIGAR calculation
    
    Returns:
        Tuple of (substitution_count, cigar_string)
    """
    sub_count = 0
    cigar_parts = []
    
    motif_len = len(motif)
    region_len = len(aligned_region)
    
    # Process repeat units
    num_units = region_len // motif_len
    
    for unit_idx in range(num_units):
        motif_start = unit_idx * motif_len
        motif_end = motif_start + motif_len
        
        for pos_in_unit in range(motif_len):
            motif_base = motif[pos_in_unit]
            region_base = aligned_region[motif_start + pos_in_unit]
            
            if motif_base != region_base:
                sub_count += 1
                cigar_parts.append("M")  # Mismatch
            else:
                cigar_parts.append("=")  # Match
    
    # Handle leftover bases if not exact multiple
    if region_len % motif_len != 0:
        leftover_start = num_units * motif_len
        for pos in range(leftover_start, region_len):
            pos_in_motif = pos % motif_len
            motif_base = motif[pos_in_motif]
            region_base = aligned_region[pos]
            if motif_base != region_base:
                sub_count += 1
                cigar_parts.append("M")
            else:
                cigar_parts.append("=")
    
    cigar = "".join(cigar_parts)
    return sub_count, cigar


def find_indels(
    motif: str,
    region: str,
    max_indel_size: int = 2
) -> Tuple[List[IndelEvent], str]:
    """Detect insertion/deletion events in a region.
    
    Uses a simple heuristic: look for length mismatches between expected
    motif repetitions and actual region. For each missing or extra base,
    record an indel event.
    
    Args:
        motif: The repeat unit
        region: The sequence region to analyze
        max_indel_size: Maximum indel size to consider
    
    Returns:
        Tuple of (list of IndelEvent, cigar_with_indels)
    """
    indels: List[IndelEvent] = []
    motif_len = len(motif)
    region_len = len(region)

    # Check if region length is not a clean multiple of motif
    remainder = region_len % motif_len

    if remainder != 0:
        # Determine if it's a deletion or insertion
        # Deletion: region is shorter than expected (missing bases)
        # Insertion: region is longer than expected (extra bases)
        
        # If remainder > motif_len/2, it's closer to the next full motif,
        # so it's a deletion (missing bases from the next unit)
        # Otherwise, it's an insertion (extra bases in current units)
        
        if remainder > motif_len // 2:
            # More likely a deletion (region is short of the next full motif)
            deletion_size = motif_len - remainder
            event = IndelEvent(
                position=region_len - remainder,
                event_type="deletion",
                length_bp=deletion_size,
                sequence=None  # Don't include sequence for deletion
            )
        else:
            # More likely an insertion (extra bases beyond complete units)
            event = IndelEvent(
                position=region_len - remainder,
                event_type="insertion",
                length_bp=remainder,
                sequence=region[region_len - remainder:region_len] if remainder <= max_indel_size else None
            )
        
        if event.length_bp <= max_indel_size:
            indels.append(event)
    
    # CIGAR for indels (simplified)
    if indels:
        cigar = "I" if indels[0].event_type == "insertion" else "D"
    else:
        cigar = ""
    
    return indels, cigar


def align_motif_to_region(
    motif: str,
    region: str,
    max_substitutions: int = 1,
    max_imperfection_pct: float = 5.0,
    start_pos: int = 0,
    max_indel_size: int = 2
) -> Optional[AlignmentResult]:
    """Attempt to align a motif to a region with tolerance for imperfections.
    
    This is the core algorithm for matching imperfect repeats. It checks if
    the region matches the motif with acceptable substitutions and indels.
    
    Args:
        motif: The repeat unit (expected pattern)
        region: The actual sequence region to match
        max_substitutions: Maximum substitutions tolerated
        max_imperfection_pct: Maximum imperfection % tolerated
        start_pos: Starting position in the sequence (for tracking)
        max_indel_size: Maximum indel size to tolerate
    
    Returns:
        AlignmentResult if alignment succeeds, None if it exceeds thresholds
    """
    if len(region) == 0 or len(motif) == 0:
        return None
    
    # Count substitutions
    subs, cigar_subs = count_substitutions_in_alignment(motif, region, start_pos)
    
    # Count indels
    indels, cigar_indels = find_indels(motif, region, max_indel_size)
    
    # Calculate metrics
    indel_count = len(indels)
    tract_length = len(region)
    imperfection_pct = calculate_imperfection_pct(subs, indel_count, tract_length)
    
    # Check thresholds
    if subs > max_substitutions:
        return None
    
    if imperfection_pct > max_imperfection_pct:
        return None
    
    if indel_count > 0 and indels[0].length_bp > max_indel_size:
        return None
    
    # Build CIGAR
    motif_reps = len(region) // len(motif)
    cigar = cigar_subs + cigar_indels if cigar_indels else cigar_subs
    
    return AlignmentResult(
        motif=motif,
        aligned_region=region,
        start=start_pos,
        end=start_pos + tract_length,
        motif_repetitions=motif_reps,
        substitutions=subs,
        indels=indels,
        imperfection_pct=imperfection_pct,
        cigar=cigar if cigar else None
    )


def generate_alignment_text(
    motif: str,
    aligned_region: str,
    substitutions: int,
    indels: List[IndelEvent],
    width: int = 50
) -> str:
    """Generate human-readable alignment text for display.
    
    Args:
        motif: Expected motif
        aligned_region: Actual sequence
        substitutions: Number of substitutions
        indels: List of indel events
        width: Line width for formatting
    
    Returns:
        Multi-line alignment text
    """
    lines = []
    lines.append(f"Motif:   {motif} (repeating)")
    lines.append(f"Aligned: {aligned_region}")
    lines.append(f"Substitutions: {substitutions}, Indels: {len(indels)}")
    
    if indels:
        for event in indels:
            lines.append(f"  - {event.event_type.capitalize()} at pos {event.position}: {event.length_bp}bp")
    
    return "\n".join(lines)
