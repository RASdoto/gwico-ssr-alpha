"""Compound SSR detection and standardization for GWICO-SSR.

Implements IMEX-compatible compound SSR chaining using dMAX gap distance
and motif standardization levels (L0, L1, L2, Full).

Algorithm:
    1. Sort SSRs by start position within each accession.
    2. For each SSR, attempt to chain it with the next SSR if gap <= dmax.
    3. Build maximal chains of connected SSRs.
    4. Apply standardization to each chain's motif sequence.
    5. Persist compound records and component mappings.

Standardization levels:
    L0:   Use raw motif from first component (no transformation).
    L1:   Canonicalize + rotation to lexicographic minimum.
    L2:   Canonicalize (rotation + reverse complement normalization) - DEFAULT.
    Full: All rotations + reverse complement representations (for comprehensive search).
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from typing import Optional, Sequence

from gwico_ssr.ssr.detector import SSRHit
from gwico_ssr.ssr.motif import (
    all_rotations,
    canonicalize_motif,
    determine_strand,
    reverse_complement,
)


class StandardizationLevel(str, Enum):
    """Standardization levels for compound motif representation.
    
    L0:   Raw motif from first component (no transformation).
    L1:   Rotation normalization (lexicographic minimum of forward rotations).
    L2:   Canonical (rotation + reverse complement normalization) - DEFAULT.
    Full: All representations (rotations + reverse complement rotations).
    """
    L0 = "L0"
    L1 = "L1"
    L2 = "L2"
    Full = "Full"


@dataclass
class CompoundHit:
    """Represents a detected compound SSR.
    
    Stores the compound chain with component SSRs, gap information,
    and standardized motif.
    """
    accession: str
    start: int                          # Start of first component
    end: int                            # End of last component
    components: list[SSRHit]            # Ordered list of component SSRs
    gaps_bp: list[int]                  # Gap sizes between consecutive components
    dmax_used: int                      # dMAX threshold applied
    standardization_level: StandardizationLevel
    standardized_motif: str             # Motif after standardization
    compound_motif: str                 # Composite sequence of all components
    total_repeat_length_bp: int         # Sum of repeat_length_bp for all components
    strand: str = "+"                   # Dominant strand (+ or -)
    
    @property
    def component_count(self) -> int:
        """Number of component SSRs in this compound."""
        return len(self.components)


def standardize_motif(
    motif: str,
    level: StandardizationLevel,
) -> str:
    """Apply standardization to a motif at the specified level.
    
    Args:
        motif: The input motif sequence (single repeat unit or composite).
        level: The standardization level to apply.
        
    Returns:
        Standardized motif string.
        
    Standardization rules:
        L0:   Return raw motif unchanged.
        L1:   Return lexicographically minimum rotation (no RC).
        L2:   Return canonical form (min rotation across forward + RC).
        Full: Return canonical form (same as L2, Full is for future expansion).
    """
    motif_upper = motif.upper()
    
    if level == StandardizationLevel.L0:
        return motif_upper
    
    if level == StandardizationLevel.L1:
        # Minimum rotation of forward strand only
        rotations = all_rotations(motif_upper)
        return min(rotations)
    
    if level in (StandardizationLevel.L2, StandardizationLevel.Full):
        # Full canonicalization: forward + RC
        return canonicalize_motif(motif_upper)
    
    # Fallback to L2 (canonical)
    return canonicalize_motif(motif_upper)


def find_compound_chains(
    hits: Sequence[SSRHit],
    dmax_bp: int,
) -> list[list[SSRHit]]:
    """Find maximal chains of SSRs where consecutive SSRs are <= dmax_bp apart.
    
    Algorithm:
        1. Sort hits by start position.
        2. For each hit, check if it can be chained with the next hit.
        3. If gap <= dmax, add to current chain.
        4. If gap > dmax, start a new chain.
        5. Return all chains with 2+ components.
        
    Args:
        hits: List of SSRHit objects, unsorted.
        dmax_bp: Maximum gap (in bp) between consecutive SSRs for chaining.
        
    Returns:
        List of chains; each chain is a list of SSRHit objects (2 or more).
    """
    if len(hits) < 2:
        return []
    
    # Sort by start position
    sorted_hits = sorted(hits, key=lambda h: h.start)
    
    chains: list[list[SSRHit]] = []
    current_chain: list[SSRHit] = [sorted_hits[0]]
    
    for i in range(1, len(sorted_hits)):
        prev_hit = sorted_hits[i - 1]
        curr_hit = sorted_hits[i]
        gap = curr_hit.start - prev_hit.end
        
        if gap <= dmax_bp:
            # Can chain: add to current chain
            current_chain.append(curr_hit)
        else:
            # Cannot chain: finalize current chain if 2+ components, start new
            if len(current_chain) >= 2:
                chains.append(current_chain)
            current_chain = [curr_hit]
    
    # Finalize last chain
    if len(current_chain) >= 2:
        chains.append(current_chain)
    
    return chains


def build_compound_hit(
    components: list[SSRHit],
    dmax_bp: int,
    standardization_level: StandardizationLevel,
) -> CompoundHit:
    """Build a compound hit from a chain of component SSRs.
    
    Args:
        components: List of SSRHit objects in order.
        dmax_bp: The dMAX threshold used for chaining.
        standardization_level: The standardization level to apply.
        
    Returns:
        CompoundHit with all fields populated.
    """
    # Ensure components are sorted by start
    sorted_components = sorted(components, key=lambda h: h.start)
    
    # Calculate gaps
    gaps: list[int] = []
    for i in range(len(sorted_components) - 1):
        gap = sorted_components[i + 1].start - sorted_components[i].end
        gaps.append(gap)
    
    # Extract region boundaries
    start = sorted_components[0].start
    end = sorted_components[-1].end
    accession = sorted_components[0].accession
    
    # Build composite motif by concatenating repeat sequences
    compound_motif = "".join(
        hit.actual_repeat if hit.actual_repeat else (hit.motif_raw * hit.repeat_units)
        for hit in sorted_components
    )
    
    # Total repeat length
    total_repeat_length_bp = sum(hit.repeat_length_bp for hit in sorted_components)
    
    # Determine dominant strand (majority strand of components)
    strand_counts = {"+" : 0, "-": 0}
    for hit in sorted_components:
        strand_counts[hit.strand] += 1
    dominant_strand = "+" if strand_counts["+"] >= strand_counts["-"] else "-"
    
    # Standardize the compound motif
    standardized = standardize_motif(compound_motif, standardization_level)
    
    return CompoundHit(
        accession=accession,
        start=start,
        end=end,
        components=sorted_components,
        gaps_bp=gaps,
        dmax_used=dmax_bp,
        standardization_level=standardization_level,
        standardized_motif=standardized,
        compound_motif=compound_motif,
        total_repeat_length_bp=total_repeat_length_bp,
        strand=dominant_strand,
    )


def detect_compound_ssrs(
    hits: Sequence[SSRHit],
    dmax_bp: int = 10,
    standardization_level: StandardizationLevel = StandardizationLevel.L2,
) -> list[CompoundHit]:
    """Detect compound SSRs from a list of individual SSR hits.
    
    Args:
        hits: List of SSRHit objects (can be from perfect and/or imperfect detection).
        dmax_bp: Maximum gap distance for chaining (default 10 bp, typical 5-20 bp).
        standardization_level: Level of motif standardization (default L2 - canonical).
        
    Returns:
        List of CompoundHit objects representing detected compounds.
    """
    if len(hits) < 2:
        return []
    
    # Group hits by accession
    by_accession: dict[str, list[SSRHit]] = {}
    for hit in hits:
        if hit.accession not in by_accession:
            by_accession[hit.accession] = []
        by_accession[hit.accession].append(hit)
    
    # Find chains within each accession
    compounds: list[CompoundHit] = []
    for accession, accession_hits in by_accession.items():
        chains = find_compound_chains(accession_hits, dmax_bp)
        for chain in chains:
            compound = build_compound_hit(
                chain,
                dmax_bp,
                standardization_level,
            )
            compounds.append(compound)
    
    return compounds


def mark_compound_components(
    all_hits: list[SSRHit],
    compounds: list[CompoundHit],
) -> list[SSRHit]:
    """Mark component SSRs as part of compound chains.
    
    Updates the repeat_class field to "compound_component" for any SSR
    that is part of a detected compound. Returns updated hit list.
    
    Args:
        all_hits: All SSRHit objects (perfect and imperfect).
        compounds: All detected CompoundHit objects.
        
    Returns:
        Updated list of SSRHit objects with repeat_class updated.
    """
    # Build set of component SSR identities (accession, start, end)
    component_ids = set()
    for compound in compounds:
        for component in compound.components:
            component_ids.add((component.accession, component.start, component.end))
    
    # Mark components
    updated_hits = []
    for hit in all_hits:
        if (hit.accession, hit.start, hit.end) in component_ids:
            # Create a modified copy
            hit_copy = SSRHit(
                start=hit.start,
                end=hit.end,
                motif_raw=hit.motif_raw,
                motif_canonical=hit.motif_canonical,
                motif_size=hit.motif_size,
                repeat_units=hit.repeat_units,
                repeat_length_bp=hit.repeat_length_bp,
                strand=hit.strand,
                actual_repeat=hit.actual_repeat,
                accession=hit.accession,
                is_perfect=hit.is_perfect,
                repeat_class="compound_component",  # Mark as component
                imperfection_pct=hit.imperfection_pct,
                num_substitutions=hit.num_substitutions,
                num_indels=hit.num_indels,
                imperfection_cigar=hit.imperfection_cigar,
            )
            updated_hits.append(hit_copy)
        else:
            updated_hits.append(hit)
    
    return updated_hits
