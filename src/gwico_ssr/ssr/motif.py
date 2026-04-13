"""Motif canonicalization utilities for GWICO-SSR.

Provides canonical motif representation using lexicographic minimum rotation
and reverse complement normalization. This ensures equivalent motifs
(e.g., AAG, AGA, GAA and their reverse complements CTT, TTC, TCT) all
map to a single canonical form.

Convention:
    canonical = min(all rotations of motif, all rotations of reverse complement)
    choosing the lexicographically smallest rotation across both strands.
"""

from __future__ import annotations

from functools import lru_cache

# Complement table for DNA bases
_COMPLEMENT = str.maketrans("ACGTacgt", "TGCAtgca")


def reverse_complement(seq: str) -> str:
    """Return the reverse complement of a DNA sequence."""
    return seq.translate(_COMPLEMENT)[::-1]


def all_rotations(motif: str) -> list[str]:
    """Return all rotations of a motif string.

    For "AAG" → ["AAG", "AGA", "GAA"].
    """
    n = len(motif)
    doubled = motif + motif
    return [doubled[i : i + n] for i in range(n)]


@lru_cache(maxsize=8192)
def canonicalize_motif(motif: str) -> str:
    """Return the canonical form of a motif.

    The canonical form is the lexicographically smallest rotation
    across both the forward motif and its reverse complement.

    Examples:
        canonicalize_motif("AAG") → "AAG"
        canonicalize_motif("AGA") → "AAG"
        canonicalize_motif("GAA") → "AAG"
        canonicalize_motif("CTT") → "AAG"  (RC of AAG)
        canonicalize_motif("TTC") → "AAG"  (RC of GAA = CTT, rotation)
    """
    motif_upper = motif.upper()
    rc = reverse_complement(motif_upper)

    forward_rotations = all_rotations(motif_upper)
    rc_rotations = all_rotations(rc)

    return min(forward_rotations + rc_rotations)


def determine_strand(motif_raw: str, motif_canonical: str) -> str:
    """Determine strand based on whether the canonical form matches
    a rotation of the raw motif (+ strand) or its reverse complement (- strand).

    Returns "+" or "-".
    """
    motif_upper = motif_raw.upper()
    if motif_upper in all_rotations(motif_canonical):
        return "+"
    return "-"


@lru_cache(maxsize=8192)
def is_sub_repeat(motif: str) -> bool:
    """Check if a motif is a sub-repeat of a shorter motif.

    E.g., "AAGAAG" is a sub-repeat of "AAG" (hexamer made of 2× trimer).
    "AAAA" is a sub-repeat of "A".
    "ATATAT" is a sub-repeat of "AT".

    We check whether the motif can be constructed by repeating a shorter
    sub-motif. Used to avoid double-counting SSRs at multiple motif sizes.
    """
    motif_upper = motif.upper()
    n = len(motif_upper)
    for sub_len in range(1, n):
        if n % sub_len == 0:
            sub = motif_upper[:sub_len]
            if sub * (n // sub_len) == motif_upper:
                return True
    return False
