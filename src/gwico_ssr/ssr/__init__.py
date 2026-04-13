"""SSR detection engine for GWICO-SSR."""

from gwico_ssr.ssr.detector import (
    DETECTOR_VERSION,
    DetectionResult,
    SSRHit,
    SSRThresholds,
    detect_ssrs,
    detect_ssrs_from_fasta,
)
from gwico_ssr.ssr.motif import (
    all_rotations,
    canonicalize_motif,
    determine_strand,
    is_sub_repeat,
    reverse_complement,
)

__all__ = [
    "DETECTOR_VERSION",
    "DetectionResult",
    "SSRHit",
    "SSRThresholds",
    "all_rotations",
    "canonicalize_motif",
    "detect_ssrs",
    "detect_ssrs_from_fasta",
    "determine_strand",
    "is_sub_repeat",
    "reverse_complement",
]