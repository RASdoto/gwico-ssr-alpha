"""Generate gold-standard fixture files for GWICO-SSR validation.

This script is kept for reproducibility - it documents exactly how
the gold_standard.fasta and gold_standard_expected.json were created.
Run it to regenerate the fixtures if needed.
"""

from __future__ import annotations

import json
import random
from pathlib import Path

from gwico_ssr.ssr.detector import SSRThresholds, detect_ssrs


def make_background(length: int, seed: int = 42) -> str:
    """Generate a random background sequence avoiding long repeats."""
    rng = random.Random(seed)
    bases: list[str] = []
    prev = ""
    for _ in range(length):
        choices = [b for b in "ACGT" if b != prev or rng.random() > 0.3]
        b = rng.choice(choices)
        bases.append(b)
        prev = b
    return "".join(bases)


def inject_ssrs(background: str, injections: list[dict]) -> str:
    """Inject SSR motifs into background sequence with flanking isolation."""
    seq_list = list(background)
    for inj in injections:
        motif = inj["motif"]
        units = inj["repeat_units"]
        pos = inj["position"]
        repeat_seq = motif * units
        # Clear flanking bases to avoid merging
        flank_start = max(0, pos - 2)
        for i in range(flank_start, pos):
            seq_list[i] = "T" if motif[0] != "T" else "C"
        for i, b in enumerate(repeat_seq):
            seq_list[pos + i] = b
        end = pos + len(repeat_seq)
        for i in range(end, min(end + 2, len(seq_list))):
            seq_list[i] = "T" if motif[-1] != "T" else "C"
    return "".join(seq_list)


def main() -> None:
    fixture_dir = Path(__file__).parent

    # Generate two gold-standard sequences
    # Sequence 1: 500bp with 11 injected SSRs covering all motif sizes
    bg1 = make_background(500, seed=42)
    injections_1 = [
        {"motif": "A", "repeat_units": 12, "position": 10},
        {"motif": "C", "repeat_units": 15, "position": 40},
        {"motif": "AT", "repeat_units": 6, "position": 70},
        {"motif": "AG", "repeat_units": 8, "position": 100},
        {"motif": "AAG", "repeat_units": 4, "position": 130},
        {"motif": "ATG", "repeat_units": 5, "position": 160},
        {"motif": "AATG", "repeat_units": 3, "position": 200},
        {"motif": "ACTG", "repeat_units": 4, "position": 230},
        {"motif": "AACGT", "repeat_units": 3, "position": 270},
        {"motif": "AACGTT", "repeat_units": 2, "position": 320},
        {"motif": "AATCGG", "repeat_units": 3, "position": 360},
    ]
    seq1 = inject_ssrs(bg1, injections_1)

    # Sequence 2: 300bp with mono/di edge cases
    bg2 = make_background(300, seed=99)
    injections_2 = [
        {"motif": "T", "repeat_units": 20, "position": 10},
        {"motif": "GC", "repeat_units": 10, "position": 50},
        {"motif": "CAG", "repeat_units": 6, "position": 90},
    ]
    seq2 = inject_ssrs(bg2, injections_2)

    # Write FASTA
    fasta_path = fixture_dir / "gold_standard.fasta"
    lines = [
        ">GOLD_001 synthetic gold-standard SSR validation sequence 1\n",
    ]
    for i in range(0, len(seq1), 70):
        lines.append(seq1[i : i + 70] + "\n")
    lines.append(">GOLD_002 synthetic gold-standard SSR validation sequence 2\n")
    for i in range(0, len(seq2), 70):
        lines.append(seq2[i : i + 70] + "\n")
    fasta_path.write_text("".join(lines))

    # Detect SSRs and record expected results
    thresholds = SSRThresholds()
    r1 = detect_ssrs(seq1, "GOLD_001", thresholds)
    r2 = detect_ssrs(seq2, "GOLD_002", thresholds)

    expected = {
        "thresholds": {
            "mono": thresholds.mono,
            "di": thresholds.di,
            "tri": thresholds.tri,
            "tetra": thresholds.tetra,
            "penta": thresholds.penta,
            "hexa": thresholds.hexa,
        },
        "sequences": {
            "GOLD_001": {
                "length": len(seq1),
                "total_hits": r1.hit_count,
                "hits": [
                    {
                        "start": h.start,
                        "end": h.end,
                        "motif_raw": h.motif_raw,
                        "motif_canonical": h.motif_canonical,
                        "motif_size": h.motif_size,
                        "repeat_units": h.repeat_units,
                        "repeat_length_bp": h.repeat_length_bp,
                        "strand": h.strand,
                        "actual_repeat": h.actual_repeat,
                    }
                    for h in r1.hits
                ],
            },
            "GOLD_002": {
                "length": len(seq2),
                "total_hits": r2.hit_count,
                "hits": [
                    {
                        "start": h.start,
                        "end": h.end,
                        "motif_raw": h.motif_raw,
                        "motif_canonical": h.motif_canonical,
                        "motif_size": h.motif_size,
                        "repeat_units": h.repeat_units,
                        "repeat_length_bp": h.repeat_length_bp,
                        "strand": h.strand,
                        "actual_repeat": h.actual_repeat,
                    }
                    for h in r2.hits
                ],
            },
        },
    }

    json_path = fixture_dir / "gold_standard_expected.json"
    json_path.write_text(json.dumps(expected, indent=2))

    print(f"FASTA written to {fasta_path}")
    print(f"Expected JSON written to {json_path}")
    print(f"GOLD_001: {r1.hit_count} hits in {len(seq1)} bp")
    print(f"GOLD_002: {r2.hit_count} hits in {len(seq2)} bp")
    for acc, result in [("GOLD_001", r1), ("GOLD_002", r2)]:
        print(f"\n{acc}:")
        for h in result.hits:
            print(
                f"  [{h.start}:{h.end}) {h.motif_canonical} "
                f"(raw={h.motif_raw}) {h.motif_size}bp×{h.repeat_units} "
                f"strand={h.strand}"
            )


if __name__ == "__main__":
    main()
