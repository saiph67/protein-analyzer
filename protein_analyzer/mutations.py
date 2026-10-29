"""Single-residue substitutions: what changes when one amino acid is swapped for another."""
from __future__ import annotations

import re

from .analyzer import AMINO_ACIDS, KYTE_DOOLITTLE, analyze, clean_sequence
from .compare import blosum62

AA_NAMES = {
    "A": "alanine", "R": "arginine", "N": "asparagine", "D": "aspartate", "C": "cysteine", "Q": "glutamine",
    "E": "glutamate", "G": "glycine", "H": "histidine", "I": "isoleucine", "L": "leucine", "K": "lysine",
    "M": "methionine", "F": "phenylalanine", "P": "proline", "S": "serine", "T": "threonine", "W": "tryptophan",
    "Y": "tyrosine", "V": "valine",
}
AA_GROUP = {**{c: "hydrophobic" for c in "AVILMFWC"}, **{c: "polar" for c in "STNQYG"},
            **{c: "positive" for c in "KRH"}, **{c: "negative" for c in "DE"}, "P": "special (proline)"}
_MUT_RE = re.compile(r"^([A-Z])(\d+)([A-Z])$")


def parse_mutation(text: str) -> tuple[str, int, str]:
    """'E7V' -> ('E', 7, 'V')."""
    m = _MUT_RE.match((text or "").strip().upper())
    if not m:
        raise ValueError("Write a mutation like E7V: original residue, position, new residue.")
    return m.group(1), int(m.group(2)), m.group(3)


def _local_hydropathy(seq: str, pos: int, half: int = 9) -> float:
    win = [KYTE_DOOLITTLE.get(c, 0.0) for c in seq[max(0, pos - 1 - half): pos + half]]
    return sum(win) / len(win) if win else 0.0


def mutate(seq: str, position: int, new: str, expected: str | None = None) -> dict:
    """Apply one substitution and report how the protein's properties change."""
    s = clean_sequence(seq)
    new = (new or "").strip().upper()
    if new not in AMINO_ACIDS:
        raise ValueError(f"'{new}' is not one of the 20 standard amino acids.")
    if not 1 <= position <= len(s):
        raise ValueError(f"Position {position} is outside the sequence (1–{len(s)}).")
    old = s[position - 1]
    if expected and expected.upper() != old:
        raise ValueError(f"Residue {position} is {old} ({AA_NAMES.get(old, old)}), not {expected.upper()}. "
                         "Positions count from 1 at the start (methionine) of this sequence.")
    if new == old:
        raise ValueError(f"Residue {position} is already {old}.")
    mutant = s[:position - 1] + new + s[position:]
    a, b = analyze(s, "wild type"), analyze(mutant, "mutant")
    keys = ["molecular_weight", "gravy", "isoelectric_point", "charge_ph7", "aromaticity", "n_tm_segments"]
    score = blosum62(old, new)
    return {
        "mutation": f"{old}{position}{new}",
        "position": position, "from": old, "to": new,
        "from_name": AA_NAMES.get(old, old), "to_name": AA_NAMES.get(new, new),
        "from_group": AA_GROUP.get(old, "other"), "to_group": AA_GROUP.get(new, "other"),
        "blosum62": score,
        "blosum_meaning": ("common, usually tolerated swap" if score > 0 else
                           "neutral / occasional swap" if score == 0 else "rare swap in evolution, often disruptive"),
        "hydropathy_change": round(KYTE_DOOLITTLE.get(new, 0) - KYTE_DOOLITTLE.get(old, 0), 2),
        "local_hydropathy": {"before": round(_local_hydropathy(s, position), 3),
                             "after": round(_local_hydropathy(mutant, position), 3)},
        "before": {k: a[k] for k in keys}, "after": {k: b[k] for k in keys},
        "delta": {k: round(b[k] - a[k], 4) for k in keys},
        "context": {"start": max(1, position - 10), "before": s[max(0, position - 11): position + 10],
                    "after": mutant[max(0, position - 11): position + 10]},
        "notes": _notes(old, new, score),
    }


def _notes(old: str, new: str, score: int) -> list[str]:
    notes = []
    g1, g2 = AA_GROUP.get(old), AA_GROUP.get(new)
    if g1 != g2:
        notes.append(f"Changes chemistry: {g1} → {g2}.")
    if {g1, g2} == {"positive", "negative"}:
        notes.append("Reverses the charge at this position (+ ↔ −).")
    if g1 in ("positive", "negative") and g2 == "hydrophobic":
        notes.append("Replaces a charged residue with a water-repelling one, which can make proteins stick "
                     "together (this is what happens in sickle-cell haemoglobin, E7V).")
    if new == "P":
        notes.append("Proline is rigid and breaks α-helices and β-strands.")
    if old == "P":
        notes.append("Removing a proline can make a loop more flexible.")
    if new == "G" or old == "G":
        notes.append("Glycine is the most flexible residue; adding or removing it changes backbone flexibility.")
    if old == "C":
        notes.append("Cysteine may form a disulfide bond; losing it can destabilise the fold.")
    if new == "C":
        notes.append("A new cysteine can form unwanted disulfide bonds.")
    if score <= -2 and not notes:
        notes.append("Very different residues; substitutions like this are rare between related proteins.")
    return notes


def missense_matrix(scores: dict[str, float], sequence: str) -> dict:
    """AlphaMissense scores as a dense [position][aa] matrix (None where missing) for a heatmap."""
    mat = [[None] * 20 for _ in sequence]
    idx = {a: i for i, a in enumerate(AMINO_ACIDS)}
    for key, v in scores.items():
        m = _MUT_RE.match(key)
        if not m:
            continue
        pos, to = int(m.group(2)), m.group(3)
        if 1 <= pos <= len(sequence) and to in idx and sequence[pos - 1] == m.group(1):
            mat[pos - 1][idx[to]] = round(v, 3)
    covered = sum(1 for row in mat for v in row if v is not None)
    return {"amino_acids": AMINO_ACIDS, "matrix": mat, "covered": covered,
            "mean_by_position": [round(sum(x for x in row if x is not None) / max(1, sum(x is not None for x in row)), 3)
                                 if any(x is not None for x in row) else None for row in mat]}
