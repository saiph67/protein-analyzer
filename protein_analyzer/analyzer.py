"""Sequence-level physicochemical analysis of proteins."""
from __future__ import annotations

import re

AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"
# Non-standard / ambiguous letters we accept (counted under "other").
EXTRA_LETTERS = "BZXJUO"
VALID_LETTERS = set(AMINO_ACIDS + EXTRA_LETTERS)
MAX_SEQUENCE_LENGTH = 35_000

# Kyte & Doolittle (1982) hydropathy scale.
KYTE_DOOLITTLE = {
    "A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5, "Q": -3.5, "E": -3.5,
    "G": -0.4, "H": -3.2, "I": 4.5, "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8,
    "P": -1.6, "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2,
}

# Average residue masses (amino acid minus water), Da.
RESIDUE_MASS = {
    "A": 71.0788, "R": 156.1875, "N": 114.1038, "D": 115.0886, "C": 103.1388,
    "E": 129.1155, "Q": 128.1307, "G": 57.0519, "H": 137.1411, "I": 113.1594,
    "L": 113.1594, "K": 128.1741, "M": 131.1926, "F": 147.1766, "P": 97.1167,
    "S": 87.0782, "T": 101.1051, "W": 186.2132, "Y": 163.1760, "V": 99.1326,
    "U": 150.0388, "O": 237.3018,
    "B": (114.1038 + 115.0886) / 2, "Z": (128.1307 + 129.1155) / 2,
    "J": 113.1594, "X": 110.0,
}
WATER = 18.01528

# pK values of Bjellqvist et al. (1993, 1994), as used by ExPASy "Compute pI/Mw" and ProtParam, so the
# app's pI and charge match the numbers students will cross-check against.
PKA_POS = {"Nterm": 7.5, "K": 10.0, "R": 12.0, "H": 5.98}
PKA_NEG = {"Cterm": 3.55, "D": 4.05, "E": 4.45, "C": 9.0, "Y": 10.0}
PKA_NTERM_BY_RESIDUE = {"A": 7.59, "M": 7.0, "S": 6.93, "P": 8.36, "T": 6.82, "V": 7.44, "E": 7.7}
# Note: ExPASy applies the residue-specific N-terminal pK but NOT a C-terminal D/E adjustment (verified against
# the live Compute pI/Mw tool, 2026-10). Biopython applies both, so it differs from ExPASy for acidic C-termini.
PI_METHOD = "Bjellqvist pK values, as ExPASy Compute pI/Mw and ProtParam"

TM_WINDOW = 19
TM_THRESHOLD = 1.6


class SequenceError(ValueError):
    """Raised for empty, too long or otherwise invalid sequences."""


def clean_sequence(seq: str) -> str:
    """Uppercase, strip whitespace/digits and a trailing stop '*'; validate letters."""
    if seq is None:
        raise SequenceError("No sequence given.")
    s = re.sub(r"[\s\d]+", "", str(seq)).upper().rstrip("*")
    if not s:
        raise SequenceError("Sequence is empty.")
    bad = sorted({c for c in s if c not in VALID_LETTERS})
    if bad:
        positions = [i + 1 for i, c in enumerate(s) if c not in VALID_LETTERS][:5]
        raise SequenceError(
            f"Invalid residue(s) {', '.join(repr(b) for b in bad)} at position(s) "
            f"{', '.join(map(str, positions))}{'...' if len(positions) == 5 else ''}. "
            f"Allowed: the 20 standard amino acids plus {EXTRA_LETTERS}."
        )
    if len(s) > MAX_SEQUENCE_LENGTH:
        raise SequenceError(f"Sequence is {len(s):,} aa; the maximum is {MAX_SEQUENCE_LENGTH:,}.")
    return s


def read_fasta(text: str) -> list[tuple[str, str]]:
    """Parse FASTA text into (name, sequence) pairs. A bare sequence (no '>') is accepted."""
    records: list[tuple[str, str]] = []
    name, chunks = None, []
    text = (text or "").replace("\r", "")
    if not text.strip():
        return records
    if not text.lstrip().startswith(">"):
        return [("sequence_1", "".join(text.split()))]
    for line in text.split("\n"):
        line = line.strip()
        if not line or line.startswith(";"):
            continue
        if line.startswith(">"):
            if name is not None:
                records.append((name, "".join(chunks)))
            name, chunks = line[1:].strip() or f"sequence_{len(records) + 1}", []
        else:
            chunks.append(line)
    if name is not None:
        records.append((name, "".join(chunks)))
    return records


def molecular_weight(seq: str) -> float:
    return sum(RESIDUE_MASS[c] for c in seq) + WATER


def net_charge(seq: str, ph: float = 7.0) -> float:
    """Net charge at a given pH (Henderson-Hasselbalch with Bjellqvist pK values)."""
    if not seq:
        return 0.0
    counts = {c: seq.count(c) for c in "KRHDECY"}
    pk_n = PKA_NTERM_BY_RESIDUE.get(seq[0], PKA_POS["Nterm"])
    pos = 1 / (1 + 10 ** (ph - pk_n))
    neg = 1 / (1 + 10 ** (PKA_NEG["Cterm"] - ph))
    for aa in "KRH":
        pos += counts[aa] / (1 + 10 ** (ph - PKA_POS[aa]))
    for aa in "DECY":
        neg += counts[aa] / (1 + 10 ** (PKA_NEG[aa] - ph))
    return pos - neg


def isoelectric_point(seq: str) -> float:
    """pH at which the net charge is zero (bisection)."""
    lo, hi = 0.0, 14.0
    for _ in range(60):
        mid = (lo + hi) / 2
        if net_charge(seq, mid) > 0:
            lo = mid
        else:
            hi = mid
    return round((lo + hi) / 2, 3)


def hydropathy_profile(seq: str, window: int = TM_WINDOW) -> list[float]:
    """Sliding-window mean KD hydropathy; value i is centred on residue i + window//2."""
    vals = [KYTE_DOOLITTLE.get(c, 0.0) for c in seq]
    if len(vals) < window:
        return []
    out, s = [], sum(vals[:window])
    out.append(s / window)
    for i in range(window, len(vals)):
        s += vals[i] - vals[i - window]
        out.append(s / window)
    return out


def tm_segments(profile: list[float], window: int = TM_WINDOW, threshold: float = TM_THRESHOLD):
    """Predicted TM segments as 1-based inclusive (start, end) residue ranges."""
    half = window // 2
    segs: list[list[int]] = []
    for i, v in enumerate(profile):
        if v >= threshold:
            centre = i + half  # 0-based residue index
            start, end = centre - half + 1, centre + half + 1
            if segs and start <= segs[-1][1] + 1:
                segs[-1][1] = max(segs[-1][1], end)
            else:
                segs.append([start, end])
    return [tuple(s) for s in segs]


def gravy_class(gravy: float) -> str:
    if gravy > 0.3:
        return "hydrophobic (membrane-like)"
    if gravy < -0.3:
        return "hydrophilic (soluble-like)"
    return "intermediate"


def analyze(seq: str, name: str = "sequence") -> dict:
    """Full analysis of one sequence. Raises SequenceError on invalid input."""
    s = clean_sequence(seq)
    n = len(s)
    counts = {aa: s.count(aa) for aa in AMINO_ACIDS}
    other = n - sum(counts.values())
    standard = [c for c in s if c in KYTE_DOOLITTLE]
    gravy = sum(KYTE_DOOLITTLE[c] for c in standard) / len(standard) if standard else 0.0
    profile = hydropathy_profile(s)
    tms = tm_segments(profile)
    return {
        "name": name,
        "length": n,
        "composition": {aa: {"count": counts[aa], "percent": round(100 * counts[aa] / n, 3)} for aa in AMINO_ACIDS},
        "other_residues": other,
        "molecular_weight": round(molecular_weight(s), 2),
        "gravy": round(gravy, 4),
        "gravy_class": gravy_class(gravy),
        "isoelectric_point": isoelectric_point(s),
        "charge_ph7": round(net_charge(s, 7.0), 3),
        "aromaticity": round(sum(counts[a] for a in "FWY") / n, 4),
        "hydropathy_window": TM_WINDOW,
        "hydropathy_profile": [round(v, 3) for v in profile],
        "tm_segments": tms,
        "n_tm_segments": len(tms),
        "sequence": s,
    }


def plot_hydrophobicity(result: dict, path: str) -> str:
    """Save a hydropathy plot (PNG) with predicted TM regions shaded."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    prof = result["hydropathy_profile"]
    half = result["hydropathy_window"] // 2
    fig, ax = plt.subplots(figsize=(9, 3.2))
    ax.plot(range(half + 1, half + 1 + len(prof)), prof, lw=1.2)
    ax.axhline(TM_THRESHOLD, ls="--", c="grey", lw=0.8)
    for a, b in result["tm_segments"]:
        ax.axvspan(a, b, color="orange", alpha=0.25)
    ax.set_xlabel("Residue")
    ax.set_ylabel(f"KD hydropathy (w={result['hydropathy_window']})")
    ax.set_title(result["name"][:80])
    fig.tight_layout()
    fig.savefig(path, dpi=120)
    plt.close(fig)
    return path
