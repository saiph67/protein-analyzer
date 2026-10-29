"""Bench tools: DNA translation and ORF finding, protease digests with peptide masses, A280 concentration."""
from __future__ import annotations

import re

from .analyzer import clean_sequence, molecular_weight

MAX_DNA = 200_000

# ------------------------------------------------------------------ genetic codes
_BASES = "TCAG"
_STANDARD = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
_VERT_MITO = "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG"
CODES = {
    1: {"name": "Standard", "table": _STANDARD, "starts": {"ATG"}},
    2: {"name": "Vertebrate mitochondrial", "table": _VERT_MITO, "starts": {"ATG", "ATA", "ATT", "ATC", "GTG"}},
    11: {"name": "Bacterial, archaeal & plant plastid", "table": _STANDARD,
         "starts": {"ATG", "GTG", "TTG", "CTG", "ATT", "ATC", "ATA"}},
}


def _codon_map(table: str) -> dict[str, str]:
    codons = [a + b + c for a in _BASES for b in _BASES for c in _BASES]
    return dict(zip(codons, table))


_MAPS = {k: _codon_map(v["table"]) for k, v in CODES.items()}
_COMP = str.maketrans("ACGTRYKMBVDHN", "TGCAYRMKVBHDN")


def clean_dna(seq: str) -> str:
    lines = [l for l in (seq or "").splitlines() if not l.startswith(">")]
    s = re.sub(r"[\s\d]+", "", "".join(lines)).upper().replace("U", "T")
    if not s:
        raise ValueError("Paste a DNA or RNA sequence.")
    bad = sorted(set(s) - set("ACGTRYKMSWBDHVN"))
    if bad:
        raise ValueError(f"Not a DNA/RNA letter: {', '.join(bad)}. (Did you paste a protein? Use Analyze for that.)")
    if len(s) > MAX_DNA:
        raise ValueError(f"Sequence is {len(s):,} nt; the limit is {MAX_DNA:,}.")
    return s


def reverse_complement(s: str) -> str:
    return s.translate(_COMP)[::-1]


def translate(dna: str, frame: int = 0, code: int = 1) -> str:
    m = _MAPS[code]
    return "".join(m.get(dna[i:i + 3], "X") for i in range(frame, len(dna) - 2, 3))


def six_frames(dna: str, code: int = 1) -> list[dict]:
    rc = reverse_complement(dna)
    return ([{"strand": "+", "frame": f + 1, "protein": translate(dna, f, code)} for f in range(3)]
            + [{"strand": "-", "frame": f + 1, "protein": translate(rc, f, code)} for f in range(3)])


def find_orfs(dna: str, min_aa: int = 50, code: int = 1, alt_starts: bool = False) -> list[dict]:
    """Open reading frames (start codon ... stop) on both strands, longest first.

    Positions are 1-based on the forward strand; `complete` is False when no stop codon was reached.
    """
    starts = CODES[code]["starts"] if alt_starts else {"ATG"}
    cmap, n, out = _MAPS[code], len(dna), []
    for strand, s in (("+", dna), ("-", reverse_complement(dna))):
        for f in range(3):
            i = f
            while i + 3 <= n:
                if s[i:i + 3] in starts:
                    j, prot = i, []
                    while j + 3 <= n:
                        aa = cmap.get(s[j:j + 3], "X")
                        if aa == "*":
                            break
                        prot.append(aa)
                        j += 3
                    complete = j + 3 <= n
                    if len(prot) >= min_aa:
                        if prot:
                            prot[0] = "M"  # alternative starts are still read as methionine
                        end = j + 3 if complete else j
                        a, b = (i + 1, end) if strand == "+" else (n - end + 1, n - i)
                        out.append({"strand": strand, "frame": f + 1, "start": a, "end": b, "length_nt": end - i,
                                    "length_aa": len(prot), "complete": complete, "start_codon": s[i:i + 3],
                                    "protein": "".join(prot)})
                    i = j + 3  # continue after this ORF (nested starts belong to the same ORF)
                else:
                    i += 3
    return sorted(out, key=lambda o: -o["length_aa"])


def dna_stats(dna: str) -> dict:
    acgt = sum(dna.count(c) for c in "ACGT")
    return {"length": len(dna), "gc_percent": round(100 * (dna.count("G") + dna.count("C")) / max(acgt, 1), 2),
            "ambiguous": len(dna) - acgt}


# ------------------------------------------------------------------ protease digest
MONO = {"G": 57.02146, "A": 71.03711, "S": 87.03203, "P": 97.05276, "V": 99.06841, "T": 101.04768,
        "C": 103.00919, "L": 113.08406, "I": 113.08406, "N": 114.04293, "D": 115.02694, "Q": 128.05858,
        "K": 128.09496, "E": 129.04259, "M": 131.04049, "H": 137.05891, "F": 147.06841, "R": 156.10111,
        "Y": 163.06333, "W": 186.07931, "U": 150.95364, "O": 237.14773}
WATER_MONO, PROTON = 18.010565, 1.007276
CARBAMIDOMETHYL = 57.021464  # monoisotopic
CAM_AVG = 57.0513  # average

# name: (description, cleave-after residues, not-before residues, cleave-before residues)
ENZYMES = {
    "trypsin": ("Trypsin: after K or R, not before P", "KR", "P", ""),
    "trypsin_p": ("Trypsin/P: after K or R, also before P", "KR", "", ""),
    "lysc": ("Lys-C: after K", "K", "", ""),
    "argc": ("Arg-C: after R, not before P", "R", "P", ""),
    "gluc": ("Glu-C (bicarbonate buffer): after E", "E", "", ""),
    "aspn": ("Asp-N: before D", "", "", "D"),
    "chymotrypsin": ("Chymotrypsin (high specificity): after F, Y or W, not before P", "FYW", "P", ""),
    "pepsin": ("Pepsin (pH 1.3): after F or L", "FL", "", ""),
    "cnbr": ("Cyanogen bromide (chemical): after M", "M", "", ""),
}


def cleavage_sites(seq: str, enzyme: str) -> list[int]:
    """0-based cut positions (a cut at k splits seq[:k] | seq[k:])."""
    _, after, not_before, before = ENZYMES[enzyme]
    cuts = []
    for k in range(1, len(seq)):
        if (seq[k - 1] in after and seq[k] not in not_before) or (before and seq[k] in before):
            cuts.append(k)
    return cuts


def peptide_mass(pep: str, cam: bool = False) -> float:
    m = sum(MONO.get(c, 0.0) for c in pep) + WATER_MONO
    if cam:
        m += CARBAMIDOMETHYL * pep.count("C")
    return m


def digest(seq: str, enzyme: str = "trypsin", missed: int = 0, min_len: int = 6, max_len: int = 50,
           min_mass: float = 500.0, max_mass: float = 5000.0, cam: bool = True) -> dict:
    if enzyme not in ENZYMES:
        raise ValueError(f"Unknown enzyme '{enzyme}'.")
    s = clean_sequence(seq)
    missed = max(0, min(int(missed), 3))
    bounds = [0] + cleavage_sites(s, enzyme) + [len(s)]
    peps, covered = [], [False] * len(s)
    for i in range(len(bounds) - 1):
        for m in range(missed + 1):
            if i + m + 1 >= len(bounds):
                break
            a, b = bounds[i], bounds[i + m + 1]
            p = s[a:b]
            mono = peptide_mass(p, cam)
            ok = min_len <= len(p) <= max_len and min_mass <= mono <= max_mass
            if ok:
                for k in range(a, b):
                    covered[k] = True
            peps.append({"start": a + 1, "end": b, "sequence": p, "missed": m, "length": len(p),
                         "mono": round(mono, 4), "average": round(molecular_weight(p) + (CAM_AVG * p.count("C") if cam else 0), 2),
                         "mz1": round(mono + PROTON, 4), "mz2": round((mono + 2 * PROTON) / 2, 4),
                         "mz3": round((mono + 3 * PROTON) / 3, 4), "in_range": ok})
    return {"enzyme": enzyme, "description": ENZYMES[enzyme][0], "sites": len(bounds) - 2, "peptides": peps,
            "n_in_range": sum(p["in_range"] for p in peps), "coverage": round(100 * sum(covered) / len(s), 1),
            "carbamidomethyl": cam, "length": len(s)}


# ------------------------------------------------------------------ A280
def extinction(seq: str) -> dict:
    """Molar extinction coefficient at 280 nm in water (Pace et al. 1995, as used by ExPASy ProtParam)."""
    s = clean_sequence(seq)
    w, y, c = s.count("W"), s.count("Y"), s.count("C")
    mw = molecular_weight(s)
    reduced = 5500 * w + 1490 * y
    oxidized = reduced + 125 * (c // 2)
    return {"W": w, "Y": y, "C": c, "molecular_weight": round(mw, 2),
            "epsilon_reduced": reduced, "epsilon_cystines": oxidized,
            "abs_0_1pct_reduced": round(reduced / mw, 3) if reduced else 0.0,
            "abs_0_1pct_cystines": round(oxidized / mw, 3) if oxidized else 0.0,
            "warning": None if (w or y) else "No W or Y: A280 can't measure this protein. Use a dye assay (Bradford/BCA) instead."}


def concentration(a280: float, epsilon: float, mw: float, path_cm: float = 1.0, dilution: float = 1.0) -> dict:
    if epsilon <= 0:
        raise ValueError("Extinction coefficient must be positive.")
    if path_cm <= 0 or dilution <= 0:
        raise ValueError("Path length and dilution factor must be positive.")
    molar = a280 / (epsilon * path_cm) * dilution
    return {"molar": molar, "micromolar": round(molar * 1e6, 3), "mg_per_ml": round(molar * mw, 4),
            "warning": "A280 above ~1.5 is outside the linear range of most spectrophotometers: dilute and re-measure."
            if a280 > 1.5 else None}
