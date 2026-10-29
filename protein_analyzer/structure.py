"""Structure-derived tracks and superposition for AlphaFold models.

Secondary structure: a C-alpha-only assignment in the spirit of P-SEA (Labesse et al. 1997). Helix and strand
are recognised from characteristic C-alpha distances (i to i+2, i+3, i+4). AlphaFold models contain no
HELIX/SHEET records, and full DSSP needs backbone hydrogen-bond geometry, so this is a good, fast approximation.

Disorder: AlphaFold pLDDT < 50 marks regions that are usually intrinsically disordered.

Superposition: residues paired by sequence alignment, then Kabsch least-squares fitting, iteratively refined
onto the core (pairs within 4 Å), reporting RMSD and an approximate TM-score.
"""
from __future__ import annotations

import numpy as np

from .compare import align

DISORDER_PLDDT = 50.0
_HELIX = {2: (5.5, 0.5), 3: (5.3, 0.5), 4: (6.4, 0.6)}
_STRAND = {2: (6.7, 0.6), 3: (9.9, 0.9), 4: (12.4, 1.1)}


def parse_ca(pdb: str) -> dict:
    """C-alpha coordinates and B-factor (= pLDDT in AlphaFold models), keyed by residue number."""
    res, xyz, b, aa = [], [], [], []
    three = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E", "GLY": "G",
             "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P", "SER": "S",
             "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V", "SEC": "U", "PYL": "O"}
    seen = set()
    for line in pdb.splitlines():
        if line.startswith("ATOM") and line[12:16].strip() == "CA":
            r = int(line[22:26])
            if r in seen:
                continue
            seen.add(r)
            res.append(r)
            xyz.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
            b.append(float(line[60:66]) if line[60:66].strip() else 0.0)
            aa.append(three.get(line[17:20], "X"))
    if not res:
        raise ValueError("No C-alpha atoms found in the structure.")
    return {"resnum": np.array(res), "xyz": np.array(xyz), "plddt": np.array(b), "sequence": "".join(aa)}


def secondary_structure(xyz: np.ndarray, plddt: np.ndarray | None = None) -> str:
    """'H' helix, 'E' strand, 'C' coil for each C-alpha (consecutive residues assumed).

    A strand also needs a partner: a C-alpha 4.2-5.6 Å away that is >= 3 residues apart in sequence (a
    neighbouring strand in a sheet). Low-confidence residues (pLDDT < 50) are left as coil, because AlphaFold
    draws disordered regions as extended ribbons that would otherwise look like strands.
    """
    n = len(xyz)
    d = {k: np.full(n, np.nan) for k in (2, 3, 4)}
    for k in (2, 3, 4):
        if n > k:
            d[k][:n - k] = np.linalg.norm(xyz[k:] - xyz[:-k], axis=1)

    def match(crit, i):
        return all(not np.isnan(d[k][i]) and abs(d[k][i] - m) <= tol for k, (m, tol) in crit.items())

    ss = ["C"] * n
    helix_starts = [i for i in range(n) if match(_HELIX, i)]
    for i in helix_starts:  # residue i..i+4 form one helical turn
        for j in range(i, min(i + 5, n)):
            ss[j] = "H"
    has_partner = np.zeros(n, bool)
    for start in range(0, n, 1000):  # chunked distance matrix keeps memory small for long proteins
        block = np.linalg.norm(xyz[start:start + 1000, None, :] - xyz[None, :, :], axis=2)
        idx = np.arange(start, min(start + 1000, n))
        far = np.abs(idx[:, None] - np.arange(n)[None, :]) >= 3
        has_partner[idx] = ((block >= 4.2) & (block <= 5.6) & far).any(axis=1)
    strand_seed = [i for i in range(n) if match(_STRAND, i)]
    for i in strand_seed:
        for j in range(i, min(i + 3, n)):
            if ss[j] != "H" and has_partner[j]:
                ss[j] = "E"
    if plddt is not None:
        ss = [c if p >= DISORDER_PLDDT else "C" for c, p in zip(ss, plddt)]
    s = "".join(ss)
    # clean up: helices shorter than 5 and strands shorter than 3 become coil
    import re
    s = re.sub(r"(?<!H)H{1,4}(?!H)", lambda m: "C" * len(m.group()), s)
    s = re.sub(r"(?<!E)E{1,2}(?!E)", lambda m: "C" * len(m.group()), s)
    return s


def segments(track: str, char: str) -> list[list[int]]:
    out, start = [], None
    for i, c in enumerate(track + "\0"):
        if c == char and start is None:
            start = i
        elif c != char and start is not None:
            out.append([start + 1, i]); start = None
    return out


SS_METHOD_DSSP = ("DSSP (Kabsch & Sander hydrogen-bond rules; PyDSSP algorithm, >97% agreement with original DSSP), "
                  "3 states; disorder = AlphaFold pLDDT < 50")
SS_METHOD_CA = "C-alpha geometry (P-SEA-like; ~78% per-residue agreement with DSSP); disorder = pLDDT < 50"


def structure_features(pdb: str) -> dict:
    from .dssp import assign_pdb

    ca = parse_ca(pdb)
    ss, method = None, SS_METHOD_CA
    try:
        dssp_ss, resnums = assign_pdb(pdb)
        if len(dssp_ss) == len(ca["resnum"]) and list(resnums) == [int(r) for r in ca["resnum"]]:
            ss, method = dssp_ss, SS_METHOD_DSSP
    except (ValueError, IndexError):
        pass
    if ss is None:  # C-alpha-only or incomplete backbone: fall back to the geometric approximation
        ss = secondary_structure(ca["xyz"], ca["plddt"])
    dis = "".join("D" if p < DISORDER_PLDDT else "-" for p in ca["plddt"])
    n = len(ss)
    first = int(ca["resnum"][0])
    shift = lambda segs: [[a + first - 1, b + first - 1] for a, b in segs]
    return {"first_residue": first, "secondary_structure": ss, "plddt": [round(float(x), 1) for x in ca["plddt"]],
            "helix_percent": round(100 * ss.count("H") / n, 1), "strand_percent": round(100 * ss.count("E") / n, 1),
            "coil_percent": round(100 * ss.count("C") / n, 1),
            "disordered_percent": round(100 * dis.count("D") / n, 1),
            "helices": shift(segments(ss, "H")), "strands": shift(segments(ss, "E")),
            "disordered": shift(segments(dis, "D")), "method": method}


# ------------------------------------------------------------------ superposition
def kabsch(P: np.ndarray, Q: np.ndarray):
    """Rotation R and translation t minimising |(P @ R.T + t) - Q|."""
    pc, qc = P.mean(axis=0), Q.mean(axis=0)
    H = (P - pc).T @ (Q - qc)
    U, _, Vt = np.linalg.svd(H)
    dsign = np.sign(np.linalg.det(Vt.T @ U.T))
    D = np.diag([1, 1, dsign])
    R = Vt.T @ D @ U.T
    return R, qc - pc @ R.T


def tm_score(dist: np.ndarray, L: int) -> float:
    d0 = max(0.5, 1.24 * (max(L, 19) - 15) ** (1 / 3) - 1.8)
    return float(np.sum(1 / (1 + (dist / d0) ** 2)) / L)


def superpose(pdb_a: str, pdb_b: str, seq_a: str | None = None, seq_b: str | None = None,
              min_plddt: float = 70.0, core_cutoff: float = 4.0) -> dict:
    """Superpose B onto A. Returns metrics, per-residue distances and B's transformed PDB text."""
    A, B = parse_ca(pdb_a), parse_ca(pdb_b)
    sa, sb = seq_a or A["sequence"], seq_b or B["sequence"]
    aln = align(sa, sb, mode="global")
    ia = {int(r): k for k, r in enumerate(A["resnum"])}
    ib = {int(r): k for k, r in enumerate(B["resnum"])}
    pairs, pa, pb = [], 0, 0
    for x, y in zip(aln["aligned_a"], aln["aligned_b"]):
        if x != "-":
            pa += 1
        if y != "-":
            pb += 1
        if x != "-" and y != "-" and pa in ia and pb in ib:
            pairs.append((pa, pb))
    if len(pairs) < 10:
        raise ValueError("Too few aligned residues with coordinates to superpose these structures.")
    ka = np.array([ia[p] for p, _ in pairs]); kb = np.array([ib[q] for _, q in pairs])
    PA, PB = A["xyz"][ka], B["xyz"][kb]
    conf = (A["plddt"][ka] >= min_plddt) & (B["plddt"][kb] >= min_plddt)
    use = conf if conf.sum() >= 10 else np.ones(len(pairs), bool)
    for _ in range(8):
        R, t = kabsch(PB[use], PA[use])
        dist = np.linalg.norm(PB @ R.T + t - PA, axis=1)
        new = (dist < core_cutoff) & (conf if conf.sum() >= 10 else True)
        if new.sum() < max(10, 0.2 * len(pairs)) or np.array_equal(new, use):
            break
        use = new
    R, t = kabsch(PB[use], PA[use])
    dist = np.linalg.norm(PB @ R.T + t - PA, axis=1)
    core = dist < core_cutoff
    return {
        "aligned_pairs": len(pairs), "core_pairs": int(core.sum()),
        "rmsd_core": round(float(np.sqrt(np.mean(dist[core] ** 2))), 2) if core.sum() >= 3 else None,
        "rmsd_all": round(float(np.sqrt(np.mean(dist ** 2))), 2),
        "tm_score_a": round(tm_score(dist, len(A["xyz"])), 3), "tm_score_b": round(tm_score(dist, len(B["xyz"])), 3),
        "sequence_identity": aln["identity"],
        "per_residue": [{"a": int(A["resnum"][i]), "b": int(B["resnum"][j]), "distance": round(float(d), 2)}
                        for i, j, d in zip(ka, kb, dist)],
        "pdb_b_transformed": _transform_pdb(pdb_b, R, t),
        "len_a": len(A["xyz"]), "len_b": len(B["xyz"]),
    }


def _transform_pdb(pdb: str, R: np.ndarray, t: np.ndarray) -> str:
    out = []
    for line in pdb.splitlines():
        if line.startswith(("ATOM", "HETATM")) and len(line) >= 54:
            x = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])]) @ R.T + t
            line = f"{line[:30]}{x[0]:8.3f}{x[1]:8.3f}{x[2]:8.3f}{line[54:]}"
        out.append(line)
    return "\n".join(out) + "\n"
