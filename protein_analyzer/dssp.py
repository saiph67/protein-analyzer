"""DSSP-style 3-state secondary structure (helix H, strand E, loop C) from backbone atoms.

A memory-efficient numpy port of PyDSSP's algorithm (https://github.com/ShintaroMinami/PyDSSP), which implements
the Kabsch & Sander (1983) hydrogen-bond energy and DSSP's helix/bridge rules; its authors report >97% per-residue
agreement with the original DSSP. Differences from PyDSSP: hydrogen-bond energies are computed in row chunks so
long proteins need little memory, there is no einops/torch dependency, and prolines are not used as H-bond
donors (as in the original DSSP).

PyDSSP is MIT-licensed:
    Copyright (c) 2022 Shintaro Minami. Permission is hereby granted, free of charge, to any person obtaining a
    copy of this software ... THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND. (Full text: MIT
    License, https://opensource.org/license/mit)
"""
from __future__ import annotations

import numpy as np

Q1Q2_F = 0.084 * 332
CUTOFF, MARGIN = -0.5, 1.0
_BACKBONE = {" N  ": 0, " CA ": 1, " C  ": 2, " O  ": 3}


def read_backbone(pdb: str) -> tuple[np.ndarray, str, list[int]]:
    """(L, 4, 3) N/CA/C/O coordinates, one-letter sequence and residue numbers; residues missing an atom are dropped."""
    three = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E", "GLY": "G",
             "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P", "SER": "S",
             "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V"}
    res: dict[tuple, dict] = {}
    order = []
    for line in pdb.splitlines():
        if not line.startswith("ATOM"):
            continue
        k = _BACKBONE.get(line[12:16])
        if k is None:
            continue
        key = (line[21], int(line[22:26]), line[26])
        if key not in res:
            res[key] = {"aa": three.get(line[17:20], "X"), "xyz": [None] * 4}
            order.append(key)
        res[key]["xyz"][k] = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
    keep = [k for k in order if all(a is not None for a in res[k]["xyz"])]
    coord = np.array([res[k]["xyz"] for k in keep], dtype=float).reshape(-1, 4, 3)
    return coord, "".join(res[k]["aa"] for k in keep), [k[1] for k in keep]


def _norm(v):
    return v / np.linalg.norm(v, axis=-1, keepdims=True)


def hbond_map(coord: np.ndarray, donor_mask: np.ndarray | None = None, chunk: int = 256) -> np.ndarray:
    """(L, L) map in [0, 1]; entry [i, j] = residue i's N-H donating to residue j's C=O."""
    L = len(coord)
    N, CA, C, O = coord[:, 0], coord[:, 1], coord[:, 2], coord[:, 3]
    H = N[1:] + 1.01 * _norm(_norm(N[1:] - C[:-1]) + _norm(N[1:] - CA[1:]))  # idealised amide hydrogen
    e = np.zeros((L, L))
    C0, O0 = C[:-1][None], O[:-1][None]
    for s in range(0, L - 1, chunk):
        n, h = N[1 + s:1 + s + chunk][:, None], H[s:s + chunk][:, None]
        d_on = np.linalg.norm(O0 - n, axis=-1); d_ch = np.linalg.norm(C0 - h, axis=-1)
        d_oh = np.linalg.norm(O0 - h, axis=-1); d_cn = np.linalg.norm(C0 - n, axis=-1)
        e[1 + s:1 + s + len(n), :L - 1] = Q1Q2_F * (1 / d_on + 1 / d_ch - 1 / d_oh - 1 / d_cn)
    hb = np.clip(CUTOFF - MARGIN - e, -MARGIN, MARGIN)
    hb = (np.sin(hb / MARGIN * np.pi / 2) + 1) / 2
    idx = np.arange(L)
    hb[(idx[:, None] - idx[None, :] >= 0) & (idx[:, None] - idx[None, :] <= 2)] = 0  # no local (i, i..i-2) pairs
    if donor_mask is not None:
        hb *= np.asarray(donor_mask, float)[:, None]
    return hb


def assign(coord: np.ndarray, donor_mask: np.ndarray | None = None) -> str:
    """3-state DSSP string ('H', 'E', 'C') for backbone coordinates (L, 4, 3)."""
    L = len(coord)
    if L < 6:
        return "C" * L
    X = hbond_map(coord, donor_mask).T > 0  # X[i, j]: C=O of i accepts from N-H of j

    def turns(k):
        t = np.diagonal(X, offset=k)  # i -> i+k
        h = np.zeros(L, bool)
        h[1:1 + len(t) - 1] = t[:-1] & t[1:]  # two consecutive turns start a minimal helix
        return h

    h3, h4, h5 = turns(3), turns(4), turns(5)
    helix4 = h4 | np.roll(h4, 1) | np.roll(h4, 2) | np.roll(h4, 3)
    h3 = h3 & ~np.roll(helix4, -1) & ~helix4
    h5 = h5 & ~np.roll(helix4, -1) & ~helix4
    helix3 = h3 | np.roll(h3, 1) | np.roll(h3, 2)
    helix5 = h5 | np.roll(h5, 1) | np.roll(h5, 2) | np.roll(h5, 3) | np.roll(h5, 4)
    # bridges between residue i+1 and j+1 (parallel p1/p2, antiparallel a1/a2)
    a, b, c = X[:-2, 1:-1], X[1:-1, 2:], X[1:-1, 1:-1]
    bridge = (a & b.T) | (a.T & b) | (c & c.T) | (X[:-2, 2:] & X[:-2, 2:].T)
    ladder = np.zeros(L, bool)
    ladder[1:-1] = bridge.any(axis=1)
    helix = helix3 | helix4 | helix5
    return "".join("H" if hh else "E" if ee else "C" for hh, ee in zip(helix, ladder))


def assign_pdb(pdb: str) -> tuple[str, list[int]]:
    coord, seq, resnums = read_backbone(pdb)
    donor = np.array([c != "P" for c in seq], float)
    return assign(coord, donor), resnums
