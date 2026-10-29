"""Fixed-length numeric features for proteins (434 per sequence).

20 amino-acid composition + 400 dipeptide composition + 8 residue-group fractions + 6 physico-chemical.
"""
from __future__ import annotations

import math

import numpy as np

from .analyzer import AMINO_ACIDS, KYTE_DOOLITTLE, clean_sequence, isoelectric_point, net_charge, hydropathy_profile, tm_segments

_IDX = {c: i for i, c in enumerate(AMINO_ACIDS)}
GROUPS = {
    "hydrophobic": "AVILMFWC",
    "polar": "STNQYHC",
    "positive": "KRH",
    "negative": "DE",
    "aromatic": "FWYH",
    "aliphatic": "AVIL",
    "tiny": "AGS",
    "pro_gly": "PG",
}
PHYSCHEM = ["log10_length", "gravy", "pI", "charge_per_residue", "aromaticity", "n_tm_segments"]
FEATURE_NAMES = (
    [f"aac_{a}" for a in AMINO_ACIDS]
    + [f"dpc_{a}{b}" for a in AMINO_ACIDS for b in AMINO_ACIDS]
    + [f"group_{g}" for g in GROUPS]
    + PHYSCHEM
)
N_FEATURES = len(FEATURE_NAMES)  # 434


def featurize(seq: str) -> np.ndarray:
    s = clean_sequence(seq)
    idx = np.array([_IDX[c] for c in s if c in _IDX], dtype=np.int64)
    n = max(len(idx), 1)
    aac = np.bincount(idx, minlength=20) / n
    if len(idx) > 1:
        dp = np.bincount(idx[:-1] * 20 + idx[1:], minlength=400) / (len(idx) - 1)
    else:
        dp = np.zeros(400)
    groups = [sum(aac[_IDX[c]] for c in members) for members in GROUPS.values()]
    gravy = float(sum(KYTE_DOOLITTLE[AMINO_ACIDS[i]] * aac[i] for i in range(20)))
    phys = [
        math.log10(len(s)),
        gravy,
        isoelectric_point(s),
        net_charge(s, 7.0) / len(s),
        float(aac[_IDX["F"]] + aac[_IDX["W"]] + aac[_IDX["Y"]]),
        len(tm_segments(hydropathy_profile(s))),
    ]
    return np.concatenate([aac, dp, groups, phys]).astype(np.float64)


def featurize_many(seqs: list[str]) -> np.ndarray:
    if not seqs:
        return np.zeros((0, N_FEATURES))
    return np.vstack([featurize(s) for s in seqs])
