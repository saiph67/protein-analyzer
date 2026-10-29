"""FAKE data for offline smoke tests only. Never report accuracy measured on this as a result."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

AA = "ACDEFGHIKLMNPQRSTVWY"
BIASES = {"alpha": "KRE", "beta": "LIVF", "gamma": "DSTN"}


def random_seq(rng, length, bias):
    p = np.ones(20)
    for c in bias:
        p[AA.index(c)] = 4
    p /= p.sum()
    return "".join(rng.choice(list(AA), size=length, p=p))


def make_dataset(path, per_class=40, holdout=10, seed=0) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    rows, n = [], 0
    for label, bias in BIASES.items():
        for split, k, date in (("train", per_class, "2010-01-01"), ("holdout", holdout, "2022-01-01")):
            for _ in range(k):
                n += 1
                rows.append({"accession": f"X{n:05d}", "label": label, "split": split, "date_created": date,
                             "length": 0, "name": f"synthetic {label}", "organism": "Synthetica fakeii",
                             "sequence": random_seq(rng, int(rng.integers(80, 300)), bias)})
    df = pd.DataFrame(rows)
    df["length"] = df.sequence.str.len()
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False)
    return df


if __name__ == "__main__":
    out = sys.argv[1] if len(sys.argv) > 1 else "synthetic/dataset.csv"
    make_dataset(out)
    print("wrote", out, "(synthetic - for smoke tests only)")
