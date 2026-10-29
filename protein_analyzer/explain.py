"""Explanations for a single prediction.

1. Region scan (both models): hide one window of the sequence at a time (replace it with X) and measure
   how much the probability of the predicted class drops. Big drops = regions the model relied on.
2. Feature contributions (hand-crafted random forest only): exact decomposition of the forest's probability
   into a starting value plus one contribution per feature, summed along each tree's decision path
   (the "tree interpreter" method). Translated into plain words.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

from .analyzer import AMINO_ACIDS, clean_sequence
from .features import FEATURE_NAMES, GROUPS

MAX_WINDOWS = 40
_AA_NAMES = {"A": "alanine", "R": "arginine", "N": "asparagine", "D": "aspartate", "C": "cysteine",
             "Q": "glutamine", "E": "glutamate", "G": "glycine", "H": "histidine", "I": "isoleucine",
             "L": "leucine", "K": "lysine", "M": "methionine", "F": "phenylalanine", "P": "proline",
             "S": "serine", "T": "threonine", "W": "tryptophan", "Y": "tyrosine", "V": "valine"}
_GROUP_WORDS = {"hydrophobic": "water-repelling residues", "polar": "polar residues",
                "positive": "positively charged residues (K, R, H)", "negative": "negatively charged residues (D, E)",
                "aromatic": "aromatic residues (F, W, Y, H)", "aliphatic": "aliphatic residues (A, V, I, L)",
                "tiny": "tiny residues (A, G, S)", "pro_gly": "proline + glycine"}
_PHYS_WORDS = {"log10_length": "protein length", "gravy": "overall hydropathy (GRAVY)", "pI": "isoelectric point",
               "charge_per_residue": "net charge per residue", "aromaticity": "aromaticity",
               "n_tm_segments": "number of predicted membrane segments"}


def feature_words(name: str) -> str:
    kind, _, rest = name.partition("_")
    if kind == "aac":
        return f"{_AA_NAMES.get(rest, rest)} content"
    if kind == "dpc":
        return f"'{rest}' residue pairs"
    if kind == "group":
        return f"share of {_GROUP_WORDS.get(rest, rest)}"
    return _PHYS_WORDS.get(name, name)


def region_scan(seq: str, predict_proba, target: int, window: int | None = None) -> dict:
    """predict_proba(list_of_seqs) -> (n, n_classes). Returns per-window probability drops."""
    s = clean_sequence(seq)
    L = len(s)
    w = window or max(10, -(-L // MAX_WINDOWS))
    starts = list(range(0, L, w))
    masked = [s[:a] + "X" * len(s[a:a + w]) + s[a + w:] for a in starts]
    probs = predict_proba([s] + masked)
    base = float(probs[0, target])
    drops = [round(base - float(p[target]), 4) for p in probs[1:]]
    return {"window": w, "base_probability": round(base, 4),
            "regions": [{"start": a + 1, "end": min(a + w, L), "drop": d} for a, d in zip(starts, drops)]}


def forest_contributions(forest, x: np.ndarray, target: int, train_mean: np.ndarray, top: int = 8) -> dict:
    """Tree-interpreter decomposition of P(target) for one sample: bias + sum(contributions)."""
    contrib = np.zeros(x.shape[0])
    bias = 0.0
    for est in forest.estimators_:
        t = est.tree_
        vals = t.value[:, 0, :]
        vals = vals / np.maximum(vals.sum(axis=1, keepdims=True), 1e-12)  # counts or fractions -> fractions
        node = 0
        bias += vals[0, target]
        while t.children_left[node] != -1:
            f = t.feature[node]
            child = t.children_left[node] if x[f] <= t.threshold[node] else t.children_right[node]
            contrib[f] += vals[child, target] - vals[node, target]
            node = child
    n = len(forest.estimators_)
    contrib /= n
    bias /= n
    order = np.argsort(-np.abs(contrib))[:top]
    rows = []
    for i in order:
        if abs(contrib[i]) < 1e-4:
            continue
        higher = x[i] > train_mean[i]
        rows.append({"feature": FEATURE_NAMES[i], "words": feature_words(FEATURE_NAMES[i]),
                     "value": round(float(x[i]), 4), "average": round(float(train_mean[i]), 4),
                     "direction": "higher" if higher else "lower",
                     "contribution": round(float(contrib[i]), 4),
                     "sentence": f"{'Higher' if higher else 'Lower'} {feature_words(FEATURE_NAMES[i])} than average "
                                 f"{'supports' if contrib[i] > 0 else 'argues against'} this class"})
    return {"bias": round(float(bias), 4), "features": rows}


def explain(seq: str, model_path: str | Path, name: str | None = None) -> dict:
    from . import ml

    bundle = ml.load_model(model_path)
    s = clean_sequence(seq)
    backend = bundle.get("backend", "handcrafted")
    model = bundle["model"]
    classes = list(model.classes_)
    cache_dir = Path(model_path).parent

    def proba(seqs, cache=True):
        if backend == "esm2":
            from .embeddings import embed
            X = embed(seqs, cache_dir=cache_dir if cache else None)
        else:
            X = ml._features(seqs, backend, cache_dir)
        return model.predict_proba(X)

    p0 = proba([s])[0]
    target = int(np.argmax(p0))
    # ESM-2 re-embeds every masked copy, so it gets fewer (wider) windows to stay responsive.
    window = max(10, -(-len(s) // (20 if backend == "esm2" else MAX_WINDOWS)))
    scan = region_scan(s, lambda ss: np.vstack([p0[None, :], proba(ss[1:], cache=False)]), target, window)
    out = {"name": name or "sequence", "length": len(s), "backend": backend, "prediction": classes[target],
           "probability": round(float(p0[target]), 4), "region_scan": scan, "features": None,
           "method_notes": [
               f"Region scan: each {scan['window']}-residue window was hidden (replaced by X) in turn, and we "
               f"measured how much the probability of '{classes[target]}' fell. Big drops show regions the model relied on."]}
    if backend == "handcrafted":
        x = ml._features([s], "handcrafted", cache_dir)[0]
        mean = bundle.get("train_mean")
        if mean is None:
            mean = ml._features(bundle["train_sequences"][:400], "handcrafted", cache_dir).mean(axis=0)
        out["features"] = forest_contributions(model, x, target, np.asarray(mean))
        out["method_notes"].append(
            "Feature contributions: the forest's probability is split exactly into a starting value (the "
            "average over training proteins) plus one push per feature, added up along every tree's decisions.")
    else:
        out["method_notes"].append(
            "Feature contributions aren't shown for ESM-2: its 480 embedding numbers have no human-readable "
            "meaning, so the region scan is the explanation.")
    return out
