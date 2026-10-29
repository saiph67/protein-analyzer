"""Function-class dataset from UniProt, training and honest time-split evaluation.

Train = reviewed UniProt entries first published on/before `cutoff`; hold-out = entries first
published after it. Classes are defined by UniProt keywords; entries carrying more than one of the
chosen keywords are excluded so labels are unambiguous. An optional "other" class holds proteins with
none of the chosen keywords, so the model can say "none of these" instead of forcing a label.

Two feature backends:
  handcrafted - 434 sequence features + random forest (no extra dependencies)
  esm2        - ESM-2 protein language-model embeddings + logistic regression (needs torch/transformers)
"""
from __future__ import annotations

import datetime as dt
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd

from .analyzer import analyze, clean_sequence
from .compare import align, nearest_neighbor_similarity
from .features import FEATURE_NAMES, featurize_many

# Keyword IDs verified against https://rest.uniprot.org/keywords/<id> (2026-10).
DEFAULT_CLASSES = {
    "kinase": "KW-0418",
    "protease": "KW-0645",
    "ion_channel": "KW-0407",
    "ribosomal_protein": "KW-0689",
    "toxin": "KW-0800",
    "dna_binding": "KW-0238",
}
OTHER = "other"
BACKENDS = {"handcrafted": "Random forest on 434 hand-crafted sequence features",
            "esm2": "Logistic regression on ESM-2 language-model embeddings"}
BASE_FILTER = "reviewed:true AND fragment:false AND length:[50 TO 2000]"
# UniProt returns well-studied (mostly human) entries first, so taking the first N results skews every class
# toward particular organisms. We sample each class proportionally across these taxonomic groups instead.
TAXON_STRATA = {
    "bacteria": "taxonomy_id:2",
    "archaea": "taxonomy_id:2157",
    "viruses": "taxonomy_id:10239",
    "fungi": "taxonomy_id:4751",
    "plants": "taxonomy_id:33090",
    "mammals": "taxonomy_id:40674",
    "other_animals": "taxonomy_id:33208 AND NOT taxonomy_id:40674",
    "other_eukaryotes": "taxonomy_id:2759 AND NOT taxonomy_id:4751 AND NOT taxonomy_id:33090 AND NOT taxonomy_id:33208",
}
DATASET_FIELDS = ["accession", "protein_name", "organism_name", "length", "sequence", "date_created"]
SIM_BINS = [0.0, 0.10, 0.20, 0.40, 1.0001]
SIM_LABELS = ["<0.10 (dissimilar)", "0.10-0.20", "0.20-0.40", ">=0.40 (near-homolog)"]
ID_BINS = [0, 30, 50, 70, 100.01]
ID_LABELS = ["<30% identity", "30-50%", "50-70%", ">=70%"]
LOW_SIMILARITY = 0.10
IDENTITY_CANDIDATES = 5     # hold-out proteins are aligned to this many most-similar training proteins
REDUNDANCY_IDENTITY = 40.0  # ... and counted as redundant at >= this % identity to any of them


def _noop(frac, msg):
    pass


def class_query(keyword: str | None, other_keywords: list[str], start: str, end: str) -> str:
    """keyword=None builds the 'other' query: proteins with none of `other_keywords`."""
    q = f"{BASE_FILTER} AND date_created:[{start} TO {end}]"
    if keyword:
        q = f"keyword:{keyword} AND " + q
    for k in other_keywords:
        q += f" AND NOT keyword:{k}"
    return q


def _allocate(counts: dict[str, int], total: int) -> dict[str, int]:
    """Split `total` across strata proportionally to `counts` (largest remainder), capped by availability."""
    avail = sum(counts.values())
    if not avail:
        return {k: 0 for k in counts}
    total = min(total, avail)
    raw = {k: total * c / avail for k, c in counts.items()}
    alloc = {k: min(int(v), counts[k]) for k, v in raw.items()}
    for k in sorted(raw, key=lambda k: raw[k] - int(raw[k]), reverse=True):
        if sum(alloc.values()) >= total:
            break
        if alloc[k] < counts[k]:
            alloc[k] += 1
    return alloc


def _stratified_pool(client, q: str, pool_n: int, progress, label: str) -> tuple[list[dict], dict]:
    counts = {s: client.count(f"{q} AND {t}") for s, t in TAXON_STRATA.items()}
    alloc = _allocate(counts, pool_n)
    pool = []
    for s, n in alloc.items():
        if n:
            progress(None, f"{label}: {s} ({n} of {counts[s]})")
            pool += list(client.iter_search(f"{q} AND {TAXON_STRATA[s]}", limit=n, fields=DATASET_FIELDS))
    return pool, counts


def build_dataset(client, classes: dict[str, str] | None = None, per_class: int = 400, holdout: int = 100,
                  cutoff: str = "2018-12-31", out_dir: str | os.PathLike = "results", pool_factor: int = 3,
                  seed: int = 0, include_other: bool = False, stratify: bool = True, progress=_noop) -> dict:
    """Download a labelled dataset from UniProt and write `dataset.csv` + `dataset_summary.json`."""
    classes = dict(classes or DEFAULT_CLASSES)
    if len(classes) + include_other < 2:
        raise ValueError("Choose at least two classes.")
    if include_other:
        classes[OTHER] = None
    cut = dt.date.fromisoformat(cutoff)
    after = (cut + dt.timedelta(days=1)).isoformat()
    out = Path(out_dir); out.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(seed)
    rows, summary, warnings = [], {}, []
    steps = len(classes) * 2
    step = 0
    for label, kw in classes.items():
        others = [k for k in classes.values() if k and k != kw]
        for split, (start, end), want in (("train", ("1986-01-01", cutoff), per_class),
                                          ("holdout", (after, "*"), holdout)):
            q = class_query(kw, others, start, end)
            frac = step / steps
            progress(frac, f"{label} / {split}: counting")
            available = client.count(q)
            pool_n = min(available, want * pool_factor)
            if stratify:
                pool, strata = _stratified_pool(client, q, pool_n,
                                                lambda _f, msg: progress(frac, msg), f"{label} / {split}")
                summary.setdefault(label, {"keyword": kw})[f"{split}_strata_available"] = strata
            else:
                progress(frac, f"{label} / {split}: downloading {pool_n} of {available} entries")
                pool = list(client.iter_search(q, limit=pool_n, fields=DATASET_FIELDS))
            if len(pool) > want:
                pool = [pool[i] for i in sorted(rng.choice(len(pool), size=want, replace=False))]
            for e in pool:
                rows.append({"accession": e["accession"], "label": label, "split": split,
                             "date_created": e["date_created"], "length": e["length"],
                             "name": e["name"], "organism": e["organism"], "sequence": e["sequence"]})
            summary.setdefault(label, {"keyword": kw})[split] = len(pool)
            summary[label][f"{split}_available"] = available
            if len(pool) < want:
                warnings.append(f"{label}: only {len(pool)} {split} entries available (wanted {want}).")
            step += 1
    df = pd.DataFrame(rows)
    before = len(df)
    df = df.drop_duplicates("accession")
    df = df.drop_duplicates("sequence")  # identical sequences across classes/splits
    train_seqs = set(df.loc[df.split == "train", "sequence"])
    df = df[~((df.split == "holdout") & df.sequence.isin(train_seqs))]
    if len(df) < before:
        warnings.append(f"Removed {before - len(df)} duplicate accessions/sequences.")
    # Sanity check of the time split
    dates = pd.to_datetime(df.date_created)
    assert (dates[df.split == "train"] <= pd.Timestamp(cut)).all()
    assert (dates[df.split == "holdout"] > pd.Timestamp(cut)).all()
    df.to_csv(out / "dataset.csv", index=False)
    counts = df.groupby(["label", "split"]).size().unstack(fill_value=0)
    info = {
        "cutoff": cutoff, "per_class": per_class, "holdout": holdout, "classes": classes,
        "base_filter": BASE_FILTER, "built": dt.datetime.now().isoformat(timespec="seconds"),
        "sampling": ("proportional across taxonomic groups (" + ", ".join(TAXON_STRATA) + ")"
                     if stratify else "first results in UniProt's default order"),
        "human_fraction": {lab: round(float((g.organism == "Homo sapiens").mean()), 3)
                           for lab, g in df.groupby("label")},
        "counts": {lab: {k: int(v) for k, v in row.items()} for lab, row in counts.iterrows()},
        "availability": summary, "warnings": warnings, "n_total": int(len(df)),
    }
    (out / "dataset_summary.json").write_text(json.dumps(info, indent=2))
    progress(1.0, f"Dataset ready: {len(df)} proteins")
    return info


def load_dataset(path) -> pd.DataFrame:
    df = pd.read_csv(path)
    df["sequence"] = df["sequence"].map(clean_sequence)
    return df


def find_patterns(df: pd.DataFrame, out_dir="results") -> dict:
    """Per-class property summaries and most enriched residues; writes a boxplot PNG."""
    out = Path(out_dir); out.mkdir(parents=True, exist_ok=True)
    props = []
    for _, r in df.iterrows():
        a = analyze(r.sequence, r.accession)
        props.append({"label": r.label, "length": a["length"], "gravy": a["gravy"], "pI": a["isoelectric_point"],
                      "charge_per_100": 100 * a["charge_ph7"] / a["length"], "aromaticity": a["aromaticity"],
                      "tm_segments": a["n_tm_segments"],
                      **{f"aa_{k}": v["percent"] for k, v in a["composition"].items()}})
    p = pd.DataFrame(props)
    overall = p.filter(like="aa_").mean()
    patterns = {}
    for lab, g in p.groupby("label"):
        enr = (g.filter(like="aa_").mean() - overall).sort_values(ascending=False)
        patterns[lab] = {
            "n": int(len(g)),
            "median_length": float(g.length.median()),
            "mean_gravy": round(float(g.gravy.mean()), 3),
            "mean_pI": round(float(g.pI.mean()), 2),
            "mean_charge_per_100": round(float(g.charge_per_100.mean()), 2),
            "mean_aromaticity": round(float(g.aromaticity.mean()), 3),
            "frac_with_tm": round(float((g.tm_segments > 0).mean()), 3),
            "enriched": [f"{k[3:]} (+{v:.1f}%)" for k, v in enr.head(3).items()],
            "depleted": [f"{k[3:]} ({v:.1f}%)" for k, v in enr.tail(3).items()],
        }
    _boxplots(p, out / "property_boxplots.png")
    return patterns


# ------------------------------------------------------------------ training

def model_path(results_dir, backend: str) -> Path:
    return Path(results_dir) / f"model_{backend}.joblib"


def available_models(results_dir) -> list[dict]:
    """Trained models in `results_dir`, best hold-out accuracy first."""
    out = []
    for b in BACKENDS:
        mp = Path(results_dir) / f"metrics_{b}.json"
        if model_path(results_dir, b).exists() and mp.exists():
            m = json.loads(mp.read_text())
            out.append({"backend": b, "label": BACKENDS[b], "accuracy": m["holdout"]["accuracy"],
                        "low_similarity_accuracy": (m.get("holdout_low_similarity") or {}).get("accuracy"),
                        "classes": m["classes"], "trained": m["trained"], "needs_esm": b == "esm2"})
    return sorted(out, key=lambda d: -d["accuracy"])


def _features(seqs, backend, cache_dir, progress=_noop, lo=0.0, hi=1.0):
    if backend == "handcrafted":
        return featurize_many(seqs)
    if backend == "esm2":
        from .embeddings import embed
        return embed(seqs, cache_dir=cache_dir,
                     progress=lambda f, msg: progress(lo + (hi - lo) * f, msg))
    raise ValueError(f"Unknown backend {backend!r}; choose one of {', '.join(BACKENDS)}.")


def train_and_evaluate(dataset_csv="results/dataset.csv", out_dir="results", n_estimators: int = 500,
                       align_identity: bool = True, seed: int = 0, backend: str = "handcrafted",
                       progress=_noop) -> dict:
    from sklearn.dummy import DummyClassifier
    from sklearn.ensemble import RandomForestClassifier
    from sklearn.linear_model import LogisticRegression
    from sklearn.metrics import accuracy_score, balanced_accuracy_score, confusion_matrix, f1_score
    from sklearn.model_selection import GridSearchCV
    from sklearn.pipeline import make_pipeline
    from sklearn.preprocessing import StandardScaler
    import joblib

    if backend not in BACKENDS:
        raise ValueError(f"Unknown backend {backend!r}.")
    out = Path(out_dir); out.mkdir(parents=True, exist_ok=True)
    df = load_dataset(dataset_csv)
    tr, ho = df[df.split == "train"].reset_index(drop=True), df[df.split == "holdout"].reset_index(drop=True)
    if tr.label.nunique() < 2 or ho.empty:
        raise ValueError("Dataset needs at least two classes in train and a non-empty hold-out set.")
    classes = sorted(tr.label.unique())
    progress(0.02, f"Hand-crafted features for {len(df)} sequences")
    Htr, Hho = featurize_many(list(tr.sequence)), featurize_many(list(ho.sequence))
    ytr, yho = tr.label.values, ho.label.values

    def scores(pred):
        return {"accuracy": round(accuracy_score(yho, pred), 4),
                "balanced_accuracy": round(balanced_accuracy_score(yho, pred), 4),
                "macro_f1": round(f1_score(yho, pred, average="macro", labels=classes, zero_division=0), 4)}

    rf = RandomForestClassifier(n_estimators=n_estimators, class_weight="balanced_subsample",
                                n_jobs=-1, random_state=seed, min_samples_leaf=1)
    baselines = {}
    if backend == "handcrafted":
        Xtr, Xho = Htr, Hho
        progress(0.25, "Training random forest (434 features)")
        model = rf.fit(Xtr, ytr)
        model_desc = f"RandomForest({n_estimators} trees, class_weight=balanced_subsample)"
    else:
        progress(0.05, "Computing ESM-2 embeddings (first run downloads the model, ~140 MB)")
        X = _features(list(df.sequence), "esm2", out, progress, 0.05, 0.40)
        Xtr, Xho = X[(df.split == "train").values], X[(df.split == "holdout").values]
        progress(0.40, "Training logistic regression on embeddings (C chosen by 5-fold CV on train only)")
        grid = GridSearchCV(make_pipeline(StandardScaler(), LogisticRegression(max_iter=5000, class_weight="balanced")),
                            {"logisticregression__C": [0.01, 0.1, 1.0]}, cv=5, scoring="balanced_accuracy", n_jobs=-1)
        model = grid.fit(Xtr, ytr).best_estimator_
        c = grid.best_params_["logisticregression__C"]
        model_desc = f"LogisticRegression(C={c}) on {Xtr.shape[1]}-d ESM-2 embeddings"
        progress(0.43, "Hand-crafted random forest for comparison")
        rf.fit(Htr, ytr)
        baselines["random forest, hand-crafted features"] = scores(rf.predict(Hho))
    pred = model.predict(Xho)
    proba = model.predict_proba(Xho)

    progress(0.45, "Baselines")
    for name, clf, cols in (
        ("majority class", DummyClassifier(strategy="most_frequent"), slice(None)),
        ("random (stratified)", DummyClassifier(strategy="stratified", random_state=seed), slice(None)),
        ("logistic regression, composition only (20 AAC)",
         make_pipeline(StandardScaler(), LogisticRegression(max_iter=3000, class_weight="balanced")), slice(0, 20)),
    ):
        clf.fit(Htr[:, cols], ytr)
        baselines[name] = scores(clf.predict(Hho[:, cols]))

    progress(0.5, "Nearest-neighbour similarity of hold-out to training set")
    nn = nearest_neighbor_similarity(list(ho.sequence), list(tr.sequence), k=3, top=IDENTITY_CANDIDATES)
    nn_idx = np.array([x[0][0] for x in nn]); nn_sim = np.array([x[0][1] for x in nn])
    baselines["1-NN k-mer homology transfer"] = scores(ytr[nn_idx])

    identity = None
    if align_identity:
        # Max local-alignment identity to the IDENTITY_CANDIDATES most k-mer-similar training proteins.
        identity = np.zeros(len(ho))
        for i, (q, cands) in enumerate(zip(ho.sequence, nn)):
            best = 0.0
            for j, _ in cands:
                r = tr.sequence[j]
                a = align(q, r, mode="local")
                best = max(best, 100 * a["identical"] / min(len(q), len(r)))
            identity[i] = best
            if i % 10 == 0:
                progress(0.5 + 0.33 * i / len(ho),
                         f"Redundancy check: aligning hold-out protein {i}/{len(ho)} to its {len(cands)} closest training proteins")

    correct = pred == yho
    by_sim = _binned(correct, nn_sim, SIM_BINS, SIM_LABELS)
    by_id = _binned(correct, identity, ID_BINS, ID_LABELS) if identity is not None else []
    low = nn_sim < LOW_SIMILARITY
    low_scores = None
    if low.sum():
        low_scores = {"n": int(low.sum()),
                      "accuracy": round(accuracy_score(yho[low], pred[low]), 4),
                      "balanced_accuracy": round(balanced_accuracy_score(yho[low], pred[low]), 4)}
    strict = filtered = None
    if identity is not None and (identity < 30).sum():
        m = identity < 30
        strict = {"n": int(m.sum()), "accuracy": round(accuracy_score(yho[m], pred[m]), 4),
                  "balanced_accuracy": round(balanced_accuracy_score(yho[m], pred[m]), 4)}
    if identity is not None and (identity < REDUNDANCY_IDENTITY).sum():
        m = identity < REDUNDANCY_IDENTITY
        filtered = {"n": int(m.sum()), "removed": int((~m).sum()), "threshold": REDUNDANCY_IDENTITY,
                    "candidates": IDENTITY_CANDIDATES,
                    "accuracy": round(accuracy_score(yho[m], pred[m]), 4),
                    "balanced_accuracy": round(balanced_accuracy_score(yho[m], pred[m]), 4),
                    "macro_f1": round(f1_score(yho[m], pred[m], average="macro", labels=classes, zero_division=0), 4),
                    "baseline_1nn": round(accuracy_score(yho[m], ytr[nn_idx][m]), 4)}

    progress(0.85, "Plots and report")
    cm = confusion_matrix(yho, pred, labels=classes)
    imp = sorted(zip(FEATURE_NAMES, rf.feature_importances_), key=lambda t: -t[1])[:25]
    per_class = {c: {"n": int((yho == c).sum()),
                     "recall": round(float((pred[yho == c] == c).mean()), 4) if (yho == c).any() else None}
                 for c in classes}
    metrics = {
        "backend": backend,
        "model": model_desc,
        "n_features": int(Xtr.shape[1]),
        "classes": classes,
        "n_train": int(len(tr)), "n_holdout": int(len(ho)),
        "train_counts": {k: int(v) for k, v in tr.label.value_counts().items()},
        "holdout_counts": {k: int(v) for k, v in ho.label.value_counts().items()},
        "holdout": scores(pred),
        "holdout_low_similarity": low_scores,
        "holdout_below_30pct_identity": strict,
        "holdout_redundancy_filtered": filtered,
        "per_class": per_class,
        "baselines": baselines,
        "confusion_matrix": cm.tolist(),
        "accuracy_by_similarity": by_sim,
        "accuracy_by_identity": by_id,
        "median_nn_similarity": round(float(np.median(nn_sim)), 3),
        "median_nn_identity": round(float(np.median(identity)), 1) if identity is not None else None,
        "top_features": [{"feature": f, "importance": round(float(v), 5)} for f, v in imp],
        "mean_confidence": round(float(proba.max(axis=1).mean()), 3),
        "trained": dt.datetime.now().isoformat(timespec="seconds"),
    }
    summary_path = Path(dataset_csv).with_name("dataset_summary.json")
    if summary_path.exists():
        metrics["dataset"] = json.loads(summary_path.read_text())

    sfx = f"_{backend}"
    metrics["figures"] = [f"confusion_matrix{sfx}.png", f"accuracy_by_similarity{sfx}.png",
                          f"feature_importance{sfx}.png", f"pca_map{sfx}.png", "property_boxplots.png"]
    _plot_confusion(cm, classes, out / f"confusion_matrix{sfx}.png")
    _plot_importance(imp, out / f"feature_importance{sfx}.png",
                     "Top feature importances (hand-crafted random forest" +
                     (", comparison model)" if backend == "esm2" else ")"))
    _plot_pca(np.vstack([Xtr, Xho]), np.concatenate([ytr, yho]), len(tr), classes, out / f"pca_map{sfx}.png",
              "ESM-2 embeddings" if backend == "esm2" else "hand-crafted features")
    _plot_by_bin(by_sim, by_id, out / f"accuracy_by_similarity{sfx}.png")
    metrics["patterns"] = find_patterns(df, out)
    progress(0.97, "Saving model")

    joblib.dump({"model": model, "backend": backend, "classes": classes, "feature_names": FEATURE_NAMES,
                 "train_mean": Htr.mean(axis=0),
                 "train_sequences": list(tr.sequence), "train_accessions": list(tr.accession),
                 "train_labels": list(ytr), "metrics": metrics}, model_path(out, backend), compress=3)
    text = json.dumps(metrics, indent=2)
    report = render_report(metrics)
    (out / f"metrics{sfx}.json").write_text(text)
    (out / f"report{sfx}.md").write_text(report, encoding="utf-8")
    (out / "report.md").write_text(report, encoding="utf-8")  # latest run
    progress(1.0, f"Done: hold-out accuracy {metrics['holdout']['accuracy']:.1%}")
    return metrics


def _binned(correct, values, edges, labels):
    rows = []
    for lo, hi, lab in zip(edges[:-1], edges[1:], labels):
        m = (values >= lo) & (values < hi)
        rows.append({"bin": lab, "n": int(m.sum()),
                     "accuracy": round(float(correct[m].mean()), 4) if m.any() else None})
    return rows


def render_report(m: dict) -> str:
    from tabulate import tabulate

    pct = lambda x: "n/a" if x is None else f"{100 * x:.1f}%"
    h = m["holdout"]
    lines = [
        "# Protein function classifier - evaluation report", "",
        f"Generated {m['trained']}. All data downloaded from the **live UniProtKB REST API**.", "",
        "## Setup", "",
        f"- Model: {m['model']} ({m['n_features']} input features). " + (
            "Features: amino-acid & dipeptide composition, residue groups, physico-chemistry."
            if m.get("backend", "handcrafted") == "handcrafted" else
            "Features: mean-pooled ESM-2 protein language-model embeddings (pre-trained on UniRef; the "
            "language model never saw our labels). Caveat: its unlabelled pre-training data (~2021) includes "
            "some sequences from our post-cutoff hold-out, so results may be slightly optimistic."),
    ]
    ds = m.get("dataset")
    if ds:
        lines += [f"- Time split: train = reviewed entries first published <= **{ds['cutoff']}**; "
                  f"hold-out = published after.",
                  f"- Filter: `{ds['base_filter']}`; entries with >1 of the class keywords excluded.",
                  "- Classes: " + ", ".join(f"{k} ({v or 'none of the other keywords'})" for k, v in ds["classes"].items())]
        for w in ds.get("warnings", []):
            lines.append(f"- Note: {w}")
    rows = [[c, m["train_counts"].get(c, 0), m["holdout_counts"].get(c, 0), pct(m["per_class"][c]["recall"])]
            for c in m["classes"]]
    lines += ["", tabulate(rows, headers=["class", "train", "hold-out", "hold-out recall"], tablefmt="github"), "",
              "## Hold-out results (time split)", "",
              tabulate([[f"**{BACKENDS[m.get('backend', 'handcrafted')]} (this model)**", pct(h["accuracy"]), pct(h["balanced_accuracy"]), pct(h["macro_f1"])]]
                       + [[k, pct(v["accuracy"]), pct(v["balanced_accuracy"]), pct(v["macro_f1"])]
                          for k, v in m["baselines"].items()],
                       headers=["method", "accuracy", "balanced acc.", "macro-F1"], tablefmt="github"), ""]
    lines += ["## Leakage check: accuracy vs. similarity to the training set", "",
              f"For each hold-out protein we find its {IDENTITY_CANDIDATES} most similar training proteins "
              "(3-mer Jaccard), align it to each (local alignment; identical residues / shorter length) and keep "
              "the highest identity. A newer entry that is a close homolog of an older one is easy; the "
              "low-similarity rows are the honest measure of generalisation to genuinely new proteins.", "",
              tabulate([[b["bin"], b["n"], pct(b["accuracy"])] for b in m["accuracy_by_similarity"]],
                       headers=["3-mer similarity", "n", "accuracy"], tablefmt="github"), ""]
    if m.get("accuracy_by_identity"):
        lines += [tabulate([[b["bin"], b["n"], pct(b["accuracy"])] for b in m["accuracy_by_identity"]],
                           headers=["identity to nearest train protein", "n", "accuracy"], tablefmt="github"), ""]
    lo, st = m.get("holdout_low_similarity"), m.get("holdout_below_30pct_identity")
    lines.append(f"- Median nearest-neighbour similarity: {m['median_nn_similarity']}"
                 + (f"; median identity: {m['median_nn_identity']}%" if m.get("median_nn_identity") is not None else ""))
    if lo:
        lines.append(f"- **Low-similarity hold-out (3-mer < 0.10, n={lo['n']}): accuracy {pct(lo['accuracy'])}, "
                     f"balanced {pct(lo['balanced_accuracy'])}**")
    if st:
        lines.append(f"- **Hold-out with < 30% identity to any training protein (n={st['n']}): accuracy "
                     f"{pct(st['accuracy'])}, balanced {pct(st['balanced_accuracy'])}**")
    rf_ = m.get("holdout_redundancy_filtered")
    if rf_:
        lines += ["", "### Redundancy-filtered hold-out", "",
                  f"Removing the {rf_['removed']} hold-out proteins that are >= {rf_['threshold']:.0f}% identical to any of "
                  f"their {rf_['candidates']} closest training proteins leaves n={rf_['n']}:", "",
                  tabulate([["this model", pct(rf_["accuracy"]), pct(rf_["balanced_accuracy"]), pct(rf_["macro_f1"])],
                            ["1-NN homology transfer", pct(rf_["baseline_1nn"]), "", ""]],
                           headers=["method", "accuracy", "balanced acc.", "macro-F1"], tablefmt="github"), "",
                  "This approximates clustering at 40% identity (as MMseqs2/CD-HIT would) for the train/hold-out "
                  "boundary. Candidates come from a 3-mer prefilter, so a rare distant homolog could be missed."]
    lines += ["", "## Confusion matrix (rows = true, columns = predicted)", "",
              tabulate([[c] + row for c, row in zip(m["classes"], m["confusion_matrix"])],
                       headers=["true \\ pred"] + m["classes"], tablefmt="github"), "",
              "## Top features" + (" (hand-crafted random forest)" if m.get("backend") == "esm2" else ""), "",
              tabulate([[f["feature"], f["importance"]] for f in m["top_features"][:15]],
                       headers=["feature", "importance"], tablefmt="github"), ""]
    if m.get("patterns"):
        lines += ["## Class patterns", "",
                  tabulate([[k, v["n"], v["median_length"], v["mean_gravy"], v["mean_pI"], v["frac_with_tm"],
                             ", ".join(v["enriched"])] for k, v in m["patterns"].items()],
                           headers=["class", "n", "median len", "GRAVY", "pI", "frac. TM", "enriched residues"],
                           tablefmt="github"), ""]
    lines += ["## Figures", ""] + [f"![{f[:-4]}]({f})" for f in m.get("figures", [])] + [""]
    return "\n".join(lines)


# ------------------------------------------------------------------ prediction

_MODEL_CACHE: dict = {}


def load_model(path="results/model_handcrafted.joblib") -> dict:
    import joblib

    p = Path(path)
    if not p.exists():
        raise FileNotFoundError("No trained model yet. Build a dataset and train on the Dataset & ML page first.")
    key = (str(p.resolve()), p.stat().st_mtime)
    if key not in _MODEL_CACHE:
        for k in [k for k in _MODEL_CACHE if k[0] == key[0]]:
            del _MODEL_CACHE[k]  # stale version of the same file
        _MODEL_CACHE[key] = joblib.load(p)
    return _MODEL_CACHE[key]


def predict(seqs: list[str], model_path="results/model_handcrafted.joblib", names: list[str] | None = None) -> list[dict]:
    bundle = load_model(model_path)
    seqs = [clean_sequence(s) for s in seqs]
    X = _features(seqs, bundle.get("backend", "handcrafted"), Path(model_path).parent)
    proba = bundle["model"].predict_proba(X)
    classes = list(bundle["model"].classes_)
    nn = nearest_neighbor_similarity(seqs, bundle["train_sequences"], k=3, top=1)
    out = []
    for i, p in enumerate(proba):
        order = np.argsort(-p)
        j, sim = nn[i][0]
        out.append({
            "name": names[i] if names else f"sequence_{i + 1}",
            "length": len(seqs[i]),
            "sequence": seqs[i],
            "prediction": classes[order[0]],
            "confidence": round(float(p[order[0]]), 4),
            "top2": [{"label": classes[k], "probability": round(float(p[k]), 4)} for k in order[:2]],
            "probabilities": {classes[k]: round(float(p[k]), 4) for k in order},
            "nn_similarity": round(sim, 4),
            "nn_accession": bundle["train_accessions"][j],
            "nn_label": bundle["train_labels"][j],
            "reliable": sim >= LOW_SIMILARITY,
        })
    return out


# ------------------------------------------------------------------ plots

def _plt():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    return plt


def _plot_confusion(cm, classes, path):
    plt = _plt()
    fig, ax = plt.subplots(figsize=(6.5, 5.5))
    norm = cm / np.maximum(cm.sum(axis=1, keepdims=True), 1)
    ax.imshow(norm, cmap="Blues", vmin=0, vmax=1)
    for i in range(len(classes)):
        for j in range(len(classes)):
            ax.text(j, i, cm[i, j], ha="center", va="center", color="white" if norm[i, j] > 0.5 else "black")
    ax.set_xticks(range(len(classes)), classes, rotation=40, ha="right")
    ax.set_yticks(range(len(classes)), classes)
    ax.set_xlabel("Predicted"); ax.set_ylabel("True"); ax.set_title("Hold-out confusion matrix")
    fig.tight_layout(); fig.savefig(path, dpi=120); plt.close(fig)


def _plot_importance(imp, path, title="Top feature importances (random forest)"):
    plt = _plt()
    fig, ax = plt.subplots(figsize=(7, 6))
    names, vals = zip(*imp[::-1])
    ax.barh(names, vals, color="#3b82f6")
    ax.set_title(title, fontsize=10)
    fig.tight_layout(); fig.savefig(path, dpi=120); plt.close(fig)


def _plot_pca(X, y, n_train, classes, path, what="features"):
    from sklearn.decomposition import PCA
    from sklearn.preprocessing import StandardScaler

    plt = _plt()
    Z = PCA(n_components=2, random_state=0).fit_transform(StandardScaler().fit_transform(X))
    fig, ax = plt.subplots(figsize=(7, 5.5))
    cmap = plt.get_cmap("tab10")
    for k, c in enumerate(classes):
        m = y == c
        tr_m, ho_m = m.copy(), m.copy()
        tr_m[n_train:] = False; ho_m[:n_train] = False
        ax.scatter(Z[tr_m, 0], Z[tr_m, 1], s=8, alpha=0.5, color=cmap(k), label=c)
        ax.scatter(Z[ho_m, 0], Z[ho_m, 1], s=18, alpha=0.9, color=cmap(k), marker="x")
    ax.set_title(f"PCA of {what} (dots = train, x = hold-out)", fontsize=10)
    ax.legend(fontsize=8, markerscale=2)
    fig.tight_layout(); fig.savefig(path, dpi=120); plt.close(fig)


def _plot_by_bin(by_sim, by_id, path):
    plt = _plt()
    panels = [("3-mer similarity to nearest training protein", by_sim)]
    if by_id:
        panels.append(("Identity to nearest training protein", by_id))
    fig, axes = plt.subplots(1, len(panels), figsize=(5.5 * len(panels), 4), squeeze=False)
    for ax, (title, rows) in zip(axes[0], panels):
        labs = [f"{r['bin']}\n(n={r['n']})" for r in rows]
        vals = [r["accuracy"] or 0 for r in rows]
        ax.bar(labs, vals, color="#10b981")
        ax.set_ylim(0, 1); ax.set_ylabel("Hold-out accuracy"); ax.set_title(title, fontsize=10)
        ax.tick_params(axis="x", labelsize=8)
    fig.tight_layout(); fig.savefig(path, dpi=120); plt.close(fig)


def _boxplots(p: pd.DataFrame, path):
    plt = _plt()
    cols = ["length", "gravy", "pI", "charge_per_100", "aromaticity", "tm_segments"]
    labels = sorted(p.label.unique())
    fig, axes = plt.subplots(2, 3, figsize=(13, 7))
    for ax, c in zip(axes.ravel(), cols):
        ax.boxplot([p.loc[p.label == l, c] for l in labels], showfliers=False)
        ax.set_xticks(range(1, len(labels) + 1), labels, rotation=35, ha="right", fontsize=8)
        ax.set_title(c)
    fig.tight_layout(); fig.savefig(path, dpi=110); plt.close(fig)
