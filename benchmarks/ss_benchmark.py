"""Benchmark the app's secondary-structure assignment against DSSP on AlphaFold models.

Reference: PyDSSP (Minami, MIT licence), a reimplementation of DSSP's Kabsch-Sander hydrogen-bond algorithm
whose authors report >97% per-residue agreement with the original DSSP in 3 states (helix H, strand E, loop).

Proteins: a stratified random sample from results/dataset.csv (all function classes, many organisms).
Run:  python benchmarks/ss_benchmark.py [--per-class 25]
Writes benchmarks/ss_benchmark.json and benchmarks/ss_benchmark.md
"""
from __future__ import annotations

import argparse
import datetime as dt
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from protein_analyzer import config  # noqa: E402
from protein_analyzer.structure import parse_ca, secondary_structure  # noqa: E402
from protein_analyzer.uniprot import UniProtClient, UniProtError  # noqa: E402

STATES = "HEC"


def dssp_states(pdb: str) -> str:
    from pydssp import read_pdbtext
    from pydssp.pydssp_numpy import assign
    coord = read_pdbtext(pdb)  # (L, 4, 3): N, CA, C, O
    onehot = assign(np.asarray(coord, dtype=float))
    return "".join("CHE"[i] for i in onehot.argmax(-1))


def scores(pred: str, ref: str, mask: np.ndarray) -> dict:
    p = np.array(list(pred)); r = np.array(list(ref))
    p, r = p[mask], r[mask]
    out = {"n": int(len(r)), "q3": round(float((p == r).mean()), 4) if len(r) else None,
           "confusion": {a: {b: int(((r == a) & (p == b)).sum()) for b in STATES} for a in STATES}}
    for s in STATES:
        tp = int(((p == s) & (r == s)).sum()); fp = int(((p == s) & (r != s)).sum()); fn = int(((p != s) & (r == s)).sum())
        prec = tp / (tp + fp) if tp + fp else None; rec = tp / (tp + fn) if tp + fn else None
        f1 = 2 * prec * rec / (prec + rec) if prec and rec else None
        out[s] = {"precision": round(prec, 4) if prec is not None else None, "recall": round(rec, 4) if rec is not None else None,
                  "f1": round(f1, 4) if f1 is not None else None, "support": int((r == s).sum())}
    return out


def run(per_class: int = 25, seed: int = 0, method=None):
    method = method or (lambda pdb, ca: secondary_structure(ca["xyz"], ca["plddt"]))
    df = pd.read_csv(ROOT / "results" / "dataset.csv")
    df = df[df.length <= 1500]
    sample = pd.concat([g.sample(min(per_class, len(g)), random_state=seed) for _, g in df.groupby("label")])
    client = UniProtClient(cache_dir=config.CACHE_DIR)
    preds, refs, plddts, per_protein, skipped = [], [], [], [], []
    for acc, label, org in zip(sample.accession, sample.label, sample.organism):
        try:
            pdb = client.fetch_alphafold(acc)["pdb"]
        except UniProtError:
            skipped.append(acc); continue
        ca = parse_ca(pdb)
        ref = dssp_states(pdb)
        pred = method(pdb, ca)
        if len(ref) != len(pred):
            skipped.append(acc); continue
        preds.append(pred); refs.append(ref); plddts.append(ca["plddt"])
        conf = ca["plddt"] >= 70
        per_protein.append({"accession": acc, "label": label, "organism": org, "length": len(ref),
                            "q3_confident": round(float(np.mean(np.array(list(pred))[conf] == np.array(list(ref))[conf])), 4) if conf.any() else None,
                            "helix_app": pred.count("H") / len(pred), "helix_dssp": ref.count("H") / len(ref),
                            "strand_app": pred.count("E") / len(pred), "strand_dssp": ref.count("E") / len(ref)})
    P, R, B = "".join(preds), "".join(refs), np.concatenate(plddts)
    pp = pd.DataFrame(per_protein)
    res = {
        "date": dt.date.today().isoformat(), "reference": "PyDSSP 0.9.1 (numpy), 3-state DSSP reimplementation",
        "n_proteins": len(per_protein), "n_skipped_no_model": len(skipped), "n_residues": len(R),
        "confident_residues": scores(P, R, B >= 70), "all_residues": scores(P, R, np.ones(len(R), bool)),
        "content_correlation": {"helix_r": round(float(np.corrcoef(pp.helix_app, pp.helix_dssp)[0, 1]), 3),
                                "strand_r": round(float(np.corrcoef(pp.strand_app, pp.strand_dssp)[0, 1]), 3)},
        "per_class_q3_confident": {k: round(float(v), 4) for k, v in pp.groupby("label").q3_confident.mean().items()},
        "proteins": per_protein,
    }
    return res


def to_markdown(res: dict, title: str) -> str:
    c, a = res["confident_residues"], res["all_residues"]
    pct = lambda x: "–" if x is None else f"{100 * x:.1f}%"
    lines = [f"# {title}", "", f"{res['date']} · reference: {res['reference']} · {res['n_proteins']} AlphaFold models "
             f"({res['n_residues']:,} residues; {res['n_skipped_no_model']} proteins without a model skipped)", "",
             "| | Q3 agreement | helix F1 | strand F1 | loop F1 |", "|---|---|---|---|---|",
             f"| Confident residues (pLDDT ≥ 70, n={c['n']:,}) | **{pct(c['q3'])}** | {pct(c['H']['f1'])} | {pct(c['E']['f1'])} | {pct(c['C']['f1'])} |",
             f"| All residues (n={a['n']:,}) | {pct(a['q3'])} | {pct(a['H']['f1'])} | {pct(a['E']['f1'])} | {pct(a['C']['f1'])} |", "",
             f"Per-protein content correlation with DSSP: helix r = {res['content_correlation']['helix_r']}, "
             f"strand r = {res['content_correlation']['strand_r']}", "",
             "Confusion (confident residues; rows = DSSP, columns = app):", "", "| DSSP \\ app | H | E | C |", "|---|---|---|---|"]
    for s in STATES:
        lines.append(f"| {s} | " + " | ".join(f"{c['confusion'][s][t]:,}" for t in STATES) + " |")
    lines += ["", "Q3 by protein class (confident residues): " +
              ", ".join(f"{k} {pct(v)}" for k, v in res["per_class_q3_confident"].items()), ""]
    return "\n".join(lines)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--per-class", type=int, default=25)
    args = ap.parse_args()
    from protein_analyzer.structure import structure_features

    out = Path(__file__).resolve().parent
    methods = {
        "app (DSSP port, proline-aware)": lambda pdb, ca: structure_features(pdb)["secondary_structure"],
        "old C-alpha method (P-SEA-like)": lambda pdb, ca: secondary_structure(ca["xyz"], ca["plddt"]),
    }
    report = ["# Secondary-structure benchmark against DSSP", "",
              "Each method is compared residue by residue with PyDSSP (3-state DSSP reimplementation, >97% agreement "
              "with the original DSSP according to its authors) on AlphaFold models of a stratified random sample of "
              "proteins from this app's dataset.", ""]
    results = {}
    for name, fn in methods.items():
        res = run(args.per_class, method=fn)
        results[name] = res
        report.append(to_markdown(res, name).replace("# ", "## ", 1))
    (out / "ss_benchmark.json").write_text(json.dumps(results, indent=1))
    (out / "ss_benchmark.md").write_text("\n".join(report), encoding="utf-8")
    print("\n".join(report))
