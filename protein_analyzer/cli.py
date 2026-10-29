"""Command-line interface: analyze, fetch, compare, dataset, train, predict."""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from . import config
from .analyzer import analyze, plot_hydrophobicity, read_fasta
from .compare import align
from .uniprot import UniProtClient, UniProtError, is_accession, save_fasta


def _client():
    return UniProtClient(cache_dir=config.CACHE_DIR)


def _resolve(arg: str, client=None) -> tuple[str, str]:
    """An accession, a FASTA file path or a raw sequence -> (name, sequence)."""
    p = Path(arg)
    if p.is_file():
        recs = read_fasta(p.read_text())
        return recs[0]
    if is_accession(arg):
        e = (client or _client()).fetch_entry(arg)
        return e["accession"], e["sequence"]
    return "sequence", arg


def _summary(a: dict) -> str:
    return (f"{a['name']}: {a['length']} aa, MW {a['molecular_weight']:.1f} Da, GRAVY {a['gravy']:+.3f} "
            f"({a['gravy_class']}), pI {a['isoelectric_point']:.2f}, charge(pH7) {a['charge_ph7']:+.2f}, "
            f"aromaticity {a['aromaticity']:.3f}, TM segments {a['tm_segments']}")


def main(argv=None):
    ap = argparse.ArgumentParser(prog="protein_analyzer")
    sub = ap.add_subparsers(dest="cmd", required=True)

    s = sub.add_parser("analyze", help="analyze a FASTA file or sequence")
    s.add_argument("input"); s.add_argument("--plot", help="directory for hydropathy PNGs"); s.add_argument("--json", action="store_true")

    s = sub.add_parser("fetch", help="fetch UniProt entries")
    s.add_argument("accessions", nargs="+"); s.add_argument("--analyze", action="store_true")
    s.add_argument("--fasta", help="save sequences to this FASTA file")

    s = sub.add_parser("compare", help="align two sequences/accessions/FASTA files")
    s.add_argument("a"); s.add_argument("b"); s.add_argument("--mode", choices=["global", "local"], default="global")

    s = sub.add_parser("dataset", help="build a labelled dataset from UniProt")
    s.add_argument("--per-class", type=int, default=400); s.add_argument("--holdout", type=int, default=100)
    s.add_argument("--cutoff", default="2018-12-31"); s.add_argument("--classes", help="comma-separated subset")
    s.add_argument("--no-other", action="store_true", help="do not add the 'other' (none of these) class")

    s = sub.add_parser("train", help="train and evaluate on results/dataset.csv")
    s.add_argument("--trees", type=int, default=500); s.add_argument("--no-align", action="store_true")
    s.add_argument("--backend", choices=["handcrafted", "esm2"], default="handcrafted")

    s = sub.add_parser("predict", help="predict function class")
    s.add_argument("input", nargs="+", help="FASTA file, sequence or accession")
    s.add_argument("--backend", choices=["handcrafted", "esm2"], help="model to use (default: best trained)")

    args = ap.parse_args(argv)
    try:
        return _run(args)
    except (UniProtError, ValueError, FileNotFoundError) as e:
        print(f"error: {e}", file=sys.stderr)
        return 1


def _run(args):
    from . import ml

    results = Path(config.RESULTS_DIR)
    progress = lambda f, msg: print(f"[{f:5.0%}] {msg}", flush=True)
    if args.cmd == "analyze":
        p = Path(args.input)
        recs = read_fasta(p.read_text()) if p.is_file() else read_fasta(args.input)
        for name, seq in recs:
            a = analyze(seq, name)
            print(json.dumps({k: v for k, v in a.items() if k not in ("hydropathy_profile", "sequence")}) if args.json else _summary(a))
            if args.plot:
                Path(args.plot).mkdir(parents=True, exist_ok=True)
                plot_hydrophobicity(a, str(Path(args.plot) / f"{name.split()[0].replace('|', '_')}.png"))
    elif args.cmd == "fetch":
        c = _client()
        entries = [c.fetch_entry(a) for a in args.accessions]
        for e in entries:
            print(f"{e['accession']} {e['id']} | {e['name']} | {e['organism']} | {e['length']} aa | "
                  f"EC {', '.join(e['ec']) or '-'} | created {e['date_created']}")
            if e["function"]:
                print("  function:", e["function"][0][:200])
            print(f"  GO: F={len(e['go']['F'])} P={len(e['go']['P'])} C={len(e['go']['C'])}; "
                  f"keywords: {', '.join(k['name'] for k in e['keywords'][:8])}")
            if args.analyze:
                print("  " + _summary(analyze(e["sequence"], e["accession"])))
        if args.fasta:
            save_fasta(entries, args.fasta)
    elif args.cmd == "compare":
        c = _client()
        (na, sa), (nb, sb) = _resolve(args.a, c), _resolve(args.b, c)
        r = align(sa, sb, mode=args.mode)
        print(f"{na} vs {nb} ({r['mode']}): score {r['score']}, identity {r['identity']}%, similarity "
              f"{r['similarity']}%, gaps {r['gaps']}, coverage {r['coverage_a']}% / {r['coverage_b']}%")
        for i in range(0, r["length"], 60):
            print(r["aligned_a"][i:i + 60]); print(r["match_line"][i:i + 60]); print(r["aligned_b"][i:i + 60]); print()
    elif args.cmd == "dataset":
        classes = ml.DEFAULT_CLASSES
        if args.classes:
            classes = {k: ml.DEFAULT_CLASSES[k] for k in args.classes.split(",")}
        info = ml.build_dataset(_client(), classes, args.per_class, args.holdout, args.cutoff, results,
                                include_other=not args.no_other, progress=progress)
        print(json.dumps(info["counts"], indent=2))
        for w in info["warnings"]:
            print("warning:", w)
    elif args.cmd == "train":
        m = ml.train_and_evaluate(results / "dataset.csv", results, n_estimators=args.trees,
                                  align_identity=not args.no_align, backend=args.backend, progress=progress)
        print(json.dumps({"holdout": m["holdout"], "low_similarity": m["holdout_low_similarity"],
                          "below_30pct_identity": m["holdout_below_30pct_identity"]}, indent=2))
        print(f"report: {results / 'report.md'}")
    elif args.cmd == "predict":
        c = _client()
        names, seqs = [], []
        for item in args.input:
            p = Path(item)
            recs = read_fasta(p.read_text()) if p.is_file() else [_resolve(item, c)]
            for n, s in recs:
                names.append(n); seqs.append(s)
        models = ml.available_models(results)
        backend = args.backend or (models[0]["backend"] if models else "handcrafted")
        for r in ml.predict(seqs, ml.model_path(results, backend), names):
            warn = "" if r["reliable"] else "  [UNRELIABLE: dissimilar to training data]"
            print(f"{r['name']}: {r['prediction']} ({r['confidence']:.0%}); 2nd {r['top2'][1]['label']} "
                  f"({r['top2'][1]['probability']:.0%}); nn-sim {r['nn_similarity']:.2f}{warn}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
