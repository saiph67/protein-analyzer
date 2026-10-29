"""Protein Analyzer web app. Start with:  python app.py   ->  http://127.0.0.1:5000"""
from __future__ import annotations

import datetime as dt
import json
import re
import threading
import time
import webbrowser
from collections import defaultdict, deque
from pathlib import Path

from flask import Flask, Response, abort, jsonify, render_template, request, send_from_directory

from protein_analyzer import config, labtools, ml, msa
from protein_analyzer.analyzer import MAX_SEQUENCE_LENGTH, SequenceError, analyze, clean_sequence, read_fasta
from protein_analyzer.compare import MAX_ALIGN_LENGTH, align, nearest_neighbor_similarity
from protein_analyzer.ebi import DATABASES as BLAST_DATABASES, BlastClient
from protein_analyzer.embeddings import EmbeddingUnavailable, esm_available
from protein_analyzer.explain import explain as explain_prediction
from protein_analyzer.jobs import JobManager
from protein_analyzer.mutations import missense_matrix, mutate, parse_mutation
from protein_analyzer import provenance
from protein_analyzer.projects import ProjectStore, project_provenance, to_markdown as project_markdown
from protein_analyzer.structure import structure_features, superpose
from protein_analyzer.uniprot import UniProtClient, UniProtError, is_accession, to_fasta

MAX_RECORDS = 200
SERVED_FILES = {"report.md", "dataset.csv", "property_boxplots.png"}
SERVED_PREFIXES = ("confusion_matrix_", "accuracy_by_similarity_", "feature_importance_", "pca_map_",
                   "report_", "metrics_")


class InputError(ValueError):
    pass


def create_app(client: UniProtClient | None = None, results_dir: str | Path | None = None,
                blast_client: BlastClient | None = None, data_dir: str | Path | None = None) -> Flask:
    app = Flask(__name__, template_folder=str(config.BUNDLE_DIR / "templates"),
                static_folder=str(config.BUNDLE_DIR / "static"))
    app.config["MAX_CONTENT_LENGTH"] = config.MAX_UPLOAD_BYTES
    client = client or UniProtClient(cache_dir=config.CACHE_DIR)
    blast = blast_client or BlastClient()
    results = Path(results_dir or config.RESULTS_DIR)
    results.mkdir(parents=True, exist_ok=True)
    jobs = JobManager(results / "jobs.json")
    projects = ProjectStore(Path(data_dir or config.DATA_DIR) / "projects.json")
    app.extensions.update(pa_client=client, pa_results=results, pa_jobs=jobs, pa_projects=projects)
    dataset_cache: dict = {}

    @app.context_processor
    def _globals():
        return {"public": config.PUBLIC, "blast_databases": BLAST_DATABASES, "app_version": provenance.__version__}

    # ------------------------------------------------------------ public-server safeguards
    # On a shared server (PA_PUBLIC=1): no dataset/training jobs, no projects, and a per-visitor rate limit
    # on endpoints that call UniProt/AlphaFold/EBI so one visitor can't exhaust the shared quota.
    hits: dict[str, deque] = defaultdict(deque)
    hits_lock = threading.Lock()
    RATE_LIMITED = ("/api/search", "/api/protein", "/api/structure", "/api/compare", "/api/mutate", "/api/blast",
                    "/api/msa", "/api/classify", "/api/explain", "/api/similar", "/api/superpose",
                    "/api/lab/digest", "/api/lab/a280")

    @app.before_request
    def _public_guard():
        if not config.PUBLIC:
            return None
        if request.path.startswith("/api/jobs/") and request.method == "POST":
            return jsonify(error="Building datasets and training are turned off on this shared server."), 403
        if request.path.startswith(RATE_LIMITED):
            ip, now = request.remote_addr or "?", time.time()
            with hits_lock:
                q = hits[ip]
                while q and q[0] < now - 60:
                    q.popleft()
                if len(q) >= config.PUBLIC_RATE_PER_MIN:
                    return jsonify(error="Too many requests from your computer; please wait a minute."), 429
                q.append(now)
        return None

    # ------------------------------------------------------------ helpers
    def get_json() -> dict:
        data = request.get_json(silent=True)
        return data if isinstance(data, dict) else {}

    def uploaded_text() -> str | None:
        f = request.files.get("file")
        if not f:
            return None
        try:
            return f.read().decode("utf-8")
        except UnicodeDecodeError:
            raise InputError("The uploaded file is not a text/FASTA file.")

    def records_from_text(text: str) -> list[tuple[str, str]]:
        if not text or not text.strip():
            raise InputError("Paste a sequence or FASTA text, or upload a file.")
        recs = read_fasta(text)
        if not recs:
            raise InputError("No sequences found.")
        if len(recs) > MAX_RECORDS:
            raise InputError(f"Too many records ({len(recs)}); the limit is {MAX_RECORDS}.")
        return recs

    def resolve_one(value: str, label: str) -> dict:
        """Accession or sequence/FASTA -> {'name', 'sequence', 'accession'?}."""
        v = (value or "").strip()
        if not v:
            raise InputError(f"{label}: enter an accession or a sequence.")
        if is_accession(v):
            e = client.fetch_entry(v)
            return {"name": f"{e['accession']} {e['name'] or ''}".strip(), "sequence": e["sequence"],
                    "accession": e["accession"]}
        name, seq = records_from_text(v)[0]
        try:
            return {"name": name, "sequence": clean_sequence(seq)}
        except SequenceError as e:
            raise InputError(f"{label}: {e}")

    def load_dataset():
        path = results / "dataset.csv"
        if not path.exists():
            return None
        key = path.stat().st_mtime
        if dataset_cache.get("key") != key:
            dataset_cache.update(key=key, df=ml.load_dataset(path))
        return dataset_cache["df"]

    def usable_models():
        models = ml.available_models(results)
        has_esm = esm_available()
        for m in models:
            m["usable"] = has_esm or not m["needs_esm"]
        return models

    # ------------------------------------------------------------ errors
    @app.errorhandler(InputError)
    @app.errorhandler(SequenceError)
    def _bad_input(e):
        return jsonify(error=str(e)), 400

    @app.errorhandler(UniProtError)
    def _uniprot(e):
        return jsonify(error=str(e), source="uniprot"), e.status

    @app.errorhandler(EmbeddingUnavailable)
    def _no_esm(e):
        return jsonify(error=str(e)), 400

    @app.errorhandler(413)
    def _too_big(e):
        return jsonify(error=f"Upload too large; the limit is {config.MAX_UPLOAD_BYTES // 1024 // 1024} MB."), 413

    @app.errorhandler(404)
    def _nf(e):
        if request.path.startswith("/api/"):
            return jsonify(error="Not found."), 404
        return render_template("error.html", message="Page not found."), 404

    # ------------------------------------------------------------ pages
    @app.get("/")
    def home_page():
        return render_template("home.html", page="home")

    @app.get("/analyze")
    def analyze_page():
        return render_template("analyze.html", page="analyze")

    @app.get("/search")
    def search_page():
        return render_template("search.html", page="search")

    @app.get("/protein/<acc>")
    def protein_page(acc):
        return render_template("protein.html", page="search", accession=acc.upper())

    @app.get("/compare")
    def compare_page():
        return render_template("compare.html", page="compare", max_align=MAX_ALIGN_LENGTH)

    @app.get("/ml")
    def ml_page():
        return render_template("ml.html", page="ml", classes=ml.DEFAULT_CLASSES, backends=ml.BACKENDS,
                               esm=esm_available(), esm_model=config.ESM_MODEL)

    @app.get("/classify")
    def classify_page():
        return render_template("classify.html", page="classify")

    # ------------------------------------------------------------ API
    @app.post("/api/analyze")
    def api_analyze():
        text = uploaded_text()
        if text is None:
            text = get_json().get("text") or request.form.get("text", "")
        out, errors = [], []
        for name, seq in records_from_text(text):
            try:
                out.append(analyze(seq, name))
            except SequenceError as e:
                errors.append({"name": name, "error": str(e)})
        if not out and errors:
            raise InputError(errors[0]["error"] if len(errors) == 1 else
                             f"No valid sequences. First problem ({errors[0]['name']}): {errors[0]['error']}")
        return jsonify(records=out, errors=errors, max_length=MAX_SEQUENCE_LENGTH)

    @app.get("/api/search")
    def api_search():
        q = (request.args.get("q") or "").strip()
        if not q:
            raise InputError("Enter a search term or choose a filter.")
        if len(q) > 1000:
            raise InputError("Query too long.")
        size = request.args.get("size", 25, type=int)
        return jsonify(client.search(q, size=max(1, min(size, 100)), cursor=request.args.get("cursor") or None))

    @app.get("/api/protein/<acc>")
    def api_protein(acc):
        e = client.fetch_entry(acc)
        e["analysis"] = analyze(e["sequence"], e["accession"]) if e.get("sequence") else None
        return jsonify(e)

    @app.get("/api/protein/<acc>/fasta")
    def api_fasta(acc):
        e = client.fetch_entry(acc)
        return Response(to_fasta([e]), mimetype="text/plain",
                        headers={"Content-Disposition": f"attachment; filename={e['accession']}.fasta"})

    @app.get("/api/structure/<acc>")
    def api_structure(acc):
        return jsonify(client.fetch_alphafold(acc))

    @app.post("/api/compare")
    def api_compare():
        d = get_json()
        mode = d.get("mode", "global")
        if mode not in ("global", "local"):
            raise InputError("Mode must be global or local.")
        a, b = resolve_one(d.get("a"), "Protein A"), resolve_one(d.get("b"), "Protein B")
        try:
            r = align(a["sequence"], b["sequence"], mode=mode)
        except SequenceError:
            raise
        except ValueError as e:
            raise InputError(str(e))
        r.update(name_a=a["name"], name_b=b["name"], accession_a=a.get("accession"), accession_b=b.get("accession"))
        return jsonify(r)

    @app.post("/api/similar")
    def api_similar():
        d = get_json()
        q = resolve_one(d.get("query"), "Query")
        df = load_dataset()
        if df is None:
            raise InputError("No dataset yet - build one on the Dataset & ML page first.")
        top = max(1, min(int(d.get("top", 10)), 50))
        hits = nearest_neighbor_similarity([q["sequence"]], list(df.sequence), k=3, top=top)[0]
        rows = [{"accession": df.accession.iat[i], "name": df.name.iat[i], "organism": df.organism.iat[i],
                 "label": df.label.iat[i], "split": df.split.iat[i], "length": int(df.length.iat[i]),
                 "similarity": round(s, 4)} for i, s in hits]
        return jsonify(query=q["name"], hits=rows, dataset_size=int(len(df)))

    @app.post("/api/classify")
    def api_classify():
        text = uploaded_text()
        d = get_json() if text is None else {"text": text, "backend": request.form.get("backend")}
        text = (d.get("text") or "").strip()
        if not text:
            raise InputError("Paste a sequence, FASTA or UniProt accession, or upload a FASTA file.")
        if is_accession(text):
            r = resolve_one(text, "Input")
            names, seqs = [r["name"]], [r["sequence"]]
        else:
            recs = records_from_text(text)
            names, seqs = [n for n, _ in recs], []
            for n, s in recs:
                try:
                    seqs.append(clean_sequence(s))
                except SequenceError as e:
                    raise InputError(f"{n}: {e}")
        models = [m for m in usable_models() if m["usable"]]
        if not models:
            raise InputError("No trained model yet. Build a dataset and train on the Dataset & ML page first.")
        backend = d.get("backend") or models[0]["backend"]
        if backend not in [m["backend"] for m in models]:
            raise InputError(f"Model '{backend}' is not trained or not available here.")
        preds = ml.predict(seqs, ml.model_path(results, backend), names)
        return jsonify(predictions=preds, threshold=ml.LOW_SIMILARITY, backend=backend,
                       other_label=ml.OTHER)

    @app.get("/api/ml/status")
    def api_ml_status():
        def read(name):
            p = results / name
            return json.loads(p.read_text()) if p.exists() else None
        metrics = {}
        for b in ml.BACKENDS:
            m = read(f"metrics_{b}.json")
            if m:
                m.pop("dataset", None)
                m["figures"] = [f for f in m.get("figures", []) if (results / f).exists()]
                metrics[b] = m
        return jsonify(dataset=read("dataset_summary.json"), metrics=metrics, models=usable_models(),
                       has_report=(results / "report.md").exists(), esm_available=esm_available(),
                       running=jobs.running(), last_job=jobs.latest())

    @app.post("/api/jobs/dataset")
    def api_job_dataset():
        d = get_json()
        chosen = d.get("classes") or list(ml.DEFAULT_CLASSES)
        include_other = bool(d.get("include_other", True))
        unknown = [c for c in chosen if c not in ml.DEFAULT_CLASSES]
        if unknown or len(chosen) + include_other < 2:
            raise InputError("Choose at least two of the listed classes (\"other\" counts as one).")
        try:
            per_class = int(d.get("per_class", 400)); holdout = int(d.get("holdout", 100))
        except (TypeError, ValueError):
            raise InputError("Sizes must be whole numbers.")
        if not (20 <= per_class <= 2000 and 5 <= holdout <= 500):
            raise InputError("Per-class must be 20-2000 and hold-out 5-500.")
        cutoff = str(d.get("cutoff", "2018-12-31"))
        try:
            dt.date.fromisoformat(cutoff)
        except ValueError:
            raise InputError("Cutoff must be a date like 2018-12-31.")
        try:
            jid = jobs.start("dataset", ml.build_dataset, client, {c: ml.DEFAULT_CLASSES[c] for c in chosen},
                             per_class, holdout, cutoff, results, include_other=include_other)
        except RuntimeError as e:
            return jsonify(error=str(e)), 409
        return jsonify(job_id=jid), 202

    @app.post("/api/jobs/train")
    def api_job_train():
        if not (results / "dataset.csv").exists():
            raise InputError("Build a dataset first.")
        d = get_json()
        backend = d.get("backend", "handcrafted")
        if backend not in ml.BACKENDS:
            raise InputError("Unknown model type.")
        if backend == "esm2" and not esm_available():
            raise EmbeddingUnavailable("ESM-2 needs PyTorch and transformers installed (see README).")
        trees = max(50, min(int(d.get("trees", 500)), 2000))
        try:
            jid = jobs.start("train", lambda progress: _train_summary(ml.train_and_evaluate(
                results / "dataset.csv", results, n_estimators=trees, align_identity=bool(d.get("align", True)),
                backend=backend, progress=progress)))
        except RuntimeError as e:
            return jsonify(error=str(e)), 409
        return jsonify(job_id=jid), 202

    @app.get("/api/jobs/<jid>")
    def api_job(jid):
        j = jobs.get(jid)
        if not j:
            abort(404)
        return jsonify(j)

    @app.post("/api/jobs/<jid>/cancel")
    def api_job_cancel(jid):
        if not jobs.get(jid):
            abort(404)
        if not jobs.cancel(jid):
            raise InputError("That job is not running.")
        return jsonify(ok=True)

    # ------------------------------------------------------------ mutations
    @app.get("/mutate")
    def mutate_page():
        return render_template("mutate.html", page="mutate")

    @app.get("/api/protein/<acc>/variants")
    def api_variants(acc):
        return jsonify(variants=client.fetch_variants(acc))

    @app.get("/api/protein/<acc>/domains")
    def api_domains(acc):
        return jsonify(domains=client.fetch_pfam(acc))

    @app.get("/api/protein/<acc>/alphamissense")
    def api_alphamissense(acc):
        e = client.fetch_entry(acc)
        am = client.fetch_alphamissense(e["accession"])
        out = missense_matrix(am["scores"], e["sequence"])
        out.update(accession=e["accession"], source=am["source"])
        return jsonify(out)

    @app.post("/api/mutate")
    def api_mutate():
        d = get_json()
        acc = (d.get("accession") or "").strip().upper()
        if acc:
            e = client.fetch_entry(acc)
            seq, acc = e["sequence"], e["accession"]
        else:
            seq = (d.get("sequence") or "").strip()
            if seq.startswith(">"):
                seq = records_from_text(seq)[0][1]
        try:
            if d.get("mutation"):
                frm, pos, to = parse_mutation(d["mutation"])
            else:
                frm, pos, to = None, int(d.get("position", 0)), str(d.get("to", ""))
            r = mutate(seq, pos, to, frm)
        except (TypeError, ValueError) as ex:
            if isinstance(ex, SequenceError):
                raise
            raise InputError(str(ex))
        r["accession"] = acc or None
        r["alphamissense"] = r["alphamissense_class"] = None
        r["known_variants"] = []
        if acc:
            try:
                am = client.fetch_alphamissense(acc)
                r["alphamissense"] = am["scores"].get(r["mutation"])
                r["alphamissense_class"] = am["classes"].get(r["mutation"])
            except UniProtError:
                pass  # not human / no model: shown as unavailable
            try:
                r["known_variants"] = [v for v in client.fetch_variants(acc) if v["position"] == r["position"]]
            except UniProtError:
                pass
        return jsonify(r)

    # ------------------------------------------------------------ BLAST (EBI) and multiple alignment
    @app.get("/blast")
    def blast_page():
        return render_template("blast.html", page="search")

    @app.post("/api/blast")
    def api_blast():
        d = get_json()
        q = resolve_one(d.get("query"), "Query")
        if len(q["sequence"]) < 10:
            raise InputError("The query is too short for a BLAST search (at least 10 residues).")
        if len(q["sequence"]) > 10_000:
            raise InputError("BLAST queries are limited to 10,000 residues here.")
        jid = blast.submit(q["sequence"], (d.get("email") or "").strip(), d.get("database", "uniprotkb_swissprot"),
                           max_hits=max(5, min(int(d.get("max_hits", 50)), 250)))
        return jsonify(job_id=jid, query=q["name"], length=len(q["sequence"])), 202

    @app.get("/api/blast/<jid>")
    def api_blast_status(jid):
        st = blast.status(jid)
        out = {"job_id": jid, "status": st}
        if st == "FINISHED":
            out["result"] = blast.result(jid)
        elif st in ("ERROR", "FAILURE", "NOT_FOUND"):
            out["error"] = "EBI reported the search as failed. Check the sequence and try again."
        return jsonify(out)

    @app.get("/align")
    def align_page():
        return render_template("align.html", page="compare", max_seqs=msa.MAX_SEQS, max_len=msa.MAX_LEN)

    @app.post("/api/msa")
    def api_msa():
        text = (get_json().get("text") or "").strip()
        if not text:
            raise InputError("Paste FASTA or list UniProt accessions (one per line).")
        if text.startswith(">"):
            recs = records_from_text(text)
        else:
            tokens = [t for t in re.split(r"[\s,;]+", text) if t]
            bad = [t for t in tokens if not is_accession(t)]
            if bad:
                raise InputError(f"Not an accession: {bad[0]}. Use FASTA for raw sequences.")
            if len(tokens) > msa.MAX_SEQS:
                raise InputError(f"At most {msa.MAX_SEQS} proteins at a time.")
            recs = []
            for t in dict.fromkeys(tok.upper() for tok in tokens):
                e = client.fetch_entry(t)
                gene = e.get("gene") or (e.get("id") or t)
                recs.append((f"{e['accession']} {gene}", e["sequence"]))
        try:
            return jsonify(msa.build(recs))
        except SequenceError:
            raise
        except ValueError as ex:
            raise InputError(str(ex))

    # ------------------------------------------------------------ provenance / methods
    def release():
        fn = getattr(client, "uniprot_release", None)
        try:
            return fn() if fn else {}
        except Exception:
            return {}

    @app.get("/api/provenance")
    def api_provenance():
        kind = request.args.get("kind") or None
        if kind and kind not in provenance.KIND_METHODS:
            raise InputError("Unknown kind.")
        return jsonify(provenance.stamp(kind, release()))

    @app.get("/methods")
    def methods_page():
        return render_template("methods.html", page="methods", prov=provenance.full(), release=release(),
                               ss_benchmark=_read_benchmark())

    # ------------------------------------------------------------ interactions, pathways, structure extras
    @app.get("/api/protein/<acc>/interactions")
    def api_interactions(acc):
        limit = max(5, min(request.args.get("limit", 15, type=int), 40))
        score = max(0.15, min(request.args.get("min_score", 0.4, type=float), 0.99))
        return jsonify(client.fetch_interactions(acc, limit=limit, min_score=score))

    @app.get("/api/protein/<acc>/pathways")
    def api_pathways(acc):
        e = client.fetch_entry(acc)
        return jsonify(pathways=client.fetch_pathways(e["accession"], e.get("organism_id")),
                       organism=e.get("organism"))

    @app.get("/api/protein/<acc>/structure-features")
    def api_structure_features(acc):
        s = client.fetch_alphafold(acc)
        try:
            return jsonify(structure_features(s["pdb"]))
        except ValueError as ex:
            raise InputError(str(ex))

    @app.get("/structures")
    def structures_page():
        return render_template("structures.html", page="compare")

    @app.post("/api/superpose")
    def api_superpose():
        d = get_json()
        accs = [(d.get(k) or "").strip().upper() for k in ("a", "b")]
        for a in accs:
            if not is_accession(a):
                raise InputError("Enter two UniProt accessions (structures come from AlphaFold DB).")
        ea, eb = client.fetch_entry(accs[0]), client.fetch_entry(accs[1])
        if max(len(ea["sequence"]), len(eb["sequence"])) > MAX_ALIGN_LENGTH:
            raise InputError(f"Superposition is limited to proteins of up to {MAX_ALIGN_LENGTH:,} residues.")
        sa, sb = client.fetch_alphafold(ea["accession"]), client.fetch_alphafold(eb["accession"])
        try:
            r = superpose(sa["pdb"], sb["pdb"], ea["sequence"], eb["sequence"])
        except ValueError as ex:
            raise InputError(str(ex))
        r.update(pdb_a=sa["pdb"], name_a=f"{ea['accession']} {ea['name'] or ''}".strip(),
                 name_b=f"{eb['accession']} {eb['name'] or ''}".strip(), accession_a=ea["accession"],
                 accession_b=eb["accession"], plddt_a=sa["mean_plddt"], plddt_b=sb["mean_plddt"],
                 features_a=structure_features(sa["pdb"]), features_b=structure_features(sb["pdb"]))
        return jsonify(r)

    # ------------------------------------------------------------ lab toolkit
    @app.get("/lab")
    def lab_page():
        return render_template("lab.html", page="lab", enzymes={k: v[0] for k, v in labtools.ENZYMES.items()},
                               codes={k: v["name"] for k, v in labtools.CODES.items()})

    def lab_protein(d) -> tuple[str, str]:
        text = (d.get("text") or "").strip()
        if not text:
            raise InputError("Enter a UniProt accession or a protein sequence.")
        r = resolve_one(text, "Protein")
        return r["name"], r["sequence"]

    @app.post("/api/lab/translate")
    def api_translate():
        d = get_json()
        try:
            dna = labtools.clean_dna(d.get("dna", ""))
            code = int(d.get("code", 1))
            if code not in labtools.CODES:
                raise ValueError("Unknown genetic code.")
            min_aa = max(10, min(int(d.get("min_aa", 50)), 5000))
        except (TypeError, ValueError) as ex:
            raise InputError(str(ex))
        orfs = labtools.find_orfs(dna, min_aa, code, bool(d.get("alt_starts")))
        return jsonify(stats=labtools.dna_stats(dna), frames=labtools.six_frames(dna, code), orfs=orfs[:100],
                       n_orfs=len(orfs), code=code, code_name=labtools.CODES[code]["name"])

    @app.post("/api/lab/digest")
    def api_digest():
        d = get_json()
        name, seq = lab_protein(d)
        try:
            r = labtools.digest(seq, d.get("enzyme", "trypsin"), int(d.get("missed", 0)),
                                int(d.get("min_len", 6)), int(d.get("max_len", 50)),
                                float(d.get("min_mass", 500)), float(d.get("max_mass", 5000)), bool(d.get("cam", True)))
        except (TypeError, ValueError) as ex:
            if isinstance(ex, SequenceError):
                raise
            raise InputError(str(ex))
        custom = str(d.get("name") or "").strip()[:100]
        r["name"] = custom if custom and not is_accession(d.get("text", "").strip()) else name
        return jsonify(r)

    @app.post("/api/lab/a280")
    def api_a280():
        d = get_json()
        name, seq = lab_protein(d)
        r = labtools.extinction(seq)
        r.update(name=name, length=len(seq))
        if d.get("a280") not in (None, ""):
            eps = r["epsilon_cystines"] if d.get("cystines", True) else r["epsilon_reduced"]
            try:
                r["concentration"] = labtools.concentration(float(d["a280"]), eps, r["molecular_weight"],
                                                            float(d.get("path_cm", 1) or 1), float(d.get("dilution", 1) or 1))
            except (TypeError, ValueError) as ex:
                raise InputError(str(ex))
        return jsonify(r)

    # ------------------------------------------------------------ explanations
    @app.post("/api/explain")
    def api_explain():
        d = get_json()
        r = resolve_one(d.get("text"), "Input")
        models = [m for m in usable_models() if m["usable"]]
        if not models:
            raise InputError("No trained model yet.")
        backend = d.get("backend") or models[0]["backend"]
        if backend not in [m["backend"] for m in models]:
            raise InputError(f"Model '{backend}' is not available.")
        return jsonify(explain_prediction(r["sequence"], ml.model_path(results, backend), r["name"]))

    # ------------------------------------------------------------ lessons + projects
    @app.get("/learn")
    def learn_page():
        return render_template("learn.html", page="learn")

    def need_projects():
        if config.PUBLIC:
            raise InputError("Projects are turned off on shared (public) servers, because everyone would "
                             "see everyone's projects. Run the app on your own computer to use them.")

    @app.get("/projects")
    def projects_page():
        return render_template("projects.html", page="projects", public=config.PUBLIC)

    @app.get("/projects/<pid>/report")
    def project_report(pid):
        need_projects()
        try:
            p = projects.get(pid)
        except KeyError:
            abort(404)
        return render_template("project_report.html", page="projects", p=p, pv=project_provenance(p),
                               when=dt.datetime.fromtimestamp)

    @app.get("/api/projects")
    def api_projects():
        need_projects()
        return jsonify(projects=projects.list())

    @app.post("/api/projects")
    def api_project_create():
        need_projects()
        try:
            return jsonify(projects.create(get_json().get("name", ""))), 201
        except ValueError as ex:
            raise InputError(str(ex))

    @app.get("/api/projects/backup")
    def api_projects_backup():
        need_projects()
        body = json.dumps(projects.backup(), indent=1)
        stamp = dt.date.today().isoformat()
        return Response(body, mimetype="application/json",
                        headers={"Content-Disposition": f"attachment; filename=protein_analyzer_backup_{stamp}.json"})

    @app.post("/api/projects/import")
    def api_projects_import():
        need_projects()
        payload = get_json()
        if not payload:
            try:
                payload = json.loads(uploaded_text() or "null")
            except ValueError:
                raise InputError("That file isn't valid JSON. Choose a project or backup file exported from this app.")
        try:
            created = projects.import_data(payload)
        except ValueError as ex:
            raise InputError(str(ex))
        return jsonify(created=created), 201

    @app.route("/api/projects/<pid>", methods=["GET", "PATCH", "DELETE"])
    def api_project(pid):
        need_projects()
        try:
            if request.method == "GET":
                return jsonify(projects.get(pid))
            if request.method == "DELETE":
                projects.delete(pid)
                return jsonify(ok=True)
            d = get_json()
            return jsonify(projects.update(pid, d.get("name"), d.get("notes")))
        except KeyError:
            abort(404)
        except ValueError as ex:
            raise InputError(str(ex))

    @app.post("/api/projects/<pid>/items")
    def api_project_add(pid):
        need_projects()
        d = get_json()
        kind = d.get("kind", "")
        try:
            return jsonify(projects.add_item(pid, kind, d.get("title", ""), d.get("data") or {}, d.get("note", ""),
                                             provenance.stamp(kind if kind in provenance.KIND_METHODS else None,
                                                              release()))), 201
        except KeyError:
            abort(404)
        except ValueError as ex:
            raise InputError(str(ex))

    @app.route("/api/projects/<pid>/items/<iid>", methods=["PATCH", "DELETE"])
    def api_project_item(pid, iid):
        need_projects()
        try:
            if request.method == "DELETE":
                projects.delete_item(pid, iid)
                return jsonify(ok=True)
            return jsonify(projects.update_item(pid, iid, get_json().get("note", "")))
        except KeyError:
            abort(404)

    @app.get("/api/projects/<pid>/export.<fmt>")
    def api_project_export(pid, fmt):
        need_projects()
        try:
            p = projects.get(pid)
        except KeyError:
            abort(404)
        safe = re.sub(r"[^A-Za-z0-9_-]+", "_", p["name"])[:50] or "project"
        if fmt == "md":
            body, mime = project_markdown(p), "text/markdown"
        elif fmt == "json":
            body, mime = json.dumps(p, indent=2), "application/json"
        else:
            abort(404)
        return Response(body, mimetype=mime, headers={"Content-Disposition": f"attachment; filename={safe}.{fmt}"})

    @app.get("/results/<path:name>")
    def result_file(name):
        ok = name in SERVED_FILES or (name.startswith(SERVED_PREFIXES) and "/" not in name and "\\" not in name
                                      and name.endswith((".png", ".md", ".json")))
        if not ok or not (results / name).exists():
            abort(404)
        return send_from_directory(results, name)

    return app


def _read_benchmark():
    p = config.BUNDLE_DIR / "benchmarks" / "ss_benchmark.json"
    if not p.exists():
        p = config.PROJECT_DIR / "benchmarks" / "ss_benchmark.json"
    try:
        data = json.loads(p.read_text())
    except (OSError, ValueError):
        return None
    return {name: {"q3": r["confident_residues"]["q3"], "n_proteins": r["n_proteins"],
                   "n": r["confident_residues"]["n"], "date": r["date"]} for name, r in data.items()}


def _train_summary(m: dict) -> dict:
    return {"backend": m["backend"], "holdout": m["holdout"], "low_similarity": m["holdout_low_similarity"]}


if __name__ == "__main__":
    application = create_app()
    url = f"http://{config.HOST}:{config.PORT}"
    print(f"Protein Analyzer running at {url}  (Ctrl+C to stop)")
    if config.OPEN_BROWSER:
        try:
            webbrowser.open(url)
        except Exception:
            pass
    application.run(host=config.HOST, port=config.PORT, debug=False, threaded=True)
