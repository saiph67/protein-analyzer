"""Offline tests for the protein_analyzer package. Run with `pytest` or `python tests/test_all.py`."""
from __future__ import annotations

import json
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / "tests"))

import pytest  # noqa: E402

from fixtures import BACR_SEQ, HBA_SEQ, HBB_SEQ, raw_entry  # noqa: E402
from protein_analyzer import ml  # noqa: E402
from protein_analyzer.analyzer import SequenceError, analyze, clean_sequence, read_fasta  # noqa: E402
from protein_analyzer.compare import MAX_ALIGN_LENGTH, align, kmer_similarity, nearest_neighbor_similarity  # noqa: E402
from protein_analyzer.features import N_FEATURES, featurize, featurize_many  # noqa: E402
from protein_analyzer.jobs import JobManager  # noqa: E402
from protein_analyzer.uniprot import (UniProtClient, UniProtBadRequest, UniProtUnavailable, is_accession,  # noqa: E402
                                      parse_entry, to_fasta)


# ------------------------------------------------------------------ analyzer
def test_analyze_hemoglobin_matches_reference():
    r = analyze(HBA_SEQ, "HBA")
    assert r["length"] == 142
    assert abs(r["molecular_weight"] - 15258) < 2          # UniProt: 15,258 Da
    assert r["isoelectric_point"] == pytest.approx(8.72, abs=0.01)  # ExPASy Compute pI/Mw: 8.72
    assert sum(v["count"] for v in r["composition"].values()) == 142
    assert r["n_tm_segments"] == 0


def test_bacteriorhodopsin_has_tm_segments():
    r = analyze(BACR_SEQ)
    assert r["gravy"] > 0.3 and r["gravy_class"].startswith("hydrophobic")
    assert r["n_tm_segments"] >= 4
    assert len(r["hydropathy_profile"]) == len(BACR_SEQ) - 18


def test_invalid_residues_are_reported_with_positions():
    with pytest.raises(SequenceError, match=r"'1'|'@'"):
        clean_sequence("MKT@L")
    with pytest.raises(SequenceError, match="position"):
        clean_sequence("MKT#L")
    with pytest.raises(SequenceError, match="empty"):
        clean_sequence("   ")
    assert clean_sequence("mk t\n12 l*") == "MKTL"
    with pytest.raises(SequenceError, match="maximum"):
        clean_sequence("A" * 35_001)


def test_read_fasta_multi_and_bare():
    recs = read_fasta(">a desc\nMKT\nLL\n\n>b\nAAA\n")
    assert recs == [("a desc", "MKTLL"), ("b", "AAA")]
    assert read_fasta("MKTL LL") == [("sequence_1", "MKTLLL")]
    assert read_fasta("") == []


# ------------------------------------------------------------------ compare
def test_align_identical_and_hemoglobins():
    r = align(HBA_SEQ, HBA_SEQ)
    assert r["identity"] == 100 and r["gaps"] == 0 and r["coverage_a"] == 100
    r = align(HBA_SEQ, HBB_SEQ, "global")
    assert 40 <= r["identity"] <= 46                         # textbook ~43%
    assert r["aligned_a"].replace("-", "") == HBA_SEQ and r["aligned_b"].replace("-", "") == HBB_SEQ


def test_align_gap_and_local():
    r = align("ACDEFGHIKLMNPQRSTVWY", "ACDEFGHKLMNPQRSTVWY")
    assert r["gaps"] == 1 and r["identical"] == 19
    r = align("WWWWWW" + "HEAGAWGHEE" + "PPPPPP", "HEAGAWGHEE", mode="local")
    assert r["aligned_a"] == "HEAGAWGHEE" and r["range_a"] == [7, 16] and r["identity"] == 100


def test_align_affine_gap_score():
    # one gap of length 3 costs 11 + 3*1 = 14 (BLAST convention)
    left, right = "MKTAYIAKQRQ", "ISFVKSHFSRQ"
    r = align(left + "WWW" + right, left + right)
    assert r["gaps"] == 3 and r["aligned_b"] == left + "---" + right
    assert r["score"] == align(left + right, left + right)["score"] - 14


def test_align_length_cap():
    with pytest.raises(ValueError, match="limited"):
        align("A" * (MAX_ALIGN_LENGTH + 1), "AAAA")


def test_kmer_similarity_and_nn():
    assert kmer_similarity(HBA_SEQ, HBA_SEQ) == 1.0
    assert kmer_similarity(HBA_SEQ, BACR_SEQ) < 0.05
    nn = nearest_neighbor_similarity([HBA_SEQ], [BACR_SEQ, HBB_SEQ, HBA_SEQ], top=2)
    assert nn[0][0] == (2, 1.0) and nn[0][1][0] == 1


# ------------------------------------------------------------------ features
def test_features_shape():
    assert N_FEATURES == 434
    x = featurize(HBA_SEQ)
    assert x.shape == (434,)
    assert abs(x[:20].sum() - 1) < 1e-9 and abs(x[20:420].sum() - 1) < 1e-9
    assert featurize_many([HBA_SEQ, BACR_SEQ]).shape == (2, 434)


# ------------------------------------------------------------------ uniprot
def test_parse_entry_live_shape():
    e = parse_entry(raw_entry(ec=["1.2.3.4"]))
    assert e["accession"] == "P69905" and e["gene"] == "HBA1" and e["organism_id"] == 9606
    assert e["function"] == ["Involved in oxygen transport."]
    assert [g["term"] for g in e["go"]["F"]] == ["oxygen binding"]
    assert len(e["go"]["P"]) == 1 and len(e["go"]["C"]) == 1
    assert e["ec"] == ["1.2.3.4"] and e["date_created"] == "1986-07-21" and e["reviewed"] is True
    assert [(f["type"], f["start"], f["end"]) for f in e["features"]] == [
        ("Signal", 1, 10), ("Transmembrane", 20, 42), ("Binding site", 88, 88)]  # unknown-position region skipped
    assert e["features"][2]["ligand"] == "heme b"
    assert to_fasta([e]).startswith(">sp|P69905|P69905_TEST Hemoglobin subunit alpha OS=Homo sapiens GN=HBA1\nMVLSP")


def test_accession_regex():
    for ok in ["P69905", "Q9H0H5", "A0A023GPI8", "P01308-2"]:
        assert is_accession(ok)
    for bad in ["kinase", "P6990", "12345", ""]:
        assert not is_accession(bad)


class _Resp:
    def __init__(self, status, body=None, headers=None):
        self.status_code, self._body, self.headers = status, body or {}, headers or {}
        self.text = json.dumps(self._body)

    def json(self):
        return self._body


class _Session:
    def __init__(self, responses):
        self.responses, self.headers, self.calls = list(responses), {}, []

    def get(self, url, params=None, timeout=None):
        self.calls.append((url, params))
        r = self.responses.pop(0)
        if isinstance(r, Exception):
            raise r
        return r


def _client(responses):
    return UniProtClient(cache_dir=None, session=_Session(responses), min_interval=0)


def test_client_pagination_count_and_search():
    nxt = 'https://rest.uniprot.org/uniprotkb/search?query=x&cursor=abc123&size=2'
    c = _client([
        _Resp(200, {"results": [raw_entry("P69905"), raw_entry("P68871")]}, {"Link": f'<{nxt}>; rel="next"', "X-Total-Results": "3"}),
        _Resp(200, {"results": [raw_entry("P02945")]}, {"X-Total-Results": "3"}),
    ])
    accs = [e["accession"] for e in c.iter_search("x", page_size=2)]
    assert accs == ["P69905", "P68871", "P02945"]
    assert c.session.calls[1] == (nxt, None)  # follows Link verbatim
    c = _client([_Resp(200, {"results": []}, {"X-Total-Results": "21273"})])
    assert c.count("keyword:KW-0418") == 21273 and c.session.calls[0][1]["size"] == 0
    c = _client([_Resp(200, {"results": [raw_entry()]}, {"Link": f'<{nxt}>; rel="next"', "X-Total-Results": "9"})])
    page = c.search("hemoglobin", size=1000)
    assert page["next_cursor"] == "abc123" and page["total"] == 9 and c.session.calls[0][1]["size"] == 500


def test_client_retries_and_errors(monkeypatch):
    import protein_analyzer.uniprot as up
    monkeypatch.setattr(up.time, "sleep", lambda s: None)
    c = _client([_Resp(429, headers={"Retry-After": "1"}), _Resp(503), _Resp(200, raw_entry())])
    assert c.fetch_entry("P69905")["accession"] == "P69905" and len(c.session.calls) == 3
    c = _client([_Resp(400, {"messages": ["'foo' is not a valid search field"]})])
    with pytest.raises(UniProtBadRequest, match="not a valid search field"):
        c.search("foo:bar")
    c = _client([up.requests.ConnectionError("down")] * 5)
    with pytest.raises(UniProtUnavailable):
        c.search("x")
    with pytest.raises(UniProtBadRequest, match="not a valid UniProt accession"):
        _client([]).fetch_entry("hello")


def test_fetch_alphafold(tmp_path):
    api = [{"uniprotAccession": "P69905-2", "pdbUrl": "https://x/iso.pdb"},
           {"uniprotAccession": "P69905", "pdbUrl": "https://x/AF-P69905.pdb", "latestVersion": 6,
            "globalMetricValue": 98.1, "uniprotStart": 1, "uniprotEnd": 142}]
    pdb = _Resp(200); pdb.text = "ATOM ...\nEND\n"
    c = UniProtClient(cache_dir=tmp_path, session=_Session([_Resp(200, api), pdb]), min_interval=0)
    s = c.fetch_alphafold("p69905")
    assert s["pdb"].startswith("ATOM") and s["version"] == 6 and s["mean_plddt"] == 98.1
    assert c.session.calls[1][0] == "https://x/AF-P69905.pdb"  # canonical model chosen over isoform
    assert c.fetch_alphafold("P69905")["pdb"] == s["pdb"] and len(c.session.calls) == 2  # cached
    c = _client([_Resp(404)])
    with pytest.raises(Exception, match="no predicted structure"):
        c.fetch_alphafold("Q99999")


def test_client_follows_merged_entry_and_caches(tmp_path):
    inactive = {"entryType": "Inactive", "primaryAccession": "P00001",
                "inactiveReason": {"inactiveReasonType": "MERGED", "mergeDemergeTo": ["P69905"]}}
    c = UniProtClient(cache_dir=tmp_path, session=_Session([_Resp(200, inactive), _Resp(200, raw_entry())]), min_interval=0)
    assert c.fetch_entry("P00001")["accession"] == "P69905"
    assert c.fetch_entry("P69905")["accession"] == "P69905"  # served from cache: no third call
    assert len(c.session.calls) == 2


# ------------------------------------------------------------------ ml
def test_class_query_syntax():
    q = ml.class_query("KW-0418", ["KW-0645"], "1986-01-01", "2018-12-31")
    assert q == ("keyword:KW-0418 AND reviewed:true AND fragment:false AND length:[50 TO 2000] "
                 "AND date_created:[1986-01-01 TO 2018-12-31] AND NOT keyword:KW-0645")


def test_class_query_other():
    q = ml.class_query(None, ["KW-0418", "KW-0645"], "2019-01-01", "*")
    assert q.startswith("reviewed:true") and q.endswith("AND NOT keyword:KW-0418 AND NOT keyword:KW-0645")


def test_allocate_proportional():
    assert ml._allocate({"a": 900, "b": 100, "c": 0}, 100) == {"a": 90, "b": 10, "c": 0}
    assert ml._allocate({"a": 5, "b": 1000}, 50) == {"a": 0, "b": 50}
    assert sum(ml._allocate({"a": 3, "b": 4, "c": 5}, 7).values()) == 7
    assert ml._allocate({"a": 2, "b": 1}, 100) == {"a": 2, "b": 1}  # capped by availability
    assert ml._allocate({"a": 0}, 10) == {"a": 0}


def test_build_dataset_offline(tmp_path):
    from fixtures import FakeUniProtClient
    info = ml.build_dataset(FakeUniProtClient(), {"kinase": "KW-0418", "protease": "KW-0645"},
                            per_class=20, holdout=5, out_dir=tmp_path, include_other=True)
    assert info["counts"]["kinase"] == {"holdout": 5, "train": 20}
    assert info["counts"]["other"] == {"holdout": 5, "train": 20}
    assert info["classes"]["other"] is None and "taxonomic" in info["sampling"]
    assert set(info["availability"]["kinase"]["train_strata_available"]) == set(ml.TAXON_STRATA)
    assert (tmp_path / "dataset.csv").exists()


def test_train_predict_smoke_on_synthetic(tmp_path):
    from make_synthetic import make_dataset
    make_dataset(tmp_path / "dataset.csv", per_class=30, holdout=8)
    m = ml.train_and_evaluate(tmp_path / "dataset.csv", tmp_path, n_estimators=40, align_identity=True)
    assert set(m["holdout"]) == {"accuracy", "balanced_accuracy", "macro_f1"}
    assert len(m["confusion_matrix"]) == 3 and m["accuracy_by_identity"] and m["backend"] == "handcrafted"
    for f in ["report.md", "report_handcrafted.md", "metrics_handcrafted.json", "model_handcrafted.joblib",
              *m["figures"]]:
        assert (tmp_path / f).exists(), f
    p = ml.predict([BACR_SEQ], ml.model_path(tmp_path, "handcrafted"), ["bR"])[0]
    assert p["prediction"] in m["classes"] and len(p["top2"]) == 2 and isinstance(p["reliable"], bool)
    models = ml.available_models(tmp_path)
    assert [x["backend"] for x in models] == ["handcrafted"]


def test_esm_backend_pipeline_with_fake_embeddings(tmp_path, monkeypatch):
    """Exercises the esm2 code path offline; the embedding function is replaced by a cheap stand-in."""
    import numpy as np
    from make_synthetic import make_dataset
    import protein_analyzer.embeddings as emb

    def fake_embed(seqs, model_name=None, cache_dir=None, progress=None, **kw):
        if progress:
            progress(1.0, "fake embeddings")
        return np.vstack([featurize(s)[:20] * 10 for s in seqs])  # composition as a stand-in

    monkeypatch.setattr(emb, "embed", fake_embed)
    make_dataset(tmp_path / "dataset.csv", per_class=30, holdout=8)
    m = ml.train_and_evaluate(tmp_path / "dataset.csv", tmp_path, n_estimators=30, align_identity=False, backend="esm2")
    assert m["backend"] == "esm2" and "LogisticRegression" in m["model"] and m["n_features"] == 20
    assert "random forest, hand-crafted features" in m["baselines"]
    assert (tmp_path / "model_esm2.joblib").exists() and "ESM-2" in (tmp_path / "report_esm2.md").read_text()
    p = ml.predict([HBA_SEQ], ml.model_path(tmp_path, "esm2"))[0]
    assert p["prediction"] in m["classes"]
    with pytest.raises(ValueError):
        ml.train_and_evaluate(tmp_path / "dataset.csv", tmp_path, backend="nope")


def _wait_job(jm, jid):
    import time
    for _ in range(200):
        if jm.get(jid)["state"] not in ("running", "cancelling"):
            return jm.get(jid)
        time.sleep(0.02)
    raise AssertionError("job did not finish")


def test_job_cancel_and_restart_persistence(tmp_path):
    import threading
    state = tmp_path / "jobs.json"
    jm = JobManager(state)
    go = threading.Event()

    def slow(progress):
        for i in range(500):
            progress(i / 500, f"step {i}")
            go.set()
            threading.Event().wait(0.01)
        return "finished"

    jid = jm.start("train", slow)
    go.wait(2)
    with pytest.raises(RuntimeError, match="already running"):
        jm.start("dataset", slow)
    assert jm.cancel(jid)
    assert _wait_job(jm, jid)["state"] == "cancelled"
    assert not jm.cancel(jid)
    # simulate the app being closed mid-job: a 'running' record on disk becomes 'interrupted'
    import json as _json
    saved = _json.loads(state.read_text())
    saved[0]["state"] = "running"
    state.write_text(_json.dumps(saved))
    jm2 = JobManager(state)
    assert jm2.get(jid)["state"] == "interrupted" and "closed" in jm2.get(jid)["error"]
    assert jm2.latest()["id"] == jid and jm2.running() is None


def test_job_manager():
    jm = JobManager()
    jid = jm.start("t", lambda progress: (progress(0.5, "half"), 42)[1])
    assert _wait_job(jm, jid)["result"] == 42 and jm.get(jid)["state"] == "done"
    jid = jm.start("t", lambda progress: 1 / 0)
    assert _wait_job(jm, jid)["state"] == "error"


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-q"]))
