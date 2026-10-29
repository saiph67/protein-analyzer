"""Flask API tests with the UniProt client mocked (fully offline)."""
from __future__ import annotations

import io
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / "tests"))

import pytest  # noqa: E402

from app import create_app  # noqa: E402
from fixtures import BACR_SEQ, HBA_SEQ, FakeUniProtClient  # noqa: E402
from protein_analyzer.uniprot import UniProtRateLimited, UniProtUnavailable  # noqa: E402


@pytest.fixture
def make(tmp_path):
    def _make(client=None):
        app = create_app(client or FakeUniProtClient(), results_dir=tmp_path)
        app.testing = True
        return app.test_client()
    return _make


@pytest.fixture
def c(make):
    return make()


@pytest.fixture(scope="module")
def trained_results(tmp_path_factory):
    from make_synthetic import make_dataset
    from protein_analyzer import ml
    d = tmp_path_factory.mktemp("results")
    make_dataset(d / "dataset.csv", per_class=25, holdout=6)
    ml.train_and_evaluate(d / "dataset.csv", d, n_estimators=30, align_identity=False)
    return d


def test_pages_render(c):
    for url in ["/", "/analyze", "/search", "/protein/P69905", "/compare", "/ml", "/classify"]:
        r = c.get(url)
        assert r.status_code == 200, url
        assert b"Protein Analyzer" in r.data
    assert c.get("/nope").status_code == 404


def test_analyze_text_fasta_and_upload(c):
    r = c.post("/api/analyze", json={"text": f">a\n{HBA_SEQ}\n>b\n{BACR_SEQ}"})
    d = r.get_json()
    assert r.status_code == 200 and [x["name"] for x in d["records"]] == ["a", "b"]
    assert d["records"][1]["n_tm_segments"] >= 4
    r = c.post("/api/analyze", data={"file": (io.BytesIO(f">up\n{HBA_SEQ}".encode()), "x.fasta")},
               content_type="multipart/form-data")
    assert r.status_code == 200 and r.get_json()["records"][0]["length"] == 142


def test_analyze_errors(c):
    r = c.post("/api/analyze", json={"text": "MKT!!L"})
    assert r.status_code == 400 and "Invalid residue" in r.get_json()["error"]
    assert c.post("/api/analyze", json={"text": ""}).status_code == 400
    r = c.post("/api/analyze", json={"text": ">ok\nMKTLL\n>bad\nMK1@"})
    assert r.status_code == 200 and len(r.get_json()["errors"]) == 1
    big = io.BytesIO(b">x\n" + b"A" * (5 * 1024 * 1024 + 10))
    r = c.post("/api/analyze", data={"file": (big, "big.fa")}, content_type="multipart/form-data")
    assert r.status_code == 413
    r = c.post("/api/analyze", data={"file": (io.BytesIO(b"\xff\xfe\x00bin"), "x.bin")}, content_type="multipart/form-data")
    assert r.status_code == 400


def test_search_and_pagination(c):
    d = c.get("/api/search?q=hemoglobin&size=2").get_json()
    assert d["total"] == 3 and len(d["results"]) == 2 and d["next_cursor"] == "2"
    d2 = c.get(f"/api/search?q=hemoglobin&size=2&cursor={d['next_cursor']}").get_json()
    assert [x["accession"] for x in d2["results"]] == ["P02945"] and d2["next_cursor"] is None
    assert c.get("/api/search?q=").status_code == 400


def test_protein_detail_and_fasta(c):
    d = c.get("/api/protein/P69905").get_json()
    assert d["gene"] == "HBA1" and d["analysis"]["length"] == 142 and d["go"]["F"]
    r = c.get("/api/protein/P69905/fasta")
    assert r.status_code == 200 and r.data.startswith(b">sp|P69905") and "attachment" in r.headers["Content-Disposition"]
    assert c.get("/api/protein/Q99999").status_code == 404
    assert c.get("/api/protein/notanacc").status_code == 400


def test_uniprot_down_and_rate_limited(make):
    c = make(FakeUniProtClient(fail=UniProtUnavailable("Could not reach UniProt")))
    r = c.get("/api/search?q=x")
    assert r.status_code == 503 and r.get_json()["source"] == "uniprot"
    c = make(FakeUniProtClient(fail=UniProtRateLimited("UniProt is rate-limiting")))
    assert "rate-limiting" in c.get("/api/protein/P69905").get_json()["error"]


def test_compare(c):
    d = c.post("/api/compare", json={"a": "P69905", "b": "P68871", "mode": "global"}).get_json()
    assert 40 <= d["identity"] <= 46 and d["accession_a"] == "P69905"
    d = c.post("/api/compare", json={"a": HBA_SEQ, "b": HBA_SEQ[20:60], "mode": "local"}).get_json()
    assert d["identity"] == 100 and d["coverage_b"] == 100
    assert c.post("/api/compare", json={"a": "MKT", "b": "A" * 3001}).status_code == 400
    assert c.post("/api/compare", json={"a": "MKT", "b": ""}).status_code == 400
    assert c.post("/api/compare", json={"a": "MKT", "b": "MKT", "mode": "x"}).status_code == 400


def test_similar_and_classify_need_dataset(c):
    assert "No dataset" in c.post("/api/similar", json={"query": HBA_SEQ}).get_json()["error"]
    assert "No trained model" in c.post("/api/classify", json={"text": HBA_SEQ}).get_json()["error"]


def test_similar_and_classify_with_model(trained_results):
    app = create_app(FakeUniProtClient(), results_dir=trained_results)
    c = app.test_client()
    d = c.post("/api/similar", json={"query": HBA_SEQ, "top": 5}).get_json()
    assert len(d["hits"]) == 5 and d["dataset_size"] == 93
    d = c.post("/api/classify", json={"text": f">x\n{BACR_SEQ}\n>y\n{HBA_SEQ}"}).get_json()
    assert len(d["predictions"]) == 2 and "confidence" in d["predictions"][0]
    d = c.post("/api/classify", json={"text": "P69905"}).get_json()
    assert d["predictions"][0]["name"].startswith("P69905") and d["backend"] == "handcrafted"
    r = c.post("/api/classify", data={"file": (io.BytesIO(f">up1\n{HBA_SEQ}\n>up2\n{BACR_SEQ}".encode()), "x.fa"),
                                      "backend": "handcrafted"}, content_type="multipart/form-data")
    assert r.status_code == 200 and [p["name"] for p in r.get_json()["predictions"]] == ["up1", "up2"]
    assert c.post("/api/classify", json={"text": HBA_SEQ, "backend": "esm2"}).status_code == 400  # not trained
    s = c.get("/api/ml/status").get_json()
    assert [m["backend"] for m in s["models"]] == ["handcrafted"] and s["models"][0]["usable"]
    m = s["metrics"]["handcrafted"]
    assert m["classes"] and "confusion_matrix_handcrafted.png" in m["figures"]
    assert c.get("/results/confusion_matrix_handcrafted.png").status_code == 200
    assert c.get("/results/report_handcrafted.md").status_code == 200
    assert c.get("/results/../app.py").status_code == 404
    assert c.get("/results/model_handcrafted.joblib").status_code == 404  # never serve the model
    assert c.get("/results/jobs.json").status_code == 404


def test_structure_endpoint(c):
    d = c.get("/api/structure/P69905").get_json()
    assert d["pdb"].startswith("ATOM") and d["mean_plddt"] == 91.5
    r = c.get("/api/structure/Q99999")
    assert r.status_code == 404 and "no predicted structure" in r.get_json()["error"]


def test_protein_features_in_detail(c):
    d = c.get("/api/protein/P69905").get_json()
    assert {f["type"] for f in d["features"]} == {"Signal", "Transmembrane", "Binding site"}


def _wait(c, jid, timeout=120):
    t = time.time()
    while time.time() - t < timeout:
        j = c.get(f"/api/jobs/{jid}").get_json()
        if j["state"] not in ("running", "cancelling"):
            return j
        time.sleep(0.1)
    raise AssertionError("job timed out")


def test_dataset_and_train_jobs(c):
    assert c.post("/api/jobs/train", json={}).status_code == 400  # no dataset yet
    assert c.post("/api/jobs/dataset", json={"cutoff": "yesterday"}).status_code == 400
    assert c.post("/api/jobs/dataset", json={"classes": ["kinase"], "include_other": False}).status_code == 400
    r = c.post("/api/jobs/dataset", json={"classes": ["kinase", "protease", "toxin"], "per_class": 20, "holdout": 6})
    assert r.status_code == 202
    jid = r.get_json()["job_id"]
    j = _wait(c, jid)
    assert j["state"] == "done", j
    s = c.get("/api/ml/status").get_json()
    assert s["dataset"]["counts"]["kinase"]["train"] == 20 and s["dataset"]["counts"]["other"]["train"] == 20
    assert s["last_job"]["id"] == jid
    assert c.post("/api/jobs/train", json={"backend": "nope"}).status_code == 400
    r = c.post("/api/jobs/train", json={"trees": 50, "align": False})
    assert r.status_code == 202
    assert c.post("/api/jobs/dataset", json={}).status_code == 409  # only one job at a time
    j = _wait(c, r.get_json()["job_id"])
    assert j["state"] == "done", j
    assert "accuracy" in j["result"]["holdout"] and j["result"]["backend"] == "handcrafted"
    assert c.get("/api/jobs/doesnotexist").status_code == 404
    assert c.post(f"/api/jobs/{j['id']}/cancel").status_code == 400  # already finished


def test_cancel_running_job(c):
    r = c.post("/api/jobs/dataset", json={"classes": ["kinase", "protease"], "per_class": 400, "holdout": 100})
    jid = r.get_json()["job_id"]
    c.post(f"/api/jobs/{jid}/cancel")
    j = _wait(c, jid)
    assert j["state"] in ("cancelled", "done")  # the fake client may finish before the cancel lands
    if j["state"] == "cancelled":
        assert "Nothing was overwritten" in j["message"]
    assert c.post("/api/jobs/nope/cancel").status_code == 404
