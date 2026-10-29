"""Offline tests for mutations, multiple alignment, explanations, BLAST client, projects and public mode."""
from __future__ import annotations

import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / "tests"))

import pytest  # noqa: E402

from fixtures import BACR_SEQ, HBA_SEQ, HBB_SEQ, FakeUniProtClient  # noqa: E402
from protein_analyzer import msa  # noqa: E402
from protein_analyzer.ebi import BlastClient, parse_blast_json  # noqa: E402
from protein_analyzer.mutations import missense_matrix, mutate, parse_mutation  # noqa: E402
from protein_analyzer.projects import ProjectStore, to_markdown  # noqa: E402
from protein_analyzer.uniprot import UniProtBadRequest  # noqa: E402

MB_SEQ = ("MGLSDGEWQLVLNVWGKVEADIPGHGQEVLIRLFKGHPETLEKFDKFKHLKSEDEMKASEDLKKHGATVLTALGGILKKKGHHEAEIKPLAQSHATKHKIPVKY"
          "LEFISECIIQVLQSKHPGDFGADAQGAMNKALELFRKDMASNYKELGFQG")


# ------------------------------------------------------------------ mutations
def test_parse_mutation():
    assert parse_mutation("e7v") == ("E", 7, "V")
    with pytest.raises(ValueError, match="like E7V"):
        parse_mutation("7V")


def test_mutate_sickle_cell():
    r = mutate(HBB_SEQ, 7, "V", expected="E")
    assert r["mutation"] == "E7V" and r["from_group"] == "negative" and r["to_group"] == "hydrophobic"
    assert r["blosum62"] == -2 and r["hydropathy_change"] == pytest.approx(4.2 - (-3.5))
    assert r["delta"]["isoelectric_point"] > 0  # losing a negative charge raises pI
    assert any("sickle" in n for n in r["notes"])
    assert r["context"]["before"][r["position"] - r["context"]["start"]] == "E"
    assert r["context"]["after"][r["position"] - r["context"]["start"]] == "V"


def test_mutate_errors():
    with pytest.raises(ValueError, match="not K"):
        mutate(HBB_SEQ, 7, "V", expected="K")
    with pytest.raises(ValueError, match="outside"):
        mutate(HBB_SEQ, 999, "V")
    with pytest.raises(ValueError, match="already"):
        mutate(HBB_SEQ, 7, "E")
    with pytest.raises(ValueError, match="standard"):
        mutate(HBB_SEQ, 7, "X")


def test_missense_matrix():
    m = missense_matrix({"M1A": 0.4, "V2L": 0.9, "Q2L": 0.1, "bad": 1.0}, "MVH")
    assert m["matrix"][0][0] == 0.4 and m["matrix"][1][9] == 0.9  # L is index 9
    assert m["covered"] == 2  # Q2L ignored: position 2 is V, not Q


# ------------------------------------------------------------------ multiple alignment + tree
def test_msa_globins_tree_and_alignment():
    r = msa.build([("HBA", HBA_SEQ), ("HBB", HBB_SEQ), ("MB", MB_SEQ), ("HBB2", HBB_SEQ[:-3])])
    assert len(set(map(len, r["rows"]))) == 1  # all rows same length
    for row, s in zip(r["rows"], [HBA_SEQ, HBB_SEQ, MB_SEQ, HBB_SEQ[:-3]]):
        assert row.replace("-", "") == s  # no residue lost
    assert r["identity_matrix"][1][3] == 100.0
    nwk = r["newick"]
    assert nwk.endswith(";") and all(n in nwk for n in ("HBA", "HBB", "MB", "HBB2"))

    # HBB and its truncated copy must be siblings in the tree
    def leaves(n):
        return [n["name"]] if not n["children"] else sum((leaves(c) for c in n["children"]), [])

    def find_parent_of(n, name):
        for c in n["children"]:
            if not c["children"] and c["name"] == name:
                return n
            p = find_parent_of(c, name)
            if p:
                return p
    assert sorted(leaves(find_parent_of(r["tree"], "HBB"))) == ["HBB", "HBB2"]


def test_msa_limits():
    with pytest.raises(ValueError, match="at least 3"):
        msa.build([("a", HBA_SEQ), ("b", HBB_SEQ)])
    with pytest.raises(ValueError, match="up to"):
        msa.build([("a", "A" * 1501), ("b", HBB_SEQ), ("c", HBA_SEQ)])


def test_neighbor_joining_known_tree():
    import numpy as np
    # Additive distances from tree ((A:1,B:2):1,(C:3,D:1)) -> NJ must recover the topology AB|CD
    d = np.array([[0, 3, 5, 3], [3, 0, 6, 4], [5, 6, 0, 4], [3, 4, 4, 0]], dtype=float)
    t = msa.neighbor_joining(d, ["A", "B", "C", "D"])
    nwk = msa.to_newick(t)
    assert ("(A:" in nwk and "B:" in nwk) and ("C:" in nwk and "D:" in nwk)
    sides = [sorted(_leaves(c)) for c in t["children"]]
    assert ["A", "B"] in sides or ["C", "D"] in sides


def _leaves(n):
    return [n["name"]] if not n["children"] else sum((_leaves(c) for c in n["children"]), [])


# ------------------------------------------------------------------ explanations
@pytest.fixture(scope="module")
def trained(tmp_path_factory):
    from make_synthetic import make_dataset
    from protein_analyzer import ml
    d = tmp_path_factory.mktemp("res")
    make_dataset(d / "dataset.csv", per_class=25, holdout=6)
    ml.train_and_evaluate(d / "dataset.csv", d, n_estimators=30, align_identity=False)
    return d


def test_explain_handcrafted(trained):
    from protein_analyzer import ml
    from protein_analyzer.explain import explain
    x = explain(BACR_SEQ, ml.model_path(trained, "handcrafted"), "bR")
    assert x["backend"] == "handcrafted" and 0 <= x["probability"] <= 1
    regs = x["region_scan"]["regions"]
    assert regs[0]["start"] == 1 and regs[-1]["end"] == len(BACR_SEQ) and len(regs) <= 40
    f = x["features"]
    assert f["features"] and all(r["sentence"] for r in f["features"])
    # tree-interpreter identity: bias + sum(all contributions) == model probability
    bundle = ml.load_model(ml.model_path(trained, "handcrafted"))
    from protein_analyzer.explain import forest_contributions
    xv = ml._features([BACR_SEQ], "handcrafted", trained)[0]
    full = forest_contributions(bundle["model"], xv, list(bundle["model"].classes_).index(x["prediction"]),
                                bundle["train_mean"], top=434)
    assert f["bias"] + sum(r["contribution"] for r in full["features"]) == pytest.approx(x["probability"], abs=0.01)


def test_feature_words():
    from protein_analyzer.explain import feature_words
    assert feature_words("aac_C") == "cysteine content"
    assert feature_words("dpc_KR") == "'KR' residue pairs"
    assert "positively charged" in feature_words("group_positive")
    assert feature_words("pI") == "isoelectric point"


# ------------------------------------------------------------------ EBI BLAST
class _R:
    def __init__(self, status, text="", js=None):
        self.status_code, self.text, self._js = status, text, js

    def json(self):
        return self._js


class _S:
    def __init__(self, rs):
        self.rs, self.calls, self.headers = list(rs), [], {}

    def request(self, method, url, timeout=None, **kw):
        self.calls.append((method, url, kw))
        return self.rs.pop(0)


BLAST_JSON = {"program": "blastp", "query_len": 142, "dbs": [{"name": "uniprotkb_swissprot"}], "hits": [
    {"hit_acc": "P69905", "hit_id": "SP:HBA_HUMAN", "hit_uni_de": "Hemoglobin subunit alpha", "hit_uni_os": "Homo sapiens",
     "hit_len": 142, "hit_hsps": [{"hsp_bit_score": 290, "hsp_identity": 100.0, "hsp_positive": 100.0, "hsp_expect": 1e-100,
                                    "hsp_align_len": 142, "hsp_query_from": 1, "hsp_query_to": 142, "hsp_hit_from": 1, "hsp_hit_to": 142}]},
    {"hit_acc": "P01942", "hit_uni_de": "Hemoglobin subunit alpha", "hit_uni_os": "Mus musculus", "hit_len": 142,
     "hit_hsps": [{"hsp_bit_score": 10, "hsp_identity": 20.0}, {"hsp_bit_score": 250, "hsp_identity": 86.0, "hsp_expect": 1e-80}]}]}


def test_blast_client_flow():
    s = _S([_R(200, "ncbiblast-R20261024-123-1-p1m\n"), _R(200, "RUNNING"), _R(200, "FINISHED"), _R(200, js=BLAST_JSON)])
    b = BlastClient(session=s, min_interval=0)
    jid = b.submit(HBA_SEQ, "student@example.org")
    assert jid == "ncbiblast-R20261024-123-1-p1m"
    sent = s.calls[0][2]["data"]
    assert sent["email"] == "student@example.org" and sent["program"] == "blastp" and sent["database"] == "uniprotkb_swissprot"
    assert b.status(jid) == "RUNNING" and b.status(jid) == "FINISHED"
    res = b.result(jid)
    assert [h["accession"] for h in res["hits"]] == ["P69905", "P01942"]
    assert res["hits"][1]["identity"] == 86.0  # best HSP chosen by bit score
    assert b.result(jid) is res and len(s.calls) == 4  # cached


def test_blast_client_validation():
    b = BlastClient(session=_S([_R(400, "<error><description>Invalid sequence</description></error>")]), min_interval=0)
    with pytest.raises(UniProtBadRequest, match="valid email"):
        b.submit(HBA_SEQ, "nope")
    with pytest.raises(UniProtBadRequest, match="Unknown database"):
        b.submit(HBA_SEQ, "a@b.co", database="evil")
    with pytest.raises(UniProtBadRequest, match="Invalid sequence"):
        b.submit(HBA_SEQ, "a@b.co")
    with pytest.raises(UniProtBadRequest, match="job id"):
        b.status("../../etc")


def test_parse_blast_json_empty():
    assert parse_blast_json({"hits": []})["hits"] == []


# ------------------------------------------------------------------ projects
def test_project_store(tmp_path):
    st = ProjectStore(tmp_path / "p.json")
    p = st.create("Haemoglobin")
    with pytest.raises(ValueError):
        st.create("  ")
    it = st.add_item(p["id"], "mutation", "HBB E7V", {"mutation": "E7V", "from_name": "glutamate", "to_name": "valine",
                                                       "blosum62": -2, "alphamissense": 0.22, "alphamissense_class": "LBen"})
    with pytest.raises(ValueError, match="Unknown item"):
        st.add_item(p["id"], "virus", "x", {})
    st.update_item(p["id"], it["id"], "note!")
    st.update(p["id"], notes="Project notes")
    assert st.list()[0]["n_items"] == 1
    md = to_markdown(st.get(p["id"]))
    assert "# Haemoglobin" in md and "E7V" in md and "AlphaMissense 0.22" in md and "> note!" in md
    st.delete_item(p["id"], it["id"])
    with pytest.raises(KeyError):
        st.delete_item(p["id"], it["id"])
    st.delete(p["id"])
    assert st.list() == []
    (tmp_path / "p.json").write_text("{corrupt")
    assert ProjectStore(tmp_path / "p.json").list() == []  # corrupt file is set aside, not fatal


# ------------------------------------------------------------------ API
@pytest.fixture
def client(tmp_path, trained):
    from app import create_app

    class FakeBlast:
        def submit(self, seq, email, database="uniprotkb_swissprot", max_hits=50):
            if "@" not in email:
                raise UniProtBadRequest("EBI requires a valid email address for every search.")
            return "ncbiblast-TEST"

        def status(self, jid):
            return "FINISHED"

        def result(self, jid):
            return parse_blast_json(BLAST_JSON)

    app = create_app(FakeUniProtClient(), results_dir=trained, blast_client=FakeBlast(), data_dir=tmp_path / "ud")
    app.testing = True
    return app.test_client()


def test_new_pages_render(client):
    for url in ["/mutate", "/blast", "/align", "/learn", "/projects", "/static/lessons.js"]:
        assert client.get(url).status_code == 200, url


def test_api_mutate_with_alphamissense_and_variants(client):
    r = client.post("/api/mutate", json={"accession": "P68871", "mutation": "E7V"}).get_json()
    assert r["mutation"] == "E7V" and r["alphamissense"] is not None and r["alphamissense_class"]
    assert [v["description"] for v in r["known_variants"] if v["to"] == "V"] == ["in SKCA; Hb S"]
    r = client.post("/api/mutate", json={"accession": "P68871", "position": 7, "to": "K"}).get_json()
    assert r["mutation"] == "E7K"
    r = client.post("/api/mutate", json={"sequence": HBB_SEQ, "mutation": "E7V"}).get_json()
    assert r["accession"] is None and r["alphamissense"] is None
    assert client.post("/api/mutate", json={"accession": "P68871", "mutation": "K7V"}).status_code == 400
    assert client.post("/api/mutate", json={"accession": "P68871", "mutation": "zzz"}).status_code == 400


def test_api_alphamissense_variants_domains(client):
    am = client.get("/api/protein/P68871/alphamissense").get_json()
    assert len(am["matrix"]) == len(HBB_SEQ) and am["covered"] == 20 * 19
    assert client.get("/api/protein/P02945/alphamissense").status_code == 404
    assert client.get("/api/protein/P68871/variants").get_json()["variants"][0]["to"] == "V"
    assert client.get("/api/protein/P69905/domains").get_json()["domains"][0]["accession"] == "PF00042"


def test_api_msa(client):
    r = client.post("/api/msa", json={"text": "P69905 P68871 P02945"})
    d = r.get_json()
    assert r.status_code == 200 and len(d["rows"]) == 3 and d["newick"].endswith(";")
    assert client.post("/api/msa", json={"text": "P69905 hello P02945"}).status_code == 400
    d = client.post("/api/msa", json={"text": f">a\n{HBA_SEQ}\n>b\n{HBB_SEQ}\n>c\n{MB_SEQ}"}).get_json()
    assert d["names"] == ["a", "b", "c"]
    assert client.post("/api/msa", json={"text": f">a\n{HBA_SEQ}\n>b\n{HBB_SEQ}"}).status_code == 400


def test_api_blast(client):
    assert client.post("/api/blast", json={"query": "P69905", "email": "x"}).status_code == 400
    assert client.post("/api/blast", json={"query": "MKT", "email": "a@b.co"}).status_code == 400  # too short
    r = client.post("/api/blast", json={"query": "P69905", "email": "a@b.co"})
    assert r.status_code == 202 and r.get_json()["job_id"] == "ncbiblast-TEST"
    d = client.get("/api/blast/ncbiblast-TEST").get_json()
    assert d["status"] == "FINISHED" and d["result"]["hits"][0]["accession"] == "P69905"


def test_api_explain(client):
    d = client.post("/api/explain", json={"text": "P02945"}).get_json()
    assert d["prediction"] and d["region_scan"]["regions"] and d["features"]["features"]
    assert client.post("/api/explain", json={"text": "P02945", "backend": "esm2"}).status_code == 400


def test_api_projects_flow(client):
    p = client.post("/api/projects", json={"name": "Course"}).get_json()
    assert client.post("/api/projects", json={"name": ""}).status_code == 400
    it = client.post(f"/api/projects/{p['id']}/items", json={"kind": "comparison", "title": "a vs b",
                     "data": {"name_a": "a", "name_b": "b", "mode": "global", "identity": 42, "similarity": 59,
                              "score": 282, "gaps": 9, "length": 3, "aligned_a": "MV-", "aligned_b": "MVH", "match_line": "|| "}}).get_json()
    assert client.patch(f"/api/projects/{p['id']}/items/{it['id']}", json={"note": "related"}).status_code == 200
    assert client.patch(f"/api/projects/{p['id']}", json={"notes": "my notes"}).get_json()["notes"] == "my notes"
    full = client.get(f"/api/projects/{p['id']}").get_json()
    assert full["items"][0]["note"] == "related"
    md = client.get(f"/api/projects/{p['id']}/export.md")
    assert md.status_code == 200 and b"identity 42%" in md.data and "attachment" in md.headers["Content-Disposition"]
    assert client.get(f"/api/projects/{p['id']}/export.json").get_json()["name"] == "Course"
    rep = client.get(f"/projects/{p['id']}/report")
    assert rep.status_code == 200 and b"42% identity" in rep.data
    assert client.get(f"/api/projects/{p['id']}/export.exe").status_code == 404
    assert client.delete(f"/api/projects/{p['id']}/items/{it['id']}").status_code == 200
    assert client.delete(f"/api/projects/{p['id']}").status_code == 200
    assert client.get(f"/api/projects/{p['id']}").status_code == 404


def test_public_mode(client, monkeypatch):
    from protein_analyzer import config
    monkeypatch.setattr(config, "PUBLIC", True)
    monkeypatch.setattr(config, "PUBLIC_RATE_PER_MIN", 3)
    assert client.post("/api/jobs/train", json={}).status_code == 403
    assert client.post("/api/jobs/dataset", json={}).status_code == 403
    assert client.get("/api/projects").status_code == 400
    assert b"Projects</a>" not in client.get("/analyze").data  # nav link hidden
    codes = [client.get("/api/search?q=x").status_code for _ in range(5)]
    assert codes[:3] == [200, 200, 200] and codes[3:] == [429, 429]
    assert client.post("/api/analyze", json={"text": HBA_SEQ}).status_code == 200  # local-only work isn't limited
