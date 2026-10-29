"""Research-readiness tests: ExPASy-compatible pI, DSSP secondary structure, provenance stamping."""
from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / "tests"))

import numpy as np  # noqa: E402
import pytest  # noqa: E402

from fixtures import HBA_SEQ, FakeUniProtClient  # noqa: E402
from protein_analyzer import provenance  # noqa: E402
from protein_analyzer.analyzer import PI_METHOD, isoelectric_point, molecular_weight  # noqa: E402
from protein_analyzer.dssp import assign, assign_pdb, read_backbone  # noqa: E402
from protein_analyzer.structure import SS_METHOD_DSSP, structure_features  # noqa: E402

DATA = ROOT / "tests" / "data"
HBA_PDB = (DATA / "AF-P69905-F1-model_v6.pdb").read_text()
GFP_PDB = (DATA / "AF-P42212-F1-model_v6.pdb").read_text()

# Reference values from the live ExPASy Compute pI/Mw tool (web.expasy.org/compute_pi), 2026-10-24.
EXPASY_PI = {
    HBA_SEQ: 8.72,
    "DDDDDDDDEEEEEEEEGGGGAD": 2.90,   # acidic C-terminal D: ExPASy applies no C-terminal correction
    "EEEEEEEEEEEEGGGGGGGGE": 3.20,
    "KKKKKKKKGGGGGGGGE": 10.48,
    "PGGGGGGGGHGGG": 7.17,           # residue-specific N-terminal pK (Pro 8.36)
    "AGGGGGGGGHGGG": 6.79,
    "SGGGGGGGGHGGG": 6.46,
    "TCYHGGGGGGGGGHGGGG": 6.61,
}


@pytest.mark.parametrize("seq,expected", list(EXPASY_PI.items()))
def test_pi_matches_expasy(seq, expected):
    assert round(isoelectric_point(seq), 2) == pytest.approx(expected, abs=0.011)


def test_mw_matches_expasy():
    assert molecular_weight(HBA_SEQ) == pytest.approx(15257.55, abs=0.02)  # ExPASy Compute pI/Mw
    assert "Bjellqvist" in PI_METHOD


# ------------------------------------------------------------------ DSSP
def test_dssp_real_models():
    hba, _ = assign_pdb(HBA_PDB)
    assert len(hba) == 142 and hba.count("E") == 0 and 95 <= hba.count("H") <= 115  # all-alpha globin
    gfp, _ = assign_pdb(GFP_PDB)
    assert gfp.count("E") >= 105  # 11-stranded beta barrel


def test_dssp_port_matches_pydssp_exactly():
    pydssp = pytest.importorskip("pydssp.pydssp_numpy")
    read = pytest.importorskip("pydssp").read_pdbtext
    for pdb in (HBA_PDB, GFP_PDB):
        ref = "".join("CHE"[i] for i in pydssp.assign(np.asarray(read(pdb), float)).argmax(-1))
        assert assign(read_backbone(pdb)[0]) == ref  # same settings (no proline mask) -> identical


def test_structure_features_use_dssp_with_fallback():
    f = structure_features(HBA_PDB)
    assert f["method"] == SS_METHOD_DSSP and f["strand_percent"] == 0.0 and f["helix_percent"] > 70
    ca_only = "\n".join(l for l in HBA_PDB.splitlines() if l[12:16] == " CA ") + "\nEND\n"
    g = structure_features(ca_only)
    assert "C-alpha" in g["method"] and g["helix_percent"] > 60  # graceful fallback without backbone


def test_dssp_short_and_long_inputs():
    assert assign(np.zeros((3, 4, 3))) == "CCC"
    coord, seq, nums = read_backbone(HBA_PDB)
    assert coord.shape == (142, 4, 3) and seq.startswith("MVLSP") and nums[0] == 1


# ------------------------------------------------------------------ provenance
def test_stamp_contents():
    s = provenance.stamp("analysis", {"release": "2026_03", "date": "02-September-2026"})
    assert s["app"].startswith("Protein Analyzer ") and s["uniprot_release"] == "2026_03"
    assert any("Bjellqvist" in m for m in s["methods"]) and s["generated"]
    m = provenance.stamp("mutation")
    assert any("AlphaMissense (CC BY 4.0)" in x for x in m["sources"])
    for k in provenance.KIND_SOURCES:
        assert all(s in provenance.SOURCES for s in provenance.KIND_SOURCES[k])
        assert all(x in provenance.METHODS for x in provenance.KIND_METHODS[k])


def test_release_captured_from_headers():
    from protein_analyzer.uniprot import UniProtClient

    class R:
        status_code, text = 200, "{}"
        headers = {"X-UniProt-Release": "2026_03", "X-UniProt-Release-Date": "02-September-2026", "X-Total-Results": "1"}

        def json(self):
            return {"results": []}

    class S:
        headers = {}

        def get(self, *a, **k):
            return R()

    c = UniProtClient(cache_dir=None, session=S(), min_interval=0)
    assert c.uniprot_release() == {"release": "2026_03", "date": "02-September-2026"}
    c2 = UniProtClient(cache_dir=None, session=S(), min_interval=0)
    c2.count("x")
    assert c2.release["release"] == "2026_03"


@pytest.fixture
def client(tmp_path):
    from app import create_app

    class Fake(FakeUniProtClient):
        def uniprot_release(self):
            return {"release": "2026_03", "date": "02-September-2026"}

    app = create_app(Fake(), results_dir=tmp_path / "r", data_dir=tmp_path / "ud")
    app.testing = True
    return app.test_client()


def test_api_provenance_and_methods_page(client):
    p = client.get("/api/provenance?kind=analysis").get_json()
    assert p["uniprot_release"] == "2026_03" and any("ExPASy" in m for m in p["methods"])
    assert client.get("/api/provenance?kind=nope").status_code == 400
    page = client.get("/methods").data.decode()
    assert "How to report your results" in page and "2026_03" in page and "CC BY 4.0" in page and "Kabsch" in page
    assert 'href="/methods"' in client.get("/analyze").data.decode()  # footer link on every page


def test_project_items_and_exports_carry_provenance(client):
    pid = client.post("/api/projects", json={"name": "Lab report"}).get_json()["id"]
    it = client.post(f"/api/projects/{pid}/items", json={"kind": "mutation", "title": "HBB E7V",
                     "data": {"mutation": "E7V", "from_name": "glutamate", "to_name": "valine", "blosum62": -2}}).get_json()
    assert it["provenance"]["uniprot_release"] == "2026_03" and it["provenance"]["methods"]
    md = client.get(f"/api/projects/{pid}/export.md").data.decode()
    assert "## Methods & data sources" in md and "UniProtKB release(s) 2026_03" in md and "Cheng J" in md
    rep = client.get(f"/projects/{pid}/report").data.decode()
    assert "Methods &amp; data sources" in rep and "AlphaMissense" in rep and "2026_03" in rep
    js = client.get(f"/api/projects/{pid}/export.json").get_json()
    assert js["items"][0]["provenance"]["app"].startswith("Protein Analyzer")


# ------------------------------------------------------------------ v1.5 QOL: home page, backup / import
def test_home_page_and_nav(client):
    html = client.get("/").data.decode()
    assert "Explore any protein" in html and 'id="homeq"' in html and 'data-theme-set="dark"' in html
    assert client.get("/analyze").status_code == 200 and 'id="menubtn"' in client.get("/analyze").data.decode()


def test_projects_backup_and_import_roundtrip(client):
    pid = client.post("/api/projects", json={"name": "Course"}).get_json()["id"]
    client.post(f"/api/projects/{pid}/items", json={"kind": "note", "title": "first note", "data": {}, "note": "hello"})
    r = client.get("/api/projects/backup")
    assert r.status_code == 200 and "attachment" in r.headers["Content-Disposition"]
    backup = r.get_json()
    assert backup["format"] == "protein-analyzer-projects" and backup["projects"][0]["name"] == "Course"
    # import the full backup -> a renamed copy; import a single exported project -> another copy
    created = client.post("/api/projects/import", json=backup).get_json()["created"]
    assert created == [{"id": created[0]["id"], "name": "Course (2)", "n_items": 1}]
    single = client.get(f"/api/projects/{pid}/export.json").get_json()
    created2 = client.post("/api/projects/import", json=single).get_json()["created"]
    assert created2[0]["name"] == "Course (3)"
    copy = client.get(f"/api/projects/{created[0]['id']}").get_json()
    assert copy["items"][0]["note"] == "hello" and copy["items"][0]["id"] != backup["projects"][0]["items"][0]["id"]
    assert len(client.get("/api/projects").get_json()["projects"]) == 3


def test_projects_import_rejects_bad_files(client):
    assert client.post("/api/projects/import", json={"hello": 1}).status_code == 400
    assert client.post("/api/projects/import", json={"projects": ["x"]}).status_code == 400
    r = client.post("/api/projects/import", json={"name": "ok", "items": [{"kind": "virus", "data": {}}, {"kind": "note", "title": "t"}]})
    assert r.status_code == 201 and r.get_json()["created"][0]["n_items"] == 1  # unknown kinds are skipped
    import io
    r = client.post("/api/projects/import", data={"file": (io.BytesIO(b"not json"), "x.json")}, content_type="multipart/form-data")
    assert r.status_code == 400 and "valid JSON" in r.get_json()["error"]


def test_page_scripts_are_valid_javascript(client, tmp_path):
    """Every inline <script> on every page must parse (catches broken string escapes in templates)."""
    import re
    import shutil
    import subprocess
    node = shutil.which("node")
    if not node:
        pytest.skip("node not installed")
    pages = ["/", "/analyze", "/search", "/blast", "/compare", "/align", "/structures", "/mutate", "/lab",
             "/classify", "/ml", "/learn", "/projects", "/methods", "/protein/P69905"]
    for page in pages:
        html = client.get(page).data.decode()
        for i, js in enumerate(re.findall(r"<script>(.*?)</script>", html, re.S)):
            f = tmp_path / "s.js"
            f.write_text(js, encoding="utf-8")
            r = subprocess.run([node, "--check", str(f)], capture_output=True, text=True)
            assert r.returncode == 0, f"{page} script {i}: {r.stderr[:300]}"
    r = subprocess.run([node, "--check", str(ROOT / "static" / "app.js")], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr[:300]
