"""Offline tests for lab tools, structure features/superposition, and the interactions/pathways endpoints."""
from __future__ import annotations

import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / "tests"))

import numpy as np  # noqa: E402
import pytest  # noqa: E402

from fixtures import HBA_SEQ, HBB_SEQ, FakeUniProtClient, ideal_pdb  # noqa: E402
from protein_analyzer import labtools as lt  # noqa: E402
from protein_analyzer.structure import parse_ca, secondary_structure, structure_features, superpose  # noqa: E402

UBIQUITIN = "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG"


# ------------------------------------------------------------------ DNA
def test_translate_and_frames():
    dna = "ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG"
    assert lt.translate(dna) == "MAIVMGR*KGAR*"
    frames = lt.six_frames(dna)
    assert [f["strand"] + str(f["frame"]) for f in frames] == ["+1", "+2", "+3", "-1", "-2", "-3"]
    assert lt.reverse_complement("ATGCN") == "NGCAT"
    assert lt.clean_dna(">x\naugc\nuu ") == "ATGCTT"  # RNA + FASTA header accepted


def test_orfs_both_strands_and_incomplete():
    gene = "ATG" + "GCT" * 60 + "TAA"
    dna = "CC" + gene + "GG"
    orfs = lt.find_orfs(dna, min_aa=50)
    assert orfs[0]["strand"] == "+" and orfs[0]["start"] == 3 and orfs[0]["end"] == 3 + len(gene) - 1
    assert orfs[0]["length_aa"] == 61 and orfs[0]["complete"] and orfs[0]["protein"] == "M" + "A" * 60
    rc = lt.find_orfs(lt.reverse_complement(dna), min_aa=50)
    assert rc[0]["strand"] == "-" and (rc[0]["start"], rc[0]["end"]) == (3, len(gene) + 2)  # same gene, other strand
    open_end = lt.find_orfs("ATG" + "GCT" * 60, min_aa=50)
    assert open_end and not open_end[0]["complete"]


def test_genetic_codes_and_alt_starts():
    assert lt.translate("TGAAGA", code=1) == "*R"  # standard: TGA stop, AGA arginine
    assert lt.translate("TGAAGA", code=2) == "W*"  # vertebrate mitochondria: TGA = Trp, AGA = stop
    dna = "GTG" + "GCT" * 60 + "TAA"
    assert lt.find_orfs(dna, min_aa=50) == []
    alt = lt.find_orfs(dna, min_aa=50, code=11, alt_starts=True)
    assert alt and alt[0]["start_codon"] == "GTG" and alt[0]["protein"].startswith("M")  # GTG start read as Met


def test_dna_validation():
    with pytest.raises(ValueError, match="protein"):
        lt.clean_dna("MKTLLEEQ")
    with pytest.raises(ValueError, match="Paste"):
        lt.clean_dna("  ")
    assert lt.dna_stats("GGCCAATT")["gc_percent"] == 50.0


# ------------------------------------------------------------------ digest
def test_peptide_masses_reference():
    assert lt.peptide_mass("DRVYIHPF") + lt.PROTON == pytest.approx(1046.5418, abs=1e-3)  # angiotensin II [M+H]+
    assert lt.peptide_mass("C", cam=True) - lt.peptide_mass("C") == pytest.approx(57.021464)


def test_digest_rules():
    s = "AKPAKAARGGDEEK"
    assert lt.cleavage_sites(s, "trypsin") == [5, 8]  # after K5 and R8, but not after K2 (followed by P)
    assert lt.cleavage_sites(s, "trypsin_p") == [2, 5, 8]
    assert lt.cleavage_sites(s, "aspn") == [10]  # before D
    assert lt.cleavage_sites(s, "lysc") == [2, 5]
    d = lt.digest(s, "trypsin", missed=1, min_len=1, min_mass=0, max_mass=1e9)
    seqs = [p["sequence"] for p in d["peptides"]]
    assert seqs == ["AKPAK", "AKPAKAAR", "AAR", "AARGGDEEK", "GGDEEK"]
    assert d["coverage"] == 100.0 and d["sites"] == 2
    p = d["peptides"][0]
    assert p["mz2"] == pytest.approx((p["mono"] + 2 * lt.PROTON) / 2, abs=1e-4)
    with pytest.raises(ValueError):
        lt.digest(s, "nope")


# ------------------------------------------------------------------ A280
def test_extinction_and_concentration():
    e = lt.extinction(UBIQUITIN)
    assert e["epsilon_reduced"] == 1490 and e["epsilon_cystines"] == 1490  # ExPASy ProtParam value
    e = lt.extinction("WWYCCC")
    assert e["epsilon_reduced"] == 2 * 5500 + 1490 and e["epsilon_cystines"] == 2 * 5500 + 1490 + 125
    assert lt.extinction("AAAK")["warning"]
    c = lt.concentration(0.5, 50000, 25000, path_cm=1, dilution=2)
    assert c["micromolar"] == pytest.approx(20.0) and c["mg_per_ml"] == pytest.approx(0.5)
    assert lt.concentration(2.0, 50000, 25000)["warning"]
    with pytest.raises(ValueError):
        lt.concentration(0.5, 0, 25000)


# ------------------------------------------------------------------ structure
def test_secondary_structure_ideal_geometry():
    helix = parse_ca(ideal_pdb("A" * 30, "helix"))
    assert secondary_structure(helix["xyz"]).count("H") >= 26
    hp = parse_ca(ideal_pdb("V" * 20, "hairpin"))
    ss = secondary_structure(hp["xyz"])
    assert ss.count("E") >= 10 and "H" not in ss
    lonely = np.array([[3.3 * i, 0.9 * (-1) ** i, 0] for i in range(10)], float)  # extended, but no partner strand
    assert "E" not in secondary_structure(lonely)
    low = parse_ca(ideal_pdb("A" * 30, "helix", plddt=30))
    f = structure_features(ideal_pdb("A" * 30, "helix", plddt=30))
    assert secondary_structure(low["xyz"], low["plddt"]) == "C" * 30 and f["disordered_percent"] == 100.0


def test_structure_features_segments():
    f = structure_features(ideal_pdb("A" * 40, "helix"))
    assert f["helix_percent"] > 85 and f["helices"][0][0] >= 1 and f["first_residue"] == 1


def test_superpose_recovers_rotation():
    seq = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQ"
    a = ideal_pdb(seq, "helix")
    th = math.radians(40)
    R = np.array([[math.cos(th), -math.sin(th), 0], [math.sin(th), math.cos(th), 0], [0, 0, 1]])
    lines = []
    for line in a.splitlines():
        if line.startswith("ATOM"):
            x = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])]) @ R.T + [5, -3, 12]
            line = f"{line[:30]}{x[0]:8.3f}{x[1]:8.3f}{x[2]:8.3f}{line[54:]}"
        lines.append(line)
    r = superpose(a, "\n".join(lines), seq, seq)
    assert r["rmsd_all"] < 0.01 and r["tm_score_a"] > 0.99 and r["core_pairs"] == len(seq)
    moved = parse_ca(r["pdb_b_transformed"])["xyz"]
    assert np.allclose(moved, parse_ca(a)["xyz"], atol=0.01)


def test_superpose_unrelated_shapes_score_low():
    r = superpose(ideal_pdb("A" * 40, "helix"), ideal_pdb("V" * 40, "hairpin"), "A" * 40, "A" * 40)
    assert r["tm_score_a"] < 0.5


# ------------------------------------------------------------------ API
@pytest.fixture
def client(tmp_path):
    from app import create_app
    app = create_app(FakeUniProtClient(), results_dir=tmp_path / "r", data_dir=tmp_path / "ud")
    app.testing = True
    return app.test_client()


def test_pages(client):
    for url in ["/lab", "/structures"]:
        assert client.get(url).status_code == 200


def test_api_lab(client):
    r = client.post("/api/lab/translate", json={"dna": "ATG" + "GCT" * 60 + "TAA", "min_aa": 50}).get_json()
    assert r["n_orfs"] == 1 and len(r["frames"]) == 6 and r["stats"]["length"] == 186
    assert client.post("/api/lab/translate", json={"dna": "MKTL"}).status_code == 400
    assert client.post("/api/lab/translate", json={"dna": "ATG", "code": 99}).status_code == 400
    d = client.post("/api/lab/digest", json={"text": "P69905", "enzyme": "trypsin", "missed": 1}).get_json()
    assert d["name"].startswith("P69905") and d["n_in_range"] > 5 and d["coverage"] > 50
    assert client.post("/api/lab/digest", json={"text": "P69905", "enzyme": "x"}).status_code == 400
    a = client.post("/api/lab/a280", json={"text": "P69905", "a280": 0.5, "dilution": 2}).get_json()
    assert a["epsilon_reduced"] == 1 * 5500 + 3 * 1490 and a["concentration"]["mg_per_ml"] > 0
    assert client.post("/api/lab/a280", json={"text": ""}).status_code == 400


def test_api_structure_interactions_pathways(client):
    f = client.get("/api/protein/P69905/structure-features").get_json()
    assert f["helix_percent"] > 80 and len(f["secondary_structure"]) == len(HBA_SEQ)
    i = client.get("/api/protein/P69905/interactions").get_json()
    assert [n["name"] for n in i["nodes"]] == ["HBA1", "HBB", "AHSP"] and len(i["edges"]) == 3
    assert client.get("/api/protein/Q99999/interactions").status_code == 404
    p = client.get("/api/protein/P69905/pathways").get_json()
    assert len(p["pathways"]) == 2 and p["pathways"][1]["disease"]
    assert client.get("/api/protein/P68871/pathways").get_json()["pathways"] == []


def test_api_superpose(client):
    r = client.post("/api/superpose", json={"a": "P69905", "b": "P68871"}).get_json()
    assert r["aligned_pairs"] > 100 and "pdb_b_transformed" in r and r["features_a"]["helix_percent"] > 80
    assert 0 <= r["tm_score_a"] <= 1
    assert client.post("/api/superpose", json={"a": "P69905", "b": "hello"}).status_code == 400
