"""Offline fixtures shaped like real UniProt REST JSON (trimmed from live responses, 2026-10)."""
from __future__ import annotations

HBA_SEQ = ("MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVADALTNAVAHVDDMPNALSALSDL"
           "HAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVLTSKYR")
HBB_SEQ = ("MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKKVLGAFSDGLAHLDNLKGTFATL"
           "SELHCDKLHVDPENFRLLGNVLVCVLAHHFGKEFTPPVQAAYQKVVAGVANALAHKYH")
BACR_SEQ = ("MLELLPTAVEGVSQAQITGRPEWIWLALGTALMGLGTLYFLVKGMGVSDPDAKKFYAITTLVPAIAFTMYLSMLLGYGLTMVPFGGEQNP"
            "IYWARYADWLFTTPLLLLDLALLVDADQGTILALVGADGIMIGTGLVGALTKVYSYRFVWWAISTAAMLYILYVLFFGFTSKAESMRPEV"
            "ASTFKVLRNVTVVLWSAYPVVWLIGSEGAGIVPLNIETLLFMVLDVSAKVGFGLILLRSRAIFGEAEAPEPSAGDGAAATSD")


def raw_entry(acc="P69905", seq=HBA_SEQ, name="Hemoglobin subunit alpha", gene="HBA1", ec=None,
              created="1986-07-21", reviewed=True):
    rec = {"fullName": {"value": name}}
    if ec:
        rec["ecNumbers"] = [{"value": e} for e in ec]
    return {
        "entryType": "UniProtKB reviewed (Swiss-Prot)" if reviewed else "UniProtKB unreviewed (TrEMBL)",
        "primaryAccession": acc,
        "uniProtkbId": f"{acc}_TEST",
        "entryAudit": {"firstPublicDate": created},
        "organism": {"scientificName": "Homo sapiens", "taxonId": 9606},
        "proteinDescription": {"recommendedName": rec},
        "genes": [{"geneName": {"value": gene}}],
        "comments": [{"commentType": "FUNCTION", "texts": [{"value": "Involved in oxygen transport."}]},
                     {"commentType": "SUBUNIT", "texts": [{"value": "Heterotetramer."}]}],
        "keywords": [{"id": "KW-0349", "category": "Ligand", "name": "Heme"}],
        "features": [
            {"type": "Signal", "location": {"start": {"value": 1}, "end": {"value": 10}}, "description": ""},
            {"type": "Transmembrane", "location": {"start": {"value": 20}, "end": {"value": 42}}, "description": "Helical"},
            {"type": "Binding site", "location": {"start": {"value": 88}, "end": {"value": 88}}, "description": "",
             "ligand": {"name": "heme b"}},
            {"type": "Region", "location": {"start": {"value": None, "modifier": "UNKNOWN"},
                                            "end": {"value": None, "modifier": "UNKNOWN"}}, "description": "Disordered"},
        ],
        "uniProtKBCrossReferences": [
            {"database": "GO", "id": "GO:0005833", "properties": [{"key": "GoTerm", "value": "C:hemoglobin complex"},
                                                                   {"key": "GoEvidenceType", "value": "IDA:UniProtKB"}]},
            {"database": "GO", "id": "GO:0019825", "properties": [{"key": "GoTerm", "value": "F:oxygen binding"}]},
            {"database": "GO", "id": "GO:0015671", "properties": [{"key": "GoTerm", "value": "P:oxygen transport"}]},
            {"database": "PDB", "id": "1A00", "properties": []},
        ],
        "sequence": {"value": seq, "length": len(seq)},
    }


def ideal_pdb(seq: str, kind: str = "helix", plddt: float = 90.0) -> str:
    """C-alpha-only PDB with ideal geometry: 'helix' (alpha helix) or 'hairpin' (two paired beta strands)."""
    import math
    three = {"A": "ALA", "R": "ARG", "N": "ASN", "D": "ASP", "C": "CYS", "Q": "GLN", "E": "GLU", "G": "GLY",
             "H": "HIS", "I": "ILE", "L": "LEU", "K": "LYS", "M": "MET", "F": "PHE", "P": "PRO", "S": "SER",
             "T": "THR", "W": "TRP", "Y": "TYR", "V": "VAL"}
    n, coords = len(seq), []
    for i in range(n):
        if kind == "helix":
            a = math.radians(100 * i)
            coords.append((2.3 * math.cos(a), 2.3 * math.sin(a), 1.5 * i))
        else:  # hairpin: first half along +x, second half back along -x, 4.8 Å away
            half = n // 2
            if i < half:
                coords.append((3.3 * i, 0.9 * (-1) ** i, 0.0))
            else:
                k = n - 1 - i
                coords.append((3.3 * k, 4.8 + 0.9 * (-1) ** k, 0.0))
    lines = [f"ATOM  {i + 1:5d}  CA  {three.get(c, 'GLY')} A{i + 1:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00{plddt:6.2f}           C"
             for i, (c, (x, y, z)) in enumerate(zip(seq, coords))]
    return "\n".join(lines) + "\nEND\n"


class FakeUniProtClient:
    """Stands in for UniProtClient in offline tests."""

    def __init__(self, entries=None, fail=None):
        self.entries = entries or {"P69905": raw_entry(), "P68871": raw_entry("P68871", HBB_SEQ, "Hemoglobin subunit beta", "HBB"),
                                   "P02945": raw_entry("P02945", BACR_SEQ, "Bacteriorhodopsin", "bop")}
        self.fail = fail
        self.calls = []

    def _check(self):
        if self.fail:
            raise self.fail

    def fetch_entry(self, acc):
        from protein_analyzer.uniprot import UniProtBadRequest, UniProtNotFound, is_accession, parse_entry
        self._check(); self.calls.append(("fetch", acc))
        acc = acc.strip().upper()
        if not is_accession(acc):
            raise UniProtBadRequest(f"'{acc}' is not a valid UniProt accession (e.g. P69905).")
        if acc not in self.entries:
            raise UniProtNotFound(f"No UniProt entry for accession {acc}.")
        return parse_entry(self.entries[acc])

    def search(self, query, size=25, cursor=None, fields=None):
        from protein_analyzer.uniprot import parse_entry
        self._check(); self.calls.append(("search", query, size, cursor))
        items = [parse_entry(e) for e in self.entries.values()]
        start = int(cursor or 0)
        page = items[start:start + size]
        nxt = str(start + size) if start + size < len(items) else None
        return {"results": page, "total": len(items), "next_cursor": nxt}

    def count(self, query):
        self._check()
        return 30

    def fetch_alphafold(self, acc):
        from protein_analyzer.uniprot import UniProtNotFound
        self._check()
        if acc.upper() not in self.entries:
            raise UniProtNotFound(f"AlphaFold DB has no predicted structure for {acc}.")
        seq = self.entries[acc.upper()]["sequence"]["value"]
        pdb = ideal_pdb(seq, "helix", 91.5)
        return {"pdb": pdb, "version": 6, "mean_plddt": 91.5, "start": 1, "end": len(seq),
                "url": f"https://alphafold.ebi.ac.uk/entry/{acc}"}

    def fetch_variants(self, acc):
        self._check()
        if acc.upper() == "P68871":
            return [{"position": 7, "from": "E", "to": "V", "description": "in SKCA; Hb S", "id": "VAR_1"},
                    {"position": 7, "from": "E", "to": "K", "description": "in Hb C", "id": "VAR_2"}]
        return []

    def fetch_alphamissense(self, acc):
        from protein_analyzer.uniprot import UniProtNotFound
        self._check()
        if acc.upper() not in ("P68871", "P69905"):
            raise UniProtNotFound(f"No AlphaMissense predictions for {acc} (they exist for human proteins only).")
        seq = self.entries[acc.upper()]["sequence"]["value"]
        scores, classes = {}, {}
        for i, c in enumerate(seq[:20], 1):
            for a in "ACDEFGHIKLMNPQRSTVWY":
                if a != c:
                    v = round(((i * 7 + ord(a)) % 100) / 100, 2)
                    scores[f"{c}{i}{a}"] = v
                    classes[f"{c}{i}{a}"] = "LPath" if v > 0.564 else "Amb" if v > 0.34 else "LBen"
        return {"scores": scores, "classes": classes, "source": "fake"}

    def fetch_interactions(self, acc, limit=15, min_score=0.4):
        from protein_analyzer.uniprot import UniProtNotFound
        self._check()
        if acc.upper() not in self.entries:
            raise UniProtNotFound(f"STRING has no entry for {acc}.")
        ev = {"experiments": 0.9, "databases": 0.5, "text mining": 0.3, "co-expression": 0, "neighbourhood": 0,
              "gene fusion": 0, "co-occurrence": 0}
        return {"query": "HBA1", "taxon": 9606, "organism": "Homo sapiens", "annotation": "",
                "nodes": [{"id": "q", "name": "HBA1", "is_query": True, "score": None},
                          {"id": "p1", "name": "HBB", "is_query": False, "score": 0.999},
                          {"id": "p2", "name": "AHSP", "is_query": False, "score": 0.98}],
                "edges": [{"a": "q", "b": "p1", "score": 0.999, "evidence": ev}, {"a": "q", "b": "p2", "score": 0.98, "evidence": ev},
                          {"a": "p1", "b": "p2", "score": 0.7, "evidence": ev}],
                "url": "https://string-db.org/network/q"}

    def fetch_pathways(self, acc, taxon=None):
        self._check()
        if acc.upper() != "P69905":
            return []
        return [{"id": "R-HSA-1237044", "name": "Erythrocytes take up oxygen", "disease": False, "inferred": False,
                 "has_diagram": True, "url": "https://reactome.org/content/detail/R-HSA-1237044"},
                {"id": "R-HSA-9999", "name": "Disease variant", "disease": True, "inferred": False,
                 "has_diagram": False, "url": "https://reactome.org/content/detail/R-HSA-9999"}]

    def fetch_pfam(self, acc):
        self._check()
        return [{"accession": "PF00042", "name": "Globin", "type": "family", "start": 2, "end": 140,
                 "url": "https://www.ebi.ac.uk/interpro/entry/pfam/PF00042/"}] if acc.upper() in self.entries else []

    def iter_search(self, query, limit=None, fields=None):
        """Deterministic fake per query (class keyword / split / taxon group), so dataset jobs run offline."""
        import zlib

        import numpy as np
        from protein_analyzer.uniprot import parse_entry
        self._check()
        first = query.split()[0]
        kw = first.split(":")[1] if first.startswith("keyword:") else "OTHER"
        holdout = "TO *]" in query
        qid = zlib.crc32(query.encode())
        rng = np.random.default_rng(qid)
        bias = {"KW-0418": "KRE", "KW-0645": "LIVF"}.get(kw, "DSTN")
        aa = "ACDEFGHIKLMNPQRSTVWY"
        p = np.array([4.0 if c in bias else 1.0 for c in aa]); p /= p.sum()
        for i in range(limit or 30):
            seq = "".join(rng.choice(list(aa), size=int(rng.integers(80, 200)), p=p))
            yield parse_entry(raw_entry(f"X{qid % 10**7:07d}{i:03d}", seq, f"fake {kw}",
                                        created="2022-05-01" if holdout else "2010-05-01"))
