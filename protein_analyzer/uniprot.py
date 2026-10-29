"""Small, polite client for the UniProtKB REST API (https://rest.uniprot.org/uniprotkb).

Verified against the live API (2026-10): field names below return HTTP 200, `size` max is 500,
pagination uses `Link: <...>; rel="next"`, totals come from the `X-Total-Results` header.
"""
from __future__ import annotations

import hashlib
import json
import os
import re
import threading
import time
from pathlib import Path

import requests

BASE_URL = "https://rest.uniprot.org/uniprotkb"
ALPHAFOLD_API = "https://alphafold.ebi.ac.uk/api/prediction"
INTERPRO_API = "https://www.ebi.ac.uk/interpro/api"
STRING_API = "https://string-db.org/api/json"
REACTOME_API = "https://reactome.org/ContentService"
MAX_PAGE_SIZE = 500
# Sequence annotations shown next to the app's own predictions (verified live 2026-10).
FEATURE_FIELDS = ["ft_transmem", "ft_intramem", "ft_signal", "ft_topo_dom", "ft_domain", "ft_region",
                  "ft_motif", "ft_zn_fing", "ft_dna_bind", "ft_act_site", "ft_binding"]
DEFAULT_FIELDS = [
    "accession", "id", "protein_name", "gene_names", "organism_name", "organism_id",
    "length", "sequence", "cc_function", "go", "ec", "keyword", "date_created", "reviewed",
] + FEATURE_FIELDS
SEARCH_FIELDS = ["accession", "id", "protein_name", "gene_names", "organism_name", "length", "reviewed"]
ACCESSION_RE = re.compile(
    r"^([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})(-\d+)?$"
)
_LINK_NEXT_RE = re.compile(r'<([^>]+)>;\s*rel="next"')


class UniProtError(Exception):
    """Base error; `.status` is the HTTP status to return to the app's user."""
    status = 502


class UniProtUnavailable(UniProtError):
    status = 503


class UniProtRateLimited(UniProtError):
    status = 503


class UniProtBadRequest(UniProtError):
    status = 400


class UniProtNotFound(UniProtError):
    status = 404


def is_accession(text: str) -> bool:
    return bool(ACCESSION_RE.match((text or "").strip().upper()))


class UniProtClient:
    def __init__(self, cache_dir: str | os.PathLike | None = None, timeout: float = 30,
                 max_retries: int = 4, min_interval: float = 0.34, cache_ttl: float = 7 * 86400,
                 session: requests.Session | None = None):
        self.cache_dir = Path(cache_dir) if cache_dir else None
        if self.cache_dir:
            self.cache_dir.mkdir(parents=True, exist_ok=True)
        self.timeout = timeout
        self.max_retries = max_retries
        self.min_interval = min_interval  # ~3 requests/sec
        self.cache_ttl = cache_ttl
        self.session = session or requests.Session()
        self.session.headers.setdefault("User-Agent", "protein-analyzer/1.0 (local research app)")
        self._lock = threading.Lock()
        self._last = 0.0

    # ------------------------------------------------------------------ low level
    def _cache_path(self, url: str, params: dict | None) -> Path | None:
        if not self.cache_dir:
            return None
        key = url + "?" + json.dumps(params or {}, sort_keys=True)
        return self.cache_dir / (hashlib.sha1(key.encode()).hexdigest() + ".json")

    def _throttle(self):
        with self._lock:
            wait = self._last + self.min_interval - time.monotonic()
            if wait > 0:
                time.sleep(wait)
            self._last = time.monotonic()

    def _get(self, url: str, params: dict | None = None, text: bool = False,
             service: str = "UniProt") -> tuple[dict | str, dict]:
        """GET JSON (or text) with caching, throttling and retries. Returns (body, headers)."""
        cp = self._cache_path(url, params)
        if cp and cp.exists() and time.time() - cp.stat().st_mtime < self.cache_ttl:
            try:
                cached = json.loads(cp.read_text(encoding="utf-8"))
                self._note_release(cached["headers"])
                return cached["body"], cached["headers"]
            except (ValueError, KeyError):
                pass
        last_exc: Exception | None = None
        for attempt in range(self.max_retries + 1):
            self._throttle()
            try:
                r = self.session.get(url, params=params, timeout=self.timeout)
            except (requests.ConnectionError, requests.Timeout) as e:
                last_exc = UniProtUnavailable(
                    f"Could not reach {service} (network error or timeout). Check your internet "
                    "connection or try again shortly.")
                last_exc.__cause__ = e
                time.sleep(min(2 ** attempt, 16))
                continue
            if r.status_code == 200:
                body = r.text if text else r.json()
                headers = {"link": r.headers.get("Link", ""), "total": r.headers.get("X-Total-Results"),
                           "release": r.headers.get("X-UniProt-Release"),
                           "release_date": r.headers.get("X-UniProt-Release-Date")}
                self._note_release(headers)
                if cp:
                    tmp = cp.with_suffix(".tmp")
                    tmp.write_text(json.dumps({"body": body, "headers": headers}), encoding="utf-8")
                    os.replace(tmp, cp)
                return body, headers
            if r.status_code == 400:
                raise UniProtBadRequest(f"{service} rejected the query: " + _messages(r))
            if r.status_code == 404:
                raise UniProtNotFound(f"Not found in {service}.")
            if r.status_code == 429 or r.status_code >= 500:
                retry_after = r.headers.get("Retry-After")
                delay = float(retry_after) if retry_after and retry_after.isdigit() else min(2 ** attempt, 16)
                last_exc = (UniProtRateLimited(f"{service} is rate-limiting requests; please wait a moment.")
                            if r.status_code == 429 else
                            UniProtUnavailable(f"{service} is having problems (HTTP {r.status_code}). Try again later."))
                time.sleep(delay)
                continue
            raise UniProtError(f"Unexpected {service} response HTTP {r.status_code}: {_messages(r)}")
        raise last_exc or UniProtUnavailable(f"{service} did not respond.")

    def _note_release(self, headers: dict):
        if headers.get("release"):
            self.release = {"release": headers["release"], "date": headers.get("release_date")}

    def uniprot_release(self) -> dict:
        """Current UniProtKB release (e.g. {'release': '2026_04', 'date': '...'}), fetched once if unknown."""
        if not getattr(self, "release", None):
            try:  # direct request: older cached responses don't carry the release header
                self._throttle()
                r = self.session.get(f"{BASE_URL}/search", params={"query": "accession:P69905", "size": 0},
                                     timeout=self.timeout)
                self._note_release({"release": r.headers.get("X-UniProt-Release"),
                                    "release_date": r.headers.get("X-UniProt-Release-Date")})
            except requests.RequestException:
                return {}
        return getattr(self, "release", None) or {}

    # ------------------------------------------------------------------ public API
    def fetch_entry(self, accession: str, _follow: int = 2) -> dict:
        acc = (accession or "").strip().upper()
        if not is_accession(acc):
            raise UniProtBadRequest(f"'{accession}' is not a valid UniProt accession (e.g. P69905).")
        try:
            body, _ = self._get(f"{BASE_URL}/{acc}", {"fields": ",".join(DEFAULT_FIELDS), "format": "json"})
        except UniProtBadRequest:
            raise UniProtNotFound(f"No UniProt entry for accession {acc}.")
        if body.get("entryType") == "Inactive":
            reason = body.get("inactiveReason", {})
            targets = reason.get("mergeDemergeTo") or []
            if targets and _follow:
                return self.fetch_entry(targets[0], _follow - 1)
            raise UniProtNotFound(f"{acc} is obsolete in UniProt ({reason.get('inactiveReasonType', 'deleted')}).")
        return parse_entry(body)

    def search(self, query: str, size: int = 25, cursor: str | None = None,
               fields: list[str] | None = None) -> dict:
        """One page of results: {'results': [...], 'total': int, 'next_cursor': str|None}."""
        if not (query or "").strip():
            raise UniProtBadRequest("Empty query.")
        size = max(1, min(int(size), MAX_PAGE_SIZE))
        params = {"query": query, "size": size, "format": "json",
                  "fields": ",".join(fields or SEARCH_FIELDS)}
        if cursor:
            params["cursor"] = cursor
        body, headers = self._get(f"{BASE_URL}/search", params)
        nxt = _next_link(headers.get("link", ""))
        next_cursor = None
        if nxt:
            m = re.search(r"[?&]cursor=([^&]+)", nxt)
            next_cursor = m.group(1) if m else None
        return {
            "results": [parse_entry(e) for e in body.get("results", [])],
            "total": int(headers["total"]) if headers.get("total") else None,
            "next_cursor": next_cursor,
        }

    def iter_search(self, query: str, limit: int | None = None, fields: list[str] | None = None,
                    page_size: int = MAX_PAGE_SIZE):
        """Yield parsed entries across pages, following Link rel=next."""
        url = f"{BASE_URL}/search"
        params: dict | None = {"query": query, "size": min(page_size, MAX_PAGE_SIZE), "format": "json",
                               "fields": ",".join(fields or DEFAULT_FIELDS)}
        n = 0
        while url:
            body, headers = self._get(url, params)
            for e in body.get("results", []):
                yield parse_entry(e)
                n += 1
                if limit is not None and n >= limit:
                    return
            url, params = _next_link(headers.get("link", "")), None

    def count(self, query: str) -> int:
        _, headers = self._get(f"{BASE_URL}/search", {"query": query, "size": 0, "format": "json"})
        return int(headers.get("total") or 0)

    def fetch_variants(self, accession: str) -> list[dict]:
        """Known natural (often disease-linked) single-residue variants from UniProt."""
        acc = (accession or "").strip().upper()
        if not is_accession(acc):
            raise UniProtBadRequest(f"'{accession}' is not a valid UniProt accession.")
        body, _ = self._get(f"{BASE_URL}/{acc}", {"fields": "accession,ft_variant", "format": "json"})
        out = []
        for f in body.get("features", []) or []:
            loc = f.get("location") or {}
            start, end = (loc.get("start") or {}).get("value"), (loc.get("end") or {}).get("value")
            alt = f.get("alternativeSequence") or {}
            orig, alts = alt.get("originalSequence"), alt.get("alternativeSequences") or []
            if start is None or start != end or not orig or len(orig) != 1:
                continue  # keep simple substitutions only
            for a in alts:
                if len(a) == 1:
                    out.append({"position": start, "from": orig, "to": a, "description": f.get("description", ""),
                                "id": f.get("featureId")})
        return out

    def _alphafold_models(self, acc: str) -> dict:
        missing = UniProtNotFound(f"AlphaFold DB has no predicted structure for {acc}.")
        try:
            models, _ = self._get(f"{ALPHAFOLD_API}/{acc}", service="AlphaFold DB")
        except (UniProtBadRequest, UniProtNotFound):
            raise missing
        models = [m for m in models or [] if m.get("pdbUrl")]
        if not models:
            raise missing
        # prefer the canonical (non-isoform) full-length model
        return next((x for x in models if x.get("uniprotAccession", "").upper() == acc), models[0])

    def fetch_alphamissense(self, accession: str) -> dict:
        """AlphaMissense pathogenicity for every single-residue substitution (human proteins only).

        Returns {'scores': {"E7V": 0.98, ...}, 'classes': {"E7V": "LPath"}, 'source': url}.
        """
        acc = (accession or "").strip().upper()
        if not is_accession(acc):
            raise UniProtBadRequest(f"'{accession}' is not a valid UniProt accession.")
        m = self._alphafold_models(acc)
        url = m.get("amAnnotationsUrl")
        if not url:
            raise UniProtNotFound(f"No AlphaMissense predictions for {acc} (they exist for human proteins only).")
        csv, _ = self._get(url, text=True, service="AlphaFold DB")
        scores, classes = {}, {}
        for line in csv.splitlines()[1:]:
            parts = line.strip().split(",")
            if len(parts) >= 3:
                try:
                    scores[parts[0]] = float(parts[1])
                except ValueError:
                    continue
                classes[parts[0]] = parts[2]
        return {"scores": scores, "classes": classes, "source": url}

    def fetch_pfam(self, accession: str) -> list[dict]:
        """Pfam domain hits from InterPro: [{'accession','name','start','end','url'}]."""
        acc = (accession or "").strip().upper()
        if not is_accession(acc):
            raise UniProtBadRequest(f"'{accession}' is not a valid UniProt accession.")
        try:
            body, _ = self._get(f"{INTERPRO_API}/entry/pfam/protein/uniprot/{acc}", {"format": "json"},
                                service="InterPro")
        except UniProtNotFound:
            return []
        out = []
        for r in (body or {}).get("results", []) if isinstance(body, dict) else []:
            md = r.get("metadata", {})
            for p in r.get("proteins", []):
                for loc in p.get("entry_protein_locations") or []:
                    for frag in loc.get("fragments", []):
                        out.append({"accession": md.get("accession"), "name": md.get("name"),
                                    "type": md.get("type"), "start": frag.get("start"), "end": frag.get("end"),
                                    "url": f"https://www.ebi.ac.uk/interpro/entry/pfam/{md.get('accession')}/"})
        return sorted(out, key=lambda d: d["start"] or 0)

    def fetch_interactions(self, accession: str, limit: int = 15, min_score: float = 0.4) -> dict:
        """Interaction partners from STRING (string-db.org): nodes + scored edges among query and partners."""
        acc = (accession or "").strip().upper()
        if not is_accession(acc):
            raise UniProtBadRequest(f"'{accession}' is not a valid UniProt accession.")
        common = {"caller_identity": "protein-analyzer"}
        try:
            ids, _ = self._get(f"{STRING_API}/get_string_ids", {"identifiers": acc, "limit": 1, **common}, service="STRING")
        except (UniProtNotFound, UniProtBadRequest):
            ids = []
        if not ids:
            raise UniProtNotFound(f"STRING has no entry for {acc}.")
        q = ids[0]
        taxon = q["ncbiTaxonId"]
        partners, _ = self._get(f"{STRING_API}/interaction_partners",
                                {"identifiers": q["stringId"], "species": taxon, "limit": limit,
                                 "required_score": int(min_score * 1000), **common}, service="STRING")
        names = {q["stringId"]: q["preferredName"]}
        for p in partners or []:
            names[p["stringId_B"]] = p["preferredName_B"]
        edges = []
        if len(names) > 1:
            net, _ = self._get(f"{STRING_API}/network",
                               {"identifiers": "\r".join(names), "species": taxon,
                                "required_score": int(min_score * 1000), **common}, service="STRING")
            seen = set()
            for e in net or []:
                key = tuple(sorted((e["stringId_A"], e["stringId_B"])))
                if key in seen or e["stringId_A"] not in names or e["stringId_B"] not in names:
                    continue
                seen.add(key)
                edges.append({"a": e["stringId_A"], "b": e["stringId_B"], "score": e["score"],
                              "evidence": {k: e.get(s, 0) for k, s in
                                           (("experiments", "escore"), ("databases", "dscore"), ("text mining", "tscore"),
                                            ("co-expression", "ascore"), ("neighbourhood", "nscore"),
                                            ("gene fusion", "fscore"), ("co-occurrence", "pscore"))}})
        score_with_query = {p["stringId_B"]: p["score"] for p in partners or []}
        nodes = [{"id": sid, "name": nm, "is_query": sid == q["stringId"],
                  "score": score_with_query.get(sid)} for sid, nm in names.items()]
        return {"query": q["preferredName"], "taxon": taxon, "organism": q.get("taxonName"),
                "annotation": q.get("annotation"), "nodes": nodes, "edges": edges,
                "url": f"https://string-db.org/network/{q['stringId']}"}

    def fetch_pathways(self, accession: str, taxon: int | None = None) -> list[dict]:
        """Reactome pathways containing this protein ([] if Reactome doesn't cover it)."""
        acc = (accession or "").strip().upper()
        if not is_accession(acc):
            raise UniProtBadRequest(f"'{accession}' is not a valid UniProt accession.")
        params = {"species": taxon} if taxon else None
        try:
            body, _ = self._get(f"{REACTOME_API}/data/mapping/UniProt/{acc}/pathways", params, service="Reactome")
        except (UniProtNotFound, UniProtBadRequest):
            return []
        return [{"id": p.get("stId"), "name": p.get("displayName"), "disease": bool(p.get("isInDisease")),
                 "inferred": bool(p.get("isInferred")), "has_diagram": bool(p.get("hasDiagram")),
                 "url": f"https://reactome.org/content/detail/{p.get('stId')}"}
                for p in body or [] if p.get("stId")]

    def fetch_alphafold(self, accession: str) -> dict:
        """AlphaFold DB predicted structure: {'pdb': text, 'version', 'mean_plddt', 'start', 'end', 'url'}."""
        acc = (accession or "").strip().upper()
        if not is_accession(acc):
            raise UniProtBadRequest(f"'{accession}' is not a valid UniProt accession.")
        m = self._alphafold_models(acc)
        pdb, _ = self._get(m["pdbUrl"], text=True, service="AlphaFold DB")
        return {"pdb": pdb, "version": m.get("latestVersion"), "mean_plddt": m.get("globalMetricValue"),
                "start": m.get("uniprotStart"), "end": m.get("uniprotEnd"),
                "url": f"https://alphafold.ebi.ac.uk/entry/{acc}"}


def _next_link(link_header: str) -> str | None:
    m = _LINK_NEXT_RE.search(link_header or "")
    return m.group(1) if m else None


def _messages(r: requests.Response) -> str:
    try:
        return "; ".join(r.json().get("messages", [])) or r.text[:200]
    except ValueError:
        return r.text[:200]


def parse_entry(e: dict) -> dict:
    """Flatten a UniProtKB JSON entry into the fields the app uses."""
    desc = e.get("proteinDescription", {}) or {}
    rec = desc.get("recommendedName") or {}
    subs = desc.get("submissionNames") or []
    name = (rec.get("fullName") or {}).get("value") or ((subs[0].get("fullName") or {}).get("value") if subs else None)
    ecs = [x.get("value") for x in rec.get("ecNumbers", []) if x.get("value")]
    for s in subs:
        ecs += [x.get("value") for x in s.get("ecNumbers", []) if x.get("value")]

    gene = None
    for g in e.get("genes", []) or []:
        for key in ("geneName",):
            if g.get(key):
                gene = g[key].get("value"); break
        if not gene:
            for key in ("orderedLocusNames", "orfNames"):
                if g.get(key):
                    gene = g[key][0].get("value"); break
        if gene:
            break

    function = []
    for c in e.get("comments", []) or []:
        if c.get("commentType") == "FUNCTION":
            function += [t.get("value") for t in c.get("texts", []) if t.get("value")]

    go = {"F": [], "P": [], "C": []}
    for x in e.get("uniProtKBCrossReferences", []) or []:
        if x.get("database") != "GO":
            continue
        props = {p.get("key"): p.get("value") for p in x.get("properties", [])}
        term = props.get("GoTerm", "")
        if len(term) > 2 and term[1] == ":" and term[0] in go:
            go[term[0]].append({"id": x.get("id"), "term": term[2:], "evidence": props.get("GoEvidenceType")})

    features = []
    for f in e.get("features", []) or []:
        loc = f.get("location") or {}
        start, end = (loc.get("start") or {}).get("value"), (loc.get("end") or {}).get("value")
        if start is None and end is None:
            continue
        features.append({"type": f.get("type"), "start": start or end, "end": end or start,
                         "description": f.get("description") or "",
                         "ligand": ((f.get("ligand") or {}).get("name") or "")})

    org = e.get("organism", {}) or {}
    seq = (e.get("sequence") or {}).get("value", "")
    acc = e.get("primaryAccession")
    entry_type = e.get("entryType", "")
    return {
        "accession": acc,
        "id": e.get("uniProtkbId"),
        "name": name,
        "gene": gene,
        "organism": org.get("scientificName"),
        "organism_id": org.get("taxonId"),
        "length": (e.get("sequence") or {}).get("length") or (len(seq) if seq else None),
        "sequence": seq,
        "function": function,
        "keywords": [{"id": k.get("id"), "name": k.get("name"), "category": k.get("category")}
                     for k in e.get("keywords", []) or []],
        "go": go,
        "ec": ecs,
        "features": features,
        "date_created": (e.get("entryAudit") or {}).get("firstPublicDate"),
        "reviewed": "unreviewed" not in entry_type.lower() if entry_type else None,
        "url": f"https://www.uniprot.org/uniprotkb/{acc}/entry" if acc else None,
    }


def save_fasta(entries: list[dict], path: str | os.PathLike) -> str:
    with open(path, "w", encoding="utf-8") as f:
        f.write(to_fasta(entries))
    return str(path)


def to_fasta(entries: list[dict], width: int = 60) -> str:
    out = []
    for e in entries:
        header = f">{e.get('accession') or ''}"
        if e.get("id"):
            header = f">sp|{e['accession']}|{e['id']}" if e.get("reviewed") else f">tr|{e['accession']}|{e['id']}"
        if e.get("name"):
            header += f" {e['name']}"
        if e.get("organism"):
            header += f" OS={e['organism']}"
        if e.get("gene"):
            header += f" GN={e['gene']}"
        seq = e.get("sequence", "")
        out.append(header)
        out += [seq[i:i + width] for i in range(0, len(seq), width)]
    return "\n".join(out) + "\n"
