"""EBI Job Dispatcher client for BLAST searches against all of UniProtKB.

EBI requires a real email address with every job (they use it only to contact users about problems).
The app never stores it: the browser keeps it locally and sends it with each search the user starts.
Docs: https://www.ebi.ac.uk/Tools/common/tools/help/
"""
from __future__ import annotations

import re
import threading
import time

import requests

from .uniprot import UniProtBadRequest, UniProtError, UniProtNotFound, UniProtUnavailable

BLAST_URL = "https://www.ebi.ac.uk/Tools/services/rest/ncbiblast"
# Names verified against .../ncbiblast/parameterdetails/database (2026-10).
DATABASES = {
    "uniprotkb_swissprot": "UniProtKB/Swiss-Prot: reviewed proteins (fast)",
    "uniref50": "UniRef50: one representative per 50%-identity cluster (broad, fast)",
    "uniprotkb_human": "UniProtKB: human only",
    "uniprotkb": "UniProtKB: everything, reviewed + unreviewed (slowest)",
}
EMAIL_RE = re.compile(r"^[^@\s]+@[^@\s]+\.[^@\s]+$")
JOB_RE = re.compile(r"^ncbiblast-[A-Za-z0-9_.-]+$")


class BlastClient:
    def __init__(self, session: requests.Session | None = None, min_interval: float = 1.0, timeout: float = 30):
        self.session = session or requests.Session()
        self.session.headers.setdefault("User-Agent", "protein-analyzer/1.0 (local research app)")
        self.min_interval, self.timeout = min_interval, timeout
        self._lock, self._last = threading.Lock(), 0.0
        self._results: dict[str, dict] = {}

    def _throttle(self):
        with self._lock:
            wait = self._last + self.min_interval - time.monotonic()
            if wait > 0:
                time.sleep(wait)
            self._last = time.monotonic()

    def _req(self, method, path, **kw):
        self._throttle()
        try:
            r = self.session.request(method, f"{BLAST_URL}/{path}", timeout=self.timeout, **kw)
        except (requests.ConnectionError, requests.Timeout):
            raise UniProtUnavailable("Could not reach EBI (network error or timeout). Try again shortly.")
        if r.status_code == 400:
            msg = re.search(r"<description>(.*?)</description>", r.text, re.S)
            raise UniProtBadRequest("EBI rejected the search: " + (msg.group(1).strip() if msg else r.text[:200]))
        if r.status_code == 404:
            raise UniProtNotFound("EBI does not know this search job (results expire after about 7 days).")
        if r.status_code == 429 or r.status_code >= 500:
            raise UniProtUnavailable(f"EBI is busy or having problems (HTTP {r.status_code}). Try again later.")
        if r.status_code != 200:
            raise UniProtError(f"Unexpected EBI response HTTP {r.status_code}.")
        return r

    def submit(self, sequence: str, email: str, database: str = "uniprotkb_swissprot",
               max_hits: int = 50, evalue: str = "1e-3") -> str:
        if not EMAIL_RE.match(email or ""):
            raise UniProtBadRequest("EBI requires a valid email address for every search.")
        if database not in DATABASES:
            raise UniProtBadRequest("Unknown database.")
        data = {"email": email, "program": "blastp", "stype": "protein", "sequence": sequence,
                "database": database, "alignments": str(max_hits), "scores": str(max_hits),
                "exp": evalue, "matrix": "BLOSUM62", "title": "protein-analyzer"}
        return self._req("POST", "run", data=data).text.strip()

    def status(self, job_id: str) -> str:
        """RUNNING, FINISHED, ERROR, FAILURE, NOT_FOUND or QUEUED."""
        _check_job(job_id)
        return self._req("GET", f"status/{job_id}").text.strip()

    def result(self, job_id: str) -> dict:
        _check_job(job_id)
        if job_id not in self._results:
            body = self._req("GET", f"result/{job_id}/json").json()
            self._results[job_id] = parse_blast_json(body)
        return self._results[job_id]


def _check_job(job_id: str):
    if not JOB_RE.match(job_id or ""):
        raise UniProtBadRequest("Not a valid EBI job id.")


def parse_blast_json(body: dict) -> dict:
    """Flatten EBI NCBI-BLAST JSON output into the fields the app shows."""
    hits = []
    for h in body.get("hits", []) or []:
        hsps = h.get("hit_hsps") or [{}]
        best = max(hsps, key=lambda s: s.get("hsp_bit_score") or 0)
        acc = h.get("hit_acc") or h.get("hit_id", "").split("|")[-1]
        hits.append({
            "accession": acc,
            "id": h.get("hit_id"),
            "name": h.get("hit_uni_de") or h.get("hit_desc") or h.get("hit_def"),
            "organism": h.get("hit_uni_os"),
            "gene": h.get("hit_uni_gn"),
            "length": h.get("hit_len"),
            "identity": best.get("hsp_identity"),
            "positives": best.get("hsp_positive"),
            "evalue": best.get("hsp_expect"),
            "bit_score": best.get("hsp_bit_score"),
            "align_len": best.get("hsp_align_len"),
            "query_range": [best.get("hsp_query_from"), best.get("hsp_query_to")],
            "hit_range": [best.get("hsp_hit_from"), best.get("hsp_hit_to")],
            "qseq": best.get("hsp_qseq"), "mseq": best.get("hsp_mseq"), "hseq": best.get("hsp_hseq"),
        })
    return {"program": body.get("program"), "database": body.get("dbs", [{}])[0].get("name") if body.get("dbs") else None,
            "query_length": body.get("query_len"), "hits": hits}
