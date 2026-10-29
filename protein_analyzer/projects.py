"""Projects: named collections of saved proteins and results with notes, stored in one local JSON file."""
from __future__ import annotations

import json
import threading
import time
import uuid
from pathlib import Path

KINDS = {"protein", "sequence", "analysis", "comparison", "alignment", "classification", "mutation",
         "blast", "superposition", "digest", "orf", "note"}
MAX_PROJECTS = 200
MAX_ITEMS = 500
MAX_ITEM_BYTES = 2_000_000


class ProjectStore:
    def __init__(self, path: str | Path):
        self.path = Path(path)
        self._lock = threading.Lock()

    def _read(self) -> dict:
        if not self.path.exists():
            return {"projects": []}
        try:
            return json.loads(self.path.read_text(encoding="utf-8"))
        except ValueError:
            backup = self.path.with_suffix(f".corrupt-{int(time.time())}.json")
            self.path.replace(backup)
            return {"projects": []}

    def _write(self, data: dict):
        self.path.parent.mkdir(parents=True, exist_ok=True)
        tmp = self.path.with_suffix(".tmp")
        tmp.write_text(json.dumps(data, indent=1), encoding="utf-8")
        tmp.replace(self.path)

    def _find(self, data, pid):
        for p in data["projects"]:
            if p["id"] == pid:
                return p
        raise KeyError("Project not found.")

    # ---------------------------------------------------------------- projects
    def list(self) -> list[dict]:
        with self._lock:
            return [{k: p[k] for k in ("id", "name", "created", "updated")} | {"n_items": len(p["items"])}
                    for p in sorted(self._read()["projects"], key=lambda p: -p["updated"])]

    def create(self, name: str) -> dict:
        name = (name or "").strip()[:100]
        if not name:
            raise ValueError("Give the project a name.")
        with self._lock:
            data = self._read()
            if len(data["projects"]) >= MAX_PROJECTS:
                raise ValueError(f"At most {MAX_PROJECTS} projects.")
            now = time.time()
            p = {"id": uuid.uuid4().hex[:10], "name": name, "notes": "", "created": now, "updated": now, "items": []}
            data["projects"].append(p)
            self._write(data)
            return p

    def get(self, pid: str) -> dict:
        with self._lock:
            return self._find(self._read(), pid)

    def update(self, pid: str, name: str | None = None, notes: str | None = None) -> dict:
        with self._lock:
            data = self._read()
            p = self._find(data, pid)
            if name is not None:
                if not name.strip():
                    raise ValueError("Name can't be empty.")
                p["name"] = name.strip()[:100]
            if notes is not None:
                p["notes"] = notes[:20000]
            p["updated"] = time.time()
            self._write(data)
            return p

    def delete(self, pid: str):
        with self._lock:
            data = self._read()
            self._find(data, pid)
            data["projects"] = [p for p in data["projects"] if p["id"] != pid]
            self._write(data)

    # ---------------------------------------------------------------- items
    def add_item(self, pid: str, kind: str, title: str, data_: dict, note: str = "",
                 provenance: dict | None = None) -> dict:
        if kind not in KINDS:
            raise ValueError(f"Unknown item type '{kind}'.")
        if len(json.dumps(data_)) > MAX_ITEM_BYTES:
            raise ValueError("That result is too large to save.")
        with self._lock:
            data = self._read()
            p = self._find(data, pid)
            if len(p["items"]) >= MAX_ITEMS:
                raise ValueError(f"A project holds at most {MAX_ITEMS} items.")
            item = {"id": uuid.uuid4().hex[:10], "kind": kind, "title": (title or kind)[:200],
                    "note": (note or "")[:5000], "data": data_, "added": time.time(), "provenance": provenance}
            p["items"].append(item)
            p["updated"] = time.time()
            self._write(data)
            return item

    def update_item(self, pid: str, iid: str, note: str) -> dict:
        with self._lock:
            data = self._read()
            p = self._find(data, pid)
            item = next((i for i in p["items"] if i["id"] == iid), None)
            if not item:
                raise KeyError("Item not found.")
            item["note"] = (note or "")[:5000]
            p["updated"] = time.time()
            self._write(data)
            return item

    def delete_item(self, pid: str, iid: str):
        with self._lock:
            data = self._read()
            p = self._find(data, pid)
            before = len(p["items"])
            p["items"] = [i for i in p["items"] if i["id"] != iid]
            if len(p["items"]) == before:
                raise KeyError("Item not found.")
            p["updated"] = time.time()
            self._write(data)


    # ---------------------------------------------------------------- backup / import
    def backup(self) -> dict:
        from . import __version__
        with self._lock:
            projects = self._read()["projects"]
        return {"format": BACKUP_FORMAT, "app_version": __version__, "exported": time.time(), "projects": projects}

    def import_data(self, payload: dict) -> list[dict]:
        """Import a single exported project (export.json) or a full backup; always creates new copies."""
        if not isinstance(payload, dict):
            raise ValueError("That file isn't a Protein Analyzer project or backup.")
        incoming = payload["projects"] if isinstance(payload.get("projects"), list) else [payload]
        cleaned = []
        for p in incoming:
            if not isinstance(p, dict) or not isinstance(p.get("name"), str) or not isinstance(p.get("items", []), list):
                raise ValueError("That file isn't a Protein Analyzer project or backup.")
            items = []
            for it in p.get("items", [])[:MAX_ITEMS]:
                if not isinstance(it, dict) or it.get("kind") not in KINDS or not isinstance(it.get("data", {}), dict):
                    continue
                if len(json.dumps(it.get("data", {}))) > MAX_ITEM_BYTES:
                    continue
                items.append({"id": uuid.uuid4().hex[:10], "kind": it["kind"], "title": str(it.get("title", it["kind"]))[:200],
                              "note": str(it.get("note", ""))[:5000], "data": it.get("data", {}),
                              "added": float(it.get("added", time.time())) if isinstance(it.get("added"), (int, float)) else time.time(),
                              "provenance": it.get("provenance") if isinstance(it.get("provenance"), dict) else None})
            cleaned.append((p["name"].strip()[:100] or "Imported project", str(p.get("notes", ""))[:20000], items))
        with self._lock:
            data = self._read()
            if len(data["projects"]) + len(cleaned) > MAX_PROJECTS:
                raise ValueError(f"Importing would exceed {MAX_PROJECTS} projects.")
            names = {q["name"] for q in data["projects"]}
            created = []
            for name, notes, items in cleaned:
                nm, k = name, 2
                while nm in names:
                    nm = f"{name} ({k})"; k += 1
                names.add(nm)
                now = time.time()
                proj = {"id": uuid.uuid4().hex[:10], "name": nm, "notes": notes, "created": now, "updated": now, "items": items}
                data["projects"].append(proj)
                created.append({"id": proj["id"], "name": nm, "n_items": len(items)})
            self._write(data)
            return created


BACKUP_FORMAT = "protein-analyzer-projects"


def to_markdown(p: dict) -> str:
    import datetime as dt

    ts = lambda t: dt.datetime.fromtimestamp(t).strftime("%Y-%m-%d %H:%M")
    out = [f"# {p['name']}", "", f"Created {ts(p['created'])} · last changed {ts(p['updated'])} · "
           f"{len(p['items'])} item(s)", ""]
    if p.get("notes"):
        out += ["## Notes", "", p["notes"], ""]
    for i, it in enumerate(p["items"], 1):
        out += [f"## {i}. {it['title']}", f"*{it['kind']}, saved {ts(it['added'])}*", ""]
        if it.get("note"):
            out += [f"> {line}" for line in it["note"].splitlines()] + [""]
        out += _summary_lines(it["kind"], it["data"]) + [""]
        pv = it.get("provenance") or {}
        if pv:
            out.append(f"<sub>{pv.get('app', '')} · {pv.get('generated', '')}"
                       + (f" · UniProtKB release {pv['uniprot_release']}" if pv.get("uniprot_release") else "") + "</sub>")
            out.append("")
    out += methods_section(p)
    return "\n".join(out)


def project_provenance(p: dict) -> dict:
    """Methods, sources (with citations/licences), app versions and UniProt releases used in a project."""
    from .provenance import KIND_METHODS, KIND_SOURCES, METHODS, SOURCES

    kinds = {it["kind"] for it in p["items"]}
    order = [it["kind"] for it in p["items"]]
    kinds = sorted(kinds, key=order.index)
    return {
        "methods": [METHODS[m] for m in dict.fromkeys(m for k in kinds for m in KIND_METHODS.get(k, []))],
        "sources": [SOURCES[s] for s in dict.fromkeys(s for k in kinds for s in KIND_SOURCES.get(k, []))],
        "releases": sorted({(it.get("provenance") or {}).get("uniprot_release") for it in p["items"]} - {None}),
        "apps": sorted({(it.get("provenance") or {}).get("app") for it in p["items"]} - {None}),
    }


def methods_section(p: dict) -> list[str]:
    """'Methods & data sources' for everything in the project, with citations and licences."""
    pv = project_provenance(p)
    if not (pv["methods"] or pv["sources"]):
        return []
    out = ["## Methods & data sources", ""]
    if pv["apps"] or pv["releases"]:
        out += [f"Generated with {', '.join(pv['apps']) or 'Protein Analyzer'}"
                + (f"; UniProtKB release(s) {', '.join(pv['releases'])}" if pv["releases"] else "") + ".", ""]
    out += [f"- {m}" for m in pv["methods"]]
    if pv["sources"]:
        out += ["", "Data sources (please cite):", ""]
        out += [f"- **{s['name']}** ({s['licence']}), {s['url']}. {s['cite']}" for s in pv["sources"]]
    out += ["", "_Check each source's own citation page for its latest reference._", ""]
    return out


def _summary_lines(kind: str, d: dict) -> list[str]:
    f = lambda x, n=2: "–" if x is None else f"{x:.{n}f}"
    if kind in ("protein", "analysis"):
        a = d.get("analysis") or d
        head = [f"- Accession: {d['accession']} ({d.get('name') or ''}, {d.get('organism') or ''})"] if d.get("accession") else []
        if "length" in a:
            head += [f"- Length {a['length']} aa · MW {f(a.get('molecular_weight', 0) / 1000)} kDa · GRAVY {f(a.get('gravy'), 3)} · "
                     f"pI {f(a.get('isoelectric_point'))} · charge {f(a.get('charge_ph7'))} · TM segments {a.get('n_tm_segments')}"]
        return head
    if kind == "comparison":
        return [f"- {d.get('name_a')} vs {d.get('name_b')} ({d.get('mode')}): identity {d.get('identity')}%, "
                f"similarity {d.get('similarity')}%, score {d.get('score')}, gaps {d.get('gaps')}"]
    if kind == "classification":
        return [f"- {r['name']}: **{r['prediction']}** ({f(100 * r['confidence'], 0)}%), 2nd {r['top2'][1]['label']}"
                f"{'' if r.get('reliable', True) else ' — unreliable (dissimilar to training data)'}"
                for r in d.get("predictions", [])]
    if kind == "mutation":
        return [f"- {d.get('mutation')}: {d.get('from_name')} → {d.get('to_name')}, BLOSUM62 {d.get('blosum62')}"
                + (f", AlphaMissense {d['alphamissense']:.2f} ({d.get('alphamissense_class')})" if d.get("alphamissense") is not None else "")]
    if kind == "alignment":
        return [f"- {len(d.get('names', []))} proteins, alignment length {d.get('length')}", f"- Tree: `{d.get('newick')}`"]
    if kind == "blast":
        return [f"- {h['accession']} {h.get('name') or ''}: {h.get('identity')}% identity, E={h.get('evalue')}"
                for h in d.get("hits", [])[:10]]
    if kind == "superposition":
        return [f"- {d.get('name_a')} vs {d.get('name_b')}: TM-score {f(d.get('tm_score'))}, core RMSD {d.get('rmsd_core')} Å "
                f"({d.get('core_pairs')} residues), sequence identity {d.get('sequence_identity')}%"]
    if kind == "sequence":
        return [f"- {len(d.get('sequence', ''))} aa: `{d.get('sequence', '')[:60]}{'…' if len(d.get('sequence', '')) > 60 else ''}`"]
    return []
