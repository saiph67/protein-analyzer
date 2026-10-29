"""Background job runner (thread + id + polling) with cancel and restart-survival.

Job state is saved to a small JSON file, so after the app restarts the UI can say that a job was
interrupted instead of silently forgetting it.
"""
from __future__ import annotations

import json
import threading
import time
import traceback
import uuid
from pathlib import Path


class JobCancelled(Exception):
    pass


class JobManager:
    def __init__(self, state_path: str | Path | None = None):
        self._jobs: dict[str, dict] = {}
        self._cancel: dict[str, threading.Event] = {}
        self._lock = threading.Lock()
        self.state_path = Path(state_path) if state_path else None
        self._load()

    # -------------------------------------------------------------- persistence
    def _load(self):
        if not (self.state_path and self.state_path.exists()):
            return
        try:
            jobs = json.loads(self.state_path.read_text())
        except ValueError:
            return
        for j in jobs:
            if j.get("state") == "running":
                j.update(state="interrupted", finished=j.get("updated"),
                         error="The app was closed while this job was running. Start it again.")
            self._jobs[j["id"]] = j

    def _save(self):  # caller holds the lock
        if not self.state_path:
            return
        recent = sorted(self._jobs.values(), key=lambda j: j["started"])[-20:]
        tmp = self.state_path.with_suffix(".tmp")
        try:
            self.state_path.parent.mkdir(parents=True, exist_ok=True)
            tmp.write_text(json.dumps(recent, default=str))
            tmp.replace(self.state_path)
        except OSError:
            pass

    # -------------------------------------------------------------- API
    def running(self, kind: str | None = None) -> dict | None:
        with self._lock:
            for j in self._jobs.values():
                if j["state"] in ("running", "cancelling") and (kind is None or j["kind"] == kind):
                    return dict(j)
        return None

    def latest(self) -> dict | None:
        with self._lock:
            if not self._jobs:
                return None
            return dict(max(self._jobs.values(), key=lambda j: j["started"]))

    def start(self, kind: str, fn, *args, **kwargs) -> str:
        """Run fn(*args, progress=cb, **kwargs) in a daemon thread. Only one job runs at a time.

        The progress callback raises JobCancelled once cancel() was requested, so long jobs stop at
        their next progress report.
        """
        with self._lock:
            busy = [j for j in self._jobs.values() if j["state"] in ("running", "cancelling")]
            if busy:
                raise RuntimeError(f"A {busy[0]['kind']} job is already running; wait for it or cancel it.")
            jid = uuid.uuid4().hex[:12]
            now = time.time()
            self._jobs[jid] = {"id": jid, "kind": kind, "state": "running", "progress": 0.0,
                               "message": "starting", "result": None, "error": None,
                               "started": now, "updated": now, "finished": None}
            self._cancel[jid] = threading.Event()
            self._save()

        def progress(frac, msg):
            if self._cancel[jid].is_set():
                raise JobCancelled()
            with self._lock:
                self._jobs[jid].update(progress=round(float(frac), 3), message=msg, updated=time.time())
                self._save()

        def run():
            try:
                result = fn(*args, progress=progress, **kwargs)
                update = {"state": "done", "progress": 1.0, "result": result}
            except JobCancelled:
                update = {"state": "cancelled", "message": "Cancelled. Nothing was overwritten."}
            except Exception as e:  # report any failure to the UI
                traceback.print_exc()
                update = {"state": "error", "error": str(e) or e.__class__.__name__}
            with self._lock:
                self._jobs[jid].update(finished=time.time(), updated=time.time(), **update)
                self._save()

        threading.Thread(target=run, daemon=True, name=f"job-{kind}-{jid}").start()
        return jid

    def cancel(self, jid: str) -> bool:
        with self._lock:
            j = self._jobs.get(jid)
            if not j or j["state"] != "running":
                return False
            self._cancel[jid].set()
            j.update(state="cancelling", message="Cancelling at the next checkpoint…")
            self._save()
            return True

    def get(self, jid: str) -> dict | None:
        with self._lock:
            j = self._jobs.get(jid)
            return dict(j) if j else None
