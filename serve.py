"""Production server for a shared/classroom deployment (see DEPLOY.md).

    PA_PUBLIC=1 PA_HOST=0.0.0.0 PA_PORT=8000 python serve.py

Uses waitress (a pure-Python production WSGI server that also runs on Windows) instead of Flask's
development server. Public mode turns off training jobs and projects and rate-limits each visitor.
"""
import os

os.environ.setdefault("PA_PUBLIC", "1")
os.environ.setdefault("PA_OPEN_BROWSER", "0")

from waitress import serve  # noqa: E402

from app import create_app  # noqa: E402
from protein_analyzer import config  # noqa: E402

if __name__ == "__main__":
    app = create_app()
    print(f"Serving Protein Analyzer on http://{config.HOST}:{config.PORT} "
          f"(public mode: {config.PUBLIC}, {config.PUBLIC_RATE_PER_MIN} requests/min per visitor)")
    serve(app, host=config.HOST, port=config.PORT, threads=int(os.environ.get("PA_THREADS", "8")))
