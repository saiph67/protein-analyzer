"""Configuration from environment variables (works from source and from a PyInstaller build)."""
import os
import sys
from pathlib import Path

if getattr(sys, "frozen", False):  # PyInstaller one-folder build
    BUNDLE_DIR = Path(getattr(sys, "_MEIPASS", Path(sys.executable).parent))  # templates/static
    PROJECT_DIR = Path(sys.executable).resolve().parent                       # writable data next to the .exe
else:
    BUNDLE_DIR = PROJECT_DIR = Path(__file__).resolve().parent.parent

CACHE_DIR = Path(os.environ.get("PA_CACHE_DIR", PROJECT_DIR / ".cache" / "uniprot"))
RESULTS_DIR = Path(os.environ.get("PA_RESULTS_DIR", PROJECT_DIR / "results"))
DATA_DIR = Path(os.environ.get("PA_DATA_DIR", PROJECT_DIR / "userdata"))  # projects
PORT = int(os.environ.get("PA_PORT", "5000"))
OPEN_BROWSER = os.environ.get("PA_OPEN_BROWSER", "1") != "0"
ESM_MODEL = os.environ.get("PA_ESM_MODEL", "facebook/esm2_t12_35M_UR50D")
# Shared-server mode (see DEPLOY.md): disables training jobs and projects, adds per-visitor rate limits.
PUBLIC = os.environ.get("PA_PUBLIC", "0") == "1"
PUBLIC_RATE_PER_MIN = int(os.environ.get("PA_RATE_PER_MIN", "30"))
# Local use binds to 127.0.0.1 only. A public server (serve.py) may bind elsewhere via PA_HOST.
HOST = os.environ.get("PA_HOST", "127.0.0.1") if PUBLIC else "127.0.0.1"
MAX_UPLOAD_BYTES = 5 * 1024 * 1024
