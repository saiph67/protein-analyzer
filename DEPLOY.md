# Putting Protein Analyzer online (for a class or group)

By default the app only listens on your own computer (`127.0.0.1`). To let others use it from a link,
run it in **public mode** on a server. Nothing in this folder puts the app online by itself; you choose where.

## What public mode changes (`PA_PUBLIC=1`)

| | Your own computer | Public mode |
|---|---|---|
| Listens on | 127.0.0.1 only | `PA_HOST` (e.g. 0.0.0.0) |
| Server | Flask development server | waitress (production) via `serve.py` |
| Build dataset / train models | ✓ | **off** (a visitor could tie up the server for minutes) |
| Projects | ✓ | **off** (all visitors would share one project list) |
| Requests that call UniProt / AlphaFold / EBI | unlimited | **30 per minute per visitor** (`PA_RATE_PER_MIN`) |
| Everything else (analyze, search, 3D, compare, mutations, classify, explain, lessons, BLAST) | ✓ | ✓ |

BLAST still asks each visitor for their own email (an EBI requirement); the server never stores it.

## Option A: Docker (any cloud host or a school server)

```bash
docker build -t protein-analyzer .
docker run -p 8000:8000 protein-analyzer
```

Open `http://<server>:8000`. Most container hosts (Render, Railway, Fly.io, Google Cloud Run, Azure Container
Apps…) can deploy straight from this Dockerfile; set the port to 8000. About 512 MB RAM is enough. ESM-2
is left out of the image (it would add about 1 GB); the hand-crafted model is used.

## Option B: a plain server with Python

```bash
pip install -r requirements.txt
PA_PUBLIC=1 PA_HOST=0.0.0.0 PA_PORT=8000 python serve.py
```

## Before you share the link

- **Put it behind HTTPS** (most hosts do this automatically; on your own server use a reverse proxy such as Caddy or nginx).
- Leave **public mode on**. Don't expose `python app.py` (the development server) to the internet.
- The app has no logins. If the class should be the only users, use your host's access controls (a password
  or single sign-on in front of the app) or keep it on the school network.
- Be fair to the free public data services: UniProt, AlphaFold DB, InterPro and EBI are shared resources. The
  built-in cache and rate limit help; for very large classes, lower `PA_RATE_PER_MIN`.
- Results and caches live inside the container and reset when it restarts. That's fine, because nothing
  user-specific is stored in public mode.
