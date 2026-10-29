"""Optional ESM-2 protein-language-model embeddings (needs `torch` + `transformers`).

Each protein becomes the mean of its per-residue ESM-2 representations. Proteins longer than the model's
1022-residue context are split into windows whose means are averaged (weighted by length).
Embeddings are cached on disk by sequence hash, so re-training and re-classifying is fast.
"""
from __future__ import annotations

import hashlib
import threading
from pathlib import Path

import numpy as np

from . import config

WINDOW = 1000
_LOCK = threading.Lock()
_MODELS: dict = {}


class EmbeddingUnavailable(RuntimeError):
    pass


def esm_available() -> bool:
    try:
        import torch  # noqa: F401
        import transformers  # noqa: F401
        return True
    except Exception:
        return False


def _load(model_name: str):
    if not esm_available():
        raise EmbeddingUnavailable(
            "ESM-2 needs PyTorch and transformers. Install with:  "
            "pip install torch --index-url https://download.pytorch.org/whl/cpu  and  pip install transformers")
    with _LOCK:
        if model_name not in _MODELS:
            import os
            os.environ.setdefault("HF_HUB_DISABLE_SYMLINKS_WARNING", "1")
            os.environ.setdefault("HF_HUB_DISABLE_PROGRESS_BARS", "1")
            import torch
            from transformers import AutoTokenizer, EsmModel
            from transformers.utils import logging as hf_logging

            hf_logging.set_verbosity_error()
            hf_logging.disable_progress_bar()

            tok = AutoTokenizer.from_pretrained(model_name)
            model = EsmModel.from_pretrained(model_name, add_pooling_layer=False).eval()
            torch.set_num_threads(max(1, torch.get_num_threads()))
            _MODELS[model_name] = (tok, model)
        return _MODELS[model_name]


def _key(seq: str) -> str:
    return hashlib.sha1(seq.encode()).hexdigest()


class _Cache:
    def __init__(self, path: Path | None):
        self.path, self.data = path, {}
        if path and path.exists():
            import joblib
            try:
                self.data = joblib.load(path)
            except Exception:
                self.data = {}

    def save(self):
        if self.path:
            import joblib
            self.path.parent.mkdir(parents=True, exist_ok=True)
            tmp = self.path.with_suffix(".tmp")
            joblib.dump(self.data, tmp)
            tmp.replace(self.path)


def embed(seqs: list[str], model_name: str | None = None, cache_dir: Path | None = None,
          batch_tokens: int = 4000, progress=None) -> np.ndarray:
    """(n, d) matrix of mean-pooled ESM-2 embeddings."""
    model_name = model_name or config.ESM_MODEL
    cache = _Cache(Path(cache_dir) / f"esm_cache_{model_name.split('/')[-1]}.joblib" if cache_dir else None)
    todo = sorted({s for s in seqs if _key(s) not in cache.data}, key=len)
    if todo:
        import torch

        tok, model = _load(model_name)
        # windows: (seq, chunk)
        chunks = [(s, s[i:i + WINDOW]) for s in todo for i in range(0, len(s), WINDOW)]
        sums: dict[str, np.ndarray] = {}
        lens: dict[str, int] = {}
        done, batch = 0, []

        def run(batch):
            enc = tok([c for _, c in batch], return_tensors="pt", padding=True)
            with torch.inference_mode():
                h = model(**enc).last_hidden_state
            mask = enc["attention_mask"].clone()
            mask[:, 0] = 0  # drop <cls>
            ends = enc["attention_mask"].sum(1) - 1
            mask[torch.arange(len(batch)), ends] = 0  # drop <eos>
            vec = (h * mask.unsqueeze(-1)).sum(1).numpy()
            for (s, c), v in zip(batch, vec):
                sums[s] = sums.get(s, 0) + v
                lens[s] = lens.get(s, 0) + len(c)

        for item in chunks:
            if batch and (len(batch) + 1) * (len(item[1]) + 2) > batch_tokens:
                run(batch); done += len(batch); batch = []
                if progress:
                    progress(done / len(chunks), f"ESM-2 embeddings {done}/{len(chunks)}")
            batch.append(item)
        if batch:
            run(batch)
        for s in todo:
            cache.data[_key(s)] = (sums[s] / lens[s]).astype(np.float32)
        cache.save()
    return np.vstack([cache.data[_key(s)] for s in seqs]) if seqs else np.zeros((0, 1))
