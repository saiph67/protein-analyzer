"""Pairwise alignment (BLOSUM62, affine gaps) and k-mer similarity."""
from __future__ import annotations

import numpy as np

from .analyzer import clean_sequence

MAX_ALIGN_LENGTH = 3000  # per sequence; DP is O(n*m)

_B62_TEXT = """
   A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  *
A  4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0 -2 -1  0 -4
R -1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3 -1  0 -1 -4
N -2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3  3  0 -1 -4
D -2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3  4  1 -1 -4
C  0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1 -3 -3 -2 -4
Q -1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2  0  3 -1 -4
E -1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4
G  0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3 -1 -2 -1 -4
H -2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3  0  0 -1 -4
I -1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3 -3 -3 -1 -4
L -1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1 -4 -3 -1 -4
K -1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2  0  1 -1 -4
M -1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1 -3 -1 -1 -4
F -2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1 -3 -3 -1 -4
P -1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2 -2 -1 -2 -4
S  1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2  0  0  0 -4
T  0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0 -1 -1  0 -4
W -3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3 -4 -3 -2 -4
Y -2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1 -3 -2 -1 -4
V  0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4 -3 -2 -1 -4
B -2 -1  3  4 -3  0  1 -1  0 -3 -4  0 -3 -3 -2  0 -1 -4 -3 -3  4  1 -1 -4
Z -1  0  0  1 -3  3  4 -2  0 -3 -3  1 -1 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4
X  0 -1 -1 -1 -2 -1 -1 -1 -1 -1 -1 -1 -1 -1 -2  0  0 -2 -1 -1 -1 -1 -1 -4
* -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4  1
"""


def _load_blosum62():
    lines = [l for l in _B62_TEXT.strip("\n").split("\n")]
    cols = lines[0].split()
    mat = np.array([[int(x) for x in l.split()[1:]] for l in lines[1:]], dtype=np.int32)
    return cols, mat


_B62_LETTERS, BLOSUM62 = _load_blosum62()
_B62_INDEX = {c: i for i, c in enumerate(_B62_LETTERS)}
# Letters outside the matrix are scored as X (U->C, O->K, J->L by convention).
_ALIAS = {"U": "C", "O": "K", "J": "L"}


def _encode(seq: str) -> np.ndarray:
    return np.array([_B62_INDEX.get(_ALIAS.get(c, c), _B62_INDEX["X"]) for c in seq], dtype=np.int64)


def blosum62(a: str, b: str) -> int:
    return int(BLOSUM62[_encode(a)[0], _encode(b)[0]])


_M, _IX, _IY, _START = 0, 1, 2, 3
_NEG = -10**9


def align(a: str, b: str, mode: str = "global", gap_open: int = 11, gap_extend: int = 1) -> dict:
    """Needleman-Wunsch (global) or Smith-Waterman (local) with affine gaps (Gotoh).

    A gap of length k costs gap_open + k*gap_extend (BLAST convention, default 11/1).
    """
    if mode not in ("global", "local"):
        raise ValueError("mode must be 'global' or 'local'")
    a, b = clean_sequence(a), clean_sequence(b)
    if len(a) > MAX_ALIGN_LENGTH or len(b) > MAX_ALIGN_LENGTH:
        raise ValueError(
            f"Alignment is limited to {MAX_ALIGN_LENGTH:,} residues per sequence "
            f"(got {len(a):,} x {len(b):,}); the pure dynamic-programming aligner is O(n*m)."
        )
    local = mode == "local"
    first, ext = gap_open + gap_extend, gap_extend
    n, m = len(a), len(b)
    ea, eb = _encode(a), _encode(b)

    best_h = np.zeros((n + 1, m + 1), dtype=np.int8)   # state achieving H
    best_hp = np.zeros((n + 1, m + 1), dtype=np.int8)  # state achieving H' = max(M, Ix)
    ix_ext = np.zeros((n + 1, m + 1), dtype=bool)
    iy_ext = np.zeros((n + 1, m + 1), dtype=bool)

    j_idx = np.arange(m + 1)
    if local:
        h_prev = np.zeros(m + 1, dtype=np.int64)
        best_h[0, :] = _START
        best_h[:, 0] = _START
    else:
        h_prev = np.where(j_idx == 0, 0, -first - ext * (j_idx - 1)).astype(np.int64)
        best_h[0, 1:] = _IY
        iy_ext[0, 2:] = True
        best_h[0, 0] = _M
    ix_prev = np.full(m + 1, _NEG, dtype=np.int64)

    best_score, best_pos = (0, (0, 0))
    for i in range(1, n + 1):
        sub = BLOSUM62[ea[i - 1], eb]
        M = np.full(m + 1, _NEG, dtype=np.int64)
        M[1:] = h_prev[:-1] + sub
        ix_open, ix_extend = h_prev - first, ix_prev - ext
        Ix = np.maximum(ix_open, ix_extend)
        ix_ext[i] = ix_extend > ix_open
        if local:
            h0 = 0
        else:
            h0 = -first - ext * (i - 1)
            Ix[0] = h0
            ix_ext[i, 0] = i > 1
        M[0] = _NEG
        Hp = np.maximum(M, Ix)
        best_hp[i] = np.where(Ix > M, _IX, _M)
        # Iy[j] = max_k<j Hp[k] - first - ext*(j-1-k)  (vectorised with a running max)
        V = Hp + ext * j_idx
        cm = np.maximum.accumulate(V)
        Iy = np.full(m + 1, _NEG, dtype=np.int64)
        Iy[1:] = cm[:-1] - first - ext * (j_idx[1:] - 1)
        open_val = np.full(m + 1, _NEG, dtype=np.int64)
        open_val[1:] = Hp[:-1] - first
        iy_ext[i, 1:] = Iy[1:] > open_val[1:]
        H = np.maximum(Hp, Iy)
        state = np.where(Iy > Hp, _IY, best_hp[i]).astype(np.int8)
        if local:
            zero = H <= 0
            H = np.where(zero, 0, H)
            state[zero] = _START
            row_best = int(H.argmax())
            if H[row_best] > best_score:
                best_score, best_pos = int(H[row_best]), (i, row_best)
        else:
            H[0] = h0
            state[0] = _IX
        best_h[i] = state
        h_prev, ix_prev = H, Ix

    if local:
        score, (i, j) = best_score, best_pos
        s = int(best_h[i, j]) if score > 0 else _START
    else:
        score, i, j = int(h_prev[m]), n, m
        s = int(best_h[i, j])
    end_i, end_j = i, j

    ra, rb = [], []
    while (i > 0 or j > 0) and s != _START:
        if s == _M:
            ra.append(a[i - 1]); rb.append(b[j - 1])
            s = int(best_h[i - 1, j - 1]); i -= 1; j -= 1
            if not local and i == 0 and j == 0:
                break
        elif s == _IX:
            ra.append(a[i - 1]); rb.append("-")
            extended = ix_ext[i, j]; i -= 1
            s = _IX if extended else int(best_h[i, j])
        else:
            ra.append("-"); rb.append(b[j - 1])
            extended = iy_ext[i, j]; j -= 1
            s = _IY if extended else int(best_hp[i, j])
    aligned_a, aligned_b = "".join(reversed(ra)), "".join(reversed(rb))
    return _alignment_stats(a, b, aligned_a, aligned_b, score, mode, (i + 1, end_i), (j + 1, end_j))


def _alignment_stats(a, b, aa, ab, score, mode, range_a, range_b) -> dict:
    ident = pos = gaps = 0
    match_line = []
    for x, y in zip(aa, ab):
        if x == "-" or y == "-":
            gaps += 1
            match_line.append(" ")
        elif x == y:
            ident += 1; pos += 1
            match_line.append("|")
        elif blosum62(x, y) > 0:
            pos += 1
            match_line.append(":")
        else:
            match_line.append(".")
    L = len(aa) or 1
    res_a = sum(c != "-" for c in aa)
    res_b = sum(c != "-" for c in ab)
    return {
        "mode": mode,
        "score": score,
        "length": len(aa),
        "identity": round(100 * ident / L, 2),
        "similarity": round(100 * pos / L, 2),
        "identical": ident,
        "positives": pos,
        "gaps": gaps,
        "gap_percent": round(100 * gaps / L, 2),
        "coverage_a": round(100 * res_a / len(a), 2),
        "coverage_b": round(100 * res_b / len(b), 2),
        "range_a": list(range_a) if res_a else [0, 0],
        "range_b": list(range_b) if res_b else [0, 0],
        "aligned_a": aa,
        "aligned_b": ab,
        "match_line": "".join(match_line),
        "len_a": len(a),
        "len_b": len(b),
    }


# ---------------------------------------------------------------- k-mer similarity

_AA_INDEX = {c: i for i, c in enumerate("ACDEFGHIKLMNPQRSTVWY")}


def kmer_set(seq: str, k: int = 3) -> set[str]:
    return {seq[i:i + k] for i in range(len(seq) - k + 1)}


def kmer_similarity(a: str, b: str, k: int = 3) -> float:
    """Jaccard similarity of k-mer sets (0 = nothing shared, 1 = identical sets)."""
    sa, sb = kmer_set(a, k), kmer_set(b, k)
    if not sa or not sb:
        return 0.0
    inter = len(sa & sb)
    return inter / (len(sa) + len(sb) - inter)


def _kmer_matrix(seqs: list[str], k: int = 3):
    from scipy.sparse import csr_matrix

    rows, cols = [], []
    for r, s in enumerate(seqs):
        ids = set()
        for i in range(len(s) - k + 1):
            code = 0
            for c in s[i:i + k]:
                v = _AA_INDEX.get(c)
                if v is None:
                    code = -1
                    break
                code = code * 20 + v
            if code >= 0:
                ids.add(code)
        rows.extend([r] * len(ids))
        cols.extend(ids)
    data = np.ones(len(rows), dtype=np.float32)
    return csr_matrix((data, (rows, cols)), shape=(len(seqs), 20 ** k))


def similarity_matrix(queries: list[str], refs: list[str], k: int = 3) -> np.ndarray:
    """Dense (len(queries) x len(refs)) matrix of k-mer Jaccard similarities."""
    Q, R = _kmer_matrix(queries, k), _kmer_matrix(refs, k)
    inter = (Q @ R.T).toarray()
    qn = np.asarray(Q.sum(axis=1)).ravel()[:, None]
    rn = np.asarray(R.sum(axis=1)).ravel()[None, :]
    union = qn + rn - inter
    with np.errstate(divide="ignore", invalid="ignore"):
        sim = np.where(union > 0, inter / union, 0.0)
    return sim


def nearest_neighbor_similarity(queries: list[str], refs: list[str], k: int = 3, top: int = 1):
    """For each query, the `top` most similar references: list of [(ref_index, similarity), ...]."""
    if not refs:
        return [[] for _ in queries]
    out = []
    # chunk queries to bound memory
    for start in range(0, len(queries), 500):
        sim = similarity_matrix(queries[start:start + 500], refs, k)
        for row in sim:
            idx = np.argsort(-row)[:top]
            out.append([(int(i), float(row[i])) for i in idx])
    return out
