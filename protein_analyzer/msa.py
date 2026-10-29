"""Multiple sequence alignment (center-star) and a neighbour-joining family tree.

Center-star: align every pair, pick the sequence most similar to all others as the "center", then merge
the center's pairwise alignments ("once a gap, always a gap"). Simple and exact enough for teaching-sized
sets (<= 20 proteins); dedicated tools (Clustal Omega, MAFFT) are better for large or distant families.
"""
from __future__ import annotations

import math
from collections import Counter

import numpy as np

from .analyzer import clean_sequence
from .compare import align

MAX_SEQS = 20
MAX_LEN = 1500


def pairwise(seqs: list[str], progress=None):
    n = len(seqs)
    scores = np.zeros((n, n))
    dist = np.zeros((n, n))
    ident = np.eye(n) * 100
    alns = {}
    total, done = n * (n - 1) // 2, 0
    for i in range(n):
        scores[i, i] = align(seqs[i], seqs[i])["score"]
        for j in range(i + 1, n):
            a = align(seqs[i], seqs[j], mode="global")
            alns[(i, j)] = (a["aligned_a"], a["aligned_b"])
            cols = sum(1 for x, y in zip(a["aligned_a"], a["aligned_b"]) if x != "-" and y != "-")
            p = 1 - (a["identical"] / cols) if cols else 1.0
            ident[i, j] = ident[j, i] = round(100 * (1 - p), 1)
            dist[i, j] = dist[j, i] = _kimura(p)
            scores[i, j] = scores[j, i] = a["score"]
            done += 1
            if progress:
                progress(done / max(total, 1), f"Aligned {done}/{total} pairs")
    return scores, dist, ident, alns


def _kimura(p: float) -> float:
    """Kimura's protein distance correction for multiple substitutions (capped for very distant pairs)."""
    x = 1 - p - 0.2 * p * p
    return round(-math.log(x), 4) if x > 0.05 else 3.0


def center_star(seqs: list[str], scores, alns) -> tuple[list[str], int]:
    n = len(seqs)
    c = int(np.argmax(scores.sum(axis=1) - scores.diagonal()))
    rows = {c: seqs[c]}
    for i in range(n):
        if i == c:
            continue
        pc, pi = alns[(c, i)] if c < i else alns[(i, c)][::-1]
        center_row = rows[c]
        x = y = 0
        new_rows = {k: [] for k in rows}
        new_i = []
        while x < len(center_row) or y < len(pc):
            mx = center_row[x] if x < len(center_row) else None
            py = pc[y] if y < len(pc) else None
            if mx == "-" and py == "-":
                for k in rows: new_rows[k].append(rows[k][x])
                new_i.append(pi[y]); x += 1; y += 1
            elif mx == "-" or py is None:
                for k in rows: new_rows[k].append(rows[k][x])
                new_i.append("-"); x += 1
            elif py == "-" or mx is None:
                for k in rows: new_rows[k].append("-")
                new_i.append(pi[y]); y += 1
            else:
                for k in rows: new_rows[k].append(rows[k][x])
                new_i.append(pi[y]); x += 1; y += 1
        rows = {k: "".join(v) for k, v in new_rows.items()}
        rows[i] = "".join(new_i)
    return [rows[i] for i in range(n)], c


def column_stats(rows: list[str]) -> tuple[str, list[float]]:
    consensus, conservation = [], []
    for col in zip(*rows):
        res = [r for r in col if r != "-"]
        if not res:
            consensus.append("-"); conservation.append(0.0); continue
        aa, k = Counter(res).most_common(1)[0]
        frac = k / len(col)
        conservation.append(round(frac, 3))
        consensus.append(aa if frac >= 0.5 else ("+" if frac >= 0.3 else " "))
    return "".join(consensus), conservation


def neighbor_joining(dist: np.ndarray, names: list[str]) -> dict:
    """Neighbour-joining tree, midpoint-rooted. Returns nested {'name', 'length', 'children'}."""
    n = len(names)
    if n == 1:
        return {"name": names[0], "length": 0.0, "children": []}
    adj: dict[int, list[tuple[int, float]]] = {i: [] for i in range(n)}
    D = {(i, j): float(dist[i, j]) for i in range(n) for j in range(n)}
    active, nxt = list(range(n)), n

    def link(a, b, w):
        adj.setdefault(a, []).append((b, max(w, 0.0)))
        adj.setdefault(b, []).append((a, max(w, 0.0)))

    while len(active) > 2:
        r = len(active)
        tot = {i: sum(D[(i, k)] for k in active) for i in active}
        best, pair = None, None
        for x in range(r):
            for y in range(x + 1, r):
                i, j = active[x], active[y]
                q = (r - 2) * D[(i, j)] - tot[i] - tot[j]
                if best is None or q < best:
                    best, pair = q, (i, j)
        i, j = pair
        li = 0.5 * D[(i, j)] + (tot[i] - tot[j]) / (2 * (r - 2))
        lj = D[(i, j)] - li
        u = nxt; nxt += 1
        link(u, i, li); link(u, j, lj)
        for k in active:
            if k not in (i, j):
                D[(u, k)] = D[(k, u)] = 0.5 * (D[(i, k)] + D[(j, k)] - D[(i, j)])
        D[(u, u)] = 0.0
        active = [k for k in active if k not in (i, j)] + [u]
    a, b = active
    link(a, b, D[(a, b)])
    return _midpoint_root(adj, names)


def _midpoint_root(adj, names):
    leaves = list(range(len(names)))

    def far(start):
        dist, parent, stack = {start: 0.0}, {start: None}, [start]
        while stack:
            v = stack.pop()
            for w, l in adj[v]:
                if w not in dist:
                    dist[w], parent[w] = dist[v] + l, v
                    stack.append(w)
        leaf = max(leaves, key=lambda x: dist[x])
        return leaf, dist, parent

    a, _, _ = far(leaves[0])
    b, dist, parent = far(a)
    half, v = dist[b] / 2, b
    while parent[v] is not None and dist[parent[v]] > half:
        v = parent[v]
    p = parent[v]
    if p is None:  # degenerate (all distances zero)
        root_children = [(w, l) for w, l in adj[a]]
        return {"name": None, "length": 0.0, "children": [_subtree(adj, w, a, l, names) for w, l in root_children]
                + [{"name": names[a], "length": 0.0, "children": []}]}
    edge = next(l for w, l in adj[v] if w == p)
    to_p = half - dist[p]  # distance from p toward v where the root sits
    return {"name": None, "length": 0.0, "children": [
        _subtree(adj, v, p, edge - to_p, names), _subtree(adj, p, v, to_p, names)]}


def _subtree(adj, v, came_from, length, names):
    kids = [_subtree(adj, w, v, l, names) for w, l in adj[v] if w != came_from]
    return {"name": names[v] if v < len(names) else None, "length": round(max(length, 0.0), 4), "children": kids}


def to_newick(node) -> str:
    def rec(nd):
        label = (nd["name"] or "").replace("(", "_").replace(")", "_").replace(",", "_").replace(":", "_").replace(";", "_").replace(" ", "_")
        if nd["children"]:
            label = "(" + ",".join(rec(c) for c in nd["children"]) + ")" + label
        return f"{label}:{nd['length']}"
    return rec(node).rsplit(":", 1)[0] + ";"


def build(records: list[tuple[str, str]], progress=None) -> dict:
    if len(records) < 3:
        raise ValueError("Enter at least 3 proteins (use Compare for two).")
    if len(records) > MAX_SEQS:
        raise ValueError(f"At most {MAX_SEQS} proteins at a time.")
    names, seqs = [], []
    for n, s in records:
        s = clean_sequence(s)
        if len(s) > MAX_LEN:
            raise ValueError(f"{n} is {len(s):,} aa; the multi-protein aligner allows up to {MAX_LEN:,} per protein.")
        names.append(n); seqs.append(s)
    seen = Counter(names)
    names = [f"{n} ({i + 1})" if seen[n] > 1 else n for i, n in enumerate(names)]
    scores, dist, ident, alns = pairwise(seqs, progress)
    rows, center = center_star(seqs, scores, alns)
    consensus, conservation = column_stats(rows)
    tree = neighbor_joining(dist, names)
    return {"names": names, "rows": rows, "center": center, "consensus": consensus,
            "conservation": conservation, "identity_matrix": ident.round(1).tolist(),
            "distance_matrix": dist.round(4).tolist(), "tree": tree, "newick": to_newick(tree),
            "length": len(rows[0]), "lengths": [len(s) for s in seqs]}
