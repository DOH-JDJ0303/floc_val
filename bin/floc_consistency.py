#!/usr/bin/env python3
"""
floc_consistency.py - helper for the floc-consistency Nextflow pipeline.
summarise uses the standard library only; tree uses numpy and plotly. Both
run in the floc container, which provides them.

  summarise  combine the final floc clusters of every batch-size run with an
             all-vs-all SNP distance matrix and write intra-/inter-cluster
             distance statistics, discrepant pairs and close pairs that floc
             placed in different clusters
  tree       neighbour-joining tree of the SNP distances (midpoint rooted, Newick)
             and an interactive plotly HTML figure: the tree with one column of
             tiles per batch size showing each sample's cluster

A run's clusters are given as files named bs<BATCH_SIZE>.clusters.csv (floc's
clusters.csv: columns sid, cluster). The reference is the run whose batch
size equals the number of samples (all samples in one batch).

The distance matrix can be
  - long form: sample_a, sample_b, distance[, ...] - e.g. `ska distance`
    output (Sample1 Sample2 Distance Mismatches Match_count ...). A column
    named Distance / snps / snp_distance is used if present, else column 3.
  - square: a header row of sample names (first cell is ignored) and one row
    per sample, first cell the sample name.
Tab, comma or whitespace separated, optionally gzipped.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import heapq
import math
import re
import sys
from array import array
from collections import Counter, defaultdict
from itertools import combinations
from pathlib import Path

STATS = ["n", "min", "max", "mean", "sd", "median"]
NAN = float("nan")
NEAREST_K = 64          # neighbours kept per sample for the nearest-other-cluster search


# --------------------------------------------------------------------------- #
# Reading clusters
# --------------------------------------------------------------------------- #
def read_runs(files: list[str]) -> dict[int, dict[str, str]]:
    runs = {}
    for f in files:
        m = re.search(r"bs(\d+)\.clusters\.csv$", Path(f).name)
        if not m:
            sys.exit(f"ERROR: cannot read the batch size from file name {f}")
        with open(f, newline="") as fh:
            rows = list(csv.DictReader(fh))
        if rows and not {"sid", "cluster"} <= set(rows[0]):
            sys.exit(f"ERROR: {f} lacks columns sid, cluster")
        runs[int(m.group(1))] = {r["sid"].strip(): r["cluster"].strip() for r in rows
                                 if r.get("sid") and r.get("cluster")}
    if not runs:
        sys.exit("ERROR: no cluster files given")
    return dict(sorted(runs.items()))


def reference_size(runs: dict[int, dict], n: int) -> int:
    if n not in runs:
        sys.exit(f"ERROR: the single-batch reference run (batch size {n}) is missing")
    return n


def members(part: dict[str, str]) -> dict[str, list[str]]:
    out = defaultdict(list)
    for s, c in part.items():
        out[c].append(s)
    return {c: sorted(m) for c, m in out.items()}


def pairs_of(mem: list[str]):
    return combinations(sorted(mem), 2)


# --------------------------------------------------------------------------- #
# Distance matrix (condensed upper triangle, NaN = missing)
# --------------------------------------------------------------------------- #
class Distances:
    def __init__(self, names: list[str]):
        self.names = names
        self.idx = {s: i for i, s in enumerate(names)}
        n = len(names)
        self.n = n
        self.d = array("d", [NAN]) * (n * (n - 1) // 2)
        self.sites = array("d", [NAN]) * (n * (n - 1) // 2)
        self.conflicts = 0

    def _k(self, i: int, j: int) -> int:
        if i > j:
            i, j = j, i
        return i * self.n - i * (i + 1) // 2 + (j - i - 1)

    def set(self, i: int, j: int, d: float, sites: float = NAN):
        k = self._k(i, j)
        old = self.d[k]
        if not math.isnan(old) and old != d:
            self.conflicts += 1
            d = min(old, d)
        self.d[k] = d
        if not math.isnan(sites):
            self.sites[k] = sites

    def get(self, a: str, b: str) -> tuple[float | None, float | None]:
        i, j = self.idx.get(a), self.idx.get(b)
        if i is None or j is None or i == j:
            return None, None
        k = self._k(i, j)
        d, s = self.d[k], self.sites[k]
        return (None if math.isnan(d) else d), (None if math.isnan(s) else s)

    def row(self, i: int):
        """(j, d) for every j != i with a distance."""
        n, d, k = self.n, self.d, self._k
        for j in range(n):
            if j != i:
                v = d[k(i, j)]
                if not math.isnan(v):
                    yield j, v


def _open(path: str):
    return gzip.open(path, "rt", newline="") if path.endswith(".gz") else open(path, newline="")


def _num(x: str) -> float | None:
    try:
        v = float(x)
    except ValueError:
        return None
    return v if not math.isnan(v) else None


def read_distances(path: str, keep: set[str] | None = None) -> Distances:
    """Read the matrix; with `keep`, samples outside it are dropped."""
    with _open(path) as fh:
        first = fh.readline()
        if "\t" in first:
            split = lambda line: line.rstrip("\r\n").split("\t")
        elif "," in first:
            split = lambda line: next(csv.reader([line.rstrip("\r\n")]))
        else:
            split = lambda line: line.split()
        rows = [split(first)] + [split(line) for line in fh if line.strip()]
    rows = [[c.strip() for c in r] for r in rows]
    if len(rows) < 2:
        sys.exit(f"ERROR: distance matrix {path} has fewer than two lines")
    header, body = rows[0], rows[1:]

    # square: header names == first column names
    hnames = header[1:]
    rnames = [r[0] for r in body]
    if len(hnames) == len(rnames) and set(hnames) == set(rnames) and len(set(rnames)) == len(rnames) \
            and all(len(r) == len(header) for r in body):
        dist = Distances(sorted(s for s in rnames if keep is None or s in keep))
        col = [dist.idx.get(s) for s in hnames]
        for r in body:
            i = dist.idx.get(r[0])
            if i is None:
                continue
            for c, v in zip(col, r[1:]):
                if c is not None and c != i:
                    x = _num(v)
                    if x is not None:
                        dist.set(i, c, x)
        fmt_name = "square"
    else:
        # long form
        has_header = len(header) >= 3 and _num(header[2]) is None
        lower = [h.lower().replace(" ", "_") for h in header] if has_header else []
        dcol = next((lower.index(c) for c in ("distance", "snps", "snp_distance", "snp")
                     if c in lower), 2)
        scol = next((lower.index(c) for c in ("match_count", "core_sites", "sites") if c in lower), None)
        data = body if has_header else rows
        if any(len(r) <= dcol for r in data):
            sys.exit(f"ERROR: distance matrix {path}: expected at least {dcol + 1} columns on every "
                     "line (sample_a, sample_b, distance) or a square matrix with a header row")
        if keep is not None:
            data = [r for r in data if r[0] in keep and r[1] in keep]
        dist = Distances(sorted({r[0] for r in data} | {r[1] for r in data}))
        bad = 0
        for r in data:
            x = _num(r[dcol])
            if x is None:
                bad += 1
                continue
            if r[0] == r[1]:
                continue
            s = _num(r[scol]) if scol is not None and len(r) > scol else None
            dist.set(dist.idx[r[0]], dist.idx[r[1]], x, NAN if s is None else s)
        if bad:
            print(f"WARNING: {bad} lines in {path} had no numeric distance", file=sys.stderr)
        fmt_name = "long"
    if dist.conflicts:
        print(f"WARNING: {dist.conflicts} pairs had two different distances in {path}; "
              "the smaller was used", file=sys.stderr)
    print(f"distance matrix: {path} ({fmt_name} form, {dist.n} samples"
          f"{' after keeping clustered samples' if keep is not None else ''})")
    return dist


# --------------------------------------------------------------------------- #
# Pairs
# --------------------------------------------------------------------------- #
def discrepant_pairs(runs: dict[int, dict], ref_bs: int) -> list[dict]:
    """Pairs grouped differently from the reference in at least one run."""
    ref = runs[ref_bs]
    together = defaultdict(set)                    # pair -> batch sizes where together
    for bs, part in runs.items():
        for mem in members(part).values():
            for p in pairs_of(mem):
                together[p].add(bs)
    ref_pairs = {p for p, b in together.items() if ref_bs in b}
    rows = []
    for p in sorted(together):
        a, b = p
        if a not in ref or b not in ref:
            continue
        both = [bs for bs, part in runs.items() if a in part and b in part]
        in_ref = p in ref_pairs
        differing = [bs for bs in both if (bs in together[p]) != in_ref]
        if not differing:
            continue
        rows.append({
            "sample_a": a, "sample_b": b,
            "reference_cluster_a": ref[a], "reference_cluster_b": ref[b],
            "pair_class": "together_in_reference" if in_ref else "apart_in_reference",
            "runs_together": len(together[p] & set(both)), "runs_with_both": len(both),
            "runs_different_from_reference": len(differing),
            "batch_sizes_different_from_reference": ",".join(map(str, sorted(differing))),
        })
    return rows


# --------------------------------------------------------------------------- #
# Statistics and output
# --------------------------------------------------------------------------- #
def describe(values) -> dict:
    return describe_counts(Counter(float(v) for v in values if v is not None))


def describe_counts(c: Counter) -> dict:
    """n, min, max, mean, sd, median of a value -> count histogram."""
    c = +c                                          # drop zero / negative counts
    n = sum(c.values())
    if not n:
        return {"n": 0, "min": "", "max": "", "mean": "", "sd": "", "median": ""}
    keys = sorted(c)
    mean = sum(k * v for k, v in c.items()) / n
    sd = math.sqrt(sum(v * (k - mean) ** 2 for k, v in c.items()) / (n - 1)) if n > 1 else ""
    lo, hi = (n - 1) // 2, n // 2                   # median positions (0-based)
    seen, mlo, mhi = 0, None, None
    for k in keys:
        seen += c[k]
        if mlo is None and seen > lo:
            mlo = k
        if seen > hi:
            mhi = k
            break
    return {"n": n, "min": keys[0], "max": keys[-1], "mean": mean, "sd": sd,
            "median": (mlo + mhi) / 2}


def fmt(v):
    if isinstance(v, bool):
        return str(v)
    if isinstance(v, float):
        return str(int(v)) if v.is_integer() and abs(v) < 1e15 else f"{v:.6g}"
    if v is None:
        return ""
    return str(v)


def write_tsv(path: str, rows: list[dict], cols: list[str], header_lines: list[str] = ()):
    with open(path, "w", newline="") as fh:
        for line in header_lines:
            fh.write(f"# {line}\n")
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(fmt(r.get(c, "")) for c in cols) + "\n")


# --------------------------------------------------------------------------- #
# summarise
# --------------------------------------------------------------------------- #
def cmd_summarise(a):
    runs = read_runs(a.clusters)
    n = a.n_samples
    ref_bs = reference_size(runs, n)
    ref = runs[ref_bs]
    clustered = set().union(*runs.values())
    dist = read_distances(a.distances, keep=clustered)    # ignore samples not in this analysis
    thr = a.snp_threshold
    below = f"snps_below_{thr:g}"

    not_in_matrix = sorted(clustered - set(dist.idx))
    if not_in_matrix:
        print(f"WARNING: {len(not_in_matrix)} clustered samples are not in the distance matrix: "
              f"{', '.join(not_in_matrix[:20])}{' ...' if len(not_in_matrix) > 20 else ''}",
              file=sys.stderr)

    # one pass over the matrix: histogram of every pair, close pairs, nearest neighbours
    all_counts: Counter = Counter()
    close: list[tuple[int, int, float]] = []
    nearest: list[list[tuple[float, int]]] = [[] for _ in range(dist.n)]
    names, N, D = dist.names, dist.n, dist.d
    k = 0
    for i in range(N):
        for j in range(i + 1, N):
            v = D[k]
            k += 1
            if math.isnan(v):
                continue
            all_counts[v] += 1
            if v < thr:
                close.append((i, j, v))
            for x, y in ((i, j), (j, i)):
                h = nearest[x]
                if len(h) < NEAREST_K:
                    heapq.heappush(h, (-v, y))
                elif -h[0][0] > v:
                    heapq.heapreplace(h, (-v, y))
    nearest_sorted = [sorted((-d, y) for d, y in h) for h in nearest]
    pairs_missing_total = N * (N - 1) // 2 - sum(all_counts.values())

    def nearest_other(i: int, part_idx: dict[int, str], own: str):
        """Nearest sample (index, distance) outside cluster `own`, among clustered samples."""
        for d, y in nearest_sorted[i]:
            if y in part_idx and part_idx[y] != own:
                return d, y
        if len(nearest_sorted[i]) < NEAREST_K:
            return None, None                       # list holds every neighbour with a distance
        best = (None, None)
        for y, d in dist.row(i):
            if y in part_idx and part_idx[y] != own and (best[0] is None or d < best[0]):
                best = (d, y)
        return best

    # discrepant pairs
    disc = discrepant_pairs(runs, ref_bs)
    for r in disc:
        r["snps"], r["sites"] = dist.get(r["sample_a"], r["sample_b"])
        r[below] = r["snps"] is not None and r["snps"] < thr

    order = [ref_bs] + [bs for bs in runs if bs != ref_bs]
    within, cstats, split_rows, metrics, summary = [], [], [], [], []
    for bs in order:
        part = runs[bs]
        part_idx = {dist.idx[s]: c for s, c in part.items() if s in dist.idx}
        is_ref = bs == ref_bs

        # every pair among this run's clustered samples (subtract pairs with unclustered samples)
        run_counts = all_counts.copy()
        outside = sorted(set(range(N)) - set(part_idx))
        out_set = set(outside)
        for i in outside:
            for j, v in dist.row(i):
                if j not in out_set or j > i:
                    run_counts[v] -= 1
        run_pairs_total = len(part_idx) * (len(part_idx) - 1) // 2

        intra_counts: Counter = Counter()
        intra_missing = 0
        mem_all = members(part)
        for c, mem in sorted(mem_all.items(), key=lambda kv: (-len(kv[1]), kv[0])):
            vals = []
            for x, y in pairs_of(mem):
                s, sites = dist.get(x, y)
                if s is None:
                    intra_missing += 1
                else:
                    intra_counts[s] += 1
                vals.append(s)
                within.append({"batch_size": bs, "is_reference": is_ref, "cluster": c,
                               "cluster_size": len(mem), "sample_a": x, "sample_b": y,
                               "snps": s, "sites": sites})
            d = describe(vals)
            # nearest sample in another cluster (inter-cluster separation)
            best_d, best_s = None, None
            for s in mem:
                i = dist.idx.get(s)
                if i is None:
                    continue
                dd, y = nearest_other(i, part_idx, c)
                if dd is not None and (best_d is None or dd < best_d):
                    best_d, best_s = dd, (s, names[y])
            cstats.append({"batch_size": bs, "is_reference": is_ref, "cluster": c,
                           "cluster_size": len(mem), "n_pairs": d["n"],
                           **{f"intra_snps_{k}": d[k] for k in STATS if k != "n"},
                           "nearest_other_cluster_snps": best_d,
                           "nearest_other_cluster": part.get(best_s[1], "") if best_s else "",
                           "nearest_other_pair": f"{best_s[0]}:{best_s[1]}" if best_s else "",
                           "members": ",".join(mem)})

        inter_counts = run_counts - intra_counts
        intra_d, inter_d = describe_counts(intra_counts), describe_counts(inter_counts)

        # close pairs placed in different clusters
        n_close, split = 0, []
        for i, j, v in close:
            if i in part_idx and j in part_idx:
                n_close += 1
                if part_idx[i] != part_idx[j]:
                    a_, b_ = names[i], names[j]
                    split.append({"batch_size": bs, "is_reference": is_ref,
                                  "sample_a": a_, "sample_b": b_, "snps": v,
                                  "sites": dist.get(a_, b_)[1],
                                  "cluster_a": part_idx[i], "cluster_b": part_idx[j],
                                  "split_in_reference": (a_ in ref and b_ in ref and ref[a_] != ref[b_])
                                  if not is_ref else True})
        split.sort(key=lambda r: (r["snps"], r["sample_a"], r["sample_b"]))
        split_rows += split
        split_samples = {s for r in split for s in (r["sample_a"], r["sample_b"])}

        sizes = [len(m) for m in mem_all.values()]
        metrics.append({
            "batch_size": bs, "is_reference": is_ref, "n_batches": math.ceil(n / bs),
            "samples_clustered": len(part), "samples_with_distances": len(part_idx),
            "clusters": len(sizes), "clusters_2plus": sum(s > 1 for s in sizes),
            "singletons": sum(s == 1 for s in sizes),
            **{f"intra_{k}": intra_d[k] for k in ("n", "min", "median", "mean", "max")},
            **{f"inter_{k}": inter_d[k] for k in ("n", "min", "median", "mean", "max")},
            f"pairs_{below}": n_close, f"pairs_{below}_split": len(split),
            f"pairs_{below}_split_frac": len(split) / n_close if n_close else "",
            "samples_in_split_pairs": len(split_samples),
            "pairs_missing_distance": run_pairs_total - sum((+run_counts).values()),
        })
        label = f"batch size {bs}" + (" (reference)" if is_ref else "")
        summary.append({"section": "intra-cluster pairs", "group": label, **intra_d})
        summary.append({"section": "inter-cluster pairs", "group": label, **inter_d})

    for label, sel in (("all", disc),
                       ("together_in_reference", [r for r in disc if r["pair_class"] == "together_in_reference"]),
                       ("apart_in_reference", [r for r in disc if r["pair_class"] == "apart_in_reference"])):
        summary.append({"section": "discrepant pairs", "group": label,
                        **describe(r["snps"] for r in sel)})

    # header / console report
    close_disc = [r for r in disc if r[below]]
    header = [f"reference: batch size {ref_bs} (all {n} samples in one batch); runs: {len(runs)}; "
              f"distances: {a.distance_source}",
              f"distance matrix: {dist.n} samples, {pairs_missing_total} pairs without a distance",
              f"pairs with < {thr:g} SNPs: {len(close)}",
              f"discrepant pairs (grouped differently from the reference): {len(disc)}, "
              f"{len(close_disc)} with < {thr:g} SNPs"]
    for m in metrics:
        tag = " (reference)" if m["is_reference"] else ""
        header.append(
            f"batch size {m['batch_size']}{tag}: {m['clusters']} clusters "
            f"({m['clusters_2plus']} with 2+ samples, {m['singletons']} singletons); "
            f"intra median {fmt(m['intra_median'])}, max {fmt(m['intra_max'])}; "
            f"inter min {fmt(m['inter_min'])}, median {fmt(m['inter_median'])}; "
            f"{m[f'pairs_{below}_split']} of {m[f'pairs_{below}']} pairs with < {thr:g} SNPs "
            f"in different clusters ({m['samples_in_split_pairs']} samples)")

    disc_cols = ["sample_a", "sample_b", "reference_cluster_a", "reference_cluster_b",
                 "pair_class", "runs_together", "runs_with_both", "runs_different_from_reference",
                 "batch_sizes_different_from_reference", "snps", "sites", below]
    disc.sort(key=lambda r: (r["snps"] is None, r["snps"] or 0, r["sample_a"], r["sample_b"]))
    write_tsv("discrepant_pairs.tsv", disc, disc_cols)
    write_tsv("within_cluster_pairs.tsv", within,
              ["batch_size", "is_reference", "cluster", "cluster_size", "sample_a", "sample_b",
               "snps", "sites"])
    write_tsv("cluster_stats.tsv", cstats,
              ["batch_size", "is_reference", "cluster", "cluster_size", "n_pairs"]
              + [f"intra_snps_{k}" for k in STATS if k != "n"]
              + ["nearest_other_cluster_snps", "nearest_other_cluster", "nearest_other_pair", "members"])
    write_tsv("split_close_pairs.tsv", split_rows,
              ["batch_size", "is_reference", "sample_a", "sample_b", "snps", "sites",
               "cluster_a", "cluster_b", "split_in_reference"],
              [f"pairs with < {thr:g} SNPs that floc placed in different clusters, per batch size"])
    write_tsv("run_metrics.tsv", metrics, list(metrics[0]))
    write_tsv("summary.tsv", summary, ["section", "group"] + STATS, header)

    print("\n".join(header))
    print("\t".join(["section", "group"] + STATS))
    for r in summary:
        print("\t".join(fmt(r[c]) for c in ["section", "group"] + STATS))


# --------------------------------------------------------------------------- #
# tree: neighbour-joining tree + cluster tiles (needs numpy and plotly, which
# the floc container provides)
# --------------------------------------------------------------------------- #
PALETTE = [  # qualitative, high contrast; cycled when there are more clusters
    "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2",
    "#000000", "#bcbd22", "#17becf", "#393b79", "#e6550d", "#31a354", "#843c39",
    "#7b4173", "#3182bd", "#fd8d3c", "#74c476", "#e7969c", "#a55194", "#6baed6",
    "#fdae6b", "#a1d99b", "#ad494a", "#de9ed6", "#637939", "#8ca252", "#bd9e39",
    "#e7ba52", "#d6616b", "#ce6dbd", "#9c9ede", "#5254a3", "#b5cf6b", "#cedb9c",
    "#e7cb94", "#c7e9c0", "#9ecae1", "#fdd0a2", "#9e0142",
]
SINGLETON = "#d9d9d9"


def neighbour_joining(D, np):
    """Neighbour joining on a full symmetric matrix. Returns (edges, n_nodes):
    edges are (u, v, length) with leaves 0..n-1; negative lengths are set to 0."""
    n = D.shape[0]
    if n == 1:
        return [], 1
    if n == 2:
        return [(0, 1, max(float(D[0, 1]), 0.0))], 2
    D = D.astype(float).copy()
    ids = list(range(n))
    edges, nxt = [], n
    while len(ids) > 2:
        m = len(ids)
        r = D.sum(axis=1)
        Q = (m - 2) * D - r[:, None] - r[None, :]
        np.fill_diagonal(Q, np.inf)
        i, j = divmod(int(np.argmin(Q)), m)
        if i > j:
            i, j = j, i
        li = 0.5 * D[i, j] + (r[i] - r[j]) / (2 * (m - 2))
        lj = D[i, j] - li
        u = nxt
        nxt += 1
        edges += [(u, ids[i], max(float(li), 0.0)), (u, ids[j], max(float(lj), 0.0))]
        du = 0.5 * (D[i] + D[j] - D[i, j])
        D[i, :] = du
        D[:, i] = du
        D[i, i] = 0.0
        D = np.delete(np.delete(D, j, axis=0), j, axis=1)
        ids[i] = u
        del ids[j]
    edges.append((ids[0], ids[1], max(float(D[0, 1]), 0.0)))
    return edges, nxt


def midpoint_root(edges, n_nodes, n_leaves):
    """Root the unrooted tree at the midpoint of its longest leaf-to-leaf path.
    Returns (root, children: dict node -> [(child, length)])."""
    adj = defaultdict(list)
    for u, v, w in edges:
        adj[u].append((v, w))
        adj[v].append((u, w))
    if n_leaves == 1:
        return 0, {0: []}

    def farthest(src):
        dist, prev, stack = {src: 0.0}, {src: None}, [src]
        while stack:
            x = stack.pop()
            for y, w in adj[x]:
                if y not in dist:
                    dist[y], prev[y] = dist[x] + w, x
                    stack.append(y)
        leaf = max((x for x in dist if x < n_leaves), key=lambda x: (dist[x], -x))
        return leaf, dist, prev

    a, _, _ = farthest(0)
    b, dist, prev = farthest(a)
    half = dist[b] / 2
    # walk from b towards a until the midpoint edge is found
    x = b
    while prev[x] is not None and dist[prev[x]] > half:
        x = prev[x]
    p = prev[x]
    root = n_nodes                                  # new root node on edge p-x
    if p is None:
        root_children = [(y, w) for y, w in adj[x]]
        root = x
    else:
        w = next(w for y, w in adj[x] if y == p)
        to_x = dist[x] - half                       # root -> x
        adj[x] = [(y, ww) for y, ww in adj[x] if y != p] + [(root, to_x)]
        adj[p] = [(y, ww) for y, ww in adj[p] if y != x] + [(root, w - to_x)]
        adj[root] = [(x, to_x), (p, w - to_x)]
    children, seen, stack = defaultdict(list), {root}, [root]
    while stack:
        u = stack.pop()
        for v, w in adj[u]:
            if v not in seen:
                seen.add(v)
                children[u].append((v, w))
                stack.append(v)
    return root, children


def layout(root, children, n_leaves):
    """x = distance from root, y = leaf order (ladderized); internal y = mean of
    first and last child. Iterative, so deep trees do not hit recursion limits."""
    size = {}
    for u in reversed(_preorder(root, children)):
        size[u] = 1 if u < n_leaves and not children.get(u) else sum(size[v] for v, _ in children[u])
    for u in children:
        children[u].sort(key=lambda c: (size[c[0]], c[0]))   # smaller clades first
    x, y, order = {root: 0.0}, {}, []
    for u in _preorder(root, children):
        for v, w in children.get(u, []):
            x[v] = x[u] + w
    for u in _preorder(root, children):
        if not children.get(u):
            y[u] = len(order)
            order.append(u)
    for u in reversed(_preorder(root, children)):
        if children.get(u):
            ys = [y[v] for v, _ in children[u]]
            y[u] = (min(ys) + max(ys)) / 2
    return x, y, order


def _preorder(root, children):
    out, stack = [], [root]
    while stack:
        u = stack.pop()
        out.append(u)
        stack.extend(v for v, _ in reversed(children.get(u, [])))
    return out


def newick(root, children, names, n_leaves) -> str:
    def label(i):
        s = names[i]
        return s if re.fullmatch(r"[A-Za-z0-9_.\-]+", s) else "'" + s.replace("'", "''") + "'"
    # iterative post-order build
    text = {}
    for u in reversed(_preorder(root, children)):
        if not children.get(u):
            text[u] = label(u) if u < n_leaves else ""
        else:
            text[u] = "(" + ",".join(f"{text[v]}:{w:.6g}" for v, w in children[u]) + ")"
    return text[root] + ";\n"


def cluster_colours(runs, ref_bs, order_names):
    """colour per (batch size, sample). Reference clusters get palette colours in
    tree order; another run's cluster takes the colour of the reference cluster it
    overlaps most (one run cluster per reference cluster, largest overlap first);
    the rest get new colours. Singletons are grey."""
    ref = runs[ref_bs]
    pal = iter(PALETTE * (1 + 4 * len(order_names) // len(PALETTE)))
    ref_mem = members(ref)
    ref_col = {}
    for s in order_names:                           # colours follow the tree top-down
        c = ref.get(s)
        if c is not None and c not in ref_col:
            ref_col[c] = SINGLETON if len(ref_mem[c]) == 1 else next(pal)
    out = {}
    for bs, part in runs.items():
        mem = members(part)
        if bs == ref_bs:
            col = dict(ref_col)
        else:
            overlap = Counter((c, ref[s]) for s, c in part.items() if s in ref)
            col, used = {}, set()
            for (c, rc), k in sorted(overlap.items(), key=lambda kv: (-kv[1], kv[0])):
                if c not in col and rc not in used and len(mem[c]) > 1 and ref_col[rc] != SINGLETON:
                    col[c] = ref_col[rc]
                    used.add(rc)
            for s in order_names:
                c = part.get(s)
                if c is not None and c not in col:
                    col[c] = SINGLETON if len(mem[c]) == 1 else next(pal)
        for s, c in part.items():
            out[bs, s] = col[c]
    return out


def cmd_tree(a):
    import numpy as np
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots

    runs = read_runs(a.clusters)
    ref_bs = reference_size(runs, a.n_samples)
    clustered = set().union(*runs.values())
    dist = read_distances(a.distances, keep=clustered)
    names, n = dist.names, dist.n
    if n == 0:
        sys.exit("ERROR: no clustered sample is in the distance matrix")
    thr = a.snp_threshold

    D = np.zeros((n, n))
    iu = np.triu_indices(n, 1)
    D[iu] = np.frombuffer(dist.d, dtype=float) if n > 1 else []
    D = D + D.T
    missing = int(np.isnan(D[iu]).sum()) if n > 1 else 0
    if missing:
        fill = float(np.nanmax(D)) if np.isfinite(np.nanmax(D)) else 0.0
        print(f"WARNING: {missing} pairs have no distance; set to the largest distance ({fill:g}) "
              "for the tree", file=sys.stderr)
        D = np.where(np.isnan(D), fill, D)
        np.fill_diagonal(D, 0.0)

    edges, n_nodes = neighbour_joining(D, np)
    root, children = midpoint_root(edges, n_nodes, n)
    x, y, order = layout(root, children, n)
    order_names = [names[i] for i in order]
    with open(a.newick, "w") as fh:
        fh.write(newick(root, children, names, n))

    # ---- tree traces: elbow segments, with None breaks between segments
    hx, hy, vx, vy = [], [], [], []
    for u, kids in children.items():
        if not kids:
            continue
        ys = [y[v] for v, _ in kids]
        vx += [x[u], x[u], None]
        vy += [min(ys), max(ys), None]
        for v, _ in kids:
            hx += [x[u], x[v], None]
            hy += [y[v], y[v], None]
    xmax = max(x.values()) or 1.0
    gx, gy = [], []                                 # faint guides from tips to the tiles
    for i in order:
        gx += [x[i], xmax, None]
        gy += [y[i], y[i], None]

    # ---- tiles: one column per batch size, reference first
    cols = [ref_bs] + [bs for bs in runs if bs != ref_bs]
    col_labels = [f"{bs} (ref)" if bs == ref_bs else str(bs) for bs in cols]
    colour = cluster_colours(runs, ref_bs, order_names)
    palette = sorted(set(colour.values()))
    cidx = {c: k for k, c in enumerate(palette)}
    K = max(len(palette), 1)
    scale = []
    for k, c in enumerate(palette):
        scale += [[k / K, c], [(k + 1) / K, c]]
    if not scale:
        scale = [[0, SINGLETON], [1, SINGLETON]]

    # close pairs (< threshold) split across clusters, per sample and run
    close = (D < thr) & ~np.eye(n, dtype=bool)
    z, hover = [], []
    for row, s in enumerate(order_names):
        zr, hr = [], []
        for bs, lab in zip(cols, col_labels):
            part = runs[bs]
            c = part.get(s)
            if c is None:
                zr.append(None)
                hr.append(f"<b>{s}</b><br>batch size {lab}<br>not clustered (failed QC)")
                continue
            mem = [t for t in part if part[t] == c]
            i = dist.idx[s]
            split = [names[j] for j in np.nonzero(close[i])[0]
                     if names[j] in part and part[names[j]] != c]
            zr.append(cidx[colour[bs, s]] + 0.5)
            hr.append(f"<b>{s}</b><br>batch size {lab}<br>cluster {c} ({len(mem)} samples)"
                      + (f"<br>{len(split)} sample(s) within {thr:g} SNPs in other clusters: "
                         + ", ".join(sorted(split)[:10]) + (" ..." if len(split) > 10 else "")
                         if split else ""))
        z.append(zr)
        hover.append(hr)

    tree_w, tile_w = 700, max(60 * len(cols), 120)
    fig = make_subplots(rows=1, cols=2, shared_yaxes=True, horizontal_spacing=0.005,
                        column_widths=[tree_w, tile_w])
    fig.add_trace(go.Scatter(x=gx, y=gy, mode="lines", hoverinfo="skip", showlegend=False,
                             line=dict(color="rgba(150,150,150,0.35)", width=1, dash="dot")), 1, 1)
    fig.add_trace(go.Scatter(x=hx + vx, y=hy + vy, mode="lines", hoverinfo="skip",
                             showlegend=False, line=dict(color="#333", width=1.2)), 1, 1)
    fig.add_trace(go.Scatter(
        x=[x[i] for i in order], y=[y[i] for i in order], mode="markers", showlegend=False,
        marker=dict(size=5, color=[colour.get((ref_bs, s), "#fff") for s in order_names],
                    line=dict(width=0.5, color="#333")),
        text=[f"<b>{s}</b><br>root distance {x[i]:.4g} SNPs<br>reference cluster "
              f"{runs[ref_bs].get(s, 'not clustered')}" for s, i in zip(order_names, order)],
        hovertemplate="%{text}<extra></extra>"), 1, 1)
    fig.add_trace(go.Heatmap(z=z, x=col_labels, y=list(range(n)), text=hover,
                             hovertemplate="%{text}<extra></extra>", colorscale=scale,
                             zmin=0, zmax=K, showscale=False, xgap=2, ygap=1), 1, 2)

    height = max(450, 16 * n + 200)
    fig.update_layout(
        height=height, width=tree_w + tile_w + 300, plot_bgcolor="white",
        margin=dict(l=40, r=200, t=150, b=60),
        title=dict(text=(f"Neighbour-joining tree of pairwise SNP distances ({n} samples) "
                         f"and floc clusters by batch size"
                         f"<br><sup>Tiles: colour = cluster, matched to the reference cluster it "
                         f"overlaps most; new colours = clusters with no reference match.<br>"
                         f"Light grey = singleton; blank = not clustered. Hover a tile for samples "
                         f"within {thr:g} SNPs placed in other clusters.</sup>"), x=0.01, y=0.98,
                    yanchor="top"))
    fig.update_xaxes(title_text="SNPs from midpoint root", showgrid=True, gridcolor="#eee",
                     zeroline=False, row=1, col=1)
    fig.update_xaxes(title_text="batch size", side="top", tickangle=0, row=1, col=2)
    fig.update_yaxes(range=[n - 0.5, -0.5], showgrid=False, zeroline=False,
                     showticklabels=False, row=1, col=1)
    fig.update_yaxes(tickmode="array", tickvals=list(range(n)), ticktext=order_names,
                     side="right", showticklabels=True, tickfont=dict(size=10 if n <= 150 else 7),
                     row=1, col=2)
    fig.write_html(a.html, include_plotlyjs=True, full_html=True,
                   config={"displaylogo": False, "toImageButtonOptions": {"format": "svg"}})
    print(f"tree: {n} samples, {len(cols)} batch sizes -> {a.html}, {a.newick}")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("summarise")
    s.add_argument("--clusters", nargs="+", required=True)
    s.add_argument("--distances", required=True, help="all-vs-all distance matrix (long or square)")
    s.add_argument("--n-samples", type=int, required=True)
    s.add_argument("--snp-threshold", type=float, default=100)
    s.add_argument("--distance-source", default="SKA2",
                   help="description of where the distances came from (for the report)")
    t = sub.add_parser("tree")
    t.add_argument("--clusters", nargs="+", required=True)
    t.add_argument("--distances", required=True, help="all-vs-all distance matrix (long or square)")
    t.add_argument("--n-samples", type=int, required=True)
    t.add_argument("--snp-threshold", type=float, default=100)
    t.add_argument("--html", default="nj_tree_clusters.html")
    t.add_argument("--newick", default="nj_tree.nwk")
    a = ap.parse_args(argv)
    {"summarise": cmd_summarise, "tree": cmd_tree}[a.cmd](a)


if __name__ == "__main__":
    main()
