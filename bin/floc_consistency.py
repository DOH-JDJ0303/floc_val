#!/usr/bin/env python3
"""
floc_consistency.py - helper for the floc-consistency Nextflow pipeline.
Standard library only, so it runs in the floc container.

  summarise  combine the final floc clusters of every batch-size run with an
             all-vs-all SNP distance matrix and write intra-/inter-cluster
             distance statistics, discrepant pairs and close pairs that floc
             placed in different clusters

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
    a = ap.parse_args(argv)
    {"summarise": cmd_summarise}[a.cmd](a)


if __name__ == "__main__":
    main()
