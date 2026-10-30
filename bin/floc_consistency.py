#!/usr/bin/env python3
"""
floc_consistency.py - helper for the floc-consistency Nextflow pipeline.
Standard library only, so it runs in the floc container.

  plan       read the final floc clusters of every batch-size run and write the
             SKA groups (group<TAB>sample) whose samples need split k-mer files
  summarise  combine the floc clusters with the SKA distances and write the
             discrepant-pair, within-cluster and summary tables

A run's clusters are given as files named bs<BATCH_SIZE>.clusters.csv (floc's
clusters.csv: columns sid, cluster). The reference is the run whose batch
size equals the number of samples (all samples in one batch).
"""
from __future__ import annotations

import argparse
import csv
import re
import statistics
import sys
from collections import defaultdict
from itertools import combinations
from pathlib import Path

STATS = ["n", "min", "max", "mean", "sd", "median"]


# --------------------------------------------------------------------------- #
# Reading inputs
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


def pairs_to_measure(runs: dict[int, dict], disc: list[dict]) -> set[tuple[str, str]]:
    pairs = {(r["sample_a"], r["sample_b"]) for r in disc}
    for part in runs.values():
        for mem in members(part).values():
            pairs.update(pairs_of(mem))
    return pairs


def components(pairs: set[tuple[str, str]]) -> list[list[str]]:
    parent: dict[str, str] = {}

    def find(x):
        parent.setdefault(x, x)
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    for a, b in pairs:
        parent[find(a)] = find(b)
    comps = defaultdict(list)
    for x in list(parent):
        comps[find(x)].append(x)
    return sorted((sorted(c) for c in comps.values()), key=lambda c: (-len(c), c[0]))


# --------------------------------------------------------------------------- #
# Statistics and output
# --------------------------------------------------------------------------- #
def describe(values) -> dict:
    x = [float(v) for v in values if v is not None]
    if not x:
        return {"n": 0, "min": "", "max": "", "mean": "", "sd": "", "median": ""}
    return {"n": len(x), "min": min(x), "max": max(x), "mean": statistics.fmean(x),
            "sd": statistics.stdev(x) if len(x) > 1 else "", "median": statistics.median(x)}


def fmt(v):
    if isinstance(v, float):
        return f"{v:.6g}"
    return str(v)


def write_tsv(path: str, rows: list[dict], cols: list[str], header_lines: list[str] = ()):
    with open(path, "w", newline="") as fh:
        for line in header_lines:
            fh.write(f"# {line}\n")
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(fmt(r.get(c, "")) for c in cols) + "\n")


def read_distances(files: list[str]) -> dict[frozenset, tuple[float, int]]:
    """ska distance output: Sample1 Sample2 Distance Mismatches Match_count ..."""
    dist = {}
    for f in files:
        with open(f, newline="") as fh:
            reader = csv.reader(fh, delimiter="\t")
            next(reader, None)
            for row in reader:
                if len(row) >= 5:
                    dist[frozenset((row[0], row[1]))] = (float(row[2]), int(float(row[4])))
    return dist


# --------------------------------------------------------------------------- #
# Commands
# --------------------------------------------------------------------------- #
def cmd_plan(a):
    runs = read_runs(a.clusters)
    ref_bs = reference_size(runs, a.n_samples)
    disc = discrepant_pairs(runs, ref_bs)
    pairs = pairs_to_measure(runs, disc)
    if a.core == "group":
        groups = components(pairs)
    else:   # global / pairwise: one group with every sample that needs measuring
        groups = [sorted({s for p in pairs for s in p})] if pairs else []
    with open(a.out, "w") as fh:
        fh.write("group\tsample\n")
        for i, g in enumerate(groups, 1):
            for s in g:
                fh.write(f"g{i:05d}\t{s}\n")
    print(f"{len(disc)} discrepant pairs; {len(pairs)} pairs to measure in "
          f"{len(groups)} SKA group(s)")


def cmd_summarise(a):
    runs = read_runs(a.clusters)
    n = a.n_samples
    ref_bs = reference_size(runs, n)
    dist = read_distances(a.distances)
    disc = discrepant_pairs(runs, ref_bs)
    thr = a.snp_threshold

    def snp(x, y):
        return dist.get(frozenset((x, y)), (None, 0))

    missing = 0
    for r in disc:
        r["core_snps"], r["core_sites"] = snp(r["sample_a"], r["sample_b"])
        missing += r["core_snps"] is None
        r[f"core_snps_below_{thr:g}"] = (r["core_snps"] is not None and r["core_snps"] < thr)

    within, cstats = [], []
    order = [ref_bs] + [bs for bs in runs if bs != ref_bs]
    for bs in order:
        for c, mem in sorted(members(runs[bs]).items(), key=lambda kv: (-len(kv[1]), kv[0])):
            vals = []
            for x, y in pairs_of(mem):
                s, sites = snp(x, y)
                missing += s is None
                vals.append(s)
                within.append({"batch_size": bs, "is_reference": bs == ref_bs, "cluster": c,
                               "cluster_size": len(mem), "sample_a": x, "sample_b": y,
                               "core_snps": s, "core_sites": sites})
            d = describe(vals)
            cstats.append({"batch_size": bs, "is_reference": bs == ref_bs, "cluster": c,
                           "cluster_size": len(mem), "n_pairs": d["n"],
                           **{f"core_snps_{k}": d[k] for k in STATS if k != "n"},
                           "members": ",".join(mem)})
    if missing:
        print(f"WARNING: {missing} pairs had no SKA distance", file=sys.stderr)

    summary = []
    for label, sel in (("all", disc),
                       ("together_in_reference", [r for r in disc if r["pair_class"] == "together_in_reference"]),
                       ("apart_in_reference", [r for r in disc if r["pair_class"] == "apart_in_reference"])):
        summary.append({"section": "discrepant pairs", "group": label,
                        **describe(r["core_snps"] for r in sel)})
    for bs in order:
        label = f"batch size {bs}" + (" (reference)" if bs == ref_bs else "")
        summary.append({"section": "within-cluster pairs", "group": label,
                        **describe(w["core_snps"] for w in within if w["batch_size"] == bs)})
    uniq = {(w["sample_a"], w["sample_b"]): w["core_snps"] for w in within}
    summary.append({"section": "within-cluster pairs", "group": "all runs (unique pairs)",
                    **describe(uniq.values())})
    summary.append({"section": "within-cluster pairs", "group": "all runs (every pair in every run)",
                    **describe(w["core_snps"] for w in within)})

    close = [r for r in disc if r[f"core_snps_below_{thr:g}"]]
    close_samples = sorted({s for r in close for s in (r["sample_a"], r["sample_b"])})
    disc_samples = sorted({s for r in disc for s in (r["sample_a"], r["sample_b"])})
    header = [f"reference: batch size {ref_bs} (all {n} samples in one batch); runs: "
              f"{len(runs)}; core SNPs: SKA2 k={a.ska_k}, core={a.core}",
              f"discrepant pairs: {len(disc)} ({len(disc_samples)} samples)",
              f"discrepant pairs with < {thr:g} core SNPs: {len(close)} "
              f"({len(close_samples)} samples: {', '.join(close_samples) or 'none'})"]
    for bs in order:
        cs = [c for c in cstats if c["batch_size"] == bs]
        header.append(f"batch size {bs}: {len(cs)} clusters, "
                      f"{sum(c['cluster_size'] > 1 for c in cs)} with 2+ samples, "
                      f"{sum(c['cluster_size'] == 1 for c in cs)} singletons")

    disc_cols = ["sample_a", "sample_b", "reference_cluster_a", "reference_cluster_b",
                 "pair_class", "runs_together", "runs_with_both", "runs_different_from_reference",
                 "batch_sizes_different_from_reference", "core_snps", "core_sites",
                 f"core_snps_below_{thr:g}"]
    disc.sort(key=lambda r: (r["core_snps"] is None, r["core_snps"] or 0, r["sample_a"], r["sample_b"]))
    write_tsv("discrepant_pairs.tsv", disc, disc_cols)
    write_tsv("within_cluster_pairs.tsv", within,
              ["batch_size", "is_reference", "cluster", "cluster_size", "sample_a", "sample_b",
               "core_snps", "core_sites"])
    write_tsv("cluster_stats.tsv", cstats,
              ["batch_size", "is_reference", "cluster", "cluster_size", "n_pairs"]
              + [f"core_snps_{k}" for k in STATS if k != "n"] + ["members"])
    write_tsv("summary.tsv", summary, ["section", "group"] + STATS, header)

    # console report (shown in the task log)
    print("\n".join(header))
    print("\t".join(["section", "group"] + STATS))
    for r in summary:
        print("\t".join(fmt(r[c]) for c in ["section", "group"] + STATS))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("plan")
    p.add_argument("--clusters", nargs="+", required=True)
    p.add_argument("--n-samples", type=int, required=True)
    p.add_argument("--core", choices=["group", "global", "pairwise"], default="group")
    p.add_argument("--out", default="ska_groups.tsv")
    s = sub.add_parser("summarise")
    s.add_argument("--clusters", nargs="+", required=True)
    s.add_argument("--distances", nargs="*", default=[])
    s.add_argument("--n-samples", type=int, required=True)
    s.add_argument("--snp-threshold", type=float, default=100)
    s.add_argument("--core", default="group")
    s.add_argument("--ska-k", type=int, default=31)
    a = ap.parse_args(argv)
    {"plan": cmd_plan, "summarise": cmd_summarise}[a.cmd](a)


if __name__ == "__main__":
    main()
