"""Snakemake step 2: perturb one time tree's branch rates, then recover ages with RED.

Reads a time tree, applies the perturbation clock named by ``--spec`` (a JSON blob from
``config['perturbations']``), computes RED on the resulting substitution branch lengths, and
writes: a one-row ``metrics.tsv`` (Pearson/Spearman r, nRMSE, realized fold-range) and a
``points.csv`` (true vs recovered age per internal node, for the scatter figure).
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from red import build_clock, derive_seed, load_tree, red_recovery, strength_of  # noqa: E402

METRIC_COLUMNS = [
    "treeset", "model", "n_tips", "tseed", "pert", "clock", "strength", "cseed", "seed",
    "total_age", "n_internal", "pearson_r", "spearman_r", "nrmse", "fold_range",
    "rate_p95_p5", "rate_logsd",
]


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--tree", required=True)
    ap.add_argument("--spec", required=True, help="JSON perturbation spec (config['perturbations'][pert])")
    ap.add_argument("--treeset", required=True)
    ap.add_argument("--model", required=True)
    ap.add_argument("--n-tips", type=int, required=True)
    ap.add_argument("--pert", required=True)
    ap.add_argument("--tseed", type=int, required=True)
    ap.add_argument("--cseed", type=int, required=True)
    ap.add_argument("--out-metrics", required=True)
    ap.add_argument("--out-points", required=True)
    args = ap.parse_args()

    spec = json.loads(args.spec)
    tree = load_tree(args.tree)
    total_age = tree.total_age
    seed = derive_seed("clock", args.treeset, args.tseed, args.pert, args.cseed)
    scaled = build_clock(spec).scale(tree, seed=seed)
    points, m = red_recovery(tree, scaled, total_age)

    os.makedirs(os.path.dirname(args.out_metrics), exist_ok=True)
    row = dict(
        treeset=args.treeset, model=args.model, n_tips=args.n_tips, tseed=args.tseed,
        pert=args.pert, clock=spec["clock"], strength=strength_of(spec), cseed=args.cseed,
        seed=seed, **m,
    )
    with open(args.out_metrics, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(METRIC_COLUMNS)
        w.writerow([_fmt(row[c]) for c in METRIC_COLUMNS])
    with open(args.out_points, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["node", "true_age", "recovered_age"])
        w.writerows(points)
    print(f"[red] {args.treeset} s{args.tseed} {args.pert} r{args.cseed}: "
          f"r={m['pearson_r']:.4f} nRMSE={100*m['nrmse']:.2f}% fold={m['fold_range']:.1f}x")


def _fmt(v):
    if isinstance(v, float):
        return f"{v:.6g}"
    return v


if __name__ == "__main__":
    main()
