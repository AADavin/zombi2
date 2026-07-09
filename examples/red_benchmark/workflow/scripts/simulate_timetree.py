"""Snakemake step 1: simulate one ultrametric time tree (known node ages) → Newick.

Run standalone or via the workflow. The tree is simulated once per (config, seed) and reused
across every perturbation × replicate, so this is the cheap, shared root of the DAG.
"""

from __future__ import annotations

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))  # import the sibling red module
from red import derive_seed, simulate_time_tree  # noqa: E402


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--model", required=True, choices=["yule", "bd"])
    ap.add_argument("--n-tips", type=int, required=True)
    ap.add_argument("--birth", type=float, default=1.0)
    ap.add_argument("--death", type=float, default=0.0)
    ap.add_argument("--treeset", required=True)   # only used to derive an independent seed
    ap.add_argument("--tseed", type=int, required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    seed = derive_seed("tree", args.treeset, args.tseed)
    tree = simulate_time_tree(args.model, args.n_tips, args.birth, args.death, seed)
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    with open(args.out, "w") as fh:
        fh.write(tree.to_newick() + "\n")
    print(f"[simulate] {args.treeset} s{args.tseed}: {len(tree.leaves())} tips, "
          f"total_age={tree.total_age:.4f}, seed={seed} -> {args.out}")


if __name__ == "__main__":
    main()
