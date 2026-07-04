#!/usr/bin/env python3.12
"""
Reusable timing harness for the ZOMBI 1 (legacy) GENE-FAMILY (genome) simulation step.

For a target extant-tip count N and a seed, it:
  (a) generates a ZOMBI 1 species tree with ~N EXTANT tips using the `T` mode,
  (b) times ONLY the `Gm` genome step (wall-clock, subprocess, hard timeout),
  (c) prints one JSON line:
        {"n_tips": <actual extant tips>, "seconds": <wall clock>,
         "n_gene_families": <count>, "status": "ok"|"timeout"|"explode"|"error"}

Rate regime (non-exploding, ZOMBI2-comparable): see GenomeParameters_bench.tsv
  D=0.2  T=0.1  L=0.25  O=0.5, initial genome size 20, REPLACEMENT_TRANSFER=1.

Usage:
  python3.12 harness.py --n 1000 --seed 1
  python3.12 harness.py --n 1000 --seed 1 --timeout 180 --keep
  python3.12 harness.py --sizes 100 300 1000 3000 10000 --seed 1   # -> results.json lines

The T-mode tree build is NOT timed (only Gm is). We use STOPPING_RULE=1 with
EXTINCTION=0 so exactly N extant lineages are produced (clean tip control, no retries).
"""

import argparse
import json
import os
import shutil
import subprocess
import sys
import time

# ---- Fixed locations -------------------------------------------------------
# ZOMBI 1 checkout: override with $ZOMBI1_DIR, else the default clone location.
ZOMBI_DIR = os.environ.get("ZOMBI1_DIR", "/Users/aadria/Desktop/Github/ZOMBI")
BENCH_DIR = os.path.dirname(os.path.abspath(__file__))
SP_TEMPLATE = os.path.join(BENCH_DIR, "SpeciesTreeParameters_template.tsv")
GENOME_PARAMS = os.path.join(BENCH_DIR, "GenomeParameters_bench.tsv")
PYTHON = sys.executable  # python3.12 that launched us
ZOMBI_PY = os.path.join(ZOMBI_DIR, "Zombi.py")


def count_extant_tips(extant_tree_file):
    """Count leaves in T/ExtantTree.nwk without importing ete3 (fast, robust)."""
    with open(extant_tree_file) as f:
        nwk = f.read().strip()
    if not nwk or nwk == ";":
        return 0
    # Leaves are the comma/paren-delimited tokens that carry a name before ':'.
    # Robust approach: a leaf is any name token that is immediately preceded by
    # '(' or ',' in the newick string. Count them.
    import re
    # names look like  n42:0.123  ; strip branch lengths and internal-node labels.
    # Internal labels appear right after ')'. Leaf labels appear after '(' or ','.
    leaves = re.findall(r'[(,]\s*([^(),:;]+):', nwk)
    return len(leaves)


def count_gene_families(genome_dir):
    """Count gene families from G/Gene_families/*_events.tsv."""
    gf_dir = os.path.join(genome_dir, "Gene_families")
    if not os.path.isdir(gf_dir):
        # fall back to profiles column count
        prof = os.path.join(genome_dir, "Profiles", "Profiles.tsv")
        if os.path.isfile(prof):
            with open(prof) as f:
                header = f.readline().strip().split("\t")
            return max(0, len(header) - 1)
        return 0
    return sum(1 for x in os.listdir(gf_dir) if x.endswith("_events.tsv"))


def build_species_tree(out_dir, n_tips, seed):
    """Run T mode to build a species tree with ~n_tips extant tips. NOT timed."""
    sp_params = os.path.join(out_dir, "sp_params.tsv")
    with open(SP_TEMPLATE) as f:
        tmpl = f.read()
    tmpl = tmpl.replace("__NTIPS__", str(int(n_tips))).replace("__SEED__", str(int(seed)))
    os.makedirs(out_dir, exist_ok=True)
    with open(sp_params, "w") as f:
        f.write(tmpl)

    proc = subprocess.run(
        [PYTHON, ZOMBI_PY, "T", sp_params, out_dir],
        cwd=ZOMBI_DIR, capture_output=True, text=True,
    )
    extant = os.path.join(out_dir, "T", "ExtantTree.nwk")
    if proc.returncode != 0 or not os.path.isfile(extant):
        return None, proc.stderr or proc.stdout
    return extant, None


def time_genome_step(out_dir, timeout_s):
    """Time ONLY the Gm genome step (subprocess, hard timeout).

    Returns (seconds, status, detail).
    status in {"ok", "timeout", "error"}.  "explode" is decided by the caller
    (a timeout under a non-exploding regime is reported as "timeout"; the
    default-parameter genome-explosion case is what 'explode' documents).
    """
    # Gm reuses / recreates the G/ folder itself.
    t0 = time.perf_counter()
    try:
        proc = subprocess.run(
            [PYTHON, ZOMBI_PY, "Gm", GENOME_PARAMS, out_dir],
            cwd=ZOMBI_DIR, capture_output=True, text=True, timeout=timeout_s,
        )
    except subprocess.TimeoutExpired:
        return timeout_s, "timeout", "Gm exceeded hard timeout"
    seconds = time.perf_counter() - t0
    if proc.returncode != 0:
        return seconds, "error", (proc.stderr or proc.stdout)[-2000:]
    return seconds, "ok", None


def run_one(n_tips, seed, timeout_s=180, keep=False, workroot=None):
    workroot = workroot or os.path.join(BENCH_DIR, "runs")
    os.makedirs(workroot, exist_ok=True)
    out_dir = os.path.join(workroot, f"n{n_tips}_s{seed}")
    if os.path.isdir(out_dir):
        shutil.rmtree(out_dir)

    result = {"target_n": int(n_tips), "seed": int(seed), "n_tips": None,
              "seconds": None, "n_gene_families": None, "status": None}

    extant, err = build_species_tree(out_dir, n_tips, seed)
    if extant is None:
        result["status"] = "error"
        result["detail"] = "species-tree build failed: " + (err or "")[-500:]
        if not keep:
            shutil.rmtree(out_dir, ignore_errors=True)
        return result

    result["n_tips"] = count_extant_tips(extant)

    seconds, status, detail = time_genome_step(out_dir, timeout_s)
    result["seconds"] = round(seconds, 3)
    result["status"] = status
    if detail:
        result["detail"] = detail

    genome_dir = os.path.join(out_dir, "G")
    result["n_gene_families"] = count_gene_families(genome_dir)

    if not keep:
        shutil.rmtree(out_dir, ignore_errors=True)
    return result


def main():
    ap = argparse.ArgumentParser(description="ZOMBI 1 Gm (genome step) timing harness.")
    ap.add_argument("--n", type=int, help="single target extant-tip count")
    ap.add_argument("--sizes", type=int, nargs="+", help="multiple target sizes")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--repeats", type=int, default=1, help="repeats per size (varies seed)")
    ap.add_argument("--timeout", type=float, default=180.0, help="hard Gm timeout (s)")
    ap.add_argument("--keep", action="store_true", help="keep output dirs (default: delete)")
    ap.add_argument("--workroot", type=str, default=None)
    args = ap.parse_args()

    sizes = args.sizes if args.sizes else ([args.n] if args.n else None)
    if not sizes:
        ap.error("provide --n or --sizes")

    for n in sizes:
        for r in range(args.repeats):
            seed = args.seed + r
            res = run_one(n, seed, timeout_s=args.timeout,
                          keep=args.keep, workroot=args.workroot)
            print(json.dumps(res), flush=True)


if __name__ == "__main__":
    main()
