#!/usr/bin/env python3
"""What can gene order constrain? An ABC recovery test on simulated truth.

Genomes are simulated at known rearrangement rates down a dated tree, then the parameters are
inferred back from gene order alone and compared to the values used to generate the data. Nothing
here is fitted to real data, so every claim can be checked against a known answer.

Three arms:

  A  inversions only        grid over (inversion rate x inversion extent)
  B  mixed                  grid over (inversion rate x translocation rate), extent fixed at truth
  C  extent misspecified    the same truth as B, with the extent fixed at twice the truth

Arm A asks which of the two inversion parameters the data pin. Arm B asks whether a mixed model's
two rates can be told apart, which is the regime that motivates ABC. Arm C asks what it costs to fix
the extent at the wrong value, which is the practical consequence of arm A's answer.

Writes ``results.json`` (committed) and ``data/tree.nwk`` (the dated tree, regenerable from
TREE_SEED). Figures come from ``figures.py``.

    python experiment.py
"""

from __future__ import annotations

import json
import multiprocessing as mp
import os
import platform
import subprocess
import sys
import time

import numpy as np

import sim
import synteny
import zombi2

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "results.json")
TREE_FILE = os.path.join(HERE, "data", "tree.nwk")

# ------------------------------------------------------------------------------------- design
MASTER_SEED = 20260817
TREE_SEED = 20260817
N_TIPS = 40
CROWN_DEPTH = 100.0          # tree time units, so patristic divergence runs 0 to 200
BIRTH, DEATH = 1.0, 0.3
N_GENES = 2000
N_CHROMOSOMES = 8
BINS = (0.0, 60.0, 120.0, 200.0)
N_REPLICATES = 3             # independent simulations per grid cell, averaged
N_OBSERVED = 3               # independent "observed" datasets, each inferred from separately

# The truth. Rates are per gene per unit of tree time; the extent is a mean, in genes.
TRUE_INVERSION = 1.0e-3
TRUE_EXTENT = 4.0
TRUE_TRANSLOCATION = 2.5e-4         # a translocation share of 0.25 / 1.25 = 0.20
MISSPECIFIED_EXTENT = 8.0           # twice the truth, used by arm C

# Grids. Each contains the true value, so exact recovery is possible.
RATE_STEP = 2.0 ** 0.5
INVERSION_GRID_A = [TRUE_INVERSION * RATE_STEP ** k for k in range(-6, 7)]
EXTENT_GRID_A = [1.0, 2.0, 4.0, 8.0, 16.0, 32.0, 64.0]
INVERSION_GRID_B = [TRUE_INVERSION * RATE_STEP ** k for k in range(-4, 5)]
TRANSLOCATION_GRID_B = [0.0, 6.25e-5, 1.25e-4, 2.5e-4, 5.0e-4, 1.0e-3]

# Set once in the parent; forked workers inherit them.
_TREE = _NAMEMAP = _DIVERGENCE = None


def cell_seed(arm: str, cell: int, replicate: int) -> int:
    """A seed fixed by (arm, cell, replicate), so the run is reproducible and does not depend on
    how many workers happen to be used."""
    entropy = [MASTER_SEED, sum(arm.encode()), cell, replicate]
    return int(np.random.default_rng(entropy).integers(1, 2**31 - 1))


def _one(job):
    """Simulate one genome set and reduce it to a summary vector. Runs in a worker."""
    inversion, extent, translocation, seed = job
    genomes = sim.simulate_signed_order(
        _TREE, _NAMEMAP, inversion=inversion, inversion_extent=extent,
        translocation=translocation, n_genes=N_GENES, n_chromosomes=N_CHROMOSOMES, seed=seed)
    return synteny.summary(genomes, _DIVERGENCE, bins=BINS)


def run_jobs(jobs, workers):
    if workers == 1:
        return [_one(j) for j in jobs]
    # fork, so workers inherit the tree and the divergence map instead of rebuilding them
    with mp.get_context("fork").Pool(workers) as pool:
        return pool.map(_one, jobs, chunksize=1)


def average(summaries, keys) -> dict:
    """Mean of each key across a cell's replicates, ignoring missing values."""
    out = {}
    for k in keys:
        vals = [s[k] for s in summaries if not np.isnan(s.get(k, np.nan))]
        out[k] = float(np.mean(vals)) if vals else np.nan
    return out


def infer(cells, targets, keys, axes):
    """Pick the closest cell to each target, and the pseudo-posterior marginals.

    ``cells`` is the list of averaged summary vectors, one per grid cell, in the order of
    ``axes["order"]`` (a list of parameter tuples). Returns one record per target.
    """
    weights = {k: 1.0 / s ** 2 for k, s in synteny.scales(cells + targets, keys).items()}
    order = axes["order"]
    out = []
    for target in targets:
        distances = np.array([synteny.distance(c, target, keys, weights) for c in cells])
        best = int(np.nanargmin(distances))
        record = {"best": {name: order[best][i] for i, name in enumerate(axes["names"])},
                  "distance_min": float(distances[best]),
                  "distances": [float(x) for x in distances]}
        # a soft credible interval: weight cells by how close they land
        epsilon = float(np.nanpercentile(distances, 10))
        w = np.exp(-(distances / max(epsilon, 1e-12)) ** 2)
        w[~np.isfinite(w)] = 0.0
        record["marginals"] = {}
        for i, name in enumerate(axes["names"]):
            values = np.array([p[i] for p in order], float)
            levels = np.array(sorted(set(values)))
            mass = np.array([w[values == v].sum() for v in levels])
            mass = mass / mass.sum() if mass.sum() > 0 else mass
            cdf = np.cumsum(mass)
            record["marginals"][name] = {
                "levels": [float(v) for v in levels],
                "mass": [float(m) for m in mass],
                "median": float(levels[int(np.searchsorted(cdf, 0.5))]),
                "q05": float(levels[int(np.searchsorted(cdf, 0.05))]),
                "q95": float(levels[min(int(np.searchsorted(cdf, 0.95)), len(levels) - 1)]),
                "best_profile": [float(np.nanmin(distances[values == v])) for v in levels],
            }
        out.append(record)
    return out


def sweep(arm, grid, names, *, extent=None, use_cross, workers, label):
    """Simulate every cell of one grid, ``N_REPLICATES`` times, and average each cell."""
    keys = synteny.fit_keys(BINS, use_cross=use_cross)
    jobs, cell_of = [], []
    for index, point in enumerate(grid):
        for replicate in range(N_REPLICATES):
            values = dict(zip(names, point))
            jobs.append((values.get("inversion"),
                         extent if extent is not None else values.get("extent"),
                         values.get("translocation", 0.0),
                         cell_seed(arm, index, replicate)))
            cell_of.append(index)
    started = time.perf_counter()
    summaries = run_jobs(jobs, workers)
    elapsed = time.perf_counter() - started
    cells = [average([summaries[j] for j in range(len(jobs)) if cell_of[j] == i],
                     synteny.all_keys(BINS))
             for i in range(len(grid))]
    print(f"  {label}: {len(grid)} cells x {N_REPLICATES} = {len(jobs)} runs "
          f"in {elapsed:.0f} s", flush=True)
    return cells, keys, elapsed


def observed(arm, *, inversion, extent, translocation, workers):
    """The datasets the study treats as observed: independent runs at the true parameters."""
    jobs = [(inversion, extent, translocation, cell_seed(arm + ":observed", 0, r))
            for r in range(N_OBSERVED)]
    return run_jobs(jobs, workers), [j[3] for j in jobs]


def git_commit():
    root = os.path.dirname(os.path.dirname(os.path.abspath(zombi2.__file__)))
    try:
        return subprocess.run(["git", "-C", root, "rev-parse", "HEAD"],
                              capture_output=True, text=True, check=True).stdout.strip()
    except Exception:
        return None


def main() -> int:
    global _TREE, _NAMEMAP, _DIVERGENCE
    workers = max(1, mp.cpu_count() - 1)
    os.makedirs(os.path.dirname(TREE_FILE), exist_ok=True)
    started = time.perf_counter()

    _TREE, _NAMEMAP, newick = sim.dated_tree(
        seed=TREE_SEED, n_tips=N_TIPS, crown_depth=CROWN_DEPTH, birth=BIRTH, death=DEATH)
    with open(TREE_FILE, "w") as fh:
        fh.write(newick + "\n")
    _DIVERGENCE = synteny.pairwise_divergence(_TREE, _NAMEMAP)
    times = np.array(list(_DIVERGENCE.values()))
    bin_counts = [int(((times >= BINS[i]) & (times < BINS[i + 1])).sum())
                  for i in range(len(BINS) - 1)]
    bin_means = [float(times[(times >= BINS[i]) & (times < BINS[i + 1])].mean())
                 if bin_counts[i] else float("nan") for i in range(len(BINS) - 1)]
    print(f"tree: {N_TIPS} tips, crown depth {CROWN_DEPTH:g}, {len(times)} pairs, "
          f"bin counts {bin_counts}", flush=True)
    print(f"workers: {workers}\n", flush=True)

    arms = {}

    # ------------------------------------------------------------------- arm A: inversions only
    print("arm A - inversions only", flush=True)
    grid_a = [(inv, ext) for ext in EXTENT_GRID_A for inv in INVERSION_GRID_A]
    obs_a, seeds_a = observed("A", inversion=TRUE_INVERSION, extent=TRUE_EXTENT,
                              translocation=0.0, workers=workers)
    cells_a, keys_a, secs_a = sweep("A", grid_a, ("inversion", "extent"), use_cross=False,
                                    workers=workers, label="grid")
    axes_a = {"names": ("inversion", "extent"), "order": grid_a}
    arms["A"] = {
        "what": "inversions only; grid over inversion rate x inversion extent",
        "truth": {"inversion": TRUE_INVERSION, "extent": TRUE_EXTENT, "translocation": 0.0},
        "grid": {"inversion": INVERSION_GRID_A, "extent": EXTENT_GRID_A},
        "fit_statistics": keys_a,
        "observed_seeds": seeds_a,
        "observed_summaries": obs_a,
        "cell_summaries": cells_a,
        "inference": infer(cells_a, obs_a, keys_a, axes_a),
        "seconds": round(secs_a, 1),
    }

    # ------------------------------------------------------------------------- arm B: mixed
    print("arm B - inversions and translocations", flush=True)
    grid_b = [(inv, tra) for tra in TRANSLOCATION_GRID_B for inv in INVERSION_GRID_B]
    obs_b, seeds_b = observed("B", inversion=TRUE_INVERSION, extent=TRUE_EXTENT,
                              translocation=TRUE_TRANSLOCATION, workers=workers)
    cells_b, keys_b, secs_b = sweep("B", grid_b, ("inversion", "translocation"),
                                    extent=TRUE_EXTENT, use_cross=True, workers=workers,
                                    label="grid")
    axes_b = {"names": ("inversion", "translocation"), "order": grid_b}
    arms["B"] = {
        "what": "inversions and translocations; grid over the two rates, extent fixed at truth",
        "truth": {"inversion": TRUE_INVERSION, "extent": TRUE_EXTENT,
                  "translocation": TRUE_TRANSLOCATION,
                  "translocation_share": TRUE_TRANSLOCATION / (TRUE_INVERSION + TRUE_TRANSLOCATION)},
        "grid": {"inversion": INVERSION_GRID_B, "translocation": TRANSLOCATION_GRID_B},
        "assumed_extent": TRUE_EXTENT,
        "fit_statistics": keys_b,
        "observed_seeds": seeds_b,
        "observed_summaries": obs_b,
        "cell_summaries": cells_b,
        "inference": infer(cells_b, obs_b, keys_b, axes_b),
        "seconds": round(secs_b, 1),
    }

    # --------------------------------------------------------------- arm C: extent misspecified
    print("arm C - the same truth, extent fixed at twice the truth", flush=True)
    cells_c, keys_c, secs_c = sweep("C", grid_b, ("inversion", "translocation"),
                                    extent=MISSPECIFIED_EXTENT, use_cross=True, workers=workers,
                                    label="grid")
    arms["C"] = {
        "what": "the same observed data as arm B, fitted with the extent fixed at twice the truth",
        "truth": arms["B"]["truth"],
        "grid": arms["B"]["grid"],
        "assumed_extent": MISSPECIFIED_EXTENT,
        "fit_statistics": keys_c,
        "observed_seeds": seeds_b,
        "cell_summaries": cells_c,
        "inference": infer(cells_c, obs_b, keys_c, axes_b),
        "seconds": round(secs_c, 1),
    }

    # ---------------------------------------------------------------------------- the verdicts
    verdicts = {}
    for name in ("A", "B", "C"):
        arm = arms[name]
        recovered = {}
        for axis in arm["inference"][0]["best"]:
            best = [rec["best"][axis] for rec in arm["inference"]]
            median = [rec["marginals"][axis]["median"] for rec in arm["inference"]]
            truth = arm["truth"].get(axis)
            entry = {"true": truth, "best_per_observed": best, "median_per_observed": median}
            if truth:
                entry["ratio_best_to_true"] = [float(b / truth) for b in best]
            # How much of the grid the credible interval covers. A constrained axis keeps a few
            # levels around the truth; an unconstrained one keeps nearly all of them. This is
            # grid-relative, so it means the same thing on both axes.
            spans = []
            for rec in arm["inference"]:
                m = rec["marginals"][axis]
                levels = np.array(m["levels"], float)
                inside = ((levels >= m["q05"]) & (levels <= m["q95"])).sum()
                spans.append(float(inside / len(levels)))
            entry["credible_span_fraction"] = spans
            entry["q05_per_observed"] = [rec["marginals"][axis]["q05"]
                                         for rec in arm["inference"]]
            entry["q95_per_observed"] = [rec["marginals"][axis]["q95"]
                                         for rec in arm["inference"]]
            recovered[axis] = entry
        if name in ("B", "C"):
            shares = []
            for rec in arm["inference"]:
                inv, tra = rec["best"]["inversion"], rec["best"]["translocation"]
                shares.append(float(tra / (inv + tra)) if inv + tra else float("nan"))
            recovered["translocation_share"] = {
                "true": arm["truth"]["translocation_share"], "best_per_observed": shares}
        verdicts[name] = recovered

    payload = {
        "what": "What can gene order constrain? ABC recovery of rearrangement parameters on "
                "simulated truth, with the nucleotide genome model.",
        "zombi2_version": zombi2.__version__,
        "git_commit": git_commit(),
        "python": sys.version.split()[0],
        "numpy": np.__version__,
        "platform": platform.platform(),
        "api_call": ("genomes.simulate_genomes_nucleotide(tree, inversion=rate*n_genes, "
                     "inversion_extent=extent_bp, translocation=rate*n_genes, "
                     "translocation_extent=extent_bp, genes=..., gene_length=1000, "
                     "chromosomes=8, root_length=..., topology='linear', seed=...)"),
        "design": {
            "master_seed": MASTER_SEED, "tree_seed": TREE_SEED, "n_tips": N_TIPS,
            "crown_depth": CROWN_DEPTH, "birth": BIRTH, "death": DEATH,
            "n_genes": N_GENES, "n_chromosomes": N_CHROMOSOMES,
            "gene_length": sim.GENE_LENGTH, "spacing": sim.SPACING,
            "bins": list(BINS), "bin_counts": bin_counts, "bin_mean_divergence": bin_means,
            "n_pairs": int(len(times)),
            "n_replicates_per_cell": N_REPLICATES, "n_observed": N_OBSERVED,
            "rate_units": "events per gene per unit of tree time",
            "extent_units": "mean extent in genes",
        },
        "statistics": {
            "conservation": "gene-order conservation, the normalised complement of the "
                            "breakpoint distance",
            "segment": "mean conserved segment length, in genes",
            "cross": "fraction of broken adjacencies whose genes lie on different chromosomes",
        },
        "arms": arms,
        "verdicts": verdicts,
        "total_seconds": round(time.perf_counter() - started, 1),
    }
    with open(OUT, "w") as fh:
        json.dump(payload, fh, indent=1)

    # ------------------------------------------------------------------------------- reporting
    print(f"\nwrote {OUT}  ({payload['total_seconds']:.0f} s total)\n")
    for name in ("A", "B", "C"):
        arm, verdict = arms[name], verdicts[name]
        print(f"arm {name} - {arm['what']}")
        for axis, entry in verdict.items():
            if entry.get("true") in (None, 0.0):
                continue
            best = entry["best_per_observed"]
            shown = ", ".join(f"{b:.4g}" for b in best)
            print(f"    {axis:<20} true {entry['true']:<10.4g} recovered {shown}"
                  + (f"   (x{np.mean(entry['ratio_best_to_true']):.2f} of truth)"
                     if "ratio_best_to_true" in entry else ""))
        for axis in arm["inference"][0]["best"]:
            span = float(np.mean(verdict[axis]["credible_span_fraction"]))
            print(f"    {axis:<20} credible interval covers {span:.0%} of the grid "
                  f"({'not constrained' if span > 0.7 else 'constrained'})")
        print()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
