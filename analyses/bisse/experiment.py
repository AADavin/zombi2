#!/usr/bin/env python3
"""A gene drives speciation: the joint runs behind the paper's BiSSE example.

Six driver factors, 200 replicates each on matched seeds, written in the current
grammar: ``joint.simulate(species.birth_death(birth=PerLineage(1.0).scaled_by(...)), ...)``.

Beyond the prevalence summaries, every sweep replicate now writes
  data/trees/f{factor}_r{rep:03d}.nwk    the extant tree (ultrametric, labels n<id>)
  data/states/f{factor}_r{rep:03d}.tsv   tip label, driver presence 0/1, control presence 0/1
which is what the diversitree BiSSE fits consume.

    python experiment.py
"""

from __future__ import annotations

import json
import os
import platform
import subprocess
import sys
import time

import numpy as np

import zombi2
from zombi2 import genomes, joint, species
from zombi2.params import PerLineage

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "results.json")
TREES = os.path.join(HERE, "data", "trees")
STATES = os.path.join(HERE, "data", "states")

# --------------------------------------------------------------------------- design
MASTER_SEED = 20260802             # same as the v0.26.0 run, for comparability
N_REPS = 200
N_EXTANT = 150
BIRTH_BASE = 1.0
DEATH = 0.3
LOSS = 0.12
DUPLICATION = 0.0
ORIGINATION = 0.0
DRIVER = "driver"
CONTROL = "control"
FACTORS = [1.0, 1.25, 1.5, 2.0, 3.0, 5.0]
HEADLINE = 3.0


def make_seeds(master: int, n: int) -> list[int]:
    return [int(s) for s in np.random.default_rng(master).integers(1, 2**31 - 1, size=n)]


def one_run(factor: float, seed: int, *, n_extant=None, total_time=None,
            dump: str | None = None) -> dict:
    """One joint run; optionally dump tree + tip states under tag `dump`."""
    t0 = time.perf_counter()
    stop = {"n_extant": n_extant} if n_extant is not None else {"total_time": total_time}
    try:
        r = joint.simulate(
            species.birth_death(
                birth=PerLineage(BIRTH_BASE).scaled_by(
                    f"genomes:{DRIVER}", {"present": factor, "absent": 1.0}),
                death=DEATH, **stop),
            genomes.genome(duplication=DUPLICATION, loss=LOSS, origination=ORIGINATION,
                           families=[genomes.family(DRIVER), genomes.family(CONTROL)]),
            seed=seed)
    except Exception as exc:            # a total_time run can die out entirely
        if total_time is None:
            raise
        return {"seed": seed, "factor": factor, "n_tips": 0, "height": None,
                "prevalence_driver": None, "prevalence_control": None,
                "extinct": True, "error": type(exc).__name__,
                "seconds": round(time.perf_counter() - t0, 4)}
    elapsed = time.perf_counter() - t0

    tree = r.extant_tree
    if tree is None or not tree.extant_leaves():     # died out before total_time
        return {"seed": seed, "factor": factor, "n_tips": 0, "height": None,
                "prevalence_driver": None, "prevalence_control": None,
                "extinct": True, "seconds": round(elapsed, 4)}
    tips = tree.extant_leaves()
    labels = tree.labels()
    ends = [tree.nodes[i].end_time for i in tips]
    drv = {i: bool(r.genome.has_family(i, DRIVER)) for i in tips}
    ctl = {i: bool(r.genome.has_family(i, CONTROL)) for i in tips}
    rec = {"seed": seed, "factor": factor, "n_tips": len(tips),
           "height": max(ends), "tip_depth_spread": max(ends) - min(ends),
           "prevalence_driver": sum(drv.values()) / len(tips),
           "prevalence_control": sum(ctl.values()) / len(tips),
           "extinct": False, "seconds": round(elapsed, 4)}
    if dump is not None:
        with open(os.path.join(TREES, dump + ".nwk"), "w") as fh:
            fh.write(tree.to_newick() + "\n")
        with open(os.path.join(STATES, dump + ".tsv"), "w") as fh:
            fh.write("tip\tdriver\tcontrol\n")
            for i in tips:
                fh.write(f"{labels[i]}\t{int(drv[i])}\t{int(ctl[i])}\n")
    return rec


def stats(values) -> dict:
    v = np.asarray([x for x in values if x is not None], dtype=float)
    if v.size == 0:
        return {"n": 0}
    return {"n": int(v.size), "mean": float(v.mean()),
            "sd": float(v.std(ddof=1)) if v.size > 1 else 0.0,
            "sem": float(v.std(ddof=1) / np.sqrt(v.size)) if v.size > 1 else 0.0,
            "median": float(np.median(v))}


def bootstrap_ci(values, reps=20000, seed=12345) -> dict:
    v = np.asarray([x for x in values if x is not None], dtype=float)
    if v.size < 2:
        return {}
    rng = np.random.default_rng(seed)
    draws = rng.choice(v, size=(reps, v.size), replace=True).mean(axis=1)
    return {"ci95_low": float(np.percentile(draws, 2.5)),
            "ci95_high": float(np.percentile(draws, 97.5))}


def signflip_p(differences, reps=200000, seed=777) -> dict:
    v = np.asarray([x for x in differences if x is not None], dtype=float)
    if v.size < 2:
        return {}
    rng = np.random.default_rng(seed)
    draws = (v * rng.choice([-1.0, 1.0], size=(reps, v.size))).mean(axis=1)
    hits = int((np.abs(draws) >= abs(v.mean())).sum())
    return {"observed_mean_difference": float(v.mean()), "n_pairs": int(v.size),
            "p_two_sided": max(hits / reps, 1.0 / reps)}


def git_commit() -> str | None:
    pkg_root = os.path.dirname(os.path.dirname(os.path.abspath(zombi2.__file__)))
    try:
        return subprocess.run(["git", "-C", pkg_root, "rev-parse", "HEAD"],
                              capture_output=True, text=True, check=True).stdout.strip()
    except Exception:
        return None


def main() -> int:
    os.makedirs(TREES, exist_ok=True)
    os.makedirs(STATES, exist_ok=True)
    seeds = make_seeds(MASTER_SEED, N_REPS)
    t_all = time.perf_counter()

    sweep: dict[str, list[dict]] = {}
    for f in FACTORS:
        t0 = time.perf_counter()
        sweep[repr(f)] = [one_run(f, s, n_extant=N_EXTANT, dump=f"f{f}_r{i:03d}")
                          for i, s in enumerate(seeds)]
        print(f"  factor {f:>4}: {time.perf_counter() - t0:6.1f} s", flush=True)

    matched = [one_run(1.0, rec["seed"], total_time=rec["height"])
               for rec in sweep[repr(HEADLINE)]]
    print("  duration-matched null: done", flush=True)

    def col(recs, key):
        return [r[key] for r in recs]

    null = sweep[repr(1.0)]
    summary = {"sweep": {}}
    for f in FACTORS:
        recs = sweep[repr(f)]
        paired = [r["prevalence_driver"] - r["prevalence_control"] for r in recs]
        vs_null = [a["prevalence_driver"] - b["prevalence_driver"]
                   for a, b in zip(recs, null)]
        h = np.asarray(col(recs, "height"), dtype=float)
        obs_ctl = np.asarray(col(recs, "prevalence_control"), dtype=float)
        d = obs_ctl - np.exp(-LOSS * h)          # control family vs its closed form
        summary["sweep"][repr(f)] = {
            "factor": f,
            "prevalence_driver": {**stats(col(recs, "prevalence_driver")),
                                  **bootstrap_ci(col(recs, "prevalence_driver"))},
            "prevalence_control": {**stats(col(recs, "prevalence_control")),
                                   **bootstrap_ci(col(recs, "prevalence_control"))},
            "driver_minus_control": {**stats(paired), **bootstrap_ci(paired),
                                     **signflip_p(paired)},
            "driver_minus_null_driver": {**stats(vs_null), **bootstrap_ci(vs_null),
                                         **signflip_p(vs_null)},
            "height": stats(col(recs, "height")),
            "control_vs_exp_minus_loss_height": {
                **stats(d), "z": float(d.mean() / (d.std(ddof=1) / np.sqrt(d.size)))},
        }

    head = sweep[repr(HEADLINE)]
    diff = [c - n for c, n in zip(col(head, "prevalence_driver"),
                                  col(null, "prevalence_driver"))]
    summary["headline"] = {
        "coupled_factor": HEADLINE,
        "coupled_mean_prevalence": stats(col(head, "prevalence_driver"))["mean"],
        "null_mean_prevalence": stats(col(null, "prevalence_driver"))["mean"],
        "paired_coupled_minus_null": {**stats(diff), **bootstrap_ci(diff),
                                      **signflip_p(diff)},
        "within_run_driver_minus_control_at_headline":
            summary["sweep"][repr(HEADLINE)]["driver_minus_control"],
        "coupled_mean_height": stats(col(head, "height"))["mean"],
        "null_mean_height": stats(col(null, "height"))["mean"],
        "duration_matched_null": {
            "prevalence_driver": {**stats(col(matched, "prevalence_driver")),
                                  **bootstrap_ci(col(matched, "prevalence_driver"))},
            "n_runs_that_died_out": int(sum(r["extinct"] for r in matched)),
        },
    }

    payload = {
        "what": "ZOMBI2 joint gene-content -> speciation case study, 0.42.x rerun, "
                "with per-replicate trees + tip states for the BiSSE extension",
        "zombi2_version": zombi2.__version__,
        "git_commit": git_commit(),
        "python": sys.version.split()[0],
        "numpy": np.__version__,
        "platform": platform.platform(),
        "api_call": ("joint.simulate(species.birth_death(birth=PerLineage(1.0)"
                     ".scaled_by('genomes:driver', {'present': f, 'absent': 1.0}), "
                     "death=0.3, n_extant=150), genomes.genome(duplication=0.0, "
                     "loss=0.12, origination=0.0, families=[family('driver'), "
                     "family('control')]), seed=...)"),
        "parameters": {"master_seed": MASTER_SEED, "n_reps": N_REPS,
                       "n_extant": N_EXTANT, "birth_base": BIRTH_BASE, "death": DEATH,
                       "loss": LOSS, "duplication": DUPLICATION,
                       "origination": ORIGINATION, "factors": FACTORS,
                       "headline_factor": HEADLINE},
        "seeds": seeds,
        "summary": summary,
        "replicates": {"sweep": sweep, "duration_matched_null": matched},
        "total_seconds": round(time.perf_counter() - t_all, 2),
    }
    with open(OUT, "w") as fh:
        json.dump(payload, fh, indent=1)

    print(f"\nwrote {OUT}  ({payload['total_seconds']} s total)\n")
    print(f"{'factor':>7} {'driver':>18} {'control':>18} {'driver-control':>18} {'height':>8}")
    for f in FACTORS:
        s = summary["sweep"][repr(f)]
        d, c, p = s["prevalence_driver"], s["prevalence_control"], s["driver_minus_control"]
        print(f"{f:>7} {d['mean']:>8.3f} +-{d['sem']:.3f} {c['mean']:>8.3f} +-{c['sem']:.3f} "
              f"{p['mean']:>8.3f} +-{p['sem']:.3f} {s['height']['mean']:>8.2f}")
    h = summary["headline"]
    print(f"\nheadline f={HEADLINE}: coupled {h['coupled_mean_prevalence']:.3f} "
          f"vs null {h['null_mean_prevalence']:.3f}; "
          f"duration-matched null "
          f"{h['duration_matched_null']['prevalence_driver'].get('mean', float('nan')):.3f}")
    print("\ncontrol vs exp(-loss*height), z by factor:")
    for f in FACTORS:
        c = summary["sweep"][repr(f)]["control_vs_exp_minus_loss_height"]
        print(f"  f={f:>4}: diff {c['mean']:+.4f} (z {c['z']:+.2f})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
