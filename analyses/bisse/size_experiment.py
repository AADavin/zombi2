#!/usr/bin/env python3
"""Tree-size arms for the BiSSE example: does power recover on bigger trees?

f in {3, 5} x n_extant in {500, 1000}, 100 replicates each, on the same first 100
seeds as the 150-tip sweep so the arms are comparable. Trees and tip states land
in data_size/, per-run summaries in results_size.json.

    python size_experiment.py
"""

from __future__ import annotations

import json
import os
import time

import numpy as np

import zombi2
from zombi2 import genomes, joint, species
from zombi2.params import PerLineage

HERE = os.path.dirname(os.path.abspath(__file__))
TREES = os.path.join(HERE, "data_size", "trees")
STATES = os.path.join(HERE, "data_size", "states")

MASTER_SEED = 20260802
N_REPS = 100
FACTORS = [3.0, 5.0]
SIZES = [500, 1000]
BIRTH_BASE, DEATH, LOSS = 1.0, 0.3, 0.12
DRIVER, CONTROL = "driver", "control"


def make_seeds(master: int, n: int) -> list[int]:
    return [int(s) for s in np.random.default_rng(master).integers(1, 2**31 - 1, size=n)]


def one_run(factor: float, n_extant: int, seed: int, dump: str) -> dict:
    t0 = time.perf_counter()
    r = joint.simulate(
        species.birth_death(
            birth=PerLineage(BIRTH_BASE).scaled_by(
                f"genomes:{DRIVER}", {"present": factor, "absent": 1.0}),
            death=DEATH, n_extant=n_extant),
        genomes.genome(duplication=0.0, loss=LOSS, origination=0.0,
                       families=[genomes.family(DRIVER), genomes.family(CONTROL)]),
        seed=seed)
    tree = r.extant_tree
    tips = tree.extant_leaves()
    labels = tree.labels()
    drv = {i: bool(r.genome.has_family(i, DRIVER)) for i in tips}
    ctl = {i: bool(r.genome.has_family(i, CONTROL)) for i in tips}
    with open(os.path.join(TREES, dump + ".nwk"), "w") as fh:
        fh.write(tree.to_newick() + "\n")
    with open(os.path.join(STATES, dump + ".tsv"), "w") as fh:
        fh.write("tip\tdriver\tcontrol\n")
        for i in tips:
            fh.write(f"{labels[i]}\t{int(drv[i])}\t{int(ctl[i])}\n")
    ends = [tree.nodes[i].end_time for i in tips]
    return {"seed": seed, "factor": factor, "n_extant": n_extant,
            "height": max(ends),
            "prevalence_driver": sum(drv.values()) / len(tips),
            "prevalence_control": sum(ctl.values()) / len(tips),
            "seconds": round(time.perf_counter() - t0, 3)}


def main() -> int:
    os.makedirs(TREES, exist_ok=True)
    os.makedirs(STATES, exist_ok=True)
    seeds = make_seeds(MASTER_SEED, 200)[:N_REPS]
    out = {"what": "tree-size arms for the BiSSE power-vs-size question",
           "zombi2_version": zombi2.__version__, "master_seed": MASTER_SEED,
           "n_reps": N_REPS, "factors": FACTORS, "sizes": SIZES, "arms": {}}
    for n in SIZES:
        for f in FACTORS:
            t0 = time.perf_counter()
            recs = [one_run(f, n, s, f"n{n}_f{f}_r{i:03d}")
                    for i, s in enumerate(seeds)]
            out["arms"][f"n{n}_f{f}"] = recs
            drv = np.mean([r["prevalence_driver"] for r in recs])
            print(f"  n={n} f={f}: {time.perf_counter() - t0:7.1f} s, "
                  f"driver prevalence {drv:.3f}", flush=True)
    with open(os.path.join(HERE, "results_size.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    print("wrote results_size.json", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
