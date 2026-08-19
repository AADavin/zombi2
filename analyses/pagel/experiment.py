#!/usr/bin/env python3
"""A gene family and a habitat that shape each other: the runs behind the paper's
correlated-evolution example.

Four arms on matched seeds, 150 replicates each, 150 extant tips. The two connections:
the habitat multiplies the loss rate of every gene family (5x in the cave), and the eye
family's absence multiplies the rate of switching to the cave (12x). The arms switch
each connection on or off:

  feedback    both on
  trait2gen   only the habitat drives the loss rate
  gen2trait   only the eye family drives the switch rate
  null        both off (the factors are 1)

Every replicate also carries a control character from an INDEPENDENT genome run on the
same tree: one family under plain loss and transfer, connected to nothing, so any signal
Pagel's test finds on it is an artifact of the tree or the method.

Each replicate writes
  data/trees/{arm}_r{rep:03d}.nwk   the extant tree
  data/states/{arm}_r{rep:03d}.tsv  tip, habitat, eye presence, control presence

    python experiment.py
"""

from __future__ import annotations

import csv
import json
import os
import platform
import time

import zombi2
from zombi2 import joint, traits
from zombi2.genomes import family, genome as genome_spec, simulate_genomes_family
from zombi2.params import PerCopy, PerLineage
from zombi2.species import simulate_species_tree

HERE = os.path.dirname(os.path.abspath(__file__))
TREES = os.path.join(HERE, "data", "trees")
STATES = os.path.join(HERE, "data", "states")

# --------------------------------------------------------------------------- design
N_REPS = 150
N_EXTANT = 150
LOSS_BASE = 0.25            # per copy; the habitat multiplies it in the cave
LOSS_FACTOR = 5.0
SWITCH_BASE = 0.08          # surface -> cave; the eye's absence multiplies it
SWITCH_FACTOR = 12.0
SWITCH_BACK = 0.10          # cave -> surface, constant
ARMS = {"feedback": (LOSS_FACTOR, SWITCH_FACTOR),
        "trait2gen": (LOSS_FACTOR, 1.0),
        "gen2trait": (1.0, SWITCH_FACTOR),
        "null": (1.0, 1.0)}
CTRL_SEED_OFFSET = 1000     # the control run's seed, relative to the replicate's


def one(arm, L, S, seed):
    ct = simulate_species_tree(birth=1.0, n_extant=N_EXTANT, seed=seed).complete_tree
    r = joint.simulate(
        genome_spec(duplication=0.05, origination=8.0, initial_families=40,
                    loss=PerCopy(LOSS_BASE).scaled_by("trait", {"cave": L, "surface": 1.0}),
                    families=[family("eye")]),
        traits.discrete(states=["surface", "cave"], start="surface",
                        switch={"surface->cave": PerLineage(SWITCH_BASE).scaled_by(
                                    "genomes:eye", {"present": 1.0, "absent": S}),
                                "cave->surface": SWITCH_BACK}),
        tree=ct, seed=seed)
    ctrl_run = simulate_genomes_family(ct, initial_families=1, duplication=0.0,
                                       origination=0.0, loss=0.2, transfer=0.25,
                                       families=[family("ctrl")],
                                       seed=seed + CTRL_SEED_OFFSET)
    lab = ct.labels()
    eye = r.genome.family_names["eye"]
    ctrl = ctrl_run.family_names["ctrl"]
    rows = []
    for n in sorted(ct.extant_leaves()):
        rows.append({"tip": lab[n],
                     "habitat": r.trait.values[lab[n]],
                     "eye": "present" if r.genome.family_counts(n)[eye] > 0 else "absent",
                     "ctrl": "present" if ctrl_run.family_counts(n)[ctrl] > 0 else "absent"})
    return r.species.extant_tree.to_newick(), rows


def main() -> int:
    os.makedirs(TREES, exist_ok=True)
    os.makedirs(STATES, exist_ok=True)
    t0 = time.time()
    for arm, (L, S) in ARMS.items():
        for rep in range(1, N_REPS + 1):
            nwk, rows = one(arm, L, S, rep)
            with open(os.path.join(TREES, f"{arm}_r{rep:03d}.nwk"), "w") as fh:
                fh.write(nwk + "\n")
            with open(os.path.join(STATES, f"{arm}_r{rep:03d}.tsv"), "w", newline="") as fh:
                w = csv.DictWriter(fh, fieldnames=["tip", "habitat", "eye", "ctrl"],
                                   delimiter="\t")
                w.writeheader()
                w.writerows(rows)
        print(f"{arm}: {N_REPS} replicates", flush=True)
    manifest = {"zombi2": zombi2.__version__, "python": platform.python_version(),
                "elapsed_s": round(time.time() - t0, 1),
                "design": {"n_reps": N_REPS, "n_extant": N_EXTANT,
                           "loss_base": LOSS_BASE, "loss_factor": LOSS_FACTOR,
                           "switch_base": SWITCH_BASE, "switch_factor": SWITCH_FACTOR,
                           "switch_back": SWITCH_BACK, "arms": {a: list(v) for a, v in ARMS.items()}}}
    with open(os.path.join(HERE, "manifest.json"), "w") as fh:
        json.dump(manifest, fh, indent=1)
    print(f"done in {manifest['elapsed_s']}s")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
