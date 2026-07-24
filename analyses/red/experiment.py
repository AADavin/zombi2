"""The RED validation experiment on the ZOMBI2 clean core.

Question: does Relative Evolutionary Divergence recover a tree's relative node ages once a realistic
molecular clock distorts branch lengths into substitutions? True ages are unknowable on a real
tree, so we borrow one model-free number from the real GTDB archaeal phylogram — root-to-tip
substitution CV = 0.2315 (``observable.py``) — reproduce that raggedness in ZOMBI2 simulations
where truth *is* known, and grade RED there.

Forward model. Simulate a dated species tree (truth); evolve it under ZOMBI2's shipped relaxed
clock — ``substitution = 1.0 * ByLineage(spread=sigma, dist=...)``, the uncorrelated lineage clock
— and read the resulting ``species_phylogram`` (branch lengths in substitutions/site). RED of the
dated tree is the ground truth (RED is exact on an ultrametric tree); RED of the ragged phylogram
is the estimate. Sweep ``sigma`` to trace accuracy (Pearson r, nRMSE) against the realized
root-to-tip CV, and read off the accuracy at the CV real archaea actually show.

Shipped-clock note. The clean core wires exactly one sequence-level clock — the uncorrelated
``ByLineage`` — with two tails, ``lognormal`` and ``gamma``. (The autocorrelated ``FromParent`` is a
species-level modifier, not wired at the sequence level; the legacy six-clock sweep is narrowed to
what ships.)
"""
from __future__ import annotations

import json
import pathlib

import numpy as np
from scipy.stats import pearsonr

from zombi2 import genomes, sequences, species
from zombi2.rates import modifiers as mod
from zombi2.sequences import substitution_models as sm
from zombi2.tree import read_newick

from red import internal_nodes
from red import relative_evolutionary_divergence as red_of

TARGET_CV = 0.2315          # measured GTDB archaeal root-to-tip substitution CV (10,122 tips)
HERE = pathlib.Path(__file__).parent
N_EXTANT = 400
REPS = 8
SPREADS = np.round(np.arange(0.0, 2.001, 0.1), 3)
DISTS = ("lognormal", "gamma")


def _preorder(tree) -> list[int]:
    order: list[int] = []
    stack = [tree.root]
    while stack:
        i = stack.pop()
        order.append(i)
        k = tree.nodes[i].children
        if k is not None:
            stack.extend(k)
    return order


def rtt_cv(tree) -> float:
    """Root-to-tip branch-length CV of a tree/phylogram."""
    depth: dict[int, float] = {}
    tips: list[float] = []
    for i in _preorder(tree):
        nd = tree.nodes[i]
        depth[i] = 0.0 if nd.parent is None else depth[nd.parent] + (nd.end_time - nd.birth_time)
        if nd.children is None:
            tips.append(depth[i])
    d = np.asarray(tips)
    return float(d.std() / d.mean())


def make_tree(n_extant: int, seed: int, birth: float = 1.0, death: float = 0.0):
    """A dated species tree (the truth) plus a small genome to carry the clock down it."""
    res = species.simulate_species_tree(birth=birth, death=death, n_extant=n_extant, seed=seed)
    g = genomes.simulate_genomes_unordered(res.complete_tree, duplication=0.01, loss=0.01,
                                           initial_families=5, seed=seed)
    return res.extant_tree, g


def grade(dated, red_true, g, dist: str, spread: float, seed: int):
    """Evolve the dated tree under a ByLineage clock; return (realized CV, r, nRMSE, (true, est))."""
    if spread == 0.0:                                   # strict clock: phylogram ∝ dated tree
        ids = internal_nodes(dated)
        t = np.array([red_true[i] for i in ids])
        return 0.0, 1.0, 0.0, (t, t.copy())
    seqres = sequences.simulate_sequences(
        g, model=sm.jc69(), length=1,
        substitution=1.0 * mod.ByLineage(spread=spread, dist=dist), seed=seed)
    phylo, _ = read_newick(seqres.species_phylogram["extant"])
    red_est = red_of(phylo)
    ids = [i for i in internal_nodes(dated) if i in red_est]
    t = np.array([red_true[i] for i in ids])
    e = np.array([red_est[i] for i in ids])
    r = float(pearsonr(t, e)[0])
    nrmse = float(np.sqrt(np.mean((e - t) ** 2)))
    return rtt_cv(phylo), r, nrmse, (t, e)


def _crossing(x: np.ndarray, y: np.ndarray, target: float) -> float:
    """Linear-interpolated x where a monotone curve y(x) hits target."""
    for i in range(1, len(y)):
        if (y[i - 1] - target) * (y[i] - target) <= 0 and y[i] != y[i - 1]:
            f = (target - y[i - 1]) / (y[i] - y[i - 1])
            return float(x[i - 1] + f * (x[i] - x[i - 1]))
    return float("nan")


def sweep() -> dict:
    trees = [make_tree(N_EXTANT, seed=100 + k) for k in range(REPS)]
    red_trues = [red_of(dated) for dated, _ in trees]
    out: dict = {"target_cv": TARGET_CV, "n_extant": N_EXTANT, "reps": REPS,
                 "spreads": [float(s) for s in SPREADS], "families": {}}
    for dist in DISTS:
        rows = []
        for si, s in enumerate(SPREADS):
            cvs, rs, nes = [], [], []
            for k, (dated, g) in enumerate(trees):
                cv, r, ne, _ = grade(dated, red_trues[k], g, dist, float(s), seed=7000 + si * 20 + k)
                cvs.append(cv); rs.append(r); nes.append(ne)
            rows.append({"spread": float(s),
                         "cv": float(np.mean(cvs)), "cv_sd": float(np.std(cvs)),
                         "r": float(np.mean(rs)), "r_sd": float(np.std(rs)),
                         "nrmse": float(np.mean(nes)), "nrmse_sd": float(np.std(nes))})
        cvarr = np.array([row["cv"] for row in rows])
        rec_spread = _crossing(SPREADS, cvarr, TARGET_CV)
        # accuracy AT the target CV: interpolate r / nRMSE as functions of CV (CV rises with sigma)
        r_at = float(np.interp(TARGET_CV, cvarr, np.array([row["r"] for row in rows])))
        nrmse_at = float(np.interp(TARGET_CV, cvarr, np.array([row["nrmse"] for row in rows])))
        out["families"][dist] = {"rows": rows, "recovered_spread": rec_spread,
                                 "r_at_target": r_at, "nrmse_at_target": nrmse_at,
                                 "cv_max": float(cvarr.max())}
        print(f"  {dist:10s}: recovered sigma={rec_spread:.3f} at CV={TARGET_CV}"
              f" -> r={r_at:.4f}, nRMSE={nrmse_at*100:.1f}%  (CV range 0..{cvarr.max():.2f})")
    return out


def scatter(sweep_rows, levels=(0.10, TARGET_CV, 0.40), dist="lognormal", n_extant=500, seed=42) -> dict:
    """True vs RED-recovered relative ages on one tree at three raggedness levels. The sigma for each
    level is read from the (stable, multi-tree) sweep's CV(sigma) map, so the realized CV lands near
    the target; each panel is labelled by its realized CV."""
    dated, g = make_tree(n_extant, seed=seed)
    red_true = red_of(dated)
    sig_grid = np.array([r["spread"] for r in sweep_rows])
    cv_grid = np.array([r["cv"] for r in sweep_rows])
    panels = []
    for lvl in levels:
        sig = float(np.interp(lvl, cv_grid, sig_grid))
        cv, r, ne, (t, e) = grade(dated, red_true, g, dist, sig, seed=1234)
        panels.append({"target_cv": lvl, "sigma": sig, "cv": cv, "r": r, "nrmse": ne,
                       "true": [float(x) for x in t], "est": [float(x) for x in e]})
        print(f"  scatter CV~{lvl:.2f}: sigma={sig:.3f} realized_CV={cv:.3f} r={r:.4f} nRMSE={ne*100:.1f}%")
    return {"dist": dist, "n_extant": n_extant, "panels": panels}


def main():
    print(f"GTDB target root-to-tip CV = {TARGET_CV}")
    print("Bridge sweep (RED accuracy vs raggedness):")
    res = sweep()
    print("Scatter (true vs recovered ages):")
    sc = scatter(res["families"]["lognormal"]["rows"])
    res["scatter"] = sc
    (HERE / "results.json").write_text(json.dumps(res, indent=2))
    print(f"wrote {HERE / 'results.json'}")


if __name__ == "__main__":
    main()
