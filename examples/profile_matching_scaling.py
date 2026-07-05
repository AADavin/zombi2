"""How profile matching improves with more data (bigger tree -> more gene families).

The identifiability ridge that biases the loss estimate on a small tree is a
*data* limitation, not a method bug. Here we hold the true rates fixed and refit on
progressively larger species trees, tracking the recovered loss (the hardest rate) and
origination (the easiest) across independent replicates. Loss bias and interval width
shrink as the tree grows.

Run: ``python examples/profile_matching_scaling.py [out.png]``  (needs the Rust engine).
"""

import sys

import numpy as np

import zombi2 as z

TRUTH = dict(duplication=0.30, transfer=0.10, loss=0.50, origination=1.50)
PRIORS = {"duplication": (0, 1.0), "transfer": (0, 0.5), "loss": (0, 1.2), "origination": (0, 3.0)}
SIZES = [15, 30, 60, 120]
N_REP = 20
OUT = sys.argv[1] if len(sys.argv) > 1 else "profile_matching_scaling.png"


def calibrate(n_tips):
    tree = z.simulate_species_tree(z.BirthDeath(0.9, 0.25), n_tips=n_tips, age=5.0, seed=0)
    rows = {p: {"err": [], "width": [], "cov": 0} for p in TRUTH}
    n_fam = []
    for r in range(N_REP):
        emp = z.simulate_genomes(tree, initial_families=15, seed=100 + r, output="profiles", **TRUTH)
        n_fam.append(emp.matrix.shape[0])
        fit = z.match_profiles(tree, emp, priors=PRIORS, n_sims=2000, accept=0.03, seed=7)
        s = fit.summary()
        for p in TRUTH:
            rows[p]["err"].append(s[p]["median"] - TRUTH[p])
            rows[p]["width"].append(s[p]["hi95"] - s[p]["lo95"])
            rows[p]["cov"] += int(s[p]["lo95"] <= TRUTH[p] <= s[p]["hi95"])
    return np.mean(n_fam), rows


def main():
    print(f"{'tips':>5}{'families':>10} | "
          + " | ".join(f"{p[:4]} bias  cov  width" for p in z.matching.RATE_PARAMS))
    results = {}
    for n in SIZES:
        fam, rows = calibrate(n)
        results[n] = rows
        cells = []
        for p in z.matching.RATE_PARAMS:
            bias = np.mean(rows[p]["err"]); width = np.mean(rows[p]["width"]); cov = rows[p]["cov"]
            cells.append(f"{bias:+.2f} {cov:2d}/{N_REP} {width:.2f}")
        print(f"{n:>5}{fam:>10.0f} | " + " | ".join(cells))

    _figure(results, OUT)
    print(f"\nsaved figure -> {OUT}")


def _figure(results, out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    for p, a in zip(("loss", "origination"), ax):
        for n in SIZES:
            err = np.array(results[n][p]["err"]) + TRUTH[p]  # recovered medians
            a.scatter(np.full_like(err, n) + np.random.default_rng(n).uniform(-2, 2, len(err)),
                      err, s=12, c="C0", alpha=0.6)
        a.axhline(TRUTH[p], color="C3", lw=1.5, label="truth")
        a.set_xlabel("species-tree tips"); a.set_ylabel(f"recovered {p} median")
        a.set_title(f"{p}: recovery vs tree size"); a.legend()
    fig.suptitle("More data (bigger tree -> more families) shrinks bias and interval width")
    fig.tight_layout()
    fig.savefig(out, dpi=90)


if __name__ == "__main__":
    main()
