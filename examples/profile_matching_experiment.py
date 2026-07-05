"""A full profile-matching experiment on a small species tree.

We (1) fix a small species tree and a set of *true* gene-family rates, (2) simulate an
empirical dataset — a whole distribution of gene families, each with its own gene tree,
summarised as a copy-number profile — and (3) run ``match_profiles`` to recover the rates
from the profile alone. Then we (4) repeat the fit over many independent empirical
datasets from the *same* truth to check calibration (does the posterior bracket the truth?).

Run: ``python examples/profile_matching_experiment.py [out.png]``  (needs the Rust fast
path for the sweep; matplotlib for the figure).
"""

import sys

import numpy as np

import zombi2 as z
from zombi2 import matching  # ABC inference is withheld from the public v1 API

TRUTH = dict(duplication=0.30, transfer=0.10, loss=0.50, origination=1.50)
PRIORS = {"duplication": (0, 1.0), "transfer": (0, 0.5),
          "loss": (0, 1.2), "origination": (0, 3.0)}
INITIAL_SIZE = 15
OUT = sys.argv[1] if len(sys.argv) > 1 else "profile_matching.png"


def main():
    # (1) a small species tree
    tree = z.simulate_species_tree(z.BirthDeath(0.9, 0.25), n_tips=18, age=5.0, seed=0)
    print(f"species tree: {len(tree.extant_leaves())} extant tips\n")

    # (2) the empirical dataset: a distribution of gene families / gene trees along the tree
    emp_genomes = z.simulate_genomes(tree, **TRUTH, initial_families=INITIAL_SIZE, seed=1)
    profile = emp_genomes.profiles
    trees = emp_genomes.gene_trees()
    surviving = {f: t for f, (c, t) in trees.items() if t}
    print(f"empirical: {profile.matrix.shape[0]} gene families "
          f"({len(surviving)} with surviving copies), "
          f"max copy number {profile.matrix.max()}")
    example = sorted(surviving)[:1]
    for f in example:
        nwk = surviving[f]
        print(f"  e.g. family {f} gene tree: "
              f"{nwk[:70]}{'...' if len(nwk) > 70 else ''}")

    # (3) recover the rates from the profile alone
    fit = matching.match_profiles(tree, profile, priors=PRIORS, n_sims=4000, accept=0.02,
                           initial_families=INITIAL_SIZE, seed=1)
    print("\nrecovered posterior (median [95% CI])  vs  truth:")
    s = fit.summary()
    for p in z.matching.RATE_PARAMS:
        cover = "in" if s[p]["lo95"] <= TRUTH[p] <= s[p]["hi95"] else "OUT"
        print(f"  {p:<12} {s[p]['median']:.2f} [{s[p]['lo95']:.2f}, {s[p]['hi95']:.2f}]"
              f"   truth {TRUTH[p]:.2f}  ({cover} 95% CI)")

    post = fit.posterior
    ridge = np.corrcoef(post["origination"], post["loss"])[0, 1]
    print(f"\nidentifiability ridge  corr(origination, loss) = {ridge:+.2f}")

    # (4) calibration sweep: refit on many independent empirical datasets from the same truth
    n_rep = 15
    print(f"\ncalibration over {n_rep} independent empirical datasets (same truth):")
    medians = {p: [] for p in TRUTH}
    covered = {p: 0 for p in TRUTH}
    for r in range(n_rep):
        emp_r = z.simulate_genomes(tree, **TRUTH, initial_families=INITIAL_SIZE, seed=100 + r,
                                   output="profiles")
        fit_r = matching.match_profiles(tree, emp_r, priors=PRIORS, n_sims=2000, accept=0.03, seed=7)
        sr = fit_r.summary()
        for p in TRUTH:
            medians[p].append(sr[p]["median"])
            covered[p] += int(sr[p]["lo95"] <= TRUTH[p] <= sr[p]["hi95"])
    for p in z.matching.RATE_PARAMS:
        med = np.array(medians[p])
        print(f"  {p:<12} median {med.mean():.2f} +/- {med.std():.2f}   truth {TRUTH[p]:.2f}"
              f"   95%-CI coverage {covered[p]}/{n_rep}")

    _figure(fit, medians, OUT)
    print(f"\nsaved figure -> {OUT}")


def _figure(fit, medians, out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    params = list(z.matching.RATE_PARAMS)
    hi = {p: PRIORS[p][1] for p in params}
    post = fit.posterior

    fig, ax = plt.subplots(2, 2, figsize=(11, 8))

    # A: frequency-spectrum posterior-predictive check
    fit.plot_spectra(ax=ax[0, 0])

    # B: the identifiability ridge (origination vs loss), all draws + accepted
    a = ax[0, 1]
    a.scatter(fit.samples[:, params.index("origination")],
              fit.samples[:, params.index("loss")], s=6, c="lightgray", label="all draws")
    a.scatter(post["origination"], post["loss"], s=14, c="C3", label="accepted")
    a.axvline(TRUTH["origination"], color="k", ls="--", lw=1)
    a.axhline(TRUTH["loss"], color="k", ls="--", lw=1)
    a.set_xlabel("origination"); a.set_ylabel("loss")
    a.set_title("Identifiability ridge"); a.legend(loc="upper right", fontsize=8)

    # C: posterior for all four rates, normalised to their prior range; truth in red
    a = ax[1, 0]
    data = [post[p] / hi[p] for p in params]
    a.boxplot(data, tick_labels=params, showfliers=False)
    a.scatter(range(1, len(params) + 1), [TRUTH[p] / hi[p] for p in params],
              c="C3", zorder=3, label="truth")
    a.set_ylim(0, 1); a.set_ylabel("value / prior upper bound")
    a.set_title("Posterior vs truth (all rates)"); a.legend(fontsize=8)

    # D: calibration — recovered medians across replicates, normalised; truth in red
    a = ax[1, 1]
    for i, p in enumerate(params, start=1):
        med = np.array(medians[p]) / hi[p]
        a.scatter(np.full_like(med, i) + np.linspace(-0.15, 0.15, len(med)), med,
                  s=12, c="C0", alpha=0.6)
        a.scatter([i], [TRUTH[p] / hi[p]], c="C3", zorder=3)
    a.set_xticks(range(1, len(params) + 1)); a.set_xticklabels(params)
    a.set_ylim(0, 1); a.set_ylabel("recovered median / prior upper bound")
    a.set_title("Calibration across replicates")

    fig.suptitle("Profile matching on a small tree: recover D/T/L/O from a copy-number profile",
                 fontsize=12)
    fig.tight_layout()
    fig.savefig(out, dpi=90)


if __name__ == "__main__":
    main()
