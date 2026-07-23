"""Figures for the inversion-rate recipe, from the forward-map ``.npz``/``.json`` + a best-fit sim.

Fig 1 observables — gene-order conservation and conserved block size vs divergence time: real yeast
                    pairs vs the best-fit simulation, one column per clade.
Fig 2 ridge       — the (inversion rate × mean length) distance surface: a diagonal ridge, because
                    rate and length trade off; the star is the best cell.
Fig 3 constraints — misfit vs rate (a clear minimum: the rate is pinned) and vs length (flat: the
                    length is not) — the identifiability statement in one figure.
"""
from __future__ import annotations

import json
import pathlib
from itertools import combinations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import sim
import stats

HERE = pathlib.Path(__file__).parent
FIG = HERE / "figures"
RES = HERE / "results"
INK = "#1a1a1a"
CLADES = {
    "lachancea": dict(name="Lachancea", tree="data/lachancea/tree_myr.nwk",
                      obs="data/lachancea/observed_signed_order.json", K=8, color="#4477AA"),
    "kluyveromyces": dict(name="Kluyveromyces + Eremothecium", tree="data/kluyveromyces/tree_myr.nwk",
                          obs="data/kluyveromyces/observed_signed_order.json", K=6, color="#EE6677"),
}
plt.rcParams.update({"font.family": "sans-serif", "font.size": 11, "axes.edgecolor": INK,
                     "axes.labelcolor": INK, "text.color": INK, "xtick.color": INK,
                     "ytick.color": INK, "svg.fonttype": "none"})


def _save(fig, name):
    FIG.mkdir(exist_ok=True)
    fig.savefig(FIG / f"{name}.png", dpi=200, bbox_inches="tight")
    fig.savefig(FIG / f"{name}.svg", bbox_inches="tight")
    plt.close(fig)
    print(f"wrote figures/{name}.png")


def per_pair(genomes, divergence, core):
    """(divergence, conservation, block size) for every species pair — for the scatter."""
    pc = stats.precompute(genomes, core)
    D, U, B = [], [], []
    for a, b in combinations(genomes, 2):
        key = frozenset((a, b))
        if key not in divergence:
            continue
        uf, _ = stats.pair_stats(pc[a], pc[b])
        blk = stats.pair_block_mean(pc[a], pc[b])
        if not np.isnan(uf):
            D.append(divergence[key]); U.append(uf); B.append(blk)
    return np.array(D), np.array(U), np.array(B)


def _best_fit_pairs(clade_key, est):
    """Real and best-fit-sim per-pair observables for one clade."""
    c = CLADES[clade_key]
    div = stats.pairwise_divergence(str(HERE / c["tree"]))
    real = stats.load_signed(str(HERE / c["obs"]))
    core = stats.core_single_copy(real, 1.0)
    Dr, Ur, Br = per_pair(real, div, core)
    tree, nm = sim.load_dated_tree(str(HERE / c["tree"]))
    from fit import markers
    # the reported model: the rate read off the ridge at the reference length (a few genes)
    g = sim.simulate_signed_order(tree, nm, inversion=est["rate_at_reference_length"],
                                  mean_length=est["reference_length"],
                                  n_total=est["G"], n_chrom=c["K"], seed=12345)
    mk = markers(est["G"], est["n_core"])
    marked = {sp: [(x, f, o) for (x, f, o) in order if f in mk] for sp, order in g.items()}
    Ds, Us, Bs = per_pair(marked, div, mk)
    return (Dr, Ur, Br), (Ds, Us, Bs)


def fig_observables(ests):
    keys = list(CLADES)
    fig, axes = plt.subplots(2, len(keys), figsize=(5.2 * len(keys), 7.4), sharex="col")
    for j, key in enumerate(keys):
        c = CLADES[key]
        (Dr, Ur, Br), (Ds, Us, Bs) = _best_fit_pairs(key, ests[key])
        au, ab = axes[0, j], axes[1, j]
        au.scatter(Ds, Us, s=42, color=c["color"], alpha=0.55, edgecolors="none", zorder=2,
                   label="best-fit simulation")
        au.scatter(Dr, Ur, s=44, facecolors="white", edgecolors=INK, linewidths=1.3, zorder=3,
                   label="real genomes")
        au.set_title(c["name"], fontsize=12.5, style="italic")
        au.set_ylim(0, 1.02); au.set_xlim(left=0)
        ab.scatter(Ds, Bs, s=42, color=c["color"], alpha=0.55, edgecolors="none", zorder=2)
        ab.scatter(Dr, Br, s=44, facecolors="white", edgecolors=INK, linewidths=1.3, zorder=3)
        ab.set_yscale("log"); ab.set_xlim(left=0)
        ab.set_xlabel("divergence time (Myr)")
        for ax in (au, ab):
            for s in ("top", "right"):
                ax.spines[s].set_visible(False)
        r = ests[key]
        au.text(0.03, 0.06, f"inversion rate\n{r['rate_at_reference_length']:.1e} /gene·Myr",
                transform=au.transAxes, fontsize=10, va="bottom")
    axes[0, 0].set_ylabel("gene-order pairwise conservation")
    axes[1, 0].set_ylabel("conserved block size (genes)")
    axes[0, -1].legend(loc="upper right", frameon=False, fontsize=10)
    _save(fig, "observables")


def fig_ridge(ests):
    keys = list(CLADES)
    fig, axes = plt.subplots(1, len(keys), figsize=(6.0 * len(keys), 4.8))
    for j, key in enumerate(keys):
        d = np.load(RES / f"{key}_forward_map.npz")
        invs, Ls, dist = d["invs"], d["Ls"], d["dist"]
        logI = np.log10(invs)
        ax = axes[j]
        pcm = ax.pcolormesh(logI, Ls, dist, shading="nearest", cmap="viridis_r")
        ax.contour(logI, Ls, dist, levels=[np.nanpercentile(dist, 8)], colors="white", linewidths=1.2)
        r = ests[key]
        # mark the reported reading: the rate off the ridge at the reference length (a few genes)
        ax.plot(np.log10(r["rate_at_reference_length"]), r["reference_length"], marker="*",
                ms=18, mfc="white", mec="black", mew=1.3, zorder=5)
        ax.set_yscale("log")
        ax.set_title(CLADES[key]["name"], fontsize=12.5, style="italic")
        ax.set_xlabel(r"$\log_{10}$ inversion rate (per gene $\cdot$ Myr)")
        if j == 0:
            ax.set_ylabel("mean inversion length (genes)")
        fig.colorbar(pcm, ax=ax, label="misfit to data", fraction=0.046, pad=0.04)
        ax.text(0.03, 0.95, f"rate = {r['rate_at_reference_length']:.1e} /gene·Myr\n"
                f"(read at length = {r['reference_length']:.0f} genes)",
                transform=ax.transAxes, va="top", fontsize=9.5,
                bbox=dict(boxstyle="round", fc="white", ec=INK, alpha=0.85))
    fig.suptitle("The rate–length ridge: rate is pinned (narrow across), length is not (runs the height)",
                 fontsize=12)
    _save(fig, "ridge")


def fig_constraints(ests):
    fig, (axR, axL) = plt.subplots(1, 2, figsize=(11, 4.4))
    for key, c in CLADES.items():
        d = np.load(RES / f"{key}_forward_map.npz")
        invs, Ls, dist = d["invs"], d["Ls"], d["dist"]
        prof_rate = np.nanmin(dist, axis=0)             # best over length, at each rate
        prof_len = np.nanmin(dist, axis=1)              # best over rate, at each length
        axR.plot(np.log10(invs), prof_rate / np.nanmin(prof_rate), lw=1.8, marker="o", ms=3,
                 mfc="white", color=c["color"], label=c["name"])
        axL.plot(Ls, prof_len / np.nanmin(prof_len), lw=1.8, marker="o", ms=3,
                 mfc="white", color=c["color"], label=c["name"])
    axR.set_title("the rate is constrained", fontsize=12)
    axR.set_xlabel(r"$\log_{10}$ inversion rate (per gene $\cdot$ Myr)")
    axR.set_ylabel("misfit to data (× best fit)"); axR.set_yscale("log")
    axL.set_title("the length is not", fontsize=12)
    axL.set_xlabel("mean inversion length (genes)"); axL.set_xscale("log"); axL.set_yscale("log")
    for ax in (axR, axL):
        ax.legend(frameon=False, fontsize=10)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
    _save(fig, "constraints")


def main():
    ests = {key: json.loads((RES / f"{key}_forward_map.json").read_text()) for key in CLADES}
    fig_observables(ests)
    fig_ridge(ests)
    fig_constraints(ests)


if __name__ == "__main__":
    main()
