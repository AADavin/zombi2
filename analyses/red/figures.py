"""Figures for the RED validation recipe, from ``results.json`` + the GTDB observable.

Fig 1 observable   — the raggedness real archaea show (GTDB root-to-tip substitution histogram).
Fig 2 clock_recovery — root-to-tip CV vs clock sigma, one curve per tail; the GTDB target crossing.
Fig 3 red_bridge   — RED accuracy (Pearson r, nRMSE) vs realized CV; the read-off at CV = 0.2315.
Fig 4 red_scatter  — true vs RED-recovered relative ages at three raggedness levels.
"""
from __future__ import annotations

import json
import pathlib

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from observable import cv as cv_of
from observable import root_to_tip_depths

HERE = pathlib.Path(__file__).parent
FIG = HERE / "figures"
COLORS = {"lognormal": "#4477AA", "gamma": "#EE6677"}   # Paul Tol 'bright'
INK = "#1a1a1a"
plt.rcParams.update({"font.family": "sans-serif", "font.size": 11, "axes.edgecolor": INK,
                     "axes.labelcolor": INK, "text.color": INK, "xtick.color": INK,
                     "ytick.color": INK, "svg.fonttype": "none"})


def _save(fig, name):
    FIG.mkdir(exist_ok=True)
    fig.savefig(FIG / f"{name}.png", dpi=200, bbox_inches="tight")
    fig.savefig(FIG / f"{name}.svg", bbox_inches="tight")
    plt.close(fig)
    print(f"wrote figures/{name}.png")


def fig_observable(res):
    d = root_to_tip_depths(HERE / "data" / "ar53.tree")
    fig, ax = plt.subplots(figsize=(6.2, 4.0))
    ax.hist(d, bins=70, color="#BBBBBB", edgecolor="white", linewidth=0.3)
    ax.axvline(d.mean(), color=INK, lw=1.5, ls="--")
    ax.set_xlabel("root-to-tip distance (substitutions/site)")
    ax.set_ylabel("archaeal genomes")
    ax.text(0.97, 0.95, f"GTDB archaea\n{len(d):,} tips\nmean {d.mean():.2f}\n"
            f"CV = {cv_of(d):.3f}", transform=ax.transAxes, ha="right", va="top",
            fontsize=11, bbox=dict(boxstyle="round", fc="white", ec=INK, alpha=0.9))
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    _save(fig, "observable")


def fig_clock_recovery(res):
    fig, ax = plt.subplots(figsize=(6.2, 4.2))
    spreads = np.array(res["spreads"])
    tgt = res["target_cv"]
    for dist, fam in res["families"].items():
        cv = np.array([r["cv"] for r in fam["rows"]])
        sd = np.array([r["cv_sd"] for r in fam["rows"]])
        ax.plot(spreads, cv, "-o", ms=3, color=COLORS[dist], label=f"{dist}")
        ax.fill_between(spreads, cv - sd, cv + sd, color=COLORS[dist], alpha=0.15, lw=0)
        rec = fam["recovered_spread"]
        ax.plot([rec], [tgt], "o", color=COLORS[dist], ms=8, mec=INK, mew=1, zorder=5)
    ax.axhline(tgt, color=INK, lw=1.2, ls="--")
    ax.text(spreads[-1], tgt, f"  GTDB CV = {tgt:.3f}", va="bottom", ha="right", fontsize=10)
    ax.set_xlabel("clock heterogeneity  σ  (ByLineage spread)")
    ax.set_ylabel("root-to-tip substitution CV")
    ax.set_title("Calibrating the clock to real raggedness", fontsize=12)
    ax.legend(loc="upper left", frameon=False, title="uncorrelated tail")
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    _save(fig, "clock_recovery")


def fig_red_bridge(res):
    fig, axes = plt.subplots(1, 2, figsize=(10.2, 4.2))
    tgt = res["target_cv"]
    for dist, fam in res["families"].items():
        cv = np.array([r["cv"] for r in fam["rows"]])
        r = np.array([r["r"] for r in fam["rows"]])
        rsd = np.array([r_["r_sd"] for r_ in fam["rows"]])
        ne = np.array([r_["nrmse"] for r_ in fam["rows"]]) * 100
        nesd = np.array([r_["nrmse_sd"] for r_ in fam["rows"]]) * 100
        axes[0].plot(cv, r, "-o", ms=3, color=COLORS[dist], label=dist)
        axes[0].fill_between(cv, r - rsd, r + rsd, color=COLORS[dist], alpha=0.15, lw=0)
        axes[1].plot(cv, ne, "-o", ms=3, color=COLORS[dist], label=dist)
        axes[1].fill_between(cv, ne - nesd, ne + nesd, color=COLORS[dist], alpha=0.15, lw=0)
        axes[0].plot([tgt], [fam["r_at_target"]], "o", color=COLORS[dist], ms=8, mec=INK, mew=1, zorder=5)
        axes[1].plot([tgt], [fam["nrmse_at_target"] * 100], "o", color=COLORS[dist], ms=8, mec=INK, mew=1, zorder=5)
    for ax in axes:
        ax.axvline(tgt, color=INK, lw=1.2, ls="--")
        ax.annotate("real archaea\n(CV = 0.232)", xy=(tgt, ax.get_ylim()[1]),
                    xytext=(tgt + 0.12, ax.get_ylim()[1]), va="top", ha="left", fontsize=9,
                    color=INK)
        ax.set_xlabel("root-to-tip substitution CV")
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
    axes[0].set_ylabel("RED–age Pearson r")
    axes[1].set_ylabel("RED–age nRMSE (% of tree depth)")
    axes[0].set_title("How faithfully RED recovers ages", fontsize=12)
    axes[1].set_title("... and by how much it errs", fontsize=12)
    axes[0].legend(loc="lower left", frameon=False, title="uncorrelated tail")
    _save(fig, "red_bridge")


def fig_red_scatter(res):
    panels = res["scatter"]["panels"]
    fig, axes = plt.subplots(1, len(panels), figsize=(4.0 * len(panels), 4.0))
    for ax, p in zip(axes, panels):
        ax.plot([0, 1], [0, 1], color=INK, lw=1, ls="--", zorder=1)
        ax.scatter(p["true"], p["est"], s=10, color=COLORS[res["scatter"]["dist"]], alpha=0.6,
                   edgecolors="none", zorder=2)
        ax.set_title(f"CV = {p['cv']:.2f}   (σ = {p['sigma']:.2f})", fontsize=11)
        ax.text(0.05, 0.95, f"r = {p['r']:.3f}\nnRMSE = {p['nrmse']*100:.1f}%",
                transform=ax.transAxes, va="top", fontsize=10)
        ax.set_xlabel("true relative age")
        ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.set_aspect("equal")
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
    axes[0].set_ylabel("RED-recovered relative age")
    fig.suptitle(f"True vs RED-recovered node ages "
                 f"({res['scatter']['dist']} clock, {res['scatter']['n_extant']} tips)", fontsize=12)
    _save(fig, "red_scatter")


def main():
    res = json.loads((HERE / "results.json").read_text())
    fig_observable(res)
    fig_clock_recovery(res)
    fig_red_bridge(res)
    fig_red_scatter(res)


if __name__ == "__main__":
    main()
