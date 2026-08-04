"""Figures for the RED validation recipe, from ``results.json`` + the GTDB observable.

Fig 1 observable   — what root-to-tip variation is (schematic), and how much of it real archaea show.
Fig 2 clock_recovery — root-to-tip CV vs clock sigma, one curve per clock; the GTDB target crossing.
Fig 3 red_bridge   — RED accuracy (Pearson r, nRMSE) vs root-to-tip CV; the read-off at CV = 0.2315.
Fig 4 red_scatter  — true vs RED-recovered relative ages below, at and above the real archaeal value.

Every figure names the quantity the same way: **root-to-tip variation**, measured as the CV of
root-to-tip distances. One name, on the axes and in the prose.
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
# Paul Tol 'bright'. The two uncorrelated tails are cool, the autocorrelated clock warm — the
# comparison the figures exist to make is uncorrelated vs autocorrelated, not tail vs tail.
COLORS = {"lognormal": "#4477AA", "gamma": "#66CCEE", "autocorrelated": "#EE6677"}
LABELS = {"lognormal": "uncorrelated, lognormal", "gamma": "uncorrelated, gamma",
          "autocorrelated": "autocorrelated"}
INK = "#1a1a1a"
plt.rcParams.update({"font.family": "sans-serif", "font.size": 11, "axes.edgecolor": INK,
                     "axes.labelcolor": INK, "text.color": INK, "xtick.color": INK,
                     "ytick.color": INK, "svg.fonttype": "none"})


#: the docs site carries this study as a worked example, so the figures are written there too rather
#: than copied by hand — one command regenerates both, and the two cannot drift apart.
DOCS = HERE.parents[1] / "docs" / "assets" / "red"


def _save(fig, name):
    FIG.mkdir(exist_ok=True)
    fig.savefig(FIG / f"{name}.png", dpi=200, bbox_inches="tight")
    fig.savefig(FIG / f"{name}.svg", bbox_inches="tight")
    if DOCS.parent.parent.is_dir():
        DOCS.mkdir(parents=True, exist_ok=True)
        fig.savefig(DOCS / f"{name}.png", dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote figures/{name}.png")


def _cladogram(ax, tips_x, y0, *, split=0.30, sub=0.62, color=INK, lw=1.5):
    """A 4-tip cladogram with topology ((0,1),(2,3)), rooted at x = 0.

    ``tips_x`` is each tip's distance from the root, so passing four equal values draws a strict
    clock and four unequal ones draws a relaxed clock. Returns the tips' y positions, top to bottom.
    """
    y = [y0 + 3, y0 + 2, y0 + 1, y0]
    ya, yb = (y[0] + y[1]) / 2, (y[2] + y[3]) / 2
    ax.plot([0, split], [(ya + yb) / 2] * 2, color=color, lw=lw, solid_capstyle="round")
    ax.plot([split] * 2, [ya, yb], color=color, lw=lw, solid_capstyle="round")
    for parent_y, pair in ((ya, (0, 1)), (yb, (2, 3))):
        ax.plot([split, sub], [parent_y] * 2, color=color, lw=lw, solid_capstyle="round")
        ax.plot([sub] * 2, [y[pair[0]], y[pair[1]]], color=color, lw=lw, solid_capstyle="round")
        for i in pair:
            ax.plot([sub, tips_x[i]], [y[i]] * 2, color=color, lw=lw, solid_capstyle="round")
            ax.plot([tips_x[i]], [y[i]], "o", ms=5, color=color, zorder=3)
    return y


def _schematic(ax):
    """What root-to-tip variation *is*, before the histogram says how much of it archaea have.

    The same four-tip tree under a strict clock (every tip the same distance from the root) and a
    relaxed one (rates differ, so the tips end at different distances). The measurement is the
    spread of those distances, and nothing about it needs a date or a rate model.
    """
    strict = [1.0] * 4
    relaxed = [0.78, 1.02, 1.24, 0.90]
    _cladogram(ax, strict, 5.2)
    tips = _cladogram(ax, relaxed, 0.0)
    rule = 1.40
    ax.plot([1.0, 1.0], [4.9, 8.5], color="#888888", lw=1.0, ls="--", zorder=1)
    for x, y in zip(relaxed, tips):                       # leaders out to a common rule
        ax.plot([x, rule], [y, y], color="#888888", lw=0.8, ls=":", zorder=1)
    ax.annotate("", xy=(rule, min(tips)), xytext=(rule, max(tips)),
                arrowprops=dict(arrowstyle="<->", color="#888888", lw=1.0))
    ax.text(0.0, 9.0, "strict clock: tips aligned", fontsize=11, va="bottom", fontweight="semibold")
    ax.text(0.0, 3.8, "relaxed clock: tips spread out", fontsize=11, va="bottom",
            fontweight="semibold")
    ax.text(rule + 0.08, sum(tips) / 4, "what we\nmeasure", fontsize=10, va="center",
            color="#555555", linespacing=1.35)
    ax.set_xlim(-0.05, 2.15)
    ax.set_ylim(-0.6, 9.9)
    ax.set_axis_off()


def fig_observable(res):
    d = root_to_tip_depths(HERE / "data" / "ar53.tree")
    fig, (axl, axr) = plt.subplots(1, 2, figsize=(11.4, 4.0),
                                   gridspec_kw={"width_ratios": [1.0, 1.25], "wspace": 0.22})
    _schematic(axl)
    axr.hist(d, bins=70, color="#BBBBBB", edgecolor="white", linewidth=0.3)
    axr.axvline(d.mean(), color=INK, lw=1.5, ls="--")
    axr.set_xlabel("root-to-tip distance (substitutions/site)")
    axr.set_ylabel("archaeal genomes")
    axr.text(0.97, 0.95, f"GTDB archaea\n{len(d):,} tips\nmean {d.mean():.2f}\n"
             f"CV = {cv_of(d):.3f}", transform=axr.transAxes, ha="right", va="top",
             fontsize=11, bbox=dict(boxstyle="round", fc="white", ec=INK, alpha=0.9))
    for s in ("top", "right"):
        axr.spines[s].set_visible(False)
    _save(fig, "observable")


def fig_clock_recovery(res):
    fig, ax = plt.subplots(figsize=(6.2, 4.2))
    tgt = res["target_cv"]
    # The crossings of the two uncorrelated clocks are ~0.05 apart in σ, so their labels are
    # staggered in depth rather than laid side by side, where they would overprint.
    for k, (dist, fam) in enumerate(res["families"].items()):
        spreads = np.array(fam["spreads"])
        cv = np.array([r["cv"] for r in fam["rows"]])
        sd = np.array([r["cv_sd"] for r in fam["rows"]])
        ax.plot(spreads, cv, "-o", ms=3, color=COLORS[dist], label=LABELS[dist])
        ax.fill_between(spreads, cv - sd, cv + sd, color=COLORS[dist], alpha=0.15, lw=0)
        rec = fam["recovered_spread"]
        ax.plot([rec], [tgt], "o", color=COLORS[dist], ms=8, mec=INK, mew=1, zorder=5)
        ax.annotate(f"σ = {rec:.2f}", xy=(rec, tgt), xytext=(rec, tgt - 0.055 - 0.045 * k),
                    ha="center", fontsize=9.5, color=COLORS[dist])
    ax.axhline(tgt, color=INK, lw=1.2, ls="--")
    # The sweeps run to σ = 2, far past every crossing; the axis stops where the read-off is legible.
    ax.set_xlim(0, 1.0)
    ax.set_ylim(0, 0.7)
    ax.text(0.02, tgt, f"real archaea, CV = {tgt:.3f}", va="bottom", ha="left", fontsize=10)
    ax.set_xlabel("clock heterogeneity  σ")
    ax.set_ylabel("root-to-tip variation  (CV)")
    ax.set_title("Which σ reproduces the real spread", fontsize=12)
    ax.legend(loc="upper left", frameon=False, title="clock")
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
        axes[0].plot(cv, r, "-o", ms=3, color=COLORS[dist], label=LABELS[dist])
        axes[0].fill_between(cv, r - rsd, r + rsd, color=COLORS[dist], alpha=0.15, lw=0)
        axes[1].plot(cv, ne, "-o", ms=3, color=COLORS[dist], label=LABELS[dist])
        axes[1].fill_between(cv, ne - nesd, ne + nesd, color=COLORS[dist], alpha=0.15, lw=0)
        axes[0].plot([tgt], [fam["r_at_target"]], "o", color=COLORS[dist], ms=8, mec=INK, mew=1, zorder=5)
        axes[1].plot([tgt], [fam["nrmse_at_target"] * 100], "o", color=COLORS[dist], ms=8, mec=INK, mew=1, zorder=5)
    axes[0].set_ylabel("RED–age Pearson r")
    axes[1].set_ylabel("RED–age nRMSE (% of tree depth)")
    axes[0].set_title("How well RED recovers node ages", fontsize=12)
    axes[1].set_title("... and by how much it is off", fontsize=12)
    axes[0].set_ylim(0.6, 1.02)
    axes[1].set_ylim(0, 22)
    for ax in axes:
        # The sweeps run to CV ~ 1.9, but past ~0.8 the autocorrelated clock's between-tree variance
        # swamps the signal at 8 replicates — the curve there is noise, not structure. The claim this
        # figure makes lives at CV = 0.232, so the axis stops where the curves are still readable.
        ax.set_xlim(0, 0.8)
        # Everything to the left of the dashed line is real archaea or milder, which is the half of
        # the figure the conclusion is read off; shading it says so without a sentence.
        ax.axvspan(0, tgt, color="#EE6677", alpha=0.07, lw=0, zorder=0)
        ax.axvline(tgt, color=INK, lw=1.2, ls="--")
        ax.annotate(f"real archaea (CV = {tgt:.3f})", xy=(tgt, ax.get_ylim()[1]),
                    xytext=(tgt + 0.03, ax.get_ylim()[1]), va="top", ha="left", fontsize=9,
                    color=INK)
        ax.set_xlabel("root-to-tip variation  (CV)")
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
    axes[0].legend(loc="lower left", frameon=False, title="clock", fontsize=9)
    _save(fig, "red_bridge")


def fig_red_scatter(res):
    panels = res["scatter"]["panels"]
    tgt = res["target_cv"]
    # Which panel is the real archaeal one: the middle panel is chosen to sit at the target CV, and
    # naming it "real archaea" rather than a number is what the reader needs to know it is the case
    # the study is actually about. Highlighting it in the accent colour says the same thing twice.
    real = min(range(len(panels)), key=lambda i: abs(panels[i]["cv"] - tgt))
    where = ["milder than real archaea", "real archaea", "more extreme than real archaea"]
    fig, axes = plt.subplots(1, len(panels), figsize=(4.0 * len(panels), 4.2))
    for i, (ax, p) in enumerate(zip(axes, panels)):
        accent = i == real
        colour = COLORS["autocorrelated"] if accent else COLORS[res["scatter"]["dist"]]
        ax.plot([0, 1], [0, 1], color=INK, lw=1, ls="--", zorder=1)
        ax.scatter(p["true"], p["est"], s=10, color=colour, alpha=0.6,
                   edgecolors="none", zorder=2)
        ax.set_title(f"{where[i]}\nCV = {p['cv']:.2f}   (σ = {p['sigma']:.2f})", fontsize=11,
                     color=colour if accent else INK,
                     fontweight="bold" if accent else "normal", linespacing=1.5)
        ax.text(0.05, 0.95, f"r = {p['r']:.3f}\nnRMSE = {p['nrmse']*100:.1f}%",
                transform=ax.transAxes, va="top", fontsize=10)
        ax.set_xlabel("true relative age")
        ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.set_aspect("equal")
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
    axes[0].set_ylabel("RED-recovered relative age")
    fig.suptitle(f"True vs RED-recovered node ages "
                 f"({LABELS[res['scatter']['dist']]} clock, {res['scatter']['n_extant']} tips)",
                 fontsize=12)
    _save(fig, "red_scatter")


def main():
    res = json.loads((HERE / "results.json").read_text())
    fig_observable(res)
    fig_clock_recovery(res)
    fig_red_bridge(res)
    fig_red_scatter(res)


if __name__ == "__main__":
    main()
