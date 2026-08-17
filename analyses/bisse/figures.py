"""The three-panel figure for the BiSSE analysis, from ``results.json`` + ``verdicts.json``.

Panel A  prevalence — driver and control prevalence among extant tips against the driver
         factor, with 95% bootstrap CIs; the gap between the curves is the within-run signal.
Panel B  decomposition — the headline factor f = 3 against its matched null as a
         three-column slopegraph, every replicate joined by its seed; the raw gap splits
         exactly into the shorter-tree part and the carriers-proliferate part.
Panel C  verdict — the fraction of replicates on which a default BiSSE fit rejects
         state-independence, for the driver (power) and the control (calibration).

Self-contained: reads the two JSON files beside it, writes ``figures/bisse.{png,pdf}``.

    python figures.py
"""
from __future__ import annotations

import json
import pathlib

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection

HERE = pathlib.Path(__file__).parent
FIG = HERE / "figures"

INK, DARK, MUTED = "#1a1a1a", "#4d4d4d", "#8a8a8a"
HEADLINE, NULL = "3.0", "1.0"

plt.rcParams.update({
    "font.size": 8.5, "axes.titlesize": 10, "axes.labelsize": 9,
    "legend.fontsize": 8, "xtick.labelsize": 8, "ytick.labelsize": 8,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.grid": True, "grid.linewidth": 0.5, "grid.alpha": 0.35,
    "lines.linewidth": 1.4, "savefig.bbox": "tight",
})


def series_style(which: str) -> dict:
    if which == "driver":
        return dict(color=INK, marker="o", markerfacecolor="white", linestyle="-")
    return dict(color=DARK, marker="s", markerfacecolor="white",
                linestyle=(0, (5, 2)))


def log_x(ax, factors) -> None:
    ax.set_xscale("log")
    ax.set_xlim(0.94, 5.5)
    ax.set_xticks(factors)
    ax.set_xticklabels([f"{f:g}" for f in factors])
    ax.xaxis.set_minor_locator(matplotlib.ticker.NullLocator())


def panel_a(ax, d) -> None:
    factors = d["parameters"]["factors"]
    xs = np.array(factors)
    for which in ("control", "driver"):
        s = [d["summary"]["sweep"][repr(f)][f"prevalence_{which}"] for f in factors]
        m = np.array([x["mean"] for x in s])
        lo = np.array([x["ci95_low"] for x in s])
        hi = np.array([x["ci95_high"] for x in s])
        ax.errorbar(xs, m, yerr=[m - lo, hi - m], capsize=2, elinewidth=0.8,
                    label=f"{which} (undriven)" if which == "control" else which,
                    **series_style(which))
    hts = [d["summary"]["sweep"][repr(f)]["height"]["mean"] for f in factors]
    for x, h in zip(xs, hts):
        ax.annotate(f"{h:.1f}", (x, 0.02), ha="center", color=MUTED, fontsize=7)
    ax.annotate("mean tree height", (0.98, 0.075), xycoords="axes fraction",
                ha="right", color=MUTED, fontsize=7)
    ax.set_ylim(0, 1.05)
    ax.set_ylabel("Prevalence among extant tips")
    ax.set_xlabel("Driver factor $f$ ($\\times$ speciation when present)")
    log_x(ax, factors)
    ax.legend(loc="upper left", frameon=False)
    ax.set_title("A   Prevalence against the driver factor", loc="left",
                 fontweight="bold")


def panel_b(ax, d) -> None:
    def rep(f, field):
        return np.array([r[field] for r in d["replicates"]["sweep"][f]])

    cols = (rep(NULL, "prevalence_driver"), rep(HEADLINE, "prevalence_control"),
            rep(HEADLINE, "prevalence_driver"))
    pos = np.array([0.0, 1.0, 2.0])
    means = np.array([c.mean() for c in cols])

    seg = np.stack([np.repeat(pos[None, :], cols[0].size, 0), np.stack(cols, 1)], -1)
    ax.add_collection(LineCollection(seg, colors=MUTED, linewidths=0.4, alpha=0.15))
    for x, v, ec in zip(pos, cols, (INK, MUTED, INK)):
        ax.boxplot([v], positions=[x], widths=0.3, patch_artist=True,
                   medianprops=dict(color=ec, linewidth=1.4),
                   boxprops=dict(facecolor="white", edgecolor=ec, alpha=0.9),
                   whiskerprops=dict(color=ec), capprops=dict(color=ec),
                   flierprops=dict(marker="none"))
    ax.plot(pos, means, "D", color=INK, markersize=4, zorder=6)

    xa = 2.55
    for lo, hi, lab in ((means[0], means[1], "shorter tree,\nless loss"),
                        (means[1], means[2], "carriers\nproliferate")):
        ax.annotate("", xy=(xa, hi), xytext=(xa, lo),
                    arrowprops=dict(arrowstyle="<->", color=DARK, linewidth=0.9))
        ax.text(xa + 0.1, (lo + hi) / 2, f"$+${hi - lo:.2f}\n{lab}",
                va="center", fontsize=7.5, color=INK)
    head = d["summary"]["headline"]["paired_coupled_minus_null"]
    ax.set_title("B   The headline factor against its matched null", loc="left",
                 fontweight="bold")
    ax.text(0.98, 0.03, f"paired gap $+${head['mean']:.2f} "
            f"(95% CI {head['ci95_low']:.2f}–{head['ci95_high']:.2f})",
            transform=ax.transAxes, ha="right", va="bottom", fontsize=7.5,
            color=DARK)
    ax.set_xticks(pos)
    ax.set_xticklabels(["driver\n$f=1$", "control\n$f=3$", "driver\n$f=3$"])
    ax.set_xlim(-0.5, 3.9)
    ax.set_ylim(0, 1.05)
    ax.set_ylabel("Prevalence among extant tips")


def panel_c(ax, d, v) -> None:
    factors = d["parameters"]["factors"]
    xs = np.array(factors)
    alpha = v["alpha"]
    ax.axhline(alpha, color=MUTED, linewidth=0.8, linestyle=(0, (2, 2)))
    ax.annotate(f"$\\alpha = {alpha:g}$", (5.2, alpha + 0.01), ha="right",
                color=MUTED, fontsize=7.5)
    for which in ("control", "driver"):
        cells = [v["per_factor"][repr(f)][which] for f in factors]
        m = np.array([c["positive_rate_among_fit"] for c in cells])
        lo = np.array([c["wilson_ci95"][0] for c in cells])
        hi = np.array([c["wilson_ci95"][1] for c in cells])
        ax.errorbar(xs, m, yerr=[m - lo, hi - m], capsize=2, elinewidth=0.8,
                    label=f"{which} (undriven)" if which == "control" else which,
                    **series_style(which))
    ax.set_ylim(0, 0.55)
    ax.set_ylabel("Fraction of fits with $p<0.05$")
    ax.set_xlabel("Driver factor $f$")
    log_x(ax, factors)
    ax.legend(loc="upper left", frameon=False)
    ax.set_title("C   What a default BiSSE fit reports", loc="left",
                 fontweight="bold")


def main() -> int:
    d = json.loads((HERE / "results.json").read_text())
    v = json.loads((HERE / "verdicts.json").read_text())
    FIG.mkdir(exist_ok=True)
    fig, axes = plt.subplots(1, 3, figsize=(12.5, 3.6))
    panel_a(axes[0], d)
    panel_b(axes[1], d)
    panel_c(axes[2], d, v)
    fig.tight_layout(w_pad=2.5)
    for ext in ("png", "pdf"):
        out = FIG / f"bisse.{ext}"
        fig.savefig(out, dpi=300)
        print(" ", out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
