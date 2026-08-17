"""The two-panel figure for the BiSSE analysis, from ``results.json`` + ``verdicts.json``.

Panel A  prevalence — driver and control prevalence among extant tips against the driver
         factor, with 95% bootstrap CIs; the gap between the curves is the within-run
         signal, and the decomposition at f = 3 is quoted in prose from the same numbers.
Panel B  verdict — the fraction of replicates on which a default BiSSE fit rejects
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

HERE = pathlib.Path(__file__).parent
FIG = HERE / "figures"

INK, DARK, MUTED = "#1a1a1a", "#4d4d4d", "#8a8a8a"

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
    ax.set_title("A   Prevalence of the driver and control families", loc="left",
                 fontweight="bold")


def panel_b(ax, d, v) -> None:
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
    ax.set_title("B   What a default BiSSE fit reports", loc="left",
                 fontweight="bold")


def panel_c(ax, vs) -> None:
    """Power against tree size, from verdicts_size.json (which folds in the
    150-tip cells of verdicts.json)."""
    sizes = [150, 500, 1000]
    alpha = 0.05
    ax.axhline(alpha, color=MUTED, linewidth=0.8, linestyle=(0, (2, 2)))
    ax.annotate(f"$\\alpha = {alpha:g}$", (137, alpha + 0.015), ha="left",
                color=MUTED, fontsize=7.5)

    def cell(n, f, char):
        return vs["cells"][f"n{n}_f{f}_{char}"]

    for f, style_kw in (("3.0", dict(color=INK, marker="o", linestyle="-")),
                        ("5.0", dict(color=DARK, marker="D",
                                     linestyle=(0, (3, 1.5))))):
        m = np.array([cell(n, f, "driver")["positive_rate_among_fit"]
                      for n in sizes])
        lo = np.array([cell(n, f, "driver")["wilson_ci95"][0] for n in sizes])
        hi = np.array([cell(n, f, "driver")["wilson_ci95"][1] for n in sizes])
        # at a rate of exactly 1.0 the Wilson centre sits below the point
        # estimate, so clamp the bar lengths at zero
        ax.errorbar(sizes, m, yerr=[np.maximum(m - lo, 0), np.maximum(hi - m, 0)],
                    capsize=2, elinewidth=0.8,
                    markerfacecolor="white", label=f"driver, $f={float(f):g}$",
                    **style_kw)
    # the control pooled over both factors at each size: the calibration line
    k = np.array([sum(cell(n, f, "control")["n_significant"]
                      for f in ("3.0", "5.0")) for n in sizes], float)
    nfit = np.array([sum(cell(n, f, "control")["n_fit"]
                         for f in ("3.0", "5.0")) for n in sizes], float)
    ax.plot(sizes, k / nfit, color=MUTED, marker="s", markerfacecolor="white",
            linestyle=(0, (5, 2)), label="control (pooled)")

    ax.set_xscale("log")
    ax.set_xticks(sizes)
    ax.set_xticklabels([str(s) for s in sizes])
    ax.xaxis.set_minor_locator(matplotlib.ticker.NullLocator())
    ax.set_xlim(130, 1150)
    ax.set_ylim(0, 1.05)
    ax.set_xlabel("Extant tips per tree")
    ax.set_ylabel("Fraction of fits with $p<0.05$")
    ax.legend(loc="center left", frameon=False)
    ax.set_title("C   Power against tree size", loc="left", fontweight="bold")


def main() -> int:
    d = json.loads((HERE / "results.json").read_text())
    v = json.loads((HERE / "verdicts.json").read_text())
    vs = json.loads((HERE / "verdicts_size.json").read_text())
    FIG.mkdir(exist_ok=True)
    fig, axes = plt.subplots(1, 3, figsize=(13.0, 3.8))
    panel_a(axes[0], d)
    panel_b(axes[1], d, v)
    panel_c(axes[2], vs)
    fig.tight_layout(w_pad=2.5)
    for ext in ("png", "pdf"):
        out = FIG / f"bisse.{ext}"
        fig.savefig(out, dpi=300)
        print(" ", out)
    plt.close(fig)
    # one image per panel too, large enough to read inline on the website
    for letter, fn, args in (("a", panel_a, (d,)), ("b", panel_b, (d, v)),
                             ("c", panel_c, (vs,))):
        fig, ax = plt.subplots(figsize=(6.4, 4.4))
        fn(ax, *args)
        fig.tight_layout()
        out = FIG / f"bisse_{letter}.png"
        fig.savefig(out, dpi=300)
        print(" ", out)
        plt.close(fig)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
