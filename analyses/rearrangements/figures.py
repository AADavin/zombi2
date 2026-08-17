"""The four-panel figure for the rearrangement recovery study, from ``results.json``.

Panel A  the ridge — the ABC misfit over arm A's (inversion rate x inversion extent) grid. The
         band of equally good fits runs the full height of the grid, so the extent is free.
Panel B  the identifiability statement — misfit against each parameter as a multiple of its true
         value, profiled over the other. The rate has a minimum at the truth; the extent is flat.
Panel C  the mixed arm — misfit over arm B's (inversion rate x translocation rate) grid. Both
         rates have a joint minimum, so a mixed model's two rates can be told apart.
Panel D  recovery — the recovered value as a multiple of the truth, for every arm, axis and
         observed replicate. Arm C repeats arm B with the extent fixed at twice the truth.

Self-contained: reads ``results.json`` beside it, writes ``figures/rearrangements.{png,pdf}``.

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
CMAP = "viridis"

plt.rcParams.update({
    "font.size": 8.5, "axes.titlesize": 10, "axes.labelsize": 9,
    "legend.fontsize": 8, "xtick.labelsize": 8, "ytick.labelsize": 8,
    "axes.spines.top": False, "axes.spines.right": False,
    "lines.linewidth": 1.4, "savefig.bbox": "tight",
})


def surface(arm, x_axis, y_axis, observed=0):
    """The misfit as a (y, x) array, with the grid's axis values."""
    xs = np.array(arm["grid"][x_axis], float)
    ys = np.array(arm["grid"][y_axis], float)
    d = np.array(arm["inference"][observed]["distances"], float)
    return xs, ys, d.reshape(len(ys), len(xs))


def draw_surface(ax, xs, ys, grid, truth_x, truth_y, best_x, best_y, *, xlabel, ylabel, title):
    # equal-width cells on a log axis: plot against the level index, label with the values
    mesh = ax.pcolormesh(np.arange(len(xs) + 1) - 0.5, np.arange(len(ys) + 1) - 0.5,
                         grid, cmap=CMAP, shading="flat")
    levels = [np.nanpercentile(grid, 5)]
    ax.contour(np.arange(len(xs)), np.arange(len(ys)), grid, levels=levels,
               colors="white", linewidths=1.0)
    ax.plot(np.argmin(np.abs(xs - truth_x)), np.argmin(np.abs(ys - truth_y)), "*",
            color="white", markersize=15, markeredgecolor=INK, markeredgewidth=0.6,
            label="truth", zorder=5)
    ax.plot(np.argmin(np.abs(xs - best_x)), np.argmin(np.abs(ys - best_y)), "o",
            color="none", markersize=10, markeredgecolor="white", markeredgewidth=1.6,
            label="best fit", zorder=5)
    step = max(1, len(xs) // 7)
    ax.set_xticks(np.arange(len(xs))[::step])
    ax.set_xticklabels([f"{v:.3g}" for v in xs[::step]], rotation=45, ha="right")
    ax.set_yticks(np.arange(len(ys)))
    ax.set_yticklabels([f"{v:g}" for v in ys])
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title, loc="left", fontweight="bold")
    ax.legend(loc="upper left", frameon=True, facecolor="white", framealpha=0.75,
              edgecolor="none", fontsize=7.5)
    return mesh


def panel_a(ax, fig, d):
    arm = d["arms"]["A"]
    xs, ys, grid = surface(arm, "inversion", "extent")
    best = arm["inference"][0]["best"]
    mesh = draw_surface(
        ax, xs, ys, grid, arm["truth"]["inversion"], arm["truth"]["extent"],
        best["inversion"], best["extent"],
        xlabel="Inversion rate (per gene per time unit)",
        ylabel="Inversion extent (genes)",
        title="A   Inversions only: the rate-extent ridge")
    fig.colorbar(mesh, ax=ax, label="ABC misfit", pad=0.02)


def panel_b(ax, d):
    arm = d["arms"]["A"]
    styles = {"inversion": dict(color=INK, marker="o", linestyle="-"),
              "extent": dict(color=DARK, marker="s", linestyle=(0, (5, 2)))}
    names = {"inversion": "inversion rate", "extent": "inversion extent"}
    for axis in ("inversion", "extent"):
        marginal = arm["inference"][0]["marginals"][axis]
        levels = np.array(marginal["levels"], float)
        profile = np.array(marginal["best_profile"], float)
        ax.plot(levels / arm["truth"][axis], profile, markerfacecolor="white",
                markersize=4, label=names[axis], **styles[axis])
    ax.axvline(1.0, color=MUTED, linewidth=0.8, linestyle=(0, (2, 2)))
    ax.annotate("truth", (1.0, 0.45), xycoords=("data", "axes fraction"),
                xytext=(4, 0), textcoords="offset points", color=MUTED, fontsize=7.5)
    ax.set_xscale("log", base=2)
    ax.set_xlabel("Value as a multiple of the truth")
    ax.set_ylabel("ABC misfit (best over the other axis)")
    ax.set_title("B   The rate has a minimum, the extent does not", loc="left",
                 fontweight="bold")
    ax.legend(loc="upper center", frameon=False)


def panel_c(ax, fig, d):
    arm = d["arms"]["B"]
    xs, ys, grid = surface(arm, "inversion", "translocation")
    best = arm["inference"][0]["best"]
    mesh = draw_surface(
        ax, xs, ys, grid, arm["truth"]["inversion"], arm["truth"]["translocation"],
        best["inversion"], best["translocation"],
        xlabel="Inversion rate (per gene per time unit)",
        ylabel="Translocation rate (per gene per time unit)",
        title="C   Mixed model: both rates are pinned")
    fig.colorbar(mesh, ax=ax, label="ABC misfit", pad=0.02)


def panel_d(ax, d):
    rows = [("A", "inversion", "A: inversion rate"), ("A", "extent", "A: inversion extent"),
            ("B", "inversion", "B: inversion rate"),
            ("B", "translocation", "B: translocation rate"),
            ("C", "inversion", "C: inversion rate\n(extent set wrong)"),
            ("C", "translocation", "C: translocation rate\n(extent set wrong)")]
    ax.axvline(1.0, color=MUTED, linewidth=0.9, linestyle=(0, (2, 2)))
    labels, spans = [], []
    for i, (arm_name, axis, label) in enumerate(rows):
        entry = d["verdicts"][arm_name][axis]
        ratios = np.array(entry["best_per_observed"], float) / entry["true"]
        y = len(rows) - 1 - i
        colour = INK if axis != "extent" else DARK
        # the three observed replicates usually agree exactly, so offset them to stay countable
        offsets = (np.arange(len(ratios)) - (len(ratios) - 1) / 2) * 0.13
        ax.plot(ratios, y + offsets, "o", color=colour, markerfacecolor="white",
                markersize=5.5, markeredgewidth=1.2, linestyle="none")
        labels.append(label)
        spans.append((y, float(np.mean(entry["credible_span_fraction"]))))
    ax.set_xscale("log", base=2)
    ax.set_ylim(-0.6, len(rows) - 0.2)
    ax.margins(x=0.25)
    # annotate only once the axis limits are settled by all the data
    right = ax.get_xlim()[1]
    for y, span in spans:
        ax.annotate(f"credible interval covers {span:.0%} of the grid", (right, y),
                    xytext=(-2, 8), textcoords="offset points", ha="right",
                    color=MUTED, fontsize=7)
    ax.set_yticks(np.arange(len(rows)))
    ax.set_yticklabels(labels[::-1])
    ax.set_xlabel("Recovered value as a multiple of the truth")
    ax.set_title("D   What each arm recovered", loc="left", fontweight="bold")
    ax.spines["left"].set_visible(False)
    ax.tick_params(axis="y", length=0)


def main() -> int:
    d = json.loads((HERE / "results.json").read_text())
    FIG.mkdir(exist_ok=True)
    fig, axes = plt.subplots(2, 2, figsize=(12.0, 8.4))
    panel_a(axes[0][0], fig, d)
    panel_b(axes[0][1], d)
    panel_c(axes[1][0], fig, d)
    panel_d(axes[1][1], d)
    fig.tight_layout(w_pad=2.5, h_pad=3.0)
    for ext in ("png", "pdf"):
        out = FIG / f"rearrangements.{ext}"
        fig.savefig(out, dpi=300)
        print(" ", out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
