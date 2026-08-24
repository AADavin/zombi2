"""The four-panel figure for the rearrangement recovery study, from ``results.json``.

Panel A  the ridge — the ABC distance over the inversions-only (rate x extent) grid. The
         band of equally good fits runs the full height of the grid, so the extent is free.
Panel B  the identifiability statement — distance against each parameter as a multiple of its true
         value, profiled over the other. The rate has a minimum at the truth; the extent is flat.
Panel C  the mixed experiment — distance over the (inversion rate x translocation rate) grid. Both
         rates have a joint minimum, so a mixed model's two rates can be told apart.
Panel D  recovery — the recovered value as a multiple of the truth, for every experiment, axis and
         observed replicate. The last experiment repeats the mixed one with the extent fixed at
         twice the truth.

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
from matplotlib.lines import Line2D

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
    """The distance as a (y, x) array, with the grid's axis values."""
    xs = np.array(arm["grid"][x_axis], float)
    ys = np.array(arm["grid"][y_axis], float)
    d = np.array(arm["inference"][observed]["distances"], float)
    return xs, ys, d.reshape(len(ys), len(xs))


def outline_cells(ax, mask, color="white", lw=1.2):
    """Outline the cells of a boolean (y, x) mask: an edge wherever inside meets outside."""
    ny, nx = mask.shape
    for j in range(ny):
        for i in range(nx):
            if not mask[j, i]:
                continue
            if i == 0 or not mask[j, i - 1]:
                ax.plot([i - .5, i - .5], [j - .5, j + .5], color=color, lw=lw, zorder=4)
            if i == nx - 1 or not mask[j, i + 1]:
                ax.plot([i + .5, i + .5], [j - .5, j + .5], color=color, lw=lw, zorder=4)
            if j == 0 or not mask[j - 1, i]:
                ax.plot([i - .5, i + .5], [j - .5, j - .5], color=color, lw=lw, zorder=4)
            if j == ny - 1 or not mask[j + 1, i]:
                ax.plot([i - .5, i + .5], [j + .5, j + .5], color=color, lw=lw, zorder=4)


def draw_surface(ax, xs, ys, grid, truth_x, truth_y, best_x, best_y, *,
                 xlabel, ylabel, xtick_step=2, ytick_labels=None):
    # equal-width cells on a log axis: plot against the level index, label with the values
    mesh = ax.pcolormesh(np.arange(len(xs) + 1) - 0.5, np.arange(len(ys) + 1) - 0.5,
                         grid, cmap=CMAP, shading="flat")
    # the best 5% of cells, outlined cell by cell (no interpolated contour)
    outline_cells(ax, grid <= np.nanpercentile(grid, 5))
    ax.plot(np.argmin(np.abs(xs - truth_x)), np.argmin(np.abs(ys - truth_y)), "*",
            color="white", markersize=15, markeredgecolor=INK, markeredgewidth=0.6,
            zorder=5)
    ax.plot(np.argmin(np.abs(xs - best_x)), np.argmin(np.abs(ys - best_y)), "o",
            color="none", markersize=10, markeredgecolor="white", markeredgewidth=1.6,
            zorder=5)
    # rates are shown in units of 10^-3, at every other grid level, so the ticks are short
    ax.set_xticks(np.arange(len(xs))[::xtick_step])
    ax.set_xticklabels([f"{v * 1e3:g}" for v in xs[::xtick_step]])
    ax.set_yticks(np.arange(len(ys)))
    ax.set_yticklabels(ytick_labels if ytick_labels is not None
                       else [f"{v:g}" for v in ys])
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    handles = [
        Line2D([], [], marker="*", linestyle="none", color="white",
               markeredgecolor=INK, markeredgewidth=0.6, markersize=11, label="truth"),
        Line2D([], [], marker="o", linestyle="none", color="none",
               markeredgecolor="white", markeredgewidth=1.6, markersize=7,
               label="best fit"),
        Line2D([], [], marker="s", linestyle="none", color="none",
               markeredgecolor="white", markeredgewidth=1.2, markersize=7,
               label="best 5% of cells"),
    ]
    ax.legend(handles=handles, loc="upper left", frameon=True, facecolor="#8a8f94",
              framealpha=0.95, edgecolor="none", fontsize=7.5, labelspacing=0.7,
              labelcolor="white")
    return mesh


def panel_a(ax, fig, d):
    arm = d["arms"]["A"]
    xs, ys, grid = surface(arm, "inversion", "extent")
    best = arm["inference"][0]["best"]
    mesh = draw_surface(
        ax, xs, ys, grid, arm["truth"]["inversion"], arm["truth"]["extent"],
        best["inversion"], best["extent"],
        xlabel=r"Inversion rate ($10^{-3}$ per gene per time unit)",
        ylabel="Inversion extent (genes)")
    fig.colorbar(mesh, ax=ax, label="ABC distance", pad=0.02)


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
    ax.set_ylabel("ABC distance (best over the other axis)")
    ax.legend(loc="upper center", frameon=False)


def panel_c(ax, fig, d):
    arm = d["arms"]["B"]
    xs, ys, grid = surface(arm, "inversion", "translocation")
    best = arm["inference"][0]["best"]
    mesh = draw_surface(
        ax, xs, ys, grid, arm["truth"]["inversion"], arm["truth"]["translocation"],
        best["inversion"], best["translocation"],
        xlabel=r"Inversion rate ($10^{-3}$ per gene per time unit)",
        ylabel=r"Translocation rate ($10^{-3}$ per gene per time unit)",
        ytick_labels=[f"{v * 1e3:g}" for v in ys])
    fig.colorbar(mesh, ax=ax, label="ABC distance", pad=0.02)


def panel_d(ax, fig, d):
    # the twin of panel C: the same grid, refit with the extent fixed at twice the truth.
    # In C the best fit sits on the star; here it slides one cell off it.
    arm = d["arms"]["C"]
    xs, ys, grid = surface(arm, "inversion", "translocation")
    best = arm["inference"][0]["best"]
    mesh = draw_surface(
        ax, xs, ys, grid, arm["truth"]["inversion"], arm["truth"]["translocation"],
        best["inversion"], best["translocation"],
        xlabel=r"Inversion rate ($10^{-3}$ per gene per time unit)",
        ylabel=r"Translocation rate ($10^{-3}$ per gene per time unit)",
        ytick_labels=[f"{v * 1e3:g}" for v in ys])
    fig.colorbar(mesh, ax=ax, label="ABC distance", pad=0.02)


def main() -> int:
    d = json.loads((HERE / "results.json").read_text())
    FIG.mkdir(exist_ok=True)
    fig, axes = plt.subplots(2, 2, figsize=(12.0, 8.4))
    panel_a(axes[0][0], fig, d)
    panel_b(axes[0][1], d)
    panel_c(axes[1][0], fig, d)
    panel_d(axes[1][1], fig, d)
    for ax, letter in zip(axes.flat, "ABCD"):
        ax.text(-0.08, 1.04, letter, transform=ax.transAxes,
                fontsize=13, fontweight="bold", va="bottom")
    fig.tight_layout(w_pad=2.5, h_pad=3.0)
    for ext in ("png", "pdf"):
        out = FIG / f"rearrangements.{ext}"
        fig.savefig(out, dpi=300)
        print(" ", out)
    plt.close(fig)
    # one image per panel too, large enough to read inline on the website
    for letter, fn, needs_fig in (("a", panel_a, True), ("b", panel_b, False),
                                  ("c", panel_c, True), ("d", panel_d, True)):
        fig, ax = plt.subplots(figsize=(6.6, 4.6))
        fn(ax, fig, d) if needs_fig else fn(ax, d)
        fig.tight_layout()
        out = FIG / f"rearrangements_{letter}.png"
        fig.savefig(out, dpi=300)
        print(" ", out)
        plt.close(fig)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
