#!/usr/bin/env python
"""Render publication figures from ``results/*.json`` (monochrome house style).

    python plot.py            # render every figure it has data for

Reads only the JSON written by ``run.py`` — never imports the simulator — so you
can restyle endlessly without re-measuring. Each benchmark gets a standalone
figure plus a combined 2x2 overview.
"""

from __future__ import annotations

import argparse

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter

from perfkit import FIGURES_DIR, load_all, one_line, style
from perfkit.style import INK, DARK, MUTED

style.apply()


# --- axis formatting ------------------------------------------------------

def _fmt_count(v, _pos=None) -> str:
    v = float(v)
    for div, suf in ((1e6, "M"), (1e3, "k")):
        if v >= div:
            q = v / div
            return f"{q:g}{suf}"
    return f"{int(v):g}" if v >= 1 else f"{v:g}"


def _fmt_time(v, _pos=None) -> str:
    v = float(v)
    if v <= 0:
        return "0"
    if v < 1e-3:
        return f"{v*1e6:g} µs"
    if v < 1:
        return f"{v*1e3:g} ms"
    return f"{v:g} s"


def _fmt_mem(v, _pos=None) -> str:
    v = float(v)
    if v >= 1000:
        return f"{v/1000:g} GB"
    return f"{v:g} MB"


def _log_axes(ax, x_label, y_label, y_fmt=_fmt_time):
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    ax.xaxis.set_major_formatter(FuncFormatter(_fmt_count))
    ax.yaxis.set_major_formatter(FuncFormatter(y_fmt))
    ax.tick_params(which="minor", length=0)


# --- reusable drawers -----------------------------------------------------

def draw_scaling(ax, result, *, series_order=None, label_map=None):
    """Median line + interquartile band + open markers, log-log (time)."""
    labels = series_order or result.series()
    label_map = label_map or {}
    for i, label in enumerate(labels):
        pts = result.by_series(label)
        if not pts:
            continue
        xs = np.array([p.x for p in pts], float)
        med = np.array([p.median for p in pts])
        lo = np.array([np.percentile(p.times, 25) for p in pts])
        hi = np.array([np.percentile(p.times, 75) for p in pts])
        st = style.style_for(label, i)
        ax.fill_between(xs, lo, hi, color=st["color"], alpha=0.12, linewidth=0)
        style.plot_series(ax, xs, med, label, i, label=label_map.get(label, label),
                          zorder=3)

    _log_axes(ax, result.x_label, "Wall-clock time")
    ax.legend(frameon=False, loc="upper left", handlelength=3.2)


def draw_single(ax, result, series, *, band=True):
    """One monochrome scaling line (median + IQR band), no legend — the model-agnostic
    'it scales linearly' panel."""
    pts = result.by_series(series)
    if not pts:
        return
    xs = np.array([p.x for p in pts], float)
    med = np.array([p.median for p in pts])
    if band:
        lo = np.array([np.percentile(p.times, 25) for p in pts])
        hi = np.array([np.percentile(p.times, 75) for p in pts])
        ax.fill_between(xs, lo, hi, color=INK, alpha=0.12, linewidth=0)
    style.plot_series(ax, xs, med, series, 0, color=INK, linestyle=style.LS_SOLID,
                      zorder=3)
    _log_axes(ax, result.x_label, "Wall-clock time")


def draw_memory(ax, result, *, series_order=None):
    """Peak RSS vs tip count, log-log."""
    labels = series_order or result.series()
    for i, label in enumerate(labels):
        pts = [p for p in result.by_series(label) if "rss_mb" in p.work]
        if not pts:
            continue
        xs = np.array([p.x for p in pts], float)
        ys = np.array([p.work["rss_mb"] for p in pts])
        style.plot_series(ax, xs, ys, label, i, label=label, zorder=3)
    _log_axes(ax, result.x_label, "Peak memory (RSS)", y_fmt=_fmt_mem)
    ax.legend(frameon=False, loc="upper left", handlelength=3.2)


def draw_parallel(ax, result):
    """Parallel speed-up t(1)/t(p) with an ideal-linear reference."""
    pts = result.by_series(result.series()[0])
    if not pts:
        return
    xs = np.array([p.x for p in pts], float)
    t = np.array([p.best for p in pts])
    base = t[np.argmin(xs)]
    speed = base / t
    lim = [1, xs.max()]
    ax.plot(lim, lim, color=style.MUTED, linestyle=(0, (4, 3)), linewidth=1.2,
            label="ideal (linear)")
    style.plot_series(ax, xs, speed, "measured", label="measured", zorder=3)
    ax.set_xlabel("Worker processes")
    ax.set_ylabel("Speed-up  t(1) / t(p)  (×)")
    peak = speed.max()
    ax.legend(frameon=False, loc="upper left", handlelength=2.4)
    return peak


# --- standalone figures ---------------------------------------------------

def _new_single():
    fig, ax = plt.subplots(figsize=(6.4, 4.8))
    fig.subplots_adjust(left=0.14, right=0.965, top=0.9, bottom=0.17)
    return fig, ax


def _caption(fig, results):
    envs = [r.env for r in results if r.env]
    if envs:
        style.caption(fig, one_line(envs[0]))


_TREE_ORDER = ["Species tree · backward", "Species tree · forward"]
_GENE_ORDER = ["Rust · full genomes", "Rust · event trace", "Rust · profiles only"]
_GENE_LABELS = {"Rust · full genomes": "full event log",
                "Rust · event trace": "event trace",
                "Rust · profiles only": "profiles (counts)"}

# Overview panels (the README figure): one species-tree line, the two extreme
# genome output modes, and the fixed-tree head-to-head box-plot.
_OVERVIEW_TREE = "Species tree · backward"        # a single generic scaling line
_OVERVIEW_GENE_ORDER = ["Rust · full genomes", "Rust · profiles only"]
_OVERVIEW_GENE_LABELS = {"Rust · full genomes": "Full gene trees",
                         "Rust · profiles only": "Only profiles"}
_MEM_ORDER = ["Species tree", "Gene families · full", "Gene families · trace",
              "Gene families · profiles"]
_VS_ORDER = ["ZOMBI2 · Rust", "ZOMBI 1 · Python"]


def draw_vs_zombi1(ax, result):
    """Two log-log curves (ZOMBI2 vs ZOMBI 1) with the speed-up annotated at the largest
    tip count both tools reached, and ZOMBI 1's practical ceiling marked."""
    draw_scaling(ax, result, series_order=_VS_ORDER)
    z2 = {p.x: p.median for p in result.by_series("ZOMBI2 · Rust")}
    z1 = {p.x: p.median for p in result.by_series("ZOMBI 1 · Python")}
    common = sorted(set(z2) & set(z1))
    if common:
        x = common[-1]
        factor = z1[x] / z2[x]
        ax.annotate(f"≈ {factor:,.0f}× faster\nat {_fmt_count(x)} tips",
                    xy=(x, (z1[x] * z2[x]) ** 0.5), xytext=(0.30, 0.72),
                    textcoords="axes fraction", fontsize=9.5, color=style.INK,
                    ha="center", va="center",
                    arrowprops=dict(arrowstyle="->", color=style.MUTED, linewidth=1.0))
    if z1:  # ZOMBI 1's largest finishing size — its practical ceiling
        ceil = max(z1)
        ax.axvline(ceil, color=style.MUTED, linestyle=(0, (1, 2)), linewidth=1.0, zorder=1)
        ax.text(ceil, ax.get_ylim()[0], "ZOMBI 1\nceiling ", rotation=0,
                fontsize=8, color=style.MUTED, ha="right", va="bottom")


_VS_BOX_ORDER = ["ZOMBI 1 · Python", "ZOMBI2 · Rust"]     # slow → fast, left → right
_VS_BOX_LABELS = {"ZOMBI 1 · Python": "ZOMBI 1\n(Python)",
                  "ZOMBI2 · Rust": "ZOMBI2\n(Rust)"}


def draw_boxplot_vs(ax, result):
    """Box-plot of per-run wall-clock times for both engines on ONE shared species tree.

    Each engine ran once per seed on the identical tree, so the boxes are honest run-to-run
    distributions at a single size. Individual runs are overlaid as open markers, and the
    median speed-up is annotated with a double-headed arrow spanning the two medians."""
    data, positions, present = [], [], []
    for i, s in enumerate(_VS_BOX_ORDER):
        pts = result.by_series(s)
        if not pts or not pts[0].times:
            continue
        data.append(np.array(pts[0].times, float))
        positions.append(i)
        present.append(s)
    if not data:
        return

    ax.set_yscale("log")
    ax.boxplot(data, positions=positions, widths=0.52, patch_artist=True,
               medianprops=dict(color=INK, linewidth=2.4),
               boxprops=dict(facecolor="white", edgecolor=INK, linewidth=1.5),
               whiskerprops=dict(color=INK, linewidth=1.3),
               capprops=dict(color=INK, linewidth=1.3),
               flierprops=dict(marker="none"), zorder=2)

    # overlay every individual run as an open marker (deterministic jitter)
    for pos, s, arr in zip(positions, present, data):
        st = style.style_for(s)
        jit = np.linspace(-0.14, 0.14, len(arr)) if len(arr) > 1 else np.array([0.0])
        ax.plot(pos + jit, arr, linestyle="none", marker=st["marker"], markersize=7,
                markerfacecolor="white", markeredgecolor=st["color"],
                markeredgewidth=1.4, alpha=0.95, zorder=4)

    ax.set_xticks(positions)
    ax.set_xticklabels([_VS_BOX_LABELS[s] for s in present])
    ax.set_xlim(-0.62, len(present) - 0.38)
    ax.set_ylabel("Wall-clock time")
    ax.yaxis.set_major_formatter(FuncFormatter(_fmt_time))
    ax.tick_params(which="minor", length=0)
    ax.margins(y=0.22)
    # controlled-comparison context (self-updating from the data)
    size = int(result.by_series(present[0])[0].x)
    nruns = max(len(a) for a in data)
    ax.set_xlabel(f"one shared {size:,}-tip tree · {nruns} runs each", fontsize=15)

    # speed-up annotation: a double-headed arrow between the two medians
    med = {s: float(np.median(a)) for s, a in zip(present, data)}
    if "ZOMBI 1 · Python" in med and "ZOMBI2 · Rust" in med:
        factor = med["ZOMBI 1 · Python"] / med["ZOMBI2 · Rust"]
        xa = 0.5  # between the two boxes
        ax.annotate("", xy=(xa, med["ZOMBI 1 · Python"]), xytext=(xa, med["ZOMBI2 · Rust"]),
                    arrowprops=dict(arrowstyle="<->", color=DARK, linewidth=1.6))
        ax.text(xa + 0.06, (med["ZOMBI 1 · Python"] * med["ZOMBI2 · Rust"]) ** 0.5,
                f"≈ {factor:,.0f}×\nfaster", ha="left", va="center", color=INK,
                fontweight="bold")


def fig_species_tree(results):
    r = results["species_tree"]
    fig, ax = _new_single()
    draw_scaling(ax, r, series_order=_TREE_ORDER)
    ax.set_title(r.title, loc="left", pad=10)
    _caption(fig, [r])
    return style.save(fig, FIGURES_DIR / "species_tree_scaling")


def fig_gene_families(results):
    r = results["gene_families"]
    fig, ax = _new_single()
    draw_scaling(ax, r, series_order=_GENE_ORDER, label_map=_GENE_LABELS)
    ax.set_title(r.title, loc="left", pad=10)
    _caption(fig, [r])
    return style.save(fig, FIGURES_DIR / "gene_family_scaling")


def fig_memory(results):
    r = results["memory_scaling"]
    fig, ax = _new_single()
    draw_memory(ax, r, series_order=_MEM_ORDER)
    ax.set_title(r.title, loc="left", pad=10)
    _caption(fig, [r])
    return style.save(fig, FIGURES_DIR / "memory_scaling")


def fig_parallel(results):
    r = results["parallel_scaling"]
    fig, ax = _new_single()
    draw_parallel(ax, r)
    ax.set_title(r.title, loc="left", pad=10)
    _caption(fig, [r])
    return style.save(fig, FIGURES_DIR / "parallel_scaling")


def fig_write(results):
    r = results["write_output"]
    fig, ax = _new_single()
    draw_scaling(ax, r)
    ax.set_title(r.title, loc="left", pad=10)
    _caption(fig, [r])
    return style.save(fig, FIGURES_DIR / "write_scaling")


def fig_vs_zombi1(results):
    r = results["vs_zombi1"]
    fig, ax = _new_single()
    draw_vs_zombi1(ax, r)
    ax.set_title(r.title, loc="left", pad=10)
    _caption(fig, [r])
    return style.save(fig, FIGURES_DIR / "vs_zombi1")


# Larger type just for the README banner: the 3-panel overview is displayed at
# column width, where the standalone figures' sizes would render tiny. Applied via
# rc_context so the single-panel figures keep their own (already-legible) sizing.
_OVERVIEW_RC = {
    "font.size": 18,
    "axes.titlesize": 22,
    "axes.labelsize": 19,
    "legend.fontsize": 18,
    "xtick.labelsize": 16,
    "ytick.labelsize": 16,
    "lines.linewidth": 2.4,
    "lines.markersize": 9.0,
    "lines.markeredgewidth": 1.7,
}


def fig_overview(results):
    """The README banner: (a) one species-tree scaling line, (b) the two extreme genome
    output modes, (c) the fixed-tree ZOMBI 1 vs ZOMBI2 box-plot. Big type for legibility
    at column width. Falls back to (a)+(b) if the head-to-head data is absent."""
    have_vs = "vs_zombi1_fixedtree" in results and results["vs_zombi1_fixedtree"].points
    with mpl.rc_context(_OVERVIEW_RC):
        ncols = 3 if have_vs else 2
        fig, axes = plt.subplots(1, ncols, figsize=(6.9 * ncols, 6.0))
        a, b = axes[0], axes[1]
        if "species_tree" in results:
            draw_single(a, results["species_tree"], _OVERVIEW_TREE)
            a.set_title("a   Species-tree simulation", loc="left", pad=12, weight="bold")
        if "gene_families" in results:
            draw_scaling(b, results["gene_families"], series_order=_OVERVIEW_GENE_ORDER,
                         label_map=_OVERVIEW_GENE_LABELS)
            b.set_title("b   Genome (unordered)", loc="left", pad=12, weight="bold")
            b.legend(frameon=False, loc="upper left", handlelength=3.0)
        if have_vs:
            draw_boxplot_vs(axes[2], results["vs_zombi1_fixedtree"])
            axes[2].set_title("c   ZOMBI 1 vs ZOMBI2", loc="left", pad=12, weight="bold")

        fig.suptitle("ZOMBI2 performance", x=0.008, ha="left", fontsize=25, weight="bold")
        fig.tight_layout(rect=(0, 0.02, 1, 0.93), w_pad=3.0)
        _caption(fig, list(results.values()))
        return style.save(fig, FIGURES_DIR / "overview")


def fig_vs_fixedtree(results):
    r = results["vs_zombi1_fixedtree"]
    if not r.points:
        return None
    fig, ax = _new_single()
    draw_boxplot_vs(ax, r)
    ax.set_title(r.title, loc="left", pad=10)
    _caption(fig, [r])
    return style.save(fig, FIGURES_DIR / "vs_zombi1_fixedtree")


# --- driver ---------------------------------------------------------------

FIGURES = [
    ("species_tree", fig_species_tree),
    ("gene_families", fig_gene_families),
    ("memory_scaling", fig_memory),
    ("parallel_scaling", fig_parallel),
    ("write_output", fig_write),
    ("vs_zombi1", fig_vs_zombi1),
    ("vs_zombi1_fixedtree", fig_vs_fixedtree),
]


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.parse_args(argv)

    results = load_all()
    if not results:
        print("No results found. Run:  python run.py")
        return 1

    made = []
    for name, fn in FIGURES:
        if name in results:
            out = fn(results)
            if out is None:
                continue
            svg, _ = out
            made.append(svg.name)
            print(f"  {svg.relative_to(FIGURES_DIR.parent)}")
    svg, _ = fig_overview(results)
    made.append(svg.name)
    print(f"  {svg.relative_to(FIGURES_DIR.parent)}")
    print(f"\nWrote {len(made)} figures to "
          f"{FIGURES_DIR.relative_to(FIGURES_DIR.parent.parent)}/")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
