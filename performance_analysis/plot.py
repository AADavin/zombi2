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
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter

from perfkit import FIGURES_DIR, load_all, one_line, style

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
_GENE_ORDER = ["Rust · full genomes", "Rust · profiles only"]
_GENE_LABELS = {"Rust · full genomes": "full gene trees",
                "Rust · profiles only": "profiles only"}
_MEM_ORDER = ["Species tree", "Gene families · full", "Gene families · profiles"]


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


def fig_overview(results):
    fig, (a, b) = plt.subplots(1, 2, figsize=(12.4, 4.8))
    if "species_tree" in results:
        draw_scaling(a, results["species_tree"], series_order=_TREE_ORDER)
        a.set_title("a  Species-tree simulation", loc="left", pad=8)
    if "gene_families" in results:
        draw_scaling(b, results["gene_families"], series_order=_GENE_ORDER,
                     label_map=_GENE_LABELS)
        b.set_title("b  Gene-family simulation", loc="left", pad=8)

    fig.suptitle("ZOMBI2 performance", x=0.01, ha="left", fontsize=16, weight="bold")
    fig.tight_layout(rect=(0, 0.015, 1, 0.94))
    _caption(fig, list(results.values()))
    return style.save(fig, FIGURES_DIR / "overview")


# --- driver ---------------------------------------------------------------

FIGURES = [
    ("species_tree", fig_species_tree),
    ("gene_families", fig_gene_families),
    ("memory_scaling", fig_memory),
    ("parallel_scaling", fig_parallel),
    ("write_output", fig_write),
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
            svg, _ = fn(results)
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
