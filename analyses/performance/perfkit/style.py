"""Publication house style — monochrome (black & white).

A print-first, journal-grade look: every series is drawn in black or grey and
told apart by **line style** (solid / dashed / dotted / dash-dot) and **open
marker shape** (white fill, dark edge), never by colour. This reproduces
perfectly in greyscale, photocopies, and single-ink print.

Consistent with the tree figures' typographic choices (Helvetica, near-black
ink, spare frame); only the colour axis is removed.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt

# --- monochrome ink ------------------------------------------------------
INK = "#1a1a1a"          # primary lines, text, axes
DARK = "#4d4d4d"         # secondary series
MUTED = "#8a8a8a"        # guide lines, captions
GRID = "#e6e6e6"         # hairline grid
PANEL = "#ffffff"

# Series → (ink, linestyle, marker). Distinguished with no colour at all:
# solid+filled-shape for the "primary" of each pair, dashed+different-shape for
# the "secondary". Markers render open (white fill) via rcParams below.
LS_SOLID = "-"
LS_DASH = (0, (5, 2))
LS_DOT = (0, (1, 2))
LS_DASHDOT = (0, (5, 2, 1, 2))

SERIES_STYLE = {
    "Species tree · backward": dict(color=INK,  linestyle=LS_SOLID,   marker="o"),
    "Species tree · forward":  dict(color=DARK, linestyle=LS_DASH,    marker="D"),
    "Rust · full genomes":     dict(color=INK,  linestyle=LS_SOLID,   marker="s"),
    "Rust · event trace":      dict(color=DARK, linestyle=LS_DASHDOT, marker="D"),
    "Rust · profiles only":    dict(color=MUTED, linestyle=LS_DOT,    marker="^"),
    # memory footprint series
    "Species tree":            dict(color=INK,  linestyle=LS_SOLID,   marker="o"),
    "Gene families · full":    dict(color=DARK, linestyle=LS_DASH,    marker="s"),
    "Gene families · trace":   dict(color=DARK, linestyle=LS_DASHDOT, marker="D"),
    "Gene families · profiles":dict(color=MUTED, linestyle=LS_DOT,    marker="^"),
    # cross-tool comparison (ZOMBI2 vs ZOMBI 1)
    "ZOMBI2 · Rust":           dict(color=INK,  linestyle=LS_SOLID,   marker="o"),
    "ZOMBI 1 · Python":        dict(color=DARK, linestyle=LS_DASH,    marker="s"),
    # writing / parallel
    "write()":                 dict(color=INK,  linestyle=LS_SOLID,   marker="v"),
    "measured":                dict(color=INK,  linestyle=LS_SOLID,   marker="o"),
    "ideal (linear)":          dict(color=MUTED, linestyle=LS_DASH,   marker=None),
}

# Fallback for unregistered series: cycle line styles + marker shapes, all black.
_FALLBACK_LS = [LS_SOLID, LS_DASH, LS_DOT, LS_DASHDOT]
_FALLBACK_MK = ["o", "s", "^", "D", "v", "P"]


def style_for(series: str, index: int = 0) -> dict:
    if series in SERIES_STYLE:
        return dict(SERIES_STYLE[series])
    return dict(color=INK, linestyle=_FALLBACK_LS[index % len(_FALLBACK_LS)],
                marker=_FALLBACK_MK[index % len(_FALLBACK_MK)])


def apply() -> None:
    """Install the monochrome house rcParams. Idempotent."""
    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
        "font.size": 11,
        "axes.titlesize": 13,
        "axes.labelsize": 12,
        "legend.fontsize": 10,
        "xtick.labelsize": 10,
        "ytick.labelsize": 10,
        "text.color": INK,
        "axes.edgecolor": INK,
        "axes.labelcolor": INK,
        "xtick.color": INK,
        "ytick.color": INK,
        "axes.linewidth": 1.0,
        "lines.linewidth": 1.8,
        "lines.markersize": 6.0,
        # open markers: white fill, dark edge — the monochrome distinguisher
        "lines.markerfacecolor": "white",
        "lines.markeredgewidth": 1.3,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": True,
        "grid.color": GRID,
        "grid.linewidth": 0.8,
        "axes.axisbelow": True,
        "figure.facecolor": PANEL,
        "axes.facecolor": PANEL,
        "savefig.facecolor": PANEL,
        "figure.dpi": 120,
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
    })


def plot_series(ax, xs, ys, series, index=0, **overrides):
    """Plot one monochrome series with its registered line style + open marker."""
    st = style_for(series, index)
    st.update(overrides)
    marker = st.pop("marker")
    return ax.plot(xs, ys, marker=marker, markeredgecolor=st["color"],
                   markerfacecolor="white", **st)


def caption(fig, text: str) -> None:
    fig.text(0.005, 0.004, text, fontsize=7.0, color=MUTED, ha="left", va="bottom")


def save(fig, stem: Path | str) -> tuple[Path, Path]:
    """Write ``<stem>.svg`` and ``<stem>.png`` (300 dpi). Returns both paths."""
    stem = Path(stem)
    stem.parent.mkdir(parents=True, exist_ok=True)
    svg, png = stem.with_suffix(".svg"), stem.with_suffix(".png")
    fig.savefig(svg)
    fig.savefig(png)
    plt.close(fig)
    return svg, png
