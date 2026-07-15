"""Figure: how rates work -- "how many clocks, how fast".

Top row: the OPPORTUNITY, i.e. how many clocks a rate is counted over -- one per lineage, one per
gene copy, or one shared by the whole process. Bottom: the consequence, drawn from real ZOMBI2 runs.
A count that TRACKS the growing quantity (per lineage / per copy) compounds -> exponential; a single
shared clock keeps a constant total rate -> linear. Same speed per clock, different NUMBER of clocks.

Run:  python figures/scripts/fig_rate_clocks.py
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

# use THIS worktree's zombi2 (has SharedBirthDeath), not the editable-installed main checkout
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

import cairosvg
import drawsvg as draw

from zombi2 import BirthDeath, SharedBirthDeath, simulate_species_tree

from zombi_style import FONT, INK, FS_TITLE, FS_LABEL, FS_ANNOT, FS_TICK

OUT = Path(__file__).resolve().parent.parent / "rate_clocks"

W, H = 1420, 940
GREY = "#8a8a8a"
VIOLET = "#5a4e86"


def clock(d, cx, cy, r, col=INK, w=2.2):
    """A small clock glyph: face + two hands (12 and 3) + centre dot."""
    d.append(draw.Circle(cx, cy, r, fill="white", stroke=col, stroke_width=w))
    d.append(draw.Line(cx, cy, cx, cy - r * 0.68, stroke=col, stroke_width=w, stroke_linecap="round"))
    d.append(draw.Line(cx, cy, cx + r * 0.52, cy, stroke=col, stroke_width=w, stroke_linecap="round"))
    d.append(draw.Circle(cx, cy, 1.7, fill=col))


def mean_ltt(model_fn, age, seeds, ngrid=160):
    """Mean lineages-through-time over `seeds` forward runs, on a shared 0..age grid."""
    grid = [age * i / (ngrid - 1) for i in range(ngrid)]
    acc = [0.0] * ngrid
    for s in seeds:
        tree = simulate_species_tree(model_fn(), age=age, direction="forward",
                                     age_type="crown", seed=s)
        spans = [(n.parent.time, n.time) for n in tree.nodes() if n.parent is not None]
        for i, t in enumerate(grid):
            acc[i] += sum(1 for a, b in spans if a < t <= b)
    return grid, [c / len(seeds) for c in acc]


def render():
    d = draw.Drawing(W, H, origin=(0, 0))
    d.append(draw.Rectangle(0, 0, W, H, fill="white"))
    d.append(draw.Text("How rates work: how many clocks, how fast", FS_TITLE, W / 2, 48,
                       font_family=FONT, text_anchor="middle", font_weight="bold", fill=INK))
    d.append(draw.Text("the opportunity -- how many clocks -- is what makes growth exponential or linear",
                       FS_ANNOT, W / 2, 84, font_family=FONT, text_anchor="middle",
                       font_style="italic", fill="#555"))

    # ---- top row: the three opportunities as three columns ----
    cy = 200
    cols = [W / 2 - 430, W / 2, W / 2 + 430]

    # per lineage: three branch bars, a clock atop each
    x0 = cols[0]
    d.append(draw.Text("per lineage", FS_LABEL, x0, 130, font_family=FONT, text_anchor="middle",
                       font_weight="bold", fill=INK))
    for k, dx in enumerate((-52, 0, 52)):
        d.append(draw.Line(x0 + dx, cy + 46, x0 + dx, cy + 96, stroke=INK, stroke_width=3,
                           stroke_linecap="round"))
        clock(d, x0 + dx, cy + 14, 18)
    d.append(draw.Text("one clock on each lineage", FS_TICK, x0, cy + 128, font_family=FONT,
                       text_anchor="middle", fill="#555"))
    d.append(draw.Text("count tracks the tree", FS_TICK, x0, cy + 156, font_family=FONT,
                       text_anchor="middle", fill="#555"))

    # per copy: a genome box with three gene squares, a clock on each
    x1 = cols[1]
    d.append(draw.Text("per copy", FS_LABEL, x1, 130, font_family=FONT, text_anchor="middle",
                       font_weight="bold", fill=INK))
    d.append(draw.Rectangle(x1 - 92, cy + 44, 184, 54, rx=10, fill="none", stroke=GREY, stroke_width=2))
    for k, dx in enumerate((-56, 0, 56)):
        d.append(draw.Rectangle(x1 + dx - 16, cy + 56, 32, 30, rx=4, fill="#eee", stroke=INK,
                                stroke_width=1.6))
        clock(d, x1 + dx, cy + 14, 18)
    d.append(draw.Text("one clock on each gene copy", FS_TICK, x1, cy + 128, font_family=FONT,
                       text_anchor="middle", fill="#555"))
    d.append(draw.Text("count tracks the family", FS_TICK, x1, cy + 156, font_family=FONT,
                       text_anchor="middle", fill="#555"))

    # shared: one clock for a bracketed whole
    x2 = cols[2]
    d.append(draw.Text("shared", FS_LABEL, x2, 130, font_family=FONT, text_anchor="middle",
                       font_weight="bold", fill=INK))
    for dx in (-52, 0, 52):
        d.append(draw.Line(x2 + dx, cy + 52, x2 + dx, cy + 96, stroke=GREY, stroke_width=3,
                           stroke_linecap="round"))
    br = draw.Path(fill="none", stroke=INK, stroke_width=2)
    br.M(x2 - 74, cy + 44).Q(x2, cy + 24, x2 + 74, cy + 44)
    d.append(br)
    clock(d, x2, cy + 6, 20)
    d.append(draw.Text("one clock for the whole thing", FS_TICK, x2, cy + 128, font_family=FONT,
                       text_anchor="middle", fill="#555"))
    d.append(draw.Text("count is fixed at one", FS_TICK, x2, cy + 156, font_family=FONT,
                       text_anchor="middle", fill="#555"))

    # ---- bottom: the consequence, from real runs ----
    seeds = list(range(1, 41))
    AGE = 5.0
    g_exp, c_exp = mean_ltt(lambda: BirthDeath(1.0, 0.2), AGE, seeds)
    g_lin, c_lin = mean_ltt(lambda: SharedBirthDeath(1.0, 0.2), AGE, seeds)
    cmax = max(max(c_exp), max(c_lin)) * 1.06

    PX0, PX1 = 150, W - 150
    PY0, PY1 = 470, 858

    def X(t):
        return PX0 + t / AGE * (PX1 - PX0)

    def Y(c):
        return PY1 - c / cmax * (PY1 - PY0)

    d.append(draw.Line(PX0, PY1, PX1, PY1, stroke="#bdbdbd", stroke_width=1.3))
    d.append(draw.Line(PX0, PY0, PX0, PY1, stroke="#bdbdbd", stroke_width=1.3))

    # curves
    def curve(grid, counts, dash=None):
        pts = []
        for t, c in zip(grid, counts):
            pts += [X(t), Y(c)]
        kw = dict(close=False, fill="none", stroke=INK, stroke_width=3.0, stroke_linejoin="round")
        if dash:
            kw["stroke_dasharray"] = dash
        d.append(draw.Lines(*pts, **kw))

    curve(g_exp, c_exp)                 # per lineage / per copy -> exponential
    curve(g_lin, c_lin, dash="9,6")     # shared -> linear

    d.append(draw.Text("one clock per lineage (or per copy)", FS_LABEL, X(1.62), PY0 + 116,
                       font_family=FONT, text_anchor="middle", font_weight="bold", fill=INK))
    d.append(draw.Text("-> EXPONENTIAL", FS_LABEL, X(1.62), PY0 + 146,
                       font_family=FONT, text_anchor="middle", fill=INK))
    d.append(draw.Text("one shared clock -> LINEAR", FS_LABEL, X(4.1), Y(c_lin[-1]) - 28,
                       font_family=FONT, text_anchor="middle", font_weight="bold", fill="#555"))

    # y ticks
    for c in (0, round(cmax / 2 / 10) * 10, round(max(c_exp) / 10) * 10):
        d.append(draw.Text(str(int(c)), FS_TICK, PX0 - 12, Y(c), font_family=FONT, text_anchor="end",
                           dominant_baseline="central", fill="#777"))
    d.append(draw.Text("lineages (count)", FS_LABEL, PX0 - 44, (PY0 + PY1) / 2, font_family=FONT,
                       text_anchor="middle", fill="#555",
                       transform=f"rotate(-90 {PX0 - 44} {(PY0 + PY1) / 2})"))
    for i in range(6):
        t = AGE * i / 5
        d.append(draw.Line(X(t), PY1, X(t), PY1 + 6, stroke=INK, stroke_width=1.2))
        d.append(draw.Text(f"{t:.0f}", FS_TICK, X(t), PY1 + 24, font_family=FONT, text_anchor="middle",
                           fill="#777"))
    d.append(draw.Text("time", FS_LABEL, (PX0 + PX1) / 2, PY1 + 50, font_family=FONT,
                       text_anchor="middle", fill="#555"))
    d.append(draw.Text("same rate per clock (birth 1.0, death 0.2) -- only the number of clocks differs",
                       FS_ANNOT, (PX0 + PX1) / 2, PY0 - 12, font_family=FONT, text_anchor="middle",
                       font_style="italic", fill="#555"))

    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "rate_clocks.svg").write_text(d.as_svg(), encoding="utf-8")
    cairosvg.svg2png(bytestring=d.as_svg().encode(), write_to=str(OUT / "rate_clocks.png"), scale=300 / 72.0)
    print(f"wrote rate_clocks  (exp peak {max(c_exp):.1f}, linear peak {max(c_lin):.1f})")


if __name__ == "__main__":
    render()
