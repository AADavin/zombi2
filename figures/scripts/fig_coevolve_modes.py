"""Figure: the couplings of the coevolve mode (S / T / G triangle).

A coevolution scenario is a set of directed edges over three processes — species diversification (S),
a trait (T), and gene-family content (G). Each ordered pair is a distinct model. Between the two
directed arrows of every pair runs a straight DOUBLE-HEADED arrow: turning on *both* edges of a pair
is that pair's joint (bidirectional) model — ClaSSE for traits<->species, and its two analogues.

The two edges that point *into* S are drawn heavy: an arrow into S makes the tree depend on the
coupled state, so the tree becomes an OUTPUT (grown jointly). The other four are overlays on a given
tree. Each directed arrow is coloured to match its label so the six labels stay legible.

House style here relaxes the B&W rule on Adrián's request: colour = arrow identity (label match).

Run:  /Users/aadria/miniconda3/bin/python figures/scripts/fig_coevolve_modes.py
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import cairosvg
import drawsvg as draw

from zombi_style import FONT, INK, MUTED, FS_TITLE, FS_LABEL, FS_ANNOT, FS_TICK

OUT_DIR = Path(__file__).resolve().parent.parent

W, H = 1180, 820
R = 48                                        # node radius

# One colour per directed arrow, matched to its label (Okabe-Ito, colour-blind safe). Warm/cool
# groups the three node-pairs: cool blues = species<->traits, warm = species<->genes, green/purple =
# traits<->genes. The bidirectional joint arrows are neutral grey (a "combination", not a 7th edge).
COL = {
    "trait-driven":     "#0072B2",   # T -> S  (into S)
    "cladogenetic":     "#56B4E9",   # S -> T
    "key innovation":   "#D55E00",   # G -> S  (into S)
    "punctuational":    "#E69F00",   # S -> G
    "trait-linked":     "#009E73",   # T -> G
    "gene-conditioned": "#8E5AA8",   # G -> T  (purple; pairs with the trait-linked green)
}
DBL = "#555555"                       # bidirectional (joint) double-headed arrows
LW_INTO_S, LW_OVERLAY = 4.6, 2.8      # weight encodes into-S (tree is an output) vs overlay


def edge_point(cx, cy, tx, ty, r):
    a = math.atan2(ty - cy, tx - cx)
    return cx + r * math.cos(a), cy + r * math.sin(a)


def _head(d, hx, hy, ax, ay, col, ah=13.0):
    """Filled arrowhead at (hx, hy) pointing away from (ax, ay)."""
    ang = math.atan2(hy - ay, hx - ax)
    d.append(draw.Lines(hx, hy, hx - ah * math.cos(ang - 0.42), hy - ah * math.sin(ang - 0.42),
                        hx - ah * math.cos(ang + 0.42), hy - ah * math.sin(ang + 0.42),
                        close=True, fill=col))


def arrow(d, c1, c2, side, bend, label, col, lw, label_gap=34):
    """Curved directed arrow c1->c2; bows to `side` (+1/-1). Coloured `col`, matching its label.

    `label_gap` is how far past the curve's control point the label sits (smaller = closer to
    the arrow); the near-horizontal T<->G edges use a small gap so their labels hug their arrows.
    """
    x1, y1 = c1
    x2, y2 = c2
    mx, my = (x1 + x2) / 2, (y1 + y2) / 2
    dx, dy = x2 - x1, y2 - y1
    L = math.hypot(dx, dy)
    px, py = -dy / L, dx / L
    if py > 0 or (abs(py) < 1e-9 and px > 0):     # canonical normal (up / left)
        px, py = -px, -py
    ox, oy = px * side, py * side
    cx, cy = mx + ox * bend, my + oy * bend
    sx, sy = edge_point(x1, y1, cx, cy, R)
    ex, ey = edge_point(x2, y2, cx, cy, R)
    p = draw.Path(fill="none", stroke=col, stroke_width=lw, stroke_linecap="round")
    p.M(sx, sy).Q(cx, cy, ex, ey)
    d.append(p)
    _head(d, ex, ey, cx, cy, col)
    lx, ly = cx + ox * label_gap, cy + oy * label_gap
    d.append(draw.Text(label, FS_TICK, lx, ly, font_family=FONT, text_anchor="middle",
                       dominant_baseline="central", fill=col, font_weight="bold"))


def biarrow(d, c1, c2, label, col=DBL, lw=3.2):
    """Straight DOUBLE-headed arrow along the c1-c2 chord = the pair's joint (bidirectional) model.
    Label sits centred on the chord over a white patch that breaks the line."""
    x1, y1 = c1
    x2, y2 = c2
    sx, sy = edge_point(x1, y1, x2, y2, R)
    ex, ey = edge_point(x2, y2, x1, y1, R)
    d.append(draw.Line(sx, sy, ex, ey, stroke=col, stroke_width=lw, stroke_linecap="round"))
    _head(d, sx, sy, ex, ey, col, ah=12.0)
    _head(d, ex, ey, sx, sy, col, ah=12.0)
    mx, my = (sx + ex) / 2, (sy + ey) / 2
    tw = 0.55 * FS_TICK * len(label) + 16
    th = FS_TICK + 12
    d.append(draw.Rectangle(mx - tw / 2, my - th / 2, tw, th, fill="white", stroke="none"))
    d.append(draw.Text(label, FS_TICK, mx, my, font_family=FONT, text_anchor="middle",
                       dominant_baseline="central", fill=col, font_weight="bold"))


def node(d, cx, cy, letter, name):
    d.append(draw.Circle(cx, cy, R, fill="white", stroke=INK, stroke_width=2.6))
    d.append(draw.Text(letter, FS_TITLE, cx, cy - 4, font_family=FONT, text_anchor="middle",
                       dominant_baseline="central", font_weight="bold", fill=INK))
    d.append(draw.Text(name, FS_TICK, cx, cy + 24, font_family=FONT, text_anchor="middle",
                       dominant_baseline="central", fill=MUTED))


def render():
    d = draw.Drawing(W, H, origin=(0, 0))
    d.append(draw.Rectangle(0, 0, W, H, fill="white"))
    d.append(draw.Text("Coevolve mode: directed edges and their joint models", FS_TITLE, W / 2, 48,
                       font_family=FONT, text_anchor="middle", font_weight="bold", fill=INK))

    S = (W / 2, 270)
    T = (W / 2 - 300, 630)
    G = (W / 2 + 300, 630)

    # directed edges (bow outward), coloured to match their labels; into-S edges heavy
    arrow(d, T, S, +1, 58, "trait-driven",     COL["trait-driven"],     LW_INTO_S)
    arrow(d, S, T, -1, 58, "cladogenetic",     COL["cladogenetic"],     LW_OVERLAY)
    arrow(d, G, S, -1, 58, "key innovation",   COL["key innovation"],   LW_INTO_S)
    arrow(d, S, G, +1, 58, "punctuational",    COL["punctuational"],    LW_OVERLAY)
    arrow(d, T, G, +1, 52, "trait-linked",     COL["trait-linked"],     LW_OVERLAY, label_gap=4)
    arrow(d, G, T, -1, 52, "gene-conditioned", COL["gene-conditioned"], LW_OVERLAY, label_gap=4)

    # the joint (bidirectional) model of each pair: both edges at once
    biarrow(d, T, S, "ClaSSE")                       # traits:species + species:traits
    biarrow(d, G, S, "co-diversification")           # genes:species + species:genes  (into S)
    biarrow(d, T, G, "trait-gene feedback")          # traits:genes + genes:traits    (overlay)

    node(d, *S, "S", "species")
    node(d, *T, "T", "traits")
    node(d, *G, "G", "genes")

    # legend, top-left: weight encodes into-S vs overlay; the double arrow = the joint model
    lx, ly = 60, 108
    d.append(draw.Line(lx, ly, lx + 44, ly, stroke=INK, stroke_width=LW_INTO_S, stroke_linecap="round"))
    d.append(draw.Text("into S: the tree is grown as an OUTPUT", FS_TICK, lx + 56, ly,
                       font_family=FONT, text_anchor="start", dominant_baseline="central", fill=INK))
    d.append(draw.Line(lx, ly + 30, lx + 44, ly + 30, stroke=INK, stroke_width=LW_OVERLAY,
                       stroke_linecap="round"))
    d.append(draw.Text("overlay: runs on a given tree", FS_TICK, lx + 56, ly + 30, font_family=FONT,
                       text_anchor="start", dominant_baseline="central", fill=INK))
    d.append(draw.Line(lx + 6, ly + 60, lx + 38, ly + 60, stroke=DBL, stroke_width=3.2,
                       stroke_linecap="round"))
    _head(d, lx + 6, ly + 60, lx + 38, ly + 60, DBL, ah=9.0)
    _head(d, lx + 38, ly + 60, lx + 6, ly + 60, DBL, ah=9.0)
    d.append(draw.Text("both edges at once: the joint (bidirectional) model", FS_TICK, lx + 56,
                       ly + 60, font_family=FONT, text_anchor="start", dominant_baseline="central",
                       fill=DBL))

    d.append(draw.Text("each arrow driver -> target = one model; a double arrow = both edges "
                       "(--couple driver:target)", FS_TICK, W / 2, H - 24, font_family=FONT,
                       text_anchor="middle", fill=MUTED))

    name = "coevolve_modes"
    out = OUT_DIR / name
    out.mkdir(parents=True, exist_ok=True)
    (out / f"{name}.svg").write_text(d.as_svg(), encoding="utf-8")
    cairosvg.svg2png(bytestring=d.as_svg().encode(), write_to=str(out / f"{name}.png"), scale=300 / 72.0)
    print(f"wrote {out}/{name}.svg / .png")


if __name__ == "__main__":
    render()
