"""Figure: the couplings of the coevolve mode (S / T / G triangle) -- colour version.

A coevolution scenario is a set of directed edges over three processes -- species diversification
(S), a trait (T), and gene-family content (G). Each ordered pair is a distinct model. Between the
two directed arrows of every pair runs a straight DOUBLE-HEADED arrow: turning on *both* edges of a
pair is that pair's joint (bidirectional) model -- ClaSSE for traits<->species, and its two analogues.

The two edges that point *into* S are drawn heavy: an arrow into S makes the tree depend on the
coupled state, so the tree becomes an OUTPUT (grown jointly). The other four are overlays.

Colour: the three circles use the ZOMBI palette (teal / sage / terracotta). The six arrows use three
(dark, light) hue-PAIRS, each in the colour of ONE endpoint node so the palette closes symmetrically
-- green (echoing T) for species<->traits, teal (echoing S) for species<->genes, terracotta (echoing
G) for traits<->genes; the label is always the dark shade so it stays legible on white.

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

W, H = 1240, 860
R = 62                                        # node radius (30% larger than the 48 mono version)

# ZOMBI palette for the three processes -- teal / sage / amber, all dark enough for white text.
NODE_COL = {
    "S": "#2f7d84",   # species -- ZOMBI teal
    "T": "#5f7d34",   # traits  -- sage / olive green
    "G": "#b5654a",   # genes   -- terracotta
}

# Each node-PAIR gets one hue as a (dark, light) pair. The arrow is drawn in its shade; the LABEL is
# always the dark shade so it reads on white.
PAIR = {
    "ST": ("#446b28", "#9ec27c"),   # species <-> traits : green      -- echoes the T (traits) node
    "SG": ("#1f5c63", "#86bdc1"),   # species <-> genes  : teal       -- echoes the S (species) node
    "TG": ("#93492e", "#dc9f80"),   # traits  <-> genes  : terracotta -- echoes the G (genes) node
}
DBL = "#6b6b6b"                       # bidirectional (joint) double-headed arrows -- neutral grey
LW_INTO_S, LW_OVERLAY = 5.2, 3.2      # weight encodes into-S (tree is an output) vs overlay
FS_MODEL = FS_TICK - 1                 # model-name labels
FS_NODE = FS_TITLE + 6                 # node letter -- bigger, to match the larger circle
FS_NAME = FS_TICK - 2                  # node name (inside the circle, below the letter)


def edge_point(cx, cy, tx, ty, r):
    a = math.atan2(ty - cy, tx - cx)
    return cx + r * math.cos(a), cy + r * math.sin(a)


def _head(d, hx, hy, ax, ay, col, ah=28.0):
    """Filled arrowhead at (hx, hy) pointing away from (ax, ay)."""
    ang = math.atan2(hy - ay, hx - ax)
    d.append(draw.Lines(hx, hy, hx - ah * math.cos(ang - 0.42), hy - ah * math.sin(ang - 0.42),
                        hx - ah * math.cos(ang + 0.42), hy - ah * math.sin(ang + 0.42),
                        close=True, fill=col))


def arrow(d, c1, c2, side, bend, label, col, label_col, lw, label_gap=34, label_dx=0.0, label_dy=0.0):
    """Curved directed arrow c1->c2 bowing to `side`; drawn in `col`, labelled in `label_col`."""
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
    lx, ly = cx + ox * label_gap + label_dx, cy + oy * label_gap + label_dy
    d.append(draw.Text(label, FS_MODEL, lx, ly, font_family=FONT, text_anchor="middle",
                       dominant_baseline="central", fill=label_col, font_weight="bold"))


def biarrow(d, c1, c2, label, col=DBL, lw=3.4):
    """Straight DOUBLE-headed arrow along the chord = the pair's joint (bidirectional) model."""
    x1, y1 = c1
    x2, y2 = c2
    sx, sy = edge_point(x1, y1, x2, y2, R)
    ex, ey = edge_point(x2, y2, x1, y1, R)
    d.append(draw.Line(sx, sy, ex, ey, stroke=col, stroke_width=lw, stroke_linecap="round"))
    _head(d, sx, sy, ex, ey, col, ah=24.0)
    _head(d, ex, ey, sx, sy, col, ah=24.0)
    mx, my = (sx + ex) / 2, (sy + ey) / 2
    tw = 0.55 * FS_MODEL * len(label) + 16
    th = FS_MODEL + 12
    d.append(draw.Rectangle(mx - tw / 2, my - th / 2, tw, th, fill="white", stroke="none"))
    d.append(draw.Text(label, FS_MODEL, mx, my, font_family=FONT, text_anchor="middle",
                       dominant_baseline="central", fill=col, font_weight="bold"))


def node(d, cx, cy, letter, name):
    col = NODE_COL[letter]
    d.append(draw.Circle(cx, cy, R, fill=col, stroke="white", stroke_width=3.2))
    d.append(draw.Circle(cx, cy, R, fill="none", stroke=INK, stroke_width=1.3))     # thin ink ring
    d.append(draw.Text(letter, FS_NODE, cx, cy - 7, font_family=FONT, text_anchor="middle",
                       dominant_baseline="central", font_weight="bold", fill="white"))
    d.append(draw.Text(name, FS_NAME, cx, cy + 28, font_family=FONT, text_anchor="middle",
                       dominant_baseline="central", fill="white"))


def render():
    d = draw.Drawing(W, H, origin=(0, 0))
    d.append(draw.Rectangle(0, 0, W, H, fill="white"))
    d.append(draw.Text("Coevolve mode: directed edges and their joint models", FS_TITLE, W / 2, 50,
                       font_family=FONT, text_anchor="middle", font_weight="bold", fill=INK))

    S = (W / 2, 288)
    T = (W / 2 - 330, 668)
    G = (W / 2 + 330, 668)

    sd, sl = PAIR["ST"]     # (dark, light)
    gd, gl = PAIR["SG"]
    td, tl = PAIR["TG"]

    # directed edges; dark shade = one direction, light = the other; label always the dark shade.
    arrow(d, T, S, +1, 62, "trait-driven",     sd, sd, LW_INTO_S)                            # T->S into S
    arrow(d, S, T, -1, 62, "cladogenetic",     sl, sd, LW_OVERLAY, label_dy=-24)             # S->T overlay
    arrow(d, G, S, -1, 62, "key innovation",   gd, gd, LW_INTO_S, label_dx=-26, label_dy=10) # G->S into S
    arrow(d, S, G, +1, 62, "punctuational",    gl, gd, LW_OVERLAY)                           # S->G overlay
    arrow(d, T, G, +1, 56, "trait-linked",     td, td, LW_OVERLAY, label_gap=4)              # T->G
    arrow(d, G, T, -1, 56, "gene-conditioned", tl, td, LW_OVERLAY, label_gap=4)              # G->T

    # the joint (bidirectional) model of each pair
    biarrow(d, T, S, "ClaSSE")
    biarrow(d, G, S, "co-diversification")
    biarrow(d, T, G, "trait-gene feedback")

    node(d, *S, "S", "species")
    node(d, *T, "T", "traits")
    node(d, *G, "G", "genes")

    # legend, top-left: weight = into-S vs overlay; the double arrow = the joint model
    lx, ly = 64, 116
    d.append(draw.Line(lx, ly, lx + 46, ly, stroke=INK, stroke_width=LW_INTO_S, stroke_linecap="round"))
    d.append(draw.Text("into S: the tree is grown as an OUTPUT", FS_TICK, lx + 60, ly,
                       font_family=FONT, text_anchor="start", dominant_baseline="central", fill=INK))
    d.append(draw.Line(lx, ly + 32, lx + 46, ly + 32, stroke=INK, stroke_width=LW_OVERLAY,
                       stroke_linecap="round"))
    d.append(draw.Text("overlay: runs on a given tree", FS_TICK, lx + 60, ly + 32, font_family=FONT,
                       text_anchor="start", dominant_baseline="central", fill=INK))
    d.append(draw.Line(lx + 7, ly + 64, lx + 39, ly + 64, stroke=DBL, stroke_width=3.4,
                       stroke_linecap="round"))
    _head(d, lx + 7, ly + 64, lx + 39, ly + 64, DBL, ah=9.0)
    _head(d, lx + 39, ly + 64, lx + 7, ly + 64, DBL, ah=9.0)
    d.append(draw.Text("both edges at once: the joint (bidirectional) model", FS_TICK, lx + 60,
                       ly + 64, font_family=FONT, text_anchor="start", dominant_baseline="central",
                       fill=DBL))

    name = "coevolve_modes"
    out = OUT_DIR / name
    out.mkdir(parents=True, exist_ok=True)
    (out / f"{name}.svg").write_text(d.as_svg(), encoding="utf-8")
    cairosvg.svg2png(bytestring=d.as_svg().encode(), write_to=str(out / f"{name}.png"), scale=300 / 72.0)
    print(f"wrote {out}/{name}.svg / .png")


if __name__ == "__main__":
    render()
