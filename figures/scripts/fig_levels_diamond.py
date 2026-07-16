"""Figure: the four levels of evolution, as a single-panel diamond.

The same diamond the coevolution part fills with coupling edges — here it just shows the four levels
and how the default *pipeline* composes them: each level is simulated *along* (conditioned on) the one
above. Species is the substrate (the timeline); traits and gene content are characters that ride the
species tree; sequences ride the gene trees below. Coupling these levels instead of layering them is
coevolution (Part VI).

Run:  python figures/scripts/fig_levels_diamond.py
"""

from __future__ import annotations

import math
from pathlib import Path

import cairosvg
import drawsvg as draw

from zombi_style import FONT, INK, FS_TITLE, FS_LABEL, FS_TICK

OUT = Path(__file__).resolve().parent.parent / "levels_diamond"

W, H = 1120, 940
R = 66

NODE_COL = {
    "S": "#2f7d84",   # species    -- teal
    "T": "#5f7d34",   # traits     -- sage
    "G": "#b5654a",   # genes      -- terracotta
    "Σ": "#5a4e86",   # sequences  -- violet
}
ARROW = "#6b6b6b"
LW = 3.2


def edge_point(cx, cy, tx, ty, r):
    a = math.atan2(ty - cy, tx - cx)
    return cx + r * math.cos(a), cy + r * math.sin(a)


def _head(d, hx, hy, ax, ay, col, ah=24.0):
    ang = math.atan2(hy - ay, hx - ax)
    d.append(draw.Lines(hx, hy, hx - ah * math.cos(ang - 0.42), hy - ah * math.sin(ang - 0.42),
                        hx - ah * math.cos(ang + 0.42), hy - ah * math.sin(ang + 0.42),
                        close=True, fill=col))


def arrow(d, c1, c2, label, dx=0.0, dy=0.0):
    x1, y1 = c1
    x2, y2 = c2
    sx, sy = edge_point(x1, y1, x2, y2, R)
    ex, ey = edge_point(x2, y2, x1, y1, R)
    d.append(draw.Line(sx, sy, ex, ey, stroke=ARROW, stroke_width=LW, stroke_linecap="round"))
    _head(d, ex, ey, sx, sy, ARROW)
    mx, my = (sx + ex) / 2, (sy + ey) / 2
    d.append(draw.Text(label, FS_TICK, mx + dx, my + dy, font_family=FONT, text_anchor="middle",
                       dominant_baseline="central", fill="#555", font_style="italic"))


def node(d, cx, cy, letter, name):
    col = NODE_COL[letter]
    d.append(draw.Circle(cx, cy, R, fill=col, stroke="white", stroke_width=3.2))
    d.append(draw.Circle(cx, cy, R + 6, fill="none", stroke=INK, stroke_width=1.6))
    d.append(draw.Text(letter, FS_TITLE + 8, cx, cy - 8, font_family=FONT, text_anchor="middle",
                       dominant_baseline="central", font_weight="bold", fill="white"))
    d.append(draw.Text(name, FS_TICK - 2, cx, cy + 28, font_family=FONT, text_anchor="middle",
                       dominant_baseline="central", fill="white"))


def render():
    d = draw.Drawing(W, H, origin=(0, 0))
    d.append(draw.Rectangle(0, 0, W, H, fill="white"))
    d.append(draw.Text("The four levels of evolution", FS_TITLE, W / 2, 48,
                       font_family=FONT, text_anchor="middle", font_weight="bold", fill=INK))
    d.append(draw.Text("each level is simulated along (conditioned on) the one above", FS_TICK,
                       W / 2, 82, font_family=FONT, text_anchor="middle", font_style="italic",
                       fill="#777"))

    S = (W / 2, 232)
    T = (W / 2 - 300, 512)
    G = (W / 2 + 300, 512)
    Q = (W / 2, 800)

    # the default pipeline: characters ride the species tree; sequences ride the gene trees
    arrow(d, S, T, "along the tree", dx=-46, dy=-6)
    arrow(d, S, G, "along the tree", dx=46, dy=-6)
    arrow(d, G, Q, "along the gene trees", dx=104, dy=0)

    node(d, *S, "S", "species")
    node(d, *T, "T", "traits")
    node(d, *G, "G", "genomes")
    node(d, *Q, "Σ", "sequences")

    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "levels_diamond.svg").write_text(d.as_svg(), encoding="utf-8")
    cairosvg.svg2png(bytestring=d.as_svg().encode(), write_to=str(OUT / "levels_diamond.png"), scale=300 / 72.0)
    print(f"wrote {OUT}/levels_diamond.svg / .png")


if __name__ == "__main__":
    render()
