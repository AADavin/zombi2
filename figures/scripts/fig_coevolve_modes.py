"""Figure: the six directed couplings of the coevolve mode (S / T / G triangle).

A coevolution scenario is a set of directed edges over three processes — species diversification (S),
a trait (T), and gene-family content (G). Each ordered pair is a distinct model. The two edges that
point *into* S are drawn heavy: an arrow into S makes the tree depend on the coupled state, so the
tree becomes an OUTPUT (forward-only). The other four are overlays on a given tree.

House style: B&W, one centered title, ASCII text.

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

W, H = 1180, 790
R = 48                                        # node radius
GREY = "#8a8a8a"


def edge_point(cx, cy, tx, ty, r):
    a = math.atan2(ty - cy, tx - cx)
    return cx + r * math.cos(a), cy + r * math.sin(a)


def arrow(d, c1, c2, side, bend, label, into_s):
    """Curved arrow c1->c2; bows to `side` (+1/-1), direction-independent. Heavy if into_s."""
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
    col = INK if into_s else GREY
    lw = 5.0 if into_s else 2.6
    p = draw.Path(fill="none", stroke=col, stroke_width=lw, stroke_linecap="round")
    p.M(sx, sy).Q(cx, cy, ex, ey)
    d.append(p)
    ang, ah = math.atan2(ey - cy, ex - cx), 13.0
    d.append(draw.Lines(ex, ey, ex - ah * math.cos(ang - 0.42), ey - ah * math.sin(ang - 0.42),
                        ex - ah * math.cos(ang + 0.42), ey - ah * math.sin(ang + 0.42),
                        close=True, fill=col))
    lx, ly = cx + ox * 34, cy + oy * 34
    d.append(draw.Text(label, FS_TICK, lx, ly, font_family=FONT, text_anchor="middle",
                       dominant_baseline="central", fill=INK if into_s else MUTED,
                       font_weight="bold" if into_s else "normal"))


def node(d, cx, cy, letter, name):
    d.append(draw.Circle(cx, cy, R, fill="white", stroke=INK, stroke_width=2.6))
    d.append(draw.Text(letter, FS_TITLE, cx, cy - 4, font_family=FONT, text_anchor="middle",
                       dominant_baseline="central", font_weight="bold", fill=INK))
    d.append(draw.Text(name, FS_TICK, cx, cy + 24, font_family=FONT, text_anchor="middle",
                       dominant_baseline="central", fill=MUTED))


def render():
    d = draw.Drawing(W, H, origin=(0, 0))
    d.append(draw.Rectangle(0, 0, W, H, fill="white"))
    d.append(draw.Text("The six couplings of coevolve mode", FS_TITLE, W / 2, 48, font_family=FONT,
                       text_anchor="middle", font_weight="bold", fill=INK))

    S = (W / 2, 250)
    T = (W / 2 - 300, 600)
    G = (W / 2 + 300, 600)

    # S <-> T (left side): T->S is INTO S (heavy)
    arrow(d, T, S, +1, 44, "SSE", into_s=True)
    arrow(d, S, T, -1, 44, "cladogenetic", into_s=False)
    # S <-> G (right side): G->S is INTO S (heavy)
    arrow(d, G, S, -1, 44, "key innovation", into_s=True)
    arrow(d, S, G, +1, 44, "punctuational", into_s=False)
    # T <-> G (bottom): both overlays
    arrow(d, T, G, +1, 38, "trait-linked", into_s=False)
    arrow(d, G, T, -1, 38, "gene-conditioned", into_s=False)

    node(d, *S, "S", "species")
    node(d, *T, "T", "traits")
    node(d, *G, "G", "genes")

    # legend, top-left
    lx, ly = 60, 110
    d.append(draw.Line(lx, ly, lx + 44, ly, stroke=INK, stroke_width=5.0, stroke_linecap="round"))
    d.append(draw.Text("into S: the tree is an output (forward-only)", FS_TICK, lx + 56, ly,
                       font_family=FONT, text_anchor="start", dominant_baseline="central", fill=INK))
    d.append(draw.Line(lx, ly + 30, lx + 44, ly + 30, stroke=GREY, stroke_width=2.6, stroke_linecap="round"))
    d.append(draw.Text("overlay: runs on a given tree", FS_TICK, lx + 56, ly + 30, font_family=FONT,
                       text_anchor="start", dominant_baseline="central", fill=MUTED))
    d.append(draw.Text("each arrow driver -> target = one model (--couple driver:target)", FS_TICK,
                       W / 2, H - 24, font_family=FONT, text_anchor="middle", fill=MUTED))

    name = "coevolve_modes"
    out = OUT_DIR / name
    out.mkdir(parents=True, exist_ok=True)
    (out / f"{name}.svg").write_text(d.as_svg(), encoding="utf-8")
    cairosvg.svg2png(bytestring=d.as_svg().encode(), write_to=str(out / f"{name}.png"), scale=300 / 72.0)
    print(f"wrote {out}/{name}.svg / .png")


if __name__ == "__main__":
    render()
