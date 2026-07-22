"""Figure: ZOMBI2's discrete-trait models as continuous-time Markov chains.

Every discrete trait in ZOMBI2 (the ``Mk`` family, correlated binary, hidden-state Mk)
is a continuous-time Markov chain over a set of states with an instantaneous rate matrix
``Q``. This figure draws the classic state-and-arrows diagram for six scenarios, so the
structure of each model (which transitions are allowed, and how the rates relate) is
visible at a glance.

Rendered in colour and B&W.  Run:  python figures/scripts/fig_trait_markov.py
"""

from __future__ import annotations

import math
from pathlib import Path

import cairosvg
import drawsvg as draw

from zombi_style import FONT, INK

OUT_DIR = Path(__file__).resolve().parent.parent

W, H = 1240, 730
NR = 25                                                   # node radius
COLOR = {0: "#4477AA", 1: "#EE6677", 2: "#228833"}
BW = {0: "#4a4a4a", 1: "#9a9a9a", 2: "#cfcfcf"}
NEUTRAL = "#e9e9e9"
MODE = "color"


def scol(i):
    return (BW if MODE == "bw" else COLOR)[i]


def _lum(h):
    return 0.299 * int(h[1:3], 16) + 0.587 * int(h[3:5], 16) + 0.114 * int(h[5:7], 16)


# --------------------------------------------------------------------------- primitives
def node(d, cx, cy, label, fill, r=NR):
    d.append(draw.Circle(cx, cy, r, fill=fill, stroke=INK, stroke_width=1.7))
    fs = 16 if len(str(label)) == 1 else 13
    d.append(draw.Text(str(label), fs, cx, cy, font_family=FONT, text_anchor="middle",
                       dominant_baseline="central", fill="white" if _lum(fill) < 150 else INK,
                       font_weight="bold"))


def _edge(cx, cy, tx, ty, r):
    a = math.atan2(ty - cy, tx - cx)
    return cx + r * math.cos(a), cy + r * math.sin(a)


def arc(d, c1, c2, bend, color, label="", r=NR, lw=2.3, lab_off=15):
    x1, y1 = c1
    x2, y2 = c2
    mx, my = (x1 + x2) / 2, (y1 + y2) / 2
    dx, dy = x2 - x1, y2 - y1
    L = math.hypot(dx, dy)
    px, py = -dy / L, dx / L
    cx, cy = mx + px * bend, my + py * bend
    sx, sy = _edge(x1, y1, cx, cy, r)
    ex, ey = _edge(x2, y2, cx, cy, r)
    p = draw.Path(fill="none", stroke=color, stroke_width=lw)
    p.M(sx, sy).Q(cx, cy, ex, ey)
    d.append(p)
    aa, ah = math.atan2(ey - cy, ex - cx), 8.5
    d.append(draw.Lines(ex, ey, ex - ah * math.cos(aa - 0.42), ey - ah * math.sin(aa - 0.42),
                        ex - ah * math.cos(aa + 0.42), ey - ah * math.sin(aa + 0.42),
                        close=True, fill=color))
    if label:
        d.append(draw.Text(label, 13.5, cx + px * lab_off, cy + py * lab_off, font_family=FONT,
                           text_anchor="middle", dominant_baseline="central", fill=INK))


def bidir(d, c1, c2, bend, color, lab1="", lab2=""):
    arc(d, c1, c2, bend, color, lab1)
    arc(d, c2, c1, bend, color, lab2)


def title(d, cx, y, text, sub):
    d.append(draw.Text(text, 16, cx, y, font_family=FONT, text_anchor="middle",
                       font_weight="bold", fill=INK))
    d.append(draw.Text(sub, 12, cx, y + 18, font_family=FONT, text_anchor="middle", fill="#777"))


# --------------------------------------------------------------------------- scenarios
def er(d, cx, cy):
    title(d, cx, cy - 112, "Equal rates (ER)", "one shared rate q")
    a, b = (cx - 66, cy), (cx + 66, cy)
    bidir(d, a, b, 30, INK, "q", "q")
    node(d, *a, "0", scol(0))
    node(d, *b, "1", scol(1))


def ard(d, cx, cy):
    title(d, cx, cy - 112, "All rates different (ARD)", "each direction its own rate")
    a, b = (cx - 66, cy), (cx + 66, cy)
    bidir(d, a, b, 30, INK, "q01", "q10")
    node(d, *a, "0", scol(0))
    node(d, *b, "1", scol(1))


def sym(d, cx, cy):
    title(d, cx, cy - 116, "Symmetric (SYM)", "each pair shares one rate")
    s0, s1, s2 = (cx, cy - 52), (cx - 66, cy + 46), (cx + 66, cy + 46)
    bidir(d, s0, s1, 16, INK, "a", "")
    bidir(d, s1, s2, 16, INK, "b", "")
    bidir(d, s0, s2, 16, INK, "", "c")
    node(d, *s0, "0", scol(0))
    node(d, *s1, "1", scol(1))
    node(d, *s2, "2", scol(2))


def ordered(d, cx, cy):
    title(d, cx, cy - 112, "Ordered / meristic", "neighbour jumps only (no 0-2)")
    s0, s1, s2 = (cx - 92, cy), (cx, cy), (cx + 92, cy)
    bidir(d, s0, s1, 26, INK, "q", "q")
    bidir(d, s1, s2, 26, INK, "q", "q")
    node(d, *s0, "0", scol(0))
    node(d, *s1, "1", scol(1))
    node(d, *s2, "2", scol(2))


def pagel(d, cx, cy):
    title(d, cx, cy - 116, "Correlated binary (Pagel)", "each trait's rate depends on the other")
    tl, tr = (cx - 66, cy - 48), (cx + 66, cy - 48)
    bl, br = (cx - 66, cy + 48), (cx + 66, cy + 48)
    bidir(d, tl, tr, 15, INK, "Y", "")
    bidir(d, bl, br, 15, INK, "", "Y")
    bidir(d, tl, bl, 15, INK, "X", "")
    bidir(d, tr, br, 15, INK, "", "X")
    for c, lab in ((tl, "00"), (tr, "01"), (bl, "10"), (br, "11")):
        node(d, *c, lab, NEUTRAL)


def hidden(d, cx, cy):
    title(d, cx, cy - 116, "Hidden-state Mk", "observed rate depends on a hidden class")
    tl, tr = (cx - 66, cy - 48), (cx + 66, cy - 48)
    bl, br = (cx - 66, cy + 48), (cx + 66, cy + 48)
    bidir(d, tl, tr, 15, INK, "slow", "")
    bidir(d, bl, br, 15, INK, "", "fast")
    bidir(d, tl, bl, 15, "#999", "", "")
    bidir(d, tr, br, 15, "#999", "", "")
    for c, lab in ((tl, "slow"), (bl, "fast")):
        d.append(draw.Text(lab, 12, cx - 118, c[1], font_family=FONT, text_anchor="middle",
                           dominant_baseline="central", font_style="italic", fill="#777"))
    for c, lab in ((tl, "0"), (tr, "1"), (bl, "0"), (br, "1")):
        node(d, *c, lab, scol(int(lab)))


# --------------------------------------------------------------------------- render
def render(mode):
    global MODE
    MODE = mode
    d = draw.Drawing(W, H, origin=(0, 0))
    d.append(draw.Rectangle(0, 0, W, H, fill="white"))
    d.append(draw.Text("Discrete-trait models as continuous-time Markov chains", 20, 40, 40,
                       font_family=FONT, text_anchor="start", font_weight="bold", fill=INK))
    d.append(draw.Text("states = circles,  allowed transitions = arrows,  instantaneous rates = labels "
                       "(the off-diagonals of Q)", 13.5, 40, 64, font_family=FONT,
                       text_anchor="start", fill="#777"))
    cols = [250, 620, 990]
    r1, r2 = 270, 560
    er(d, cols[0], r1)
    sym(d, cols[1], r1)
    ordered(d, cols[2], r1)
    ard(d, cols[0], r2)
    pagel(d, cols[1], r2)
    hidden(d, cols[2], r2)

    name = "trait_markov" if mode == "color" else "trait_markov_bw"
    out = OUT_DIR / name
    out.mkdir(parents=True, exist_ok=True)
    (out / f"{name}.svg").write_text(d.as_svg(), encoding="utf-8")
    cairosvg.svg2png(bytestring=d.as_svg().encode(), write_to=str(out / f"{name}.png"), scale=300 / 72.0)


def main():
    for mode in ("color", "bw"):
        render(mode)
    print("wrote trait_markov (+_bw)")


if __name__ == "__main__":
    main()
