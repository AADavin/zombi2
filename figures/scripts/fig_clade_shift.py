"""Figure: the clade-shift birth-death model — one clade switches rate regime.

``CladeShiftBirthDeath(birth, death, clade_shifts=[(age, birth, death), ...])``: the tree runs
at the background rates until a scheduled age, when one random then-alive lineage AND all its
descendants adopt a new ``(birth, death)`` regime — the discrete, hand-specified version of
clade rate heterogeneity (a key innovation that sparks a radiation, a lineage entering a
slow-down). Here the shifted clade speeds up and radiates.

ZOMBI2 picks the shifted lineage at random and doesn't tag it, so the tree here is grown by a
faithful re-implementation of the documented process so the shifted clade can be highlighted.

Monochrome (species-tree house style).  Run:  python figures/scripts/fig_clade_shift.py
"""

from __future__ import annotations

import math
from pathlib import Path

import cairosvg
import drawsvg as draw
import numpy as np

from zombi_style import FONT, INK

OUT = Path(__file__).resolve().parent.parent / "clade_shift"

W, H = 1180, 720
XL, XR = 96, 980
TREE_TOP, TREE_H = 130, 500
DASH = "6,5"
GREY = "#9a9a9a"

B0, D0 = 0.45, 0.08              # background (slow)
B1, D1 = 1.7, 0.10              # shifted regime (fast radiation)
T, T_SHIFT = 6.0, 2.1


class Node:
    __slots__ = ("time", "parent", "children", "is_extant", "shifted", "shift_here", "surv")

    def __init__(self, time, parent=None, shifted=False):
        self.time, self.parent, self.children = time, parent, []
        self.is_extant, self.shifted, self.shift_here = False, shifted, False


def grow(rng):
    while True:
        root = Node(0.0)
        a, b = Node(0.0, root), Node(0.0, root)
        root.children = [a, b]
        live, t, applied, sn = [a, b], 0.0, False, None
        ok = True
        while live:
            rate = [(B1 + D1) if n.shifted else (B0 + D0) for n in live]
            total = math.fsum(rate)
            dt = rng.exponential(1.0 / total)
            if not applied and t + dt >= T_SHIFT:
                sn = live[int(rng.integers(len(live)))]
                sn.shifted, sn.shift_here, applied = True, True, True
                t = T_SHIFT
                continue
            if t + dt >= T:
                for n in live:
                    n.time, n.is_extant = T, True
                break
            t += dt
            i = int(rng.choice(len(live), p=np.array(rate) / total))
            n = live[i]
            n.time = t
            bi, di = (B1, D1) if n.shifted else (B0, D0)
            if rng.random() < bi / (bi + di):
                c1, c2 = Node(t, n, n.shifted), Node(t, n, n.shifted)
                n.children = [c1, c2]
                live[i] = c1
                live.append(c2)
            else:
                live.pop(i)
            if len(live) > 500:
                ok = False
                break
        if not ok or sn is None:
            continue
        _mark(root)
        shifted_tips = sum(1 for n in _leaves(root) if n.shifted)
        if 26 <= _nleaf(root) <= 44 and shifted_tips >= 12 and root.children and sn.parent is not None:
            return root, sn


def _leaves(n):
    return [n] if not n.children else [x for c in n.children for x in _leaves(c)]


def _nleaf(n):
    return 1 if not n.children else sum(_nleaf(c) for c in n.children)


def _mark(n):
    n.surv = n.is_extant if not n.children else any([_mark(c) for c in n.children])
    return n.surv


def nodes(root):
    out, st = [], [root]
    while st:
        n = st.pop(); out.append(n); st.extend(n.children)
    return out


def layout(root):
    order, ys = [0], {}

    def rec(n):
        if not n.children:
            ys[n] = order[0]; order[0] += 1
        else:
            for c in sorted(n.children, key=_nleaf):
                rec(c)
            ys[n] = sum(ys[c] for c in n.children) / len(n.children)
    rec(root)
    return ys, order[0]


def render():
    rng = np.random.default_rng(2)
    root, sn = grow(rng)
    ys, nleaf = layout(root)
    dy = min(11.5, TREE_H / max(1, nleaf - 1))
    top = TREE_TOP + (TREE_H - dy * (nleaf - 1)) / 2
    present = T

    def X(t):
        return XL + t / present * (XR - XL)

    def Y(n):
        return top + ys[n] * dy

    d = draw.Drawing(W, H, origin=(0, 0))
    d.append(draw.Rectangle(0, 0, W, H, fill="white"))
    d.append(draw.Text("The clade-shift model — one clade switches to a new rate regime", 20, 40, 40,
                       font_family=FONT, text_anchor="start", font_weight="bold", fill=INK))
    d.append(draw.Text("at a set age, a random lineage and all its descendants adopt new "
                       "(birth, death) rates — here a burst of speciation.  ZOMBI2: "
                       "CladeShiftBirthDeath(b, d, clade_shifts=[(age, b, d)])", 13, 40, 64,
                       font_family=FONT, text_anchor="start", fill="#777"))

    alln = nodes(root)
    shifted_tip_ys = [ys[n] for n in _leaves(root) if n.shifted]
    y0, y1 = top + min(shifted_tip_ys) * dy, top + max(shifted_tip_ys) * dy
    d.append(draw.Rectangle(X(T_SHIFT) - 4, y0 - 10, XR - X(T_SHIFT) + 12, y1 - y0 + 20,
                            fill="#f0f0f0"))                          # shifted-clade backdrop

    def seg(x1, y1_, x2, y2, n):
        col = INK if n.shifted else GREY
        wid = 2.4 if n.shifted else 1.6
        kw = dict(stroke=col, stroke_width=wid, stroke_linecap="round")
        if not n.surv:
            kw = dict(stroke=col, stroke_width=1.4, stroke_dasharray=DASH, stroke_linecap="butt")
        d.append(draw.Line(x1, y1_, x2, y2, **kw))

    for n in alln:
        px = X(0) - 12 if n.parent is None else X(n.parent.time)
        if n.parent is not None:
            seg(px, Y(n), X(n.time), Y(n), n)
        for c in n.children:
            seg(X(n.time), Y(n), X(n.time), Y(c), c)

    # shift marker
    sy = Y(sn)
    d.append(draw.Line(X(T_SHIFT), y0 - 14, X(T_SHIFT), y1 + 14, stroke=INK, stroke_width=1.4,
                       stroke_dasharray="5,4"))
    d.append(draw.Circle(X(T_SHIFT), sy, 6.5, fill="white", stroke=INK, stroke_width=2.2))
    d.append(draw.Lines(X(T_SHIFT) - 3, sy, X(T_SHIFT) + 2, sy - 3.5, X(T_SHIFT) + 2, sy + 3.5,
                        close=True, fill=INK))                        # a little "play" mark = speeds up
    d.append(draw.Text("rate shift", 13, X(T_SHIFT), y0 - 22, font_family=FONT, text_anchor="middle",
                       font_weight="bold", fill=INK))
    d.append(draw.Text("this clade now diversifies fast", 12, X(T_SHIFT) + 14, y1 + 26,
                       font_family=FONT, text_anchor="start", fill="#555", font_style="italic"))

    # legend
    lx, ly = XL + 6, TREE_TOP - 44
    d.append(draw.Line(lx, ly, lx + 26, ly, stroke=GREY, stroke_width=1.6, stroke_linecap="round"))
    d.append(draw.Text("background rate", 12, lx + 32, ly, font_family=FONT, text_anchor="start",
                       dominant_baseline="central", fill=INK))
    d.append(draw.Line(lx + 180, ly, lx + 206, ly, stroke=INK, stroke_width=2.4, stroke_linecap="round"))
    d.append(draw.Text("shifted clade", 12, lx + 212, ly, font_family=FONT, text_anchor="start",
                       dominant_baseline="central", fill=INK))
    d.append(draw.Line(lx + 350, ly, lx + 376, ly, stroke=GREY, stroke_width=1.4, stroke_dasharray=DASH))
    d.append(draw.Text("extinct", 12, lx + 382, ly, font_family=FONT, text_anchor="start",
                       dominant_baseline="central", fill=INK))

    ya = top + dy * (nleaf - 1) + 26
    for i in range(6):
        t = present * i / 5
        d.append(draw.Line(X(t), ya, X(t), ya + 5, stroke=INK, stroke_width=1.2))
        d.append(draw.Text(f"{t:.1f}", 11, X(t), ya + 18, font_family=FONT, text_anchor="middle", fill="#777"))
    d.append(draw.Text("time (root -> present)", 12.5, (XL + XR) / 2, ya + 38, font_family=FONT,
                       text_anchor="middle", fill="#777"))

    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "clade_shift.svg").write_text(d.as_svg(), encoding="utf-8")
    cairosvg.svg2png(bytestring=d.as_svg().encode(), write_to=str(OUT / "clade_shift.png"), scale=300 / 72.0)
    print(f"wrote clade_shift  ({nleaf} tips)")


if __name__ == "__main__":
    render()
