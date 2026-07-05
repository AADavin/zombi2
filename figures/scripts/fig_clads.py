"""Figure: the ClaDS model -- every lineage carries its own, drifting speciation rate.

ClaDS (Maliet, Hartig & Morlon 2019): at each speciation the two daughters inherit the
parent's rate times an independent lognormal jump (a per-lineage lognormal trend plus
jump spread; extinction follows a constant turnover). Rates therefore drift lineage by
lineage down the tree, so different clades diversify at very different tempos.

ZOMBI2 doesn't expose the per-lineage rates on the output tree, so here the tree is grown by
a faithful re-implementation of the documented ClaDS process (same formulas) so each branch
can be painted by its rate. Branch shade AND width encode the rate (dark+thick = fast).

Monochrome (species-tree house style).  Run:  python figures/scripts/fig_clads.py
"""

from __future__ import annotations

import math
from pathlib import Path

import cairosvg
import drawsvg as draw
import numpy as np

from zombi_style import FONT, INK, FS_TITLE, FS_LABEL, FS_ANNOT, FS_TICK

OUT = Path(__file__).resolve().parent.parent / "clads"

W, H = 1400, 760
XL, XR = 100, 1300
TREE_TOP, TREE_H = 176, 496

LAM0, ALPHA, SIGMA, TURN = 1.0, 0.9, 0.55, 0.1


class Node:
    __slots__ = ("time", "rate", "parent", "children", "is_extant")

    def __init__(self, time, rate, parent=None):
        self.time, self.rate, self.parent, self.children, self.is_extant = time, rate, parent, [], False


def grow_clads(rng, T, target=(46, 72)):
    """Forward ClaDS: crown start, per-lineage λ with lognormal jumps at each split."""
    while True:
        root = Node(0.0, LAM0)
        a, b = Node(0.0, LAM0, root), Node(0.0, LAM0, root)
        root.children = [a, b]
        live, t = [a, b], 0.0
        while live:
            tot = [n.rate * (1.0 + TURN) for n in live]
            total = math.fsum(tot)
            t += rng.exponential(1.0 / total)
            if t >= T:
                for n in live:
                    n.time, n.is_extant = T, True
                break
            i = int(rng.choice(len(live), p=np.array(tot) / total))
            n = live[i]
            n.time = t
            lam = n.rate
            if rng.random() < 1.0 / (1.0 + TURN):                 # speciation
                jf = lambda: lam * math.exp(rng.normal(math.log(ALPHA), SIGMA))
                c1, c2 = Node(t, jf(), n), Node(t, jf(), n)
                n.children = [c1, c2]
                live[i] = c1
                live.append(c2)
            else:                                                 # extinction
                live.pop(i)
            if len(live) > 500:
                break
        nl = _count_leaves(root)
        if target[0] <= nl <= target[1]:
            return root


def _count_leaves(n):
    return 1 if not n.children else sum(_count_leaves(c) for c in n.children)


def nodes(root):
    out, stack = [], [root]
    while stack:
        n = stack.pop()
        out.append(n)
        stack.extend(n.children)
    return out


def layout(root):
    order, ys = [0], {}

    def rec(n):
        if not n.children:
            ys[n] = order[0]; order[0] += 1
        else:
            for c in sorted(n.children, key=_count_leaves):
                rec(c)
            ys[n] = sum(ys[c] for c in n.children) / len(n.children)
    rec(root)
    return ys, order[0]


def render():
    rng = np.random.default_rng(45)
    T = 6.0
    root = grow_clads(rng, T)
    ys, nleaf = layout(root)
    alln = nodes(root)
    present = max(n.time for n in alln)
    logr = [math.log(n.rate) for n in alln if n.parent is not None]
    lo, hi = min(logr), max(logr)

    dy = min(11.5, TREE_H / max(1, nleaf - 1))
    top = TREE_TOP + (TREE_H - dy * (nleaf - 1)) / 2

    def X(t):
        return XL + t / present * (XR - XL)

    def Y(n):
        return top + ys[n] * dy

    def rate_style(rate):
        f = (math.log(rate) - lo) / (hi - lo) if hi > lo else 0.5
        g = int(round(200 - 182 * f))                            # light (slow) -> dark (fast)
        return "#%02x%02x%02x" % (g, g, g), 1.3 + 3.4 * f

    d = draw.Drawing(W, H, origin=(0, 0))
    d.append(draw.Rectangle(0, 0, W, H, fill="white"))
    d.append(draw.Text("Per-lineage rates", FS_TITLE, W / 2, 48,
                       font_family=FONT, text_anchor="middle", font_weight="bold", fill=INK))

    for n in alln:
        col, wid = rate_style(n.rate)
        px = X(0) - 12 if n.parent is None else X(n.parent.time)
        if n.parent is not None:
            d.append(draw.Line(px, Y(n), X(n.time), Y(n), stroke=col, stroke_width=wid, stroke_linecap="round"))
        for c in n.children:
            cc, cw = rate_style(c.rate)
            d.append(draw.Line(X(n.time), Y(n), X(n.time), Y(c), stroke=cc, stroke_width=cw))

    # rate colour/width bar -- top-left, close to the tree, never over the branches
    bw, bh = 250, 18
    bx, by = XL + 20, 128
    for k in range(60):
        f = k / 59
        g = int(round(200 - 182 * f))
        d.append(draw.Rectangle(bx + f * bw, by, bw / 60 + 0.6, bh, fill="#%02x%02x%02x" % (g, g, g)))
    d.append(draw.Rectangle(bx, by, bw, bh, fill="none", stroke=INK, stroke_width=0.9))
    d.append(draw.Text("speciation rate", FS_LABEL, bx, by - 12, font_family=FONT,
                       text_anchor="start", font_weight="bold", fill=INK))
    d.append(draw.Text("slow", FS_TICK, bx - 10, by + bh / 2, font_family=FONT, text_anchor="end",
                       dominant_baseline="central", fill="#555"))
    d.append(draw.Text("fast", FS_TICK, bx + bw + 10, by + bh / 2, font_family=FONT, text_anchor="start",
                       dominant_baseline="central", fill="#555"))

    # time axis -- baseline + ticks + labels, matching fig_diversity_dependent
    ya = top + dy * (nleaf - 1) + 28
    d.append(draw.Line(XL, ya, XR, ya, stroke="#bdbdbd", stroke_width=1.2))
    for i in range(6):
        t = present * i / 5
        d.append(draw.Line(X(t), ya, X(t), ya + 6, stroke=INK, stroke_width=1.2))
        d.append(draw.Text(f"{t:.1f}", FS_TICK, X(t), ya + 24, font_family=FONT, text_anchor="middle", fill="#777"))
    d.append(draw.Text("time (root to present)", FS_LABEL, (XL + XR) / 2, ya + 52, font_family=FONT,
                       text_anchor="middle", fill="#555"))

    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "clads.svg").write_text(d.as_svg(), encoding="utf-8")
    cairosvg.svg2png(bytestring=d.as_svg().encode(), write_to=str(OUT / "clads.png"), scale=300 / 72.0)
    print(f"wrote clads  ({nleaf} tips, log-rate range {lo:.2f}..{hi:.2f})")


if __name__ == "__main__":
    render()
