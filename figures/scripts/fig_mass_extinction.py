"""Figure: the mass-extinction model — a tree-wide survival pulse.

A ZOMBI2 species model may carry ``mass_extinctions=[(age, fraction), ...]``: at each age
before the present, every standing lineage independently dies with probability ``fraction``
— an instantaneous cataclysm (forward simulation only). This figure draws the COMPLETE tree
(surviving skeleton solid, extinct lineages dashed) with the pulse as a vertical wall where
a whole cohort of lineages terminates at once, and an aligned lineages-through-time curve
that shows the diversity crash and the recovery after it.

Monochrome (species-tree house style).  Run:  python figures/scripts/fig_mass_extinction.py
"""

from __future__ import annotations

from pathlib import Path

import cairosvg
import drawsvg as draw

from zombi2 import BirthDeath, simulate_species_tree

from zombi_style import FONT, INK, FS_TITLE, FS_LABEL, FS_ANNOT, FS_TICK

OUT = Path(__file__).resolve().parent.parent / "mass_extinction"

W, H = 1180, 780
XL, XR = 100, 1000
TREE_TOP, TREE_H = 150, 372
LTT_TOP, LTT_H = 572, 128
DASH = "6,5"
BAND = "#e9e9e9"
GREY = "#9a9a9a"
DATA = "#2b6cb0"          # one accent colour: the lineages-through-time curve

BIRTH, DEATH = 1.0, 0.15
AGE, PULSE_AGE, FRAC, SEED = 5.0, 2.5, 0.8, 9


def subtree_leaves(n):
    return 1 if n.is_leaf() else sum(subtree_leaves(c) for c in n.children)


def layout(tree):
    """Assign each node a y (leaf order, ladderised) and ``surv`` (has an extant descendant)."""
    order = [0]
    ys, surv = {}, {}

    def rec(n):
        if n.is_leaf():
            ys[n] = order[0]
            order[0] += 1
            surv[n] = n.is_extant
        else:
            kids = sorted(n.children, key=subtree_leaves)
            for c in kids:
                rec(c)
            ys[n] = sum(ys[c] for c in kids) / len(kids)
            surv[n] = any(surv[c] for c in kids)
    rec(tree.root)
    return ys, surv, order[0]


def render():
    tree = simulate_species_tree(BirthDeath(BIRTH, DEATH, mass_extinctions=[(PULSE_AGE, FRAC)]),
                                 age=AGE, direction="forward", age_type="crown", seed=SEED)
    present = tree.total_age
    pulse_t = present - PULSE_AGE
    ys, surv, nleaf = layout(tree)
    dy = min(13.0, TREE_H / max(1, nleaf - 1))
    top = TREE_TOP + (TREE_H - dy * (nleaf - 1)) / 2

    def X(t):
        return XL + t / present * (XR - XL)

    def Y(n):
        return top + ys[n] * dy

    d = draw.Drawing(W, H, origin=(0, 0))
    d.append(draw.Rectangle(0, 0, W, H, fill="white"))
    d.append(draw.Text("The mass-extinction model", FS_TITLE, W / 2, 48,
                       font_family=FONT, text_anchor="middle", font_weight="bold", fill=INK))

    # legend, lifted clear of the tree (a row under the title, left of centre)
    lx = 40
    ly = 96
    d.append(draw.Line(lx, ly, lx + 34, ly, stroke=INK, stroke_width=2.2, stroke_linecap="round"))
    d.append(draw.Text("surviving lineage", FS_LABEL, lx + 44, ly, font_family=FONT,
                       text_anchor="start", dominant_baseline="central", fill=INK))
    d.append(draw.Line(lx + 250, ly, lx + 284, ly, stroke=GREY, stroke_width=1.6, stroke_dasharray=DASH))
    d.append(draw.Text("extinct lineage", FS_LABEL, lx + 294, ly, font_family=FONT,
                       text_anchor="start", dominant_baseline="central", fill=INK))

    # present reference + pulse band (spans tree and LTT)
    d.append(draw.Line(X(present), TREE_TOP - 16, X(present), LTT_TOP + LTT_H,
                       stroke="#cccccc", stroke_width=1.0, stroke_dasharray="2,4"))
    px = X(pulse_t)
    d.append(draw.Rectangle(px - 7, TREE_TOP - 16, 14, (LTT_TOP + LTT_H) - (TREE_TOP - 16), fill=BAND))
    d.append(draw.Line(px, TREE_TOP - 16, px, LTT_TOP + LTT_H, stroke=INK, stroke_width=1.6, stroke_dasharray="5,4"))

    # tree skeleton
    def seg(x1, y1, x2, y2, dashed):
        kw = dict(stroke=INK, stroke_width=2.0)
        if dashed:
            kw.update(stroke=GREY, stroke_width=1.5, stroke_dasharray=DASH, stroke_linecap="butt")
        else:
            kw.update(stroke_linecap="round")
        d.append(draw.Line(x1, y1, x2, y2, **kw))

    victims = []
    for n in tree.nodes():
        px0 = X(0) - 14 if n.parent is None else X(n.parent.time)
        seg(px0, Y(n), X(n.time), Y(n), not surv[n])
        for c in n.children:
            seg(X(n.time), Y(n), X(n.time), Y(c), not surv[c])
        if n.is_leaf() and not n.is_extant and abs(n.time - pulse_t) < 1e-6:
            victims.append(n)
    for n in victims:                                                # cohort killed by the pulse
        d.append(draw.Circle(X(n.time), Y(n), 3.2, fill=INK))

    n_ext = sum(1 for n in tree.leaves() if n.is_extant)

    # --- lineages-through-time (aligned) ---
    grid = [present * i / 700 for i in range(701)]

    def alive(t):
        return sum(1 for n in tree.nodes() if n.parent is not None and n.parent.time < t <= n.time)
    counts = [alive(t) for t in grid]
    cmax = max(counts)

    def LY(c):
        return LTT_TOP + LTT_H - c / cmax * LTT_H
    d.append(draw.Line(XL, LTT_TOP + LTT_H, XR, LTT_TOP + LTT_H, stroke="#bdbdbd", stroke_width=1.2))
    d.append(draw.Line(XL, LTT_TOP, XL, LTT_TOP + LTT_H, stroke="#bdbdbd", stroke_width=1.2))
    pts = []
    for t, c in zip(grid, counts):
        pts += [X(t), LY(c)]
    d.append(draw.Lines(*pts, close=False, fill="none", stroke=INK, stroke_width=2.8,
                        stroke_linejoin="round"))
    for c in (0, cmax):
        d.append(draw.Text(str(c), FS_TICK, XL - 10, LY(c), font_family=FONT, text_anchor="end",
                           dominant_baseline="central", fill="#555"))
    d.append(draw.Text("lineages", FS_LABEL, XL - 40, LTT_TOP + LTT_H / 2, font_family=FONT,
                       text_anchor="middle", fill="#555",
                       transform=f"rotate(-90 {XL - 40} {LTT_TOP + LTT_H / 2})"))
    d.append(draw.Text("crash, then recovery", FS_ANNOT,
                       X(pulse_t) + 16, LTT_TOP + 22, font_family=FONT,
                       text_anchor="start", fill="#555", font_style="italic"))

    # time axis
    ya = LTT_TOP + LTT_H
    for i in range(6):
        t = present * i / 5
        d.append(draw.Line(X(t), ya, X(t), ya + 6, stroke=INK, stroke_width=1.2))
        d.append(draw.Text(f"{t:.0f}" if t == int(t) else f"{t:.1f}", FS_TICK, X(t), ya + 24,
                           font_family=FONT, text_anchor="middle", fill="#555"))
    d.append(draw.Text("time (root to present)", FS_LABEL, (XL + XR) / 2, ya + 52, font_family=FONT,
                       text_anchor="middle", fill="#555"))

    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "mass_extinction.svg").write_text(d.as_svg(), encoding="utf-8")
    cairosvg.svg2png(bytestring=d.as_svg().encode(), write_to=str(OUT / "mass_extinction.png"), scale=300 / 72.0)
    print(f"wrote mass_extinction  ({nleaf} leaves, {len(victims)} pulse victims, {n_ext} survivors)")


if __name__ == "__main__":
    render()
