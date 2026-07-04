"""Figure: the diversity-dependent birth-death model.

``DiversityDependent(lambda_0, death, carrying_capacity=K)``: the speciation rate declines
linearly with the number of standing lineages, ``λ(n) = λ0·(1 − n/K)``, so the tree grows
fast when small and saturates as it fills its carrying capacity K. This figure draws the
complete tree (survivors solid, extinct dashed) and the aligned lineages-through-time curve,
which climbs and then plateaus near the equilibrium — the density-dependent signature.

Monochrome (species-tree house style).  Run:  python figures/scripts/fig_diversity_dependent.py
"""

from __future__ import annotations

from pathlib import Path

import cairosvg
import drawsvg as draw

from zombi2 import DiversityDependent, simulate_species_tree

from zombi_style import FONT, INK

OUT = Path(__file__).resolve().parent.parent / "diversity_dependent"

W, H = 1180, 760
XL, XR = 100, 1000
TREE_TOP, TREE_H = 104, 388
LTT_TOP, LTT_H = 566, 150
DASH = "6,5"
GREY = "#9a9a9a"

LAM0, DEATH, K, AGE, SEED = 1.5, 0.12, 24, 7.0, 3


def subtree_leaves(n):
    return 1 if n.is_leaf() else sum(subtree_leaves(c) for c in n.children)


def layout(tree):
    order = [0]
    ys, surv = {}, {}

    def rec(n):
        if n.is_leaf():
            ys[n] = order[0]; order[0] += 1; surv[n] = n.is_extant
        else:
            kids = sorted(n.children, key=subtree_leaves)
            for c in kids:
                rec(c)
            ys[n] = sum(ys[c] for c in kids) / len(kids)
            surv[n] = any(surv[c] for c in kids)
    rec(tree.root)
    return ys, surv, order[0]


def render():
    tree = simulate_species_tree(DiversityDependent(LAM0, DEATH, carrying_capacity=K),
                                 age=AGE, direction="forward", age_type="crown", seed=SEED)
    present = tree.total_age
    ys, surv, nleaf = layout(tree)
    dy = min(12.5, TREE_H / max(1, nleaf - 1))
    top = TREE_TOP + (TREE_H - dy * (nleaf - 1)) / 2

    def X(t):
        return XL + t / present * (XR - XL)

    def Y(n):
        return top + ys[n] * dy

    d = draw.Drawing(W, H, origin=(0, 0))
    d.append(draw.Rectangle(0, 0, W, H, fill="white"))
    d.append(draw.Text("The diversity-dependent model — diversification slows as the tree fills up", 20,
                       40, 40, font_family=FONT, text_anchor="start", font_weight="bold", fill=INK))
    d.append(draw.Text("speciation rate falls with standing diversity: λ(n) = λ0·(1 − n/K).  ZOMBI2: "
                       "DiversityDependent(λ0, death, carrying_capacity=K)", 13, 40, 64,
                       font_family=FONT, text_anchor="start", fill="#777"))

    d.append(draw.Line(X(present), TREE_TOP - 12, X(present), LTT_TOP + LTT_H,
                       stroke="#cccccc", stroke_width=1.0, stroke_dasharray="2,4"))

    def seg(x1, y1, x2, y2, dashed):
        kw = dict(stroke=INK, stroke_width=2.0, stroke_linecap="round")
        if dashed:
            kw = dict(stroke=GREY, stroke_width=1.5, stroke_dasharray=DASH, stroke_linecap="butt")
        d.append(draw.Line(x1, y1, x2, y2, **kw))

    for n in tree.nodes():
        px0 = X(0) - 14 if n.parent is None else X(n.parent.time)
        seg(px0, Y(n), X(n.time), Y(n), not surv[n])
        for c in n.children:
            seg(X(n.time), Y(n), X(n.time), Y(c), not surv[c])

    n_ext = sum(1 for n in tree.leaves() if n.is_extant)

    # --- LTT ---
    grid = [present * i / 720 for i in range(721)]

    def alive(t):
        return sum(1 for n in tree.nodes() if n.parent is not None and n.parent.time < t <= n.time)
    counts = [alive(t) for t in grid]
    cmax = max(K, max(counts))

    def LY(c):
        return LTT_TOP + LTT_H - c / cmax * LTT_H
    d.append(draw.Line(XL, LTT_TOP + LTT_H, XR, LTT_TOP + LTT_H, stroke="#bdbdbd", stroke_width=1.2))
    d.append(draw.Line(XL, LTT_TOP, XL, LTT_TOP + LTT_H, stroke="#bdbdbd", stroke_width=1.2))
    # carrying capacity K
    d.append(draw.Line(XL, LY(K), XR, LY(K), stroke=INK, stroke_width=1.4, stroke_dasharray="7,4"))
    d.append(draw.Text("carrying capacity K", 12, XR - 4, LY(K) - 7, font_family=FONT,
                       text_anchor="end", font_weight="bold", fill=INK))
    pts = []
    for t, c in zip(grid, counts):
        pts += [X(t), LY(c)]
    d.append(draw.Lines(*pts, close=False, fill="none", stroke=INK, stroke_width=2.0, stroke_linejoin="round"))
    for c in (0, K):
        d.append(draw.Text(str(c), 11, XL - 8, LY(c), font_family=FONT, text_anchor="end",
                           dominant_baseline="central", fill="#777"))
    d.append(draw.Text("lineages", 12.5, XL - 32, LTT_TOP + LTT_H / 2, font_family=FONT,
                       text_anchor="middle", fill="#777",
                       transform=f"rotate(-90 {XL - 32} {LTT_TOP + LTT_H / 2})"))
    d.append(draw.Text("fast growth while small, then a plateau near K", 12.5, X(present * 0.55),
                       LTT_TOP + LTT_H - 16, font_family=FONT, text_anchor="middle", fill="#555",
                       font_style="italic"))

    ya = LTT_TOP + LTT_H
    for i in range(6):
        t = present * i / 5
        d.append(draw.Line(X(t), ya, X(t), ya + 5, stroke=INK, stroke_width=1.2))
        d.append(draw.Text(f"{t:.0f}" if t == int(t) else f"{t:.1f}", 11, X(t), ya + 18,
                           font_family=FONT, text_anchor="middle", fill="#777"))
    d.append(draw.Text("time (root -> present)", 12.5, (XL + XR) / 2, ya + 38, font_family=FONT,
                       text_anchor="middle", fill="#777"))

    lx, ly = XR - 250, TREE_TOP + 4
    d.append(draw.Line(lx, ly, lx + 26, ly, stroke=INK, stroke_width=2.0, stroke_linecap="round"))
    d.append(draw.Text("surviving lineage", 12, lx + 32, ly, font_family=FONT, text_anchor="start",
                       dominant_baseline="central", fill=INK))
    d.append(draw.Line(lx, ly + 20, lx + 26, ly + 20, stroke=GREY, stroke_width=1.5, stroke_dasharray=DASH))
    d.append(draw.Text("extinct lineage", 12, lx + 32, ly + 20, font_family=FONT, text_anchor="start",
                       dominant_baseline="central", fill=INK))

    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "diversity_dependent.svg").write_text(d.as_svg(), encoding="utf-8")
    cairosvg.svg2png(bytestring=d.as_svg().encode(), write_to=str(OUT / "diversity_dependent.png"), scale=300 / 72.0)
    print(f"wrote diversity_dependent  ({nleaf} leaves, {n_ext} extant, K={K})")


if __name__ == "__main__":
    render()
