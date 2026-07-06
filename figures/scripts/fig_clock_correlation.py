"""Figure (Ch16): uncorrelated vs autocorrelated clocks -- the central split.

The one distinction that organises the whole family: does a branch's rate tell you
anything about its neighbours'? Two painted trees (the same tree, branches coloured by
rate) sit above the *proof*: a scatter of every parent branch's rate against its child
branch's rate.

  * Left (uncorrelated lognormal): each branch draws its rate independently. The colours
    are salt-and-pepper and the scatter is a shapeless cloud -- knowing a branch's rate
    tells you nothing about its child's (correlation near 0).
  * Right (autocorrelated lognormal): a branch inherits its parent's rate and drifts from
    there. The colour changes smoothly down the tree and the scatter hugs the diagonal --
    child rate tracks parent rate (strong positive correlation).

House style: painted trees, viridis rate scale shared with the family figure, ASCII text.

Run:  /Users/aadria/miniconda3/bin/python figures/scripts/fig_clock_correlation.py
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import cairosvg
import drawsvg as draw

import clock_common as C
from zombi_style import FONT, INK, MUTED, FS_TITLE, FS_LABEL, FS_ANNOT, FS_TICK

OUT_DIR = Path(__file__).resolve().parent.parent

W, H = 1300, 940
N_TIPS, TREE_SEED = 30, 6
SIGMA = 0.55
BRANCH_W = 2.8


def pearson_logr(pairs):
    xs = [math.log(a) for a, b in pairs]
    ys = [math.log(b) for a, b in pairs]
    n = len(xs)
    mx, my = sum(xs) / n, sum(ys) / n
    sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    sxx = sum((x - mx) ** 2 for x in xs)
    syy = sum((y - my) ** 2 for y in ys)
    return sxy / math.sqrt(sxx * syy) if sxx > 0 and syy > 0 else 0.0


def parent_child_pairs(tree, rate):
    """(parent-branch rate, child-branch rate) for every branch whose parent is
    itself a real branch (i.e. the parent is not the root)."""
    out = []
    for n in tree.traverse():
        if n.is_root() or n.up.is_root():
            continue
        out.append((rate[n.up.name], rate[n.name]))
    return out


def draw_tree(d, ox, oy, pw, ph, tree, rate, sublen, title, tag):
    ys, nleaf = C.leaf_ys(tree)
    dist = C.subst_dist_to_root(tree, sublen)
    maxd = max(dist[l.name] for l in tree.get_leaves())
    x_at = lambda v: ox + 6 + (v / maxd) * (pw - 12)          # noqa: E731
    y_at = lambda k: oy + 40 + (k / max(1, nleaf - 1)) * (ph - 56)  # noqa: E731

    d.append(draw.Text(title, FS_LABEL, ox + pw / 2, oy + 2, font_family=FONT,
                       text_anchor="middle", font_weight="bold", fill=INK))
    d.append(draw.Text(tag, FS_TICK, ox + pw / 2, oy + 24, font_family=FONT,
                       text_anchor="middle", fill=MUTED, font_style="italic"))
    for n in tree.traverse():
        if n.is_root():
            continue
        y = y_at(ys[n.name])
        d.append(draw.Line(x_at(dist[n.up.name]), y, x_at(dist[n.name]), y,
                           stroke=C.rate_hex(rate[n.name]), stroke_width=BRANCH_W, stroke_linecap="round"))
    for n in tree.traverse("postorder"):
        if not n.is_leaf():
            x = x_at(dist[n.name])
            yy = [y_at(ys[c.name]) for c in n.children]
            d.append(draw.Line(x, min(yy), x, max(yy), stroke=INK, stroke_width=1.3))


def draw_scatter(d, ox, oy, s, pairs, r):
    """rate(parent) on x, rate(child) on y, on a shared log axis; y=x diagonal drawn."""
    lo, hi = C.RATE_LO, C.RATE_HI
    lx = lambda v: ox + (math.log(max(lo, min(hi, v))) - math.log(lo)) / (math.log(hi) - math.log(lo)) * s  # noqa: E731
    ly = lambda v: oy + s - (math.log(max(lo, min(hi, v))) - math.log(lo)) / (math.log(hi) - math.log(lo)) * s  # noqa: E731
    # frame
    d.append(draw.Rectangle(ox, oy, s, s, fill="white", stroke=MUTED, stroke_width=1.0))
    # y = x diagonal
    d.append(draw.Line(lx(lo), ly(lo), lx(hi), ly(hi), stroke=MUTED, stroke_width=1.2,
                       stroke_dasharray="4,4"))
    # 1x guides
    for gx in (1.0,):
        d.append(draw.Line(lx(gx), oy, lx(gx), oy + s, stroke="#e2e2e2", stroke_width=1.0))
        d.append(draw.Line(ox, ly(gx), ox + s, ly(gx), stroke="#e2e2e2", stroke_width=1.0))
    for a, b in pairs:
        d.append(draw.Circle(lx(a), ly(b), 4.2, fill=C.rate_hex((a * b) ** 0.5),
                             stroke=INK, stroke_width=0.6, fill_opacity=0.9))
    # axis labels
    d.append(draw.Text("parent branch rate", FS_TICK, ox + s / 2, oy + s + 42, font_family=FONT,
                       text_anchor="middle", fill=INK))
    d.append(draw.Text("child branch rate", FS_TICK, ox - 26, oy + s / 2, font_family=FONT,
                       text_anchor="middle", fill=INK, transform=f"rotate(-90,{ox - 26},{oy + s / 2})"))
    fmt = {lo: "1/3", 1.0: "1", hi: "3"}
    for v in (lo, 1.0, hi):
        anc = "start" if v == lo else ("end" if v == hi else "middle")
        d.append(draw.Text(fmt[v], FS_TICK, lx(v), oy + s + 20, font_family=FONT,
                           text_anchor=anc, fill=MUTED))
        d.append(draw.Text(fmt[v], FS_TICK, ox - 8, ly(v) + 4, font_family=FONT,
                           text_anchor="end", fill=MUTED))
    # correlation readout
    d.append(draw.Text(f"correlation r = {r:+.2f}", FS_ANNOT, ox + s / 2, oy - 12, font_family=FONT,
                       text_anchor="middle", font_weight="bold", fill=INK))


def main():
    tree = C.build_tree(n_tips=N_TIPS, seed=TREE_SEED)
    Ru, Su = C.uncorrelated_lognormal(tree, SIGMA, seed=3)
    Ra, Sa = C.autocorrelated_lognormal(tree, SIGMA, seed=4)

    d = draw.Drawing(W, H, origin=(0, 0))
    d.append(draw.Rectangle(0, 0, W, H, fill="white"))
    d.append(draw.Text("Uncorrelated vs autocorrelated clocks", FS_TITLE, W / 2, 46,
                       font_family=FONT, text_anchor="middle", font_weight="bold", fill=INK))
    _rate_bar(d, W / 2 - 175, 84, 350, 16)

    # trees (top), scatters (bottom), each half its column
    colw = (W - 120) / 2
    tx0, tx1 = 60, 60 + colw + 20
    draw_tree(d, tx0, 150, colw - 20, 420, tree, Ru, Su,
              "Uncorrelated lognormal", "each branch drawn independently")
    draw_tree(d, tx1, 150, colw - 20, 420, tree, Ra, Sa,
              "Autocorrelated lognormal", "each branch inherits its parent's rate")

    s = 250
    sc_y = 640
    draw_scatter(d, tx0 + (colw - 20 - s) / 2, sc_y, s, parent_child_pairs(tree, Ru),
                 pearson_logr(parent_child_pairs(tree, Ru)))
    draw_scatter(d, tx1 + (colw - 20 - s) / 2, sc_y, s, parent_child_pairs(tree, Ra),
                 pearson_logr(parent_child_pairs(tree, Ra)))

    name = "clock_correlation"
    out = OUT_DIR / name
    out.mkdir(parents=True, exist_ok=True)
    (out / f"{name}.svg").write_text(d.as_svg(), encoding="utf-8")
    cairosvg.svg2png(bytestring=d.as_svg().encode(), write_to=str(out / f"{name}.png"), scale=300 / 72.0)
    print(f"wrote {out}/{name}.svg / .png")


def _rate_bar(d, x, y, w, h):
    grad = draw.LinearGradient(x, y, x + w, y)
    for i in range(21):
        t = i / 20.0
        grad.add_stop(t, C.rate_hex(C.RATE_LO * (C.RATE_HI / C.RATE_LO) ** t))
    d.append(grad)
    d.append(draw.Text("branch rate multiplier", FS_TICK, x + w / 2, y - 8, font_family=FONT,
                       text_anchor="middle", font_weight="bold", fill=INK))
    d.append(draw.Rectangle(x, y, w, h, fill=grad, stroke=INK, stroke_width=0.8))
    for tx, lab, anc in ((x, "slow", "start"), (x + w / 2, "1x", "middle"), (x + w, "fast", "end")):
        d.append(draw.Text(lab, FS_TICK, tx, y + h + 16, font_family=FONT, text_anchor=anc, fill=MUTED))


if __name__ == "__main__":
    main()
