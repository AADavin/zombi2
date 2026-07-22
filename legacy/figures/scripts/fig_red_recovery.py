"""Figure (Tools/RED): RED recovers a timetree from a rate-distorted phylogram.

One canonical tree, shown three ways, same tip order in every panel:

  A. the TRUE timetree (a chronogram; ultrametric, tips aligned at the present),
  B. the same tree PERTURBED into a phylogram by a strong relaxed clock (branch
     lengths are substitutions; rate variation makes the tips ragged),
  C. the tree RECOVERED by Relative Evolutionary Divergence — each node placed at
     its RED value, which pulls the tips back into line and lands on panel A.

Each panel's x-axis is relative depth (root 0 -> present 1), so A and C compare
directly. House style: near-black/red/green trees, one centred title, ASCII text.

Run:  python figures/scripts/fig_red_recovery.py
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import cairosvg
import drawsvg as draw

import clock_common as C
from zombi_style import FONT, INK, MUTED, FS_TITLE, FS_LABEL, FS_TICK, AVOID, COOCCUR  # noqa: F401

OUT_DIR = Path(__file__).resolve().parent.parent

TREE = dict(n_tips=16, age=1.0, seed=3)
SIGMA = 1.3             # autocorrelated-lognormal clock strength (strong, for a vivid B)
CLOCK_SEED = 9

# panels: (key, title, subtitle, colour)
TRUE_COL = "#2f6f8f"
PANELS = [
    ("true", "A  ·  True timetree", "simulated · ultrametric (tips aligned)", TRUE_COL),
    ("pert", "B  ·  Perturbed", "strong rate variation · tips ragged", AVOID),
    ("red",  "C  ·  Recovered by RED", "node depth = RED · ultrametric again", COOCCUR),
]

W, H = 1180, 560
MARGIN_X = 30
TOP = 92
ROW_AXIS = 84
PANEL_W = (W - 2 * MARGIN_X) / 3
PAD_L, PAD_R, PAD_T = 16, 22, 64
INNER_W = PANEL_W - PAD_L - PAD_R
INNER_H = H - TOP - PAD_T - ROW_AXIS
BRANCH_W = 3.0
AXIS_MAX = 1.0                     # every panel is relative depth in [0, 1]


def red_values(tree, subst_len: dict) -> dict:
    """RED (Parks et al. 2018) on the phylogram, keyed by node name. Root 0, leaves 1.
    (The library tool is ``zombi2.tools.relative_evolutionary_divergence``; recomputed
    here on the figure's ete3 tree, as the sibling clock figures do for their clocks.)"""
    mtd, ntip = {}, {}
    for n in tree.traverse("postorder"):
        if n.is_leaf():
            mtd[n.name], ntip[n.name] = 0.0, 1
        else:
            tot = sum(ntip[c.name] * (subst_len[c.name] + mtd[c.name]) for c in n.children)
            k = sum(ntip[c.name] for c in n.children)
            mtd[n.name], ntip[n.name] = tot / k, k
    red = {}
    for n in tree.traverse("preorder"):
        if n.is_root():
            red[n.name] = 0.0
        else:
            a, b, rp = subst_len[n.name], mtd[n.name], red[n.up.name]
            red[n.name] = rp + (a / (a + b)) * (1.0 - rp) if (a + b) > 0 else rp
    return red


def main():
    tree = C.build_tree(**TREE)
    tfo, present = C.node_times(tree)
    ys, nleaf = C.leaf_ys(tree)

    # perturb: autocorrelated-lognormal clock -> substitution branch lengths
    _rate, subst = C.autocorrelated_lognormal(tree, sigma=SIGMA, seed=CLOCK_SEED)
    dr = C.subst_dist_to_root(tree, subst)
    smax = max(dr[l.name] for l in tree.get_leaves())
    red = red_values(tree, subst)

    # relative depth (root 0 -> present 1) per panel
    depth = {
        "true": {name: tfo[name] / present for name in tfo},
        "pert": {name: dr[name] / smax for name in dr},
        "red":  red,
    }

    d = draw.Drawing(W, H, origin=(0, 0))
    d.append(draw.Rectangle(0, 0, W, H, fill="white"))
    d.append(draw.Text("RED recovers the timetree from a rate-distorted phylogram", FS_TITLE,
                       W / 2, 46, font_family=FONT, text_anchor="middle", font_weight="bold", fill=INK))

    for idx, (key, title, sub, col) in enumerate(PANELS):
        ox = MARGIN_X + idx * PANEL_W
        _panel(d, ox, tree, ys, nleaf, depth[key], col, title, sub,
               aligned=(key in ("true", "red")))

    name = "red_recovery"
    out = OUT_DIR / name
    out.mkdir(parents=True, exist_ok=True)
    (out / f"{name}.svg").write_text(d.as_svg(), encoding="utf-8")
    cairosvg.svg2png(bytestring=d.as_svg().encode(), write_to=str(out / f"{name}.png"), scale=300 / 72.0)

    # a couple of numbers for the caption
    import numpy as np
    rag = smax / min(dr[l.name] for l in tree.get_leaves())
    internal = [n.name for n in tree.traverse() if not n.is_leaf()]
    corr = np.corrcoef([depth["true"][k] for k in internal], [red[k] for k in internal])[0, 1]
    print(f"wrote {out}/{name}.svg / .png  (B raggedness {rag:.1f}x, A-vs-C corr {corr:.4f})")


def _panel(d, ox, tree, ys, nleaf, depth, col, title, subtitle, aligned):
    oy = TOP
    x_at = lambda v: ox + PAD_L + v * (INNER_W / AXIS_MAX)                 # noqa: E731
    y_at = lambda k: oy + PAD_T + (k / max(1, nleaf - 1)) * INNER_H        # noqa: E731

    d.append(draw.Text(title, FS_LABEL, ox + PANEL_W / 2, oy + 20, font_family=FONT,
                       text_anchor="middle", font_weight="bold", fill=INK))
    d.append(draw.Text(subtitle, FS_TICK, ox + PANEL_W / 2, oy + 42, font_family=FONT,
                       text_anchor="middle", fill=MUTED, font_style="italic"))

    # a light 'present' reference at depth 1 (tips reach it when ultrametric)
    d.append(draw.Line(x_at(1.0), y_at(0) - 7, x_at(1.0), y_at(nleaf - 1) + 7,
                       stroke=MUTED, stroke_width=1.0, stroke_dasharray="3,3"))

    for n in tree.traverse():                                  # horizontal branches
        if n.is_root():
            continue
        y = y_at(ys[n.name])
        d.append(draw.Line(x_at(depth[n.up.name]), y, x_at(depth[n.name]), y,
                           stroke=col, stroke_width=BRANCH_W, stroke_linecap="round"))
    for n in tree.traverse("postorder"):                       # vertical connectors
        if not n.is_leaf():
            x = x_at(depth[n.name])
            yy = [y_at(ys[ch.name]) for ch in n.children]
            d.append(draw.Line(x, min(yy), x, max(yy), stroke=col, stroke_width=BRANCH_W,
                               stroke_linecap="round"))
    for lf in tree.get_leaves():                               # tip dots
        d.append(draw.Circle(x_at(depth[lf.name]), y_at(ys[lf.name]), 3.0, fill=col,
                             stroke="white", stroke_width=0.7))

    base = oy + PAD_T + INNER_H + 30                           # axis
    d.append(draw.Line(x_at(0.0), base, x_at(1.0), base, stroke=INK, stroke_width=1.4))
    for k in range(3):
        v = 0.5 * k
        d.append(draw.Line(x_at(v), base, x_at(v), base + 5, stroke=INK, stroke_width=1.4))
        d.append(draw.Text(f"{v:g}", FS_TICK, x_at(v), base + 20, font_family=FONT,
                           text_anchor="middle", fill=MUTED))
    d.append(draw.Text("relative depth  (root to present)", FS_TICK, (x_at(0) + x_at(1)) / 2,
                       base + 40, font_family=FONT, text_anchor="middle", fill=MUTED,
                       font_style="italic"))


if __name__ == "__main__":
    main()
