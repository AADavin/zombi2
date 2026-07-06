"""Trait figure: a continuous trait evolving by Brownian motion on the tree.

A single continuous character drifts down the time-calibrated species tree under Brownian motion
(variance rate sigma^2). Each branch is painted as a gradient from its parent's value to its own,
so the whole tree is coloured by the trait; a tip strip and a colour bar complete the read. Drawn
in the same painted-tree house style as the Ornstein-Uhlenbeck figure, so the two read as a pair.

This module also exports the shared viridis colormap (VIRIDIS / viridis / hexc) used by the other
continuous-trait figures.

Run:  /Users/aadria/miniconda3/bin/python figures/scripts/fig_trait_bm.py
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))   # local: zombi_style, model_common

import drawsvg as draw

import phylustrator as ph
from zombi2.species import BirthDeath, simulate_species_tree
from zombi2.traits import BrownianMotion, simulate_traits

from model_common import zombi_to_ete3
from zombi_style import INK, species_style, FS_TITLE, FS_LABEL, FS_TICK

OUT_STEM = Path(__file__).resolve().parent.parent / "trait_bm" / "trait_bm"

N_TIPS, AGE, TREE_SEED = 28, 1.0, 3
SIGMA2, X0, TRAIT_SEED = 1.0, 0.0, 12
BRANCH_W = 4.0

# viridis anchor stops (t, (r,g,b))
VIRIDIS = [(0.0, (68, 1, 84)), (0.25, (59, 82, 139)), (0.5, (33, 145, 140)),
           (0.75, (94, 201, 98)), (1.0, (253, 231, 37))]


def viridis(t):
    t = max(0.0, min(1.0, t))
    for (t0, c0), (t1, c1) in zip(VIRIDIS, VIRIDIS[1:]):
        if t <= t1:
            f = (t - t0) / (t1 - t0) if t1 > t0 else 0.0
            return tuple(int(round(c0[i] + (c1[i] - c0[i]) * f)) for i in range(3))
    return VIRIDIS[-1][1]


def hexc(rgb):
    return "#%02x%02x%02x" % tuple(rgb)


def color_bar(d, x, y, w, h, vmin, vmax, title):
    """House-style viridis colour bar with readable (FS_*) labels."""
    fam = d.style.font_family
    grad = draw.LinearGradient(x, y, x + w, y)
    for t, c in VIRIDIS:
        grad.add_stop(t, hexc(c))
    d.drawing.append(grad)
    d.drawing.append(draw.Text(title, FS_LABEL, x, y - 14, font_weight="bold", font_family=fam,
                               text_anchor="start", fill=INK))
    d.drawing.append(draw.Rectangle(x, y, w, h, fill=grad, stroke=INK, stroke_width=0.8))
    for tx, val, anchor in ((x, vmin, "start"), (x + w, vmax, "end")):
        d.drawing.append(draw.Text(f"{val:+.2f}", FS_TICK, tx, y + h + 22, font_family=fam,
                                   text_anchor=anchor, fill="#555"))


def pick_seed(ztree):
    """A trait seed whose tip values spread evenly (avoids one outlier washing out the gradient)."""
    import statistics
    best = None
    for s in range(1, 40):
        res = simulate_traits(ztree, BrownianMotion(sigma2=SIGMA2, x0=X0), seed=s)
        vals = [float(v) for n, v in res.node_values.items()]
        lo, hi = min(vals), max(vals)
        if hi <= lo:
            continue
        tips = [float(v) for n, v in res.node_values.items() if n.is_leaf()]
        # score = how uniformly the tips fill their range (std / range; higher = better spread)
        score = statistics.pstdev(tips) / (hi - lo)
        if best is None or score > best[0]:
            best = (score, s, res, lo, hi)
    return best[1], best[2], best[3], best[4]


def main():
    ztree = simulate_species_tree(BirthDeath(1.0, 0.3), n_tips=N_TIPS, age=AGE,
                                  direction="backward", seed=TREE_SEED)
    seed, res, vmin, vmax = pick_seed(ztree)
    name2val = {n.name: float(v) for n, v in res.node_values.items()}
    norm = lambda v: (v - vmin) / (vmax - vmin) if vmax > vmin else 0.5   # noqa: E731
    node_to_rgb = {name: viridis(norm(v)) for name, v in name2val.items()}

    tree = zombi_to_ete3(ztree)
    style = species_style(width=1180, height=740, margin=120, font_size=FS_TICK)
    d = ph.VerticalTreeDrawer(tree, style=style)
    d._calculate_layout()
    d.plot_continuous_variable(node_to_rgb, stroke_width=BRANCH_W)

    tip_x = max(l.coordinates[0] for l in tree.get_leaves()) + 22
    for lf in tree.get_leaves():
        _, y = lf.coordinates
        d._draw_shape_at(tip_x, y, "square", hexc(node_to_rgb[lf.name]), r=6, stroke=INK, stroke_width=1.0)

    d.add_time_axis(ticks=[round(AGE * i / 4, 6) for i in range(5)],
                    tick_labels=[f"{AGE * i / 4:.2f}" for i in range(5)],
                    label="Time (root to present)", tick_size=6.0, padding=16.0, stroke_width=1.6)

    d.drawing.append(draw.Text("A trait under Brownian motion", FS_TITLE, 0, -style.height / 2 + 42,
                               font_weight="bold", font_family=style.font_family,
                               text_anchor="middle", dominant_baseline="central", fill=INK))
    color_bar(d, x=-style.width / 2 + 32, y=-style.height / 2 + 96, w=220, h=18,
              vmin=vmin, vmax=vmax, title="Trait value")

    d.save_svg(f"{OUT_STEM}.svg")
    d.save_png(f"{OUT_STEM}.png", dpi=300)
    print(f"wrote {OUT_STEM}.svg / .png  ({N_TIPS} tips, trait seed {seed}, "
          f"range [{vmin:+.2f}, {vmax:+.2f}])")


if __name__ == "__main__":
    main()
