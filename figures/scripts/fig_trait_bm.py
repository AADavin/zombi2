"""Trait figure: a continuous trait evolving by Brownian motion on the tree.

A single continuous character drifts up the time-calibrated species tree under
Brownian motion (variance rate sigma^2). Each branch is coloured as a gradient
from its parent's value to its own, so the whole tree is 'painted' by the trait;
a strip of tip values and a colour bar complete the read.

Uses the trait build in the sibling worktree (../ZOMBI2-traits, branch `traits`).
Run:  python figures/scripts/fig_trait_bm.py
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, "/Users/aadria/Desktop/CLAUDE/ZOMBI2-traits")   # traits-enabled zombi2

import drawsvg as draw

import phylustrator as ph
from zombi2 import BirthDeath, BrownianMotion, simulate_species_tree, simulate_traits

from model_common import zombi_to_ete3
from zombi_style import INK, PANEL, species_style

OUT_STEM = Path("/Users/aadria/Desktop/CLAUDE/ZOMBI2/figures/trait_bm/trait_bm")

N_TIPS, AGE, TREE_SEED = 14, 1.0, 3
SIGMA2, TRAIT_SEED = 2.0, 7
BRANCH_W = 5.0

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


def color_bar(d, x, y, w=150, h=16, vmin=0.0, vmax=1.0, title="Trait value"):
    grad = draw.LinearGradient(x, y, x + w, y)
    for t, c in VIRIDIS:
        grad.add_stop(t, hexc(c))
    d.drawing.append(grad)
    fam = d.style.font_family
    d.drawing.append(draw.Text(title, 16, x, y - 12, font_weight="bold", font_family=fam,
                               text_anchor="start", fill=INK))
    d.drawing.append(draw.Rectangle(x, y, w, h, fill=grad, stroke=INK, stroke_width=0.6))
    for tx, val, anchor in ((x, vmin, "start"), (x + w, vmax, "end")):
        d.drawing.append(draw.Text(f"{val:+.2f}", 13, tx, y + h + 15, font_family=fam,
                                   text_anchor=anchor, fill=INK))


def main():
    ztree = simulate_species_tree(BirthDeath(1.0, 0.3), n_tips=N_TIPS, age=AGE,
                                  direction="backward", seed=TREE_SEED)
    res = simulate_traits(ztree, BrownianMotion(sigma2=SIGMA2), seed=TRAIT_SEED)
    # node_values covers internal (ancestral) states too — values is leaves only,
    # which would leave every internal branch uncoloured.
    name2val = {n.name: float(v) for n, v in res.node_values.items()}
    vmin, vmax = min(name2val.values()), max(name2val.values())
    norm = lambda v: (v - vmin) / (vmax - vmin) if vmax > vmin else 0.5   # noqa: E731
    node_to_rgb = {name: viridis(norm(v)) for name, v in name2val.items()}

    tree = zombi_to_ete3(ztree)
    style = species_style(height=max(680, 44 * N_TIPS + 180), font_size=14)
    d = ph.VerticalTreeDrawer(tree, style=style)
    d._calculate_layout()
    d.plot_continuous_variable(node_to_rgb, stroke_width=BRANCH_W)

    # tip value strip: a colour chip per tip, just past the (aligned) tips
    tip_x = max(l.coordinates[0] for l in tree.get_leaves()) + 22
    for lf in tree.get_leaves():
        _, y = lf.coordinates
        d._draw_shape_at(tip_x, y, "square", hexc(node_to_rgb[lf.name]), r=8,
                         stroke=INK, stroke_width=1.0)

    d.add_time_axis(ticks=[round(AGE * i / 4, 6) for i in range(5)],
                    tick_labels=[f"{AGE * i / 4:.2f}" for i in range(5)],
                    label="Time (root to present)", tick_size=6.0, padding=14.0, stroke_width=1.6)

    color_bar(d, x=-style.width / 2 + 34, y=-style.height / 2 + 44,
              vmin=vmin, vmax=vmax, title="Trait value (Brownian motion)")

    d.save_svg(f"{OUT_STEM}.svg")
    d.save_png(f"{OUT_STEM}.png", dpi=300)
    print(f"wrote {OUT_STEM}.svg / .png  ({N_TIPS} tips, trait range [{vmin:+.2f}, {vmax:+.2f}])")


if __name__ == "__main__":
    main()
