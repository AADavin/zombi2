"""Trait figure: a continuous trait under Ornstein-Uhlenbeck (pull to an optimum).

Same painted-tree style as the Brownian-motion figure, on the SAME tree, so the
two read as a pair: BM drifts freely and spreads out, whereas OU is pulled toward
an optimum theta (here from an ancestral value x0 = 0 up to theta = 2). The trait
therefore climbs from the root's value and then clusters near theta at the tips;
the optimum is marked on the colour bar.

Run:  python figures/scripts/fig_trait_ou.py
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, "/Users/aadria/Desktop/CLAUDE/ZOMBI2-traits")   # traits-enabled zombi2

import drawsvg as draw

import phylustrator as ph
from zombi2 import BirthDeath, OrnsteinUhlenbeck, simulate_species_tree, simulate_traits

from fig_trait_bm import VIRIDIS, hexc, viridis   # shared colormap
from model_common import zombi_to_ete3
from zombi_style import INK, species_style, FS_TITLE, FS_LABEL, FS_ANNOT, FS_TICK

OUT_STEM = Path("/Users/aadria/Desktop/CLAUDE/ZOMBI2/figures/trait_ou/trait_ou")


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

N_TIPS, AGE, TREE_SEED = 28, 1.0, 3          # ~2x the tips of fig_trait_bm
SIGMA2, ALPHA, THETA, X0 = 1.2, 4.0, 2.0, 0.0
TRAIT_SEED = 7
BRANCH_W = 4.0


def main():
    ztree = simulate_species_tree(BirthDeath(1.0, 0.3), n_tips=N_TIPS, age=AGE,
                                  direction="backward", seed=TREE_SEED)
    res = simulate_traits(ztree, OrnsteinUhlenbeck(sigma2=SIGMA2, alpha=ALPHA, theta=THETA, x0=X0),
                          seed=TRAIT_SEED)
    name2val = {n.name: float(v) for n, v in res.node_values.items()}
    vmin, vmax = min(name2val.values()), max(name2val.values())
    norm = lambda v: (v - vmin) / (vmax - vmin) if vmax > vmin else 0.5   # noqa: E731
    node_to_rgb = {name: viridis(norm(v)) for name, v in name2val.items()}

    tree = zombi_to_ete3(ztree)
    # wide, low canvas: generous top margin keeps the header band (title + colour bar)
    # clear of the tree, while the extra width and reduced height make the tree
    # broad and short rather than tall.
    style = species_style(width=1180, height=740, margin=120, font_size=FS_TICK)
    d = ph.VerticalTreeDrawer(tree, style=style)
    d._calculate_layout()
    d.plot_continuous_variable(node_to_rgb, stroke_width=BRANCH_W)

    tip_x = max(l.coordinates[0] for l in tree.get_leaves()) + 22
    for lf in tree.get_leaves():
        _, y = lf.coordinates
        d._draw_shape_at(tip_x, y, "square", hexc(node_to_rgb[lf.name]), r=6,
                         stroke=INK, stroke_width=1.0)

    d.add_time_axis(ticks=[round(AGE * i / 4, 6) for i in range(5)],
                    tick_labels=[f"{AGE * i / 4:.2f}" for i in range(5)],
                    label="Time (root to present)", tick_size=6.0, padding=16.0, stroke_width=1.6)

    # title centered at the top; colour bar (the "data") in the top-left -- both clear of the tree
    left = -style.width / 2 + 30
    d.drawing.append(draw.Text("A trait under Ornstein-Uhlenbeck", FS_TITLE, 0,
                               -style.height / 2 + 42, font_weight="bold",
                               font_family=style.font_family, text_anchor="middle",
                               dominant_baseline="central", fill=INK))
    bar_x, bar_y, bar_w, bar_h = left + 2, -style.height / 2 + 96, 220, 18
    color_bar(d, x=bar_x, y=bar_y, w=bar_w, h=bar_h, vmin=vmin, vmax=vmax,
              title="Trait value")
    # a mark on the bar where the optimum sits, plus a caption to its right
    tx = bar_x + max(0.0, min(1.0, (THETA - vmin) / (vmax - vmin))) * bar_w
    d.drawing.append(draw.Line(tx, bar_y - 3, tx, bar_y + bar_h + 3, stroke=INK, stroke_width=2.4))
    d.add_text(f"tips pulled to optimum (theta = {THETA:g})", x=bar_x + bar_w + 40,
               y=bar_y + bar_h / 2, font_size=FS_ANNOT, color=INK)

    d.save_svg(f"{OUT_STEM}.svg")
    d.save_png(f"{OUT_STEM}.png", dpi=300)
    print(f"wrote {OUT_STEM}.svg / .png  ({N_TIPS} tips, trait range [{vmin:+.2f}, {vmax:+.2f}], theta={THETA})")


if __name__ == "__main__":
    main()
