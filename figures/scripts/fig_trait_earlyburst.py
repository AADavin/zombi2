"""Trait figure: an early-burst (ACDC) continuous trait.

Brownian motion whose rate decays through time: sigma^2(t) = sigma2 * e^{rate*t}
with rate < 0. Most disparity accumulates early, so the deep branches carry big
trait changes (clades diverge early) while recent branches barely change (tips
within a young clade stay similar) — the opposite of BM's even spread.

A grey strip along the time axis shows the rate decaying from fast (dark, early)
to slow (light, late). Same tree as the BM / OU figures for comparison.

Run:  python figures/scripts/fig_trait_earlyburst.py
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, "/Users/aadria/Desktop/CLAUDE/ZOMBI2-traits")   # traits-enabled zombi2

import drawsvg as draw

import phylustrator as ph
from zombi2 import BirthDeath, EarlyBurst, simulate_species_tree, simulate_traits

from fig_trait_bm import color_bar, hexc, viridis
from model_common import zombi_to_ete3
from zombi_style import INK, species_style

OUT_STEM = Path("/Users/aadria/Desktop/CLAUDE/ZOMBI2/figures/trait_earlyburst/trait_earlyburst")

N_TIPS, AGE, TREE_SEED = 14, 1.0, 3          # same tree as fig_trait_bm / fig_trait_ou
SIGMA2, RATE, TRAIT_SEED = 8.0, -5.0, 7      # rate < 0 -> early burst
BRANCH_W = 5.0


def main():
    ztree = simulate_species_tree(BirthDeath(1.0, 0.3), n_tips=N_TIPS, age=AGE,
                                  direction="backward", seed=TREE_SEED)
    res = simulate_traits(ztree, EarlyBurst(sigma2=SIGMA2, rate=RATE, x0=0.0), seed=TRAIT_SEED)
    name2val = {n.name: float(v) for n, v in res.node_values.items()}
    vmin, vmax = min(name2val.values()), max(name2val.values())
    norm = lambda v: (v - vmin) / (vmax - vmin) if vmax > vmin else 0.5   # noqa: E731
    node_to_rgb = {name: viridis(norm(v)) for name, v in name2val.items()}

    tree = zombi_to_ete3(ztree)
    style = species_style(height=max(680, 44 * N_TIPS + 180), font_size=14)
    d = ph.VerticalTreeDrawer(tree, style=style)
    d._calculate_layout()

    # a thin grey strip along the time axis: rate sigma^2(t) fast (dark, early) -> slow (light)
    ys = [l.y_coord for l in tree.get_leaves()]
    strip_y = max(ys) + 18
    x0, x1 = d.root_x, d.root_x + AGE * d.sf
    grad = draw.LinearGradient(x0, 0, x1, 0, id="rate_decay")
    grad.add_stop(0.0, "#3a3a3a")
    grad.add_stop(1.0, "#ededed")
    d.drawing.append(grad)
    d.drawing.append(draw.Rectangle(x0, strip_y, x1 - x0, 11, fill=grad, stroke=INK, stroke_width=0.5))
    d.drawing.append(draw.Text("rate σ²(t):  fast", 12, x0 + 4, strip_y + 8, font_family=style.font_family,
                               text_anchor="start", fill="white"))
    d.drawing.append(draw.Text("slow", 12, x1 - 4, strip_y + 8, font_family=style.font_family,
                               text_anchor="end", fill=INK))

    d.plot_continuous_variable(node_to_rgb, stroke_width=BRANCH_W)

    tip_x = max(l.coordinates[0] for l in tree.get_leaves()) + 22
    for lf in tree.get_leaves():
        d._draw_shape_at(tip_x, lf.coordinates[1], "square", hexc(node_to_rgb[lf.name]), r=8,
                         stroke=INK, stroke_width=1.0)

    d.add_time_axis(ticks=[round(AGE * i / 4, 6) for i in range(5)],
                    tick_labels=[f"{AGE * i / 4:.2f}" for i in range(5)],
                    label="Time (root to present)", tick_size=6.0, padding=46.0, stroke_width=1.6)
    color_bar(d, x=-style.width / 2 + 34, y=-style.height / 2 + 44, vmin=vmin, vmax=vmax,
              title="Trait value (early burst)")

    d.save_svg(f"{OUT_STEM}.svg")
    d.save_png(f"{OUT_STEM}.png", dpi=300)
    print(f"wrote {OUT_STEM}.svg / .png  ({N_TIPS} tips, range [{vmin:+.2f}, {vmax:+.2f}])")


if __name__ == "__main__":
    main()
