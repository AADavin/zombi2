"""Trait figure: multi-optimum OU — a trait chasing a different optimum per regime.

A discrete character paints *regimes* on the tree (an Mk stochastic map); a continuous trait then
evolves under Ornstein-Uhlenbeck with a **different optimum per regime**. Each lineage is pulled
toward whichever optimum its current regime sits under, so the tips split into two clusters.

  * branches painted by the continuous trait value (viridis);
  * the two optima (theta) marked on the colour bar;
  * a regime bar at each tip (dark = regime A, light = regime B) — regime-A tips sit near the low
    optimum, regime-B tips near the high one.

House style: painted tree, one centered title, ASCII text.

Run:  /Users/aadria/miniconda3/bin/python figures/scripts/fig_trait_multioptimum.py
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import drawsvg as draw

import phylustrator as ph
from zombi2.species import BirthDeath, simulate_species_tree
from zombi2.traits import Mk, MultiOptimumOU, simulate_traits

from fig_trait_bm import VIRIDIS, hexc, viridis, color_bar
from model_common import zombi_to_ete3
from zombi_style import INK, MUTED, species_style, FS_TITLE, FS_LABEL, FS_ANNOT, FS_TICK

OUT_STEM = Path(__file__).resolve().parent.parent / "trait_multioptimum" / "trait_multioptimum"
OUT_STEM.parent.mkdir(parents=True, exist_ok=True)

N_TIPS, AGE, TREE_SEED = 28, 1.0, 3
THETA = [-5.0, 5.0]
ALPHA, SIGMA2 = 4.0, 0.5
REGIME_SHADE = {0: "#4a4a4a", 1: "#cfcfcf"}          # regime A (low opt) / regime B (high opt)
BRANCH_W = 4.0


def pick_seeds(ztree):
    """Regime + trait seeds giving both regimes present and a clean two-cluster split."""
    best = None
    for rs in range(1, 25):
        reg = simulate_traits(ztree, Mk.equal_rates(2, 0.7), seed=rs)
        tip_reg = [reg.node_values[n] for n in ztree.extant_leaves()]
        if len(set(tip_reg)) < 2:
            continue
        bal = -abs(sum(tip_reg) - len(tip_reg) / 2)
        mou = MultiOptimumOU(reg, theta=THETA, alpha=ALPHA, sigma2=SIGMA2)
        res = simulate_traits(ztree, mou, seed=rs + 1)
        if best is None or bal > best[0]:
            best = (bal, reg, res)
    return best[1], best[2]


def main():
    ztree = simulate_species_tree(BirthDeath(1.0, 0.3), n_tips=N_TIPS, age=AGE,
                                  direction="backward", seed=TREE_SEED)
    reg, res = pick_seeds(ztree)
    val = {n.name: float(v) for n, v in res.node_values.items()}
    regime = {n.name: int(r) for n, r in reg.node_values.items()}
    # Include both optima in the colour range so their markers fall inside the colour bar
    # (a tip need not reach theta exactly); branches use the same range, so bar and tree agree.
    # A small margin beyond the extremes keeps either theta marker off the bar's border.
    lo = min(min(val.values()), *THETA)
    hi = max(max(val.values()), *THETA)
    pad = 0.05 * (hi - lo)
    vmin, vmax = lo - pad, hi + pad
    norm = lambda v: (v - vmin) / (vmax - vmin) if vmax > vmin else 0.5     # noqa: E731
    node_to_rgb = {name: viridis(norm(v)) for name, v in val.items()}

    tree = zombi_to_ete3(ztree)
    style = species_style(width=1180, height=740, margin=120, font_size=FS_TICK)
    d = ph.VerticalTreeDrawer(tree, style=style)
    d._calculate_layout()

    # Compress the tree vertically a touch and drop its top edge, so the header band
    # (colour bar + the two-optima caption) sits clear of the top branches instead of
    # crowding them. Remap every node's y into a slightly shorter band that starts lower.
    ys = [n.coordinates[1] for n in tree.traverse()]
    y_top, y_bot = min(ys), max(ys)
    NEW_TOP, NEW_BOT = y_top + 55.0, y_bot - 15.0          # start lower, end a touch higher
    span = (y_bot - y_top) or 1.0
    for n in tree.traverse():
        x, y = n.coordinates
        n.coordinates = (x, NEW_TOP + (y - y_top) / span * (NEW_BOT - NEW_TOP))
        n.y_coord = n.coordinates[1]

    d.plot_continuous_variable(node_to_rgb, stroke_width=BRANCH_W)

    # regime bar per tip (the discrete map that sets which optimum a lineage chases)
    tip_x = max(l.coordinates[0] for l in tree.get_leaves()) + 24
    for lf in tree.get_leaves():
        _, y = lf.coordinates
        d._draw_shape_at(tip_x, y, "square", REGIME_SHADE[regime[lf.name]], r=6,
                         stroke=INK, stroke_width=1.1)

    d.add_time_axis(ticks=[round(AGE * i / 4, 6) for i in range(5)],
                    tick_labels=[f"{AGE * i / 4:.2f}" for i in range(5)],
                    label="Time (root to present)", tick_size=6.0, padding=16.0, stroke_width=1.6)

    fam = style.font_family
    d.drawing.append(draw.Text("Multi-optimum OU: a trait chases its regime's optimum", FS_TITLE, 0,
                               -style.height / 2 + 42, font_weight="bold", font_family=fam,
                               text_anchor="middle", dominant_baseline="central", fill=INK))
    bar_x, bar_y, bar_w, bar_h = -style.width / 2 + 32, -style.height / 2 + 96, 220, 18
    color_bar(d, x=bar_x, y=bar_y, w=bar_w, h=bar_h, vmin=vmin, vmax=vmax, title="Trait value")
    # mark the two optima (theta) on the colour bar; the ticks sit INSIDE the bar boundary
    for th in THETA:
        tx = bar_x + max(0.0, min(1.0, (th - vmin) / (vmax - vmin))) * bar_w
        d.drawing.append(draw.Line(tx, bar_y + 1, tx, bar_y + bar_h - 1, stroke=INK, stroke_width=2.4))
    d.drawing.append(draw.Text(f"two optima: theta = {THETA[0]:g} and {THETA[1]:g}", FS_TICK,
                               bar_x, bar_y + bar_h + 46, font_family=fam, text_anchor="start",
                               fill=MUTED))

    # regime legend, bottom-left
    lx, ly = bar_x, style.height / 2 - 168
    d.drawing.append(draw.Rectangle(lx, ly, 16, 16, fill=REGIME_SHADE[0], stroke=INK, stroke_width=1.2))
    d.drawing.append(draw.Text("regime A (low optimum)", FS_ANNOT, lx + 26, ly + 8, font_family=fam,
                               text_anchor="start", dominant_baseline="central", fill=INK))
    d.drawing.append(draw.Rectangle(lx, ly + 30, 16, 16, fill=REGIME_SHADE[1], stroke=INK, stroke_width=1.2))
    d.drawing.append(draw.Text("regime B (high optimum)", FS_ANNOT, lx + 26, ly + 38, font_family=fam,
                               text_anchor="start", dominant_baseline="central", fill=INK))

    d.save_svg(f"{OUT_STEM}.svg")
    d.save_png(f"{OUT_STEM}.png", dpi=300)
    print(f"wrote {OUT_STEM}.svg / .png  ({N_TIPS} tips, range [{vmin:+.2f}, {vmax:+.2f}])")


if __name__ == "__main__":
    main()
