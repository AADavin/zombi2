"""Trait figure: the threshold model — a latent liability read out as a discrete state.

Felsenstein's threshold model bridges the continuous and discrete trait families: an unobserved
**liability** evolves by Brownian motion, and the observed discrete state is simply which side of
a threshold the liability is on. Drawn in the same painted-tree style as the Brownian-motion and
Ornstein-Uhlenbeck figures, so the family reads together:

  * branches are painted by the *continuous* liability (viridis), the latent evolving value;
  * the threshold is marked on the colour bar — colours left of it are below (state 0), right
    of it are above (state 1);
  * each tip carries the *discrete* observed state it ends in: open square = 0, filled = 1.

The continuous gradient flowing into a binary tip pattern is exactly the model's point.

Run:  /Users/aadria/miniconda3/bin/python figures/scripts/fig_trait_threshold.py
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))   # local: zombi_style, model_common, fig_trait_bm

import drawsvg as draw

import phylustrator as ph
from zombi2.species import BirthDeath, simulate_species_tree
from zombi2.traits import ThresholdModel, simulate_traits

from fig_trait_bm import VIRIDIS, hexc, viridis       # shared colormap
from model_common import zombi_to_ete3
from zombi_style import INK, species_style, FS_TITLE, FS_LABEL, FS_ANNOT, FS_TICK

OUT_STEM = Path("/Users/aadria/Desktop/CLAUDE/ZOMBI2/figures/trait_threshold/trait_threshold")

N_TIPS, AGE, TREE_SEED = 28, 1.0, 3
SIGMA2, X0, THRESHOLD = 1.2, 0.0, 0.0
BRANCH_W = 4.0


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


def pick_seed(ztree, model):
    """A trait seed giving a readable mix of both states (near-balanced, both present)."""
    best = None
    for s in range(1, 40):
        res = simulate_traits(ztree, model, seed=s)
        states = list(res.labeled_values().values())
        n1 = sum(states)
        if 0 < n1 < len(states):
            bal = -abs(n1 - len(states) / 2)          # prefer a balanced split
            if best is None or bal > best[0]:
                best = (bal, s, res)
    return best[1], best[2]


def main():
    ztree = simulate_species_tree(BirthDeath(1.0, 0.3), n_tips=N_TIPS, age=AGE,
                                  direction="backward", seed=TREE_SEED)
    model = ThresholdModel(thresholds=[THRESHOLD], sigma2=SIGMA2, x0=X0)
    seed, res = pick_seed(ztree, model)
    liab = {n.name: float(v) for n, v in res.node_values.items()}          # continuous liability
    disc = {n.name: int(s) for n, s in res.labeled_values().items()}        # observed 0/1 (tips)
    vmin, vmax = min(liab.values()), max(liab.values())
    norm = lambda v: (v - vmin) / (vmax - vmin) if vmax > vmin else 0.5     # noqa: E731
    node_to_rgb = {name: viridis(norm(v)) for name, v in liab.items()}

    tree = zombi_to_ete3(ztree)
    style = species_style(width=1180, height=740, margin=120, font_size=FS_TICK)
    d = ph.VerticalTreeDrawer(tree, style=style)
    d._calculate_layout()
    d.plot_continuous_variable(node_to_rgb, stroke_width=BRANCH_W)

    # tips carry the DISCRETE observed state: open = below (0), filled = above (1)
    tip_x = max(l.coordinates[0] for l in tree.get_leaves()) + 22
    for lf in tree.get_leaves():
        _, y = lf.coordinates
        on = disc[lf.name] == 1
        d._draw_shape_at(tip_x, y, "square", INK if on else "white", r=6, stroke=INK, stroke_width=1.3)

    d.add_time_axis(ticks=[round(AGE * i / 4, 6) for i in range(5)],
                    tick_labels=[f"{AGE * i / 4:.2f}" for i in range(5)],
                    label="Time (root to present)", tick_size=6.0, padding=16.0, stroke_width=1.6)

    # title centered; colour bar (the "data") top-left, both clear of the tree
    left = -style.width / 2 + 30
    d.drawing.append(draw.Text("The threshold model", FS_TITLE, 0, -style.height / 2 + 42,
                               font_weight="bold", font_family=style.font_family,
                               text_anchor="middle", dominant_baseline="central", fill=INK))
    bar_x, bar_y, bar_w, bar_h = left + 2, -style.height / 2 + 96, 220, 18
    color_bar(d, x=bar_x, y=bar_y, w=bar_w, h=bar_h, vmin=vmin, vmax=vmax, title="Liability")
    # threshold mark on the bar (label below, in the gap between the two range labels)
    tx = bar_x + max(0.0, min(1.0, (THRESHOLD - vmin) / (vmax - vmin))) * bar_w
    d.drawing.append(draw.Line(tx, bar_y - 4, tx, bar_y + bar_h + 4, stroke=INK, stroke_width=2.6))
    d.drawing.append(draw.Text("threshold", FS_TICK, tx, bar_y + bar_h + 22, font_family=style.font_family,
                               text_anchor="middle", fill=INK))

    # tip-state legend (discrete readout), clear of the tree
    lx, ly = bar_x, bar_y + bar_h + 58
    fam = style.font_family
    d.drawing.append(draw.Rectangle(lx, ly, 16, 16, fill=INK, stroke=INK, stroke_width=1.3))
    d.drawing.append(draw.Text("state 1  (above)", FS_ANNOT, lx + 26, ly + 8, font_family=fam,
                               text_anchor="start", dominant_baseline="central", fill=INK))
    d.drawing.append(draw.Rectangle(lx + 190, ly, 16, 16, fill="white", stroke=INK, stroke_width=1.3))
    d.drawing.append(draw.Text("state 0  (below)", FS_ANNOT, lx + 216, ly + 8, font_family=fam,
                               text_anchor="start", dominant_baseline="central", fill=INK))

    d.save_svg(f"{OUT_STEM}.svg")
    d.save_png(f"{OUT_STEM}.png", dpi=300)
    print(f"wrote {OUT_STEM}.svg / .png  ({N_TIPS} tips, trait seed {seed}, "
          f"liability [{vmin:+.2f}, {vmax:+.2f}])")


if __name__ == "__main__":
    main()
