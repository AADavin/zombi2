"""Trait figure: an early-burst (ACDC) continuous trait.

Brownian motion whose rate decays through time: sigma^2(t) = sigma2 * exp(rate*t) with rate < 0.
Most disparity accumulates early, so the deep branches carry big trait changes (clades diverge
early) while recent branches barely change (tips within a young clade stay similar) — the opposite
of Brownian motion's even spread. A grey strip under the tree shows the rate decaying from fast
(dark, early) to slow (light, late). Same painted-tree style as the BM / OU figures.

Run:  python figures/scripts/fig_trait_earlyburst.py
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import drawsvg as draw

import phylustrator as ph
from zombi2.species import BirthDeath, simulate_species_tree
from zombi2.traits import EarlyBurst, simulate_traits

from fig_trait_bm import color_bar, hexc, viridis
from model_common import zombi_to_ete3
from zombi_style import INK, MUTED, species_style, FS_TITLE, FS_LABEL, FS_TICK

OUT_STEM = Path(__file__).resolve().parent.parent / "trait_earlyburst" / "trait_earlyburst"
OUT_STEM.parent.mkdir(parents=True, exist_ok=True)

N_TIPS, AGE, TREE_SEED = 28, 1.0, 3
SIGMA2, RATE, X0 = 6.0, -4.0, 0.0            # rate < 0 -> early burst
BRANCH_W = 4.0


def pick_seed(ztree):
    """A trait seed whose tip values spread evenly across the range (good gradient contrast)."""
    import statistics
    best = None
    for s in range(1, 40):
        res = simulate_traits(ztree, EarlyBurst(sigma2=SIGMA2, rate=RATE, x0=X0), seed=s)
        vals = [float(v) for n, v in res.node_values.items()]
        lo, hi = min(vals), max(vals)
        if hi <= lo:
            continue
        tips = [float(v) for n, v in res.node_values.items() if n.is_leaf()]
        score = statistics.pstdev(tips) / (hi - lo)
        if best is None or score > best[0]:
            best = (score, res, lo, hi)
    return best[1], best[2], best[3]


def main():
    ztree = simulate_species_tree(BirthDeath(1.0, 0.3), n_tips=N_TIPS, age=AGE,
                                  direction="backward", seed=TREE_SEED)
    res, vmin, vmax = pick_seed(ztree)
    name2val = {n.name: float(v) for n, v in res.node_values.items()}
    norm = lambda v: (v - vmin) / (vmax - vmin) if vmax > vmin else 0.5   # noqa: E731
    node_to_rgb = {name: viridis(norm(v)) for name, v in name2val.items()}

    tree = zombi_to_ete3(ztree)
    style = species_style(width=1180, height=740, margin=120, font_size=FS_TICK)
    d = ph.VerticalTreeDrawer(tree, style=style)
    d._calculate_layout()

    # rate strip under the tree: sigma^2(t) fast (dark, early) -> slow (light, late).
    # The strip is made a good deal TALLER than the label height so the "rate: fast"
    # and "slow" labels sit comfortably inside it (vertically centred) instead of
    # spilling above/below the band.
    strip_y = max(l.coordinates[1] for l in tree.get_leaves()) + 20
    strip_h = 30
    strip_pad = 8
    x0, x1 = d.root_x - strip_pad, d.root_x + AGE * d.sf + strip_pad
    grad = draw.LinearGradient(x0, 0, x1, 0, id="rate_decay")
    grad.add_stop(0.0, "#3a3a3a")
    grad.add_stop(1.0, "#ededed")
    d.drawing.append(grad)
    d.drawing.append(draw.Rectangle(x0, strip_y, x1 - x0, strip_h, fill=grad, stroke=INK, stroke_width=0.5))
    # Baseline is placed explicitly (renderer-independent) so the text's visual centre
    # (~0.35*font above the baseline) lands on the strip's mid-height.
    label_y = strip_y + strip_h / 2 + FS_TICK * 0.35
    d.drawing.append(draw.Text("rate: fast", FS_TICK, x0 + 10, label_y,
                               font_family=style.font_family,
                               text_anchor="start", fill="white"))
    d.drawing.append(draw.Text("slow", FS_TICK, x1 - 10, label_y,
                               font_family=style.font_family,
                               text_anchor="end", fill=INK))

    d.plot_continuous_variable(node_to_rgb, stroke_width=BRANCH_W)

    tip_x = max(l.coordinates[0] for l in tree.get_leaves()) + 22
    for lf in tree.get_leaves():
        d._draw_shape_at(tip_x, lf.coordinates[1], "square", hexc(node_to_rgb[lf.name]), r=6,
                         stroke=INK, stroke_width=1.0)

    d.add_time_axis(ticks=[round(AGE * i / 4, 6) for i in range(5)],
                    tick_labels=[f"{AGE * i / 4:.2f}" for i in range(5)],
                    label="Time (root to present)", tick_size=6.0, padding=48.0, stroke_width=1.6)

    fam = style.font_family
    d.drawing.append(draw.Text("A trait under early burst (ACDC)", FS_TITLE, 0, -style.height / 2 + 42,
                               font_weight="bold", font_family=fam, text_anchor="middle",
                               dominant_baseline="central", fill=INK))
    color_bar(d, x=-style.width / 2 + 32, y=-style.height / 2 + 96, w=220, h=18,
              vmin=vmin, vmax=vmax, title="Trait value")

    d.save_svg(f"{OUT_STEM}.svg")
    d.save_png(f"{OUT_STEM}.png", dpi=300)
    print(f"wrote {OUT_STEM}.svg / .png  ({N_TIPS} tips, range [{vmin:+.2f}, {vmax:+.2f}])")


if __name__ == "__main__":
    main()
