"""Trait figure: the threshold model on the tree (double-track, like Hidden-state Mk).

Same idea as the hidden-state Mk figure, but the two tracks are the two levels of
the threshold model:

  * thick main line = the OBSERVED discrete state (below / above the threshold),
    a simmap that switches colour exactly where the liability crosses the threshold
  * thin track beneath = the underlying continuous LIABILITY (viridis), the thing
    that actually evolves by Brownian motion

The threshold value is marked on the liability colour bar; the discrete state flips
precisely where the thin (liability) track passes that value.

Run:  python figures/scripts/fig_trait_threshold_tree.py
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, "/Users/aadria/Desktop/CLAUDE/ZOMBI2-traits")   # traits-enabled zombi2

import drawsvg as draw

import phylustrator as ph
from zombi2 import BirthDeath, ThresholdModel, simulate_species_tree, simulate_traits

from fig_trait_bm import hexc, viridis
from model_common import zombi_to_ete3
from zombi_style import INK, species_style
from phylustrator.utils import generate_id

OUT_STEM = Path("/Users/aadria/Desktop/CLAUDE/ZOMBI2/figures/trait_threshold_tree/trait_threshold_tree")

N_TIPS, AGE, TREE_SEED = 12, 1.0, 2
SIGMA2, THRESH, TRAIT_SEED = 1.6, 0.0, 3
STATES = ["below", "above"]
STATE_PALETTE = {"below": "#4477AA", "above": "#EE6677"}
MAIN_W, TRACK_W, TRACK_OFF = 5.5, 3.5, 8.0


def state_of(v):
    return STATES[1] if v >= THRESH else STATES[0]


def main():
    ztree = simulate_species_tree(BirthDeath(1.0, 0.3), n_tips=N_TIPS, age=AGE,
                                  direction="backward", seed=TREE_SEED)
    res = simulate_traits(ztree, ThresholdModel(thresholds=[THRESH], sigma2=SIGMA2, x0=0.0,
                                                states=STATES), seed=TRAIT_SEED)
    liab = {n.name: float(v) for n, v in res.node_values.items()}
    vmin, vmax = min(liab.values()), max(liab.values())
    norm = lambda v: (v - vmin) / (vmax - vmin) if vmax > vmin else 0.5   # noqa: E731

    tree = zombi_to_ete3(ztree)
    for n in tree.traverse("preorder"):
        n.add_feature("tfo", 0.0 if n.is_root() else n.up.tfo + n.dist)

    style = species_style(height=max(680, 46 * N_TIPS + 180), font_size=14)
    d = ph.VerticalTreeDrawer(tree, style=style)
    d._calculate_layout()
    xt = lambda t: d.root_x + t * d.sf                       # noqa: E731

    for n in tree.traverse():
        if n.is_root():
            continue
        y = n.coordinates[1]
        lp, ln = liab[n.up.name], liab[n.name]
        # thin liability track (continuous viridis gradient)
        gid = generate_id("liab")
        g = draw.LinearGradient(xt(n.up.tfo), 0, xt(n.tfo), 0, id=gid)
        g.add_stop(0, hexc(viridis(norm(lp))))
        g.add_stop(1, hexc(viridis(norm(ln))))
        d.drawing.append(g)
        d.drawing.append(draw.Line(xt(n.up.tfo), y + TRACK_OFF, xt(n.tfo), y + TRACK_OFF,
                                   stroke=g, stroke_width=TRACK_W, stroke_linecap="butt"))
        # thick discrete-state simmap (split at the threshold crossing)
        sp, sn = state_of(lp), state_of(ln)
        if sp == sn:
            d.drawing.append(draw.Line(xt(n.up.tfo), y, xt(n.tfo), y, stroke=STATE_PALETTE[sn],
                                       stroke_width=MAIN_W, stroke_linecap="butt"))
        else:
            tc = n.up.tfo + (n.tfo - n.up.tfo) * (THRESH - lp) / (ln - lp)
            d.drawing.append(draw.Line(xt(n.up.tfo), y, xt(tc), y, stroke=STATE_PALETTE[sp],
                                       stroke_width=MAIN_W, stroke_linecap="butt"))
            d.drawing.append(draw.Line(xt(tc), y, xt(n.tfo), y, stroke=STATE_PALETTE[sn],
                                       stroke_width=MAIN_W, stroke_linecap="butt"))

    for n in tree.traverse("postorder"):    # vertical connectors by discrete state
        if not n.is_leaf():
            x0, _ = n.coordinates
            ys = [c.coordinates[1] for c in n.children]
            d.drawing.append(draw.Line(x0, min(ys), x0, max(ys),
                                       stroke=STATE_PALETTE[state_of(liab[n.name])],
                                       stroke_width=MAIN_W, stroke_linecap="round"))

    tip_x = max(l.coordinates[0] for l in tree.get_leaves()) + 22
    for lf in tree.get_leaves():
        d._draw_shape_at(tip_x, lf.coordinates[1], "square", STATE_PALETTE[state_of(liab[lf.name])],
                         r=8, stroke=INK, stroke_width=1.0)

    d.add_time_axis(ticks=[round(AGE * i / 4, 6) for i in range(5)],
                    tick_labels=[f"{AGE * i / 4:.2f}" for i in range(5)],
                    label="Time (root to present)", tick_size=6.0, padding=14.0, stroke_width=1.6)
    add_legend(d, style, x=-style.width / 2 + 34, y=-style.height / 2 + 24, norm=norm,
               vmin=vmin, vmax=vmax)

    d.save_svg(f"{OUT_STEM}.svg")
    d.save_png(f"{OUT_STEM}.png", dpi=300)
    n_cross = sum(1 for n in tree.traverse() if not n.is_root()
                  and state_of(liab[n.name]) != state_of(liab[n.up.name]))
    print(f"wrote {OUT_STEM}.svg / .png  ({N_TIPS} tips, {n_cross} threshold crossings)")


def add_legend(d, style, x, y, norm, vmin, vmax, font_size=14):
    fam = style.font_family
    d.drawing.append(draw.Text("Threshold model", font_size + 3, x, y, font_weight="bold",
                               font_family=fam, text_anchor="start"))
    d.drawing.append(draw.Text("observed state (thick) + liability (thin)", font_size - 2, x, y + 16,
                               font_family=fam, text_anchor="start", fill=INK))
    # discrete state swatches
    ry, L = y + 38, 28
    for i, s in enumerate(STATES):
        cy = ry + i * 20
        d.drawing.append(draw.Line(x, cy, x + L, cy, stroke=STATE_PALETTE[s], stroke_width=MAIN_W,
                                   stroke_linecap="round"))
        d.drawing.append(draw.Text(f"state: {s}", font_size, x + L + 10, cy, font_family=fam,
                                   text_anchor="start", dominant_baseline="middle"))
    # liability colour bar with a threshold tick
    bx, by, bw, bh = x + 180, y + 34, 130, 13
    gid = generate_id("liabbar")
    g = draw.LinearGradient(bx, 0, bx + bw, 0, id=gid)
    for t in (0.0, 0.25, 0.5, 0.75, 1.0):
        r, gg, b = viridis(t)
        g.add_stop(t, "#%02x%02x%02x" % (r, gg, b))
    d.drawing.append(g)
    d.drawing.append(draw.Text("Liability", font_size - 1, bx, by - 6, font_family=fam,
                               text_anchor="start", fill=INK))
    d.drawing.append(draw.Rectangle(bx, by, bw, bh, fill=g, stroke=INK, stroke_width=0.5))
    d.drawing.append(draw.Text(f"{vmin:+.1f}", 12, bx, by + bh + 13, font_family=fam, text_anchor="start", fill=INK))
    d.drawing.append(draw.Text(f"{vmax:+.1f}", 12, bx + bw, by + bh + 13, font_family=fam, text_anchor="end", fill=INK))
    txp = bx + max(0.0, min(1.0, norm(THRESH))) * bw          # threshold tick on the bar
    d.drawing.append(draw.Line(txp, by - 3, txp, by + bh + 3, stroke=INK, stroke_width=1.6))
    d.drawing.append(draw.Text("threshold", 12, txp, by + bh + 26, font_family=fam, text_anchor="middle", fill=INK))


if __name__ == "__main__":
    main()
