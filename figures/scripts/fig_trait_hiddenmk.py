"""Trait figure: a hidden-state Mk model (hidden rate classes).

A binary observed character evolves along the tree, but its transition rate
depends on an unobserved 'hidden' class that itself switches between slow and
fast. We draw a double-track stochastic map: the thick branch is coloured by the
OBSERVED state, and a thin track beneath it by the HIDDEN rate class. Observed-
state switches cluster on the 'fast' (dark-track) branches and are rare on the
'slow' (light-track) ones — the whole point of a hidden-rates model.

Run:  python figures/scripts/fig_trait_hiddenmk.py
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, "/Users/aadria/Desktop/CLAUDE/ZOMBI2-traits")   # traits-enabled zombi2

import drawsvg as draw

import phylustrator as ph
from zombi2 import BirthDeath, HiddenStateMk, simulate_species_tree, simulate_traits

from model_common import zombi_to_ete3
from zombi_style import INK, species_style

OUT_STEM = Path("/Users/aadria/Desktop/CLAUDE/ZOMBI2/figures/trait_hiddenmk/trait_hiddenmk")

N_TIPS, AGE, TREE_SEED = 12, 1.0, 1
OBS_STATES, HID_STATES = ["0", "1"], ["slow", "fast"]
TRAIT_SEED = 9
MAIN_W, TRACK_W, TRACK_OFF = 5.5, 4.0, 8.0
OBS_PALETTE = {"0": "#4477AA", "1": "#EE6677"}          # observed character: blue / red
# hidden rate class: two DISCRETE states, both clearly visible (not a faint gradient)
HID_SHADE = {"slow": "#999999", "fast": "#000000"}      # grey / black


def draw_double_simmap(d, tree, full_by_name, trans_by_name):
    xt = lambda t: d.root_x + t * d.sf                       # noqa: E731

    for n in tree.traverse():
        if n.is_root():
            continue
        y = n.coordinates[1]
        cur_t, cur_full = n.up.tfo, full_by_name[n.up.name]
        segs = []
        for tt, _frm, to in sorted(trans_by_name.get(n.name, []), key=lambda c: c[0]):
            segs.append((cur_t, tt, cur_full))
            cur_t, cur_full = tt, to
        segs.append((cur_t, n.tfo, cur_full))
        for a, b, (obs, hid) in segs:
            d.drawing.append(draw.Line(xt(a), y + TRACK_OFF, xt(b), y + TRACK_OFF,
                                       stroke=HID_SHADE[hid], stroke_width=TRACK_W, stroke_linecap="butt"))
            d.drawing.append(draw.Line(xt(a), y, xt(b), y, stroke=OBS_PALETTE[obs],
                                       stroke_width=MAIN_W, stroke_linecap="butt"))

    for n in tree.traverse("postorder"):        # vertical connectors, coloured by observed state
        if not n.is_leaf():
            x0, _ = n.coordinates
            ys = [c.coordinates[1] for c in n.children]
            d.drawing.append(draw.Line(x0, min(ys), x0, max(ys), stroke=OBS_PALETTE[full_by_name[n.name][0]],
                                       stroke_width=MAIN_W, stroke_linecap="round"))


def add_legend(d, x, y, font_size=14):
    """Compact 2-column key (observed character | hidden rate class) that fits in the
    top strip above the tips."""
    fam = d.style.font_family
    L, c2 = 28, x + 210
    d.drawing.append(draw.Text("Hidden-state Mk", font_size + 3, x, y, font_weight="bold",
                               font_family=fam, text_anchor="start"))
    d.drawing.append(draw.Text("observed character (thick) + hidden rate class (thin)", font_size - 2,
                               x, y + 16, font_family=fam, text_anchor="start", fill=INK))

    def swatch(sx, sy, label, color, w):
        d.drawing.append(draw.Line(sx, sy, sx + L, sy, stroke=color, stroke_width=w, stroke_linecap="round"))
        d.drawing.append(draw.Text(label, font_size, sx + L + 10, sy, font_family=fam,
                                   text_anchor="start", dominant_baseline="middle"))

    ry = y + 34
    swatch(x, ry, "Observed 0", OBS_PALETTE["0"], MAIN_W)
    swatch(c2, ry, "Hidden: slow", HID_SHADE["slow"], TRACK_W)
    ry += 20
    swatch(x, ry, "Observed 1", OBS_PALETTE["1"], MAIN_W)
    swatch(c2, ry, "Hidden: fast", HID_SHADE["fast"], TRACK_W)


def main():
    ztree = simulate_species_tree(BirthDeath(1.0, 0.3), n_tips=N_TIPS, age=AGE,
                                  direction="backward", seed=TREE_SEED)
    Qslow = [[-0.4, 0.4], [0.4, -0.4]]
    Qfast = [[-3.5, 3.5], [3.5, -3.5]]
    model = HiddenStateMk(observed_rates=[Qslow, Qfast], hidden_rate=0.7,
                          observed_states=OBS_STATES, hidden_states=HID_STATES)
    res = simulate_traits(ztree, model, seed=TRAIT_SEED)
    full_by_name = {n.name: res.full_label(i) for n, i in res.node_values.items()}
    trans_by_name = {}
    for node, t, frm, to in res.changes():
        trans_by_name.setdefault(node.name, []).append((t, frm, to))

    tree = zombi_to_ete3(ztree)
    for n in tree.traverse("preorder"):
        n.add_feature("tfo", 0.0 if n.is_root() else n.up.tfo + n.dist)

    style = species_style(height=max(680, 46 * N_TIPS + 180), font_size=14)
    d = ph.VerticalTreeDrawer(tree, style=style)
    d._calculate_layout()
    draw_double_simmap(d, tree, full_by_name, trans_by_name)

    tip_x = max(l.coordinates[0] for l in tree.get_leaves()) + 22
    for lf in tree.get_leaves():
        d._draw_shape_at(tip_x, lf.coordinates[1], "square", OBS_PALETTE[full_by_name[lf.name][0]],
                         r=8, stroke=INK, stroke_width=1.0)

    d.add_time_axis(ticks=[round(AGE * i / 4, 6) for i in range(5)],
                    tick_labels=[f"{AGE * i / 4:.2f}" for i in range(5)],
                    label="Time (root to present)", tick_size=6.0, padding=14.0, stroke_width=1.6)
    add_legend(d, x=-style.width / 2 + 34, y=-style.height / 2 + 22)

    n_obs = sum(1 for c in res.changes() if c[2][0] != c[3][0])
    d.save_svg(f"{OUT_STEM}.svg")
    d.save_png(f"{OUT_STEM}.png", dpi=300)
    print(f"wrote {OUT_STEM}.svg / .png  ({N_TIPS} tips, {n_obs} observed-state changes)")


if __name__ == "__main__":
    main()
