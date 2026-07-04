"""Trait figure: a discrete character under the Mk model (stochastic map).

A 3-state discrete trait evolves along the tree as a continuous-time Markov chain
(Mk). Because ZOMBI2 records every transition and its time, we draw a true
stochastic character map: each branch is painted in segments, switching colour
sharply at the exact time a state change happened.

  * each state -> a colour (categorical)
  * a state change on a branch -> a sharp colour switch at its time
  * internal nodes / tips -> coloured by their (ancestral / observed) state

Run:  python figures/scripts/fig_trait_mk.py
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, "/Users/aadria/Desktop/CLAUDE/ZOMBI2-traits")   # traits-enabled zombi2

import drawsvg as draw

import phylustrator as ph
from zombi2 import BirthDeath, Mk, simulate_species_tree, simulate_traits

from model_common import zombi_to_ete3
from zombi_style import INK, species_style

OUT_STEM = Path("/Users/aadria/Desktop/CLAUDE/ZOMBI2/figures/trait_mk/trait_mk")

N_TIPS, AGE, TREE_SEED = 12, 1.0, 1
STATES = ["A", "B", "C"]
Q_RATE, TRAIT_SEED = 0.55, 2
BRANCH_W = 5.0
# colour-blind-friendly categorical palette (Paul Tol 'bright')
PALETTE = {"A": "#4477AA", "B": "#EE6677", "C": "#228833"}


def build_Q(q, k):
    return [[-(k - 1) * q if i == j else q for j in range(k)] for i in range(k)]


def draw_simmap(d, tree, state_of, trans_by_name):
    sw = BRANCH_W
    xt = lambda t: d.root_x + t * d.sf                       # noqa: E731  time -> x

    # root stub in the root's state
    rx, ry = tree.coordinates
    d.drawing.append(draw.Line(rx - d.style.root_stub_length, ry, rx, ry,
                               stroke=PALETTE[state_of[tree.name]], stroke_width=sw, stroke_linecap="round"))

    for n in tree.traverse():
        if n.is_root():
            continue
        y = n.coordinates[1]
        cur_t, cur_state = n.up.tfo, state_of[n.up.name]
        for tt, _frm, to in sorted(trans_by_name.get(n.name, []), key=lambda c: c[0]):
            d.drawing.append(draw.Line(xt(cur_t), y, xt(tt), y, stroke=PALETTE[cur_state],
                                       stroke_width=sw, stroke_linecap="butt"))
            cur_t, cur_state = tt, to
        d.drawing.append(draw.Line(xt(cur_t), y, xt(n.tfo), y, stroke=PALETTE[cur_state],
                                   stroke_width=sw, stroke_linecap="butt"))

    # vertical connectors, coloured by the node's (ancestral) state
    for n in tree.traverse("postorder"):
        if not n.is_leaf():
            x0, _ = n.coordinates
            ys = [c.coordinates[1] for c in n.children]
            d.drawing.append(draw.Line(x0, min(ys), x0, max(ys), stroke=PALETTE[state_of[n.name]],
                                       stroke_width=sw, stroke_linecap="round"))


def main():
    ztree = simulate_species_tree(BirthDeath(1.0, 0.3), n_tips=N_TIPS, age=AGE,
                                  direction="backward", seed=TREE_SEED)
    model = Mk(build_Q(Q_RATE, len(STATES)), states=STATES)
    res = simulate_traits(ztree, model, seed=TRAIT_SEED)
    # node_values are integer state indices; map to labels
    state_of = {n.name: STATES[i] for n, i in res.node_values.items()}
    trans_by_name = {}
    for node, t, frm, to in res.changes():
        trans_by_name.setdefault(node.name, []).append((t, frm, to))

    tree = zombi_to_ete3(ztree)
    present = 0.0
    for n in tree.traverse("preorder"):
        n.add_feature("tfo", 0.0 if n.is_root() else n.up.tfo + n.dist)
        present = max(present, n.tfo)

    style = species_style(height=max(680, 44 * N_TIPS + 180), font_size=14)
    d = ph.VerticalTreeDrawer(tree, style=style)
    d._calculate_layout()
    draw_simmap(d, tree, state_of, trans_by_name)

    # tip chips coloured by observed state
    tip_x = max(l.coordinates[0] for l in tree.get_leaves()) + 22
    for lf in tree.get_leaves():
        d._draw_shape_at(tip_x, lf.coordinates[1], "square", PALETTE[state_of[lf.name]], r=8,
                         stroke=INK, stroke_width=1.0)

    d.add_time_axis(ticks=[round(AGE * i / 4, 6) for i in range(5)],
                    tick_labels=[f"{AGE * i / 4:.2f}" for i in range(5)],
                    label="Time (root to present)", tick_size=6.0, padding=14.0, stroke_width=1.6)
    d.add_categorical_legend({f"State {s}": PALETTE[s] for s in STATES}, title="Mk trait (3 states)",
                             x=-style.width / 2 + 34, y=-style.height / 2 + 40, r=7)

    d.save_svg(f"{OUT_STEM}.svg")
    d.save_png(f"{OUT_STEM}.png", dpi=300)
    print(f"wrote {OUT_STEM}.svg / .png  ({N_TIPS} tips, {len(res.changes())} transitions)")


if __name__ == "__main__":
    main()
