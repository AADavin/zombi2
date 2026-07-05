"""Figure: the Mk model of discrete-character evolution — two panels.

A discrete character evolves as a continuous-time Markov chain over k states with a rate matrix Q.

  * Panel A (the model): the state graph — here equal-rates over three states, every transition at
    one shared rate q. The structure lives entirely in Q (equal_rates, symmetric, ordered, ARD).
  * Panel B (a realization): the exact stochastic character map — each branch is painted in the
    state it is in, switching at the instant a transition happens; tip chips give the observed
    state each lineage ends in.

House style: B&W. Three categorical states are shown as three greys (light / mid / dark), the one
readable monochrome device for a small categorical set.

Run:  /Users/aadria/miniconda3/bin/python figures/scripts/fig_trait_mk.py
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import cairosvg
import drawsvg as draw

from zombi2.species import BirthDeath, simulate_species_tree
from zombi2.traits import Mk, simulate_traits

from fig_trait_pagel import curved_arrow, _layout
from model_common import zombi_to_ete3
from zombi_style import FONT, INK, MUTED, FS_TITLE, FS_LABEL, FS_ANNOT, FS_TICK

OUT_DIR = Path(__file__).resolve().parent.parent

W, H = 1200, 660
STATES = ["A", "B", "C"]
SHADE = {"A": "#1f1f1f", "B": "#8a8a8a", "C": "#cccccc"}     # dark / mid / light
Q_RATE, N_TIPS, AGE, TREE_SEED, TRAIT_SEED = 0.7, 12, 1.0, 4, 5
NR = 34


def _lum(h):
    return 0.299 * int(h[1:3], 16) + 0.587 * int(h[3:5], 16) + 0.114 * int(h[5:7], 16)


def mk_node(d, cx, cy, label):
    fill = SHADE[label]
    d.append(draw.Circle(cx, cy, NR, fill=fill, stroke=INK, stroke_width=2.2))
    d.append(draw.Text(label, FS_LABEL, cx, cy, font_family=FONT, text_anchor="middle",
                       dominant_baseline="central", font_weight="bold",
                       fill="white" if _lum(fill) < 140 else INK))


def state_chip(d, cx, cy, label, s=12):
    d.append(draw.Rectangle(cx - s, cy - s, 2 * s, 2 * s, fill=SHADE[label], stroke=INK, stroke_width=1.5))


# --------------------------------------------------------------------------- panel A: the model
def panel_model(d, cx, cy):
    d.append(draw.Text("A   the model", FS_LABEL, cx - 150, 126, font_family=FONT,
                       text_anchor="start", fill=INK, font_weight="bold"))
    P = {"A": (cx, cy - 118), "B": (cx - 118, cy + 92), "C": (cx + 118, cy + 92)}
    # equal-rates: every pair bidirectional at rate q (label once per edge, on the outer arc)
    from math import inf  # noqa: F401  (keep import list tidy across figures)
    edges = [("A", "B"), ("B", "C"), ("A", "C")]
    for a, b in edges:
        curved_arrow(d, P[a], P[b], +1, 16, 3.4, "q")
        curved_arrow(d, P[b], P[a], -1, 16, 3.4, "")
    for lab, (x, y) in P.items():
        mk_node(d, x, y, lab)
    d.append(draw.Text("equal-rates Q (one shared rate q)", FS_TICK, cx, cy + 150, font_family=FONT,
                       text_anchor="middle", fill=MUTED))


# --------------------------------------------------------------------------- panel B: a realization
def panel_realization(d, ox, oy, pw, ph):
    ztree = simulate_species_tree(BirthDeath(1.0, 0.3), n_tips=N_TIPS, age=AGE,
                                  direction="backward", seed=TREE_SEED)
    Q = [[-(len(STATES) - 1) * Q_RATE if i == j else Q_RATE for j in range(len(STATES))]
         for i in range(len(STATES))]
    res = simulate_traits(ztree, Mk(Q, states=STATES), seed=TRAIT_SEED)
    state_of = {n.name: STATES[i] for n, i in res.node_values.items()}
    norm = lambda v: v if isinstance(v, str) else STATES[v]      # noqa: E731  (label or index)
    chg = {}
    for node, t, frm, to in res.changes():
        chg.setdefault(node.name, []).append((t, norm(to)))

    tree = zombi_to_ete3(ztree)
    tfo, present, ys, nleaf = _layout(tree)
    x_at = lambda t: ox + 50 + (t / present) * (pw - 150)        # noqa: E731
    y_at = lambda k: oy + 40 + (k / max(1, nleaf - 1)) * (ph - 90)

    d.append(draw.Text("B   a simulated realization", FS_LABEL, ox, oy - 24, font_family=FONT,
                       text_anchor="start", fill=INK, font_weight="bold"))

    for n in tree.traverse():
        if n.is_root():
            continue
        y = y_at(ys[n.name])
        t0, cur = tfo[n.up.name], state_of[n.up.name]
        for tt, to in sorted(chg.get(n.name, [])):
            d.append(draw.Line(x_at(t0), y, x_at(tt), y, stroke=SHADE[cur], stroke_width=5.0,
                               stroke_linecap="butt"))
            t0, cur = tt, to
        d.append(draw.Line(x_at(t0), y, x_at(tfo[n.name]), y, stroke=SHADE[cur], stroke_width=5.0,
                           stroke_linecap="butt"))
    for n in tree.traverse("postorder"):
        if not n.is_leaf():
            x = x_at(tfo[n.name])
            yy = [y_at(ys[c.name]) for c in n.children]
            d.append(draw.Line(x, min(yy), x, max(yy), stroke=SHADE[state_of[n.name]], stroke_width=5.0))

    colc = ox + pw - 70
    for n in tree.get_leaves():
        y = y_at(ys[n.name])
        state_chip(d, colc, y, state_of[n.name])

    base = oy + ph - 14
    d.append(draw.Line(x_at(0), base, x_at(present), base, stroke=INK, stroke_width=1.6))
    for k in range(5):
        tv = present * k / 4
        xx = x_at(tv)
        d.append(draw.Line(xx, base, xx, base + 6, stroke=INK, stroke_width=1.6))
        d.append(draw.Text(f"{tv:.2f}", FS_TICK, xx, base + 22, font_family=FONT,
                           text_anchor="middle", fill=INK))
    d.append(draw.Text("time (root to present)", FS_TICK, (x_at(0) + x_at(present)) / 2, base + 44,
                       font_family=FONT, text_anchor="middle", fill=MUTED))


# --------------------------------------------------------------------------- render
def render():
    d = draw.Drawing(W, H, origin=(0, 0))
    d.append(draw.Rectangle(0, 0, W, H, fill="white"))
    d.append(draw.Text("A discrete character under the Mk model", FS_TITLE, W / 2, 46,
                       font_family=FONT, text_anchor="middle", font_weight="bold", fill=INK))
    # categorical legend (three greys)
    lx = W / 2 - 150
    for i, s in enumerate(STATES):
        x = lx + i * 110
        d.append(draw.Rectangle(x, 74, 18, 18, fill=SHADE[s], stroke=INK, stroke_width=1.4))
        d.append(draw.Text(f"state {s}", FS_TICK, x + 26, 83, font_family=FONT, text_anchor="start",
                           dominant_baseline="central", fill=INK))

    panel_model(d, 300, 360)
    panel_realization(d, 560, 150, 580, 470)

    name = "trait_mk"
    out = OUT_DIR / name
    out.mkdir(parents=True, exist_ok=True)
    (out / f"{name}.svg").write_text(d.as_svg(), encoding="utf-8")
    cairosvg.svg2png(bytestring=d.as_svg().encode(), write_to=str(out / f"{name}.png"), scale=300 / 72.0)
    print(f"wrote {out}/{name}.svg / .png")


if __name__ == "__main__":
    render()
