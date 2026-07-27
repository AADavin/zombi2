"""Trait-evolution examples: a trait evolving down the tree, branches coloured by its value."""

from __future__ import annotations

import helpers as h
from helpers import Example

import phylustrator as ph
from zombi2.species import simulate_species_tree
from zombi2.traits import simulate_continuous, simulate_discrete

SEED = 42
N = 100

# two binary characters X, Y as compound Mk states; only single-bit flips are named (all else 0).
# the dependence: gaining Y is slow while X is absent (00->01) but fast once X is present (10->11).
_DEP_STATES = ("00", "01", "10", "11")
_DEP_RATES = {"00->01": 0.05, "01->00": 2.0, "10->11": 2.0, "11->10": 0.05,   # flip Y (rate ~ X)
              "00->10": 0.5, "10->00": 0.5, "01->11": 0.5, "11->01": 0.5}      # flip X
# one binary palette per trait (shared "absent"); X and Y are shown as two lanes on each branch.
_ABSENT = "#c2cac8"
_XPAL = {"0": _ABSENT, "1": "#3C8D6E"}      # X present = green
_YPAL = {"0": _ABSENT, "1": "#8B6B9E"}      # Y present = muted purple


def _tree():
    sp = simulate_species_tree(birth=1.0, n_extant=N, seed=SEED)
    return sp.complete_tree, ph.trees.loads(sp.complete_tree.to_newick())


def _style():
    return ph.Style(width=1200, height=1000, margin=80, branch_width=1.4)


def brownian_motion(out):
    ct, tree = _tree()
    res = simulate_continuous(ct, start=0.0, rate=1.0, seed=SEED)
    (ph.trees.plot(tree, style=_style())
     + ph.trees.color_branches(h.node_values(res))
     + ph.trees.colorbar("trait value", width=240, height=16, size=20)
     + ph.trees.time_axis("time", tick_size=22, label_size=28)).save(out)


def ornstein_uhlenbeck(out):
    ct, tree = _tree()
    res = simulate_continuous(ct, start=6.0, rate=1.0, reverts_to=0.0, pull=2.5, seed=SEED)
    (ph.trees.plot(tree, style=_style())
     + ph.trees.color_branches(h.node_values(res))
     + ph.trees.colorbar("trait value  —  OU: optimum 0, pull 2.5, σ 1, start 6",
                   width=240, height=16, size=20)
     + ph.trees.time_axis("time", tick_size=22, label_size=28)).save(out)


def discrete_states(out):
    ct, tree = _tree()
    res = simulate_discrete(ct, states=("marine", "terrestrial"), switch=0.5, start="marine", seed=5)
    history = {f"n{i}": segs for i, segs in res.history.items()}   # per-branch (state, duration)
    (ph.trees.plot(tree, style=_style())
     + ph.trees.color_history(history, palette={"marine": "#4C9AA6", "terrestrial": "#E08A3C"})
     + ph.trees.legend("habitat")
     + ph.trees.time_axis("time", tick_size=22, label_size=28)).save(out)


def correlated(out):
    ct, _ = _tree()
    res = simulate_continuous(ct, start={"x": 0.0, "y": 0.0}, rate={"x": 1.0, "y": 1.0},
                              correlation={("x", "y"): 0.9}, seed=SEED)
    vx = {f"n{i}": v["x"] for i, v in res.node_values.items()}
    vy = {f"n{i}": v["y"] for i, v in res.node_values.items()}
    px, py = out.replace(".png", "_x.png"), out.replace(".png", "_y.png")
    for png, vals in ((px, vx), (py, vy)):
        h.render_tree_for_composite(ph.trees.loads(ct.to_newick()), png, vals,
                                    width=1000, height=520, branch_width=1.2)
    tips = list(ct.extant())
    xs = [res.node_values[n.id]["x"] for n in tips]
    ys = [res.node_values[n.id]["y"] for n in tips]
    h.composite_two_trees_scatter(px, py, xs, ys, out)


def dependent_characters(out):
    sp = simulate_species_tree(birth=1.0, n_extant=35, seed=7)             # Yule, all lineages extant
    ct = sp.complete_tree
    tree = ph.trees.loads(ct.to_newick())
    res = simulate_discrete(ct, states=_DEP_STATES, switch=_DEP_RATES, start="00", seed=3)
    raw = res.history                                                      # {id: [(compound, dur), …]}
    xh = {f"n{i}": [(s[0], d) for s, d in segs] for i, segs in raw.items()}   # project onto X …
    yh = {f"n{i}": [(s[1], d) for s, d in segs] for i, segs in raw.items()}   # … and onto Y
    tips = list(ct.extant())
    M = ph.genomes.Matrix(rows=[f"n{n.id}" for n in tips], cols=["X", "Y"],
                          values=[list(res.node_values[n.id]) for n in tips])   # "01" -> ["0","1"]
    realization = out.replace(".png", "_real.png")
    # no base skeleton — the lanes paint each branch as a two-tone band and carry the coloured joints
    fig = (ph.trees.plot(tree, style=ph.Style(width=820, height=1000, margin=62, branch_width=6.5),
                         skeleton=False)
           + ph.trees.color_lanes([(xh, _XPAL), (yh, _YPAL)], gap=1.0)
           + ph.trees.time_axis("time", tick_size=20, label_size=26))
    ph.beside(fig, ph.genomes.states(M, col_palettes=[_XPAL, _YPAL], legend=False),
              width=1080, tree_fraction=0.82, footer=24).save(realization)
    h.composite_model_realization(realization, out,
                                  lambda ax: h.draw_grid_markov(ax, _DEP_RATES, _XPAL, _YPAL))


_C_BM = '''\
### simulate  —  Brownian motion down a 100-tip tree
from zombi2.species import simulate_species_tree
from zombi2.traits import simulate_continuous

sp = simulate_species_tree(birth=1.0, n_extant=100, seed=42)
res = simulate_continuous(sp.complete_tree, start=0.0, rate=1.0, seed=42)

### plot  —  branches coloured by trait value, with a colour key
import phylustrator as ph

tree = ph.trees.loads(sp.complete_tree.to_newick())
values = {f"n{i}": v for i, v in res.node_values.items()}
(ph.trees.plot(tree)
 + ph.trees.color_branches(values)
 + ph.trees.colorbar("trait value", width=240, height=16, size=20)
 + ph.trees.time_axis("time", tick_size=22, label_size=28)).save("bm.png")'''

_C_OU = '''\
### simulate  —  Ornstein–Uhlenbeck (pulled to an optimum), 100 tips
from zombi2.species import simulate_species_tree
from zombi2.traits import simulate_continuous

sp = simulate_species_tree(birth=1.0, n_extant=100, seed=42)
res = simulate_continuous(sp.complete_tree, start=6.0, rate=1.0, reverts_to=0.0, pull=2.5, seed=42)

### plot  —  the colour key names the OU parameters
import phylustrator as ph

tree = ph.trees.loads(sp.complete_tree.to_newick())
values = {f"n{i}": v for i, v in res.node_values.items()}
(ph.trees.plot(tree)
 + ph.trees.color_branches(values)
 + ph.trees.colorbar("trait value  —  OU: optimum 0, pull 2.5, σ 1, start 6", width=240, height=16, size=20)
 + ph.trees.time_axis("time", tick_size=22, label_size=28)).save("ou.png")'''

_C_DISCRETE = '''\
### simulate  —  a two-state trait hopping between habitats, 100 tips
from zombi2.species import simulate_species_tree
from zombi2.traits import simulate_discrete

sp = simulate_species_tree(birth=1.0, n_extant=100, seed=42)
res = simulate_discrete(sp.complete_tree, states=("marine", "terrestrial"),
                        switch=0.5, start="marine", seed=5)

### plot  —  color_history paints each branch as segments (states change *along* a branch)
import phylustrator as ph

tree = ph.trees.loads(sp.complete_tree.to_newick())
history = {f"n{i}": segs for i, segs in res.history.items()}   # [(state, duration), ...] per branch
(ph.trees.plot(tree)
 + ph.trees.color_history(history, palette={"marine": "#4C9AA6", "terrestrial": "#E08A3C"})
 + ph.trees.legend("habitat")
 + ph.trees.time_axis("time", tick_size=22, label_size=28)).save("discrete.png")'''

_C_CORRELATED = '''\
### simulate  —  two traits evolving together (correlation 0.9)
from zombi2.species import simulate_species_tree
from zombi2.traits import simulate_continuous

sp = simulate_species_tree(birth=1.0, n_extant=100, seed=42)
ct = sp.complete_tree
res = simulate_continuous(ct, start={"x": 0.0, "y": 0.0}, rate={"x": 1.0, "y": 1.0},
                          correlation={("x", "y"): 0.9}, seed=42)

### plot  —  two trees (coloured by x, then y) stacked left; a tip x-vs-y scatter right
import phylustrator as ph
import helpers as h                          # gallery helper: two trees + scatter

vx = {f"n{i}": v["x"] for i, v in res.node_values.items()}
vy = {f"n{i}": v["y"] for i, v in res.node_values.items()}
h.render_tree_for_composite(ph.trees.loads(ct.to_newick()), "tree_x.png", vx)
h.render_tree_for_composite(ph.trees.loads(ct.to_newick()), "tree_y.png", vy)
tips = list(ct.extant())
xs = [res.node_values[n.id]["x"] for n in tips]
ys = [res.node_values[n.id]["y"] for n in tips]
h.composite_two_trees_scatter("tree_x.png", "tree_y.png", xs, ys, "correlated.png")'''


_C_DEPENDENT = '''\
### simulate  —  two binary characters X, Y where Y's flip rate depends on X's state
from zombi2.species import simulate_species_tree
from zombi2.traits import simulate_discrete

# compound Mk states 00/01/10/11; only single-bit flips are named (everything else stays 0).
# gaining Y is slow while X is absent (00->01) but fast once X is present (10->11): the dependence.
rates = {"00->01": 0.05, "01->00": 2.0, "10->11": 2.0, "11->10": 0.05,   # flip Y (rate ~ X)
         "00->10": 0.5,  "10->00": 0.5, "01->11": 0.5, "11->01": 0.5}     # flip X
sp = simulate_species_tree(birth=1.0, n_extant=35, seed=7)
ct = sp.complete_tree
res = simulate_discrete(ct, states=("00", "01", "10", "11"), switch=rates, start="00", seed=3)

### plot  —  the two traits as two lanes on each branch, beside a matching tip matrix
import phylustrator as ph

xpal = {"0": "#c2cac8", "1": "#3C8D6E"}       # X: absent / present (green)
ypal = {"0": "#c2cac8", "1": "#8B6B9E"}       # Y: absent / present (muted purple)
tree = ph.trees.loads(ct.to_newick())
raw = res.history                             # {id: [(compound_state, dur), …]}
xh = {f"n{i}": [(s[0], d) for s, d in segs] for i, segs in raw.items()}   # project onto X
yh = {f"n{i}": [(s[1], d) for s, d in segs] for i, segs in raw.items()}   # project onto Y
tips = list(ct.extant())
M = ph.genomes.Matrix(rows=[f"n{n.id}" for n in tips], cols=["X", "Y"],
                      values=[list(res.node_values[n.id]) for n in tips])   # "01" -> ["0","1"]
# no base skeleton: color_lanes paints each branch as a two-tone band (X | Y) and the coloured joints
fig = (ph.trees.plot(tree, skeleton=False)
       + ph.trees.color_lanes([(xh, xpal), (yh, ypal)], gap=1.0))
ph.beside(fig, ph.genomes.states(M, col_palettes=[xpal, ypal], legend=False)).save("dependent.png")'''


EXAMPLES = [
    Example("bm", "Brownian motion", "Free diffusion — sister lineages drift apart with time.",
            "continuous", brownian_motion, code=_C_BM),
    Example("ou", "Ornstein–Uhlenbeck", "Pulled to an optimum: a high start (yellow) converges to blue.",
            "continuous", ornstein_uhlenbeck, code=_C_OU),
    Example("discrete", "Discrete states",
            "A two-state trait hops between habitats; each branch is painted by its state history.",
            "discrete · Mk", discrete_states, code=_C_DISCRETE),
    Example("correlated", "Dependent continuous traits",
            "Two traits evolve together (r&nbsp;=&nbsp;0.9) — two trees, coloured by each trait, and the "
            "tip scatter.",
            "continuous · +&nbsp;scatter", correlated, code=_C_CORRELATED),
    Example("dependent", "Dependent discrete traits",
            "Two binary characters where one's flip rate depends on the other's state — the model as a "
            "2×2 chain (arrow width = rate), and a run: the tree painted by compound state, beside the "
            "presence/absence tips. <code>simulate_discrete(states=(&quot;00&quot;,…),&nbsp;switch={…})</code>.",
            "discrete · dependent", dependent_characters, code=_C_DEPENDENT),
]
