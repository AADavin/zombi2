"""Trait-evolution examples: a trait evolving down the tree, branches coloured by its value."""

from __future__ import annotations

import helpers as h
from helpers import Example

import phylustrator as ph
from zombi2.rates import modifiers as mod
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
# the driven example: a discrete habitat (the driver) and the body size it speeds up (the target)
_HABITAT = {"stable": "#3A7CA5", "fluctuating": "#E4572E"}
_SWITCH, _BASE, _FACTOR = 0.4, 0.25, 20.0


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
    tips = list(ct.extant_leaves())
    xs = [res.node_values[n.id]["x"] for n in tips]
    ys = [res.node_values[n.id]["y"] for n in tips]
    h.composite_two_trees_scatter(px, py, xs, ys, out)


def dependent_characters(out):
    sp = simulate_species_tree(birth=1.0, n_extant=55, seed=7)             # Yule, all lineages extant
    ct = sp.complete_tree
    res = simulate_discrete(ct, states=_DEP_STATES, switch=_DEP_RATES, start="00", seed=3)
    raw = res.history                                                      # {id: [(compound, dur), …]}
    xh = {f"n{i}": [(s[0], d) for s, d in segs] for i, segs in raw.items()}   # project onto X …
    yh = {f"n{i}": [(s[1], d) for s, d in segs] for i, segs in raw.items()}   # … and onto Y
    # two trees, coloured by each character's history (like the continuous "correlated" figure) — you
    # can see Y (purple) tends to be present where X (green) is; the model panel says why.
    px, py = out.replace(".png", "_x.png"), out.replace(".png", "_y.png")
    for png, hist, pal in ((px, xh, _XPAL), (py, yh, _YPAL)):
        (ph.trees.plot(ph.trees.loads(ct.to_newick()),
                       style=ph.Style(width=1000, height=520, margin=45, branch_width=2.6),
                       skeleton=False)
         + ph.trees.color_history(hist, palette=pal)).save(png)
    h.composite_two_trees_panel(px, py, lambda ax: h.draw_grid_markov(ax, _DEP_RATES, _XPAL, _YPAL),
                                out, x_label="character X", y_label="character Y")


def driven_trait(out):
    """A trait driving a trait, with the same ``DrivenBy`` a trait uses to drive a genome rate. The
    habitat is grown first and held fixed; body size then diffuses down the same tree, twenty times
    faster on the stretches of branch where the habitat fluctuates.

    The figure is the tree painted by the habitat, with the driver·mapping·target diagram above it.
    What the driving *did* to body size is the Conditioning section's "One trait drives another",
    which paints the same tree twice; repeating it here only made this figure taller."""
    ct = simulate_species_tree(birth=1.0, n_extant=50, seed=7).complete_tree
    hab = simulate_discrete(ct, states=["stable", "fluctuating"], switch=_SWITCH, start="stable", seed=5)
    simulate_continuous(ct, start=0.0, seed=9,           # the target the diagram names
            rate=_BASE * mod.DrivenBy(hab, {"fluctuating": _FACTOR, "stable": 1.0}))
    lab = ct.labels()                      # {id: 'n<id>'} — the tree's own names, never built by hand
    tree_png = out.replace(".png", "_tree.png")
    (ph.trees.plot(ph.trees.loads(ct.to_newick()), skeleton=False,
                   # margin clears the time axis and its label: tighter and "time" is cut in half
                   style=ph.Style(width=1320, height=640, margin=96, branch_width=2.6))
     + ph.trees.color_history({lab[i]: segs for i, segs in hab.history.items()},
                              palette=_HABITAT)
     + ph.trees.time_axis("time", tick_size=20, label_size=26, bold=False)).save(tree_png)
    diag = h.conditioning_png(out.replace(".png", "_diag.png"), driver="habitat",
        states=["stable", "fluctuating"],
        switch={"stable->fluctuating": _SWITCH, "fluctuating->stable": _SWITCH},
        mapping={"fluctuating": _FACTOR, "stable": 1}, target="body size", target_base=_BASE,
        target_sub="variance-rate", state_colors=_HABITAT)
    h.composite_under_diagram(out, diag, [(tree_png, "habitat")])


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
tips = list(ct.extant_leaves())
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

### plot  —  two trees (coloured by each character's history) + the 2x2 model chain, like "correlated"
import phylustrator as ph
import helpers as h                           # gallery helper: two trees + a side panel

xpal = {"0": "#c2cac8", "1": "#3C8D6E"}       # X: absent / present (green)
ypal = {"0": "#c2cac8", "1": "#8B6B9E"}       # Y: absent / present (muted purple)
raw = res.history                             # {id: [(compound_state, dur), …]}
xh = {f"n{i}": [(s[0], d) for s, d in segs] for i, segs in raw.items()}   # project onto X
yh = {f"n{i}": [(s[1], d) for s, d in segs] for i, segs in raw.items()}   # project onto Y
for png, hist, pal in (("tree_x.png", xh, xpal), ("tree_y.png", yh, ypal)):
    (ph.trees.plot(ph.trees.loads(ct.to_newick()), skeleton=False)
     + ph.trees.color_history(hist, palette=pal)).save(png)
# right panel = the 2x2 CTMC (the model), kept as the key; Y present tends to track X present
h.composite_two_trees_panel("tree_x.png", "tree_y.png",
                            lambda ax: h.draw_grid_markov(ax, rates, xpal, ypal),
                            "dependent.png", x_label="character X", y_label="character Y")'''


_C_DRIVEN = '''\
### simulate  —  a discrete trait drives a continuous one: the driver first, then the target
from zombi2.species import simulate_species_tree
from zombi2.traits import simulate_continuous, simulate_discrete
from zombi2.rates import modifiers as mod

ct = simulate_species_tree(birth=1.0, n_extant=50, seed=7).complete_tree
hab = simulate_discrete(ct, states=["stable", "fluctuating"], switch=0.4, start="stable", seed=5)
# the same DrivenBy that lets a trait drive a genome rate lets it drive another trait: body size
# diffuses 20x faster on the stretches of branch where the habitat fluctuates
size = simulate_continuous(ct, start=0.0, seed=9,
        rate=0.25 * mod.DrivenBy(hab, {"fluctuating": 20.0, "stable": 1.0}))

### plot  —  the tree coloured by the habitat, the driver that sets body size's diffusion rate
import phylustrator as ph

pal = {"stable": "#3A7CA5", "fluctuating": "#E4572E"}
lab = ct.labels()                                # {id: 'n<id>'}
(ph.trees.plot(ph.trees.loads(ct.to_newick()), skeleton=False)
 + ph.trees.color_history({lab[i]: segs for i, segs in hab.history.items()}, palette=pal)
 + ph.trees.time_axis("time", bold=False)).save("tree.png")
# the figure then composites the driver->mapping->target diagram (habitat -> body size) on top'''


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
            "Two binary characters where one's flip rate depends on the other's state. X in green, Y in "
            "purple, so you can see Y is present where X is; the 2×2 chain is the model.",
            "discrete · dependent", dependent_characters, code=_C_DEPENDENT),
    Example("driven", "A trait driving a trait",
            "A habitat trait sets how fast body size diffuses: twenty times faster where the habitat "
            "fluctuates. <code>rate&nbsp;=&nbsp;0.25&nbsp;*&nbsp;mod.DrivenBy(habitat,&nbsp;{…})</code>.",
            "trait → trait", driven_trait, code=_C_DRIVEN),
]
