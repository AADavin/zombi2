"""Joining examples: a trait drives diversification (BiSSE / MuSSE / SSE), so the tree's *shape*
reflects the state — the tree and its driver grow together. Each figure carries the state Markov
chain (the model) as an inset, bottom-left."""

from __future__ import annotations

import helpers as h
from helpers import Example

import phylustrator as ph
from zombi2 import joint, traits
from zombi2.rates import modifiers as mod

_BISSE = {"fast": "#E4572E", "slow": "#3A7CA5"}
_SSE = {"doomed": "#B0413E", "safe": "#2A9D8F"}
_MUSSE = {"slow": "#3A7CA5", "medium": "#F2A541", "fast": "#E4572E"}


def _style():
    return ph.Style(width=1300, height=1000, margin=82, branch_width=3.4)


def _history(r):
    return {f"n{i}": segs for i, segs in r.trait.history.items()}    # per-branch (state, duration)


def bisse(out):
    r = joint.simulate_joint(
        birth=1.0 * mod.DrivenBy("trait", {"fast": 2.6, "slow": 0.7}),
        trait=traits.discrete(states=["fast", "slow"], switch=0.35),
        n_extant=70, seed=3)
    tree_png = out.replace(".png", "_tree.png")
    (ph.trees.plot(ph.trees.loads(r.complete_tree.to_newick()), style=_style(), skeleton=False)
     + ph.trees.color_history(_history(r), palette=_BISSE)
     + ph.trees.time_axis("time", tick_size=22, label_size=28)).save(tree_png)
    h.composite_markov(tree_png, out, lambda ax: h.draw_markov(
        ax, ["fast", "slow"], _BISSE, 0.35, {"fast": 2.6, "slow": 0.7}, symbol="λ"),
        loc=(0.02, 0.04, 0.34, 0.30))


def state_extinction(out):
    r = joint.simulate_joint(
        birth=1.0,
        death=1.0 * mod.DrivenBy("trait", {"doomed": 0.75, "safe": 0.05}),
        trait=traits.discrete(states=["doomed", "safe"], switch=0.3),
        n_extant=35, seed=1)
    ct = r.complete_tree
    tree = ph.trees.loads(ct.to_newick())
    dashed = h.dashed_extinct(tree, {f"n{n.id}" for n in ct.extinct()})
    tree_png = out.replace(".png", "_tree.png")
    (ph.trees.plot(tree, style=_style(), skeleton=False)
     + ph.trees.color_history(_history(r), palette=_SSE, dashed=dashed)     # dead lineages: dashed + coloured
     + ph.trees.time_axis("time", tick_size=22, label_size=28)).save(tree_png)
    h.composite_markov(tree_png, out, lambda ax: h.draw_markov(
        ax, ["doomed", "safe"], _SSE, 0.3, {"doomed": 0.75, "safe": 0.05}, symbol="μ"),
        loc=(0.02, 0.05, 0.34, 0.30))


def musse(out):
    r = joint.simulate_joint(
        birth=1.0 * mod.DrivenBy("trait", {"slow": 0.6, "medium": 1.3, "fast": 2.6}),
        death=0.4,                                                    # so some lineages die (dashed)
        trait=traits.discrete(states=["slow", "medium", "fast"], switch=0.3),
        n_extant=50, seed=2)
    ct = r.complete_tree
    tree = ph.trees.loads(ct.to_newick())
    dashed = h.dashed_extinct(tree, {f"n{n.id}" for n in ct.extinct()})
    tree_png = out.replace(".png", "_tree.png")
    (ph.trees.plot(tree, style=_style(), skeleton=False)
     + ph.trees.color_history(_history(r), palette=_MUSSE, dashed=dashed)
     + ph.trees.time_axis("time", tick_size=22, label_size=28)).save(tree_png)
    h.composite_markov(tree_png, out, lambda ax: h.draw_markov(
        ax, ["slow", "medium", "fast"], _MUSSE, 0.3, {"slow": 0.6, "medium": 1.3, "fast": 2.6},
        symbol="λ"), loc=(0.02, 0.075, 0.35, 0.40))


_C_BISSE = '''\
### simulate  —  a 2-state trait drives speciation (BiSSE)
from zombi2 import joint, traits
from zombi2.rates import modifiers as mod

r = joint.simulate_joint(
    birth=1.0 * mod.DrivenBy("trait", {"fast": 2.6, "slow": 0.7}),
    trait=traits.discrete(states=["fast", "slow"], switch=0.35),
    n_extant=70, seed=3)

### plot  —  color_history paints each branch by its state history; a Markov-chain inset shows the model
import phylustrator as ph
import helpers as h

palette = {"fast": "#E4572E", "slow": "#3A7CA5"}
tree = ph.trees.loads(r.complete_tree.to_newick())
history = {f"n{i}": segs for i, segs in r.trait.history.items()}
(ph.trees.plot(tree, skeleton=False)
 + ph.trees.color_history(history, palette=palette)
 + ph.trees.time_axis("time", tick_size=22, label_size=28)).save("tree.png")
h.composite_markov("tree.png", "bisse.png", lambda ax: h.draw_markov(
    ax, ["fast", "slow"], palette, 0.35, {"fast": 2.6, "slow": 0.7}, symbol="λ"))'''

_C_STATE = '''\
### simulate  —  one state dies far faster (state-dependent extinction)
from zombi2 import joint, traits
from zombi2.rates import modifiers as mod

r = joint.simulate_joint(
    birth=1.0,
    death=1.0 * mod.DrivenBy("trait", {"doomed": 0.75, "safe": 0.05}),
    trait=traits.discrete(states=["doomed", "safe"], switch=0.3),
    n_extant=35, seed=1)
ct = r.complete_tree

### plot  —  coloured by state history, extinct lineages dashed, model inset
import phylustrator as ph
import helpers as h

palette = {"doomed": "#B0413E", "safe": "#2A9D8F"}
tree = ph.trees.loads(ct.to_newick())
dashed = {f"n{n.id}" for n in ct.extinct()}                 # (+ their all-extinct ancestors)
history = {f"n{i}": segs for i, segs in r.trait.history.items()}
(ph.trees.plot(tree, skeleton=False)
 + ph.trees.color_history(history, palette=palette, dashed=dashed)
 + ph.trees.time_axis("time", tick_size=22, label_size=28)).save("tree.png")
h.composite_markov("tree.png", "sse.png", lambda ax: h.draw_markov(
    ax, ["doomed", "safe"], palette, 0.3, {"doomed": 0.75, "safe": 0.05}, symbol="μ"))'''

_C_MUSSE = '''\
### simulate  —  three graded speciation rates + constant death (MuSSE)
from zombi2 import joint, traits
from zombi2.rates import modifiers as mod

r = joint.simulate_joint(
    birth=1.0 * mod.DrivenBy("trait", {"slow": 0.6, "medium": 1.3, "fast": 2.6}),
    death=0.4,
    trait=traits.discrete(states=["slow", "medium", "fast"], switch=0.3),
    n_extant=50, seed=2)
ct = r.complete_tree

### plot  —  state history + dashed extinct lineages + a 3-state Markov chain
import phylustrator as ph
import helpers as h

palette = {"slow": "#3A7CA5", "medium": "#F2A541", "fast": "#E4572E"}
tree = ph.trees.loads(ct.to_newick())
dashed = {f"n{n.id}" for n in ct.extinct()}
history = {f"n{i}": segs for i, segs in r.trait.history.items()}
(ph.trees.plot(tree, skeleton=False)
 + ph.trees.color_history(history, palette=palette, dashed=dashed)
 + ph.trees.time_axis("time", tick_size=22, label_size=28)).save("tree.png")
h.composite_markov("tree.png", "musse.png", lambda ax: h.draw_markov(
    ax, ["slow", "medium", "fast"], palette, 0.3, {"slow": 0.6, "medium": 1.3, "fast": 2.6},
    symbol="λ"))'''


EXAMPLES = [
    Example("bisse", "BiSSE",
            "A two-state trait drives speciation — the fast state's clades take over; the inset is the "
            "state Markov chain.",
            "trait → speciation", bisse, code=_C_BISSE),
    Example("state_extinction", "State-dependent extinction",
            "One state dies far faster; the doomed lineages (dashed) drop out.",
            "trait → extinction", state_extinction, code=_C_STATE),
    Example("musse", "MuSSE",
            "Three graded speciation rates with constant death — the fastest state fills the tree, "
            "extinct lineages dashed.",
            "trait → speciation", musse, code=_C_MUSSE),
]
