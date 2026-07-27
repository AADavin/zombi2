"""Joining and conditioning: one level drives another through the same ``DrivenBy`` mechanism. A trait
drives diversification (BiSSE / MuSSE / SSE), so the tree's *shape* reflects the state; or a trait
conditions the genome, driving gene loss so genome *size* reflects the state. The SSE figures carry the
state Markov chain (the model) as an inset, bottom-left."""

from __future__ import annotations

import helpers as h
from helpers import Example

import phylustrator as ph
from zombi2 import joint, traits
from zombi2.genomes import simulate_genomes_family
from zombi2.rates import modifiers as mod
from zombi2.species import simulate_species_tree
from zombi2.traits import simulate_discrete

_BISSE = {"fast": "#E4572E", "slow": "#3A7CA5"}
_SSE = {"doomed": "#B0413E", "safe": "#2A9D8F"}
_MUSSE = {"slow": "#3A7CA5", "medium": "#F2A541", "fast": "#E4572E"}
_HAB = {"free-living": "#2E8B6F", "endosymbiont": "#C25A3C"}


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


def genome_reduction(out):
    sp = simulate_species_tree(birth=1.0, n_extant=36, seed=4)
    ct = sp.complete_tree
    # an irreversible lifestyle: once a lineage turns endosymbiont it stays (Dollo-ish)
    hab = simulate_discrete(ct, states=["free-living", "endosymbiont"], start="free-living", seed=8,
                            switch={"free-living->endosymbiont": 0.09, "endosymbiont->free-living": 0.0})
    # the trait CONDITIONS the genome: endosymbionts shed genes fast and gain almost none
    g = simulate_genomes_family(ct, initial_families=55, duplication=0.1,
            origination=3.0 * mod.DrivenBy(hab, {"endosymbiont": 0.1, "free-living": 1.0}),
            loss=0.08 * mod.DrivenBy(hab, {"endosymbiont": 12.0, "free-living": 1.0}), seed=9)
    tree = ph.trees.loads(ct.to_newick())
    history = {f"n{i}": segs for i, segs in hab.history.items()}
    tips = list(ct.extant())
    sizes = {f"n{n.id}": len(g.genomes[n.id]) for n in tips}
    tipcol = {f"n{n.id}": _HAB[hab.values[n.id]] for n in tips}
    fig = (ph.trees.plot(tree, style=ph.Style(width=900, height=1000, margin=98, branch_width=3.0),
                         skeleton=False)
           + ph.trees.color_history(history, palette=_HAB)
           + ph.trees.legend("lifestyle", size=25)
           + ph.trees.time_axis("time", tick_size=22, label_size=28, bold=False))
    ph.beside(fig, ph.genomes.bars(sizes, colors=tipcol, label="genome size (genes)",
                                   tick_size=22, label_size=28),
              width=1150, tree_fraction=0.58, footer=36).save(out)


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

_C_REDUCTION = '''\
### simulate  —  a trait CONDITIONS the genome: grow the tree + trait, then the genome
from zombi2.species import simulate_species_tree
from zombi2.traits import simulate_discrete
from zombi2.genomes import simulate_genomes_family
from zombi2.rates import modifiers as mod

sp = simulate_species_tree(birth=1.0, n_extant=36, seed=4)
ct = sp.complete_tree
# an irreversible lifestyle: free-living -> endosymbiont, never back
hab = simulate_discrete(ct, states=["free-living", "endosymbiont"], start="free-living", seed=8,
                        switch={"free-living->endosymbiont": 0.09, "endosymbiont->free-living": 0.0})
# the SAME DrivenBy that couples a trait to speciation couples it to the genome: endosymbionts
# shed genes fast (loss x12) and gain almost none (origination x0.1)
g = simulate_genomes_family(ct, initial_families=55, duplication=0.1,
        origination=3.0 * mod.DrivenBy(hab, {"endosymbiont": 0.1, "free-living": 1.0}),
        loss=0.08 * mod.DrivenBy(hab, {"endosymbiont": 12.0, "free-living": 1.0}), seed=9)

### plot  —  tree coloured by lifestyle, beside per-tip genome-size bars (aligned axes)
import phylustrator as ph

pal = {"free-living": "#2E8B6F", "endosymbiont": "#C25A3C"}
tree = ph.trees.loads(ct.to_newick())
history = {f"n{i}": segs for i, segs in hab.history.items()}
tips = list(ct.extant())
sizes  = {f"n{n.id}": len(g.genomes[n.id]) for n in tips}          # gene count per tip
colors = {f"n{n.id}": pal[hab.values[n.id]] for n in tips}         # bar colour = lifestyle
fig = (ph.trees.plot(tree, skeleton=False)
       + ph.trees.color_history(history, palette=pal)
       + ph.trees.legend("lifestyle")
       + ph.trees.time_axis("time", bold=False))
ph.beside(fig, ph.genomes.bars(sizes, colors=colors, label="genome size (genes)")).save("reduction.png")'''


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
    Example("genome_reduction", "Genome reduction",
            "The same coupling, aimed at the genome instead of the tree: an irreversible endosymbiont "
            "lifestyle drives fast gene loss and near-zero gene gain, so those lineages' genomes "
            "collapse. Tree coloured by lifestyle; bars are per-tip genome size.",
            "trait → genome", genome_reduction, code=_C_REDUCTION),
]
