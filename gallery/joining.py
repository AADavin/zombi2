"""Two ways one level reads another, through the same ``DrivenBy`` mechanism.

**Conditioning** grows the driver first and holds it fixed: a trait is simulated on the tree, then a
genome run reads it, so the trait's state sets a genome rate. **Joining** grows both at once, because
the trait drives speciation or extinction and so shapes the tree it is evolving on
(``joint.simulate_joint``). The two lists below feed two gallery sections.
"""

from __future__ import annotations

import matplotlib.image as mpimg
import matplotlib.pyplot as plt
from matplotlib import cm, colors

import helpers as h
from helpers import Example

import phylustrator as ph
from zombi2 import joint, traits
from zombi2.genomes import simulate_genomes_family
from zombi2.rates import modifiers as mod
from zombi2.rates.mapping import Curve
from zombi2.species import simulate_species_tree
from zombi2.traits import simulate_continuous, simulate_discrete

_BISSE = {"fast": "#E4572E", "slow": "#3A7CA5"}
_SSE = {"doomed": "#B0413E", "safe": "#2A9D8F"}
_MUSSE = {"slow": "#3A7CA5", "medium": "#F2A541", "fast": "#E4572E"}
_HAB = {"free-living": "#2E8B6F", "endosymbiont": "#C25A3C"}
_SEL = {"purifying": "#3A7CA5", "relaxed": "#C25A3C"}       # relaxed selection → duplicates accumulate
_COMP = {"quiet": "#8f99a3", "competent": "#2E8B6F"}        # competent → takes up more DNA


def _style():
    return ph.Style(width=1300, height=1000, margin=82, branch_width=3.4)


def _state_history(ct, trait):
    lab = ct.labels()                   # 'n<id>', or 'e<id>' for a lineage that went extinct
    return {lab[i]: segs for i, segs in trait.history.items()}     # per-branch (state, duration)


def _history(r):
    return _state_history(r.complete_tree, r.trait)


# --- joining: the trait and the tree grow together ------------------------------------------

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
        ax, ["fast", "slow"], _BISSE, {"fast": 2.6, "slow": 0.7}, symbol="λ"),
        loc=(0.02, 0.04, 0.34, 0.30))


def state_extinction(out):
    r = joint.simulate_joint(
        birth=1.0,
        death=1.0 * mod.DrivenBy("trait", {"doomed": 0.75, "safe": 0.05}),
        trait=traits.discrete(states=["doomed", "safe"], switch=0.3),
        n_extant=35, seed=1)
    ct = r.complete_tree
    tree = ph.trees.loads(ct.to_newick())
    dashed = h.dashed_extinct(tree, ct)
    tree_png = out.replace(".png", "_tree.png")
    (ph.trees.plot(tree, style=_style(), skeleton=False)
     + ph.trees.color_history(_history(r), palette=_SSE, dashed=dashed)     # dead lineages: dashed + coloured
     + ph.trees.time_axis("time", tick_size=22, label_size=28)).save(tree_png)
    h.composite_markov(tree_png, out, lambda ax: h.draw_markov(
        ax, ["doomed", "safe"], _SSE, {"doomed": 0.75, "safe": 0.05}, symbol="μ"),
        loc=(0.02, 0.05, 0.34, 0.30))


def musse(out):
    r = joint.simulate_joint(
        birth=1.0 * mod.DrivenBy("trait", {"slow": 0.6, "medium": 1.3, "fast": 2.6}),
        death=0.4,                                                    # so some lineages die (dashed)
        trait=traits.discrete(states=["slow", "medium", "fast"], switch=0.3),
        n_extant=50, seed=2)
    ct = r.complete_tree
    tree = ph.trees.loads(ct.to_newick())
    dashed = h.dashed_extinct(tree, ct)
    tree_png = out.replace(".png", "_tree.png")
    (ph.trees.plot(tree, style=_style(), skeleton=False)
     + ph.trees.color_history(_history(r), palette=_MUSSE, dashed=dashed)
     + ph.trees.time_axis("time", tick_size=22, label_size=28)).save(tree_png)
    h.composite_markov(tree_png, out, lambda ax: h.draw_markov(
        ax, ["slow", "medium", "fast"], _MUSSE, {"slow": 0.6, "medium": 1.3, "fast": 2.6},
        symbol="λ"), loc=(0.02, 0.075, 0.35, 0.40))


# --- conditioning: the driver is grown first, then the genome reads it -----------------------

def _conditioned_genome(out, ct, layers, sizes, tipcol, diagram):
    """The shared conditioning-figure layout: the tree painted by the driver, beside per-tip genome-size
    bars, with the driver·modifier·target diagram small on top. ``layers`` are the Phylustrator layers
    that colour the tree; ``diagram`` is the kwargs for :func:`helpers.conditioning_png`."""
    fig = ph.trees.plot(ph.trees.loads(ct.to_newick()), skeleton=False,
                        style=ph.Style(width=900, height=900, margin=92, branch_width=3.0))
    for layer in layers:                      # no legend on the tree — the diagram is the key
        fig = fig + layer
    fig = fig + ph.trees.time_axis("time", tick_size=20, label_size=26, bold=False)
    real = out.replace(".png", "_real.png")
    ph.beside(fig, ph.genomes.bars(sizes, colors=tipcol, label="genome size (genes)",
                                   tick_size=20, label_size=26),
              width=1150, tree_fraction=0.58, footer=36).save(real)
    diag = h.conditioning_png(out.replace(".png", "_diag.png"), **diagram)
    fig2 = plt.figure(figsize=(12, 9.6))
    axr = fig2.add_axes([0.0, 0.0, 1.0, 0.80])
    axr.imshow(mpimg.imread(real))
    axr.set_axis_off()
    axd = fig2.add_axes([0.30, 0.80, 0.40, 0.185])
    axd.imshow(mpimg.imread(diag))
    axd.set_axis_off()
    fig2.savefig(out, dpi=140, bbox_inches="tight")
    plt.close(fig2)


def _sizes(ct, g, trait, palette):
    lab = ct.labels()
    tips = list(ct.extant_leaves())
    sizes = {lab[n.id]: len(g.genomes[n.id]) for n in tips}
    tipcol = {lab[n.id]: palette[trait.values[n.id]] for n in tips}
    return sizes, tipcol


def genome_reduction(out):
    ct = simulate_species_tree(birth=1.0, n_extant=36, seed=4).complete_tree
    # an irreversible lifestyle: once a lineage turns endosymbiont it stays (Dollo-ish)
    hab = simulate_discrete(ct, states=["free-living", "endosymbiont"], start="free-living", seed=8,
                            switch={"free-living->endosymbiont": 0.09, "endosymbiont->free-living": 0.0})
    # the trait CONDITIONS the genome: endosymbionts shed genes fast and gain almost none
    g = simulate_genomes_family(ct, initial_families=55, duplication=0.1,
            origination=3.0 * mod.DrivenBy(hab, {"endosymbiont": 0.1, "free-living": 1.0}),
            loss=0.08 * mod.DrivenBy(hab, {"endosymbiont": 12.0, "free-living": 1.0}), seed=9)
    sizes, tipcol = _sizes(ct, g, hab, _HAB)
    _conditioned_genome(out, ct, [ph.trees.color_history(_state_history(ct, hab), palette=_HAB)],
                        sizes, tipcol, dict(
        driver="lifestyle", states=["free-living", "endosymbiont"],
        switch={"free-living->endosymbiont": 0.09},               # irreversible → one arrow
        mapping={"endosymbiont": 12, "free-living": 1}, target="loss", target_base=0.08,
        state_colors=_HAB))


def genome_expansion(out):
    ct = simulate_species_tree(birth=1.0, n_extant=32, seed=4).complete_tree
    # under relaxed selection (irreversible here) duplicates pile up
    sel = simulate_discrete(ct, states=["purifying", "relaxed"], start="purifying", seed=8,
                            switch={"purifying->relaxed": 0.09, "relaxed->purifying": 0.0})
    g = simulate_genomes_family(ct, initial_families=25, loss=0.07,
            duplication=0.05 * mod.DrivenBy(sel, {"relaxed": 11.0, "purifying": 1.0}), seed=9)
    sizes, tipcol = _sizes(ct, g, sel, _SEL)
    _conditioned_genome(out, ct, [ph.trees.color_history(_state_history(ct, sel), palette=_SEL)],
                        sizes, tipcol, dict(
        driver="selection", states=["purifying", "relaxed"], switch={"purifying->relaxed": 0.09},
        mapping={"relaxed": 11, "purifying": 1}, target="duplication", target_base=0.05,
        state_colors=_SEL))


def hgt_uptake(out):
    ct = simulate_species_tree(birth=1.0, n_extant=30, seed=4).complete_tree
    comp = simulate_discrete(ct, states=["quiet", "competent"], start="quiet", seed=8,
                             switch={"quiet->competent": 0.12, "competent->quiet": 0.05})
    # competence conditions WHO RECEIVES a transfer (the choice slot) — competent lineages take up more
    g = simulate_genomes_family(ct, initial_families=35, transfer=0.5, loss=0.05, duplication=0.03,
            transfer_to=mod.DrivenBy(comp, {"competent": 8.0, "quiet": 1.0}), seed=3)
    sizes, tipcol = _sizes(ct, g, comp, _COMP)
    _conditioned_genome(out, ct, [ph.trees.color_history(_state_history(ct, comp), palette=_COMP)],
                        sizes, tipcol, dict(
        driver="competence", states=["quiet", "competent"],
        switch={"quiet->competent": 0.12, "competent->quiet": 0.05},
        mapping={"competent": 8, "quiet": 1}, target="transfer\nuptake", target_base=None,
        target_sub="who receives a transfer", state_colors=_COMP))


def continuous_conditioning(out):
    """A CONTINUOUS trait conditions a genome rate. A diffusing "activity" trait drives gene gain
    through a Curve (high activity → more originations), so genome size tracks the trait. Same layout
    as the discrete conditioning examples; the diagram's modifier column plots the curve, because a
    continuous driver has no per-state multiplier to list."""
    ct = simulate_species_tree(birth=1.0, n_extant=50, seed=4).complete_tree
    act = simulate_continuous(ct, start=0.0, rate=1.8, seed=3)
    factor = (lambda v: 2.0 ** v)                                # value → factor, the whole mapping
    g = simulate_genomes_family(ct, initial_families=12, loss=0.05,
            origination=0.6 * mod.DrivenBy(act, Curve(factor)), seed=9)
    lab = ct.labels()
    vals = {lab[i]: act.node_values[i] for i in ct.nodes}         # the continuous trait, per node
    tips = list(ct.extant_leaves())
    cmap, norm = cm.viridis, colors.Normalize(min(vals.values()), max(vals.values()))
    sizes = {lab[n.id]: len(g.genomes[n.id]) for n in tips}
    tipcol = {lab[n.id]: colors.to_hex(cmap(norm(act.node_values[n.id]))) for n in tips}
    _conditioned_genome(out, ct, [ph.trees.color_branches(vals, cmap="viridis")],
                        sizes, tipcol, dict(
        draw=h.draw_conditioning_curve, driver="activity", curve=factor,
        vrange=(min(vals.values()), max(vals.values())), value_label="activity",
        target="origination", target_base=0.6))


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
    ax, ["fast", "slow"], palette, {"fast": 2.6, "slow": 0.7}, symbol="λ"))'''

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
lab = ct.labels()                           # {id: 'n<id>'} — or 'e<id>' where the lineage went extinct
dashed = {lab[n.id] for n in ct.extinct_leaves()}                  # (+ their all-extinct ancestors)
history = {lab[i]: segs for i, segs in r.trait.history.items()}
(ph.trees.plot(tree, skeleton=False)
 + ph.trees.color_history(history, palette=palette, dashed=dashed)
 + ph.trees.time_axis("time", tick_size=22, label_size=28)).save("tree.png")
h.composite_markov("tree.png", "sse.png", lambda ax: h.draw_markov(
    ax, ["doomed", "safe"], palette, {"doomed": 0.75, "safe": 0.05}, symbol="μ"))'''

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
lab = ct.labels()                           # {id: 'n<id>'} — or 'e<id>' where the lineage went extinct
dashed = {lab[n.id] for n in ct.extinct_leaves()}
history = {lab[i]: segs for i, segs in r.trait.history.items()}
(ph.trees.plot(tree, skeleton=False)
 + ph.trees.color_history(history, palette=palette, dashed=dashed)
 + ph.trees.time_axis("time", tick_size=22, label_size=28)).save("tree.png")
h.composite_markov("tree.png", "musse.png", lambda ax: h.draw_markov(
    ax, ["slow", "medium", "fast"], palette, {"slow": 0.6, "medium": 1.3, "fast": 2.6},
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
# the SAME DrivenBy that drives speciation with a trait drives a genome rate: endosymbionts
# shed genes fast (loss x12) and gain almost none (origination x0.1)
g = simulate_genomes_family(ct, initial_families=55, duplication=0.1,
        origination=3.0 * mod.DrivenBy(hab, {"endosymbiont": 0.1, "free-living": 1.0}),
        loss=0.08 * mod.DrivenBy(hab, {"endosymbiont": 12.0, "free-living": 1.0}), seed=9)

### plot  —  tree coloured by lifestyle, beside per-tip genome-size bars (aligned axes)
import phylustrator as ph

pal = {"free-living": "#2E8B6F", "endosymbiont": "#C25A3C"}
lab = ct.labels()                                                  # {id: 'n<id>'}
tree = ph.trees.loads(ct.to_newick())
history = {lab[i]: segs for i, segs in hab.history.items()}
tips = list(ct.extant_leaves())
sizes  = {lab[n.id]: len(g.genomes[n.id]) for n in tips}           # gene count per tip
colors = {lab[n.id]: pal[hab.values[n.id]] for n in tips}          # bar colour = lifestyle
fig = (ph.trees.plot(tree, skeleton=False)
       + ph.trees.color_history(history, palette=pal)
       + ph.trees.time_axis("time", bold=False))
ph.beside(fig, ph.genomes.bars(sizes, colors=colors, label="genome size (genes)")).save("reduction.png")
# the figure then composites the driver->modifier->target diagram (lifestyle -> loss) on top'''

_C_EXPANSION = '''\
### simulate  —  a trait conditions DUPLICATION, so genomes grow
from zombi2.species import simulate_species_tree
from zombi2.traits import simulate_discrete
from zombi2.genomes import simulate_genomes_family
from zombi2.rates import modifiers as mod

ct = simulate_species_tree(birth=1.0, n_extant=32, seed=4).complete_tree
sel = simulate_discrete(ct, states=["purifying", "relaxed"], start="purifying", seed=8,
                        switch={"purifying->relaxed": 0.09, "relaxed->purifying": 0.0})  # irreversible
# under relaxed selection, duplicates pile up: DrivenBy on the duplication rate
g = simulate_genomes_family(ct, initial_families=25, loss=0.07,
        duplication=0.05 * mod.DrivenBy(sel, {"relaxed": 11.0, "purifying": 1.0}), seed=9)
### plot  —  tree coloured by selection, beside per-tip genome-size bars (relaxed clades grow)
import phylustrator as ph

pal = {"purifying": "#3A7CA5", "relaxed": "#C25A3C"}
lab = ct.labels()                                                  # {id: 'n<id>'}
tree = ph.trees.loads(ct.to_newick())
history = {lab[i]: segs for i, segs in sel.history.items()}
tips = list(ct.extant_leaves())
sizes  = {lab[n.id]: len(g.genomes[n.id]) for n in tips}
colors = {lab[n.id]: pal[sel.values[n.id]] for n in tips}
fig = (ph.trees.plot(tree, skeleton=False)
       + ph.trees.color_history(history, palette=pal)
       + ph.trees.time_axis("time", bold=False))
ph.beside(fig, ph.genomes.bars(sizes, colors=colors, label="genome size (genes)")).save("expansion.png")
# the figure then composites the driver->modifier->target diagram (selection -> duplication) on top'''

_C_UPTAKE = '''\
### simulate  —  competence conditions WHO RECEIVES a transfer (uptake), not a rate
from zombi2.species import simulate_species_tree
from zombi2.traits import simulate_discrete
from zombi2.genomes import simulate_genomes_family
from zombi2.rates import modifiers as mod

ct = simulate_species_tree(birth=1.0, n_extant=30, seed=4).complete_tree
comp = simulate_discrete(ct, states=["quiet", "competent"], start="quiet", seed=8,
                         switch={"quiet->competent": 0.12, "competent->quiet": 0.05})
# DrivenBy on transfer_to (the choice slot) makes competent lineages likelier recipients, so
# competent genomes take up more DNA
g = simulate_genomes_family(ct, initial_families=35, transfer=0.5, loss=0.05, duplication=0.03,
        transfer_to=mod.DrivenBy(comp, {"competent": 8.0, "quiet": 1.0}), seed=3)
### plot  —  tree coloured by competence, beside per-tip genome-size bars (competent take up more)
import phylustrator as ph

pal = {"quiet": "#8f99a3", "competent": "#2E8B6F"}
lab = ct.labels()                                                  # {id: 'n<id>'}
tree = ph.trees.loads(ct.to_newick())
history = {lab[i]: segs for i, segs in comp.history.items()}
tips = list(ct.extant_leaves())
sizes  = {lab[n.id]: len(g.genomes[n.id]) for n in tips}
colors = {lab[n.id]: pal[comp.values[n.id]] for n in tips}
fig = (ph.trees.plot(tree, skeleton=False)
       + ph.trees.color_history(history, palette=pal)
       + ph.trees.time_axis("time", bold=False))
ph.beside(fig, ph.genomes.bars(sizes, colors=colors, label="genome size (genes)")).save("uptake.png")
# the figure then composites the driver->modifier->target diagram (competence -> uptake) on top'''

_C_CONTINUOUS = '''\
### simulate  —  a CONTINUOUS trait conditions a genome rate (via a Curve, not a state table)
from zombi2.species import simulate_species_tree
from zombi2.traits import simulate_continuous
from zombi2.genomes import simulate_genomes_family
from zombi2.rates import modifiers as mod
from zombi2.rates.mapping import Curve

ct = simulate_species_tree(birth=1.0, n_extant=50, seed=4).complete_tree
act = simulate_continuous(ct, start=0.0, rate=1.8, seed=3)          # a diffusing "activity" trait
# a continuous driver maps its VALUE to a factor with a Curve; here high activity -> more gene gain.
# (Each branch is cut into constant sub-steps internally, so the same engine consumes it.)
factor = lambda v: 2.0 ** v
g = simulate_genomes_family(ct, initial_families=12, loss=0.05,
        origination=0.6 * mod.DrivenBy(act, Curve(factor)), seed=9)

### plot  —  tree coloured by the continuous trait, beside per-tip genome-size bars
import phylustrator as ph
from matplotlib import cm, colors as mcolors

lab = ct.labels()                                                   # {id: 'n<id>'}
tree = ph.trees.loads(ct.to_newick())
vals = {lab[i]: act.node_values[i] for i in ct.nodes}               # the continuous value, per node
tips = list(ct.extant_leaves())
sizes = {lab[n.id]: len(g.genomes[n.id]) for n in tips}
norm = mcolors.Normalize(min(vals.values()), max(vals.values()))
bar_c = {lab[n.id]: mcolors.to_hex(cm.viridis(norm(act.node_values[n.id]))) for n in tips}
fig = (ph.trees.plot(tree, skeleton=False)
       + ph.trees.color_branches(vals, cmap="viridis")
       + ph.trees.time_axis("time", bold=False))
ph.beside(fig, ph.genomes.bars(sizes, colors=bar_c, label="genome size (genes)")).save("cont.png")
# on top goes the same driver->modifier->target diagram, its middle column plotting value -> factor'''


CONDITIONING = [
    Example("genome_reduction", "Genome reduction",
            "A driver (a trait for the lifestyle) modifies the rate of loss (the target). "
            "Endosymbionts also gain genes more slowly, so their genomes shrink. The tree is coloured "
            "by the lifestyle and the bars are genome size at each tip.",
            "trait → loss", genome_reduction, code=_C_REDUCTION),
    Example("genome_expansion", "Genome expansion",
            "A driver (a trait for the strength of selection) modifies the rate of duplication (the "
            "target). Under relaxed selection duplicates accumulate and the genomes grow. The tree is "
            "coloured by the selection regime and the bars are genome size at each tip.",
            "trait → duplication", genome_expansion, code=_C_EXPANSION),
    Example("hgt_uptake", "HGT uptake by competence",
            "A driver (a trait for competence) modifies <b>who receives</b> a transfer (the target), "
            "not a rate. Competent lineages take up DNA more often and their genomes grow. The tree is "
            "coloured by competence and the bars are genome size at each tip.",
            "trait → transfer uptake", hgt_uptake, code=_C_UPTAKE),
    Example("continuous_conditioning", "A continuous driver",
            "A driver (a diffusing continuous trait) modifies the rate of origination (the target). A "
            "<code>Curve</code> turns each value into a factor, so genome size follows the trait. The "
            "tree is coloured by the trait value and the bars are genome size at each tip.",
            "continuous trait → origination", continuous_conditioning, code=_C_CONTINUOUS),
]

JOINING = [
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

EXAMPLES = CONDITIONING + JOINING        # the module's full list; build.py takes the two separately
