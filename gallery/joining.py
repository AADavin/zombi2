"""Two ways one level reads another, through the same ``DrivenBy`` mechanism.

**Conditioning** grows the driver first and holds it fixed: a trait is simulated on the tree, then a
genome run reads it, so the trait's state sets a genome rate. **Joining** grows both at once, because
the trait drives speciation or extinction and so shapes the tree it is evolving on
(``joint.simulate_joint``). The two lists below feed two gallery sections.
"""

from __future__ import annotations

import math

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
_TOX = {"absent": "#b9bec4", "present": "#2E8B6F"}          # a gene family, present or not
_DISEASE = {"harmless": "#8f99a3", "pathogenic": "#C2453C"}
_METAB = {"anaerobic": "#6b5b95", "aerobic": "#2E8B6F"}


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
    tipcol = {lab[n.id]: palette[trait.values[lab[n.id]]] for n in tips}
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


def _continuous_figure(out, factor, *, driver, value_label, target, base, layer_cmap="viridis",
                       rate="origination", trait_rate=1.2, trait_seed=3, genome_seed=9):
    """The continuous-driver figure with the curve left free. `continuous_conditioning` is the
    exponential case; the two below are the same run with a different value → factor function."""
    ct = simulate_species_tree(birth=1.0, n_extant=50, seed=4).complete_tree
    tr = simulate_continuous(ct, start=0.0, rate=trait_rate, seed=trait_seed)
    g = simulate_genomes_family(ct, initial_families=12, loss=0.05,
            **{rate: base * mod.DrivenBy(tr, Curve(factor))}, seed=genome_seed)
    lab = ct.labels()
    vals = {lab[i]: tr.node_values[i] for i in ct.nodes}
    tips = list(ct.extant_leaves())
    cmap, norm = plt.get_cmap(layer_cmap), colors.Normalize(min(vals.values()), max(vals.values()))
    sizes = {lab[n.id]: len(g.genomes[n.id]) for n in tips}
    tipcol = {lab[n.id]: colors.to_hex(cmap(norm(tr.node_values[n.id]))) for n in tips}
    _conditioned_genome(out, ct, [ph.trees.color_branches(vals, cmap=layer_cmap)],
                        sizes, tipcol, dict(
        draw=h.draw_conditioning_curve, driver=driver, curve=factor, cmap=layer_cmap,
        vrange=(min(vals.values()), max(vals.values())), value_label=value_label,
        target=target, target_base=base))


def curve_saturating(out):
    """A SATURATING curve. Gene gain rises with a "resource" trait and then levels off at a ceiling:
    once resources are plentiful, more of them buys nothing. Unlike the exponential, this response is
    bounded — the factor can never exceed 5."""
    factor = (lambda v: 0.2 + 5.8 / (1.0 + math.exp(-1.6 * v)))
    _continuous_figure(out, factor, driver="resources", value_label="resources",
                       target="origination", base=0.6)


def curve_optimum(out):
    """A HUMPED curve. Gene gain is fastest at an intermediate temperature and falls away on both
    sides. A table of per-state multipliers can express neither this nor the saturating curve: the
    response is not monotone, so the tallest bars sit in the middle of the colour ramp, not at
    its top."""
    factor = (lambda v: 0.2 + 5.8 * math.exp(-((v - 1.2) ** 2) / (2 * 0.8 ** 2)))
    _continuous_figure(out, factor, driver="temperature", value_label="temperature",
                       target="origination", base=0.6)


def trait_drives_trait(out):
    """The driver and the target are BOTH traits, on one tree. A continuous "temperature" trait is
    grown first; the rate at which a second trait — body size — diffuses then reads it through the
    same saturating curve. Warm lineages redesign their size quickly and their descendants spread
    apart; cold ones barely move from where their ancestor left them. Nothing here is a new
    mechanism: the driver is grown first and held fixed, which is conditioning, spelled the way
    every other conditioned rate is spelled."""
    ct = simulate_species_tree(birth=1.0, n_extant=40, seed=4).complete_tree
    temp = simulate_continuous(ct, start=0.0, rate=1.2, seed=6)
    # the threshold sits between the two halves of this tree, so one is switched on
    # (~6x) and the other is all but switched off (~0.08x) — a 75-fold contrast
    factor = (lambda v: 0.05 + 6.0 / (1.0 + math.exp(-1.6 * (v - 2.0))))
    size = simulate_continuous(ct, start=0.0, rate=0.35 * mod.DrivenBy(temp, Curve(factor)), seed=11)
    lab = ct.labels()
    driver = {lab[i]: temp.node_values[i] for i in ct.nodes}
    driven = {lab[i]: size.node_values[i] for i in ct.nodes}
    # the driven trait is drawn on a DIVERGING scale centred on where it started, so white reads
    # "has not moved" and the saturated ends read "has moved a long way" — which is what the driver
    # actually controls. On its own min-to-max range the same colours would say nothing about that.
    reach = max(abs(v) for v in driven.values())
    pngs = []
    for k, (vals, cmap, lim) in enumerate(((driver, "viridis", None),
                                           (driven, "coolwarm", (-reach, reach)))):
        png = out.replace(".png", f"_t{k}.png")
        (ph.trees.plot(ph.trees.loads(ct.to_newick()), skeleton=False,
                       style=ph.Style(width=1000, height=560, margin=70, branch_width=3.4))
         + ph.trees.color_branches(vals, cmap=cmap, limits=lim)
         + ph.trees.time_axis("time", tick_size=20, label_size=26, bold=False)).save(png)
        pngs.append(png)
    diag = h.conditioning_png(out.replace(".png", "_diag.png"),
                              draw=h.draw_conditioning_curve, driver="temperature", curve=factor,
                              vrange=(min(driver.values()), max(driver.values())),
                              value_label="temperature", target="size diffusion", target_base=0.35)
    fig = plt.figure(figsize=(12, 11.4))
    fig.add_axes([0.30, 0.845, 0.40, 0.155]).imshow(mpimg.imread(diag))
    rows = ((pngs[0], "temperature — grown first, then held fixed", 0.44),
            (pngs[1], "body size — its diffusion rate reads the temperature", 0.02))
    for png, name, y in rows:
        fig.add_axes([0.0, y, 1.0, 0.375]).imshow(mpimg.imread(png))
        fig.text(0.045, y + 0.385, name, fontsize=15, ha="left", va="bottom")
    for ax in fig.axes:
        ax.set_axis_off()
    fig.savefig(out, dpi=140, bbox_inches="tight")
    plt.close(fig)


def gene_drives_trait(out):
    """A GENE FAMILY is the driver and a trait is the target — the other direction of a relation that
    only ran one way before. A toxin family is grown first, gained and lost down the tree; a
    pathogenicity trait then *becomes* pathogenic forty times faster in the lineages that carry it —
    one direction only, since a toxin makes a lineage dangerous rather than helping it recover. The
    same tree is painted twice, so the answer is in the alignment of the two panels: 80% of the tips
    carrying the gene end up pathogenic against 37% of those without.

    Presence is exact and changes mid-branch, at the instant the last copy actually went."""
    ct = simulate_species_tree(birth=1.0, death=0.2, n_extant=45, seed=4).complete_tree
    # loss well above duplication, so the family is genuinely lost in whole clades: 62% of the
    # tree's branch length carries it and 38% does not, which is what gives the two panels something
    # to disagree about
    g = simulate_genomes_family(ct, initial_families=20, family_names=["tox"],
                                duplication=0.05, loss=0.3, seed=9)
    tox = g.presence("tox")
    # the gene drives ONE direction: carrying a toxin makes a lineage *become* pathogenic, it does
    # not make it revert faster. `switch` takes a rate per transition, and only one of them reads the
    # driver — so the signal is which tips end up pathogenic, not how much they flicker.
    disease = simulate_discrete(
        ct, states=["harmless", "pathogenic"], start="harmless", seed=2,
        switch={"harmless->pathogenic": 0.02 * mod.DrivenBy(tox, {"present": 40.0, "absent": 1.0}),
                "pathogenic->harmless": 0.6})

    lab = ct.labels()
    top = {lab[i]: segs for i, segs in tox.history(ct).items()}
    bottom = _state_history(ct, disease)
    style = ph.Style(width=1000, height=560, margin=70, branch_width=3.4)
    pngs = []
    for k, (hist, palette) in enumerate(((top, _TOX), (bottom, _DISEASE))):
        png = out.replace(".png", f"_g{k}.png")
        (ph.trees.plot(ph.trees.loads(ct.to_newick()), skeleton=False, style=style)
         + ph.trees.color_history(hist, palette=palette)
         + ph.trees.time_axis("time", tick_size=20, label_size=26, bold=False)).save(png)
        pngs.append(png)
    diag = h.conditioning_png(out.replace(".png", "_diag.png"),
                              driver="tox", states=["absent", "present"],
                              switch={"present->absent": 0.3},    # the family is lost, not regained
                              mapping={"present": 40, "absent": 1},
                              target="→ pathogenic",
                              target_base=0.02, state_colors=_TOX,
                              # the target is itself a discrete trait, so it gets its own chain:
                              # what the gene drives is the rate of *these* two arrows
                              target_states=["harmless", "pathogenic"],
                              target_switch={"harmless->pathogenic": 0.8,
                                             "pathogenic->harmless": 0.6},
                              target_colors=_DISEASE,
                              target_driven="harmless->pathogenic")
    fig = plt.figure(figsize=(12, 11.4))
    fig.add_axes([0.30, 0.845, 0.40, 0.155]).imshow(mpimg.imread(diag))
    rows = ((pngs[0], "the toxin family — present or absent", 0.44),
            (pngs[1], "pathogenicity — the gene drives becoming it, not reverting", 0.02))
    for png, name, y in rows:
        fig.add_axes([0.0, y, 1.0, 0.375]).imshow(mpimg.imread(png))
        fig.text(0.045, y + 0.385, name, fontsize=15, ha="left", va="bottom")
    for ax in fig.axes:
        ax.set_axis_off()
    fig.savefig(out, dpi=140, bbox_inches="tight")
    plt.close(fig)


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
colors = {lab[n.id]: pal[hab.values[lab[n.id]]] for n in tips}     # bar colour = lifestyle
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
colors = {lab[n.id]: pal[sel.values[lab[n.id]]] for n in tips}
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
colors = {lab[n.id]: pal[comp.values[lab[n.id]]] for n in tips}
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


_C_SATURATING = '''\
### simulate  —  the same run, a SATURATING curve instead of an exponential one
import math
from zombi2.species import simulate_species_tree
from zombi2.traits import simulate_continuous
from zombi2.genomes import simulate_genomes_family
from zombi2.rates import modifiers as mod
from zombi2.rates.mapping import Curve

ct = simulate_species_tree(birth=1.0, n_extant=50, seed=4).complete_tree
res = simulate_continuous(ct, start=0.0, rate=1.2, seed=3)          # a diffusing "resources" trait
# a logistic response: gene gain switches on as resources rise and then levels off at a ceiling.
# The factor is bounded (0.2 to 6.0), which an exponential curve never is.
factor = lambda v: 0.2 + 5.8 / (1.0 + math.exp(-1.6 * v))
g = simulate_genomes_family(ct, initial_families=12, loss=0.05,
        origination=0.6 * mod.DrivenBy(res, Curve(factor)), seed=9)

### plot  —  identical to the exponential example: tree by trait value, bars by genome size
# (see "A continuous driver" above; only `factor` changed)'''


_C_OPTIMUM = '''\
### simulate  —  the same run again, this time a HUMPED curve
import math
from zombi2.species import simulate_species_tree
from zombi2.traits import simulate_continuous
from zombi2.genomes import simulate_genomes_family
from zombi2.rates import modifiers as mod
from zombi2.rates.mapping import Curve

ct = simulate_species_tree(birth=1.0, n_extant=50, seed=4).complete_tree
temp = simulate_continuous(ct, start=0.0, rate=1.2, seed=3)         # a diffusing "temperature" trait
# gene gain is fastest at an intermediate temperature and falls away on both sides. A Table of
# per-state multipliers cannot express this: the response is not monotone in the driver.
factor = lambda v: 0.2 + 5.8 * math.exp(-((v - 1.2) ** 2) / (2 * 0.8 ** 2))
g = simulate_genomes_family(ct, initial_families=12, loss=0.05,
        origination=0.6 * mod.DrivenBy(temp, Curve(factor)), seed=9)

### plot  —  identical again; the tallest bars now sit in the MIDDLE of the colour ramp
# (see "A continuous driver" above; only `factor` changed)'''


_C_TRAIT_TRAIT = '''\
### simulate  —  the driver and the target are both TRAITS, on one tree
import math
from zombi2.species import simulate_species_tree
from zombi2.traits import simulate_continuous
from zombi2.rates import modifiers as mod
from zombi2.rates.mapping import Curve

ct = simulate_species_tree(birth=1.0, n_extant=40, seed=4).complete_tree
temp = simulate_continuous(ct, start=0.0, rate=1.2, seed=6)         # grown first, then held fixed
# the threshold sits between the two halves of this tree: one is switched on (~6x), the other
# all but switched off (~0.08x). The rate a trait diffuses AT is itself driven.
factor = lambda v: 0.05 + 6.0 / (1.0 + math.exp(-1.6 * (v - 2.0)))
size = simulate_continuous(ct, start=0.0, rate=0.35 * mod.DrivenBy(temp, Curve(factor)), seed=11)

### plot  —  the same tree twice: painted by the driver, then by what it drove
import phylustrator as ph

lab = ct.labels()
tree = ph.trees.loads(ct.to_newick())
driver = {lab[i]: temp.node_values[i] for i in ct.nodes}
driven = {lab[i]: size.node_values[i] for i in ct.nodes}
reach = max(abs(v) for v in driven.values())            # centre the diverging scale on the start
(ph.trees.plot(tree, skeleton=False)
 + ph.trees.color_branches(driver, cmap="viridis")
 + ph.trees.time_axis("time", bold=False)).save("driver.png")
(ph.trees.plot(tree, skeleton=False)
 + ph.trees.color_branches(driven, cmap="coolwarm", limits=(-reach, reach))
 + ph.trees.time_axis("time", bold=False)).save("driven.png")
# white = has not moved from where its ancestor left it; the cold half of the tree stays white'''


def module_drives_metabolism(out):
    """A MODULE drives a trait, through a **discontinuous** response. Four families make up aerobic
    respiration; a lineage that keeps more than half of them turns aerobic, one that drops below
    turns back. The response is a step, not a slope — the Curve is
    ``lambda f: 20.0 if f > 0.5 else 1.0`` — which is how a threshold is written in ZOMBI2: in the
    mapping, where every other response shape lives, rather than as a separate kind of driver.

    Both directions read the module, oppositely, so the trait tracks the gene content rather than
    just accumulating: 97% of the tree's branch length has the trait on the side of the threshold its
    completion is."""
    ct = simulate_species_tree(birth=1.0, death=0.2, n_extant=45, seed=4).complete_tree
    nuo = [f"nuo{c}" for c in "ABCD"]
    # loss just under duplication: the module survives in about half the tree and decays in the rest,
    # which is what puts branches on both sides of the threshold
    g = simulate_genomes_family(ct, initial_families=20, family_names=nuo,
                                modules={"aerobic": nuo}, duplication=0.08, loss=0.06, seed=3)
    comp = g.completion("aerobic")
    step = (lambda f: 20.0 if f > 0.5 else 1.0)
    back = (lambda f: 1.0 if f > 0.5 else 20.0)
    metabolism = simulate_discrete(
        ct, states=["anaerobic", "aerobic"], start="anaerobic", seed=2,
        switch={"anaerobic->aerobic": 0.3 * mod.DrivenBy(comp, Curve(step)),
                "aerobic->anaerobic": 0.3 * mod.DrivenBy(comp, Curve(back))})

    lab = ct.labels()
    levels = sorted({f for segs in comp.history(ct).values() for f, _ in segs})
    ramp = {f: colors.to_hex(cm.viridis(f)) for f in levels}
    top = {lab[i]: segs for i, segs in comp.history(ct).items()}
    bottom = _state_history(ct, metabolism)
    style = ph.Style(width=1000, height=560, margin=70, branch_width=3.4)
    pngs = []
    for k, (hist, palette) in enumerate(((top, ramp), (bottom, _METAB))):
        png = out.replace(".png", f"_m{k}.png")
        (ph.trees.plot(ph.trees.loads(ct.to_newick()), skeleton=False, style=style)
         + ph.trees.color_history(hist, palette=palette)
         + ph.trees.time_axis("time", tick_size=20, label_size=26, bold=False)).save(png)
        pngs.append(png)
    diag = h.conditioning_png(out.replace(".png", "_diag.png"),
                              draw=h.draw_conditioning_curve, driver="aerobic module", curve=step,
                              vrange=(0.0, 1.0), value_label="fraction of the module present",
                              target="→ aerobic", target_base=0.3,
                              target_states=["anaerobic", "aerobic"],
                              target_switch={"anaerobic->aerobic": 6.0,
                                             "aerobic->anaerobic": 0.3},
                              target_colors=_METAB, target_driven="anaerobic->aerobic")
    fig = plt.figure(figsize=(12, 11.4))
    fig.add_axes([0.30, 0.845, 0.40, 0.155]).imshow(mpimg.imread(diag))
    rows = ((pngs[0], "aerobic respiration — how much of the module a lineage keeps", 0.44),
            (pngs[1], "metabolism — aerobic above half the module, anaerobic below", 0.02))
    for png, name, y in rows:
        fig.add_axes([0.0, y, 1.0, 0.375]).imshow(mpimg.imread(png))
        fig.text(0.045, y + 0.385, name, fontsize=15, ha="left", va="bottom")
    for ax in fig.axes:
        ax.set_axis_off()
    fig.savefig(out, dpi=140, bbox_inches="tight")
    plt.close(fig)


_C_GENE_TRAIT = '''\
### simulate  —  a GENE FAMILY drives a trait (the other direction of the same relation)
from zombi2.species import simulate_species_tree
from zombi2.genomes import simulate_genomes_family
from zombi2.traits import simulate_discrete
from zombi2.rates import modifiers as mod

ct = simulate_species_tree(birth=1.0, death=0.2, n_extant=45, seed=4).complete_tree

# 1. the driver: genomes with ONE family named, so it can be referred to. Loss well above
#    duplication, so the family is genuinely lost in whole clades rather than kept everywhere.
g = simulate_genomes_family(ct, initial_families=20, family_names=["tox"],
                            duplication=0.05, loss=0.3, seed=9)

# 2. the target: a trait whose switch rate reads whether that family is there. `presence` is a
#    driver like a grown trait, so the mapping is an ordinary table over its two states.
# `switch` takes a rate per transition, and only ONE of them reads the driver: a toxin makes a
# lineage become pathogenic, it does not help it revert. So the signal is which tips END UP
# pathogenic — 80% of those carrying the gene against 37% without — not how much they flicker.
disease = simulate_discrete(
    ct, states=["harmless", "pathogenic"], start="harmless", seed=2,
    switch={"harmless->pathogenic": 0.02 * mod.DrivenBy(g.presence("tox"),
                                                        {"present": 40.0, "absent": 1.0}),
            "pathogenic->harmless": 0.6})

### plot  —  the same tree painted twice: by the gene, then by what the gene drove
import phylustrator as ph

lab = ct.labels()
tree = ph.trees.loads(ct.to_newick())
gene = {lab[i]: segs for i, segs in g.presence("tox").history(ct).items()}
trait = {lab[i]: segs for i, segs in disease.history.items()}
(ph.trees.plot(tree, skeleton=False)
 + ph.trees.color_history(gene, palette={"absent": "#b9bec4", "present": "#2E8B6F"})
 + ph.trees.time_axis("time", bold=False)).save("gene.png")
(ph.trees.plot(tree, skeleton=False)
 + ph.trees.color_history(trait, palette={"harmless": "#8f99a3", "pathogenic": "#C2453C"})
 + ph.trees.time_axis("time", bold=False)).save("trait.png")
# presence changes MID-BRANCH, at the instant the last copy actually went'''


_C_MODULE = '''\
### simulate  —  a MODULE of gene families drives a trait, through a step
from zombi2.species import simulate_species_tree
from zombi2.genomes import simulate_genomes_family
from zombi2.traits import simulate_discrete
from zombi2.rates import modifiers as mod
from zombi2.rates.mapping import Curve

ct = simulate_species_tree(birth=1.0, death=0.2, n_extant=45, seed=4).complete_tree
nuo = [f"nuo{c}" for c in "ABCD"]

# 1. the driver: four families grouped into one module. Loss just under duplication, so the
#    module survives in about half the tree and decays in the rest.
g = simulate_genomes_family(ct, initial_families=20, family_names=nuo,
                            modules={"aerobic": nuo}, duplication=0.08, loss=0.06, seed=3)

# 2. the target: a DISCONTINUOUS response. `completion` is a number in [0, 1] and the Curve is a
#    step — more than half the module and the lineage turns aerobic, below and it turns back.
#    A threshold is written here, in the mapping, not as a separate kind of driver.
comp = g.completion("aerobic")
metabolism = simulate_discrete(
    ct, states=["anaerobic", "aerobic"], start="anaerobic", seed=2,
    switch={"anaerobic->aerobic": 0.3 * mod.DrivenBy(comp, Curve(lambda f: 20.0 if f > 0.5 else 1.0)),
            "aerobic->anaerobic": 0.3 * mod.DrivenBy(comp, Curve(lambda f: 1.0 if f > 0.5 else 20.0))})

### plot  —  the tree by module completion, then by the metabolism it decides
import phylustrator as ph
from matplotlib import cm, colors

lab = ct.labels()
tree = ph.trees.loads(ct.to_newick())
kept = {lab[i]: segs for i, segs in comp.history(ct).items()}
ramp = {f: colors.to_hex(cm.viridis(f)) for f in {f for s in kept.values() for f, _ in s}}
(ph.trees.plot(tree, skeleton=False)
 + ph.trees.color_history(kept, palette=ramp)
 + ph.trees.time_axis("time", bold=False)).save("module.png")
(ph.trees.plot(tree, skeleton=False)
 + ph.trees.color_history({lab[i]: s for i, s in metabolism.history.items()},
                          palette={"anaerobic": "#6b5b95", "aerobic": "#2E8B6F"})
 + ph.trees.time_axis("time", bold=False)).save("metabolism.png")
# 97% of the tree's branch length has the trait on the side of the threshold its completion is'''


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
    Example("curve_saturating", "A saturating curve",
            "The same run with a different <code>Curve</code>. Gene gain switches on as the trait "
            "rises and then levels off at a ceiling, so the factor is <b>bounded</b> — once resources "
            "are plentiful, more of them buys nothing.",
            "continuous trait → origination", curve_saturating, code=_C_SATURATING),
    Example("curve_optimum", "A humped curve",
            "The same run again, with the response fastest at an <b>intermediate</b> value and falling "
            "away on both sides. A table of per-state multipliers can express neither this nor the "
            "saturating curve: the tallest bars sit in the middle of the colour ramp, not at its top.",
            "continuous trait → origination", curve_optimum, code=_C_OPTIMUM),
    Example("trait_drives_trait", "One trait drives another",
            "Driver and target are both traits on one tree. A temperature trait is grown first; the "
            "rate at which body size diffuses then reads it. The same tree is painted twice — by the "
            "driver, then by what it drove, on a scale centred where the trait started, so white "
            "means <b>has not moved</b>. The cold half of the tree stays white.",
            "trait → trait", trait_drives_trait, code=_C_TRAIT_TRAIT),
    Example("gene_drives_trait", "A gene drives a trait",
            "The other direction of the same relation. A driver (a named gene family, present or "
            "absent) modifies <b>one direction</b> of a trait's switch (the target): carrying a toxin "
            "makes a lineage <i>become</i> pathogenic forty times faster, but does not help it revert. "
            "The same tree is painted twice — by the gene, then by what the gene drove — so the answer "
            "is in the alignment of the two: 80% of the tips carrying the gene end up pathogenic "
            "against 37% of those without. Presence changes <b>mid-branch</b>, at the instant the "
            "last copy went.",
            "gene → trait", gene_drives_trait, code=_C_GENE_TRAIT),
    Example("module_drives_metabolism", "A module, through a step",
            "Genes rarely act alone. Four families make up aerobic respiration, and a driver "
            "(<code>completion</code>, the fraction of the module a lineage keeps) modifies the trait "
            "that decides its metabolism. The response is <b>discontinuous</b> — "
            "<code>lambda f: 20.0 if f > 0.5 else 1.0</code> — so more than half the module makes a "
            "lineage aerobic and less makes it revert. That is where a threshold belongs in ZOMBI2: "
            "in the mapping, alongside every other response shape, rather than as its own kind of "
            "driver. 97% of the tree's branch length has the trait on the side of the threshold its "
            "completion is.",
            "module → trait", module_drives_metabolism, code=_C_MODULE),
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
