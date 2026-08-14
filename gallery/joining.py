"""Two ways one evolved value drives another, through the same driver mechanism.

**Conditioning** grows the driver first and holds it fixed: a trait is simulated on the tree, then a
second run reads it, driving a genome rate, another trait's rate, or which lineage receives a
transfer. **Joining** grows both at once, because
the trait drives speciation or extinction and so shapes the tree it is evolving on
(``joint.simulate_joint``). The two lists below feed two gallery sections.
"""

from __future__ import annotations

import math

import matplotlib.pyplot as plt
from matplotlib import cm, colors

import helpers as h
from helpers import Example

import phylustrator as ph
from zombi2.params import Curve, PerCopy, PerLineage, Recipients
from zombi2 import joint, traits
from zombi2.genomes import family as genomes_spec
from zombi2.genomes import simulate_genomes_family
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
_KEY = {"absent": "#b9bec4", "present": "#2E8B6F"}          # the gene that drives the splitting


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
        birth=PerLineage(1.0).scaled_by("trait", {"fast": 2.6, "slow": 0.7}),
        trait=traits.discrete(states=["fast", "slow"], switch=0.35),
        n_extant=70, seed=3)
    tree_png = out.replace(".png", "_tree.png")
    (ph.trees.plot(ph.trees.loads(r.complete_tree.to_newick()), style=_style(), skeleton=False)
     + ph.trees.color_history(_history(r), palette=_BISSE)
     + ph.trees.time_axis("time", tick_size=22, label_size=28)).save(tree_png)
    h.composite_markov(tree_png, out, lambda ax: h.draw_markov(
        ax, ["fast", "slow"], _BISSE, {"fast": 2.6, "slow": 0.7}, symbol="λ"),
        loc=(0.02, 0.04, 0.34, 0.30))


def key_innovation(out):
    """GENE CONTENT drives speciation — the joint model that is not about a trait. A lineage carrying
    the family splits five times as often, and the tree is an **output** of that.

    The family is losing ground on its own: loss runs at twice origination, and on a tree that does
    not feel it that leaves it on about 9% of the tips (eight seeds, 0-36%). Here it reaches 83%,
    because the lineages that carry it split faster than the ones that shed it. Nothing selects for
    the gene directly — what selects for it is the tree it is shaping, which is what makes this
    joint rather than a genome run on a tree handed to it."""
    r = joint.simulate_joint(
        birth=PerLineage(0.7).scaled_by("genomes:toxin", {"present": 5.0, "absent": 1.0}),
        death=0.15,
        genome=genomes_spec(origination=0.15, loss=0.3, family_names=["toxin"],
                            initial_families=5),
        n_extant=70, seed=3)
    ct = r.complete_tree
    lab = ct.labels()
    history = {lab[i]: segs for i, segs in r.genome.presence("toxin").history(ct).items()}
    tree_png = out.replace(".png", "_tree.png")
    (ph.trees.plot(ph.trees.loads(ct.to_newick()), style=_style(), skeleton=False)
     + ph.trees.color_history(history, palette=_KEY)
     + ph.trees.time_axis("time", tick_size=22, label_size=28)).save(tree_png)
    h.composite_markov(tree_png, out, lambda ax: h.draw_markov(
        ax, ["absent", "present"], _KEY, {"absent": 0.7, "present": 3.5}, symbol="λ"),
        loc=(0.02, 0.04, 0.34, 0.30))


def state_extinction(out):
    r = joint.simulate_joint(
        birth=1.0,
        death=PerLineage(1.0).scaled_by("trait", {"doomed": 0.75, "safe": 0.05}),
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
        birth=PerLineage(1.0).scaled_by("trait", {"slow": 0.6, "medium": 1.3, "fast": 2.6}),
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
    """The conditioning-figure layout with genome size on the bars — `helpers.conditioned_figure`,
    which the cross-level examples share at other targets."""
    h.conditioned_figure(out, ct, layers, sizes, tipcol, diagram)


def _panel_style():
    """The tree style the three paint-the-tree-twice figures share.

    Wide and short on purpose: each panel is laid out at its own aspect ratio, so a squarer tree
    makes the whole figure taller than a gallery card can show. The margin has to clear the time
    axis *and* its label — at 64 the word "time" was cut in half inside the PNG, before the
    compositing ever saw it.
    """
    return ph.Style(width=1400, height=660, margin=96, branch_width=3.0)


def _sizes(ct, g, trait, palette):
    lab = ct.labels()
    tips = list(ct.extant_leaves())
    sizes = {name: len(gen) for name, gen in g.genomes.items()}
    tipcol = {lab[n]: palette[trait.values[lab[n]]] for n in tips}
    return sizes, tipcol


# The lifestyle / selection chains below run in BOTH directions. An irreversible switch saturates:
# the stem alone is ~4 time units, so one draw there paints the whole tree one colour and the figure
# has nothing left to show. A reversible chain gives the mosaic these figures exist for, and
# ``tests/test_gallery_api.py`` pins it — both states must reach the extant tips.
_LIFESTYLE = {"free-living->endosymbiont": 0.20, "endosymbiont->free-living": 0.08}
_SELECTION = {"purifying->relaxed": 0.20, "relaxed->purifying": 0.08}


def genome_reduction(out):
    ct = simulate_species_tree(birth=1.0, n_extant=36, seed=4).complete_tree
    hab = simulate_discrete(ct, states=["free-living", "endosymbiont"], start="free-living", seed=6,
                            switch=_LIFESTYLE)
    # the trait CONDITIONS the genome: endosymbionts shed genes fast and gain few
    g = simulate_genomes_family(ct, initial_families=200, duplication=0.1,
            origination=PerLineage(3.0).scaled_by(hab, {"endosymbiont": 0.3, "free-living": 1.0}),
            loss=PerCopy(0.08).scaled_by(hab, {"endosymbiont": 6.0, "free-living": 1.0}), seed=9)
    sizes, tipcol = _sizes(ct, g, hab, _HAB)
    _conditioned_genome(out, ct, [ph.trees.color_history(_state_history(ct, hab), palette=_HAB)],
                        sizes, tipcol, dict(
        driver=("traits", "lifestyle", "two states"),
        connection=("scaled_by", "table"),
        target_level="genomes",
        targets=[("loss", "rate · per copy", "endosymbiont × 6"),
                 ("origination", "rate · per lineage", "endosymbiont × 0.3")],
        chain=(("free-living", "endosymbiont"),
               [("0.20", "0.08")],
               (_HAB["free-living"], _HAB["endosymbiont"]))))


def genome_expansion(out):
    ct = simulate_species_tree(birth=1.0, n_extant=32, seed=4).complete_tree
    sel = simulate_discrete(ct, states=["purifying", "relaxed"], start="purifying", seed=6,
                            switch=_SELECTION)
    # under relaxed selection duplicates pile up
    g = simulate_genomes_family(ct, initial_families=120, loss=0.07,
            duplication=PerCopy(0.05).scaled_by(sel, {"relaxed": 11.0, "purifying": 1.0}), seed=9)
    sizes, tipcol = _sizes(ct, g, sel, _SEL)
    _conditioned_genome(out, ct, [ph.trees.color_history(_state_history(ct, sel), palette=_SEL)],
                        sizes, tipcol, dict(
        driver=("traits", "selection", "two states"),
        connection=("scaled_by", "table"),
        target_level="genomes",
        targets=[("duplication", "rate · per copy", "relaxed × 11    purifying × 1")],
        chain=(("purifying", "relaxed"),
               [("0.20", "0.08")],
               (_SEL["purifying"], _SEL["relaxed"]))))


def hgt_uptake(out):
    ct = simulate_species_tree(birth=1.0, n_extant=30, seed=4).complete_tree
    comp = simulate_discrete(ct, states=["quiet", "competent"], start="quiet", seed=8,
                             switch={"quiet->competent": 0.12, "competent->quiet": 0.05})
    # competence conditions WHO RECEIVES a transfer (the choice). Competent lineages take up more
    g = simulate_genomes_family(ct, initial_families=35, transfer=0.5, loss=0.05, duplication=0.03,
            transfer_to=Recipients().weighted_by(comp, {"competent": 8.0, "quiet": 1.0}), seed=3)
    sizes, tipcol = _sizes(ct, g, comp, _COMP)
    _conditioned_genome(out, ct, [ph.trees.color_history(_state_history(ct, comp), palette=_COMP)],
                        sizes, tipcol, dict(
        driver=("traits", "competence", "two states"),
        connection=("weighted_by", "table"),
        target_level="genomes",
        targets=[("transfer_to", "choice · a weight", "competent × 8    quiet × 1")],
        chain=(("quiet", "competent"),
               [("0.12", "0.05")],
               (_COMP["quiet"], _COMP["competent"]))))


def continuous_conditioning(out):
    """A CONTINUOUS trait conditions a genome rate. A diffusing "activity" trait drives gene gain
    through a Curve (high activity → more originations), so genome size tracks the trait. Same layout
    as the discrete conditioning examples; the diagram's modifier column plots the curve, because a
    continuous driver has no per-state multiplier to list."""
    ct = simulate_species_tree(birth=1.0, n_extant=50, seed=4).complete_tree
    act = simulate_continuous(ct, start=0.0, rate=1.8, seed=3)
    factor = (lambda v: 2.0 ** v)                                # value → factor, the whole mapping
    g = simulate_genomes_family(ct, initial_families=12, loss=0.05,
            origination=PerLineage(0.6).scaled_by(act, Curve(factor)), seed=9)
    lab = ct.labels()
    vals = {lab[i]: act.node_values[i] for i in ct.nodes}         # the continuous trait, per node
    tips = list(ct.extant_leaves())
    cmap, norm = cm.viridis, colors.Normalize(min(vals.values()), max(vals.values()))
    sizes = {name: len(gen) for name, gen in g.genomes.items()}
    tipcol = {lab[n]: colors.to_hex(cmap(norm(act.node_values[n]))) for n in tips}
    _conditioned_genome(out, ct, [ph.trees.color_branches(vals, cmap="viridis")],
                        sizes, tipcol, dict(
        driver=("traits", "activity", "a number"),
        connection=("scaled_by", "curve"),
        target_level="genomes",
        targets=[("origination", "rate · per lineage", "")],
        curve=(factor, "activity", (min(vals.values()), max(vals.values())))))


def _continuous_figure(out, factor, *, driver, value_label, target, base, layer_cmap="viridis",
                       rate="origination", trait_rate=1.2, trait_seed=3, genome_seed=9):
    """The continuous-driver figure with the curve left free. `continuous_conditioning` is the
    exponential case; the two below are the same run with a different value → factor function."""
    ct = simulate_species_tree(birth=1.0, n_extant=50, seed=4).complete_tree
    tr = simulate_continuous(ct, start=0.0, rate=trait_rate, seed=trait_seed)
    g = simulate_genomes_family(ct, initial_families=12, loss=0.05,
            **{rate: PerLineage(base).scaled_by(tr, Curve(factor))}, seed=genome_seed)
    lab = ct.labels()
    vals = {lab[i]: tr.node_values[i] for i in ct.nodes}
    tips = list(ct.extant_leaves())
    cmap, norm = plt.get_cmap(layer_cmap), colors.Normalize(min(vals.values()), max(vals.values()))
    sizes = {name: len(gen) for name, gen in g.genomes.items()}
    tipcol = {lab[n]: colors.to_hex(cmap(norm(tr.node_values[n]))) for n in tips}
    _conditioned_genome(out, ct, [ph.trees.color_branches(vals, cmap=layer_cmap)],
                        sizes, tipcol, dict(
        driver=("traits", driver, "a number"),
        connection=("scaled_by", "curve"),
        target_level="genomes",
        targets=[(target, "rate · per lineage", "")],
        curve=(factor, value_label, (min(vals.values()), max(vals.values())))))


def curve_saturating(out):
    """A SATURATING curve. Gene gain rises with a "resource" trait and then levels off at a ceiling:
    once resources are plentiful, more of them buys nothing. Unlike the exponential, this response is
    bounded — the factor can never exceed 6."""
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
    size = simulate_continuous(ct, start=0.0, seed=11,
                               rate=PerLineage(0.35).scaled_by(temp, Curve(factor)))
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
                       style=_panel_style())
         + ph.trees.color_branches(vals, cmap=cmap, limits=lim)
         + ph.trees.time_axis("time", tick_size=20, label_size=26, bold=False)).save(png)
        pngs.append(png)
    diag = h.conditioning_png(out.replace(".png", "_diag.png"),
                              driver=("traits", "temperature", "a number"),
                              connection=("scaled_by", "curve"),
                              target_level="traits",
                              targets=[("size diffusion", "rate · per lineage", "")],
                              curve=(factor, "temperature",
                                     (min(driver.values()), max(driver.values()))))
    h.composite_under_diagram(out, diag,
                              [(pngs[0], "temperature", ("viridis", "cold", "warm")),
                               (pngs[1], "body size",
                                ("coolwarm", "smaller", "has not moved", "larger"))])


def gene_drives_trait(out):
    """A GENE FAMILY is the driver and a trait is the target — the other direction of a relation that
    only ran one way before. A toxin family is grown first, gained and lost down the tree; a
    pathogenicity trait then *becomes* pathogenic forty times faster in the lineages that carry it —
    one direction only, since a toxin makes a lineage dangerous rather than helping it recover. The
    same tree is painted twice, so the answer is in the alignment of the two panels: 74% of the tips
    carrying the gene end up pathogenic against none of those without.

    Presence is exact and changes mid-branch, at the instant the last copy actually went."""
    ct = simulate_species_tree(birth=1.0, death=0.2, n_extant=45, seed=4).complete_tree
    # Loss a little above duplication, so the family survives in some clades and decays in others:
    # 59% of the tree's branch length carries it, and 23 tips do against 22 that do not. At the old
    # loss=0.3 the family died on the stem and both panels came out blank.
    g = simulate_genomes_family(ct, initial_families=20, family_names=["tox"],
                                duplication=0.1, loss=0.13, seed=5)
    tox = g.presence("tox")
    # the gene drives ONE direction: carrying a toxin makes a lineage *become* pathogenic, it does
    # not make it revert faster. `switch` takes a rate per transition, and only one of them reads the
    # driver — so the signal is which tips end up pathogenic, not how much they flicker.
    disease = simulate_discrete(
        ct, states=["harmless", "pathogenic"], start="harmless", seed=2,
        switch={"harmless->pathogenic":
                    PerLineage(0.02).scaled_by(tox, {"present": 40.0, "absent": 1.0}),
                "pathogenic->harmless": 0.6})

    lab = ct.labels()
    top = {lab[i]: segs for i, segs in tox.history(ct).items()}
    bottom = _state_history(ct, disease)
    style = _panel_style()
    pngs = []
    for k, (hist, palette) in enumerate(((top, _TOX), (bottom, _DISEASE))):
        png = out.replace(".png", f"_g{k}.png")
        (ph.trees.plot(ph.trees.loads(ct.to_newick()), skeleton=False, style=style)
         + ph.trees.color_history(hist, palette=palette)
         + ph.trees.time_axis("time", tick_size=20, label_size=26, bold=False)).save(png)
        pngs.append(png)
    diag = h.conditioning_png(
        out.replace(".png", "_diag.png"),
        driver=("genomes", "tox", "present or absent"),
        connection=("scaled_by", "table"),
        target_level="traits",
        targets=[("pathogenic", "rate · per lineage", "present × 40    absent × 1")],
        # the family is lost and not regained, so one arc; the TARGET is itself a trait with states,
        # and what the gene drives is the rate of one of its two arrows
        chain=(("present", "absent"), [("0.13", None)],
               (_TOX["present"], _TOX["absent"])),
        target_chain=(("harmless", "pathogenic"), [("driven", "0.6")],
                      (_DISEASE["harmless"], _DISEASE["pathogenic"])))
    h.composite_under_diagram(out, diag, [(pngs[0], "the toxin family", _TOX),
                                          (pngs[1], "pathogenicity", _DISEASE)])


_C_BISSE = '''\
### simulate  —  a 2-state trait drives speciation (BiSSE)
from zombi2 import joint, traits
from zombi2.params import PerLineage

r = joint.simulate_joint(
    birth=PerLineage(1.0).scaled_by("trait", {"fast": 2.6, "slow": 0.7}),
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
from zombi2.params import PerLineage

r = joint.simulate_joint(
    birth=1.0,
    death=PerLineage(1.0).scaled_by("trait", {"doomed": 0.75, "safe": 0.05}),
    trait=traits.discrete(states=["doomed", "safe"], switch=0.3),
    n_extant=35, seed=1)
ct = r.complete_tree

### plot  —  coloured by state history, extinct lineages dashed, model inset
import phylustrator as ph
import helpers as h

palette = {"doomed": "#B0413E", "safe": "#2A9D8F"}
tree = ph.trees.loads(ct.to_newick())
lab = ct.labels()                           # {id: 'n<id>'} — or 'e<id>' where the lineage went extinct
dashed = {lab[n.id] for n in (ct.nodes[_i] for _i in ct.extinct_leaves())}                  # (+ their all-extinct ancestors)
history = {lab[i]: segs for i, segs in r.trait.history.items()}
(ph.trees.plot(tree, skeleton=False)
 + ph.trees.color_history(history, palette=palette, dashed=dashed)
 + ph.trees.time_axis("time", tick_size=22, label_size=28)).save("tree.png")
h.composite_markov("tree.png", "sse.png", lambda ax: h.draw_markov(
    ax, ["doomed", "safe"], palette, {"doomed": 0.75, "safe": 0.05}, symbol="μ"))'''

_C_MUSSE = '''\
### simulate  —  three graded speciation rates + constant death (MuSSE)
from zombi2 import joint, traits
from zombi2.params import PerLineage

r = joint.simulate_joint(
    birth=PerLineage(1.0).scaled_by("trait", {"slow": 0.6, "medium": 1.3, "fast": 2.6}),
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
dashed = {lab[n.id] for n in (ct.nodes[_i] for _i in ct.extinct_leaves())}
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
from zombi2.params import PerCopy, PerLineage

sp = simulate_species_tree(birth=1.0, n_extant=36, seed=4)
ct = sp.complete_tree
# the lifestyle switches both ways, so the tree ends up a mosaic rather than one colour
hab = simulate_discrete(ct, states=["free-living", "endosymbiont"], start="free-living", seed=6,
                        switch={"free-living->endosymbiont": 0.20, "endosymbiont->free-living": 0.08})
# the SAME driver that drives speciation with a trait drives a genome rate: endosymbionts
# shed genes fast (loss x6) and gain few (origination x0.3)
g = simulate_genomes_family(ct, initial_families=200, duplication=0.1,
        origination=PerLineage(3.0).scaled_by(hab, {"endosymbiont": 0.3, "free-living": 1.0}),
        loss=PerCopy(0.08).scaled_by(hab, {"endosymbiont": 6.0, "free-living": 1.0}), seed=9)
# 241 genes at the median free-living tip against 48 at the median endosymbiont one

### plot  —  tree coloured by lifestyle, beside per-tip genome-size bars (aligned axes)
import phylustrator as ph

pal = {"free-living": "#2E8B6F", "endosymbiont": "#C25A3C"}
lab = ct.labels()                                                  # {id: 'n<id>'}
tree = ph.trees.loads(ct.to_newick())
history = {lab[i]: segs for i, segs in hab.history.items()}
tips = list(ct.extant_leaves())
sizes  = {name: len(gen) for name, gen in g.genomes.items()}  # gene count per tip
colors = {lab[n]: pal[hab.values[lab[n]]] for n in tips}     # bar colour = lifestyle
fig = (ph.trees.plot(tree, skeleton=False)
       + ph.trees.color_history(history, palette=pal)
       + ph.trees.time_axis("time", bold=False))
ph.beside(fig, ph.genomes.bars(sizes, colors=colors, label="genome size (genes)")).save("reduction.png")
# the figure then composites the driver->mapping->target diagram (lifestyle -> loss) on top'''

_C_EXPANSION = '''\
### simulate  —  a trait conditions DUPLICATION, so genomes grow
from zombi2.species import simulate_species_tree
from zombi2.traits import simulate_discrete
from zombi2.genomes import simulate_genomes_family
from zombi2.params import PerCopy

ct = simulate_species_tree(birth=1.0, n_extant=32, seed=4).complete_tree
sel = simulate_discrete(ct, states=["purifying", "relaxed"], start="purifying", seed=6,
                        switch={"purifying->relaxed": 0.20, "relaxed->purifying": 0.08})
# under relaxed selection, duplicates pile up: scaled_by on the duplication rate
g = simulate_genomes_family(ct, initial_families=120, loss=0.07,
        duplication=PerCopy(0.05).scaled_by(sel, {"relaxed": 11.0, "purifying": 1.0}), seed=9)
# 96 genes at the median purifying tip against 496 at the median relaxed one
### plot  —  tree coloured by selection, beside per-tip genome-size bars (relaxed clades grow)
import phylustrator as ph

pal = {"purifying": "#3A7CA5", "relaxed": "#C25A3C"}
lab = ct.labels()                                                  # {id: 'n<id>'}
tree = ph.trees.loads(ct.to_newick())
history = {lab[i]: segs for i, segs in sel.history.items()}
tips = list(ct.extant_leaves())
sizes  = {name: len(gen) for name, gen in g.genomes.items()}
colors = {lab[n]: pal[sel.values[lab[n]]] for n in tips}
fig = (ph.trees.plot(tree, skeleton=False)
       + ph.trees.color_history(history, palette=pal)
       + ph.trees.time_axis("time", bold=False))
ph.beside(fig, ph.genomes.bars(sizes, colors=colors, label="genome size (genes)")).save("expansion.png")
# the figure then composites the driver->mapping->target diagram (selection -> duplication) on top'''

_C_UPTAKE = '''\
### simulate  —  competence conditions WHO RECEIVES a transfer (uptake), not a rate
from zombi2.species import simulate_species_tree
from zombi2.traits import simulate_discrete
from zombi2.genomes import simulate_genomes_family
from zombi2.params import Recipients

ct = simulate_species_tree(birth=1.0, n_extant=30, seed=4).complete_tree
comp = simulate_discrete(ct, states=["quiet", "competent"], start="quiet", seed=8,
                         switch={"quiet->competent": 0.12, "competent->quiet": 0.05})
# weighted_by on transfer_to (the choice) makes competent lineages likelier recipients, so
# competent genomes take up more DNA
g = simulate_genomes_family(ct, initial_families=35, transfer=0.5, loss=0.05, duplication=0.03,
        transfer_to=Recipients().weighted_by(comp, {"competent": 8.0, "quiet": 1.0}), seed=3)
### plot  —  tree coloured by competence, beside per-tip genome-size bars (competent take up more)
import phylustrator as ph

pal = {"quiet": "#8f99a3", "competent": "#2E8B6F"}
lab = ct.labels()                                                  # {id: 'n<id>'}
tree = ph.trees.loads(ct.to_newick())
history = {lab[i]: segs for i, segs in comp.history.items()}
tips = list(ct.extant_leaves())
sizes  = {name: len(gen) for name, gen in g.genomes.items()}
colors = {lab[n]: pal[comp.values[lab[n]]] for n in tips}
fig = (ph.trees.plot(tree, skeleton=False)
       + ph.trees.color_history(history, palette=pal)
       + ph.trees.time_axis("time", bold=False))
ph.beside(fig, ph.genomes.bars(sizes, colors=colors, label="genome size (genes)")).save("uptake.png")
# the figure then composites the driver->mapping->target diagram (competence -> uptake) on top'''

_C_CONTINUOUS = '''\
### simulate  —  a CONTINUOUS trait conditions a genome rate (via a Curve, not a state table)
from zombi2.species import simulate_species_tree
from zombi2.traits import simulate_continuous
from zombi2.genomes import simulate_genomes_family
from zombi2.params import Curve, PerLineage

ct = simulate_species_tree(birth=1.0, n_extant=50, seed=4).complete_tree
act = simulate_continuous(ct, start=0.0, rate=1.8, seed=3)          # a diffusing "activity" trait
# a continuous driver maps its VALUE to a factor with a Curve; here high activity -> more gene gain.
# (Each branch is cut into constant sub-steps internally, so the same engine consumes it.)
factor = lambda v: 2.0 ** v
g = simulate_genomes_family(ct, initial_families=12, loss=0.05,
        origination=PerLineage(0.6).scaled_by(act, Curve(factor)), seed=9)

### plot  —  tree coloured by the continuous trait, beside per-tip genome-size bars
import phylustrator as ph
from matplotlib import cm, colors as mcolors

lab = ct.labels()                                                   # {id: 'n<id>'}
tree = ph.trees.loads(ct.to_newick())
vals = {lab[i]: act.node_values[i] for i in ct.nodes}               # the continuous value, per node
tips = list(ct.extant_leaves())
sizes = {name: len(gen) for name, gen in g.genomes.items()}
norm = mcolors.Normalize(min(vals.values()), max(vals.values()))
bar_c = {lab[n]: mcolors.to_hex(cm.viridis(norm(act.node_values[n]))) for n in tips}
fig = (ph.trees.plot(tree, skeleton=False)
       + ph.trees.color_branches(vals, cmap="viridis")
       + ph.trees.time_axis("time", bold=False))
ph.beside(fig, ph.genomes.bars(sizes, colors=bar_c, label="genome size (genes)")).save("cont.png")
# on top goes the same driver->mapping->target diagram, its middle column plotting value -> factor'''


_C_SATURATING = '''\
### simulate  —  the same run, a SATURATING curve instead of an exponential one
import math
from zombi2.species import simulate_species_tree
from zombi2.traits import simulate_continuous
from zombi2.genomes import simulate_genomes_family
from zombi2.params import Curve, PerLineage

ct = simulate_species_tree(birth=1.0, n_extant=50, seed=4).complete_tree
res = simulate_continuous(ct, start=0.0, rate=1.2, seed=3)          # a diffusing "resources" trait
# a logistic response: gene gain switches on as resources rise and then levels off at a ceiling.
# The factor is bounded (0.2 to 6.0), which an exponential curve never is.
factor = lambda v: 0.2 + 5.8 / (1.0 + math.exp(-1.6 * v))
g = simulate_genomes_family(ct, initial_families=12, loss=0.05,
        origination=PerLineage(0.6).scaled_by(res, Curve(factor)), seed=9)

### plot  —  identical to the exponential example: tree by trait value, bars by genome size
# (see "A continuous driver" above; only `factor` changed)'''


_C_OPTIMUM = '''\
### simulate  —  the same run again, with a curve of your own
import math
from zombi2.species import simulate_species_tree
from zombi2.traits import simulate_continuous
from zombi2.genomes import simulate_genomes_family
from zombi2.params import Curve, PerLineage

ct = simulate_species_tree(birth=1.0, n_extant=50, seed=4).complete_tree
temp = simulate_continuous(ct, start=0.0, rate=1.2, seed=3)         # a diffusing "temperature" trait
# gene gain is fastest at an intermediate temperature and falls away on both sides. A Table of
# per-state multipliers cannot express this: the response is not monotone in the driver.
factor = lambda v: 0.2 + 5.8 * math.exp(-((v - 1.2) ** 2) / (2 * 0.8 ** 2))
g = simulate_genomes_family(ct, initial_families=12, loss=0.05,
        origination=PerLineage(0.6).scaled_by(temp, Curve(factor)), seed=9)

### plot  —  identical again; the tallest bars now sit in the MIDDLE of the colour ramp
# (see "A continuous driver" above; only `factor` changed)'''


_C_TRAIT_TRAIT = '''\
### simulate  —  the driver and the target are both TRAITS, on one tree
import math
from zombi2.species import simulate_species_tree
from zombi2.traits import simulate_continuous
from zombi2.params import Curve, PerLineage

ct = simulate_species_tree(birth=1.0, n_extant=40, seed=4).complete_tree
temp = simulate_continuous(ct, start=0.0, rate=1.2, seed=6)         # grown first, then held fixed
# the threshold sits between the two halves of this tree: one is switched on (~6x), the other
# all but switched off (~0.08x). The rate a trait diffuses AT is itself driven.
factor = lambda v: 0.05 + 6.0 / (1.0 + math.exp(-1.6 * (v - 2.0)))
size = simulate_continuous(ct, start=0.0, rate=PerLineage(0.35).scaled_by(temp, Curve(factor)), seed=11)

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
    just accumulating: every tip that keeps more than half the module ends up aerobic, and every tip
    that does not ends up anaerobic."""
    ct = simulate_species_tree(birth=1.0, death=0.2, n_extant=45, seed=4).complete_tree
    nuo = [f"nuo{c}" for c in "ABCD"]
    # Loss matched to duplication: the module stays above half on 51% of the tree's branch length and
    # decays below it on the rest, which is what puts branches on both sides of the threshold. At the
    # old duplication=0.08 / loss=0.06 it survived nearly everywhere and the figure was one colour.
    g = simulate_genomes_family(ct, initial_families=20, family_names=nuo,
                                modules={"aerobic": nuo}, duplication=0.08, loss=0.08, seed=11)
    comp = g.completion("aerobic")
    step = (lambda f: 20.0 if f > 0.5 else 1.0)
    back = (lambda f: 1.0 if f > 0.5 else 20.0)
    metabolism = simulate_discrete(
        ct, states=["anaerobic", "aerobic"], start="anaerobic", seed=2,
        switch={"anaerobic->aerobic": PerLineage(0.3).scaled_by(comp, Curve(step)),
                "aerobic->anaerobic": PerLineage(0.3).scaled_by(comp, Curve(back))})

    lab = ct.labels()
    levels = sorted({f for segs in comp.history(ct).values() for f, _ in segs})
    ramp = {f: colors.to_hex(cm.viridis(f)) for f in levels}
    top = {lab[i]: segs for i, segs in comp.history(ct).items()}
    bottom = _state_history(ct, metabolism)
    style = _panel_style()
    pngs = []
    for k, (hist, palette) in enumerate(((top, ramp), (bottom, _METAB))):
        png = out.replace(".png", f"_m{k}.png")
        (ph.trees.plot(ph.trees.loads(ct.to_newick()), skeleton=False, style=style)
         + ph.trees.color_history(hist, palette=palette)
         + ph.trees.time_axis("time", tick_size=20, label_size=26, bold=False)).save(png)
        pngs.append(png)
    diag = h.conditioning_png(
        out.replace(".png", "_diag.png"),
        driver=("genomes", "aerobic module", "a fraction, 0–1"),
        connection=("scaled_by", "curve"),
        target_level="traits",
        targets=[("aerobic", "rate · per lineage", "")],
        curve=(step, "fraction of the module present", (0.0, 1.0)),
        target_chain=(("anaerobic", "aerobic"), [("driven", "0.3")],
                      (_METAB["anaerobic"], _METAB["aerobic"])))
    h.composite_under_diagram(out, diag,
                              [(pngs[0], "the aerobic module",
                                ("viridis", "none of it", "all of it")),
                               (pngs[1], "metabolism", _METAB)])


_C_GENE_TRAIT = '''\
### simulate  —  a GENE FAMILY drives a trait (the other direction of the same relation)
from zombi2.species import simulate_species_tree
from zombi2.genomes import simulate_genomes_family
from zombi2.traits import simulate_discrete
from zombi2.params import PerLineage

ct = simulate_species_tree(birth=1.0, death=0.2, n_extant=45, seed=4).complete_tree

# 1. the driver: genomes with ONE family named, so it can be referred to. Loss a little above
#    duplication, so the family survives in some clades and decays in others.
g = simulate_genomes_family(ct, initial_families=20, family_names=["tox"],
                            duplication=0.1, loss=0.13, seed=5)

# 2. the target: a trait whose switch rate reads whether that family is there. `presence` is a
#    driver like a grown trait, so the mapping is an ordinary table over its two states.
# `switch` takes a rate per transition, and only ONE of them reads the driver: a toxin makes a
# lineage become pathogenic, it does not help it revert. So the signal is which tips END UP
# pathogenic — 74% of the 23 tips carrying the gene, against none of the 22 without.
disease = simulate_discrete(
    ct, states=["harmless", "pathogenic"], start="harmless", seed=2,
    switch={"harmless->pathogenic": PerLineage(0.02).scaled_by(g.presence("tox"),
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
from zombi2.params import Curve, PerLineage

ct = simulate_species_tree(birth=1.0, death=0.2, n_extant=45, seed=4).complete_tree
nuo = [f"nuo{c}" for c in "ABCD"]

# 1. the driver: four families grouped into one module. Loss matched to duplication, so the
#    module survives in about half the tree and decays in the rest.
g = simulate_genomes_family(ct, initial_families=20, family_names=nuo,
                            modules={"aerobic": nuo}, duplication=0.08, loss=0.08, seed=11)

# 2. the target: a DISCONTINUOUS response. `completion` is a number in [0, 1] and the Curve is a
#    step — more than half the module and the lineage turns aerobic, below and it turns back.
#    A threshold is written here, in the mapping, not as a separate kind of driver.
comp = g.completion("aerobic")
metabolism = simulate_discrete(
    ct, states=["anaerobic", "aerobic"], start="anaerobic", seed=2,
    switch={"anaerobic->aerobic": PerLineage(0.3).scaled_by(comp, Curve(lambda f: 20.0 if f > 0.5 else 1.0)),
            "aerobic->anaerobic": PerLineage(0.3).scaled_by(comp, Curve(lambda f: 1.0 if f > 0.5 else 20.0))})

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
# every tip above half the module ends aerobic, every tip below it anaerobic'''


CONDITIONING = [
    Example("genome_reduction", "Genome reduction",
            "A lifestyle trait drives gene loss. Endosymbionts shed genes faster and gain fewer, so "
            "their genomes end up a fifth the size.",
            "trait → loss", genome_reduction, code=_C_REDUCTION),
    Example("genome_expansion", "Genome expansion",
            "A selection trait drives duplication. Under relaxed selection duplicates accumulate and "
            "the genomes grow fivefold.",
            "trait → duplication", genome_expansion, code=_C_EXPANSION),
    Example("hgt_uptake", "HGT uptake by competence",
            "A competence trait drives <b>who receives</b> a transfer rather than a rate. Competent "
            "lineages take up DNA more often.",
            "trait → transfer uptake", hgt_uptake, code=_C_UPTAKE),
    Example("continuous_conditioning", "A continuous driver",
            "A diffusing trait drives gene gain. A <code>Curve</code> turns each value into a factor, "
            "so genome size follows the trait.",
            "continuous trait → origination", continuous_conditioning, code=_C_CONTINUOUS),
    Example("curve_saturating", "A saturating curve",
            "The same run with a different <code>Curve</code>. Gene gain rises with the trait and "
            "then levels off, so the factor is bounded.",
            "continuous trait → origination", curve_saturating, code=_C_SATURATING),
    Example("curve_optimum", "A curve of your own",
            "The response is any function you write. Here it peaks at an intermediate value, which no "
            "table of per-state multipliers can express.",
            "continuous trait → origination", curve_optimum, code=_C_OPTIMUM),
    Example("trait_drives_trait", "One trait drives another",
            "A temperature trait is grown first; body size then diffuses at a rate that reads it. The "
            "scale is centred on where it started, so white means it has not moved.",
            "trait → trait", trait_drives_trait, code=_C_TRAIT_TRAIT),
    Example("gene_drives_trait", "A gene drives a trait",
            "Carrying a toxin family makes a lineage <i>become</i> pathogenic forty times faster. 74% of "
            "the tips with the gene end up pathogenic, against none of those without.",
            "gene → trait", gene_drives_trait, code=_C_GENE_TRAIT),
    Example("module_drives_metabolism", "A module, through a step",
            "How much of a four-gene module a lineage keeps decides its metabolism, through a step: "
            "<code>lambda f: 20.0 if f > 0.5 else 1.0</code>.",
            "module → trait", module_drives_metabolism, code=_C_MODULE),
]

_C_KEY = """### gene content drives speciation — the tree is an OUTPUT
from zombi2 import genomes, joint
from zombi2.params import PerLineage

r = joint.simulate_joint(
    birth=PerLineage(0.7).scaled_by("genomes:toxin", {"present": 5.0, "absent": 1.0}),
    death=0.15,
    genome=genomes.family(origination=0.15, loss=0.3,
                          family_names=["toxin"], initial_families=5),
    n_extant=70, seed=3)

# loss is twice origination, so the family is losing ground on its own — on a tree
# that does not feel it, it ends up on ~9% of tips. Here it reaches 83%, because the
# lineages carrying it split five times as often.
r.species        # the grown tree
r.genome         # the gene content grown with it"""


JOINING = [
    Example("key_innovation", "A gene drives the splitting",
            "Gene content drives speciation, so the tree is an <b>output</b>. Loss runs at twice "
            "origination, yet the family reaches 83% of tips — its carriers split five times as often.",
            "gene content → speciation", key_innovation, code=_C_KEY),
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
