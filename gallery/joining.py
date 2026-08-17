"""Two ways one evolved value drives another, through the same driver mechanism.

**Conditioning** grows the driver first and holds it fixed: a trait is simulated on the tree, then a
second run reads it, driving a genome rate, another trait's rate, or which lineage receives a
transfer. **Joining** grows both at once, because
the trait drives speciation or extinction and so shapes the tree it is evolving on
(``joint.simulate``). The two lists below feed two gallery sections.
"""

from __future__ import annotations

import collections
import math

import matplotlib.pyplot as plt
from matplotlib import cm, colors

import helpers as h
from helpers import Example

import phylustrator as ph
from zombi2 import species
from zombi2.params import Curve, PerCopy, PerLineage, Recipients
from zombi2 import joint, traits
from zombi2.genomes import family, genome as genomes_spec
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
_IS1J = {"absent": "#b9bec4", "present": "#C25A3C"}         # the element that drives its own genome
_SIZE = {"small": "#b9bec4", "large": "#8C5E8B"}            # the other half of the trait loop
_CAVE = {"surface": "#2E8B6F", "cave": "#4A4A6A"}           # a habitat and a genome, each driving the other
_CLIMATE = {"cold": "#3A7CA5", "hot": "#E4572E"}            # a habitat and a gene's sequence, likewise

#: The two colours the joint figures mark events in, on the tree and in the diagram's key alike: a
#: gene family's presence is drawn in ink, a loss in red. They are constants because the same value
#: has to reach both, and a glyph that does not match its key is worse than no key.
_INK, _LOSS = "#1a1a1a", "#c1443c"
_BAR = "#6a9bd8"                                            # Phylustrator's default bar colour


def _switch_events(ct, trait):
    """A discrete trait's switches as branch events: ``to <state>``, on the branch, at its own time.

    ``initial`` is dropped — it is the root's starting state, not a change — and so is any lineage
    the tree does not name. The kind carries the state switched **to**, which is what colours the
    triangle: a marker in the colour a branch is about to become reads as the moment it turned.
    """
    lab = ct.labels()
    return [{"kind": f"to {c.to_state}", "node": lab[c.lineage], "x": float(c.time)}
            for c in trait.events if c.kind != "initial" and lab.get(c.lineage)]


def _switch_styles(palette):
    return {f"to {state}": ("triangle_right", colour) for state, colour in palette.items()}


def _switch_key(palette, order):
    """The diagram key that matches :func:`_switch_styles` — one triangle per destination state.

    ``order`` is given rather than taken from the palette because a key reads best with the state the
    model is *about* first: the cave before the surface, hot before cold."""
    return [("triangle", palette[s], f"to {s}") for s in order]


def _presence_events(ct, genome, name):
    """A named family's presence at the INTERNAL nodes: a filled dot where it is, a ring where it is not.

    Ancestors only. Every tip would be marked too, and fifty markers stacked against the tip labels
    say nothing the bars beside them do not already say."""
    lab, tips = ct.labels(), set(ct.extant_leaves())
    return [{"kind": f"{name} present" if genome.has_family(n, name) else f"{name} absent",
             "node": lab[n], "x": ct.nodes[n].end_time}
            for n in ct.nodes if n not in tips and lab.get(n) is not None]


def _loss_events(ct, genome, name):
    """One named family's losses as branch events, for the red cross the diagram's key names."""
    lab, fam = ct.labels(), genome.family_names[name]
    return [{"kind": "loss", "node": lab[e.lineage], "x": float(e.time)}
            for e in genome.edges
            if getattr(e, "kind", None) == "loss" and e.family == fam and lab.get(e.lineage)]


def _style():
    return ph.Style(width=1300, height=1000, margin=82, branch_width=3.4)


def _state_history(ct, trait):
    lab = ct.labels()                   # 'n<id>', or 'e<id>' for a lineage that went extinct
    return {lab[i]: segs for i, segs in trait.history.items()}     # per-branch (state, duration)


def _history(r):
    return _state_history(r.complete_tree, r.trait)


# --- joining: the trait and the tree grow together ------------------------------------------

def bisse(out):
    r = joint.simulate(species.birth_death(birth=PerLineage(1.0).scaled_by("trait", {"fast": 2.6, "slow": 0.7}), n_extant=70), traits.discrete(states=["fast", "slow"], switch=0.35), seed=3)
    tree_png = out.replace(".png", "_tree.png")
    (ph.trees.plot(ph.trees.loads(r.complete_tree.to_newick()), style=_style(), skeleton=False)
     + ph.trees.color_history(_history(r), palette=_BISSE)
     + ph.trees.time_axis("time", tick_size=22, label_size=28)).save(tree_png)
    h.composite_markov(tree_png, out, lambda ax: h.draw_markov(
        ax, ["fast", "slow"], _BISSE, {"fast": 2.6, "slow": 0.7}, symbol="λ"),
        loc=(0.02, 0.04, 0.34, 0.30))


_QUASSE_SLOPE = 0.5      # λ = 0.4·e^{0.5·size}: the size doubles the rate every 1.4 units
_QUASSE_BASE = 0.4
_QUASSE_STEP = 0.05


def _draw_response(ax, lo, hi):
    """The rate the colours buy: λ against body size, painted in the tree's own colours.

    A continuous driver has no states to draw a Markov chain between, which is what the other joint
    cards put here. What replaces it is the response curve, in the same colours as the branches, so
    it doubles as the key: find a branch's colour on the curve and read across to its rate."""
    from matplotlib.collections import LineCollection
    from matplotlib.ticker import MaxNLocator

    ax.set_facecolor("#ffffff")
    xs = [lo + (hi - lo) * k / 240 for k in range(241)]
    ys = [_QUASSE_BASE * math.exp(_QUASSE_SLOPE * x) for x in xs]
    ax.add_collection(LineCollection(
        [[(xs[k], ys[k]), (xs[k + 1], ys[k + 1])] for k in range(len(xs) - 1)],
        cmap="viridis", norm=colors.Normalize(lo, hi), array=xs[:-1], linewidth=3.4,
        capstyle="round"))
    ax.set_xlim(lo, hi)
    ax.set_ylim(0.0, ys[-1] * 1.08)
    ax.set_xlabel("body size", fontsize=10.5, labelpad=2)
    ax.set_ylabel("speciation rate", fontsize=10.5, labelpad=2)
    ax.tick_params(labelsize=9, length=3, pad=2)
    ax.xaxis.set_major_locator(MaxNLocator(4))
    ax.yaxis.set_major_locator(MaxNLocator(4))
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)


def quasse(out):
    """A **continuously diffusing** trait driving speciation — QuaSSE, and the one joint model here
    that does not race exactly.

    Body size diffuses down every lineage, and a lineage splits at a rate that rises with it. Neither
    can be simulated first: the size needs a tree to diffuse on, and the tree's shape is what the size
    decides. So one run does both, and the tree is an output.

    The tree is the whole picture, and what to look at is which lineages left descendants. The
    yellow clade in the middle is a radiation: once a lineage was large it split fast, and its
    daughters started large. The blue and purple lineages above and below it are the ones the trait
    never favoured, and they sit on long branches with almost no daughters. That bias is what makes
    this joint rather than a trait painted on a tree it did not shape — a diffusion has no direction
    of its own, so on a tree it does not touch the tips would sit around the root value.

    **This one slices.** Every other driver in this section changes only at events, and an event ends
    the Gillespie step, so the rate is constant in between. A diffusion moves at every instant. The
    run holds each lineage's size fixed across a step of 0.05 and releases it at the boundary, where
    the exact transition law applies — so the trait is exact and only its grip on speciation is
    approximated. Halve the step, rerun the same seed, and see whether the answer moves.
    """
    r = joint.simulate(
        species.birth_death(
            birth=PerLineage(_QUASSE_BASE).scaled_by(
                "trait", Curve(lambda x: math.exp(_QUASSE_SLOPE * x)), step=_QUASSE_STEP),
            death=0.05, n_extant=70),
        traits.continuous(start=0.0, rate=1.0), seed=20)
    ct = r.complete_tree
    lab = ct.labels()
    vals = {lab[i]: r.trait.node_values[i] for i in ct.nodes}
    limits = (math.floor(min(vals.values())), math.ceil(max(vals.values())))
    tree_png = out.replace(".png", "_tree.png")
    (ph.trees.plot(ph.trees.loads(ct.to_newick()), style=_style(), skeleton=False)
     # the same limits the key spans, so the two say the same thing about a colour
     + ph.trees.color_branches(vals, cmap="viridis", limits=limits)
     + ph.trees.time_axis("time", tick_size=22, label_size=28)).save(tree_png)
    # the key goes in the tree's empty upper-left, on a card: it is a plot with axes of its own, and
    # unframed beside the tree's time axis the two read as one set of coordinates
    h.composite_markov(tree_png, out, lambda ax: _draw_response(ax, *limits),
                       loc=(0.091, 0.787, 0.147, 0.12), keep_axes=True,
                       panel=(0.041, 0.047, 0.016, 0.019))


def classe(out):
    """The state drives the split, and the split changes the state — change at the fork, not along it.

    BiSSE's trait wanders down the branches and the tree feels it. Here almost nothing happens along a
    branch (`switch=0.08`): the state changes **at** speciations, one daughter taking a new one, which
    is what the field calls a cladogenetic model. The squares mark the splits where it happened, and
    they are on the tree the same trait was shaping — which is why this is joint and not a trait
    painted on a tree that was finished first."""
    r = joint.simulate(
        species.birth_death(birth=PerLineage(1.0).scaled_by("trait", {"fast": 2.6, "slow": 0.7}),
                            n_extant=70),
        traits.discrete(states=["fast", "slow"], switch=0.08, at_speciation=0.15), seed=3)
    label = r.complete_tree.labels()
    at_split = [{"kind": "change at the split", "node": label[e.lineage], "x": e.time}
                for e in r.trait.events if e.kind == "on_speciation"]
    tree_png = out.replace(".png", "_tree.png")
    (ph.trees.plot(ph.trees.loads(r.complete_tree.to_newick()), style=_style(), skeleton=False)
     + ph.trees.color_history(_history(r), palette=_BISSE)
     + ph.trees.branch_events(at_split, size=5.5, legend_title="", legend_loc="top-right",
                              legend_size=20,
                              styles={"change at the split": ("square", "#111111")})
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
    r = joint.simulate(species.birth_death(birth=PerLineage(0.7).scaled_by("genomes:toxin", {"present": 5.0, "absent": 1.0}), death=0.15, n_extant=70), genomes_spec(origination=0.15, loss=0.3, families=[family("toxin")],
                            initial_families=5), seed=3)
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


def mobile_element_joint(out):
    """A gene family driving the rest of its OWN genome, in one run — the level joined to itself.

    An insertion sequence is planted on one branch. Where it is present the whole genome donates
    genes thirty times as often, and the element also moves itself, so it jumps into lineages that
    never inherited it and makes them donors too. That is the loop: the element decides the transfer
    rate, and the transfer rate is how the element spreads.

    Nothing here can be simulated first. In the conditioning version of this figure the element lives
    in a *separate* genome run, finished before the run it drives — so it could never be moved by the
    transfer it caused, and the genes it mobilised were not in the same organism. One run removes
    both of those, and the driver is a live name rather than a finished result.

    The bars count transfers **given**, since a driven ``transfer`` says how often a lineage donates.
    The 13 carriers give away 4.5 genes each against 0.35 for the other 17.
    """
    ct = simulate_species_tree(birth=1.0, n_extant=30, seed=4).complete_tree
    g = simulate_genomes_family(
        ct, initial_families=25, duplication=0.05, loss=0.12, seed=7, joint=True, max_family_size=8,
        families=[family("IS1", origin=(11, None), transfer=PerCopy(0.30), loss=0.08)],
        transfer=PerCopy(0.025).scaled_by("genomes:IS1", {"present": 30.0, "absent": 1.0}))

    is1 = g.family_names["IS1"]
    lab = ct.labels()
    # every transfer writes two edges; take the recipient-side one so each is counted once, and group
    # it by the branch that gave the copy away. The element's own moves are left out — the claim is
    # about the genes it mobilises, not about itself.
    donated = collections.Counter(e.donor for e in g.edges
                                  if e.kind == "transfer" and e.recipient is not None
                                  and e.family != is1)
    tips = list(ct.extant_leaves())
    given = {lab[n]: donated.get(n, 0) for n in sorted(tips)}
    carriers = {lab[n] for n in tips if g.family_counts(n)[is1] > 0}
    tipcol = {n: _IS1J["present" if n in carriers else "absent"] for n in given}
    history = {lab[i]: segs for i, segs in g.presence("IS1").history(ct).items()}
    fig = (ph.trees.plot(ph.trees.loads(ct.to_newick()), skeleton=False,
                         style=ph.Style(width=900, height=900, margin=92, branch_width=3.0))
           + ph.trees.color_history(history, palette=_IS1J)
           + ph.trees.time_axis("time", tick_size=20, label_size=26, bold=False))
    real = out.replace(".png", "_real.png")
    ph.beside(fig, ph.genomes.bars(given, colors=tipcol, label="transfers given",
                                   tick_size=20, label_size=26),
              width=1150, tree_fraction=0.58, footer=36).save(real)
    # The same inset the other joining cards use: the two states a lineage can be in, and the rate
    # each one transfers at. The arrows are deliberately unlabelled — a lineage leaves `present` by
    # losing its last copy and enters it by receiving one from elsewhere, and neither of those is a
    # rate this lineage carries, so putting a number on them would be inventing one.
    h.composite_markov(real, out, lambda ax: h.draw_markov(
        ax, ["present", "absent"], _IS1J, {"present": 0.75, "absent": 0.025}, symbol="τ"),
        loc=(0.005, 0.665, 0.33, 0.30))


def cave_genomes(out):
    """A trait and a genome driving **each other**, on a tree the run is handed.

    Two statements point opposite ways. A lineage in the cave loses genes five times as fast, so the
    habitat drives the genome. And a lineage that has lost its eye gene turns to the cave twenty-five
    times as readily, so the genome drives the habitat. Neither can be simulated first: to grow the
    habitat you would need to know whether the eye is still there, and whether the eye survives
    depends on the habitat.

    This is the first joint model whose tree is an **input**. The two that drive speciation produce
    the tree; here it is handed over like any other input, and what comes back is two levels' results
    from one run.

    Twenty-one of the fifty tips end up underground, and their genomes are half the size: twenty-seven
    genes at the median cave tip against fifty-three on the surface. The eye is gone from every cave
    tip — the twenty-three tips that still carry it are all on the surface. With the returning arrow
    taken out — the habitat grown on its own, which is the only way the model can be written as two
    runs in order — a median of one tip in fifty goes underground, over twenty seeds.

    The tree carries the whole story, and the diagram above it is the key. Dots and rings mark the
    ancestors that still had the eye and the ones that had lost it; the six red crosses are where it
    went; the triangles are the ten switches to the cave, and nine of the ten hang off a branch whose
    ancestor is already a ring.
    """
    import matplotlib.image as mpimg

    ct = simulate_species_tree(birth=1.0, n_extant=50, seed=2).complete_tree
    r = joint.simulate(
        genomes_spec(duplication=0.05, origination=12.0, initial_families=60,
                     loss=PerCopy(0.30).scaled_by("trait", {"cave": 5.0, "surface": 1.0}),
                     families=[family("eye")]),
        traits.discrete(states=["surface", "cave"], start="surface",
                        switch={"surface->cave": PerLineage(0.02).scaled_by(
                                    "genomes:eye", {"present": 1.0, "absent": 25.0}),
                                "cave->surface": 0.10}),
        tree=ct, seed=2)

    lab, tips = ct.labels(), sorted(ct.extant_leaves())
    sizes = {lab[n]: len(r.genome.genomes[lab[n]]) for n in tips}
    tipcol = {lab[n]: _CAVE[r.trait.values[lab[n]]] for n in tips}
    # three annotation layers, no legend on any of them: the diagram's boxes carry the key
    fig = (ph.trees.plot(ph.trees.loads(ct.to_newick()), skeleton=False,
                         style=ph.Style(width=900, height=900, margin=92, branch_width=3.0))
           + ph.trees.color_history(_state_history(ct, r.trait), palette=_CAVE)
           + ph.trees.branch_events(_presence_events(ct, r.genome, "eye"),
                                    styles={"eye present": ("circle", _INK),
                                            "eye absent": ("ring", _INK)},
                                    size=3.8, legend=False)
           + ph.trees.branch_events(_loss_events(ct, r.genome, "eye"),
                                    styles={"loss": ("cross", _LOSS)}, size=7.5, legend=False)
           + ph.trees.branch_events(_switch_events(ct, r.trait), styles=_switch_styles(_CAVE),
                                    size=6.0, legend=False)
           + ph.trees.time_axis("time", tick_size=20, label_size=26, bold=False))
    real = out.replace(".png", "_real.png")
    ph.beside(fig, ph.genomes.bars(sizes, colors=tipcol, label="genome size (genes)",
                                   tick_size=20, label_size=26),
              width=1150, tree_fraction=0.58, footer=36).save(real)

    # one row per statement, each row reading driver · sentence · the rate it drives
    diag = h.joint_png(out.replace(".png", "_diag.png"), [
        (("traits", "habitat", [("swatch", _CAVE["surface"], "surface"),
                                ("swatch", _CAVE["cave"], "cave")]),
         ("genomes", "loss rate", [("cross", _LOSS, "a gene is lost")]),
         "in the cave, genes are lost 5× faster"),
        (("genomes", "the eye family", [("circle", _INK, "present"), ("ring", _INK, "absent")]),
         ("traits", "switch rate", _switch_key(_CAVE, ["cave", "surface"])),
         "with no eye, turns cave 25× faster"),
    ])
    fig2 = plt.figure(figsize=(12, 9.1))
    axr = fig2.add_axes([0.0, 0.0, 1.0, 0.79]); axr.imshow(mpimg.imread(real)); axr.set_axis_off()
    axd = fig2.add_axes([0.14, 0.762, 0.72, 0.185]); axd.imshow(mpimg.imread(diag)); axd.set_axis_off()
    fig2.savefig(out, dpi=140, bbox_inches="tight")
    plt.close(fig2)


_SEQ_LOOP_STEP = 0.05
_SEQ_LOOP_RESPONSE = (0.10, 25.0)      # factor = 0.10 + 25·(the other gene's KR share)


def _kr_poor():
    """LG's chemistry over KR-depleted frequencies.

    Both genes are founded here rather than at LG's own equilibrium, and that is the whole reason the
    figure has anything to show: a gene founded at its model's equilibrium sits there, so its
    composition never moves and a rate reading it reads a constant. Founded away from it, a gene
    **ameliorates** toward it, which is a real thing a horizontally acquired gene does — and while it
    does, its composition is a quantity with somewhere to go."""
    import numpy as np
    from zombi2.sequences.substitution_models import AMINO_ACIDS, lg, reversible

    m = lg()
    S = m.Q / m.stationary[None, :]
    S = (S + S.T) / 2.0
    np.fill_diagonal(S, 0.0)
    pi = m.stationary.copy()
    pi[[AMINO_ACIDS.index(c) for c in "KR"]] *= 0.12
    return reversible(S, pi / pi.sum(), name="KR-poor LG", alphabet=AMINO_ACIDS)


def _kr_by_node(result, fam, letters="KR"):
    """``{node label: the share of that copy's sequence that is ``letters``}``, over every node of one
    family's gene tree, pooled onto the species lineage each copy sits on."""
    hits, total = collections.Counter(), collections.Counter()
    for table in (result.alignments[fam], result.ancestral[fam]):
        for label, seq in table.items():
            node = label.split("_", 1)[0]
            hits[node] += sum(seq.count(c) for c in letters)
            total[node] += len(seq)
    return {node: hits[node] / total[node] for node in total}


def genome_and_sequence(out):
    """A GENOME and a gene's SEQUENCE, each driving the other — the last cell of the map.

    `hisA` arrives AT-rich, as a horizontally acquired gene does, and ameliorates toward the
    composition its host's model settles on. While it is still AT-rich its lineage loses genes fast.
    So how far one gene has got decides how much of the whole genome survives, and the genome decides
    which copies of that gene exist to be counted.

    Neither can be simulated first. To evolve the sequences you need to know which copies exist and
    when. To know which copies exist you need the loss rate, which reads the sequences.

    The tree is the only thing handed to the run. The gene trees are not: they come out of the genome
    participant, which is what makes this different from every conditioned run at this level.

    One panel: how far `hisA` has ameliorated, and beside it how many genes each lineage has left.
    The pale clades are the ameliorated ones, and they are the full genomes. In this run `hisA`
    survives at all thirty tips, so no branch is left grey for `absent=0.35` to answer for, and the
    twenty-two losses the driven rate fired are all of other families.
    """
    from zombi2.genomes import genome as genomes_spec_
    from zombi2.sequences import composition as offers_composition
    from zombi2.sequences import gene as gene_spec
    from zombi2.sequences import hky85

    ct = simulate_species_tree(birth=1.0, n_extant=30, seed=1).complete_tree
    at_rich = hky85(2.0, frequencies=(0.40, 0.10, 0.10, 0.40))
    r = joint.simulate(
        genomes_spec_(duplication=0.15, origination=0.05, initial_families=25,
                      loss=PerCopy(0.15).scaled_by(
                          "sequences:hisA", Curve(lambda gc: 30.0 ** ((0.35 - gc) / 0.2)),
                          step=0.05),
                      families=[family("hisA")]),
        gene_spec(name="hisA", model=hky85(2.0), length=250, start=at_rich, substitution=0.8,
                  offers=offers_composition("GC", absent=0.35)),
        tree=ct, seed=3)

    lab = ct.labels()
    fam = r.genome.family_names["hisA"]
    gc = _kr_by_node(r.sequences, fam, letters="GC")
    span = (min(gc.values()), max(gc.values()))
    sizes = {lab[n]: len(r.genome.genomes[lab[n]]) for n in sorted(ct.extant_leaves())}
    # a lineage that carries no copy of hisA keeps the default colour, so the grey branches are
    # exactly the ones `absent=` had to answer for
    style = _panel_style()
    style = ph.Style(width=style.width, height=style.height, margin=style.margin,
                     branch_width=style.branch_width, branch_color="#c3c8cc")
    # no loss crosses here, unlike cave_genomes: `hisA` itself is never lost in this run, and the
    # twenty-two losses of everything else pile up half a dozen deep on the stem, where they read as
    # one red smudge rather than as events. The bars are what the loss rate did, and they say it
    fig = (ph.trees.plot(ph.trees.loads(ct.to_newick()), skeleton=False, style=style)
           + ph.trees.color_branches(gc, cmap="magma_dark", limits=span)
           + ph.trees.time_axis("time", tick_size=20, label_size=26, bold=False))
    png = out.replace(".png", "_panel.png")
    ph.beside(fig, ph.genomes.bars(sizes, label="genes left", tick_size=20, label_size=26),
              width=1150, tree_fraction=0.58, footer=36).save(png)
    ramp = [("gradient", "magma_dark", ("AT-rich", "at the model's own"))]
    diag = h.joint_png(out.replace(".png", "_diag.png"), [
        (("genomes", "the genome", [("swatch", _BAR, "genes left")]),
         ("sequences", "hisA", ramp),
         "the genome says which copies of hisA exist"),
        (("sequences", "hisA", ramp),
         ("genomes", "loss rate", [("cross", _LOSS, "a gene is lost")]),
         "while hisA is AT-rich, genes are lost 30× faster"),
    ])
    h.composite_under_diagram(out, diag, [(png, "hisA, and what each tip has left")],
                              diagram_frac=0.72)


def sequence_loop(out):
    """Two genes, each one's substitution rate reading the other's composition — the **sequence**
    level joined to itself, and the last of the three loops.

    Both genes arrive KR-poor and ameliorate toward LG's own composition. What each one's rate reads
    is how far the *other* has got: a gene beside a KR-poor partner runs at a tenth of its base rate,
    and one beside a partner at LG's equilibrium runs at three times it. So the pair is slow while
    both are still poor and accelerates as neither is — a loop that pushes rather than settles.

    Neither can be simulated first. To evolve the first gene you would need the second's composition
    at every instant of every branch, and that composition is what the first gene's own rate decides.

    **This one slices.** A composition moves with every substitution, so there is no interval where
    either rate holds still. The run holds both fixed across a step of 0.05 of species time and
    releases them at the boundary; inside a slice the transition matrix is the ordinary one. Against
    the same run with the response flattened to 1.0, the loop carries the pair measurably further —
    every seed, not on average.

    The same tree twice, painted by each gene's share of K and R. What the loop produces is that the
    two panels move together: neither gene ameliorates far while the other is still behind.
    """
    from zombi2.sequences import composition as offers_composition
    from zombi2.sequences import gene as gene_spec
    from zombi2.sequences import lg, simulate_sequences
    from zombi2.params import PerSite

    ct = simulate_species_tree(birth=1.0, n_extant=40, seed=1).complete_tree
    g = simulate_genomes_family(ct, initial_families=4, duplication=0.0, loss=0.0, origination=0.0,
                                seed=2, families=[family("hisA"), family("hisF")])
    a, b = _SEQ_LOOP_RESPONSE
    start = _kr_poor()
    r = simulate_sequences(g, joint=True, seed=3, genes=[
        gene_spec(name=name, model=lg(), length=250, start=start,
                  offers=offers_composition("KR", absent=0.02),
                  substitution=PerSite(0.5).scaled_by(f"sequences:{other}",
                                                      Curve(lambda x: a + b * x),
                                                      step=_SEQ_LOOP_STEP))
        for name, other in (("hisA", "hisF"), ("hisF", "hisA"))])

    shares = {name: _kr_by_node(r, g.family_names[name]) for name in ("hisA", "hisF")}
    span = (min(min(v.values()) for v in shares.values()),
            max(max(v.values()) for v in shares.values()))
    style = _panel_style()
    pngs = []
    for name in ("hisA", "hisF"):
        png = out.replace(".png", f"_{name}.png")
        (ph.trees.plot(ph.trees.loads(ct.to_newick()), skeleton=False, style=style)
         + ph.trees.color_branches(shares[name], cmap="magma_dark", limits=span)
         + ph.trees.time_axis("time", tick_size=20, label_size=26, bold=False)).save(png)
        pngs.append(png)
    ramp = [("gradient", "magma_dark", ("KR-poor", "at LG's own"))]
    diag = h.joint_png(out.replace(".png", "_diag.png"), [
        (("sequences", "hisA", ramp),
         ("sequences", "hisF's rate", [("word", None, "a site changes")]),
         "the richer hisA is, the faster hisF evolves"),
        (("sequences", "hisF", ramp),
         ("sequences", "hisA's rate", [("word", None, "a site changes")]),
         "the richer hisF is, the faster hisA evolves"),
    ])
    # the ramp is named once per gene in the diagram, so the rows below carry no key of their own
    h.composite_under_diagram(out, diag, [(pngs[0], "hisA"), (pngs[1], "hisF")],
                              diagram_frac=0.72)


def trait_and_sequence(out):
    """A trait and a gene's SEQUENCE, each driving the other — the cross-level join whose two ends
    are the furthest apart, and the last cell of the map to be built.

    Both directions already ran as conditioning, one card each. A trait drives a gene's substitution
    rate (Co8), and that gene's composition drives how fast a trait switches (Co19). Written at once
    they are a cycle: to evolve the sequence you would need the trait on every stretch of every
    branch, and the trait's own rate is what the sequence decides.

    The model here closes that loop in one direction. `rpoB` arrives KR-poor and ameliorates toward
    LG's own composition; the richer it gets, the more readily its lineage turns hot; and a hot
    lineage substitutes four times faster, which is what carries `rpoB` further. Against the same run
    with the trait's response flattened, four-fifths of the tree ends hot rather than three-fifths,
    over eight seeds.

    Species time is sliced. Inside a slice the trait's own Gillespie is exact — switches mid-slice
    included — and the gene's branch length is the trait's factor **integrated** across those
    switches, drawn once, which is exact too. What is approximated is only the composition the trait
    reads: it belongs to the top of the slice rather than to each instant.

    The tree twice: the habitat the trait walked, and the same tree painted by how far `rpoB` has
    got. The hot clades are the ameliorated ones, in both directions at once.
    """
    from zombi2.sequences import composition as offers_composition
    from zombi2.sequences import gene as gene_spec
    from zombi2.sequences import lg
    from zombi2.params import PerSite

    ct = simulate_species_tree(birth=1.0, n_extant=30, seed=1).complete_tree
    g = simulate_genomes_family(ct, initial_families=3, duplication=0.0, loss=0.0, origination=0.0,
                                seed=2, families=[family("rpoB")])
    r = joint.simulate(
        traits.discrete(name="habitat", states=["cold", "hot"], start="cold",
                        switch={"cold->hot": PerLineage(0.5).scaled_by(
                                    "sequences:rpoB", Curve(lambda x: 0.05 + 30.0 * x), step=0.05),
                                "hot->cold": 0.3}),
        gene_spec(name="rpoB", model=lg(), length=250, start=_kr_poor(),
                  offers=offers_composition("KR", absent=0.02),
                  substitution=PerSite(0.6).scaled_by("trait", {"hot": 4.0, "cold": 1.0})),
        genomes=g, seed=5)

    lab = ct.labels()
    fam = g.family_names["rpoB"]
    share = _kr_by_node(r.sequences, fam)
    span = (min(share.values()), max(share.values()))
    style = _panel_style()
    habitat = out.replace(".png", "_habitat.png")
    (ph.trees.plot(ph.trees.loads(ct.to_newick()), skeleton=False, style=style)
     + ph.trees.color_history({lab[i]: segs for i, segs in r.trait.history.items()},
                              palette=_CLIMATE)
     + ph.trees.branch_events(_switch_events(ct, r.trait), styles=_switch_styles(_CLIMATE),
                              size=6.0, legend=False)
     + ph.trees.time_axis("time", tick_size=20, label_size=26, bold=False)).save(habitat)
    gene_png = out.replace(".png", "_rpoB.png")
    (ph.trees.plot(ph.trees.loads(ct.to_newick()), skeleton=False, style=style)
     + ph.trees.color_branches(share, cmap="magma_dark", limits=span)
     + ph.trees.time_axis("time", tick_size=20, label_size=26, bold=False)).save(gene_png)
    diag = h.joint_png(out.replace(".png", "_diag.png"), [
        (("traits", "habitat", [("swatch", _CLIMATE["cold"], "cold"),
                                ("swatch", _CLIMATE["hot"], "hot")]),
         ("sequences", "rpoB's rate", [("word", None, "a site changes")]),
         "hot lineages substitute 4× faster"),
        (("sequences", "rpoB", [("gradient", "magma_dark", ("KR-poor", "at LG's own"))]),
         ("traits", "switch rate", _switch_key(_CLIMATE, ["hot", "cold"])),
         "the richer rpoB is, the readier the switch to hot"),
    ])
    h.composite_under_diagram(out, diag, [(habitat, "habitat"), (gene_png, "rpoB")],
                              diagram_frac=0.72)


def trait_loop(out):
    """Two traits, each reading the other, on one tree — the trait level joined to itself.

    Body size sets how readily a lineage goes underground; living underground sets how readily it
    grows. Neither can be simulated first, so one run does both.

    This one is **exact**. Two traits that read each other are one Markov chain over the pairs of
    their states, so the pair has an ordinary generator and the same branch walk a single trait takes
    runs it — nothing thinned, nothing approximated.

    The same tree twice, painted by each trait, with a triangle at every switch in the colour of the
    state it turned into. What the loop produces is the alignment between the two panels, and it is a
    statement neither trait makes on its own: twelve of the twenty-one cave tips are large, against
    one of the nineteen on the surface.
    """
    ct = simulate_species_tree(birth=1.0, n_extant=40, seed=4).complete_tree
    r = traits.simulate_traits(ct, [
        traits.discrete(name="habitat", states=["surface", "cave"], start="surface",
                        switch={"surface->cave": PerLineage(0.05).scaled_by(
                                    "traits:size", {"small": 1.0, "large": 8.0}),
                                "cave->surface": 0.1}),
        traits.discrete(name="size", states=["small", "large"], start="small",
                        switch={"small->large": PerLineage(0.05).scaled_by(
                                    "traits:habitat", {"surface": 1.0, "cave": 6.0}),
                                "large->small": 0.1})],
        joint=True, seed=16)

    style = _panel_style()
    pngs = []
    for k, (name, palette) in enumerate((("habitat", _CAVE), ("size", _SIZE))):
        png = out.replace(".png", f"_t{k}.png")
        (ph.trees.plot(ph.trees.loads(ct.to_newick()), skeleton=False, style=style)
         + ph.trees.color_history(_state_history(ct, r[name]), palette=palette)
         + ph.trees.branch_events(_switch_events(ct, r[name]), styles=_switch_styles(palette),
                                  size=6.0, legend=False)
         + ph.trees.time_axis("time", tick_size=20, label_size=26, bold=False)).save(png)
        pngs.append(png)
    diag = h.joint_png(out.replace(".png", "_diag.png"), [
        (("traits", "body size", [("swatch", _SIZE["small"], "small"),
                                  ("swatch", _SIZE["large"], "large")]),
         ("traits", "habitat switch", _switch_key(_CAVE, ["cave", "surface"])),
         "when large, turns cave 8× more readily"),
        (("traits", "habitat", [("swatch", _CAVE["surface"], "surface"),
                                ("swatch", _CAVE["cave"], "cave")]),
         ("traits", "size switch", _switch_key(_SIZE, ["large", "small"])),
         "in the cave, grows 6× more readily"),
    ])
    # no colour key under either row label: the diagram's boxes carry both palettes
    h.composite_under_diagram(out, diag, [(pngs[0], "habitat"), (pngs[1], "body size")],
                              diagram_frac=0.72)


def state_extinction(out):
    r = joint.simulate(species.birth_death(birth=1.0, death=PerLineage(1.0).scaled_by("trait", {"doomed": 0.75, "safe": 0.05}), n_extant=35), traits.discrete(states=["doomed", "safe"], switch=0.3), seed=1)
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
    r = joint.simulate(species.birth_death(birth=PerLineage(1.0).scaled_by("trait", {"slow": 0.6, "medium": 1.3, "fast": 2.6}), death=0.4, n_extant=50), traits.discrete(states=["slow", "medium", "fast"], switch=0.3), seed=2)
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
    g = simulate_genomes_family(ct, initial_families=20, families=[family("tox")],
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
from zombi2 import joint, species, traits
from zombi2.params import PerLineage

r = joint.simulate(
    species.birth_death(birth=PerLineage(1.0).scaled_by("trait", {"fast": 2.6, "slow": 0.7}),
                        n_extant=70),
    traits.discrete(states=["fast", "slow"], switch=0.35), seed=3)

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
from zombi2 import joint, species, traits
from zombi2.params import PerLineage

r = joint.simulate(
    species.birth_death(birth=1.0,
                        death=PerLineage(1.0).scaled_by("trait", {"doomed": 0.75, "safe": 0.05}),
                        n_extant=35),
    traits.discrete(states=["doomed", "safe"], switch=0.3), seed=1)
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
from zombi2 import joint, species, traits
from zombi2.params import PerLineage

r = joint.simulate(
    species.birth_death(
        birth=PerLineage(1.0).scaled_by("trait", {"slow": 0.6, "medium": 1.3, "fast": 2.6}),
        death=0.4, n_extant=50),
    traits.discrete(states=["slow", "medium", "fast"], switch=0.3), seed=2)
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
from zombi2.genomes import family, simulate_genomes_family
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
from zombi2.genomes import family, simulate_genomes_family
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
from zombi2.genomes import family, simulate_genomes_family
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
from zombi2.genomes import family, simulate_genomes_family
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
from zombi2.genomes import family, simulate_genomes_family
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
from zombi2.genomes import family, simulate_genomes_family
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
    g = simulate_genomes_family(ct, initial_families=20, duplication=0.08, loss=0.08, seed=11,
                                families=[family(n, module="aerobic") for n in nuo])
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
from zombi2.genomes import family, simulate_genomes_family
from zombi2.traits import simulate_discrete
from zombi2.params import PerLineage

ct = simulate_species_tree(birth=1.0, death=0.2, n_extant=45, seed=4).complete_tree

# 1. the driver: genomes with ONE family named, so it can be referred to. Loss a little above
#    duplication, so the family survives in some clades and decays in others.
g = simulate_genomes_family(ct, initial_families=20, families=[family("tox")],
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
from zombi2.genomes import family, simulate_genomes_family
from zombi2.traits import simulate_discrete
from zombi2.params import Curve, PerLineage

ct = simulate_species_tree(birth=1.0, death=0.2, n_extant=45, seed=4).complete_tree
nuo = [f"nuo{c}" for c in "ABCD"]

# 1. the driver: four families grouped into one module. Loss matched to duplication, so the
#    module survives in about half the tree and decays in the rest.
g = simulate_genomes_family(ct, initial_families=20, duplication=0.08, loss=0.08, seed=11,
                            families=[family(n, module="aerobic") for n in nuo])

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


_C_CLASSE = '''\
### simulate  —  the state drives the split, and the split changes the state
from zombi2 import joint, species, traits
from zombi2.params import PerLineage

r = joint.simulate(
    species.birth_death(birth=PerLineage(1.0).scaled_by("trait", {"fast": 2.6, "slow": 0.7}),
                        n_extant=70),
    traits.discrete(states=["fast", "slow"],
                    switch=0.08,          # almost nothing happens along a branch
                    at_speciation=0.15),  # it happens at the fork instead
    seed=3)

# the squares are the splits where the state changed: r.trait.events, kind "on_speciation"'''


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
from zombi2 import genomes, joint, species
from zombi2.genomes import family
from zombi2.params import PerLineage

r = joint.simulate(
    species.birth_death(
        birth=PerLineage(0.7).scaled_by("genomes:toxin", {"present": 5.0, "absent": 1.0}),
        death=0.15, n_extant=70),
    genomes.genome(origination=0.15, loss=0.3, families=[family("toxin")], initial_families=5),
    seed=3)

# loss is twice origination, so the family is losing ground on its own — on a tree
# that does not feel it, it ends up on ~9% of tips. Here it reaches 83%, because the
# lineages carrying it split five times as often.
r.species        # the tree that came out
r.genome         # the gene content simulated with it"""


_C_QUASSE = """### a CONTINUOUSLY diffusing trait drives speciation — QuaSSE, and the one model that slices
import math
from zombi2 import joint, species, traits
from zombi2.params import Curve, PerLineage

# the driver is a number, so the mapping is a Curve (value -> factor), never a table.
# `step=` is the stretch of time the size is held fixed across: a diffusion moves at
# every instant, so there is no interval where the birth rate holds still on its own.
r = joint.simulate(
    species.birth_death(
        birth=PerLineage(0.4).scaled_by("trait", Curve(lambda x: math.exp(0.5 * x)), step=0.05),
        death=0.05, n_extant=70),
    traits.continuous(start=0.0, rate=1.0), seed=11)

r.species        # the tree that came out
r.trait          # the diffusion simulated with it — an ordinary continuous TraitsResult

### plot  —  the tree painted by the value the birth rate was reading
import phylustrator as ph

ct, lab = r.complete_tree, r.complete_tree.labels()
(ph.trees.plot(ph.trees.loads(ct.to_newick()))
 + ph.trees.color_branches({lab[i]: r.trait.node_values[i] for i in ct.nodes}, cmap="viridis")
 + ph.trees.time_axis("time")).save("quasse.png")
# the check the manual asks for: halve step, rerun the same seed, see whether it moves"""


_C_MOBILE_JOINT = """### a gene family drives the rest of its OWN genome — one run, the level joined to itself
from zombi2.species import simulate_species_tree
from zombi2.genomes import family, simulate_genomes_family
from zombi2.params import PerCopy

ct = simulate_species_tree(birth=1.0, n_extant=30, seed=4).complete_tree

# `joint=True` says the two things drive each other. The driver is a live NAME —
# "genomes:IS1" — not a finished run, because the element is what this run is producing.
g = simulate_genomes_family(
    ct, initial_families=25, duplication=0.05, loss=0.12, seed=7, joint=True, max_family_size=8,
    # the element moves itself, and is planted on one branch rather than at the origin
    families=[family("IS1", origin=(11, None), transfer=PerCopy(0.30), loss=0.08)],
    # ...and where it is present, the whole genome donates 30x as often
    transfer=PerCopy(0.025).scaled_by("genomes:IS1", {"present": 30.0, "absent": 1.0}))

g.presence("IS1")                # where it went, exact and mid-branch
# 13 carriers give away 4.5 genes each; the other 17 give away 0.35

### plot  —  the tree by where the element is, beside transfers GIVEN by each tip
import phylustrator as ph

lab = ct.labels()
history = {lab[i]: segs for i, segs in g.presence("IS1").history(ct).items()}
(ph.trees.plot(ph.trees.loads(ct.to_newick()), skeleton=False)
 + ph.trees.color_history(history, palette={"absent": "#b9bec4", "present": "#C25A3C"})
 + ph.trees.time_axis("time", bold=False)).save("tree.png")
# a driven `transfer` says how often a lineage DONATES, so the bars count transfers given"""


_C_CAVE = """### a trait and a genome driving EACH OTHER, on a tree the run is handed
from zombi2 import genomes, joint, traits
from zombi2.genomes import family
from zombi2.params import PerCopy, PerLineage
from zombi2.species import simulate_species_tree

ct = simulate_species_tree(birth=1.0, n_extant=50, seed=2).complete_tree

# no species.birth_death among the participants, so the tree is an INPUT: pass it with tree=
r = joint.simulate(
    # the cave costs genes — the habitat drives the genome
    genomes.genome(duplication=0.05, origination=12.0, initial_families=60,
                   loss=PerCopy(0.30).scaled_by("trait", {"cave": 5.0, "surface": 1.0}),
                   families=[family("eye")]),
    # ...and losing the eye commits a lineage to the cave — the genome drives the habitat
    traits.discrete(states=["surface", "cave"], start="surface",
                    switch={"surface->cave": PerLineage(0.02).scaled_by(
                                "genomes:eye", {"present": 1.0, "absent": 25.0}),
                            "cave->surface": 0.10}),
    tree=ct, seed=2)

r.trait          # the habitat, as the traits level would have written it
r.genome         # the genome, as the genomes level would have written it
# 27 genes at the median cave tip against 53 on the surface

### plot  —  the tree by habitat, beside genome size at each tip
import phylustrator as ph

lab = ct.labels()
history = {lab[i]: segs for i, segs in r.trait.history.items()}
pal = {"surface": "#2E8B6F", "cave": "#4A4A6A"}
# every switch as a triangle in the colour of the state it turned INTO
switches = [{"kind": f"to {c.to_state}", "node": lab[c.lineage], "x": c.time}
            for c in r.trait.events if c.kind != "initial"]
fig = (ph.trees.plot(ph.trees.loads(ct.to_newick()), skeleton=False)
       + ph.trees.color_history(history, palette=pal)
       + ph.trees.branch_events(switches, legend=False,
                                styles={f"to {s}": ("triangle_right", c) for s, c in pal.items()})
       + ph.trees.time_axis("time", bold=False))
sizes = {lab[n]: len(r.genome.genomes[lab[n]]) for n in sorted(ct.extant_leaves())}
ph.beside(fig, ph.genomes.bars(sizes, label="genome size (genes)")).save("cave.png")"""


_C_TRAIT_LOOP = """### two traits, each reading the other — the trait level joined to itself
from zombi2 import traits
from zombi2.params import PerLineage
from zombi2.species import simulate_species_tree

ct = simulate_species_tree(birth=1.0, n_extant=40, seed=4).complete_tree

# `joint=True`, and each switch rate reads the OTHER trait by name. Neither can be
# simulated first: to grow the habitat you would need the size, and the size needs
# the habitat. One run does both, exactly — the pair is a Markov chain over the
# pairs of their states, so there is nothing to thin and nothing to approximate.
r = traits.simulate_traits(ct, [
    traits.discrete(name="habitat", states=["surface", "cave"], start="surface",
                    switch={"surface->cave": PerLineage(0.05).scaled_by(
                                "traits:size", {"small": 1.0, "large": 8.0}),
                            "cave->surface": 0.1}),
    traits.discrete(name="size", states=["small", "large"], start="small",
                    switch={"small->large": PerLineage(0.05).scaled_by(
                                "traits:habitat", {"surface": 1.0, "cave": 6.0}),
                            "large->small": 0.1})],
    joint=True, seed=16)

r["habitat"], r["size"]    # one ordinary TraitsResult each, keyed by name

### plot  —  the same tree painted by each trait, one above the other
import phylustrator as ph

lab = ct.labels()
for name, palette in (("habitat", {"surface": "#2E8B6F", "cave": "#4A4A6A"}),
                      ("size", {"small": "#b9bec4", "large": "#8C5E8B"})):
    history = {lab[i]: segs for i, segs in r[name].history.items()}
    # a triangle at each switch, in the colour of the state it turned INTO
    switches = [{"kind": f"to {c.to_state}", "node": lab[c.lineage], "x": c.time}
                for c in r[name].events if c.kind != "initial"]
    (ph.trees.plot(ph.trees.loads(ct.to_newick()), skeleton=False)
     + ph.trees.color_history(history, palette=palette)
     + ph.trees.branch_events(switches, legend=False,
                              styles={f"to {s}": ("triangle_right", c) for s, c in palette.items()})
     + ph.trees.time_axis("time", bold=False)).save(f"{name}.png")
# the loop shows in the alignment: cave lineages are far likelier to be large"""


_C_SEQ_LOOP = """### two genes, each one's rate reading the other's composition — the SEQUENCE level joined to itself
from zombi2.genomes import family, simulate_genomes_family
from zombi2.params import Curve, PerSite
from zombi2.sequences import composition, gene, lg, simulate_sequences

# both families out of ONE genome run: a joint run walks their gene trees together
g = simulate_genomes_family(ct, initial_families=4, duplication=0.0, loss=0.0, seed=2,
                            families=[family("hisA"), family("hisF")])

# `start=` founds each gene AWAY from its model's equilibrium, so it ameliorates toward
# it — which is the one way a composition has anywhere to go. `offers=` is what a gene
# publishes; `absent=` answers for a lineage carrying none of it.
#
# Each rate reads the OTHER gene by name, with a step=: a composition moves with every
# substitution, so there is no interval where either rate holds still, and the run slices.
r = simulate_sequences(g, joint=True, seed=3, genes=[
    gene(name=name, model=lg(), length=250, start=kr_poor,
         offers=composition("KR", absent=0.02),
         substitution=PerSite(0.5).scaled_by(f"sequences:{other}",
                                             Curve(lambda x: 0.10 + 25.0 * x), step=0.05))
    for name, other in (("hisA", "hisF"), ("hisF", "hisA"))])

r.alignments        # one entry per family, as any sequence run returns
r.phylograms        # branch lengths accumulated slice by slice, not one sample of a rate"""


_C_TRAIT_SEQ = """### a TRAIT and a GENE's SEQUENCE, each driving the other — the two ends furthest apart
from zombi2 import joint, traits
from zombi2.params import Curve, PerLineage, PerSite
from zombi2.sequences import composition, gene, lg

# the genome run is handed over, and carries the tree its gene trees sit on:
# a sequence lives on a gene tree, which lives on the species tree
r = joint.simulate(
    traits.discrete(name="habitat", states=["cold", "hot"], start="cold",
                    switch={"cold->hot": PerLineage(0.5).scaled_by(
                                "sequences:rpoB", Curve(lambda x: 0.05 + 30.0 * x), step=0.05),
                            "hot->cold": 0.3}),
    gene(name="rpoB", model=lg(), length=250, start=kr_poor,
         offers=composition("KR", absent=0.02),
         substitution=PerSite(0.6).scaled_by("trait", {"hot": 4.0, "cold": 1.0})),
    genomes=g, seed=5)

r.trait          # the habitat, an ordinary TraitsResult
r.sequences      # the gene, an ordinary SequencesResult holding the one family"""


_C_GENOME_SEQ = """### a GENOME and a gene's SEQUENCE, each driving the other — the last cell of the map
from zombi2 import joint
from zombi2.genomes import family, genome
from zombi2.params import Curve, PerCopy
from zombi2.sequences import composition, gene, hky85

# BOTH levels are participants, so both come out. The tree is the one thing handed
# over — the gene trees are not, because the genome participant produces them.
r = joint.simulate(
    genome(duplication=0.15, origination=0.05, initial_families=25,
           families=[family("hisA")],
           loss=PerCopy(0.15).scaled_by("sequences:hisA",
                                        Curve(lambda gc: 30.0 ** ((0.35 - gc) / 0.2)),
                                        step=0.05)),
    gene(name="hisA", model=hky85(2.0), length=250, start=at_rich, substitution=0.8,
         offers=composition("GC", absent=0.35)),
    tree=ct, seed=3)

r.genome        # the gene content, an ordinary FamilyGenomesResult
r.sequences     # hisA's alignment, along the gene tree this run produced"""


JOINING = [
    Example("genome_and_sequence", "A genome and a gene's sequence",
            "The last cell of the map. <code>hisA</code> arrives AT-rich and ameliorates; while it "
            "is AT-rich its lineage loses genes fast. So how far one gene has got decides how much "
            "of the whole genome survives — and the genome decides which copies of it exist.",
            "gene content \u2194 sequence", genome_and_sequence, code=_C_GENOME_SEQ),
    Example("trait_and_sequence", "A trait and a gene's sequence",
            "The two ends of the map, joined. <code>rpoB</code> arrives KR-poor and ameliorates; "
            "the richer it gets the readier its lineage turns hot, and a hot lineage substitutes "
            "four times faster — which is what carries it further.",
            "trait \u2194 sequence", trait_and_sequence, code=_C_TRAIT_SEQ),
    Example("sequence_loop", "Two genes, each other's driver",
            "Two genes arrive KR-poor and ameliorate together: each one's substitution rate reads "
            "how far the other has got, so the pair is slow while both are behind and accelerates "
            "as neither is. The <b>sequence</b> level joined to itself.",
            "sequence \u2194 sequence", sequence_loop, code=_C_SEQ_LOOP),
    Example("trait_loop", "Two traits, each other's driver",
            "Body size decides how readily a lineage goes underground, and underground decides how "
            "readily it grows. One run, and exact: the pair is a Markov chain over their states.",
            "trait ↔ trait", trait_loop, code=_C_TRAIT_LOOP),
    Example("cave_genomes", "A trait and a genome, each other's driver",
            "The cave costs genes, and losing the eye commits a lineage to the cave. Neither can be "
            "simulated first, and the tree is an <b>input</b> here rather than an output.",
            "trait ↔ gene content", cave_genomes, code=_C_CAVE),
    Example("mobile_element_joint", "A gene drives its own genome",
            "An insertion sequence makes the genome it sits in donate genes 30x as often — and moves "
            "itself, so it spreads into lineages that never inherited it. One run, because neither "
            "half can be finished first.",
            "gene content → transfer", mobile_element_joint, code=_C_MOBILE_JOINT),
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
    Example("classe", "Change at the split",
            "The state drives how fast a lineage splits, and the split is where the state changes — "
            "the squares. Along the branches almost nothing happens "
            "(<code>switch=0.08</code>, <code>at_speciation=0.15</code>).",
            "at_speciation · joint", classe, code=_C_CLASSE),
    Example("quasse", "QuaSSE",
            "A <b>diffusing</b> body size drives speciation. Brownian motion has no direction, but "
            "the big lineages leave more descendants, so the tips end up well above the root. The "
            "one model here that slices rather than racing exactly.",
            "continuous trait → speciation", quasse, code=_C_QUASSE),
]

EXAMPLES = CONDITIONING + JOINING        # the module's full list; build.py takes the two separately
