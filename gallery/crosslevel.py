"""Conditioning across the pairs the other modules do not reach.

``joining.py`` shows a trait driving a genome and a gene driving a trait. The map in Chapter 9 has
eleven pairs, and the rest of them are here: the ones that reach the **sequence** level, at either
end, and the ones whose driver is an **ordered** genome, where a gene has neighbours and an event has
an extent.

Every figure is the same shape as the conditioning ones beside it — the driver painted on the tree,
what it did in the bars — so the pairs can be read against each other.
"""

from __future__ import annotations

import collections

import helpers as h
from helpers import Example

import phylustrator as ph
import matplotlib.pyplot as plt
from matplotlib import colors

from zombi2.genomes import family, simulate_genomes_family, simulate_genomes_ordered
from zombi2.params import Curve, PerCopy, PerLineage, PerSite
from zombi2.sequences import jc69, simulate_sequences
from zombi2.species import simulate_species_tree
from zombi2.traits import simulate_continuous, simulate_discrete

_REPAIR = {"absent": "#C2453C", "present": "#2E8B6F"}     # lose the repair gene → evolve faster
_IS1 = {"absent": "#b9bec4", "present": "#6b5b95"}        # a mobile element, present or not
_CLIMATE = {"cold": "#3A7CA5", "hot": "#E4572E"}
_OPERON = "viridis_dark"                                       # a module's completion, 0 to 1


def _tree(n=45, seed=4):
    return simulate_species_tree(birth=1.0, death=0.2, n_extant=n, seed=seed).complete_tree


def _carrier_run(ct, name, *, seed=5):
    """The driver run: one named family, gained and lost down the tree, so presence is patchy.

    Loss a little above duplication is what makes the figure worth drawing — the family survives in
    some clades and decays in others, instead of blanketing the tree or dying on the stem."""
    return simulate_genomes_family(ct, initial_families=20, families=[family(name)],
                                   duplication=0.1, loss=0.13, seed=seed)


def _genes_to_evolve(ct, seed=2):
    """A small, quiet genome run for the sequence level to evolve down. Its own D/L are low: this run
    supplies gene trees, it is not the thing under study."""
    return simulate_genomes_family(ct, initial_families=4, duplication=0.02, loss=0.02, seed=seed)


def _presence_history(ct, driver):
    return {ct.labels()[i]: segs for i, segs in driver.history(ct).items()}


def _state_history(ct, trait):
    lab = ct.labels()
    return {lab[i]: segs for i, segs in trait.history.items()}


# --- a gene family drives the sequence level ---------------------------------------------------

def _in_substitutions(hist, mapping):
    """A per-branch history retimed from **time** into **substitutions**.

    A segment's share of a phylogram branch is its duration times the factor its state carries, since
    that is what the rate was multiplied by while the lineage sat in it. Colouring the phylogram with
    the time history instead would put every switch in the wrong place: the fast half of a branch
    holds most of its substitutions, so it should take up most of the branch."""
    return {name: [(state, duration * mapping[state]) for state, duration in segments]
            for name, segments in hist.items()}


def _clock_panels(out, ct, seqs, hist, mapping, palette, diagram, *, labels):
    """The two-panel sequence figure: the same tree in time, then in substitutions.

    A driven ``substitution`` rate shows in the tree's **shape**, so that is what the figure shows —
    the dated tree above, the phylogram below, the driver's colour on both. Bars of root-to-tip
    distance were the first attempt and read as two flat blocks: the factor is one of two numbers and
    root-to-tip integrates over a long shared path, so every tip landed in one of two piles.

    Both panels carry the exact history, switches mid-branch included; the lower one is retimed by
    `_in_substitutions` so each segment takes the share of the branch it actually earned."""
    style = ph.Style(width=1400, height=660, margin=96, branch_width=3.0)
    pngs = []
    for k, (newick, segments, axis) in enumerate((
            (ct.to_newick(), hist, "time"),
            (seqs.species_phylogram["complete"], _in_substitutions(hist, mapping),
             "substitutions/site"))):
        png = out.replace(".png", f"_p{k}.png")
        (ph.trees.plot(ph.trees.loads(newick), skeleton=False, style=style)
         + ph.trees.color_history(segments, palette=palette)
         + ph.trees.time_axis(axis, tick_size=20, label_size=26, bold=False)).save(png)
        pngs.append(png)
    diag = h.conditioning_png(out.replace(".png", "_diag.png"), **diagram)
    h.composite_under_diagram(out, diag, [(pngs[0], labels[0], palette), (pngs[1], labels[1])])


def repair_gene(out):
    """A GENE FAMILY drives the SUBSTITUTION rate. A mismatch-repair family is grown first; the
    sequences then evolve four times faster on the branches that have lost it.

    The same tree twice — in time, then in substitutions per site. Every tip is the same age, so the
    second panel's shape is the driver's doing and nothing else: the clades without the family run
    out to several times the divergence of the ones that kept it."""
    ct = _tree()
    repair = _carrier_run(ct, "mutS")
    genes = _genes_to_evolve(ct)
    seqs = simulate_sequences(
        genes, model=jc69(), length=60, seed=1,
        substitution=PerSite(0.15).scaled_by(repair.presence("mutS"),
                                             {"present": 1.0, "absent": 4.0}))
    _clock_panels(out, ct, seqs, _presence_history(ct, repair.presence("mutS")),
                  {"present": 1.0, "absent": 4.0}, _REPAIR,
                  dict(driver=("genomes", "mutS", "present or absent"),
                       connection=("scaled_by", "table"),
                       target_level="sequences",
                       targets=[("substitution", "rate · per site", "absent × 4    present × 1")],
                       chain=(("present", "absent"), [("0.13", None)],
                              (_REPAIR["present"], _REPAIR["absent"]))),
                  labels=("the tree in time", "the same tree in substitutions"))


def climate_substitution(out):
    """A TRAIT drives the SUBSTITUTION rate — the same target as the example above, reached from the
    other end of the map. The sequences evolve four times faster in the hot state.

    Where the gene family above is mostly held for a whole branch, a trait switches along one: a
    fifth of the branches here change state partway down. That is what the second panel is for. A
    driver is read wherever it changes, so a hot stretch takes four times the share of a branch in
    substitutions that it takes in time, and the lower tree shows exactly that."""
    ct = _tree()
    climate = simulate_discrete(ct, states=["cold", "hot"], switch=0.35, seed=6)
    genes = _genes_to_evolve(ct)
    seqs = simulate_sequences(genes, model=jc69(), length=60, seed=1,
                              substitution=PerSite(0.15).scaled_by(climate,
                                                                   {"hot": 4.0, "cold": 1.0}))
    _clock_panels(out, ct, seqs, _state_history(ct, climate), {"hot": 4.0, "cold": 1.0}, _CLIMATE,
                  dict(driver=("traits", "climate", "two states"),
                       connection=("scaled_by", "table"),
                       target_level="sequences",
                       targets=[("substitution", "rate · per site", "hot × 4    cold × 1")],
                       chain=(("cold", "hot"), [("0.35", "0.35")],
                              (_CLIMATE["cold"], _CLIMATE["hot"]))),
                  labels=("the tree in time", "the same tree in substitutions"))


# --- a gene family drives another gene family --------------------------------------------------

def mobile_element(out):
    """A GENE FAMILY drives ANOTHER — the diagonal of the map, one level conditioning itself. A
    mobile element is grown first; a second genome run then transfers twenty-five times more often
    on the branches carrying it.

    A driven ``transfer`` says how often a lineage **donates**, so the bars count transfers *given*,
    not received. Nothing here is cyclic — the element is finished before the second run starts —
    which is what keeps this conditioning rather than a joint model."""
    ct = _tree()
    element = _carrier_run(ct, "IS1")
    is1 = element.presence("IS1")
    genome = simulate_genomes_family(
        ct, initial_families=25, loss=0.15, seed=7,
        transfer=PerCopy(0.02).scaled_by(is1, {"present": 25.0, "absent": 1.0}))

    lab = ct.labels()
    # every transfer writes two edges, donor-side and recipient-side; take the recipient-side one so
    # each transfer is counted once, and group it by the branch that gave the copy away
    donated = collections.Counter(lab[e.donor] for e in genome.edges
                                  if e.kind == "transfer" and e.recipient is not None)
    given = {lab[i]: donated.get(lab[i], 0) for i in sorted(ct.extant_leaves())}
    carriers = {n for n in given if element.has_family(int(n[1:]), "IS1")}
    tipcol = {n: _IS1["present" if n in carriers else "absent"] for n in given}
    h.conditioned_figure(
        out, ct, [ph.trees.color_history(_presence_history(ct, is1), palette=_IS1)],
        given, tipcol, dict(driver=("genomes", "IS1", "present or absent"),
                            connection=("scaled_by", "table"),
                            target_level="genomes",
                            targets=[("transfer", "rate · per copy", "present × 25    absent × 1")],
                            chain=(("present", "absent"), [("0.13", None)],
                                   (_IS1["present"], _IS1["absent"]))),
        label="transfers donated")


# --- a trait drives an ordered genome, rate and extent together --------------------------------

def climate_inversions(out):
    """A TRAIT drives an ORDERED genome — and drives the **extent** as well as the rate. Hot lineages
    invert fifteen times as often, and each inversion takes a longer run of genes.

    Only the ordered and nucleotide resolutions have this pair: an extent is how much an event takes,
    which needs genes to have positions. Rate and extent are separate targets and multiply, so the
    hot branches rearrange both more often and more coarsely."""
    from zombi2.params import Extent
    ct = _tree()
    climate = simulate_discrete(ct, states=["cold", "hot"], switch=0.35, seed=6)
    genome = simulate_genomes_ordered(
        ct, initial_families=30, duplication=0.05, loss=0.05, seed=3,
        inversion=PerCopy(0.05).scaled_by(climate, {"hot": 15.0, "cold": 1.0}),
        inversion_extent=Extent(3).scaled_by(climate, {"hot": 3.0, "cold": 1.0}))

    lab = ct.labels()
    counted = collections.Counter()
    for r in genome.rearrangements:
        branch = getattr(r, "lineage", None)
        if branch is not None:
            counted[lab.get(branch)] += 1
    inversions = {lab[i]: counted.get(lab[i], 0) for i in sorted(ct.extant_leaves())}
    tipcol = {n: _CLIMATE[climate.values[n]] for n in inversions if n in climate.values}
    h.conditioned_figure(
        out, ct, [ph.trees.color_history(_state_history(ct, climate), palette=_CLIMATE)],
        inversions, tipcol, dict(driver=("traits", "climate", "two states"),
                                 connection=("scaled_by", "table"),
                                 target_level="genomes",
                                 targets=[("inversion", "rate · per copy", "hot × 15"),
                                          ("inversion_extent", "extent · in genes", "hot × 3")],
                                 chain=(("cold", "hot"), [("0.35", "0.35")],
                                        (_CLIMATE["cold"], _CLIMATE["hot"]))),
        label="inversions on the tip branch")


def _gc_by_node(seqs) -> dict:
    """``{node name: GC fraction}`` at every node of the tree, extant and not.

    The sequence level names each sequence ``n<species>_g<copy>``, so the species falls out of the
    label and a node's GC is taken over every gene it carries — which is what the ``gc()`` driver
    reads as the run walks the tree."""
    hits, total = collections.Counter(), collections.Counter()
    for table in (seqs.alignments, seqs.ancestral):
        for by_label in table.values():
            for label, seq in by_label.items():
                node = label.split("_", 1)[0]
                hits[node] += sum(seq.count(c) for c in "GC")
                total[node] += len(seq)
    return {n: hits[n] / total[n] for n in total if total[n]}


# --- an ordered genome drives, where a gene has neighbours ------------------------------------

def _operon_run(ct, names, *, seed=2):
    """The driver run at the ORDERED resolution: an operon dealt across the chromosome, losing genes
    in runs rather than one at a time. ``loss`` here takes a run of *consecutive* genes, so the
    operon decays in blocks — the reason this pair is a different model from the family one, and not
    just a different spelling of it."""
    return simulate_genomes_ordered(ct, initial_families=24,
                                    families=[family(n, module="repair") for n in names],
                                    duplication=0.08, loss=0.055, loss_extent=1.5, seed=seed)


_REPAIR_OPERON = ["mutS", "mutL", "mutH", "uvrA", "uvrB", "uvrC"]


def operon_substitution(out):
    """An ORDERED genome drives the SUBSTITUTION rate, and the driver is a **fraction** rather than a
    yes/no. How much of a four-gene repair operon a lineage still holds sets how fast its sequences
    evolve, through a `Curve`.

    Coordinates are what makes this its own model. At the ordered resolution a loss takes a *run* of
    consecutive genes, so an operon goes in blocks and completion drops in steps — a whole operon can
    disappear in one event because its genes are neighbours."""
    ct = _tree()
    operon = _operon_run(ct, _REPAIR_OPERON)
    factor = (lambda kept: 8.0 - 7.0 * kept)          # a full operon repairs; an empty one does not
    genes = _genes_to_evolve(ct)
    seqs = simulate_sequences(genes, model=jc69(), length=60, seed=1,
                              substitution=PerSite(0.15).scaled_by(operon.completion("repair"),
                                                                   Curve(factor)))
    dist = h.root_to_tip(seqs)
    lab = ct.labels()
    frac = {lab[i]: sum(operon.has_family(i, g) for g in _REPAIR_OPERON) / len(_REPAIR_OPERON)
            for i in ct.nodes}
    cmap, norm = plt.get_cmap(_OPERON), colors.Normalize(0.0, 1.0)
    tipcol = {n: colors.to_hex(cmap(norm(frac[n]))) for n in dist if n in frac}
    h.conditioned_figure(
        out, ct, [ph.trees.color_branches(frac, cmap=_OPERON)],
        dist, tipcol, dict(driver=("genomes", "repair operon", "a fraction, 0–1"),
                           connection=("scaled_by", "curve"),
                           target_level="sequences",
                           targets=[("substitution", "rate · per site", "")],
                           curve=(factor, "fraction of the operon kept", (0.0, 1.0))),
        label="root-to-tip substitutions/site")


def operon_trait(out):
    """An ORDERED genome drives a TRAIT. The same repair operon, read the same way, now sets how fast
    a continuous trait diffuses instead of how fast sequences evolve — the driver is unchanged and
    only the target moved.

    Read against `operon_substitution` beside it, this is the point of the map: one driver reaches
    every level that sits on the same tree. The colour steps *along* a branch at the instant each
    gene goes, because a completion is a number that changes mid-branch and the losses are what
    change it."""
    ct = _tree()
    operon = _operon_run(ct, _REPAIR_OPERON)
    factor = (lambda kept: 8.0 - 7.0 * kept)
    drift = simulate_continuous(ct, start=0.0, seed=4,
                                rate=PerLineage(0.25).scaled_by(operon.completion("repair"),
                                                                Curve(factor)))
    lab = ct.labels()
    hist = {lab[i]: segs for i, segs in operon.completion("repair").history(ct).items()}
    spread = {lab[i]: abs(drift.node_values[i]) for i in sorted(ct.extant_leaves())}
    cmap, norm = plt.get_cmap(_OPERON), colors.Normalize(0.0, 1.0)
    tipcol = {n: colors.to_hex(cmap(norm(hist[n][-1][0]))) for n in spread if n in hist}
    h.conditioned_figure(
        out, ct, [ph.trees.color_history(hist, cmap=_OPERON, limits=(0.0, 1.0))],
        spread, tipcol, dict(driver=("genomes", "repair operon", "a fraction, 0–1"),
                             connection=("scaled_by", "curve"),
                             target_level="traits",
                             targets=[("rate", "rate · per lineage", "")],
                             curve=(factor, "fraction of the operon kept", (0.0, 1.0))),
        label="distance travelled from the start")


# --- a sequence drives, at either end -----------------------------------------------------------

_AT_RICH = (0.40, 0.10, 0.10, 0.40)     # A+T = 0.80 at equilibrium, so GC settles near 0.20
_GC = "magma_dark"


def _composition_run(ct):
    """The driver run: one clade evolves under an AT-rich model and the rest under an even one, so
    base composition genuinely differs between lineages instead of drifting around one value.

    A composition driver needs something to read. Under a single model GC wanders within a few
    percent of its equilibrium and most of that is sampling, which would make the figure below a
    picture of noise."""
    from zombi2.params import Clade
    from zombi2.sequences import Models, hky85
    genes = simulate_genomes_family(ct, initial_families=6, duplication=0.0, loss=0.0, seed=2)
    clade = Clade({"at": ["n13", "n33"]})
    seqs = simulate_sequences(
        genes, length=300, substitution=0.55, seed=2,
        model=Models().set_by(clade, {"at": hky85(kappa=2.0, frequencies=_AT_RICH),
                                      "rest": hky85(kappa=2.0)}))
    return genes, seqs


def _composition_tree():
    """The tree `sequences.clade_own_model` uses, because the clade it names is defined on it."""
    return simulate_species_tree(birth=1.0, n_extant=34, seed=2).complete_tree


def gc_drives_sequence(out):
    """A SEQUENCE drives a SEQUENCE — the second diagonal of the map. One set of genes is grown
    first; its **base composition** then sets how fast a second set evolves, an AT-rich lineage
    running four times faster than a GC-rich one.

    A composition is not a mechanism — no count of bases reaches out and sets another gene's rate —
    so what it is good for is standing in for a property of the *lineage*. A bacterium's GC content
    is largely the balance between an AT-biased mutation spectrum and what counteracts it, and the
    lineages that let it slide toward AT are the ones whose mismatch repair has weakened, which is
    also what makes them substitute faster everywhere. So this run reads GC as an index of the
    mutational regime, not as a cause.

    A sequence cannot drive the gene it grows inside — that would condition a run on its own output —
    but it can drive a different one, run after it. That is the whole of the restriction. Here the
    two sets come from two genome runs; the card after this one takes both families out of **one**,
    which is what makes the driver's gene tree the one the target's run also saw.

    Coloured with `color_branches` rather than the stepped history the operon examples use: GC drifts
    continuously with every substitution, so there is no event to step at, and a gradient down the
    branch is what the quantity actually does."""
    ct = _composition_tree()
    _, first = _composition_run(ct)
    factor = (lambda gc: 4.0 ** ((0.5 - gc) / 0.3))       # AT-rich → fast, GC-rich → slow
    second = simulate_genomes_family(ct, initial_families=4, duplication=0.02, loss=0.02, seed=8)
    seqs = simulate_sequences(second, model=jc69(), length=60, seed=5,
                              substitution=PerSite(0.2).scaled_by(first.gc(), Curve(factor)))
    gc = _gc_by_node(first)
    span = (min(gc.values()), max(gc.values()))
    style = ph.Style(width=1400, height=660, margin=96, branch_width=3.0)
    pngs = []
    for k, (newick, axis) in enumerate(((ct.to_newick(), "time"),
                                        (seqs.species_phylogram["complete"], "substitutions/site"))):
        png = out.replace(".png", f"_p{k}.png")
        (ph.trees.plot(ph.trees.loads(newick), skeleton=False, style=style)
         + ph.trees.color_branches(gc, cmap=_GC, limits=span)
         + ph.trees.time_axis(axis, tick_size=20, label_size=26, bold=False)).save(png)
        pngs.append(png)
    diag = h.conditioning_png(out.replace(".png", "_diag.png"),
                              driver=("sequences", "GC", "a number"),
                              connection=("scaled_by", "curve"),
                              target_level="sequences",
                              targets=[("substitution", "rate · per site", "")],
                              curve=(factor, "GC content", span))
    h.composite_under_diagram(out, diag, [(pngs[0], "the tree in time"),
                                          (pngs[1], "the same tree in substitutions")])


_IVYWREL = "IVYWREL"          # the thermophily signature: Zeldovich, Berezovsky & Shakhnovich 2007
_MARKER_ABSENT = 0.40         # no marker → the driver reads LG's own IVYWREL share


def _ivywrel_model(scale: float, name: str):
    """LG's exchangeabilities over frequencies with the seven residues re-weighted by ``scale``.

    A thermophile's proteins are not a different chemistry — which residues swap easily for which is
    the same — they are the same chemistry at a different composition. So the exchangeabilities stay
    LG's and only ``π`` moves, which is what `~zombi2.sequences.substitution_models.reversible` is
    the door for. The equilibria this reaches, 0.32 and 0.52, straddle the real mesophile-to-
    thermophile range (about 0.38 to 0.48) and are pulled apart a little so the gradient reads on a
    tree this size."""
    import numpy as np
    from zombi2.sequences.substitution_models import AMINO_ACIDS, lg, reversible

    base = lg()
    S = base.Q / base.stationary[None, :]
    S = (S + S.T) / 2.0
    np.fill_diagonal(S, 0.0)
    pi = base.stationary.copy()
    pi[[AMINO_ACIDS.index(c) for c in _IVYWREL]] *= scale
    return reversible(S, pi / pi.sum(), name=name, alphabet=AMINO_ACIDS)


def named_family_drives_sequence(out):
    """One **named family's** composition driving another named family's rate, out of a single
    genome run — and the thing that used to make this unwritable.

    A composition is pooled over whatever a run evolved, so on a whole genome it is the lineage's
    whole complement and belongs to no family in particular. Restricting the run to one family with
    ``families=["marker"]`` makes the pooled statistic that family's. Both runs then read the
    **same** genome run, so the driver's history is the one on disk; two separate genome runs would
    have put the driver on gene trees the target never saw.

    The model. A composition is not a mechanism — nothing about a protein's residue counts reaches
    out and sets another gene's rate — so what it is good for is standing in for a property of the
    *lineage*. The share of a protein that is I, V, Y, W, R, E or L rises with the temperature its
    owner lives at, which is the one compositional signature of thermophily that holds across the
    tree of life (Zeldovich, Berezovsky & Shakhnovich 2007). So the marker's IVYWREL share is read
    here as how hot the lineage runs, and this run says hot lineages substitute four times faster.
    One clade is thermophilic, and its marker settles near 0.52 against 0.32 for the rest.

    What that costs is the branches with no marker at all — grey here, and a fifth of them. A
    driver has to answer for every branch the target walks, so the run declares what they read:
    ``absent=0.40``, the background share. Leaving it out is an error rather than a guess, because
    carrying the parent's value forward would drive those branches as though the marker were still
    there.

    The second panel is the same tree measured in substitutions, so the branches the driver sped up
    are the long ones.
    """
    from zombi2.params import Clade
    from zombi2.sequences import Models, lg

    ct = _composition_tree()
    g = simulate_genomes_family(ct, initial_families=6, duplication=0.0, loss=0.0, origination=0.0,
                                seed=2, families=[family("marker", loss=0.30), family("ribosomal")])
    marker = simulate_sequences(
        g, families=["marker"], length=300, substitution=0.55, seed=2,
        model=Models().set_by(Clade({"hot": ["n13", "n33"]}),
                              {"hot": _ivywrel_model(1.6, "LG+thermophile"),
                               "rest": _ivywrel_model(0.7, "LG+mesophile")}))
    factor = (lambda share: 4.0 ** ((share - 0.32) / 0.20))    # hot → fast, cool → slow
    driver = marker.composition(_IVYWREL, absent=_MARKER_ABSENT)
    target = simulate_sequences(g, families=["ribosomal"], model=lg(), length=200, seed=5,
                                substitution=PerSite(0.2).scaled_by(driver, Curve(factor)))

    lab = ct.labels()
    read = driver._node_values(ct)
    # only the branches that HAVE the family are painted, so the grey ones are exactly the branches
    # `absent=` had to answer for — a partial mapping, which Phylustrator leaves in its own colour
    share = {lab[i]: v for i, v in read.items() if v != _MARKER_ABSENT}
    span = (min(share.values()), max(share.values()))
    # the branch colour a lineage keeps when the map has no value for it, set light on purpose: the
    # low end of this ramp is near-black, and against the default dark it was the one thing on the
    # figure a reader could not tell from another
    style = ph.Style(width=1400, height=660, margin=96, branch_width=3.0, branch_color="#c3c8cc")
    pngs = []
    for k, (newick, axis) in enumerate(((ct.to_newick(), "time"),
                                        (target.species_phylogram["complete"],
                                         "substitutions/site"))):
        png = out.replace(".png", f"_p{k}.png")
        (ph.trees.plot(ph.trees.loads(newick), skeleton=False, style=style)
         + ph.trees.color_branches(share, cmap=_GC, limits=span)
         + ph.trees.time_axis(axis, tick_size=20, label_size=26, bold=False)).save(png)
        pngs.append(png)
    diag = h.conditioning_png(out.replace(".png", "_diag.png"),
                              driver=("sequences", "marker IVYWREL", "a number"),
                              connection=("scaled_by", "curve"),
                              target_level="sequences",
                              targets=[("ribosomal substitution", "rate · per site", "")],
                              curve=(factor, "the marker's IVYWREL share", span))
    h.composite_under_diagram(
        out, diag,
        [(pngs[0], "the marker's IVYWREL share  (grey: no marker, so the driver reads absent=0.40)",
          (_GC, "cool", "hot")),
         (pngs[1], "the ribosomal family, in substitutions")])


def gc_drives_trait(out):
    """A SEQUENCE drives a TRAIT — the last pair on the map, and the only one whose driver is grown at
    the sequence level and whose target is not. Base composition sets how fast a discrete character
    switches, so the AT-rich clade flickers between states while the rest stays put.

    ``gc()`` is a number, so the connection carries a `Curve` rather than a table of per-state
    multipliers: there are no states to enumerate on a driver that is continuous."""
    ct = _composition_tree()
    _, first = _composition_run(ct)
    factor = (lambda gc: 4.0 ** ((0.5 - gc) / 0.3))
    switching = simulate_discrete(ct, states=["A", "B"], start="A", seed=11,
                                  switch=PerLineage(0.35).scaled_by(first.gc(), Curve(factor)))
    gc = _gc_by_node(first)
    lab = ct.labels()
    flips = collections.Counter(c.lineage for c in switching.events if c.kind == "on_branch")
    switches = {lab[i]: flips.get(i, 0) for i in sorted(ct.extant_leaves())}
    norm = colors.Normalize(min(gc.values()), max(gc.values()))
    tipcol = {n: colors.to_hex(plt.get_cmap(_GC)(norm(gc[n]))) for n in switches if n in gc}
    h.conditioned_figure(
        out, ct, [ph.trees.color_branches(gc, cmap=_GC)],
        switches, tipcol, dict(driver=("sequences", "GC", "a number"),
                               connection=("scaled_by", "curve"),
                               target_level="traits",
                               targets=[("switch", "rate · per lineage", "")],
                               curve=(factor, "GC content",
                                      (min(gc.values()), max(gc.values())))),
        label="switches on the tip branch")

# --- the snippets shown on the detail views ----------------------------------------------------

_C_REPAIR = """### the driver: a repair family, grown first and held fixed
from zombi2.genomes import family, simulate_genomes_family
from zombi2.params import PerSite
from zombi2.sequences import jc69, simulate_sequences
from zombi2.species import simulate_species_tree

ct = simulate_species_tree(birth=1.0, death=0.2, n_extant=45, seed=4).complete_tree
repair = simulate_genomes_family(ct, initial_families=20, families=[family('mutS')],
                                 duplication=0.1, loss=0.13, seed=5)

### the target: sequences evolve 8x faster where the family has been lost
genes = simulate_genomes_family(ct, initial_families=4, duplication=0.02, loss=0.02, seed=2)
seqs = simulate_sequences(
    genes, model=jc69(), length=60, seed=1,
    substitution=PerSite(0.15).scaled_by(repair.presence("mutS"),
                                         {"present": 1.0, "absent": 4.0}))"""

_C_CLIMATE_SUB = """### a trait reaches the sequence level, the same target from the other end
from zombi2.traits import simulate_discrete

climate = simulate_discrete(ct, states=["cold", "hot"], switch=0.35, seed=6)
seqs = simulate_sequences(genes, model=jc69(), length=60, seed=1,
                          substitution=PerSite(0.15).scaled_by(climate,
                                                               {"hot": 8.0, "cold": 1.0}))"""

_C_MOBILE = """### one level conditioning itself: an element makes its genome donate more
from zombi2.params import PerCopy

element = simulate_genomes_family(ct, initial_families=families=[family('IS1')]S1"],
                                  duplication=0.1, loss=0.13, seed=5)
genome = simulate_genomes_family(
    ct, initial_families=25, loss=0.15, seed=7,
    transfer=PerCopy(0.02).scaled_by(element.presence("IS1"),
                                     {"present": 25.0, "absent": 1.0}))
# a driven transfer says how often a lineage DONATES, so count the donor side
donated = [e for e in genome.edges if e.kind == "transfer" and e.recipient is not None]"""

_C_INVERSION = """### a trait drives an ordered genome — the rate AND the extent
from zombi2.genomes import family, simulate_genomes_ordered
from zombi2.params import Extent, PerCopy

genome = simulate_genomes_ordered(
    ct, initial_families=30, duplication=0.05, loss=0.05, seed=3,
    inversion=PerCopy(0.05).scaled_by(climate, {"hot": 15.0, "cold": 1.0}),
    inversion_extent=Extent(3).scaled_by(climate, {"hot": 3.0, "cold": 1.0}))
# rate and extent are separate targets and multiply: hot lineages rearrange
# more often AND more coarsely"""

_C_OPERON_SUB = """### the driver is a FRACTION: how much of a six-gene operon survives
from zombi2.genomes import family, simulate_genomes_ordered
from zombi2.params import Curve, PerSite

operon = ["mutS", "mutL", "mutH", "uvrA", "uvrB", "uvrC"]
# at the ordered resolution a loss takes a RUN of consecutive genes, so an
# operon goes in blocks and completion drops in steps
repair = simulate_genomes_ordered(ct, initial_families=24,
                                  families=[family(n, module="repair") for n in operon],
                                  duplication=0.08, loss=0.055, loss_extent=1.5, seed=2)

seqs = simulate_sequences(genes, model=jc69(), length=60, seed=1,
                          substitution=PerSite(0.15).scaled_by(
                              repair.completion("repair"), Curve(lambda kept: 8.0 - 7.0 * kept)))"""

_C_OPERON_TRAIT = """### the same driver, a different level: it now sets how fast a trait diffuses
from zombi2.params import Curve, PerLineage
from zombi2.traits import simulate_continuous

drift = simulate_continuous(ct, start=0.0, seed=4,
                            rate=PerLineage(0.25).scaled_by(
                                repair.completion("repair"),
                                Curve(lambda kept: 8.0 - 7.0 * kept)))"""

_C_GC_SEQ = """### one set of genes sets how fast another evolves — compensatory evolution
from zombi2.params import Clade, Curve, PerSite
from zombi2.sequences import Models, hky85, jc69, simulate_sequences

# one clade under an AT-rich model, so composition really differs
first = simulate_sequences(
    genes, length=300, substitution=0.55, seed=2,
    model=Models().set_by(Clade({"at": ["n13", "n33"]}),
                          {"at": hky85(kappa=2.0, frequencies=(0.40, 0.10, 0.10, 0.40)),
                           "rest": hky85(kappa=2.0)}))

second = simulate_genomes_family(ct, initial_families=4, duplication=0.02, loss=0.02, seed=8)
seqs = simulate_sequences(second, model=jc69(), length=60, seed=5,
                          substitution=PerSite(0.2).scaled_by(
                              first.gc(), Curve(lambda gc: 4.0 ** ((0.5 - gc) / 0.3))))"""

_C_NAMED_SEQ = """### ONE NAMED FAMILY's composition drives another's rate — one genome run, two sequence runs
from zombi2.genomes import family, simulate_genomes_family
from zombi2.params import Clade, Curve, PerSite
from zombi2.sequences import Models, lg, simulate_sequences

# both families come out of ONE genome run, so the driver's gene tree is the one
# the target's run also saw. Per-family rates: the marker is lost sometimes.
g = simulate_genomes_family(ct, initial_families=6, duplication=0.0, loss=0.0, seed=2,
                            families=[family("marker", loss=0.30), family("ribosomal")])

# one clade is thermophilic: LG's chemistry over IVYWREL-enriched frequencies, so the
# marker's composition really differs (`hot`/`cool` are built with reversible(S, pi))
marker = simulate_sequences(
    g, families=["marker"], length=300, substitution=0.55, seed=2,
    model=Models().set_by(Clade({"hot": ["n13", "n33"]}), {"hot": hot, "rest": cool}))

# `families=` restricts the run, so its POOLED composition is that family's. The share
# of a protein that is I,V,Y,W,R,E or L rises with growth temperature (Zeldovich et al.
# 2007), so it is read here as how hot the lineage runs — a proxy, not a mechanism.
#
# absent= is what a branch WITHOUT the marker reads. Required, not guessed: a driver
# has to answer for every branch the target walks, and carrying the parent's value
# forward would say the family was still there.
ribosomal = simulate_sequences(
    g, families=["ribosomal"], model=lg(), length=200, seed=5,
    substitution=PerSite(0.2).scaled_by(marker.composition("IVYWREL", absent=0.40),
                                        Curve(lambda s: 4.0 ** ((s - 0.32) / 0.20))))"""

_C_GC_TRAIT = """### base composition sets how fast a character switches
from zombi2.params import Curve, PerLineage
from zombi2.traits import simulate_discrete

switching = simulate_discrete(ct, states=["A", "B"], start="A", seed=11,
                              switch=PerLineage(0.35).scaled_by(
                                  first.gc(), Curve(lambda gc: 4.0 ** ((0.5 - gc) / 0.3))))
# gc() is a number, so the connection carries a Curve: a continuous driver has
# no states to enumerate in a table"""


EXAMPLES = [
    Example("repair_gene", "Lose the repair gene, evolve faster",
            "A gene family drives the <b>substitution rate</b>. The same tree in time, then in "
            "substitutions: branches that have lost the mismatch-repair family run four times longer.",
            "gene → substitution", repair_gene, code=_C_REPAIR),
    Example("climate_substitution", "A trait sets the substitution rate",
            "The same target reached from the other end of the map. A driver is read wherever it "
            "changes, so a lineage that switches halfway down accumulates at two rates.",
            "trait → substitution", climate_substitution, code=_C_CLIMATE_SUB),
    Example("mobile_element", "A mobile element spreads its genome",
            "One level conditioning itself. Carrying the element makes a lineage <b>donate</b> "
            "twenty-five times more often — a driven transfer is about giving, not receiving.",
            "gene → transfer", mobile_element, code=_C_MOBILE),
    Example("climate_inversions", "A trait drives rate and extent",
            "Hot lineages invert fifteen times as often, and each inversion takes a longer run of "
            "genes. Rate and extent are separate targets and multiply.",
            "trait → inversion", climate_inversions, code=_C_INVERSION),
    Example("operon_substitution", "How much of an operon is left",
            "The driver is a <b>fraction</b>, not a yes/no. At the ordered resolution a loss takes a "
            "run of neighbours, so an operon goes in blocks and completion drops in steps.",
            "module → substitution", operon_substitution, code=_C_OPERON_SUB),
    Example("operon_trait", "The same driver, another level",
            "The repair operon again, now setting how fast a trait diffuses instead of how fast "
            "sequences evolve. One driver reaches every level on the tree.",
            "module → trait rate", operon_trait, code=_C_OPERON_TRAIT),
    Example("gc_drives_sequence", "One gene’s composition, another’s rate",
            "GC is read as an <b>index of the mutational regime</b> rather than as a cause: the "
            "lineages that slide AT-ward are the ones whose repair has weakened, and they "
            "substitute faster everywhere.",
            "GC → substitution", gc_drives_sequence, code=_C_GC_SEQ),
    Example("named_family_drives_sequence", "One named gene drives another",
            "A marker’s <b>IVYWREL</b> share stands in for how hot the lineage runs — the one "
            "compositional signature of thermophily. Out of <b>one</b> genome run, with "
            "<code>absent=</code> answering for the branches without the marker (grey).",
            "one family’s composition → substitution", named_family_drives_sequence,
            code=_C_NAMED_SEQ),
    Example("gc_drives_trait", "Composition drives a character",
            "The only pair whose driver is grown at the sequence level and whose target is not. "
            "<code>gc()</code> is a number, so the connection carries a <code>Curve</code>.",
            "GC → switch", gc_drives_trait, code=_C_GC_TRAIT),
]
