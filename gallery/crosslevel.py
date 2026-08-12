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

from zombi2.genomes import simulate_genomes_family, simulate_genomes_ordered
from zombi2.params import Curve, PerCopy, PerLineage, PerSite
from zombi2.sequences import jc69, simulate_sequences
from zombi2.species import simulate_species_tree
from zombi2.traits import simulate_continuous, simulate_discrete

_REPAIR = {"absent": "#C2453C", "present": "#2E8B6F"}     # lose the repair gene → evolve faster
_IS1 = {"absent": "#b9bec4", "present": "#6b5b95"}        # a mobile element, present or not
_CLIMATE = {"cold": "#3A7CA5", "hot": "#E4572E"}
_OPERON = "viridis"                                       # a module's completion, 0 to 1


def _tree(n=45, seed=4):
    return simulate_species_tree(birth=1.0, death=0.2, n_extant=n, seed=seed).complete_tree


def _carrier_run(ct, name, *, seed=5):
    """The driver run: one named family, gained and lost down the tree, so presence is patchy.

    Loss a little above duplication is what makes the figure worth drawing — the family survives in
    some clades and decays in others, instead of blanketing the tree or dying on the stem."""
    return simulate_genomes_family(ct, initial_families=20, family_names=[name],
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

def repair_gene(out):
    """A GENE FAMILY drives the SUBSTITUTION rate. A mismatch-repair family is grown first; the
    sequences then evolve eight times faster on the branches that have lost it.

    The bars are root-to-tip substitutions per site. Every tip sits the same amount of *time* from
    the root, so the spread in the bars is the driver's doing and nothing else — which is what makes
    this readable at all."""
    ct = _tree()
    repair = _carrier_run(ct, "mutS")
    genes = _genes_to_evolve(ct)
    seqs = simulate_sequences(
        genes, model=jc69(), length=60, seed=1,
        substitution=PerSite(0.15).scaled_by(repair.presence("mutS"),
                                             {"present": 1.0, "absent": 8.0}))
    dist = h.root_to_tip(seqs)
    carriers = {n for n in genes.genomes if repair.has_family(int(n[1:]), "mutS")}
    tipcol = {n: _REPAIR["present" if n in carriers else "absent"] for n in dist}
    h.conditioned_figure(
        out, ct, [ph.trees.color_history(_presence_history(ct, repair.presence("mutS")),
                                         palette=_REPAIR)],
        dist, tipcol, dict(driver="mutS", states=["absent", "present"],
                           switch={"present->absent": 0.13},
                           mapping={"absent": 8, "present": 1},
                           target="substitution", target_base=0.15,
                           target_scope="PerSite", state_colors=_REPAIR),
        label="root-to-tip substitutions/site")


def climate_substitution(out):
    """A TRAIT drives the SUBSTITUTION rate — the same target as the example above, reached from the
    other end of the map. A climate trait switches along the tree and the sequences evolve eight
    times faster in the hot state.

    A driver is read wherever it changes, not once per branch, so a lineage that switches halfway
    down accumulates substitutions at one rate before the switch and another after it. The bars are
    the sum of that."""
    ct = _tree()
    climate = simulate_discrete(ct, states=["cold", "hot"], switch=0.35, seed=6)
    genes = _genes_to_evolve(ct)
    seqs = simulate_sequences(genes, model=jc69(), length=60, seed=1,
                              substitution=PerSite(0.15).scaled_by(climate,
                                                                   {"hot": 8.0, "cold": 1.0}))
    dist = h.root_to_tip(seqs)
    tipcol = {n: _CLIMATE[climate.values[n]] for n in dist if n in climate.values}
    h.conditioned_figure(
        out, ct, [ph.trees.color_history(_state_history(ct, climate), palette=_CLIMATE)],
        dist, tipcol, dict(driver="climate", states=["cold", "hot"],
                           switch={"cold->hot": 0.35, "hot->cold": 0.35},
                           mapping={"hot": 8, "cold": 1},
                           target="substitution", target_base=0.15,
                           target_scope="PerSite", state_colors=_CLIMATE),
        label="root-to-tip substitutions/site")


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
        given, tipcol, dict(driver="IS1", states=["absent", "present"],
                            switch={"present->absent": 0.13},
                            mapping={"present": 25, "absent": 1},
                            target="transfer", target_base=0.02,
                            target_sub="how often a lineage donates", state_colors=_IS1),
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
        inversions, tipcol, dict(driver="climate", states=["cold", "hot"],
                                 switch={"cold->hot": 0.35, "hot->cold": 0.35},
                                 mapping={"hot": 15, "cold": 1},
                                 target="inversion", target_base=0.05,
                                 target_sub="rate and extent", state_colors=_CLIMATE),
        label="inversions on the tip branch")


def _gc_by_node(seqs) -> dict:
    """``{node name: GC fraction}`` at every node of the tree, extant and not.

    The sequence level names each sequence ``n<species>_g<copy>``, so the species falls out of the
    label and a node's GC is taken over every gene it carries — which is what the ``gc()`` driver
    reads as the run walks the tree."""
    hits, total = collections.Counter(), collections.Counter()
    for table in (seqs.alignments, seqs.ancestral):
        for family in table.values():
            for label, seq in family.items():
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
    return simulate_genomes_ordered(ct, initial_families=24, family_names=names,
                                    modules={"repair": tuple(names)},
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
        dist, tipcol, dict(draw=h.draw_conditioning_curve, driver="repair operon", curve=factor,
                           vrange=(0.0, 1.0), value_label="fraction of the operon kept",
                           cmap=_OPERON, target="substitution", target_base=0.15,
                           target_scope="PerSite"),
        label="root-to-tip substitutions/site")


def operon_trait(out):
    """An ORDERED genome drives a TRAIT. The same repair operon, read the same way, now sets how fast
    a continuous trait diffuses instead of how fast sequences evolve — the driver is unchanged and
    only the target moved.

    Read against `operon_substitution` beside it, this is the point of the map: one driver reaches
    every level that sits on the same tree."""
    ct = _tree()
    operon = _operon_run(ct, _REPAIR_OPERON)
    factor = (lambda kept: 8.0 - 7.0 * kept)
    drift = simulate_continuous(ct, start=0.0, seed=4,
                                rate=PerLineage(0.25).scaled_by(operon.completion("repair"),
                                                                Curve(factor)))
    lab = ct.labels()
    frac = {lab[i]: sum(operon.has_family(i, g) for g in _REPAIR_OPERON) / len(_REPAIR_OPERON)
            for i in ct.nodes}
    spread = {lab[i]: abs(drift.node_values[i]) for i in ct.extant_leaves()}
    cmap, norm = plt.get_cmap(_OPERON), colors.Normalize(0.0, 1.0)
    tipcol = {n: colors.to_hex(cmap(norm(frac[n]))) for n in spread if n in frac}
    h.conditioned_figure(
        out, ct, [ph.trees.color_branches(frac, cmap=_OPERON)],
        spread, tipcol, dict(draw=h.draw_conditioning_curve, driver="repair operon", curve=factor,
                             vrange=(0.0, 1.0), value_label="fraction of the operon kept",
                             cmap=_OPERON, target="rate", target_base=0.25,
                             target_scope="PerLineage"),
        label="distance travelled from the start")


# --- a sequence drives, at either end -----------------------------------------------------------

_AT_RICH = (0.40, 0.10, 0.10, 0.40)     # A+T = 0.80 at equilibrium, so GC settles near 0.20
_GC = "cividis"


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
    """A SEQUENCE drives a SEQUENCE — compensatory evolution, and the second diagonal of the map. One
    set of genes is grown first; its **base composition** then sets how fast a second set evolves, an
    AT-rich lineage running four times faster than a GC-rich one.

    A sequence cannot drive the gene it grows inside — that would condition a run on its own output —
    but it can drive a different one, run after it. That is the whole of the restriction."""
    ct = _composition_tree()
    _, first = _composition_run(ct)
    factor = (lambda gc: 4.0 ** ((0.5 - gc) / 0.3))       # AT-rich → fast, GC-rich → slow
    second = simulate_genomes_family(ct, initial_families=4, duplication=0.02, loss=0.02, seed=8)
    seqs = simulate_sequences(second, model=jc69(), length=60, seed=5,
                              substitution=PerSite(0.2).scaled_by(first.gc(), Curve(factor)))
    dist = h.root_to_tip(seqs)
    gc = _gc_by_node(first)
    norm = colors.Normalize(min(gc.values()), max(gc.values()))
    tipcol = {n: colors.to_hex(plt.get_cmap(_GC)(norm(gc[n]))) for n in dist if n in gc}
    h.conditioned_figure(
        out, ct, [ph.trees.color_branches(gc, cmap=_GC)],
        dist, tipcol, dict(draw=h.draw_conditioning_curve, driver="GC",
                           curve=factor, vrange=(min(gc.values()), max(gc.values())),
                           value_label="GC content", cmap=_GC,
                           target="substitution", target_base=0.2, target_scope="PerSite"),
        label="root-to-tip substitutions/site")


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
        switches, tipcol, dict(draw=h.draw_conditioning_curve, driver="GC",
                               curve=factor, vrange=(min(gc.values()), max(gc.values())),
                               value_label="GC content", cmap=_GC,
                               target="switch", target_base=0.35, target_scope="PerLineage"),
        label="switches on the tip branch")

# --- the snippets shown on the detail views ----------------------------------------------------

_C_REPAIR = """### the driver: a repair family, grown first and held fixed
from zombi2.genomes import simulate_genomes_family
from zombi2.params import PerSite
from zombi2.sequences import jc69, simulate_sequences
from zombi2.species import simulate_species_tree

ct = simulate_species_tree(birth=1.0, death=0.2, n_extant=45, seed=4).complete_tree
repair = simulate_genomes_family(ct, initial_families=20, family_names=["mutS"],
                                 duplication=0.1, loss=0.13, seed=5)

### the target: sequences evolve 8x faster where the family has been lost
genes = simulate_genomes_family(ct, initial_families=4, duplication=0.02, loss=0.02, seed=2)
seqs = simulate_sequences(
    genes, model=jc69(), length=60, seed=1,
    substitution=PerSite(0.15).scaled_by(repair.presence("mutS"),
                                         {"present": 1.0, "absent": 8.0}))"""

_C_CLIMATE_SUB = """### a trait reaches the sequence level, the same target from the other end
from zombi2.traits import simulate_discrete

climate = simulate_discrete(ct, states=["cold", "hot"], switch=0.35, seed=6)
seqs = simulate_sequences(genes, model=jc69(), length=60, seed=1,
                          substitution=PerSite(0.15).scaled_by(climate,
                                                               {"hot": 8.0, "cold": 1.0}))"""

_C_MOBILE = """### one level conditioning itself: an element makes its genome donate more
from zombi2.params import PerCopy

element = simulate_genomes_family(ct, initial_families=20, family_names=["IS1"],
                                  duplication=0.1, loss=0.13, seed=5)
genome = simulate_genomes_family(
    ct, initial_families=25, loss=0.15, seed=7,
    transfer=PerCopy(0.02).scaled_by(element.presence("IS1"),
                                     {"present": 25.0, "absent": 1.0}))
# a driven transfer says how often a lineage DONATES, so count the donor side
donated = [e for e in genome.edges if e.kind == "transfer" and e.recipient is not None]"""

_C_INVERSION = """### a trait drives an ordered genome — the rate AND the extent
from zombi2.genomes import simulate_genomes_ordered
from zombi2.params import Extent, PerCopy

genome = simulate_genomes_ordered(
    ct, initial_families=30, duplication=0.05, loss=0.05, seed=3,
    inversion=PerCopy(0.05).scaled_by(climate, {"hot": 15.0, "cold": 1.0}),
    inversion_extent=Extent(3).scaled_by(climate, {"hot": 3.0, "cold": 1.0}))
# rate and extent are separate targets and multiply: hot lineages rearrange
# more often AND more coarsely"""

_C_OPERON_SUB = """### the driver is a FRACTION: how much of a six-gene operon survives
from zombi2.genomes import simulate_genomes_ordered
from zombi2.params import Curve, PerSite

operon = ["mutS", "mutL", "mutH", "uvrA", "uvrB", "uvrC"]
# at the ordered resolution a loss takes a RUN of consecutive genes, so an
# operon goes in blocks and completion drops in steps
repair = simulate_genomes_ordered(ct, initial_families=24, family_names=operon,
                                  modules={"repair": tuple(operon)},
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
            "A gene family drives the <b>substitution rate</b>. Branches that have lost the "
            "mismatch-repair family evolve eight times faster, and the bars are what came out.",
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
            "A sequence drives a sequence. It cannot drive the gene it grows inside — that "
            "would read a run’s own output — but it can drive a different one, run after it.",
            "GC → substitution", gc_drives_sequence, code=_C_GC_SEQ),
    Example("gc_drives_trait", "Composition drives a character",
            "The only pair whose driver is grown at the sequence level and whose target is not. "
            "<code>gc()</code> is a number, so the connection carries a <code>Curve</code>.",
            "GC → switch", gc_drives_trait, code=_C_GC_TRAIT),
]
