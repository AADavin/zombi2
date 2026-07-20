"""Genomes I — the unordered D/T/L/O gene-family core.

A genome is a multiset of gene families that evolves along the species tree by four events:
**origination** (a new family arises in a lineage — per lineage), **duplication** (a gene copy
duplicates, +1 copy in its family — per copy), **loss** (a gene copy is lost — per copy), and
**transfer** (a copy is donated to a *contemporaneous* lineage — per copy). Rates are the
cross-level ``scope(base) × modifiers`` grammar (``SPEC §5``); the defaults are the natural
"per what?" for each event.

This reads as the genome twin of :mod:`zombi2.species_tree`: one forward Gillespie over the
**complete** tree, plain/frozen dataclasses, an event log as the source of truth (per-family gene
trees are derived from it later), ``as_rate``/``.effective`` for every rate. Because a transfer at
time ``t`` couples two lineages alive at ``t``, the engine evolves **all lineages alive at once**
along one global clock — exactly like ``species_tree._grow`` over its ``alive`` list, except the
species tree is a fixed input (its ``end_time``s form a schedule that decides who is alive), so
there is no birth-death race, no survival conditioning. Speciations and extinctions from that
schedule enter/retire lineages; between them one Gillespie fires D/T/L/O. ``transfer=0`` is the
special case where the lineages are independent — same law as evolving each segment alone.

Still to come: per-family heterogeneity (``ByFamily`` + ``Speed``), the sparse profiles and lazy
gene-tree views behind the ``record=`` memory dial, and the Rust core. This lives here for now so
the legacy ``zombi2/genomes`` package is untouched.
"""

from __future__ import annotations

import collections
import pathlib
from dataclasses import dataclass
from functools import cached_property

import numpy as np

from ..rates.modifiers import OnTime
from ..rates.rate import as_rate
from ..rates.scope import PerCopy, PerLineage
from ..species import SpeciesResult, Tree
from ._live import enter, retire
from ._transfer import Distance, mean_root_to_tip, recipient_index
from .events import Event, events_tsv
from .gene_trees import GeneNode, GeneTree, gene_trees_from_events
from .ordered import (
    Chromosome,
    ChromosomeEvent,
    Gene,
    Inversion,
    OrderedGenomesResult,
    Translocation,
    Transposition,
    simulate_genomes_ordered,
)
from .profiles import Profiles, profiles_from_genomes


@dataclass(frozen=True)
class GeneCopy:
    """One gene copy: a member of family ``family``, identified by a globally-unique ``id``. Its
    birth/death times and parentage live in the event log (the source of truth); the copy carries
    only what a genome snapshot needs to be self-describing — who it is and which family it is in. A
    genome may hold several copies sharing a ``family`` (that family's copy count)."""

    id: int
    family: int


@dataclass
class GenomesResult:
    """What ``simulate_genomes_unordered`` returns: the ``complete_tree`` it ran on, the final
    ``genomes`` at **every** node (extant and extinct), the ``events`` log (the compact source of
    truth), and the ``seed``. The observed genomes are the extant tips —
    ``{n.id: genomes[n.id] for n in complete_tree.extant()}``. The phyletic ``profiles`` are derived
    from those tips on access, and ``write`` materialises the chosen outputs to disk. (Lazy gene
    trees and the ``record=`` scale dial are later slices.)"""

    complete_tree: Tree
    genomes: dict[int, tuple[GeneCopy, ...]]
    events: list[Event]
    seed: int | None

    def family_counts(self, node_id: int) -> collections.Counter:
        """A multiset view of one node's genome: ``family id → copy count``."""
        return collections.Counter(c.family for c in self.genomes[node_id])

    @cached_property
    def profiles(self) -> Profiles:
        """The phyletic profiles — each gene family's copy count in each extant species — derived
        from the observed genomes (the classic comparative-genomics matrix). See :mod:`.profiles`."""
        extant = [n.id for n in self.complete_tree.extant()]
        return profiles_from_genomes(self.genomes, extant)

    @cached_property
    def gene_trees(self) -> dict[int, GeneTree]:
        """``{family id: GeneTree}`` — each family's true genealogy inside the complete tree,
        derived from the event log. Each ``GeneTree`` exposes ``.complete`` and ``.extant``. See
        :mod:`.gene_trees`."""
        return gene_trees_from_events(self.events, self.complete_tree)

    def write(self, directory, outputs=("events", "profiles")) -> None:
        """Materialise chosen ``outputs`` to ``directory`` (created if needed):

        - ``"events"`` → ``genome_events.tsv``, the event log (the source of truth).
        - ``"profiles"`` → ``profiles.tsv``, the family × extant-species copy-count matrix.
        """
        d = pathlib.Path(directory)
        d.mkdir(parents=True, exist_ok=True)
        if "events" in outputs:
            (d / "genome_events.tsv").write_text(events_tsv(self.events))
        if "profiles" in outputs:
            (d / "profiles.tsv").write_text(self.profiles.to_tsv())


# --- the live genomes: parallel arrays under swap-remove, the ``species_tree._grow`` shape --------

def _pick_copy(rng, gen, total_copies) -> tuple[int, int]:
    """A uniform global copy pick → ``(lineage index k, copy index j in gen[k])``. Realises
    per-copy scope across the whole pool: every copy, in any lineage, is equally likely."""
    j = int(rng.integers(total_copies))
    for k, g in enumerate(gen):
        if j < len(g):
            return k, j
        j -= len(g)
    raise AssertionError("total_copies out of sync with the genomes")  # unreachable


# --- the D/T/L/O mutators (each records to the event log; ids from the minters) -------------------

def _originate(genome, node, t, events, new_copy, new_family) -> None:
    """A new gene family arises: mint a founding copy in a fresh family and record it."""
    c = new_copy(new_family())
    genome.append(c)
    events.append(Event(t, "origination", node.id, c.family, c.id))


def _duplicate(genome, j, node, t, events, new_copy) -> None:
    """The gene at index ``j`` duplicates. In the ZOMBI1 per-segment model every event re-ids: the
    gene *ends* and **two** fresh copies descend from it, so both carry new ids (and the id in any
    node is that node's own)."""
    old = genome[j]
    cont, dup = new_copy(old.family), new_copy(old.family)
    genome[j] = cont                                   # the continuing lineage (a fresh id)
    genome.append(dup)                                 # the new copy (a fresh id)
    events.append(Event(t, "duplication", node.id, old.family, cont.id, parent=old.id))
    events.append(Event(t, "duplication", node.id, old.family, dup.id, parent=old.id))


def _lose_at(genome, j, node, t, events) -> None:
    """The copy at index ``j`` is lost (swap-remove — the genome is an order-agnostic multiset)."""
    lost = genome[j]
    genome[j] = genome[-1]
    genome.pop()
    events.append(Event(t, "loss", node.id, lost.family, lost.id))


def _do_transfer(rng, tree, alive, gen, total_copies, t, events, new_copy,
                 transfer_to, replacement, self_transfer, depth) -> int:
    """A gene transfers from a donor copy to a contemporaneous recipient lineage. Returns the change
    in total copy count: +1 additive, 0 replacement (the arriving copy displaces a resident)."""
    kd, jd = _pick_copy(rng, gen, total_copies)
    donor = alive[kd]
    src = gen[kd][jd]
    fam = src.family
    cand = [k for k in range(len(alive)) if self_transfer or k != kd]
    kr = recipient_index(rng, tree, alive, cand, donor, t, transfer_to, depth)
    recipient = alive[kr]
    rg = gen[kr]
    # the donor gene ends; two fresh copies descend from it (ZOMBI1 re-id): the continuation on the
    # donor branch and the transferred copy on the recipient branch — a horizontal edge in the gene tree.
    cont, xfer = new_copy(fam), new_copy(fam)
    gen[kd][jd] = cont
    delta = 1
    if replacement:
        residents = [p for p, c in enumerate(rg) if c.family == fam and c.id != cont.id]
        if residents:  # homologous overwrite; empty ⇒ additive fallback (the gene still arrives)
            p = residents[int(rng.integers(len(residents)))]
            victim = rg[p]
            rg[p] = rg[-1]
            rg.pop()
            events.append(Event(t, "loss", recipient, fam, victim.id))
            delta = 0
    rg.append(xfer)
    events.append(Event(t, "transfer", donor, fam, cont.id, parent=src.id))
    events.append(Event(t, "transfer", recipient, fam, xfer.id, parent=src.id, recipient=recipient))
    return delta


def simulate_genomes_unordered(tree, *, duplication=0.0, transfer=0.0, loss=0.0, origination=0.0,
                               transfer_to="uniform", replacement=False, self_transfer=False,
                               initial_families=0, seed=None) -> GenomesResult:
    """Evolve a multiset of gene families along a species tree by duplication, transfer, loss, and
    origination.

    ``tree`` is the **complete** species tree (a :class:`~zombi2.species_tree.Tree`, or a
    :class:`~zombi2.species_tree.SpeciesResult` whose ``complete_tree`` is used). Genomes evolve on
    **every** lineage, extant and extinct alike, so the true gene-tree history is complete and a
    transfer can arrive "from the dead"; the observed genomes are the extant tips.

    Rates (each a ``scope(base) × modifiers`` spec): ``duplication``/``transfer``/``loss`` default
    **per copy**, ``origination`` **per lineage**. When a transfer fires it moves a copy from a
    uniformly-chosen donor copy to a recipient lineage alive at that instant, chosen by
    ``transfer_to`` — ``"uniform"`` (any other contemporaneous lineage) or ``"distance"`` /
    ``Distance(decay=)`` (closer relatives likelier). ``replacement=True`` overwrites a homologous
    copy in the recipient (additive fallback if it has none); ``self_transfer=True`` lets a lineage
    donate to itself. The root starts with ``initial_families`` families of one copy each, recorded
    as originations at the crown. Deterministic given ``seed``.
    """
    tree = tree.complete_tree if isinstance(tree, SpeciesResult) else tree
    dup = as_rate(duplication, default_scope=PerCopy)
    tra = as_rate(transfer, default_scope=PerCopy)
    los = as_rate(loss, default_scope=PerCopy)
    org = as_rate(origination, default_scope=PerLineage)
    # this slice wires only the default scope (D/T/L per copy, origination per lineage) and OnTime
    # (skyline) modifiers. A non-default scope would set the *total* rate one way while the engine
    # still picks the affected copy/lineage the default way — a silent mismatch (e.g. a PerCopy
    # origination is base×0 copies, a no-op) — so reject it, as we reject non-OnTime modifiers.
    for label, rate, want in (("duplication", dup, PerCopy), ("transfer", tra, PerCopy),
                              ("loss", los, PerCopy), ("origination", org, PerLineage)):
        if not isinstance(rate.scope, want):
            raise ValueError(
                f"{label} has a {type(rate.scope).__name__} scope, but the unordered genome engine "
                f"wires only {want.__name__} for {label} this slice — scope overrides are a later slice."
            )
        for m in rate.modifiers:
            if not isinstance(m, OnTime):
                raise ValueError(
                    f"{label} carries {type(m).__name__}, which the unordered genome engine does not "
                    f"support yet — only OnTime (skyline) is wired. Per-family heterogeneity (ByFamily, "
                    f"Speed) and clade drift are later slices."
                )
    if transfer_to == "distance":
        transfer_to = Distance()
    if transfer_to != "uniform" and not isinstance(transfer_to, Distance):
        raise ValueError(f"transfer_to must be 'uniform', 'distance', or Distance(decay=), got {transfer_to!r}")
    if isinstance(initial_families, bool) or not isinstance(initial_families, int) or initial_families < 0:
        raise ValueError(f"initial_families must be a non-negative integer, got {initial_families!r}")

    rng = np.random.default_rng(seed)
    copy_counter = 0
    family_counter = 0

    def new_copy(family: int) -> GeneCopy:
        nonlocal copy_counter
        c = GeneCopy(copy_counter, family)
        copy_counter += 1
        return c

    def new_family() -> int:
        nonlocal family_counter
        f = family_counter
        family_counter += 1
        return f

    depth = mean_root_to_tip(tree)  # timescale for Distance weighting (unused by "uniform")
    schedule = sorted((tree.nodes[i].end_time, i) for i in tree.nodes)  # (end_time, node_id)

    root = tree.nodes[tree.root]
    t = root.birth_time
    alive: list[int] = []
    gen: list[list[GeneCopy]] = []
    pos: dict[int, int] = {}
    genomes: dict[int, tuple[GeneCopy, ...]] = {}
    events: list[Event] = []
    enter(alive, gen, pos, root.id, [])
    for _ in range(initial_families):  # seed the crown as originations at t = root.birth_time
        _originate(gen[0], root, t, events, new_copy, new_family)
    total_copies = len(gen[0])

    si = 0
    while si < len(schedule):
        n = total_copies
        k_alive = len(alive)
        ctx = {"copies": n, "lineages": k_alive, "time": t}
        r_dup = dup.effective(**ctx) if n else 0.0      # per copy → base × N × f(t), pooled over all lineages
        r_los = los.effective(**ctx) if n else 0.0      # per copy
        r_org = org.effective(**ctx)                     # per lineage → base × K × f(t)
        can_xfer = n > 0 and (k_alive >= 2 or self_transfer)  # a recipient must be able to exist
        r_tra = tra.effective(**ctx) if can_xfer else 0.0    # per copy
        total = r_dup + r_los + r_org + r_tra

        next_species = schedule[si][0]  # the tree's own next event: who is alive changes only here
        horizon = min(next_species, dup.next_change(t), los.next_change(t),
                      org.next_change(t), tra.next_change(t))

        if total > 0.0:
            t_ev = t + float(rng.exponential(1.0 / total))
            if t_ev < horizon:  # a genome event fires before the alive set or the rate changes
                t = t_ev
                r = float(rng.random()) * total
                if r < r_dup:
                    k, j = _pick_copy(rng, gen, n)
                    _duplicate(gen[k], j, tree.nodes[alive[k]], t, events, new_copy)
                    total_copies += 1
                elif r < r_dup + r_los:
                    k, j = _pick_copy(rng, gen, n)
                    _lose_at(gen[k], j, tree.nodes[alive[k]], t, events)
                    total_copies -= 1
                elif r < r_dup + r_los + r_org:
                    k = int(rng.integers(k_alive))  # origination is per lineage: a uniform lineage
                    _originate(gen[k], tree.nodes[alive[k]], t, events, new_copy, new_family)
                    total_copies += 1
                else:
                    total_copies += _do_transfer(rng, tree, alive, gen, n, t, events, new_copy,
                                                 transfer_to, replacement, self_transfer, depth)
                continue

        if horizon == next_species:  # advance to the tree's next event(s); process the whole tie-batch
            t = next_species
            while si < len(schedule) and schedule[si][0] == t:
                i = schedule[si][1]
                g = gen[pos[i]]
                genomes[i] = tuple(g)  # finalise this lineage (extant, extinct, or unsampled)
                total_copies -= len(g)
                retire(alive, gen, pos, pos[i])
                node = tree.nodes[i]
                if node.children is not None:  # a speciation: each gene re-ids into each daughter
                    for c in node.children:
                        child_genome = []
                        for old in g:  # ZOMBI1: the gene ends here and continues under a fresh id
                            nc = new_copy(old.family)
                            child_genome.append(nc)
                            events.append(Event(t, "speciation", c, old.family, nc.id, parent=old.id))
                        enter(alive, gen, pos, c, child_genome)
                        total_copies += len(child_genome)
                si += 1
        else:
            t = horizon  # a skyline breakpoint: advance and re-evaluate the (now changed) rate

    return GenomesResult(tree, genomes, events, seed)


__all__ = ["simulate_genomes_unordered", "GenomesResult", "Event", "GeneCopy", "Distance",
           "Profiles", "GeneTree", "GeneNode",
           "simulate_genomes_ordered", "OrderedGenomesResult", "Gene", "Chromosome",
           "ChromosomeEvent", "Inversion", "Transposition", "Translocation"]
