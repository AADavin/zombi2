"""Genomes I — the unordered D/L/O gene-family core.

A genome is a multiset of gene families that evolves along the species tree by three events:
**origination** (a new family arises in a lineage — per lineage), **duplication** (a gene copy
duplicates, +1 copy in its family — per copy), and **loss** (a gene copy is lost — per copy). Rates
are the cross-level ``scope(base) × modifiers`` grammar (``SPEC §5``); the defaults are the natural
"per what?" for each event.

This is the genome twin of :mod:`zombi2.species_tree`, and it reads the same way — a forward walk
over the **complete** tree, plain/frozen dataclasses, an event log as the source of truth (per-family
gene trees are derived from it later), ``as_rate``/``.effective`` for every rate. The one structural
difference: with **no transfer**, lineages are conditionally independent given the parent genome, so
the global birth-death race of :func:`~zombi2.species_tree._grow` collapses to a per-segment Gillespie
inside a plain pre-order tree walk. Transfer (the next slice) reintroduces cross-lineage coupling and
grows the walk back toward a global timeline — a deliberate seam.

Still to come: transfer, per-family heterogeneity (``ByFamily`` + ``Speed``), the sparse profiles and
lazy gene-tree views behind the ``record=`` memory dial, and the Rust core. This lives here for now
so the legacy ``zombi2/genomes`` package is untouched.
"""

from __future__ import annotations

import collections
from dataclasses import dataclass

import numpy as np

from .modifiers import Time
from .rate import as_rate
from .scope import PerCopy, PerLineage
from .species_tree import SpeciesResult, Tree


@dataclass(frozen=True)
class GeneCopy:
    """One gene copy: a member of family ``family``, identified by a globally-unique ``id``. Its
    birth/death times and parentage live in the event log (the source of truth); the copy carries
    only what a genome snapshot needs to be self-describing — who it is and which family it is in. A
    genome may hold several copies sharing a ``family`` (that family's copy count)."""

    id: int
    family: int


@dataclass(frozen=True)
class Event:
    """A recorded genome event — the true history every per-family gene tree is later derived from.
    ``lineage`` is the species-tree node id it fired on; ``time`` is when (crown-forward, the same
    clock as the species tree). By kind:

    - ``"origination"`` — ``copy`` is the family's founding copy (``parent`` is ``None``).
    - ``"duplication"`` — ``copy`` is the NEW copy; ``parent`` is the copy that duplicated and survives.
    - ``"loss"`` — ``copy`` is the copy removed (``parent`` is ``None``).
    """

    time: float
    kind: str  # "origination" | "duplication" | "loss"
    lineage: int  # the species-tree node id where it fired
    family: int
    copy: int  # the copy born (origination / duplication) or removed (loss)
    parent: int | None = None  # duplication only: the copy that spawned `copy`


@dataclass
class GenomesResult:
    """What ``simulate_unordered`` returns: the ``complete_tree`` it ran on, the final ``genomes`` at
    **every** node (extant and extinct), the ``events`` log (the compact source of truth), and the
    ``seed``. (Sparse profiles, lazy gene trees, the ``record=`` dial and ``write`` are later slices.)"""

    complete_tree: Tree
    genomes: dict[int, tuple[GeneCopy, ...]]
    events: list[Event]
    seed: int | None

    @property
    def extant_genomes(self) -> dict[int, tuple[GeneCopy, ...]]:
        """The observed genomes — the content at the extant tips only."""
        return {n.id: self.genomes[n.id] for n in self.complete_tree.extant()}

    def family_counts(self, node_id: int) -> collections.Counter:
        """A multiset view of one node's genome: ``family id → copy count``."""
        return collections.Counter(c.family for c in self.genomes[node_id])


def _originate(genome, node, t, events, new_copy, new_family) -> None:
    """A new gene family arises: mint a founding copy in a fresh family and record it."""
    c = new_copy(new_family())
    genome.append(c)
    events.append(Event(t, "origination", node.id, c.family, c.id))


def _duplicate(rng, genome, node, t, events, new_copy) -> None:
    """A uniformly-chosen copy duplicates: +1 copy in its family; the parent copy survives."""
    p = genome[int(rng.integers(len(genome)))]  # uniform: every copy carries equal per-copy weight
    c = new_copy(p.family)
    genome.append(c)
    events.append(Event(t, "duplication", node.id, p.family, c.id, parent=p.id))


def _lose(rng, genome, node, t, events) -> None:
    """A uniformly-chosen copy is lost (swap-remove — the genome is an order-agnostic multiset)."""
    j = int(rng.integers(len(genome)))
    lost = genome[j]
    genome[j] = genome[-1]
    genome.pop()
    events.append(Event(t, "loss", node.id, lost.family, lost.id))


def _evolve_segment(rng, genome, node, dup, los, org, events, new_copy, new_family) -> None:
    """Play D/L/O forward along one lineage segment (``birth_time..end_time`` on the complete tree)
    by an exact Gillespie: the total rate holds until the next skyline breakpoint or the segment end,
    so draw a waiting time, pick the event in proportion to its rate, apply it, and re-draw."""
    t, end = node.birth_time, node.end_time
    while True:
        n = len(genome)
        ctx = {"copies": n, "lineages": 1, "time": t}  # one lineage per genome
        r_dup = dup.effective(**ctx) if n else 0.0  # per copy → nothing to duplicate on a void
        r_los = los.effective(**ctx) if n else 0.0  # per copy → nothing to lose on a void
        r_org = org.effective(**ctx)  # per lineage → size-independent (can populate an empty genome)
        total = r_dup + r_los + r_org
        # the total is constant until the next skyline breakpoint (or the segment end)
        horizon = min(end, dup.next_change(t), los.next_change(t), org.next_change(t))

        if total > 0.0:
            t_event = t + float(rng.exponential(1.0 / total))
            if t_event < horizon:  # an event fires inside this rate regime
                t = t_event
                r = float(rng.random()) * total
                if r < r_dup:
                    _duplicate(rng, genome, node, t, events, new_copy)
                elif r < r_dup + r_los:
                    _lose(rng, genome, node, t, events)
                else:
                    _originate(genome, node, t, events, new_copy, new_family)
                continue

        if horizon >= end:  # reached the segment end / the present
            break
        t = horizon  # a skyline breakpoint: advance and re-evaluate the (now changed) rate


def simulate_unordered(tree, *, duplication=0.0, loss=0.0, origination=0.0,
                       initial_families=0, seed=None) -> GenomesResult:
    """Evolve a multiset of gene families along a species tree by duplication, loss, and origination.

    ``tree`` is the **complete** species tree (a :class:`~zombi2.species_tree.Tree`, or a
    :class:`~zombi2.species_tree.SpeciesResult` whose ``complete_tree`` is used). Genomes evolve on
    **every** lineage, extant and extinct alike, so the true gene-tree history is complete; the
    observed genomes are the extant tips (``result.extant_genomes``).

    ``duplication`` and ``loss`` default to **per copy**, ``origination`` to **per lineage** (each is a
    ``scope(base) × modifiers`` rate spec — a number, a scope wrapper, or a product). The root starts
    with ``initial_families`` families of one copy each, recorded as originations at the crown so the
    event log is the whole story. Deterministic given ``seed``.
    """
    tree = tree.complete_tree if isinstance(tree, SpeciesResult) else tree
    dup = as_rate(duplication, default_scope=PerCopy)
    los = as_rate(loss, default_scope=PerCopy)
    org = as_rate(origination, default_scope=PerLineage)
    # only Time (skyline) is wired for genome rates this slice: the walk supplies `time`, but not the
    # per-lineage/per-family threading a drift or heterogeneity modifier needs. Reject the rest loudly
    # rather than silently drop them (a drift modifier would otherwise no-op) or crash on missing context.
    for label, rate in (("duplication", dup), ("loss", los), ("origination", org)):
        for m in rate.modifiers:
            if not isinstance(m, Time):
                raise ValueError(
                    f"{label} carries {type(m).__name__}, which the unordered genome engine does not "
                    f"support yet — only Time (skyline) is wired. Per-family heterogeneity (ByFamily, "
                    f"Speed) and clade drift are later slices."
                )
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

    genomes: dict[int, tuple[GeneCopy, ...]] = {}
    events: list[Event] = []

    for i in sorted(tree.nodes):  # ascending id = parents before children (children mint id > parent)
        node = tree.nodes[i]
        if node.parent is None:  # the root: seed the initial families as originations at the crown
            genome: list[GeneCopy] = []
            for _ in range(initial_families):
                _originate(genome, node, node.birth_time, events, new_copy, new_family)
        else:  # both children inherit the parent's final content; frozen copies → independent list
            genome = list(genomes[node.parent])
        _evolve_segment(rng, genome, node, dup, los, org, events, new_copy, new_family)
        genomes[i] = tuple(genome)

    return GenomesResult(tree, genomes, events, seed)


__all__ = ["simulate_unordered", "GenomesResult", "Event", "GeneCopy"]
