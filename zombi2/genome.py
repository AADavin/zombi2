"""Genome representations.

The simulator programs against the abstract :class:`Genome` interface; v1 ships one
trivial implementation, :class:`UnorderedGenome` (a multiset of gene families with copy
numbers). Future representations (ordered chromosomes, length-aware genomes) are new
subclasses that implement the same interface — the simulator, rate model, event log and
profile matrix never change.

Design notes tying back to the spec's "six load-bearing signatures":

* ``draw_target`` returns a :class:`~zombi2.events.Selection` (segment-shaped even in
  v1, where it holds a single gene) and accepts a ``params`` object it currently ignores.
* Transfer is a three-method handoff (``extract_segment`` / ``choose_insertion_point`` /
  ``insert_segment``) so a copy can move between two genomes and land at a
  representation-appropriate position.
* ``total_length`` lives on the base interface (v1 returns ``size``) so length-dependent
  rate models type-check against ``Genome``.
* ``supported_events`` gates the simulator loop, so adding inversion/transposition later
  needs no loop change.

Gene identity: speciation does **not** mint new gene ids (a gene lineage keeps its id
across the species tree; a copy is distinguished by the branch it lives on). New ids are
minted only when a gene is *created* — origination, duplication, or the recipient copy
of a transfer — via a shared :class:`IdManager`.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass

import numpy as np

from .events import (
    EventType,
    GeneOp,
    InsertionPoint,
    Selection,
    TargetParams,
    TransferSegment,
    STOCHASTIC_EVENTS,
)


class IdManager:
    """Mints globally unique gene ids and family ids for one simulation.

    Shared across every genome in a run (propagated through :meth:`Genome.clone`), so
    ids are unique and reproducible for a given seed and call order.
    """

    def __init__(self) -> None:
        self._gene = 0
        self._family = 0

    def new_gene(self) -> str:
        self._gene += 1
        return f"g{self._gene}"

    def new_family(self) -> str:
        self._family += 1
        return str(self._family)


@dataclass
class Gene:
    """A single gene copy. v1 carries only ``gid`` and ``family``.

    The ordered-genome extension adds ``orientation``; the length extension adds
    ``length`` — additive fields on a subclass, not a type change.
    """

    gid: str
    family: str


class Genome(ABC):
    """Abstract genome interface (the contract the simulator depends on)."""

    ids: IdManager

    # --- queries -----------------------------------------------------------
    @abstractmethod
    def families(self) -> list[str]:
        """Families with at least one copy present."""

    @abstractmethod
    def copy_number(self, family: str) -> int:
        ...

    @abstractmethod
    def size(self) -> int:
        """Total number of gene copies."""

    @abstractmethod
    def total_length(self) -> float:
        """Total length (v1: equals :meth:`size`, every gene of length 1)."""

    @abstractmethod
    def genes(self) -> list[Gene]:
        ...

    @abstractmethod
    def presence_vector(self, family_order) -> np.ndarray:
        """Binary presence vector ``σ`` over ``family_order`` (the future Potts state)."""

    @abstractmethod
    def supported_events(self) -> frozenset[EventType]:
        """The stochastic events this representation can undergo."""

    # --- mutation ----------------------------------------------------------
    @abstractmethod
    def draw_target(self, event: EventType, rng, params: TargetParams) -> Selection:
        """Choose what a D / L / T event acts on."""

    @abstractmethod
    def apply(self, event: EventType, selection: Selection, rng, params: TargetParams) -> list[GeneOp]:
        """Apply a D or L event in place; return the affected gene rows."""

    @abstractmethod
    def originate(self, rng, params: TargetParams) -> list[GeneOp]:
        """Create a brand-new family (one copy) in place; return its gene row."""

    # --- transfer handoff (three methods) ---------------------------------
    @abstractmethod
    def extract_segment(self, selection: Selection, rng, *, keep_copy: bool = True) -> TransferSegment:
        """Build a portable copy of the selected genes (fresh ids) for a transfer."""

    @abstractmethod
    def choose_insertion_point(self, segment: TransferSegment, rng) -> InsertionPoint:
        """Where an incoming transferred segment should land in this (recipient) genome."""

    @abstractmethod
    def insert_segment(self, segment: TransferSegment, at, rng) -> list[GeneOp]:
        """Insert an incoming transferred segment in place; return its gene rows."""

    # --- speciation --------------------------------------------------------
    @abstractmethod
    def clone(self) -> "Genome":
        """A deep copy for a child branch at a speciation node (ids preserved)."""


class UnorderedGenome(Genome):
    """v1 genome: an order-free multiset of families -> live gene copies."""

    def __init__(self, ids: IdManager):
        self.ids = ids
        self._genes: dict[str, list[Gene]] = {}

    # --- internal helpers --------------------------------------------------
    def _add(self, gene: Gene) -> None:
        self._genes.setdefault(gene.family, []).append(gene)

    def _remove(self, gene: Gene) -> None:
        lst = self._genes[gene.family]
        lst.remove(gene)
        if not lst:
            del self._genes[gene.family]

    # --- queries -----------------------------------------------------------
    def families(self) -> list[str]:
        return list(self._genes.keys())

    def copy_number(self, family: str) -> int:
        return len(self._genes.get(family, ()))

    def size(self) -> int:
        return sum(len(lst) for lst in self._genes.values())

    def total_length(self) -> float:
        return float(self.size())

    def genes(self) -> list[Gene]:
        return [g for lst in self._genes.values() for g in lst]

    def presence_vector(self, family_order) -> np.ndarray:
        return np.fromiter(
            (1 if self.copy_number(f) > 0 else 0 for f in family_order),
            dtype=np.int8,
            count=len(family_order),
        )

    def supported_events(self) -> frozenset[EventType]:
        return frozenset(STOCHASTIC_EVENTS)

    # --- mutation ----------------------------------------------------------
    def draw_target(self, event: EventType, rng, params: TargetParams) -> Selection:
        all_genes = self.genes()  # deterministic order (dict insertion order)
        idx = int(rng.integers(len(all_genes)))
        return Selection(genes=(all_genes[idx],))

    def apply(self, event: EventType, selection: Selection, rng, params: TargetParams) -> list[GeneOp]:
        if event is EventType.DUPLICATION:
            parent = selection.genes[0]
            copy = Gene(self.ids.new_gene(), parent.family)
            self._add(copy)
            return [
                GeneOp(parent.gid, parent.family, "parent"),
                GeneOp(copy.gid, copy.family, "copy"),
            ]
        if event is EventType.LOSS:
            ops = []
            for g in selection.genes:
                self._remove(g)
                ops.append(GeneOp(g.gid, g.family, "lost"))
            return ops
        raise ValueError(f"apply() does not handle {event!r}")

    def originate(self, rng, params: TargetParams) -> list[GeneOp]:
        family = self.ids.new_family()
        gene = Gene(self.ids.new_gene(), family)
        self._add(gene)
        return [GeneOp(gene.gid, family, "origin")]

    # --- transfer handoff --------------------------------------------------
    def extract_segment(self, selection: Selection, rng, *, keep_copy: bool = True) -> TransferSegment:
        new_genes = tuple(Gene(self.ids.new_gene(), g.family) for g in selection.genes)
        if not keep_copy:  # replacement transfer (not used in v1)
            for g in selection.genes:
                self._remove(g)
        return TransferSegment(family=selection.family, genes=new_genes)

    def choose_insertion_point(self, segment: TransferSegment, rng) -> InsertionPoint:
        return InsertionPoint.ANYWHERE  # a multiset has no meaningful position

    def insert_segment(self, segment: TransferSegment, at, rng) -> list[GeneOp]:
        ops = []
        for g in segment.genes:
            self._add(g)
            ops.append(GeneOp(g.gid, g.family, "transferred"))
        return ops

    # --- speciation --------------------------------------------------------
    def clone(self) -> "UnorderedGenome":
        new = UnorderedGenome(self.ids)
        for family, lst in self._genes.items():
            new._genes[family] = [Gene(g.gid, g.family) for g in lst]
        return new
