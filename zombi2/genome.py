"""Genome representations.

The simulator programs against the abstract :class:`Genome` interface; v1 ships one
trivial implementation, :class:`UnorderedGenome` (a multiset of gene families with copy
numbers). Future representations (ordered chromosomes, length-aware genomes) are new
subclasses implementing the same interface.

Gene-lineage identity (for gene-tree reconstruction): a gene copy's ``gid`` is the id of
its **current lineage segment**. Every event that changes a lineage — duplication,
transfer, speciation — *terminates* the incoming segment(s) and *opens* fresh ones with
new ids, so each ``gid`` is born once and dies once. The event log therefore records a
full genealogy (from-id → to-ids) that :mod:`zombi2.reconciliation` turns into gene trees.
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
    """Mints globally unique gene-lineage ids and family ids for one simulation."""

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
    """A single gene copy. v1 carries only ``gid`` (current lineage segment) and ``family``.

    The ordered-genome extension adds ``orientation``; the length extension adds
    ``length`` — additive fields on a subclass.
    """

    gid: str
    family: str


class Genome(ABC):
    """Abstract genome interface (the contract the simulator depends on)."""

    ids: IdManager

    # --- queries -----------------------------------------------------------
    @abstractmethod
    def families(self) -> list[str]: ...
    @abstractmethod
    def copy_number(self, family: str) -> int: ...
    @abstractmethod
    def size(self) -> int: ...
    @abstractmethod
    def total_length(self) -> float: ...
    @abstractmethod
    def genes(self) -> list[Gene]: ...
    @abstractmethod
    def presence_vector(self, family_order) -> np.ndarray: ...
    @abstractmethod
    def supported_events(self) -> frozenset[EventType]: ...

    # --- mutation ----------------------------------------------------------
    @abstractmethod
    def draw_target(self, event: EventType, rng, params: TargetParams, family: str | None = None) -> Selection:
        """Choose what a D / L / T event acts on (optionally restricted to one family)."""

    @abstractmethod
    def apply(self, event: EventType, selection: Selection, rng, params: TargetParams) -> list[GeneOp]:
        """Apply a D or L event in place; return the affected gene rows (from, to...)."""

    @abstractmethod
    def originate(self, rng, params: TargetParams) -> list[GeneOp]:
        """Create a brand-new family (one copy) in place; return its origin row."""

    # --- transfer handoff --------------------------------------------------
    @abstractmethod
    def extract_segment(self, selection: Selection, rng) -> TransferSegment:
        """Bifurcate the donor lineage: re-mint its continuation and build the copy to send."""

    @abstractmethod
    def choose_insertion_point(self, segment: TransferSegment, rng) -> InsertionPoint: ...

    @abstractmethod
    def insert_segment(self, segment: TransferSegment, at, rng) -> list[GeneOp]:
        """Insert an incoming transferred segment in place; return its gene rows."""

    # --- speciation --------------------------------------------------------
    @abstractmethod
    def clone_reminting(self) -> tuple["Genome", list[tuple[str, str, str]]]:
        """A child copy for a speciation, with re-minted ids.

        Returns the child genome and a list of ``(old_gid, new_gid, family)`` so the
        driver can log the speciation as a lineage bifurcation.
        """


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
            dtype=np.int8, count=len(family_order),
        )

    def supported_events(self) -> frozenset[EventType]:
        return frozenset(STOCHASTIC_EVENTS)

    # --- mutation ----------------------------------------------------------
    def draw_target(self, event: EventType, rng, params: TargetParams, family: str | None = None) -> Selection:
        pool = self.genes() if family is None else self._genes[family]
        idx = int(rng.integers(len(pool)))
        return Selection(genes=(pool[idx],))

    def apply(self, event: EventType, selection: Selection, rng, params: TargetParams) -> list[GeneOp]:
        if event is EventType.DUPLICATION:
            parent = selection.genes[0]
            fam = parent.family
            self._remove(parent)  # the ancestral lineage terminates here
            left = Gene(self.ids.new_gene(), fam)   # continuation
            right = Gene(self.ids.new_gene(), fam)  # new copy
            self._add(left)
            self._add(right)
            return [
                GeneOp(parent.gid, fam, "parent"),
                GeneOp(left.gid, fam, "left"),
                GeneOp(right.gid, fam, "right"),
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
    def extract_segment(self, selection: Selection, rng) -> TransferSegment:
        old_gids, cont_gids, transferred = [], [], []
        for g in selection.genes:
            fam = g.family
            old_gids.append(g.gid)
            self._remove(g)
            cont = Gene(self.ids.new_gene(), fam)  # donor lineage continues (new segment)
            self._add(cont)
            cont_gids.append(cont.gid)
            transferred.append(Gene(self.ids.new_gene(), fam))
        return TransferSegment(
            family=selection.family, genes=tuple(transferred),
            donor_old_gids=old_gids, donor_cont_gids=cont_gids,
        )

    def choose_insertion_point(self, segment: TransferSegment, rng) -> InsertionPoint:
        return InsertionPoint.ANYWHERE  # a multiset has no meaningful position

    def insert_segment(self, segment: TransferSegment, at, rng) -> list[GeneOp]:
        ops = []
        for g in segment.genes:
            self._add(g)
            ops.append(GeneOp(g.gid, g.family, "transfer_copy"))
        return ops

    # --- speciation --------------------------------------------------------
    def clone_reminting(self) -> tuple["UnorderedGenome", list[tuple[str, str, str]]]:
        new = type(self)(self.ids)  # same representation (works for subclasses)
        mapping = []
        for family, lst in self._genes.items():
            for g in lst:
                ng = Gene(self.ids.new_gene(), family)
                new._add(ng)
                mapping.append((g.gid, ng.gid, family))
        return new, mapping
