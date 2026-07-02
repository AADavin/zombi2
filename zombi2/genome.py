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
    Region,
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


@dataclass(slots=True)
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
    def apply(self, event: EventType, selection: Selection, rng, params: TargetParams) -> list[list[GeneOp]]:
        """Apply a D/L/I/P event in place; return one gene-op *group* per event to log.

        A single-gene event yields one group; a segment event (ordered genomes) yields one
        group per gene it touches (each an independent genealogical edge). Each group's
        first row is the ``from`` lineage, the rest its ``to`` lineages.
        """

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
        self._size = 0  # cached total copy count (kept O(1) by _add/_remove)

    # --- internal helpers --------------------------------------------------
    def _add(self, gene: Gene) -> None:
        self._genes.setdefault(gene.family, []).append(gene)
        self._size += 1

    def _remove(self, gene: Gene) -> None:
        lst = self._genes[gene.family]
        lst.remove(gene)
        self._size -= 1
        if not lst:
            del self._genes[gene.family]

    # --- queries -----------------------------------------------------------
    def families(self) -> list[str]:
        return list(self._genes.keys())

    def copy_number(self, family: str) -> int:
        return len(self._genes.get(family, ()))

    def size(self) -> int:
        return self._size

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
        if family is not None:
            pool = self._genes[family]
            return Selection(genes=(pool[int(rng.integers(len(pool)))],))
        # uniform over all copies, without materialising the full gene list
        r = int(rng.integers(self._size))
        for lst in self._genes.values():
            if r < len(lst):
                return Selection(genes=(lst[r],))
            r -= len(lst)
        raise RuntimeError("draw_target on an empty genome")  # guarded by the rate model

    def apply(self, event: EventType, selection: Selection, rng, params: TargetParams) -> list[list[GeneOp]]:
        if event is EventType.DUPLICATION:
            parent = selection.genes[0]
            fam = parent.family
            self._remove(parent)  # the ancestral lineage terminates here
            left = Gene(self.ids.new_gene(), fam)   # continuation
            right = Gene(self.ids.new_gene(), fam)  # new copy
            self._add(left)
            self._add(right)
            return [[
                GeneOp(parent.gid, fam, "parent"),
                GeneOp(left.gid, fam, "left"),
                GeneOp(right.gid, fam, "right"),
            ]]
        if event is EventType.LOSS:
            groups = []
            for g in selection.genes:
                self._remove(g)
                groups.append([GeneOp(g.gid, g.family, "lost")])
            return groups
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


@dataclass(slots=True)
class OrderedGene(Gene):
    """A gene copy on an ordered chromosome: adds strand ``orientation`` (+1 / -1)."""

    orientation: int = 1


class OrderedGenome(Genome):
    """ZOMBI-1-style ordered (circular) chromosome, genes carrying orientation, no
    intergenes.

    Events act on a **contiguous segment** whose length is drawn from ``extension``
    (a per-step continuation probability; ``None`` -> single-gene events). Supported:
    duplication (tandem), loss, transfer, **inversion** (reverse + flip strands), and
    **transposition** (cut and paste elsewhere). The chromosome is circular, so a segment
    may wrap; we handle that by rotating the ring to bring the segment to the front.

    The ``extension`` can also be supplied per event via ``TargetParams.extension`` from a
    rate model; the per-genome value is the fallback. Construct via a factory, e.g.
    ``genome_factory=lambda ids: OrderedGenome(ids, extension=0.5)``.
    """

    def __init__(self, ids: IdManager, extension: float | None = None):
        self.ids = ids
        self.extension = extension
        self.chromosome: list[OrderedGene] = []

    # --- queries -----------------------------------------------------------
    def families(self) -> list[str]:
        return list(dict.fromkeys(g.family for g in self.chromosome))

    def copy_number(self, family: str) -> int:
        return sum(1 for g in self.chromosome if g.family == family)

    def size(self) -> int:
        return len(self.chromosome)

    def total_length(self) -> float:
        return float(len(self.chromosome))  # no intergenes: length == gene count

    def genes(self) -> list[OrderedGene]:
        return list(self.chromosome)

    def presence_vector(self, family_order) -> np.ndarray:
        present = {g.family for g in self.chromosome}
        return np.fromiter((1 if f in present else 0 for f in family_order),
                           dtype=np.int8, count=len(family_order))

    def supported_events(self) -> frozenset[EventType]:
        return frozenset(
            STOCHASTIC_EVENTS + (EventType.INVERSION, EventType.TRANSPOSITION)
        )

    # --- helpers -----------------------------------------------------------
    def _segment_length(self, rng, params: TargetParams) -> int:
        ext = params.extension if params.extension is not None else self.extension
        n = len(self.chromosome)
        if ext is None or n <= 1:
            return 1
        length = 1
        while length < n and rng.random() < ext:
            length += 1
        return length

    def _rotate_to(self, start: int) -> None:
        if start:
            self.chromosome = self.chromosome[start:] + self.chromosome[:start]

    # --- mutation ----------------------------------------------------------
    def draw_target(self, event, rng, params, family=None) -> Selection:
        n = len(self.chromosome)
        if family is None:
            start = int(rng.integers(n))
        else:
            positions = [i for i, g in enumerate(self.chromosome) if g.family == family]
            start = positions[int(rng.integers(len(positions)))]
        length = self._segment_length(rng, params)
        genes = tuple(self.chromosome[(start + i) % n] for i in range(length))
        return Selection(genes=genes, region=Region(chromosome=0, start=start, length=length))

    def apply(self, event, selection, rng, params) -> list[list[GeneOp]]:
        if event is EventType.LOSS:  # may be a non-contiguous single gene (replacement)
            groups = []
            for g in selection.genes:
                self.chromosome.remove(g)
                groups.append([GeneOp(g.gid, g.family, "lost")])
            return groups

        # contiguous segment operations: bring the segment to the front of the ring
        length = selection.region.length
        self._rotate_to(selection.region.start)
        segment = self.chromosome[:length]

        if event is EventType.DUPLICATION:
            groups, copies = [], []
            for i, g in enumerate(segment):
                cont = OrderedGene(self.ids.new_gene(), g.family, g.orientation)
                self.chromosome[i] = cont  # ancestral lineage continues (re-minted)
                copy = OrderedGene(self.ids.new_gene(), g.family, g.orientation)
                copies.append(copy)
                groups.append([GeneOp(g.gid, g.family, "parent"),
                               GeneOp(cont.gid, g.family, "left"),
                               GeneOp(copy.gid, g.family, "right")])
            self.chromosome[length:length] = copies  # tandem block after the segment
            return groups

        if event is EventType.INVERSION:
            for g in segment:
                g.orientation = -g.orientation
            self.chromosome[:length] = list(reversed(segment))
            return [[GeneOp(g.gid, g.family, "inverted") for g in segment]]

        if event is EventType.TRANSPOSITION:
            block = list(segment)
            del self.chromosome[:length]
            j = int(rng.integers(len(self.chromosome) + 1))
            self.chromosome[j:j] = block
            return [[GeneOp(g.gid, g.family, "transposed") for g in block]]

        raise ValueError(f"apply() does not handle {event!r}")

    def originate(self, rng, params) -> list[GeneOp]:
        family = self.ids.new_family()
        gene = OrderedGene(self.ids.new_gene(), family, 1 if rng.random() < 0.5 else -1)
        j = int(rng.integers(len(self.chromosome) + 1))
        self.chromosome[j:j] = [gene]
        return [GeneOp(gene.gid, family, "origin")]

    # --- transfer handoff --------------------------------------------------
    def extract_segment(self, selection, rng) -> TransferSegment:
        length = selection.region.length
        self._rotate_to(selection.region.start)
        segment = self.chromosome[:length]
        old_gids, cont_gids, transferred = [], [], []
        for i, g in enumerate(segment):
            old_gids.append(g.gid)
            cont = OrderedGene(self.ids.new_gene(), g.family, g.orientation)
            self.chromosome[i] = cont
            cont_gids.append(cont.gid)
            transferred.append(OrderedGene(self.ids.new_gene(), g.family, g.orientation))
        return TransferSegment(family=selection.family, genes=tuple(transferred),
                               donor_old_gids=old_gids, donor_cont_gids=cont_gids)

    def choose_insertion_point(self, segment, rng) -> int:
        return int(rng.integers(len(self.chromosome) + 1))

    def insert_segment(self, segment, at, rng) -> list[GeneOp]:
        j = at if isinstance(at, int) else int(rng.integers(len(self.chromosome) + 1))
        self.chromosome[j:j] = list(segment.genes)
        return [GeneOp(g.gid, g.family, "transfer_copy") for g in segment.genes]

    # --- speciation --------------------------------------------------------
    def clone_reminting(self) -> tuple["OrderedGenome", list[tuple[str, str, str]]]:
        new = type(self)(self.ids, self.extension)
        mapping = []
        for g in self.chromosome:
            ng = OrderedGene(self.ids.new_gene(), g.family, g.orientation)
            new.chromosome.append(ng)
            mapping.append((g.gid, ng.gid, g.family))
        return new, mapping
