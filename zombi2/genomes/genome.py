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

import copy
from abc import ABC, abstractmethod
from dataclasses import dataclass, field, replace

import numpy as np

from zombi2.genomes.events import (
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
        self._order = 0
        self._chromosome = 0

    def new_gene(self) -> str:
        self._gene += 1
        return f"g{self._gene}"

    def new_family(self) -> str:
        self._family += 1
        return str(self._family)

    def new_chromosome(self) -> int:
        """A fresh, globally-unique chromosome id.

        A **separate id space** from genes / families / order: minting a chromosome never advances
        the gene, family or order counters, so introducing the chromosome tier cannot shift a gene
        id and a single-chromosome run stays byte-identical. Chromosomes are re-minted at speciation
        (like gene copies), so every chromosome instance has a unique id across the whole run."""
        self._chromosome += 1
        return self._chromosome

    def new_order(self) -> int:
        """A fresh, monotonically increasing seniority stamp for a genuinely new gene copy.

        Used by :attr:`Gene.origin_order`: a lineage's *continuation* inherits its parent's stamp,
        while every truly new copy (origination, duplication copy, transferred-in copy, converted
        copy) gets a larger one. So the founder lineage of a family always holds the minimum stamp
        among its live copies — the anchor for biased (directional) gene conversion."""
        self._order += 1
        return self._order


@dataclass(slots=True)
class Gene:
    """A single gene copy. v1 carries ``gid`` (current lineage segment) and ``family``.

    The ordered-genome extension adds ``orientation``; the length extension adds
    ``length`` — additive fields on a subclass.

    ``origin_order`` is an inherited **seniority** stamp (lower = older lineage): a lineage's
    continuation keeps its parent's stamp, every truly new copy gets a larger one (see
    :meth:`IdManager.new_order`), so the founder copy of a family always has the minimum. It is
    keyword-only (it never participates in positional construction, so subclasses like
    :class:`OrderedGene` are unaffected) and is used only to bias directional gene conversion; it
    is pure bookkeeping and appears in no output, so it leaves every existing run byte-identical.
    """

    gid: str
    family: str
    origin_order: int = field(default=0, kw_only=True)


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

    def snapshot(self) -> "Genome":
        """A frozen, independent copy at the current state, keeping the SAME gene ids (no
        re-minting). Records the ancestral genome at a degree-two (sampled-ancestor) node while
        the live lineage keeps mutating. Subclasses override with an efficient same-id copy that
        shares the IdManager; the default deep-copies."""
        return copy.deepcopy(self)


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
        return frozenset(STOCHASTIC_EVENTS + (EventType.CONVERSION,))

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
            left = Gene(self.ids.new_gene(), fam, origin_order=parent.origin_order)  # continuation
            right = Gene(self.ids.new_gene(), fam, origin_order=self.ids.new_order())  # new copy
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
        gene = Gene(self.ids.new_gene(), family, origin_order=self.ids.new_order())
        self._add(gene)
        return [GeneOp(gene.gid, family, "origin")]

    # --- intra-genome gene conversion --------------------------------------
    def convert(self, family: str, rng, bias: float = 0.0) -> tuple[list[GeneOp], list[GeneOp]]:
        """One copy of ``family`` overwrites ("converts") another copy of the same family in place.

        Non-reciprocal and copy-number-neutral: the **donor** lineage bifurcates (exactly like a
        duplication of the donor — its old segment ends, a continuation and a fresh converted copy
        open), and the **recipient** copy's old lineage is removed. So the converted copy descends
        from the donor and the two coalesce at the conversion time, while the recipient's ancestry
        above the event is gone — the concerted-evolution signal. The recipient is chosen uniformly;
        ``bias`` in [0, 1] tilts the donor toward the family's founder (smallest ``origin_order``):
        ``0`` = uniform donor, ``1`` = always the oldest candidate. ``bias`` is inert when only two
        copies exist (one donor candidate). Requires ``copy_number(family) >= 2`` (the rate model
        only fires it then). Returns ``(donor_bifurcation_group, recipient_loss_group)``.
        """
        pool = self._genes[family]
        n = len(pool)
        recipient = pool[int(rng.integers(n))]           # the copy being overwritten
        others = [g for g in pool if g is not recipient]
        donor = self._choose_donor(others, rng, bias)
        # donor bifurcates: old segment ends, continuation + converted copy open (a duplication)
        self._remove(donor)
        donor_cont = Gene(self.ids.new_gene(), family, origin_order=donor.origin_order)
        converted = Gene(self.ids.new_gene(), family, origin_order=self.ids.new_order())
        self._add(donor_cont)
        self._add(converted)
        # the recipient's old lineage is overwritten -> it ends here (a loss of that ancestry)
        self._remove(recipient)
        return (
            [GeneOp(donor.gid, family, "parent"),
             GeneOp(donor_cont.gid, family, "donor_copy"),
             GeneOp(converted.gid, family, "converted_copy")],
            [GeneOp(recipient.gid, family, "converted_out")],
        )

    @staticmethod
    def _choose_donor(others, rng, bias):
        """Pick the donor among ``others``: uniform (``bias=0`` or a single candidate), else with
        probability ``bias`` the oldest lineage (smallest ``origin_order``). The ``bias<=0`` /
        single-candidate path draws exactly one uniform integer, so unbiased runs stay
        byte-identical to a bias-free implementation."""
        if bias > 0.0 and len(others) > 1 and rng.random() < bias:
            return min(others, key=lambda g: g.origin_order)
        return others[int(rng.integers(len(others)))]

    # --- transfer handoff --------------------------------------------------
    def extract_segment(self, selection: Selection, rng) -> TransferSegment:
        old_gids, cont_gids, transferred = [], [], []
        for g in selection.genes:
            fam = g.family
            old_gids.append(g.gid)
            self._remove(g)
            # donor lineage continues (new segment, keeps its seniority)
            cont = Gene(self.ids.new_gene(), fam, origin_order=g.origin_order)
            self._add(cont)
            cont_gids.append(cont.gid)
            # the transferred copy is a newcomer in the recipient -> a fresh seniority stamp
            transferred.append(Gene(self.ids.new_gene(), fam, origin_order=self.ids.new_order()))
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
                ng = Gene(self.ids.new_gene(), family, origin_order=g.origin_order)  # keeps seniority
                new._add(ng)
                mapping.append((g.gid, ng.gid, family))
        return new, mapping

    def snapshot(self) -> "UnorderedGenome":
        new = copy.copy(self)                       # shares the IdManager (a snapshot never mints)
        new._genes = {fam: [replace(g) for g in lst] for fam, lst in self._genes.items()}
        return new


@dataclass(slots=True)
class OrderedGene(Gene):
    """A gene copy on an ordered chromosome: adds strand ``orientation`` (+1 / -1)."""

    orientation: int = 1


@dataclass(slots=True)
class Chromosome:
    """One chromosome: an identified, topology-aware ordered run of elements.

    The **shared container** for the chromosome tier: an identity (``chrom_id``), a per-chromosome
    ``circular`` topology (circular wraps at the origin — bacteria; linear has ends — eukaryotes),
    and an ordered ``elements`` list. ``elements`` are :class:`OrderedGene`s for the ordered model
    and nucleotide ``Segment``s for the nucleotide model — the container is the same, and each model
    brings its own coordinate arithmetic ("where is position X"). ``genes`` is a backward-compatible
    alias for ``elements`` (the ordered case, where elements are genes).

    ``chrom_id`` comes from :meth:`IdManager.new_chromosome` — a *separate* id space from genes and
    families, so introducing chromosomes never shifts a gene id.

    The coordinate-agnostic list moves (``bring_to_front`` / ``insert`` / ``remove``) are shared by
    both models; ``segment`` / ``invert_block`` are the ordered model's index-based helpers (the
    nucleotide model, working in base pairs with segment-splitting, supplies its own).
    """

    chrom_id: int
    circular: bool
    elements: list = field(default_factory=list)  # OrderedGene (ordered) or Segment (nucleotide)

    @property
    def genes(self) -> list:
        """Backward-compatible alias for :attr:`elements` (the ordered model's elements are genes)."""
        return self.elements

    @genes.setter
    def genes(self, value: list) -> None:
        self.elements = value

    def __len__(self) -> int:
        return len(self.elements)

    # --- list moves (coordinate-agnostic; shared by both models) -----------
    def bring_to_front(self, start: int) -> int:
        """Rotate the ring so the run at index ``start`` begins at index 0, returning the new start.

        A circular chromosome has no privileged origin, so rotating it is a structural no-op that
        makes a (possibly wrapped) run contiguous for slicing, and it returns 0. A linear chromosome
        is left untouched and ``start`` is returned unchanged."""
        if self.circular:
            if start:
                self.elements = self.elements[start:] + self.elements[:start]
            return 0
        return start

    def insert(self, pos: int, items) -> None:
        """Insert ``items`` in order so the first lands at index ``pos``."""
        self.elements[pos:pos] = list(items)

    def remove(self, item) -> bool:
        """Remove ``item`` by identity if present; return whether it was found."""
        for i, x in enumerate(self.elements):
            if x is item:
                del self.elements[i]
                return True
        return False

    # --- ordered index-based helpers (the ordered model's "find / cut") ----
    def segment(self, start: int, length: int) -> tuple[tuple, int]:
        """Read the contiguous ``length``-element run at ``start``, as ``(elements, length)``.

        A circular chromosome follows the ring, so the run may wrap the origin; a linear one clamps
        it to the end (the returned length is the clamped one). Read-only. Index-based — the ordered
        model's coordinate read (the nucleotide model reads by base-pair arc instead)."""
        n = len(self.elements)
        if self.circular:
            return tuple(self.elements[(start + i) % n] for i in range(length)), length
        length = min(length, n - start)
        return tuple(self.elements[start + i] for i in range(length)), length

    @staticmethod
    def invert_block(items: list) -> list:
        """Invert a block of oriented elements: flip every ``orientation`` (in place) and return it
        reversed — what an inversion does to an ordered segment, and what a reverse-complemented
        transposition reuses."""
        for x in items:
            x.orientation = -x.orientation
        return list(reversed(items))


class OrderedGenome(Genome):
    """ZOMBI-1-style ordered chromosome(s), genes carrying orientation, no intergenes.

    State is a **dict of chromosomes** keyed by identity, ``self.chromosomes: dict[int,
    Chromosome]`` (insertion-ordered); a single-chromosome genome is just ``len == 1``. Each
    :class:`Chromosome` carries its own ``chrom_id`` (minted from a *separate* id space, so it
    never shifts a gene id) and its own ``circular`` topology. ``self.circular`` is the
    genome-wide seed used when constructing / cloning (this stage is topology-homogeneous).

    Events act on a **contiguous segment** whose length is drawn from ``extension``
    (a per-step continuation probability; ``None`` -> single-gene events). Supported:
    duplication (tandem), loss, transfer, **inversion** (reverse + flip strands), and
    **transposition** (cut and paste elsewhere). These gene events stay **within one chromosome**
    (whole-chromosome fission / fusion / origination / loss are the separate chromosome tier —
    :meth:`fission`, :meth:`fusion`, :meth:`originate_chromosome`, :meth:`lose_chromosome`):
    ``draw_target`` first picks a chromosome (size-weighted, returning its ``chrom_id``), then a
    segment within it, and ``apply`` looks that chromosome up by id. A circular chromosome's segment
    may wrap (we rotate the ring to the front); a linear chromosome's segment is clamped to the end.

    The ``extension`` can also be supplied per event via ``TargetParams.extension`` from a
    rate model; the per-genome value is the fallback. Construct via a factory, e.g.
    ``genome_factory=lambda ids: OrderedGenome(ids, extension=0.5, n_chromosomes=8,
    circular=False)`` for eight linear chromosomes.

    ``transposition_flip`` in [0, 1] is the probability that a transposed segment reinserts
    in **reversed** orientation (gene order reversed and every strand flipped), modelling a
    reverse-complemented reinsertion. The default ``0.0`` always preserves orientation, so
    it leaves every existing run byte-identical.

    **Backward compatibility.** ``self.chromosome`` (singular) is a read-only flattened view of
    every gene, so external readers keep working; the per-chromosome structure lives on
    ``self.chromosomes`` (the dict). With ``n_chromosomes=1, circular=True`` (the defaults) every
    RNG draw matches the pre-multichromosome engine, so a run is byte-identical (the chromosome
    choice short-circuits without drawing when there is only one chromosome, and chromosome ids
    come from their own counter).
    """

    def __init__(self, ids: IdManager, extension: float | None = None,
                 transposition_flip: float = 0.0, n_chromosomes: int = 1,
                 circular: bool = True):
        if n_chromosomes < 1:
            raise ValueError("n_chromosomes must be >= 1")
        self.ids = ids
        self.extension = extension
        self.transposition_flip = transposition_flip
        self.circular = circular  # genome-wide seed; Chromosome.circular is the per-chromosome authority
        self.chromosomes: dict[int, Chromosome] = {}
        for _ in range(n_chromosomes):
            cid = ids.new_chromosome()
            self.chromosomes[cid] = Chromosome(cid, circular)

    @property
    def chromosome(self) -> list[OrderedGene]:
        """Read-only flattened view of all genes across every chromosome (backward compat).

        A single-chromosome genome returns its one inner gene list directly; otherwise a fresh
        concatenation. Mutating state goes through ``self.chromosomes`` (only this class does).
        """
        if len(self.chromosomes) == 1:
            return next(iter(self.chromosomes.values())).genes
        return [g for chrom in self.chromosomes.values() for g in chrom.genes]

    # --- queries -----------------------------------------------------------
    def families(self) -> list[str]:
        return list(dict.fromkeys(g.family for chrom in self.chromosomes.values()
                                  for g in chrom.genes))

    def copy_number(self, family: str) -> int:
        return sum(1 for chrom in self.chromosomes.values() for g in chrom.genes
                   if g.family == family)

    def size(self) -> int:
        return sum(len(chrom.genes) for chrom in self.chromosomes.values())

    def total_length(self) -> float:
        return float(self.size())  # no intergenes: length == gene count

    def genes(self) -> list[OrderedGene]:
        return [g for chrom in self.chromosomes.values() for g in chrom.genes]

    def presence_vector(self, family_order) -> np.ndarray:
        present = {g.family for chrom in self.chromosomes.values() for g in chrom.genes}
        return np.fromiter((1 if f in present else 0 for f in family_order),
                           dtype=np.int8, count=len(family_order))

    def supported_events(self) -> frozenset[EventType]:
        return frozenset(STOCHASTIC_EVENTS + (
            EventType.INVERSION, EventType.TRANSPOSITION,
            EventType.CHROMOSOME_ORIGINATION, EventType.CHROMOSOME_LOSS,
            EventType.FISSION, EventType.FUSION,
        ))

    # --- helpers -----------------------------------------------------------
    def _segment_length(self, rng, chrom: "Chromosome", params: TargetParams) -> int:
        ext = params.extension if params.extension is not None else self.extension
        n = len(chrom.genes)
        if ext is None or n <= 1:
            return 1
        length = 1
        while length < n and rng.random() < ext:
            length += 1
        return length

    def _choose_chromosome_weighted(self, rng) -> int:
        """Pick a chromosome (return its ``chrom_id``), weighted by its gene count.

        Short-circuits without drawing when there is only one chromosome, so a single-chromosome
        genome consumes exactly the RNG draws of the pre-multichromosome engine (byte-identity)."""
        if len(self.chromosomes) == 1:
            return next(iter(self.chromosomes))
        chroms = list(self.chromosomes.values())
        sizes = [len(c.genes) for c in chroms]
        r = int(rng.integers(sum(sizes)))
        for c, s in zip(chroms, sizes):
            if r < s:
                return c.chrom_id
            r -= s
        return chroms[-1].chrom_id  # pragma: no cover -- unreachable while total > 0

    def _choose_chromosome_uniform(self, rng) -> int:
        """Pick a chromosome uniformly (return its ``chrom_id``), for origination / transfer
        insertion where empty chromosomes are valid targets. Draws nothing when there is one."""
        if len(self.chromosomes) == 1:
            return next(iter(self.chromosomes))
        cids = list(self.chromosomes)
        return cids[int(rng.integers(len(cids)))]

    def _remove_gene(self, gene: OrderedGene) -> None:
        """Remove ``gene`` (by identity) from whichever chromosome holds it."""
        for chrom in self.chromosomes.values():
            if chrom.remove(gene):
                return
        raise ValueError("gene not present in any chromosome")  # pragma: no cover

    # --- mutation ----------------------------------------------------------
    def draw_target(self, event, rng, params, family=None) -> Selection:
        if family is None:
            cid = self._choose_chromosome_weighted(rng)  # size-weighted; no draw if 1 chromosome
            chrom = self.chromosomes[cid]
            start = int(rng.integers(len(chrom.genes)))
        else:
            # a uniform occurrence of the family across all chromosomes (implicitly weights a
            # chromosome by how many copies it carries); one integer draw, as before
            occ = [(cix, i) for cix, chrom in self.chromosomes.items()
                   for i, g in enumerate(chrom.genes) if g.family == family]
            cid, start = occ[int(rng.integers(len(occ)))]
            chrom = self.chromosomes[cid]
        length = self._segment_length(rng, chrom, params)
        genes, length = chrom.segment(start, length)  # topology-aware read (wrap / clamp)
        return Selection(genes=genes, region=Region(chromosome=cid, start=start, length=length))

    def apply(self, event, selection, rng, params) -> list[list[GeneOp]]:
        if event is EventType.LOSS:  # may be a non-contiguous single gene (replacement)
            groups = []
            for g in selection.genes:
                self._remove_gene(g)
                groups.append([GeneOp(g.gid, g.family, "lost")])
            return groups

        # contiguous segment operations within one chromosome
        region = selection.region
        chrom = self.chromosomes[region.chromosome]
        length = region.length
        start = chrom.bring_to_front(region.start)  # circular: rotate to front (=0); linear: unchanged
        end = start + length
        segment = chrom.genes[start:end]

        if event is EventType.DUPLICATION:
            groups, copies = [], []
            for i, g in enumerate(segment):
                cont = OrderedGene(self.ids.new_gene(), g.family, g.orientation)
                chrom.genes[start + i] = cont  # ancestral lineage continues (re-minted)
                copy = OrderedGene(self.ids.new_gene(), g.family, g.orientation)
                copies.append(copy)
                groups.append([GeneOp(g.gid, g.family, "parent"),
                               GeneOp(cont.gid, g.family, "left"),
                               GeneOp(copy.gid, g.family, "right")])
            chrom.insert(end, copies)  # tandem block right after the segment
            return groups

        if event is EventType.INVERSION:
            chrom.genes[start:end] = Chromosome.invert_block(segment)  # reverse + flip strands
            return [[GeneOp(g.gid, g.family, "inverted") for g in segment]]

        if event is EventType.TRANSPOSITION:
            block = list(segment)
            del chrom.genes[start:end]
            # With probability ``transposition_flip`` the segment reinserts reverse-complemented
            # (gene order reversed and every strand flipped). The ``and`` guards the draw so
            # ``transposition_flip == 0`` never consumes an ``rng`` value — an orientation-preserving
            # run stays byte-identical to a flip-free implementation (the same pattern as biased gene
            # conversion in :meth:`UnorderedGenome._choose_donor`).
            if self.transposition_flip and rng.random() < self.transposition_flip:
                block = Chromosome.invert_block(block)
            j = int(rng.integers(len(chrom.genes) + 1))  # destination within the SAME chromosome
            chrom.insert(j, block)
            return [[GeneOp(g.gid, g.family, "transposed") for g in block]]

        raise ValueError(f"apply() does not handle {event!r}")

    def originate(self, rng, params) -> list[GeneOp]:
        family = self.ids.new_family()
        chrom = self.chromosomes[self._choose_chromosome_uniform(rng)]  # no draw if 1 chromosome
        gene = OrderedGene(self.ids.new_gene(), family, 1 if rng.random() < 0.5 else -1)
        chrom.insert(int(rng.integers(len(chrom.genes) + 1)), [gene])
        return [GeneOp(gene.gid, family, "origin")]

    # --- transfer handoff --------------------------------------------------
    def extract_segment(self, selection, rng) -> TransferSegment:
        region = selection.region
        chrom = self.chromosomes[region.chromosome]
        length = region.length
        start = chrom.bring_to_front(region.start)  # circular: rotate to front (=0); linear: unchanged
        segment = chrom.genes[start:start + length]
        old_gids, cont_gids, transferred = [], [], []
        for i, g in enumerate(segment):
            old_gids.append(g.gid)
            cont = OrderedGene(self.ids.new_gene(), g.family, g.orientation)
            chrom.genes[start + i] = cont
            cont_gids.append(cont.gid)
            transferred.append(OrderedGene(self.ids.new_gene(), g.family, g.orientation))
        return TransferSegment(family=selection.family, genes=tuple(transferred),
                               donor_old_gids=old_gids, donor_cont_gids=cont_gids)

    def choose_insertion_point(self, segment, rng) -> tuple[int, int]:
        cid = self._choose_chromosome_uniform(rng)  # recipient chromosome; no draw if 1 chromosome
        return (cid, int(rng.integers(len(self.chromosomes[cid].genes) + 1)))

    def insert_segment(self, segment, at, rng) -> list[GeneOp]:
        if isinstance(at, tuple):
            cid, j = at
        else:  # fallback: no precomputed point -> choose one now
            cid = self._choose_chromosome_uniform(rng)
            j = int(rng.integers(len(self.chromosomes[cid].genes) + 1))
        self.chromosomes[cid].insert(j, segment.genes)
        return [GeneOp(g.gid, g.family, "transfer_copy") for g in segment.genes]

    # --- speciation --------------------------------------------------------
    def clone_reminting(self) -> tuple["OrderedGenome", list[tuple[str, str, str]]]:
        new = type(self)(self.ids, self.extension, self.transposition_flip,
                         n_chromosomes=len(self.chromosomes), circular=self.circular)
        mapping = []
        for pchrom, nchrom in zip(self.chromosomes.values(), new.chromosomes.values()):
            nchrom.circular = pchrom.circular  # per-chromosome topology (homogeneous today)
            for g in pchrom.genes:
                ng = OrderedGene(self.ids.new_gene(), g.family, g.orientation)
                nchrom.genes.append(ng)
                mapping.append((g.gid, ng.gid, g.family))
        return new, mapping

    def snapshot(self) -> "OrderedGenome":
        new = copy.copy(self)                       # shares the IdManager and config
        new.chromosomes = {cid: replace(chrom, elements=[replace(g) for g in chrom.elements])
                           for cid, chrom in self.chromosomes.items()}
        return new

    # --- chromosome-tier events (origination / loss / fission / fusion) -----
    #
    # These act on whole chromosomes, one tier above the gene events. Fission / fusion only move
    # genes between chromosomes, so gene lineages (gids) are untouched and gene-tree reconstruction
    # is unaffected; only chromosome LOSS ends gene lineages, so it reports them for the log.
    def originate_chromosome(self, rng, params) -> int:
        """Originate a new, empty chromosome — a de-novo replicon (a *plasmid* in the bacterial
        case), with a fresh id and the genome's topology. Genes reach it later via origination or
        transfer (both can target any chromosome). Returns the new chrom_id. Draws no rng."""
        cid = self.ids.new_chromosome()
        self.chromosomes[cid] = Chromosome(cid, self.circular)
        return cid

    def lose_chromosome(self, rng) -> tuple[int, list[list[GeneOp]]]:
        """Lose an entire chromosome (chosen uniformly) and everything on it, returning its chrom_id
        and one LOSS group per gene so the caller can log the gene deaths for gene-tree
        reconstruction. The caller guarantees >= 2 chromosomes (a genome keeps at least one)."""
        cid = self._choose_chromosome_uniform(rng)
        chrom = self.chromosomes.pop(cid)
        groups = [[GeneOp(g.gid, g.family, "lost")] for g in chrom.genes]
        return cid, groups

    def fission(self, rng, params) -> tuple[int, int]:
        """Split one chromosome into two (genes keep their ids — fission only reorganises which
        chromosome holds them). A **linear** chromosome is cut at ONE breakpoint (the prefix stays,
        the suffix becomes a new chromosome); a **circular** one at TWO breakpoints, excising the arc
        between them into a new circular replicon (the remainder stays circular). The source is
        size-weighted, so bigger chromosomes fission more. Returns (source_cid, new_cid)."""
        cid = self._choose_chromosome_weighted(rng)
        chrom = self.chromosomes[cid]
        n = len(chrom.genes)
        new_cid = self.ids.new_chromosome()
        if chrom.circular:
            i, j = int(rng.integers(n + 1)), int(rng.integers(n + 1))
            if i > j:
                i, j = j, i
            arc = chrom.genes[i:j]
            chrom.genes = chrom.genes[:i] + chrom.genes[j:]
            self.chromosomes[new_cid] = Chromosome(new_cid, True, arc)
        else:
            k = int(rng.integers(n + 1))
            self.chromosomes[new_cid] = Chromosome(new_cid, False, chrom.genes[k:])
            chrom.genes = chrom.genes[:k]
        return cid, new_cid

    def fusion(self, rng, params) -> tuple[int, int]:
        """Fuse two chromosomes into one: append the second's genes onto the first and drop the
        second (genes keep their ids). The two are chosen uniformly; the caller guarantees >= 2
        chromosomes. Genomes are topology-homogeneous today, so the same-topology rule is met
        automatically. Returns (kept_cid, absorbed_cid)."""
        cids = list(self.chromosomes)
        m = len(cids)
        i = int(rng.integers(m))
        j = int(rng.integers(m - 1))
        if j >= i:  # a uniform partner distinct from i
            j += 1
        keep_cid, absorb_cid = cids[i], cids[j]
        self.chromosomes[keep_cid].genes.extend(self.chromosomes.pop(absorb_cid).genes)
        return keep_cid, absorb_cid
