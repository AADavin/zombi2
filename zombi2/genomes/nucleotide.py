"""Genomes III — nucleotide: a genome as a coordinate space of base pairs.

The third and hardest resolution. Where an ordered genome is a *list of gene tokens*, a nucleotide
genome is a **karyotype of chromosomes**, each a coordinate axis of base pairs represented as an
ordered list of **blocks**. The vocabulary (fixed here, once):

- **block** — the persistent unit of both the representation *and* the ancestry: a **maximal** run of
  one unbroken ancestry, the half-open interval ``[start, end)`` on an ancestral ``source``, read
  forward (``strand`` ``+1``) or reverse-complemented (``-1``). *Maximal*: two adjacent collinear
  blocks are one lineage, so they are merged. A block carries (later) one gene tree.
- **gene / intergene** — a *classification* of blocks (genic mode, later): a gene is a declared,
  indivisible block; an intergene fragments freely into many blocks.
- **segment** — *not* an object: the **extent** a single event acts on, the arc ``[start,
  start+length)``. "An inversion inverts a segment of a chromosome."

This module is grown **slice by slice** (the legacy engine is ~1900 lines, so we build it very
slowly). Built so far:

- **The representation + inversion** — split a chromosome at a coordinate, invert an arc (which may
  wrap the origin on a circular chromosome), keep blocks maximal, and trace every nucleotide back to
  its ancestral origin.
- **The tree wiring + the multi-chromosome container** (``simulate_genomes_nucleotide``): a genome is
  a karyotype of chromosomes (heterogeneous **sizes and shapes**), each **identity-bearing** (a
  chromosome id re-minted at every speciation, the edge recorded in the shared chromosome network);
  inversions act within a length-weighted chromosome; the whole karyotype is inherited at speciation.
- **The number-changing chromosome tier**: ``fission`` (a bifurcation — one chromosome into two) and
  ``fusion`` (the reticulation — two same-topology chromosomes into one), both re-minting their
  children into the chromosome network. Every event conserves total length but fission/fusion change
  the chromosome *count*, so a branch is a small Gillespie.
- **Translocation** — an arc moved from one chromosome to a different one (optionally inverted, with
  probability ``inversion_probability``). Both chromosomes persist, so it is a *rearrangement*, not a
  network edge.

Those are all **ancestry-neutral** (the strong invariant: every node carries the whole root sequence,
permuted). The first ancestry-**changing** event is now here too:

- **Loss** — an arc deleted (per nucleotide, never emptying a chromosome). A death, in the ``events``
  (genealogy) log. On its own it weakens the invariant to *subset*: each ancestral position at most
  once, monotonically down every path.
- **Duplication** — an arc copied in **tandem** (per nucleotide). The first *birth*: the copied
  material now has an extra copy (same source coordinates), so a position can appear more than once.

Threaded through all of the above is the **copy lineage** (``Block.copy``): the persistent identity
that a split preserves, a duplication mints fresh (parent → child), and a speciation re-mints into
each daughter. The genealogy ``events`` — ``origination`` / ``loss`` / ``duplication`` / ``speciation``
— carry these copy ids, and the last piece here reads them:

- **The gene-tree recovery** (``result.root_blocks`` / ``result.gene_trees``) — the **root partition**
  (cut each source at the union of extant-leaf breakpoints, keep the intervals that survive in a leaf:
  each is a maximal never-cut interval = one gene family) and, per block, the copy-lineage log replayed
  into one gene tree (a duplication is a ladder, a speciation a bifurcation), reusing the shared
  per-segment builder. Validated end to end: the recovered extant tips equal the copies actually
  present in every extant leaf.

Deferred to later slices: transfer (→ the global timeline), transposition, origination, and chromosome
origination / loss; indels; declared genes / intergenes (pseudogenization, replacement); GFF input and
BED / FASTA output. (The recovery above builds the *extant* tree exactly; the *complete* tree's dead
partial-overlap lineages need the finer all-node partition — a later refinement.)
"""

from __future__ import annotations

import collections
from dataclasses import dataclass

import numpy as np

from ..species import SpeciesResult, Tree
from .chromosomes import ChromosomeEvent
from .gene_trees import GeneTree, gene_trees_from_events


@dataclass
class Block:
    """A run of one unbroken ancestry: the half-open interval ``[start, end)`` on ``source`` this run
    descends from (``start < end`` always), read forward (``strand`` ``+1``) or reverse-complemented
    (``-1``). Its physical length is ``end - start``. ``copy`` names the **copy lineage** it belongs to
    — the persistent identity that threads through splits and moves (a split preserves it, a
    duplication mints a fresh child), re-minted into daughters at speciation; the gene-tree recovery
    reads it to resolve *which copy begat which*. (A gene / intergene classification joins later.)

    Blocks are **not** kept maximal during the run — a rearrangement leaves collinear neighbours split
    rather than merging them (merging would muddy the copy-lineage bookkeeping). Maximality is instead
    a property of the *recovered* root-blocks: the coarsest partition that survives in the leaves."""

    source: int
    start: int
    end: int
    strand: int
    copy: int = 0

    @property
    def length(self) -> int:
        return self.end - self.start


def _split_block(b: Block, o: int) -> list[Block]:
    """Split ``b`` after ``o`` physical positions into ``[left, right]``, strand-aware. For a forward
    block the cut falls at ``start + o``; for a reversed one the first ``o`` physical positions are
    the **high** end of the source, so the cut falls at ``end - o``. Both pieces keep ``b``'s copy
    lineage — a split is not a birth."""
    if b.strand == 1:
        return [Block(b.source, b.start, b.start + o, 1, b.copy),
                Block(b.source, b.start + o, b.end, 1, b.copy)]
    return [Block(b.source, b.end - o, b.end, -1, b.copy),
            Block(b.source, b.start, b.end - o, -1, b.copy)]


@dataclass
class Chromosome:
    """One replicon: an ordered list of :class:`Block`\\ s over a nucleotide coordinate axis, with an
    ``id`` and a ``topology`` — ``"circular"`` (a ring, where coordinate ``length`` wraps to ``0`` and
    an arc may cross the origin) or ``"linear"`` (two ends, no wrap). Owns the per-chromosome
    operations. Blocks are **not** merged after an event (see :class:`Block`): an event only ever
    *splits*, so its endpoints stay as breakpoints and the copy lineages are never fused."""

    id: int
    topology: str
    blocks: list[Block]

    def __post_init__(self) -> None:
        if self.topology not in ("circular", "linear"):
            raise ValueError(f"topology must be 'circular' or 'linear', got {self.topology!r}")

    @property
    def length(self) -> int:
        return sum(b.length for b in self.blocks)

    def _split_at(self, c: int) -> None:
        """Ensure a block boundary at physical coordinate ``c`` (``0 <= c <= length``). A no-op at the
        ends or an existing boundary; otherwise the block straddling ``c`` is split."""
        if c <= 0:
            return
        pos = 0
        for i, b in enumerate(self.blocks):
            if pos == c:
                return
            if pos < c < pos + b.length:
                self.blocks[i:i + 1] = _split_block(b, c - pos)
                return
            pos += b.length

    def _index_at(self, phys: int) -> int:
        pos = 0
        for i, b in enumerate(self.blocks):
            if pos == phys:
                return i
            pos += b.length
        if pos == phys:
            return len(self.blocks)
        raise ValueError(f"no block boundary at physical {phys}")

    def _arc_range(self, start: int, length: int) -> tuple[int, int] | None:
        """Prepare the arc ``[start, start+length)`` and return the block index range ``[i, j)`` it
        occupies, or ``None`` for an empty arc. Splits at both ends first. A **linear** chromosome
        clamps the arc to its ends; a **circular** one may wrap the origin, rotating the ring to bring
        the arc to the front (the origin drifts to a real breakpoint — harmless on a ring)."""
        total = self.length
        if total == 0:
            return None
        if self.topology == "linear":
            start = max(0, min(start, total))
            end = min(start + max(0, length), total)
            if end <= start:
                return None
            self._split_at(start)
            self._split_at(end)
            return self._index_at(start), self._index_at(end)
        ell = max(1, min(length, total))
        s = start % total
        self._split_at(s)
        self._split_at((s + ell) % total)
        if s + ell <= total:
            j = self._index_at(s + ell) if s + ell < total else len(self.blocks)
            return self._index_at(s), j
        i = self._index_at(s)
        self.blocks[:] = self.blocks[i:] + self.blocks[:i]
        j, acc = 0, 0
        while acc < ell:
            acc += self.blocks[j].length
            j += 1
        return 0, j

    def invert(self, start: int, length: int) -> None:
        """Invert the arc ``[start, start+length)``: reverse the block order and flip each strand.
        Ancestry (and every copy lineage) is unchanged; the arc endpoints stay as breakpoints."""
        span = self._arc_range(start, length)
        if span is None:
            return
        i, j = span
        arc = self.blocks[i:j]
        arc.reverse()
        for b in arc:
            b.strand = -b.strand
        self.blocks[i:j] = arc

    def mosaic(self) -> list[tuple[int, int, int, int]]:
        """The chromosome as ordered blocks — ``[(source, start, end, strand), ...]`` in physical
        order (already maximal)."""
        return [(b.source, b.start, b.end, b.strand) for b in self.blocks]

    def trace_back(self) -> list[tuple[int, int, int]]:
        """Every nucleotide's ancestral origin, left to right: ``(source, source_position, strand)``.
        A forward block reads its source low→high, a reversed one high→low."""
        out: list[tuple[int, int, int]] = []
        for b in self.blocks:
            if b.strand == 1:
                out.extend((b.source, p, 1) for p in range(b.start, b.end))
            else:
                out.extend((b.source, p, -1) for p in range(b.end - 1, b.start - 1, -1))
        return out


@dataclass
class NucleotideGenome:
    """A **karyotype**: an ordered list of :class:`Chromosome`\\ s. Its :attr:`length` is the total
    over all chromosomes; a length-scaled event lands on a chromosome in proportion to its bp."""

    chromosomes: list[Chromosome]

    @property
    def length(self) -> int:
        return sum(c.length for c in self.chromosomes)

    def _pick_position(self, rng) -> tuple[Chromosome, int]:
        """A uniform nucleotide pick → ``(chromosome, physical position)`` — realises a per-nucleotide
        (length-weighted) choice of chromosome."""
        m = int(rng.integers(self.length))
        for c in self.chromosomes:
            if m < c.length:
                return c, m
            m -= c.length
        raise AssertionError("length out of sync with the chromosomes")  # unreachable

    def mosaic(self) -> dict[int, list[tuple[int, int, int, int]]]:
        """``{chromosome id: its block mosaic}``."""
        return {c.id: c.mosaic() for c in self.chromosomes}

    def trace_back(self) -> dict[int, list[tuple[int, int, int]]]:
        """``{chromosome id: its per-nucleotide trace-back}``."""
        return {c.id: c.trace_back() for c in self.chromosomes}

    def ancestry(self) -> list[tuple[int, int]]:
        """The sorted multiset of ancestral ``(source, position)`` over the whole genome — orientation
        and chromosome agnostic. Conserved by inversion, so it is the strong invariant."""
        return sorted((src, pos) for c in self.chromosomes for (src, pos, _s) in c.trace_back())


# --- the tree wiring (inversions along the species tree; a multi-chromosome, identity-bearing karyotype)

@dataclass(frozen=True)
class Inversion:
    """A recorded nucleotide inversion: on species branch ``lineage`` at ``time``, the arc
    ``[start, start+length)`` of chromosome ``chromosome`` was reversed. Ancestry is unchanged, so it
    is a rearrangement record, not a gene-genealogy event."""

    time: float
    lineage: int
    chromosome: int
    start: int
    length: int


@dataclass(frozen=True)
class Translocation:
    """A recorded nucleotide translocation: on branch ``lineage`` at ``time``, the arc
    ``[start, start+length)`` of chromosome ``source`` was moved to chromosome ``dest`` (a different
    chromosome of the same genome), landing ``flipped`` (reversed) or not. Ancestry is unchanged — the
    blocks keep their source coordinates — and both chromosomes persist, so it is a rearrangement, not
    a chromosome-network edge."""

    time: float
    lineage: int
    source: int
    dest: int
    start: int
    length: int
    flipped: bool


@dataclass(frozen=True)
class Origination:
    """A recorded **birth of a seed copy lineage** — the root of a gene tree. At ``time`` on branch
    ``lineage`` a seed replicon was laid down on chromosome ``chromosome`` as copy lineage ``copy``,
    covering the ancestral interval ``[start, end)`` on ``source``. One per seed replicon; the
    gene-tree recovery reads these as the roots."""

    time: float
    lineage: int
    chromosome: int
    copy: int
    source: int
    start: int
    end: int


@dataclass(frozen=True)
class Loss:
    """A recorded nucleotide loss — an ancestry-**changing** event (a death), so it belongs to the
    genealogy log, not the rearrangements. On branch ``lineage`` at ``time`` an arc of chromosome
    ``chromosome`` was deleted; ``lost`` names the material that died on this lineage as
    ``(copy, source, start, end)`` rows — copy lineage ``copy`` lost ``[start, end)`` of ``source``.
    The per-block gene trees read these deaths."""

    time: float
    lineage: int
    chromosome: int
    lost: tuple[tuple[int, int, int, int], ...]


@dataclass(frozen=True)
class Duplication:
    """A recorded nucleotide duplication — an ancestry-**changing** event (a birth). On branch
    ``lineage`` at ``time`` an arc of chromosome ``chromosome`` was copied in **tandem**; ``copied``
    names the parentage as ``(parent_copy, child_copy, source, start, end)`` rows — over ``[start,
    end)`` of ``source``, copy lineage ``parent_copy`` begat the fresh ``child_copy``. It is a branch
    point in a block's gene tree; the copy ids resolve *which* copy begat which."""

    time: float
    lineage: int
    chromosome: int
    copied: tuple[tuple[int, int, int, int, int], ...]


@dataclass(frozen=True)
class Speciation:
    """A recorded **copy-lineage speciation** — parallel to the chromosome network's speciation edge
    but at gene-tree granularity. At ``time`` on branch ``lineage`` (the parent node) copy lineage
    ``parent`` was re-minted into one fresh ``children`` copy per daughter species (in the daughters'
    order). Chromosome-agnostic: a copy may span chromosomes, and the whole lineage is inherited."""

    time: float
    lineage: int
    parent: int
    children: tuple[int, ...]


@dataclass
class NucleotideGenomesResult:
    """What :func:`simulate_genomes_nucleotide` returns: the ``complete_tree`` it ran on, the final
    nucleotide ``genomes`` (karyotypes) at **every** node, the **copy-lineage genealogy** ``events``
    (``origination``, ``loss``, ``duplication``, ``speciation`` — carrying the copy ids the gene-tree
    recovery reads), the ancestry-neutral ``rearrangements`` (inversion, translocation), the
    ``chromosome_events`` (the chromosome network), and the ``seed``. ``mosaic`` / ``trace_back`` /
    ``ancestry`` read a node's genome."""

    complete_tree: Tree
    genomes: dict[int, NucleotideGenome]
    events: list[Origination | Loss | Duplication | Speciation]
    rearrangements: list[Inversion | Translocation]
    chromosome_events: list[ChromosomeEvent]
    seed: int | None

    def mosaic(self, node_id: int) -> dict[int, list[tuple[int, int, int, int]]]:
        return self.genomes[node_id].mosaic()

    def trace_back(self, node_id: int) -> dict[int, list[tuple[int, int, int]]]:
        return self.genomes[node_id].trace_back()

    def ancestry(self, node_id: int) -> list[tuple[int, int]]:
        return self.genomes[node_id].ancestry()

    def _recover(self) -> tuple[list[tuple[int, int, int]], dict[int, GeneTree]]:
        if not hasattr(self, "_recovered"):
            self._recovered = _recover_gene_trees(self)
        return self._recovered

    @property
    def root_blocks(self) -> list[tuple[int, int, int]]:
        """The recovered **root partition**: ``(source, start, end)`` for each maximal never-cut
        interval that survives in an extant leaf — one per :attr:`gene_trees` family (by index)."""
        return self._recover()[0]

    @property
    def gene_trees(self) -> dict[int, GeneTree]:
        """``{family: GeneTree}`` — one gene tree per root-block (family id = the block's index in
        :attr:`root_blocks`), recovered from the copy-lineage genealogy."""
        return self._recover()[1]


def _valid_length(length) -> int:
    if isinstance(length, bool) or not isinstance(length, int) or length < 1:
        raise ValueError(f"a chromosome length must be a positive integer, got {length!r}")
    return length


def _replicon_specs(chromosomes, root_length, topology) -> list[tuple[int, str]]:
    """Resolve the ``chromosomes`` argument to a list of ``(length, topology)`` replicon specs. An int
    ``N`` seeds ``N`` equal replicons of ``root_length`` and ``topology``; a list gives heterogeneous
    replicons of different **sizes and shapes**, e.g. ``[(1000, "circular"), (50, "linear")]``."""
    if isinstance(chromosomes, bool) or isinstance(chromosomes, int):
        if isinstance(chromosomes, bool) or chromosomes < 1:
            raise ValueError(f"chromosomes must be a positive integer or a list of specs, got {chromosomes!r}")
        return [(_valid_length(root_length), topology)] * chromosomes
    specs = [(_valid_length(length), top) for (length, top) in chromosomes]
    if not specs:
        raise ValueError("chromosomes must have at least one replicon")
    return specs


def _copy_chromosome(c: Chromosome, cid: int, copy_map: dict[int, int]) -> Chromosome:
    """A daughter chromosome: a fresh id, a deep copy of the blocks (fresh :class:`Block` objects, so
    a daughter's inversions never mutate the parent's genome), with every block's copy lineage
    re-minted through ``copy_map`` (parent copy id → this daughter's fresh copy id)."""
    return Chromosome(cid, c.topology,
                      [Block(b.source, b.start, b.end, b.strand, copy_map[b.copy]) for b in c.blocks])


@dataclass(frozen=True)
class _Rates:
    inversion: float
    translocation: float
    loss: float
    duplication: float
    fission: float
    fusion: float
    inversion_length: float
    translocation_length: float
    loss_length: float
    duplication_length: float
    inversion_probability: float


def _do_duplication(g, node_id, t, duplication_length, rng, events, new_copy) -> None:
    """Copy a geometric-length arc of a length-weighted chromosome **in tandem** (the copy inserted
    right after the arc). An ancestry-changing *birth*: the copied material now has an extra copy (same
    source coordinates). Each distinct copy lineage in the arc begets one fresh child lineage, so the
    tandem copy is a new copy of that material; the parentage is recorded as a :class:`Duplication`."""
    chrom, start = g._pick_position(rng)
    ell = min(chrom.length, max(1, int(rng.geometric(1.0 / duplication_length))))
    span = chrom._arc_range(start, ell)
    if span is None:
        return
    i, j = span
    arc = chrom.blocks[i:j]
    child_of: dict[int, int] = {}
    for b in arc:
        if b.copy not in child_of:
            child_of[b.copy] = new_copy()
    copied = tuple((b.copy, child_of[b.copy], b.source, b.start, b.end) for b in arc)
    chrom.blocks[j:j] = [Block(b.source, b.start, b.end, b.strand, child_of[b.copy]) for b in arc]
    events.append(Duplication(t, node_id, chrom.id, copied))


def _do_loss(g, node_id, t, loss_length, rng, events) -> None:
    """Delete a geometric-length arc from a length-weighted chromosome — an ancestry-changing event (a
    death). Never empties a chromosome (leaves at least one nucleotide; whole-chromosome loss is a
    deferred tier event). Records the deleted material — which copy lineage lost which arc — as a
    :class:`Loss`."""
    chrom, start = g._pick_position(rng)
    if chrom.length < 2:
        return
    ell = min(chrom.length - 1, max(1, int(rng.geometric(1.0 / loss_length))))
    span = chrom._arc_range(start, ell)
    if span is None:
        return
    i, j = span
    lost = tuple((b.copy, b.source, b.start, b.end) for b in chrom.blocks[i:j])
    chrom.blocks = chrom.blocks[:i] + chrom.blocks[j:]
    events.append(Loss(t, node_id, chrom.id, lost))


def _do_inversion(g, node_id, t, inversion_length, rng, rearrangements) -> None:
    """Invert a geometric-length arc of a length-weighted chromosome."""
    chrom, pos = g._pick_position(rng)
    length = min(chrom.length, int(rng.geometric(1.0 / inversion_length)))
    chrom.invert(pos, length)
    rearrangements.append(Inversion(t, node_id, chrom.id, pos, length))


def _do_translocation(g, node_id, t, translocation_length, inversion_probability, rng,
                      rearrangements) -> None:
    """Move a geometric-length arc from a length-weighted source chromosome to a uniformly-chosen
    **different** chromosome, landing inverted with probability ``inversion_probability``.
    Ancestry-neutral (blocks keep their source coordinates); both chromosomes persist. No-op if the
    genome has one chromosome or the source is below two nucleotides (an arc never empties a
    chromosome — it leaves at least one)."""
    if len(g.chromosomes) < 2:
        return
    source, start = g._pick_position(rng)
    if source.length < 2:
        return
    ell = min(source.length - 1, max(1, int(rng.geometric(1.0 / translocation_length))))
    span = source._arc_range(start, ell)
    if span is None:
        return
    i, j = span
    arc = source.blocks[i:j]
    source.blocks = source.blocks[:i] + source.blocks[j:]
    flipped = bool(rng.random() < inversion_probability)
    if flipped:
        arc = [Block(b.source, b.start, b.end, -b.strand, b.copy) for b in reversed(arc)]
    others = [c for c in g.chromosomes if c is not source]
    dest = others[int(rng.integers(len(others)))]
    pos = int(rng.integers(dest.length + 1))
    dest._split_at(pos)
    k = dest._index_at(pos)
    dest.blocks[k:k] = arc
    rearrangements.append(Translocation(t, node_id, source.id, dest.id, start, ell, flipped))


def _do_fission(g, node_id, t, rng, chromosome_events, new_chrom_id) -> None:
    """Split a uniformly-chosen chromosome into two re-minted children — a **bifurcation**. A linear
    chromosome cuts at one breakpoint (prefix + suffix); a circular one at two (the arc between them
    a new ring, the remainder another). Blocks keep their ancestry; no-op below two nucleotides."""
    ci = int(rng.integers(len(g.chromosomes)))
    chrom = g.chromosomes[ci]
    if chrom.length < 2:
        return
    if chrom.topology == "linear":
        c = int(rng.integers(1, chrom.length))
        chrom._split_at(c)
        i = chrom._index_at(c)
        a = Chromosome(new_chrom_id(), "linear", chrom.blocks[:i])
        b = Chromosome(new_chrom_id(), "linear", chrom.blocks[i:])
    else:
        c1, c2 = sorted(int(x) for x in rng.choice(chrom.length, size=2, replace=False))
        chrom._split_at(c1)
        chrom._split_at(c2)
        i, j = chrom._index_at(c1), chrom._index_at(c2)
        a = Chromosome(new_chrom_id(), "circular", chrom.blocks[i:j])
        b = Chromosome(new_chrom_id(), "circular", chrom.blocks[:i] + chrom.blocks[j:])
    g.chromosomes[ci:ci + 1] = [a, b]
    chromosome_events.append(ChromosomeEvent(t, "fission", node_id, (chrom.id,), (a.id, b.id)))


def _do_fusion(g, node_id, t, rng, chromosome_events, new_chrom_id) -> None:
    """Merge a uniformly-chosen chromosome with another of the **same topology** — a **reticulation**
    (two parents, one child): concatenate their blocks, re-mint the child. No-op if the genome has no
    same-topology partner for the chosen chromosome."""
    ci = int(rng.integers(len(g.chromosomes)))
    a = g.chromosomes[ci]
    partners = [k for k in range(len(g.chromosomes))
                if k != ci and g.chromosomes[k].topology == a.topology]
    if not partners:
        return
    cj = partners[int(rng.integers(len(partners)))]
    b = g.chromosomes[cj]
    fused = Chromosome(new_chrom_id(), a.topology, a.blocks + b.blocks)
    g.chromosomes[:] = [c for k, c in enumerate(g.chromosomes) if k not in (ci, cj)] + [fused]
    chromosome_events.append(ChromosomeEvent(t, "fusion", node_id, (a.id, b.id), (fused.id,)))


def _evolve_branch(g, node_id, t0, t1, rates, rng, events, rearrangements, chromosome_events,
                   new_chrom_id, new_copy) -> None:
    """A Gillespie over the branch ``[t0, t1]``. The per-nucleotide rates (inversion, translocation,
    loss) scale with the current length, and fission/fusion with the current chromosome count — both
    of which change as events fire — so the loop recomputes the rates after each. All within-lineage
    (no coupling yet); the global timeline waits for transfer."""
    if t1 <= t0:
        return
    t = t0
    while True:
        nc = len(g.chromosomes)
        length = g.length
        r_inv = rates.inversion * length
        r_trl = rates.translocation * length if nc >= 2 else 0.0
        r_los = rates.loss * length
        r_dup = rates.duplication * length
        r_fis = rates.fission * nc
        r_fus = rates.fusion * nc if nc >= 2 else 0.0
        total = r_inv + r_trl + r_los + r_dup + r_fis + r_fus
        if total <= 0:
            return
        t += float(rng.exponential(1.0 / total))
        if t >= t1:
            return
        r = float(rng.random()) * total
        if r < r_inv:
            _do_inversion(g, node_id, t, rates.inversion_length, rng, rearrangements)
        elif r < r_inv + r_trl:
            _do_translocation(g, node_id, t, rates.translocation_length, rates.inversion_probability,
                              rng, rearrangements)
        elif r < r_inv + r_trl + r_los:
            _do_loss(g, node_id, t, rates.loss_length, rng, events)
        elif r < r_inv + r_trl + r_los + r_dup:
            _do_duplication(g, node_id, t, rates.duplication_length, rng, events, new_copy)
        elif r < r_inv + r_trl + r_los + r_dup + r_fis:
            _do_fission(g, node_id, t, rng, chromosome_events, new_chrom_id)
        else:
            _do_fusion(g, node_id, t, rng, chromosome_events, new_chrom_id)


def simulate_genomes_nucleotide(tree, *, inversion=0.0, inversion_length=50.0, translocation=0.0,
                                translocation_length=50.0, inversion_probability=0.0, loss=0.0,
                                loss_length=50.0, duplication=0.0, duplication_length=50.0, fission=0.0,
                                fusion=0.0, chromosomes=1, root_length=1000, topology="circular",
                                seed=None) -> NucleotideGenomesResult:
    """Evolve a nucleotide genome along a species tree by inversion, translocation, **loss**, and the
    number-changing chromosome tier. The root is seeded with a **karyotype** — ``chromosomes``
    replicons, each its own source: an int ``N`` gives ``N`` equal replicons of
    ``root_length``/``topology``, or pass a list of ``(length, topology)`` for heterogeneous **sizes
    and shapes**. Each lineage inherits a copy of its parent's karyotype at speciation, with **every
    chromosome re-minted** (recorded in ``chromosome_events``, the chromosome network), then evolves
    along its branch:

    - ``inversion`` (**per nucleotide**) reverses a geometric-length (mean ``inversion_length``) arc
      of a length-weighted chromosome.
    - ``translocation`` (**per nucleotide**) moves a geometric-length (mean ``translocation_length``)
      arc to a different chromosome, landing inverted with probability ``inversion_probability``. Both
      chromosomes persist — it is a rearrangement, not a network edge.
    - ``loss`` (**per nucleotide**) deletes a geometric-length (mean ``loss_length``) arc — the first
      ancestry-**changing** event (a death), recorded in ``events``. Never empties a chromosome.
    - ``fission`` (**per chromosome**) splits a chromosome in two (a **bifurcation**); ``fusion``
      (**per chromosome**) merges two chromosomes of the same topology (the **reticulation**). Both
      re-mint their children and record a network edge.

    The branch is a small Gillespie recomputing its rates as the length (which loss shrinks) and the
    chromosome count (which fission/fusion change) evolve; all within-lineage (the global timeline
    waits for transfer). With loss, the strong invariant weakens: every node carries a **subset** of
    the root sequence (each ancestral position at most once, and monotonically down every path).
    Deterministic given ``seed``."""
    tree = tree.complete_tree if isinstance(tree, SpeciesResult) else tree
    for label, rate in (("inversion", inversion), ("translocation", translocation), ("loss", loss),
                        ("duplication", duplication), ("fission", fission), ("fusion", fusion)):
        if rate < 0:
            raise ValueError(f"{label} must be >= 0, got {rate}")
    for label, mean in (("inversion_length", inversion_length),
                        ("translocation_length", translocation_length), ("loss_length", loss_length),
                        ("duplication_length", duplication_length)):
        if mean <= 0:
            raise ValueError(f"{label} must be > 0, got {mean}")
    if not 0.0 <= inversion_probability <= 1.0:
        raise ValueError(f"inversion_probability must be in [0, 1], got {inversion_probability}")
    specs = _replicon_specs(chromosomes, root_length, topology)
    rates = _Rates(inversion, translocation, loss, duplication, fission, fusion, inversion_length,
                   translocation_length, loss_length, duplication_length, inversion_probability)

    rng = np.random.default_rng(seed)
    chrom_counter = 0
    copy_counter = 0

    def new_chrom_id() -> int:
        nonlocal chrom_counter
        cid = chrom_counter
        chrom_counter += 1
        return cid

    def new_copy() -> int:
        nonlocal copy_counter
        copy_counter += 1
        return copy_counter                             # copy ids start at 1 (0 = the unset sentinel)

    genomes: dict[int, NucleotideGenome] = {}
    events: list[Origination | Loss | Duplication | Speciation] = []
    rearrangements: list[Inversion | Translocation] = []
    chromosome_events: list[ChromosomeEvent] = []
    root = tree.nodes[tree.root]

    root_chroms = []
    for source, (length, top) in enumerate(specs):     # one source per seed replicon; each a network root
        cid = new_chrom_id()
        cp = new_copy()                                 # ...and one seed copy lineage per replicon
        root_chroms.append(Chromosome(cid, top, [Block(source, 0, length, 1, cp)]))
        chromosome_events.append(ChromosomeEvent(root.birth_time, "origination", root.id, (), (cid,)))
        events.append(Origination(root.birth_time, root.id, cid, cp, source, 0, length))
    start_genome = {root.id: NucleotideGenome(root_chroms)}

    stack = [root.id]
    while stack:  # DFS from the root: a parent's final genome exists before its children inherit it
        node_id = stack.pop()
        node = tree.nodes[node_id]
        g = start_genome.pop(node_id)
        _evolve_branch(g, node_id, node.birth_time, node.end_time, rates, rng, events,
                       rearrangements, chromosome_events, new_chrom_id, new_copy)
        genomes[node_id] = g
        if node.children is not None:  # speciation: re-mint every chromosome AND copy lineage per daughter
            copy_maps: dict[int, dict[int, int]] = {c: {} for c in node.children}
            for pc in sorted({b.copy for chrom in g.chromosomes for b in chrom.blocks}):
                dcs = []
                for c in node.children:
                    dc = new_copy()
                    copy_maps[c][pc] = dc
                    dcs.append(dc)
                events.append(Speciation(node.end_time, node.id, pc, tuple(dcs)))
            starts: dict[int, list[Chromosome]] = {c: [] for c in node.children}
            for pchrom in g.chromosomes:
                dcids = []
                for c in node.children:
                    dcid = new_chrom_id()
                    dcids.append(dcid)
                    starts[c].append(_copy_chromosome(pchrom, dcid, copy_maps[c]))
                chromosome_events.append(
                    ChromosomeEvent(node.end_time, "speciation", node.id, (pchrom.id,), tuple(dcids)))
            for c in node.children:
                start_genome[c] = NucleotideGenome(starts[c])
                stack.append(c)
    return NucleotideGenomesResult(tree, genomes, events, rearrangements, chromosome_events, seed)


# --- the gene-tree recovery: root partition -> per-block genealogy -> one tree per block ----------
#
# A gene tree per *root-block*: the coarsest interval of a source that is never cut in any extant leaf
# and survives there. Within such a block every copy is un-cut in every leaf that carries it (a cut
# would be an observable breakpoint, which would have split the block), so all its copies share one
# genealogy. The recovery has two moves:
#
#   1. the **root partition** — cut each source at the union of extant-leaf breakpoints and keep the
#      intervals some extant leaf still covers (the surviving material);
#   2. per block, replay the copy-lineage log restricted to that block into the **per-segment** model
#      the shared ``gene_trees_from_events`` reads (every event ends a segment and starts fresh ids),
#      so a duplication is a ladder (parent continues + child) and a speciation a bifurcation.
#
# Because the partition is by *surviving* breakpoints, every event that reaches an extant copy is
# atomic on it (covers the whole block or none) — an event that only partially covers a block made a
# breakpoint that did not survive, so its child is dead and irrelevant to the extant tree.


@dataclass(frozen=True)
class _SegEvent:
    """One per-segment event in the model :func:`gene_trees_from_events` reads: ``kind`` is
    ``origination`` (a root) / ``duplication`` / ``speciation`` (a parent segment ends → a child
    begins on species branch ``lineage``) / ``loss`` (a dead leaf). ``copy`` is the fresh segment id,
    ``parent`` the segment it descends from (``None`` for a root or a loss)."""

    kind: str
    family: int
    lineage: int
    time: float
    copy: int
    parent: int | None


def _root_block_partition(result) -> list[tuple[int, int, int]]:
    """The root partition: per source, the maximal intervals bounded by extant-leaf breakpoints that
    survive in some extant leaf. Returns a sorted list of ``(source, start, end)`` root-blocks."""
    leaves = [n.id for n in result.complete_tree.extant()]
    bounds: dict[int, set[int]] = collections.defaultdict(set)
    spans: dict[int, set[tuple[int, int]]] = collections.defaultdict(set)
    for lid in leaves:
        for chrom in result.genomes[lid].chromosomes:
            for b in chrom.blocks:
                bounds[b.source].update((b.start, b.end))
                spans[b.source].add((b.start, b.end))
    blocks: list[tuple[int, int, int]] = []
    for source in bounds:
        cuts = sorted(bounds[source])
        covers = spans[source]
        for a, c in zip(cuts, cuts[1:]):
            if any(x <= a and c <= y for (x, y) in covers):   # some leaf still carries [a, c)
                blocks.append((source, a, c))
    return sorted(blocks)


def _emit_block_events(fam, s, a, b, tree, origs, dups, losses, specs, new_seg, out) -> None:
    """Replay the copy-lineage log for one root-block ``(s, [a, b))`` into per-segment events on
    ``out``. A copy lineage is a *block-copy* when it covers ``[a, b)`` in full; a duplication that
    covers it in full begets a block-copy child (a ladder rung on its parent), a speciation re-mints
    it, and any loss overlapping it ends it."""
    def covers(x, y):
        return x <= a and b <= y

    def overlaps(x, y):
        return x < b and a < y

    root = tree.nodes[tree.root]
    root_copies = [e.copy for e in origs if e.source == s and covers(e.start, e.end)]
    if not root_copies:
        return

    dup_child: dict[int, list[tuple[float, int]]] = collections.defaultdict(list)   # parent -> [(t, child)]
    species: dict[int, int] = {}                                                    # copy -> species branch
    for e in dups:
        for (pc, cc, src, x, y) in e.copied:
            if src == s and covers(x, y):
                dup_child[pc].append((e.time, cc))
                species[cc] = e.lineage
    loss_of: dict[int, float] = {}                                                  # copy -> earliest loss time
    for e in losses:
        for (cp, src, x, y) in e.lost:
            if src == s and overlaps(x, y) and (cp not in loss_of or e.time < loss_of[cp]):
                loss_of[cp] = e.time
    for rc in root_copies:
        species[rc] = root.id

    block_copies: set[int] = set()
    order: list[int] = []
    stack = list(root_copies)
    while stack:                                            # BFS the block-copy forest (single-parent)
        c = stack.pop()
        if c in block_copies:
            continue
        block_copies.add(c)
        order.append(c)
        for (_t, cc) in dup_child.get(c, ()):
            stack.append(cc)
        if c not in loss_of and c in specs:                # survived to its node's end -> re-minted
            pnode = tree.nodes[specs[c].lineage]
            for i, d in enumerate(specs[c].children):
                species[d] = pnode.children[i]
                stack.append(d)

    seg_in = {c: new_seg() for c in order}
    for rc in root_copies:
        out.append(_SegEvent("origination", fam, root.id, root.birth_time, seg_in[rc], None))
    for c in order:
        prev = seg_in[c]
        for (t, cc) in sorted(dup_child.get(c, ())):       # ladder: each duplication rung is a bifurcation
            nxt = new_seg()
            out.append(_SegEvent("duplication", fam, species[c], t, nxt, prev))
            out.append(_SegEvent("duplication", fam, species[cc], t, seg_in[cc], prev))
            prev = nxt
        if c in loss_of:                                   # a death (dead leaf)
            out.append(_SegEvent("loss", fam, species[c], loss_of[c], prev, None))
        elif c in specs:                                   # a bifurcation into the daughter species
            pnode = tree.nodes[specs[c].lineage]
            for i, d in enumerate(specs[c].children):
                if d in block_copies:
                    out.append(_SegEvent("speciation", fam, pnode.children[i], specs[c].time,
                                         seg_in[d], prev))
        # else: prev survives to an extant/extinct leaf — gene_trees_from_events tags it by species fate


def _recover_gene_trees(result) -> tuple[list[tuple[int, int, int]], dict[int, GeneTree]]:
    """The full recovery: the root partition, and one :class:`GeneTree` per block (its family id is the
    block's index in the partition). Reuses the shared per-segment tree builder."""
    tree = result.complete_tree
    blocks = _root_block_partition(result)
    origs = [e for e in result.events if isinstance(e, Origination)]
    dups = [e for e in result.events if isinstance(e, Duplication)]
    losses = [e for e in result.events if isinstance(e, Loss)]
    specs = {e.parent: e for e in result.events if isinstance(e, Speciation)}
    counter = [0]

    def new_seg():
        counter[0] += 1
        return counter[0]

    seg_events: list[_SegEvent] = []
    for fam, (s, a, b) in enumerate(blocks):
        _emit_block_events(fam, s, a, b, tree, origs, dups, losses, specs, new_seg, seg_events)
    return blocks, gene_trees_from_events(seg_events, tree)


__all__ = ["Block", "Chromosome", "NucleotideGenome", "NucleotideGenomesResult",
           "Origination", "Loss", "Duplication", "Speciation", "Inversion", "Translocation",
           "simulate_genomes_nucleotide"]
