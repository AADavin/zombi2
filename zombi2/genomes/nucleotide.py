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
  Because inversion conserves ancestry, every node still carries the whole root sequence, merely
  permuted across its chromosomes — the strong invariant.

Deferred to later slices: loss / duplication / transfer / transposition / origination (the events that
birth and kill ancestry, creating per-block gene trees), the chromosome **tier** (fission, fusion,
translocation — the chromosomes *moving*) and chromosome origination / loss, indels, declared
genes / intergenes (pseudogenization, replacement), GFF input and BED / FASTA output.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..species import SpeciesResult, Tree
from .chromosomes import ChromosomeEvent


@dataclass
class Block:
    """A **maximal run of one unbroken ancestry**: the half-open interval ``[start, end)`` on
    ``source`` this run descends from (``start < end`` always), read forward (``strand`` ``+1``) or
    reverse-complemented (``-1``). Its physical length is ``end - start``. (A gene / intergene
    classification joins later.)"""

    source: int
    start: int
    end: int
    strand: int

    @property
    def length(self) -> int:
        return self.end - self.start


def _split_block(b: Block, o: int) -> list[Block]:
    """Split ``b`` after ``o`` physical positions into ``[left, right]``, strand-aware. For a forward
    block the cut falls at ``start + o``; for a reversed one the first ``o`` physical positions are
    the **high** end of the source, so the cut falls at ``end - o``."""
    if b.strand == 1:
        return [Block(b.source, b.start, b.start + o, 1), Block(b.source, b.start + o, b.end, 1)]
    return [Block(b.source, b.end - o, b.end, -1), Block(b.source, b.start, b.end - o, -1)]


@dataclass
class Chromosome:
    """One replicon: an ordered list of :class:`Block`\\ s over a nucleotide coordinate axis, with an
    ``id`` and a ``topology`` — ``"circular"`` (a ring, where coordinate ``length`` wraps to ``0`` and
    an arc may cross the origin) or ``"linear"`` (two ends, no wrap). Owns the per-chromosome
    operations. Blocks are kept **maximal**: after an event any collinear neighbours are merged."""

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
        """Invert the arc ``[start, start+length)``: reverse the block order and flip each strand, then
        re-canonicalise. Ancestry is unchanged."""
        span = self._arc_range(start, length)
        if span is None:
            return
        i, j = span
        arc = self.blocks[i:j]
        arc.reverse()
        for b in arc:
            b.strand = -b.strand
        self.blocks[i:j] = arc
        self._canonicalize()

    def _canonicalize(self) -> None:
        """Merge adjacent collinear blocks so every block is maximal (one unbroken ancestry). Two
        blocks of the same source and strand are collinear when their source coordinates are
        contiguous in reading order (forward: ``a.end == b.start``; reversed: ``a.start == b.end``)."""
        merged: list[Block] = []
        for b in self.blocks:
            if merged and merged[-1].source == b.source and merged[-1].strand == b.strand:
                a = merged[-1]
                if a.strand == 1 and a.end == b.start:
                    merged[-1] = Block(a.source, a.start, b.end, 1)
                    continue
                if a.strand == -1 and a.start == b.end:
                    merged[-1] = Block(a.source, b.start, a.end, -1)
                    continue
            merged.append(b)
        self.blocks = merged

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


@dataclass
class NucleotideGenomesResult:
    """What :func:`simulate_genomes_nucleotide` returns: the ``complete_tree`` it ran on, the final
    nucleotide ``genomes`` (karyotypes) at **every** node, the ``rearrangements`` (inversion) log, the
    ``chromosome_events`` (the chromosome network — tree-shaped until the tier reticulates it), and
    the ``seed``. ``mosaic`` / ``trace_back`` / ``ancestry`` read a node's genome."""

    complete_tree: Tree
    genomes: dict[int, NucleotideGenome]
    rearrangements: list[Inversion]
    chromosome_events: list[ChromosomeEvent]
    seed: int | None

    def mosaic(self, node_id: int) -> dict[int, list[tuple[int, int, int, int]]]:
        return self.genomes[node_id].mosaic()

    def trace_back(self, node_id: int) -> dict[int, list[tuple[int, int, int]]]:
        return self.genomes[node_id].trace_back()

    def ancestry(self, node_id: int) -> list[tuple[int, int]]:
        return self.genomes[node_id].ancestry()


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


def _copy_chromosome(c: Chromosome, cid: int) -> Chromosome:
    """A daughter chromosome: a fresh id, a deep copy of the blocks (fresh :class:`Block` objects, so
    a daughter's inversions never mutate the parent's genome)."""
    return Chromosome(cid, c.topology, [Block(b.source, b.start, b.end, b.strand) for b in c.blocks])


def _evolve_branch(g, node_id, t0, t1, inversion, inversion_length, rng, rearrangements) -> None:
    """Fire a homogeneous Poisson process of inversions along the branch ``[t0, t1]``. An inversion
    conserves length, so the total per-nucleotide rate ``inversion × length`` is constant on the
    branch; each inversion lands on a length-weighted chromosome, spanning a geometric run of mean
    ``inversion_length`` nucleotides at a uniform start within that chromosome."""
    total_rate = inversion * g.length
    if total_rate <= 0 or t1 <= t0:
        return
    t = t0 + float(rng.exponential(1.0 / total_rate))
    while t < t1:
        chrom, pos = g._pick_position(rng)
        length = min(chrom.length, int(rng.geometric(1.0 / inversion_length)))
        chrom.invert(pos, length)
        rearrangements.append(Inversion(t, node_id, chrom.id, pos, length))
        t += float(rng.exponential(1.0 / total_rate))


def simulate_genomes_nucleotide(tree, *, inversion=0.0, inversion_length=50.0, chromosomes=1,
                                root_length=1000, topology="circular", seed=None
                                ) -> NucleotideGenomesResult:
    """Evolve a nucleotide genome along a species tree by **inversion only** (the tree-wiring +
    multi-chromosome step). The root is seeded with a **karyotype** — ``chromosomes`` replicons, each
    its own source: an int ``N`` gives ``N`` equal replicons of ``root_length``/``topology``, or pass
    a list of ``(length, topology)`` for heterogeneous **sizes and shapes**. Each lineage inherits a
    copy of its parent's karyotype at speciation, with **every chromosome re-minted** (the edge
    recorded in ``chromosome_events``, the chromosome network), then accumulates inversions along its
    branch: ``inversion`` is a **per-nucleotide** rate (total ``inversion × length``, constant since
    an inversion conserves length), an inversion landing on a length-weighted chromosome and spanning
    a geometric run of mean ``inversion_length``. Deterministic given ``seed``.

    No coupling yet ⇒ the tree is walked branch by branch; no ancestry-changing event yet ⇒ **every**
    node carries the whole root sequence, merely permuted across its chromosomes."""
    tree = tree.complete_tree if isinstance(tree, SpeciesResult) else tree
    if inversion < 0:
        raise ValueError(f"inversion must be >= 0, got {inversion}")
    if inversion_length <= 0:
        raise ValueError(f"inversion_length must be > 0, got {inversion_length}")
    specs = _replicon_specs(chromosomes, root_length, topology)

    rng = np.random.default_rng(seed)
    chrom_counter = 0

    def new_chrom_id() -> int:
        nonlocal chrom_counter
        cid = chrom_counter
        chrom_counter += 1
        return cid

    genomes: dict[int, NucleotideGenome] = {}
    rearrangements: list[Inversion] = []
    chromosome_events: list[ChromosomeEvent] = []
    root = tree.nodes[tree.root]

    root_chroms = []
    for source, (length, top) in enumerate(specs):     # one source per seed replicon; each a network root
        cid = new_chrom_id()
        root_chroms.append(Chromosome(cid, top, [Block(source, 0, length, 1)]))
        chromosome_events.append(ChromosomeEvent(root.birth_time, "origination", root.id, (), (cid,)))
    start_genome = {root.id: NucleotideGenome(root_chroms)}

    stack = [root.id]
    while stack:  # DFS from the root: a parent's final genome exists before its children inherit it
        node_id = stack.pop()
        node = tree.nodes[node_id]
        g = start_genome.pop(node_id)
        _evolve_branch(g, node_id, node.birth_time, node.end_time, inversion, inversion_length,
                       rng, rearrangements)
        genomes[node_id] = g
        if node.children is not None:  # speciation: re-mint every chromosome into each daughter
            starts: dict[int, list[Chromosome]] = {c: [] for c in node.children}
            for pchrom in g.chromosomes:
                dcids = []
                for c in node.children:
                    dcid = new_chrom_id()
                    dcids.append(dcid)
                    starts[c].append(_copy_chromosome(pchrom, dcid))
                chromosome_events.append(
                    ChromosomeEvent(node.end_time, "speciation", node.id, (pchrom.id,), tuple(dcids)))
            for c in node.children:
                start_genome[c] = NucleotideGenome(starts[c])
                stack.append(c)
    return NucleotideGenomesResult(tree, genomes, rearrangements, chromosome_events, seed)


__all__ = ["Block", "Chromosome", "NucleotideGenome", "NucleotideGenomesResult",
           "simulate_genomes_nucleotide"]
