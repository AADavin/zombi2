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

All of the above are **ancestry-neutral**, so every node still carries the whole root sequence, merely
permuted across a changing number of chromosomes — the strong invariant carried through the scary
boundary-crossing code. Deferred to later slices: loss / duplication / transfer / transposition /
origination (the events that birth and kill ancestry, creating per-block gene trees) and chromosome
origination / loss; indels; declared genes / intergenes (pseudogenization, replacement); GFF input and
BED / FASTA output.
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


@dataclass
class NucleotideGenomesResult:
    """What :func:`simulate_genomes_nucleotide` returns: the ``complete_tree`` it ran on, the final
    nucleotide ``genomes`` (karyotypes) at **every** node, the ``rearrangements`` (inversion) log, the
    ``chromosome_events`` (the chromosome network — tree-shaped until the tier reticulates it), and
    the ``seed``. ``mosaic`` / ``trace_back`` / ``ancestry`` read a node's genome."""

    complete_tree: Tree
    genomes: dict[int, NucleotideGenome]
    rearrangements: list[Inversion | Translocation]
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


@dataclass(frozen=True)
class _Rates:
    inversion: float
    translocation: float
    fission: float
    fusion: float
    inversion_length: float
    translocation_length: float
    inversion_probability: float


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
    source._canonicalize()
    flipped = bool(rng.random() < inversion_probability)
    if flipped:
        arc = [Block(b.source, b.start, b.end, -b.strand) for b in reversed(arc)]
    others = [c for c in g.chromosomes if c is not source]
    dest = others[int(rng.integers(len(others)))]
    pos = int(rng.integers(dest.length + 1))
    dest._split_at(pos)
    k = dest._index_at(pos)
    dest.blocks[k:k] = arc
    dest._canonicalize()
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
    a._canonicalize()
    b._canonicalize()
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
    fused._canonicalize()
    g.chromosomes[:] = [c for k, c in enumerate(g.chromosomes) if k not in (ci, cj)] + [fused]
    chromosome_events.append(ChromosomeEvent(t, "fusion", node_id, (a.id, b.id), (fused.id,)))


def _evolve_branch(g, node_id, t0, t1, rates, rng, rearrangements, chromosome_events,
                   new_chrom_id) -> None:
    """A Gillespie over the branch ``[t0, t1]``. Every event conserves total length, so the
    per-nucleotide inversion rate is constant, but fission and fusion change the chromosome *count*,
    so the per-chromosome rates vary — the loop recomputes after each. All within-lineage (no
    coupling yet); the global timeline waits for transfer."""
    if t1 <= t0:
        return
    t = t0
    while True:
        nc = len(g.chromosomes)
        length = g.length
        r_inv = rates.inversion * length
        r_trl = rates.translocation * length if nc >= 2 else 0.0
        r_fis = rates.fission * nc
        r_fus = rates.fusion * nc if nc >= 2 else 0.0
        total = r_inv + r_trl + r_fis + r_fus
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
        elif r < r_inv + r_trl + r_fis:
            _do_fission(g, node_id, t, rng, chromosome_events, new_chrom_id)
        else:
            _do_fusion(g, node_id, t, rng, chromosome_events, new_chrom_id)


def simulate_genomes_nucleotide(tree, *, inversion=0.0, inversion_length=50.0, translocation=0.0,
                                translocation_length=50.0, inversion_probability=0.0, fission=0.0,
                                fusion=0.0, chromosomes=1, root_length=1000, topology="circular",
                                seed=None) -> NucleotideGenomesResult:
    """Evolve a nucleotide genome along a species tree by inversion, translocation, and the
    **number-changing chromosome tier**. The root is seeded with a **karyotype** — ``chromosomes``
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
    - ``fission`` (**per chromosome**) splits a chromosome in two (a **bifurcation**); ``fusion``
      (**per chromosome**) merges two chromosomes of the same topology (the **reticulation**). Both
      re-mint their children and record a network edge.

    Every event conserves total length, but fission/fusion change the chromosome count, so the branch
    is a small Gillespie recomputing its rates as the karyotype changes (still within-lineage — the
    global timeline waits for transfer). No ancestry-changing event yet ⇒ **every** node still carries
    the whole root sequence, merely permuted across its chromosomes. Deterministic given ``seed``."""
    tree = tree.complete_tree if isinstance(tree, SpeciesResult) else tree
    for label, rate in (("inversion", inversion), ("translocation", translocation),
                        ("fission", fission), ("fusion", fusion)):
        if rate < 0:
            raise ValueError(f"{label} must be >= 0, got {rate}")
    for label, mean in (("inversion_length", inversion_length),
                        ("translocation_length", translocation_length)):
        if mean <= 0:
            raise ValueError(f"{label} must be > 0, got {mean}")
    if not 0.0 <= inversion_probability <= 1.0:
        raise ValueError(f"inversion_probability must be in [0, 1], got {inversion_probability}")
    specs = _replicon_specs(chromosomes, root_length, topology)
    rates = _Rates(inversion, translocation, fission, fusion, inversion_length, translocation_length,
                   inversion_probability)

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
        _evolve_branch(g, node_id, node.birth_time, node.end_time, rates, rng, rearrangements,
                       chromosome_events, new_chrom_id)
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
