"""Genomes III — nucleotide: a genome as a coordinate space of base pairs.

The third and hardest resolution. Where an ordered genome is a *list of gene tokens*, a nucleotide
genome is a **karyotype of chromosomes**, each a coordinate axis of base pairs represented as an
ordered list of **blocks**. The vocabulary (fixed here, once):

- **block** — the persistent unit of both the representation *and* the ancestry: a run of one unbroken
  ancestry, the half-open interval ``[start, end)`` on an ancestral ``source``, read forward (``strand``
  ``+1``) or reverse-complemented (``-1``). Blocks are **not** merged during the run — an event only
  ever *splits* — so maximality is a property of the *recovered* root-blocks, not the working ones.
- **gene / intergene** — a *classification* of blocks: a **gene** is a *declared, indivisible* block
  (one family, one id, **never split**) carrying one gene tree; an **intergene** is the spacer, which
  fragments freely into many blocks.
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
- **Translocation / transposition** — an arc moved to a *different* chromosome (translocation) or
  *within* its own (transposition), optionally inverted (probability ``inversion_probability``). Both
  keep source coordinates, so they are *rearrangements*, not network edges.

Those are all **ancestry-neutral** (the strong invariant: every node carries the whole root sequence,
permuted). The ancestry-**changing** events are here too:

- **Loss** — an arc deleted (per lineage, never emptying a chromosome). A death, in the ``events``
  (genealogy) log. On its own it weakens the invariant to *subset*: each ancestral position at most
  once, monotonically down every path.
- **Duplication** — an arc copied in **tandem** (per lineage). A *birth*: the copied material now
  has an extra copy (same source coordinates), so a position can appear more than once.
- **Transfer** — an arc copied into a **contemporaneous recipient** lineage (per lineage,
  additive). A *horizontal* birth that couples two lineages — the reason the whole engine runs on a
  **global timeline** (all lineages alive at once, one clock) rather than branch-by-branch.
- **Origination** — a **new gene** laid down de novo on a fresh source (per lineage): its own family,
  indivisible from birth, rooting its own gene tree at its own branch. Origination mints a gene, never
  plain spacer.

Threaded through all of the above is the **copy lineage** (``Block.copy``): the persistent identity
that a split preserves, a duplication / transfer mints fresh (parent → child), and a speciation
re-mints into each daughter. The genealogy ``events`` — ``origination`` / ``loss`` / ``duplication`` /
``transfer`` / ``speciation`` — carry these copy ids, and the recovery reads them:

- **The gene-tree recovery** (``result.root_blocks`` / ``result.gene_trees``) — the **root partition**
  (cut each source at the union of extant-leaf breakpoints, keep the intervals that survive in a leaf:
  each is a maximal never-cut interval = one gene family) and, per block, the copy-lineage log replayed
  into one gene tree (a duplication is a ladder, a transfer a horizontal edge, a speciation a
  bifurcation), reusing the shared per-segment builder. Validated end to end: the recovered extant tips
  equal the copies actually present in every extant leaf.

Finally, the **genic layer** (``genes`` / ``gene_length``) declares some blocks to be **genes**, the
rest **intergenes**. There is exactly **one** way an event picks its target — an *extension* — and it is
drawn **directly from the legal arcs**, never guessed and retried: the event **nucleates at a uniform
intergenic position**, and its far end is chosen among the positions where a breakpoint is legal (never
strictly inside a gene), **weighted by** ``exp(-d / mean)`` — the extent distribution you asked for,
restricted to the ends that exist. Landing spots and chromosome cuts are picked the same way. Nothing is
clipped, nothing is snapped (snapping would inflate an extent to "whatever it takes to clear the gene"),
and nothing is silently dropped for running out of retries. So a gene is only ever engulfed **whole**,
never split, and each recovers as exactly one root-block with one gene tree.

**Two consequences the guide must state plainly.** (1) Gene turnover is *emergent* and size-dependent —
a gene changes copy number only when an event engulfs it whole, so large genes are rarely lost or
duplicated. (2) The **realised extent is shorter than the mean you ask for**, the more so the denser the
genome: on a 94%-genic bacterial genome, asking for 3 000 bp yields ~1 300. A long stretch that both
begins *and* ends in a spacer often simply does not exist, and an event that would cut a gene does not
happen. That conditioning is the model; what it must not do — and no longer does — is quietly eat the
event *rate* along with it.

Deferred to later slices: homologous *replacement* transfer (only additive for now); indels;
pseudogenization (``gene → intergene``); declaring genes from a GFF / distribution file; BED / FASTA
output; and the opt-in per-copy dial (size-blind, settable per-gene turnover — a *second* selection
method, deliberately kept out). (The recovery above builds the *extant* tree exactly; the *complete*
tree's dead partial-overlap lineages need the finer all-node partition — a later refinement.)
"""

from __future__ import annotations

import collections
import math
from dataclasses import dataclass, field

import numpy as np

from ..species import SpeciesResult, Tree
from ._live import enter, retire
from ._transfer import Distance, mean_root_to_tip, recipient_index
from .chromosomes import ChromosomeEvent
from .gene_trees import GeneTree, gene_trees_from_events
from .gff import read_gff


@dataclass
class Block:
    """A run of one unbroken ancestry: the half-open interval ``[start, end)`` on ``source`` this run
    descends from (``start < end`` always), read forward (``strand`` ``+1``) or reverse-complemented
    (``-1``). Its physical length is ``end - start``. ``copy`` names the **copy lineage** it belongs to
    — the persistent identity that threads through splits and moves (a split preserves it, a
    duplication mints a fresh child), re-minted into daughters at speciation; the gene-tree recovery
    reads it to resolve *which copy begat which*.

    ``gene`` is the **genic classification**: ``0`` for an **intergene** (the free-fragmenting spacer),
    or the **gene family id** for a **gene** — a *declared, indivisible* block that is **never split**
    (an event that would cut one redraws instead). A duplicated gene keeps its family and gets a fresh
    ``copy``.

    Blocks are **not** kept maximal during the run — a rearrangement leaves collinear neighbours split
    rather than merging them (merging would muddy the copy-lineage bookkeeping). Maximality is instead
    a property of the *recovered* root-blocks: the coarsest partition that survives in the leaves."""

    source: int
    start: int
    end: int
    strand: int
    copy: int = 0
    gene: int = 0

    @property
    def length(self) -> int:
        return self.end - self.start

    @property
    def is_gene(self) -> bool:
        """Whether this block is a declared gene (indivisible) rather than intergenic spacer."""
        return self.gene != 0


class _CutsGene(Exception):
    """Raised when a breakpoint would fall **strictly inside a gene** — a **guard, not a control flow**.

    Every breakpoint is now drawn from the legal set to begin with (:meth:`Chromosome._pick_legal_cut`,
    :meth:`Chromosome._pick_arc_extent`), so this should never fire; if it does, a caller invented a
    position instead of picking one, and failing loudly beats silently mangling a gene."""


def _split_block(b: Block, o: int) -> list[Block]:
    """Split ``b`` after ``o`` physical positions into ``[left, right]``, strand-aware. For a forward
    block the cut falls at ``start + o``; for a reversed one the first ``o`` physical positions are
    the **high** end of the source, so the cut falls at ``end - o``. Both pieces keep ``b``'s copy
    lineage and genic tag — a split is not a birth."""
    if b.strand == 1:
        return [Block(b.source, b.start, b.start + o, 1, b.copy, b.gene),
                Block(b.source, b.start + o, b.end, 1, b.copy, b.gene)]
    return [Block(b.source, b.end - o, b.end, -1, b.copy, b.gene),
            Block(b.source, b.start, b.end - o, -1, b.copy, b.gene)]


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

    def _check_cut(self, c: int) -> None:
        """Raise :class:`_CutsGene` if a breakpoint at ``c`` would fall **strictly inside a gene**.
        Pure — mutates nothing — so a caller can test both ends of an arc *before* splitting either."""
        if c <= 0:
            return
        pos = 0
        for b in self.blocks:
            if pos == c:
                return
            if pos < c < pos + b.length:
                if b.is_gene:
                    raise _CutsGene
                return
            pos += b.length

    def _split_at(self, c: int) -> None:
        """Ensure a block boundary at physical coordinate ``c`` (``0 <= c <= length``). A no-op at the
        ends or an existing boundary; otherwise the block straddling ``c`` is split — unless that block
        is a **gene**, which is indivisible, in which case :class:`_CutsGene` is raised (and nothing is
        mutated) so the event redraws."""
        if c <= 0:
            return
        pos = 0
        for i, b in enumerate(self.blocks):
            if pos == c:
                return
            if pos < c < pos + b.length:
                if b.is_gene:
                    raise _CutsGene
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

    def _legal_cut_spans(self) -> list[tuple[int, int]]:
        """The physical windows where a breakpoint is legal — the intergenic blocks, endpoints included
        (a gene's own edge is legal; only its interior is not). ``[(lo, hi), …]`` in block order."""
        out, pos = [], 0
        for b in self.blocks:
            if not b.is_gene:
                out.append((pos, pos + b.length))
            pos += b.length
        return out

    def _pick_legal_cut(self, rng) -> int | None:
        """A **uniform** position at which a breakpoint is legal — where an arc may land, or a
        chromosome be cut, without splitting a gene. ``None`` when the chromosome is all gene."""
        spans = self._legal_cut_spans()
        total = sum(hi - lo for lo, hi in spans)
        if total == 0:
            return None
        m = int(rng.integers(total))
        for lo, hi in spans:
            if m < hi - lo:
                return lo + m
            m -= hi - lo
        raise AssertionError("legal-cut length out of sync with the blocks")  # unreachable

    def _pick_arc_extent(self, start: int, mean: float, rng) -> int | None:
        """Choose the arc's far end, forward from ``start``: an extent ``d >= 1`` whose breakpoint at
        ``start + d`` is **legal** (never strictly inside a gene), drawn with weight ``exp(-d/mean)``.

        This is the extent distribution you asked for, **restricted to the ends that exist** — the same
        law that draw-then-reject would converge to, but sampled directly, so nothing is ever wasted and
        no event silently vanishes because its retries ran out. On a gene-dense genome the realised
        extent therefore comes out *shorter* than ``mean``: a long stretch that both begins and ends in a
        spacer often simply does not exist. That conditioning is the model (an event that would cut a
        gene does not happen); what it must not do is quietly eat the event *rate* as well.

        ``None`` when no legal end exists."""
        total = self.length
        if total < 2:
            return None
        circular = self.topology == "circular"
        limit = total - 1 if circular else total - start     # the longest arc that still fits
        if limit < 1:
            return None
        spans = self._legal_cut_spans()
        if not spans:
            return None
        if circular:                                         # the ring: the same windows one lap on
            spans = spans + [(lo + total, hi + total) for lo, hi in spans]
        windows = []
        for lo, hi in spans:
            d_lo, d_hi = max(lo - start, 1), min(hi - start, limit)
            if d_lo <= d_hi:
                windows.append((d_lo, d_hi))
        if not windows:
            return None
        q = math.exp(-1.0 / mean)                            # a geometric extent, summed per window
        weights = [max(q ** d_lo - q ** (d_hi + 1), 0.0) for d_lo, d_hi in windows]
        bulk = sum(weights)
        if bulk <= 0.0:                                      # every window is far beyond `mean`
            return min(windows, key=lambda w: w[0])[0]       # ...so take the nearest legal end
        r = float(rng.random()) * bulk
        for (d_lo, d_hi), w in zip(windows, weights):
            if r < w:
                head, tail = q ** d_lo, q ** (d_hi + 1)      # inverse-CDF inside the window
                u = float(rng.random())
                target = head - u * (head - tail)
                d = int(math.log(target) / math.log(q)) - 1 if target > 0 else d_hi
                return min(max(d, d_lo), d_hi)
            r -= w
        return windows[-1][1]                                # floating-point guard

    def _arc_range(self, start: int, length: int) -> tuple[int, int] | None:
        """Prepare the arc ``[start, start+length)`` and return the block index range ``[i, j)`` it
        occupies, or ``None`` for an empty arc. Splits at both ends first. A **linear** chromosome
        clamps the arc to its ends; a **circular** one may wrap the origin, rotating the ring to bring
        the arc to the front (the origin drifts to a real breakpoint — harmless on a ring).

        Both ends are checked **before** either is split, so an illegal arc leaves the chromosome
        untouched rather than half-cut. With the ends drawn from the legal set this never triggers; it
        is a guard."""
        total = self.length
        if total == 0:
            return None
        if self.topology == "linear":
            start = max(0, min(start, total))
            end = min(start + max(0, length), total)
            if end <= start:
                return None
            self._check_cut(start)
            self._check_cut(end)
            self._split_at(start)
            self._split_at(end)
            return self._index_at(start), self._index_at(end)
        ell = max(1, min(length, total))
        s = start % total
        self._check_cut(s)
        self._check_cut((s + ell) % total)
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

    def _pick_intergenic_position(self, rng) -> tuple[Chromosome, int] | None:
        """A uniform pick over the genome's **intergenic** nucleotides → ``(chromosome, physical
        position)``. An event nucleates in the spacer, never inside a gene — genes are indivisible, so
        starting inside one could only ever be redrawn away. With no genes declared every nucleotide is
        intergenic and this is just a uniform pick. ``None`` when the genome is all gene (nowhere to
        start)."""
        total = sum(b.length for c in self.chromosomes for b in c.blocks if not b.is_gene)
        if total == 0:
            return None
        m = int(rng.integers(total))
        for c in self.chromosomes:
            pos = 0
            for b in c.blocks:
                if not b.is_gene:
                    if m < b.length:
                        return c, pos + m
                    m -= b.length
                pos += b.length
        raise AssertionError("intergenic length out of sync with the blocks")  # unreachable

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
class Transposition:
    """A recorded nucleotide transposition: on branch ``lineage`` at ``time``, the arc
    ``[start, start+length)`` of chromosome ``chromosome`` was excised and reinserted **elsewhere on
    the same chromosome** (at physical position ``dest`` of the remainder), landing ``flipped``
    (reversed) or not. Ancestry is unchanged — the blocks keep their source coordinates — so it is a
    rearrangement, not a gene-genealogy event."""

    time: float
    lineage: int
    chromosome: int
    start: int
    length: int
    dest: int
    flipped: bool


@dataclass(frozen=True)
class Origination:
    """A recorded **birth of a copy lineage** — the root of a gene tree. At ``time`` on branch
    ``lineage`` new material was laid down on chromosome ``chromosome`` as copy lineage ``copy``,
    covering the ancestral interval ``[start, end)`` on ``source``. Either a **seed** replicon (one per
    root replicon, on the root branch) or a **de-novo** origination (a fresh source arising mid-tree);
    the gene-tree recovery reads each as the root of its family."""

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
class Transfer:
    """A recorded nucleotide transfer — an ancestry-**changing** *birth* that couples two
    contemporaneous lineages, so it is a **horizontal** edge in a block's gene tree. At ``time`` an arc
    on the donor branch ``lineage`` was copied into the recipient branch ``recipient`` (additive — the
    donor keeps its copy). ``transferred`` names the parentage as ``(parent_copy, child_copy, source,
    start, end)`` rows — the donor's ``parent_copy`` begat the recipient's fresh ``child_copy`` over
    ``[start, end)`` of ``source``."""

    time: float
    lineage: int
    recipient: int
    transferred: tuple[tuple[int, int, int, int, int], ...]


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
    (``origination``, ``loss``, ``duplication``, ``transfer``, ``speciation`` — carrying the copy ids
    the gene-tree recovery reads), the ancestry-neutral ``rearrangements`` (inversion, translocation), the
    ``chromosome_events`` (the chromosome network), and the ``seed``. ``mosaic`` / ``trace_back`` /
    ``ancestry`` read a node's genome."""

    complete_tree: Tree
    genomes: dict[int, NucleotideGenome]
    events: list[Origination | Loss | Duplication | Transfer | Speciation]
    rearrangements: list[Inversion | Translocation | Transposition]
    chromosome_events: list[ChromosomeEvent]
    seed: int | None
    #: ``{gene family id: (source, start, end)}`` — where each **declared gene** sits in the root
    #: coordinate space. Empty when no genes were declared. A gene is never split, so its span is
    #: fixed for the whole run and is exactly the root-block that carries its gene tree.
    gene_spans: dict[int, tuple[int, int, int]] = field(default_factory=dict)
    #: ``{name: gene family id}`` for genes declared with a name (a GFF ``ID`` / ``Name``) — the handle
    #: to look a named gene up in :attr:`gene_spans` / :attr:`gene_trees`. Empty for the even layout.
    gene_names: dict[str, int] = field(default_factory=dict)
    #: ``{gene family id: +1 / -1}`` — each declared gene's **coding** strand (which strand carries the
    #: ORF), as given by the GFF. This is annotation, *not* ancestry: it is fixed for the family and is
    #: unrelated to :attr:`Block.strand`, which records whether a stretch has been inverted since the
    #: root. The even layout declares every gene on ``+1``.
    gene_strands: dict[int, int] = field(default_factory=dict)

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
        """``{family: GeneTree}`` — the recovered gene trees.

        With **genes declared**, one tree per gene, keyed by its **gene family id** (see
        :attr:`gene_spans`); the intergenic root-blocks keep their block ancestry in the log but are not
        built into trees. With **no genes declared** the whole genome is one big intergene, so every
        recovered root-block is a family in its own right and the key is its index in
        :attr:`root_blocks`. A declared gene that survives in no extant leaf has no root-block and so no
        tree."""
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


def _even_gene_intervals(length, genes, gene_length) -> list[tuple[int, int, int, str | None]]:
    """The evenly-spaced layout: ``genes`` genes of ``gene_length``, with the leftover spread as evenly
    as possible over the intergenes that precede them. Returns ``(start, end, strand, name)`` intervals,
    0-based half-open, in order."""
    if genes == 0:
        return []
    base, extra = divmod(length - genes * gene_length, genes)   # `extra` intergenes are 1 bp longer
    out, at = [], 0
    for i in range(genes):
        at += base + (1 if i < extra else 0)
        out.append((at, at + gene_length, 1, None))
        at += gene_length
    return out


def _seed_blocks(source, length, cp, intervals, new_family, gene_spans, gene_names,
                 gene_strands) -> list[Block]:
    """Lay one seed replicon of ``length`` down as its blocks: the declared genes at ``intervals``
    (0-based half-open, sorted, non-overlapping) and **intergene** everywhere else — the alternating
    chain. Every block shares the replicon's seed copy lineage ``cp`` (they are one copy of one
    replicon); each **gene** additionally gets a fresh **family** id, which is what makes it indivisible
    and gives it a gene tree, and is recorded in ``gene_spans``, ``gene_strands`` (its **coding**
    strand) and ``gene_names`` (when it is named).
    With no intervals the replicon is a single intergenic block — today's uniform sequence."""
    blocks, at = [], 0
    for (start, end, strand, name) in intervals:
        if start > at:
            blocks.append(Block(source, at, start, 1, cp))       # intergene before this gene
        fam = new_family()
        # NB: every seed block is strand +1. `Block.strand` is orientation *relative to the ancestral
        # source*, and at the root the genome IS its own source — nothing is inverted yet. A gene's
        # coding strand from the GFF is a different thing entirely (which strand carries the ORF), a
        # constant property of the family, so it is recorded separately.
        blocks.append(Block(source, start, end, 1, cp, fam))
        gene_spans[fam] = (source, start, end)
        gene_strands[fam] = strand
        if name is not None:
            gene_names[name] = fam
        at = end
    if at < length:
        blocks.append(Block(source, at, length, 1, cp))           # trailing intergene
    return blocks or [Block(source, 0, length, 1, cp)]


def _copy_chromosome(c: Chromosome, cid: int, copy_map: dict[int, int]) -> Chromosome:
    """A daughter chromosome: a fresh id, a deep copy of the blocks (fresh :class:`Block` objects, so
    a daughter's inversions never mutate the parent's genome), with every block's copy lineage
    re-minted through ``copy_map`` (parent copy id → this daughter's fresh copy id)."""
    return Chromosome(cid, c.topology,
                      [Block(b.source, b.start, b.end, b.strand, copy_map[b.copy], b.gene) for b in c.blocks])


@dataclass(frozen=True)
class _Rates:
    inversion: float
    translocation: float
    transposition: float
    loss: float
    duplication: float
    transfer: float
    origination: float
    fission: float
    fusion: float
    chromosome_origination: float
    chromosome_loss: float
    inversion_length: float
    translocation_length: float
    transposition_length: float
    loss_length: float
    duplication_length: float
    transfer_length: float
    origination_length: float
    inversion_probability: float


def _do_duplication(g, node_id, t, duplication_length, rng, events, new_copy) -> int:
    """Copy a geometric-length arc of a length-weighted chromosome **in tandem** (the copy inserted
    right after the arc). An ancestry-changing *birth*: the copied material now has an extra copy (same
    source coordinates). Each distinct copy lineage in the arc begets one fresh child lineage, so the
    tandem copy is a new copy of that material; the parentage is recorded as a :class:`Duplication`.
    Returns the length added (0 on a no-op)."""
    spot = g._pick_intergenic_position(rng)
    if spot is None:
        return 0
    chrom, start = spot
    ell = chrom._pick_arc_extent(start, duplication_length, rng)
    if ell is None:
        return 0
    span = chrom._arc_range(start, ell)
    if span is None:
        return 0
    i, j = span
    arc = chrom.blocks[i:j]
    child_of: dict[int, int] = {}
    for b in arc:
        if b.copy not in child_of:
            child_of[b.copy] = new_copy()
    copied = tuple((b.copy, child_of[b.copy], b.source, b.start, b.end) for b in arc)
    chrom.blocks[j:j] = [Block(b.source, b.start, b.end, b.strand, child_of[b.copy], b.gene)
                         for b in arc]
    events.append(Duplication(t, node_id, chrom.id, copied))
    return sum(b.length for b in arc)


def _do_transfer(rng, tree, alive, gen, kd, t, transfer_length, transfer_to, self_transfer, depth,
                 events, new_copy) -> int:
    """Copy a geometric-length arc of the donor lineage ``alive[kd]`` into a **contemporaneous
    recipient** (chosen by ``transfer_to``: uniform, or a :class:`Distance` weighting): the arc's copy
    lineages beget fresh children that arrive as a block at a random spot on a random recipient
    chromosome (strands travel with them). A horizontal edge in each block's gene tree. **Additive**
    — the donor keeps its copy — so it returns the recipient's length gain (0 on a no-op)."""
    donor_g = gen[kd]
    spot = donor_g._pick_intergenic_position(rng)
    if spot is None:
        return 0
    chrom, start = spot
    ell = chrom._pick_arc_extent(start, transfer_length, rng)
    if ell is None:
        return 0
    span = chrom._arc_range(start, ell)
    if span is None:
        return 0
    i, j = span
    arc = chrom.blocks[i:j]                              # additive: the arc stays on the donor
    donor = alive[kd]
    cand = [k for k in range(len(alive)) if self_transfer or k != kd]
    if not cand:
        return 0
    kr = recipient_index(rng, tree, alive, cand, donor, t, transfer_to, depth)
    recipient = alive[kr]
    rgenome = gen[kr]
    child_of: dict[int, int] = {}                        # each donor copy lineage begets one fresh child
    for b in arc:
        if b.copy not in child_of:
            child_of[b.copy] = new_copy()
    transferred = tuple((b.copy, child_of[b.copy], b.source, b.start, b.end) for b in arc)
    xfers = [Block(b.source, b.start, b.end, b.strand, child_of[b.copy], b.gene) for b in arc]
    rchrom = rgenome.chromosomes[int(rng.integers(len(rgenome.chromosomes)))]
    p = rchrom._pick_legal_cut(rng)                     # arrive at a legal spot on the recipient
    if p is None:
        return 0
    rchrom._split_at(p)
    q = rchrom._index_at(p)
    rchrom.blocks[q:q] = xfers
    events.append(Transfer(t, donor, recipient, transferred))
    return sum(b.length for b in arc)                   # the recipient's length gain


def _do_origination(g, node_id, t, origination_length, rng, events, new_source, new_copy,
                    new_family, gene_spans, gene_strands) -> int:
    """A **new gene** arises de novo: a fresh source (a geometric-length stretch, mean
    ``origination_length``), its own copy lineage and its own **gene family**, inserted at a random spot
    on a uniformly-chosen chromosome. Per lineage — a family is born once. Origination mints a *gene*,
    never plain spacer, so the new block is indivisible from birth and gets its own gene tree. Returns
    the length added."""
    chrom = g.chromosomes[int(rng.integers(len(g.chromosomes)))]
    length = max(1, int(rng.geometric(1.0 / origination_length)))
    src, cp, fam = new_source(), new_copy(), new_family()
    p = chrom._pick_legal_cut(rng)
    if p is None:
        return 0
    chrom._split_at(p)
    k = chrom._index_at(p)
    chrom.blocks[k:k] = [Block(src, 0, length, 1, cp, fam)]
    gene_spans[fam] = (src, 0, length)                   # a de-novo gene, tracked like a declared one
    gene_strands[fam] = 1
    events.append(Origination(t, node_id, chrom.id, cp, src, 0, length))
    return length


def _do_loss(g, node_id, t, loss_length, rng, events) -> int:
    """Delete a geometric-length arc from a length-weighted chromosome — an ancestry-changing event (a
    death). Never empties a chromosome (leaves at least one nucleotide; whole-chromosome loss is a
    deferred tier event). Records the deleted material — which copy lineage lost which arc — as a
    :class:`Loss`. Returns the length removed as a **negative** delta (0 on a no-op)."""
    spot = g._pick_intergenic_position(rng)
    if spot is None:
        return 0
    chrom, start = spot
    if chrom.length < 2:
        return 0
    ell = chrom._pick_arc_extent(start, loss_length, rng)
    if ell is None:
        return 0
    span = chrom._arc_range(start, ell)
    if span is None:
        return 0
    i, j = span
    gone = chrom.blocks[i:j]
    lost = tuple((b.copy, b.source, b.start, b.end) for b in gone)
    chrom.blocks = chrom.blocks[:i] + chrom.blocks[j:]
    events.append(Loss(t, node_id, chrom.id, lost))
    return -sum(b.length for b in gone)


def _do_inversion(g, node_id, t, inversion_length, rng, rearrangements) -> None:
    """Invert a geometric-length arc of a length-weighted chromosome (length-neutral)."""
    spot = g._pick_intergenic_position(rng)
    if spot is None:
        return
    chrom, pos = spot
    length = chrom._pick_arc_extent(pos, inversion_length, rng)
    if length is None:
        return
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
    spot = g._pick_intergenic_position(rng)
    if spot is None:
        return
    source, start = spot
    if source.length < 2:
        return
    ell = source._pick_arc_extent(start, translocation_length, rng)
    if ell is None:
        return
    span = source._arc_range(start, ell)
    if span is None:
        return
    i, j = span
    arc = source.blocks[i:j]
    intact = source.blocks                               # keep for rollback if the landing cuts a gene
    source.blocks = source.blocks[:i] + source.blocks[j:]
    flipped = bool(rng.random() < inversion_probability)
    if flipped:
        arc = [Block(b.source, b.start, b.end, -b.strand, b.copy, b.gene) for b in reversed(arc)]
    others = [c for c in g.chromosomes if c is not source]
    dest = others[int(rng.integers(len(others)))]
    pos = dest._pick_legal_cut(rng)                      # land where it will not split a gene
    if pos is None:
        source.blocks = intact                           # nowhere legal to land: undo the excision
        return
    dest._split_at(pos)
    k = dest._index_at(pos)
    dest.blocks[k:k] = arc
    rearrangements.append(Translocation(t, node_id, source.id, dest.id, start, ell, flipped))


def _do_transposition(g, node_id, t, transposition_length, inversion_probability, rng,
                      rearrangements) -> None:
    """Excise a geometric-length arc of a length-weighted chromosome and reinsert it **elsewhere on
    the same chromosome**, landing inverted with probability ``inversion_probability``. Ancestry-neutral
    (blocks keep their source coordinates). No-op below two nucleotides (an arc leaves at least one, so
    there is a landing spot)."""
    spot = g._pick_intergenic_position(rng)
    if spot is None:
        return
    chrom, start = spot
    if chrom.length < 2:
        return
    ell = chrom._pick_arc_extent(start, transposition_length, rng)
    if ell is None:
        return
    span = chrom._arc_range(start, ell)
    if span is None:
        return
    i, j = span
    arc = chrom.blocks[i:j]
    intact = chrom.blocks                                # keep for rollback if the landing cuts a gene
    chrom.blocks = chrom.blocks[:i] + chrom.blocks[j:]
    flipped = bool(rng.random() < inversion_probability)
    if flipped:
        arc = [Block(b.source, b.start, b.end, -b.strand, b.copy, b.gene) for b in reversed(arc)]
    dest = chrom._pick_legal_cut(rng)                   # a legal spot on the remainder
    if dest is None:
        chrom.blocks = intact                            # nowhere legal to land: undo the excision
        return
    chrom._split_at(dest)
    k = chrom._index_at(dest)
    chrom.blocks[k:k] = arc
    rearrangements.append(Transposition(t, node_id, chrom.id, start, ell, dest, flipped))


def _do_fission(g, node_id, t, rng, chromosome_events, new_chrom_id) -> int:
    """Split a uniformly-chosen chromosome into two re-minted children — a **bifurcation**. A linear
    chromosome cuts at one breakpoint (prefix + suffix); a circular one at two (the arc between them
    a new ring, the remainder another). Blocks keep their ancestry; no-op below two nucleotides.
    Returns the chromosome-count delta (``+1``, or ``0`` on a no-op)."""
    ci = int(rng.integers(len(g.chromosomes)))
    chrom = g.chromosomes[ci]
    if chrom.length < 2:
        return 0
    if chrom.topology == "linear":
        c = chrom._pick_legal_cut(rng)
        if c is None or c == 0:
            return 0
        chrom._split_at(c)
        i = chrom._index_at(c)
        a = Chromosome(new_chrom_id(), "linear", chrom.blocks[:i])
        b = Chromosome(new_chrom_id(), "linear", chrom.blocks[i:])
    else:
        cuts = {chrom._pick_legal_cut(rng), chrom._pick_legal_cut(rng)} - {None}
        if len(cuts) < 2:
            return 0
        c1, c2 = sorted(cuts)
        chrom._split_at(c1)
        chrom._split_at(c2)
        i, j = chrom._index_at(c1), chrom._index_at(c2)
        a = Chromosome(new_chrom_id(), "circular", chrom.blocks[i:j])
        b = Chromosome(new_chrom_id(), "circular", chrom.blocks[:i] + chrom.blocks[j:])
    g.chromosomes[ci:ci + 1] = [a, b]
    chromosome_events.append(ChromosomeEvent(t, "fission", node_id, (chrom.id,), (a.id, b.id)))
    return 1


def _do_fusion(g, node_id, t, rng, chromosome_events, new_chrom_id) -> int:
    """Merge a uniformly-chosen chromosome with another of the **same topology** — a **reticulation**
    (two parents, one child): concatenate their blocks, re-mint the child. No-op if the genome has no
    same-topology partner for the chosen chromosome. Returns the chromosome-count delta (``-1``, or
    ``0`` on a no-op)."""
    ci = int(rng.integers(len(g.chromosomes)))
    a = g.chromosomes[ci]
    partners = [k for k in range(len(g.chromosomes))
                if k != ci and g.chromosomes[k].topology == a.topology]
    if not partners:
        return 0
    cj = partners[int(rng.integers(len(partners)))]
    b = g.chromosomes[cj]
    fused = Chromosome(new_chrom_id(), a.topology, a.blocks + b.blocks)
    g.chromosomes[:] = [c for k, c in enumerate(g.chromosomes) if k not in (ci, cj)] + [fused]
    chromosome_events.append(ChromosomeEvent(t, "fusion", node_id, (a.id, b.id), (fused.id,)))
    return -1


def _do_chromosome_origination(g, node_id, t, chromosome_events, new_chrom_id) -> int:
    """A de-novo replicon (a plasmid): a fresh **empty** circular chromosome — a **root** of the
    chromosome network (no parent). It carries no sequence yet (length 0); material arrives later by
    origination / transfer / translocation. Returns the chromosome-count delta (``+1``)."""
    cid = new_chrom_id()
    g.chromosomes.append(Chromosome(cid, "circular", []))
    chromosome_events.append(ChromosomeEvent(t, "origination", node_id, (), (cid,)))
    return 1


def _do_chromosome_loss(g, node_id, t, rng, events, chromosome_events) -> tuple[int, int]:
    """A whole chromosome dies — a **leaf** of the chromosome network (no child): its material dies as
    one :class:`Loss` (each copy lineage on it ends). No-op if it is the genome's last chromosome (a
    lineage never loses its whole genome this way). Returns ``(chromosome delta, length delta)``."""
    if len(g.chromosomes) < 2:
        return (0, 0)
    ci = int(rng.integers(len(g.chromosomes)))
    lost = g.chromosomes[ci]
    if lost.blocks:                                     # an empty de-novo replicon dies with no material
        rows = tuple((b.copy, b.source, b.start, b.end) for b in lost.blocks)
        events.append(Loss(t, node_id, lost.id, rows))
    del g.chromosomes[ci]
    chromosome_events.append(ChromosomeEvent(t, "loss", node_id, (lost.id,), ()))
    return (-1, -lost.length)


def _pick_lineage_by_chromosomes(rng, gen, total_chromosomes) -> int:
    """Pick a lineage index proportional to its chromosome count — a global per-chromosome pick."""
    m = int(rng.integers(total_chromosomes))
    for k, g in enumerate(gen):
        if m < len(g.chromosomes):
            return k
        m -= len(g.chromosomes)
    raise AssertionError("total_chromosomes out of sync with the alive set")  # unreachable


def _speciate(node, g, new_chrom_id, new_copy, events, chromosome_events):
    """Re-mint a parent's final genome ``g`` into one fresh child karyotype per daughter species:
    every chromosome id and every copy lineage is re-minted (recorded as ``ChromosomeEvent`` /
    :class:`Speciation` edges). Returns ``{child species id: NucleotideGenome}``."""
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
    return {c: NucleotideGenome(starts[c]) for c in node.children}


def simulate_genomes_nucleotide(tree, *, inversion=0.0, inversion_length=50.0, translocation=0.0,
                                translocation_length=50.0, transposition=0.0, transposition_length=50.0,
                                inversion_probability=0.0, loss=0.0, loss_length=50.0, duplication=0.0,
                                duplication_length=50.0, transfer=0.0, transfer_length=50.0,
                                transfer_to="uniform", self_transfer=False, origination=0.0,
                                origination_length=50.0, fission=0.0, fusion=0.0,
                                chromosome_origination=0.0, chromosome_loss=0.0, chromosomes=1,
                                root_length=1000, topology="circular", genes=0, gene_length=100,
                                gff=None, trim_overlaps=False, seed=None) -> NucleotideGenomesResult:
    """Evolve a nucleotide genome along a species tree by inversion, translocation, transposition,
    **loss**, **duplication**, **transfer**, **origination**, and the number-changing chromosome tier.
    The root is seeded with a **karyotype** — ``chromosomes`` replicons, each its own source: an int
    ``N`` gives ``N`` equal replicons of ``root_length``/``topology``, or pass a list of ``(length,
    topology)`` for heterogeneous **sizes and shapes**. Each lineage inherits a copy of its parent's
    karyotype at speciation, with **every chromosome re-minted** (the chromosome network), and evolves:

    - ``inversion`` (**per lineage**) reverses a geometric-length (mean ``inversion_length``) arc
      of a length-weighted chromosome.
    - ``translocation`` (**per lineage**) moves a geometric-length (mean ``translocation_length``)
      arc to a **different** chromosome; ``transposition`` (**per lineage**, mean
      ``transposition_length``) moves one **within** its chromosome. Both land inverted with
      probability ``inversion_probability``, keep source coordinates, and are rearrangements, not edges.
    - ``loss`` (**per lineage**) deletes a geometric-length (mean ``loss_length``) arc — an
      ancestry-**changing** event (a death), recorded in ``events``. Never empties a chromosome.
    - ``duplication`` (**per lineage**) copies a geometric-length (mean ``duplication_length``) arc
      in tandem — an ancestry-**changing** *birth*, recorded in ``events``.
    - ``transfer`` (**per lineage**) copies a geometric-length (mean ``transfer_length``) arc into a
      **contemporaneous recipient** (``transfer_to``: ``"uniform"`` or ``"distance"`` / a
      :class:`Distance`; ``self_transfer`` allows the donor itself) — a horizontal *birth*, additive
      (the donor keeps its copy). This is what needs the global timeline.
    - ``origination`` (**per lineage**) lays down a **new gene** on a fresh source (geometric length,
      mean ``origination_length``) — a *birth* of a wholly new family, indivisible from birth.
    - ``fission`` (**per chromosome**) splits a chromosome in two (a **bifurcation**); ``fusion``
      (**per chromosome**) merges two chromosomes of the same topology (the **reticulation**).
      ``chromosome_origination`` (**per lineage**) adds a de-novo empty circular replicon (a plasmid, a
      network **root**); ``chromosome_loss`` (**per chromosome**) kills a whole chromosome (its material
      dies as a loss; never the last one) — a network **leaf**. All record a chromosome-network edge.

    The engine runs a **global-timeline** Gillespie: all lineages alive at once evolve along one clock
    (every extension event is **per lineage** — the rate says how often a lineage does it, the extent how
    much it touches, so a bigger genome does not get proportionally more events; the chromosome tier is
    per chromosome), so a transfer
    couples two contemporaries. With loss, the strong invariant weakens: every node carries a **subset**
    of the root sequence (each ancestral position at most once, monotonically down every path);
    origination further adds fresh sources beyond the root. Deterministic given ``seed``. (Transfer is
    additive for now; homologous *replacement* transfer is a later refinement.)"""
    tree = tree.complete_tree if isinstance(tree, SpeciesResult) else tree
    for label, rate in (("inversion", inversion), ("translocation", translocation),
                        ("transposition", transposition), ("loss", loss), ("duplication", duplication),
                        ("transfer", transfer), ("origination", origination), ("fission", fission),
                        ("fusion", fusion), ("chromosome_origination", chromosome_origination),
                        ("chromosome_loss", chromosome_loss)):
        if rate < 0:
            raise ValueError(f"{label} must be >= 0, got {rate}")
    for label, mean in (("inversion_length", inversion_length),
                        ("translocation_length", translocation_length),
                        ("transposition_length", transposition_length), ("loss_length", loss_length),
                        ("duplication_length", duplication_length), ("transfer_length", transfer_length),
                        ("origination_length", origination_length)):
        if mean <= 0:
            raise ValueError(f"{label} must be > 0, got {mean}")
    if not 0.0 <= inversion_probability <= 1.0:
        raise ValueError(f"inversion_probability must be in [0, 1], got {inversion_probability}")
    if transfer_to == "distance":
        transfer_to = Distance()
    if transfer_to != "uniform" and not isinstance(transfer_to, Distance):
        raise ValueError(f"transfer_to must be 'uniform', 'distance', or Distance(decay=), "
                         f"got {transfer_to!r}")
    if isinstance(genes, bool) or not isinstance(genes, int) or genes < 0:
        raise ValueError(f"genes must be a non-negative integer, got {genes!r}")
    if genes and (isinstance(gene_length, bool) or not isinstance(gene_length, int) or gene_length < 1):
        raise ValueError(f"gene_length must be a positive integer, got {gene_length!r}")
    if gff is not None:                              # declared from a GFF: exact coordinates and names
        if genes:
            raise ValueError("pass either gff= or genes=, not both — a GFF already declares the genes")
        lengths, gff_genes = read_gff(gff, trim_overlaps=trim_overlaps)
        seqids = sorted(lengths)                     # a deterministic replicon order
        by_seqid: dict[str, list] = {sq: [] for sq in seqids}
        for gene in gff_genes:
            by_seqid[gene.seqid].append((gene.start, gene.end, gene.strand, gene.name))
        specs = [(_valid_length(lengths[sq]), topology) for sq in seqids]
        layouts = [by_seqid[sq] for sq in seqids]
    else:
        specs = _replicon_specs(chromosomes, root_length, topology)
        for _length, _top in specs:                  # every gene must fit, and leave intergenic room
            if genes and genes * gene_length >= _length:
                raise ValueError(f"{genes} genes of {gene_length} bp do not fit in a {_length} bp "
                                 f"replicon with room left for intergenes")
        layouts = [_even_gene_intervals(length, genes, gene_length) for (length, _t) in specs]
    rates = _Rates(inversion, translocation, transposition, loss, duplication, transfer, origination,
                   fission, fusion, chromosome_origination, chromosome_loss, inversion_length,
                   translocation_length, transposition_length, loss_length, duplication_length,
                   transfer_length, origination_length, inversion_probability)
    depth = mean_root_to_tip(tree)                       # timescale for Distance weighting

    rng = np.random.default_rng(seed)
    chrom_counter = 0
    copy_counter = 0
    source_counter = len(specs)                          # de-novo sources continue past the seed sources

    def new_chrom_id() -> int:
        nonlocal chrom_counter
        cid = chrom_counter
        chrom_counter += 1
        return cid

    def new_copy() -> int:
        nonlocal copy_counter
        copy_counter += 1
        return copy_counter                             # copy ids start at 1 (0 = the unset sentinel)

    family_counter = 0

    def new_family() -> int:
        nonlocal family_counter
        family_counter += 1
        return family_counter                           # gene family ids start at 1 (0 = intergene)

    def new_source() -> int:
        nonlocal source_counter
        src = source_counter
        source_counter += 1
        return src

    genomes: dict[int, NucleotideGenome] = {}
    events: list[Origination | Loss | Duplication | Transfer | Speciation] = []
    rearrangements: list[Inversion | Translocation | Transposition] = []
    chromosome_events: list[ChromosomeEvent] = []
    root = tree.nodes[tree.root]
    schedule = sorted((tree.nodes[i].end_time, i) for i in tree.nodes)   # (end_time, node) in time order

    root_chroms = []
    gene_spans: dict[int, tuple[int, int, int]] = {}
    gene_names: dict[str, int] = {}
    gene_strands: dict[int, int] = {}
    for source, ((length, top), intervals) in enumerate(zip(specs, layouts)):  # one source per replicon
        cid = new_chrom_id()
        cp = new_copy()                                 # ...and one seed copy lineage per replicon
        root_chroms.append(Chromosome(cid, top, _seed_blocks(source, length, cp, intervals, new_family,
                                                             gene_spans, gene_names,
                                                             gene_strands)))
        chromosome_events.append(ChromosomeEvent(root.birth_time, "origination", root.id, (), (cid,)))
        events.append(Origination(root.birth_time, root.id, cid, cp, source, 0, length))

    t = root.birth_time
    alive: list[int] = []                               # the live-lineage set (species._grow shape)
    gen: list[NucleotideGenome] = []
    pos: dict[int, int] = {}
    enter(alive, gen, pos, root.id, NucleotideGenome(root_chroms))
    total_length = sum(c.length for c in root_chroms)
    total_chromosomes = len(root_chroms)

    si = 0
    while si < len(schedule):
        length, count, nlin = total_length, total_chromosomes, len(alive)
        can_xfer = nlin >= 2 or self_transfer
        r_inv = rates.inversion * nlin                  # every extension event is PER LINEAGE:
        r_trl = rates.translocation * nlin              # the rate says how often a lineage does this,
        r_trp = rates.transposition * nlin              # and the extent says how much it touches — so
        r_los = rates.loss * nlin                       # a bigger genome does NOT get more events
        r_dup = rates.duplication * nlin                # (that would double-count size and explode).
        r_tra = rates.transfer * nlin if can_xfer else 0.0
        r_org = rates.origination * nlin
        r_fis = rates.fission * count
        r_fus = rates.fusion * count
        r_cor = rates.chromosome_origination * nlin     # per lineage (de-novo replicon)
        r_clo = rates.chromosome_loss * count           # per chromosome
        total = (r_inv + r_trl + r_trp + r_los + r_dup + r_tra + r_org + r_fis + r_fus + r_cor + r_clo)
        next_species = schedule[si][0]
        if total > 0.0:
            t_ev = t + float(rng.exponential(1.0 / total))
            if t_ev < next_species:                     # a genome event fires before the next species event
                t = t_ev
                r = float(rng.random()) * total
                b_trl = r_inv + r_trl
                b_trp = b_trl + r_trp
                b_los = b_trp + r_los
                b_dup = b_los + r_dup
                b_tra = b_dup + r_tra
                b_org = b_tra + r_org
                b_fis = b_org + r_fis
                b_fus = b_fis + r_fus
                b_cor = b_fus + r_cor
                if r < r_inv:
                    k = int(rng.integers(nlin))
                    _do_inversion(gen[k], alive[k], t, rates.inversion_length, rng, rearrangements)
                elif r < b_trl:
                    k = int(rng.integers(nlin))
                    _do_translocation(gen[k], alive[k], t, rates.translocation_length,
                                      rates.inversion_probability, rng, rearrangements)
                elif r < b_trp:
                    k = int(rng.integers(nlin))
                    _do_transposition(gen[k], alive[k], t, rates.transposition_length,
                                      rates.inversion_probability, rng, rearrangements)
                elif r < b_los:
                    k = int(rng.integers(nlin))
                    total_length += _do_loss(gen[k], alive[k], t, rates.loss_length, rng, events)
                elif r < b_dup:
                    k = int(rng.integers(nlin))
                    total_length += _do_duplication(gen[k], alive[k], t, rates.duplication_length,
                                                    rng, events, new_copy)
                elif r < b_tra:
                    kd = int(rng.integers(nlin))
                    total_length += _do_transfer(rng, tree, alive, gen, kd, t, rates.transfer_length,
                                                 transfer_to, self_transfer, depth, events, new_copy)
                elif r < b_org:
                    k = int(rng.integers(nlin))         # origination is per lineage: a uniform lineage
                    total_length += _do_origination(gen[k], alive[k], t, rates.origination_length,
                                                    rng, events, new_source, new_copy,
                                                    new_family, gene_spans, gene_strands)
                elif r < b_fis:
                    k = _pick_lineage_by_chromosomes(rng, gen, count)
                    total_chromosomes += _do_fission(gen[k], alive[k], t, rng, chromosome_events,
                                                     new_chrom_id)
                elif r < b_fus:
                    k = _pick_lineage_by_chromosomes(rng, gen, count)
                    total_chromosomes += _do_fusion(gen[k], alive[k], t, rng, chromosome_events,
                                                    new_chrom_id)
                elif r < b_cor:
                    k = int(rng.integers(nlin))         # chromosome origination is per lineage
                    total_chromosomes += _do_chromosome_origination(gen[k], alive[k], t,
                                                                    chromosome_events, new_chrom_id)
                else:
                    k = _pick_lineage_by_chromosomes(rng, gen, count)
                    dc, dl = _do_chromosome_loss(gen[k], alive[k], t, rng, events, chromosome_events)
                    total_chromosomes += dc
                    total_length += dl
                continue

        t = next_species                                # advance to the next species event(s)
        while si < len(schedule) and schedule[si][0] == t:   # process the whole tie-batch
            i = schedule[si][1]
            g = gen[pos[i]]
            genomes[i] = g                              # freeze: the lineage retires, never mutated again
            total_length -= g.length
            total_chromosomes -= len(g.chromosomes)
            retire(alive, gen, pos, pos[i])
            node = tree.nodes[i]
            if node.children is not None:              # a speciation: re-mint into the daughters
                for c, cg in _speciate(node, g, new_chrom_id, new_copy, events,
                                       chromosome_events).items():
                    enter(alive, gen, pos, c, cg)
                    total_length += cg.length
                    total_chromosomes += len(cg.chromosomes)
            si += 1
    return NucleotideGenomesResult(tree, genomes, events, rearrangements, chromosome_events, seed,
                                  gene_spans, gene_names, gene_strands)


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
    ``origination`` (a root) / ``duplication`` / ``transfer`` / ``speciation`` (a parent segment ends →
    a child begins on species branch ``lineage``) / ``loss`` (a dead leaf). ``copy`` is the fresh
    segment id, ``parent`` the segment it descends from (``None`` for a root or a loss)."""

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


def _emit_block_events(fam, s, a, b, tree, origs, dups, transfers, losses, specs, new_seg, out) -> None:
    """Replay the copy-lineage log for one root-block ``(s, [a, b))`` into per-segment events on
    ``out``. A copy lineage is a *block-copy* when it covers ``[a, b)`` in full; a duplication or a
    transfer that covers it in full begets a block-copy child (a ladder rung on its parent — the
    transfer's child on the recipient branch, a horizontal edge), a speciation re-mints it, and any
    loss overlapping it ends it."""
    def covers(x, y):
        return x <= a and b <= y

    def overlaps(x, y):
        return x < b and a < y

    root_origs = [e for e in origs if e.source == s and covers(e.start, e.end)]   # 1 per source
    if not root_origs:
        return

    spawns: dict[int, list[tuple[float, int, str]]] = collections.defaultdict(list)  # parent -> [(t, child, kind)]
    species: dict[int, int] = {}                                                     # copy -> species branch
    for e in dups:
        for (pc, cc, src, x, y) in e.copied:
            if src == s and covers(x, y):
                spawns[pc].append((e.time, cc, "duplication"))
                species[cc] = e.lineage                    # the tandem copy stays on the donor branch
    for e in transfers:
        for (pc, cc, src, x, y) in e.transferred:
            if src == s and covers(x, y):
                spawns[pc].append((e.time, cc, "transfer"))
                species[cc] = e.recipient                  # the transferred copy lands on the recipient
    loss_of: dict[int, float] = {}                                                   # copy -> earliest loss time
    for e in losses:
        for (cp, src, x, y) in e.lost:
            if src == s and overlaps(x, y) and (cp not in loss_of or e.time < loss_of[cp]):
                loss_of[cp] = e.time
    for e in root_origs:
        species[e.copy] = e.lineage                    # a de-novo origination roots at its own branch

    block_copies: set[int] = set()
    order: list[int] = []
    stack = [e.copy for e in root_origs]
    while stack:                                            # BFS the block-copy forest (single-parent)
        c = stack.pop()
        if c in block_copies:
            continue
        block_copies.add(c)
        order.append(c)
        for (_t, cc, _kind) in spawns.get(c, ()):
            stack.append(cc)
        if c not in loss_of and c in specs:                # survived to its node's end -> re-minted
            pnode = tree.nodes[specs[c].lineage]
            for i, d in enumerate(specs[c].children):
                species[d] = pnode.children[i]
                stack.append(d)

    seg_in = {c: new_seg() for c in order}
    for e in root_origs:
        out.append(_SegEvent("origination", fam, e.lineage, e.time, seg_in[e.copy], None))
    for c in order:
        prev = seg_in[c]
        for (t, cc, kind) in sorted(spawns.get(c, ())):    # ladder: each rung a bifurcation (dup or transfer)
            nxt = new_seg()
            out.append(_SegEvent(kind, fam, species[c], t, nxt, prev))       # continuation, on c's branch
            out.append(_SegEvent(kind, fam, species[cc], t, seg_in[cc], prev))  # the new copy
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
    """The full recovery: the root partition, and the gene trees.

    With genes declared we build a tree **only** for the root-blocks that are declared genes, keyed by
    gene family id — the intergenic blocks keep their genealogy in the log but are not worth a tree, and
    skipping them is most of the work. With none declared, every root-block is a family (keyed by index).
    Reuses the shared per-segment tree builder."""
    tree = result.complete_tree
    blocks = _root_block_partition(result)
    origs = [e for e in result.events if isinstance(e, Origination)]
    dups = [e for e in result.events if isinstance(e, Duplication)]
    transfers = [e for e in result.events if isinstance(e, Transfer)]
    losses = [e for e in result.events if isinstance(e, Loss)]
    specs = {e.parent: e for e in result.events if isinstance(e, Speciation)}
    counter = [0]

    def new_seg():
        counter[0] += 1
        return counter[0]

    if result.gene_spans:                                # genic: one family per surviving declared gene
        family_of = {span: fam for fam, span in result.gene_spans.items()}
        targets = [(family_of[iv], iv) for iv in blocks if iv in family_of]
    else:                                                # uniform: every root-block is its own family
        targets = list(enumerate(blocks))

    seg_events: list[_SegEvent] = []
    for fam, (s, a, b) in targets:
        _emit_block_events(fam, s, a, b, tree, origs, dups, transfers, losses, specs, new_seg, seg_events)
    return blocks, gene_trees_from_events(seg_events, tree)


__all__ = ["Block", "Chromosome", "NucleotideGenome", "NucleotideGenomesResult",
           "Origination", "Loss", "Duplication", "Transfer", "Speciation", "Inversion",
           "Translocation", "Transposition", "Distance", "simulate_genomes_nucleotide"]
