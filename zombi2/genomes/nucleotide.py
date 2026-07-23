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

Those are all **ancestry-neutral** (the strong invariant: every node carries the whole initial sequence,
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
  (cut each source at the union of **every node's** breakpoints, keep the intervals some node still
  carries: each is a maximal never-cut interval = one gene family) and, per block, the copy-lineage log replayed
  into one gene tree (a duplication is a ladder, a transfer a horizontal edge, a speciation a
  bifurcation), reusing the shared per-segment builder. Validated end to end: the recovered extant tips
  equal the copies actually present in every extant leaf.

Finally, the **genic layer** (``genes`` / ``gene_length``) declares some blocks to be **genes**, the
rest **intergenes**. One rule governs all of it: **a breakpoint may never fall strictly inside a gene.**
That set of legal positions is written down once (``Chromosome._legal_cuts``) and everything reads it —
where an event starts, where its far end lands, where an arc arrives, where a chromosome is cut. A gene
contributes its own edge; an intergene contributes its whole interior. Every breakpoint is drawn
**directly from that set**, never guessed and retried, with the far end weighted by ``exp(-d / mean)``:
the extent distribution you asked for, restricted to the ends that exist. Nothing is clipped, nothing is
snapped (snapping would inflate an extent to "whatever it takes to clear the gene"), and nothing is
silently dropped for running out of retries. So a gene is only ever engulfed **whole**, never split, and
each recovers as exactly one root-block with one gene tree.

A consequence worth stating, because the code got it wrong for a while: a genome with **no spacer at
all** — ten 100 bp genes in 1000 bp — is not frozen. Its gene boundaries are legal cuts like any other,
so it inverts, duplicates, loses and transfers whole genes. The sampling used to say "nucleate in the
spacer", which silently ate every rate when there was none.

**Two consequences the guide must state plainly.** (1) Gene turnover is *emergent* and size-dependent —
a gene changes copy number only when an event engulfs it whole, so large genes are rarely lost or
duplicated. (2) The **realised extent is shorter than the mean you ask for**, the more so the denser the
genome: on a 94%-genic bacterial genome, asking for 3 000 bp yields ~1 300. A stretch ending exactly
where you asked often does not exist, so the arc stops at the nearest legal breakpoint. That
conditioning is the model; what it must not do — and no longer does — is quietly eat the event *rate*
along with it.

Because every node votes on the partition, the recovery covers the **complete** tree, not only its
extant tips: any node's genome can be put back together (``NucleotideGenomesResult.assembly``), an
extinct lineage and the root as readily as a survivor.

Deferred to later slices: homologous *replacement* transfer (only additive for now); indels;
pseudogenization (``gene → intergene``); BED output; and the opt-in per-copy dial (size-blind,
settable per-gene turnover — a *second* selection method, deliberately kept out).
"""

from __future__ import annotations

import bisect
import collections
import math
import pathlib
from dataclasses import dataclass, field

import numpy as np

from ..species import SpeciesResult, Tree
from ._live import enter, retire, without_cyclic_gc
from ._transfer import Distance, mean_root_to_tip, recipient_index
from .chromosomes import ChromosomeEvent, chromosome_events_tsv
from ..progress import progress_bar
from .events import node_from_label, node_label
from .gene_trees import GeneTree, gene_trees_from_events, write_gene_trees
from .gff import read_fasta, read_gff


@dataclass(slots=True)
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


#: sentinel for "no entry", so ``None`` can mean the distinct "the copy died here"
_MISSING = object()


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

    @property
    def n_genes(self) -> int:
        """How many gene copies this chromosome carries. A chromosome may never be left with none:
        a replicon is born with a gene, and an event that would strip its last one does not happen."""
        return sum(1 for b in self.blocks if b.is_gene)

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

    def _legal_cuts(self) -> list[tuple[int, int]]:
        """**The** legal cut set of this chromosome: every physical position where a breakpoint may
        fall, as inclusive ranges ``[(lo, hi), …]`` in position order.

        There is one rule — *a cut may never fall strictly inside a gene* — and this is the one place
        it is written down. Everything that needs a breakpoint reads it: where an event starts, where
        its far end lands, where an arc arrives, where a chromosome is cut. It used to be spelled out
        separately for each, and the copies disagreed: the ones choosing an event said "in the spacer",
        which silently froze a genome with no spacer, while the ones landing one said "not inside a
        gene" and were right.

        So a **gene** contributes its own leading edge — one position, no more — and an **intergene**
        contributes its whole interior. The ranges partition the set, so counting is
        ``hi - lo + 1`` apiece. Both extremes then take care of themselves: an **empty** chromosome
        (a de-novo replicon) still offers position 0, so material can arrive on it, and a **fully
        genic** one still offers every boundary between its genes, so whole genes are moved about
        rather than nothing happening at all."""
        out, pos = [], 0
        for b in self.blocks:
            out.append((pos, pos + b.length - 1) if not b.is_gene else (pos, pos))
            pos += b.length
        if self.topology == "linear":
            out.append((pos, pos))                              # the far end is a legal cut too
        return out

    def _pick_legal_cut(self, rng) -> int:
        """A **uniform** position from :meth:`_legal_cuts`."""
        cuts = self._legal_cuts()
        total = sum(hi - lo + 1 for lo, hi in cuts)
        if not total:
            return 0                                            # an empty ring: position 0 it is
        m = int(rng.integers(total))
        for lo, hi in cuts:
            if m <= hi - lo:
                return lo + m
            m -= hi - lo + 1
        raise AssertionError("legal-cut count out of sync with the blocks")  # unreachable

    def _pick_arc_extent(self, start: int, mean: float, rng) -> int | None:
        """Choose the arc's far end, forward from ``start``: an extent ``d >= 1`` whose breakpoint at
        ``start + d`` is **legal** (never strictly inside a gene), drawn with weight ``exp(-d/mean)``.

        This is the extent distribution you asked for, **restricted to the ends that exist** — sampled
        directly from :meth:`_legal_cuts` rather than drawn and rejected, so nothing is ever wasted and
        no event silently vanishes because its retries ran out. On a gene-dense genome the realised
        extent therefore comes out *shorter* than ``mean``: a long stretch ending exactly where you
        asked often simply does not exist, so the arc stops at the nearest legal breakpoint instead.
        That conditioning is the model; what it must not do is quietly eat the event *rate* as well.

        ``None`` when no legal end exists."""
        total = self.length
        if total < 2:
            return None
        circular = self.topology == "circular"
        limit = total - 1 if circular else total - start     # the longest arc that still fits
        if limit < 1:
            return None
        spans = self._legal_cuts()
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

    # --- the mutators, by explicit coordinates ----------------------------------------------------
    # Each takes the arc it acts on rather than drawing it, so the ``_do_*`` events above can pick
    # with the rng and then call one of these, and a test (or a replay of a written event log) can
    # call the same one with coordinates in hand. There is one implementation, not two: the engine
    # runs exactly the code a scripted event runs.

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

    def duplicate(self, start: int, length: int, new_copy) -> tuple | None:
        """Copy the arc ``[start, start+length)`` **in tandem**, the copy landing immediately after it.

        Each distinct copy lineage in the arc begets one fresh child (minted by ``new_copy``), so the
        tandem copy is new material of the same ancestry. Returns the ``(parent copy, child copy,
        source, start, end)`` record per block — what a :class:`Duplication` carries — or ``None`` if
        the arc is empty."""
        span = self._arc_range(start, length)
        if span is None:
            return None
        i, j = span
        arc = self.blocks[i:j]
        child_of: dict[int, int] = {}
        for b in arc:
            if b.copy not in child_of:
                child_of[b.copy] = new_copy()
        copied = tuple((b.copy, child_of[b.copy], b.source, b.start, b.end) for b in arc)
        self.blocks[j:j] = [Block(b.source, b.start, b.end, b.strand, child_of[b.copy], b.gene)
                            for b in arc]
        return copied

    def delete(self, start: int, length: int) -> tuple | None:
        """Remove the arc ``[start, start+length)``. Returns the ``(copy, source, start, end)`` record
        per block removed — what a :class:`Loss` carries — or ``None`` when the deletion does not
        happen: an empty arc, a chromosome of under two nucleotides, or a cut that would strip the
        chromosome of its last gene (a chromosome never exists without one)."""
        if self.length < 2:
            return None
        span = self._arc_range(start, length)
        if span is None:
            return None
        i, j = span
        gone = self.blocks[i:j]
        if self.n_genes and sum(1 for b in gone if b.is_gene) == self.n_genes:
            return None
        lost = tuple((b.copy, b.source, b.start, b.end) for b in gone)
        self.blocks = self.blocks[:i] + self.blocks[j:]
        return lost

    def originate(self, at: int, length: int, source: int, copy: int, family: int) -> None:
        """Lay down a new ``length``-nucleotide **gene** of its own fresh ``source`` at position
        ``at``. Indivisible from birth, like a declared gene. Raises :class:`_CutsGene` if ``at``
        falls inside an existing gene."""
        self._split_at(at)
        k = self._index_at(at)
        self.blocks[k:k] = [Block(source, 0, length, 1, copy, family)]

    def excise(self, start: int, length: int) -> list | None:
        """Lift the arc ``[start, start+length)`` out of the chromosome and return its blocks (or
        ``None`` for an empty arc). The other half of a transposition — see :meth:`place`."""
        span = self._arc_range(start, length)
        if span is None:
            return None
        i, j = span
        arc = self.blocks[i:j]
        self.blocks = self.blocks[:i] + self.blocks[j:]
        return arc

    def place(self, arc: list, at: int, flipped: bool = False) -> None:
        """Insert ``arc`` at position ``at``, reversed and strand-flipped if ``flipped``.

        ``at`` is a position in the chromosome **as it stands now** — for a transposition, that is
        the genome *after* the arc was excised, which is shorter. Raises :class:`_CutsGene` if ``at``
        falls inside a gene."""
        if flipped:
            arc = [Block(b.source, b.start, b.end, -b.strand, b.copy, b.gene) for b in reversed(arc)]
        self._split_at(at)
        k = self._index_at(at)
        self.blocks[k:k] = arc

    def transpose(self, start: int, length: int, dest: int, flipped: bool = False) -> bool:
        """Move the arc ``[start, start+length)`` to ``dest``, flipped or not — :meth:`excise` then
        :meth:`place`, with the chromosome restored if the landing is not a legal cut.

        ``dest`` is a position in the **remainder**, after the arc is lifted out (the engine picks it
        that way round, so a scripted call means the same thing an engine-drawn one does). Returns
        whether the move happened."""
        intact = self.blocks
        arc = self.excise(start, length)
        if arc is None:
            return False
        try:
            self.place(arc, dest, flipped)
        except _CutsGene:
            self.blocks = intact                        # nowhere legal to land: undo the excision
            raise
        return True

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

    def _pick_legal_cut(self, rng) -> tuple[Chromosome, int] | None:
        """A uniform pick over the **whole genome's** legal breakpoints → ``(chromosome, physical
        position)``. :meth:`Chromosome._legal_cuts` one scope up: it is where an event *starts*, as
        against where one lands, and both are the same set of positions.

        With no genes declared every position is legal and this is a plain uniform pick. ``None`` only
        when there is nowhere at all — a genome with no chromosomes."""
        cuts = [c._legal_cuts() for c in self.chromosomes]
        counts = [sum(hi - lo + 1 for lo, hi in cc) for cc in cuts]
        total = sum(counts)
        if total == 0:
            return None
        m = int(rng.integers(total))
        for chrom, cc, n in zip(self.chromosomes, cuts, counts):
            if m < n:
                for lo, hi in cc:
                    if m <= hi - lo:
                        return chrom, lo + m
                    m -= hi - lo + 1
            m -= n
        raise AssertionError("legal-cut count out of sync with the blocks")  # unreachable

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

@dataclass(frozen=True, slots=True)
class Inversion:
    """A recorded nucleotide inversion: on species branch ``lineage`` at ``time``, the arc
    ``[start, start+length)`` of chromosome ``chromosome`` was reversed. Ancestry is unchanged, so it
    is a rearrangement record, not a gene-genealogy event."""

    time: float
    lineage: int
    chromosome: int
    start: int
    length: int


@dataclass(frozen=True, slots=True)
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


@dataclass(frozen=True, slots=True)
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


@dataclass(frozen=True, slots=True)
class Origination:
    """A recorded **birth of a copy lineage** — the root of a gene tree. At ``time`` on branch
    ``lineage`` new material was laid down on chromosome ``chromosome`` as copy lineage ``copy``,
    covering the ancestral interval ``[start, end)`` on ``source``.

    ``initial`` tells the two roots apart. An **initial** origination (``initial=True``) lays down the
    initial genome at time 0 — one per initial replicon — and is what the run *starts* with, not
    something it *did*; it is written with kind ``"initial"`` so counting ``"origination"`` in the log
    gives the de-novo births alone (what the ``origination`` rate controls). A **de-novo** origination
    (``initial=False``) is a fresh source arising mid-tree. The gene-tree recovery reads either as the
    root of its family."""

    time: float
    lineage: int
    chromosome: int
    copy: int
    source: int
    start: int
    end: int
    initial: bool = False


@dataclass(frozen=True, slots=True)
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


@dataclass(frozen=True, slots=True)
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


@dataclass(frozen=True, slots=True)
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


@dataclass(frozen=True, slots=True)
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
    #: The genome the run **started** with — the initial karyotype, laid down at the start of the root branch,
    #: before any event. It is not in :attr:`genomes`, which holds a genome per *node*, and a node
    #: sits at the **end** of its branch: the root branch is real simulated time, so
    #: ``genomes[root]`` is this genome plus whatever happened along the stem. It votes on the root
    #: partition like every other genome, so it can be reconstructed too — see
    #: :meth:`initial_assembly`.
    initial_genome: NucleotideGenome = field(default_factory=lambda: NucleotideGenome([]))
    #: ``{source: DNA}`` — the **initial sequence** the run was given, one entry per initial replicon,
    #: from ``fasta=`` paired with ``gff=``. Empty when no FASTA was given (then the sequence level
    #: draws the founding sequence from the model instead). A *de-novo* originated source is never
    #: here: it arose mid-run, so nothing was supplied for it. The letters live here, not in
    #: :attr:`genomes` (which is pure ancestry) — the sequence level reads them as each block's
    #: founding sequence, and an assembled genome then descends from exactly this input.
    initial_sequence: dict[int, str] = field(default_factory=dict)

    def mosaic(self, node_id: int) -> dict[int, list[tuple[int, int, int, int]]]:
        return self.genomes[node_id].mosaic()

    def trace_back(self, node_id: int) -> dict[int, list[tuple[int, int, int]]]:
        return self.genomes[node_id].trace_back()

    def ancestry(self, node_id: int) -> list[tuple[int, int]]:
        return self.genomes[node_id].ancestry()

    def _recover(self):
        if not hasattr(self, "_recovered"):
            self._recovered = _recover_gene_trees(self)
        return self._recovered

    def _recover_blocks(self):
        """The every-block recovery, cached: :attr:`block_trees` and :meth:`assembly` are two reads of
        the same replay, and it is far too expensive to run twice."""
        if not hasattr(self, "_recovered_blocks"):
            self._recovered_blocks = _recover_gene_trees(self, every_block=True)
        return self._recovered_blocks

    @property
    def root_blocks(self) -> list[tuple[int, int, int]]:
        """The recovered **root partition**: ``(source, start, end)`` for each maximal never-cut
        interval that some node still carries — one per :attr:`block_trees` entry (by index).

        Cut at the breakpoints of **every** node's genome, not only the extant leaves', which is what
        lets any node be reconstructed (:meth:`assembly`) rather than the survivors alone."""
        return self._recover()[0]

    @property
    def block_trees(self) -> dict[int, GeneTree]:
        """``{root-block index: GeneTree}`` — a tree for **every** recovered root block, spacer as
        well as genes, keyed by its index in :attr:`root_blocks`.

        :attr:`gene_trees` covers the declared genes; this covers the whole genome. A block never
        splits, so its size is fixed and its genealogy is in the event log just as a gene's is — the
        recovery is the same one, pointed at every block instead of a chosen few. That is what makes
        an ancestral genome reconstructable at any node rather than only at the loci you declared.

        A gene's tree here has the **same topology and branch lengths** as its :attr:`gene_trees` one,
        but not the same ``g<id>`` leaf labels: segment ids are handed out as the recovery walks its
        targets, and walking every block numbers them differently from walking three. Use one accessor
        or the other within a piece of analysis — they are the same genealogy under different names."""
        return self._recover_blocks()[1]

    def block_of(self, family: int) -> int:
        """The index in :attr:`root_blocks` of the block a declared **gene family** occupies — the join
        between the two numbering schemes this resolution has.

        They are both plain ints over overlapping ranges, so mixing them up is silent: ``gene_spans``
        and ``gene_trees`` are keyed by **gene family id**, while ``root_blocks``, ``block_trees`` and
        everything a sequence run produces here are keyed by **block index** (every block evolves, and
        spacer has no family). ``block_of`` is how you get from a gene to its sequences::

            r.alignments[g.block_of(g.gene_names["dnaA"])]     # that gene's alignment

        Raises ``KeyError`` for a family that was never declared, and ``LookupError`` for one declared
        but surviving nowhere at all — it has no recovered block, so there is nothing to point at.
        The reverse lookup is one line: ``{span: fam for fam, span in g.gene_spans.items()}`` read at
        ``root_blocks[i]``."""
        span = self.gene_spans[family]                    # KeyError: never declared
        try:
            return self.root_blocks.index(span)
        except ValueError:
            raise LookupError(
                f"gene family {family} spans {span} but has no recovered root block — no node in the "
                "tree still carries it, so nothing was reconstructed for it. gene_trees leaves such a "
                "family out for the same reason.") from None

    def assembly(self, node_id: int) -> dict[int, list[tuple[int, int, int]]]:
        """How this node's genome is built out of the recovered root blocks:
        ``{chromosome id: [(block, gene, strand), …]}`` in **physical order**, where ``block`` indexes
        :attr:`root_blocks`, ``gene`` is the gene id that block's tree gives this node's copy (the
        ``g<id>`` label in :attr:`block_trees`), and ``strand`` is ``+1`` read forward or ``-1``
        reverse-complemented.

        To reconstruct a genome: pair each piece with its block's evolved sequence, flip the
        ``-1``\\ s, and concatenate. The sequence level does exactly that; nothing here knows about
        letters. **Every** node works — an extinct leaf and the root as readily as a surviving tip —
        which is what makes the whole history recoverable rather than only its leaves.
        :meth:`initial_assembly` does the same for the genome the run started with.

        A piece is always a **whole** block, never part of one, because every node votes on where the
        partition is cut (see :func:`_root_block_partition`): this node's own breakpoints are all in
        it, so each of its blocks is a whole number of root blocks. What a block *is* cut into is one
        piece per root block it spans — and on a reversed block those come out in descending
        coordinate order, since physical order runs *down* the source."""
        tips = self._recover_blocks()[2]
        blocks = self.root_blocks
        what = node_label(node_id)
        out = {}
        for cid, pieces in self._pieces(self.genomes[node_id], what).items():
            named = []
            for (i, copy, strand) in pieces:
                gene = tips.get((i, copy), _MISSING)
                if gene is _MISSING or gene is None:
                    raise AssertionError(                            # a guard — see the class docstring
                        f"{what} carries {blocks[i]} under copy lineage {copy}, but that block's "
                        + ("genealogy has no such copy" if gene is _MISSING
                           else "genealogy ends that copy in a loss")
                        + " — the event log and the genomes disagree")
                named.append((i, gene, strand))
            out[cid] = named
        return out

    def initial_assembly(self) -> dict[int, list[tuple[int, int]]]:
        """:meth:`assembly` for :attr:`initial_genome`: ``{chromosome id: [(block, strand), …]}``.

        No gene id here, unlike :meth:`assembly`, and that is the honest shape rather than a saving.
        The initial genome sits at the **start** of the root branch, before any event, so each of its
        blocks has exactly one sequence — the founding draw the sequence level records as
        ``founding[block]`` — and there is no copy to disambiguate. A gene id would in fact be *wrong*
        here: the one :meth:`assembly` gives is the **last** gene a copy held, and for an initial copy
        that is at the far end of the stem. A loss on the stem can even end it, which is the same
        thing said louder."""
        return {cid: [(i, strand) for (i, _copy, strand) in pieces]
                for cid, pieces in self._pieces(self.initial_genome, "the initial genome").items()}

    def _pieces(self, genome: NucleotideGenome, what: str
                ) -> dict[int, list[tuple[int, int, int]]]:
        """The walk both assemblies share: ``{chromosome id: [(block, copy, strand), …]}`` in physical
        order. Cutting each of the genome's blocks at the partition, which is at least as fine."""
        blocks = self.root_blocks
        index = self._block_index()
        out: dict[int, list[tuple[int, int, int]]] = {}
        for chrom in genome.chromosomes:
            pieces: list[tuple[int, int, int]] = []
            for b in chrom.blocks:
                cut = []
                starts, idx = index.get(b.source, ([], []))
                k = bisect.bisect_left(starts, b.start)   # the partition starts a block exactly here
                at = b.start
                while k < len(idx) and blocks[idx[k]][1] == at < b.end:
                    cut.append((idx[k], b.copy, b.strand))
                    at = blocks[idx[k]][2]
                    k += 1
                if at != b.end:
                    raise AssertionError(                            # a guard — see the class docstring
                        f"{what} carries [{at}, {b.end}) of source {b.source}, which the root "
                        "partition does not cover, though every genome in the run votes on it")
                pieces.extend(reversed(cut) if b.strand == -1 else cut)
            out[chrom.id] = pieces
        return out

    def _block_index(self) -> dict[int, tuple[list[int], list[int]]]:
        """``{source: ([block start, …], [block index, …])}`` — the root partition indexed by source
        for lookup. It comes back sorted, so each source's starts are already ascending."""
        if not hasattr(self, "_block_ix"):
            ix: dict[int, tuple[list[int], list[int]]] = {}
            for i, (src, a, _b) in enumerate(self.root_blocks):
                starts, idx = ix.setdefault(src, ([], []))
                starts.append(a)
                idx.append(i)
            self._block_ix = ix
        return self._block_ix

    @property
    def gene_trees(self) -> dict[int, GeneTree]:
        """``{family: GeneTree}`` — the recovered gene trees.

        With **genes declared**, one tree per gene, keyed by its **gene family id** (see
        :attr:`gene_spans`); the intergenic root-blocks keep their block ancestry in the log but are not
        built into trees. With **no genes declared** the whole genome is one big intergene, so every
        recovered root-block is a family in its own right and the key is its index in
        :attr:`root_blocks`. Every node votes on the partition, so a gene surviving only in lineages
        that died still gets a tree — a complete one, with no extant tree to go with it. Only a gene
        lost from *every* node has no root-block and no tree."""
        return self._recover()[1]

    def write(self, directory,
              outputs=("events", "genes", "blocks", "initial_genome", "initial_sequence",
                       "gene_trees", "chromosome_events", "gff", "bed")) -> None:
        """Materialise chosen ``outputs`` to ``directory`` (created if needed):

        - ``"events"`` → ``genome_events.tsv``, the run's whole history in one time-ordered table:
          the copy-lineage genealogy and the ancestry-neutral rearrangements. One row per
          **ancestral interval** an event touched, so an event that spanned several blocks writes
          several rows sharing a ``time`` and ``kind``.
        - ``"blocks"`` → ``blocks.tsv``, every node's genome as its block mosaic (ancestors
          included, as for the ordered resolution's ``gene_order``). The one big file here: blocks
          are not kept maximal during a run, so a rearrangement-heavy genome carries far more of
          them than it has distinct ancestral runs, and this grows with their number × every node.
        - ``"genes"`` → ``genes.tsv``, the declared genes and where they sit in root coordinates.
          Header-only for a run that declared none.
        - ``"initial_genome"`` → ``initial_genome.tsv``, the block mosaic the run started with. Its
          own file, not a row in ``blocks.tsv``, because it belongs to no node: it sits at the start
          of the root branch, and every ``lineage`` in that table is a node at the end of one.
        - ``"chromosome_events"`` → ``chromosome_events.tsv``, the chromosome network's edges. The
          one log kept apart: it is a network over chromosome **ids**, with list-valued parents and
          children, joined on a different key from everything above.
        - ``"gene_trees"`` → ``gene_tree_fam<family>_{complete,extant}.nwk``, one recovered
          genealogy per family some node still carries; the ``_extant`` file only where the family
          has a surviving copy.
        - ``"initial_sequence"`` → ``initial_sequence.fasta``, the initial DNA the run was given (``fasta=``),
          one ``>source<n>`` record per replicon. Written only when a FASTA was supplied — it is what
          lets a separate ``zombi2 sequences`` run found its blocks from the real sequence.
        - ``"gff"`` → ``genome_<lineage>.gff``, that genome's **genes**, in its own coordinates: the
          annotation to read beside the sequence level's ``genome_<lineage>.fasta``.
        - ``"bed"`` → ``genome_<lineage>.bed``, that genome's **blocks** — every piece, spacer
          included, named by the ancestral interval it descends from. The ancestry as a browser track.

        Both name their sequences ``<lineage>_chr<c>``, exactly as the FASTA records are named, so a
        genome and its annotation join without renaming anything. Written for every node and for the
        initial genome, so there are two files per genome.
        """
        d = pathlib.Path(directory)
        d.mkdir(parents=True, exist_ok=True)
        if "events" in outputs:
            (d / "genome_events.tsv").write_text(
                _nucleotide_events_tsv(self.events, self.rearrangements))
        if "blocks" in outputs:
            (d / "blocks.tsv").write_text(self._blocks_tsv())
        if "genes" in outputs:
            (d / "genes.tsv").write_text(self._genes_tsv())
        if "initial_genome" in outputs:
            (d / "initial_genome.tsv").write_text(self._initial_genome_tsv())
        if "chromosome_events" in outputs:
            (d / "chromosome_events.tsv").write_text(chromosome_events_tsv(self.chromosome_events))
        if "gene_trees" in outputs:
            write_gene_trees(self.gene_trees, d)
        if "initial_sequence" in outputs and self.initial_sequence:
            (d / "initial_sequence.fasta").write_text(
                "".join(f">source{src}\n{self.initial_sequence[src]}\n"
                        for src in sorted(self.initial_sequence)))
        for token, ext, render in (("gff", "gff", self._gff), ("bed", "bed", self._bed)):
            if token in outputs:
                for label, genome in self._every_genome():
                    (d / f"genome_{label}.{ext}").write_text(render(label, genome))

    def _every_genome(self):
        """``(label, genome)`` for every genome the run holds — each node, and the initial one. The
        labels are the ones the sequence level writes its FASTA under, so the files pair up."""
        for node_id in sorted(self.genomes):
            yield node_label(node_id), self.genomes[node_id]
        yield "initial", self.initial_genome

    def _laid_out(self, genome: NucleotideGenome):
        """Walk one genome's blocks in physical order, yielding ``(chromosome, at, block)`` — ``at``
        being the block's 0-based offset along its chromosome. The one place that walk is written."""
        for chrom in genome.chromosomes:
            at = 0
            for b in chrom.blocks:
                yield chrom, at, b
                at += b.length

    def _gff(self, label: str, genome: NucleotideGenome) -> str:
        """One genome's **genes** as GFF3, in that genome's own coordinates.

        GFF is 1-based inclusive and blocks are 0-based half-open, so a gene at ``[at, at+len)``
        is written ``at+1 .. at+len``. The strand is the one the gene reads on **here**: its coding
        strand from the annotation, flipped if the block carrying it has been inverted since the root.
        That product is the point of the file — a gene that an inversion turned over says so."""
        name_of = {fam: name for name, fam in self.gene_names.items()}
        out = ["##gff-version 3"]
        for chrom in genome.chromosomes:
            out.append(f"##sequence-region {label}_chr{chrom.id} 1 {chrom.length}")
        for n, (chrom, at, b) in enumerate(
                (x for x in self._laid_out(genome) if x[2].is_gene), start=1):
            strand = "+" if self.gene_strands.get(b.gene, 1) * b.strand == 1 else "-"
            # ID is a plain per-file counter, because GFF wants one unique handle and nothing else
            # here is: `copy` is the *block's* copy lineage, and a whole replicon is laid down as one, so
            # every gene on it shares it until something duplicates them apart. What joins this row to
            # the rest of the run is `family` — its gene tree, and its alignment via block_of(family).
            named = f"Name={name_of[b.gene]};" if b.gene in name_of else ""
            attrs = (f"ID=gene{n};{named}family={b.gene};copy={b.copy};"
                     f"source={b.source}:{b.start}-{b.end}")
            out.append(f"{label}_chr{chrom.id}\tZOMBI2\tgene\t{at + 1}\t{at + b.length}\t.\t"
                       f"{strand}\t.\t{attrs}")
        return "\n".join(out) + "\n"

    def _bed(self, label: str, genome: NucleotideGenome) -> str:
        """One genome's **blocks** as BED — every piece, spacer included. BED is 0-based half-open,
        which is what a block already is, so the coordinates go out unchanged.

        The name column is the **ancestral interval** the block descends from, which is the thing this
        format is worth having for: laid over the genome it says which stretch came from where, and
        two genomes' tracks line up on shared ancestry. ``strand`` is the block's orientation relative
        to that ancestor, not a gene's coding strand (the GFF carries that)."""
        out = []
        for chrom, at, b in self._laid_out(genome):
            out.append(f"{label}_chr{chrom.id}\t{at}\t{at + b.length}\t"
                       f"{b.source}:{b.start}-{b.end}\t0\t{'+' if b.strand == 1 else '-'}")
        return "\n".join(out) + "\n"

    def _blocks_tsv(self) -> str:
        """Every node's genome, block by block. ``position`` is the block's physical offset along its
        chromosome, so the rows of one chromosome tile it end to end from 0."""
        cols = ("lineage", "chromosome", "position", "source", "start", "end", "strand", "copy", "gene")
        rows = []
        for s in sorted(self.genomes):
            for c in self.genomes[s].chromosomes:
                at = 0
                for b in c.blocks:
                    rows.append(f"{node_label(s)}\t{c.id}\t{at}\t{b.source}\t{b.start}\t{b.end}\t{b.strand}\t"
                                f"{b.copy}\t{b.gene}")
                    at += b.length
        return "\n".join(["\t".join(cols), *rows]) + "\n"

    def _initial_genome_tsv(self) -> str:
        """The mosaic the run started with — ``blocks.tsv``'s columns without ``lineage``, which is
        the whole point: it belongs to the start of the root branch, not to a node."""
        cols = ("chromosome", "position", "source", "start", "end", "strand", "copy", "gene")
        rows = []
        for c in self.initial_genome.chromosomes:
            at = 0
            for b in c.blocks:
                rows.append(f"{c.id}\t{at}\t{b.source}\t{b.start}\t{b.end}\t{b.strand}\t"
                            f"{b.copy}\t{b.gene}")
                at += b.length
        return "\n".join(["\t".join(cols), *rows]) + "\n"

    def _genes_tsv(self) -> str:
        """The declared genes: family id, name (blank when unnamed), the root interval the gene
        occupies, and the **coding** strand the GFF gave it."""
        cols = ("family", "name", "source", "start", "end", "strand")
        name_of = {fam: name for name, fam in self.gene_names.items()}
        rows = [f"{fam}\t{name_of.get(fam, '')}\t{src}\t{start}\t{end}\t{self.gene_strands.get(fam, 1)}"
                for fam, (src, start, end) in sorted(self.gene_spans.items())]
        return "\n".join(["\t".join(cols), *rows]) + "\n"


#: One table for the whole history of a run: the copy-lineage genealogy and the ancestry-neutral
#: rearrangements, interleaved by time. ``source`` / ``start`` / ``end`` are **ancestral** coordinates
#: — which stretch of which source — while ``position`` / ``length`` are **physical** ones on
#: the chromosome named by ``chromosome``, as ``blocks.tsv`` numbers it. They are different frames,
#: so they are different columns.
_NUCLEOTIDE_EVENT_COLS = ("time", "kind", "lineage", "chromosome", "copy", "parent", "recipient",
                          "source", "start", "end",
                          "position", "length", "dest_chromosome", "dest_position", "flipped")


def _nucleotide_events_tsv(events, rearrangements=()) -> str:
    """The run's whole history as one time-ordered table (see :data:`_NUCLEOTIDE_EVENT_COLS`).

    An event here can span several blocks at once (a loss deletes an arc covering many), and each
    carries its own copy lineage and ancestral interval, so a flat table needs one row apiece. Rows
    of the same event share ``time``, ``kind`` and ``lineage``. Empty cells for the fields a kind
    does not use; a speciation re-mints a copy lineage without touching sequence, so it names only
    ``parent`` and ``copy``.

    The **rearrangements** are here too, as their own kinds. They end no gene lineage, which is why
    they used to be a separate file, but they are events on the same branches at the same clock and a
    reader replaying one has to interleave them anyway.
    """
    rows = []
    # the columns holding a species-tree node, labelled n<id> like every other table
    node_at = {i for i, c in enumerate(_NUCLEOTIDE_EVENT_COLS) if c in ("lineage", "recipient")}

    def row(*cells):
        cells = cells + ("",) * (len(_NUCLEOTIDE_EVENT_COLS) - len(cells))
        rows.append((cells[0], "\t".join(
            "" if c is None else (node_label(c) if i in node_at else str(c))
            for i, c in enumerate(cells))))

    for e in events:
        if isinstance(e, Origination):
            row(e.time, "initial" if e.initial else "origination", e.lineage, e.chromosome, e.copy,
                None, None, e.source, e.start, e.end)
        elif isinstance(e, Loss):
            for (copy, source, start, end) in e.lost:
                row(e.time, "loss", e.lineage, e.chromosome, copy, None, None, source, start, end)
        elif isinstance(e, Duplication):
            for (parent, child, source, start, end) in e.copied:
                row(e.time, "duplication", e.lineage, e.chromosome, child, parent, None,
                    source, start, end)
        elif isinstance(e, Transfer):
            for (parent, child, source, start, end) in e.transferred:
                row(e.time, "transfer", e.lineage, None, child, parent, e.recipient,
                    source, start, end)
        elif isinstance(e, Speciation):
            for child in e.children:
                row(e.time, "speciation", e.lineage, None, child, e.parent, None, None, None, None)
        else:
            raise AssertionError(f"unhandled event {type(e).__name__}")
    for r in rearrangements:                      # ancestry-neutral: no copy, no ancestral interval
        if isinstance(r, Inversion):
            row(r.time, "inversion", r.lineage, r.chromosome, None, None, None, None, None, None,
                r.start, r.length, None, None, None)
        elif isinstance(r, Transposition):
            row(r.time, "transposition", r.lineage, r.chromosome, None, None, None, None, None, None,
                r.start, r.length, None, r.dest, int(r.flipped))
        else:
            row(r.time, "translocation", r.lineage, r.source, None, None, None, None, None, None,
                r.start, r.length, r.dest, None, int(r.flipped))
    rows.sort(key=lambda tr: tr[0])               # one stream, in the order it happened
    return "\n".join(["\t".join(_NUCLEOTIDE_EVENT_COLS), *[r for _t, r in rows]]) + "\n"


def _rows(text: str, cols: tuple[str, ...], what: str):
    """Parse a TSV this module wrote, checking the header. Yields one cell-list per row."""
    lines = text.splitlines()
    if not lines:
        raise ValueError(f"empty {what} — is the file empty?")
    header = tuple(lines[0].split("\t"))
    if header != cols:
        raise ValueError(f"unexpected {what} columns {list(header)}; expected {list(cols)}")
    for lineno, raw in enumerate(lines[1:], 2):
        if not raw:                                      # tolerate a trailing blank line
            continue
        cells = raw.split("\t")
        if len(cells) != len(cols):
            raise ValueError(f"{what} line {lineno}: expected {len(cols)} columns, got {len(cells)}")
        yield cells


def _karyotype(rows) -> NucleotideGenome:
    """Group ``(chromosome id, topology-less block fields)`` rows into one karyotype, in file order."""
    chroms: dict[int, list[Block]] = {}
    for (cid, source, start, end, strand, copy, gene) in rows:
        chroms.setdefault(cid, []).append(Block(source, start, end, strand, copy, gene))
    # Topology is not in blocks.tsv and the recovery never reads it — every event that needs it has
    # already happened. "circular" is the honest placeholder; a replayed run is for reading, not
    # for evolving further.
    return NucleotideGenome([Chromosome(cid, "circular", blocks) for cid, blocks in chroms.items()])


def _blocks_from_tsv(text: str) -> dict[int, NucleotideGenome]:
    """``blocks.tsv`` → ``{node id: NucleotideGenome}``. ``position`` is redundant (the rows tile each
    chromosome in order) and is checked rather than used, so a hand-edited file fails loudly."""
    per_node: dict[int, list] = {}
    at: dict[tuple[int, int], int] = {}
    cols = ("lineage", "chromosome", "position", "source", "start", "end", "strand", "copy", "gene")
    for cells in _rows(text, cols, "blocks.tsv"):
        node = node_from_label(cells[0])
        cid, position = int(cells[1]), int(cells[2])
        source, start, end, strand, copy, gene = (int(c) for c in cells[3:])
        if at.get((node, cid), 0) != position:
            raise ValueError(f"blocks.tsv: chromosome {cid} of {cells[0]} does not tile — expected "
                             f"position {at.get((node, cid), 0)}, got {position}")
        at[(node, cid)] = position + (end - start)
        per_node.setdefault(node, []).append((cid, source, start, end, strand, copy, gene))
    return {node: _karyotype(rows) for node, rows in per_node.items()}


def _initial_genome_from_tsv(text: str) -> NucleotideGenome:
    """``initial_genome.tsv`` → the genome the run started with. No ``lineage`` column: it is not a
    node's (see :attr:`NucleotideGenomesResult.initial_genome`)."""
    cols = ("chromosome", "position", "source", "start", "end", "strand", "copy", "gene")
    return _karyotype([(int(c[0]), *(int(x) for x in c[2:])) for c in
                       _rows(text, cols, "initial_genome.tsv")])


def _genes_from_tsv(text: str):
    """``genes.tsv`` → ``(gene_spans, gene_names, gene_strands)``."""
    spans, names, strands = {}, {}, {}
    for (fam, name, source, start, end, strand) in _rows(
            text, ("family", "name", "source", "start", "end", "strand"), "genes.tsv"):
        f = int(fam)
        spans[f] = (int(source), int(start), int(end))
        strands[f] = int(strand)
        if name:
            names[name] = f
    return spans, names, strands


def _events_from_tsv(text: str) -> tuple[list, list]:
    """The nucleotide ``genome_events.tsv`` → ``(genealogy, rearrangements)``, the inverse of
    :func:`_nucleotide_events_tsv`.

    An event that spanned several ancestral intervals was written as several rows, so the rows have to
    be regrouped. What identifies one event differs by kind, and getting it wrong merges two events or
    splits one: a **speciation** writes one row per daughter and every copy re-minted at a node shares
    its time and lineage, so the *parent copy* is what tells two apart; a **loss** / **duplication** /
    **transfer** writes one row per interval, and those rows differ in parent copy but share the
    event's time, lineage, chromosome and recipient; an **origination** is always a single row."""
    def num(cell):
        return int(cell) if cell else None

    events: list = []
    rearrangements: list = []
    pending: list = []
    key = None

    def flush():
        nonlocal key
        if not pending:
            return
        kind, time, lineage, chrom, recipient = pending[0][:5]
        if kind in ("origination", "initial"):
            (_k, _t, _l, _c, _r, copy, _p, src, start, end) = pending[0]
            events.append(Origination(time, lineage, chrom, copy, src, start, end, initial=kind == "initial"))
        elif kind == "loss":
            events.append(Loss(time, lineage, chrom,
                               tuple((c, s, a, b) for (*_h, c, _p, s, a, b) in pending)))
        elif kind == "duplication":
            events.append(Duplication(time, lineage, chrom,
                                      tuple((p, c, s, a, b) for (*_h, c, p, s, a, b) in pending)))
        elif kind == "transfer":
            events.append(Transfer(time, lineage, recipient,
                                   tuple((p, c, s, a, b) for (*_h, c, p, s, a, b) in pending)))
        elif kind == "speciation":
            events.append(Speciation(time, lineage, pending[0][6],
                                     tuple(row[5] for row in pending)))
        else:
            raise ValueError(f"genome_events.tsv: unknown event kind {kind!r}")
        pending.clear()
        key = None

    for cells in _rows(text, _NUCLEOTIDE_EVENT_COLS, "genome_events.tsv"):
        (time, kind, lineage, chrom, copy, parent, recipient, source, start, end,
         *_physical) = cells
        if kind not in ("origination", "initial", "loss", "duplication", "transfer", "speciation"):
            flush()                                  # a rearrangement: it ends no copy lineage
            t, ln = float(time), node_from_label(lineage)
            at, ell, dc, dp, fl = (num(c) for c in _physical)
            if kind == "inversion":
                rearrangements.append(Inversion(t, ln, num(chrom), at, ell))
            elif kind == "transposition":
                rearrangements.append(Transposition(t, ln, num(chrom), at, ell, dp, bool(fl)))
            elif kind == "translocation":
                rearrangements.append(Translocation(t, ln, num(chrom), dc, at, ell, bool(fl)))
            else:
                raise ValueError(f"genome_events.tsv: unknown event kind {kind!r}")
            continue
        row = (kind, float(time), node_from_label(lineage), num(chrom),
               node_from_label(recipient) if recipient else None,
               num(copy), num(parent), num(source), num(start), num(end))
        # what makes this row part of the *same* event as the last one
        row_key = (None if kind in ("origination", "initial") else
                   (*row[:3], row[6]) if kind == "speciation" else row[:5])
        if pending and row_key != key:
            flush()
        pending.append(row)
        key = row_key
        if kind in ("origination", "initial"):
            flush()
    flush()
    return events, rearrangements


def read_nucleotide_genomes(directory, tree) -> NucleotideGenomesResult:
    """Rebuild a :class:`NucleotideGenomesResult` from the files a run wrote, so a later level can
    replay it from disk. ``tree`` is the species tree it ran on.

    Reads ``blocks.tsv``, ``initial_genome.tsv``, ``genome_events.tsv`` and ``genes.tsv`` — the four
    the recovery needs — and ``initial_sequence.fasta`` if present, so a run given real DNA still
    founds its blocks from it. The rearrangements come back too, since they share the event table now.
    ``chromosome_events.tsv`` is not read: it records how the karyotype got its shape, and the shape
    itself is already in ``blocks.tsv``. What comes back reconstructs and writes exactly as the
    original did; it is not for evolving further."""
    d = pathlib.Path(directory)

    def read(name):
        try:
            return (d / name).read_text()
        except FileNotFoundError:
            raise FileNotFoundError(
                f"{d / name} not found — re-run 'zombi2 genomes --resolution nucleotide' with "
                f"{name.split('.')[0]!r} in --write, or leave --write off to get everything"
            ) from None

    spans, names, strands = _genes_from_tsv(read("genes.tsv"))
    events, rearrangements = _events_from_tsv(read("genome_events.tsv"))
    initial_sequence: dict[int, str] = {}
    fpath = d / "initial_sequence.fasta"
    if fpath.exists():                               # a run given real DNA; keyed by source id
        for sq, seq in read_fasta(fpath).items():
            initial_sequence[int(sq[len("source"):] if sq.startswith("source") else sq)] = seq
    return NucleotideGenomesResult(
        tree, _blocks_from_tsv(read("blocks.tsv")), events, rearrangements,
        [], None, spans, names, strands, _initial_genome_from_tsv(read("initial_genome.tsv")),
        initial_sequence)


def _valid_length(length) -> int:
    if isinstance(length, bool) or not isinstance(length, int) or length < 1:
        raise ValueError(f"a chromosome length must be a positive integer, got {length!r}")
    return length


def _replicon_specs(chromosomes, root_length, topology) -> list[tuple[int, str]]:
    """Resolve the ``chromosomes`` argument to a list of ``(length, topology)`` replicon specs. An int
    ``N`` gives ``N`` equal replicons of ``root_length`` and ``topology``; a list gives heterogeneous
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


def _initial_blocks(source, length, cp, intervals, new_family, gene_spans, gene_names,
                 gene_strands) -> list[Block]:
    """Lay one initial replicon of ``length`` down as its blocks: the declared genes at ``intervals``
    (0-based half-open, sorted, non-overlapping) and **intergene** everywhere else — the alternating
    chain. Every block shares the replicon's initial copy lineage ``cp`` (they are one copy of one
    replicon); each **gene** additionally gets a fresh **family** id, which is what makes it indivisible
    and gives it a gene tree, and is recorded in ``gene_spans``, ``gene_strands`` (its **coding**
    strand) and ``gene_names`` (when it is named).
    With no intervals the replicon is a single intergenic block — today's uniform sequence."""
    blocks, at = [], 0
    for (start, end, strand, name) in intervals:
        if start > at:
            blocks.append(Block(source, at, start, 1, cp))       # intergene before this gene
        fam = new_family()
        # NB: every initial block is strand +1. `Block.strand` is orientation *relative to the ancestral
        # source*, and at the start the genome IS its own source — nothing is inverted yet. A gene's
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
    spot = g._pick_legal_cut(rng)
    if spot is None:
        return 0
    chrom, start = spot
    ell = chrom._pick_arc_extent(start, duplication_length, rng)
    if ell is None:
        return 0
    copied = chrom.duplicate(start, ell, new_copy)
    if copied is None:
        return 0
    events.append(Duplication(t, node_id, chrom.id, copied))
    return sum(end - beg for (_par, _child, _src, beg, end) in copied)


def _do_transfer(rng, tree, alive, gen, kd, t, transfer_length, transfer_to, self_transfer, depth,
                 events, new_copy) -> int:
    """Copy a geometric-length arc of the donor lineage ``alive[kd]`` into a **contemporaneous
    recipient** (chosen by ``transfer_to``: uniform, or a :class:`Distance` weighting): the arc's copy
    lineages beget fresh children that arrive as a block at a random spot on a random recipient
    chromosome (strands travel with them). A horizontal edge in each block's gene tree. **Additive**
    — the donor keeps its copy — so it returns the recipient's length gain (0 on a no-op)."""
    donor_g = gen[kd]
    spot = donor_g._pick_legal_cut(rng)
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
    chrom.originate(p, length, src, cp, fam)
    gene_spans[fam] = (src, 0, length)                   # a de-novo gene, tracked like a declared one
    gene_strands[fam] = 1
    events.append(Origination(t, node_id, chrom.id, cp, src, 0, length))
    return length


def _do_loss(g, node_id, t, loss_length, rng, events) -> int:
    """Delete a geometric-length arc from a length-weighted chromosome — an ancestry-changing event (a
    death). Never empties a chromosome (leaves at least one nucleotide; whole-chromosome loss is a
    deferred tier event). Records the deleted material — which copy lineage lost which arc — as a
    :class:`Loss`. Returns the length removed as a **negative** delta (0 on a no-op)."""
    spot = g._pick_legal_cut(rng)
    if spot is None:
        return 0
    chrom, start = spot
    if chrom.length < 2:
        return 0
    ell = chrom._pick_arc_extent(start, loss_length, rng)
    if ell is None:
        return 0
    lost = chrom.delete(start, ell)
    if lost is None:
        return 0
    events.append(Loss(t, node_id, chrom.id, lost))
    return -sum(end - beg for (_cp, _src, beg, end) in lost)


def _do_inversion(g, node_id, t, inversion_length, rng, rearrangements) -> None:
    """Invert a geometric-length arc of a length-weighted chromosome (length-neutral)."""
    spot = g._pick_legal_cut(rng)
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
    spot = g._pick_legal_cut(rng)
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
    if source.n_genes and sum(1 for b in arc if b.is_gene) == source.n_genes:
        return                                           # would leave the donor geneless: no-op
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
    spot = g._pick_legal_cut(rng)
    if spot is None:
        return
    chrom, start = spot
    if chrom.length < 2:
        return
    ell = chrom._pick_arc_extent(start, transposition_length, rng)
    if ell is None:
        return
    intact = chrom.blocks                                # keep for rollback if there is nowhere to land
    arc = chrom.excise(start, ell)
    if arc is None:
        return
    flipped = bool(rng.random() < inversion_probability)
    dest = chrom._pick_legal_cut(rng)                   # a legal spot on the remainder
    if dest is None:
        chrom.blocks = intact                            # nowhere legal to land: undo the excision
        return
    chrom.place(arc, dest, flipped)
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
        if chrom.n_genes and not (any(x.is_gene for x in chrom.blocks[:i])
                                  and any(x.is_gene for x in chrom.blocks[i:])):
            return 0                                     # a half without a gene: the fission fails
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
        if chrom.n_genes and not (any(x.is_gene for x in chrom.blocks[i:j])
                                  and any(x.is_gene for x in chrom.blocks[:i] + chrom.blocks[j:])):
            return 0                                     # a half without a gene: the fission fails
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


def _do_chromosome_origination(g, node_id, t, origination_length, rng, events, chromosome_events,
                               new_chrom_id, new_source, new_copy, new_family,
                               gene_spans, gene_strands) -> tuple[int, int]:
    """A de-novo replicon (a plasmid): a fresh circular chromosome — a **root** of the chromosome
    network (no parent) — **carrying one new gene**. A chromosome never exists without a gene, so the
    replicon is born with one rather than as an empty shell: its own source, copy lineage and family,
    exactly like a de-novo :func:`_do_origination`. Returns ``(chromosome delta, length delta)``."""
    cid, src, cp, fam = new_chrom_id(), new_source(), new_copy(), new_family()
    length = max(1, int(rng.geometric(1.0 / origination_length)))
    g.chromosomes.append(Chromosome(cid, "circular", [Block(src, 0, length, 1, cp, fam)]))
    gene_spans[fam] = (src, 0, length)
    gene_strands[fam] = 1
    chromosome_events.append(ChromosomeEvent(t, "origination", node_id, (), (cid,)))
    events.append(Origination(t, node_id, cid, cp, src, 0, length))
    return (1, length)


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


@without_cyclic_gc
def simulate_genomes_nucleotide(tree, *, inversion=0.0, inversion_length=50.0, translocation=0.0,
                                translocation_length=50.0, transposition=0.0, transposition_length=50.0,
                                inversion_probability=0.0, loss=0.0, loss_length=50.0, duplication=0.0,
                                duplication_length=50.0, transfer=0.0, transfer_length=50.0,
                                transfer_to="uniform", self_transfer=False, origination=0.0,
                                origination_length=50.0, fission=0.0, fusion=0.0,
                                chromosome_origination=0.0, chromosome_loss=0.0, chromosomes=1,
                                root_length=1000, topology="circular", genes=0, gene_length=100,
                                gff=None, fasta=None, trim_overlaps=False, seed=None,
                                progress=False) -> NucleotideGenomesResult:
    """Evolve a nucleotide genome along a species tree by inversion, translocation, transposition,
    **loss**, **duplication**, **transfer**, **origination**, and the number-changing chromosome tier.
    The run starts from a **karyotype** — ``chromosomes`` replicons, each its own source: an int
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
      ``chromosome_origination`` (**per lineage**) adds a de-novo circular replicon (a plasmid, a
      network **root**) **carrying one new gene**; ``chromosome_loss`` (**per chromosome**) kills a
      whole chromosome (its material dies as a loss; never the last one) — a network **leaf**. All
      record a chromosome-network edge.

    **A chromosome never exists without a gene.** A replicon is born with one, and any event that
    would strip a chromosome of its last gene — a loss, a translocation carrying it away, a fission
    splitting off a geneless half — simply does not happen. (Vacuous when no genes are declared.)

    The engine runs a **global-timeline** Gillespie: all lineages alive at once evolve along one clock
    (every extension event is **per lineage** — the rate says how often a lineage does it, the extent how
    much it touches, so a bigger genome does not get proportionally more events; the chromosome tier is
    per chromosome), so a transfer
    couples two contemporaries. With loss, the strong invariant weakens: every node carries a **subset**
    of the initial sequence (each ancestral position at most once, monotonically down every path);
    origination further adds fresh sources beyond the root. Deterministic given ``seed``. (Transfer is
    additive for now; homologous *replacement* transfer is a later refinement.)"""
    tree = tree.complete_tree if isinstance(tree, SpeciesResult) else tree
    for label, rate in (("inversion", inversion), ("translocation", translocation),
                        ("transposition", transposition), ("loss", loss), ("duplication", duplication),
                        ("transfer", transfer), ("origination", origination), ("fission", fission),
                        ("fusion", fusion), ("chromosome_origination", chromosome_origination),
                        ("chromosome_loss", chromosome_loss)):
        if isinstance(rate, bool) or not isinstance(rate, (int, float)):
            # a Rate expression (OnTime / DrivenBy / ByFamily). The nucleotide engine holds rates
            # constant, so it cannot honour one — say so, rather than crashing on `rate < 0` below.
            raise ValueError(
                f"{label} carries a rate modifier, but the nucleotide genome engine takes constant "
                "rates only — a skyline (OnTime) or a conditioned/driven rate (DrivenBy) is not wired "
                "here. Driving a nucleotide rate with a trait is a later slice (see the conditioning "
                "chapter). Pass a plain number.")
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
    initial_sequence: dict[int, str] = {}               # {source: initial DNA}, empty unless a FASTA is given
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
        if fasta is not None:                        # the initial DNA, one record per replicon, by seqid
            seqs = read_fasta(fasta)
            if set(seqs) != set(seqids):
                raise ValueError(
                    f"the FASTA's records {sorted(seqs)} do not match the GFF's replicons "
                    f"{seqids} — every ##sequence-region needs exactly one > record, same id")
            for i, sq in enumerate(seqids):
                if len(seqs[sq]) != lengths[sq]:
                    raise ValueError(
                        f"replicon {sq!r} is {lengths[sq]} bp in the GFF but {len(seqs[sq])} bp in "
                        "the FASTA — the sequence must be exactly as long as its sequence-region")
                initial_sequence[i] = seqs[sq]
    else:
        if fasta is not None:
            raise ValueError("fasta= needs gff=: the FASTA's records are matched to the GFF's "
                             "replicons by id, so there is nothing to lay down without one")
        specs = _replicon_specs(chromosomes, root_length, topology)
        for _length, _top in specs:                  # the genes must fit; they need not leave a gap
            if genes and genes * gene_length > _length:
                raise ValueError(f"{genes} genes of {gene_length} bp do not fit in a {_length} bp "
                                 f"replicon")
        layouts = [_even_gene_intervals(length, genes, gene_length) for (length, _t) in specs]
    rates = _Rates(inversion, translocation, transposition, loss, duplication, transfer, origination,
                   fission, fusion, chromosome_origination, chromosome_loss, inversion_length,
                   translocation_length, transposition_length, loss_length, duplication_length,
                   transfer_length, origination_length, inversion_probability)
    depth = mean_root_to_tip(tree)                       # timescale for Distance weighting

    rng = np.random.default_rng(seed)
    chrom_counter = 0
    copy_counter = 0
    source_counter = len(specs)                          # de-novo sources continue past the initial sources

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

    initial_chroms = []
    gene_spans: dict[int, tuple[int, int, int]] = {}
    gene_names: dict[str, int] = {}
    gene_strands: dict[int, int] = {}
    for source, ((length, top), intervals) in enumerate(zip(specs, layouts)):  # one source per replicon
        cid = new_chrom_id()
        cp = new_copy()                                 # ...and one initial copy lineage per replicon
        initial_chroms.append(Chromosome(cid, top, _initial_blocks(source, length, cp, intervals, new_family,
                                                             gene_spans, gene_names,
                                                             gene_strands)))
        chromosome_events.append(ChromosomeEvent(root.birth_time, "origination", root.id, (), (cid,)))
        events.append(Origination(root.birth_time, root.id, cid, cp, source, 0, length, initial=True))

    # the run's starting genome: a deep snapshot, so the live genome's events never reach it
    initial_genome = NucleotideGenome(
        [Chromosome(c.id, c.topology, [Block(b.source, b.start, b.end, b.strand, b.copy, b.gene)
                                       for b in c.blocks]) for c in initial_chroms])

    t = root.birth_time
    alive: list[int] = []                               # the live-lineage set (species._grow shape)
    gen: list[NucleotideGenome] = []
    pos: dict[int, int] = {}
    enter(alive, gen, pos, root.id, NucleotideGenome(initial_chroms))
    total_length = sum(c.length for c in initial_chroms)
    total_chromosomes = len(initial_chroms)

    bar = progress_bar(len(schedule), "genomes", unit="branch", enabled=progress)
    si = 0
    while si < len(schedule):
        bar.to(si)
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
                    dc, dl = _do_chromosome_origination(
                        gen[k], alive[k], t, rates.origination_length, rng, events, chromosome_events,
                        new_chrom_id, new_source, new_copy, new_family, gene_spans, gene_strands)
                    total_chromosomes += dc
                    total_length += dl
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
    bar.close()
    return NucleotideGenomesResult(tree, genomes, events, rearrangements, chromosome_events, seed,
                                  gene_spans, gene_names, gene_strands, initial_genome, initial_sequence)


# --- the gene-tree recovery: root partition -> per-block genealogy -> one tree per block ----------
#
# A gene tree per *root-block*: the coarsest interval of a source that is never cut in any node
# and survives there. Within such a block every copy is un-cut in every leaf that carries it (a cut
# would be an observable breakpoint, which would have split the block), so all its copies share one
# genealogy. The recovery has two moves:
#
#   1. the **root partition** — cut each source at the union of every node's breakpoints and keep the
#      intervals some node still covers (the material that exists anywhere);
#   2. per block, replay the copy-lineage log restricted to that block into the **per-segment** model
#      the shared ``gene_trees_from_events`` reads (every event ends a segment and starts fresh ids),
#      so a duplication is a ladder (parent continues + child) and a speciation a bifurcation.
#
# Because the partition is by *surviving* breakpoints, every event that reaches an extant copy is
# atomic on it (covers the whole block or none) — an event that only partially covers a block made a
# breakpoint that did not survive, so its child is dead and irrelevant to the extant tree.


@dataclass(frozen=True, slots=True)
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
    """The root partition: per source, the maximal intervals bounded by the breakpoints of **every**
    node's genome, that some node still carries. Returns a sorted ``(source, start, end)`` list.

    Every genome the run holds votes — every node, and the initial genome — not only the extant
    leaves, and that is what makes the whole tree reconstructable. Cutting at the survivors alone leaves two kinds of hole: material no survivor
    kept has no block at all, and an ancestor can hold a *fragment* of a block whose genealogy — being
    the survivors' — has no lineage for it. Counting every node closes both at once, because a node's
    own breakpoints are then all in the partition, so its every block is a whole number of root blocks.

    Reading the **final** genome at each node is enough: a breakpoint matters only where material
    survives on one side of it and not the other, and the surviving side carries that breakpoint to
    the end of its branch. A boundary that vanishes took its material with it."""
    bounds: dict[int, set[int]] = collections.defaultdict(set)
    spans: dict[int, set[tuple[int, int]]] = collections.defaultdict(set)
    for genome in [*result.genomes.values(), result.initial_genome]:
        for chrom in genome.chromosomes:
            for b in chrom.blocks:
                bounds[b.source].update((b.start, b.end))
                spans[b.source].add((b.start, b.end))
    blocks: list[tuple[int, int, int]] = []
    for source in bounds:
        cuts = sorted(bounds[source])
        covers = spans[source]
        for a, c in zip(cuts, cuts[1:]):
            if any(x <= a and c <= y for (x, y) in covers):   # some node still carries [a, c)
                blocks.append((source, a, c))
    return sorted(blocks)


def _emit_block_events(fam, s, a, b, tree, origs, dups, transfers, losses, specs, new_seg, out,
                       tip_of) -> None:
    """Replay the copy-lineage log for one root-block ``(s, [a, b))`` into per-segment events on
    ``out``. A copy lineage is a *block-copy* when it covers ``[a, b)`` in full; a duplication or a
    transfer that covers it in full begets a block-copy child (a ladder rung on its parent — the
    transfer's child on the recipient branch, a horizontal edge), a speciation re-mints it, and any
    loss overlapping it ends it.

    ``tip_of[(fam, copy)]`` collects the **last** gene id each copy lineage held: an event ends a gene
    and starts a fresh id, so a copy that duplicated twice is three genes in a row, and the one a
    node's genome still carries is the last rung of that ladder. This is the join between the blocks
    a genome is made of and the sequences evolved down their trees."""
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
        # The gene a genome still carrying c holds. ``None`` when a loss ended c over this block —
        # and then no node carries that block under c at all, since a copy's blocks are disjoint in
        # source coordinates and every node's breakpoints are in the partition. Recorded rather than
        # dropped so the assembly's guard can tell "died here" from "never existed".
        tip_of[(fam, c)] = None if c in loss_of else prev
        if c in loss_of:                                   # a death (dead leaf)
            out.append(_SegEvent("loss", fam, species[c], loss_of[c], prev, None))
        elif c in specs:                                   # a bifurcation into the daughter species
            pnode = tree.nodes[specs[c].lineage]
            for i, d in enumerate(specs[c].children):
                if d in block_copies:
                    out.append(_SegEvent("speciation", fam, pnode.children[i], specs[c].time,
                                         seg_in[d], prev))
        # else: prev survives to an extant/extinct leaf — gene_trees_from_events tags it by species fate


def _recover_gene_trees(result, *, every_block: bool = False
                        ) -> tuple[list[tuple[int, int, int]], dict[int, GeneTree],
                                   dict[tuple[int, int], int | None]]:
    """The full recovery: the root partition, a tree per family, and ``{(family, copy): gene id}`` —
    the last gene each copy lineage held, which is what a genome still carrying that copy is made of
    (``None`` where a loss ended that copy over that block, leaving nothing to read a fragment from).

    With genes declared we build a tree for the root-blocks that are declared genes, keyed by gene
    family id; with none declared, every root-block is a family (keyed by index). Either way the
    per-block builder below is the same — which block it is *pointed at* is the only difference.

    ``every_block`` points it at all of them: a block never splits, so its size is fixed and its whole
    genealogy is already in the log, exactly as a gene's is. The trees come back keyed by the block's
    index in the partition, so intergenic spacer is reconstructed on the same footing as a gene — which
    is what lets a sequence run rebuild a whole ancestral genome rather than a handful of loci."""
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

    if every_block:                                      # every root-block, genic or spacer, by index
        targets = list(enumerate(blocks))
    elif result.gene_spans:                              # genic: one family per surviving declared gene
        family_of = {span: fam for fam, span in result.gene_spans.items()}
        targets = [(family_of[iv], iv) for iv in blocks if iv in family_of]
    else:                                                # uniform: every root-block is its own family
        targets = list(enumerate(blocks))

    seg_events: list[_SegEvent] = []
    tip_of: dict[tuple[int, int], int] = {}
    for fam, (s, a, b) in targets:
        _emit_block_events(fam, s, a, b, tree, origs, dups, transfers, losses, specs, new_seg,
                           seg_events, tip_of)
    return blocks, gene_trees_from_events(seg_events, tree), tip_of


__all__ = ["Block", "Chromosome", "NucleotideGenome", "NucleotideGenomesResult",
           "Origination", "Loss", "Duplication", "Transfer", "Speciation", "Inversion",
           "Translocation", "Transposition", "Distance", "simulate_genomes_nucleotide"]
