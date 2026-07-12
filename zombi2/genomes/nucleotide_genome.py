"""Nucleotide-level circular genome — M1 (inversions only).

A genome is an ordered, **circular** sequence of :class:`Segment`s. Each segment is
*root-anchored*: it names the interval ``[src_start, src_end)`` on the **source**
(ancestral) genome it descends from, plus a ``strand`` (+1 forward, -1 reverse-
complement). Because the ancestor is never renumbered, tracing any present-day
position back to its origin is O(1) — the segment already carries its ancestral
coordinates (see :meth:`NucleotideGenome.to_cells`).

Events act on a variable-length **arc** of the circle, which may wrap the origin. The
one genuinely new primitive is :meth:`NucleotideGenome._split_at`: introduce a
breakpoint by cutting whatever segment straddles a coordinate into two collinear
pieces. Every interval event is then "split both ends, act on the whole segments now
covered". M1 implements INVERSION (reverse the arc's segment order + flip each strand);
duplication / loss / transfer land in later milestones through the same primitive.

This is the nucleotide-granular sibling of :class:`~zombi2.genome.OrderedGenome`, whose
genes are atomic and are never subdivided. Here sequence is continuous, so segments are
*emergent* — created by breakpoints, not declared up front.

A shared ``registry`` (``seg_id -> (source, src_start, src_end)``) records the provenance
of every segment ever minted, including the transient ones an event re-mints away; the
trace-back post-processing reads it to attribute events to ancestral coordinates. All
clones produced at speciation share the one registry (and the one :class:`IdManager`).
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from zombi2.genomes.events import EventType, GeneOp, Region, Selection, TransferSegment
from zombi2.genomes.genome import Chromosome, Gene, Genome, IdManager

#: Smallest a chromosome may shrink to under deletion — the min-genome floor. A deletion is
#: clamped so the genome never drops below this many nucleotides (guards against emptying the
#: chromosome, which would stall the size-proportional event process). One base, effectively
#: "never delete the last nucleotide"; insertions have no floor.
MIN_GENOME_LENGTH = 1


@dataclass(frozen=True)
class GeneInterval:
    """A gene as a half-open ancestral interval ``[start, end)`` on one ``source``.

    Genes are *indivisible*: no event breakpoint ever falls strictly inside one, so a gene
    is exactly one block wherever it survives — one genealogy per gene. ``gene_id`` is the
    immutable identity (kept even after pseudogenization).
    """

    gene_id: str
    source: str
    start: int
    end: int

    @property
    def length(self) -> int:
        return self.end - self.start


class SegmentRegistry:
    """Provenance shared across all genomes of one simulation (clones share one instance).

    ``provenance`` maps every segment id ever minted to its ancestral interval
    ``(source, src_start, src_end)`` — enough to attribute an event to blocks. ``split_parent``
    maps a segment born from a :meth:`NucleotideGenome._split_segment` to the segment it was
    cut from; those births are the only ones *not* in the event log (a split is a degree-2
    node), so recording them here keeps the segment genealogy connected for gene-tree
    reconstruction.

    ``genes`` (``source -> [GeneInterval]``, sorted by start) is the gene annotation, shared
    across clones so gene identity is a global property of *source coordinates* while
    functional state stays per-:class:`Segment`. ``_pending_genes`` holds the user's
    ``(start, end, name)`` intervals until the root chromosome is seeded and its source id is
    known (see :meth:`consume_pending_genes`).
    """

    __slots__ = ("provenance", "split_parent", "genes", "_pending_genes")

    def __init__(self, pending_genes=None) -> None:
        self.provenance: dict[str, tuple[str, int, int]] = {}
        self.split_parent: dict[str, str] = {}
        self.genes: dict[str, list[GeneInterval]] = {}
        self._pending_genes: list[tuple[int, int, str | None]] = list(pending_genes or [])

    def has_genes(self) -> bool:
        return bool(self.genes) or bool(self._pending_genes)

    def register_gene(self, gi: GeneInterval) -> None:
        lst = self.genes.setdefault(gi.source, [])
        lst.append(gi)
        lst.sort(key=lambda g: g.start)

    def consume_pending_genes(self, source: str) -> list[GeneInterval]:
        """Bind the user's pending intervals to the freshly-minted seed ``source`` and register.

        Called once per seed chromosome by :meth:`NucleotideGenome.originate` (a genome may seed
        several initial chromosomes, each an independent copy under its own source). The template
        is *not* consumed — each call re-registers the same intervals under ``source`` and returns
        them sorted by start (also now in ``self.genes[source]``).
        """
        out: list[GeneInterval] = []
        for k, (a, b, name) in enumerate(sorted(self._pending_genes)):
            gi = GeneInterval(name if name else f"gene{k + 1}", source, a, b)
            self.register_gene(gi)
            out.append(gi)
        return out

    def gene_id_at(self, source: str, start: int, end: int) -> str | None:
        """The gene whose interval contains ``[start, end)`` (else ``None``) — block classifier."""
        for gi in self.genes.get(source, ()):  # sorted by start
            if gi.start <= start and end <= gi.end:
                return gi.gene_id
            if gi.start >= end:
                break
        return None


@dataclass(slots=True)
class Segment:
    """A maximal run of sequence that is contiguous *and* has one shared history.

    ``[src_start, src_end)`` is the half-open interval on ``source`` this run descends
    from (``src_start < src_end`` always); ``strand`` says whether it is read forward
    (+1) or reverse-complemented (-1) relative to that source.

    ``gene_id`` is ``None`` for intergene, else the ancestral gene this block belongs to;
    ``is_gene`` is the *functional* state on this lineage (``True`` = functional gene,
    ``False`` = pseudogene — a gene demoted to intergene by pseudogenization, sequence
    retained). A block with ``gene_id is not None`` is never cut internally (functional or
    pseudogenized), which is what keeps a gene to exactly one block.
    """

    seg_id: str
    source: str
    src_start: int
    src_end: int
    strand: int = 1
    gene_id: str | None = None
    is_gene: bool = True
    length: int = field(init=False)  # cached (src coords never change after construction)

    def __post_init__(self):
        self.length = self.src_end - self.src_start

    # A Segment also plays the driver's "gene" role during a transfer handoff
    # (TransferSegment.genes), where the loop reads ``.gid`` / ``.family``.
    @property
    def gid(self) -> str:
        return self.seg_id

    @property
    def family(self) -> str:
        return self.source


class NucleotideGenome(Genome):
    """Circular nucleotide genome supporting variable-length **inversions** (M1).

    Construct via a factory so the driver can clone it, e.g.
    ``genome_factory=lambda ids: NucleotideGenome(ids, root_length=1000, extension=0.99)``.
    ``extension`` is the per-position continuation probability of the geometric length
    model (mean inversion length ``1/(1-extension)``), mirroring ``OrderedGenome``.
    """

    def __init__(self, ids: IdManager, root_length: int = 1000,
                 extension: float | None = 0.99, registry: SegmentRegistry | None = None,
                 pseudogenization: float = 0.0, replacement: float = 0.0,
                 indel_mean_length: float = 10.0, initial_chromosomes: int = 1):
        self.ids = ids
        self.root_length = int(root_length)
        self.extension = extension
        self._registry: SegmentRegistry = registry if registry is not None else SegmentRegistry()
        # The shared Chromosome container (bp coordinates): the genome is a dict of circular
        # chromosomes, each holding ``Segment`` elements. ``_segments`` is a view over the main
        # chromosome, so the single-chromosome engine is unchanged. The main chromosome takes the
        # canonical id 0 — matching the historical nucleotide ``Region(chromosome=0)`` convention;
        # additional chromosomes (initial copies, fission / fusion / plasmids) mint distinct ids
        # from :meth:`IdManager.new_chromosome`.
        self._cid = 0
        self.chromosomes: dict[int, Chromosome] = {self._cid: Chromosome(self._cid, True)}
        # number of root chromosomes to lay down at seeding (default 1); each an independent copy of
        # the root chromosome under its own source. ``_seeded`` flips true on the first origination
        # so seeding happens exactly once, then every later origination adds a novel gene.
        self._initial_chromosomes = max(1, int(initial_chromosomes))
        self._seeded = False
        # genic-model parameters (0 in the plain nucleotide model); inherited by clones
        self.pseudogenization = pseudogenization  # P(a loss on genes demotes instead of deletes)
        self.replacement = replacement            # P(a transfer is a homologous replacement)
        # intergenic indels: mean length of an insertion/deletion run (geometric, own parameter —
        # independent of `extension`, which drives the structural events). Inherited by clones.
        self.indel_mean_length = float(indel_mean_length)
        self._last_replaced: list[Segment] | None = None  # recipient homolog removed by a transfer

    # --- single-chromosome container (parallel-change: the nucleotide engine still works in one
    # flat circular sequence; ``_segments`` is a view over the one chromosome's elements) --------
    @property
    def _segments(self) -> list:
        return self.chromosomes[self._cid].elements

    @_segments.setter
    def _segments(self, value: list) -> None:
        self.chromosomes[self._cid].elements = value

    @property
    def _length(self) -> int:
        """Total nucleotide length across all chromosomes (computed from the segments, not
        maintained). Per-chromosome length is :meth:`_chrom_length`."""
        return self.size()

    # --- minting ------------------------------------------------------------
    def _new_segment(self, source: str, a: int, b: int, strand: int,
                     gene_id: str | None = None, is_gene: bool = True) -> Segment:
        seg = Segment(self.ids.new_gene(), source, a, b, strand, gene_id, is_gene)
        self._registry.provenance[seg.seg_id] = (source, a, b)
        return seg

    # --- queries ------------------------------------------------------------
    def _iter_segments(self):
        """Every segment across all chromosomes, in chromosome (insertion) then physical order."""
        for chrom in self.chromosomes.values():
            yield from chrom.elements

    def size(self) -> int:
        # nucleotide length -> inversion rate scales with genome size (summed over chromosomes)
        return sum(self._chrom_length(c) for c in self.chromosomes.values())

    def total_length(self) -> float:
        return float(self.size())

    def n_segments(self) -> int:
        return sum(len(c.elements) for c in self.chromosomes.values())

    def families(self) -> list[str]:
        return list(dict.fromkeys(s.source for s in self._iter_segments()))

    def copy_number(self, family: str) -> int:
        return sum(1 for s in self._iter_segments() if s.source == family)

    def genes(self) -> list[Gene]:
        return [Gene(s.seg_id, s.source) for s in self._iter_segments()]

    def presence_vector(self, family_order) -> np.ndarray:
        present = {s.source for s in self._iter_segments()}
        return np.fromiter((1 if f in present else 0 for f in family_order),
                           dtype=np.int8, count=len(family_order))

    def supported_events(self) -> frozenset[EventType]:
        return frozenset((EventType.ORIGINATION, EventType.INVERSION, EventType.LOSS,
                          EventType.DUPLICATION, EventType.TRANSFER, EventType.TRANSPOSITION,
                          EventType.INSERTION, EventType.DELETION,
                          EventType.CHROMOSOME_ORIGINATION, EventType.CHROMOSOME_LOSS,
                          EventType.FISSION, EventType.FUSION))

    # --- trace-back view: one cell per nucleotide, in physical order --------
    def to_cells(self) -> list[tuple[str, int, int]]:
        """Expand to ``[(source, src_pos, strand), ...]`` of length ``size()``.

        This *is* the position-level trace-back: cell ``p`` names the ancestral
        nucleotide the present-day position ``p`` descends from. Used by the tests as
        the ground truth and by the post-processor.
        """
        out: list[tuple[str, int, int]] = []
        for seg in self._iter_segments():
            if seg.strand == 1:
                out.extend((seg.source, p, 1) for p in range(seg.src_start, seg.src_end))
            else:
                out.extend((seg.source, p, -1) for p in range(seg.src_end - 1, seg.src_start - 1, -1))
        return out

    # --- seeding / origination ---------------------------------------------
    def originate(self, rng, params) -> list[GeneOp]:
        """Create brand-new sequence under a fresh ``source`` namespace.

        The first call (the seed) lays down the ``initial_chromosomes`` root chromosomes — each an
        independent full-length copy (tiled into gene / intergene segments when a gene annotation
        was supplied), under its own source. Later calls insert a *novel* gene of geometric length
        at a random position — it descends from no ancestral position, so it opens its own source
        and its own blocks (in genic mode it is a whole gene), with the gene tree rooted at this
        origination time.
        """
        if not self._seeded:                        # seed: the initial root chromosome(s)
            self._seeded = True
            return self._seed_root_chromosomes()

        source = self.ids.new_family()
        chrom = self.chromosomes[self._cid]
        # novel gene inserted somewhere
        length = self._draw_length(rng, params, chrom)
        L = self._chrom_length(chrom)
        if self._registry.has_genes():              # a novel gene is one whole, indivisible gene
            self._registry.register_gene(GeneInterval(source, source, 0, length))
            seg = self._new_segment(source, 0, length, 1, gene_id=source, is_gene=True)
            at = self._snap(int(rng.integers(L + 1)), chrom=chrom)
        else:
            seg = self._new_segment(source, 0, length, 1)
            at = int(rng.integers(L + 1))
        if at >= L:
            chrom.elements.append(seg)
        else:
            self._split_at(at, chrom)
            chrom.elements.insert(self._index_at(at, chrom), seg)
        return [GeneOp(seg.seg_id, source, "origin")]

    def _new_root_chromosome(self) -> Chromosome:
        """Mint an additional empty circular chromosome (a distinct id from the gene counter)."""
        cid = self.ids.new_chromosome()
        chrom = Chromosome(cid, True)
        self.chromosomes[cid] = chrom
        return chrom

    def _seed_root_chromosomes(self) -> list[GeneOp]:
        """Lay down the ``initial_chromosomes`` root chromosomes, returning every origin op.

        Chromosome 0 is the one the constructor already created; each further chromosome is an
        independent copy of the root — same length and (in genic mode) the same gene layout, but
        under its own source namespace, so the copies are distinct lineages that evolve apart.
        """
        ops: list[GeneOp] = []
        for k in range(self._initial_chromosomes):
            chrom = self.chromosomes[self._cid] if k == 0 else self._new_root_chromosome()
            source = self.ids.new_family()
            pending = self._registry.consume_pending_genes(source)
            if pending:
                ops.extend(self._seed_with_genes(source, pending, chrom))
            else:
                seg = self._new_segment(source, 0, self.root_length, 1)
                chrom.elements.append(seg)
                ops.append(GeneOp(seg.seg_id, source, "origin"))
        return ops

    def _seed_with_genes(self, source: str, pending: list[GeneInterval], chrom) -> list[GeneOp]:
        """Tile ``chrom`` into intergene / gene segments (genes get their own block).

        One ORIGINATION row per seed segment, so each becomes the root of the blocks it covers.
        """
        ops: list[GeneOp] = []
        cursor = 0
        for gi in pending:                          # sorted, non-overlapping
            if gi.start > cursor:
                ig = self._new_segment(source, cursor, gi.start, 1)
                chrom.elements.append(ig)
                ops.append(GeneOp(ig.seg_id, source, "origin"))
            g = self._new_segment(source, gi.start, gi.end, 1, gene_id=gi.gene_id, is_gene=True)
            chrom.elements.append(g)
            ops.append(GeneOp(g.seg_id, source, "origin"))
            cursor = gi.end
        if cursor < self.root_length:
            ig = self._new_segment(source, cursor, self.root_length, 1)
            chrom.elements.append(ig)
            ops.append(GeneOp(ig.seg_id, source, "origin"))
        return ops

    # --- the split primitive (the whole novelty vs OrderedGenome) ----------
    def _split_segment(self, seg: Segment, o: int) -> tuple[Segment, Segment]:
        """Split ``seg`` after ``o`` physical positions into (left, right), strand-aware."""
        if seg.strand == 1:
            left = self._new_segment(seg.source, seg.src_start, seg.src_start + o, 1,
                                     seg.gene_id, seg.is_gene)
            right = self._new_segment(seg.source, seg.src_start + o, seg.src_end, 1,
                                      seg.gene_id, seg.is_gene)
        else:  # reversed: first o physical positions are the HIGH end of the source
            left = self._new_segment(seg.source, seg.src_end - o, seg.src_end, -1,
                                     seg.gene_id, seg.is_gene)
            right = self._new_segment(seg.source, seg.src_start, seg.src_end - o, -1,
                                      seg.gene_id, seg.is_gene)
        # a split is a degree-2 birth (unlogged); record it so the genealogy stays connected
        self._registry.split_parent[left.seg_id] = seg.seg_id
        self._registry.split_parent[right.seg_id] = seg.seg_id
        return left, right

    def _chrom_length(self, chrom: Chromosome) -> int:
        """Total nucleotide length of ``chrom`` (sum of its segment lengths)."""
        return sum(seg.length for seg in chrom.elements)

    def _choose_chromosome_weighted(self, rng) -> Chromosome:
        """Pick a chromosome, weighted by its nucleotide length (where an event lands).

        Short-circuits without drawing when there is only one chromosome, so a single-chromosome
        genome consumes exactly the RNG draws of the pre-multichromosome engine (byte-identity)."""
        if len(self.chromosomes) == 1:
            return next(iter(self.chromosomes.values()))
        chroms = list(self.chromosomes.values())
        sizes = [self._chrom_length(c) for c in chroms]
        r = int(rng.integers(sum(sizes)))
        for c, s in zip(chroms, sizes):
            if r < s:
                return c
            r -= s
        return chroms[-1]  # pragma: no cover -- unreachable while total > 0

    def _choose_chromosome_uniform(self, rng) -> int:
        """Pick a chromosome uniformly (return its ``chrom_id``), for transfer insertion — every
        chromosome is an equally likely recipient regardless of length (the Q4 behaviour, matching
        ``OrderedGenome``). Draws nothing when there is a single chromosome (byte-identity)."""
        if len(self.chromosomes) == 1:
            return next(iter(self.chromosomes))
        cids = list(self.chromosomes)
        return cids[int(rng.integers(len(cids)))]

    def _split_at(self, c: int, chrom: Chromosome | None = None) -> None:
        """Ensure a segment boundary at physical coordinate ``c`` on ``chrom`` (``0`` is always one).

        ``chrom`` defaults to the genome's single chromosome (the migration is still one chromosome
        end to end). The genic invariant — ``c`` never lands strictly inside a gene block — is enforced
        upstream by :meth:`_snap`; the assertion here is the cheap safety net (a no-op for intergene).
        """
        if c == 0:
            return
        els = (chrom if chrom is not None else self.chromosomes[self._cid]).elements
        pos = 0
        for i, seg in enumerate(els):
            if pos == c:
                return  # already a boundary
            seglen = seg.length
            if pos < c < pos + seglen:
                assert seg.gene_id is None, f"illegal split inside gene {seg.gene_id!r} at {c}"
                els[i:i + 1] = list(self._split_segment(seg, c - pos))
                return
            pos += seglen
        # c == chromosome length: a boundary (wrap point); no-op

    def _snap(self, c: int, direction: int = 0, chrom: Chromosome | None = None) -> int:
        """Move ``c`` out of any gene interior to a legal breakpoint on ``chrom`` (no-op in intergene).

        ``direction`` picks the boundary when ``c`` is inside a gene: ``-1`` = the gene's
        physical start, ``+1`` = its end, ``0`` = the nearer of the two. Positions in intergene
        (or on an existing boundary) are returned unchanged — cutting there is legal and is how
        intergene blocks form. ``chrom`` defaults to the genome's single chromosome. Walks the
        segment list (O(n), like :meth:`_split_at`), so it stays correct mid-mutation.
        """
        pos = 0
        for seg in (chrom if chrom is not None else self.chromosomes[self._cid]).elements:
            nxt = pos + seg.length
            if c == pos:
                return c                       # already a boundary
            if c < nxt:                        # pos < c < nxt: strictly inside this segment
                if seg.gene_id is None:
                    return c                   # intergene interior: a legal cut
                if direction < 0:
                    return pos
                if direction > 0:
                    return nxt
                return pos if (c - pos) <= (nxt - c) else nxt
            pos = nxt
        return c                               # c == total length: the wrap / append boundary

    def _index_at(self, phys: int, chrom: Chromosome | None = None) -> int:
        els = (chrom if chrom is not None else self.chromosomes[self._cid]).elements
        pos = 0
        for i, seg in enumerate(els):
            if pos == phys:
                return i
            pos += seg.length
        raise ValueError(f"no segment boundary at physical {phys}")

    # --- inversion ----------------------------------------------------------
    def _apply_inversion(self, s: int, ell: int, chrom: Chromosome | None = None) -> list[Segment]:
        """Invert the circular arc ``[s, s+ell)`` on ``chrom`` and return the reversed segments.

        A non-wrapping arc is reversed in place with the physical origin fixed. A
        **wrapping** arc has no origin-fixed representation without inventing a spurious
        breakpoint at position 0, so the ring is first rotated to bring the arc to the
        front — the origin drifts to a real breakpoint. This is harmless: a circular
        genome has no privileged origin, blocks live in origin-independent *source*
        coordinates, and the invariants are rotation-invariant.
        """
        chrom = chrom if chrom is not None else self.chromosomes[self._cid]
        els = chrom.elements
        L = self._chrom_length(chrom)
        if L == 0:
            return []
        ell = max(1, min(ell, L))
        s %= L
        e = (s + ell) % L
        self._split_at(s, chrom)
        self._split_at(e, chrom)  # no-op when e == 0 (arc ends at the wrap) or e == s (whole genome)

        if s + ell <= L:  # --- non-wrapping: origin fixed ---
            i_s = self._index_at(s, chrom)
            i_e = self._index_at(s + ell, chrom) if s + ell < L else len(els)
        else:             # --- wrapping: rotate the arc to the front, origin drifts ---
            i_s = self._index_at(s, chrom)
            els[:] = els[i_s:] + els[:i_s]
            i_s, i_e, acc = 0, 0, 0
            while acc < ell:
                acc += els[i_e].length
                i_e += 1

        arc = els[i_s:i_e]
        arc.reverse()                     # segment order flips...
        for seg in arc:
            seg.strand = -seg.strand       # ...and each segment's strand flips
        els[i_s:i_e] = arc
        return arc

    # --- deletion -----------------------------------------------------------
    def _apply_loss(self, s: int, ell: int, chrom: Chromosome | None = None) -> list[Segment]:
        """Delete the circular arc ``[s, s+ell)`` on ``chrom`` and return the removed segments.

        Boundaries at ``s`` and the arc end are split first, so the arc is a whole
        number of segments; the survivors keep their circular order (the physical origin
        follows the surviving material when the deletion swallows position 0).
        """
        chrom = chrom if chrom is not None else self.chromosomes[self._cid]
        els = chrom.elements
        L = self._chrom_length(chrom)
        if L == 0:
            return []
        ell = min(ell, L)
        s %= L
        if ell == L:                       # deletes the whole chromosome's content
            removed = els[:]
            els[:] = []
            return removed
        e = (s + ell) % L
        self._split_at(s, chrom)
        self._split_at(e, chrom)

        wrapping = s + ell > L
        removed, kept, pos = [], [], 0
        for seg in els:
            start = pos
            pos += seg.length
            in_arc = (start >= s or pos <= e) if wrapping else (start >= s and pos <= s + ell)
            (removed if in_arc else kept).append(seg)
        els[:] = kept
        return removed

    def _apply_loss_or_pseudogenize(self, s: int, ell: int, rng, chrom: Chromosome | None = None):
        """Resolve a LOSS: either delete the arc, or *pseudogenize* the genes it contains.

        Returns ``("loss", removed_segments)`` or ``("pseudo", [(old, cont), ...])``. With
        probability ``pseudogenization`` — and only when the arc holds a functional gene — the
        gene segments are re-minted in place with ``is_gene=False`` (same ``gene_id`` / ancestral
        coords, all sequence retained, length unchanged); every other arc segment stays. This is
        a logged ``parent -> continuation`` edge on the gene lineage (a state change), not a
        split. Otherwise the arc is deleted as usual.
        """
        chrom = chrom if chrom is not None else self.chromosomes[self._cid]
        els = chrom.elements
        if (self.pseudogenization > 0.0 and self._registry.has_genes()
                and self._chrom_length(chrom) > 0):
            i_s, i_e = self._arc_range(s, ell, chrom)   # splits the (snapped, legal) ends
            arc = els[i_s:i_e]
            if (any(sg.gene_id is not None and sg.is_gene for sg in arc)
                    and rng.random() < self.pseudogenization):
                demoted, new_arc = [], []
                for sg in arc:
                    if sg.gene_id is not None and sg.is_gene:
                        cont = self._new_segment(sg.source, sg.src_start, sg.src_end, sg.strand,
                                                 gene_id=sg.gene_id, is_gene=False)
                        demoted.append((sg, cont))
                        new_arc.append(cont)
                    else:
                        new_arc.append(sg)
                els[i_s:i_e] = new_arc      # length unchanged (sequence retained)
                return "pseudo", demoted
            # coin said delete (or no functional gene): remove the already-split arc
            removed = els[i_s:i_e]
            del els[i_s:i_e]
            return "loss", removed
        return "loss", self._apply_loss(s, ell, chrom)

    # --- duplication --------------------------------------------------------
    def _arc_range(self, s: int, ell: int, chrom: Chromosome | None = None) -> tuple[int, int]:
        """Split the ends of arc ``[s, s+ell)`` on ``chrom``, return its ``[i_s, i_e)`` segment range.

        A wrapping arc is rotated to the front first (the origin drifts to a real breakpoint), so
        the arc is always a contiguous slice ``chrom.elements[i_s:i_e]``.
        """
        chrom = chrom if chrom is not None else self.chromosomes[self._cid]
        els = chrom.elements
        L = self._chrom_length(chrom)
        ell = max(1, min(ell, L))
        s %= L
        e = (s + ell) % L
        self._split_at(s, chrom)
        self._split_at(e, chrom)
        if s + ell <= L:
            i_s = self._index_at(s, chrom)
            i_e = self._index_at(s + ell, chrom) if s + ell < L else len(els)
            return i_s, i_e
        i_s = self._index_at(s, chrom)
        els[:] = els[i_s:] + els[:i_s]
        i_e, acc = 0, 0
        while acc < ell:
            acc += els[i_e].length
            i_e += 1
        return 0, i_e

    def _apply_duplication(self, s: int, ell: int,
                           chrom: Chromosome | None = None) -> list[list[GeneOp]]:
        """Duplicate the arc ``[s, s+ell)`` in tandem; return one op-group per segment.

        Each arc segment forks into a re-minted *continuation* (in place) and a *copy*
        (the paralog, inserted immediately after the arc). Both carry the segment's source
        interval, so the two are paralogs that coalesce at this duplication.
        """
        chrom = chrom if chrom is not None else self.chromosomes[self._cid]
        els = chrom.elements
        L = self._chrom_length(chrom)
        if L == 0:
            return []
        ell = max(1, min(ell, L))
        s %= L
        i_s, i_e = self._arc_range(s, ell, chrom)
        arc = els[i_s:i_e]
        groups, conts, copies = [], [], []
        for a in arc:
            cont = self._new_segment(a.source, a.src_start, a.src_end, a.strand, a.gene_id, a.is_gene)
            copy = self._new_segment(a.source, a.src_start, a.src_end, a.strand, a.gene_id, a.is_gene)
            conts.append(cont)
            copies.append(copy)
            groups.append([GeneOp(a.seg_id, a.source, "parent"),
                           GeneOp(cont.seg_id, a.source, "left"),
                           GeneOp(copy.seg_id, a.source, "right")])
        els[i_s:i_e] = conts          # continuations replace the originals...
        els[i_e:i_e] = copies          # ...and the copy block lands in tandem
        return groups

    # --- transposition (cut a segment, paste it elsewhere) -----------------
    def _apply_transposition(self, s: int, ell: int, dest: int,
                             chrom: Chromosome | None = None) -> list[Segment]:
        """Cut the arc ``[s, s+ell)`` on ``chrom`` and splice it back in at physical ``dest``.

        Content- and length-preserving and genealogically neutral (segments keep their
        ids), so it only permutes the mosaic — like inversion, no lineage re-mint.
        """
        chrom = chrom if chrom is not None else self.chromosomes[self._cid]
        els = chrom.elements
        L = self._chrom_length(chrom)
        if L <= 1:
            return []
        ell = max(1, min(ell, L - 1))
        i_s, i_e = self._arc_range(s, ell, chrom)
        arc = els[i_s:i_e]
        del els[i_s:i_e]
        rem = L - ell
        dest %= rem + 1
        if self._registry.has_genes() and dest < rem:  # never paste into a gene interior
            dest = self._snap(dest, chrom=chrom)        # (may snap up to the wrap boundary == rem)
        if dest >= rem:
            els.extend(arc)
        else:
            self._split_at(dest, chrom)
            idx = self._index_at(dest, chrom)
            els[idx:idx] = arc
        return arc

    def _draw_length(self, rng, params, chrom: Chromosome | None = None) -> int:
        ext = params.extension if params.extension is not None else self.extension
        n = self._chrom_length(chrom if chrom is not None else self.chromosomes[self._cid])
        if ext is None or n <= 1:
            return 1
        if ext >= 1.0:
            return n
        # truncated geometric (mean 1/(1-ext)); drawn in one shot instead of a Python loop
        return min(int(rng.geometric(1.0 - ext)), n)

    # --- intergenic indels --------------------------------------------------
    def _draw_indel_length(self, rng) -> int:
        """Geometric indel length with mean ``indel_mean_length`` (own parameter, not ``extension``).

        ``rng.geometric(p)`` has mean ``1/p``, so ``p = 1/mean`` gives the requested mean and a
        per-step continuation probability ``1 - 1/mean``. Always at least 1 nucleotide.
        """
        m = self.indel_mean_length
        if not (m > 1.0):                       # mean <= 1 -> always a single nucleotide
            return 1
        return max(1, int(rng.geometric(1.0 / m)))

    def _intergene_run_at(self, phys: int,
                          chrom: Chromosome | None = None) -> tuple[int, int] | None:
        """Physical bounds ``[lo, hi)`` of the maximal intergene run covering ``phys`` on ``chrom``.

        A "run" is a contiguous stretch whose segments are all intergene (``gene_id is None``),
        bounded by genes (or the chromosome ends — runs do **not** wrap the origin, so an indel
        stays within one physical stretch). Returns ``None`` when ``phys`` lands inside a gene.
        Used to clamp an intergenic deletion so it never reaches into a neighbouring gene.
        """
        pos, run_lo, cur = 0, 0, None
        for seg in (chrom if chrom is not None else self.chromosomes[self._cid]).elements:
            nxt = pos + seg.length
            if seg.gene_id is None:
                if cur is None:
                    run_lo = pos            # a new intergene run opens here
                cur = nxt                   # ...and extends to here
            else:
                if cur is not None and run_lo <= phys < cur:
                    return (run_lo, cur)    # phys was in the run just closed by this gene
                cur = None
            if pos <= phys < nxt and seg.gene_id is not None:
                return None                 # phys is strictly inside a gene
            pos = nxt
        if cur is not None and run_lo <= phys < cur:
            return (run_lo, cur)            # trailing intergene run to the chromosome end
        return None

    def _apply_insertion(self, rng, chrom: Chromosome | None = None) -> list[GeneOp]:
        """Insert a run of novel nucleotides (a fresh source) at an intergene position on ``chrom``.

        The insertion lengthens the intergene stretch it lands in — it is its own block, with an
        independent genealogy rooted at this event (like an intergenic origination of new
        sequence). In genic mode the insertion point is snapped out of any gene interior, so a
        gene is never split; with no genes it may land anywhere. Returns the log op-group.
        """
        chrom = chrom if chrom is not None else self.chromosomes[self._cid]
        els = chrom.elements
        source = self.ids.new_family()
        length = self._draw_indel_length(rng)
        seg = self._new_segment(source, 0, length, 1)   # intergene (gene_id None) novel sequence
        L = self._chrom_length(chrom)
        at = int(rng.integers(L + 1))
        if self._registry.has_genes():
            at = self._snap(at, chrom=chrom)            # never inside a gene interior
        if at >= L:
            els.append(seg)
        else:
            self._split_at(at, chrom)
            els.insert(self._index_at(at, chrom), seg)
        return [GeneOp(seg.seg_id, source, "origin")]

    def _apply_deletion(self, rng, chrom: Chromosome | None = None) -> list[Segment]:
        """Delete a run of nucleotides from *within a single intergene* stretch of ``chrom``.

        Picks a physical start, snaps it out of any gene, then clamps the arc to the intergene run
        containing it (so the deletion can never reach into, span, or remove a gene) and to the
        min-genome floor (:data:`MIN_GENOME_LENGTH`). Returns the removed segments (possibly empty
        — e.g. no room to delete without touching a gene / hitting the floor). Never wraps.
        """
        chrom = chrom if chrom is not None else self.chromosomes[self._cid]
        L = self._chrom_length(chrom)
        if L <= MIN_GENOME_LENGTH:
            return []
        s = int(rng.integers(L))
        ell = self._draw_indel_length(rng)
        if self._registry.has_genes():
            s = self._snap(s, chrom=chrom)              # move a gene-interior start to a boundary
            run = self._intergene_run_at(s, chrom)
            if run is None:                             # snapped onto a gene start/boundary with
                run = self._intergene_run_at(s - 1, chrom) if s > 0 else None  # gene right: look left
            if run is None:
                return []                               # nowhere legal to delete near s
            lo, hi = run
            s = min(max(s, lo), hi)                     # keep the start inside the run
            ell = min(ell, hi - s)                      # ...and the arc entirely inside it
        # honour the min-genome floor: never shrink below MIN_GENOME_LENGTH
        ell = min(ell, L - MIN_GENOME_LENGTH)
        if ell <= 0:
            return []
        return self._apply_loss(s, ell, chrom)

    _TARGETABLE = frozenset((EventType.INVERSION, EventType.LOSS, EventType.DUPLICATION,
                             EventType.TRANSFER, EventType.TRANSPOSITION,
                             EventType.INSERTION, EventType.DELETION))

    def draw_target(self, event, rng, params, family=None) -> Selection:
        if event not in self._TARGETABLE:
            raise ValueError(f"NucleotideGenome does not target {event!r}")
        chrom = self._choose_chromosome_weighted(rng)   # bp-length-weighted; no draw if one chromosome
        if event in (EventType.INSERTION, EventType.DELETION):
            # indels self-draw their position + length inside apply() (intergene-clamped), so the
            # selection is a placeholder — no arc is chosen up front and no RNG is consumed here.
            return Selection(genes=(), region=Region(chromosome=chrom.chrom_id, start=0, length=0))
        L = self._chrom_length(chrom)
        s = int(rng.integers(L))
        ell = self._draw_length(rng, params, chrom)
        if self._registry.has_genes():
            if ell >= L:                       # whole genome: keep it whole, just legalise the
                s = self._snap(s, -1, chrom)   # single split point (inversion splits at s)
                ell = L
            else:                              # snap start down, end up: genes stay whole; a
                s2 = self._snap(s, -1, chrom)  # sub-gene arc is promoted to the whole gene
                e2 = self._snap((s + ell) % L, +1, chrom)
                ell = (e2 - s2) % L or L
                s = s2
        return Selection(genes=(), region=Region(chromosome=chrom.chrom_id, start=s, length=ell))

    def apply(self, event, selection, rng, params) -> list[list[GeneOp]]:
        region = selection.region
        chrom = self.chromosomes[region.chromosome]
        if event is EventType.INVERSION:
            arc = self._apply_inversion(region.start, region.length, chrom)
            return [[GeneOp(seg.seg_id, seg.source, "inverted", orientation=seg.strand)
                     for seg in arc]]
        if event is EventType.LOSS:
            kind, payload = self._apply_loss_or_pseudogenize(region.start, region.length, rng, chrom)
            if kind == "loss":
                return [[GeneOp(seg.seg_id, seg.source, "lost")] for seg in payload]
            # pseudogenization: parent -> continuation (is_gene=False), one group per gene segment.
            # The driver logs these under the LOSS it fired; the role lets the post-processor
            # rewrite them to PSEUDOGENIZATION (a state change, not a terminal loss).
            return [[GeneOp(old.seg_id, old.source, "parent"),
                     GeneOp(cont.seg_id, cont.source, "pseudogenized")]
                    for (old, cont) in payload]
        if event is EventType.DUPLICATION:
            return self._apply_duplication(region.start, region.length, chrom)
        if event is EventType.TRANSPOSITION:
            L = self._chrom_length(chrom)
            if L <= 1:
                return []
            dest = int(rng.integers(L))
            arc = self._apply_transposition(region.start, region.length, dest, chrom)
            return [[GeneOp(seg.seg_id, seg.source, "transposed") for seg in arc]]
        if event is EventType.INSERTION:
            return [self._apply_insertion(rng, chrom)]  # one op-group: the novel intergene block
        if event is EventType.DELETION:
            removed = self._apply_deletion(rng, chrom)
            return [[GeneOp(seg.seg_id, seg.source, "lost")] for seg in removed]
        raise ValueError(f"NucleotideGenome does not handle {event!r}")

    # --- transfer handoff ---------------------------------------------------
    def extract_segment(self, selection, rng) -> TransferSegment:
        """Donor side: fork the arc — re-mint a continuation in place, build a copy to send.

        Each arc segment ``g`` forks into a continuation (stays in the donor) and a copy
        (travels to the recipient). The driver logs ``g -> [continuation, copy]`` as a
        TRANSFER, a bifurcation the block's gene tree reads as a transfer node. The copy
        keeps the arc's ancestral coordinates, so it is a xenolog of the same block.
        """
        region = selection.region
        chrom = self.chromosomes[region.chromosome]
        els = chrom.elements
        i_s, i_e = self._arc_range(region.start, region.length, chrom)
        arc = els[i_s:i_e]
        old_gids, cont_gids, copies, conts, arc_sources = [], [], [], [], []
        for a in arc:
            old_gids.append(a.seg_id)
            cont = self._new_segment(a.source, a.src_start, a.src_end, a.strand, a.gene_id, a.is_gene)
            conts.append(cont)
            cont_gids.append(cont.seg_id)
            copies.append(self._new_segment(a.source, a.src_start, a.src_end, a.strand,
                                            a.gene_id, a.is_gene))
            arc_sources.append((a.source, a.src_start, a.src_end))
        els[i_s:i_e] = conts          # donor length unchanged (arc -> continuations)
        # genic model: decide a homologous replacement and record the flanking genes (the
        # homology anchors) so the recipient can locate its syntenic copy to replace.
        replacement = self._registry.has_genes() and rng.random() < self.replacement
        left_flank = self._nearest_flank_gene(i_s - 1, -1, chrom) if replacement else None
        right_flank = self._nearest_flank_gene(i_e, +1, chrom) if replacement else None
        return TransferSegment(family=None, genes=tuple(copies),
                               donor_old_gids=old_gids, donor_cont_gids=cont_gids,
                               replacement=replacement, left_flank=left_flank,
                               right_flank=right_flank, arc_sources=tuple(arc_sources))

    def _nearest_flank_gene(self, start_idx: int, step: int,
                            chrom: Chromosome | None = None) -> tuple | None:
        """``(source, gene_id, strand)`` of the nearest gene outward from ``start_idx`` on ``chrom``
        (no wrap), else None. The strand is part of the homology anchor: a homologous locus in the
        recipient must carry the flank gene in the **same orientation**, so a flank that has since
        been inverted no longer anchors the replacement (it falls back to additive insertion)."""
        els = (chrom if chrom is not None else self.chromosomes[self._cid]).elements
        i = start_idx
        while 0 <= i < len(els):
            seg = els[i]
            if seg.gene_id is not None:
                return (seg.source, seg.gene_id, seg.strand)
            i += step
        return None

    def _find_homologous_span(self, segment,
                              chrom: Chromosome | None = None) -> tuple[int, int] | None:
        """Recipient span ``[i_s, i_e)`` on ``chrom`` between copies of the donor's flank genes, else
        None.

        "Search on the sides": locate the left-flank gene, then the next right-flank gene after
        it; the segments strictly between them are the homologous locus to replace (an empty span
        means the flanks are adjacent — the copy is inserted homologously with no removal). Both
        anchors must match the donor's ``(source, gene_id, strand)`` — same gene **and** same
        orientation — so an inverted flank breaks synteny and no homolog is found.
        """
        lf, rf = segment.left_flank, segment.right_flank
        if lf is None or rf is None:
            return None
        els = (chrom if chrom is not None else self.chromosomes[self._cid]).elements
        li = ri = None
        for idx, seg in enumerate(els):
            if seg.gene_id is None:
                continue
            if li is None and (seg.source, seg.gene_id, seg.strand) == lf:
                li = idx
            elif li is not None and (seg.source, seg.gene_id, seg.strand) == rf:
                ri = idx
                break
        if li is None or ri is None or ri <= li:
            return None
        return (li + 1, ri)

    def choose_insertion_point(self, segment, rng):
        """Recipient side: which chromosome, and where on it, the transferred copy lands.

        Returns ``(chrom_id, inner)`` where ``inner`` is ``("homolog", (i_s, i_e))`` when a
        replacement transfer finds its syntenic locus, else an integer physical position (snapped
        off any gene interior). The recipient chromosome is chosen **uniformly** (every chromosome
        is an equally likely target — the Q4 behaviour); a replacement homolog, when present, is
        honoured on whichever chromosome carries it. No homolog → additive.
        """
        if getattr(segment, "replacement", False) and self._registry.has_genes():
            # a self-transfer (recipient is the donor) still holds the arc's continuation
            # segments — homologous replacement would delete them, so fall back to additive
            cont_ids = set(segment.donor_cont_gids or ())
            is_self = any(seg.seg_id in cont_ids for seg in self._iter_segments())
            if not is_self:
                for chrom in self.chromosomes.values():       # find the homolog wherever it sits
                    span = self._find_homologous_span(segment, chrom)
                    if span is not None:
                        return (chrom.chrom_id, ("homolog", span))
        cid = self._choose_chromosome_uniform(rng)            # uniform recipient chromosome
        chrom = self.chromosomes[cid]
        at = int(rng.integers(self._chrom_length(chrom) + 1))  # any physical position, 0..length
        return (cid, self._snap(at, chrom=chrom) if self._registry.has_genes() else at)

    def insert_segment(self, segment, at, rng) -> list[GeneOp]:
        """Recipient side: splice the transferred copy block in — homologously, or additively.

        Sets ``_last_replaced`` to the recipient segments a homologous swap removed (``[]`` for a
        plain additive insertion in genic mode; ``None`` in the plain nucleotide model, which lets
        the driver fall back to its own random-removal replacement). :meth:`pop_replaced_segments`
        hands the removed list to the driver, which logs them as recipient losses.
        """
        cid, inner = at                        # (recipient chrom_id, position | homolog span)
        chrom = self.chromosomes[cid]
        els = chrom.elements
        copies = list(segment.genes)
        if isinstance(inner, tuple) and inner and inner[0] == "homolog":
            i_s, i_e = inner[1]
            removed = els[i_s:i_e]
            els[i_s:i_e] = copies
            self._last_replaced = removed
            return [GeneOp(seg.seg_id, seg.source, "transfer_copy") for seg in copies]
        if inner is None or inner >= self._chrom_length(chrom):
            idx = len(els)
        else:
            self._split_at(inner, chrom)
            idx = self._index_at(inner, chrom)
        els[idx:idx] = copies
        self._last_replaced = [] if self._registry.has_genes() else None
        return [GeneOp(seg.seg_id, seg.source, "transfer_copy") for seg in copies]

    def pop_replaced_segments(self):
        """Recipient segments a homologous transfer removed (for the driver to log as losses).

        ``None`` when this is the plain nucleotide model (no genes) — the signal for the driver to
        run its own random-removal replacement instead. A list (possibly empty) means the genome
        already handled replacement homologously. Resets after each read.
        """
        r = self._last_replaced
        self._last_replaced = None
        return r

    # --- speciation ---------------------------------------------------------
    def clone_reminting(self) -> tuple["NucleotideGenome", list[tuple[str, str, str]]]:
        child = NucleotideGenome(self.ids, self.root_length, self.extension, self._registry,
                                 self.pseudogenization, self.replacement,
                                 self.indel_mean_length, self._initial_chromosomes)
        child._seeded = True                 # a clone is seeded by copying, never by originate()
        mapping: list[tuple[str, str, str]] = []
        for i, pchrom in enumerate(self.chromosomes.values()):
            if i == 0:                       # the child already has its seed (main) chromosome, id 0
                cchrom = child.chromosomes[child._cid]
            else:                            # additional chromosomes get fresh ids
                cid = self.ids.new_chromosome()
                cchrom = Chromosome(cid, pchrom.circular)
                child.chromosomes[cid] = cchrom
            for seg in pchrom.elements:
                ns = child._new_segment(seg.source, seg.src_start, seg.src_end, seg.strand,
                                        seg.gene_id, seg.is_gene)
                cchrom.elements.append(ns)
                mapping.append((seg.seg_id, ns.seg_id, seg.source))
        return child, mapping

    # --- chromosome-tier events (origination / loss / fission / fusion) ------
    #
    # One tier above the segment events: they act on whole chromosomes. Fission / fusion only move
    # segments between chromosomes (segment ids untouched), so gene-tree reconstruction is unaffected;
    # only chromosome LOSS ends the lineages it carries. All nucleotide chromosomes are circular.
    def originate_chromosome(self, rng, params) -> int:
        """Originate a new, empty circular chromosome — a de-novo replicon (a plasmid). Returns the
        new chrom_id. Draws no rng."""
        cid = self.ids.new_chromosome()
        self.chromosomes[cid] = Chromosome(cid, True)
        return cid

    def lose_chromosome(self, rng) -> tuple[int, list[list[GeneOp]]]:
        """Lose a whole chromosome (chosen uniformly) and every segment on it, returning its chrom_id
        and one LOSS group per segment (so the caller can log the deaths). The caller guarantees
        >= 2 chromosomes (a genome keeps at least one)."""
        cids = list(self.chromosomes)
        cid = cids[int(rng.integers(len(cids)))]
        chrom = self.chromosomes.pop(cid)
        groups = [[GeneOp(seg.seg_id, seg.source, "lost")] for seg in chrom.elements]
        return cid, groups

    def fission(self, rng, params) -> tuple[int, int]:
        """Split a circular chromosome into two by excising the arc between two breakpoints into a
        new circular replicon (segments keep their ids; genes stay whole via snapping). The source
        is bp-length-weighted. Returns (source_cid, new_cid)."""
        chrom = self._choose_chromosome_weighted(rng)
        new_cid = self.ids.new_chromosome()
        L = self._chrom_length(chrom)
        if L <= 1:                            # nothing meaningful to split
            self.chromosomes[new_cid] = Chromosome(new_cid, True)
            return chrom.chrom_id, new_cid
        s = int(rng.integers(L))
        ell = int(rng.integers(1, L))         # arc length in [1, L-1]
        if self._registry.has_genes():        # snap so genes stay whole (start down, end up)
            s2 = self._snap(s, -1, chrom)
            e2 = self._snap((s + ell) % L, +1, chrom)
            ell = (e2 - s2) % L or L
            s = s2
        i_s, i_e = self._arc_range(s, ell, chrom)   # splits the ends; handles a wrapping arc
        arc = chrom.elements[i_s:i_e]
        del chrom.elements[i_s:i_e]
        self.chromosomes[new_cid] = Chromosome(new_cid, True, arc)
        return chrom.chrom_id, new_cid

    def fusion(self, rng, params) -> tuple[int, int]:
        """Fuse two chromosomes into one: append the second's segments onto the first and drop the
        second (segment ids untouched). The two are chosen uniformly; the caller guarantees >= 2
        chromosomes. Returns (kept_cid, absorbed_cid)."""
        cids = list(self.chromosomes)
        m = len(cids)
        i = int(rng.integers(m))
        j = int(rng.integers(m - 1))
        if j >= i:                            # a uniform partner distinct from i
            j += 1
        keep_cid, absorb_cid = cids[i], cids[j]
        self.chromosomes[keep_cid].elements.extend(self.chromosomes.pop(absorb_cid).elements)
        return keep_cid, absorb_cid
