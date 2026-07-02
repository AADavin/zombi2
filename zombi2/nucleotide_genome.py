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

from dataclasses import dataclass

import numpy as np

from .events import EventType, GeneOp, Region, Selection, TransferSegment
from .genome import Gene, Genome, IdManager


class SegmentRegistry:
    """Provenance shared across all genomes of one simulation (clones share one instance).

    ``provenance`` maps every segment id ever minted to its ancestral interval
    ``(source, src_start, src_end)`` — enough to attribute an event to atoms. ``split_parent``
    maps a segment born from a :meth:`NucleotideGenome._split_segment` to the segment it was
    cut from; those births are the only ones *not* in the event log (a split is a degree-2
    node), so recording them here keeps the segment genealogy connected for gene-tree
    reconstruction.
    """

    __slots__ = ("provenance", "split_parent")

    def __init__(self) -> None:
        self.provenance: dict[str, tuple[str, int, int]] = {}
        self.split_parent: dict[str, str] = {}


@dataclass(slots=True)
class Segment:
    """A maximal run of sequence that is contiguous *and* has one shared history.

    ``[src_start, src_end)`` is the half-open interval on ``source`` this run descends
    from (``src_start < src_end`` always); ``strand`` says whether it is read forward
    (+1) or reverse-complemented (-1) relative to that source.
    """

    seg_id: str
    source: str
    src_start: int
    src_end: int
    strand: int = 1

    @property
    def length(self) -> int:
        return self.src_end - self.src_start

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
                 extension: float | None = 0.99, registry: SegmentRegistry | None = None):
        self.ids = ids
        self.root_length = int(root_length)
        self.extension = extension
        self._registry: SegmentRegistry = registry if registry is not None else SegmentRegistry()
        self._segments: list[Segment] = []
        self._length = 0  # total nucleotide length; O(1), maintained on every event

    # --- minting ------------------------------------------------------------
    def _new_segment(self, source: str, a: int, b: int, strand: int) -> Segment:
        seg = Segment(self.ids.new_gene(), source, a, b, strand)
        self._registry.provenance[seg.seg_id] = (source, a, b)
        return seg

    # --- queries ------------------------------------------------------------
    def size(self) -> int:
        return self._length  # nucleotide length -> inversion rate scales with genome size

    def total_length(self) -> float:
        return float(self._length)

    def n_segments(self) -> int:
        return len(self._segments)

    def families(self) -> list[str]:
        return list(dict.fromkeys(s.source for s in self._segments))

    def copy_number(self, family: str) -> int:
        return sum(1 for s in self._segments if s.source == family)

    def genes(self) -> list[Gene]:
        return [Gene(s.seg_id, s.source) for s in self._segments]

    def presence_vector(self, family_order) -> np.ndarray:
        present = {s.source for s in self._segments}
        return np.fromiter((1 if f in present else 0 for f in family_order),
                           dtype=np.int8, count=len(family_order))

    def supported_events(self) -> frozenset[EventType]:
        return frozenset((EventType.ORIGINATION, EventType.INVERSION, EventType.LOSS,
                          EventType.DUPLICATION, EventType.TRANSFER, EventType.TRANSPOSITION))

    # --- trace-back view: one cell per nucleotide, in physical order --------
    def to_cells(self) -> list[tuple[str, int, int]]:
        """Expand to ``[(source, src_pos, strand), ...]`` of length ``size()``.

        This *is* the position-level trace-back: cell ``p`` names the ancestral
        nucleotide the present-day position ``p`` descends from. Used by the tests as
        the ground truth and by the post-processor.
        """
        out: list[tuple[str, int, int]] = []
        for seg in self._segments:
            if seg.strand == 1:
                out.extend((seg.source, p, 1) for p in range(seg.src_start, seg.src_end))
            else:
                out.extend((seg.source, p, -1) for p in range(seg.src_end - 1, seg.src_start - 1, -1))
        return out

    # --- seeding / origination ---------------------------------------------
    def originate(self, rng, params) -> list[GeneOp]:
        """Create brand-new sequence under a fresh ``source`` namespace.

        The first call (empty genome, the seed) lays down the full-length root chromosome.
        Later calls insert a *novel* gene of geometric length at a random position — it
        descends from no ancestral position, so it opens its own source and its own atoms,
        with the gene tree rooted at this origination time.
        """
        source = self.ids.new_family()
        if not self._segments:                      # seed: the root chromosome
            seg = self._new_segment(source, 0, self.root_length, 1)
            self._segments.append(seg)
        else:                                       # novel gene inserted somewhere
            length = self._draw_length(rng, params)
            seg = self._new_segment(source, 0, length, 1)
            at = int(rng.integers(self._length + 1))
            if at >= self._length:
                self._segments.append(seg)
            else:
                self._split_at(at)
                self._segments.insert(self._index_at(at), seg)
        self._length += seg.length
        return [GeneOp(seg.seg_id, source, "origin")]

    # --- the split primitive (the whole novelty vs OrderedGenome) ----------
    def _split_segment(self, seg: Segment, o: int) -> tuple[Segment, Segment]:
        """Split ``seg`` after ``o`` physical positions into (left, right), strand-aware."""
        if seg.strand == 1:
            left = self._new_segment(seg.source, seg.src_start, seg.src_start + o, 1)
            right = self._new_segment(seg.source, seg.src_start + o, seg.src_end, 1)
        else:  # reversed: first o physical positions are the HIGH end of the source
            left = self._new_segment(seg.source, seg.src_end - o, seg.src_end, -1)
            right = self._new_segment(seg.source, seg.src_start, seg.src_end - o, -1)
        # a split is a degree-2 birth (unlogged); record it so the genealogy stays connected
        self._registry.split_parent[left.seg_id] = seg.seg_id
        self._registry.split_parent[right.seg_id] = seg.seg_id
        return left, right

    def _split_at(self, c: int) -> None:
        """Ensure a segment boundary at physical coordinate ``c`` (``0`` is always one)."""
        if c == 0:
            return
        pos = 0
        for i, seg in enumerate(self._segments):
            if pos == c:
                return  # already a boundary
            seglen = seg.length
            if pos < c < pos + seglen:
                self._segments[i:i + 1] = list(self._split_segment(seg, c - pos))
                return
            pos += seglen
        # c == self._length: a boundary (wrap point); no-op

    def _index_at(self, phys: int) -> int:
        pos = 0
        for i, seg in enumerate(self._segments):
            if pos == phys:
                return i
            pos += seg.length
        raise ValueError(f"no segment boundary at physical {phys}")

    # --- inversion ----------------------------------------------------------
    def _apply_inversion(self, s: int, ell: int) -> list[Segment]:
        """Invert the circular arc ``[s, s+ell)`` and return the reversed segments.

        A non-wrapping arc is reversed in place with the physical origin fixed. A
        **wrapping** arc has no origin-fixed representation without inventing a spurious
        breakpoint at position 0, so the ring is first rotated to bring the arc to the
        front — the origin drifts to a real breakpoint. This is harmless: a circular
        genome has no privileged origin, atoms live in origin-independent *source*
        coordinates, and the invariants are rotation-invariant.
        """
        L = self._length
        if L == 0:
            return []
        ell = max(1, min(ell, L))
        s %= L
        e = (s + ell) % L
        self._split_at(s)
        self._split_at(e)  # no-op when e == 0 (arc ends at the wrap) or e == s (whole genome)

        if s + ell <= L:  # --- non-wrapping: origin fixed ---
            i_s = self._index_at(s)
            i_e = self._index_at(s + ell) if s + ell < L else len(self._segments)
        else:             # --- wrapping: rotate the arc to the front, origin drifts ---
            i_s = self._index_at(s)
            self._segments = self._segments[i_s:] + self._segments[:i_s]
            i_s, i_e, acc = 0, 0, 0
            while acc < ell:
                acc += self._segments[i_e].length
                i_e += 1

        arc = self._segments[i_s:i_e]
        arc.reverse()                     # segment order flips...
        for seg in arc:
            seg.strand = -seg.strand       # ...and each segment's strand flips
        self._segments[i_s:i_e] = arc
        return arc

    # --- deletion -----------------------------------------------------------
    def _apply_loss(self, s: int, ell: int) -> list[Segment]:
        """Delete the circular arc ``[s, s+ell)`` and return the removed segments.

        Boundaries at ``s`` and the arc end are split first, so the arc is a whole
        number of segments; the survivors keep their circular order (the physical origin
        follows the surviving material when the deletion swallows position 0).
        """
        L = self._length
        if L == 0:
            return []
        ell = min(ell, L)
        s %= L
        if ell == L:                       # deletes the whole genome
            removed = self._segments
            self._segments, self._length = [], 0
            return removed
        e = (s + ell) % L
        self._split_at(s)
        self._split_at(e)

        wrapping = s + ell > L
        removed, kept, pos = [], [], 0
        for seg in self._segments:
            start = pos
            pos += seg.length
            in_arc = (start >= s or pos <= e) if wrapping else (start >= s and pos <= s + ell)
            (removed if in_arc else kept).append(seg)
        self._segments = kept
        self._length -= ell
        return removed

    # --- duplication --------------------------------------------------------
    def _arc_range(self, s: int, ell: int) -> tuple[int, int]:
        """Split the ends of arc ``[s, s+ell)`` and return its ``[i_s, i_e)`` segment range.

        A wrapping arc is rotated to the front first (the origin drifts to a real
        breakpoint), so the arc is always a contiguous slice ``self._segments[i_s:i_e]``.
        """
        L = self._length
        ell = max(1, min(ell, L))
        s %= L
        e = (s + ell) % L
        self._split_at(s)
        self._split_at(e)
        if s + ell <= L:
            i_s = self._index_at(s)
            i_e = self._index_at(s + ell) if s + ell < L else len(self._segments)
            return i_s, i_e
        i_s = self._index_at(s)
        self._segments = self._segments[i_s:] + self._segments[:i_s]
        i_e, acc = 0, 0
        while acc < ell:
            acc += self._segments[i_e].length
            i_e += 1
        return 0, i_e

    def _apply_duplication(self, s: int, ell: int) -> list[list[GeneOp]]:
        """Duplicate the arc ``[s, s+ell)`` in tandem; return one op-group per segment.

        Each arc segment forks into a re-minted *continuation* (in place) and a *copy*
        (the paralog, inserted immediately after the arc). Both carry the segment's source
        interval, so the two are paralogs that coalesce at this duplication.
        """
        L = self._length
        if L == 0:
            return []
        ell = max(1, min(ell, L))
        s %= L
        i_s, i_e = self._arc_range(s, ell)
        arc = self._segments[i_s:i_e]
        groups, conts, copies = [], [], []
        for a in arc:
            cont = self._new_segment(a.source, a.src_start, a.src_end, a.strand)
            copy = self._new_segment(a.source, a.src_start, a.src_end, a.strand)
            conts.append(cont)
            copies.append(copy)
            groups.append([GeneOp(a.seg_id, a.source, "parent"),
                           GeneOp(cont.seg_id, a.source, "left"),
                           GeneOp(copy.seg_id, a.source, "right")])
        self._segments[i_s:i_e] = conts          # continuations replace the originals...
        self._segments[i_e:i_e] = copies          # ...and the copy block lands in tandem
        self._length += ell
        return groups

    # --- transposition (cut a segment, paste it elsewhere) -----------------
    def _apply_transposition(self, s: int, ell: int, dest: int) -> list[Segment]:
        """Cut the arc ``[s, s+ell)`` and splice it back in at physical ``dest``.

        Content- and length-preserving and genealogically neutral (segments keep their
        ids), so it only permutes the mosaic — like inversion, no lineage re-mint.
        """
        L = self._length
        if L <= 1:
            return []
        ell = max(1, min(ell, L - 1))
        i_s, i_e = self._arc_range(s, ell)
        arc = self._segments[i_s:i_e]
        del self._segments[i_s:i_e]
        rem = L - ell
        dest %= rem + 1
        if dest >= rem:
            self._segments.extend(arc)
        else:
            self._split_at(dest)
            idx = self._index_at(dest)
            self._segments[idx:idx] = arc
        return arc

    def _draw_length(self, rng, params) -> int:
        ext = params.extension if params.extension is not None else self.extension
        n = self._length
        if ext is None or n <= 1:
            return 1
        if ext >= 1.0:
            return n
        # truncated geometric (mean 1/(1-ext)); drawn in one shot instead of a Python loop
        return min(int(rng.geometric(1.0 - ext)), n)

    _TARGETABLE = frozenset((EventType.INVERSION, EventType.LOSS, EventType.DUPLICATION,
                             EventType.TRANSFER, EventType.TRANSPOSITION))

    def draw_target(self, event, rng, params, family=None) -> Selection:
        if event not in self._TARGETABLE:
            raise ValueError(f"NucleotideGenome does not target {event!r}")
        s = int(rng.integers(self._length))
        ell = self._draw_length(rng, params)
        return Selection(genes=(), region=Region(chromosome=0, start=s, length=ell))

    def apply(self, event, selection, rng, params) -> list[list[GeneOp]]:
        region = selection.region
        if event is EventType.INVERSION:
            arc = self._apply_inversion(region.start, region.length)
            return [[GeneOp(seg.seg_id, seg.source, "inverted", orientation=seg.strand)
                     for seg in arc]]
        if event is EventType.LOSS:
            removed = self._apply_loss(region.start, region.length)
            return [[GeneOp(seg.seg_id, seg.source, "lost")] for seg in removed]
        if event is EventType.DUPLICATION:
            return self._apply_duplication(region.start, region.length)
        if event is EventType.TRANSPOSITION:
            if self._length <= 1:
                return []
            dest = int(rng.integers(self._length))
            arc = self._apply_transposition(region.start, region.length, dest)
            return [[GeneOp(seg.seg_id, seg.source, "transposed") for seg in arc]]
        raise ValueError(f"NucleotideGenome does not handle {event!r}")

    # --- transfer handoff ---------------------------------------------------
    def extract_segment(self, selection, rng) -> TransferSegment:
        """Donor side: fork the arc — re-mint a continuation in place, build a copy to send.

        Each arc segment ``g`` forks into a continuation (stays in the donor) and a copy
        (travels to the recipient). The driver logs ``g -> [continuation, copy]`` as a
        TRANSFER, a bifurcation the atom's gene tree reads as a transfer node. The copy
        keeps the arc's ancestral coordinates, so it is a xenolog of the same atom.
        """
        region = selection.region
        i_s, i_e = self._arc_range(region.start, region.length)
        arc = self._segments[i_s:i_e]
        old_gids, cont_gids, copies, conts = [], [], [], []
        for a in arc:
            old_gids.append(a.seg_id)
            cont = self._new_segment(a.source, a.src_start, a.src_end, a.strand)
            conts.append(cont)
            cont_gids.append(cont.seg_id)
            copies.append(self._new_segment(a.source, a.src_start, a.src_end, a.strand))
        self._segments[i_s:i_e] = conts          # donor length unchanged (arc -> continuations)
        return TransferSegment(family=None, genes=tuple(copies),
                               donor_old_gids=old_gids, donor_cont_gids=cont_gids)

    def choose_insertion_point(self, segment, rng) -> int:
        return int(rng.integers(self._length + 1))  # any physical position, 0..length

    def insert_segment(self, segment, at, rng) -> list[GeneOp]:
        """Recipient side: splice the transferred copy block in at physical position ``at``."""
        copies = list(segment.genes)
        if at is None or at >= self._length:
            idx = len(self._segments)
        else:
            self._split_at(at)
            idx = self._index_at(at)
        self._segments[idx:idx] = copies
        self._length += sum(seg.length for seg in copies)
        return [GeneOp(seg.seg_id, seg.source, "transfer_copy") for seg in copies]

    # --- speciation ---------------------------------------------------------
    def clone_reminting(self) -> tuple["NucleotideGenome", list[tuple[str, str, str]]]:
        child = NucleotideGenome(self.ids, self.root_length, self.extension, self._registry)
        child._length = self._length
        mapping: list[tuple[str, str, str]] = []
        for seg in self._segments:
            ns = child._new_segment(seg.source, seg.src_start, seg.src_end, seg.strand)
            child._segments.append(ns)
            mapping.append((seg.seg_id, ns.seg_id, seg.source))
        return child, mapping
