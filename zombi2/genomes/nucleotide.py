"""Genomes III — nucleotide: a genome as a coordinate space of base pairs.

The third and hardest resolution. Where an ordered genome is a *list of gene tokens*, a nucleotide
genome is a **coordinate space**: a chromosome of *N* base pairs, represented as an ordered list of
**segments**, each a *signed interval into an ancestral source* — ``Segment(source, start, end,
strand)`` mapping a physical stretch to the half-open interval ``[start, end)`` it descends from on
``source``, read forward (``+1``) or reverse-complemented (``-1``). Physical position is the running
sum of segment lengths; every structural event cuts at nucleotide coordinates and copies / moves /
reverses / deletes **arcs**.

This module is grown **slice by slice** (the level is ~1900 lines of coordinate arithmetic in the
legacy engine, so we build it very slowly). **Slice 1 — this file — is the representation and the one
event that changes no ancestry: inversion**, on a single chromosome that may be **circular** (the
canonical bacterial ring, where an inversion arc can wrap the origin) or **linear**. It stands up the
hard machinery — segment splitting at a coordinate, wrap-around arcs, and the ancestral trace-back —
via the safest event, with no genes, no genealogy, and no tree simulation yet. Deliberately deferred
to later slices: the tree/speciation wiring and a ``simulate_genomes_nucleotide``,
loss/duplication/transfer/transposition/origination (the events that create per-block gene trees),
translocation and the chromosome tier, indels, declared genes/intergenes (pseudogenization,
replacement), GFF input and BED/FASTA output.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass
class Segment:
    """A run of nucleotides with one ancestral origin: the half-open interval ``[start, end)`` on
    ``source`` this run descends from (``start < end`` always), read forward (``strand`` ``+1``) or
    reverse-complemented (``-1``) relative to that source. Its physical length is ``end - start``."""

    source: int
    start: int
    end: int
    strand: int

    @property
    def length(self) -> int:
        return self.end - self.start


def _split_segment(seg: Segment, o: int) -> list[Segment]:
    """Split ``seg`` after ``o`` physical positions into ``[left, right]``, strand-aware. For a
    forward segment the cut falls at ``start + o``; for a reversed one the first ``o`` physical
    positions are the **high** end of the source, so the cut falls at ``end - o``."""
    if seg.strand == 1:
        return [Segment(seg.source, seg.start, seg.start + o, 1),
                Segment(seg.source, seg.start + o, seg.end, 1)]
    return [Segment(seg.source, seg.end - o, seg.end, -1),
            Segment(seg.source, seg.start, seg.end - o, -1)]


@dataclass
class NucleotideGenome:
    """A single chromosome as an ordered list of :class:`Segment`\\ s over a nucleotide coordinate
    axis. Physical coordinate 0 is the origin; the total :attr:`length` is the sum of the segment
    lengths. ``topology`` is ``"circular"`` (a ring — the default, and the canonical bacterial case —
    where coordinate ``length`` wraps to ``0`` and an arc may cross the origin) or ``"linear"`` (two
    ends, no wrap). Slice 1 supports :meth:`split_at` (ensure a boundary) and :meth:`invert` (reverse
    an arc) plus the readers :meth:`mosaic` and :meth:`trace_back`."""

    segments: list[Segment]
    topology: str = "circular"

    def __post_init__(self) -> None:
        if self.topology not in ("circular", "linear"):
            raise ValueError(f"topology must be 'circular' or 'linear', got {self.topology!r}")

    @classmethod
    def seed(cls, root_length: int, topology: str = "circular", source: int = 0) -> "NucleotideGenome":
        """A fresh genome of ``root_length`` nucleotides: one segment ``[0, root_length)`` of
        ``source``, read forward, on a ``topology`` chromosome (circular by default)."""
        if root_length < 1:
            raise ValueError(f"root_length must be >= 1, got {root_length}")
        return cls([Segment(source, 0, root_length, 1)], topology)

    @property
    def length(self) -> int:
        return sum(s.length for s in self.segments)

    def split_at(self, c: int) -> None:
        """Ensure a segment boundary at physical coordinate ``c`` (``0 <= c <= length``). A no-op at
        ``0``, at ``length``, or where a boundary already exists; otherwise the segment straddling
        ``c`` is split in two."""
        if c <= 0:
            return
        pos = 0
        for i, seg in enumerate(self.segments):
            if pos == c:
                return                                   # already a boundary
            if pos < c < pos + seg.length:
                self.segments[i:i + 1] = _split_segment(seg, c - pos)
                return
            pos += seg.length
        # c == length (the right end): a boundary, no-op

    def _index_at(self, phys: int) -> int:
        """The segment index that starts at physical coordinate ``phys`` (which must be a boundary)."""
        pos = 0
        for i, seg in enumerate(self.segments):
            if pos == phys:
                return i
            pos += seg.length
        if pos == phys:
            return len(self.segments)                    # phys == length: past the last segment
        raise ValueError(f"no segment boundary at physical {phys}")

    def _arc_range(self, start: int, length: int) -> tuple[int, int] | None:
        """Prepare the arc ``[start, start+length)`` and return the segment index range ``[i, j)`` it
        occupies, or ``None`` for an empty arc. Splits at both ends first. On a **linear** chromosome
        the arc is clamped to the ends (no wrap). On a **circular** one it may cross the origin: there
        is no origin-fixed representation of a wrapping arc without a spurious breakpoint at 0, so the
        ring is rotated to bring the arc to the front — the origin drifts to a real breakpoint, which
        is harmless (a ring has no privileged origin and blocks live in origin-independent source
        coordinates)."""
        length_total = self.length
        if length_total == 0:
            return None
        if self.topology == "linear":
            start = max(0, min(start, length_total))
            end = min(start + max(0, length), length_total)
            if end <= start:
                return None
            self.split_at(start)
            self.split_at(end)
            return self._index_at(start), self._index_at(end)
        # circular
        ell = max(1, min(length, length_total))
        s = start % length_total
        self.split_at(s)
        self.split_at((s + ell) % length_total)   # no-op when the arc ends at the wrap or is the whole ring
        if s + ell <= length_total:                # non-wrapping: origin stays fixed
            j = self._index_at(s + ell) if s + ell < length_total else len(self.segments)
            return self._index_at(s), j
        i = self._index_at(s)                       # wrapping: rotate the arc to the front
        self.segments[:] = self.segments[i:] + self.segments[:i]
        j, acc = 0, 0
        while acc < ell:
            acc += self.segments[j].length
            j += 1
        return 0, j

    def invert(self, start: int, length: int) -> None:
        """Invert the arc ``[start, start+length)``: reverse the order of the segments in the arc and
        flip each one's strand. Ancestry is unchanged — an inversion only reshapes order and
        orientation — so the block genealogy (a later slice) is untouched. On a circular chromosome
        the arc may wrap the origin (see :meth:`_arc_range`)."""
        span = self._arc_range(start, length)
        if span is None:
            return
        i, j = span
        arc = self.segments[i:j]
        arc.reverse()
        for seg in arc:
            seg.strand = -seg.strand
        self.segments[i:j] = arc

    def mosaic(self) -> list[tuple[int, int, int, int]]:
        """The genome as ordered signed segments — ``[(source, start, end, strand), ...]`` in
        physical order (the coarse view; consecutive same-source runs are *not* merged)."""
        return [(s.source, s.start, s.end, s.strand) for s in self.segments]

    def trace_back(self) -> list[tuple[int, int, int]]:
        """Every nucleotide's ancestral origin, one entry per physical position left to right:
        ``(source, source_position, strand)``. A forward segment reads its source low→high; a
        reversed one reads it high→low. This is the exact provenance an ancestral-sequence
        reconstruction reads."""
        out: list[tuple[int, int, int]] = []
        for s in self.segments:
            if s.strand == 1:
                out.extend((s.source, p, 1) for p in range(s.start, s.end))
            else:
                out.extend((s.source, p, -1) for p in range(s.end - 1, s.start - 1, -1))
        return out


__all__ = ["NucleotideGenome", "Segment"]
