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
event that changes no ancestry: inversion.** It stands up the hard machinery (segment splitting at a
coordinate, and the ancestral trace-back) via the safest event, on a **single linear chromosome**,
with no genes, no genealogy, and no tree simulation yet. Deliberately deferred to later slices:
circular topology (wrap-around arcs), the tree/speciation wiring and a ``simulate_genomes_nucleotide``,
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
    """A single **linear** chromosome as an ordered list of :class:`Segment`\\ s over a nucleotide
    coordinate axis. Physical coordinate 0 is the left end; the total :attr:`length` is the sum of
    the segment lengths. Slice 1 supports two operations — :meth:`split_at` (ensure a boundary) and
    :meth:`invert` (reverse an arc) — plus two readers, :meth:`mosaic` and :meth:`trace_back`."""

    segments: list[Segment]

    @classmethod
    def seed(cls, root_length: int, source: int = 0) -> "NucleotideGenome":
        """A fresh genome of ``root_length`` nucleotides: one segment ``[0, root_length)`` of
        ``source``, read forward."""
        if root_length < 1:
            raise ValueError(f"root_length must be >= 1, got {root_length}")
        return cls([Segment(source, 0, root_length, 1)])

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

    def invert(self, start: int, length: int) -> None:
        """Invert the arc ``[start, start+length)``: split at both ends, reverse the order of the
        segments in the arc, and flip each one's strand. Ancestry is unchanged — an inversion only
        reshapes order and orientation — so the block genealogy (a later slice) is untouched. The arc
        is clamped to the chromosome (no wrap: this is a linear chromosome for now)."""
        start = max(0, min(start, self.length))
        end = min(start + max(0, length), self.length)
        if end <= start:
            return
        self.split_at(start)
        self.split_at(end)
        i, j = self._index_at(start), self._index_at(end)
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
