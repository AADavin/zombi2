"""Reading a **GFF** — declaring the genes of a seed genome.

The nucleotide engine can be seeded from a GFF3 instead of an evenly-spaced layout: each ``gene``
feature becomes a **declared, indivisible gene** at exactly those coordinates, and whatever lies
between genes becomes **intergene**. That is the "start from a real genome" path.

Only what the seeding needs is read:

- ``##sequence-region <seqid> <start> <end>`` — one **replicon** and its extent. When a GFF omits it,
  the replicon is taken to end at its last gene.
- feature lines whose **type** is ``gene`` — the declared genes.

GFF is **1-based inclusive**; blocks are **0-based half-open**, so coordinates are converted on the way
in (``start-1``, ``end``). ``strand`` becomes ``+1`` / ``-1`` (``.`` reads as ``+1``), so a gene declared
on the minus strand is seeded reverse-complemented. The gene's name is its ``ID`` attribute (else
``Name``, else a generated ``seqid:start-end``).
"""

from __future__ import annotations

import pathlib
from dataclasses import dataclass


@dataclass(frozen=True)
class GffGene:
    """One declared gene: ``[start, end)`` **0-based half-open** on replicon ``seqid``, read forward
    (``strand`` ``+1``) or reverse-complemented (``-1``), under the name ``name``."""

    seqid: str
    start: int
    end: int
    strand: int
    name: str


def _attribute(attrs: str, key: str) -> str | None:
    for field in attrs.split(";"):
        field = field.strip()
        if field.startswith(f"{key}="):
            return field[len(key) + 1:].strip() or None
    return None


def read_gff(source) -> tuple[dict[str, int], list[GffGene]]:
    """Read ``source`` (a path or an iterable of lines) and return ``({seqid: length}, [GffGene, …])``,
    the genes sorted by replicon and start.

    Raises :class:`ValueError` on a malformed line, a gene outside its replicon, or two genes that
    **overlap** — genes are indivisible blocks laid end to end, so they may touch but never overlap."""
    if isinstance(source, (str, pathlib.Path)):
        lines = pathlib.Path(source).read_text().splitlines()
    else:
        lines = list(source)

    lengths: dict[str, int] = {}
    genes: list[GffGene] = []
    for n, raw in enumerate(lines, 1):
        line = raw.strip()
        if not line:
            continue
        if line.startswith("#"):
            parts = line.lstrip("#").split()
            if parts and parts[0] == "sequence-region":
                if len(parts) != 4:
                    raise ValueError(f"line {n}: ##sequence-region needs <seqid> <start> <end>, got {line!r}")
                seqid, start, end = parts[1], int(parts[2]), int(parts[3])
                lengths[seqid] = end - start + 1
            continue
        cols = line.split("\t")
        if len(cols) < 8:
            raise ValueError(f"line {n}: a GFF feature needs at least 8 tab-separated columns, got {len(cols)}")
        seqid, _src, kind, start, end, _score, strand = cols[0], cols[1], cols[2], cols[3], cols[4], cols[5], cols[6]
        if kind.lower() != "gene":                       # only genes are declared; ignore other features
            continue
        try:
            start, end = int(start), int(end)
        except ValueError:
            raise ValueError(f"line {n}: start/end must be integers, got {cols[3]!r} {cols[4]!r}") from None
        if start < 1 or end < start:
            raise ValueError(f"line {n}: need 1 <= start <= end, got start={start} end={end}")
        attrs = cols[8] if len(cols) > 8 else ""
        name = _attribute(attrs, "ID") or _attribute(attrs, "Name") or f"{seqid}:{start}-{end}"
        genes.append(GffGene(seqid, start - 1, end, -1 if strand == "-" else 1, name))

    genes.sort(key=lambda g: (g.seqid, g.start))
    for a, b in zip(genes, genes[1:]):                   # laid end to end: they may touch, never overlap
        if a.seqid == b.seqid and b.start < a.end:
            raise ValueError(f"genes {a.name!r} and {b.name!r} overlap on {a.seqid!r} "
                             f"([{a.start}, {a.end}) and [{b.start}, {b.end}))")
    declared = set(lengths)                              # replicons given by ##sequence-region
    for g in genes:                                      # any other replicon ends at its last gene
        if g.seqid not in declared:
            lengths[g.seqid] = max(lengths.get(g.seqid, 0), g.end)
    for g in genes:
        if g.end > lengths[g.seqid]:
            raise ValueError(f"gene {g.name!r} ends at {g.end} beyond replicon {g.seqid!r} "
                             f"({lengths[g.seqid]} bp)")
    if not lengths:
        raise ValueError("the GFF declares no replicon and no gene")
    return lengths, genes


__all__ = ["GffGene", "read_gff"]
