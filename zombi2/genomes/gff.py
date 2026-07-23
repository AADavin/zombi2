"""Reading a **GFF** — declaring the genes of an initial genome.

The nucleotide engine can be declared from a GFF3 instead of an evenly-spaced layout: each ``gene``
feature becomes a **declared, indivisible gene** at exactly those coordinates, and whatever lies
between genes becomes **intergene**. That is the "start from a real genome" path.

Only what declaring the genes needs is read:

- ``##sequence-region <seqid> <start> <end>`` — one **replicon** and its extent. When a GFF omits it,
  the replicon is taken to end at its last gene.
- feature lines whose **type** is ``gene`` — the declared genes.

GFF is **1-based inclusive**; blocks are **0-based half-open**, so coordinates are converted on the way
in (``start-1``, ``end``). ``strand`` becomes ``+1`` / ``-1`` (``.`` reads as ``+1``), so a gene declared
on the minus strand is laid down reverse-complemented. The gene's name is its ``ID`` attribute (else
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


def trim_overlapping_genes(genes: list[GffGene]) -> tuple[list[GffGene], int, list[str]]:
    """Shorten genes so that none overlap, and report what that cost.

    Real bacterial annotations **do** overlap — usually by a base or two where genes abut in an operon,
    occasionally by more. Here a gene is an *indivisible block* laid end to end, so overlaps have to go.
    Walking each replicon in coordinate order, a gene starting before the previous one ends has its
    **start pushed** to that end; one swallowed whole is dropped. Returns ``(genes, trimmed, dropped)``
    — the resolved genes, how many were shortened, and the names of any dropped."""
    out: list[GffGene] = []
    trimmed, dropped = 0, []
    end_of: dict[str, int] = {}
    for g in sorted(genes, key=lambda x: (x.seqid, x.start)):
        previous_end = end_of.get(g.seqid, 0)
        if g.start < previous_end:
            if g.end <= previous_end:                    # swallowed whole by its neighbour
                dropped.append(g.name)
                continue
            g = GffGene(g.seqid, previous_end, g.end, g.strand, g.name)
            trimmed += 1
        out.append(g)
        end_of[g.seqid] = g.end
    return out, trimmed, dropped


def read_gff(source, *, trim_overlaps: bool = False) -> tuple[dict[str, int], list[GffGene]]:
    """Read ``source`` (a path or an iterable of lines) and return ``({seqid: length}, [GffGene, …])``,
    the genes sorted by replicon and start.

    Genes are indivisible blocks laid end to end, so they may **touch but never overlap**. Real
    annotations do overlap, so pass ``trim_overlaps=True`` to shorten them instead of raising (see
    :func:`trim_overlapping_genes`).

    Raises :class:`ValueError` on a malformed line, a gene outside its replicon, or — unless
    ``trim_overlaps`` — two genes that overlap."""
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
    if trim_overlaps:
        genes, _trimmed, _dropped = trim_overlapping_genes(genes)
    else:
        for a, b in zip(genes, genes[1:]):               # laid end to end: may touch, never overlap
            if a.seqid == b.seqid and b.start < a.end:
                raise ValueError(f"genes {a.name!r} and {b.name!r} overlap on {a.seqid!r} "
                                 f"([{a.start}, {a.end}) and [{b.start}, {b.end})) — pass "
                                 f"trim_overlaps=True to shorten them instead")
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


def read_fasta(source) -> dict[str, str]:
    """Read ``source`` (a path or an iterable of lines) into ``{seqid: sequence}`` — the root DNA a
    nucleotide genome run starts from, paired with the GFF that lays its genes out.

    A record is a ``>seqid`` header (the id is its first whitespace-delimited token, matching a GFF
    ``##sequence-region``) followed by sequence lines, which are concatenated and upper-cased. Bases
    are validated as ``ACGT`` here — the one place letters enter, so a stray character fails loudly
    rather than surfacing as an evolved-sequence bug later."""
    if isinstance(source, (str, pathlib.Path)):
        lines = pathlib.Path(source).read_text().splitlines()
    else:
        lines = list(source)
    records: dict[str, list[str]] = {}
    seqid = None
    for n, raw in enumerate(lines, 1):
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            seqid = line[1:].split()[0]
            if seqid in records:
                raise ValueError(f"FASTA: sequence-region {seqid!r} appears twice")
            records[seqid] = []
        elif seqid is None:
            raise ValueError(f"FASTA line {n}: sequence before any '>' header: {raw!r}")
        else:
            records[seqid].append(line.upper())
    if not records:
        raise ValueError("the FASTA has no records")
    out = {}
    for sq, parts in records.items():
        seq = "".join(parts)
        bad = set(seq) - set("ACGT")
        if bad:
            raise ValueError(f"FASTA {sq!r} has non-ACGT characters {sorted(bad)} — a nucleotide "
                             "genome is given DNA (ambiguity codes and gaps are not supported)")
        out[sq] = seq
    return out


__all__ = ["GffGene", "read_gff", "read_fasta", "trim_overlapping_genes"]
