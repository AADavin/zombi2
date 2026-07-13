"""Read a genome's size and gene coordinates from a GFF3 annotation.

The nucleotide genic model (:func:`~zombi2.simulate_nucleotide_genomes`) takes a chromosome
length and a list of non-overlapping gene intervals. :func:`read_gff` extracts both from a
real GFF3 file — e.g. a RefSeq bacterial chromosome — so a simulation can start from an actual
genome's architecture (its length, its genes, and the intergenes between them).

GFF is 1-based inclusive; ZOMBI2 uses 0-based half-open, so a GFF gene ``[start, end]`` becomes
``[start-1, end)``. Bacterial genes sometimes overlap (shared start/stop codons, nested ORFs);
because the genic model forbids breakpoints inside a gene, overlaps are removed by trimming: the
genes are sorted, each gene's start is clipped up to the running end of the kept genes, and a gene
swallowed whole is dropped. The number trimmed/dropped is reported on the result.
"""

from __future__ import annotations

import gzip
from dataclasses import dataclass


@dataclass
class GffGenome:
    """A chromosome's architecture read from a GFF: its length and non-overlapping genes.

    ``genes`` is ``[(start, end, name), ...]`` in 0-based half-open coordinates, sorted and
    guaranteed non-overlapping — ready to pass as ``gene_intervals`` (with ``root_length=length``)
    to :func:`~zombi2.simulate_nucleotide_genomes`. ``n_features`` is how many gene features were
    read for this sequence; ``n_trimmed`` / ``n_dropped`` record the overlap cleanup.
    """

    seqid: str
    length: int
    circular: bool
    genes: list[tuple[int, int, str]]
    n_features: int
    n_trimmed: int
    n_dropped: int


def _open(path: str):
    """Open a plain or gzipped GFF as text."""
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path)


def _attrs(field: str) -> dict[str, str]:
    return dict(kv.split("=", 1) for kv in field.split(";") if "=" in kv)


def _gene_name(attrs: dict[str, str], used: set[str], fallback: int) -> str:
    """A stable, unique gene id: locus_tag > Name > ID (minus a ``gene-`` prefix) > geneN."""
    name = attrs.get("locus_tag") or attrs.get("Name") or attrs.get("gene")
    if not name:
        gid = attrs.get("ID", "")
        name = gid[5:] if gid.startswith("gene-") else gid
    if not name:
        name = f"gene{fallback}"
    if name in used:                      # disambiguate a repeated annotation name
        k = 2
        while f"{name}.{k}" in used:
            k += 1
        name = f"{name}.{k}"
    used.add(name)
    return name


def _parse_gff(path, types: set):
    """Parse a GFF3 into per-sequence raw features. Returns ``(region_len, feat_len, circular,
    genes)`` — the last a ``seqid -> [(start, end, attr), ...]`` map (1-based, as in the file)."""
    region_len: dict[str, int] = {}       # from ##sequence-region pragmas
    feat_len: dict[str, int] = {}         # from `region` features
    circular: dict[str, bool] = {}
    genes: dict[str, list[tuple[int, int, str]]] = {}

    with _open(path) as fh:
        for line in fh:
            if line.startswith("##FASTA"):
                break                     # sequence section — no more features
            if line.startswith("##sequence-region"):
                p = line.split()
                if len(p) >= 4:
                    region_len[p[1]] = int(p[3])
                continue
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            sid, _src, typ, start, end, _score, _strand, _phase, attr = f[:9]
            if typ == "region":
                feat_len[sid] = max(feat_len.get(sid, 0), int(end))
                circular[sid] = _attrs(attr).get("Is_circular", "").lower() == "true"
            elif typ in types:
                genes.setdefault(sid, []).append((int(start), int(end), attr))

    if not genes:
        raise ValueError(f"no {'/'.join(sorted(types))} features found in {path!r}")
    return region_len, feat_len, circular, genes


def _build_genome(seqid, gene_list, region_len, feat_len, circular) -> GffGenome:
    """Assemble one :class:`GffGenome` from a single sequence's raw gene features."""
    length = region_len.get(seqid) or feat_len.get(seqid) or max(e for _s, e, _a in gene_list)

    # 0-based half-open, with a stable unique name; validate against the chromosome length
    used: set[str] = set()
    raw: list[tuple[int, int, str]] = []
    for i, (s, e, attr) in enumerate(sorted(gene_list)):
        a, b = s - 1, e
        if not (0 <= a < b <= length):
            raise ValueError(f"gene [{s},{e}] on {seqid} falls outside [1, {length}]")
        raw.append((a, b, _gene_name(_attrs(attr), used, i)))

    # trim overlaps: clip each start up to the running end of kept genes; drop the swallowed
    kept: list[tuple[int, int, str]] = []
    max_end = 0
    n_trimmed = n_dropped = 0
    for a, b, name in raw:
        if a < max_end:
            if b <= max_end:
                n_dropped += 1
                continue
            a = max_end
            n_trimmed += 1
        kept.append((a, b, name))
        max_end = b

    return GffGenome(seqid=seqid, length=length, circular=circular.get(seqid, False),
                     genes=kept, n_features=len(gene_list),
                     n_trimmed=n_trimmed, n_dropped=n_dropped)


def read_gff(path, *, seqid: str | None = None, feature_types=("gene",)) -> GffGenome:
    """Read one chromosome's length + non-overlapping gene intervals from a GFF3 file.

    Parameters
    ----------
    path:
        A GFF3 file (optionally gzipped, ``.gz``).
    seqid:
        Which sequence to read. Default: the sequence carrying the most gene features (the
        chromosome of a single-chromosome bacterium; a plasmid is skipped). Raises if the named
        ``seqid`` is absent. Use :func:`read_gff_all` to read *every* sequence as a chromosome.
    feature_types:
        GFF feature types treated as genes (default ``("gene",)``; e.g. add ``"pseudogene"`` to
        include annotated pseudogenes as blocks too).

    Returns a :class:`GffGenome`. The length comes from the ``##sequence-region`` pragma, else the
    ``region`` feature, else the largest gene end. Overlapping genes are trimmed (later start
    clipped) or dropped (swallowed whole); the counts are on the result.
    """
    region_len, feat_len, circular, genes = _parse_gff(path, set(feature_types))
    if seqid is None:
        seqid = max(genes, key=lambda s: len(genes[s]))   # the most-annotated sequence
    elif seqid not in genes:
        raise ValueError(f"seqid {seqid!r} not found (have: {', '.join(sorted(genes))})")
    return _build_genome(seqid, genes[seqid], region_len, feat_len, circular)


def read_gff_all(path, *, feature_types=("gene",)) -> list[GffGenome]:
    """Read **every** sequence in a GFF3 as its own :class:`GffGenome`, most-annotated first.

    A real genome file usually holds several sequences — a main chromosome plus plasmids or
    secondary chromosomes/contigs. This returns one :class:`GffGenome` per sequence (ordered by
    gene count, so the main chromosome leads), ready to seed a **multi-chromosome** nucleotide
    genome (pass as ``root_chromosomes`` to :func:`~zombi2.simulate_nucleotide_genomes`).
    """
    region_len, feat_len, circular, genes = _parse_gff(path, set(feature_types))
    order = sorted(genes, key=lambda s: (-len(genes[s]), s))   # main chromosome first, then plasmids
    return [_build_genome(s, genes[s], region_len, feat_len, circular) for s in order]
