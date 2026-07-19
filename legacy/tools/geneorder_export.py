"""Gene-order export — derive gene-order study formats from a ``zombi2 genomes`` nucleotide run.

The complement of the fork's ``zombiExporter``. It consumes a genomes output directory produced by
``zombi2 genomes --genome-model nucleotide`` and turns the reconstructed genomes / event log into
the formats gene-order studies want. Inputs, by format:

- ``BED/<node>.bed`` (written by ``--write bed``) — each node's genes in genome order, with
  orientation. This *is* the reconstructed gene order at every node, so the order-based formats
  read it directly rather than re-simulating.
- ``geneorder_events.tsv`` (written by ``--write geneorder``) — the structural-event log with
  breakpoints, used for the event-count formats.
- ``species_tree.nwk`` — the tree, to walk parent→child edges.

Implemented: ``breakpoints`` (adjacencies broken on each tree edge), ``gff`` (every node's genes as
one GFF3) and ``posortho`` (positional ortholog sets over the leaves) — all from the per-node BED
gene orders. ``ffgc`` and a gene-level ``dupinfo`` are planned (see docs/design/geneorder-export.md:
dupinfo needs ``reconciliation_events.tsv`` to carry the user gene name rather than the internal
segment id, a small genomes-side fix).
"""

from __future__ import annotations

import os

from zombi2.tree import read_newick

BED_DIR = "bed"
SPECIES_TREE = "species_tree.nwk"


# --------------------------------------------------------------------------- #
# Reading the reconstructed per-node gene orders
# --------------------------------------------------------------------------- #
def read_node_beds(genomes_dir: str) -> dict[str, list[tuple[int, int, str, int]]]:
    """``{node: [(start, end, gene, strand), ...]}`` sorted by start, from ``BED/<node>.bed``.

    Each BED row is ``chrom  start  end  name  score  strand`` (BED6, 0-based half-open);
    ``strand`` is normalised to ``+1`` / ``-1``.
    """
    bed_dir = os.path.join(genomes_dir, BED_DIR)
    if not os.path.isdir(bed_dir):
        raise FileNotFoundError(
            f"no {BED_DIR}/ directory in {genomes_dir} — re-run 'zombi2 genomes' with 'bed' in "
            "--write (e.g. --write bed geneorder) so the per-node gene orders are available")
    beds: dict[str, list[tuple[int, int, str, int]]] = {}
    for fn in os.listdir(bed_dir):
        if not fn.endswith(".bed"):
            continue
        rows: list[tuple[int, int, str, int]] = []
        with open(os.path.join(bed_dir, fn)) as f:
            for line in f:
                if not line.strip():
                    continue
                p = line.rstrip("\n").split("\t")
                rows.append((int(p[1]), int(p[2]), p[3], 1 if p[5] == "+" else -1))
        rows.sort()
        beds[fn[:-4]] = rows
    return beds


def read_node_orders(genomes_dir: str) -> dict[str, list[tuple[str, int]]]:
    """``{node: [(gene, strand), ...]}`` in genome (position) order — the gene order at each node."""
    return {node: [(gene, strand) for _s, _e, gene, strand in rows]
            for node, rows in read_node_beds(genomes_dir).items()}


# --------------------------------------------------------------------------- #
# Breakpoints: adjacencies broken on each tree edge
# --------------------------------------------------------------------------- #
def _adjacencies(order: list[tuple[str, int]]) -> set[frozenset]:
    """The set of adjacencies of a circular signed gene order.

    Each gene has two extremities, head ``h`` and tail ``t``; read forward (``+``) a gene presents
    ``t`` then ``h``, reversed (``-``) it presents ``h`` then ``t``. An adjacency is the unordered
    pair of the extremities that meet between two consecutive genes (circular, so the last meets
    the first). Two genomes with the same gene content differ exactly in the adjacencies a
    rearrangement broke.
    """
    n = len(order)
    if n < 2:
        return set()

    def right(gene: str, strand: int) -> str:
        return f"{gene}_h" if strand == 1 else f"{gene}_t"

    def left(gene: str, strand: int) -> str:
        return f"{gene}_t" if strand == 1 else f"{gene}_h"

    adj = set()
    for i in range(n):
        g, s = order[i]
        g2, s2 = order[(i + 1) % n]
        adj.add(frozenset((right(g, s), left(g2, s2))))
    return adj


def broken_adjacencies(genomes_dir: str) -> list[tuple[str, str, list[str]]]:
    """``[(parent, child, [adjacency, ...]), ...]`` — adjacencies present in the parent genome but
    not the child's, i.e. the adjacencies broken by the rearrangements on that edge.

    Exact for content-conserving rearrangements (inversion / transposition). With duplication or
    loss the two genomes differ in gene content, so a differing adjacency may reflect a
    gained/lost gene rather than a broken adjacency — interpret those edges accordingly.
    """
    with open(os.path.join(genomes_dir, SPECIES_TREE)) as f:
        tree = read_newick(f.read())
    orders = read_node_orders(genomes_dir)
    out: list[tuple[str, str, list[str]]] = []
    for node in tree.nodes_preorder():
        if node.parent is None:
            continue
        parent, child = node.parent.name, node.name
        if parent not in orders or child not in orders:
            continue
        broken = _adjacencies(orders[parent]) - _adjacencies(orders[child])
        out.append((parent, child, sorted("|".join(sorted(a)) for a in broken)))
    return out


def breakpoints_tsv(genomes_dir: str) -> str:
    """``breakpoints.tsv`` text: one row per broken adjacency, ``parent  child  adjacency``."""
    rows = ["parent\tchild\tadjacency"]
    for parent, child, broken in broken_adjacencies(genomes_dir):
        for adj in broken:
            rows.append(f"{parent}\t{child}\t{adj}")
    return "\n".join(rows) + "\n"


# --------------------------------------------------------------------------- #
# gff: every node's genes as a single GFF3
# --------------------------------------------------------------------------- #
def gff_text(genomes_dir: str) -> str:
    """A single GFF3 of every node's genes — each node a seqid; BED 0-based → GFF 1-based inclusive.

    The gene family is the ``Name`` attribute; ``ID`` is per-node-unique (a family may recur after
    a duplication).
    """
    beds = read_node_beds(genomes_dir)
    lines = ["##gff-version 3"]
    for node in sorted(beds):
        rows = beds[node]
        if rows:
            lines.append(f"##sequence-region {node} 1 {max(e for _s, e, _g, _st in rows)}")
        for i, (start, end, gene, strand) in enumerate(rows):
            lines.append(f"{node}\tzombi2\tgene\t{start + 1}\t{end}\t.\t"
                         f"{'+' if strand == 1 else '-'}\t.\tID={node}.{i};Name={gene}")
    return "\n".join(lines) + "\n"


# --------------------------------------------------------------------------- #
# posortho: positional ortholog sets (genes grouped by ancestral family)
# --------------------------------------------------------------------------- #
def positional_orthologs(genomes_dir: str) -> dict[str, list[tuple[str, int, int]]]:
    """``{family: [(leaf, strand, start), ...]}`` over the extant leaves.

    Genes keep their ancestral identity across the tree, so every occurrence of one family is an
    ortholog (all descend from the same root gene). Exact positional orthology for
    content-conserving runs; under duplication a family has several copies per leaf and they are
    grouped together (family-level).
    """
    with open(os.path.join(genomes_dir, SPECIES_TREE)) as f:
        tree = read_newick(f.read())
    leaves = {n.name for n in tree.leaves()}
    sets: dict[str, list[tuple[str, int, int]]] = {}
    for node, rows in read_node_beds(genomes_dir).items():
        if node not in leaves:
            continue
        for start, _end, gene, strand in rows:
            sets.setdefault(gene, []).append((node, strand, start))
    for members in sets.values():
        members.sort()
    return sets


def posortho_tsv(genomes_dir: str) -> str:
    """``positional_orthologs.tsv``: ``family  leaf  strand  start`` — one row per gene occurrence."""
    sets = positional_orthologs(genomes_dir)
    rows = ["family\tleaf\tstrand\tstart"]
    for family in sorted(sets):
        for leaf, strand, start in sets[family]:
            rows.append(f"{family}\t{leaf}\t{'+' if strand == 1 else '-'}\t{start}")
    return "\n".join(rows) + "\n"
