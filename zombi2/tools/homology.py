"""Homology classification — the true ortholog / paralog / xenolog relation of every gene pair.

ZOMBI *simulated* each gene tree's embedding in the species tree, so the event at every internal node
is recorded, not inferred (see :mod:`zombi2.genomes.gene_trees`). The relation between two genes is
then read straight off the event at their most-recent common ancestor:

- a **speciation** at the MRCA — the two genes descend from one gene a speciation split between
  daughter species, so they are **orthologs** (``O``);
- a **duplication** — they descend from a gene a duplication copied within a genome, so they are
  **paralogs** (``P``);
- a **transfer** — they descend from a gene a transfer carried between lineages, so they are
  **xenologs** (``X``).

This is the standard event-at-the-LCA scheme, and it is exact here because the history is known rather
than reconstructed. The result is one n×n table per family (n the extant leaves): a symmetric grid
whose off-diagonal cells are ``O`` / ``P`` / ``X`` and whose diagonal (a gene against itself) is ``-``.
"""
from __future__ import annotations

import pathlib

from zombi2.genomes.events import node_label
from zombi2.genomes.gene_trees import GeneNode, GeneTree

#: the internal-node event kind → the one-letter relation of the leaf pairs it is the MRCA of.
_RELATION = {"speciation": "O", "duplication": "P", "transfer": "X"}


def _leaf_label(leaf: GeneNode) -> str:
    """A leaf's row/column header: ``n<species>|g<id>`` — the species branch the gene sits on and its
    gene id. Both tokens are exactly the ones the gene-tree Newick and the event log write, so the
    table joins to either without translation."""
    return f"{node_label(leaf.species)}|g{leaf.copy}"


def homology_table(root: GeneNode) -> tuple[list[str], list[list[str]]]:
    """Classify every pair of leaves under ``root``. Return ``(labels, matrix)``: ``labels`` the leaf
    headers left-to-right, ``matrix`` the n×n grid of ``O`` / ``P`` / ``X`` with ``-`` on the diagonal.

    Each internal node is the MRCA of exactly those leaf pairs that first meet there — one leaf in one
    of its child subtrees, the other in a different one — so a single bottom-up pass, carrying each
    node's set of descendant leaves, fills every cell with the relation its MRCA's event implies.
    """
    # a pre-order (parent before children); pushing children reversed makes pop() visit them
    # left-to-right, so the leaves come out in Newick order.
    order: list[GeneNode] = []
    stack = [root]
    while stack:
        n = stack.pop()
        order.append(n)
        stack.extend(reversed(n.children))

    index: dict[int, int] = {}
    leaves: list[GeneNode] = []
    for n in order:
        if n.is_leaf:
            index[id(n)] = len(leaves)
            leaves.append(n)
    labels = [_leaf_label(leaf) for leaf in leaves]
    matrix = [["-"] * len(leaves) for _ in leaves]

    under: dict[int, list[int]] = {}                    # node -> the leaf indices beneath it
    for n in reversed(order):                           # children before parents
        if n.is_leaf:
            under[id(n)] = [index[id(n)]]
            continue
        rel = _RELATION[n.kind]                          # extant-tree internals are S / D / T only
        child_sets = [under[id(c)] for c in n.children]
        for a in range(len(child_sets)):
            for b in range(a + 1, len(child_sets)):
                for i in child_sets[a]:
                    for j in child_sets[b]:
                        matrix[i][j] = matrix[j][i] = rel
        merged: list[int] = []
        for s in child_sets:
            merged.extend(s)
        under[id(n)] = merged
    return labels, matrix


def homology_tsv(root: GeneNode) -> str:
    """The homology table of the gene tree under ``root`` as TSV: a leading empty corner cell, then the
    leaf headers; one row per leaf, its label followed by the ``O`` / ``P`` / ``X`` / ``-`` cells. The
    top-left blank keeps it a square matrix a reader can load with the first column as the index."""
    labels, matrix = homology_table(root)
    header = "\t".join(["", *labels])
    rows = ["\t".join([labels[i], *matrix[i]]) for i in range(len(labels))]
    return "\n".join([header, *rows]) + "\n"


def write_homology(gene_trees: dict[int, GeneTree], directory) -> None:
    """Write ``homology_fam<family>.tsv`` — one n×n O/P/X table per family — into ``directory``.

    The table is over the **extant** gene tree, the leaves a real dataset would hold, so it mirrors
    :func:`zombi2.genomes.gene_trees.write_gene_trees`: a family with no surviving copy has no extant
    leaves to relate and so writes no table."""
    d = pathlib.Path(directory)
    d.mkdir(parents=True, exist_ok=True)
    for fam, gt in sorted(gene_trees.items()):
        root = gt.extant
        if root is not None:
            (d / f"homology_fam{fam}.tsv").write_text(homology_tsv(root))
