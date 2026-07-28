"""Homology classification — the true ortholog / paralog / xenolog relation of every gene pair.

ZOMBI *simulated* each gene tree's embedding in the species tree, so the event at every internal node
is recorded, not inferred (see `zombi2.genomes.gene_trees`). Two genes are then related along
**two independent axes**, and a cell carries both:

**How they diverged** — the event at their most-recent common ancestor:

- a **speciation** — one gene, and the species carrying it split in two, so the two genes are
  **orthologs** (``O``);
- a **copying** — one gene became two, so they are **paralogs** (``P``). Either a duplication, which
  copies within a genome, or a transfer, which copies into another lineage; both leave two genes
  where there was one, which is what the letter says.

**Whether horizontal transfer is in their history** — an ``x`` suffix when a transfer sits anywhere on
the path from that common ancestor down to *either* gene, making them **xenologs** as well.

So a cell is ``O``, ``P``, ``Ox`` or ``Px``, and the diagonal (a gene against itself) is ``-``.

The two axes have to be separate, because they answer different questions and the answers are not
exclusive. Divergence is a property of one *node*; xenology is a property of the whole *path*, which
is how Fitch defined it — and a pair can perfectly well have diverged at a speciation and then had
one of its two genes carried somewhere else by a transfer. Reading xenology off the MRCA event alone,
as this table used to, catches only the pair whose divergence **is** the transfer — the copy left
behind and the copy that left — and calls every other transfer-affected pair a plain ortholog. On an
ordinary run that is a sixfold undercount, and it hides the case that matters most to anyone
benchmarking against these tables: two genes in the **same genome**, one of them an arrival from a
relative, whose MRCA is the speciation that separated the two species. That pair used to read ``O``,
which no orthology method can ever reproduce, so the answer key scored the method wrong. It now reads
``Ox``.

The result is one n×n symmetric table per family, n the extant leaves.
"""
from __future__ import annotations

import pathlib

from zombi2.genomes.events import copy_label
from zombi2.genomes.gene_trees import GeneNode, GeneTree

#: the MRCA's event → how the pair **diverged**. A transfer copies a gene into another lineage, so
#: like a duplication it leaves two genes where there was one: both are ``P``. What tells them apart
#: is the ``x`` suffix, which those pairs always carry (the transfer is on their path by definition).
_DIVERGENCE = {"speciation": "O", "duplication": "P", "transfer": "P"}


def _arrived(node: GeneNode) -> GeneNode | None:
    """Of a transfer's two children, the copy that **moved** — or ``None`` for any other node.

    A transfer ends one gene and starts two: the continuation on the donor's own branch, which never
    went anywhere, and the copy that landed somewhere else. Only the second has horizontal transfer in
    its history, so only its descendants are xenologous to the rest of the tree; a gene descending
    from the continuation is related to everything else by ordinary vertical descent, and marking it
    would make almost every pair in a transfer-rich family a xenolog.

    The one that moved is the one now on another branch. A self-transfer puts both on the same branch,
    and there it is the second child — the engine logs the donor's continuation first."""
    if node.kind != "transfer" or not node.children:
        return None
    return next((c for c in node.children if c.species != node.species), node.children[-1])


def _leaf_label(leaf: GeneNode) -> str:
    """A leaf's row/column header: ``n<species>_g<copy>`` — the species branch the gene sits on and its
    gene id. Both tokens are exactly the ones the gene-tree Newick and the event log write, so the
    table joins to either without translation."""
    return copy_label(leaf.species, leaf.copy)


def homology_table(root: GeneNode) -> tuple[list[str], list[list[str]]]:
    """Classify every pair of **extant** leaves under ``root``. Return ``(labels, matrix)``: ``labels``
    the leaf headers left-to-right, ``matrix`` the n×n grid of ``O`` / ``P`` / ``Ox`` / ``Px`` with
    ``-`` on the diagonal.

    Give it the family's **complete** root, not its extant one. The pair of leaves is the same either
    way and so is the divergence letter, but the ``x`` suffix is a fact about the *path*, and pruning
    to the survivors erases part of it: a transfer whose donor-side continuation left no extant
    descendant becomes a degree-two node and is suppressed, taking the record of the transfer with it.
    Measured at a fifth of all cells on an ordinary run.

    Each internal node is the MRCA of exactly those leaf pairs that first meet there — one leaf in one
    of its child subtrees, the other in a different one — so one bottom-up pass, carrying each node's
    set of descendant leaves, reaches every cell exactly once and reads its divergence off that node.

    The ``x`` suffix needs the whole path rather than the one node, but not a second walk: count the
    transfers **above** each node on the way down, and a pair has one on its path exactly when the two
    leaves' counts exceed twice their MRCA's. (A transfer at the MRCA itself is counted on both sides,
    so it lands in the total too — which is right, since that pair's divergence *is* a transfer.)
    """
    # a pre-order (parent before children); pushing children reversed makes pop() visit them
    # left-to-right, so the leaves come out in Newick order.
    order: list[GeneNode] = []
    stack = [root]
    while stack:
        n = stack.pop()
        order.append(n)
        stack.extend(reversed(n.children))

    above: dict[int, int] = {id(root): 0}        # node -> transfers strictly above it
    index: dict[int, int] = {}
    leaves: list[GeneNode] = []
    for n in order:                              # pre-order, so a node's own count is set already
        moved = _arrived(n)                      # only the copy that MOVED carries the transfer
        for c in n.children:
            above[id(c)] = above[id(n)] + (c is moved)
        if n.kind == "extant":                   # the observable leaves; a loss or a dead species
            index[id(n)] = len(leaves)           # is real history but nothing a dataset would hold
            leaves.append(n)
    labels = [_leaf_label(leaf) for leaf in leaves]
    t_leaf = [above[id(leaf)] for leaf in leaves]
    matrix = [["-"] * len(leaves) for _ in leaves]

    under: dict[int, list[int]] = {}                    # node -> the extant leaf indices beneath it
    for n in reversed(order):                           # children before parents
        if n.is_leaf:
            under[id(n)] = [index[id(n)]] if n.kind == "extant" else []
            continue
        how = _DIVERGENCE[n.kind]        # an internal node is a speciation, duplication or transfer
        here = 2 * above[id(n)]
        child_sets = [under[id(c)] for c in n.children]
        for a in range(len(child_sets)):
            for b in range(a + 1, len(child_sets)):
                for i in child_sets[a]:
                    for j in child_sets[b]:
                        rel = how + "x" if t_leaf[i] + t_leaf[j] > here else how
                        matrix[i][j] = matrix[j][i] = rel
        merged: list[int] = []
        for s in child_sets:
            merged.extend(s)
        under[id(n)] = merged
    return labels, matrix


def homology_tsv(root: GeneNode) -> str:
    """The homology table of the **complete** gene tree under ``root`` as TSV: a leading empty corner
    cell, then the leaf headers; one row per leaf, its label followed by the ``O`` / ``P`` / ``Ox`` / ``Px`` / ``-``
    cells. The top-left blank keeps it a square matrix a reader can load with the first column as the
    index."""
    labels, matrix = homology_table(root)
    header = "\t".join(["", *labels])
    rows = ["\t".join([labels[i], *matrix[i]]) for i in range(len(labels))]
    return "\n".join([header, *rows]) + "\n"


def write_homology(gene_trees: dict[int, GeneTree], tree, directory) -> str:
    """Write ``homology_fam<family>.tsv`` — one n×n O/P/Ox/Px table per family — into ``directory``.

    The rows are the **extant** leaves, the genes a real dataset would hold, so it mirrors
    `zombi2.genomes.gene_trees.write_gene_trees()`: a family with no surviving copy has no extant
    leaves to relate and so writes no table. It is read off the family's *complete* tree, though — see
    `homology_table()` for why the pruned one is not enough. (``tree``, the species tree, is part
    of the writer contract every ``--format`` shares; this table needs only the gene tree.)"""
    d = pathlib.Path(directory)
    d.mkdir(parents=True, exist_ok=True)
    n = 0
    for fam, gt in sorted(gene_trees.items()):
        if gt.extant is not None:                    # i.e. the family left at least one survivor
            (d / f"homology_fam{fam}.tsv").write_text(homology_tsv(gt.complete))
            n += 1
    return f"{n} table(s)"
