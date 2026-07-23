"""Tests for :mod:`zombi2.tools.homology` — the ortholog / paralog / xenolog classifier.

The relation of two genes is the event at their most-recent common ancestor: a speciation makes them
orthologs (``O``), a duplication paralogs (``P``), a transfer xenologs (``X``). Because ZOMBI records
that event on every gene-tree node, the table is exact, so these tests pin the mapping on hand-built
trees and then cross-check the whole matrix against a naive pairwise-LCA computation on real runs."""

import pytest

from zombi2.genomes.gene_trees import GeneNode, gene_trees_from_events
from zombi2.genomes import simulate_genomes_unordered
from zombi2.species import simulate_species_tree
from zombi2.tools.homology import homology_table, homology_tsv


def _leaf(species: int, copy: int) -> GeneNode:
    return GeneNode("extant", species, 1.0, copy)


def _internal(kind: str, species: int, children: list[GeneNode]) -> GeneNode:
    n = GeneNode(kind, species, 0.5, -1)
    n.children = children
    return n


def test_mrca_event_maps_to_the_relation():
    # a tree whose three internal nodes are one of each kind, so every relation appears once:
    #   root=duplication( speciation(g10,g11) , transfer(g20,g21) )
    left = _internal("speciation", 2, [_leaf(3, 10), _leaf(4, 11)])
    right = _internal("transfer", 5, [_leaf(5, 20), _leaf(6, 21)])
    root = _internal("duplication", 1, [left, right])

    labels, m = homology_table(root)
    assert labels == ["n3|g10", "n4|g11", "n5|g20", "n6|g21"]     # left-to-right, Newick order
    idx = {lab: i for i, lab in enumerate(labels)}

    def rel(a, b):
        return m[idx[a]][idx[b]]

    assert rel("n3|g10", "n4|g11") == "O"                         # MRCA is the speciation
    assert rel("n5|g20", "n6|g21") == "X"                         # MRCA is the transfer
    for a in ("n3|g10", "n4|g11"):                                # MRCA of any cross pair is the root
        for b in ("n5|g20", "n6|g21"):
            assert rel(a, b) == "P"                               # ...a duplication


def test_matrix_is_symmetric_with_a_dashed_diagonal():
    root = _internal("duplication", 1,
                     [_internal("speciation", 2, [_leaf(3, 10), _leaf(4, 11)]), _leaf(5, 20)])
    _, m = homology_table(root)
    n = len(m)
    assert all(m[i][i] == "-" for i in range(n))
    assert all(m[i][j] == m[j][i] for i in range(n) for j in range(n))


def test_a_single_leaf_is_a_one_by_one_table():
    labels, m = homology_table(_leaf(7, 99))
    assert labels == ["n7|g99"]
    assert m == [["-"]]


def test_tsv_is_a_square_grid_with_a_blank_corner():
    root = _internal("speciation", 1, [_leaf(2, 10), _leaf(3, 11)])
    text = homology_tsv(root)
    lines = text.rstrip("\n").split("\n")
    assert lines[0] == "\tn2|g10\tn3|g11"                         # blank corner, then the headers
    assert lines[1] == "n2|g10\t-\tO"
    assert lines[2] == "n3|g11\tO\t-"


def _naive_matrix(root: GeneNode):
    """A deliberately simple oracle: leaf order + parent pointers, then walk both leaves to the root
    and take the first shared ancestor. Slower and obviously correct — the check the fast set-based
    :func:`homology_table` must agree with."""
    relation = {"speciation": "O", "duplication": "P", "transfer": "X"}
    leaves, parent, stack = [], {}, [root]
    while stack:
        n = stack.pop()
        if n.is_leaf:
            leaves.append(n)
        for c in reversed(n.children):
            parent[id(c)] = n
            stack.append(c)

    def ancestors(n):
        chain = [n]
        while id(chain[-1]) in parent:
            chain.append(parent[id(chain[-1])])
        return chain

    m = [["-"] * len(leaves) for _ in leaves]
    for a in range(len(leaves)):
        seen = {id(x) for x in ancestors(leaves[a])}
        for b in range(a + 1, len(leaves)):
            lca = next(x for x in ancestors(leaves[b]) if id(x) in seen)
            m[a][b] = m[b][a] = relation[lca.kind]
    return m


@pytest.mark.parametrize("seed", range(1, 8))
def test_matches_a_naive_pairwise_lca_on_real_runs(seed):
    sp = simulate_species_tree(birth=1, death=0.3, n_extant=10, seed=seed)
    g = simulate_genomes_unordered(sp.complete_tree, duplication=0.4, transfer=0.3,
                                   loss=0.25, origination=0.6, seed=seed * 3)
    trees = gene_trees_from_events(g.events, g.complete_tree)
    checked = 0
    for gt in trees.values():
        root = gt.extant
        if root is None:
            continue
        _, fast = homology_table(root)
        assert fast == _naive_matrix(root)
        checked += 1
    assert checked > 0                                           # the run actually had extant families
