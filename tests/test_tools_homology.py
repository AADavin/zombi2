"""Tests for :mod:`zombi2.tools.homology` — the event at every gene pair's common ancestor.

Two genes are related on two axes. **How they diverged** is the event at their most-recent common
ancestor: ``S`` speciation, ``D`` duplication, ``T`` transfer. **Whether transfer is in their history
since** is a fact about the whole path: an ``x`` suffix when a transfer sits between the common
ancestor and either gene — never on a ``T``, whose divergence *is* a transfer. The table states the
event rather than reading it as ortholog or paralog, because the readings disagree; see the Tools
appendix. Because ZOMBI records every event on the tree, both axes are exact, so these tests pin them
on hand-built trees and then cross-check the whole matrix against a naive pairwise-LCA oracle."""

import pytest

from zombi2.genomes.gene_trees import GeneNode, gene_trees_from_events
from zombi2.genomes import simulate_genomes_family
from zombi2.species import simulate_species_tree
from zombi2.tools.homology import homology_table, homology_tsv


@pytest.fixture(scope="module")
def seeded_run():
    """A run with enough transfer and loss to exercise both axes."""
    sp = simulate_species_tree(birth=1.0, death=0.3, n_extant=14, seed=1)
    return simulate_genomes_family(sp, duplication=0.3, transfer=0.3, loss=0.3, origination=0.5,
                                   initial_families=20, seed=1)


def _leaf(species: int, copy: int) -> GeneNode:
    return GeneNode("extant", species, 1.0, copy)


def _internal(kind: str, species: int, children: list[GeneNode]) -> GeneNode:
    n = GeneNode(kind, species, 0.5, -1)
    n.children = children
    return n


def test_mrca_event_maps_to_the_relation():
    # a tree whose three internal nodes are one of each kind, so every relation appears once:
    #   root=duplication( speciation(n3_g10,n4_g11) , transfer(n5_g20,n6_g21) )
    left = _internal("speciation", 2, [_leaf(3, 10), _leaf(4, 11)])
    right = _internal("transfer", 5, [_leaf(5, 20), _leaf(6, 21)])
    root = _internal("duplication", 1, [left, right])

    labels, m = homology_table(root)
    assert labels == ["n3_g10", "n4_g11", "n5_g20", "n6_g21"]     # left-to-right, Newick order
    idx = {lab: i for i, lab in enumerate(labels)}

    def rel(a, b):
        return m[idx[a]][idx[b]]

    assert rel("n3_g10", "n4_g11") == "S"                         # MRCA is the speciation
    # the transfer's own two children — the divergence IS the transfer, so no x for it
    assert rel("n5_g20", "n6_g21") == "T"
    for a in ("n3_g10", "n4_g11"):                                # MRCA of any cross pair is the root
        for b in ("n5_g20", "n6_g21"):
            # ...a duplication — and n6_g21 is the copy that moved, so pairs with it carry the x
            assert rel(a, b) == ("Dx" if b == "n6_g21" else "D")


def test_a_transferred_gene_is_a_xenolog_of_everything_it_meets():
    # The case the old table got wrong. Species ((a,b),c); c donates a copy of its gene into a, so a
    # ends up holding two genes. Every pair involving the ARRIVAL has transfer in its history; no
    # other pair does — including the copy c kept, which never went anywhere.
    a1 = _leaf(3, 1)                                  # a's own gene (species n3)
    b1 = _leaf(4, 2)                                  # b's gene       (species n4)
    c1 = _leaf(2, 3)                                  # the copy c KEPT (species n2)
    a2 = _leaf(3, 4)                                  # the copy c SENT to a
    root = _internal("speciation", 0, [_internal("speciation", 1, [a1, b1]),
                                       _internal("transfer", 2, [c1, a2])])
    labels, m = homology_table(root)
    rel = {(labels[i], labels[j]): m[i][j] for i in range(4) for j in range(4)}

    assert rel[("n3_g1", "n4_g2")] == "S"             # both vertical, diverged at a speciation
    assert rel[("n3_g1", "n2_g3")] == "S"             # the kept copy never moved, so this is too
    assert rel[("n4_g2", "n2_g3")] == "S"
    assert rel[("n2_g3", "n3_g4")] == "T"             # kept vs sent: the transfer itself
    assert rel[("n3_g1", "n3_g4")] == "Sx"            # SAME GENOME — used to read a plain O
    assert rel[("n4_g2", "n3_g4")] == "Sx"


def test_two_genes_in_one_genome_never_read_as_a_plain_speciation(seeded_run):
    # they are in one genome, so a duplication or a transfer put them there — the divergence may
    # still be a speciation (a gene that left and came back), but then transfer is in the history
    plain_O = 0
    for gt in seeded_run.gene_trees.values():
        if gt.extant is None:
            continue
        labels, m = homology_table(gt.complete)
        for i in range(len(labels)):
            for j in range(i + 1, len(labels)):
                if labels[i].split("_")[0] == labels[j].split("_")[0]:
                    assert m[i][j] != "S", (labels[i], labels[j])
                    plain_O += m[i][j] == "Sx"
    assert plain_O > 0                                # the run really did transfer within a genome


def test_pruning_would_lose_transfers_so_the_complete_tree_is_used(seeded_run):
    # a transfer whose donor-side copy left no extant descendant is a degree-two node in the pruned
    # tree and is suppressed, taking the record of the transfer with it. Measured at a fifth of all
    # cells, which is why write_homology reads the complete tree.
    differs = same = 0
    for gt in seeded_run.gene_trees.values():
        if gt.extant is None:
            continue
        full_labels, full = homology_table(gt.complete)
        cut_labels, cut = homology_table(gt.extant)
        assert full_labels == cut_labels              # the same genes either way...
        for i in range(len(full_labels)):
            for j in range(i + 1, len(full_labels)):
                differs += full[i][j] != cut[i][j]    # ...but not the same history
                same += full[i][j] == cut[i][j]
    assert differs > 0.05 * (differs + same), f"only {differs}/{differs + same} cells differ"


def test_matrix_is_symmetric_with_a_dashed_diagonal():
    root = _internal("duplication", 1,
                     [_internal("speciation", 2, [_leaf(3, 10), _leaf(4, 11)]), _leaf(5, 20)])
    _, m = homology_table(root)
    n = len(m)
    assert all(m[i][i] == "-" for i in range(n))
    assert all(m[i][j] == m[j][i] for i in range(n) for j in range(n))


def test_a_single_leaf_is_a_one_by_one_table():
    labels, m = homology_table(_leaf(7, 99))
    assert labels == ["n7_g99"]
    assert m == [["-"]]


def test_tsv_is_a_square_grid_with_a_blank_corner():
    root = _internal("speciation", 1, [_leaf(2, 10), _leaf(3, 11)])
    text = homology_tsv(root)
    lines = text.rstrip("\n").split("\n")
    assert lines[0] == "\tn2_g10\tn3_g11"                         # blank corner, then the headers
    assert lines[1] == "n2_g10\t-\tS"
    assert lines[2] == "n3_g11\tS\t-"


def _naive_matrix(root: GeneNode):
    """A deliberately simple oracle: leaf order + parent pointers, then walk both leaves to the root
    and take the first shared ancestor. Slower and obviously correct — the check the fast set-based
    :func:`homology_table` must agree with."""
    relation = {"speciation": "S", "duplication": "D", "transfer": "T"}
    leaves, parent, moved, stack = [], {}, set(), [root]
    while stack:
        n = stack.pop()
        if n.kind == "extant":
            leaves.append(n)
        if n.kind == "transfer" and n.children:      # which child is the copy that went somewhere
            arrived = next((c for c in n.children if c.species != n.species), n.children[-1])
            moved.add(id(arrived))
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
        chain_a = ancestors(leaves[a])
        seen = {id(x) for x in chain_a}
        for b in range(a + 1, len(leaves)):
            chain_b = ancestors(leaves[b])
            lca = next(x for x in chain_b if id(x) in seen)
            # walk each leaf up to the LCA, counting the transfers it descends through
            # a transfer AT the divergence is what T already says, so its own arrival child does
            # not count; the suffix marks a transfer SINCE
            skip = {id(c) for c in lca.children} if lca.kind == "transfer" else set()
            hgt = any(id(x) in moved and id(x) not in skip
                      for chain in (chain_a, chain_b) for x in chain[:chain.index(lca)])
            m[a][b] = m[b][a] = relation[lca.kind] + ("x" if hgt else "")
    return m


@pytest.mark.parametrize("seed", range(1, 8))
def test_matches_a_naive_pairwise_lca_on_real_runs(seed):
    sp = simulate_species_tree(birth=1, death=0.3, n_extant=10, seed=seed)
    g = simulate_genomes_family(sp.complete_tree, duplication=0.4, transfer=0.3,
                                loss=0.25, origination=0.6, seed=seed * 3)
    trees = gene_trees_from_events(g.events, g.complete_tree)
    checked = 0
    for gt in trees.values():
        if gt.extant is None:
            continue
        _, fast = homology_table(gt.complete)      # the COMPLETE tree: pruning erases transfers
        assert fast == _naive_matrix(gt.complete)
        checked += 1
    assert checked > 0                                           # the run actually had extant families
