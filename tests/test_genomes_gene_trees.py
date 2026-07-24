"""Tests for gene trees — each family's true genealogy inside the complete tree (zombi2.genomes)."""

import pytest

from zombi2.species import simulate_species_tree
from zombi2.tree import Node, Tree
from zombi2.genomes import GeneTree, simulate_genomes_unordered

LEAF_KINDS = {"extant", "extinct", "unsampled", "loss"}
INTERNAL_KINDS = {"duplication", "transfer", "speciation"}   # ZOMBI1: a node's kind is what ended the gene


def _leaves(node, kinds=None):
    if node is None:
        return []
    if node.is_leaf:
        return [node] if (kinds is None or node.kind in kinds) else []
    return [lf for c in node.children for lf in _leaves(c, kinds)]


def _all_nodes(node):
    if node is None:
        return []
    return [node] + [n for c in node.children for n in _all_nodes(c)]


def _run(seed=1, n_extant=15, death=0.4, **kw):
    sp = simulate_species_tree(birth=1.0, death=death, n_extant=n_extant, seed=seed)
    params = dict(duplication=0.4, transfer=0.2, loss=0.3, origination=0.7, initial_families=6, seed=seed)
    params.update(kw)                                # caller overrides (e.g. transfer=0.0, loss=1.2)
    return sp, simulate_genomes_unordered(sp, **params)


# --- the tree exists and is well-formed ------------------------------------

def test_one_gene_tree_per_family():
    sp, g = _run(seed=2)
    assert set(g.gene_trees) == {e.family for e in g.events if e.kind == "origination"}
    assert all(isinstance(t, GeneTree) for t in g.gene_trees.values())


def test_every_node_is_annotated_with_a_known_kind():
    _, g = _run(seed=3)
    for tree in g.gene_trees.values():
        for n in _all_nodes(tree.complete):
            assert n.kind in (LEAF_KINDS | INTERNAL_KINDS)
            assert (n.kind in LEAF_KINDS) == n.is_leaf     # leaf kinds are leaves, internal kinds branch
            assert isinstance(n.species, int) and isinstance(n.copy, int)


def test_root_is_the_founding_gene():
    # ZOMBI1: no separate origination node — the root IS the founding gene (its kind is what ended it)
    _, g = _run(seed=4)
    for fam, tree in g.gene_trees.items():
        origin_copy = next(e.copy for e in g.events if e.kind == "origination" and e.family == fam)
        assert tree.complete.copy == origin_copy


def test_origination_is_the_exact_event_time_and_starts_the_root_branch():
    # A GeneNode records when it *ended*, so the founding gene's start is the one time the tree
    # cannot derive. It comes straight off the origination event — bit-for-bit, not rounded — and is
    # where the root's branch begins, so the stem is not silently dropped from the Newick.
    _, g = _run(seed=4)
    for fam, tree in g.gene_trees.items():
        event_time = next(e.time for e in g.events if e.kind == "origination" and e.family == fam)
        assert tree.origination == event_time                    # exact, no tolerance

        stem = tree.complete.time - tree.origination
        assert stem >= 0.0
        written = float(tree.to_newick("complete").rsplit(":", 1)[1].rstrip(";"))
        assert written == pytest.approx(stem, rel=1e-5)

        # the extant root may have had ancestors suppressed; its branch still starts at origination
        if tree.extant is not None:
            extant_stem = float(tree.to_newick("extant").rsplit(":", 1)[1].rstrip(";"))
            assert extant_stem == pytest.approx(tree.extant.time - tree.origination, rel=1e-5)


def test_a_family_that_never_splits_is_one_node_with_its_lifespan():
    # the degenerate tree: a lone gene, written as its own branch rather than a bare bald label
    _, g = _run(seed=4)
    singles = [t for t in g.gene_trees.values() if t.complete.is_leaf]
    if not singles:                                              # seed-dependent; skip if none arose
        pytest.skip("no single-gene family under this seed")
    t = singles[0]
    assert t.to_newick("complete") == f"g{t.complete.copy}:{t.complete.time - t.origination:.6g};"


# --- the key cross-check: gene tree agrees with the profiles ---------------

def test_extant_leaves_equal_the_extant_copy_total():
    # the strongest invariant: a family's surviving gene-tree leaves == its copies across extant tips
    sp, g = _run(seed=5, death=0.5)
    extant_sp = {n.id for n in sp.complete_tree.extant()}
    for fam, tree in g.gene_trees.items():
        n_leaves = len(_leaves(tree.extant, {"extant"}))
        copies = sum(g.profiles.counts.get((fam, s), 0) for s in extant_sp)
        assert n_leaves == copies


def test_every_extant_leaf_sits_on_an_extant_species():
    sp, g = _run(seed=6)
    extant_sp = {n.id for n in sp.complete_tree.extant()}
    for tree in g.gene_trees.values():
        assert all(lf.species in extant_sp for lf in _leaves(tree.extant, {"extant"}))


def test_extant_is_none_iff_the_family_left_no_extant_copy():
    sp, g = _run(seed=3, death=0.7)                 # high death -> some families fully extinct
    extant_sp = {n.id for n in sp.complete_tree.extant()}
    saw_none = False
    for fam, tree in g.gene_trees.items():
        copies = sum(g.profiles.counts.get((fam, s), 0) for s in extant_sp)
        assert (tree.extant is None) == (copies == 0)
        saw_none |= tree.extant is None
    assert saw_none                                  # the high-death run really does strand some families


# --- events place correctly in the tree ------------------------------------

def test_transfer_children_cross_species_branches():
    # a transfer's two children (donor continuation + recipient copy) sit on different species branches
    _, g = _run(seed=7, self_transfer=False)
    xfers = [n for tree in g.gene_trees.values() for n in _all_nodes(tree.complete) if n.kind == "transfer"]
    assert xfers                                     # the run actually transferred
    assert all(a.species != b.species for n in xfers for a, b in [n.children])


def test_self_transfer_is_allowed_on_the_same_branch_and_still_balances():
    # self-transfer lands on the donor's own branch, so a transfer node CAN have same-branch children;
    # the extant-leaf-count invariant must still hold
    sp, g = _run(seed=2, self_transfer=True)
    extant_sp = {n.id for n in sp.complete_tree.extant()}
    for fam, tree in g.gene_trees.items():
        assert len(_leaves(tree.extant, {"extant"})) == sum(g.profiles.counts.get((fam, s), 0) for s in extant_sp)


def test_no_transfer_nodes_when_transfer_is_zero():
    _, g = _run(seed=8, transfer=0.0)
    kinds = {n.kind for tree in g.gene_trees.values() for n in _all_nodes(tree.complete)}
    assert "transfer" not in kinds
    assert "duplication" in kinds                    # but the other events still happen


def test_complete_tree_keeps_lost_and_dead_lineages():
    _, g = _run(seed=6, loss=1.2)                    # heavy loss
    leaf_kinds = {lf.kind for tree in g.gene_trees.values() for lf in _leaves(tree.complete)}
    assert "loss" in leaf_kinds                      # the complete tree records losses as leaves


# --- serialisation + determinism -------------------------------------------

def test_newick_is_balanced_and_both_trees_serialise():
    _, g = _run(seed=9)
    fam = next(f for f, t in g.gene_trees.items() if t.extant is not None)
    for which in ("complete", "extant"):
        nw = g.gene_trees[fam].to_newick(which)
        assert nw.endswith(";") and nw.count("(") == nw.count(")")


def test_deterministic_given_seed():
    sp, g = _run(seed=3)
    g2 = simulate_genomes_unordered(sp, duplication=0.4, transfer=0.2, loss=0.3, origination=0.7,
                                    initial_families=6, seed=3)
    assert all(g.gene_trees[f].to_newick("complete") == g2.gene_trees[f].to_newick("complete")
               for f in g.gene_trees)


def test_empty_run_has_no_gene_trees():
    sp = simulate_species_tree(birth=1.0, death=0.3, n_extant=10, seed=1)
    g = simulate_genomes_unordered(sp, initial_families=0, seed=1)       # no families
    assert g.gene_trees == {}


def test_deep_tree_serialises_without_recursion_error():
    # a 700-node pectinate spine -> a ~700-deep gene tree; the serializer must be iterative, since
    # CPython's C-stack recursion guard (~500 levels) can't be lifted by setrecursionlimit
    L = 700
    nodes = {}
    for k in range(L):
        s, leaf, nxt = 2 * k, 2 * k + 1, 2 * (k + 1)
        nodes[s] = Node(s, (2 * (k - 1) if k > 0 else None), float(k), float(k + 1), (leaf, nxt), "speciation")
        nodes[leaf] = Node(leaf, s, float(k + 1), float(L + 1), None, "extant")
    nodes[2 * L] = Node(2 * L, 2 * (L - 1), float(L), float(L + 1), None, "extant")
    g = simulate_genomes_unordered(Tree(nodes, 0), initial_families=1, seed=0)
    for which in ("complete", "extant"):
        nw = g.gene_trees[0].to_newick(which)         # must not raise RecursionError
        assert nw.endswith(";") and nw.count("(") == nw.count(")")
