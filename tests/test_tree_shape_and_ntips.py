"""Regression tests for the tree-shape contract and the n_tips floor (from the 2026-07 audit).

Two structural fixes:

* **tree-shape contract** — ``simulate_genomes`` rejects polytomous species trees cleanly (the Rust
  log/trace kernel used to panic on them); degree-two (FBD sampled-ancestor) nodes are handled
  uniformly (root pass-through instead of a crash; a *frozen* ancestral snapshot instead of an alias
  to the live genome; no cladogenetic jump); and cladogenesis works on multivariate continuous traits.
* **n_tips floor** — forward birth-death, gene-diversification, and SSE place the present strictly
  *after* the N-th lineage appears (last event + ``Exp(total rate)``), so no ``n_tips``-mode tree has a
  zero-length pendant edge, and ``n_tips == 2`` is no longer a degenerate zero-age tree.
"""

import numpy as np
import pytest

import zombi2 as z
from zombi2.tree import Tree, TreeNode, read_newick


# ============================================================ Pass 1: tree-shape contract

@pytest.mark.parametrize("output", ["profiles", "trace", "genomes"])
def test_polytomy_species_tree_rejected_cleanly(output):
    # ``P`` has three children — a hard polytomy. Every output mode must raise a clear ValueError
    # rather than crash the Rust kernel with a panic.
    tree = read_newick("((a:1,b:1,c:1)P:1,d:2)R:0;")
    with pytest.raises(ValueError, match="bifurcating"):
        z.simulate_genomes(tree, output=output, initial_families=5,
                           duplication=0.1, loss=0.1, seed=1)


def test_bifurcating_tree_still_simulates():
    tree = read_newick("((a:1,b:1)i:1,(c:1,d:1)j:1)R:0;")
    prof = z.simulate_genomes(tree, output="profiles", initial_families=5,
                              duplication=0.1, loss=0.1, seed=1)
    assert len(prof.families) == 5


def test_degree_two_root_does_not_crash():
    # a degree-two root (a sampled ancestor on the stem lineage) used to unpack-crash in _speciate
    tree = read_newick("((A:1,B:1)mid:0.5)root:0;")
    g = z.simulate_genomes(tree, duplication=0.1, loss=0.1, initial_families=3, seed=1)
    assert {leaf.name for leaf in g.leaf_genomes} == {"A", "B"}


def test_degree_two_ancestral_snapshot_is_frozen():
    # retain_internal at a degree-two node must store an independent, same-id snapshot rather than
    # alias the live genome (which keeps mutating along the child branch).
    from zombi2.genomes.genome_sim import GenomeSimulator

    root = TreeNode("r", 0.0)
    sa, leaf_c = TreeNode("sa", 0.3), TreeNode("C", 1.0)   # sa is degree-two (one child)
    root.add_child(sa)
    root.add_child(leaf_c)
    mid = TreeNode("mid", 0.6)
    sa.add_child(mid)
    mid.add_child(TreeNode("A", 1.0))
    mid.add_child(TreeNode("B", 1.0))
    tree = Tree(root, 1.0)

    res = GenomeSimulator().simulate(
        tree, z.SharedRates(duplication=0.8, transfer=0.0, loss=0.0, origination=0.0),
        np.random.default_rng(1), initial_size=5, retain_internal=True,
    )
    assert sa in res.node_genomes and mid in res.node_genomes
    # the sampled-ancestor snapshot is a distinct object from the continuing lineage's genome
    assert res.node_genomes[sa] is not res.node_genomes[mid]


def test_cladogenesis_handles_multivariate_continuous():
    tree = read_newick("((A:1,B:1):1,C:2);")
    res = z.simulate_traits(tree, z.MultivariateBrownian(np.eye(2)),
                            cladogenesis=z.Cladogenesis(jump_sigma2=0.5), seed=1)
    # no crash, and every node carries a length-2 vector (the jump is applied per dimension)
    assert all(np.asarray(v).shape == (2,) for v in res.node_values.values())


def test_cladogenesis_skips_degree_two_nodes():
    # S is a degree-two node; with ~zero anagenetic rate its single child C must inherit S's state
    # unchanged (a cladogenetic jump fires only at a real branching event).
    tree = read_newick("(((A:1,B:1)C:1)S:1,D:3)root:0;")
    res = z.simulate_traits(tree, z.Mk.equal_rates(4, 1e-9), root_state=0,
                            cladogenesis=z.Cladogenesis(shift=1.0), seed=1)
    nv = {n.name: res.node_values[n] for n in tree.nodes_preorder()}
    assert nv["C"] == nv["S"]


# ============================================================ Pass 2: n_tips floor

def _forward(model, **kw):
    return z.simulate_species_tree(model, direction="forward", **kw)


@pytest.mark.parametrize("model", [z.Yule(1.0), z.ClaDS(1.0)])   # thinning path + gillespie path
@pytest.mark.parametrize("seed", range(4))
def test_forward_ntips_has_positive_pendant_edges(model, seed):
    t = _forward(model, n_tips=8, seed=seed)
    assert len(t.extant_leaves()) == 8
    assert t.total_age > 0.0
    assert all(leaf.branch_length() > 0.0 for leaf in t.leaves())


@pytest.mark.parametrize("model", [z.Yule(1.0), z.ClaDS(1.0)])
def test_forward_ntips_two_is_nondegenerate(model):
    t = _forward(model, n_tips=2, seed=1)
    assert len(t.extant_leaves()) == 2
    assert t.total_age > 0.0
    assert all(leaf.branch_length() > 0.0 for leaf in t.leaves())


def test_gene_diversification_ntips_two_nondegenerate():
    m = z.GeneDiversification(2, driver_speciation=1.0, transfer=0.5, root_drivers=1)
    res = z.simulate_gene_diversification(m, n_tips=2, seed=1)
    assert len(res.tree.extant_leaves()) == 2
    assert res.tree.total_age > 0.0


def test_sse_ntips_two_nondegenerate():
    res = z.simulate_sse(z.BiSSE(1.5, 1.5, 0.2, 0.2, 0.1, 0.1), n_tips=2, seed=99)
    assert res.tree.total_age > 0.0
    assert all(leaf.branch_length() > 0.0 for leaf in res.tree.leaves())
