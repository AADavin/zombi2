"""RateVariation (autocorrelated nearest-neighbour relaxed clock) and the Newick reader."""

import numpy as np

from zombi2 import (
    BirthDeath,
    RateVariation,
    UniformRates,
    read_newick,
    simulate_genomes,
    simulate_species_tree,
)


def test_strict_clock_scales_every_branch_uniformly():
    tree = simulate_species_tree(BirthDeath(1.0, 0.2), n_tips=10, age=4.0, seed=1)
    scaled = RateVariation(bins=[2.0], switch_rate=0.0).scale(tree, seed=1)
    for node in tree.nodes_preorder():
        if node.parent is not None:
            assert abs(scaled.branch_lengths[node] - 2.0 * node.branch_length()) < 1e-9


def test_single_unit_bin_preserves_lengths():
    tree = simulate_species_tree(BirthDeath(1.0, 0.2), n_tips=8, age=3.0, seed=2)
    scaled = RateVariation(bins=[1.0], switch_rate=0.0).scale(tree, seed=2)
    assert scaled.to_newick() == tree.to_newick()


def test_transitions_are_only_to_adjacent_bins():
    """The defining property: within a branch, the bin only ever steps by +/-1."""
    tree = simulate_species_tree(BirthDeath(1.0, 0.1), n_tips=40, age=8.0, seed=3)
    scaled = RateVariation(bins=[0.25, 0.5, 1.0, 2.0, 4.0], switch_rate=4.0).scale(tree, seed=3)
    switches = 0
    for node in tree.nodes_preorder():
        segs = scaled.segments.get(node, [])
        for (b0, _), (b1, _) in zip(segs, segs[1:]):
            assert abs(b1 - b0) == 1  # never jumps across bins
            switches += 1
    assert switches > 0  # some switching actually happened


def test_symmetric_walk_mean_rate_matches_bin_mean():
    """A symmetric nearest-neighbour walk has a uniform stationary distribution."""
    tree = simulate_species_tree(BirthDeath(1.0, 0.1), n_tips=200, age=8.0, seed=4)
    bins = [0.5, 1.0, 2.0]
    scaled = RateVariation(bins=bins, switch_rate=6.0).scale(tree, seed=4)
    total_time = sum(n.branch_length() for n in tree.nodes_preorder() if n.parent)
    total_subs = sum(scaled.branch_lengths[n] for n in tree.nodes_preorder() if n.parent)
    assert abs(total_subs / total_time - np.mean(bins)) < 0.1 * np.mean(bins)


def test_up_bias_drives_the_chain_upward():
    tree = simulate_species_tree(BirthDeath(1.0, 0.1), n_tips=60, age=8.0, seed=5)
    top = 4
    scaled = RateVariation(bins=[0.25, 0.5, 1.0, 2.0, 4.0], switch_rate=3.0,
                           up_bias=1.0, start=0).scale(tree, seed=5)
    # always stepping up (from the bottom) -> most branches end near the top bin
    ends = [scaled.end_bin[n] for n in tree.nodes_preorder() if n.parent]
    assert np.mean(ends) > top / 2


def test_reproducible():
    tree = simulate_species_tree(BirthDeath(1.0, 0.2), n_tips=12, age=4.0, seed=6)
    rv = RateVariation(bins=[0.5, 1.0, 2.0], switch_rate=1.5)
    assert rv.scale(tree, seed=7).to_newick() == rv.scale(tree, seed=7).to_newick()


def test_newick_round_trip():
    tree = simulate_species_tree(BirthDeath(1.0, 0.3), n_tips=15, age=5.0, seed=8)
    nwk = tree.to_newick()
    back = read_newick(nwk)
    assert len(back.leaves()) == len(tree.leaves())
    assert {l.name for l in back.leaves()} == {l.name for l in tree.leaves()}
    assert back.to_newick() == nwk


def test_rate_variation_applies_to_gene_trees():
    tree = simulate_species_tree(BirthDeath(1.0, 0.2), n_tips=12, age=4.0, seed=9)
    g = simulate_genomes(tree, UniformRates(duplication=0.2, transfer=0.1, loss=0.2,
                                            origination=0.5), initial_families=12, seed=9)
    fam = max(g.profiles.families, key=lambda f: g.profiles.matrix[g.profiles.families.index(f)].sum())
    _, extant = g.gene_trees()[fam]
    gene_tree = read_newick(extant)                        # gene tree -> Tree
    scaled = RateVariation(bins=[0.5, 1.0, 2.0], switch_rate=1.0).scale(gene_tree, seed=1)
    assert len(scaled.branch_lengths) == len(gene_tree.nodes())
    assert read_newick(scaled.to_newick()).leaves()        # phylogram is valid Newick
