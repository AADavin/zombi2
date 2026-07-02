"""RateVariation (Markov-modulated relaxed clock) and the Newick reader."""

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


def test_mean_rate_matches_weighted_bins():
    tree = simulate_species_tree(BirthDeath(1.0, 0.1), n_tips=200, age=8.0, seed=3)
    bins, weights = [0.5, 1.0, 2.0], [0.5, 0.3, 0.2]
    rv = RateVariation(bins=bins, switch_rate=2.0, weights=weights)
    scaled = rv.scale(tree, seed=3)
    total_time = sum(n.branch_length() for n in tree.nodes_preorder() if n.parent)
    total_subs = sum(scaled.branch_lengths[n] for n in tree.nodes_preorder() if n.parent)
    expected = sum(b * w for b, w in zip(bins, weights))  # stationary mean rate
    assert abs(total_subs / total_time - expected) < 0.05 * expected


def test_switching_splits_branches_into_segments():
    tree = simulate_species_tree(BirthDeath(1.0, 0.1), n_tips=30, age=6.0, seed=4)
    scaled = RateVariation(bins=[0.2, 5.0], switch_rate=5.0).scale(tree, seed=4)
    multi = [n for n in tree.nodes_preorder() if n.parent and len(scaled.segments[n]) > 1]
    assert multi  # at least some branches were split across bins


def test_reproducible():
    tree = simulate_species_tree(BirthDeath(1.0, 0.2), n_tips=12, age=4.0, seed=5)
    rv = RateVariation(bins=[0.5, 2.0], switch_rate=1.5)
    a = rv.scale(tree, seed=7).to_newick()
    b = rv.scale(tree, seed=7).to_newick()
    assert a == b


def test_newick_round_trip():
    tree = simulate_species_tree(BirthDeath(1.0, 0.3), n_tips=15, age=5.0, seed=6)
    nwk = tree.to_newick()
    back = read_newick(nwk)
    assert len(back.leaves()) == len(tree.leaves())
    assert {l.name for l in back.leaves()} == {l.name for l in tree.leaves()}
    assert back.to_newick() == nwk


def test_rate_variation_applies_to_gene_trees():
    tree = simulate_species_tree(BirthDeath(1.0, 0.2), n_tips=12, age=4.0, seed=8)
    g = simulate_genomes(tree, UniformRates(duplication=0.2, transfer=0.1, loss=0.2,
                                            origination=0.5), initial_size=12, seed=8)
    fam = max(g.profiles.families, key=lambda f: g.profiles.matrix[g.profiles.families.index(f)].sum())
    _, extant = g.gene_trees()[fam]
    gene_tree = read_newick(extant)                      # gene tree -> Tree
    scaled = RateVariation(bins=[0.5, 2.0], switch_rate=1.0).scale(gene_tree, seed=1)
    assert len(scaled.branch_lengths) == len(gene_tree.nodes())
    assert read_newick(scaled.to_newick()).leaves()      # phylogram is valid Newick
