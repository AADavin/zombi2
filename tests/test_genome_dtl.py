"""Verification for Algorithm 2 (forward D/T/L/O gene families)."""

import math

import numpy as np

from zombi2 import (
    BirthDeath,
    GenomeSimulator,
    UniformRates,
    simulate_genomes,
    simulate_species_tree,
)
from zombi2.events import EventType


def _name_to_node(tree):
    return {n.name: n for n in tree.nodes_preorder()}


def _full_run(seed):
    tree = simulate_species_tree(BirthDeath(1.0, 0.2), n_tips=12, age=3.0, seed=seed)
    genomes = simulate_genomes(
        tree, duplication=0.15, transfer=0.1, loss=0.2, origination=0.4,
        initial_families=8, seed=seed,
    )
    return tree, genomes


def test_transfer_time_consistency():
    tree = simulate_species_tree(BirthDeath(1.0, 0.2), n_tips=15, age=4.0, seed=5)
    genomes = simulate_genomes(
        tree, duplication=0.1, transfer=0.3, loss=0.15, origination=0.3,
        initial_families=15, seed=5,
    )
    name2node = _name_to_node(tree)
    n_transfers = 0
    for r in genomes.event_log:
        if r.event is not EventType.TRANSFER:
            continue
        n_transfers += 1
        donor, recipient = name2node[r.donor], name2node[r.recipient]
        assert r.donor != r.recipient
        for node in (donor, recipient):
            assert node.parent is not None
            assert node.parent.time < r.time <= node.time + 1e-9
    assert n_transfers > 0


def test_profile_shape_and_binarization():
    tree = simulate_species_tree(BirthDeath(1.0, 0.2), n_tips=10, age=3.0, seed=1)
    genomes = simulate_genomes(tree, duplication=0.1, transfer=0.05, loss=0.15,
                               origination=0.3, initial_families=10, seed=1)
    P = genomes.profiles
    assert P.matrix.shape == (len(P.families), len(P.species))
    assert len(P.species) == 10
    assert (P.matrix >= 0).all()
    assert (P.presence() == (P.matrix > 0)).all()


def test_reproducible_full_run():
    ta, ga = _full_run(99)
    tb, gb = _full_run(99)
    assert ta.to_newick() == tb.to_newick()
    assert np.array_equal(ga.profiles.matrix, gb.profiles.matrix)
    assert len(ga.event_log) == len(gb.event_log)


def test_rate_model_or_shorthand_not_both():
    tree = simulate_species_tree(BirthDeath(1.0, 0.2), n_tips=6, age=2.0, seed=0)
    try:
        simulate_genomes(tree, UniformRates(loss=0.1), loss=0.2)
    except ValueError:
        pass
    else:
        raise AssertionError("expected ValueError when passing both a model and shorthand")


def test_conservation_at_speciation():
    """With only seed families (no events), every extant leaf inherits them exactly."""
    tree = simulate_species_tree(BirthDeath(1.0, 0.2), n_tips=6, age=2.0, seed=3)
    genomes = simulate_genomes(tree, initial_families=7, seed=3)  # no rates -> only seeds
    seeds = [r for r in genomes.event_log
             if r.event is EventType.ORIGINATION and r.branch == "root"]
    assert len(seeds) == 7
    for genome in genomes.leaf_genomes.values():
        assert genome.size() == 7
        assert all(genome.copy_number(str(i)) == 1 for i in range(1, 8))


def test_dl_mean_copy_number():
    """D/L-only: mean copy number at a leaf matches exp((d-l)*path_length)."""
    d, l, A = 0.1, 0.25, 2.0
    tree = simulate_species_tree(BirthDeath(1.0, 0.2), n_tips=12, age=A, seed=0)
    rates = UniformRates(duplication=d, transfer=0.0, loss=l, origination=0.0)
    gs = GenomeSimulator()

    counts = []
    for rep in range(800):
        gr = gs.simulate(tree, rates, np.random.default_rng(1000 + rep), initial_size=1)
        counts.extend(genome.copy_number("1") for genome in gr.leaf_genomes.values())

    mean = float(np.mean(counts))
    expected = math.exp((d - l) * A)  # root at time 0, path length == total_age
    assert abs(mean - expected) < 0.05, (mean, expected)
