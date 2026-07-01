"""Verification for Algorithm 2 (forward D/T/L/O gene families)."""

import math

import numpy as np

from zombi2 import (
    EventRates,
    GenomeSimulator,
    RateModel,
    Simulation,
    SpeciesTreeModel,
    SpeciesTreeSimulator,
)
from zombi2.events import EventType


def _name_to_node(tree):
    return {n.name: n for n in tree.nodes_preorder()}


def test_transfer_time_consistency():
    """Every transfer's donor and recipient co-exist at the transfer time."""
    sp = SpeciesTreeModel(1.0, 0.2, 15, age=4.0)
    rates = RateModel(EventRates(duplication=0.1, transfer=0.3, loss=0.15, origination=0.3))
    res = Simulation(sp, rates, seed=5, initial_size=15).run()
    name2node = _name_to_node(res.species_tree)

    n_transfers = 0
    for r in res.event_log:
        if r.event is not EventType.TRANSFER:
            continue
        n_transfers += 1
        donor, recipient = name2node[r.donor], name2node[r.recipient]
        assert r.donor != r.recipient
        for node in (donor, recipient):
            assert node.parent is not None
            assert node.parent.time < r.time <= node.time + 1e-9
    assert n_transfers > 0  # confirm the invariant was actually exercised


def test_profile_shape_and_binarization():
    sp = SpeciesTreeModel(1.0, 0.2, 10, age=3.0)
    rates = RateModel(EventRates(0.1, 0.05, 0.15, 0.3))
    res = Simulation(sp, rates, seed=1, initial_size=10).run()
    P = res.profiles
    assert P.matrix.shape == (len(P.families), len(P.species))
    assert len(P.species) == 10
    assert (P.matrix >= 0).all()
    assert (P.presence() == (P.matrix > 0)).all()


def test_reproducible_full_run():
    sp = SpeciesTreeModel(1.0, 0.2, 12, age=3.0)
    rates = RateModel(EventRates(0.15, 0.1, 0.2, 0.4))
    a = Simulation(sp, rates, seed=99, initial_size=8).run()
    b = Simulation(sp, rates, seed=99, initial_size=8).run()
    assert a.species_tree.to_newick() == b.species_tree.to_newick()
    assert np.array_equal(a.profiles.matrix, b.profiles.matrix)
    assert len(a.event_log) == len(b.event_log)


def test_conservation_at_speciation():
    """With only seed families (no events), every extant leaf inherits them exactly."""
    sp = SpeciesTreeModel(1.0, 0.2, 6, age=2.0)
    rates = RateModel(EventRates(0.0, 0.0, 0.0, 0.0))
    res = Simulation(sp, rates, seed=3, initial_size=7).run()
    seeds = [r for r in res.event_log if r.event is EventType.ORIGINATION and r.branch == "root"]
    assert len(seeds) == 7
    for genome in res.genomes.leaf_genomes.values():
        assert genome.size() == 7
        assert all(genome.copy_number(str(i)) == 1 for i in range(1, 8))


def test_dl_mean_copy_number():
    """D/L-only: mean copy number at a leaf matches exp((d-l)*path_length)."""
    d, l, A = 0.1, 0.25, 2.0
    sp = SpeciesTreeModel(1.0, 0.2, 12, age=A)
    tree = SpeciesTreeSimulator().simulate(sp, np.random.default_rng(0))  # one fixed tree
    rates = RateModel(EventRates(duplication=d, transfer=0.0, loss=l, origination=0.0))
    gs = GenomeSimulator()

    counts = []
    for rep in range(800):
        gr = gs.simulate(tree, rates, np.random.default_rng(1000 + rep), initial_size=1)
        counts.extend(genome.copy_number("1") for genome in gr.leaf_genomes.values())

    mean = float(np.mean(counts))
    expected = math.exp((d - l) * A)  # root at time 0, so path length == total_age
    assert abs(mean - expected) < 0.05, (mean, expected)
