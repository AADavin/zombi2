"""BranchRates — per-branch scaling of a base rate model (explicit / i.i.d. / autocorrelated)."""

import numpy as np

from zombi2 import (
    BirthDeath,
    BranchRates,
    LogNormal,
    UniformRates,
    simulate_genomes,
    simulate_species_tree,
)
from zombi2.events import EventType


def _base():
    return UniformRates(duplication=0.2, transfer=0.1, loss=0.2, origination=0.4)


def test_explicit_map_concentrates_events_on_boosted_branch():
    tree = simulate_species_tree(BirthDeath(1.0, 0.2), n_tips=12, age=4.0, seed=1)
    # pick an internal branch and give it a 10x rate boost
    boosted = tree.internal_nodes()[2].name
    rates = BranchRates(_base(), factors={boosted: 10.0})
    g = simulate_genomes(tree, rates, initial_size=10, seed=2)
    dtl = [EventType.DUPLICATION, EventType.TRANSFER, EventType.LOSS]
    from collections import Counter
    per_branch = Counter(r.branch for r in g.event_log if r.event in dtl)
    # the boosted branch should carry an outsized share of D/T/L events
    assert per_branch[boosted] > 0
    others = [c for b, c in per_branch.items() if b != boosted]
    assert per_branch[boosted] > (max(others) if others else 0)


def test_autocorr_zero_sigma_equals_base():
    tree = simulate_species_tree(BirthDeath(1.0, 0.2), n_tips=10, age=3.0, seed=3)
    plain = simulate_genomes(tree, _base(), initial_size=10, seed=4)
    flat = simulate_genomes(tree, BranchRates(_base(), autocorr_sigma=0.0), initial_size=10, seed=4)
    # sigma=0 -> every factor is 1 -> identical dynamics
    assert np.array_equal(plain.profiles.matrix, flat.profiles.matrix)
    assert len(plain.event_log) == len(flat.event_log)


def test_autocorrelated_runs_and_correlates_relatives():
    tree = simulate_species_tree(BirthDeath(1.0, 0.2), n_tips=20, age=4.0, seed=5)
    rates = BranchRates(_base(), autocorr_sigma=0.6)
    g = simulate_genomes(tree, rates, initial_size=15, seed=6)
    assert g.profiles.matrix.shape[1] == 20
    # a child branch's factor is its parent's times one lognormal step -> positively correlated
    parent_f, child_f = [], []
    for node in tree.nodes_preorder():
        if node.parent is not None and node.parent.parent is not None:
            parent_f.append(rates._factor[node.parent.name])
            child_f.append(rates._factor[node.name])
    assert np.corrcoef(parent_f, child_f)[0, 1] > 0.3


def test_iid_per_branch_runs():
    tree = simulate_species_tree(BirthDeath(1.0, 0.2), n_tips=12, age=3.0, seed=7)
    rates = BranchRates(_base(), per_branch=LogNormal(0.0, 0.5))
    g = simulate_genomes(tree, rates, initial_size=10, seed=8)
    assert g.profiles.matrix.shape[1] == 12
    assert len(set(rates._factor.values())) > 1  # branches really got different factors


def test_branch_rates_reproducible():
    tree = simulate_species_tree(BirthDeath(1.0, 0.2), n_tips=10, age=3.0, seed=9)
    make = lambda: BranchRates(_base(), autocorr_sigma=0.5)
    a = simulate_genomes(tree, make(), initial_size=8, seed=10)
    b = simulate_genomes(tree, make(), initial_size=8, seed=10)
    assert np.array_equal(a.profiles.matrix, b.profiles.matrix)


def test_bad_source_count_rejected():
    import pytest
    with pytest.raises(ValueError):
        BranchRates(_base())  # no source
    with pytest.raises(ValueError):
        BranchRates(_base(), autocorr_sigma=0.5, factors={})  # two sources
