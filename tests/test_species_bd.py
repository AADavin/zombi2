"""Verification for Algorithm 1 (backward reconstructed birth-death)."""

import math

import numpy as np
import pytest
from scipy import stats

from zombi2 import SpeciesTreeModel, SpeciesTreeSimulator, simulate_species_tree


def _cdf(a, lam, mu, A):
    """Analytic reconstructed-process age CDF F(a) on (0, A). Vectorized in ``a``."""
    a = np.asarray(a, dtype=float)
    r = lam - mu
    if mu == 0:  # Yule
        return (1 - np.exp(-lam * a)) / (1 - math.exp(-lam * A))
    if abs(r) < 1e-12:  # critical
        return (lam * a / (1 + lam * a)) * ((1 + lam * A) / (lam * A))
    g = lambda x: (1 - np.exp(-r * x)) / (lam - mu * np.exp(-r * x))
    return g(a) / g(A)


@pytest.mark.parametrize("lam,mu", [(1.0, 0.0), (1.0, 0.3), (1.0, 0.9), (0.7, 0.7)])
def test_sample_age_matches_cdf(lam, mu):
    """The inverse-CDF sampler reproduces the analytic age distribution (KS test)."""
    rng = np.random.default_rng(0)
    A = 5.0
    sim = SpeciesTreeSimulator()
    samples = np.array([sim._sample_age(rng.random(), lam, mu, A) for _ in range(20000)])
    assert samples.min() > 0.0 and samples.max() < A
    _, p = stats.kstest(samples, lambda a: _cdf(a, lam, mu, A))
    assert p > 1e-3, (lam, mu, p)


@pytest.mark.parametrize("N", [2, 5, 20, 100])
def test_tip_count_and_structure(N):
    rng = np.random.default_rng(1)
    tree = simulate_species_tree(SpeciesTreeModel(1.0, 0.3, N, age=4.0), rng)
    assert len(tree.leaves()) == N
    assert len(tree.internal_nodes()) == N - 1
    assert all(len(n.children) == 2 for n in tree.internal_nodes())
    assert all(abs(leaf.time - tree.total_age) < 1e-9 for leaf in tree.leaves())
    assert abs(tree.root.time) < 1e-9  # crown: root at time 0


def test_stem_conditioning_root_above_zero():
    rng = np.random.default_rng(2)
    tree = simulate_species_tree(
        SpeciesTreeModel(1.0, 0.2, 10, age=4.0, age_type="stem"), rng
    )
    assert 0.0 < tree.root.time < tree.total_age


def test_reproducible():
    model = SpeciesTreeModel(1.0, 0.3, 15, age=3.0)
    a = simulate_species_tree(model, np.random.default_rng(7)).to_newick()
    b = simulate_species_tree(model, np.random.default_rng(7)).to_newick()
    assert a == b


def test_pull_of_the_present():
    """More internal nodes near the present than near the root."""
    rng = np.random.default_rng(3)
    model = SpeciesTreeModel(1.0, 0.2, 60, age=5.0)
    ages = []
    for _ in range(30):
        tree = simulate_species_tree(model, rng)
        ages.extend(tree.total_age - n.time for n in tree.internal_nodes())
    ages = np.array(ages)
    assert (ages < model.age / 2).mean() > 0.5
