"""Verification for Algorithm 1 (backward reconstructed birth-death)."""

import math

import numpy as np
import pytest

from zombi2 import BirthDeath, Yule, simulate_species_tree


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
    """The model's inverse-CDF sampler reproduces the analytic age distribution (KS)."""
    stats = pytest.importorskip("scipy.stats")  # scipy is a [dev]-only dependency
    rng = np.random.default_rng(0)
    A = 5.0
    model = BirthDeath(lam, mu)
    samples = np.array([model.sample_internal_age(rng.random(), A) for _ in range(20000)])
    assert samples.min() > 0.0 and samples.max() < A
    _, p = stats.kstest(samples, lambda a: _cdf(a, lam, mu, A))
    assert p > 1e-3, (lam, mu, p)


@pytest.mark.parametrize("N", [2, 5, 20, 100])
def test_tip_count_and_structure(N):
    tree = simulate_species_tree(BirthDeath(1.0, 0.3), n_tips=N, age=4.0, seed=1)
    assert len(tree.leaves()) == N
    assert len(tree.internal_nodes()) == N - 1
    assert all(len(n.children) == 2 for n in tree.internal_nodes())
    assert all(abs(leaf.time - tree.total_age) < 1e-9 for leaf in tree.leaves())
    assert abs(tree.root.time) < 1e-9  # crown: root at time 0


def test_yule_is_birth_death_without_extinction():
    a = simulate_species_tree(Yule(1.0), n_tips=12, age=3.0, seed=5).to_newick()
    b = simulate_species_tree(BirthDeath(1.0, 0.0), n_tips=12, age=3.0, seed=5).to_newick()
    assert a == b


def test_stem_conditioning_root_above_zero():
    tree = simulate_species_tree(BirthDeath(1.0, 0.2), n_tips=10, age=4.0, age_type="stem", seed=2)
    assert 0.0 < tree.root.time < tree.total_age


def test_reproducible():
    a = simulate_species_tree(BirthDeath(1.0, 0.3), n_tips=15, age=3.0, seed=7).to_newick()
    b = simulate_species_tree(BirthDeath(1.0, 0.3), n_tips=15, age=3.0, seed=7).to_newick()
    assert a == b


def test_pull_of_the_present():
    ages = []
    for s in range(30):
        tree = simulate_species_tree(BirthDeath(1.0, 0.2), n_tips=60, age=5.0, seed=s)
        ages.extend(tree.total_age - n.time for n in tree.internal_nodes())
    ages = np.array(ages)
    assert (ages < 5.0 / 2).mean() > 0.5
