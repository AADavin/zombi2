"""EpisodicBirthDeath (skyline) and GenomeWiseRates."""

import numpy as np
import pytest

from zombi2 import (
    BirthDeath,
    EpisodicBirthDeath,
    GenomeWiseRates,
    UniformRates,
    simulate_genomes,
    simulate_species_tree,
)


def test_episodic_single_epoch_matches_analytic_cdf():
    """One epoch, complete sampling: the numerical CDF must equal the analytic one."""
    lam, mu, A = 1.0, 0.3, 5.0
    m = EpisodicBirthDeath(birth=[lam], death=[mu], shifts=[])
    m._prepare(A)
    a = m._ages
    r = lam - mu
    g = (1 - np.exp(-r * a)) / (lam - mu * np.exp(-r * a))
    analytic = g / g[-1]
    assert np.max(np.abs(m._cdf - analytic)) < 1e-4


@pytest.mark.parametrize("rho", [1.0, 0.5])
def test_episodic_single_epoch_matches_birth_death_sampler(rho):
    """Feeding the same uniforms to episodic and (for rho=1) BirthDeath gives ~same ages."""
    epi = EpisodicBirthDeath(birth=[1.0], death=[0.3], shifts=[], sampling_fraction=rho)
    u = np.random.default_rng(0).random(3000)
    ages = np.array([epi.sample_internal_age(x, 5.0) for x in u])
    assert ages.min() > 0 and ages.max() < 5.0
    if rho == 1.0:
        bd = np.array([BirthDeath(1.0, 0.3).sample_internal_age(x, 5.0) for x in u])
        assert np.max(np.abs(ages - bd)) < 1e-2


def test_episodic_tree_structure_and_reproducible():
    def make():
        return EpisodicBirthDeath(birth=[1.0, 1.0], death=[0.2, 3.0], shifts=[1.0])  # mass extinction
    a = simulate_species_tree(make(), n_tips=25, age=4.0, seed=1)
    assert len(a.leaves()) == 25
    assert all(abs(l.time - a.total_age) < 1e-9 for l in a.leaves())
    b = simulate_species_tree(make(), n_tips=25, age=4.0, seed=1)
    assert a.to_newick() == b.to_newick()


def test_sampling_fraction_produces_valid_tree():
    t = simulate_species_tree(
        EpisodicBirthDeath([1.0], [0.3], [], sampling_fraction=0.3),
        n_tips=12, age=5.0, seed=1,
    )
    assert len(t.leaves()) == 12


def test_episodic_validation():
    with pytest.raises(ValueError):
        EpisodicBirthDeath([1.0, 1.0], [0.2], [1.0]).validate()      # death wrong length
    with pytest.raises(ValueError):
        EpisodicBirthDeath([1.0, 1.0], [0.2, 0.3], []).validate()    # shifts wrong length
    with pytest.raises(ValueError):
        EpisodicBirthDeath([1.0], [0.3], [], sampling_fraction=0).validate()


def test_genome_wise_rates_grow_slower_than_gene_wise():
    tree = simulate_species_tree(BirthDeath(1.0, 0.1), n_tips=12, age=3.0, seed=1)
    gw = simulate_genomes(tree, GenomeWiseRates(duplication=1.0, loss=0.2, origination=0.3),
                          initial_size=10, seed=2)
    uw = simulate_genomes(tree, UniformRates(duplication=0.3, loss=0.2, origination=0.3),
                          initial_size=10, seed=2)
    assert gw.profiles.matrix.shape[1] == 12
    # genome-wise growth is linear, so families stay smaller than the exponential gene-wise case
    assert gw.profiles.matrix.max() <= uw.profiles.matrix.max()
