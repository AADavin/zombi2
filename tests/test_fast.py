"""Optional Rust fast-path engine (profiles-only). Skipped if zombi2_core isn't built."""

import numpy as np
import pytest

import zombi2 as z

pytestmark = pytest.mark.skipif(not z.rust_available(),
                                reason="zombi2_core (Rust extension) not built")

RATES = dict(duplication=0.15, transfer=0.1, loss=0.2, origination=0.5,
             initial_size=20, max_family_size=0.5)


def _tree(n=30, seed=1):
    return z.simulate_species_tree(z.BirthDeath(1.0, 0.3), n_tips=n, age=5.0, seed=seed)


def test_reproducible_same_seed():
    tree = _tree()
    a = z.simulate_profiles_fast(tree, seed=7, **RATES)
    b = z.simulate_profiles_fast(tree, seed=7, **RATES)
    assert np.array_equal(a.matrix, b.matrix)


def test_shape_and_species_are_extant_leaves():
    tree = _tree(n=30)
    pm = z.simulate_profiles_fast(tree, seed=3, **RATES)
    assert pm.matrix.shape[1] == len(tree.extant_leaves())
    assert set(pm.species) == {n.name for n in tree.extant_leaves()}
    assert (pm.matrix >= 0).all()


def test_hard_cap_respected():
    tree = _tree(n=40)
    cap = 5
    pm = z.simulate_profiles_fast(tree, duplication=0.6, transfer=0.2, loss=0.05,
                                  origination=0.4, initial_size=20, max_family_size=cap,
                                  seed=11)
    assert pm.matrix.max() <= cap


def test_accepts_uniform_rates_object():
    tree = _tree()
    obj = z.simulate_profiles_fast(tree, z.UniformRates(0.15, 0.1, 0.2, 0.5),
                                   initial_size=20, max_family_size=0.5, seed=7)
    kw = z.simulate_profiles_fast(tree, seed=7, **RATES)
    assert np.array_equal(obj.matrix, kw.matrix)


def test_rejects_unsupported_models():
    tree = _tree()
    with pytest.raises(ValueError):
        z.simulate_profiles_fast(tree, z.UniformRates(0.2, 0, 0.1, 0.3, carrying_capacity=10),
                                 seed=1)
    with pytest.raises(TypeError):
        z.simulate_profiles_fast(tree, z.GenomeWiseRates(0.2, 0.1, 0.2, 0.5), seed=1)


def test_statistically_matches_python_engine():
    # mean copy-number over the matrix should agree within Monte-Carlo error
    tree = _tree(n=60, seed=2)
    r = np.mean([z.simulate_profiles_fast(tree, duplication=0.2, loss=0.1, origination=0.4,
                                          initial_size=30, max_family_size=0.5,
                                          seed=1000 + s).matrix.mean() for s in range(15)])
    p = np.mean([z.simulate_genomes(tree, duplication=0.2, loss=0.1, origination=0.4,
                                    initial_size=30, max_family_size=0.5,
                                    seed=2000 + s).profiles.matrix.mean() for s in range(15)])
    assert abs(r - p) / p < 0.15
