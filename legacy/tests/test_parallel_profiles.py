"""Parallel counts-only path — ``simulate_genomes(..., output="profiles", threads=N)`` / ``--threads``.

The built-in model's profile can be simulated on several cores by Poisson-thinning the gene
families across ``N`` independent copies of the engine and summing the profiles. Because families
are independent and a Poisson process splits, this is *distributionally identical* to one serial
run (a different but equivalent realization). ``threads=1`` is the exact serial engine.
"""

import numpy as np
import pytest

import zombi2 as z

pytestmark = pytest.mark.skipif(not z.rust_available(),
                                reason="zombi2_core (Rust extension) not built")

RATES = dict(duplication=0.2, transfer=0.15, loss=0.3, origination=0.6,
             initial_families=20, max_family_size=0.5)


def _tree(n=300, seed=1):
    return z.simulate_species_tree(z.BirthDeath(1.0, 0.3), n_tips=n, age=2.0, seed=seed)


def _key(pm):
    r, c, d = pm.coo
    return (pm.families, pm.species, sorted(zip(r.tolist(), c.tolist(), d.tolist())))


def test_threads1_is_exactly_serial():
    """threads=1 must be byte-identical to the default serial engine (same seed)."""
    tree = _tree()
    a = z.simulate_genomes(tree, output="profiles", seed=7, **RATES)
    b = z.simulate_genomes(tree, output="profiles", seed=7, threads=1, **RATES)
    assert _key(a) == _key(b)


def test_parallel_is_deterministic():
    """Same seed + same threads → identical profile (reproducible)."""
    tree = _tree()
    a = z.simulate_genomes(tree, output="profiles", seed=3, threads=4, **RATES)
    b = z.simulate_genomes(tree, output="profiles", seed=3, threads=4, **RATES)
    assert _key(a) == _key(b)


def test_parallel_profile_is_well_formed():
    """A parallel run yields a valid profile over the same species, with a plausible family count
    and no empty rows/cols — distributionally like the serial run, not identical to it."""
    tree = _tree()
    ser = z.simulate_genomes(tree, output="profiles", seed=5, **RATES)
    par = z.simulate_genomes(tree, output="profiles", seed=5, threads=4, **RATES)

    assert par.species == ser.species                       # same extant leaves
    assert par.nnz > 0 and par.presence_per_family().min() >= 1   # every family present somewhere
    # family count is a different realization but the same ballpark (well within 3x)
    assert 0.5 < len(par.families) / len(ser.families) < 2.0
    # copy numbers are positive integers
    _, _, data = par.coo
    assert data.min() >= 1


def test_parallel_matches_serial_in_distribution():
    """Over several seeds the parallel total copy-number matches the serial mean within a few
    standard errors (statistical, not bit-identical)."""
    tree = _tree(n=250, seed=11)
    kw = dict(duplication=0.2, transfer=0.1, loss=0.3, origination=0.5,
              initial_families=20, max_family_size=0.5)

    def totals(threads):
        return np.array([z.simulate_genomes(tree, output="profiles", seed=2000 + s,
                                            threads=threads, **kw).copies_per_species().sum()
                         for s in range(40)], float)

    ser, par = totals(1), totals(4)
    se = np.sqrt(ser.var(ddof=1) / len(ser) + par.var(ddof=1) / len(par))
    assert abs(ser.mean() - par.mean()) < 4 * se                # within 4 SE (very lenient)


def test_threads_requires_profiles_output():
    with pytest.raises(ValueError, match="only supported for output='profiles'"):
        z.simulate_genomes(_tree(), output="genomes", seed=1, threads=4, **RATES)


def test_threads_requires_builtin_model():
    """The parallel path is the Rust built-in engine only; a Python-engine model rejects threads>1."""
    with pytest.raises(ValueError, match="requires the built-in model"):
        z.simulate_genomes(_tree(), rates=z.PerGenomeRates(0.2, 0.0, 0.3, 0.5),
                           output="profiles", seed=1, threads=4)
