"""SharedBirthDeath — one tree-wide diversification budget (a shared clock).

The defining property: the *total* speciation rate is fixed at ``birth`` regardless of how many
lineages stand, so diversity grows **linearly** (``E[n] = n_crown + birth·age``), not exponentially
like the per-lineage :class:`BirthDeath`.
"""
from __future__ import annotations

import numpy as np
import pytest

import zombi2 as z
from zombi2 import BirthDeath, SharedBirthDeath, simulate_species_tree


def _n_extant(tree) -> int:
    return sum(1 for leaf in tree.leaves() if leaf.is_extant)


def test_shared_grows_linearly_not_exponentially():
    # pure-birth shared clock: births are a Poisson process of constant rate `birth`, so
    # E[n(age)] = crown(2) + birth·age. Here 2 + 1·10 = 12 — vs 2·e^10 ≈ 44000 exponential.
    birth, age, reps = 1.0, 10.0, 300
    counts = [
        _n_extant(simulate_species_tree(SharedBirthDeath(birth), age=age,
                                        direction="forward", seed=s))
        for s in range(reps)
    ]
    mean = np.mean(counts)
    assert 9.0 < mean < 15.0, f"expected linear ~12, got {mean:.1f}"


def test_constant_rate_explodes_over_the_same_age():
    # the same birth rate under the standard per-lineage model is dramatically larger — the whole
    # point of the shared clock. (A handful of seeds; each should already dwarf the linear mean.)
    for s in range(3):
        n = _n_extant(simulate_species_tree(BirthDeath(1.0, 0.0), age=10.0,
                                            direction="forward", seed=s, max_lineages=10_000_000))
        assert n > 200, n


def test_shared_supports_n_tips_mode():
    tree = simulate_species_tree(SharedBirthDeath(2.0), n_tips=30, direction="forward", seed=3)
    assert _n_extant(tree) == 30


def test_shared_with_extinction_is_net_linear():
    # birth 2, death 0.5 (both shared) → net drift 1.5/time; E[n] ≈ 2 + 1.5·8 = 14 (still linear).
    counts = [
        _n_extant(simulate_species_tree(SharedBirthDeath(2.0, 0.5), age=8.0,
                                        direction="forward", seed=s))
        for s in range(300)
    ]
    assert 9.0 < np.mean(counts) < 20.0


def test_shared_in_public_api_and_namespace():
    assert "SharedBirthDeath" in z.__all__
    assert z.SharedBirthDeath is z.species.SharedBirthDeath


def test_shared_validate_rejects_nonpositive_birth():
    with pytest.raises(ValueError):
        simulate_species_tree(SharedBirthDeath(0.0), age=5.0, direction="forward", seed=1)
