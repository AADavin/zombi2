"""The shared-clock opportunity for birth–death: ``BirthDeath(per="shared")``.

One tree-wide diversification budget — the *total* speciation rate is fixed at ``birth`` regardless
of how many lineages stand, so diversity grows **linearly** (``E[n] = n_crown + birth·age``), not
exponentially like the per-lineage default ``BirthDeath(per="lineage")``. Also covers the deprecated
``SharedBirthDeath`` preset (kept working, warns, byte-identical to the knob).
"""
from __future__ import annotations

import warnings

import numpy as np
import pytest

import zombi2 as z
from zombi2 import BirthDeath, simulate_species_tree


def _n_extant(tree) -> int:
    return sum(1 for leaf in tree.leaves() if leaf.is_extant)


def _shared(birth, death=0.0):
    return BirthDeath(birth, death, per="shared")


def test_shared_grows_linearly_not_exponentially():
    # pure-birth shared clock: births are a Poisson process of constant rate `birth`, so
    # E[n(age)] = crown(2) + birth·age. Here 2 + 1·10 = 12 — vs 2·e^10 ≈ 44000 exponential.
    birth, age, reps = 1.0, 10.0, 300
    counts = [
        _n_extant(simulate_species_tree(_shared(birth), age=age, direction="forward", seed=s))
        for s in range(reps)
    ]
    mean = np.mean(counts)
    assert 9.0 < mean < 15.0, f"expected linear ~12, got {mean:.1f}"


def test_per_lineage_default_explodes_over_the_same_age():
    # the same birth rate under the (default) per-lineage clock is dramatically larger — the whole
    # point of the shared clock. (A handful of seeds; each should already dwarf the linear mean.)
    for s in range(3):
        n = _n_extant(simulate_species_tree(BirthDeath(1.0, 0.0), age=10.0,
                                            direction="forward", seed=s, max_lineages=10_000_000))
        assert n > 200, n
    assert BirthDeath(1.0).per == "lineage"   # per="lineage" is the default


def test_shared_supports_n_tips_mode():
    tree = simulate_species_tree(_shared(2.0), n_tips=30, direction="forward", seed=3)
    assert _n_extant(tree) == 30


def test_shared_with_extinction_is_net_linear():
    # birth 2, death 0.5 (both shared) → net drift 1.5/time; E[n] ≈ 2 + 1.5·8 = 14 (still linear).
    counts = [
        _n_extant(simulate_species_tree(_shared(2.0, 0.5), age=8.0, direction="forward", seed=s))
        for s in range(300)
    ]
    assert 9.0 < np.mean(counts) < 20.0


def test_shared_is_forward_only():
    # no closed-form reconstructed CDF → backward direction is rejected
    with pytest.raises(ValueError, match="forward"):
        simulate_species_tree(_shared(1.0), n_tips=10, age=5.0, direction="backward", seed=1)


def test_shared_validate_rejects_nonpositive_birth_and_bad_per():
    with pytest.raises(ValueError):
        simulate_species_tree(_shared(0.0), age=5.0, direction="forward", seed=1)
    with pytest.raises(ValueError, match="per"):
        simulate_species_tree(BirthDeath(1.0, per="bogus"), age=5.0, direction="forward", seed=1)


def test_shared_forbids_fossilization_removal():
    with pytest.raises(ValueError, match="per='shared'"):
        simulate_species_tree(BirthDeath(1.0, per="shared", fossilization=0.1),
                              age=5.0, direction="forward", seed=1)


# --- backwards compatibility: the deprecated SharedBirthDeath preset ---------------------------

def test_sharedbirthdeath_preset_warns_and_is_byte_identical():
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        old = z.SharedBirthDeath(1.5, 0.3)
    assert any(issubclass(w.category, DeprecationWarning) for w in caught)
    # same process → same tree for the same seed as the knob spelling
    a = simulate_species_tree(old, age=7.0, direction="forward", seed=11).to_newick()
    b = simulate_species_tree(BirthDeath(1.5, 0.3, per="shared"),
                              age=7.0, direction="forward", seed=11).to_newick()
    assert a == b


def test_sharedbirthdeath_left_public_all_but_still_importable():
    assert "SharedBirthDeath" not in z.__all__          # deprecated → out of the catalog
    assert z.SharedBirthDeath is z.species.SharedBirthDeath   # still resolvable
    assert issubclass(z.SharedBirthDeath, BirthDeath)         # a preset over the knob
