"""Guard-rail tests for `simulate_species_tree` on degenerate / extreme inputs.

These assert that clearly-invalid inputs raise a clean `ValueError` (rather than silently
producing NaN/inf branch lengths, hanging, or dumping a raw traceback), while valid edge
cases still simulate normally.
"""

import math

import numpy as np
import pytest

from zombi2.species.model import (
    BirthDeath, CladeShiftBirthDeath, ClaDS, DiversityDependent, EpisodicBirthDeath, Yule,
)
from zombi2.species.sim import simulate_species_tree


def _finite_tree(tree):
    """A tree with all-finite, non-negative branch lengths."""
    for node in tree.nodes_preorder():
        if node is tree.root:
            continue
        bl = node.time - node.parent.time
        assert math.isfinite(bl), f"non-finite branch length {bl}"
        assert bl >= -1e-9, f"negative branch length {bl}"
    return tree


# --- non-finite rates: previously produced NaN/inf branch lengths SILENTLY -----------

@pytest.mark.parametrize("bad", [float("nan"), float("inf"), float("-inf")])
def test_backward_nonfinite_birth_rejected(bad):
    with pytest.raises(ValueError, match="finite"):
        simulate_species_tree(BirthDeath(bad, 0.3), n_tips=10, age=1.0, seed=1)


@pytest.mark.parametrize("bad", [float("nan"), float("inf")])
def test_backward_nonfinite_death_rejected(bad):
    with pytest.raises(ValueError, match="finite"):
        simulate_species_tree(BirthDeath(1.0, bad), n_tips=10, age=1.0, seed=1)


@pytest.mark.parametrize("bad", [float("nan"), float("inf")])
def test_forward_nonfinite_birth_rejected(bad):
    with pytest.raises(ValueError, match="finite"):
        simulate_species_tree(BirthDeath(bad, 0.3), age=1.0, direction="forward",
                              seed=1, max_attempts=5)


# --- non-finite age: previously NaN age -> NaN/hang; inf age -> hang ------------------

@pytest.mark.parametrize("bad", [float("nan"), float("inf")])
def test_backward_nonfinite_age_rejected(bad):
    with pytest.raises(ValueError, match="finite"):
        simulate_species_tree(BirthDeath(1.0, 0.3), n_tips=10, age=bad, seed=1)


@pytest.mark.parametrize("bad", [float("nan"), float("inf")])
def test_forward_nonfinite_age_rejected(bad):
    # inf previously hung forever; nan previously mis-triggered the max_lineages bound
    with pytest.raises(ValueError, match="finite"):
        simulate_species_tree(BirthDeath(1.0, 0.3), age=bad, direction="forward",
                              seed=1, max_attempts=5)


# --- n_tips: type / range / lineage cap ----------------------------------------------

def test_backward_huge_n_tips_rejected():
    # 1e7 tips previously assembled slowly (a >100 s effective hang); now bounded by max_lineages
    with pytest.raises(ValueError, match="max_lineages"):
        simulate_species_tree(BirthDeath(1.0, 0.3), n_tips=10_000_000, age=1.0, seed=1)


def test_backward_n_tips_over_custom_max_lineages_rejected():
    with pytest.raises(ValueError, match="max_lineages"):
        simulate_species_tree(BirthDeath(1.0, 0.3), n_tips=5000, age=1.0,
                              max_lineages=1000, seed=1)


def test_forward_huge_n_tips_rejected():
    with pytest.raises(ValueError, match="max_lineages"):
        simulate_species_tree(BirthDeath(2.0, 0.3), n_tips=10_000_000, direction="forward",
                              max_lineages=50_000, seed=1)


def test_backward_fractional_n_tips_rejected():
    with pytest.raises(ValueError, match="whole number"):
        simulate_species_tree(BirthDeath(1.0, 0.3), n_tips=2.5, age=1.0, seed=1)


@pytest.mark.parametrize("bad", [0, 1, -5])
def test_backward_small_n_tips_rejected(bad):
    with pytest.raises(ValueError, match=">= 2"):
        simulate_species_tree(BirthDeath(1.0, 0.3), n_tips=bad, age=1.0, seed=1)


def test_backward_integer_valued_float_n_tips_accepted():
    tree = simulate_species_tree(BirthDeath(1.0, 0.3), n_tips=8.0, age=1.0, seed=1)
    assert sum(1 for n in tree.leaves() if n.is_extant) == 8


def test_backward_numpy_integer_n_tips_accepted():
    tree = simulate_species_tree(BirthDeath(1.0, 0.3), n_tips=np.int64(6), age=1.0, seed=1)
    assert sum(1 for n in tree.leaves() if n.is_extant) == 6


# --- non-finite parameters on the other models ---------------------------------------

def test_clads_nonfinite_alpha_rejected():
    with pytest.raises(ValueError, match="finite"):
        simulate_species_tree(ClaDS(1.0, alpha=float("inf")), age=1.0,
                              direction="forward", seed=1)


def test_diversity_dependent_nonfinite_K_rejected():
    with pytest.raises(ValueError, match="finite"):
        simulate_species_tree(DiversityDependent(1.0, carrying_capacity=float("nan")),
                              age=1.0, direction="forward", seed=1)


def test_episodic_nonfinite_rate_rejected():
    with pytest.raises(ValueError, match="finite"):
        simulate_species_tree(EpisodicBirthDeath([1.0, float("nan")], [0.3, 0.4], [0.5]),
                              n_tips=10, age=1.0, seed=1)


def test_clade_shift_nonfinite_age_rejected():
    with pytest.raises(ValueError, match="finite"):
        simulate_species_tree(CladeShiftBirthDeath(1.0, 0.3, clade_shifts=[(float("inf"), 2.0, 0.1)]),
                              age=2.0, direction="forward", seed=1)


def test_mass_extinction_nonfinite_fraction_rejected():
    with pytest.raises(ValueError, match="finite"):
        simulate_species_tree(BirthDeath(1.0, 0.3, mass_extinctions=[(0.5, float("nan"))]),
                              age=1.0, direction="forward", seed=1)


# --- valid edge cases must STILL work (behavior must not change) ----------------------

def test_valid_death_greater_than_birth_forward():
    # death > birth is a valid (usually-extinct) forward process
    tree = simulate_species_tree(BirthDeath(0.5, 2.0), age=1.0, direction="forward",
                                 seed=1, max_attempts=500)
    _finite_tree(tree)
    assert sum(1 for n in tree.leaves() if n.is_extant) >= 2


def test_valid_critical_birth_equals_death_backward():
    tree = simulate_species_tree(BirthDeath(1.0, 1.0), n_tips=10, age=1.0, seed=1)
    _finite_tree(tree)


def test_valid_yule_backward():
    tree = simulate_species_tree(Yule(1.0), n_tips=10, age=1.0, seed=1)
    _finite_tree(tree)
    assert sum(1 for n in tree.leaves() if n.is_extant) == 10


def test_valid_huge_but_finite_birth_backward():
    tree = simulate_species_tree(BirthDeath(1e6, 0.3), n_tips=10, age=1.0, seed=1)
    _finite_tree(tree)


def test_valid_n_tips_2():
    tree = simulate_species_tree(BirthDeath(1.0, 0.3), n_tips=2, age=1.0, seed=1)
    assert sum(1 for n in tree.leaves() if n.is_extant) == 2


def test_valid_high_turnover_clads():
    # turnover close to 1 is valid (high extinction/speciation ratio)
    tree = simulate_species_tree(ClaDS(1.0, turnover=0.95), age=0.5, direction="forward",
                                 seed=1, max_attempts=2000)
    _finite_tree(tree)


# --------------------------------------------------------------------------- #
# SpeciesCaps capability contract — an unregistered model must fail LOUDLY at
# dispatch (the silent-misroute failure mode the isinstance ladders allowed),
# and every shipped model must route to the right growth engine.
# --------------------------------------------------------------------------- #
def test_unregistered_model_raises_loudly_not_silent_misroute():
    class Bogus:                       # a duck-typed "model" with no SpeciesCaps
        birth = 1.0

        def validate(self):
            pass

    with pytest.raises(TypeError, match="no SpeciesCaps"):
        simulate_species_tree(Bogus(), age=1.0, direction="forward")


def test_species_caps_growth_and_backward():
    from zombi2.species._caps import GrowthEngine, species_caps
    thinning_backward = [BirthDeath(1.0, 0.3), Yule(1.0), EpisodicBirthDeath([1.0], [0.3], [])]
    for m in thinning_backward:
        caps = species_caps(m)         # Yule inherits BirthDeath's caps via the MRO
        assert caps.growth is GrowthEngine.THINNING
        assert caps.supports_backward
    for m in [ClaDS(1.0), DiversityDependent(1.0, carrying_capacity=50), CladeShiftBirthDeath(1.0, clade_shifts=[(0.5, 2.0, 0.1)])]:
        caps = species_caps(m)
        assert caps.growth is GrowthEngine.GILLESPIE
        assert not caps.supports_backward


def test_episodic_allows_incomplete_sampling_backward_but_constant_bd_does_not():
    # capability split the two arity-heuristic classes need distinct caps for
    simulate_species_tree(EpisodicBirthDeath([1.0], [0.3], [], sampling_fraction=0.5),
                          n_tips=8, age=2.0, seed=1)          # allowed
    with pytest.raises(ValueError, match="complete sampling"):
        simulate_species_tree(BirthDeath(1.0, 0.3, sampling_fraction=0.5), n_tips=8, age=2.0)
