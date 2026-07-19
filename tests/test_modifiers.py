"""Tests for zombi2.modifiers — the deterministic rate modifiers (SPEC §5)."""

import pytest

from zombi2 import modifiers as mod


# --- Time -----------------------------------------------------------------

def test_time_piecewise_constant():
    t = mod.Time({0: 1.0, 3: 0.3})
    assert t.factor(time=0.0) == 1.0
    assert t.factor(time=2.9) == 1.0
    assert t.factor(time=3.0) == pytest.approx(0.3)   # inclusive at the breakpoint
    assert t.factor(time=10.0) == pytest.approx(0.3)


def test_time_before_first_breakpoint_uses_earliest():
    t = mod.Time({2: 0.5, 5: 0.1})
    assert t.factor(time=0.0) == 0.5   # before the earliest key, earliest factor applies


def test_time_single_entry_is_constant():
    t = mod.Time({0: 2.0})
    assert t.factor(time=0.0) == 2.0
    assert t.factor(time=99.0) == 2.0


def test_time_ignores_extra_context():
    assert mod.Time({0: 1.0, 4: 0.5}).factor(time=5.0, diversity=10, branch="x") == pytest.approx(0.5)


def test_time_missing_context_raises():
    with pytest.raises(TypeError):
        mod.Time({0: 1.0}).factor(diversity=3)   # no 'time'


def test_time_validation():
    with pytest.raises(ValueError):
        mod.Time({})                      # empty
    with pytest.raises(ValueError):
        mod.Time({0: -1.0})               # negative factor
    with pytest.raises(ValueError):
        mod.Time({0: float("inf")})       # non-finite factor


def test_time_equality_and_repr():
    assert mod.Time({0: 1.0, 3: 0.3}) == mod.Time({3: 0.3, 0: 1.0})   # order-independent
    assert mod.Time({0: 1.0}) != mod.Time({0: 2.0})
    assert "Time(" in repr(mod.Time({0: 1.0, 3: 0.3}))
    assert hash(mod.Time({0: 1.0})) == hash(mod.Time({0: 1.0}))


# --- Diversity ------------------------------------------------------------

def test_diversity_linear_falloff():
    d = mod.Diversity(cap=100)
    assert d.factor(diversity=0) == 1.0
    assert d.factor(diversity=50) == pytest.approx(0.5)
    assert d.factor(diversity=100) == 0.0


def test_diversity_clamps_at_zero_beyond_cap():
    assert mod.Diversity(cap=10).factor(diversity=25) == 0.0   # never negative


def test_diversity_ignores_extra_context():
    assert mod.Diversity(cap=100).factor(diversity=25, time=3.0) == pytest.approx(0.75)


def test_diversity_missing_context_raises():
    with pytest.raises(TypeError):
        mod.Diversity(cap=100).factor(time=1.0)   # no 'diversity'


def test_diversity_validation():
    with pytest.raises(ValueError):
        mod.Diversity(cap=0)
    with pytest.raises(ValueError):
        mod.Diversity(cap=-5)
    with pytest.raises(ValueError):
        mod.Diversity(cap=float("nan"))
    with pytest.raises(TypeError):
        mod.Diversity(cap="big")          # type: ignore[arg-type]


def test_diversity_frozen_and_equal():
    d = mod.Diversity(cap=100)
    with pytest.raises(Exception):
        d.cap = 50                        # type: ignore[misc]
    assert mod.Diversity(cap=100) == mod.Diversity(cap=100)
    assert mod.Diversity(cap=100) != mod.Diversity(cap=50)


# --- the base / the module surface ---------------------------------------

def test_base_modifier_is_abstract():
    with pytest.raises(NotImplementedError):
        mod.Modifier().factor(time=1.0)


def test_only_deterministic_modifiers_here():
    # the stochastic ones live in the next module; assert we did not leak them yet
    for stochastic in ("Inherited", "ByBranch", "ByFamily", "Speed", "Markov", "DrivenBy"):
        assert not hasattr(mod, stochastic), f"{stochastic} should not be in this module yet"


def test_time_next_change():
    t = mod.Time({0: 1.0, 3: 0.3, 7: 0.1})
    assert t.next_change(0.0) == 3
    assert t.next_change(3.0) == 7        # strictly after the current time
    assert t.next_change(5.0) == 7
    assert t.next_change(7.0) == float("inf")  # nothing after the last breakpoint


def test_diversity_never_changes_with_time():
    assert mod.Diversity(cap=100).next_change(3.0) == float("inf")
