"""Tests for zombi2.rate — the internal * composition and effective-rate evaluation (SPEC §5)."""

import pytest

from zombi2 import modifiers as mod
from zombi2 import scope
from zombi2.rate import Rate, as_rate


# --- composition with * ---------------------------------------------------

def test_number_times_modifier():
    r = 1.0 * mod.Time({0: 1.0, 3: 0.3})
    assert isinstance(r, Rate)
    assert r.base == 1.0 and r.scope is None
    assert r.modifiers == (mod.Time({0: 1.0, 3: 0.3}),)


def test_scope_times_modifier():
    r = scope.PerLineage(0.25) * mod.Diversity(cap=100)
    assert isinstance(r, Rate)
    assert r.scope == scope.PerLineage(0.25)
    assert r.modifiers == (mod.Diversity(cap=100),)


def test_chaining_accumulates_modifiers_in_order():
    r = 1.0 * mod.Time({0: 1.0, 3: 0.3}) * mod.Diversity(cap=100)
    assert r.base == 1.0 and r.scope is None
    assert r.modifiers == (mod.Time({0: 1.0, 3: 0.3}), mod.Diversity(cap=100))


def test_modifier_times_modifier_implies_unit_base():
    r = mod.Time({0: 2.0}) * mod.Diversity(cap=50)
    assert isinstance(r, Rate)
    assert r.base == 1.0 and r.scope is None
    assert len(r.modifiers) == 2


def test_number_on_either_side():
    a = 2.0 * mod.Diversity(cap=100)
    b = mod.Diversity(cap=100) * 2.0
    assert a.base == b.base == 2.0
    assert a.modifiers == b.modifiers


# --- effective evaluation -------------------------------------------------

def test_effective_scope_times_modifiers():
    # PerLineage(0.25).total(lineages=4)=1.0 ; Diversity(100).factor(diversity=50)=0.5 -> 0.5
    r = scope.PerLineage(0.25) * mod.Diversity(cap=100)
    assert r.effective(lineages=4, diversity=50) == pytest.approx(0.5)


def test_effective_with_default_scope_and_time():
    # base 1.0 per lineage: total(lineages=3)=3 ; Time factor at t=4 is 0.3 -> 0.9
    r = (1.0 * mod.Time({0: 1.0, 3: 0.3})).with_default_scope(scope.PerLineage)
    assert r.effective(lineages=3, time=4.0) == pytest.approx(0.9)


def test_effective_plain_number():
    r = as_rate(2.0, default_scope=scope.PerLineage)
    assert r.effective(lineages=5) == pytest.approx(10.0)


def test_effective_global_is_constant():
    r = as_rate(scope.Global(1.5), default_scope=scope.PerLineage)
    assert r.effective(lineages=100, diversity=3) == pytest.approx(1.5)


def test_effective_requires_scope():
    with pytest.raises(ValueError):
        (1.0 * mod.Time({0: 1.0})).effective(time=1.0)  # scope never resolved


# --- as_rate coercion -----------------------------------------------------

def test_as_rate_number_gets_default_scope():
    r = as_rate(0.25, default_scope=scope.PerCopy)
    assert r.scope == scope.PerCopy(0.25)


def test_as_rate_scope_kept():
    r = as_rate(scope.PerLineage(0.5), default_scope=scope.PerCopy)
    assert r.scope == scope.PerLineage(0.5)  # explicit scope wins over the default


def test_as_rate_modifier_gets_unit_base_and_default_scope():
    r = as_rate(mod.Diversity(cap=100), default_scope=scope.PerLineage)
    assert r.scope == scope.PerLineage(1.0)
    assert r.modifiers == (mod.Diversity(cap=100),)


def test_as_rate_existing_rate_resolved():
    r = as_rate(1.0 * mod.Time({0: 1.0}), default_scope=scope.PerLineage)
    assert r.scope == scope.PerLineage(1.0)


def test_as_rate_rejects_junk():
    with pytest.raises(TypeError):
        as_rate("fast", default_scope=scope.PerLineage)
    with pytest.raises(TypeError):
        as_rate(True, default_scope=scope.PerLineage)  # a bool is not a rate


def test_with_default_scope_is_noop_when_scope_set():
    r = (scope.Global(2.0) * mod.Diversity(cap=100)).with_default_scope(scope.PerLineage)
    assert r.scope == scope.Global(2.0)  # unchanged


# --- Rate is internal plumbing, not user-facing ---------------------------

def test_rate_is_frozen():
    r = 1.0 * mod.Diversity(cap=100)
    with pytest.raises(Exception):
        r.base = 2.0  # type: ignore[misc]


def test_rate_next_change_is_earliest_breakpoint():
    r = 1.0 * mod.Time({0: 1.0, 5: 0.2}) * mod.Diversity(cap=100)
    assert r.next_change(0.0) == 5
    assert r.next_change(5.0) == float("inf")
    assert (1.0 * mod.Diversity(cap=100)).next_change(0.0) == float("inf")  # no time-varying part
