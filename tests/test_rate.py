"""Tests for zombi2.rates.rate — how the verbs compose a rate, and how it evaluates (SPEC §5)."""

import pytest

from zombi2.rates import LogNormal, TotalDiversity, modifiers as mod
from zombi2.rates import scope
from zombi2.rates.rate import Rate, RateCompositionError, as_rate


# --- composition: a scope, with verbs chained onto it ---------------------

def test_a_scope_with_a_verb_carries_the_modifier_that_verb_builds():
    r = scope.PerLineage(1.0).changing_at({0: 1.0, 3: 0.3})
    assert isinstance(r, Rate)
    assert r.base == 1.0 and r.scope is scope.PerLineage
    assert r.modifiers == (mod.OnTime({0: 1.0, 3: 0.3}),)


def test_scaled_by_a_driver_keeps_the_scope():
    r = scope.PerLineage(0.25).scaled_by(TotalDiversity(cap=100))
    assert r.scope is scope.PerLineage and r.base == 0.25
    assert r.modifiers == (mod.OnTotalDiversity(cap=100),)


def test_chaining_accumulates_modifiers_in_order():
    r = scope.PerLineage(1.0).changing_at({0: 1.0, 3: 0.3}).scaled_by(TotalDiversity(cap=100))
    assert r.base == 1.0 and r.scope is scope.PerLineage
    assert r.modifiers == (mod.OnTime({0: 1.0, 3: 0.3}), mod.OnTotalDiversity(cap=100))


def test_every_verb_returns_a_new_rate_so_a_chain_cannot_alias():
    """A parameter sweep built by chaining shares nothing: the rate a verb was called on is
    unchanged, which is the one thing a mutable model could not promise."""
    a = scope.PerCopy(0.25)
    b = a.changing_at({0: 1.0, 5: 0.2})
    assert a.modifiers == () and len(b.modifiers) == 1
    assert a is not b


def test_star_no_longer_composes_and_names_the_verbs():
    with pytest.raises(RateCompositionError, match=r"scaled_by"):
        scope.PerLineage(0.25) * mod.OnTotalDiversity(cap=100)
    with pytest.raises(RateCompositionError, match=r"varying_among"):
        0.25 * mod.OnTotalDiversity(cap=100)
    with pytest.raises(RateCompositionError, match=r"changing_at"):
        mod.OnTime({0: 2.0}) * mod.OnTotalDiversity(cap=50)


def test_the_star_refusal_is_a_type_error_so_old_handlers_still_catch_it():
    assert issubclass(RateCompositionError, TypeError)


# --- set_by is written first, on a bare scope -----------------------------

def test_set_by_is_written_from_the_bare_scope():
    r = scope.PerCopy().set_by("habitat.tsv", {"aquatic": 1.0, "terrestrial": 0.25})
    assert r.base is None and r.scope is scope.PerCopy
    assert isinstance(r.modifiers[0], mod.SetBy)


def test_set_by_after_a_base_is_refused_because_it_would_discard_it():
    with pytest.raises(RateCompositionError, match="silently discard"):
        scope.PerCopy(0.25).set_by("habitat.tsv", {"aquatic": 1.0})


def test_set_by_after_another_verb_is_refused_for_the_same_reason():
    with pytest.raises(RateCompositionError, match="silently discard"):
        scope.PerCopy().changing_at({0: 1.0}).set_by("habitat.tsv", {"aquatic": 1.0})


def test_a_replaced_base_may_still_be_scaled():
    r = (scope.PerCopy().set_by("habitat.tsv", {"aquatic": 1.0})
         .scaled_by("size.tsv", {"big": 2.0}))
    assert len(r.modifiers) == 2


def test_two_set_by_verbs_are_refused_rather_than_letting_the_last_one_win():
    r = scope.PerCopy().set_by("a.tsv", {"x": 1.0})._and(
        mod.SetBy("b.tsv", {"x": 2.0}, verb=mod.SET_BY))
    with pytest.raises(ValueError, match="replaced only once|can only be replaced once"):
        r.check_one_base("loss")


def test_a_scope_with_no_number_and_no_set_by_says_what_is_missing():
    with pytest.raises(ValueError, match="says per what but not how fast"):
        as_rate(scope.PerCopy(), default_scope=scope.PerCopy, label="loss")


# --- weighted_by belongs to a choice, not to a rate -----------------------

def test_weighted_by_on_a_rate_names_scaled_by():
    with pytest.raises(RateCompositionError, match="scaled_by"):
        scope.PerCopy(0.25).weighted_by("competence.tsv", {"competent": 3.0})


def test_a_rate_carrying_a_weight_is_refused_when_it_is_coerced():
    from zombi2.rates import verbs
    r = scope.PerCopy(0.25)._and(verbs.weighted_by("competence.tsv", {"competent": 3.0}))
    with pytest.raises(ValueError, match="weights the candidates of a choice"):
        r.check_one_base("loss")


# --- effective evaluation -------------------------------------------------

def test_effective_scope_times_modifiers():
    # PerLineage.total_of(0.25, lineages=4)=1.0 ; OnTotalDiversity(100) at diversity 50 = 0.5 -> 0.5
    r = scope.PerLineage(0.25).scaled_by(TotalDiversity(cap=100))
    assert r.effective(lineages=4, diversity=50) == pytest.approx(0.5)


def test_effective_with_default_scope_and_time():
    # base 1.0 per lineage: total_of(1.0, lineages=3)=3 ; OnTime factor at t=4 is 0.3 -> 0.9
    r = as_rate(1.0, default_scope=scope.PerLineage).changing_at({0: 1.0, 3: 0.3})
    assert r.effective(lineages=3, time=4.0) == pytest.approx(0.9)


def test_effective_plain_number():
    r = as_rate(2.0, default_scope=scope.PerLineage)
    assert r.effective(lineages=5) == pytest.approx(10.0)


def test_effective_global_is_constant():
    r = as_rate(scope.Global(1.5), default_scope=scope.PerLineage)
    assert r.effective(lineages=100, diversity=3) == pytest.approx(1.5)


def test_effective_requires_scope():
    with pytest.raises(ValueError, match="no scope"):
        Rate(1.0).changing_at({0: 1.0}).effective(time=1.0)  # scope never resolved


def test_effective_skips_a_carried_modifier_and_takes_the_engine_s_factor():
    """A drawn value's number is not in the context: the engine drew it when the unit was born and
    hands it back, so the loop must skip the modifier rather than ask it for a factor."""
    r = scope.PerCopy(0.25).varying_among("families", LogNormal(0.0, 0.5))
    assert r.effective(copies=4, carried_factor=2.0) == pytest.approx(2.0)


def test_a_set_by_supplies_the_base_and_the_scope_still_multiplies_it():
    r = scope.PerCopy().set_by("habitat.tsv", {"aquatic": 1.0, "terrestrial": 0.25})
    assert r.effective(copies=4, drivers={"habitat.tsv": "aquatic"}) == pytest.approx(4.0)
    assert r.effective(copies=4, drivers={"habitat.tsv": "terrestrial"}) == pytest.approx(1.0)


# --- as_rate coercion -----------------------------------------------------

def test_as_rate_number_gets_default_scope():
    r = as_rate(0.25, default_scope=scope.PerCopy)
    assert r.scope is scope.PerCopy and r.base == 0.25


def test_as_rate_scope_kept():
    r = as_rate(scope.PerLineage(0.5), default_scope=scope.PerCopy)
    assert r.scope is scope.PerLineage  # explicit scope wins over the default


def test_as_rate_takes_a_number_or_a_rate_and_nothing_else():
    """There is no third case: a scope constructor returns a `Rate` and so does every verb, so a
    bare modifier is a value rather than a parameter."""
    with pytest.raises(TypeError, match="a rate is a number, or a scope"):
        as_rate(mod.OnTotalDiversity(cap=100), default_scope=scope.PerLineage)


def test_as_rate_existing_rate_resolved():
    r = as_rate(Rate(1.0).changing_at({0: 1.0}), default_scope=scope.PerLineage)
    assert r.scope is scope.PerLineage


def test_as_rate_rejects_junk():
    with pytest.raises(TypeError):
        as_rate("fast", default_scope=scope.PerLineage)
    with pytest.raises(TypeError):
        as_rate(True, default_scope=scope.PerLineage)  # a bool is not a rate


def test_with_default_scope_is_noop_when_scope_set():
    r = scope.Global(2.0).scaled_by(TotalDiversity(cap=100)).with_default_scope(scope.PerLineage)
    assert r.scope is scope.Global  # unchanged


# --- the written form -----------------------------------------------------

def test_repr_is_the_expression_that_was_written():
    """One renderer for the log, the errors and the debugger, so a run's record is the thing that
    was typed rather than a dataclass dump."""
    assert repr(scope.PerLineage(0.5)) == "PerLineage(0.5)"
    assert repr(scope.PerLineage(0.5).changing_at({0: 1.0, 3: 0.3})) \
        == "PerLineage(0.5).changing_at({0.0: 1.0, 3.0: 0.3})"
    assert repr(scope.PerCopy().set_by("h.tsv", {"aquatic": 1.0})) \
        == "PerCopy().set_by('h.tsv', Table({'aquatic': 1.0}))"


def test_a_set_by_is_written_first_whatever_order_it_was_attached_in():
    r = scope.PerCopy().set_by("h.tsv", {"a": 1.0}).changing_at({0: 1.0, 3: 0.3})
    assert repr(r).startswith("PerCopy().set_by(")


# --- Rate is internal plumbing, not user-facing ---------------------------

def test_rate_is_frozen():
    r = scope.PerLineage(1.0).scaled_by(TotalDiversity(cap=100))
    with pytest.raises(Exception):
        r.base = 2.0  # type: ignore[misc]


def test_rate_next_change_is_earliest_breakpoint():
    r = scope.PerLineage(1.0).changing_at({0: 1.0, 5: 0.2}).scaled_by(TotalDiversity(cap=100))
    assert r.next_change(0.0) == 5
    assert r.next_change(5.0) == float("inf")
    # no time-varying part
    assert scope.PerLineage(1.0).scaled_by(TotalDiversity(cap=100)).next_change(0.0) == float("inf")
