"""Tests for zombi2.params.scope — the "per what?" of a rate (SPEC §5)."""

import pytest

from zombi2.params import scope
from zombi2.params.rate import Rate


# --- calling a scope IS the rate ------------------------------------------

def test_calling_a_scope_builds_a_rate():
    r = scope.PerCopy(0.25)
    assert isinstance(r, Rate)
    assert r.base == 0.25


def test_the_rate_holds_the_scope_CLASS_not_an_instance():
    """A scope is a marker for the unit, with no fields — so a `Scope` instance never exists and
    `Rate.scope` is the class itself. Engines ask ``rate.scope is PerLineage``; an ``isinstance``
    check against the same name is False for every rate that exists, and reports the wrong scope
    rather than raising, which is why it is worth pinning here."""
    assert scope.PerCopy(0.25).scope is scope.PerCopy
    assert not isinstance(scope.PerCopy(0.25).scope, scope.PerCopy)


def test_a_scope_with_no_number_is_what_set_by_is_written_from():
    r = scope.PerCopy()
    assert r.base is None and r.scope is scope.PerCopy


def test_the_abstract_scope_names_the_ones_to_write():
    with pytest.raises(TypeError, match="PerLineage"):
        scope.Scope(1.0)


# --- total_of: the base, times the count this scope reads -----------------

def test_per_lineage_scales_with_lineages():
    assert scope.PerLineage.total_of(1.0, lineages=4) == 4.0


def test_per_copy_scales_with_copies():
    assert scope.PerCopy.total_of(0.25, copies=8) == 2.0


def test_global_is_constant():
    assert scope.Global.total_of(1.5, lineages=100, copies=50) == 1.5


def test_per_site_and_per_chromosome():
    assert scope.PerSite.total_of(0.1, sites=10) == pytest.approx(1.0)
    assert scope.PerChromosome.total_of(0.02, chromosomes=3) == pytest.approx(0.06)


def test_extra_counts_are_ignored():
    assert scope.PerLineage.total_of(1.0, lineages=2, copies=99, sites=5) == 2.0


def test_missing_count_raises():
    with pytest.raises(KeyError):
        scope.PerCopy.total_of(0.25, lineages=4)  # no 'copies' supplied


def test_the_units_are_the_plural_count_keywords():
    assert scope.Global.unit is None
    assert [s.unit for s in (scope.PerLineage, scope.PerCopy, scope.PerSite, scope.PerChromosome)] \
        == ["lineages", "copies", "sites", "chromosomes"]


# --- the base a scope is written with -------------------------------------

def test_negative_base_rejected():
    with pytest.raises(ValueError):
        scope.PerLineage(-1.0)


def test_nonfinite_base_rejected():
    with pytest.raises(ValueError):
        scope.Global(float("inf"))
    with pytest.raises(ValueError):
        scope.PerCopy(float("nan"))


def test_nonnumeric_base_rejected():
    with pytest.raises(TypeError):
        scope.PerLineage("fast")  # type: ignore[arg-type]
    with pytest.raises(TypeError):
        scope.Global(True)  # a bool is not a rate


def test_frozen():
    r = scope.PerLineage(1.0)
    with pytest.raises(Exception):
        r.base = 2.0  # type: ignore[misc]


def test_equality_is_by_unit_and_base():
    assert scope.PerLineage(1.0) == scope.PerLineage(1.0)
    assert scope.PerLineage(1.0) != scope.PerCopy(1.0)  # different unit → different rate
    assert scope.PerCopy(1.0) != scope.PerCopy(2.0)


def test_no_per_genome():
    # dropped: one genome per lineage, so "per genome" is PerLineage
    assert not hasattr(scope, "PerGenome")


def test_zero_base_is_allowed():
    # a zero rate (e.g. death=0 = Yule) is legal
    assert scope.PerLineage.total_of(0.0, lineages=10) == 0.0
