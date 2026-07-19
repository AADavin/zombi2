"""Tests for zombi2.rates.scope — the count wrappers (SPEC §5)."""

import pytest

from zombi2.rates import scope


def test_per_lineage_scales_with_lineages():
    assert scope.PerLineage(1.0).total(lineages=4) == 4.0


def test_per_copy_scales_with_copies():
    assert scope.PerCopy(0.25).total(copies=8) == 2.0


def test_global_is_constant():
    assert scope.Global(1.5).total(lineages=100, copies=50) == 1.5


def test_per_site_and_per_chromosome():
    assert scope.PerSite(0.1).total(sites=10) == pytest.approx(1.0)
    assert scope.PerChromosome(0.02).total(chromosomes=3) == pytest.approx(0.06)


def test_extra_counts_are_ignored():
    assert scope.PerLineage(1.0).total(lineages=2, copies=99, sites=5) == 2.0


def test_missing_count_raises():
    with pytest.raises(KeyError):
        scope.PerCopy(0.25).total(lineages=4)  # no 'copies' supplied


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
    assert scope.PerLineage(0.0).total(lineages=10) == 0.0
