"""`Modifier.reads` and `Rate.carried` — which value a modifier reads, and on what unit.

These are the two halves of the grammar written where the code can see them: a modifier is a
*reading* of a **value**, and a value's kind decides who produces the number. A measured value the
modifier computes for itself; a drawn or inherited one the engine has to draw when a unit is born
and hand back afterwards. `Rate.carried` is how a level asks for the second sort without knowing
which modifier classes exist.
"""

from __future__ import annotations

import pytest

from zombi2.rates import modifiers as mod
from zombi2.rates import scope
from zombi2.rates.modifiers import CARRIED_KINDS, DRAWN, DRIVEN, INHERITED, MEASURED
from zombi2.rates.rate import as_rate


def _rate(spec):
    return as_rate(spec, default_scope=scope.PerLineage)


@pytest.mark.parametrize("modifier, expected", [
    (mod.OnTime({0: 1.0, 3: 0.3}), (MEASURED, "run")),
    (mod.OnTotalDiversity(cap=100), (MEASURED, "run")),
    (mod.FromParent(spread=0.2), (INHERITED, "lineage")),
    (mod.ByLineage(spread=0.3), (DRAWN, "lineage")),
    (mod.ByFamily(spread=0.5), (DRAWN, "family")),
    (mod.DrivenBy("habitat.tsv", {"cave": 2.0}), (DRIVEN, "lineage")),
])
def test_every_modifier_declares_what_it_reads(modifier, expected):
    assert modifier.reads == expected


def test_measured_and_driven_are_not_carried():
    """The engine draws nothing for these: a measured value it already has, and a driven one it
    resolves per lineage into ``drivers``."""
    rate = _rate(1.0 * mod.OnTime({0: 1.0}) * mod.OnTotalDiversity(cap=50)
                 * mod.DrivenBy("habitat.tsv", {"cave": 2.0}))
    assert rate.carried() == ()


def test_a_carried_modifier_is_reported_with_its_unit():
    drift = mod.FromParent(spread=0.2)
    assert _rate(1.0 * drift).carried() == ((drift, "lineage"),)

    speed = mod.ByFamily(spread=0.5)
    assert _rate(0.25 * speed).carried() == ((speed, "family"),)


def test_all_of_them_are_kept_not_just_the_first():
    """The point of the query. Each engine used to hunt for *the* per-unit modifier and take the
    first match, so a second one was silently dropped and the run was not the model asked for."""
    a, b = mod.ByFamily(spread=0.5), mod.ByFamily(spread=2.0)
    assert _rate(0.25 * a * b).carried() == ((a, "family"), (b, "family"))


def test_written_order_is_kept():
    first, second = mod.ByLineage(spread=0.3), mod.ByFamily(spread=0.5)
    assert [m for m, _ in _rate(0.25 * first * second).carried()] == [first, second]


def test_unit_narrows_the_answer():
    per_lineage, per_family = mod.ByLineage(spread=0.3), mod.ByFamily(spread=0.5)
    rate = _rate(0.25 * per_lineage * per_family)
    assert rate.carried(unit="lineage") == ((per_lineage, "lineage"),)
    assert rate.carried(unit="family") == ((per_family, "family"),)
    assert rate.carried(unit="site") == ()


def test_a_bare_rate_carries_nothing():
    assert _rate(1.0).carried() == ()


def test_a_third_party_modifier_without_reads_is_not_carried():
    """`reads` defaults to None, so a modifier written before this existed — or one of a user's
    own — is treated as computing its own factor, which is what it was already doing."""
    class OnLogTime(mod.Modifier):
        implemented_for = ("species",)

        def factor(self, *, time: float = 0.0, **_):
            return 1.0 / (1.0 + time)

    assert OnLogTime.reads is None
    assert _rate(1.0 * OnLogTime()).carried() == ()


def test_carried_kinds_are_the_ones_needing_a_generator():
    assert set(CARRIED_KINDS) == {DRAWN, INHERITED}
    assert MEASURED not in CARRIED_KINDS and DRIVEN not in CARRIED_KINDS


def test_effective_skips_carried_modifiers_and_takes_their_product_instead():
    """A carried modifier's number does not come from the context, so `effective` must not ask it
    for one — the engine holds it and passes it in as ``carried``, already multiplied out."""
    rate = _rate(2.0 * mod.ByLineage(spread=0.5))
    assert rate.effective(lineages=1) == pytest.approx(2.0)            # no factor supplied yet
    assert rate.effective(lineages=1, carried=3.0) == pytest.approx(6.0)


def test_a_measured_modifier_still_computes_its_own_factor():
    rate = _rate(2.0 * mod.OnTotalDiversity(cap=100))
    assert rate.effective(lineages=1, diversity=50) == pytest.approx(1.0)


def test_the_species_engine_now_applies_every_per_lineage_modifier():
    """The regression this whole query exists for. Two per-lineage modifiers on one rate used to
    give the first one's factor, because the engine took the first match and stopped; the run then
    reported a model it had not simulated."""
    import numpy as np

    from zombi2.species import _born, _per_lineage, _product

    a, b = mod.ByLineage(spread=0.4), mod.ByLineage(spread=0.9)
    rate = _rate(1.0 * a * b)

    assert _per_lineage(rate) == (a, b)

    drawn = _born(_per_lineage(rate), np.random.default_rng(1))
    assert len(drawn) == 2                                   # both drew, not just the first
    effective = rate.effective(lineages=1, carried=_product(drawn))
    assert effective == pytest.approx(drawn[0] * drawn[1])   # and both reached the rate
