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


def test_the_genome_engines_now_apply_every_per_family_modifier():
    """The same regression on the other side. `draw_product` is what the three genome engines use to
    give a new family its factor, and it draws from every modifier the rate carries."""
    import numpy as np

    a, b = mod.ByFamily(spread=0.4), mod.ByFamily(spread=0.9)
    rate = as_rate(0.25 * a * b, default_scope=scope.PerCopy)
    carried = tuple(m for m, _ in rate.carried(unit="family"))
    assert carried == (a, b)

    seen = mod.draw_product(carried, np.random.default_rng(7))
    rng = np.random.default_rng(7)
    assert seen == pytest.approx(a.draw(rng) * b.draw(rng))   # both draws, in written order


class TestOneObjectIsOneDraw:
    """Sharing is by **object**, not by value. Writing one modifier and reading it from two rates
    says "these move together"; building two says "these are independent". Nothing else is needed to
    tell the two models apart, and no extra argument."""

    def test_one_object_on_two_rates_draws_once(self):
        import numpy as np

        speed = mod.ByFamily(spread=0.5)
        shared: dict[int, float] = {}       # one unit's cache, passed to each of its rates
        rng = np.random.default_rng(11)

        first = mod.draw_product((speed,), rng, shared)
        second = mod.draw_product((speed,), rng, shared)
        assert first == second                                  # the same number both times

        expected = mod.ByFamily(spread=0.5).draw(np.random.default_rng(11))
        assert first == pytest.approx(expected)                 # and the generator moved once

    def test_two_objects_of_the_same_spread_draw_separately(self):
        import numpy as np

        a, b = mod.ByFamily(spread=0.5), mod.ByFamily(spread=0.5)
        shared: dict[int, float] = {}
        rng = np.random.default_rng(11)

        assert mod.draw_product((a,), rng, shared) != mod.draw_product((b,), rng, shared)

    def test_the_cache_is_per_unit_so_a_later_unit_draws_afresh(self):
        import numpy as np

        speed = mod.ByFamily(spread=0.5)
        rng = np.random.default_rng(11)
        one = mod.draw_product((speed,), rng, {})               # family one
        two = mod.draw_product((speed,), rng, {})               # family two
        assert one != two

    def test_without_a_cache_every_modifier_draws_for_itself(self):
        import numpy as np

        speed = mod.ByFamily(spread=0.5)
        rng = np.random.default_rng(11)
        assert mod.draw_product((speed, speed), rng) != pytest.approx(
            mod.draw_product((speed,), np.random.default_rng(11), {}) ** 2)

    def test_species_shares_one_object_between_birth_and_death(self):
        """A lineage's cache spans both rates, so `1.0 * s` on birth and `0.5 * s` on death is one
        number per lineage — the two rates rise and fall together."""
        import numpy as np

        from zombi2.species import _born, _per_lineage

        s = mod.ByLineage(spread=0.4)
        birth = _per_lineage(_rate(1.0 * s))
        death = _per_lineage(_rate(0.5 * s))

        lineage: dict[int, float] = {}
        rng = np.random.default_rng(4)
        (b,) = _born(birth, rng, lineage)
        (d,) = _born(death, rng, lineage)
        assert b == d

    def test_species_keeps_two_objects_independent(self):
        import numpy as np

        from zombi2.species import _born, _per_lineage

        birth = _per_lineage(_rate(1.0 * mod.ByLineage(spread=0.4)))
        death = _per_lineage(_rate(0.5 * mod.ByLineage(spread=0.4)))

        lineage: dict[int, float] = {}
        rng = np.random.default_rng(4)
        (b,) = _born(birth, rng, lineage)
        (d,) = _born(death, rng, lineage)
        assert b != d


def test_a_modifier_that_does_not_draw_says_so():
    """`FromParent` starts from its parent rather than from nothing, so it has no `draw`. The base
    raises with the reason instead of handing back a plausible 1.0."""
    with pytest.raises(NotImplementedError, match="does not draw a value per unit"):
        mod.FromParent(spread=0.2).draw(object())

    with pytest.raises(NotImplementedError):
        mod.OnTime({0: 1.0}).draw(object())
