"""`Modifier.reads` and `Rate.carried_modifiers` — which value a modifier reads, and on what unit.

These are the two halves of the grammar written where the code can see them: a modifier is a
*reading* of a **value**, and a value's kind decides who produces the number. A measured value the
modifier computes for itself; a drawn or inherited one the engine has to draw when a unit is born
and hand back afterwards. `Rate.carried_modifiers` is how a level asks for the second sort without knowing
which modifier classes exist.
"""

from __future__ import annotations

import math

import pytest

from zombi2.params import Drift, LogNormal, PerCopy, PerLineage, Random
from zombi2.params import driver as drv
from zombi2.params import evaluate as ev
from zombi2.params import law as law
from zombi2.params import scope
from zombi2.params.evaluate import CARRIED_KINDS, DRAWN, DRIVEN, INHERITED, MEASURED
from zombi2.params.parameter import as_rate


def _rate(spec):
    return as_rate(spec, default_scope=scope.PerLineage)


@pytest.mark.parametrize("modifier, expected", [
    (drv.OnTime({0: 1.0, 3: 0.3}), (MEASURED, "run")),
    (drv.OnTotalDiversity(cap=100), (MEASURED, "run")),
    (Random('lineages', Drift(LogNormal(0.0, 0.2))), (INHERITED, "lineages")),
    (Random('lineages', LogNormal(0.0, 0.3)), (DRAWN, "lineages")),
    (Random('families', LogNormal(0.0, 0.5)), (DRAWN, "families")),
    (PerCopy(0.25).scaled_by("habitat.tsv", {"cave": 2.0}).modifiers[0], (DRIVEN, "lineages")),
])
def test_every_modifier_declares_what_it_reads(modifier, expected):
    assert modifier.reads == expected


def test_measured_and_driven_are_not_carried():
    """The engine draws nothing for these: a measured value it already has, and a driven one it
    resolves per lineage into ``drivers``."""
    rate = _rate(PerLineage(1.0).changing_at({0: 1.0})
                 .scaled_by(drv.TotalDiversity(cap=50))
                 .scaled_by("habitat.tsv", {"cave": 2.0}))
    assert rate.carried_modifiers() == ()


def test_a_carried_modifier_is_reported_with_its_unit():
    drift = Random('lineages', Drift(LogNormal(0.0, 0.2)))
    assert _rate(PerLineage(1.0).varying_among(drift)).carried_modifiers() == ((drift, "lineages"),)

    speed = Random('families', LogNormal(0.0, 0.5))
    assert _rate(PerCopy(0.25).varying_among(speed)).carried_modifiers() == ((speed, "families"),)


def test_all_of_them_are_kept_not_just_the_first():
    """The point of the query. Each engine used to hunt for *the* per-unit modifier and take the
    first match, so a second one was silently dropped and the run was not the model asked for."""
    a, b = Random('families', LogNormal(0.0, 0.5)), Random('families', LogNormal(0.0, 2.0))
    rate = _rate(PerCopy(0.25).varying_among(a).varying_among(b))
    assert rate.carried_modifiers() == ((a, "families"), (b, "families"))


def test_written_order_is_kept():
    first, second = Random('lineages', LogNormal(0.0, 0.3)), Random('families', LogNormal(0.0, 0.5))
    rate = _rate(PerCopy(0.25).varying_among(first).varying_among(second))
    assert [m for m, _ in rate.carried_modifiers()] == [first, second]


def test_unit_narrows_the_answer():
    """The unit names are **plural**. A level asking for ``'lineage'`` matches nothing and its draw
    is dropped in silence — a relaxed clock quietly becoming a strict one — which is why the two
    spellings had to differ."""
    per_lineage = Random('lineages', LogNormal(0.0, 0.3))
    per_family = Random('families', LogNormal(0.0, 0.5))
    rate = _rate(PerCopy(0.25).varying_among(per_lineage).varying_among(per_family))
    assert rate.carried_modifiers(unit="lineages") == ((per_lineage, "lineages"),)
    assert rate.carried_modifiers(unit="families") == ((per_family, "families"),)
    assert rate.carried_modifiers(unit="sites") == ()
    assert rate.carried_modifiers(unit="lineage") == ()      # the retired singular matches nothing


def test_a_bare_rate_carries_nothing():
    assert _rate(1.0).carried_modifiers() == ()


def test_a_third_party_modifier_without_reads_is_not_carried():
    """`reads` defaults to None, so a modifier written before this existed — or one of a user's
    own — is treated as computing its own factor, which is what it was already doing."""
    class OnLogTime(ev.Modifier):
        implemented_for = ("species",)

        def factor(self, *, time: float = 0.0, **_):
            return 1.0 / (1.0 + time)

    assert OnLogTime.reads is None
    assert _rate(PerLineage(1.0)._and(OnLogTime())).carried_modifiers() == ()


def test_carried_kinds_are_the_ones_needing_a_generator():
    assert set(CARRIED_KINDS) == {DRAWN, INHERITED}
    assert MEASURED not in CARRIED_KINDS and DRIVEN not in CARRIED_KINDS


def test_effective_skips_carried_modifiers_and_takes_their_product_instead():
    """A carried modifier's number does not come from the context, so `effective` must not ask it
    for one — the engine holds it and passes it in as ``carried_factor``, already multiplied out."""
    rate = _rate(PerLineage(2.0).varying_among('lineages', LogNormal(0.0, 0.5)))
    assert rate.effective(lineages=1) == pytest.approx(2.0)            # no factor supplied yet
    assert rate.effective(lineages=1, carried_factor=3.0) == pytest.approx(6.0)


def test_a_measured_modifier_still_computes_its_own_factor():
    rate = _rate(PerLineage(2.0).scaled_by(drv.TotalDiversity(cap=100)))
    assert rate.effective(lineages=1, diversity=50) == pytest.approx(1.0)


def test_the_species_engine_now_applies_every_per_lineage_modifier():
    """The regression this whole query exists for. Two per-lineage modifiers on one rate used to
    give the first one's factor, because the engine took the first match and stopped; the run then
    reported a model it had not simulated."""
    import numpy as np

    from zombi2.species import _per_lineage

    a, b = Random('lineages', LogNormal(0.0, 0.4)), Random('lineages', LogNormal(0.0, 0.9))
    rate = _rate(PerLineage(1.0).varying_among(a).varying_among(b))

    assert _per_lineage(rate) == (a, b)

    drawn = ev.values_at_birth(_per_lineage(rate), np.random.default_rng(1))
    assert len(drawn) == 2                                   # both drew, not just the first
    effective = rate.effective(lineages=1, carried_factor=math.prod(drawn))
    assert effective == pytest.approx(drawn[0] * drawn[1])   # and both reached the rate


def test_the_genome_engines_now_apply_every_per_family_modifier():
    """The same regression on the other side. `values_at_birth` is what the three genome engines use
    to give a new family its values, and it draws from every modifier the rate carries — they take
    `math.prod` of the result, because a family never splits and only the combined factor matters."""
    import numpy as np

    a, b = Random('families', LogNormal(0.0, 0.4)), Random('families', LogNormal(0.0, 0.9))
    rate = as_rate(PerCopy(0.25).varying_among(a).varying_among(b), default_scope=scope.PerCopy)
    carried = tuple(m for m, _ in rate.carried_modifiers(unit="families"))
    assert carried == (a, b)

    seen = ev.values_at_birth(carried, np.random.default_rng(7))
    rng = np.random.default_rng(7)
    assert seen == pytest.approx((a.draw(rng), b.draw(rng)))   # both draws, in written order
    assert math.prod(seen) == pytest.approx(seen[0] * seen[1])


class TestOneObjectIsOneDraw:
    """Sharing is by **object**, not by value. Writing one modifier and reading it from two rates
    says "these move together"; building two says "these are independent". Nothing else is needed to
    tell the two models apart, and no extra argument."""

    def test_one_object_on_two_rates_draws_once(self):
        import numpy as np

        speed = Random('families', LogNormal(0.0, 0.5))
        shared: dict[int, float] = {}       # one unit's cache, passed to each of its rates
        rng = np.random.default_rng(11)

        first = ev.values_at_birth((speed,), rng, shared)
        second = ev.values_at_birth((speed,), rng, shared)
        assert first == second                                  # the same number both times

        expected = Random('families', LogNormal(0.0, 0.5)).draw(np.random.default_rng(11))
        assert first == pytest.approx((expected,))              # and the generator moved once

    def test_two_objects_of_the_same_spread_draw_separately(self):
        import numpy as np

        a = Random('families', LogNormal(0.0, 0.5))
        b = Random('families', LogNormal(0.0, 0.5))
        shared: dict[int, float] = {}
        rng = np.random.default_rng(11)

        assert ev.values_at_birth((a,), rng, shared) != ev.values_at_birth((b,), rng, shared)

    def test_the_cache_is_per_unit_so_a_later_unit_draws_afresh(self):
        import numpy as np

        speed = Random('families', LogNormal(0.0, 0.5))
        rng = np.random.default_rng(11)
        one = ev.values_at_birth((speed,), rng, {})               # family one
        two = ev.values_at_birth((speed,), rng, {})               # family two
        assert one != two

    def test_without_a_cache_every_modifier_draws_for_itself(self):
        import numpy as np

        speed = Random('families', LogNormal(0.0, 0.5))
        one, two = ev.values_at_birth((speed, speed), np.random.default_rng(11))
        assert one != two          # the same object, twice, but no cache — so two draws

    def test_species_shares_one_object_between_birth_and_death(self):
        """A lineage's cache spans both rates, so one `Random` read by birth and by death is one
        number per lineage — the two rates rise and fall together."""
        import numpy as np

        from zombi2.species import _per_lineage

        s = Random('lineages', LogNormal(0.0, 0.4))
        birth = _per_lineage(_rate(PerLineage(1.0).varying_among(s)))
        death = _per_lineage(_rate(PerLineage(0.5).varying_among(s)))

        lineage: dict[int, float] = {}
        rng = np.random.default_rng(4)
        (b,) = ev.values_at_birth(birth, rng, lineage)
        (d,) = ev.values_at_birth(death, rng, lineage)
        assert b == d

    def test_species_keeps_two_objects_independent(self):
        import numpy as np

        from zombi2.species import _per_lineage

        birth = _per_lineage(_rate(PerLineage(1.0).varying_among('lineages', LogNormal(0.0, 0.4))))
        death = _per_lineage(_rate(PerLineage(0.5).varying_among('lineages', LogNormal(0.0, 0.4))))

        lineage: dict[int, float] = {}
        rng = np.random.default_rng(4)
        (b,) = ev.values_at_birth(birth, rng, lineage)
        (d,) = ev.values_at_birth(death, rng, lineage)
        assert b != d


def test_a_modifier_that_does_not_draw_says_so():
    """an inherited value starts from its parent rather than from nothing, so it has no `draw`. The base
    raises with the reason instead of handing back a plausible 1.0."""
    with pytest.raises(NotImplementedError, match="does not draw a value per unit"):
        Random('lineages', Drift(LogNormal(0.0, 0.2))).draw(object())

    with pytest.raises(NotImplementedError):
        drv.OnTime({0: 1.0}).draw(object())


class TestTheEscapeHatchCannotVouchForACarriedValue:
    """`implemented_for` lets a modifier of your own declare which engines it works with. It can
    promise that for a factor it *computes* — that is a promise only the modifier has to keep. It
    cannot promise it for a carried value, because that number is drawn by the engine when a unit is
    born and handed back, and an engine can only do that for the units it declares.

    Admitting one anyway drops it twice: `Rate.carried_modifiers` never returns it (wrong unit, so nothing
    draws it) and `Rate.effective` skips it because its kind is carried, so `factor` is never called
    either. The rate then runs undriven in silence — the exact failure this gate exists to prevent,
    and newly possible because `reads` did not exist before.
    """

    def _third_party(self, reads, engine):
        class Mine(ev.Modifier):
            implemented_for = (engine,)

            def factor(self, **_):
                return 2.0

            def draw(self, rng):
                return 2.0

        Mine.reads = reads
        return Mine()

    def test_it_vouches_for_a_computed_factor(self):
        for kind in (MEASURED, DRIVEN):
            m = self._third_party((kind, "run"), "species")
            assert ev.is_implemented(m, (), "species")

    def test_it_does_not_vouch_for_a_drawn_or_inherited_one(self):
        for kind in (DRAWN, INHERITED):
            m = self._third_party((kind, "sites"), "genomes.family")
            assert not ev.is_implemented(m, (), "genomes.family")

    def test_a_modifier_with_no_reads_is_unaffected(self):
        """The pre-existing case: a modifier written before `reads` existed computes its own factor,
        which is what it always did, so the hatch keeps working for it."""
        m = self._third_party(None, "species")
        assert ev.is_implemented(m, (), "species")

    def test_the_engine_refuses_it_rather_than_running_undriven(self):
        from zombi2 import genomes
        from zombi2.species import simulate_species_tree

        tree = simulate_species_tree(birth=1.0, death=0.2, n_extant=6, seed=1).complete_tree
        m = self._third_party((DRAWN, "sites"), "genomes.family")
        with pytest.raises(ValueError, match="does not support"):
            genomes.simulate_genomes_family(tree, loss=PerCopy(0.2)._and(m),
                                            initial_families=3, seed=1)
