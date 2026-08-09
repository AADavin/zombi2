"""Tests for zombi2.params.evaluate — what the verbs build, and what each one reads (SPEC §5)."""

import pytest

from zombi2.params import Drift, Gamma, LogNormal, PerCopy, PerLineage, Time
from zombi2.params import driver as drv
from zombi2.params import evaluate as ev
from zombi2.params import law as law
# --- OnTime -----------------------------------------------------------------

def test_time_piecewise_constant():
    t = drv.OnTime({0: 1.0, 3: 0.3})
    assert t.factor(time=0.0) == 1.0
    assert t.factor(time=2.9) == 1.0
    assert t.factor(time=3.0) == pytest.approx(0.3)   # inclusive at the breakpoint
    assert t.factor(time=10.0) == pytest.approx(0.3)


def test_time_before_first_breakpoint_uses_earliest():
    t = drv.OnTime({2: 0.5, 5: 0.1})
    assert t.factor(time=0.0) == 0.5   # before the earliest key, earliest factor applies


def test_time_single_entry_is_constant():
    t = drv.OnTime({0: 2.0})
    assert t.factor(time=0.0) == 2.0
    assert t.factor(time=99.0) == 2.0


def test_time_ignores_extra_context():
    assert drv.OnTime({0: 1.0, 4: 0.5}).factor(time=5.0, diversity=10, branch="x") == pytest.approx(0.5)


def test_time_missing_context_raises():
    with pytest.raises(TypeError):
        drv.OnTime({0: 1.0}).factor(diversity=3)   # no 'time'


def test_time_validation():
    with pytest.raises(ValueError):
        drv.OnTime({})                      # empty
    with pytest.raises(ValueError):
        drv.OnTime({0: -1.0})               # negative factor
    with pytest.raises(ValueError):
        drv.OnTime({0: float("inf")})       # non-finite factor


def test_time_equality_and_repr():
    assert drv.OnTime({0: 1.0, 3: 0.3}) == drv.OnTime({3: 0.3, 0: 1.0})   # order-independent
    assert drv.OnTime({0: 1.0}) != drv.OnTime({0: 2.0})
    assert repr(drv.OnTime({0: 1.0, 3: 0.3})) == "changing_at({0.0: 1.0, 3.0: 0.3})"
    assert hash(drv.OnTime({0: 1.0})) == hash(drv.OnTime({0: 1.0}))


def test_the_two_readings_of_the_clock_are_one_object_that_records_which_was_written():
    """``changing_at`` holds factors and ``set_by(Time(), ...)`` holds the rates themselves, and a
    schedule on a base of 1.0 *is* the rate — so they run the same model and compare equal. Only the
    written form differs, which is why the verb is recorded."""
    factors = PerLineage(0.5).changing_at({0: 1.0, 3: 0.3}).modifiers[0]
    rates = PerLineage().set_by(Time(), {0: 0.5, 3: 0.15}).modifiers[0]
    assert factors.verb == ev.CHANGING_AT and rates.verb == ev.SET_BY
    assert drv.OnTime({0: 0.5}) == drv.OnTime({0: 0.5}, verb=ev.SET_BY)   # same model
    assert repr(rates) == "set_by(Time(), {0.0: 0.5, 3.0: 0.15})"


# --- OnTotalDiversity ------------------------------------------------------------

def test_diversity_linear_falloff():
    d = drv.OnTotalDiversity(cap=100)
    assert d.factor(diversity=0) == 1.0
    assert d.factor(diversity=50) == pytest.approx(0.5)
    assert d.factor(diversity=100) == 0.0


def test_diversity_clamps_at_zero_beyond_cap():
    assert drv.OnTotalDiversity(cap=10).factor(diversity=25) == 0.0   # never negative


def test_diversity_ignores_extra_context():
    assert drv.OnTotalDiversity(cap=100).factor(diversity=25, time=3.0) == pytest.approx(0.75)


def test_diversity_missing_context_raises():
    with pytest.raises(TypeError):
        drv.OnTotalDiversity(cap=100).factor(time=1.0)   # no 'diversity'


def test_diversity_validation():
    with pytest.raises(ValueError):
        drv.OnTotalDiversity(cap=0)
    with pytest.raises(ValueError):
        drv.OnTotalDiversity(cap=-5)
    with pytest.raises(ValueError):
        drv.OnTotalDiversity(cap=float("nan"))
    with pytest.raises(TypeError):
        drv.OnTotalDiversity(cap="big")          # type: ignore[arg-type]


def test_diversity_frozen_and_equal():
    d = drv.OnTotalDiversity(cap=100)
    with pytest.raises(Exception):
        d.cap = 50                        # type: ignore[misc]
    assert drv.OnTotalDiversity(cap=100) == drv.OnTotalDiversity(cap=100)
    assert drv.OnTotalDiversity(cap=100) != drv.OnTotalDiversity(cap=50)


def test_the_driver_carries_its_cap_and_refuses_to_be_written_without_one():
    with pytest.raises(ValueError, match="needs its carrying capacity"):
        drv.TotalDiversity()


# --- the base / the module surface ---------------------------------------

def test_base_modifier_is_abstract():
    with pytest.raises(NotImplementedError):
        ev.Modifier().factor(time=1.0)


def test_a_modifier_of_your_own_reprs_rather_than_recursing():
    """`written_call` builds its placeholder from the class name and never from ``repr(self)``:
    `__repr__` calls `written_call`, so a subclass overriding neither used to send the two into each
    other, and every log line and error message that named the rate died of a RecursionError while
    the run itself carried on fine."""
    class OnLogTime(ev.Modifier):
        implemented_for = ("species",)

        def factor(self, *, time: float = 0.0, **_):
            return 1.0 / (1.0 + time)

    r = PerLineage(0.5)._and(OnLogTime())
    assert repr(r) == "PerLineage(0.5).<OnLogTime>"
    assert r.effective(time=1.0, lineages=3) == pytest.approx(0.75)


def test_stochastic_status_built_vs_deferred():
    # each lives with the thing that writes it: a law with the law that chooses it, a connection
    # with the verb that builds it
    from zombi2.params import connection as conn
    for built in ("Drawn", "Inherited"):
        assert hasattr(law, built), f"{built} should be built"
    for built in ("Driven", "SetBy"):
        assert hasattr(conn, built), f"{built} should be built"
    for later in ("Markov",):
        assert not any(hasattr(m, later) for m in (law, conn, drv)), f"{later} is not built yet"


# --- Random with a Drift law (autocorrelated): the mean-corrected drift -----

def test_inherited_initial_is_one():
    assert law.Inherited('lineages', LogNormal(0.0, 0.3)).initial() == 1.0


def test_inherited_descend_is_mean_corrected():
    import numpy as np
    rng = np.random.default_rng(0)
    inh = law.Inherited('lineages', LogNormal(0.0, 0.5))
    vals = [inh.descend(1.0, rng) for _ in range(50000)]
    # E[factor] = 1 exactly (the -σ²/2 correction); the buggy version gives E ≈ e^{σ²/2} = 1.13
    assert abs(sum(vals) / len(vals) - 1.0) < 0.02


def test_inherited_no_inflation_over_a_chain():
    import numpy as np
    rng = np.random.default_rng(1)
    inh = law.Inherited('lineages', LogNormal(0.0, 0.4))
    ends = []
    for _ in range(20000):
        v = 1.0
        for _ in range(10):
            v = inh.descend(v, rng)
        ends.append(v)
    # 10 corrected steps still average ~1; the buggy version drifts to e^{10·σ²/2} ≈ 2.2
    assert abs(sum(ends) / len(ends) - 1.0) < 0.2


def test_inherited_deterministic():
    import numpy as np
    a = law.Inherited('lineages', LogNormal(0.0, 0.3)).descend(1.0, np.random.default_rng(7))
    b = law.Inherited('lineages', LogNormal(0.0, 0.3)).descend(1.0, np.random.default_rng(7))
    assert a == b


def test_binned_drift_stays_on_its_ladder_and_averages_one():
    """The discrete-bin clock: the same inherit-and-perturb model, in steps. A daughter takes its
    parent's rung or one either side, the ends reflect, and the ladder is scaled so the mean is 1 —
    which holds because a reflecting nearest-neighbour walk is uniform at stationarity."""
    import numpy as np
    inh = law.Inherited('lineages', LogNormal(0.0, 0.35), bins=6)
    rungs = {round(r, 9) for r in inh._ladder()}
    rng = np.random.default_rng(0)
    v, seen = inh.initial(), []
    for _ in range(60000):
        v = inh.descend(v, rng)
        seen.append(v)
    assert all(round(x, 9) in rungs for x in seen)      # never lands between two rungs
    assert abs(sum(seen) / len(seen) - 1.0) < 0.02      # and does not inflate the rate
    assert round(min(seen), 9) == min(rungs) and round(max(seen), 9) == max(rungs)  # ends reached


def test_binned_drift_leaves_the_continuous_form_alone():
    # `bins` defaults to None, so a run written before it existed draws exactly as it did
    import numpy as np
    assert law.Inherited('lineages', LogNormal(0.0, 0.3)).bins is None
    a = law.Inherited('lineages', LogNormal(0.0, 0.3)).descend(1.0, np.random.default_rng(7))
    b = law.Inherited('lineages', LogNormal(0.0, 0.3), bins=None).descend(1.0, np.random.default_rng(7))
    assert a == b


@pytest.mark.parametrize("bad", [1, 0, -3])
def test_binned_drift_needs_at_least_two_bins(bad):
    with pytest.raises(ValueError, match="at least 2"):
        law.Inherited('lineages', LogNormal(0.0, 0.3), bins=bad)


def test_a_carried_modifier_has_no_factor_to_give():
    """Its number never came from the context: the engine draws it when a unit is born, keeps it, and
    hands it back through `Rate.effective`'s ``carried_factor``. Asking for a factor without that is
    meaningless, so it raises with the reason rather than returning a plausible 1.0."""
    for m in (law.Inherited('lineages', LogNormal(0.0, 0.3)),
              law.Drawn('lineages', LogNormal(0.0, 0.3)),
              law.Drawn('families', LogNormal(0.0, 0.3))):
        with pytest.raises(NotImplementedError, match="carried"):
            m.factor()


# --- Random with a bare distribution (uncorrelated): i.i.d. corrected draws ---

def test_a_zero_spread_draw_is_a_strict_clock():
    import numpy as np
    rng = np.random.default_rng(0)
    byl = law.Drawn('lineages', LogNormal(0.0, 0.0))
    assert all(byl.draw(rng) == 1.0 for _ in range(100))


def test_a_draw_is_a_mean_corrected_lognormal():
    import numpy as np
    rng = np.random.default_rng(0)
    byl = law.Drawn('lineages', LogNormal(0.0, 0.5))
    vals = [byl.draw(rng) for _ in range(100000)]
    # E[factor] = 1 (the -σ²/2 correction); the buggy uncorrected draw gives E ≈ e^{σ²/2} = 1.13
    assert abs(sum(vals) / len(vals) - 1.0) < 0.02


def test_a_gamma_draw_is_normalised_to_mean_one():
    """Any distribution works, and a drawn value divides by its mean — so what a distribution
    contributes is its *shape*, and the base keeps meaning the average rate."""
    import numpy as np
    rng = np.random.default_rng(1)
    # shape k, scale θ has CV = 1/√k, so k = 4 is a coefficient of variation of 0.5
    byl = law.Drawn('lineages', Gamma(shape=4.0, scale=0.25))
    vals = [byl.draw(rng) for _ in range(100000)]
    assert abs(sum(vals) / len(vals) - 1.0) < 0.02
    var = sum((v - 1.0) ** 2 for v in vals) / len(vals)
    assert abs(var - 0.5 ** 2) < 0.02                        # CV = 0.5, whatever the scale was


def test_draws_are_independent_no_memory():
    import numpy as np
    rng = np.random.default_rng(2)
    byl = law.Drawn('lineages', LogNormal(0.0, 0.6))
    a = [byl.draw(rng) for _ in range(2000)]
    # i.i.d.: successive draws are uncorrelated (unlike a Drift, whose draws depend on the parent)
    lag1 = sum((a[i] - 1) * (a[i + 1] - 1) for i in range(len(a) - 1)) / (len(a) - 1)
    assert abs(lag1) < 0.05


def test_a_draw_is_deterministic():
    import numpy as np
    a = law.Drawn('lineages', LogNormal(0.0, 0.3)).draw(np.random.default_rng(7))
    b = law.Drawn('lineages', LogNormal(0.0, 0.3)).draw(np.random.default_rng(7))
    assert a == b


def test_drawn_validates_its_arguments():
    with pytest.raises(ValueError, match="LogNormal sigma must be >= 0"):
        law.Drawn('lineages', LogNormal(0.0, -0.1))
    with pytest.raises(ValueError, match="needs the law its value follows"):
        law.Drawn('lineages')                     # no law given
    with pytest.raises(ValueError, match="unknown unit"):
        law.Drawn('genomes', LogNormal(0.0, 0.3))


def test_inherited_validates_its_arguments():
    with pytest.raises(ValueError, match="LogNormal sigma must be >= 0"):
        law.Inherited('lineages', LogNormal(0.0, -0.1))
    with pytest.raises(ValueError, match="needs the distribution of its per-split step"):
        law.Inherited('lineages')
    with pytest.raises(ValueError, match="takes a LogNormal step"):
        law.Inherited('lineages', Gamma(shape=4.0, scale=0.25), bins=8)


# --- Random: the one value a user can name --------------------------------

def test_random_builds_a_draw_or_a_drift_from_its_law():
    assert isinstance(drv.Random('families', LogNormal(0.0, 0.5)), law.Drawn)
    assert isinstance(drv.Random('lineages', Drift(LogNormal(0.0, 0.2))), law.Inherited)
    assert drv.Random('lineages', Drift(LogNormal(0.0, 0.2), bins=8)).bins == 8


def test_reads_reports_the_kind_and_the_PLURAL_unit():
    """The unit an engine dispatches on. It is plural because a value varies *among* families
    rather than being counted *per* family — 'per' is the scope word (SPEC §5) — and the two
    strings differ, so a level still declaring the singular matches nothing and says so loudly."""
    assert drv.Random('families', LogNormal(0.0, 0.5)).reads == (ev.DRAWN, 'families')
    assert drv.Random('lineages', Drift(LogNormal(0.0, 0.2))).reads == (ev.INHERITED, 'lineages')
    assert ev.UNITS == ("run", "lineages", "chromosomes", "families", "copies", "sites")


def test_a_random_written_by_name_is_the_written_form_of_a_carried_value():
    assert repr(drv.Random('families', LogNormal(0.0, 0.5))) \
        == "Random('families', LogNormal(0.0, 0.5))"
    assert repr(drv.Random('lineages', Drift(LogNormal(0.0, 0.2)))) \
        == "Random('lineages', Drift(LogNormal(0.0, 0.2)))"
    # on a rate it is the verb call, because that is how it is typed there
    assert repr(PerCopy(0.25).varying_among('families', LogNormal(0.0, 0.5))) \
        == "PerCopy(0.25).varying_among('families', LogNormal(0.0, 0.5))"


def test_random_needs_its_unit():
    with pytest.raises(TypeError, match="needs the plural unit"):
        drv.Random()


def test_the_retired_spread_names_its_replacement():
    """`spread` meant the drawn value in one class and the per-split step in the other, so it is
    gone rather than aliased — and the error has to say what to write instead, including the
    distribution it stood for, or the reader has to guess one."""
    with pytest.raises(TypeError, match=r"varying_among\('families', LogNormal\(0.0, spread\)\)"):
        drv.Random('families', spread=0.5)          # type: ignore[call-arg]
    with pytest.raises(TypeError, match=r"Drift\(LogNormal\(0.0, spread\)\)"):
        PerCopy(0.25).varying_among('lineages', spread=0.2)   # type: ignore[call-arg]


def test_the_retired_per_keyword_names_the_plural_first_argument():
    with pytest.raises(TypeError, match="units are plural"):
        drv.Random(per='family', dist=LogNormal(0.0, 0.5))    # type: ignore[call-arg]
    with pytest.raises(TypeError, match="units are plural"):
        PerCopy(0.25).varying_among(per='families')           # type: ignore[call-arg]


def test_an_unknown_keyword_is_still_an_unknown_keyword():
    """The catch-all that answers ``per=`` must not swallow a typo: a binned drift written
    ``bns=8`` would otherwise run unbinned and say nothing."""
    with pytest.raises(TypeError, match="unexpected keyword argument 'bns'"):
        drv.Random('lineages', LogNormal(0.0, 0.3), bns=8)    # type: ignore[call-arg]


def test_time_next_change():
    t = drv.OnTime({0: 1.0, 3: 0.3, 7: 0.1})
    assert t.next_change(0.0) == 3
    assert t.next_change(3.0) == 7        # strictly after the current time
    assert t.next_change(5.0) == 7
    assert t.next_change(7.0) == float("inf")  # nothing after the last breakpoint


def test_diversity_never_changes_with_time():
    assert drv.OnTotalDiversity(cap=100).next_change(3.0) == float("inf")
