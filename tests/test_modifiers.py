"""Tests for zombi2.rates.modifiers — the deterministic rate modifiers (SPEC §5)."""

import pytest

from zombi2.rates import modifiers as mod


# --- OnTime -----------------------------------------------------------------

def test_time_piecewise_constant():
    t = mod.OnTime({0: 1.0, 3: 0.3})
    assert t.factor(time=0.0) == 1.0
    assert t.factor(time=2.9) == 1.0
    assert t.factor(time=3.0) == pytest.approx(0.3)   # inclusive at the breakpoint
    assert t.factor(time=10.0) == pytest.approx(0.3)


def test_time_before_first_breakpoint_uses_earliest():
    t = mod.OnTime({2: 0.5, 5: 0.1})
    assert t.factor(time=0.0) == 0.5   # before the earliest key, earliest factor applies


def test_time_single_entry_is_constant():
    t = mod.OnTime({0: 2.0})
    assert t.factor(time=0.0) == 2.0
    assert t.factor(time=99.0) == 2.0


def test_time_ignores_extra_context():
    assert mod.OnTime({0: 1.0, 4: 0.5}).factor(time=5.0, diversity=10, branch="x") == pytest.approx(0.5)


def test_time_missing_context_raises():
    with pytest.raises(TypeError):
        mod.OnTime({0: 1.0}).factor(diversity=3)   # no 'time'


def test_time_validation():
    with pytest.raises(ValueError):
        mod.OnTime({})                      # empty
    with pytest.raises(ValueError):
        mod.OnTime({0: -1.0})               # negative factor
    with pytest.raises(ValueError):
        mod.OnTime({0: float("inf")})       # non-finite factor


def test_time_equality_and_repr():
    assert mod.OnTime({0: 1.0, 3: 0.3}) == mod.OnTime({3: 0.3, 0: 1.0})   # order-independent
    assert mod.OnTime({0: 1.0}) != mod.OnTime({0: 2.0})
    assert "OnTime(" in repr(mod.OnTime({0: 1.0, 3: 0.3}))
    assert hash(mod.OnTime({0: 1.0})) == hash(mod.OnTime({0: 1.0}))


# --- OnTotalDiversity ------------------------------------------------------------

def test_diversity_linear_falloff():
    d = mod.OnTotalDiversity(cap=100)
    assert d.factor(diversity=0) == 1.0
    assert d.factor(diversity=50) == pytest.approx(0.5)
    assert d.factor(diversity=100) == 0.0


def test_diversity_clamps_at_zero_beyond_cap():
    assert mod.OnTotalDiversity(cap=10).factor(diversity=25) == 0.0   # never negative


def test_diversity_ignores_extra_context():
    assert mod.OnTotalDiversity(cap=100).factor(diversity=25, time=3.0) == pytest.approx(0.75)


def test_diversity_missing_context_raises():
    with pytest.raises(TypeError):
        mod.OnTotalDiversity(cap=100).factor(time=1.0)   # no 'diversity'


def test_diversity_validation():
    with pytest.raises(ValueError):
        mod.OnTotalDiversity(cap=0)
    with pytest.raises(ValueError):
        mod.OnTotalDiversity(cap=-5)
    with pytest.raises(ValueError):
        mod.OnTotalDiversity(cap=float("nan"))
    with pytest.raises(TypeError):
        mod.OnTotalDiversity(cap="big")          # type: ignore[arg-type]


def test_diversity_frozen_and_equal():
    d = mod.OnTotalDiversity(cap=100)
    with pytest.raises(Exception):
        d.cap = 50                        # type: ignore[misc]
    assert mod.OnTotalDiversity(cap=100) == mod.OnTotalDiversity(cap=100)
    assert mod.OnTotalDiversity(cap=100) != mod.OnTotalDiversity(cap=50)


# --- the base / the module surface ---------------------------------------

def test_base_modifier_is_abstract():
    with pytest.raises(NotImplementedError):
        mod.Modifier().factor(time=1.0)


def test_stochastic_status_built_vs_deferred():
    for built in ("FromParent", "ByLineage", "ByFamily", "DrivenBy"):
        assert hasattr(mod, built), f"{built} should be built"
    for later in ("Markov",):
        assert not hasattr(mod, later), f"{later} is not built yet"


# --- FromParent (clade drift): the mean-corrected drift ---------------------

def test_inherited_initial_is_one():
    assert mod.FromParent(spread=0.3).initial() == 1.0


def test_inherited_descend_is_mean_corrected():
    import numpy as np
    rng = np.random.default_rng(0)
    inh = mod.FromParent(spread=0.5)
    vals = [inh.descend(1.0, rng) for _ in range(50000)]
    # E[factor] = 1 exactly (the -σ²/2 correction); the buggy version gives E ≈ e^{σ²/2} = 1.13
    assert abs(sum(vals) / len(vals) - 1.0) < 0.02


def test_inherited_no_inflation_over_a_chain():
    import numpy as np
    rng = np.random.default_rng(1)
    inh = mod.FromParent(spread=0.4)
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
    a = mod.FromParent(spread=0.3).descend(1.0, np.random.default_rng(7))
    b = mod.FromParent(spread=0.3).descend(1.0, np.random.default_rng(7))
    assert a == b


def test_binned_drift_stays_on_its_ladder_and_averages_one():
    """The discrete-bin clock: the same inherit-and-perturb model, in steps. A daughter takes its
    parent's rung or one either side, the ends reflect, and the ladder is scaled so the mean is 1 —
    which holds because a reflecting nearest-neighbour walk is uniform at stationarity."""
    import numpy as np
    inh = mod.FromParent(spread=0.35, bins=6)
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
    assert mod.FromParent(spread=0.3).bins is None
    a = mod.FromParent(spread=0.3).descend(1.0, np.random.default_rng(7))
    b = mod.FromParent(spread=0.3, bins=None).descend(1.0, np.random.default_rng(7))
    assert a == b


@pytest.mark.parametrize("bad", [1, 0, -3])
def test_binned_drift_needs_at_least_two_bins(bad):
    with pytest.raises(ValueError, match="at least 2"):
        mod.FromParent(spread=0.3, bins=bad)


def test_a_carried_modifier_has_no_factor_to_give():
    """Its number never came from the context: the engine draws it when a unit is born, keeps it, and
    hands it back through `Rate.effective`'s ``carried``. Asking for a factor without that is
    meaningless, so it raises with the reason rather than returning a plausible 1.0."""
    for m in (mod.FromParent(spread=0.3), mod.ByLineage(spread=0.3), mod.ByFamily(spread=0.3)):
        with pytest.raises(NotImplementedError, match="carried"):
            m.factor()


# --- ByLineage (the uncorrelated / relaxed clock): i.i.d. mean-corrected draws ---

def test_bylineage_zero_spread_is_a_strict_clock():
    import numpy as np
    rng = np.random.default_rng(0)
    byl = mod.ByLineage(spread=0.0)
    assert all(byl.draw(rng) == 1.0 for _ in range(100))


def test_bylineage_draw_is_mean_corrected_lognormal():
    import numpy as np
    rng = np.random.default_rng(0)
    byl = mod.ByLineage(spread=0.5)  # default dist = lognormal
    vals = [byl.draw(rng) for _ in range(100000)]
    # E[factor] = 1 (the -σ²/2 correction); the buggy uncorrected draw gives E ≈ e^{σ²/2} = 1.13
    assert abs(sum(vals) / len(vals) - 1.0) < 0.02


def test_bylineage_draw_is_mean_corrected_gamma():
    import numpy as np
    rng = np.random.default_rng(1)
    byl = mod.ByLineage(spread=0.5, dist="gamma")
    vals = [byl.draw(rng) for _ in range(100000)]
    assert abs(sum(vals) / len(vals) - 1.0) < 0.02          # mean-1 gamma
    var = sum((v - 1.0) ** 2 for v in vals) / len(vals)
    assert abs(var - 0.5 ** 2) < 0.02                        # variance = spread² (CV = spread)


def test_bylineage_draws_are_independent_no_memory():
    import numpy as np
    rng = np.random.default_rng(2)
    byl = mod.ByLineage(spread=0.6)
    a = [byl.draw(rng) for _ in range(2000)]
    # i.i.d.: successive draws are uncorrelated (unlike FromParent, whose draws depend on the parent)
    lag1 = sum((a[i] - 1) * (a[i + 1] - 1) for i in range(len(a) - 1)) / (len(a) - 1)
    assert abs(lag1) < 0.05


def test_bylineage_deterministic():
    import numpy as np
    a = mod.ByLineage(spread=0.3).draw(np.random.default_rng(7))
    b = mod.ByLineage(spread=0.3).draw(np.random.default_rng(7))
    assert a == b


def test_bylineage_validates_its_arguments():
    for bad in (-0.1, float("inf"), float("nan"), True):
        with pytest.raises((ValueError, TypeError)):
            mod.ByLineage(spread=bad)
    with pytest.raises(ValueError):
        mod.ByLineage(spread=0.3, dist="weibull")


def test_inherited_validation():
    with pytest.raises(ValueError):
        mod.FromParent(spread=-0.1)
    with pytest.raises(ValueError):
        mod.FromParent(spread=float("inf"))
    with pytest.raises(TypeError):
        mod.FromParent(spread="wide")  # type: ignore[arg-type]


def test_time_next_change():
    t = mod.OnTime({0: 1.0, 3: 0.3, 7: 0.1})
    assert t.next_change(0.0) == 3
    assert t.next_change(3.0) == 7        # strictly after the current time
    assert t.next_change(5.0) == 7
    assert t.next_change(7.0) == float("inf")  # nothing after the last breakpoint


def test_diversity_never_changes_with_time():
    assert mod.OnTotalDiversity(cap=100).next_change(3.0) == float("inf")
