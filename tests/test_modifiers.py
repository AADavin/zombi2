"""Tests for zombi2.rates.modifiers — the deterministic rate modifiers (SPEC §5)."""

import pytest

from zombi2.rates import modifiers as mod


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


def test_stochastic_status_built_vs_deferred():
    for built in ("Inherited", "ByLineage"):
        assert hasattr(mod, built), f"{built} should be built"
    for later in ("ByFamily", "Markov", "DrivenBy"):
        assert not hasattr(mod, later), f"{later} is not built yet"


# --- Inherited (clade drift): the mean-corrected drift ---------------------

def test_inherited_initial_is_one():
    assert mod.Inherited(spread=0.3).initial() == 1.0


def test_inherited_descend_is_mean_corrected():
    import numpy as np
    rng = np.random.default_rng(0)
    inh = mod.Inherited(spread=0.5)
    vals = [inh.descend(1.0, rng) for _ in range(50000)]
    # E[factor] = 1 exactly (the -σ²/2 correction); the buggy version gives E ≈ e^{σ²/2} = 1.13
    assert abs(sum(vals) / len(vals) - 1.0) < 0.02


def test_inherited_no_inflation_over_a_chain():
    import numpy as np
    rng = np.random.default_rng(1)
    inh = mod.Inherited(spread=0.4)
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
    a = mod.Inherited(spread=0.3).descend(1.0, np.random.default_rng(7))
    b = mod.Inherited(spread=0.3).descend(1.0, np.random.default_rng(7))
    assert a == b


def test_inherited_factor_reads_lineage_multiplier():
    inh = mod.Inherited(spread=0.3)
    assert inh.factor(inherited=2.5, time=1.0) == 2.5
    assert inh.factor() == 1.0  # default: no drift


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
    # i.i.d.: successive draws are uncorrelated (unlike Inherited, whose draws depend on the parent)
    lag1 = sum((a[i] - 1) * (a[i + 1] - 1) for i in range(len(a) - 1)) / (len(a) - 1)
    assert abs(lag1) < 0.05


def test_bylineage_deterministic():
    import numpy as np
    a = mod.ByLineage(spread=0.3).draw(np.random.default_rng(7))
    b = mod.ByLineage(spread=0.3).draw(np.random.default_rng(7))
    assert a == b


def test_bylineage_factor_reads_lineage_multiplier():
    byl = mod.ByLineage(spread=0.3)
    assert byl.factor(bylineage=2.5) == 2.5
    assert byl.factor() == 1.0  # default: no clock


def test_bylineage_validates_its_arguments():
    for bad in (-0.1, float("inf"), float("nan"), True):
        with pytest.raises((ValueError, TypeError)):
            mod.ByLineage(spread=bad)
    with pytest.raises(ValueError):
        mod.ByLineage(spread=0.3, dist="weibull")


def test_inherited_validation():
    with pytest.raises(ValueError):
        mod.Inherited(spread=-0.1)
    with pytest.raises(ValueError):
        mod.Inherited(spread=float("inf"))
    with pytest.raises(TypeError):
        mod.Inherited(spread="wide")  # type: ignore[arg-type]


def test_time_next_change():
    t = mod.Time({0: 1.0, 3: 0.3, 7: 0.1})
    assert t.next_change(0.0) == 3
    assert t.next_change(3.0) == 7        # strictly after the current time
    assert t.next_change(5.0) == 7
    assert t.next_change(7.0) == float("inf")  # nothing after the last breakpoint


def test_diversity_never_changes_with_time():
    assert mod.Diversity(cap=100).next_change(3.0) == float("inf")
