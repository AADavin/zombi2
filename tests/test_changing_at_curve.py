"""A curve of time read in steps: ``changing_at(Curve(f), step=...)`` (issue #502).

The one claim under test is that it **is** the schedule a reader would have built by hand — one
breakpoint per step, each interval holding the curve at its midpoint — so every check compares the
two, down to the same seed giving the same run at each level that takes a schedule.
"""

import math
import pickle

import pytest

from zombi2 import genomes, species, traits
from zombi2.params import Curve, Extent, PerCopy, PerLineage, Time
from zombi2.params import driver as drv

K = 2.0
STEP = 0.01


def decay(t):
    return math.exp(-K * t)


def by_hand(f=decay, step=STEP, intervals=3000):
    """The schedule issue #502 asks for, written the way a user writes it today."""
    return {i * step: f((i + 0.5) * step) for i in range(intervals)}


# --- the modifier, against the schedule it stands for ---------------------------------------------

def test_each_interval_holds_the_curve_at_its_midpoint():
    m = drv.OnTimeCurve(Curve(decay), 0.1)
    assert m.factor(time=0.0) == decay(0.05)
    assert m.factor(time=0.0999) == decay(0.05)
    assert m.factor(time=0.1) == decay(0.15000000000000002)     # (1 + 0.5) * 0.1, as written by hand
    assert m.factor(time=0.25) == decay(0.25)
    # the value at an interval's start or end is a different number, so reading there would fail
    assert m.factor(time=0.25) != decay(0.2) and m.factor(time=0.25) != decay(0.30000000000000004)


def test_the_step_sets_the_breakpoints():
    assert drv.OnTimeCurve(Curve(decay), 0.1).next_change(0.25) == pytest.approx(0.3)
    assert drv.OnTimeCurve(Curve(decay), 0.5).next_change(0.25) == 0.5
    assert drv.OnTimeCurve(Curve(decay), 0.1).factor(time=0.32) == decay(0.35000000000000003)
    assert drv.OnTimeCurve(Curve(decay), 0.5).factor(time=0.32) == decay(0.25)


@pytest.mark.parametrize("step", [0.01, 0.1, 0.07, 1 / 3, 2.5])
def test_it_answers_every_instant_as_the_hand_built_schedule_does(step):
    """Both questions an engine asks, at the breakpoints themselves and between them. The
    breakpoints are the hard instants: ``0.29 / 0.01`` is ``28.999999999999996``, and an index taken
    from that division alone would answer from the interval before."""
    curve = drv.OnTimeCurve(Curve(decay), step)
    hand = drv.OnTime(by_hand(step=step, intervals=1200))
    instants = [i * step for i in range(1000)] + [i * step + 0.37 * step for i in range(1000)]
    for t in instants:
        assert curve.factor(time=t) == hand.factor(time=t), t
        assert curve.next_change(t) == hand.next_change(t), t


def test_before_time_zero_it_behaves_as_a_schedule_does():
    curve, hand = drv.OnTimeCurve(Curve(decay), STEP), drv.OnTime(by_hand())
    for t in (-3.0, -1e-12):
        assert curve.factor(time=t) == hand.factor(time=t) == decay(STEP / 2)
        assert curve.next_change(t) == hand.next_change(t) == 0.0
    assert curve.next_change(-math.inf) == 0.0      # how an engine asks "does it change at all?"


def test_a_curve_is_read_once_per_interval():
    calls = []

    def counted(t):
        calls.append(t)
        return 1.0

    m = drv.OnTimeCurve(Curve(counted), 0.1)
    for t in (0.01, 0.02, 0.03, 0.11, 0.12):
        m.factor(time=t)
    assert calls == [0.05, 0.15000000000000002]


def test_a_bound_on_the_curve_is_kept():
    m = drv.OnTimeCurve(Curve(lambda t: 10.0, bound=4.0), 0.1)
    assert m.factor(time=0.5) == 4.0


# --- the verbs ------------------------------------------------------------------------------------

def test_changing_at_takes_a_curve_or_a_bare_function():
    for written in (Curve(decay), decay):
        rate = PerLineage(2.0).changing_at(written, step=STEP)
        assert rate.effective(time=0.295, lineages=1) == 2.0 * decay(0.295)
        assert rate.next_change(0.295) == 30 * STEP


def test_set_by_time_takes_a_curve_and_holds_the_rate_itself():
    rate = PerLineage().set_by(Time(), lambda t: 2.0 * decay(t), step=STEP)
    assert rate.effective(time=0.295, lineages=1) == 2.0 * decay(0.295)


def test_the_written_form_records_the_curve_and_its_step():
    assert repr(PerLineage(2.0).changing_at(Curve(decay), step=0.01)) \
        == "PerLineage(2.0).changing_at(Curve(decay), step=0.01)"
    assert repr(PerLineage().set_by(Time(), decay, step=0.01)) \
        == "PerLineage().set_by(Time(), Curve(decay), step=0.01)"
    assert repr(Extent(500).changing_at(decay, step=0.5)) \
        == "Extent(Geometric(500.0)).changing_at(Curve(decay), step=0.5)"


def test_two_rates_are_equal_when_function_and_step_are():
    a = PerLineage(2.0).changing_at(decay, step=0.01)
    assert a == PerLineage(2.0).changing_at(Curve(decay), step=0.01)
    assert a != PerLineage(2.0).changing_at(decay, step=0.02)
    assert a != PerLineage(2.0).changing_at(lambda t: decay(t), step=0.01)
    assert a != PerLineage(2.0).changing_at(by_hand())
    assert hash(a.modifiers[0]) == hash(PerLineage(1.0).changing_at(decay, step=0.01).modifiers[0])


def test_a_curve_with_no_step_is_refused():
    with pytest.raises(ValueError, match="without a step"):
        PerLineage(2.0).changing_at(Curve(decay))
    with pytest.raises(ValueError, match=r"set_by received .* Write set_by\(Time\(\), Curve"):
        PerLineage().set_by(Time(), decay)


def test_a_schedule_with_a_step_is_refused():
    with pytest.raises(ValueError, match="a schedule together with a step"):
        PerLineage(2.0).changing_at({0: 1.0, 3: 0.3}, step=0.1)


@pytest.mark.parametrize("step", [0, -0.1, math.inf, math.nan, True, "0.1"])
def test_a_step_that_is_not_a_positive_number_is_refused(step):
    with pytest.raises(ValueError, match="step"):
        PerLineage(2.0).changing_at(decay, step=step)


def test_something_that_is_neither_is_refused():
    with pytest.raises(ValueError, match="neither a schedule nor a curve of time"):
        PerLineage(2.0).changing_at(3.0)


def test_scaled_by_time_still_names_the_verb_and_now_the_curve():
    with pytest.raises(ValueError, match=r"changing_at\(Curve"):
        PerLineage(2.0).scaled_by(Time(), decay, step=0.01)


def test_a_curve_that_goes_negative_is_refused_when_it_is_read():
    rate = PerLineage(2.0).changing_at(lambda t: 1.0 - t, step=0.5)
    assert rate.effective(time=0.2, lineages=1) == 1.5
    with pytest.raises(ValueError, match="finite non-negative"):
        rate.effective(time=1.2, lineages=1)


# --- written out for a worker process -------------------------------------------------------------

def test_written_out_it_is_the_same_schedule_up_to_the_end():
    curve = drv.OnTimeCurve(Curve(decay), STEP)
    out = curve.as_schedule(1.234)
    assert type(out) is drv.OnTime
    for t in [i * 0.001 for i in range(1235)]:
        assert out.factor(time=t) == curve.factor(time=t)
        assert out.next_change(t) == curve.next_change(t)


def test_a_rate_is_written_out_only_where_it_carries_a_curve():
    plain = PerCopy(0.3).changing_at({0: 1.0, 1: 0.5})
    assert plain.with_schedules(2.0) is plain
    out = PerCopy(0.3).changing_at(lambda t: decay(t), step=STEP).with_schedules(2.0)
    assert all(type(m) is drv.OnTime for m in out.modifiers)
    pickle.dumps(out)                                   # a lambda would not survive this


# --- the same seed, the same run ------------------------------------------------------------------

@pytest.fixture(scope="module")
def tree():
    return species.simulate_species_tree(birth=1.0, death=0.3, n_extant=20, seed=1).complete_tree


def test_species_tree_is_the_one_the_schedule_grows():
    def slow(t):
        return math.exp(-0.2 * t)

    a = species.simulate_species_tree(birth=PerLineage(3.0).changing_at(by_hand(slow)), death=0.3,
                                      n_extant=30, seed=3)
    b = species.simulate_species_tree(birth=PerLineage(3.0).changing_at(Curve(slow), step=STEP),
                                      death=0.3, n_extant=30, seed=3)
    assert a.complete_tree.to_newick() == b.complete_tree.to_newick()


def _family_run(tree, origination, loss, **more):
    r = genomes.simulate_genomes_family(tree, origination=origination, loss=loss, duplication=0.3,
                                        transfer=0.2, seed=5, **more)
    return r.events


def test_family_genomes_are_the_ones_the_schedule_gives(tree):
    """The issue's own case: origination strongest at the root, decaying as exp(-k t)."""
    hand = _family_run(tree, PerLineage(2.0).changing_at(by_hand()), PerCopy(0.4).changing_at(by_hand()))
    curve = _family_run(tree, PerLineage(2.0).changing_at(Curve(lambda t: math.exp(-K * t)), step=STEP),
                        PerCopy(0.4).changing_at(lambda t: math.exp(-K * t), step=STEP))
    assert len(hand) > 50 and curve == hand


def test_family_genomes_differ_when_the_step_does(tree):
    fine = _family_run(tree, PerLineage(2.0).changing_at(decay, step=STEP), 0.4)
    coarse = _family_run(tree, PerLineage(2.0).changing_at(decay, step=1.0), 0.4)
    assert fine != coarse


def test_parallel_family_run_takes_a_lambda_and_matches_the_schedule(tree):
    """A lambda cannot be sent to a worker, so the engine sends the schedule it stands for."""
    hand = _family_run(tree, PerLineage(2.0).changing_at(by_hand()),
                       PerCopy(0.4).changing_at(by_hand()), parallel=2)
    curve = _family_run(tree, PerLineage(2.0).changing_at(lambda t: math.exp(-K * t), step=STEP),
                        PerCopy(0.4).changing_at(lambda t: math.exp(-K * t), step=STEP), parallel=2)
    assert len(hand) > 50 and curve == hand


def test_ordered_genomes_are_the_ones_the_schedule_gives(tree):
    def run(loss, extent):
        r = genomes.simulate_genomes_ordered(tree, loss=loss, duplication=0.3, inversion=0.5,
                                             inversion_extent=extent, initial_families=20, seed=5)
        return r.events

    hand = run(PerCopy(0.4).changing_at(by_hand()), Extent(3).changing_at(by_hand()))
    curve = run(PerCopy(0.4).changing_at(decay, step=STEP), Extent(3).changing_at(decay, step=STEP))
    assert len(hand) > 20 and curve == hand


def test_nucleotide_genomes_are_the_ones_the_schedule_gives(tree):
    def run(loss, inversion, extent):
        r = genomes.simulate_genomes_nucleotide(tree, loss=loss, inversion=inversion,
                                                inversion_extent=extent, seed=2)
        return [repr(e) for e in r.events]

    hand = run(PerLineage(0.5).changing_at(by_hand()), PerLineage(2.0).changing_at(by_hand()),
               Extent(200).changing_at(by_hand()))
    curve = run(PerLineage(0.5).changing_at(decay, step=STEP), PerLineage(2.0).changing_at(decay, step=STEP),
                Extent(200).changing_at(decay, step=STEP))
    assert len(hand) > 5 and curve == hand


def test_an_early_burst_trait_is_the_one_the_schedule_gives(tree):
    hand = traits.simulate_continuous(tree, start=0.0, rate=PerLineage(1.0).changing_at(by_hand()), seed=1)
    curve = traits.simulate_continuous(tree, start=0.0, rate=PerLineage(1.0).changing_at(decay, step=STEP),
                                       seed=1)
    assert curve.values == hand.values
