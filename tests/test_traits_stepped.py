"""Two traits feeding back into each other with no exact solution — `zombi2.traits.stepped`.

The correctness-critical test is the **exact agreement**: where both optima are straight lines the
slice map is affine, so its mean and covariance can be iterated in closed form with no Monte Carlo
at all, and both that iteration and the engine's own draws must land on the law the exact
drift-matrix engine gives. Everything else here — the order of the error, the draw order, the
refusals — hangs off that one anchor.

A straight optimum is refused by `simulate_traits`, on purpose: it has an exact answer. The tests
that need one call `stepped.simulate` with ``check_affine=False``, because it is the only case where
the stepped engine and the exact engine describe the same law and so the only case where they can
be compared at all.
"""

import math

import numpy as np
import pytest

from zombi2 import species, traits
from zombi2.params import PerLineage, Scalar
from zombi2.params.connection import set_by
from zombi2.traits import stepped
from zombi2.traits._shared import _ou_transition


def _tree(n=6, seed=11):
    return species.simulate_species_tree(birth=1.0, death=0.0, n_extant=n, seed=seed).complete_tree


# --- the affine case, where the exact answer is known ----------------------------------------------

#: One affine pair, written once: pulls, the slopes and intercepts of the two optimum lines, and σ².
#: ``θ_x(y) = IX + BX·y`` and ``θ_y(x) = IY + BY·x``.
AX, AY, BX, BY, IX, IY, S2 = 1.0, 0.5, 0.5, -0.3, 1.0, 0.0, 0.1
X0 = np.array([0.3, -0.2])


def _exact_law(dt):
    """The mean and covariance the drift-matrix engine gives over a branch of length ``dt``.

    Writing ``θ_x(y) = a_x + b_x·y`` into ``dx = −α_x(x − θ_x(y))dt`` gives ``P_xx = α_x`` and
    ``P_xy = −α_x·b_x``, and the pair's own optima are where the two lines cross. That is the
    algebra the affine refusal prints, so this pins the refusal's arithmetic too."""
    P = np.array([[AX, -AX * BX], [-AY * BY, AY]])
    det = 1.0 - BX * BY
    theta = np.array([(IX + BX * IY) / det, (IY + BY * IX) / det])
    decay, cov = _ou_transition(P, np.diag([S2, S2]), dt)
    return theta + decay @ (X0 - theta), cov


def _iterated_law(dt, step):
    """The stepped scheme's own mean and covariance over one branch, in **closed form**.

    With straight optima the slice map is affine — ``z' = M·z + c + noise`` — so the mean and the
    covariance propagate through it exactly, with no sampling. The map is the engine's: each
    optimum read from the other trait's value at the **start** of the slice, then one exact OU
    transition across the slice. The boundaries come from the engine, so a change to how a branch is
    cut moves this test with it."""
    m, C = X0.copy(), np.zeros((2, 2))
    bounds = stepped._boundaries(0.0, dt, step)
    for a, b in zip(bounds, bounds[1:]):
        s = b - a
        ex, ey = math.exp(-AX * s), math.exp(-AY * s)
        M = np.array([[ex, BX * (1.0 - ex)], [BY * (1.0 - ey), ey]])
        c = np.array([IX * (1.0 - ex), IY * (1.0 - ey)])
        V = np.diag([S2 / (2.0 * AX) * (1.0 - ex * ex), S2 / (2.0 * AY) * (1.0 - ey * ey)])
        m = M @ m + c
        C = M @ C @ M.T + V
    return m, C


def _gap(dt, step):
    """How far the stepped scheme's law sits from the exact one — the worst entry of either."""
    m, C = _iterated_law(dt, step)
    em, eC = _exact_law(dt)
    return max(float(np.abs(m - em).max()), float(np.abs(C - eC).max()))


def _affine_specs(step):
    return [traits.continuous(name="x", start=float(X0[0]), rate=S2, pull=AX,
                              reverts_to=set_by("traits:y", lambda v: IX + BX * v, step=step)),
            traits.continuous(name="y", start=float(X0[1]), rate=S2, pull=AY,
                              reverts_to=set_by("traits:x", lambda v: IY + BY * v, step=step))]


def test_the_slice_map_iterates_to_the_exact_law():
    """**The important one.** At a fine step the scheme's own law is the exact law, to 1e-4."""
    assert _gap(1.3, 1e-3) < 1e-4


def test_the_error_is_first_order_in_the_step():
    """Halving the step halves the error. Two consecutive halvings, so a scheme that happened to
    look first-order at one pair of steps does not pass."""
    for step in (0.08, 0.04):
        ratio = _gap(1.3, step) / _gap(1.3, step / 2.0)
        assert 1.8 <= ratio <= 2.2, (step, ratio)


def test_the_stepped_run_has_the_exact_law_where_the_optimum_is_straight():
    """The engine's own draws, not just its scheme: across many seeds the root branch's end vector
    has the mean and covariance the exact engine gives. The root branch is used because nothing
    comes before it, so its law is the branch transition with no tree behind it."""
    tree = _tree(n=2, seed=3)
    root = tree.root
    dt = tree.nodes[root].end_time - tree.nodes[root].birth_time
    specs = _affine_specs(0.01)
    ends = np.array([[r["x"].node_values[root], r["y"].node_values[root]] for r in (
        stepped.simulate(tree, specs, seed=s, check_affine=False) for s in range(3000))])
    mean, cov = _exact_law(dt)
    se = np.sqrt(np.diag(cov) / len(ends))
    assert np.all(np.abs(ends.mean(axis=0) - mean) < 4 * se)
    assert np.allclose(np.cov(ends.T), cov, rtol=0.12, atol=0.002)


# --- the two models that ship ----------------------------------------------------------------------

def _curved(step=0.01):
    """Case A: a curved optimum, each trait's in the other. The tanh keeps every optimum bounded,
    so the pair settles rather than running away."""
    return [traits.continuous(name="brain", start=0.0, rate=0.1, pull=1.0,
                              reverts_to=set_by("traits:group", lambda v: 2.0 * math.tanh(v),
                                                step=step)),
            traits.continuous(name="group", start=0.0, rate=0.1, pull=0.5,
                              reverts_to=set_by("traits:brain",
                                                lambda v: 1.5 * math.tanh(0.7 * v), step=step))]


def _with_discrete(step=0.02, switch_reads=True):
    """Case C: a discrete trait painting a continuous trait's optimum, and the continuous trait
    driving the discrete trait's switch rate. ``switch_reads=False`` cuts the second direction,
    which leaves the pair conditioned — and then the run takes no step and is exact."""
    switch = (PerLineage(0.6).scaled_by("traits:size", Scalar(0.4), step=step) if switch_reads
              else 0.6)
    return [traits.discrete(name="habitat", states=["surface", "cave"], start="surface",
                            switch=switch),
            traits.continuous(name="size", start=0.0, rate=0.4, pull=1.0,
                              regimes="traits:habitat",
                              reverts_to={"surface": 0.0, "cave": 2.0})]


def test_both_models_run_and_give_one_ordinary_result_per_trait():
    tree = _tree()
    for specs in (_curved(), _with_discrete()):
        out = traits.simulate_traits(tree, specs, joint=True, seed=1)
        assert sorted(out) == sorted(s.name for s in specs)
        for name, res in out.items():
            assert set(res.node_values) == set(tree.nodes)
            first = res.events[0]
            assert (first.kind, first.lineage, first.from_state) == ("initial", tree.root, None)
            if res.kind == "discrete":      # the map covers every branch, as a lone run's does
                for i, node in tree.nodes.items():
                    assert sum(d for _s, d in res.history[i]) == pytest.approx(
                        node.end_time - node.birth_time)


def test_it_is_deterministic():
    tree = _tree()
    a = traits.simulate_traits(tree, _with_discrete(), joint=True, seed=7)
    b = traits.simulate_traits(tree, _with_discrete(), joint=True, seed=7)
    other = traits.simulate_traits(tree, _with_discrete(), joint=True, seed=8)
    for name in a:
        assert a[name].node_values == b[name].node_values
        assert [(e.time, e.kind, e.lineage, e.from_state, e.to_state) for e in a[name].events] \
            == [(e.time, e.kind, e.lineage, e.from_state, e.to_state) for e in b[name].events]
        assert a[name].node_values != other[name].node_values


def test_halving_the_step_changes_the_numbers():
    """By design, and written down so nobody later "fixes" it. A finer step takes more draws, so
    the same seed lands somewhere else: the two runs are samples from two nearby laws, not one
    answer refined. Convergence is checked across seeds, never within one."""
    tree = _tree()
    coarse = traits.simulate_traits(tree, _curved(step=0.02), joint=True, seed=3)
    fine = traits.simulate_traits(tree, _curved(step=0.01), joint=True, seed=3)
    assert coarse["brain"].node_values != fine["brain"].node_values


def test_the_draw_order_holds_when_a_still_trait_joins():
    """Adding a discrete trait that never switches must not move the continuous traits' values.

    That is what fixes the draw order: the full ``k``-vector every segment, the discrete walks
    before it, and a generator of all zeros that draws nothing. The engine is called directly
    because it is linear in the number of traits while the front door keeps a two-trait scope."""
    tree = _tree()
    specs = _curved()
    still = traits.discrete(name="frozen", states=["a", "b"], start="a", switch=0.0)
    without = stepped.simulate(tree, specs, seed=5)
    with_it = stepped.simulate(tree, [*specs, still], seed=5)
    for name in ("brain", "group"):
        assert without[name].node_values == with_it[name].node_values
    assert set(with_it["frozen"].values.values()) == {"a"}


# --- case C against the runs it has to agree with --------------------------------------------------

def _tip_summaries(out, tips):
    return ([out["size"].values_by_id[t] for t in tips],
            [out["habitat"].values_by_id[t] == "cave" for t in tips])


def test_the_discrete_side_is_exact_when_its_switch_reads_nothing():
    """With the switch rate not reading the continuous trait, the discrete trait's own marginal is
    the plain Mk law: the generator is constant, so its branch walk is the one `simulate_discrete`
    takes. Checked distributionally, because the draw orders differ."""
    tree, n = _tree(), 1200
    tips = sorted(tree.extant_leaves())
    joint = np.array([[traits.simulate_traits(tree, _with_discrete(switch_reads=False),
                                              joint=True, seed=s)["habitat"].values_by_id[t]
                       == "cave" for t in tips] for s in range(n)])
    alone = np.array([[traits.simulate_discrete(tree, states=["surface", "cave"], start="surface",
                                                switch=0.6, seed=s).values_by_id[t] == "cave"
                       for t in tips] for s in range(n)])
    assert joint.mean() == pytest.approx(alone.mean(), abs=0.04)


def test_case_c_matches_conditioning_when_only_one_direction_is_written():
    """One direction is conditioning, and the joint run has to agree with it: grow the discrete
    trait alone, paint the continuous trait with it, and the pair grown together gives the same
    law. Distributional, for the same reason as above."""
    tree, n = _tree(), 1200
    tips = sorted(tree.extant_leaves())
    theta = {"surface": 0.0, "cave": 2.0}
    together = np.array([_tip_summaries(
        traits.simulate_traits(tree, _with_discrete(switch_reads=False), joint=True, seed=s), tips
    )[0] for s in range(n)])
    apart = []
    for s in range(n):
        painted = traits.simulate_discrete(tree, states=["surface", "cave"], start="surface",
                                           switch=0.6, seed=s)
        apart.append([traits.simulate_continuous(tree, start=0.0, rate=0.4, pull=1.0,
                                                 regimes=painted, reverts_to=theta,
                                                 seed=s + 100000).values_by_id[t] for t in tips])
    apart = np.array(apart)
    assert together.mean() == pytest.approx(apart.mean(), abs=0.08)
    assert together.var() == pytest.approx(apart.var(), rel=0.12)


# --- what is refused, each by name -----------------------------------------------------------------

def test_a_straight_optimum_is_refused_and_names_the_exact_call():
    """Not routed to the exact engine: five probes of a callable are a guess about the whole line,
    and a silent switch would change the draw order and what step means. The numbers it prints are
    the ones `_exact_law` builds, so the arithmetic is pinned by the tests above."""
    with pytest.raises(ValueError) as e:
        traits.simulate_traits(_tree(), _affine_specs(0.01), joint=True, seed=1)
    message = str(e.value)
    assert "straight line" in message and "simulate_continuous(tree" in message
    det = 1.0 - BX * BY
    assert f"'x': {(IX + BX * IY) / det:.10g}" in message      # the crossing point of the two lines
    assert f"('x', 'y'): {-AX * BX:.10g}" in message           # P_xy = −α_x·b_x
    assert f"('y', 'x'): {-AY * BY:.10g}" in message


def test_joint_true_needs_something_to_read():
    with pytest.raises(ValueError, match="none reads another"):
        traits.simulate_traits(_tree(), [
            traits.continuous(name="x", start=0.0, rate=0.1, pull=1.0, reverts_to=0.0),
            traits.continuous(name="y", start=0.0, rate=0.1, pull=1.0, reverts_to=0.0)],
            joint=True, seed=1)


def test_reading_another_trait_needs_joint_true():
    with pytest.raises(ValueError, match="joint=True"):
        traits.simulate_traits(_tree(), _curved(), seed=1)


def test_two_steps_that_disagree_are_refused():
    with pytest.raises(ValueError, match="two resolutions"):
        traits.simulate_traits(_tree(), [
            traits.continuous(name="brain", start=0.0, rate=0.1, pull=1.0,
                              reverts_to=set_by("traits:group", lambda v: 2.0 * math.tanh(v),
                                                step=0.01)),
            traits.continuous(name="group", start=0.0, rate=0.1, pull=0.5,
                              reverts_to=set_by("traits:brain", lambda v: math.tanh(v),
                                                step=0.02))], joint=True, seed=1)


def test_a_missing_step_is_refused():
    with pytest.raises(ValueError, match="needs a step="):
        traits.simulate_traits(_tree(), _curved(step=None), joint=True, seed=1)


def test_three_traits_say_why_they_are_refused():
    """The old limit was the product generator's cost, which grows exponentially. The stepper is
    linear, so the limit is now a scope choice, and the message has to say so."""
    with pytest.raises(NotImplementedError, match="scope choice"):
        traits.simulate_traits(_tree(), [*_curved(),
                                         traits.continuous(name="z", start=0.0, rate=0.1)],
                               joint=True, seed=1)


def test_every_trait_needs_a_name():
    with pytest.raises(ValueError, match="needs a name"):
        traits.simulate_traits(_tree(), [
            traits.continuous(start=0.0, rate=0.1),
            traits.continuous(name="y", start=0.0, rate=0.1)], joint=True, seed=1)


def test_reading_a_name_that_is_not_in_the_run_is_refused():
    with pytest.raises(ValueError, match="not a trait in this run"):
        traits.simulate_traits(_tree(), [
            traits.continuous(name="x", start=0.0, rate=0.1, pull=1.0,
                              reverts_to=set_by("traits:nope", lambda v: v, step=0.01)),
            traits.continuous(name="y", start=0.0, rate=0.1, pull=1.0,
                              reverts_to=set_by("traits:x", lambda v: math.tanh(v), step=0.01))],
            joint=True, seed=1)


def test_a_variance_rate_reading_a_live_trait_says_it_is_not_built():
    """The third case of issue #456, and the one Adrián set aside: two traits setting each other's
    σ². Not meaningless — not built, and the message has to say which."""
    with pytest.raises(ValueError, match="not built"):
        traits.simulate_traits(_tree(), [
            traits.continuous(name="x", start=0.0,
                              rate=PerLineage(0.1).scaled_by("traits:y", Scalar(0.3), step=0.01)),
            traits.continuous(name="y", start=0.0, rate=0.1, pull=1.0,
                              reverts_to=set_by("traits:x", lambda v: math.tanh(v), step=0.01))],
            joint=True, seed=1)


def test_regimes_painted_by_a_continuous_trait_is_refused():
    with pytest.raises(ValueError, match="Regimes are painted by a discrete trait"):
        traits.simulate_traits(_tree(), [
            traits.continuous(name="x", start=0.0, rate=0.1, pull=1.0, regimes="traits:y",
                              reverts_to={"a": 0.0}),
            traits.continuous(name="y", start=0.0, rate=0.1, pull=1.0,
                              reverts_to=set_by("traits:x", lambda v: math.tanh(v), step=0.01))],
            joint=True, seed=1)


def test_an_optimum_missing_a_regime_state_is_refused():
    with pytest.raises(ValueError, match="missing an optimum"):
        traits.simulate_traits(_tree(), [
            traits.discrete(name="habitat", states=["surface", "cave", "canopy"],
                            switch=PerLineage(0.3).scaled_by("traits:size", Scalar(0.2),
                                                             step=0.02)),
            traits.continuous(name="size", start=0.0, rate=0.4, pull=1.0,
                              regimes="traits:habitat",
                              reverts_to={"surface": 0.0, "cave": 2.0})], joint=True, seed=1)


def test_a_driven_optimum_is_refused_where_the_tree_is_still_growing():
    """The spec is a bundle; each runner declares what it takes. A tree being grown has no trait
    beside it to read, so `joint.simulate` refuses what `simulate_traits` takes."""
    from zombi2 import joint

    with pytest.raises(ValueError, match="reverts_to=set_by"):
        joint.simulate(species.birth_death(
            birth=PerLineage(0.6).scaled_by("trait", Scalar(0.2), step=0.05), n_extant=10),
            traits.continuous(start=0.0, rate=0.1, pull=1.0,
                              reverts_to=set_by("traits:other", lambda v: v, step=0.01)), seed=1)
    with pytest.raises(ValueError, match="regimes="):
        joint.simulate(species.birth_death(
            birth=PerLineage(0.6).scaled_by("trait", Scalar(0.2), step=0.05), n_extant=10),
            traits.continuous(start=0.0, rate=0.1, pull=1.0, regimes="traits:other",
                              reverts_to={"a": 1.0}), seed=1)
