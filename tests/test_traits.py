"""Tests for the continuous trait core — zombi2.traits.simulate_continuous (Brownian motion).

The correctness-critical test is the **exact BM tip law** (Felsenstein 1985): the node-by-node
preorder walk must reproduce the multivariate-normal law over the extant tips, so across replicates
each tip has variance σ²·(root-to-tip depth) and each tip pair has covariance σ²·(shared path
length). Both are checked against the tree geometry, with fixed seeds so the statistics are
deterministic, not flaky.
"""

import math

import numpy as np
import pytest

from zombi2.params import Drift, Global, LogNormal, PerLineage, TotalDiversity
from zombi2.params import law as law
from zombi2.params import scope
from zombi2 import species
from zombi2.species import simulate_species_tree
from zombi2.traits import Change, TraitsResult, simulate_continuous, simulate_discrete


def _tree(seed=1, n_extant=12, death=0.3):
    return simulate_species_tree(birth=1.0, death=death, n_extant=n_extant, seed=seed)


def _mrca_split_time(tree, a, b):
    """The time (crown-forward) at which lineages ``a`` and ``b`` last shared an ancestor — the
    end_time of their MRCA node, i.e. the shared root-to-MRCA path length."""
    anc_a = []
    p = tree.nodes[a].parent
    while p is not None:
        anc_a.append(p)
        p = tree.nodes[p].parent
    seen = set(anc_a)
    p = tree.nodes[b].parent
    while p is not None:
        if p in seen:
            return tree.nodes[p].end_time
        p = tree.nodes[p].parent
    raise AssertionError("no common ancestor")  # unreachable on one connected tree


# --- determinism & the trivial laws ---------------------------------------------

def test_deterministic_given_seed():
    sp = _tree(seed=2)
    a = simulate_continuous(sp, start=0.0, rate=1.0, seed=9)
    b = simulate_continuous(sp, start=0.0, rate=1.0, seed=9)
    assert a.node_values == b.node_values


def test_different_seeds_differ():
    sp = _tree(seed=2)
    a = simulate_continuous(sp, rate=1.0, seed=1)
    b = simulate_continuous(sp, rate=1.0, seed=2)
    assert a.node_values != b.node_values


def test_root_branch_diffuses():
    # convention B: `start` is the value at t=0 and the root lineage diffuses over its own branch,
    # so node_values[root] is start + a diffusion (not `start` itself) and averages back to `start`.
    sp = _tree(seed=5, n_extant=6, death=0.0)
    root = sp.complete_tree.root
    vals = [simulate_continuous(sp, start=3.5, rate=1.0, seed=s).node_values[root] for s in range(3000)]
    assert not all(v == 3.5 for v in vals)               # it diffused (not pinned to `start`)
    assert abs(float(np.mean(vals)) - 3.5) < 0.12        # unbiased: averages back to `start`


def test_zero_rate_is_constant_trait():
    # σ² = 0: nothing diffuses, so every node keeps the root value exactly (inheritance in isolation)
    sp = _tree(seed=4)
    r = simulate_continuous(sp, start=2.0, rate=0.0, seed=1)
    assert all(v == 2.0 for v in r.node_values.values())


def test_every_node_valued_including_extinct():
    sp = _tree(seed=3, death=0.6)
    r = simulate_continuous(sp, rate=0.5, seed=1)
    assert set(r.node_values) == set(sp.complete_tree.nodes)  # every node has a value
    extinct = set(sp.complete_tree.extinct_leaves())
    assert extinct and extinct <= set(r.node_values)          # extinct lineages included


def test_accepts_a_result_or_a_bare_tree():
    sp = _tree(seed=7)
    a = simulate_continuous(sp, rate=0.7, seed=1)
    b = simulate_continuous(sp.complete_tree, rate=0.7, seed=1)
    assert a.node_values == b.node_values


# --- the result bundle ----------------------------------------------------------

def test_values_are_the_extant_tips():
    sp = _tree(seed=8)
    r = simulate_continuous(sp, rate=0.5, seed=1)
    extant = set(sp.complete_tree.extant_leaves())
    assert set(r.values_by_id) == extant
    assert all(r.values_by_id[i] == r.node_values[i] for i in extant)


def test_continuous_events_are_only_the_initial_marker():
    sp = _tree(seed=8)
    r = simulate_continuous(sp, rate=0.5, seed=1)
    # a diffusion has no along-branch events; the one row is the t=0 marker (the initial state),
    # which every trait log now carries so it defines where the run started
    assert [e.kind for e in r.events] == ["initial"]
    assert r.events[0].lineage == sp.complete_tree.root and r.events[0].from_state is None
    assert r.kind == "continuous"
    assert r.history is None


def test_write_values_tsv(tmp_path):
    sp = _tree(seed=8)
    r = simulate_continuous(sp, rate=0.5, seed=1)
    r.write(tmp_path, outputs=["values"])
    text = (tmp_path / "trait_values.tsv").read_text(encoding="utf-8")
    lines = text.splitlines()
    assert lines[0] == "node\tkind\ttrait"
    assert len(lines) - 1 == len(r.node_values)               # one row per node: tips, extinct, internal
    ids = {int(line.split("\t")[0][1:]) for line in lines[1:]}  # strip the n/e prefix
    assert ids == set(r.node_values)
    kinds = [line.split("\t")[1] for line in lines[1:]]        # extant/extinct[/unsampled]/ancestor
    assert set(kinds) <= {"extant", "extinct", "unsampled", "ancestor"} and "ancestor" in kinds
    n_leaves = len([n for n in sp.complete_tree.nodes.values() if n.is_leaf])
    assert len(kinds) - kinds.count("ancestor") == n_leaves   # every non-ancestor row is a tip
    n_extant = len(sp.complete_tree.extant_leaves())          # and `kind == "extant"` isolates the
    assert kinds.count("extant") == n_extant                  # observed tips a comparative method wants
    assert "extinct" in kinds                                 # this tree (death=0.3) has extinct tips too


def test_write_rejects_unknown_output(tmp_path):
    sp = _tree(seed=8)
    r = simulate_continuous(sp, rate=0.5, seed=1)
    with pytest.raises(ValueError, match="unknown write outputs"):
        r.write(tmp_path, outputs=["history"])


# --- input validation (what the engine deliberately does not wire) --------------

def test_rejects_an_unknown_modifier():
    # changing_at (early burst), an inherited law (variable-rates BM) and TotalDiversity
    # (diversity-dependent) are wired;
    # any other Modifier is rejected loudly
    from zombi2.params.evaluate import Modifier

    class _Bogus(Modifier):
        def factor(self, **_):
            return 1.0

    sp = _tree(seed=1)
    with pytest.raises(ValueError, match="does not support"):
        simulate_continuous(sp, rate=PerLineage(1.0)._and(_Bogus()), seed=1)


def test_rejects_a_non_default_scope():
    sp = _tree(seed=1)
    with pytest.raises(ValueError, match="scope"):
        simulate_continuous(sp, rate=scope.Global(1.0), seed=1)


def test_rejects_bad_start():
    sp = _tree(seed=1)
    with pytest.raises(ValueError, match="start"):
        simulate_continuous(sp, start="big", rate=1.0, seed=1)


def test_rejects_negative_rate():
    sp = _tree(seed=1)
    with pytest.raises(ValueError):  # the scope base rejects a negative variance-rate
        simulate_continuous(sp, rate=-1.0, seed=1)


# --- the exact BM tip law (Felsenstein 1985): the correctness-critical invariant --

def test_bm_tip_law_variance_and_covariance():
    # a small fixed tree; over many replicates the extant tips are jointly normal with
    # Var(tip) = σ²·(root-to-tip depth) and Cov(tip_a, tip_b) = σ²·(shared root-to-MRCA path).
    sp = _tree(seed=11, n_extant=6, death=0.0)  # Yule → clean ultrametric extant tips
    tree = sp.complete_tree
    tips = sorted(tree.extant_leaves())
    depth = tree.nodes[tips[0]].end_time        # ultrametric: every extant tip shares this depth
    assert all(np.isclose(tree.nodes[i].end_time, depth) for i in tips)
    sigma2 = 2.0

    n_rep = 6000
    data = np.array([
        [simulate_continuous(tree, start=0.0, rate=sigma2, seed=s).node_values[i] for i in tips]
        for s in range(n_rep)
    ])  # (n_rep, n_tips)

    # each tip: mean ≈ start (0), variance ≈ σ²·depth (a sanity check on the marginal law)
    means = data.mean(axis=0)
    variances = data.var(axis=0)
    assert np.allclose(means, 0.0, atol=0.15)
    assert np.allclose(variances, sigma2 * depth, rtol=0.1)

    # each pair: covariance ≈ σ²·(shared root-to-MRCA path). This is what pins the *shared* history
    # (and convention B): the covariance sampling error is set by the tip variances (~0.06 here), not
    # by the covariance size, so an absolute tolerance is the right model — and it stays well under
    # the σ²·stem ≈ 0.4 shift a stem-less engine (convention A) would put on every pair.
    cov = np.cov(data, rowvar=False)
    for a in range(len(tips)):
        for b in range(a + 1, len(tips)):
            expected = sigma2 * _mrca_split_time(tree, tips[a], tips[b])
            assert cov[a, b] == pytest.approx(expected, abs=0.22)


def test_bm_trend_free_mean_is_flat():
    # with no trend the tip mean stays at `start` regardless of depth — the walk is unbiased
    sp = _tree(seed=13, n_extant=5, death=0.0)
    tips = sorted(sp.complete_tree.extant_leaves())
    n_rep = 4000
    vals = np.array([
        [simulate_continuous(sp, start=1.5, rate=1.0, seed=s).node_values[i] for i in tips]
        for s in range(n_rep)
    ])
    assert np.allclose(vals.mean(axis=0), 1.5, atol=0.12)


def test_returns_a_traits_result():
    sp = _tree(seed=1)
    assert isinstance(simulate_continuous(sp, rate=1.0, seed=1), TraitsResult)


# --- OU (reverts_to / pull): the exact mean-reverting transition law ------------

def test_ou_tip_law():
    # OU is a Gauss–Markov process whose per-branch transition composes exactly along the path, so
    # an extant tip at depth T (from t=0, convention B) is Normal(θ + (start−θ)·e^{−αT},
    # σ²/(2α)·(1−e^{−2αT})). This is the correctness-critical check on the OU transition + composition.
    sp = _tree(seed=11, n_extant=6, death=0.0)
    tree = sp.complete_tree
    tips = sorted(tree.extant_leaves())
    T = tree.nodes[tips[0]].end_time
    theta, alpha, sigma2, start = 5.0, 1.2, 2.0, 0.0

    n_rep = 6000
    data = np.array([
        [simulate_continuous(tree, start=start, rate=sigma2, reverts_to=theta, pull=alpha,
                             seed=s).node_values[i] for i in tips]
        for s in range(n_rep)
    ])
    e = math.exp(-alpha * T)
    assert np.allclose(data.mean(axis=0), theta + (start - theta) * e, atol=0.1)
    assert np.allclose(data.var(axis=0), sigma2 / (2 * alpha) * (1 - e * e), rtol=0.12)


def test_ou_reverts_toward_the_optimum():
    # a strong pull on a deep tree drives the tips near θ, far from `start` — the qualitative OU
    # signature that a BM walk (no reversion) would not show.
    sp = _tree(seed=7, n_extant=8, death=0.0)
    theta = 10.0
    tips = simulate_continuous(sp, start=0.0, rate=1.0, reverts_to=theta, pull=3.0, seed=1).values
    m = float(np.mean(list(tips.values())))
    assert abs(m - theta) < 2.0 and m > 5.0     # clustered near θ, pulled well away from start=0


def test_ou_deterministic():
    sp = _tree(seed=2)
    kw = dict(rate=1.0, reverts_to=1.5, pull=0.7, seed=3)
    assert simulate_continuous(sp, **kw).node_values == simulate_continuous(sp, **kw).node_values


def test_ou_needs_both_knobs():
    sp = _tree(seed=1)
    with pytest.raises(ValueError, match="both"):
        simulate_continuous(sp, rate=1.0, reverts_to=2.0, seed=1)     # pull missing
    with pytest.raises(ValueError, match="both"):
        simulate_continuous(sp, rate=1.0, pull=0.5, seed=1)           # reverts_to missing


def test_ou_pull_must_be_positive():
    sp = _tree(seed=1)
    with pytest.raises(ValueError, match="pull"):
        simulate_continuous(sp, rate=1.0, reverts_to=2.0, pull=0.0, seed=1)


def _ou_segment(sigma2, a, b, t1, alpha):
    """The variance an OU process with constant σ² over ``[a, b)`` still carries at ``t1`` — the
    exact ``∫_a^b e^{−2α(t1−s)} σ² ds``. Written out here from the transition law, not taken from
    the engine, so the tests below check the model rather than restate the implementation."""
    return sigma2 * math.exp(-2 * alpha * (t1 - b)) * (1 - math.exp(-2 * alpha * (b - a))) / (2 * alpha)


def test_ou_with_a_time_varying_variance_matches_the_weighted_integral():
    # THE crux of OU × a σ² modifier. Under OU the noise injected at time s has decayed by
    # e^{−α(t1−s)} when the branch ends, so it enters the end-of-branch variance with weight
    # e^{−2α(t1−s)}: the per-branch variance is ∫ e^{−2α(t1−s)}σ²(s)ds, NOT the Brownian ∫σ²(s)ds.
    # One branch [0, 2], σ² = 4 then 0.25 at τ = 0.8, α = 1.3. Reusing the Brownian integral here
    # would give 3.5 — more than twenty times the right answer, and still look like a normal,
    # mean-reverting trait, which is why this test states both numbers.
    T, tau, alpha, theta = 2.0, 0.8, 1.3, 2.0
    tree = _one_branch(T).complete_tree
    rate = PerLineage(4.0).changing_at({0.0: 1.0, tau: 0.0625})          # σ² = 4.0, then 0.25
    vals = np.array([simulate_continuous(tree, start=0.0, rate=rate, reverts_to=theta, pull=alpha,
                                         seed=s).node_values[0] for s in range(6000)])

    expected_var = _ou_segment(4.0, 0.0, tau, T, alpha) + _ou_segment(0.25, tau, T, T, alpha)
    assert expected_var == pytest.approx(0.15136, rel=1e-3)   # the number, pinned
    assert vals.mean() == pytest.approx(theta + (0.0 - theta) * math.exp(-alpha * T), abs=0.06)
    assert vals.var() == pytest.approx(expected_var, rel=0.12)

    brownian = 4.0 * tau + 0.25 * (T - tau)                   # what the unweighted integral gives
    assert brownian == pytest.approx(3.5) and brownian > 20 * expected_var


# --- early burst (a skyline on rate): the exact ∫σ²(t)dt over each branch --

def test_eb_constant_schedule_equals_bm():
    # a single-step schedule with factor 1.0 everywhere is σ² constant → byte-identical to bare BM
    # (same one draw per branch), which pins the integral's constant-rate special case.
    sp = _tree(seed=2)
    a = simulate_continuous(sp, rate=3.0, seed=5)
    b = simulate_continuous(sp, rate=PerLineage(3.0).changing_at({0.0: 1.0}), seed=5)
    assert a.node_values == b.node_values


def test_eb_tip_variance_and_covariance_match_the_integral():
    # early burst: σ²(t) drops from `base` to `base·c` at time τ. The per-branch variance is the
    # exact integral, so a tip at depth T has variance base·(τ + c·(T−τ)) and a pair covariance
    # base·∫_0^s over the shared path to their MRCA split s. This pins the skyline integral.
    sp = _tree(seed=11, n_extant=6, death=0.0)
    tree = sp.complete_tree
    tips = sorted(tree.extant_leaves())
    T = tree.nodes[tips[0]].end_time
    base, c = 2.0, 0.25
    tau = 0.4 * T                                  # guaranteed inside (0, T) so the branch crosses it
    sched = PerLineage(base).changing_at({0.0: 1.0, tau: c})

    n_rep = 6000
    data = np.array([
        [simulate_continuous(tree, start=0.0, rate=sched, seed=s).node_values[i] for i in tips]
        for s in range(n_rep)
    ])

    def integral(upto):                            # ∫_0^upto σ²(t) dt for this two-step skyline
        return base * (1.0 * min(upto, tau) + c * max(0.0, upto - tau))

    assert np.allclose(data.var(axis=0), integral(T), rtol=0.1)   # < base·T: the burst decayed
    assert integral(T) < base * T
    cov = np.cov(data, rowvar=False)
    for a in range(len(tips)):
        for b in range(a + 1, len(tips)):
            expected = integral(_mrca_split_time(tree, tips[a], tips[b]))
            assert cov[a, b] == pytest.approx(expected, abs=0.22)


def test_eb_deterministic():
    sp = _tree(seed=2)
    sched = PerLineage(1.5).changing_at({0.0: 1.0, 2.0: 0.3})
    a = simulate_continuous(sp, rate=sched, seed=4)
    b = simulate_continuous(sp, rate=sched, seed=4)
    assert a.node_values == b.node_values


# --- variable-rates BM (an inherited rate): σ² drifts branch-to-branch ----------

def _kurtosis(col):
    """Pearson kurtosis of a sample (3.0 = Gaussian); computed with numpy, no scipy dependency."""
    x = np.asarray(col, float)
    return float(((x - x.mean()) ** 4).mean() / x.var() ** 2)


def _vrbm_tips(spread, n_rep=2500):
    """Extant-tip values over `n_rep` variable-rates-BM replicates on a fixed 8-tip Yule tree."""
    tree = simulate_species_tree(birth=1.0, death=0.0, n_extant=8, seed=11).complete_tree
    tips = sorted(tree.extant_leaves())
    depth = tree.nodes[tips[0]].end_time
    data = np.array([
        [simulate_continuous(tree, start=0.0, rate=PerLineage(2.0).varying_among('lineages', Drift(LogNormal(0.0, spread))),
                             seed=s).node_values[i] for i in tips] for s in range(n_rep)
    ])
    return data, depth


def test_variable_rates_bm_is_mean_corrected():
    # the correctness-critical property: an inherited law is mean-corrected (E[factor]=1), so a drifting σ²
    # does NOT inflate down the tree — E[tip variance] stays σ²·depth, exactly as plain BM. (A missing
    # mean-correction — a real historical bug elsewhere in the codebase — would blow the variance up.)
    data, depth = _vrbm_tips(0.6)
    assert np.allclose(data.var(axis=0), 2.0 * depth, rtol=0.08)


def test_variable_rates_bm_is_heterogeneous():
    # the drift makes σ² vary branch-to-branch, so a tip is a scale-mixture of Gaussians —
    # leptokurtic (kurtosis > 3). Plain BM (dist=LogNormal(0.0, 0)) is Gaussian (≈ 3). This is what tells the two
    # apart, since the mean-correction keeps their variances equal.
    flat, _ = _vrbm_tips(0.0)
    drift, _ = _vrbm_tips(1.2)
    assert np.mean([_kurtosis(flat[:, j]) for j in range(flat.shape[1])]) < 3.3    # BM: Gaussian
    assert np.mean([_kurtosis(drift[:, j]) for j in range(drift.shape[1])]) > 5.0  # drift: heavy-tailed


def test_variable_rates_composes_with_time():
    # an inherited law ∘ a skyline: the drift factor (E=1) rides on top of the early-burst integral, so
    # E[tip variance] equals the plain EB integral ∫σ²(t)dt.
    tree = simulate_species_tree(birth=1.0, death=0.0, n_extant=8, seed=11).complete_tree
    tips = sorted(tree.extant_leaves())
    T = tree.nodes[tips[0]].end_time
    base, c, tau = 2.0, 0.25, 0.4 * T
    rate = PerLineage(base).changing_at({0.0: 1.0, tau: c}).varying_among('lineages', Drift(LogNormal(0.0, 0.8)))
    data = np.array([
        [simulate_continuous(tree, start=0.0, rate=rate, seed=s).node_values[i] for i in tips]
        for s in range(2500)
    ])
    assert np.allclose(data.var(axis=0), base * (1.0 * tau + c * (T - tau)), rtol=0.1)


def test_variable_rates_deterministic():
    sp = _tree(seed=2)
    rate = PerLineage(1.0).varying_among('lineages', Drift(LogNormal(0.0, 0.5)))
    assert simulate_continuous(sp, rate=rate, seed=4).node_values == \
        simulate_continuous(sp, rate=rate, seed=4).node_values


def test_ou_with_an_inherited_rate_keeps_the_ou_variance_and_gets_heavy_tails():
    # an inherited factor is constant along a branch and mean-corrected (E[factor] = 1), so it factors
    # straight out of the pull-weighted integral: E[tip variance] is the plain-OU tip variance
    # σ²(1−e^{−2αT})/(2α), unchanged. The variance therefore cannot tell the two models apart — the
    # drifting σ² makes a tip a scale-mixture of Gaussians instead, so what does tell them apart is
    # the kurtosis. Checking only the variance is the trap the plain vrBM tests were written around.
    tree = _corr_tree()
    tips = sorted(tree.extant_leaves())
    T = tree.nodes[tips[0]].end_time
    sigma2, alpha = 2.0, 0.9

    def tips_under(rate, n_rep=2500):
        return np.array([[simulate_continuous(tree, start=0.0, rate=rate, reverts_to=0.0, pull=alpha,
                                              seed=s).node_values[i] for i in tips]
                         for s in range(n_rep)])

    drift = tips_under(PerLineage(sigma2).varying_among('lineages', Drift(LogNormal(0.0, 0.8))))
    flat = tips_under(sigma2)
    expected = sigma2 / (2 * alpha) * (1 - math.exp(-2 * alpha * T))
    assert np.allclose(drift.var(axis=0), expected, rtol=0.12)
    assert np.allclose(flat.var(axis=0), expected, rtol=0.12)
    assert np.mean([_kurtosis(flat[:, j]) for j in range(flat.shape[1])]) < 3.3     # OU: Gaussian
    assert np.mean([_kurtosis(drift[:, j]) for j in range(drift.shape[1])]) > 4.0   # drift: heavy tails


def test_two_inherited_drifts_compose_rather_than_being_refused():
    """SPEC §5's rule is *one memory structure per axis*, not *one modifier*: two inherited factors
    are the same structure and multiply, exactly as any two modifiers do. (Two lognormal drifts are
    one lognormal drift of the combined spread, so this is redundant rather than useful — but it is
    a composition, and refusing it was a stricter rule than the spec has.) Mixing an inherited factor
    with a drawn one is what raises, and that check now lives in one place for every level."""
    sp = _tree(seed=1)
    res = simulate_continuous(
        sp, rate=PerLineage(1.0).varying_among('lineages', Drift(LogNormal(0.0, 0.2))).varying_among('lineages', Drift(LogNormal(0.0, 0.3))), seed=1)
    assert len(res.node_values) == len(sp.complete_tree.nodes)

    one = simulate_continuous(sp, rate=PerLineage(1.0).varying_among('lineages', Drift(LogNormal(0.0, 0.2))), seed=1)
    assert res.node_values != one.node_values          # the second drift really is in the run


def test_bm_unchanged_by_the_inherited_wiring():
    # a bare rate carries no so it must draw no extra rng and stay byte-identical to slice 1
    sp = _tree(seed=3, death=0.4)
    a = simulate_continuous(sp, start=0.0, rate=1.5, seed=1)
    # reproduced from an independent run — the plain-BM path is untouched by the drift threading
    b = simulate_continuous(sp.complete_tree, start=0.0, rate=1.5, seed=1)
    assert a.node_values == b.node_values


# --- diversity-dependent BM (TotalDiversity on rate): σ² slows as the clade fills ----

def _ltt_integral(tree, cap, upto, pull=0.0):
    """∫_0^{upto} max(0, 1 − LTT(t)/cap) dt, LTT = lineages alive at t — computed independently of the
    engine (a re-derivation, so it validates the engine rather than restating it).

    ``pull`` (the OU strength α) puts the mean-reversion weight e^{−2α(upto−t)} inside the integral,
    which is the whole difference between the Brownian and the OU accrual. The same event sweep
    serves both, since the standing diversity is piecewise-constant either way."""
    ev = []
    for n in tree.nodes.values():
        ev.append((n.birth_time, 1))
        ev.append((n.end_time, -1))
    ev.sort()
    total, div, t_prev = 0.0, 0, 0.0
    for t, d in ev:
        hi = min(t, upto)
        if hi > t_prev:
            factor = max(0.0, 1.0 - div / cap)
            if pull > 0.0:
                total += _ou_segment(factor, t_prev, hi, upto, pull)
            else:
                total += factor * (hi - t_prev)
        div += d
        t_prev = t
        if t_prev >= upto:
            break
    return total


def test_diversity_dependence_matches_the_ltt_integral():
    # the correctness-critical check: σ² is scaled by (1 − LTT(t)/cap), so a tip's variance is the
    # exact integral base·∫(1−LTT/cap)dt over its path, and a pair's covariance the integral to their
    # MRCA. Verified against an independent LTT integrator; suppressed below plain BM's σ²·depth.
    sp = _tree(seed=11, n_extant=8, death=0.0)
    tree = sp.complete_tree
    tips = sorted(tree.extant_leaves())
    T = tree.nodes[tips[0]].end_time
    base, cap = 2.0, 6.0
    data = np.array([
        [simulate_continuous(tree, start=0.0, rate=PerLineage(base).scaled_by(TotalDiversity(cap=cap)), seed=s).node_values[i]
         for i in tips] for s in range(2500)
    ])
    assert np.allclose(data.var(axis=0), base * _ltt_integral(tree, cap, T), rtol=0.1)
    assert base * _ltt_integral(tree, cap, T) < base * T          # diversity-dependence suppresses σ²
    cov = np.cov(data, rowvar=False)
    for a in range(len(tips)):
        for b in range(a + 1, len(tips)):
            s = _mrca_split_time(tree, tips[a], tips[b])
            assert cov[a, b] == pytest.approx(base * _ltt_integral(tree, cap, s), abs=0.22)


def test_diversity_dependence_freezes_at_a_small_cap():
    # a small cap → σ² hits 0 once the standing diversity reaches it → far more suppression than a
    # cap the tree never approaches.
    sp = _tree(seed=11, n_extant=8, death=0.0)
    tree = sp.complete_tree
    tips = sorted(tree.extant_leaves())

    def tip_var(cap):
        d = np.array([[simulate_continuous(tree, rate=PerLineage(1.0).scaled_by(TotalDiversity(cap=cap)), seed=s).node_values[i]
                       for i in tips] for s in range(1000)])
        return d.var(axis=0).mean()

    assert tip_var(3.0) < 0.5 * tip_var(200.0)     # a tight cap chokes evolution


def test_diversity_composes_with_inherited():
    # TotalDiversity ∘ an inherited law: the drift factor (E=1) rides on the diversity-scaled integral, so
    # E[tip variance] equals the plain diversity integral.
    sp = _tree(seed=11, n_extant=8, death=0.0)
    tree = sp.complete_tree
    tips = sorted(tree.extant_leaves())
    T = tree.nodes[tips[0]].end_time
    base, cap = 2.0, 6.0
    rate = PerLineage(base).scaled_by(TotalDiversity(cap=cap)).varying_among('lineages', Drift(LogNormal(0.0, 0.6)))
    data = np.array([
        [simulate_continuous(tree, start=0.0, rate=rate, seed=s).node_values[i] for i in tips]
        for s in range(2000)
    ])
    assert np.allclose(data.var(axis=0), base * _ltt_integral(tree, cap, T), rtol=0.1)


def test_diversity_deterministic():
    sp = _tree(seed=2)
    rate = PerLineage(1.0).scaled_by(TotalDiversity(cap=10.0))
    assert simulate_continuous(sp, rate=rate, seed=4).node_values == \
        simulate_continuous(sp, rate=rate, seed=4).node_values


def test_ou_with_diversity_dependence_matches_the_weighted_ltt_integral():
    # OU × TotalDiversity: σ² is scaled by (1 − LTT(t)/cap) *and* the pull discounts what was
    # accrued early, so a tip's variance is base·∫e^{−2α(T−t)}(1−LTT(t)/cap)dt. Both effects
    # suppress it, and they are not interchangeable — the unweighted diversity integral is asserted
    # here too, an order of magnitude out, so a run that forgot the weight cannot pass.
    tree = _corr_tree()
    tips = sorted(tree.extant_leaves())
    T = tree.nodes[tips[0]].end_time
    base, cap, alpha = 2.0, 6.0, 0.8
    data = np.array([
        [simulate_continuous(tree, start=0.0, rate=PerLineage(base).scaled_by(TotalDiversity(cap=cap)),
                             reverts_to=0.0, pull=alpha, seed=s).node_values[i] for i in tips]
        for s in range(2500)
    ])
    expected = base * _ltt_integral(tree, cap, T, pull=alpha)
    assert np.allclose(data.var(axis=0), expected, rtol=0.12)
    assert base * _ltt_integral(tree, cap, T) > 10 * expected     # the unweighted integral is not it


# --- discrete traits (Mk): the exact CTMC + stochastic character map ------------

def _one_branch(total_time, seed=0):
    """A single lineage of length `total_time` (birth=0 → never splits): a clean single branch on
    which to check the exact transition law."""
    return simulate_species_tree(birth=0.0, total_time=total_time, seed=seed)


def test_discrete_deterministic():
    sp = _tree(seed=2)
    kw = dict(states=["a", "b", "c"], switch=0.4, start="a", seed=5)
    a, b = simulate_discrete(sp, **kw), simulate_discrete(sp, **kw)
    assert a.node_values == b.node_values and a.history == b.history


def test_discrete_zero_switch_is_constant():
    sp = _tree(seed=3, death=0.5)
    r = simulate_discrete(sp, states=["x", "y"], switch=0.0, start="x", seed=1)
    assert all(v == "x" for v in r.node_values.values())   # no rate → never leaves the start state
    assert [e.kind for e in r.events] == ["initial"]       # only the t=0 marker, no transitions
    assert r.events[0].to_state == "x"
    assert all(len(segs) == 1 for segs in r.history.values())


def test_discrete_result_shape():
    sp = _tree(seed=8)
    r = simulate_discrete(sp, states=["marine", "terrestrial"], switch=0.3, start="marine", seed=1)
    assert r.kind == "discrete"
    assert set(r.values_by_id) == set(sp.complete_tree.extant_leaves())
    assert set(r.values.values()) <= {"marine", "terrestrial"}     # labels, not indices
    assert set(r.history) == set(sp.complete_tree.nodes)            # a branch history for every node


def test_discrete_history_segments_sum_to_branch_length():
    sp = _tree(seed=4)
    r = simulate_discrete(sp, states=["a", "b", "c"], switch=0.7, start="a", seed=2)
    for i, segs in r.history.items():
        node = sp.complete_tree.nodes[i]
        assert np.isclose(sum(d for _s, d in segs), node.end_time - node.birth_time)
        assert r.node_values[i] == segs[-1][0]                     # end value = last segment's state


def test_discrete_events_track_the_stochastic_map():
    sp = _tree(seed=6)
    r = simulate_discrete(sp, states=["a", "b", "c"], switch=0.8, start="a", seed=3)
    # one event per jump between consecutive segments, across all branches (excluding the t=0 row)
    n_jumps = sum(len(segs) - 1 for segs in r.history.values())
    switches = [e for e in r.events if e.kind != "initial"]
    assert len(switches) == n_jumps and n_jumps > 0
    assert all(isinstance(e, Change) and e.from_state != e.to_state for e in switches)
    assert r.events == sorted(r.events, key=lambda e: e.time)      # time-ordered, initial row first
    assert r.events[0].kind == "initial"
    for e in switches:                                            # each change sits on its branch
        node = sp.complete_tree.nodes[e.lineage]
        assert node.birth_time <= e.time <= node.end_time + 1e-9


def test_discrete_transition_law_two_state_er():
    # the correctness-critical check: on a single branch of length T, a symmetric 2-state chain at
    # rate q ends in the other state with probability (1 − e^{−2qT})/2 (the exact CTMC law). If the
    # Gillespie is right, the empirical switch frequency matches — fixed seeds, so deterministic.
    T, q = 1.5, 0.4
    sp = _one_branch(T, seed=0)
    root = sp.complete_tree.root
    n_rep = 8000
    switched = sum(simulate_discrete(sp, states=["A", "B"], switch=q, start="A",
                                     seed=s).node_values[root] == "B" for s in range(n_rep))
    expected = (1.0 - math.exp(-2 * q * T)) / 2.0
    assert abs(switched / n_rep - expected) < 0.02


def test_discrete_asymmetric_reaches_stationary():
    # asymmetric gain/loss (a dict of rates): on a deep branch the state distribution forgets the
    # start and reaches the stationary π(present) = gain / (gain + loss).
    gain, loss, T = 0.3, 0.1, 20.0
    sp = _one_branch(T, seed=1)
    root = sp.complete_tree.root
    n_rep = 6000
    present = sum(simulate_discrete(sp, states=["absent", "present"],
                                    switch={"absent->present": gain, "present->absent": loss},
                                    start="absent", seed=s).node_values[root] == "present"
                  for s in range(n_rep))
    assert abs(present / n_rep - gain / (gain + loss)) < 0.02


def test_discrete_three_forms_agree():
    # the symmetric scalar, the dict, and the matrix build the SAME Q, so with one seed they give
    # byte-identical histories — one chain, three ways to spell it.
    sp = _tree(seed=7)
    q = 0.5
    scalar = simulate_discrete(sp, states=["A", "B"], switch=q, start="A", seed=9)
    asdict = simulate_discrete(sp, states=["A", "B"],
                               switch={"A->B": q, "B->A": q}, start="A", seed=9)
    matrix = simulate_discrete(sp, states=["A", "B"], switch=[[0.0, q], [q, 0.0]], start="A", seed=9)
    assert scalar.node_values == asdict.node_values == matrix.node_values


def test_a_switch_rate_may_be_written_from_its_scope():
    # every rate in the library is written from its scope, and `switch` is no exception: the two
    # spellings are the same rate, so one seed gives byte-identical histories. Both rate positions
    # take it — the symmetric rate and a {'from->to': rate} entry — which is what the engine already
    # allowed for a rate carrying a verb (`PerLineage(0.4).scaled_by(...)`), so accepting the bare
    # one costs nothing and removes the case where chaining a verb made a spec legal.
    sp = _tree(seed=7)
    events = lambda r: [(c.time, c.kind, c.lineage, c.from_state, c.to_state) for c in r.events]
    bare = simulate_discrete(sp, states=["A", "B"], switch=0.5, start="A", seed=9)
    scoped = simulate_discrete(sp, states=["A", "B"], switch=PerLineage(0.5), start="A", seed=9)
    assert bare.node_values == scoped.node_values and events(bare) == events(scoped)
    d_bare = simulate_discrete(sp, states=["A", "B"], start="A", seed=9,
                               switch={"A->B": 0.5, "B->A": 0.2})
    d_scoped = simulate_discrete(sp, states=["A", "B"], start="A", seed=9,
                                 switch={"A->B": PerLineage(0.5), "B->A": 0.2})
    assert d_bare.node_values == d_scoped.node_values and events(d_bare) == events(d_scoped)


def test_a_switch_rate_the_engine_cannot_mean_is_refused_by_what_it_says():
    # each refusal names the thing that is wrong rather than the shape list, and every spelling a
    # message offers is one that runs — the defect being guarded against is a message that tells you
    # to write what it has just turned away.
    sp = _tree(seed=1)
    with pytest.raises(ValueError, match="switches per lineage"):
        simulate_discrete(sp, states=["a", "b"], switch=Global(0.4), seed=1)
    with pytest.raises(ValueError, match="switches per lineage"):
        simulate_discrete(sp, states=["a", "b"], switch={"a->b": Global(0.4)}, seed=1)
    # a matrix cell is a plain number: every entry is per lineage, so a cell has no scope to state,
    # and the dict is the shape that does take a rate object
    with pytest.raises(ValueError, match="plain numbers"):
        simulate_discrete(sp, states=["a", "b"],
                          switch=[[0.0, PerLineage(0.4)], [0.1, 0.0]], seed=1)
    with pytest.raises(ValueError, match="PerLineage"):
        simulate_discrete(sp, states=["a", "b"], switch="fast", seed=1)


def test_a_joint_run_takes_the_same_switch_spellings():
    # `DiscreteTrait` is the same trait model, bundled unexecuted, so it answers a spec the way
    # `simulate_discrete` does — including the driven one it cannot grow, which it now refuses by
    # saying why rather than by listing shapes.
    from zombi2 import joint, traits

    grow = lambda sw: joint.simulate(species.birth_death(birth=PerLineage(1.0).scaled_by("trait", {"a": 2.0, "b": 1.0}), death=0.1, n_extant=20), traits.discrete(states=["a", "b"], switch=sw), seed=1)
    bare, scoped = grow(0.4), grow(PerLineage(0.4))
    assert bare.complete_tree.to_newick() == scoped.complete_tree.to_newick()
    assert bare.trait.node_values == scoped.trait.node_values
    with pytest.raises(ValueError, match="switches per lineage"):
        grow(Global(0.4))
    with pytest.raises(ValueError, match="one constant matrix"):
        grow(PerLineage(0.4).scaled_by("habitat.tsv", {"a": 2.0}))


def test_discrete_start_none_draws_uniformly():
    # over replicates a None start is uniform over states (checked on a zero-rate chain, so the root
    # state is exactly the drawn start with nothing overwriting it)
    sp = _one_branch(1.0, seed=2)
    root = sp.complete_tree.root
    n_rep = 4000
    a = sum(simulate_discrete(sp, states=["A", "B"], switch=0.0, seed=s).node_values[root] == "A"
            for s in range(n_rep))
    assert abs(a / n_rep - 0.5) < 0.05


def test_discrete_write(tmp_path):
    sp = _tree(seed=8)
    r = simulate_discrete(sp, states=["lo", "hi"], switch=0.6, start="lo", seed=1)
    r.write(tmp_path, outputs=["values", "events"])
    vals = (tmp_path / "trait_values.tsv").read_text(encoding="utf-8").splitlines()
    assert vals[0] == "node\tkind\ttrait"
    assert set(line.split("\t")[2] for line in vals[1:]) <= {"lo", "hi"}   # the trait is the 3rd column
    assert len(vals) - 1 == len(r.node_values)                    # every node, not only the extant tips
    ev = (tmp_path / "trait_events.tsv").read_text(encoding="utf-8").splitlines()
    assert ev[0] == "time\tkind\tlineage\tfrom\tto" and len(ev) - 1 == len(r.events)
    assert ev[1].split("\t")[1] == "initial"                       # the t=0 row comes first


def test_discrete_validation():
    sp = _tree(seed=1)
    with pytest.raises(ValueError, match="at least 2 states"):
        simulate_discrete(sp, states=["only"], switch=0.1, seed=1)
    with pytest.raises(ValueError, match="unique"):
        simulate_discrete(sp, states=["a", "a"], switch=0.1, seed=1)
    with pytest.raises(ValueError, match="switch"):
        simulate_discrete(sp, states=["a", "b"], seed=1)              # switch omitted
    with pytest.raises(ValueError, match="non-negative"):
        simulate_discrete(sp, states=["a", "b"], switch=-0.1, seed=1)
    with pytest.raises(ValueError, match="start"):
        simulate_discrete(sp, states=["a", "b"], switch=0.1, start="z", seed=1)
    with pytest.raises(ValueError, match="from->to"):
        simulate_discrete(sp, states=["a", "b"], switch={"a=>b": 0.1}, seed=1)
    with pytest.raises(ValueError, match="not in states"):
        simulate_discrete(sp, states=["a", "b"], switch={"a->z": 0.1}, seed=1)
    with pytest.raises(ValueError, match="3×3|3x3|shape"):
        simulate_discrete(sp, states=["a", "b", "c"], switch=[[0.0, 0.1], [0.1, 0.0]], seed=1)


# --- correlated continuous traits (the correlation= overlay) --------------------

def _corr_tree():
    return simulate_species_tree(birth=1.0, death=0.0, n_extant=8, seed=11).complete_tree


def test_correlated_tip_correlation_matches_rho():
    # the correctness-critical invariant: per-trait rates + a correlation overlay ρ give a tip-level
    # trait correlation of exactly ρ (independent of the rates and the tree), and each trait's
    # marginal variance is σ²_i·depth (plain BM). ρ=0 recovers independence.
    tree = _corr_tree()
    tip = sorted(tree.extant_leaves())[0]
    depth = tree.nodes[tip].end_time
    for rho in (0.6, -0.5, 0.0):
        a = np.empty(4000)
        b = np.empty(4000)
        for s in range(4000):
            v = simulate_continuous(tree, start={"a": 0.0, "b": 0.0}, rate={"a": 2.0, "b": 0.5},
                                    correlation={("a", "b"): rho}, seed=s).values_by_id[tip]
            a[s], b[s] = v["a"], v["b"]
        assert abs(float(np.corrcoef(a, b)[0, 1]) - rho) < 0.04
        assert np.isclose(a.var(), 2.0 * depth, rtol=0.08)
        assert np.isclose(b.var(), 0.5 * depth, rtol=0.08)


def test_correlated_deterministic():
    tree = _corr_tree()
    kw = dict(start={"x": 0.0, "y": 1.0}, rate={"x": 1.0, "y": 2.0},
              correlation={("x", "y"): 0.3}, seed=7)
    assert simulate_continuous(tree, **kw).node_values == simulate_continuous(tree, **kw).node_values


def test_correlated_result_shape():
    tree = _corr_tree()
    r = simulate_continuous(tree, start={"x": 0.0, "y": 0.0}, rate={"x": 1.0, "y": 1.0},
                            correlation={("x", "y"): 0.5}, seed=1)
    assert set(r.node_values) == set(tree.nodes)                       # every node valued
    assert all(set(v) == {"x", "y"} for v in r.node_values.values())   # each value a per-trait dict
    assert set(r.values_by_id) == set(tree.extant_leaves())
    # the log carries the run's discrete moments (here only the initial row — no jumps were asked
    # for); `history` stays None, because a stochastic character map is a discrete-trait thing
    assert r.kind == "continuous" and r.history is None
    assert [e.kind for e in r.events] == ["initial"]


def test_correlated_write(tmp_path):
    tree = _corr_tree()
    r = simulate_continuous(tree, start={"size": 0.0, "limb": 0.0}, rate={"size": 1.0, "limb": 1.0},
                            correlation={("size", "limb"): 0.4}, seed=1)
    r.write(tmp_path, outputs=["values"])
    lines = (tmp_path / "trait_values.tsv").read_text(encoding="utf-8").splitlines()
    assert lines[0] == "node\tkind\tsize\tlimb"                        # kind, then one column per trait
    assert len(lines) - 1 == len(r.node_values)                        # every node


def test_correlated_needs_psd():
    # three traits cannot all be strongly negatively correlated (the matrix is not PSD)
    tree = _corr_tree()
    with pytest.raises(ValueError, match="positive-semidefinite"):
        simulate_continuous(tree, start={"a": 0.0, "b": 0.0, "c": 0.0},
                            rate={"a": 1.0, "b": 1.0, "c": 1.0},
                            correlation={("a", "b"): -0.9, ("a", "c"): -0.9, ("b", "c"): -0.9}, seed=1)


def test_correlated_validation():
    tree = _corr_tree()
    good = dict(start={"a": 0.0, "b": 0.0}, rate={"a": 1.0, "b": 1.0})
    with pytest.raises(ValueError, match="must be a number in"):
        simulate_continuous(tree, **good, correlation={("a", "b"): 1.5}, seed=1)
    with pytest.raises(ValueError, match="not in"):
        simulate_continuous(tree, **good, correlation={("a", "z"): 0.5}, seed=1)
    with pytest.raises(ValueError, match="self-correlation"):
        simulate_continuous(tree, **good, correlation={("a", "a"): 0.5}, seed=1)
    with pytest.raises(ValueError, match="same traits"):
        simulate_continuous(tree, start={"a": 0.0, "b": 0.0}, rate={"a": 1.0}, seed=1)
    with pytest.raises(ValueError, match="one trait"):
        simulate_continuous(tree, start={"a": 0.0}, rate={"a": 1.0}, seed=1)
    with pytest.raises(ValueError, match="dicts keyed by trait"):
        simulate_continuous(tree, start=0.0, rate=1.0, correlation={("a", "b"): 0.5}, seed=1)
    with pytest.raises(ValueError, match="not implemented yet"):
        simulate_continuous(tree, start={"a": 0.0, "b": 0.0},
                            rate={"a": 1.0, "b": PerLineage(1.0).changing_at({0: 1.0, 3: 0.2})},
                            correlation={("a", "b"): 0.5}, seed=1)
    with pytest.raises(ValueError, match="both"):                      # OU needs the pair, here too
        simulate_continuous(tree, **good, correlation={("a", "b"): 0.5}, reverts_to=1.0, seed=1)
    with pytest.raises(ValueError, match="same traits as start"):
        simulate_continuous(tree, **good, correlation={("a", "b"): 0.5},
                            reverts_to={"a": 1.0}, pull=0.5, seed=1)
    with pytest.raises(ValueError, match="finite positive"):
        simulate_continuous(tree, **good, correlation={("a", "b"): 0.5},
                            reverts_to=1.0, pull={"a": 0.5, "b": 0.0}, seed=1)


# --- multivariate OU: the diagonal-drift restriction ----------------------------

def test_multivariate_ou_transition_law():
    # each trait reverts to its own optimum, the correlation stays in the diffusion. With ONE shared
    # α the OU factor is the same for every entry of the covariance, so it cancels out of the tip
    # correlation: Corr(a, b) is exactly ρ, as under plain correlated BM. That is the cleanest
    # statement that the noise matrix is right — the marginals below pin the reversion itself.
    tree = _corr_tree()
    tip = sorted(tree.extant_leaves())[0]
    T = tree.nodes[tip].end_time
    rho, alpha = 0.6, 0.9
    s2 = {"a": 2.0, "b": 0.5}
    theta = {"a": 3.0, "b": -1.0}
    va, vb = np.empty(4000), np.empty(4000)
    for s in range(4000):
        v = simulate_continuous(tree, start={"a": 0.0, "b": 0.0}, rate=s2,
                                correlation={("a", "b"): rho}, reverts_to=theta, pull=alpha,
                                seed=s).values_by_id[tip]
        va[s], vb[s] = v["a"], v["b"]
    e = math.exp(-alpha * T)
    for got, name in ((va, "a"), (vb, "b")):
        assert got.mean() == pytest.approx(theta[name] * (1 - e), abs=0.1)
        assert got.var() == pytest.approx(s2[name] / (2 * alpha) * (1 - e * e), rel=0.08)
    assert float(np.corrcoef(va, vb)[0, 1]) == pytest.approx(rho, abs=0.04)


def test_multivariate_ou_with_a_per_trait_pull():
    # one strength per trait. The branch covariance is Σ_ij·(1−e^{−(α_i+α_j)dt})/(α_i+α_j), so with
    # unequal α the tip correlation is NOT ρ any more — the two traits forget their shared history
    # at different speeds. An implementation that lazily scaled Σ by a single scalar would still
    # pass the marginal checks, so the last assertion (corr ≠ ρ) is the one that catches it.
    tree = _corr_tree()
    tip = sorted(tree.extant_leaves())[0]
    T = tree.nodes[tip].end_time
    rho = 0.6
    s2 = {"a": 2.0, "b": 0.5}
    alpha = {"a": 1.5, "b": 0.4}
    va, vb = np.empty(6000), np.empty(6000)
    for s in range(6000):
        v = simulate_continuous(tree, start={"a": 0.0, "b": 0.0}, rate=s2,
                                correlation={("a", "b"): rho}, reverts_to={"a": 0.0, "b": 0.0},
                                pull=alpha, seed=s).values_by_id[tip]
        va[s], vb[s] = v["a"], v["b"]
    for got, name in ((va, "a"), (vb, "b")):
        a = alpha[name]
        assert got.var() == pytest.approx(s2[name] / (2 * a) * (1 - math.exp(-2 * a * T)), rel=0.1)
    asum = alpha["a"] + alpha["b"]
    expected_cov = rho * math.sqrt(s2["a"] * s2["b"]) * (1 - math.exp(-asum * T)) / asum
    assert float(np.cov(va, vb)[0, 1]) == pytest.approx(expected_cov, abs=0.05)
    assert abs(float(np.corrcoef(va, vb)[0, 1]) - rho) > 0.02


def test_multivariate_ou_rejects_a_drift_matrix():
    # a full drift matrix — one trait's deviation pulling another — is a different model, not this
    # one written differently, so it is named and refused rather than quietly read as its diagonal.
    tree = _corr_tree()
    with pytest.raises(ValueError, match="full drift matrix"):
        simulate_continuous(tree, start={"a": 0.0, "b": 0.0}, rate={"a": 1.0, "b": 1.0},
                            correlation={("a", "b"): 0.5}, reverts_to={"a": 0.0, "b": 0.0},
                            pull=[[1.0, 0.2], [0.0, 0.5]], seed=1)
    with pytest.raises(ValueError, match="not implemented yet"):
        simulate_continuous(tree, start={"a": 0.0, "b": 0.0}, rate={"a": 1.0, "b": 1.0},
                            correlation={("a", "b"): 0.5}, reverts_to={"a": 0.0, "b": 0.0},
                            pull=[[1.0, 0.2], [0.0, 0.5]], seed=1)


def test_correlated_jumps_at_speciation_add_the_jump_covariance():
    # the jump rides the SAME correlation overlay the diffusion does. A tip's variance is
    # σ²_i·depth + v_i·(splits on its path), and its cross-trait covariance ρ·(σ_aσ_b·depth +
    # √(v_a v_b)·splits) — the second is what pins that the jumps are correlated rather than
    # independent, since the variances alone are the same either way.
    tree = _corr_tree()
    tip = sorted(tree.extant_leaves())[0]
    depth = tree.nodes[tip].end_time
    splits = _n_splits(tree, tip)
    rho, jump = 0.6, {"a": 1.5, "b": 0.5}
    va, vb = np.empty(5000), np.empty(5000)
    for s in range(5000):
        v = simulate_continuous(tree, start={"a": 0.0, "b": 0.0}, rate={"a": 1.0, "b": 1.0},
                                correlation={("a", "b"): rho}, at_speciation=jump,
                                seed=s).values_by_id[tip]
        va[s], vb[s] = v["a"], v["b"]
    assert va.var() == pytest.approx(depth + jump["a"] * splits, rel=0.1)
    assert vb.var() == pytest.approx(depth + jump["b"] * splits, rel=0.1)
    expected_cov = rho * (depth + math.sqrt(jump["a"] * jump["b"]) * splits)
    assert float(np.cov(va, vb)[0, 1]) == pytest.approx(expected_cov, abs=0.25)


def test_correlated_jump_of_zero_leaves_the_run_untouched():
    # at_speciation=0.0 must draw nothing, so a correlated run asking for a zero jump is the run
    # that never asked — the same guarantee the single-trait path gives.
    tree = _corr_tree()
    kw = dict(start={"a": 0.0, "b": 0.0}, rate={"a": 1.0, "b": 1.0},
              correlation={("a", "b"): 0.4}, seed=5)
    assert simulate_continuous(tree, **kw).node_values == \
        simulate_continuous(tree, at_speciation=0.0, **kw).node_values


# --- threshold traits (discrete from a continuous liability) + correlated discrete ---

def _phi(z):
    return 0.5 * (1.0 + math.erf(z / math.sqrt(2.0)))


def test_threshold_state_frequency_law():
    # a 2-state threshold trait: the liability ~ Normal(start, σ²·depth), so the fraction of tips in
    # the upper state is Φ((start − cut)/√(σ²·depth)) — the exact Wright–Felsenstein law.
    tree = _corr_tree()
    tip = sorted(tree.extant_leaves())[0]
    depth = tree.nodes[tip].end_time
    for start, cut, s2 in [(0.0, 0.0, 1.0), (0.5, 0.0, 1.0), (0.0, -0.4, 2.0)]:
        n = 6000
        p = sum(simulate_discrete(tree, states=["absent", "present"], liability=s2, threshold=cut,
                                  start=start, seed=z).values_by_id[tip] == "present" for z in range(n)) / n
        assert abs(p - _phi((start - cut) / math.sqrt(s2 * depth))) < 0.03


def test_threshold_shape_and_multistate():
    tree = _corr_tree()
    r = simulate_discrete(tree, states=["absent", "present"], liability=1.0, threshold=0.0, seed=1)
    assert r.kind == "threshold" and r.history is None and r.events == []   # no timed event log/map
    assert set(r.node_values) == set(tree.nodes)
    assert set(r.values.values()) <= {"absent", "present"}
    m = simulate_discrete(tree, states=["low", "mid", "high"], liability=1.0, threshold=[-1.0, 1.0], seed=1)
    assert set(m.values.values()) <= {"low", "mid", "high"}                # 3 states, 2 increasing cuts


def test_threshold_deterministic():
    tree = _corr_tree()
    kw = dict(states=["a", "b"], liability=1.5, threshold=0.2, start=0.1, seed=4)
    assert simulate_discrete(tree, **kw).node_values == simulate_discrete(tree, **kw).node_values


def test_correlated_discrete_agreement_matches_tetrachoric():
    # correlated liabilities (symmetric, cut at 0): the fraction of tips where the two states agree is
    # 1/2 + asin(ρ)/π (the tetrachoric law) — 0.5 when independent, → 1 as ρ → 1, < 0.5 for ρ < 0.
    tree = _corr_tree()
    tip = sorted(tree.extant_leaves())[0]
    for rho in (0.0, 0.8, -0.6):
        n = 5000
        agree = sum(
            (lambda v: v["w"] == v["f"])(
                simulate_discrete(tree, states=["absent", "present"], liability={"w": 1.0, "f": 1.0},
                                  correlation={("w", "f"): rho}, threshold=0.0, seed=z).values_by_id[tip])
            for z in range(n))
        assert abs(agree / n - (0.5 + math.asin(rho) / math.pi)) < 0.03


def test_correlated_discrete_shape_and_write(tmp_path):
    tree = _corr_tree()
    r = simulate_discrete(tree, states=["absent", "present"], liability={"wings": 1.0, "flight": 1.0},
                          correlation={("wings", "flight"): 0.7}, threshold=0.0, seed=1)
    assert all(set(v) == {"wings", "flight"} for v in r.node_values.values())
    assert r.history is None and r.events == []
    r.write(tmp_path, outputs=["values"])
    assert (tmp_path / "trait_values.tsv").read_text(encoding="utf-8").splitlines()[0] == "node\tkind\twings\tflight"


def test_threshold_validation():
    tree = _corr_tree()
    with pytest.raises(ValueError, match="not both"):
        simulate_discrete(tree, states=["a", "b"], switch=0.1, liability=1.0, threshold=0.0, seed=1)
    with pytest.raises(ValueError, match="cut point"):
        simulate_discrete(tree, states=["a", "b", "c"], liability=1.0, threshold=0.0, seed=1)
    with pytest.raises(ValueError, match="strictly increasing"):
        simulate_discrete(tree, states=["a", "b", "c"], liability=1.0, threshold=[1.0, 0.0], seed=1)
    with pytest.raises(ValueError, match="needs a dict liability"):
        simulate_discrete(tree, states=["a", "b"], liability=1.0, threshold=0.0,
                          correlation={("a", "b"): 0.5}, seed=1)
    with pytest.raises(ValueError, match="threshold model"):
        simulate_discrete(tree, states=["a", "b"], switch=0.1, correlation={("a", "b"): 0.5}, seed=1)
    with pytest.raises(ValueError, match="needs threshold"):
        simulate_discrete(tree, states=["a", "b"], liability=1.0, seed=1)


# --- at_speciation: on-speciation jumps (continuous) / shifts (discrete) ----------

def _n_splits(tree, i):
    c = 0
    p = tree.nodes[i].parent
    while p is not None:
        c += 1
        p = tree.nodes[p].parent
    return c


def test_at_speciation_continuous_adds_jump_variance():
    # a Normal(0, at_speciation) jump at each speciation: a tip's variance is σ²·depth (anagenesis)
    # plus at_speciation·(number of splits on its path) — the punctuational contribution, verified.
    tree = _corr_tree()
    tips = sorted(tree.extant_leaves())
    depth = tree.nodes[tips[0]].end_time
    for jump in (0.0, 1.5):
        vals = np.array([[simulate_continuous(tree, start=0.0, rate=1.0, at_speciation=jump,
                                              seed=s).node_values[i] for i in tips] for s in range(5000)])
        for j, i in enumerate(tips):
            assert np.isclose(vals[:, j].var(), depth + jump * _n_splits(tree, i), rtol=0.1)


def test_at_speciation_off_is_unchanged():
    tree = _corr_tree()
    a = simulate_continuous(tree, rate=1.0, seed=5)
    b = simulate_continuous(tree, rate=1.0, at_speciation=0.0, seed=5)   # jump_sd 0 → no draw
    assert a.node_values == b.node_values


def test_at_speciation_discrete_flips_at_every_split():
    # 2 states, no anagenesis (switch=0), a certain flip at every speciation (at_speciation=1.0): a
    # tip is `start` flipped once per split, so it equals start iff its path has an even split count.
    tree = _corr_tree()
    r = simulate_discrete(tree, states=["a", "b"], switch=0.0, at_speciation=1.0, start="a", seed=1)
    for i in sorted(tree.extant_leaves()):
        assert r.values_by_id[i] == ("a" if _n_splits(tree, i) % 2 == 0 else "b")


def test_at_speciation_deterministic():
    tree = _corr_tree()
    assert simulate_continuous(tree, rate=1.0, at_speciation=0.5, seed=3).node_values == \
        simulate_continuous(tree, rate=1.0, at_speciation=0.5, seed=3).node_values
    assert simulate_discrete(tree, states=["x", "y"], switch=0.3, at_speciation=0.4, seed=3).node_values == \
        simulate_discrete(tree, states=["x", "y"], switch=0.3, at_speciation=0.4, seed=3).node_values


def test_at_speciation_validation():
    tree = _corr_tree()
    with pytest.raises(ValueError, match="non-negative"):
        simulate_continuous(tree, rate=1.0, at_speciation=-1.0, seed=1)
    with pytest.raises(ValueError, match="non-negative"):        # per trait, under correlation
        simulate_continuous(tree, start={"a": 0.0, "b": 0.0}, rate={"a": 1.0, "b": 1.0},
                            correlation={("a", "b"): 0.5}, at_speciation={"a": 1.0, "b": -1.0}, seed=1)
    with pytest.raises(ValueError, match="probability in"):
        simulate_discrete(tree, states=["a", "b"], switch=0.1, at_speciation=1.5, seed=1)
    with pytest.raises(ValueError, match="not implemented for threshold"):
        simulate_discrete(tree, states=["a", "b"], liability=1.0, threshold=0.0, at_speciation=0.5, seed=1)


# --- regimes: multi-optimum OU (optima painted by a discrete stochastic map) -----

def test_regimes_constant_equals_plain_ou():
    # a regime map that never switches (switch=0) is one regime everywhere, so multi-optimum OU
    # collapses to plain OU toward that regime's optimum — byte-identical (one OU draw per branch).
    # This is also the canary on the plain-OU fast path: the main engine now reaches its OU variance
    # through the pull-weighted integral, and this test fails on the last bits unless a bare σ² still
    # takes the closed form σ²/(2α)·(1−e^{−2α·dt}) the regimes path spells out below.
    tree = _corr_tree()
    const = simulate_discrete(tree, states=["a", "b"], switch=0.0, start="a", seed=99)
    got = simulate_continuous(tree, rate=1.0, pull=2.0, reverts_to={"a": 5.0, "b": 0.0}, regimes=const, seed=7)
    plain = simulate_continuous(tree, rate=1.0, pull=2.0, reverts_to=5.0, seed=7)
    assert got.node_values == plain.node_values


def test_regimes_track_their_optima():
    # strong pull → a tip sits near the optimum of the regime it ends in
    tree = simulate_species_tree(birth=1.0, death=0.0, n_extant=12, seed=5).complete_tree
    regime = simulate_discrete(tree, states=["lo", "hi"], switch=0.8, seed=1)
    r = simulate_continuous(tree, rate=0.3, pull=6.0, reverts_to={"lo": 0.0, "hi": 10.0}, regimes=regime, seed=2)
    tips = sorted(tree.extant_leaves())
    lo = [r.values_by_id[i] for i in tips if regime.values_by_id[i] == "lo"]
    hi = [r.values_by_id[i] for i in tips if regime.values_by_id[i] == "hi"]
    if lo:
        assert abs(float(np.mean(lo))) < 3.0
    if hi:
        assert abs(float(np.mean(hi)) - 10.0) < 3.5
    assert (np.mean(hi) if hi else 10) > (np.mean(lo) if lo else 0)


def test_regimes_deterministic():
    tree = _corr_tree()
    regime = simulate_discrete(tree, states=["a", "b"], switch=0.5, seed=1)
    kw = dict(rate=1.0, pull=2.0, reverts_to={"a": 0.0, "b": 3.0}, regimes=regime, seed=4)
    assert simulate_continuous(tree, **kw).node_values == simulate_continuous(tree, **kw).node_values


def test_regimes_take_a_jump_at_speciation():
    # (a) with a regime map that never switches, a multi-optimum run with jumps must be the plain-OU
    # run with jumps, value for value — the jump is drawn in the same place in the same order.
    tree = _corr_tree()
    const = simulate_discrete(tree, states=["a", "b"], switch=0.0, start="a", seed=99)
    got = simulate_continuous(tree, rate=1.0, pull=2.0, reverts_to={"a": 5.0, "b": 0.0},
                              regimes=const, at_speciation=1.2, seed=7)
    plain = simulate_continuous(tree, rate=1.0, pull=2.0, reverts_to=5.0, at_speciation=1.2, seed=7)
    assert got.node_values == plain.node_values

    # (b) on a real two-regime map the jump adds variance at each split on a tip's path — but the
    # pull works on a jump exactly as it works on anagenesis, so a jump taken at time t is worth
    # e^{−2α(T−t)} of its variance by the time the tip is reached. Old jumps count for almost
    # nothing. Measured as the difference between the same run with and without the jump, so
    # everything the regimes contribute cancels.
    alpha = 0.5
    regime = simulate_discrete(tree, states=["lo", "hi"], switch=0.8, seed=1)
    tips = sorted(tree.extant_leaves())
    depth = tree.nodes[tips[0]].end_time
    kw = dict(rate=1.0, pull=alpha, reverts_to={"lo": 0.0, "hi": 4.0}, regimes=regime)
    jump = 1.3

    def tip_vars(**extra):
        d = np.array([[simulate_continuous(tree, seed=s, **kw, **extra).node_values[i] for i in tips]
                      for s in range(4000)])
        return d.var(axis=0)

    def decayed_jumps(i):
        total, node = 0.0, tree.nodes[i]
        while node.parent is not None:                 # one jump per split above this tip
            total += math.exp(-2 * alpha * (depth - node.birth_time))
            node = tree.nodes[node.parent]
        return jump * total

    rise = tip_vars(at_speciation=jump) - tip_vars()
    assert np.allclose(rise, [decayed_jumps(i) for i in tips], rtol=0.2, atol=0.12)


def test_regimes_carry_the_same_event_log_every_continuous_run_carries():
    # a regimes run is a continuous run: the t=0 initial row, plus one on_speciation row per jump.
    # It used to return an empty log, which made Appendix B's Events row untrue for this path.
    tree = _corr_tree()
    regime = simulate_discrete(tree, states=["a", "b"], switch=0.5, seed=1)
    kw = dict(rate=1.0, pull=2.0, reverts_to={"a": 0.0, "b": 3.0}, regimes=regime, seed=4)
    plain = simulate_continuous(tree, **kw)
    assert [e.kind for e in plain.events] == ["initial"]
    assert plain.events[0].lineage == tree.root and plain.events[0].from_state is None
    jumped = simulate_continuous(tree, at_speciation=0.5, **kw)
    jumps = [e for e in jumped.events if e.kind == "on_speciation"]
    assert len(jumps) == len(tree.nodes) - 1                      # one per non-root node, i.e. per split
    assert all(isinstance(e.from_state, float) and isinstance(e.to_state, float) for e in jumps)


def test_regimes_reject_a_per_regime_jump_size():
    tree = _corr_tree()
    regime = simulate_discrete(tree, states=["lo", "hi"], switch=0.5, seed=1)
    with pytest.raises(ValueError, match="not implemented yet"):
        simulate_continuous(tree, rate=1.0, pull=2.0, reverts_to={"lo": 0.0, "hi": 1.0},
                            regimes=regime, at_speciation={"lo": 1.0, "hi": 2.0}, seed=1)


def test_regimes_reject_a_modified_variance_rate():
    # the regimes path walks the map's (state, duration) segments, not absolute time, so it has no
    # place to read a schedule from. The message says that plainly and points at the path that does.
    tree = _corr_tree()
    regime = simulate_discrete(tree, states=["lo", "hi"], switch=0.5, seed=1)
    with pytest.raises(ValueError, match="changing_at.*not implemented yet"):
        simulate_continuous(tree, rate=PerLineage(1.0).changing_at({0: 1.0, 3: 0.2}), pull=2.0,
                            reverts_to={"lo": 0.0, "hi": 1.0}, regimes=regime, seed=1)


def test_regimes_validation():
    tree = _corr_tree()
    regime = simulate_discrete(tree, states=["a", "b"], switch=0.5, seed=1)
    with pytest.raises(ValueError, match="reverts_to is a dict"):
        simulate_continuous(tree, rate=1.0, pull=2.0, reverts_to=1.0, regimes=regime, seed=1)
    with pytest.raises(ValueError, match="needs pull"):
        simulate_continuous(tree, rate=1.0, reverts_to={"a": 0.0, "b": 1.0}, regimes=regime, seed=1)
    with pytest.raises(ValueError, match="missing an optimum"):
        simulate_continuous(tree, rate=1.0, pull=2.0, reverts_to={"a": 0.0}, regimes=regime, seed=1)
    with pytest.raises(ValueError, match="discrete TraitsResult"):
        simulate_continuous(tree, rate=1.0, pull=2.0, reverts_to={"a": 0.0},
                            regimes=simulate_continuous(tree, rate=1.0, seed=1), seed=1)
    with pytest.raises(ValueError, match="SAME tree"):
        other = simulate_discrete(
            simulate_species_tree(birth=1.0, death=0.0, n_extant=6, seed=1).complete_tree,
            states=["a", "b"], switch=0.5, seed=1)
        simulate_continuous(tree, rate=1.0, pull=2.0, reverts_to={"a": 0.0, "b": 1.0}, regimes=other, seed=1)


def test_write_trait_tree(tmp_path):
    # the "tree" output is a Newick with every node annotated [&trait=…] (a trait tree carrying the
    # exact ancestral values) — for continuous (float), discrete (label), and correlated (per-trait).
    import re
    sp = _tree(seed=8, n_extant=6)
    cases = [
        (simulate_continuous(sp, rate=1.0, seed=1), "[&trait="),
        (simulate_discrete(sp, states=["a", "b"], switch=0.4, seed=1), "[&trait="),
        (simulate_continuous(sp, start={"x": 0.0, "y": 0.0}, rate={"x": 1.0, "y": 1.0},
                             correlation={("x", "y"): 0.5}, seed=1), "[&x="),
    ]
    for r, marker in cases:
        r.write(tmp_path, outputs=["tree"])
        nwk = (tmp_path / "trait_tree.nwk").read_text(encoding="utf-8").strip()
        assert nwk.count("(") == nwk.count(")") and nwk.endswith(";")   # structurally valid Newick
        assert marker in nwk                                            # the right annotation shape
        ids = {int(m) for m in re.findall(r"n(\d+)\[", nwk)}
        assert ids == set(r.node_values)                               # every node annotated


# --- the event log (mirrors genomes) + history derived from it ------------------

def test_events_log_records_speciation_changes():
    # on-speciation jumps/shifts now live in the event log (kind="on_speciation"); a plain run has none
    tree = _corr_tree()
    plain = [e for e in simulate_continuous(tree, rate=1.0, seed=1).events if e.kind != "initial"]
    assert plain == []                                                         # pure BM: no switches
    jumps = [e for e in simulate_continuous(tree, rate=1.0, at_speciation=0.5, seed=1).events
             if e.kind != "initial"]
    assert jumps and all(e.kind == "on_speciation" for e in jumps)
    assert isinstance(jumps[0].from_state, float) and isinstance(jumps[0].to_state, float)
    d = [e for e in simulate_discrete(tree, states=["a", "b", "c"], switch=1.5, at_speciation=0.4,
                                      seed=2).events if e.kind != "initial"]
    assert {e.kind for e in d} == {"on_branch", "on_speciation"}                 # both in the log
    assert all(e.from_state != e.to_state for e in d)


def test_events_are_time_sorted():
    d = simulate_discrete(_corr_tree(), states=["a", "b"], switch=2.0, seed=3)
    times = [e.time for e in d.events]
    assert times == sorted(times)


def test_history_is_derived_from_events():
    # history is DERIVED from the event log (discrete only): each branch's segments sum to the branch
    # length, end at node_values, and their transitions reproduce exactly the log's on-branch events.
    tree = _corr_tree()
    d = simulate_discrete(tree, states=["a", "b", "c"], switch=1.5, at_speciation=0.3, seed=2)
    for i, segs in d.history.items():
        node = tree.nodes[i]
        assert abs(sum(dur for _, dur in segs) - (node.end_time - node.birth_time)) < 1e-9
        assert segs[-1][0] == d.node_values[i]
    on_branch = sorted((e.lineage, e.from_state, e.to_state) for e in d.events if e.kind == "on_branch")
    from_hist = sorted((i, segs[k][0], segs[k + 1][0])
                       for i, segs in d.history.items() for k in range(len(segs) - 1))
    assert on_branch == from_hist


# --- Driven: one trait driving another, on the same tree ----------------------

def _write_driver(tmp_path, rows, name="driver.tsv"):
    """A trait event log written by hand — the file a conditioned ``Driven`` names. Rows are
    ``(time, kind, lineage, from, to)``, exactly what ``TraitsResult.write(outputs=("events",))``
    puts on disk, so a test can place a driver switch at an exact time."""
    path = tmp_path / name
    text = "time\tkind\tlineage\tfrom\tto\n" + "".join(
        f"{t}\t{kind}\t{lin}\t{frm}\t{to}\n" for t, kind, lin, frm, to in rows)
    path.write_text(text, encoding="utf-8")
    return str(path)


def test_a_driven_variance_integrates_across_the_drivers_mid_branch_switch(tmp_path):
    # THE crux of a driven rate. One branch, [0, 2]: the driver is 'slow' over its first half and
    # 'fast' over its second, so the accrued variance is the integral across the two pieces —
    # 1.0×1 + 9.0×1 = 10 — not the 1.0×2 = 2 a single sample at the branch start would give (nor the
    # 9.0×2 = 18 a sample at the end would). A per-branch sample is not merely coarse here: it is a
    # different model, and on this branch it is off by a factor of five.
    from zombi2.params.conditioned import resolve_driver
    from zombi2.params.parameter import as_rate
    from zombi2.params.scope import PerLineage
    from zombi2.traits.continuous import _accrued_variance

    tree = _one_branch(2.0).complete_tree
    src = _write_driver(tmp_path, [(0.0, "initial", "n0", "", "slow"),
                                   (1.0, "on_branch", "n0", "slow", "fast")])
    rate = PerLineage(1.0).scaled_by(src, {"slow": 1.0, "fast": 9.0})
    r = as_rate(rate, default_scope=PerLineage)
    trajs = {rate.modifiers[0].key: resolve_driver(src, tree)}

    stepped = _accrued_variance(r, 0.0, 2.0, trajs=trajs, node_id=0)
    assert stepped == pytest.approx(1.0 * 1.0 + 9.0 * 1.0)
    at_start = r.effective(lineages=1, time=0.0, drivers={rate.modifiers[0].key: "slow"}) * 2.0
    assert at_start == pytest.approx(2.0) and abs(stepped - at_start) > 1.0   # the sample misses it

    # and the engine realises exactly that variance over the branch
    drawn = np.array([simulate_continuous(tree, rate=rate, seed=s).node_values[0] for s in range(1500)])
    assert drawn.std() == pytest.approx(math.sqrt(10.0), rel=0.06)


def test_ou_with_a_driven_variance_steps_where_the_driver_steps(tmp_path):
    # a driven σ² under an OU pull. One branch [0, 2], the driver switching slow→fast at t=1 with
    # factors 1.0 and 9.0, α = 1.0. Where the fast stretch sits matters twice over: it contributes
    # more variance, AND it sits at the end of the branch, where the pull has had least time to undo
    # it. The variance is 1.0·∫_0^1 e^{−2(2−s)}ds + 9.0·∫_1^2 e^{−2(2−s)}ds ≈ 3.949; a driver
    # sampled once at the branch start would give 1.0·(1−e^{−4})/2 ≈ 0.49, wrong by a factor of 8.
    alpha, T = 1.0, 2.0
    tree = _one_branch(T).complete_tree
    src = _write_driver(tmp_path, [(0.0, "initial", "n0", "", "slow"),
                                   (1.0, "on_branch", "n0", "slow", "fast")])
    rate = PerLineage(1.0).scaled_by(src, {"slow": 1.0, "fast": 9.0})

    expected = _ou_segment(1.0, 0.0, 1.0, T, alpha) + _ou_segment(9.0, 1.0, T, T, alpha)
    assert expected == pytest.approx(3.9495, rel=1e-3)
    drawn = np.array([simulate_continuous(tree, start=0.0, rate=rate, reverts_to=0.0, pull=alpha,
                                          seed=s).node_values[0] for s in range(2000)])
    assert drawn.std() == pytest.approx(math.sqrt(expected), rel=0.07)
    at_start = _ou_segment(1.0, 0.0, T, T, alpha)
    assert at_start == pytest.approx(0.4908, rel=1e-3) and expected > 8 * at_start


def test_a_driven_switch_rate_breaks_where_the_driver_breaks(tmp_path):
    # the discrete twin of the crux: the driver turns the switch rate ON half-way along a single
    # branch (factor 0 then 5). Every switch must therefore fall in the branch's second half — and a
    # generator sampled once at the branch start would be the zero matrix, so nothing would ever
    # switch at all.
    tree = _one_branch(2.0).complete_tree
    src = _write_driver(tmp_path, [(0.0, "initial", "n0", "", "off"),
                                   (1.0, "on_branch", "n0", "off", "on")])
    switch = PerLineage(1.0).scaled_by(src, {"off": 0.0, "on": 5.0})
    times = [e.time for s in range(60)
             for e in simulate_discrete(tree, states=["x", "y"], switch=switch, start="x",
                                        seed=s).events if e.kind == "on_branch"]
    assert times, "the driver switched the rate on, so the trait must switch somewhere"
    assert min(times) > 1.0, f"a switch fired while the driver held the rate at zero: t={min(times)}"
    # the realised rate over the live half is the 5.0 asked for (60 runs × 1 time unit of exposure)
    assert len(times) / 60.0 == pytest.approx(5.0, rel=0.2)


def test_a_driver_that_never_switches_leaves_a_run_byte_identical(tmp_path):
    # the undriven guarantee, from the inside: a driver that holds one state everywhere, mapped to
    # 1.0, adds no breakpoint and no draw — so both engines return exactly the numbers they return
    # with no driver at all, value for value.
    tree = _tree(seed=5, n_extant=14)
    flat = simulate_discrete(tree, states=["one", "two"], switch=0.0, start="one", seed=3)
    table = {"one": 1.0, "two": 1.0}

    plain = simulate_continuous(tree, rate=1.0, seed=17)
    driven = simulate_continuous(tree, rate=PerLineage(1.0).scaled_by(flat, table), seed=17)
    assert driven.node_values == plain.node_values

    plain_d = simulate_discrete(tree, states=["a", "b", "c"], switch=0.4, start="a", seed=19)
    driven_d = simulate_discrete(tree, states=["a", "b", "c"], start="a", seed=19,
                                 switch=PerLineage(0.4).scaled_by(flat, table))
    assert driven_d.node_values == plain_d.node_values
    assert driven_d.history == plain_d.history
    assert [(e.time, e.lineage, e.to_state) for e in driven_d.events] == \
           [(e.time, e.lineage, e.to_state) for e in plain_d.events]


def test_a_trait_diffuses_faster_where_its_driver_says_so():
    # the model, end to end: grow a two-state trait A, then a continuous trait B whose variance-rate
    # is 6× on the lineages where A is 'fast'. Measured on the branches that lie WHOLLY in one state,
    # the per-branch step sd must be √6 times larger in 'fast' than in 'slow'.
    tree = simulate_species_tree(birth=1.0, death=0.2, n_extant=300, seed=4).complete_tree
    A = simulate_discrete(tree, states=["slow", "fast"], switch=0.5, start="slow", seed=7)
    B = simulate_continuous(tree, rate=PerLineage(1.0).scaled_by(A, {"slow": 1.0, "fast": 6.0}), seed=11)

    steps = {"slow": [], "fast": []}
    for i, node in tree.nodes.items():
        dt = node.end_time - node.birth_time
        if dt <= 0 or len(A.history[i]) != 1:      # only branches spent wholly in one state
            continue
        parent = 0.0 if node.parent is None else B.node_values[node.parent]
        steps[A.history[i][0][0]].append((B.node_values[i] - parent) / math.sqrt(dt))

    assert len(steps["slow"]) > 100 and len(steps["fast"]) > 100
    ratio = np.std(steps["fast"]) / np.std(steps["slow"])
    assert ratio == pytest.approx(math.sqrt(6.0), rel=0.1), (
        f"σ² was told to be 6× on 'fast' branches, so the step sd should be √6 = {math.sqrt(6.0):.3f}× "
        f"larger there; measured {ratio:.3f}")


def test_a_trait_switches_faster_where_its_driver_says_so():
    # the discrete twin: B's switch rate is 8× on the stretches where A is 'fast'. Counted against
    # the time B actually spends under each of A's states, the realised rates are the asked-for ones.
    tree = simulate_species_tree(birth=1.0, death=0.2, n_extant=400, seed=4).complete_tree
    A = simulate_discrete(tree, states=["slow", "fast"], switch=0.5, start="slow", seed=7)
    B = simulate_discrete(tree, states=["x", "y"], start="x", seed=13,
                          switch=PerLineage(0.2).scaled_by(A, {"slow": 1.0, "fast": 8.0}))

    exposure = {"slow": 0.0, "fast": 0.0}
    switches = {"slow": 0, "fast": 0}
    for i, node in tree.nodes.items():
        spans, t = [], node.birth_time
        for state, dur in A.history[i]:
            spans.append((t, t + dur, state))
            t += dur
        t = node.birth_time
        for k, (_state, dur) in enumerate(B.history[i]):
            for lo, hi, a in spans:
                overlap = min(t + dur, hi) - max(t, lo)
                if overlap > 0:
                    exposure[a] += overlap
            t += dur
            if k + 1 < len(B.history[i]):                       # B switched at time t
                switches[next(a for lo, hi, a in spans if lo <= t < hi)] += 1

    slow = switches["slow"] / exposure["slow"]
    fast = switches["fast"] / exposure["fast"]
    assert slow == pytest.approx(0.2, rel=0.25), f"B switched at {slow:.4f} where A is slow (asked 0.2)"
    assert fast == pytest.approx(1.6, rel=0.15), f"B switched at {fast:.4f} where A is fast (asked 1.6)"


def test_a_grown_result_and_its_written_log_drive_identically(tmp_path):
    # the two conditioned sources are the same conditioning: the in-memory result is the file's
    # shortcut, so driving by either must give the same run, node for node.
    tree = _tree(seed=6, n_extant=20)
    A = simulate_discrete(tree, states=["dry", "wet"], switch=0.6, start="dry", seed=2)
    A.write(tmp_path, outputs=["events"])
    table = {"dry": 1.0, "wet": 4.0}
    by_object = simulate_continuous(tree, rate=PerLineage(1.0).scaled_by(A, table), seed=8)
    by_file = simulate_continuous(
        tree, rate=PerLineage(1.0).scaled_by(str(tmp_path / "trait_events.tsv"), table), seed=8)
    assert by_file.node_values == by_object.node_values


def test_driven_trait_validation(tmp_path):
    tree = _tree(seed=6, n_extant=12)
    A = simulate_discrete(tree, states=["dry", "wet"], switch=0.6, start="dry", seed=2)

    # a mapping that names none of the driver's states would leave every lineage at the default
    # factor — the undriven model wearing a driven rate — so it is refused, not run
    with pytest.raises(ValueError, match="match none of"):
        simulate_continuous(tree, rate=PerLineage(1.0).scaled_by(A, {"marine": 3.0}), seed=1)
    with pytest.raises(ValueError, match="match none of"):
        simulate_discrete(tree, states=["x", "y"], seed=1,
                          switch=PerLineage(0.5).scaled_by(A, {"marine": 3.0}))

    # a driver grown on a different tree: the join key is the species node id, so a lineage the
    # driver never saw has no value to read and the run stops instead of inventing one
    other = simulate_species_tree(birth=1.0, death=0.0, n_extant=3, seed=9).complete_tree
    B = simulate_discrete(other, states=["dry", "wet"], switch=0.6, start="dry", seed=2)
    with pytest.raises(KeyError, match="SAME|node ids"):
        simulate_continuous(tree, rate=PerLineage(1.0).scaled_by(B, {"dry": 1.0, "wet": 2.0}), seed=1)

    # scaled_by is the only verb a switch rate takes; anything else on it would be read by nothing
    with pytest.raises(ValueError, match="switch rate carries changing_at"):
        simulate_discrete(tree, states=["x", "y"], switch=PerLineage(0.5).changing_at({0: 1.0, 1: 0.5}), seed=1)
    with pytest.raises(ValueError, match="per lineage"):
        simulate_discrete(tree, states=["x", "y"], seed=1,
                          switch=Global(0.5).scaled_by(A, {"dry": 1.0, "wet": 2.0}))


def test_a_driven_switch_rate_can_be_written_per_transition(tmp_path):
    # the {'from->to': rate} form takes a driven rate on any subset of its transitions: here only
    # x→y is driven, so the driver makes the trait fall into y and stay there where it is 'wet'
    tree = _one_branch(4.0).complete_tree
    src = _write_driver(tmp_path, [(0.0, "initial", "n0", "", "dry"),
                                   (2.0, "on_branch", "n0", "dry", "wet")])
    switch = {"x->y": PerLineage(1.0).scaled_by(src, {"dry": 0.0, "wet": 3.0}), "y->x": 0.05}
    times = [e.time for s in range(40)
             for e in simulate_discrete(tree, states=["x", "y"], switch=switch, start="x",
                                        seed=s).events
             if e.kind == "on_branch" and e.from_state == "x"]
    assert times and min(times) > 2.0    # x→y is off until the driver turns it on at t=2


def test_a_correlated_run_carries_the_event_log_every_other_continuous_run_carries():
    """It used to return no events at all — not even the ``initial`` row — so `trait_events.tsv` came
    out header-only, which reads as data having gone missing rather than as a diffusion being
    unreconstructable (which it is, here as for one trait).

    The rows hold the whole **vector**: a correlated jump moves every trait at once, so it is one
    event, and the table widens to ``from:<trait>`` / ``to:<trait>`` exactly as ``trait_values.tsv``
    widens. Checked against the values file: a node's post-jump value is where its branch starts."""
    tree = simulate_species_tree(birth=1.0, n_extant=6, seed=1)
    r = simulate_continuous(tree, start={"a": 0.0, "b": 0.0}, rate={"a": 1.0, "b": 0.8},
                            correlation={("a", "b"): 0.6}, at_speciation=0.5, seed=1)
    kinds = [e.kind for e in r.events]
    assert kinds[0] == "initial" and set(kinds[1:]) == {"on_speciation"}
    assert r.events[0].from_state is None
    assert r.events[0].to_state == {"a": 0.0, "b": 0.0}
    # one jump per non-root node, and every payload is the full vector
    assert len(r.events) - 1 == sum(1 for n in tree.complete_tree.nodes.values() if n.parent is not None)
    assert all(set(e.to_state) == {"a", "b"} for e in r.events)


def test_a_correlated_event_log_widens_one_column_pair_per_trait(tmp_path):
    tree = simulate_species_tree(birth=1.0, n_extant=4, seed=2)
    r = simulate_continuous(tree, start={"a": 0.0, "b": 0.0}, rate={"a": 1.0, "b": 1.0},
                            correlation={("a", "b"): 0.3}, at_speciation=0.4, seed=2)
    r.write(tmp_path, outputs=("events",))
    lines = (tmp_path / "trait_events.tsv").read_text(encoding="utf-8").splitlines()
    assert lines[0].split("\t") == ["time", "kind", "lineage", "from:a", "from:b", "to:a", "to:b"]
    assert all(len(ln.split("\t")) == 7 for ln in lines)      # the initial row keeps its empty cells


def test_a_correlated_run_without_jumps_logs_only_the_initial_row():
    tree = simulate_species_tree(birth=1.0, n_extant=5, seed=3)
    r = simulate_continuous(tree, start={"a": 0.0, "b": 0.0}, rate={"a": 1.0, "b": 1.0},
                            correlation={("a", "b"): 0.5}, seed=3)
    assert [e.kind for e in r.events] == ["initial"]          # a diffusion has no other moments


def test_the_dataset_keys_are_the_tree_tip_names():
    """The comparative vector and the tree it belongs to must join.

    `.values` used to be keyed by bare node ids while every Newick label and `trait_values.tsv` row
    says ``n5`` — so the two objects a user is meant to line up shared **no keys at all**, and
    nothing said so. This is the check that was failing silently on the Python path; the written
    files always agreed with each other."""
    import re

    tree = simulate_species_tree(birth=1.0, death=0.3, n_extant=12, seed=1)
    r = simulate_continuous(tree, rate=1.0, seed=1)

    extant_newick = tree.extant_tree.to_newick()
    tips = set(re.findall(r"[(,]([ne]\d+):", extant_newick))
    assert tips, "the extant tree should have named tips"
    assert set(r.values) == tips, (
        f"dataset keys {sorted(r.values)[:4]}… do not join the tree's tips {sorted(tips)[:4]}…")

    # and the id-keyed view still exists, with the same values behind the other key
    labels = tree.complete_tree.labels()
    assert {labels[i]: v for i, v in r.values_by_id.items()} == r.values


def test_correlation_beside_regimes_is_refused():
    """Regression. ``regimes=`` dispatches before the correlated engine and threads no correlation,
    so a correlation passed alongside it was read by nothing: the run was silently the uncorrelated
    model. Two runs differing only in a bogus correlation came back byte-identical."""
    from zombi2 import species, traits
    tree = species.simulate_species_tree(birth=1.0, death=0.3, n_extant=20, seed=1)
    reg = traits.simulate_discrete(tree, states=["r1", "r2"], switch=0.3, seed=2)
    with pytest.raises(ValueError, match="correlation= with regimes="):
        traits.simulate_continuous(tree, start=0.0, rate=1.0, pull=1.0,
                                   reverts_to={"r1": -2.0, "r2": 2.0}, regimes=reg,
                                   correlation={("x", "y"): 0.9}, seed=5)


def test_summary_of_a_correlated_run_is_per_trait():
    """Regression. A correlated run holds one value per trait at every node, and ``summary()``
    assumed one number: the continuous case raised ``TypeError: float() argument must be … not
    'dict'`` (so ``write(..., "summary")`` failed outright), and the threshold case counted whole
    dicts as states, giving keys like "{'a': 'b', 'b': 'b'}"."""
    from zombi2 import species, traits
    tree = species.simulate_species_tree(birth=1.0, death=0.3, n_extant=20, seed=1)

    c = traits.simulate_continuous(tree, start={"x": 0.0, "y": 0.0}, rate={"x": 1.0, "y": 1.0},
                                   correlation={("x", "y"): 0.7}, seed=6)
    s = c.summary()
    assert s["traits"] == ["x", "y"]
    assert set(s["values"]) == {"x", "y"} and "mean" in s["values"]["x"]
    assert set(s["value_at_root_node"]) == {"x", "y"}

    d = traits.simulate_discrete(tree, states=["a", "b"], liability={"p": 1.0, "q": 1.0},
                                 threshold=0.0, correlation={("p", "q"): 0.6}, seed=6)
    s2 = d.summary()
    assert s2["traits"] == ["p", "q"]
    assert set(s2["states"]) == {"p", "q"}
    assert set(s2["states"]["p"]) <= {"a", "b"}          # states, not stringified dicts


def test_summary_of_a_single_trait_run_is_unchanged():
    """The per-trait branch must not capture the ordinary one-trait case."""
    from zombi2 import species, traits
    tree = species.simulate_species_tree(birth=1.0, death=0.3, n_extant=20, seed=1)
    c = traits.simulate_continuous(tree, start=0.0, rate=1.0, seed=6).summary()
    assert "traits" not in c and "mean" in c["values"]
    d = traits.simulate_discrete(tree, states=["a", "b"], switch=0.3, seed=6).summary()
    assert "traits" not in d and set(d["states"]) <= {"a", "b"}



def test_a_discrete_trait_steps_at_a_rate_that_changes_on_its_own_clock():
    """`traits.discrete` built its stretch horizon from the drivers alone and never asked
    `Rate.next_change`, so a switch rate that changes with time was read once at the start of a
    stretch and held for the rest of the branch — the run silently not the model asked for, which
    is exactly what the modifier gate exists to prevent. Every other engine steps to it.

    Only a modifier of your own can reach this: the engine's own declaration takes `Driven`
    alone, and a driver's switches were already stepped to. Appendix A publishes `traits.discrete`
    as an engine an `implemented_for` may name, so this is a path the manual invites."""
    import math

    from zombi2 import species, traits
    from zombi2.params.evaluate import Modifier

    class StopsAtOne(Modifier):
        """Full rate until time 1, nothing after — and it says so, as a modifier must."""

        implemented_for = ("traits.discrete",)

        def factor(self, *, time=0.0, **_):
            return 1.0 if time < 1.0 else 0.0

        def next_change(self, time):
            return 1.0 if time < 1.0 else math.inf

    tree = species.simulate_species_tree(birth=1.0, death=0.0, total_time=3.0, seed=5).complete_tree
    result = traits.simulate_discrete(tree, states=["a", "b"], seed=3,
                                      switch=PerLineage(3.0)._and(StopsAtOne()))

    late = [c.time for c in result.events if c.time > 1.0 + 1e-9]
    assert not late, f"switched after the rate went to zero: {late[:5]}"
    assert result.events, "the test proves nothing if nothing switched before time 1 either"


def test_the_continuous_process_spec_takes_the_package_attribute():
    """`traits.discrete` is a process spec for a joint run, and `traits.continuous` is what a reader
    reaches for next: the diffusing twin, for QuaSSE. The name has to be the function rather than the
    *module* of the same name — that reach used to fail with "'module' object is not callable", which
    names neither traits nor speciation nor what to write instead."""
    from zombi2 import traits
    from zombi2.traits import ContinuousTrait

    assert callable(traits.continuous), "the module is shadowing the spec again"
    spec = traits.continuous(rate=1.0, start=2.0, name="size")
    assert isinstance(spec, ContinuousTrait)
    assert (spec.rate, spec.start, spec.name) == (1.0, 2.0, "size")


def test_shadowing_the_module_leaves_both_ways_in_working():
    """The refusal takes the package attribute; the module keeps its import path, and the standalone
    continuous engine is untouched."""
    from zombi2 import species, traits
    from zombi2.traits.continuous import simulate_continuous as by_module_path

    tree = species.simulate_species_tree(birth=1.0, n_extant=6, seed=1)
    assert by_module_path is traits.simulate_continuous
    assert traits.simulate_continuous(tree, rate=1.0, seed=1).kind == "continuous"
