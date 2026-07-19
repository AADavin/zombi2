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

from zombi2.rates import modifiers as mod
from zombi2.rates import scope
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
    extinct = {n.id for n in sp.complete_tree.extinct()}
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
    extant = {n.id for n in sp.complete_tree.extant()}
    assert set(r.values) == extant
    assert all(r.values[i] == r.node_values[i] for i in extant)


def test_continuous_events_are_empty():
    sp = _tree(seed=8)
    r = simulate_continuous(sp, rate=0.5, seed=1)
    assert r.events == []           # a continuous trait has no instantaneous events
    assert r.kind == "continuous"
    assert r.history is None


def test_write_values_tsv(tmp_path):
    sp = _tree(seed=8)
    r = simulate_continuous(sp, rate=0.5, seed=1)
    r.write(tmp_path, outputs=["values"])
    text = (tmp_path / "trait_values.tsv").read_text()
    lines = text.splitlines()
    assert lines[0] == "node\ttrait"
    assert len(lines) - 1 == len(r.values)                    # one row per extant tip
    ids = {int(line.split("\t")[0][1:]) for line in lines[1:]}  # strip the "n" prefix
    assert ids == set(r.values)


def test_write_rejects_unknown_output(tmp_path):
    sp = _tree(seed=8)
    r = simulate_continuous(sp, rate=0.5, seed=1)
    with pytest.raises(ValueError, match="unknown write outputs"):
        r.write(tmp_path, outputs=["history"])


# --- input validation (what this slice deliberately does not wire) --------------

def test_rejects_a_non_time_modifier():
    # Time (early burst) is wired as of slice 2; other modifiers are still later slices
    sp = _tree(seed=1)
    with pytest.raises(ValueError, match="Diversity"):
        simulate_continuous(sp, rate=1.0 * mod.Diversity(cap=100), seed=1)


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
    tips = sorted(n.id for n in tree.extant())
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
    tips = sorted(n.id for n in sp.complete_tree.extant())
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
    tips = sorted(n.id for n in tree.extant())
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


def test_ou_with_time_is_deferred():
    sp = _tree(seed=1)
    with pytest.raises(ValueError, match="not wired yet"):
        simulate_continuous(sp, rate=1.0 * mod.Time({0: 1.0, 3: 0.2}),
                            reverts_to=2.0, pull=0.5, seed=1)


# --- early burst (a Time skyline on rate): the exact ∫σ²(t)dt over each branch --

def test_eb_constant_schedule_equals_bm():
    # a single-step schedule with factor 1.0 everywhere is σ² constant → byte-identical to bare BM
    # (same one draw per branch), which pins the integral's constant-rate special case.
    sp = _tree(seed=2)
    a = simulate_continuous(sp, rate=3.0, seed=5)
    b = simulate_continuous(sp, rate=3.0 * mod.Time({0.0: 1.0}), seed=5)
    assert a.node_values == b.node_values


def test_eb_tip_variance_and_covariance_match_the_integral():
    # early burst: σ²(t) drops from `base` to `base·c` at time τ. The per-branch variance is the
    # exact integral, so a tip at depth T has variance base·(τ + c·(T−τ)) and a pair covariance
    # base·∫_0^s over the shared path to their MRCA split s. This pins the skyline integral.
    sp = _tree(seed=11, n_extant=6, death=0.0)
    tree = sp.complete_tree
    tips = sorted(n.id for n in tree.extant())
    T = tree.nodes[tips[0]].end_time
    base, c = 2.0, 0.25
    tau = 0.4 * T                                  # guaranteed inside (0, T) so the branch crosses it
    sched = base * mod.Time({0.0: 1.0, tau: c})

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
    sched = 1.5 * mod.Time({0.0: 1.0, 2.0: 0.3})
    a = simulate_continuous(sp, rate=sched, seed=4)
    b = simulate_continuous(sp, rate=sched, seed=4)
    assert a.node_values == b.node_values


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
    assert r.events == []                                  # and no transitions
    assert all(len(segs) == 1 for segs in r.history.values())


def test_discrete_result_shape():
    sp = _tree(seed=8)
    r = simulate_discrete(sp, states=["marine", "terrestrial"], switch=0.3, start="marine", seed=1)
    assert r.kind == "discrete"
    assert set(r.values) == {n.id for n in sp.complete_tree.extant()}
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
    # one event per jump between consecutive segments, across all branches
    n_jumps = sum(len(segs) - 1 for segs in r.history.values())
    assert len(r.events) == n_jumps and n_jumps > 0
    assert all(isinstance(e, Change) and e.from_state != e.to_state for e in r.events)
    assert r.events == sorted(r.events, key=lambda e: e.time)      # time-ordered
    for e in r.events:                                             # each change sits on its branch
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
    r.write(tmp_path, outputs=["values", "changes"])
    vals = (tmp_path / "trait_values.tsv").read_text().splitlines()
    assert vals[0] == "node\ttrait" and set(line.split("\t")[1] for line in vals[1:]) <= {"lo", "hi"}
    changes = (tmp_path / "trait_changes.tsv").read_text().splitlines()
    assert changes[0] == "time\tlineage\tfrom\tto" and len(changes) - 1 == len(r.events)


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


def test_threshold_is_deferred():
    sp = _tree(seed=1)
    with pytest.raises(ValueError, match="later slice"):
        simulate_discrete(sp, states=["absent", "present"], liability=1.0, threshold=0.0, seed=1)
