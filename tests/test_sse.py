"""Tests for state-dependent speciation/extinction (``zombi2.sse``: BiSSE / MuSSE).

The stochastic joint tree+trait process is checked against known reductions and the expected
direction of its signal:

* state-independent rates with ``mu=0`` reduce to a Yule tree (mean extant ``2·e^{λ·age}``);
* symmetric dynamics give ~50/50 tip states;
* a higher speciation rate in one state biases the standing tips toward it, and a higher
  extinction rate biases them away — the core BiSSE effect;
* the tree is well-formed and the per-branch character map tiles every branch.
"""

import numpy as np
import pytest

import zombi2 as z


def _tip_state_fraction(model, *, age, reps, seed=0, state=1):
    """Aggregate fraction of extant tips in ``state`` over many independent trees."""
    rng = np.random.default_rng(seed)
    n_state = n_total = 0
    for _ in range(reps):
        res = z.simulate_sse(model, age=age, rng=rng)
        vals = list(res.values.values())
        n_total += len(vals)
        n_state += sum(1 for v in vals if v == state)
    return n_state / n_total


# --------------------------------------------------------------------------- basics / validation
def test_sse_returns_traitresult_with_tree():
    res = z.simulate_sse(z.BiSSE(1, 1, 0.1, 0.1, 0.2, 0.2), age=3.0, seed=1)
    assert isinstance(res, z.TraitResult)
    assert isinstance(res.tree, z.Tree)
    assert res.kind == "discrete"
    assert set(res.values) == set(res.tree.extant_leaves())
    assert all(v in (0, 1) for v in res.values.values())


def test_sse_reproducible():
    m = z.BiSSE(1, 2, 0.2, 0.2, 0.15, 0.15)
    a = z.simulate_sse(m, age=3.0, seed=7)
    b = z.simulate_sse(m, age=3.0, seed=7)
    assert a.to_newick() == b.to_newick()


def test_sse_validation():
    with pytest.raises(ValueError):
        z.BiSSE(-1, 1, 0, 0, 0.1, 0.1)                 # negative speciation
    with pytest.raises(ValueError):
        z.BiSSE(1, 1, 0, 0, -0.1, 0.1)                 # negative transition
    with pytest.raises(ValueError):
        z.simulate_sse(z.BiSSE(1, 1, 0, 0, 0.1, 0.1))  # neither age nor n_tips
    with pytest.raises(ValueError):
        z.simulate_sse(z.BiSSE(1, 1, 0, 0, 0.1, 0.1), age=3, n_tips=5)  # both
    with pytest.raises(ValueError):
        z.simulate_sse(z.BiSSE(1, 1, 0, 0, 0.1, 0.1), age=3, root_state=5)  # bad state
    with pytest.raises(TypeError):
        z.simulate_sse(z.BirthDeath(1, 0.3), age=3)    # not an SSE model


def test_sse_no_transitions_requires_root_state():
    with pytest.raises(ValueError):
        z.simulate_sse(z.BiSSE(1, 1, 0, 0, 0, 0), age=3)  # Q=0 -> stationary undefined
    # ... but works with an explicit root state
    res = z.simulate_sse(z.BiSSE(1, 1, 0, 0, 0, 0), age=3, root_state=1, seed=1)
    assert all(v == 1 for v in res.node_values.values())  # no transitions -> all state 1


# --------------------------------------------------------------------------- tree structure
def test_sse_tree_is_wellformed():
    res = z.simulate_sse(z.BiSSE(1.5, 1.5, 0.4, 0.4, 0.2, 0.2), age=3.0, seed=3)
    tree = res.tree
    assert tree.root.time == 0.0 and len(tree.root.children) == 2
    for node in tree.internal_nodes():
        assert len(node.children) == 2                 # bifurcating
    for leaf in tree.leaves():
        if leaf.is_extant:
            assert abs(leaf.time - tree.total_age) < 1e-9
        else:
            assert leaf.time < tree.total_age + 1e-9   # extinct before the present


def test_sse_no_extinction_has_no_extinct_leaves():
    res = z.simulate_sse(z.BiSSE(1, 1, 0, 0, 0.3, 0.3), age=2.5, seed=2)
    assert all(leaf.is_extant for leaf in res.tree.leaves())


def test_sse_n_tips_mode_is_exact():
    res = z.simulate_sse(z.BiSSE(1, 1, 0.2, 0.2, 0.2, 0.2), n_tips=12, seed=5)
    assert len(res.tree.extant_leaves()) == 12


# --------------------------------------------------------------------------- stochastic map
def test_sse_history_tiles_branches():
    res = z.simulate_sse(z.BiSSE(1, 1, 0.2, 0.2, 0.6, 0.6), age=3.0, seed=4)
    tree = res.tree
    for node in tree.nodes():
        if node.parent is None:
            continue
        segs = res.history[node]
        assert abs(sum(d for _, d in segs) - node.branch_length()) < 1e-9
        assert segs[-1][0] == res.node_values[node]
        assert segs[0][0] == res.node_values[node.parent]
    for node, time, frm, to in res.changes():
        assert frm != to


# --------------------------------------------------------------------------- reductions & signal
def test_sse_reduces_to_yule_mean_count():
    """State-independent, mu=0 -> a Yule tree with mean extant count 2·e^{λ·age}."""
    age, lam = 1.5, 1.0
    rng = np.random.default_rng(0)
    m = z.BiSSE(lam, lam, 0.0, 0.0, 0.3, 0.3)
    counts = [len(z.simulate_sse(m, age=age, rng=rng).tree.extant_leaves()) for _ in range(600)]
    assert abs(np.mean(counts) - 2 * np.exp(lam * age)) < 0.9


def test_sse_symmetric_is_balanced():
    frac1 = _tip_state_fraction(z.BiSSE(1, 1, 0.0, 0.0, 0.5, 0.5), age=3.0, reps=300)
    assert abs(frac1 - 0.5) < 0.05


def test_sse_faster_speciation_biases_tips():
    """State 1 speciates 3x faster -> standing tips are strongly biased toward state 1."""
    frac1 = _tip_state_fraction(z.BiSSE(1, 3, 0.3, 0.3, 0.1, 0.1), age=1.8, reps=150)
    assert frac1 > 0.7


def test_sse_faster_extinction_biases_tips_away():
    """State 1 goes extinct much faster -> standing tips are biased away from state 1."""
    frac1 = _tip_state_fraction(z.BiSSE(2, 2, 0.1, 1.5, 0.2, 0.2), age=2.5, reps=300)
    assert frac1 < 0.45


# --------------------------------------------------------------------------- MuSSE (general k)
def test_musse_three_states_runs():
    Q = [[0, 0.2, 0.2], [0.2, 0, 0.2], [0.2, 0.2, 0]]
    m = z.MuSSE(birth=[1, 1, 1], death=[0.1, 0.1, 0.1], Q=Q, states=["a", "b", "c"])
    res = z.simulate_sse(m, age=3.0, seed=1)
    seen = {res.label(v) for v in res.values.values()}
    assert seen <= {"a", "b", "c"} and len(res.values) >= 2


def test_musse_stationary_distribution():
    Q = [[0, 1, 0], [1, 0, 1], [0, 1, 0]]
    m = z.MuSSE(birth=[1, 1, 1], death=[0, 0, 0], Q=Q)
    pi = m.stationary_distribution()
    assert np.allclose(pi.sum(), 1.0) and np.all(pi >= 0)
    assert np.allclose(m.Q.T @ pi, 0.0, atol=1e-9)     # πQ = 0
