"""Tests for trait evolution (``zombi2.traits``).

The core models are validated against their analytic laws:

* Brownian motion — tip mean ``x0 + trend·depth`` and tip covariance ``sigma2 · C`` (``C`` =
  shared-path-length / MRCA-time matrix);
* Mk — the closed-form equal-rates transition probabilities ``P_ii(t) = 1/k + (1-1/k)e^{-kqt}``,
  a uniform stationary distribution for symmetric ``Q``, and the exactness of the simulated
  stochastic map (segment durations tile the branch; end state == node value).
"""

import numpy as np
import pytest

import zombi2 as z


# --------------------------------------------------------------------------- helpers
def _fixed_tree(n_tips=6, seed=7):
    return z.simulate_species_tree(z.Yule(1.0), n_tips=n_tips, age=2.0, seed=seed)


def _ancestors(node):
    out = []
    while node is not None:
        out.append(node)
        node = node.parent
    return out


def _mrca_time(a, b):
    aset = {id(n): n for n in _ancestors(a)}
    for n in _ancestors(b):  # first (deepest) shared ancestor walking up from b
        if id(n) in aset:
            return n.time
    raise AssertionError("no common ancestor")


# --------------------------------------------------------------------------- Brownian motion
def test_bm_reproducible():
    tree = _fixed_tree()
    a = z.simulate_traits(tree, z.BrownianMotion(0.5), seed=3).values
    b = z.simulate_traits(tree, z.BrownianMotion(0.5), seed=3).values
    assert {k.name: v for k, v in a.items()} == {k.name: v for k, v in b.items()}


def test_bm_zero_sigma_is_constant():
    tree = _fixed_tree()
    res = z.simulate_traits(tree, z.BrownianMotion(0.0, x0=3.0), seed=1)
    assert all(v == 3.0 for v in res.node_values.values())


def test_bm_root_state_override():
    tree = _fixed_tree()
    res = z.simulate_traits(tree, z.BrownianMotion(0.0, x0=1.0), root_state=9.0, seed=1)
    assert res.node_values[tree.root] == 9.0
    assert all(v == 9.0 for v in res.values.values())  # sigma2=0 keeps it flat


def test_bm_continuous_has_no_history():
    tree = _fixed_tree()
    res = z.simulate_traits(tree, z.BrownianMotion(0.5), seed=1)
    assert res.kind == "continuous"
    assert res.history is None


def test_bm_tip_moments_match_theory():
    """Empirical tip mean and covariance match ``x0 + trend·depth`` and ``sigma2·C``."""
    tree = _fixed_tree(n_tips=5, seed=11)
    leaves = tree.extant_leaves()
    sigma2, x0, trend = 0.8, 1.0, 0.3

    rng = np.random.default_rng(0)
    reps = 6000
    data = np.empty((reps, len(leaves)))
    model = z.BrownianMotion(sigma2, x0=x0, trend=trend)
    for r in range(reps):
        vals = z.simulate_traits(tree, model, rng=rng).values
        data[r] = [vals[leaf] for leaf in leaves]

    depths = np.array([leaf.time for leaf in leaves])
    exp_mean = x0 + trend * depths
    assert np.allclose(data.mean(axis=0), exp_mean, atol=0.06)

    C = np.array([[_mrca_time(a, b) for b in leaves] for a in leaves])
    emp_cov = np.cov(data, rowvar=False)
    assert np.allclose(emp_cov, sigma2 * C, atol=0.12)


# --------------------------------------------------------------------------- Mk: construction
def test_mk_diagonal_recomputed():
    m = z.Mk([[0, 1, 2], [3, 0, 1], [1, 1, 0]])
    assert np.allclose(m.Q.sum(axis=1), 0.0)
    assert m.Q[0, 0] == -3 and m.Q[1, 1] == -4


def test_mk_equal_rates_structure():
    m = z.Mk.equal_rates(4, 0.5)
    off = m.Q.copy()
    np.fill_diagonal(off, np.nan)
    assert np.nanmin(off) == 0.5 and np.nanmax(off) == 0.5
    assert np.allclose(np.diag(m.Q), -1.5)  # (k-1)*rate


def test_mk_symmetric_is_symmetrized():
    m = z.Mk.symmetric([[0, 2, 0], [1, 0, 3], [0, 3, 0]])
    off = m.Q.copy()
    np.fill_diagonal(off, 0.0)
    assert np.allclose(off, off.T)
    assert off[0, 1] == 1.5  # (2+1)/2


def test_mk_invalid_inputs():
    with pytest.raises(ValueError):
        z.Mk([[0, 1], [1, 0], [0, 0]])          # not square
    with pytest.raises(ValueError):
        z.Mk([[0, -1], [1, 0]])                 # negative off-diagonal
    with pytest.raises(ValueError):
        z.Mk.equal_rates(1)                     # < 2 states
    with pytest.raises(ValueError):
        z.Mk([[0, 1], [1, 0]], root="bogus")


# --------------------------------------------------------------------------- Mk: analytic laws
def test_mk_equal_rates_transition_closed_form():
    """P(t) for ER matches P_ii = 1/k + (1-1/k)e^{-kqt}, P_ij = (1/k)(1-e^{-kqt})."""
    k, q, t = 4, 0.7, 1.3
    m = z.Mk.equal_rates(k, q)
    P = m.transition_matrix(t)
    e = np.exp(-k * q * t)
    diag = 1.0 / k + (1.0 - 1.0 / k) * e
    off = (1.0 / k) * (1.0 - e)
    expected = np.full((k, k), off)
    np.fill_diagonal(expected, diag)
    assert np.allclose(P, expected, atol=1e-9)
    assert np.allclose(P.sum(axis=1), 1.0)


def test_mk_transition_semigroup():
    m = z.Mk([[0, 1, 0.5], [0.2, 0, 0.8], [0.3, 0.3, 0]])
    P1, P2 = m.transition_matrix(0.4), m.transition_matrix(0.9)
    assert np.allclose(m.transition_matrix(1.3), P1 @ P2, atol=1e-9)
    assert np.allclose(m.transition_matrix(0.0), np.eye(3), atol=1e-12)


def test_mk_symmetric_stationary_is_uniform():
    m = z.Mk.symmetric([[0, 2, 1], [2, 0, 3], [1, 3, 0]])
    assert np.allclose(m.stationary_distribution(), 1.0 / 3, atol=1e-9)


def test_mk_simulation_matches_transition_matrix():
    """Empirical end-state distribution over a single branch matches P(t)[start]."""
    from zombi2.tree import Tree, TreeNode
    t = 1.1
    root = TreeNode("r", 0.0)
    tip = TreeNode("a", t)
    root.add_child(tip)
    tree = Tree(root, t)

    k, start = 4, 0
    m = z.Mk([[0, 1.0, 0.2, 0.5], [0.3, 0, 0.7, 0.1],
              [0.4, 0.2, 0, 0.6], [0.1, 0.9, 0.3, 0]], root=start)
    rng = np.random.default_rng(1)
    counts = np.zeros(k)
    reps = 20000
    for _ in range(reps):
        v = z.simulate_traits(tree, m, rng=rng).node_values[tip]
        counts[v] += 1
    emp = counts / reps
    assert np.allclose(emp, m.transition_matrix(t)[start], atol=0.02)


# --------------------------------------------------------------------------- Mk: stochastic map
def test_mk_history_tiles_branches_and_matches_values():
    tree = _fixed_tree()
    res = z.simulate_traits(tree, z.Mk.equal_rates(3, 0.8), seed=5)
    assert res.kind == "discrete"
    for node in tree.nodes():
        if node.parent is None:
            continue
        segs = res.history[node]
        # durations tile the branch
        assert abs(sum(d for _, d in segs) - node.branch_length()) < 1e-9
        # end state == node value; start state == parent's end value (continuity)
        assert segs[-1][0] == res.node_values[node]
        assert segs[0][0] == res.node_values[node.parent]


def test_mk_changes_are_real_transitions():
    tree = _fixed_tree(n_tips=10, seed=2)
    res = z.simulate_traits(tree, z.Mk.equal_rates(3, 1.5), seed=9)
    for node, time, frm, to in res.changes():
        assert frm != to
        assert node.parent.time <= time <= node.time


def test_mk_root_policies():
    tree = _fixed_tree()
    # fixed index
    assert z.simulate_traits(tree, z.Mk.equal_rates(3, 0.5, root=2), seed=1).node_values[tree.root] == 2
    # probability vector pinning state 1
    r = z.simulate_traits(tree, z.Mk.equal_rates(3, 0.5, root=[0, 1, 0]), seed=1)
    assert r.node_values[tree.root] == 1
    # stationary sampling stays a valid state
    r2 = z.simulate_traits(tree, z.Mk.symmetric([[0, 1, 1], [1, 0, 1], [1, 1, 0]],
                                                 root="stationary"), seed=1)
    assert r2.node_values[tree.root] in (0, 1, 2)


def test_mk_state_labels():
    tree = _fixed_tree()
    m = z.Mk.equal_rates(3, 0.6, states=["marine", "brackish", "freshwater"])
    res = z.simulate_traits(tree, m, seed=4)
    for leaf in tree.extant_leaves():
        assert res.label(res.node_values[leaf]) in ("marine", "brackish", "freshwater")


# --------------------------------------------------------------------------- driver / result / IO
def test_values_are_extant_leaves_only():
    tree = z.simulate_species_tree(z.BirthDeath(1.0, 0.6), age=4.0,
                                   direction="forward", seed=3)  # complete tree w/ extinct leaves
    res = z.simulate_traits(tree, z.BrownianMotion(0.4), seed=1)
    extant = set(tree.extant_leaves())
    assert set(res.values) == extant
    assert len(res.leaf_values()) >= len(res.values)  # extinct tips included there


def test_works_on_read_newick_tree():
    tree = z.read_newick("((a:1,b:1)x:1,(c:1.5,d:0.5)y:0.5)root;")
    res = z.simulate_traits(tree, z.Mk.equal_rates(2, 0.5), seed=1)
    assert {leaf.name for leaf in tree.extant_leaves()} == {"a", "b", "c", "d"}
    assert set(res.values) == set(tree.extant_leaves())


def test_to_tsv_and_newick():
    tree = _fixed_tree()
    res = z.simulate_traits(tree, z.Mk.equal_rates(3, 0.5, states=["a", "b", "c"]), seed=1)
    tsv = res.to_tsv()
    assert tsv.splitlines()[0] == "node\ttrait"
    assert len(tsv.splitlines()) == 1 + len(tree.extant_leaves())
    assert res.to_tsv(nodes="all").count("\n") == 1 + len(tree.nodes())
    nwk = res.to_newick()
    assert nwk.endswith(";") and "[&trait=" in nwk


def test_ancestral_states_cover_internal_nodes():
    tree = _fixed_tree()
    res = z.simulate_traits(tree, z.BrownianMotion(0.5), seed=1)
    assert set(res.ancestral_states()) == set(tree.internal_nodes())
