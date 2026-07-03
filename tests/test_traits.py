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


# --------------------------------------------------------------------------- Ornstein-Uhlenbeck
def _single_branch(t):
    from zombi2.tree import Tree, TreeNode
    root = TreeNode("r", 0.0)
    tip = TreeNode("a", t)
    root.add_child(tip)
    return Tree(root, t), tip


def test_ou_defaults_and_reproducible():
    tree = _fixed_tree()
    m = z.OrnsteinUhlenbeck(0.5, 1.0, theta=2.0)
    assert m.x0 == 2.0  # default root is the optimum
    a = z.simulate_traits(tree, m, seed=4).values
    b = z.simulate_traits(tree, m, seed=4).values
    assert {k.name: v for k, v in a.items()} == {k.name: v for k, v in b.items()}


def test_ou_validation():
    with pytest.raises(ValueError):
        z.OrnsteinUhlenbeck(0.5, 0.0, 1.0)   # alpha must be > 0
    with pytest.raises(ValueError):
        z.OrnsteinUhlenbeck(-1.0, 1.0, 1.0)  # sigma2 must be >= 0


def test_ou_transition_moments_match_theory():
    """Single-branch endpoint: mean reverts toward theta; variance matches the OU formula."""
    t = 0.8
    tree, tip = _single_branch(t)
    sigma2, alpha, theta, x0 = 0.6, 1.5, 2.0, -1.0
    m = z.OrnsteinUhlenbeck(sigma2, alpha, theta, x0=x0)
    rng = np.random.default_rng(0)
    vals = np.array([z.simulate_traits(tree, m, rng=rng).node_values[tip] for _ in range(20000)])
    e = np.exp(-alpha * t)
    assert abs(vals.mean() - (theta + (x0 - theta) * e)) < 0.03
    assert abs(vals.var() - sigma2 / (2 * alpha) * (1 - e * e)) < 0.03


def test_ou_pulls_toward_optimum():
    """Deep tips started far from theta concentrate near theta (stationary mean = theta)."""
    tree = _fixed_tree(n_tips=8, seed=1)
    m = z.OrnsteinUhlenbeck(0.3, 3.0, theta=5.0, x0=-5.0)
    tips = np.array(list(z.simulate_traits(tree, m, seed=2).values.values()))
    assert abs(tips.mean() - 5.0) < 1.0  # far from the x0=-5 start


# --------------------------------------------------------------------------- multivariate Brownian
def test_mvbm_reproducible_and_vector_valued():
    tree = _fixed_tree(n_tips=4)
    m = z.MultivariateBrownian([[1.0, 0.3], [0.3, 0.8]])
    res = z.simulate_traits(tree, m, seed=1)
    leaf = tree.extant_leaves()[0]
    assert np.asarray(res.values[leaf]).shape == (2,)
    again = z.simulate_traits(tree, m, seed=1).values
    assert np.allclose(res.values[leaf], again[leaf])


def test_mvbm_validation():
    with pytest.raises(ValueError):
        z.MultivariateBrownian([[1.0, 0.0]])                  # not square
    with pytest.raises(ValueError):
        z.MultivariateBrownian([[1.0, 2.0], [2.0, 1.0]])      # not PSD
    with pytest.raises(ValueError):
        z.MultivariateBrownian([[1.0, 0.0], [0.0, 1.0]], x0=[0, 0, 0])  # wrong length


def test_mvbm_tip_covariance_matches_R_times_C():
    """Per-tip covariance is R·depth; cross-tip same-dim covariance is R[a,a]·MRCA-time."""
    tree = _fixed_tree(n_tips=4, seed=13)
    leaves = tree.extant_leaves()
    R = np.array([[1.0, 0.6], [0.6, 0.8]])
    m = z.MultivariateBrownian(R)
    rng = np.random.default_rng(0)
    reps = 6000
    data = np.empty((reps, len(leaves), 2))
    for r in range(reps):
        vals = z.simulate_traits(tree, m, rng=rng).values
        for i, leaf in enumerate(leaves):
            data[r, i] = vals[leaf]

    depth = leaves[0].time
    assert np.allclose(np.cov(data[:, 0], rowvar=False), R * depth, atol=0.12)
    mrca = _mrca_time(leaves[0], leaves[1])
    cross = np.cov(data[:, 0, 0], data[:, 1, 0])[0, 1]
    assert abs(cross - R[0, 0] * mrca) < 0.12


# --------------------------------------------------------------------------- multivariate OU
def test_mvou_isotropic_stationary_is_R_over_2alpha():
    R = np.array([[1.0, 0.5], [0.5, 2.0]])
    alpha = 1.3
    m = z.MultivariateOU(R, alpha, theta=[0.0, 0.0])
    assert np.allclose(m.V, R / (2 * alpha))


def test_mvou_lyapunov_holds_for_full_A():
    R = np.array([[1.0, 0.3], [0.3, 0.7]])
    A = np.array([[1.2, 0.2], [0.0, 0.9]])
    m = z.MultivariateOU(R, A, theta=[0.0, 0.0])
    assert np.allclose(m.A @ m.V + m.V @ m.A.T, R, atol=1e-9)


def test_mvou_reduces_to_bm_as_alpha_small():
    R = np.array([[1.0, 0.2], [0.2, 0.9]])
    t = 0.5
    m = z.MultivariateOU(R, 1e-4, theta=[0.0, 0.0])
    E = m._E(t)
    cov = m.V - E @ m.V @ E.T
    assert np.allclose(cov, R * t, atol=1e-2)


def test_mvou_single_branch_moments():
    t = 0.7
    tree, tip = _single_branch(t)
    R = np.array([[1.0, 0.4], [0.4, 1.2]])
    theta = np.array([0.5, -0.5])
    x0 = np.array([3.0, -3.0])
    m = z.MultivariateOU(R, [1.0, 1.6], theta=theta, x0=x0)
    rng = np.random.default_rng(1)
    reps = 12000
    data = np.empty((reps, 2))
    for r in range(reps):
        data[r] = z.simulate_traits(tree, m, rng=rng).node_values[tip]
    E = m._E(t)
    assert np.allclose(data.mean(0), theta + E @ (x0 - theta), atol=0.05)
    assert np.allclose(np.cov(data, rowvar=False), m.V - E @ m.V @ E.T, atol=0.08)


def test_mvou_validation():
    with pytest.raises(ValueError):
        z.MultivariateOU([[1, 0], [0, 1]], -1.0, theta=[0, 0])   # not mean-reverting
    with pytest.raises(ValueError):
        z.MultivariateOU([[1, 0], [0, 1]], 1.0, theta=[0, 0, 0])  # theta wrong length


def test_vector_trait_formatting():
    tree = _fixed_tree(n_tips=4)
    res = z.simulate_traits(tree, z.MultivariateBrownian([[1.0, 0.0], [0.0, 1.0]]), seed=1)
    assert "{" in res.to_tsv()
    assert "[&trait={" in res.to_newick()
