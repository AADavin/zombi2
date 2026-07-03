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


# --------------------------------------------------------------------------- correlated binary (Pagel 1994)
def _pooled_corr(tree, model, reps, seed=0):
    """corr(X, Y) over extant tips pooled across independent replicates on a fixed tree."""
    rng = np.random.default_rng(seed)
    xs, ys = [], []
    for _ in range(reps):
        res = z.simulate_traits(tree, model, rng=rng)
        for v in res.values.values():
            x, y = res.label(v)
            xs.append(x)
            ys.append(y)
    return np.corrcoef(xs, ys)[0, 1]


def test_correlated_binary_Q_structure():
    m = z.CorrelatedBinary(x_gain_y0=1, x_gain_y1=2, x_loss_y0=3, x_loss_y1=4,
                           y_gain_x0=5, y_gain_x1=6, y_loss_x0=7, y_loss_x1=8)
    Q = m.Q
    # no simultaneous double changes: (0,0)<->(1,1) and (0,1)<->(1,0)
    assert Q[0, 3] == 0 and Q[3, 0] == 0 and Q[1, 2] == 0 and Q[2, 1] == 0
    assert np.allclose(Q.sum(axis=1), 0.0)
    assert Q[0, 1] == 5 and Q[0, 2] == 1        # y_gain_x0, x_gain_y0
    assert Q[3, 1] == 4 and Q[3, 2] == 8        # x_loss_y1, y_loss_x1
    assert m.states == [(0, 0), (0, 1), (1, 0), (1, 1)]


def test_correlated_binary_independent_factorizes():
    m = z.CorrelatedBinary.independent(x_gain=0.7, x_loss=0.4, y_gain=0.9, y_loss=0.5)
    assert m.Q[0, 2] == m.Q[1, 3]   # x_gain independent of Y
    assert m.Q[2, 0] == m.Q[3, 1]   # x_loss independent of Y
    assert m.Q[0, 1] == m.Q[2, 3]   # y_gain independent of X


def test_correlated_binary_validation():
    with pytest.raises(ValueError):
        z.CorrelatedBinary(-1, 1, 1, 1, 1, 1, 1, 1)


def test_correlated_binary_root_tuple_and_labels():
    tree = _fixed_tree()
    m = z.CorrelatedBinary.independent(0, 0, 0, 0, root=(1, 0))  # frozen at (1,0)
    res = z.simulate_traits(tree, m, seed=1)
    assert res.node_values[tree.root] == 2                       # index of (1,0)
    assert all(lab == (1, 0) for lab in res.labeled_values().values())


def test_correlated_binary_independent_has_no_tip_association():
    tree = z.simulate_species_tree(z.Yule(1.0), n_tips=30, age=3.0, seed=7)
    m = z.CorrelatedBinary.independent(x_gain=0.5, x_loss=0.5, y_gain=0.5, y_loss=0.5)
    assert abs(_pooled_corr(tree, m, reps=250, seed=1)) < 0.1


def test_correlated_binary_dependent_induces_association():
    """Y tracks X (gains fast when X=1, lost fast when X=0) -> positive tip association."""
    tree = z.simulate_species_tree(z.Yule(1.0), n_tips=30, age=3.0, seed=7)
    m = z.CorrelatedBinary(x_gain_y0=0.5, x_gain_y1=0.5, x_loss_y0=0.5, x_loss_y1=0.5,
                           y_gain_x0=0.05, y_gain_x1=2.0, y_loss_x0=2.0, y_loss_x1=0.05)
    assert _pooled_corr(tree, m, reps=250, seed=1) > 0.3


# --------------------------------------------------------------------------- multi-optimum OU
def test_multi_optimum_ou_tracks_local_optima():
    """Tips in each regime concentrate near that regime's optimum (strong alpha)."""
    tree = z.simulate_species_tree(z.Yule(1.0), n_tips=40, age=3.0, seed=7)
    regimes = z.simulate_traits(tree, z.Mk.equal_rates(2, 0.4), seed=1)
    tip_reg = {leaf: regimes.node_values[leaf] for leaf in tree.extant_leaves()}
    assert set(tip_reg.values()) == {0, 1}          # both regimes present at tips
    mou = z.MultiOptimumOU(regimes, theta=[-5.0, 5.0], alpha=4.0, sigma2=0.4)
    rng = np.random.default_rng(3)
    acc = {leaf: [] for leaf in tree.extant_leaves()}
    for _ in range(30):
        res = z.simulate_traits(tree, mou, rng=rng)
        for leaf in acc:
            acc[leaf].append(res.node_values[leaf])
    m0 = np.mean([np.mean(acc[leaf]) for leaf in acc if tip_reg[leaf] == 0])
    m1 = np.mean([np.mean(acc[leaf]) for leaf in acc if tip_reg[leaf] == 1])
    assert m0 < -2.0 and m1 > 2.0


def test_multi_optimum_ou_equal_thetas_reduces_to_single_ou():
    tree = _fixed_tree(n_tips=20, seed=5)
    regimes = z.simulate_traits(tree, z.Mk.equal_rates(2, 0.5), seed=1)
    mou = z.MultiOptimumOU(regimes, theta=[3.0, 3.0], alpha=2.0, sigma2=0.3)
    rng = np.random.default_rng(0)
    vals = []
    for _ in range(200):
        vals += list(z.simulate_traits(tree, mou, rng=rng).values.values())
    assert abs(np.mean(vals) - 3.0) < 0.2         # single shared optimum


def test_multi_optimum_ou_x0_defaults_to_root_regime_optimum():
    tree = _fixed_tree()
    regimes = z.simulate_traits(tree, z.Mk.equal_rates(2, 0.5, root=1), seed=1)
    mou = z.MultiOptimumOU(regimes, theta=[-2.0, 7.0], alpha=1.0, sigma2=0.1)
    assert mou.x0 == 7.0                            # root regime is 1


def test_multi_optimum_ou_per_regime_params_and_reproducible():
    tree = _fixed_tree()
    regimes = z.simulate_traits(tree, z.Mk.equal_rates(2, 0.5), seed=1)
    mou = z.MultiOptimumOU(regimes, theta=[0.0, 1.0], alpha=[1.0, 3.0], sigma2=[0.1, 0.5])
    a = z.simulate_traits(tree, mou, seed=9).values
    b = z.simulate_traits(tree, mou, seed=9).values
    assert {k.name: v for k, v in a.items()} == {k.name: v for k, v in b.items()}


def test_multi_optimum_ou_validation():
    tree = _fixed_tree()
    regimes = z.simulate_traits(tree, z.Mk.equal_rates(2, 0.5), seed=1)
    with pytest.raises(ValueError):
        z.MultiOptimumOU(regimes, theta=[1.0], alpha=1.0, sigma2=0.1)      # wrong theta length
    with pytest.raises(ValueError):
        z.MultiOptimumOU(regimes, theta=[1.0, 2.0], alpha=0.0, sigma2=0.1)  # alpha > 0
    cont = z.simulate_traits(tree, z.BrownianMotion(0.5), seed=1)
    with pytest.raises(ValueError):
        z.MultiOptimumOU(cont, theta=[1.0, 2.0], alpha=1.0, sigma2=0.1)     # regimes not discrete


# --------------------------------------------------------------------------- threshold model
def test_threshold_liability_continuous_state_derived():
    tree = _fixed_tree()
    res = z.simulate_traits(tree, z.ThresholdModel(thresholds=[0.0]), seed=1)
    assert res.kind == "continuous"
    assert isinstance(list(res.values.values())[0], float)   # values are liabilities
    for leaf in tree.extant_leaves():
        liability = res.node_values[leaf]
        assert res.labeled_values()[leaf] == (1 if liability > 0 else 0)


def test_threshold_binary_symmetric_is_balanced():
    tree = z.simulate_species_tree(z.Yule(1.0), n_tips=40, age=3.0, seed=7)
    m = z.ThresholdModel(thresholds=[0.0], x0=0.0)
    rng = np.random.default_rng(0)
    states = []
    for _ in range(300):
        res = z.simulate_traits(tree, m, rng=rng)
        states += [res.labeled_values()[leaf] for leaf in tree.extant_leaves()]
    assert abs(np.mean(states) - 0.5) < 0.05


def test_threshold_positive_x0_biases_toward_high_state():
    tree = z.simulate_species_tree(z.Yule(1.0), n_tips=30, age=1.0, seed=7)
    m = z.ThresholdModel(thresholds=[0.0], sigma2=1.0, x0=3.0)  # liability starts well above 0
    rng = np.random.default_rng(0)
    states = []
    for _ in range(200):
        res = z.simulate_traits(tree, m, rng=rng)
        states += [res.labeled_values()[leaf] for leaf in tree.extant_leaves()]
    assert np.mean(states) > 0.8


def test_threshold_ordered_multistate_discretize():
    m = z.ThresholdModel(thresholds=[-1.0, 1.0], states=["low", "mid", "high"])
    assert m.discretize(-2.0) == "low"
    assert m.discretize(0.0) == "mid"
    assert m.discretize(5.0) == "high"


def test_threshold_validation():
    with pytest.raises(ValueError):
        z.ThresholdModel(thresholds=[1.0, 0.0])       # not strictly increasing
    with pytest.raises(ValueError):
        z.ThresholdModel(thresholds=[])               # empty
    with pytest.raises(ValueError):
        z.ThresholdModel(thresholds=[0.0], sigma2=-1.0)


def test_threshold_tsv_reports_discrete_state():
    tree = _fixed_tree()
    res = z.simulate_traits(tree, z.ThresholdModel([0.0], states=["a", "b"]), seed=1)
    tsv = res.to_tsv()
    assert "'a'" in tsv or "'b'" in tsv           # discrete labels, not liabilities


# --------------------------------------------------------------------------- early burst / ACDC
def test_early_burst_reduces_to_bm_at_rate_zero():
    tree, tip = _single_branch(1.0)
    m = z.EarlyBurst(sigma2=0.7, rate=0.0)
    rng = np.random.default_rng(0)
    v = np.array([z.simulate_traits(tree, m, rng=rng).node_values[tip] for _ in range(20000)])
    assert abs(v.var() - 0.7 * 1.0) < 0.03        # var = sigma2 * dt


def test_early_burst_variance_matches_integral():
    t = 2.0
    tree, tip = _single_branch(t)
    for rate in (-1.5, 0.8):
        m = z.EarlyBurst(sigma2=0.5, rate=rate)
        rng = np.random.default_rng(1)
        v = np.array([z.simulate_traits(tree, m, rng=rng).node_values[tip] for _ in range(20000)])
        exp_var = 0.5 * (np.exp(rate * t) - 1.0) / rate
        assert abs(v.var() - exp_var) / exp_var < 0.06


def test_early_burst_deceleration_accumulates_less():
    """rate<0 accrues less total variance than BM (rate=0), rate>0 more (same root sigma2)."""
    t = 2.0
    tree, tip = _single_branch(t)

    def tip_var(rate):
        m = z.EarlyBurst(0.5, rate)
        rng = np.random.default_rng(2)
        return np.array([z.simulate_traits(tree, m, rng=rng).node_values[tip]
                         for _ in range(15000)]).var()

    assert tip_var(-2.0) < tip_var(0.0) < tip_var(2.0)


def test_early_burst_validation():
    with pytest.raises(ValueError):
        z.EarlyBurst(-1.0, 0.0)
    assert "EarlyBurst" in repr(z.EarlyBurst(1.0, -0.5))


# --------------------------------------------------------------------------- Pagel tree transforms
def _times_by_name(tree):
    return {n.name: n.time for n in tree.nodes()}


def _times_equal(a, b):
    ta, tb = _times_by_name(a), _times_by_name(b)
    return ta.keys() == tb.keys() and all(abs(ta[k] - tb[k]) < 1e-9 for k in ta)


def test_pagel_lambda_identity_and_star():
    tree = _fixed_tree(n_tips=6, seed=11)
    assert _times_equal(z.pagel_lambda(tree, 1.0), tree)          # lam=1 -> unchanged
    star = z.pagel_lambda(tree, 0.0)
    for n in star.internal_nodes():
        if n.parent is not None:
            assert n.time == 0.0                                  # internals collapse -> star


def test_pagel_lambda_scales_shared_paths():
    tree = _fixed_tree(n_tips=6, seed=11)
    lam = 0.4
    t2 = z.pagel_lambda(tree, lam)
    old = {n.name: n for n in tree.nodes()}
    new = {n.name: n for n in t2.nodes()}
    for name, o in old.items():
        if o.is_leaf():
            assert abs(new[name].time - o.time) < 1e-9            # tip depths preserved
        elif o.parent is not None:
            assert abs(new[name].time - lam * o.time) < 1e-9      # internal depths scaled
    assert t2.root.time == 0.0


def test_pagel_delta_preserves_root_and_tips():
    tree = _fixed_tree(n_tips=6, seed=3)
    T = tree.total_age
    assert _times_equal(z.pagel_delta(tree, 1.0), tree)           # delta=1 -> unchanged
    t2 = z.pagel_delta(tree, 2.0)
    assert t2.root.time == 0.0
    for leaf in t2.leaves():
        assert abs(leaf.time - T) < 1e-9                          # tips preserved
    old = {n.name: n for n in tree.nodes()}
    for name, n in ((n.name, n) for n in t2.internal_nodes()):
        o = old[name]
        if o.parent is not None:
            assert abs(n.time - T * (o.time / T) ** 2.0) < 1e-9


def test_pagel_kappa_speciational_unit_branches():
    tree = _fixed_tree(n_tips=6, seed=5)
    assert _times_equal(z.pagel_kappa(tree, 1.0), tree)           # kappa=1 -> unchanged
    t2 = z.pagel_kappa(tree, 0.0)
    for n in t2.nodes():
        if n.parent is not None:
            assert abs(n.branch_length() - 1.0) < 1e-9            # every branch length 1


def test_pagel_validation():
    tree = _fixed_tree()
    with pytest.raises(ValueError):
        z.pagel_lambda(tree, 1.5)
    with pytest.raises(ValueError):
        z.pagel_delta(tree, 0.0)
    with pytest.raises(ValueError):
        z.pagel_kappa(tree, -1.0)


def test_pagel_transform_feeds_simulation():
    tree = _fixed_tree(n_tips=8, seed=7)
    t2 = z.pagel_lambda(tree, 0.5)
    res = z.simulate_traits(t2, z.BrownianMotion(0.5), seed=1)
    assert set(res.values) == set(t2.extant_leaves())


# --------------------------------------------------------------------------- hidden-state Mk (corHMM)
def _slow_fast_hmm(hidden_rate=0.5):
    slow = [[0.0, 0.1], [0.1, 0.0]]
    fast = [[0.0, 3.0], [3.0, 0.0]]
    return z.HiddenStateMk(observed_rates=[slow, fast], hidden_rate=hidden_rate,
                           observed_states=[0, 1], hidden_states=["slow", "fast"])


def test_hidden_state_mk_collapses_to_observed():
    tree = _fixed_tree(n_tips=20, seed=2)
    res = z.simulate_traits(tree, _slow_fast_hmm(), seed=1)
    for leaf in tree.extant_leaves():
        obs = res.labeled_values()[leaf]
        full = res.full_label(res.node_values[leaf])
        assert obs in (0, 1)
        assert full[0] == obs and full[1] in ("slow", "fast")     # full = (observed, hidden)


def test_hidden_state_mk_discretize():
    m = _slow_fast_hmm()
    # states order: (0,slow)=0, (0,fast)=1, (1,slow)=2, (1,fast)=3
    assert m.discretize(0) == 0 and m.discretize(1) == 0
    assert m.discretize(2) == 1 and m.discretize(3) == 1


def test_hidden_state_mk_changes_show_hidden_switches():
    tree = _fixed_tree(n_tips=30, seed=4)
    res = z.simulate_traits(tree, _slow_fast_hmm(hidden_rate=1.0), seed=2)
    # at least one change keeps the observed state but switches hidden class
    assert any(frm[0] == to[0] and frm[1] != to[1] for _, _, frm, to in res.changes())


def test_hidden_state_mk_same_rates_are_observed_symmetric():
    tree = z.simulate_species_tree(z.Yule(1.0), n_tips=30, age=3.0, seed=7)
    same = [[0.0, 0.8], [0.8, 0.0]]
    m = z.HiddenStateMk(observed_rates=[same, same], hidden_rate=0.5)  # hidden irrelevant to observed
    rng = np.random.default_rng(0)
    obs = []
    for _ in range(100):
        res = z.simulate_traits(tree, m, rng=rng)
        obs += [res.labeled_values()[leaf] for leaf in tree.extant_leaves()]
    assert abs(np.mean(obs) - 0.5) < 0.08


def test_hidden_state_mk_validation():
    with pytest.raises(ValueError):
        z.HiddenStateMk([[[0, 1], [1, 0]], [[0, 1, 0], [0, 0, 1], [1, 0, 0]]], 0.5)  # mismatched O
    with pytest.raises(ValueError):
        z.HiddenStateMk([[[0, 1], [1, 0]]], hidden_rate=[[0, 1], [1, 0]])            # H=1 vs 2x2


def test_hidden_state_mk_reproducible():
    tree = _fixed_tree(seed=1)
    m = _slow_fast_hmm()
    a = z.simulate_traits(tree, m, seed=5).labeled_values()
    b = z.simulate_traits(tree, m, seed=5).labeled_values()
    assert {k.name: v for k, v in a.items()} == {k.name: v for k, v in b.items()}


# --------------------------------------------------------------------------- replicates
def test_replicate_traits_returns_n_independent_results():
    tree = _fixed_tree(n_tips=6, seed=1)
    reps = z.replicate_traits(tree, z.BrownianMotion(0.5), 4, seed=1)
    assert len(reps) == 4
    leaf = tree.extant_leaves()[0]
    values = [r.values[leaf] for r in reps]
    assert len(set(values)) == 4                     # independent draws differ


def test_replicate_traits_reproducible():
    tree = _fixed_tree(seed=1)
    a = z.replicate_traits(tree, z.BrownianMotion(0.5), 3, seed=7)
    b = z.replicate_traits(tree, z.BrownianMotion(0.5), 3, seed=7)
    for ra, rb in zip(a, b):
        assert {k.name: v for k, v in ra.values.items()} == \
               {k.name: v for k, v in rb.values.items()}


def test_replicate_traits_validation():
    with pytest.raises(ValueError):
        z.replicate_traits(_fixed_tree(), z.BrownianMotion(0.5), 0)
