"""Pagel correlated evolution generalized to k binary characters (CorrelatedBinaryK).

These tests pin the k-trait Pagel model (2**k states, one-flip-at-a-time) against the
trusted 2-trait :class:`CorrelatedBinary`, verify the Kronecker-sum structure of the
independent null, the single-change invariant, the state<->tuple codec, the convenience
constructors, and end-to-end simulation + coupling recovery via :func:`simulate_traits`.
"""

import numpy as np
import pytest

import zombi2 as z
from zombi2.traits import CorrelatedBinaryK


# --------------------------------------------------------------------------- helpers
def _kron_sum(mats):
    """Kronecker sum of the given square matrices (outer index first)."""
    n = int(np.prod([m.shape[0] for m in mats]))
    Q = np.zeros((n, n))
    eyes = [np.eye(m.shape[0]) for m in mats]
    for idx, m in enumerate(mats):
        factors = [eyes[j] if j != idx else m for j in range(len(mats))]
        term = factors[0]
        for f in factors[1:]:
            term = np.kron(term, f)
        Q += term
    return Q


def _two_state_chain(gain, loss):
    """CTMC of a lone binary trait: state 0 = absent, 1 = present."""
    return np.array([[-gain, gain], [loss, -loss]])


def _fixed_tree():
    return z.simulate_species_tree(z.Yule(1.0), n_tips=6, age=2.0, seed=3)


# --------------------------------------------------------------------------- codec
def test_state_tuple_codec_roundtrip():
    for k in (2, 3, 4):
        n = 1 << k
        for i in range(n):
            t = CorrelatedBinaryK.index_to_tuple(i, k)
            assert len(t) == k
            assert set(t) <= {0, 1}
            assert CorrelatedBinaryK.tuple_to_index(t) == i


def test_state_encoding_is_big_endian():
    # trait 0 is the most-significant bit (matches CorrelatedBinary: index = 2*X + Y)
    assert CorrelatedBinaryK.index_to_tuple(0, 3) == (0, 0, 0)
    assert CorrelatedBinaryK.index_to_tuple(4, 3) == (1, 0, 0)  # only trait 0 set
    assert CorrelatedBinaryK.index_to_tuple(1, 3) == (0, 0, 1)  # only trait 2 set
    assert CorrelatedBinaryK.tuple_to_index((1, 0, 1)) == 5


def test_tuple_to_index_rejects_non_binary():
    with pytest.raises(ValueError):
        CorrelatedBinaryK.tuple_to_index((0, 2))


# --------------------------------------------------------------------------- k=2 cross-check
def test_k2_reproduces_correlated_binary_Q():
    """For k=2 the general class rebuilds the exact CorrelatedBinary rate matrix."""
    cb = z.CorrelatedBinary(x_gain_y0=1, x_gain_y1=2, x_loss_y0=3, x_loss_y1=4,
                            y_gain_x0=5, y_gain_x1=6, y_loss_x0=7, y_loss_x1=8)

    def rate_fn(trait, direction, others):
        if trait == 0:                 # X, others = (Y,)
            y = others[0]
            if direction == "gain":
                return 2 if y else 1   # x_gain_y1 / x_gain_y0
            return 4 if y else 3       # x_loss_y1 / x_loss_y0
        x = others[0]                  # Y, others = (X,)
        if direction == "gain":
            return 6 if x else 5       # y_gain_x1 / y_gain_x0
        return 8 if x else 7           # y_loss_x1 / y_loss_x0

    cbk = CorrelatedBinaryK(2, rate_fn)
    assert cbk.states == cb.states               # same state ordering
    assert np.allclose(cbk.Q, cb.Q)              # identical Q


def test_k2_independent_matches_correlated_binary_independent():
    cb = z.CorrelatedBinary.independent(x_gain=0.7, x_loss=0.4, y_gain=0.9, y_loss=0.5)
    cbk = CorrelatedBinaryK.independent(gains=[0.7, 0.9], losses=[0.4, 0.5])
    assert np.allclose(cbk.Q, cb.Q)


# --------------------------------------------------------------------------- independent = Kronecker sum
@pytest.mark.parametrize("k", [2, 3, 4])
def test_independent_is_kronecker_sum(k):
    rng = np.random.default_rng(k)
    gains = rng.uniform(0.1, 1.5, k)
    losses = rng.uniform(0.1, 1.5, k)
    m = CorrelatedBinaryK.independent(gains=gains, losses=losses)
    chains = [_two_state_chain(gains[i], losses[i]) for i in range(k)]
    assert np.allclose(m.Q, _kron_sum(chains))


def test_independent_scalar_broadcasts_against_a_vector():
    # a scalar gain/loss is broadcast to the length of the other (vector) argument
    m = CorrelatedBinaryK.independent(gains=0.6, losses=[0.3, 0.3, 0.3])
    assert m.n_traits == 3
    chains = [_two_state_chain(0.6, 0.3)] * 3
    assert np.allclose(m.Q, _kron_sum(chains))


def test_independent_two_scalars_cannot_infer_k():
    # two bare scalars give k=1, which is below the 2-trait minimum -> error.
    # (use equal_rates(k, ...) when you only have scalar rates.)
    with pytest.raises(ValueError):
        CorrelatedBinaryK.independent(gains=0.6, losses=0.3)


def test_equal_rates_is_symmetric_independent():
    k = 3
    m = CorrelatedBinaryK.equal_rates(k, gain=0.5)          # loss defaults to gain
    chains = [_two_state_chain(0.5, 0.5)] * k
    assert np.allclose(m.Q, _kron_sum(chains))
    m2 = CorrelatedBinaryK.equal_rates(k, gain=0.8, loss=0.2)
    chains2 = [_two_state_chain(0.8, 0.2)] * k
    assert np.allclose(m2.Q, _kron_sum(chains2))


# --------------------------------------------------------------------------- one-flip invariant
@pytest.mark.parametrize("k", [2, 3, 4, 5])
def test_only_single_bit_transitions_are_nonzero(k):
    """Pagel constraint: every multi-bit off-diagonal of Q is exactly 0."""
    rng = np.random.default_rng(100 + k)
    # a fully-general random rate table so no coincidental zeros hide a leak
    table = {}
    for trait in range(k):
        for direction in ("gain", "loss"):
            for cfg in range(1 << (k - 1)):
                others = tuple((cfg >> b) & 1 for b in range(k - 1))
                table[(trait, direction, others)] = rng.uniform(0.1, 2.0)
    m = CorrelatedBinaryK.from_table(k, table)
    n = 1 << k
    for i in range(n):
        for j in range(n):
            if i != j and bin(i ^ j).count("1") > 1:
                assert m.Q[i, j] == 0.0
    assert np.allclose(m.Q.sum(axis=1), 0.0)  # rows sum to zero


def test_rows_sum_to_zero_and_offdiag_nonneg():
    m = CorrelatedBinaryK.equal_rates(3, gain=0.5, loss=0.7)
    off = m.Q.copy()
    np.fill_diagonal(off, 0.0)
    assert np.all(off >= 0)
    assert np.allclose(m.Q.sum(axis=1), 0.0)


# --------------------------------------------------------------------------- marginals of independent
def test_independent_marginal_matches_lone_mk_via_Pt():
    """Each trait's marginal transition law equals a lone 2-state Mk's P(t) = exp(Qt)."""
    k = 3
    gains = [0.3, 0.9, 0.5]
    losses = [0.6, 0.2, 0.8]
    m = CorrelatedBinaryK.independent(gains=gains, losses=losses)
    t = 0.7
    P = m.transition_matrix(t)
    states = m.states
    for trait in range(k):
        lone = z.Mk(_two_state_chain(gains[trait], losses[trait]))
        Plone = lone.transition_matrix(t)
        # marginal P over `trait`, starting from the all-zero joint state
        start = 0  # (0,0,...,0)
        marg = np.zeros(2)
        for j, s in enumerate(states):
            marg[s[trait]] += P[start, j]
        assert np.allclose(marg, Plone[0], atol=1e-9)


def test_independent_marginal_matches_lone_mk_by_simulation():
    """Statistical check: simulated tip marginals of a trait under independent(...) match a
    lone 2-state Mk run on the same tree (fixed seed, loose tolerance)."""
    tree = z.simulate_species_tree(z.Yule(1.0), n_tips=25, age=3.0, seed=7)
    gains = [0.5, 0.5, 0.5]
    losses = [0.5, 0.5, 0.5]
    m = CorrelatedBinaryK.independent(gains=gains, losses=losses, root="stationary")
    lone = z.Mk(_two_state_chain(0.5, 0.5), root="stationary")

    rng = np.random.default_rng(1)
    reps = 200
    joint_freq = np.zeros(3)  # per-trait fraction present, trait 0/1/2
    joint_tips = 0
    for _ in range(reps):
        res = z.simulate_traits(tree, m, rng=rng)
        for v in res.values.values():
            cfg = res.label(v)
            joint_freq += np.array(cfg)
            joint_tips += 1
    joint_freq /= joint_tips

    rng2 = np.random.default_rng(1)
    lone_present = 0
    lone_tips = 0
    for _ in range(reps):
        res = z.simulate_traits(tree, lone, rng=rng2)
        for v in res.values.values():
            lone_present += res.label(v)
            lone_tips += 1
    lone_freq = lone_present / lone_tips

    # each trait's marginal frequency should sit near the lone-Mk frequency (~0.5)
    for f in joint_freq:
        assert abs(f - lone_freq) < 0.08


# --------------------------------------------------------------------------- simulation shape / integration
def test_simulate_returns_ktuples_and_integrates():
    tree = _fixed_tree()
    m = CorrelatedBinaryK.equal_rates(3, gain=0.6, loss=0.4)
    res = z.simulate_traits(tree, m, seed=2)
    lv = res.labeled_values()
    assert lv, "expected some extant tips"
    for tup in lv.values():
        assert isinstance(tup, tuple) and len(tup) == 3
        assert set(tup) <= {0, 1}
    # ancestral states are k-tuples too
    for node, v in res.ancestral_states().items():
        assert len(res.full_label(v)) == 3
    # changes() yields transitions between full k-tuple labels differing in one bit
    for (_node, _t, frm, to) in res.changes():
        assert len(frm) == 3 and len(to) == 3
        assert sum(a != b for a, b in zip(frm, to)) == 1


def test_root_tuple_pins_configuration():
    tree = _fixed_tree()
    # frozen (all rates 0) at (1,0,1) -> every node stays (1,0,1)
    m = CorrelatedBinaryK.equal_rates(3, gain=0.0, loss=0.0, root=(1, 0, 1))
    res = z.simulate_traits(tree, m, seed=5)
    assert res.node_values[tree.root] == CorrelatedBinaryK.tuple_to_index((1, 0, 1))
    assert all(lab == (1, 0, 1) for lab in res.labeled_values().values())


# --------------------------------------------------------------------------- partner coupling
def test_partner_coupling_reduces_to_independent_without_partners():
    m = CorrelatedBinaryK.partner_coupling(
        gains=[0.4, 0.5, 0.6], losses=[0.3, 0.2, 0.1],
        partners=[None, None, None])
    ind = CorrelatedBinaryK.independent(gains=[0.4, 0.5, 0.6], losses=[0.3, 0.2, 0.1])
    assert np.allclose(m.Q, ind.Q)


def test_partner_coupling_structure():
    # trait 1's gain is boosted 10x when its partner trait 0 is present
    m = CorrelatedBinaryK.partner_coupling(
        gains=[0.5, 0.1, 0.5], losses=[0.5, 0.5, 0.5],
        partners=[None, 0, None], boost_gain=[1, 10, 1])
    # trait-1 gain when trait0=0 -> 0.1 ; when trait0=1 -> 1.0
    # state (0,0,0) idx0 -> flip trait1 -> (0,1,0)
    i_000 = CorrelatedBinaryK.tuple_to_index((0, 0, 0))
    j_010 = CorrelatedBinaryK.tuple_to_index((0, 1, 0))
    assert np.isclose(m.Q[i_000, j_010], 0.1)
    i_100 = CorrelatedBinaryK.tuple_to_index((1, 0, 0))
    j_110 = CorrelatedBinaryK.tuple_to_index((1, 1, 0))
    assert np.isclose(m.Q[i_100, j_110], 1.0)


def _pair_cooccurrence(tree, model, pair, reps, seed):
    """Fraction of tips (pooled over reps) where both traits in `pair` are present."""
    rng = np.random.default_rng(seed)
    both = tips = 0
    a, b = pair
    for _ in range(reps):
        res = z.simulate_traits(tree, model, rng=rng)
        for v in res.values.values():
            cfg = res.label(v)
            both += 1 if (cfg[a] == 1 and cfg[b] == 1) else 0
            tips += 1
    return both / tips


def _pair_corr(tree, model, pair, reps, seed):
    rng = np.random.default_rng(seed)
    xs, ys = [], []
    a, b = pair
    for _ in range(reps):
        res = z.simulate_traits(tree, model, rng=rng)
        for v in res.values.values():
            cfg = res.label(v)
            xs.append(cfg[a])
            ys.append(cfg[b])
    if np.std(xs) == 0 or np.std(ys) == 0:
        return 0.0
    return float(np.corrcoef(xs, ys)[0, 1])


def test_strong_coupling_recovers_association():
    """A coupled pair (trait 2 tracks trait 1) co-occurs far more than under independent,
    while an uncoupled trait stays unassociated."""
    tree = z.simulate_species_tree(z.Yule(1.0), n_tips=30, age=3.0, seed=7)

    # trait 2 strongly tracks trait 1: gains fast when partner=1, lost fast when partner=0.
    def rate_fn(trait, direction, others):
        if trait == 2:                      # others = (t0, t1)
            partner = others[1]             # trait 1 in position (2->1 shift)
            if direction == "gain":
                return 2.0 if partner == 1 else 0.05
            return 0.05 if partner == 1 else 2.0
        # traits 0 and 1 evolve on their own at a moderate symmetric rate
        return 0.5

    coupled = CorrelatedBinaryK(3, rate_fn)
    independent = CorrelatedBinaryK.independent(gains=[0.5, 0.5, 0.5],
                                                losses=[0.5, 0.5, 0.5])

    corr_coupled = _pair_corr(tree, coupled, (1, 2), reps=200, seed=1)
    corr_indep = _pair_corr(tree, independent, (1, 2), reps=200, seed=1)
    assert corr_coupled > 0.3
    assert abs(corr_indep) < 0.1

    # the uncoupled pair (0, 2) shows no strong association even under the coupled model
    corr_uncoupled = _pair_corr(tree, coupled, (0, 2), reps=200, seed=2)
    assert abs(corr_uncoupled) < 0.15


def test_partner_coupling_induces_cooccurrence():
    """The ergonomic partner_coupling constructor also produces tip co-occurrence."""
    tree = z.simulate_species_tree(z.Yule(1.0), n_tips=30, age=3.0, seed=7)
    # trait 1 partners trait 0: strongly boosted gain + strongly damped loss when partner present
    coupled = CorrelatedBinaryK.partner_coupling(
        gains=[0.5, 0.05, 0.5], losses=[0.5, 2.0, 0.5],
        partners=[None, 0, None], boost_gain=[1, 40, 1], boost_loss=[1, 0.02, 1])
    independent = CorrelatedBinaryK.independent(gains=[0.5, 0.05, 0.5],
                                                losses=[0.5, 2.0, 0.5])
    co_coupled = _pair_cooccurrence(tree, coupled, (0, 1), reps=200, seed=1)
    co_indep = _pair_cooccurrence(tree, independent, (0, 1), reps=200, seed=1)
    assert co_coupled > co_indep + 0.05


# --------------------------------------------------------------------------- validation
def test_negative_rate_rejected():
    def bad(trait, direction, others):
        return -1.0
    with pytest.raises(ValueError):
        CorrelatedBinaryK(2, bad)


def test_k_too_small_rejected():
    with pytest.raises(ValueError):
        CorrelatedBinaryK(1, lambda t, d, o: 1.0)


def test_partner_index_out_of_range_rejected():
    with pytest.raises(ValueError):
        CorrelatedBinaryK.partner_coupling(gains=[1, 1], losses=[1, 1], partners=[5, None])


def test_repr():
    assert "k=3" in repr(CorrelatedBinaryK.equal_rates(3))


# --------------------------------------------------------------------------- re-export contract
def test_reexported_at_top_level_and_namespace():
    assert z.CorrelatedBinaryK is CorrelatedBinaryK
    import zombi2.traits as traits_ns
    assert traits_ns.CorrelatedBinaryK is CorrelatedBinaryK
    assert "CorrelatedBinaryK" in z.__all__
