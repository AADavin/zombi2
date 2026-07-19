"""Tests for the continuous trait core — zombi2.traits.simulate_continuous (Brownian motion).

The correctness-critical test is the **exact BM tip law** (Felsenstein 1985): the node-by-node
preorder walk must reproduce the multivariate-normal law over the extant tips, so across replicates
each tip has variance σ²·(root-to-tip depth) and each tip pair has covariance σ²·(shared path
length). Both are checked against the tree geometry, with fixed seeds so the statistics are
deterministic, not flaky.
"""

import numpy as np
import pytest

from zombi2.rates import modifiers as mod
from zombi2.rates import scope
from zombi2.species import simulate_species_tree
from zombi2.traits import TraitsResult, simulate_continuous


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

def test_rejects_a_modifier():
    sp = _tree(seed=1)
    with pytest.raises(ValueError, match="Time"):
        simulate_continuous(sp, rate=1.0 * mod.Time({0: 1.0, 3: 0.2}), seed=1)


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
