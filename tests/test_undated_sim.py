"""The forward undated simulator — the generative twin of the undated-DTL likelihood.

Correctness is pinned three ways: (1) the joint likelihood of the simulated survivors, evaluated
under the undated model, peaks at the *generating* odds on every axis — the defining property of a
generative dual; (2) the same ``pS^(2k-1)`` closed form the likelihood oracle uses is recovered
here as an empirical event frequency; (3) each ground-truth reconciliation scores 1.0 against
itself through ``recon-accuracy``, proving the emitted annotation is well-formed.
"""

import math

import pytest

from zombi2 import read_newick
from zombi2.tools.reconciliation import (
    GeneTree,
    UndatedDTL,
    simulate_undated,
    undated_joint_loglik,
)
from zombi2.tools.recon_accuracy import reconciliation_accuracy

SP4 = "((A:1,B:1)X:1,(C:1,D:1)Y:1)root:1;"
CLADO = "((A,B)X,(C,D)Y)root;"  # no branch lengths


def test_deterministic():
    m = UndatedDTL(0.3, 0.2, 0.4)
    a = simulate_undated(read_newick(SP4), m, n_families=200, seed=42)
    b = simulate_undated(read_newick(SP4), m, n_families=200, seed=42)
    assert [r.extant for r in a.reconciliations] == [r.extant for r in b.reconciliations]
    c = simulate_undated(read_newick(SP4), m, n_families=200, seed=43)
    assert [r.extant for r in a.reconciliations] != [r.extant for r in c.reconciliations]


def test_cladogram_runs_without_dates():
    # a length-less cladogram must simulate (unit branches assumed), not silently do nothing
    res = simulate_undated(read_newick(CLADO), UndatedDTL(0.2, 0.1, 0.3), n_families=100, seed=1)
    assert res.n_surviving > 0
    assert res.event_counts["S"] > 0


def test_no_dup_transfer_when_odds_zero():
    res = simulate_undated(read_newick(SP4), UndatedDTL(0.0, 0.0, 0.5), n_families=300, seed=2)
    assert res.event_counts.get("D", 0) == 0
    assert res.event_counts.get("T", 0) == 0
    assert res.event_counts["S"] > 0 and res.n_surviving > 0


def test_roundtrip_score_is_finite():
    m = UndatedDTL(0.2, 0.12, 0.5)
    res = simulate_undated(read_newick(SP4), m, n_families=400, seed=3)
    ll = undated_joint_loglik(res.gene_trees(), res.species_tree, m, origination="root",
                              transfers="global", n_extinct=res.n_extinct, backend="python")
    assert math.isfinite(ll)


def test_likelihood_peaks_at_truth():
    """A generative dual: the joint undated loglik at the true odds beats clearly-wrong odds
    (2x and 0.5x) on every axis."""
    true = dict(dup=0.25, transfer=0.15, loss=0.55)
    res = simulate_undated(read_newick(SP4), UndatedDTL(**true), n_families=1500, seed=4)
    gts = res.gene_trees()

    def jll(**kw):
        return undated_joint_loglik(gts, res.species_tree, UndatedDTL(**kw), origination="root",
                                    transfers="global", n_extinct=res.n_extinct, backend="python")

    base = jll(**true)
    for p in ("dup", "transfer", "loss"):
        hi = {**true, p: true[p] * 2.0}
        lo = {**true, p: true[p] * 0.5}
        assert base > jll(**hi), f"{p}: truth not above 2x"
        assert base > jll(**lo), f"{p}: truth not above 0.5x"


def test_matches_undated_oracle_frequency():
    """d=t=0 on (A,B): a family survives in BOTH tips with probability pS^3 (root speciates, each
    tip is sampled) — the pS^(2k-1) oracle for k=2. The sampler must reproduce it as a frequency."""
    loss = 1.0
    pS = 1.0 / (1.0 + loss)
    res = simulate_undated(read_newick("(A:1,B:1)root:1;"), UndatedDTL(0.0, 0.0, loss),
                           n_families=20000, seed=5)
    both = sum(1 for r in res.reconciliations
               if r.extant is not None and GeneTree.from_newick(r.extant).species_set() == {"A", "B"})
    freq = both / res.n_families
    assert abs(freq - pS ** 3) < 0.01, (freq, pS ** 3)


def test_reldated_requires_dates():
    with pytest.raises(ValueError):
        simulate_undated(read_newick(CLADO), UndatedDTL(0.2, 0.1, 0.3), n_families=5,
                         transfers="dated")
    res = simulate_undated(read_newick(SP4), UndatedDTL(0.2, 0.1, 0.3), n_families=50,
                           transfers="dated", seed=6)
    assert res.n_families == 50


def test_max_events_guard():
    # supercritical odds (dup+transfer >> loss) make families explode; the guard must fire
    with pytest.raises(RuntimeError):
        simulate_undated(read_newick(SP4), UndatedDTL(5.0, 5.0, 0.01), n_families=100, seed=7,
                         max_events=5000)


def test_self_recon_accuracy_is_perfect():
    res = simulate_undated(read_newick(SP4), UndatedDTL(0.3, 0.2, 0.3), n_families=200, seed=8)
    checked = 0
    for r in res.reconciliations:
        if r.extant is None:
            continue
        acc = reconciliation_accuracy(r, r)
        if acc.n_nodes > 0:  # single-survivor families have no internal node to score
            assert acc.event_accuracy == 1.0
            assert acc.mapping_accuracy == 1.0
            checked += 1
    assert checked > 0


def test_profiles_match_extant_tips():
    res = simulate_undated(read_newick(SP4), UndatedDTL(0.3, 0.2, 0.3), n_families=150, seed=9)
    rows = dict(res.profile_rows())
    assert len(rows) == res.n_surviving           # one row per surviving family
    assert res.leaf_names == ["A", "B", "C", "D"]  # columns = sorted species
    for i, r in enumerate(res.reconciliations, 1):
        if r.extant is None:
            assert str(i) not in rows              # extinct families are omitted
            continue
        n_tips = sum(g.is_leaf for g in GeneTree.from_newick(r.extant).nodes)
        assert sum(rows[str(i)]) == n_tips         # total copies == extant gene-tree tips
