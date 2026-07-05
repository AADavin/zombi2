"""Tests for the traits<->genes JOINT model — "trait-gene feedback"
(``zombi2.trait_gene_feedback``).

A continuous trait and a family panel are coupled *both* ways: the trait sets the panel's retention
(``traits:genes``) and the panel sets the trait's OU optimum (``genes:traits``). The key check is
**emergent**: with the loop closed the tips end up with the trait and the panel correlated, whereas a
decoupled control (neither arrow) shows no association — the feedback *writes* the signal.
"""

import numpy as np
import pytest

import zombi2 as z
from zombi2.coevolve import TraitGeneFeedback, simulate_trait_gene_feedback


def _tree(n_tips=120, age=8.0, seed=1):
    return z.simulate_species_tree(z.BirthDeath(1.0, 0.3), n_tips=n_tips, age=age, seed=seed)


def _mean_corr(model, tree, seeds=range(10)):
    vals = [simulate_trait_gene_feedback(tree, model, seed=s).trait_gene_correlation() for s in seeds]
    return float(np.nanmean(vals))


# --------------------------------------------------------------------------- validation
def test_tgf_validation():
    with pytest.raises(ValueError):
        TraitGeneFeedback(n_families=0)
    with pytest.raises(ValueError):
        TraitGeneFeedback(effect_loss=-1.0)
    with pytest.raises(ValueError):
        TraitGeneFeedback(root_fraction=1.5)
    with pytest.raises(ValueError):
        TraitGeneFeedback(steps=0)


def test_tgf_optimum_is_monotone_in_panel():
    m = TraitGeneFeedback(n_families=10, theta_low=-2.0, theta_high=4.0)
    assert m.optimum(0) == pytest.approx(-2.0)
    assert m.optimum(10) == pytest.approx(4.0)
    assert m.optimum(5) == pytest.approx(1.0)


def test_tgf_reproducible():
    tree = _tree()
    m = TraitGeneFeedback(n_families=16, effect_loss=1.5, theta_low=-3, theta_high=3)
    a = simulate_trait_gene_feedback(tree, m, seed=7)
    b = simulate_trait_gene_feedback(tree, m, seed=7)
    assert a.node_trait == b.node_trait
    assert a.profiles.to_tsv() == b.profiles.to_tsv()


def test_tgf_structure():
    tree = _tree(n_tips=60)
    m = TraitGeneFeedback(n_families=12)
    res = simulate_trait_gene_feedback(tree, m, seed=3)
    assert set(res.node_trait) == set(tree.nodes())
    assert res.profiles.shape == (12, len(tree.extant_leaves()))
    # panel occupancy is a fraction in [0, 1]
    assert all(0.0 <= v <= 1.0 for v in res.panel_occupancy().values())


# --------------------------------------------------------------------------- the emergent signal
def test_feedback_writes_a_trait_gene_correlation():
    """Full feedback -> tips correlate trait and panel; a decoupled control shows ~no association."""
    tree = _tree()
    coupled = TraitGeneFeedback(n_families=24, effect_loss=1.5, base_loss=1.0, gain=1.0,
                                theta_low=-3.0, theta_high=3.0, alpha=1.0, sigma2=0.5, steps=24)
    decoupled = TraitGeneFeedback(n_families=24, effect_loss=0.0, base_loss=1.0, gain=1.0,
                                  theta_low=0.0, theta_high=0.0, alpha=1.0, sigma2=0.5, steps=24)
    r_on = _mean_corr(coupled, tree)
    r_off = _mean_corr(decoupled, tree)
    assert r_on > 0.35                         # the loop induces a clear association
    assert abs(r_off) < 0.2                    # no coupling -> no association
    assert r_on > r_off + 0.25


def test_single_edge_limits_run():
    """Each single edge is a limit of the joint model and still simulates."""
    tree = _tree(n_tips=50)
    # effect_loss = 0 -> pure genes:traits (panel ignores the trait)
    only_g2t = simulate_trait_gene_feedback(
        tree, TraitGeneFeedback(n_families=12, effect_loss=0.0), seed=1)
    # theta_high == theta_low -> pure traits:genes (trait ignores the panel)
    only_t2g = simulate_trait_gene_feedback(
        tree, TraitGeneFeedback(n_families=12, theta_low=1.0, theta_high=1.0), seed=1)
    assert only_g2t.profiles.shape[0] == 12
    assert only_t2g.profiles.shape[0] == 12
