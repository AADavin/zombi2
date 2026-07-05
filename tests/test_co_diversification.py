"""Tests for the species<->genes JOINT model — "co-diversification" (``genes:species`` +
``species:genes``), the genomic analogue of ClaSSE.

The same driver panel both **drives diversification** (``genes:species``) and is **reshuffled by a
cladogenetic burst at each speciation** (``species:genes``). The checks:

* the burst genuinely differentiates sister lineages (the punctuational ``species:genes`` signal),
  which the plain ``genes:species`` process cannot do (daughters are born identical);
* the drivers still drive the tree under the joint model (a key innovation stays over-represented);
* with both burst probabilities 0 the model reduces *exactly* to :func:`simulate_gene_diversification`;
* validation and deterministic corners.
"""

import numpy as np
import pytest

import zombi2 as z
from zombi2.coevolve import simulate_co_diversification


def _pure_punctuational(**kw):
    """A co-diversification model whose drivers change ONLY at cladogenetic bursts (anagenetic off)."""
    base = dict(loss=0.0, origination=0.0, transfer=0.0)
    base.update(kw)
    return z.GeneDiversification(**base)


def _sister_diff_fraction(res):
    """Fraction of bifurcating nodes whose two daughters carry *different* driver sets."""
    nd, diff, tot = res.node_drivers, 0, 0
    for node in res.tree.nodes():
        if len(node.children) == 2:
            tot += 1
            if nd[node.children[0]] != nd[node.children[1]]:
                diff += 1
    return diff / tot if tot else 0.0


# --------------------------------------------------------------------------- validation
def test_codiv_validation():
    with pytest.raises(ValueError):
        z.GeneDiversification(2, cladogenetic_loss=1.5)               # prob out of [0, 1]
    with pytest.raises(ValueError):
        z.GeneDiversification(2, cladogenetic_gain=-0.1)              # prob out of [0, 1]
    # simulate_co_diversification refuses a model with no burst (that is plain genes:species)
    with pytest.raises(ValueError):
        simulate_co_diversification(z.GeneDiversification(2, root_drivers=1), n_tips=10, seed=1)
    assert not z.GeneDiversification(2).is_co_diversification
    assert z.GeneDiversification(2, cladogenetic_loss=0.1).is_co_diversification


def test_codiv_reproducible():
    m = _pure_punctuational(n_drivers=3, driver_speciation=1.0, root_drivers=1,
                            cladogenetic_loss=0.15, cladogenetic_gain=0.2)
    a = simulate_co_diversification(m, n_tips=60, seed=11)
    b = simulate_co_diversification(m, n_tips=60, seed=11)
    assert a.tree.to_newick() == b.tree.to_newick()
    assert a.to_tsv() == b.to_tsv()


def test_codiv_reduces_to_genes_species_when_no_burst():
    """cladogenetic_loss = cladogenetic_gain = 0 -> identical stream to plain genes:species."""
    common = dict(n_drivers=2, driver_speciation=1.0, transfer=0.6, root_drivers=1)
    plain = z.simulate_gene_diversification(z.GeneDiversification(**common), n_tips=50, seed=9)
    zero_burst = z.simulate_gene_diversification(
        z.GeneDiversification(cladogenetic_loss=0.0, cladogenetic_gain=0.0, **common),
        n_tips=50, seed=9)
    assert plain.tree.to_newick() == zero_burst.tree.to_newick()
    assert plain.to_tsv() == zero_burst.to_tsv()


# --------------------------------------------------------------------------- the joint signature
def test_burst_differentiates_sisters():
    """The species:genes burst makes sister lineages differ; without it they are born identical."""
    kw = dict(n_drivers=4, driver_speciation=0.6, root_drivers=2)
    burst = simulate_co_diversification(
        _pure_punctuational(cladogenetic_loss=0.2, cladogenetic_gain=0.15, **kw), n_tips=150, seed=4)
    none = z.simulate_gene_diversification(_pure_punctuational(**kw), n_tips=150, seed=4)
    assert _sister_diff_fraction(burst) > 0.2       # bursts split sisters at many nodes
    assert _sister_diff_fraction(none) == 0.0       # anagenetic-off, no burst -> sisters identical


def test_driver_still_drives_under_the_joint_model():
    """A key innovation stays over-represented among the tips even with the burst active."""
    def mean_prev(driver_speciation, reps=30, seed=0):
        rng = np.random.default_rng(seed)
        vals = []
        for _ in range(reps):
            m = z.GeneDiversification(1, lambda0=1.0, mu0=0.2, driver_speciation=driver_speciation,
                                      loss=0.0, origination=0.0, transfer=0.0, root_drivers=1,
                                      cladogenetic_loss=0.1, cladogenetic_gain=0.1)
            vals.append(simulate_co_diversification(m, n_tips=70, rng=rng).tip_prevalence()[0])
        return float(np.mean(vals))
    assert mean_prev(1.2) > mean_prev(0.0) + 0.1


def test_codiv_no_channel_never_appears():
    """No gain of any kind (clado_gain = origination = transfer = 0, none at root) -> prevalence 0,
    even though a loss burst is active."""
    m = z.GeneDiversification(2, loss=0.0, origination=0.0, transfer=0.0, root_drivers=0,
                              cladogenetic_loss=0.3, cladogenetic_gain=0.0)
    res = simulate_co_diversification(m, n_tips=40, seed=2)
    assert res.tip_prevalence() == [0.0, 0.0]


def test_codiv_tree_wellformed():
    # n_tips mode bounds the run (a fixed driver with loss=0 would otherwise run away under age mode)
    res = simulate_co_diversification(
        _pure_punctuational(n_drivers=2, driver_speciation=0.8, root_drivers=1,
                            cladogenetic_loss=0.1, cladogenetic_gain=0.1), n_tips=60, seed=3)
    tree = res.tree
    assert tree.root.time == 0.0 and len(tree.root.children) == 2
    assert len(tree.extant_leaves()) == 60
    for node in tree.internal_nodes():
        assert len(node.children) == 2
    for leaf in tree.extant_leaves():
        assert abs(leaf.time - tree.total_age) < 1e-9
