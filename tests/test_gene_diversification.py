"""Tests for gene-content-dependent diversification (``zombi2.gene_diversification``).

The joint tree+gene-content forward process (the ``genes:species`` edge) is checked against its
expected signal:

* a **speciation** driver biases the standing tips *toward* carrying it; an **extinction** driver
  biases them *away*;
* **transfer** spreads a driver (frequency-dependent gain), so it reaches higher prevalence than
  with no HGT;
* deterministic corners: no gain channel -> a driver never appears; no loss + present at root ->
  it is fixed;
* the tree is well-formed and the run is reproducible.
"""

import numpy as np
import pytest

import zombi2 as z


def _mean_prevalence(reps=40, seed=0, n_tips=80, driver=0, **kw):
    """Mean fraction of extant tips carrying ``driver`` over independent trees."""
    rng = np.random.default_rng(seed)
    vals = []
    for _ in range(reps):
        m = z.GeneDiversification(1, **kw)
        vals.append(z.simulate_gene_diversification(m, n_tips=n_tips, rng=rng).tip_prevalence()[driver])
    return float(np.mean(vals))


# --------------------------------------------------------------------------- validation
def test_gene_div_validation():
    with pytest.raises(ValueError):
        z.GeneDiversification(0)                                  # need >= 1 driver
    with pytest.raises(ValueError):
        z.GeneDiversification(2, lambda0=-1)                      # negative base rate
    with pytest.raises(ValueError):
        z.GeneDiversification(2, transfer=-0.1)                   # negative transfer
    with pytest.raises(ValueError):
        z.GeneDiversification(2, root_drivers=3)                  # more root drivers than families
    with pytest.raises(ValueError):
        z.GeneDiversification(2, driver_speciation=[1.0, 2.0, 3.0])   # wrong-length effect vector
    # a scalar or a correctly-sized vector are both fine (scalar broadcasts to all K)
    assert z.GeneDiversification(2, driver_speciation=1.0).beta_lambda.tolist() == [1.0, 1.0]
    assert z.GeneDiversification(2, driver_speciation=[0.5, 1.5]).beta_lambda.tolist() == [0.5, 1.5]
    with pytest.raises(ValueError):
        z.simulate_gene_diversification(z.GeneDiversification(2))              # no stop condition
    with pytest.raises(ValueError):
        z.simulate_gene_diversification(z.GeneDiversification(2), age=3, n_tips=5)  # both


def test_gene_div_reproducible():
    m = z.GeneDiversification(2, driver_speciation=1.0, transfer=0.6, root_drivers=1)
    a = z.simulate_gene_diversification(m, n_tips=50, seed=7)
    b = z.simulate_gene_diversification(m, n_tips=50, seed=7)
    assert a.tree.to_newick() == b.tree.to_newick()
    assert a.to_tsv() == b.to_tsv()


# --------------------------------------------------------------------------- tree structure / stopping
def test_gene_div_tree_wellformed():
    res = z.simulate_gene_diversification(
        z.GeneDiversification(2, driver_speciation=0.8, root_drivers=1), age=3.0, seed=3)
    tree = res.tree
    assert tree.root.time == 0.0 and len(tree.root.children) == 2
    for node in tree.internal_nodes():
        assert len(node.children) == 2                           # bifurcating
    for leaf in tree.extant_leaves():
        assert abs(leaf.time - tree.total_age) < 1e-9            # extant tips at the present


def test_gene_div_n_tips_is_exact():
    res = z.simulate_gene_diversification(
        z.GeneDiversification(2, driver_speciation=0.5, root_drivers=1), n_tips=25, seed=5)
    assert len(res.tree.extant_leaves()) == 25


def test_gene_div_root_seeding():
    res = z.simulate_gene_diversification(z.GeneDiversification(3, root_drivers=2), n_tips=20, seed=1)
    assert set(res.node_drivers[res.tree.root]) == {0, 1}        # first two present at the root


# --------------------------------------------------------------------------- deterministic corners
def test_gene_div_no_gain_channel_never_appears():
    """No origination and no transfer, none at the root -> a driver can never appear."""
    m = z.GeneDiversification(2, origination=0.0, transfer=0.0, root_drivers=0)
    res = z.simulate_gene_diversification(m, n_tips=40, seed=2)
    assert res.tip_prevalence() == [0.0, 0.0]


def test_gene_div_no_loss_root_driver_is_fixed():
    """Present at the root with no loss and no gain -> fixed in every lineage (prevalence 1)."""
    m = z.GeneDiversification(1, loss=0.0, origination=0.0, transfer=0.0, root_drivers=1)
    res = z.simulate_gene_diversification(m, n_tips=40, seed=2)
    assert res.tip_prevalence() == [1.0]
    assert all(0 in res.node_drivers[n] for n in res.tree.nodes())


# --------------------------------------------------------------------------- the diversification signal
_BASE = dict(lambda0=1.0, mu0=0.2, loss=0.3, origination=0.05, transfer=0.6, root_drivers=1)


def test_speciation_driver_biases_tips_toward_it():
    """A key innovation (positive speciation effect) is over-represented among the standing tips."""
    p_pos = _mean_prevalence(driver_speciation=1.2, driver_extinction=0.0, **_BASE)
    p_neu = _mean_prevalence(driver_speciation=0.0, driver_extinction=0.0, **_BASE)
    assert p_pos > p_neu + 0.1


def test_extinction_driver_biases_tips_away():
    """A driver that raises extinction is under-represented among the standing tips."""
    p_ext = _mean_prevalence(driver_speciation=0.0, driver_extinction=1.2, **_BASE)
    p_neu = _mean_prevalence(driver_speciation=0.0, driver_extinction=0.0, **_BASE)
    assert p_ext < p_neu - 0.1


def test_transfer_spreads_a_driver():
    """Frequency-dependent HGT raises a neutral driver's prevalence vs no transfer."""
    common = dict(lambda0=1.0, mu0=0.2, loss=0.3, origination=0.05,
                  driver_speciation=0.0, root_drivers=1)
    p_tr = _mean_prevalence(transfer=1.5, **common)
    p_notr = _mean_prevalence(transfer=0.0, **common)
    assert p_tr > p_notr + 0.1


# --------------------------------------------------------------------------- result views
def test_gene_div_to_tsv_shape():
    res = z.simulate_gene_diversification(
        z.GeneDiversification(3, root_drivers=1), n_tips=15, seed=4)
    lines = res.to_tsv(nodes="all").strip().splitlines()
    assert lines[0] == "node\tD0\tD1\tD2"
    assert len(lines) == 1 + len(res.tree.nodes())
    for row in lines[1:]:
        cells = row.split("\t")[1:]
        assert all(c in ("0", "1") for c in cells)
