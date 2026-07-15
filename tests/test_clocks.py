"""The relaxed-molecular-clock family (``zombi2.clocks``).

Every clock turns a chronogram (timetree) into a phylogram (branch lengths in expected
substitutions per site) by drawing per-branch rate multipliers, sharing one :class:`Clock`
interface. This suite pins the defining property of each model — strict vs. relaxed,
uncorrelated vs. autocorrelated, and the branch-length dependence of white noise — plus the
common contract (reproducibility, valid Newick, integration with ``SequenceEvolution``).
"""

import numpy as np
import pytest

import zombi2 as z
from zombi2 import (
    AutocorrelatedLogNormalClock,
    CIRClock,
    Clock,
    RateVariation,
    StrictClock,
    UncorrelatedGammaClock,
    UncorrelatedLogNormalClock,
    WhiteNoiseClock,
)


def _tree(n_tips=400, age=10.0, seed=7):
    return z.simulate_species_tree(z.BirthDeath(1.0, 0.1), n_tips=n_tips, age=age, seed=seed)


def _branches(tree):
    return [n for n in tree.nodes_preorder() if n.parent is not None]


def _mean_rate(scaled, tree):
    branches = _branches(tree)
    tot_time = sum(n.branch_length() for n in branches)
    tot_subs = sum(scaled.branch_lengths[n] for n in branches)
    return tot_subs / tot_time


def _parent_child_corr(scaled, tree):
    par, chi = [], []
    for n in _branches(tree):
        if n.parent is not None and n.parent.parent is not None:
            par.append(scaled.branch_rate[n.parent])
            chi.append(scaled.branch_rate[n])
    return np.corrcoef(par, chi)[0, 1]


# every clock, with mean-1 (or unit-rate) settings, for the shared-contract tests
ALL_CLOCKS = {
    "strict": StrictClock(1.0),
    "ucln": UncorrelatedLogNormalClock(0.5),
    "ugam": UncorrelatedGammaClock(3.0),
    "white_noise": WhiteNoiseClock(0.5),
    "autocorr_lognormal": AutocorrelatedLogNormalClock(0.3),
    "cir": CIRClock(theta=1.0, sigma=0.4, mean=1.0),
    "rate_variation": RateVariation(bins=[0.5, 1.0, 2.0], switch_rate=1.0),
}
# the clocks built to have mean rate 1 (autocorrelated-lognormal drifts; strict is exactly 1)
UNIT_MEAN_CLOCKS = ["ucln", "ugam", "white_noise", "cir"]


# --- shared Clock contract ------------------------------------------------------------------

@pytest.mark.parametrize("name", sorted(ALL_CLOCKS))
def test_scale_is_a_valid_phylogram(name):
    """Every clock rescales all branches and emits parseable Newick."""
    tree = _tree(n_tips=40, age=6.0, seed=1)
    scaled = ALL_CLOCKS[name].scale(tree, seed=2)
    assert set(scaled.branch_lengths) == set(tree.nodes_preorder())
    assert all(v >= 0 for v in scaled.branch_lengths.values())
    back = z.read_newick(scaled.to_newick())
    assert {l.name for l in back.leaves()} == {l.name for l in tree.leaves()}


@pytest.mark.parametrize("name", sorted(ALL_CLOCKS))
def test_reproducible_given_seed(name):
    tree = _tree(n_tips=40, age=6.0, seed=1)
    clock = ALL_CLOCKS[name]
    assert clock.scale(tree, seed=99).to_newick() == clock.scale(tree, seed=99).to_newick()


@pytest.mark.parametrize("name", sorted(ALL_CLOCKS))
def test_applies_to_gene_trees_too(name):
    """A clock scales any Tree — including a gene tree loaded from Newick."""
    tree = _tree(n_tips=12, age=4.0, seed=9)
    g = z.simulate_genomes(tree, z.Rates(duplication=0.2, transfer=0.1, loss=0.2,
                                               origination=0.5), initial_families=12, seed=9)
    fam = next(f for f, (_, e) in g.gene_trees().items() if e)
    gene_tree = z.read_newick(g.gene_trees()[fam][1])
    scaled = ALL_CLOCKS[name].scale(gene_tree, seed=1)
    assert len(scaled.branch_lengths) == len(gene_tree.nodes())
    assert z.read_newick(scaled.to_newick()).leaves()


@pytest.mark.parametrize("name", UNIT_MEAN_CLOCKS)
def test_unit_mean_clocks_average_to_one(name):
    """A mean-1 relaxed clock leaves the tree's total length ≈ unchanged (∑ rate·time ≈ ∑ time)."""
    tree = _tree()
    assert _mean_rate(ALL_CLOCKS[name].scale(tree, seed=3), tree) == pytest.approx(1.0, abs=0.2)


def test_clock_base_class_is_abstract():
    """The base Clock has no rate rule of its own."""
    with pytest.raises(NotImplementedError):
        Clock().scale(_tree(n_tips=8, age=3.0, seed=1), seed=1)


# --- strict clock ---------------------------------------------------------------------------

def test_strict_clock_of_rate_one_is_the_chronogram():
    tree = _tree(n_tips=30, age=5.0, seed=4)
    assert StrictClock(1.0).scale(tree, seed=1).to_newick() == tree.to_newick()


def test_strict_clock_uniformly_stretches():
    tree = _tree(n_tips=30, age=5.0, seed=4)
    scaled = StrictClock(2.5).scale(tree, seed=1)
    for n in _branches(tree):
        assert scaled.branch_lengths[n] == pytest.approx(2.5 * n.branch_length())
    assert all(r == 2.5 for r in scaled.branch_rate.values())


# --- uncorrelated vs. autocorrelated (the defining split) -----------------------------------

@pytest.mark.parametrize("clock", [
    UncorrelatedLogNormalClock(0.5),
    UncorrelatedGammaClock(3.0),
    WhiteNoiseClock(0.5),
])
def test_uncorrelated_clocks_have_near_zero_parent_child_correlation(clock):
    tree = _tree()
    assert abs(_parent_child_corr(clock.scale(tree, seed=3), tree)) < 0.2


@pytest.mark.parametrize("clock", [
    AutocorrelatedLogNormalClock(0.3),
    CIRClock(theta=1.0, sigma=0.4, mean=1.0),
])
def test_autocorrelated_clocks_have_positive_parent_child_correlation(clock):
    tree = _tree()
    assert _parent_child_corr(clock.scale(tree, seed=3), tree) > 0.4


def test_relaxed_clocks_actually_vary_rates():
    tree = _tree(n_tips=60, age=8.0, seed=2)
    for clock in [UncorrelatedLogNormalClock(0.5), UncorrelatedGammaClock(2.0),
                  WhiteNoiseClock(0.5), AutocorrelatedLogNormalClock(0.4),
                  CIRClock(theta=1.0, sigma=0.5)]:
        rates = clock.scale(tree, seed=1).branch_rate.values()
        assert len({round(r, 6) for r in rates}) > 1


# --- lognormal / gamma: sigma=0 (or huge shape) collapses to the strict clock ----------------

def test_uncorrelated_lognormal_sigma_zero_is_strict():
    tree = _tree(n_tips=30, age=5.0, seed=4)
    assert UncorrelatedLogNormalClock(0.0, mean=1.3).scale(tree, seed=1).to_newick() == \
        StrictClock(1.3).scale(tree, seed=1).to_newick()


def test_positive_rates_for_lognormal_gamma_cir():
    tree = _tree()
    for clock in [UncorrelatedLogNormalClock(0.8), UncorrelatedGammaClock(1.5),
                  AutocorrelatedLogNormalClock(0.5), CIRClock(theta=0.5, sigma=0.5)]:
        assert min(clock.scale(tree, seed=3).branch_rate.values()) > 0.0


# --- white noise: variance of the branch rate falls off with branch length ------------------

def test_white_noise_variance_scales_inversely_with_branch_length():
    tree = _tree()
    branches = _branches(tree)
    scaled = WhiteNoiseClock(0.6).scale(tree, seed=5)
    durs = np.array([n.branch_length() for n in branches])
    rates = np.array([scaled.branch_rate[n] for n in branches])
    med = np.median(durs)
    assert rates[durs < med].var() > rates[durs >= med].var()


# --- autocorrelated lognormal: the branch_sigma clock, now a first-class model ---------------

def test_autocorrelated_lognormal_sigma_zero_is_strict_at_root_rate():
    tree = _tree(n_tips=30, age=5.0, seed=4)
    scaled = AutocorrelatedLogNormalClock(0.0, root_rate=1.7).scale(tree, seed=1)
    assert all(r == 1.7 for r in scaled.branch_rate.values())


# --- CIR: mean reversion and within-branch variation ----------------------------------------

def test_cir_reverts_to_its_long_run_mean():
    tree = _tree()
    assert _mean_rate(CIRClock(theta=2.0, sigma=0.4, mean=1.5).scale(tree, seed=3), tree) == \
        pytest.approx(1.5, abs=0.25)


def test_cir_varies_rate_within_a_branch():
    """CIR is a diffusion — a branch longer than max_step is split into several rate segments."""
    tree = _tree(n_tips=30, age=8.0, seed=4)
    scaled = CIRClock(theta=1.0, sigma=0.4, max_step=0.05).scale(tree, seed=5)
    assert max(len(scaled.segments[n]) for n in _branches(tree)) > 1


def test_cir_sigma_zero_started_at_mean_is_strict():
    tree = _tree(n_tips=30, age=5.0, seed=4)
    scaled = CIRClock(theta=1.0, sigma=0.0, mean=1.0).scale(tree, seed=1)
    for n in _branches(tree):
        assert scaled.branch_rate[n] == pytest.approx(1.0)


# --- integration with SequenceEvolution (the shared lineage clock) --------------------------

def _genomes():
    tree = z.simulate_species_tree(z.BirthDeath(1.0, 0.2), n_tips=12, age=3.0, seed=1)
    g = z.simulate_genomes(tree, z.PerGenomeRates(duplication=0.3, transfer=0.15, loss=0.2,
                                                  origination=0.6), initial_families=12, seed=2)
    return tree, g


@pytest.mark.parametrize("clock", [
    StrictClock(1.2),
    UncorrelatedLogNormalClock(0.5),
    UncorrelatedGammaClock(2.0),
    WhiteNoiseClock(0.5),
    CIRClock(theta=1.0, sigma=0.4),
])
def test_any_clock_drives_the_lineage_clock(clock):
    """SequenceEvolution accepts any Clock as the shared lineage clock."""
    _, g = _genomes()
    ph = z.SequenceEvolution(lineage=clock, family_speed=z.LogNormal(0.0, 0.4)).scale(g, seed=5)
    extant = [e for e in ph.extant.values() if e]
    assert extant and all(e.endswith(";") for e in extant)


def test_branch_sigma_equals_autocorrelated_lognormal_lineage():
    """The branch_sigma convenience is exactly lineage=AutocorrelatedLogNormalClock(sigma)."""
    _, g = _genomes()
    a = z.SequenceEvolution(branch_sigma=0.5, family_speed=z.LogNormal(0.0, 0.4)).scale(g, seed=9)
    b = z.SequenceEvolution(lineage=AutocorrelatedLogNormalClock(0.5, root_rate=1.0),
                            family_speed=z.LogNormal(0.0, 0.4)).scale(g, seed=9)
    assert a.extant == b.extant
    assert a.branch_rate == b.branch_rate


def test_strict_lineage_clock_reproduces_the_chronogram():
    """A unit strict lineage clock with unit family speed leaves gene trees in time units."""
    _, g = _genomes()
    ph = z.SequenceEvolution(lineage=StrictClock(1.0), family_speed=1.0).scale(g, seed=3)
    chrono = g.gene_trees()
    checked = 0
    for fam, (_, extant) in chrono.items():
        if extant is None:
            continue
        checked += 1
        assert ph.extant[fam] == extant
    assert checked > 0


def test_non_clock_lineage_rejected():
    with pytest.raises(TypeError):
        z.SequenceEvolution(lineage="not-a-clock")


# --- bad parameters -------------------------------------------------------------------------

@pytest.mark.parametrize("factory", [
    lambda: StrictClock(0.0),
    lambda: StrictClock(-1.0),
    lambda: UncorrelatedLogNormalClock(-0.1),
    lambda: UncorrelatedLogNormalClock(0.5, mean=0.0),
    lambda: UncorrelatedGammaClock(0.0),
    lambda: WhiteNoiseClock(-0.1),
    lambda: AutocorrelatedLogNormalClock(-0.1),
    lambda: AutocorrelatedLogNormalClock(0.5, root_rate=0.0),
    lambda: CIRClock(theta=-1.0, sigma=0.4),
    lambda: CIRClock(theta=1.0, sigma=-0.4),
    lambda: CIRClock(theta=1.0, sigma=0.4, max_step=0.0),
])
def test_bad_parameters_rejected(factory):
    with pytest.raises(ValueError):
        factory()
