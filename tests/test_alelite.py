"""ALElite undated-DTL likelihood: closed-form oracles + ZOMBI2 end-to-end.

The undated model has exact limits that pin its speciation / speciation-loss / loss /
extinction machinery, so we check the DP against hand-derived probabilities rather than
against a reference binary (which may not be installed):

* a gene tree perfectly matching a k-tip species subtree, with d=t=0, has probability
  ``pS^(2k-1)`` — one slot per speciation and per tip sampling, nothing to lose;
* a gene present in only one of two sister species has probability ``l/(1+l)^3`` — root
  speciates, the present copy is sampled, the absent copy is lost.
"""

import math

import pytest

from zombi2 import Yule, read_newick, simulate_genomes, simulate_species_tree
from zombi2.alelite import (
    DatedDTL,
    GeneTree,
    SpeciesTree,
    UndatedDTL,
    dated_extinction,
    dated_joint_loglik,
    dated_loglik,
    undated_loglik,
)


def _ps(loss: float) -> float:
    return 1.0 / (1.0 + loss)  # d = t = 0 ⇒ denom = 1 + l


def _bd_u(s: float, d: float, lo: float) -> float:
    """Birth-death extinction probability over duration s (complete sampling, tau=0)."""
    if abs(d - lo) < 1e-12:
        return lo * s / (1 + lo * s)
    e = math.exp((d - lo) * s)
    return lo * (e - 1) / (d * e - lo)


def _bd_p1(s: float, d: float, lo: float) -> float:
    """Birth-death probability of exactly one sampled descendant over duration s (tau=0)."""
    if abs(d - lo) < 1e-12:
        return 1.0 / (1 + lo * s) ** 2
    e = math.exp(-(d - lo) * s)
    return (d - lo) ** 2 * e / (d - lo * e) ** 2


# --------------------------------------------------------------------------- oracles

@pytest.mark.parametrize("loss", [0.0, 0.3, 1.0, 2.5])
@pytest.mark.parametrize(
    "sp_newick, gt_newick, k",
    [
        ("(A:1,B:1)root;", "(A|1,B|2)r;", 2),
        ("((A:1,B:1)i1:1,C:2)root;", "((A|1,B|2)x,C|3)r;", 3),
        ("((A:1,B:1)i1:1,(C:1,D:1)i2:1)root;", "((A|1,B|2)x,(C|3,D|4)y)r;", 4),
    ],
)
def test_matching_tree_is_pS_to_the_2k_minus_1(sp_newick, gt_newick, k, loss):
    """d=t=0, gene tree == species tree ⇒ P = pS^(2k-1)."""
    sp = SpeciesTree.from_tree(read_newick(sp_newick))
    gt = GeneTree.from_newick(gt_newick)
    ll = undated_loglik(gt, sp, UndatedDTL(dup=0.0, transfer=0.0, loss=loss), origination="root")
    assert math.isclose(ll, (2 * k - 1) * math.log(_ps(loss)), rel_tol=1e-9, abs_tol=1e-12)


@pytest.mark.parametrize("loss", [0.2, 0.7, 1.5])
def test_gene_present_in_one_of_two_species(loss):
    """d=t=0, species ((A,B)) but the gene is only in A ⇒ P = l/(1+l)^3."""
    sp = SpeciesTree.from_tree(read_newick("(A:1,B:1)root;"))
    gt = GeneTree.from_newick("A|1;")  # single tip in species A
    ll = undated_loglik(gt, sp, UndatedDTL(dup=0.0, transfer=0.0, loss=loss), origination="root")
    expected = loss / (1.0 + loss) ** 3
    assert math.isclose(ll, math.log(expected), rel_tol=1e-9, abs_tol=1e-12)


def test_zero_rates_matching_tree_is_certain():
    """With no DTL at all a perfectly matching gene tree has probability 1 (log-lik 0)."""
    sp = SpeciesTree.from_tree(read_newick("((A:1,B:1)i1:1,C:2)root;"))
    gt = GeneTree.from_newick("((A|1,B|2)x,C|3)r;")
    ll = undated_loglik(gt, sp, UndatedDTL(0.0, 0.0, 0.0), origination="root")
    assert math.isclose(ll, 0.0, abs_tol=1e-12)


# --------------------------------------------------------------------------- properties

def test_it_is_a_normalised_probability_not_a_density():
    """Every reconciliation-summed likelihood is a probability in (0, 1]."""
    sp = SpeciesTree.from_tree(read_newick("((A:1,B:1)i1:1,(C:1,D:1)i2:1)root;"))
    gt = GeneTree.from_newick("((A|1,B|2)x,(C|3,D|4)y)r;")
    for model in [UndatedDTL(0.1, 0.05, 0.2), UndatedDTL(0.5, 0.5, 0.5), UndatedDTL(0.0, 0.3, 0.1)]:
        ll = undated_loglik(gt, sp, model)
        assert ll <= 1e-12  # log P <= 0
        assert math.isfinite(ll)


def test_transfer_makes_a_discordant_gene_possible():
    """A gene tree that groups A with C (crossing the species split) is impossible without
    transfer (log-lik -inf) but has positive probability once tau > 0."""
    sp = SpeciesTree.from_tree(read_newick("((A:1,B:1)i1:1,(C:1,D:1)i2:1)root;"))
    gt = GeneTree.from_newick("((A|1,C|2)x,(B|3,D|4)y)r;")  # A+C, B+D: discordant with species
    no_t = undated_loglik(gt, sp, UndatedDTL(dup=0.1, transfer=0.0, loss=0.1), origination="uniform")
    with_t = undated_loglik(gt, sp, UndatedDTL(dup=0.1, transfer=0.3, loss=0.1), origination="uniform")
    assert no_t == -math.inf or no_t < with_t
    assert math.isfinite(with_t)


# --------------------------------------------------------------------------- ZOMBI2 end-to-end

def test_zombi_reconciled_gene_tree_scores_finite():
    """Every extant ZOMBI2 family reconciled onto its species tree gets a finite log-lik <= 0."""
    tree = simulate_species_tree(Yule(1.0), n_tips=6, age=2.0, seed=7)
    g = simulate_genomes(tree, duplication=0.1, transfer=0.05, loss=0.15,
                         origination=0.0, initial_size=8, seed=7)
    sp = SpeciesTree.from_tree(tree)
    model = UndatedDTL(dup=0.1, transfer=0.05, loss=0.15)
    scored = 0
    for recon in g.reconciliations().values():
        if recon.extant is None:
            continue
        gt = GeneTree.from_reconciliation(recon)  # consumes .extant — the observable tree
        ll = undated_loglik(gt, sp, model, origination="root")
        assert math.isfinite(ll) and ll <= 1e-9
        scored += 1
    assert scored > 0


def test_complete_tree_is_rejected():
    """Feeding the complete tree (which carries LOSS tips) is refused with a clear error."""
    tree = simulate_species_tree(Yule(1.0), n_tips=5, age=2.0, seed=1)
    g = simulate_genomes(tree, duplication=0.1, transfer=0.0, loss=0.4,
                         origination=0.0, initial_size=8, seed=1)
    for recon in g.reconciliations().values():
        if recon.complete and "LOSS" in recon.complete:
            with pytest.raises(ValueError, match="EXTANT"):
                GeneTree.from_newick(recon.complete)
            break
    else:
        pytest.skip("no losses were simulated in this scenario")


def test_from_reconciliation_matches_extant_newick():
    """The convenience constructor uses exactly the extant tree."""
    tree = simulate_species_tree(Yule(1.0), n_tips=6, age=2.0, seed=7)
    g = simulate_genomes(tree, duplication=0.1, transfer=0.05, loss=0.15,
                         origination=0.0, initial_size=8, seed=7)
    for recon in g.reconciliations().values():
        if recon.extant is not None:
            a = GeneTree.from_reconciliation(recon)
            b = GeneTree.from_newick(recon.extant)
            assert a.n == b.n and a.species_set() == b.species_set()
            return
    pytest.fail("no extant family produced")


# ==================================================================== DATED engine

def test_dated_extinction_matches_birth_death():
    """With tau=0 the coupled extinction ODE reduces to the birth-death closed form."""
    sp = SpeciesTree.from_tree(read_newick("(A:1,B:1)root;"))
    d, lo = 0.6, 0.4
    E = dated_extinction(sp, DatedDTL(dup=d, transfer=0.0, loss=lo), n_steps=1600)
    assert math.isclose(E[sp.leaf_index["A"]], _bd_u(1.0, d, lo), abs_tol=1e-4)


def test_dated_zero_rates_matching_tree_is_certain():
    sp = SpeciesTree.from_tree(read_newick("((A:1,B:1)i1:1,C:2)root;"))
    gt = GeneTree.from_newick("((A|1,B|2)x,C|3)r;")
    ll = dated_loglik(gt, sp, DatedDTL(0.0, 0.0, 0.0), origination="root", n_steps=16)
    assert math.isclose(ll, 0.0, abs_tol=1e-12)


def test_dated_single_gene_matches_birth_death_closed_form():
    """Gene present in only one of two sisters, tau=0 ⇒ P = p1(s)·E_sister (birth-death)."""
    sp = SpeciesTree.from_tree(read_newick("(A:1,B:1)root;"))
    gt = GeneTree.from_newick("A|1;")
    d, lo = 0.5, 0.3
    expected = _bd_p1(1.0, d, lo) * _bd_u(1.0, d, lo)
    ll = dated_loglik(gt, sp, DatedDTL(dup=d, transfer=0.0, loss=lo), origination="root", n_steps=3200)
    assert math.isclose(math.exp(ll), expected, rel_tol=2e-3)


def test_dated_matching_pair_matches_birth_death_closed_form():
    """Two sisters each with one gene, tau=0 ⇒ P = p1(s)^2 (root speciation, both survive)."""
    sp = SpeciesTree.from_tree(read_newick("(A:1,B:1)root;"))
    gt = GeneTree.from_newick("(A|1,B|2)r;")
    d, lo = 0.5, 0.3
    ll = dated_loglik(gt, sp, DatedDTL(dup=d, transfer=0.0, loss=lo), origination="root", n_steps=3200)
    assert math.isclose(math.exp(ll), _bd_p1(1.0, d, lo) ** 2, rel_tol=2e-3)


def test_dated_likelihood_converges_in_resolution():
    """Refining the time grid changes the log-lik by a vanishing amount (first-order)."""
    sp = SpeciesTree.from_tree(read_newick("((A:1,B:1)i1:1,(C:1,D:1)i2:1)root;"))
    gt = GeneTree.from_newick("((A|1,B|2)x,(C|3,D|4)y)r;")
    model = DatedDTL(dup=0.2, transfer=0.15, loss=0.25)
    lls = [dated_loglik(gt, sp, model, n_steps=ns) for ns in (100, 400, 1600)]
    assert abs(lls[1] - lls[0]) > abs(lls[2] - lls[1])  # converging
    assert abs(lls[2] - lls[1]) < 1e-2


def test_dated_transfer_enables_discordant_gene():
    sp = SpeciesTree.from_tree(read_newick("((A:1,B:1)i1:1,(C:1,D:1)i2:1)root;"))
    gt = GeneTree.from_newick("((A|1,C|2)x,(B|3,D|4)y)r;")  # A+C, B+D crosses the species split
    no_t = dated_loglik(gt, sp, DatedDTL(dup=0.1, transfer=0.0, loss=0.1),
                        origination="uniform", n_steps=200)
    with_t = dated_loglik(gt, sp, DatedDTL(dup=0.1, transfer=0.4, loss=0.1),
                          origination="uniform", n_steps=200)
    assert no_t == -math.inf or no_t < with_t
    assert math.isfinite(with_t)


def test_dated_inject_recover_prefers_true_rates():
    """Inject-recover: the dated joint log-lik of ZOMBI2 families is higher at the true rates
    than at rates off by ~2.5x either way — the model is faithful to the generative process.
    (A Yule tree has no species extinction, so ZOMBI2's contemporaneous transfers match the
    dated model exactly.)"""
    true = dict(dup=0.30, transfer=0.20, loss=0.45)
    tree = simulate_species_tree(Yule(2.0), n_tips=7, age=1.0, seed=11)
    g = simulate_genomes(tree, duplication=true["dup"], transfer=true["transfer"],
                         loss=true["loss"], origination=0.0, initial_size=30, seed=11)
    sp = SpeciesTree.from_tree(tree)
    trees, n_extinct = [], 0
    for r in g.reconciliations().values():
        if r.extant is None:
            n_extinct += 1
        else:
            trees.append(GeneTree.from_reconciliation(r))
    assert len(trees) >= 10

    def joint(**rates):
        return dated_joint_loglik(trees, sp, DatedDTL(**rates), n_extinct=n_extinct, n_steps=20)

    base = joint(**true)
    for which in ("dup", "transfer", "loss"):
        low, high = dict(true), dict(true)
        low[which] = true[which] * 0.4
        high[which] = true[which] * 2.5
        assert base > joint(**low), f"{which}: truth not preferred over low"
        assert base > joint(**high), f"{which}: truth not preferred over high"


def test_dated_rust_matches_python():
    """The compiled Rust kernel reproduces the Python reference (bit-for-bit up to summation)."""
    from zombi2.alelite import _rust
    if not _rust.available():
        pytest.skip("zombi2_core extension not built")
    tree = simulate_species_tree(Yule(2.0), n_tips=8, age=1.0, seed=3)
    g = simulate_genomes(tree, duplication=0.3, transfer=0.2, loss=0.4,
                         origination=0.0, initial_size=25, seed=3)
    sp = SpeciesTree.from_tree(tree)
    trees, n_extinct = [], 0
    for r in g.reconciliations().values():
        if r.extant is None:
            n_extinct += 1
        else:
            trees.append(GeneTree.from_reconciliation(r))
    for origination in ("root", "uniform"):
        for rates in [(0.3, 0.2, 0.4), (0.1, 0.05, 0.6), (0.5, 0.4, 0.3)]:
            m = DatedDTL(*rates)
            py = dated_joint_loglik(trees, sp, m, origination=origination,
                                    n_extinct=n_extinct, n_steps=30, backend="python")
            ru = dated_joint_loglik(trees, sp, m, origination=origination,
                                    n_extinct=n_extinct, n_steps=30, backend="rust")
            assert abs(py - ru) < 1e-7 * max(1.0, abs(py)), (origination, rates, py, ru)


def test_dated_zombi_families_score_finite():
    tree = simulate_species_tree(Yule(1.0), n_tips=6, age=2.0, seed=7)
    g = simulate_genomes(tree, duplication=0.1, transfer=0.05, loss=0.15,
                         origination=0.0, initial_size=8, seed=7)
    sp = SpeciesTree.from_tree(tree)
    model = DatedDTL(dup=0.1, transfer=0.05, loss=0.15)
    scored = 0
    for recon in g.reconciliations().values():
        if recon.extant is None:
            continue
        gt = GeneTree.from_reconciliation(recon)
        ll = dated_loglik(gt, sp, model, origination="root", n_steps=60)
        assert math.isfinite(ll) and ll <= 1e-6
        scored += 1
    assert scored > 0
