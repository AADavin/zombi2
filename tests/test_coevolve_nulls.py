"""Null models of coevolution — the decoupled nulls that cut a coupling's arrow.

Every ``zombi2 coevolve`` edge exposes a matched null via ``.null(kind=...)``. The tests below
are the mirror image of the coevolution validation suite: instead of *injecting and recovering* a
coupling, they inject a coupling, cut its arrow, and confirm the signal is **gone** while the
target keeps its marginal variance (the honest-null property). This file covers the
``traits:species`` flagship — the ``CID`` (Character-Independent Diversification) and ``neutral``
nulls. See ``docs/guide/coevolution_nulls.md``.
"""

import numpy as np
import pytest

import zombi2 as z
from zombi2.coevolve import (
    CID, GeneDiversification, GeneConditionedTrait, CladogeneticGenome, TraitGeneCoupling,
    simulate_gene_diversification, simulate_cladogenetic_genome,
    simulate_gene_conditioned_trait, simulate_trait_linked_genomes,
)
from zombi2.coevolve.cladogenetic_genome import _branch_count_and_length


# --------------------------------------------------------------------------- helpers
def _tip_state_fraction(model, *, age, reps, seed=0, state=1):
    """Fraction of extant tips in the observed ``state``, pooled over ``reps`` trees."""
    rng = np.random.default_rng(seed)
    hits = tot = 0
    for _ in range(reps):
        res = z.simulate_sse(model, age=age, rng=rng)
        vals = list(res.labeled_values().values())
        hits += sum(1 for v in vals if v == state)
        tot += len(vals)
    return hits / tot


def _observed_and_hidden(model, *, age, reps, seed):
    """Pooled (observed state, is-fast-hidden) over the tips of ``reps`` CID trees."""
    rng = np.random.default_rng(seed)
    obs, fast = [], []
    for _ in range(reps):
        res = z.simulate_sse(model, age=age, rng=rng)
        for leaf in res.tree.extant_leaves():
            obs.append(res.labeled_values()[leaf])
            fast.append(1 if res.full_label(res.node_values[leaf])[1] == 1 else 0)
    return np.array(obs), np.array(fast)


# --------------------------------------------------------------------------- CID construction
def test_cid_two_carries_the_bisse_rate_pairs_exactly():
    """CID-2 built from a BiSSE keeps its two observed rate pairs, reassigned to hidden classes."""
    alt = z.BiSSE(lambda0=1, lambda1=3, mu0=0.2, mu1=0.4, q01=0.1, q10=0.2)
    cid = alt.null(kind="cid", n_hidden=2)
    assert isinstance(cid, CID)
    assert cid._class_rates == [(1.0, 0.2), (3.0, 0.4)]     # (lambda, mu) per hidden class
    assert cid._q == (0.1, 0.2)                             # observed transitions preserved
    # within each hidden class the two OBSERVED states share the rate (character-independent)
    assert np.isclose(cid.lambdas[0], cid.lambdas[2])       # obs 0 vs obs 1, hidden 0
    assert np.isclose(cid.lambdas[1], cid.lambdas[3])       # obs 0 vs obs 1, hidden 1


def test_cid_four_spreads_four_regimes_over_the_rate_range():
    alt = z.BiSSE(lambda0=0.5, lambda1=2.5, mu0=0.1, mu1=0.1, q01=0.2, q10=0.2)
    cid = alt.null(kind="cid", n_hidden=4)
    assert isinstance(cid, CID) and cid._H == 4
    lams = sorted(r[0] for r in cid._class_rates)
    assert lams[0] == pytest.approx(0.5) and lams[-1] == pytest.approx(2.5)   # brackets the range


def test_cid_factories():
    two = CID.two(lambda_slow=0.5, lambda_fast=2.0, mu=0.2, switch=0.15)
    four = CID.four([0.4, 0.8, 1.6, 2.4], switch=0.1)
    assert two._H == 2 and two.k == 4
    assert four._H == 4 and four.k == 8
    assert repr(two) == "CID(classes=2)"


def test_cid_requires_an_evolving_observed_character():
    """A frozen observed character (q01=q10=0) is trivially independent — reject it."""
    with pytest.raises(ValueError, match="still evolve"):
        CID([(1.0, 0.2), (3.0, 0.2)], hidden_transition=0.1, q01=0.0, q10=0.0)


def test_cid_four_needs_four_classes():
    with pytest.raises(ValueError, match="exactly four"):
        CID.four([0.4, 0.8, 1.6], switch=0.1)


# --------------------------------------------------------------------------- neutral null
def test_neutral_null_makes_rates_state_independent():
    """The neutral null averages the per-state rates and leaves the anagenetic Q untouched."""
    alt = z.BiSSE(lambda0=1, lambda1=3, mu0=0.2, mu1=0.4, q01=0.1, q10=0.3)
    null = alt.null(kind="neutral")
    assert np.allclose(null.lambdas, 2.0)          # (1 + 3) / 2
    assert np.allclose(null.mus, 0.3)              # (0.2 + 0.4) / 2
    assert np.allclose(null.Q, alt.Q)              # the character still evolves as before


def test_neutral_null_generalises_to_musse():
    Q = [[0, 0.2, 0.2], [0.2, 0, 0.2], [0.2, 0.2, 0]]
    alt = z.MuSSE(birth=[1, 1, 3], death=[0.2, 0.2, 0.2], Q=Q)
    null = alt.null(kind="neutral")
    assert np.allclose(null.lambdas, 5 / 3)        # mean of [1, 1, 3]
    assert np.allclose(null.mus, 0.2)


# --------------------------------------------------------------------------- the honest-null property
def test_neutral_null_removes_the_diversification_bias():
    """Coupled BiSSE biases the standing tips toward the fast state; its neutral null does not."""
    alt = z.BiSSE(1, 3, 0.3, 0.3, 0.3, 0.3)                 # state 1 speciates 3x faster
    frac_alt = _tip_state_fraction(alt, age=1.8, reps=120)
    frac_null = _tip_state_fraction(alt.null(kind="neutral"), age=1.8, reps=120)
    assert frac_alt > 0.7                                   # the coupling shows up
    assert abs(frac_null - 0.5) < 0.1                       # ... and is gone in the null


def test_cid_null_hidden_drives_diversification_not_the_observed_character():
    """Under the CID null the tree gets real fast/slow clades from a HIDDEN class, while the
    observed character — which a raw BiSSE fit would test — stays neutral. Mirrors the HiSSE
    validation test (fast observed mixing q=0.3 keeps the observed marginal near 0.5)."""
    alt = z.BiSSE(1, 3, 0.2, 0.2, 0.3, 0.3)
    cid = alt.null(kind="cid", n_hidden=2)
    obs, fast = _observed_and_hidden(cid, age=1.5, reps=100, seed=0)
    assert abs(obs.mean() - 0.5) < 0.1                      # observed character is not biased
    assert fast.mean() > 0.6                                # the fast hidden class dominates tips


# --------------------------------------------------------------------------- plumbing
def test_cid_null_is_reproducible():
    alt = z.BiSSE(1, 3, 0.2, 0.2, 0.3, 0.3)
    cid = alt.null(kind="cid")
    a = z.simulate_sse(cid, age=1.5, seed=7).tree.to_newick()
    b = z.simulate_sse(cid, age=1.5, seed=7).tree.to_newick()
    assert a == b


def test_traits_species_has_no_timing_null():
    with pytest.raises(ValueError, match="no 'timing' null"):
        z.BiSSE(1, 3, 0.2, 0.2, 0.1, 0.1).null(kind="timing")


def test_unknown_null_kind_rejected():
    with pytest.raises(ValueError, match="unknown null kind"):
        z.BiSSE(1, 3, 0.2, 0.2, 0.1, 0.1).null(kind="bogus")


# =========================================================================== other edges
def _tree(seed=1, tips=40):
    return z.simulate_species_tree(z.BirthDeath(1.0, 0.3), n_tips=tips, age=5.0, seed=seed)


# --------------------------------------------------------------------------- genes:species
def test_gene_diversification_neutral_zeroes_the_rate_effect():
    gd = GeneDiversification(2, lambda0=1, mu0=0.2, driver_speciation=1.5, transfer=0.8,
                             root_drivers=1)
    n = gd.null("neutral")
    assert np.allclose(n.beta_lambda, 0) and np.allclose(n.beta_mu, 0)
    assert n.rates({0, 1}) == (1.0, 0.2) == n.rates(set())         # constant-rate, any gene content
    assert gd.rates({0, 1})[0] > 1.0                               # the coupling did scale λ


def _driver0_prevalence(model, *, reps=30, n_tips=70, seed=0):
    rng = np.random.default_rng(seed)
    ps = [simulate_gene_diversification(model, n_tips=n_tips, rng=rng).tip_prevalence()[0]
          for _ in range(reps)]
    return float(np.mean(ps))


def test_gene_diversification_neutral_removes_the_tip_bias():
    """A strong speciation driver is over-represented among tips under coupling; cutting the arrow
    (neutral) leaves it to plain HGT dynamics, so it is markedly less over-represented."""
    gd = GeneDiversification(1, lambda0=1, mu0=0.2, driver_speciation=1.6, transfer=0.6,
                             loss=0.3, root_drivers=1)
    prev_coupled = _driver0_prevalence(gd, seed=0)
    prev_null = _driver0_prevalence(gd.null("neutral"), seed=0)
    assert prev_coupled > prev_null + 0.1


def test_genes_species_cid_is_a_workflow_not_a_method():
    with pytest.raises(TypeError, match="workflow"):
        GeneDiversification(2, root_drivers=1).null("cid")


# --------------------------------------------------------------------------- genes:traits
def _carrier_gap(model, tree, *, reps=8, seed=0):
    """|mean(trait | modifier present) − mean(trait | modifier absent)|, pooled over reps."""
    rng = np.random.default_rng(seed)
    car, non = [], []
    for _ in range(reps):
        res = simulate_gene_conditioned_trait(tree, model, rng=rng)
        gp = res.gene_presence()
        for leaf, v in res.trait_values().items():
            (car if gp[leaf] else non).append(v)
    if not car or not non:
        return 0.0
    return abs(np.mean(car) - np.mean(non))


def test_gene_conditioned_trait_neutral_flattens_the_optimum():
    gc = GeneConditionedTrait(gene_gain=0.8, gene_loss=0.8, theta_absent=0, theta_present=6,
                              alpha=3.0)
    n = gc.null("neutral")
    assert n.theta_present == n.theta_absent == 0.0
    tree = _tree()
    gap_coupled = _carrier_gap(gc, tree, seed=1)
    gap_null = _carrier_gap(n, tree, seed=1)
    assert gap_coupled > 2.0            # carriers sit near 6, non-carriers near 0
    assert gap_null < 1.0               # ... and the gene tells you nothing about the trait now


# --------------------------------------------------------------------------- species:genes (timing)
def test_cladogenetic_neutral_drops_the_burst():
    cg = CladogeneticGenome(30, cladogenetic_loss=0.15, cladogenetic_gain=3.0)
    n = cg.null("neutral")
    assert n.cladogenetic_loss == 0.0 and n.cladogenetic_gain == 0.0
    assert n.loss == cg.loss and n.origination == cg.origination


def test_cladogenetic_timing_match_is_analytic():
    tree = _tree()
    cg = CladogeneticGenome(30, loss=0.05, origination=0.1, cladogenetic_loss=0.2,
                            cladogenetic_gain=2.0)
    tim = cg.null("timing", tree=tree)
    nb, L = _branch_count_and_length(tree)
    assert tim.cladogenetic_loss == 0.0 and tim.cladogenetic_gain == 0.0
    assert np.isclose(tim.loss, 0.05 + 0.2 * nb / L)
    assert np.isclose(tim.origination, 0.1 + 2.0 * nb / L)
    with pytest.raises(ValueError, match="needs the tree"):
        cg.null("timing")


def _mean_tip_genome(tree, model, *, reps=6, seed=0):
    rng = np.random.default_rng(seed)
    means = []
    for _ in range(reps):
        pres = simulate_cladogenetic_genome(tree, model, rng=rng).profile_matrix().presence()
        means.append(float(np.asarray(pres).sum(axis=0).mean()))
    return float(np.mean(means))


def test_cladogenetic_timing_preserves_the_marginal_turnover():
    """The timing null spreads the burst along branches — the punctuation *pattern* changes, but
    the marginal amount of gene content is matched in expectation."""
    tree = _tree()
    cg = CladogeneticGenome(30, loss=0.0, origination=0.0, cladogenetic_loss=0.12,
                            cladogenetic_gain=1.5)
    size_burst = _mean_tip_genome(tree, cg, seed=0)
    size_timing = _mean_tip_genome(tree, cg.null("timing", tree=tree), seed=0)
    assert abs(size_burst - size_timing) / size_burst < 0.25       # same marginal, ~in expectation


def test_species_genes_has_no_cid_null():
    with pytest.raises(ValueError, match="no 'cid' null"):
        CladogeneticGenome(30).null("cid")


# --------------------------------------------------------------------------- traits:genes
def test_trait_gene_coupling_neutral_cuts_both_channels():
    coup = TraitGeneCoupling.build(30, 0.4, weight=1.0, effect_loss=3.0, effect_gain=0.0,
                                   base_loss=0.5, transfer=1.0, state_values=[-1.0, 1.0], seed=1)
    n = coup.null("neutral")
    assert n.effect_loss == 0.0 and n.effect_gain == 0.0
    # still simulates and yields the full panel
    res = simulate_trait_linked_genomes(_tree(), z.Mk.equal_rates(2, 0.4), n, seed=2)
    assert res.profiles.presence().shape[0] == 30


def test_traits_genes_cid_is_a_workflow_not_a_method():
    coup = TraitGeneCoupling.build(10, 0.4, effect_loss=2.0)
    with pytest.raises(TypeError, match="workflow"):
        coup.null("cid")


# =========================================================================== CLI --null
from zombi2.cli import main


def _tree_file(tmp_path, seed=1, tips=25, age=5.0):
    p = tmp_path / "tree.nwk"
    p.write_text(z.prune(z.simulate_species_tree(z.BirthDeath(1.0, 0.3), n_tips=tips, age=age,
                                                 seed=seed)).to_newick() + "\n")
    return str(p)


def test_cli_traits_species_cid_writes_manifest(tmp_path):
    out = tmp_path / "o"
    rc = main(["coevolve", "--couple", "traits:species", "--lambda0", "1", "--lambda1", "3",
               "--tips", "80", "--seed", "1", "--null", "cid", "--hidden", "2", "-o", str(out)])
    assert rc == 0
    man = (out / "null_manifest.tsv").read_text()
    assert "null\tcid" in man and "hidden_classes\t2" in man


def test_cli_species_genes_timing_writes_output_and_manifest(tmp_path):
    out = tmp_path / "o"
    rc = main(["coevolve", "--couple", "species:genes", "-t", _tree_file(tmp_path),
               "--genome-size", "25", "--clado-gene-loss", "0.15", "--clado-gene-gain", "3",
               "--null", "timing", "-o", str(out)])
    assert rc == 0
    assert (out / "profiles.tsv").exists()
    assert "null\ttiming" in (out / "null_manifest.tsv").read_text()


def test_cli_species_traits_neutral(tmp_path):
    out = tmp_path / "o"
    rc = main(["coevolve", "--couple", "species:traits", "-t", _tree_file(tmp_path),
               "--sse-model", "bisse", "--q01", "0.1", "--q10", "0.1", "--clado-shift", "0.4",
               "--null", "neutral", "-o", str(out)])
    assert rc == 0
    assert (out / "null_manifest.tsv").exists()


def test_cli_null_declines_joint_model(tmp_path):
    with pytest.raises(SystemExit):
        main(["coevolve", "--couple", "traits:species", "--couple", "species:traits",
              "--tips", "40", "--null", "cid", "-o", str(tmp_path / "o")])


# --------------------------------------------------------------------------- CID workflows (CLI)
def test_cli_genes_species_cid_workflow(tmp_path):
    """genes:species --null cid: a driver-shaped tree + a NEUTRAL overlay genome + the drivers
    withheld as ground-truth."""
    out = tmp_path / "o"
    rc = main(["coevolve", "--couple", "genes:species", "--drivers", "2", "--root-drivers", "1",
               "--driver-speciation", "1.4", "--tips", "40", "--seed", "1", "--null", "cid",
               "--genome-size", "25", "-o", str(out)])
    assert rc == 0
    assert (out / "profiles.tsv").exists()                       # the observed neutral genome
    assert (out / "drivers_ground_truth.tsv").exists()           # the hidden drivers
    assert "null\tcid" in (out / "null_manifest.tsv").read_text()


def test_cli_genes_traits_cid_workflow(tmp_path):
    out = tmp_path / "o"
    rc = main(["coevolve", "--couple", "genes:traits", "-t", _tree_file(tmp_path),
               "--modifier-gain", "0.6", "--modifier-loss", "0.6", "--theta-present", "5",
               "--seed", "2", "--null", "cid", "--genome-size", "20", "-o", str(out)])
    assert rc == 0
    assert (out / "profiles.tsv").exists() and (out / "traits.tsv").exists()
    assert (out / "modifier_ground_truth.tsv").exists()          # the hidden modifier


def _read_presence(path):
    lines = path.read_text().splitlines()
    tips = lines[0].split("\t")[1:]
    tot = {t: 0 for t in tips}
    nfam = 0
    for row in lines[1:]:
        nfam += 1
        for t, v in zip(tips, row.split("\t")[1:]):
            tot[t] += int(v)
    return {t: tot[t] / nfam for t in tips}                      # per-tip mean panel occupancy


def _read_trait(path):
    return {r.split("\t")[0]: float(r.split("\t")[1]) for r in path.read_text().splitlines()[1:]}


def _state_gap(occ, trait):
    """|mean occupancy in state 1 − mean occupancy in state 0| over the tips."""
    hi = [occ[t] for t in occ if t in trait and trait[t] >= 0.5]
    lo = [occ[t] for t in occ if t in trait and trait[t] < 0.5]
    return abs(np.mean(hi) - np.mean(lo)) if hi and lo else 0.0


def test_cli_traits_genes_cid_observed_trait_is_decoupled(tmp_path):
    """The honest-null property: under the CID null the panel tracks the HIDDEN trait, while the
    OBSERVED trait a user would test is decoupled from it."""
    out = tmp_path / "o"
    rc = main(["coevolve", "--couple", "traits:genes", "-t", _tree_file(tmp_path, tips=80, age=6.0),
               "--trait-model", "mk", "--states", "2", "--trait-center", "--rate", "0.1",
               "--panel", "50", "--responsive", "0.8", "--effect-loss", "5", "--loss", "0.5",
               "--trans", "1.0", "--seed", "5", "--null", "cid", "--write", "profiles", "-o", str(out)])
    assert rc == 0
    occ = _read_presence(out / "presence.tsv")
    gap_hidden = _state_gap(occ, _read_trait(out / "trait_ground_truth.tsv"))
    gap_observed = _state_gap(occ, _read_trait(out / "traits.tsv"))
    # the panel strongly tracks the HIDDEN driver; the OBSERVED trait keeps only the faint
    # shared-tree imprint (the confound a real detector must see through) — never the causal signal
    assert gap_hidden > 0.3
    assert gap_hidden > gap_observed + 0.15


# =========================================================================== MuSSE / QuaSSE / HiSSE nulls
def test_musse_kstate_cid_is_character_independent_and_hides_the_hidden_class():
    Q = [[0, 0.2, 0.2], [0.2, 0, 0.2], [0.2, 0.2, 0]]
    alt = z.MuSSE(birth=[1, 1, 3], death=[0.2, 0.2, 0.2], Q=Q)   # observed state 2 fast
    cid = alt.null("cid", n_hidden=2)
    assert cid.k == 6                                            # 3 observed x 2 hidden
    b = cid.lambdas                                              # rate set by hidden class only
    assert b[0] == b[2] == b[4] and b[1] == b[3] == b[5]
    assert cid.discretize(0) == cid.discretize(1) == 0          # collapses to observed
    assert cid.discretize(4) == cid.discretize(5) == 2
    res = z.simulate_sse(cid, age=1.4, seed=1)
    assert set(res.labeled_values().values()) <= {0, 1, 2}     # output shows observed only


def test_musse_cid_needs_two_hidden():
    alt = z.MuSSE(birth=[1, 2], death=[0.2, 0.2], Q=[[0, 0.2], [0.2, 0]])
    with pytest.raises(ValueError, match="n_hidden >= 2"):
        alt.null("cid", n_hidden=1)


def test_quasse_neutral_is_constant_speciation():
    spec = z.QuaSSE.sigmoid(0.4, 3.0, center=0.0, slope=3.0)
    alt = z.QuaSSE(spec, lambda x: 0.2, sigma2=0.5, rate_bound=3.4, x0=-1.5)
    null = alt.null("neutral")
    assert null.speciation(-10) == null.speciation(10)          # trait-independent now
    assert null.speciation(0.0) == pytest.approx(spec(-1.5))    # constant = rate at x0
    assert alt.null("neutral", rate=1.0).speciation(99) == 1.0  # overridable


def test_quasse_neutral_removes_the_trait_bias():
    spec = z.QuaSSE.sigmoid(0.4, 3.0, center=0.0, slope=3.0)
    alt = z.QuaSSE(spec, lambda x: 0.2, sigma2=0.5, rate_bound=3.4, x0=-1.5)

    def mean_x(m, reps=40, seed=0):
        rng = np.random.default_rng(seed)
        return float(np.mean([v for _ in range(reps)
                              for v in z.simulate_sse(m, age=2.0, rng=rng).labeled_values().values()]))
    assert mean_x(alt) > mean_x(alt.null("neutral"))            # coupling pulls tips high; null doesn't


def test_quasse_has_no_cid_null():
    alt = z.QuaSSE(z.QuaSSE.sigmoid(0.4, 3, 0, 3), lambda x: 0.2, sigma2=0.5, rate_bound=3.4)
    with pytest.raises(TypeError, match="no 'cid'"):
        alt.null("cid")


def test_hisse_neutral_collapses_to_a_constant_rate_bisse():
    fast, slow = z.BiSSE(2.5, 2.5, 0.2, 0.2, 0.3, 0.3), z.BiSSE(0.4, 0.4, 0.2, 0.2, 0.3, 0.3)
    hisse = z.HiSSE([fast, slow], hidden_transition=0.15)
    null = hisse.null("neutral")
    assert isinstance(null, z.BiSSE)
    assert null.lambdas[0] == null.lambdas[1]                   # constant-rate on the observed trait
    assert np.isclose(null.Q[0, 1], 0.3) and np.isclose(null.Q[1, 0], 0.3)   # observed transitions kept
    assert set(z.simulate_sse(null, age=1.5, seed=1).labeled_values().values()) <= {0, 1}


def test_hisse_has_no_cid_null():
    hisse = z.HiSSE([z.BiSSE(2, 2, .2, .2, .3, .3), z.BiSSE(.4, .4, .2, .2, .3, .3)], 0.15)
    with pytest.raises(ValueError, match="already a hidden-state"):
        hisse.null("cid")


def test_cli_sse_model_hisse_runs(tmp_path):
    out = tmp_path / "o"
    rc = main(["coevolve", "--couple", "traits:species", "--sse-model", "hisse",
               "--hidden-classes", "2", "--hidden-scale", "4", "--hidden-switch", "0.15",
               "--lambda0", "0.6", "--lambda1", "0.6", "--q01", "0.3", "--q10", "0.3",
               "--tips", "100", "--seed", "1", "-o", str(out)])
    assert rc == 0 and (out / "species_tree.nwk").exists()


def test_cli_musse_cid_writes_observed_states_only(tmp_path):
    q = tmp_path / "q3.txt"
    q.write_text("0 0.2 0.2\n0.2 0 0.2\n0.2 0.2 0\n")
    out = tmp_path / "o"
    rc = main(["coevolve", "--couple", "traits:species", "--sse-model", "musse",
               "--birth", "1", "1", "3", "--death", "0.2", "0.2", "0.2", "--q-matrix", str(q),
               "--tips", "100", "--seed", "1", "--null", "cid", "-o", str(out)])
    assert rc == 0
    states = {r.split("\t")[1] for r in (out / "traits.tsv").read_text().splitlines()[1:]}
    assert states <= {"0", "1", "2"}                            # observed only, no (o, h) leak


def test_cli_hisse_and_quasse_cid_rejected(tmp_path):
    with pytest.raises(SystemExit):
        main(["coevolve", "--couple", "traits:species", "--sse-model", "hisse", "--tips", "40",
              "--null", "cid", "-o", str(tmp_path / "a")])
    with pytest.raises(SystemExit):
        main(["coevolve", "--couple", "traits:species", "--sse-model", "quasse", "--tips", "40",
              "--null", "cid", "-o", str(tmp_path / "b")])
