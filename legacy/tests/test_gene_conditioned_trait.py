"""Tests for gene-conditioned trait evolution (``zombi2.gene_conditioned_trait``: genes:traits).

A modifier gene switches a continuous trait's OU optimum, so the model is checked by its signal:

* tips carrying the modifier sit near ``theta_present``; those without near ``theta_absent``;
* deterministic corners: the modifier never appears -> the trait stays in the absent regime; the
  modifier is fixed present -> the trait is pulled to ``theta_present``;
* the run is reproducible; the model rejects an unreachable modifier.
"""

import numpy as np
import pytest

import zombi2 as z
from zombi2.cli import main


def _tree(seed=1, tips=80):
    return z.simulate_species_tree(z.BirthDeath(1, 0.3), n_tips=tips, age=6, seed=seed)


# --------------------------------------------------------------------------- validation
def test_gene_conditioned_trait_validation():
    with pytest.raises(ValueError):
        z.GeneConditionedTrait(gene_gain=-0.1)
    with pytest.raises(ValueError):
        z.GeneConditionedTrait(alpha=-1)
    with pytest.raises(ValueError):
        z.GeneConditionedTrait(sigma2=-1)
    with pytest.raises(ValueError):                              # can never appear
        z.GeneConditionedTrait(gene_gain=0, gene_loss=0, root_gene=False)


def test_gene_conditioned_trait_reproducible():
    tree = _tree()
    m = z.GeneConditionedTrait(gene_gain=0.6, gene_loss=0.6, theta_present=5.0, alpha=2.0)
    assert z.simulate_gene_conditioned_trait(tree, m, seed=9).to_tsv() == \
           z.simulate_gene_conditioned_trait(tree, m, seed=9).to_tsv()


# --------------------------------------------------------------------------- deterministic corners
def test_modifier_never_appears_stays_in_absent_regime():
    tree = _tree()
    m = z.GeneConditionedTrait(gene_gain=0.0, gene_loss=0.5, root_gene=False,
                               theta_absent=0.0, theta_present=5.0, alpha=3.0, sigma2=0.2)
    res = z.simulate_gene_conditioned_trait(tree, m, seed=1)
    assert set(res.gene_presence().values()) == {0}             # modifier never gained
    assert abs(np.mean(list(res.trait_values().values()))) < 1.0   # trait near theta_absent = 0


def test_modifier_fixed_present_pulls_to_present_optimum():
    tree = _tree()
    m = z.GeneConditionedTrait(gene_gain=0.0, gene_loss=0.0, root_gene=True,
                               theta_absent=0.0, theta_present=5.0, alpha=3.0, sigma2=0.2)
    res = z.simulate_gene_conditioned_trait(tree, m, seed=1)
    assert set(res.gene_presence().values()) == {1}             # modifier fixed present
    assert abs(np.mean(list(res.trait_values().values())) - 5.0) < 1.0   # trait near theta_present


def test_brownian_limit_alpha_zero_runs():
    """alpha=0 is the Brownian-motion limit (no pull); the run still works."""
    tree = _tree(tips=30)
    m = z.GeneConditionedTrait(gene_gain=0.5, gene_loss=0.5, alpha=0.0, sigma2=0.5)
    res = z.simulate_gene_conditioned_trait(tree, m, seed=2)
    assert all(isinstance(v, float) for v in res.trait_values().values())


# --------------------------------------------------------------------------- the coupling signal
def test_carriers_track_present_optimum():
    """Aggregated over trees: modifier carriers sit well above non-carriers (toward theta_present)."""
    tree = _tree()
    m = z.GeneConditionedTrait(gene_gain=0.6, gene_loss=0.6, theta_absent=0.0, theta_present=5.0,
                               alpha=2.5, sigma2=0.4)
    rng = np.random.default_rng(0)
    car, non = [], []
    for _ in range(25):
        res = z.simulate_gene_conditioned_trait(tree, m, rng=rng)
        tv, gp = res.trait_values(), res.gene_presence()
        for leaf in tv:
            (car if gp[leaf] else non).append(tv[leaf])
    assert np.mean(car) > np.mean(non) + 2.0                     # clear phenotypic separation


# --------------------------------------------------------------------------- CLI (genes:traits edge)
def test_cli_genes_traits(tmp_path):
    sp = tmp_path / "sp"
    main(["species", "--tips", "40", "--age", "5", "--seed", "1", "-o", str(sp)])
    out = tmp_path / "gt"
    rc = main(["coevolve", "--couple", "genes:traits", "-t", str(sp / "species_tree.nwk"),
               "--modifier-gain", "0.6", "--modifier-loss", "0.6", "--theta-present", "5",
               "--trait-alpha", "2.5", "--seed", "2", "-o", str(out)])
    assert rc == 0
    assert (out / "species_tree.nwk").exists() and (out / "trait_tree.nwk").exists()
    header = (out / "traits.tsv").read_text().splitlines()[0]
    assert header == "node\tmodifier\ttrait"
    for row in (out / "traits.tsv").read_text().splitlines()[1:]:
        assert row.split("\t")[1] in ("0", "1")                 # modifier column is 0/1


def test_cli_genes_traits_needs_tree(tmp_path):
    with pytest.raises(SystemExit):
        main(["coevolve", "--couple", "genes:traits", "-o", str(tmp_path / "a")])


def test_cli_genes_traits_rejects_age(tmp_path):
    sp = tmp_path / "sp"
    main(["species", "--tips", "20", "--seed", "1", "-o", str(sp)])
    with pytest.raises(SystemExit):
        main(["coevolve", "--couple", "genes:traits", "-t", str(sp / "species_tree.nwk"),
              "--age", "3", "-o", str(tmp_path / "a")])
