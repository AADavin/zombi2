"""FamilyModifier — per-family rate heterogeneity as a composable overlay (the emission seam).

Contrast :class:`FamilySampledRates` (per-family rates baked into a base): a ``FamilyModifier``
stacks on *any* base, so per-genome × per-family and family × branch become expressible. The
composer expands a base's aggregate ``family=None`` weights into per-family weights so the family
factor has a family to attach to.
"""

import numpy as np
import pytest

from zombi2 import (
    BirthDeath,
    BranchModifier,
    FamilyModifier,
    LogNormal,
    ModifiedRates,
    Rates,
    PerGenomeRates,
    simulate_genomes,
    simulate_species_tree,
)


def _tree():
    return simulate_species_tree(BirthDeath(1.0, 0.2), n_tips=10, age=3.0, seed=1)


def _family_copies(g, fam):
    """Total copies of family ``fam`` summed across all leaf species."""
    fams = g.profiles.families
    return int(g.profiles.matrix[fams.index(fam)].sum()) if fam in fams else 0


def test_explicit_family_factor_boosts_that_family():
    tree = _tree()

    def copies(factor, seeds=range(6)):
        # Loss-only (pure death) keeps genomes bounded; scaling *only* family "1"'s loss makes it
        # disappear from more leaves. origination=0 keeps the family set to the seeded 1..10.
        rates = ModifiedRates(Rates(loss=0.25, origination=0.0),
                              [FamilyModifier(factors={"1": factor}, events=("loss",))])
        return sum(_family_copies(simulate_genomes(tree, rates, initial_families=10, seed=s), "1")
                   for s in seeds)

    # family "1" with 6x loss survives in far fewer leaves than unscaled
    assert copies(6.0) < copies(1.0)


def test_per_genome_times_per_family_composes():
    # per-genome base (constant totals, linear growth) × per-family overlay — the combination the old
    # --rate-model enum could not name. Per-genome rates keep it intrinsically bounded.
    tree = _tree()
    fam_mod = FamilyModifier(per_family=LogNormal(0.0, 0.7))
    rates = ModifiedRates(
        PerGenomeRates(duplication=1.0, transfer=0.2, loss=0.5, origination=0.4), [fam_mod])
    g = simulate_genomes(tree, rates, initial_families=10, seed=3)
    assert g.profiles.matrix.shape[1] == 10
    assert len(g.event_log) > 0
    assert len(set(fam_mod._factor.values())) > 1  # families really drew different factors


def test_family_and_branch_compose():
    # stack a family modifier and a branch modifier on one base: family × branch heterogeneity.
    # dup < loss keeps every family a shrinking (bounded) process whatever factors it draws.
    tree = _tree()
    rates = ModifiedRates(
        Rates(duplication=0.2, loss=0.3, origination=0.3),
        [FamilyModifier(per_family=LogNormal(0.0, 0.5)), BranchModifier(autocorr_sigma=0.4)])
    g = simulate_genomes(tree, rates, initial_families=10, seed=5, max_family_size=40)
    assert g.profiles.matrix.shape[1] == 10
    assert len(g.event_log) > 0


def test_reproducible_given_seed():
    tree = _tree()
    make = lambda: ModifiedRates(Rates(duplication=0.2, loss=0.2, origination=0.3),
                                 [FamilyModifier(per_family=LogNormal(0.0, 0.5))])
    a = simulate_genomes(tree, make(), initial_families=10, seed=11, max_family_size=40)
    b = simulate_genomes(tree, make(), initial_families=10, seed=11, max_family_size=40)
    assert np.array_equal(a.profiles.matrix, b.profiles.matrix)


def test_family_modifier_requires_exactly_one_source():
    with pytest.raises(ValueError):
        FamilyModifier()  # no source
    with pytest.raises(ValueError):
        FamilyModifier(factors={"1": 2.0}, per_family=LogNormal(0.0, 0.5))  # two sources
