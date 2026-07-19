"""Explicit per-family substitution speed (``SequenceEvolution.family_factors``).

The sequence-level analogue of :class:`~zombi2.FamilyModifier`: a *named* per-family multiplier on
the substitution rate that composes with (multiplies) the random ``family_speed`` draw and the
per-branch lineage clock — "make a given family evolve faster, on top of branch and random effects".
"""

import pytest

from zombi2 import (
    BirthDeath,
    LogNormal,
    Rates,
    SequenceEvolution,
    simulate_genomes,
    simulate_species_tree,
)


def _genomes():
    tree = simulate_species_tree(BirthDeath(1.0, 0.2), n_tips=8, age=3.0, seed=1)
    return simulate_genomes(tree, Rates(duplication=0.2, loss=0.1, origination=0.5),
                            initial_families=8, seed=1)


def _total_subst(node_trees, fam):
    entry = node_trees[fam]["complete"]
    return sum(entry[1].values()) if entry else 0.0


def test_family_factor_scales_substitution_lengths():
    g = _genomes()
    # strict clock + fixed family_speed=1.0: the ONLY difference is the explicit factor
    base = SequenceEvolution(family_speed=1.0)
    scaled = SequenceEvolution(family_speed=1.0, family_factors={"1": 3.0})
    args = (g.species_tree, g.gene_families, g._gid_to_species())
    _, tb = base.scale_families_trees(*args, seed=5)
    _, ts = scaled.scale_families_trees(*args, seed=5)
    assert _total_subst(tb, "1") > 0
    assert abs(_total_subst(ts, "1") - 3.0 * _total_subst(tb, "1")) < 1e-9   # family 1 = 3x
    assert abs(_total_subst(ts, "2") - _total_subst(tb, "2")) < 1e-9         # family 2 unchanged


def test_family_factor_composes_with_random_speed():
    g = _genomes()
    base = SequenceEvolution(family_speed=LogNormal(0.0, 0.5))
    scaled = SequenceEvolution(family_speed=LogNormal(0.0, 0.5), family_factors={"1": 10.0})
    pb = base.scale(g, seed=7)
    ps = scaled.scale(g, seed=7)
    # family 1's reported speed = 10x its random draw; every other family's draw is untouched
    assert abs(ps.family_speed["1"] - 10.0 * pb.family_speed["1"]) < 1e-9
    for fam, s in pb.family_speed.items():
        if fam != "1":
            assert abs(ps.family_speed[fam] - s) < 1e-9


def test_family_factor_composes_with_branch_clock():
    g = _genomes()
    base = SequenceEvolution(branch_sigma=0.5, family_speed=1.0)
    scaled = SequenceEvolution(branch_sigma=0.5, family_speed=1.0, family_factors={"1": 2.0})
    args = (g.species_tree, g.gene_families, g._gid_to_species())
    _, tb = base.scale_families_trees(*args, seed=9)
    _, ts = scaled.scale_families_trees(*args, seed=9)
    # under a non-trivial clock, family 1 is scaled 2x uniformly (family x branch = product)
    assert _total_subst(tb, "1") > 0
    assert abs(_total_subst(ts, "1") - 2.0 * _total_subst(tb, "1")) < 1e-9


def test_family_factors_validation_and_noop_default():
    with pytest.raises(ValueError):
        SequenceEvolution(family_factors={"1": -2.0})
    # an empty map is a no-op: identical to not passing one (same RNG stream, same speeds)
    g = _genomes()
    a = SequenceEvolution(family_speed=LogNormal(0.0, 0.4)).scale(g, seed=3)
    b = SequenceEvolution(family_speed=LogNormal(0.0, 0.4), family_factors={}).scale(g, seed=3)
    assert a.family_speed == b.family_speed
