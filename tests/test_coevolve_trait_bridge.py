"""Tests for the optimum-target bridge (:mod:`zombi2.coevolve.trait_bridge`).

The bridge walks an OU trait whose optimum is set at each point by a grammar :class:`Response`
applied to a :class:`DriverSignal` — realizing a state-target coupling on ``traits.optimum``. Checked
by its signal: the trait is pulled to the driver-selected optimum, the walk is reproducible, and the
Brownian (``alpha=0``) limit runs. (Byte-identity to the pre-reframe ``genes:traits`` edge is
verified in ``tests/test_gene_conditioned_trait.py``, which now runs through this bridge.)
"""

import numpy as np
import pytest

import zombi2 as z
from zombi2.coevolve.grammar import Table
from zombi2.coevolve.trait_bridge import walk_optimum_coupled_trait


class _ConstDriver:
    """A trivial DriverSignal: one constant value, no change points."""

    def __init__(self, value):
        self._v = float(value)

    def value(self, lineage, time):
        return self._v

    def refresh_times(self, t0, t1):
        return []


def _tree(seed=1, tips=60):
    return z.simulate_species_tree(z.BirthDeath(1, 0.3), n_tips=tips, age=6, seed=seed)


def test_trait_is_pulled_to_the_driver_selected_optimum():
    tree = _tree()
    response = Table({0: 0.0, 1: 5.0})                          # absent → 0, present → 5
    rng = np.random.default_rng(0)
    present = walk_optimum_coupled_trait(tree, _ConstDriver(1.0), response,
                                         alpha=3.0, sigma2=0.2, x0=0.0, rng=rng)
    absent = walk_optimum_coupled_trait(tree, _ConstDriver(0.0), response,
                                        alpha=3.0, sigma2=0.2, x0=0.0, rng=rng)
    tips_present = [present[l] for l in tree.extant_leaves()]
    tips_absent = [absent[l] for l in tree.extant_leaves()]
    assert abs(np.mean(tips_present) - 5.0) < 1.0              # pulled to θ_present
    assert abs(np.mean(tips_absent) - 0.0) < 1.0              # pulled to θ_absent


def test_walk_is_reproducible_under_the_same_rng():
    tree = _tree(tips=40)
    response = Table({0: 0.0, 1: 4.0})
    a = walk_optimum_coupled_trait(tree, _ConstDriver(1.0), response,
                                   alpha=2.0, sigma2=0.3, x0=1.0, rng=np.random.default_rng(5))
    b = walk_optimum_coupled_trait(tree, _ConstDriver(1.0), response,
                                   alpha=2.0, sigma2=0.3, x0=1.0, rng=np.random.default_rng(5))
    assert a == b


def test_root_value_is_x0_and_every_node_is_assigned():
    tree = _tree(tips=30)
    nt = walk_optimum_coupled_trait(tree, _ConstDriver(0.0), Table({0: 0.0, 1: 5.0}),
                                    alpha=1.0, sigma2=0.1, x0=2.5, rng=np.random.default_rng(1))
    assert nt[tree.root] == 2.5
    assert set(nt) == set(tree.nodes())                       # a value at every node
    assert all(isinstance(v, float) for v in nt.values())


def test_brownian_limit_alpha_zero_runs():
    tree = _tree(tips=25)
    nt = walk_optimum_coupled_trait(tree, _ConstDriver(1.0), Table({0: 0.0, 1: 5.0}),
                                    alpha=0.0, sigma2=0.5, x0=0.0, rng=np.random.default_rng(2))
    assert all(isinstance(v, float) for v in nt.values())


def test_stronger_pull_lands_closer_to_the_optimum():
    tree = _tree(tips=80)
    response = Table({0: 0.0, 1: 5.0})
    strong = walk_optimum_coupled_trait(tree, _ConstDriver(1.0), response,
                                        alpha=5.0, sigma2=0.2, x0=0.0, rng=np.random.default_rng(3))
    weak = walk_optimum_coupled_trait(tree, _ConstDriver(1.0), response,
                                      alpha=0.3, sigma2=0.2, x0=0.0, rng=np.random.default_rng(3))
    tips_strong = np.mean([abs(strong[l] - 5.0) for l in tree.extant_leaves()])
    tips_weak = np.mean([abs(weak[l] - 5.0) for l in tree.extant_leaves()])
    assert tips_strong < tips_weak                            # stronger α → tighter around θ
