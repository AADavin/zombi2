"""Tests for the sequence-tier bridge (:mod:`zombi2.coevolve.sequence_bridge`).

:class:`DriverClock` sets a lineage's substitution rate from a grammar coupling on
``sequences.substitution_speed``. Checked against the :class:`~zombi2.sequences.clocks.Clock`
contract: it reduces to a strict clock under a null (or constant) driver, scales per-lineage with the
driver, tracks a within-branch change, and reports the correct time-averaged rate.
"""

import math

import pytest

import zombi2 as z
from zombi2.coevolve.grammar import Scalar
from zombi2.coevolve.sequence_bridge import DriverClock
from zombi2.sequences.clocks import StrictClock


class _ConstDriver:
    def __init__(self, value):
        self._v = float(value)

    def value(self, lineage, time):
        return self._v

    def refresh_times(self, t0, t1):
        return []


class _MapDriver:
    """A per-lineage constant driver: value depends only on the branch name."""

    def __init__(self, by_branch, default=0.0):
        self._m = dict(by_branch)
        self._d = float(default)

    def value(self, lineage, time):
        return self._m.get(lineage, self._d)

    def refresh_times(self, t0, t1):
        return []


class _SwitchDriver:
    """0 before ``t_switch`` on ``branch``, 1 after — to exercise a within-branch change."""

    def __init__(self, branch, t_switch):
        self._b = branch
        self._t = float(t_switch)

    def value(self, lineage, time):
        return 1.0 if (lineage == self._b and time >= self._t) else 0.0

    def refresh_times(self, t0, t1):
        return [(self._t, self._b)] if t0 < self._t < t1 else []


def _tree(seed=1, tips=40):
    return z.simulate_species_tree(z.BirthDeath(1, 0.3), n_tips=tips, age=6, seed=seed)


def test_null_driver_reduces_to_a_strict_clock():
    tree = _tree()
    driven = DriverClock(_ConstDriver(9.0), Scalar(0.0), base_rate=2.0)   # null: multiplier 1
    strict = StrictClock(2.0)
    dsegs, davg = driven.lineage_segments(tree, None)
    ssegs, savg = strict.lineage_segments(tree, None)
    assert davg == savg
    assert dsegs == ssegs


def test_constant_driver_is_a_uniform_rescale():
    tree = _tree()
    clock = DriverClock(_ConstDriver(1.0), Scalar(0.5), base_rate=2.0)
    segs, avg = clock.lineage_segments(tree, None)
    expected = 2.0 * math.exp(0.5)
    for node in tree.nodes_preorder():
        if node.parent is None:
            continue
        (rate, t0, t1), = segs[node.name]                     # one segment per branch
        assert rate == pytest.approx(expected)
        assert avg[node.name] == pytest.approx(expected)


def test_rate_scales_per_lineage_with_the_driver():
    tree = _tree(tips=30)
    leaves = list(tree.extant_leaves())
    hot, cold = leaves[0].name, leaves[1].name
    clock = DriverClock(_MapDriver({hot: 1.0, cold: -1.0}, default=0.0), Scalar(1.0), base_rate=1.0)
    _, avg = clock.lineage_segments(tree, None)
    assert avg[hot] == pytest.approx(math.exp(1.0))          # driver +1 → faster
    assert avg[cold] == pytest.approx(math.exp(-1.0))        # driver −1 → slower
    assert avg[hot] > avg[cold]


def test_within_branch_change_splits_into_two_segments():
    tree = _tree(tips=20)
    leaf = list(tree.extant_leaves())[0]
    b0, b1 = leaf.parent.time, leaf.time
    mid = 0.5 * (b0 + b1)
    clock = DriverClock(_SwitchDriver(leaf.name, mid), Scalar(1.0), base_rate=1.0)
    segs, avg = clock.lineage_segments(tree, None)
    pieces = segs[leaf.name]
    assert len(pieces) == 2                                   # split at the change point
    assert pieces[0][0] == pytest.approx(1.0)                # driver 0 → exp(0) = 1
    assert pieces[1][0] == pytest.approx(math.exp(1.0))      # driver 1 → exp(1)
    # time-averaged rate over the equal halves
    assert avg[leaf.name] == pytest.approx(0.5 * (1.0 + math.exp(1.0)))


def test_scale_produces_a_phylogram():
    tree = _tree(tips=25)
    clock = DriverClock(_ConstDriver(1.0), Scalar(0.3), base_rate=1.5)
    scaled = clock.scale(tree)                                # exercises the full Clock machinery
    assert scaled.to_newick().endswith(";")


def test_base_rate_must_be_positive():
    with pytest.raises(ValueError, match="base_rate"):
        DriverClock(_ConstDriver(1.0), Scalar(0.0), base_rate=0.0)


# ── selection (ω): OmegaSelector (T→Σ) and GeneEventOmega (G→Σ) ────────────────
import numpy as np

from zombi2.coevolve.sequence_bridge import GeneEventOmega, OmegaSelector
from zombi2.genomes.events import EventType
from zombi2.sequences.models import GammaRates, evolve_on_tree


class _FakeNode:
    """A minimal reconciliation-node stand-in for the ω selectors / evolve_on_tree."""

    def __init__(self, gid="g", *, branch=None, species=None, birth=0.0, end=1.0, kind=None,
                 children=None):
        self.gid = gid
        self.branch = branch
        self.species = species
        self.birth = birth
        self.end = end
        self.kind = kind
        self.children = children or []


def test_omega_selector_scales_omega_with_the_trait():
    node = _FakeNode(branch="b1", birth=0.0, end=2.0)
    sel = OmegaSelector(_ConstDriver(1.0), Scalar(0.5), base_omega=0.2)
    assert sel.omega_for(node) == pytest.approx(0.2 * math.exp(0.5))     # driver 1.0
    null = OmegaSelector(_ConstDriver(9.0), Scalar(0.0), base_omega=0.2)
    assert null.omega_for(node) == pytest.approx(0.2)                    # null → uniform ω


def test_omega_class_cache_reuses_one_model_per_class():
    sel = OmegaSelector(_MapDriver({"hot": 1.0, "cold": -1.0}), Scalar(1.0),
                        base_omega=0.3, resolution=0.02)
    m_hot = sel.model_for(_FakeNode(branch="hot", birth=0.0, end=1.0))
    m_hot_again = sel.model_for(_FakeNode(branch="hot", birth=0.5, end=1.5))
    m_cold = sel.model_for(_FakeNode(branch="cold"))
    assert m_hot is m_hot_again                                          # same ω class → cached
    assert m_cold is not m_hot                                          # different ω class
    assert m_hot.k == 61                                                # a 61-state codon model


def test_gene_event_omega_relaxes_selection_on_the_event_branch():
    geo = GeneEventOmega(Scalar(0.7), base_omega=0.2, events=(EventType.DUPLICATION,))
    dup = _FakeNode(kind=EventType.DUPLICATION)
    spec = _FakeNode(kind=EventType.SPECIATION)
    assert geo.omega_for(dup) == pytest.approx(0.2 * math.exp(0.7))     # relaxed after duplication
    assert geo.omega_for(spec) == pytest.approx(0.2)                    # base elsewhere
    assert geo.model_for(dup) is not geo.model_for(spec)


def test_omega_base_model_is_a_codon_model_and_validation():
    assert OmegaSelector(_ConstDriver(0.0), Scalar(1.0)).base_model.k == 61
    with pytest.raises(ValueError, match="base_omega"):
        OmegaSelector(_ConstDriver(1.0), Scalar(0.0), base_omega=-1.0)
    with pytest.raises(ValueError, match="resolution"):
        GeneEventOmega(Scalar(0.0), resolution=0.0)


def test_evolve_on_tree_uses_the_per_branch_omega_model():
    sel = OmegaSelector(_MapDriver({"root_sp": 0.0, "child_sp": 2.0}), Scalar(1.0), base_omega=0.2)
    root = _FakeNode("r", branch="root_sp", birth=0.0, end=1.0)
    child = _FakeNode("c", branch="child_sp", birth=1.0, end=3.0)
    root.children = [child]
    out = evolve_on_tree(root, {root: 0.0, child: 0.5}, sel.base_model,
                         np.random.default_rng(0), length=6, model_for=sel.model_for)
    assert set(out) == {"r", "c"}
    assert len(out["c"]) == 18                                          # 6 codon sites × 3 nt
    assert sel.model_for(child) is not sel.base_model                   # child ran under its own ω


def test_model_for_is_mutually_exclusive_with_gamma():
    sel = OmegaSelector(_ConstDriver(1.0), Scalar(0.5))
    root = _FakeNode("r")
    with pytest.raises(ValueError, match="mutually exclusive"):
        evolve_on_tree(root, {root: 0.0}, sel.base_model, np.random.default_rng(0),
                       length=3, model_for=sel.model_for, gamma=GammaRates(0.5))
