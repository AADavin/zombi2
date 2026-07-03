"""Tests for the gene-family coupling model (:mod:`zombi2.coupling`).

Three layers:

1. **Unit** — the coupling spec constructors and the exact loss-rate formula
   ``base_loss·exp(-β·f_i)`` with ``f_i = h_i + Σ_j J_ij σ_j`` (partners only).
2. **Ground-truth recovery** — inject a known ``J`` and confirm the generated profiles show
   the prescribed co-occurrence structure: positive ``J`` → co-occurrence, negative ``J`` →
   avoidance, zero ``J`` → no structure, and uncoupled families stay uncorrelated.
3. **Driver** — :func:`simulate_coupled` shape/panel/reproducibility contract.

The recovery tests run on a *near-star* tree (all lineages split just below the root, then
evolve independently for the whole age) which isolates the injected coupling from the
phylogenetic confounding that inflates co-occurrence on a normal birth–death tree — the very
"shared ancestry" trap the design note (``docs/non_independence.tex``) warns about, and the
reason real inference (Fukunaga & Iwasaki 2022) corrects for the tree.
"""

import itertools
import math

import numpy as np
import pytest

from zombi2 import BirthDeath, simulate_species_tree
from zombi2.coupling import (
    CouplingSpec,
    PottsRates,
    pathway_blocks,
    simulate_coupled,
)
from zombi2.events import EventType
from zombi2.genome import Gene, IdManager, UnorderedGenome
from zombi2.tree import Tree, TreeNode


# ── fixtures ─────────────────────────────────────────────────────────────────
def near_star_tree(k: int, age: float, delta: float = 1e-3) -> Tree:
    """A balanced binary tree of ``2**k`` tips whose internal nodes all sit within
    ``[0, delta·k]`` of the root — so every terminal branch is ~``age`` long and the tips
    are near-independent draws (minimal shared ancestry)."""
    counter = itertools.count()
    root = TreeNode(name="root", time=0.0)

    def split(node, level):
        if level == k:
            return
        for _ in range(2):
            child = TreeNode(name=f"i{next(counter)}", time=min(age, delta * (level + 1)))
            node.add_child(child)
            split(child, level + 1)

    split(root, 0)
    stack, leaves = [root], []
    while stack:
        n = stack.pop()
        leaves.append(n) if not n.children else stack.extend(n.children)
    for idx, leaf in enumerate(leaves):
        leaf.time, leaf.name, leaf.is_extant = age, f"n{idx}", True
    return Tree(root, age)


def _corr(a, b) -> float:
    a, b = a.astype(float), b.astype(float)
    if a.std() == 0 or b.std() == 0:
        return 0.0
    return float(np.corrcoef(a, b)[0, 1])


def _make_genome(families):
    ids = IdManager()
    g = UnorderedGenome(ids)
    for fam in families:
        g._add(Gene(ids.new_gene(), fam))
    return g


# canonical recovery-test regime (see the parameter sweep in the design work): moderate
# occupancy (~0.78) so both co-occurrence and avoidance register.
_REGIME = dict(base_loss=1.0, transfer=0.2, beta=1.0, h=2.0)


# ── 1. unit: spec constructors ───────────────────────────────────────────────
def test_from_dense_and_from_edges_agree():
    J = np.array([[0.0, 1.2, -0.4],
                  [1.2, 0.0, 0.0],
                  [-0.4, 0.0, 0.0]])
    a = CouplingSpec.from_dense(J)
    b = CouplingSpec.from_edges(3, {(0, 1): 1.2, (0, 2): -0.4})
    assert np.allclose(a.dense_J(), J)
    assert np.allclose(b.dense_J(), J)          # edges are symmetrised
    assert np.allclose(a.dense_J(), b.dense_J())


def test_from_edges_rejects_self_and_out_of_range():
    with pytest.raises(ValueError):
        CouplingSpec.from_edges(3, {(1, 1): 1.0})          # self-coupling
    with pytest.raises(ValueError):
        CouplingSpec.from_edges(3, {(0, 5): 1.0})          # out of range


def test_pathway_blocks_structure():
    spec = pathway_blocks([2, 2], within=3.0, between=-1.0)
    J = spec.dense_J()
    assert np.all(np.diag(J) == 0.0)                       # zero diagonal
    assert J[0, 1] == 3.0 and J[2, 3] == 3.0               # within-block
    assert J[0, 2] == -1.0 and J[1, 3] == -1.0             # between-block
    assert np.allclose(J, J.T)                             # symmetric
    assert spec.panel_ids == ["F0", "F1", "F2", "F3"]


def test_bad_spec_rejected():
    with pytest.raises(ValueError):
        CouplingSpec.from_dense(np.zeros((3, 3)), base_loss=-1.0)
    with pytest.raises(ValueError):
        CouplingSpec.from_dense(np.zeros((3, 3)), h=np.zeros(2))   # wrong h length


# ── 1. unit: the loss-rate formula ───────────────────────────────────────────
def test_loss_rate_matches_field_formula():
    spec = CouplingSpec.from_edges(
        3, {(0, 1): 1.2, (0, 2): -0.4},
        h=[0.5, -0.3, 0.0], base_loss=1.0, beta=0.5, transfer=0.7,
    )
    rates = PottsRates(spec)
    g = _make_genome(["F0", "F1", "F2"])            # all three present
    ws = rates.event_weights(g, "b", 0.0)
    loss = {e.family: e.rate for e in ws if e.event is EventType.LOSS}

    # f0 = 0.5 + J01·1 + J02·1 = 0.5 + 1.2 - 0.4 = 1.3
    assert loss["F0"] == pytest.approx(math.exp(-0.5 * 1.3))
    # f1 = -0.3 + J01·1 = 0.9
    assert loss["F1"] == pytest.approx(math.exp(-0.5 * 0.9))
    # f2 =  0.0 + J02·1 = -0.4   → field < 0 raises loss above base
    assert loss["F2"] == pytest.approx(math.exp(-0.5 * -0.4))
    assert loss["F2"] > spec.base_loss

    # transfer is the field-blind gain channel: per-copy rate × genome size
    trans = [e for e in ws if e.event is EventType.TRANSFER]
    assert len(trans) == 1 and trans[0].family is None
    assert trans[0].rate == pytest.approx(0.7 * g.size())


def test_present_partner_changes_the_field():
    spec = CouplingSpec.from_edges(3, {(0, 1): 1.2, (0, 2): -0.4},
                                   h=[0.5, -0.3, 0.0], base_loss=1.0, beta=0.5)
    rates = PottsRates(spec)
    # drop F2: now f0 = 0.5 + 1.2 = 1.7 (only the F1 partner remains)
    g = _make_genome(["F0", "F1"])
    loss = {e.family: e.rate for e in rates.event_weights(g, "b", 0.0)
            if e.event is EventType.LOSS}
    assert loss["F0"] == pytest.approx(math.exp(-0.5 * 1.7))


def test_positive_partner_protects_negative_partner_exposes():
    spec = CouplingSpec.from_edges(3, {(0, 1): 2.0, (0, 2): -2.0}, base_loss=1.0, beta=1.0)
    rates = PottsRates(spec)
    solo = {e.family: e.rate for e in rates.event_weights(_make_genome(["F0"]), "b", 0)
            if e.event is EventType.LOSS}["F0"]
    with_pos = {e.family: e.rate for e in rates.event_weights(_make_genome(["F0", "F1"]), "b", 0)
                if e.event is EventType.LOSS}["F0"]
    with_neg = {e.family: e.rate for e in rates.event_weights(_make_genome(["F0", "F2"]), "b", 0)
                if e.event is EventType.LOSS}["F0"]
    assert with_pos < solo < with_neg                     # protection vs exposure
    assert solo == pytest.approx(1.0)                     # field 0 → base_loss


def test_non_panel_family_is_uncoupled():
    spec = CouplingSpec.from_edges(2, {(0, 1): 3.0}, base_loss=0.7, beta=1.0)
    rates = PottsRates(spec)
    g = _make_genome(["F0", "F1", "X"])                   # X is not in the panel
    loss = {e.family: e.rate for e in rates.event_weights(g, "b", 0.0)
            if e.event is EventType.LOSS}
    assert loss["X"] == pytest.approx(0.7)               # base_loss, no field


def test_origination_channel_optional():
    spec = CouplingSpec.from_edges(2, {(0, 1): 1.0}, origination=0.3)
    ws = PottsRates(spec).event_weights(_make_genome(["F0"]), "b", 0.0)
    orig = [e for e in ws if e.event is EventType.ORIGINATION]
    assert len(orig) == 1 and orig[0].rate == pytest.approx(0.3)
    # default: no origination channel (closed panel)
    ws0 = PottsRates(CouplingSpec.from_edges(2, {(0, 1): 1.0})).event_weights(
        _make_genome(["F0"]), "b", 0.0)
    assert not [e for e in ws0 if e.event is EventType.ORIGINATION]


# ── 2. ground-truth recovery ─────────────────────────────────────────────────
def test_positive_coupling_creates_cooccurrence():
    """Injected +J → coupled pair co-occurs, while an uncoupled pair does not."""
    tree = near_star_tree(8, age=6.0)
    spec = pathway_blocks([2, 1, 1], within=3.0, between=0.0, **_REGIME)  # (0,1) coupled
    P = simulate_coupled(tree, spec, seed=3).profiles.presence()
    coupled = _corr(P[0], P[1])
    uncoupled = _corr(P[2], P[3])
    assert coupled > 0.25
    assert coupled > uncoupled + 0.2


def test_zero_coupling_is_null():
    """J = 0 → the same panel positions show no injected co-occurrence."""
    tree = near_star_tree(8, age=6.0)
    spec = pathway_blocks([2, 1, 1], within=0.0, between=0.0, **_REGIME)
    P = simulate_coupled(tree, spec, seed=3).profiles.presence()
    assert abs(_corr(P[0], P[1])) < 0.25
    assert abs(_corr(P[2], P[3])) < 0.25


def test_negative_coupling_creates_avoidance():
    """Injected -J → coupled pair avoids each other (anti-correlated presence)."""
    tree = near_star_tree(8, age=6.0)
    spec = pathway_blocks([2, 1, 1], within=-3.0, between=0.0, **_REGIME)
    P = simulate_coupled(tree, spec, seed=3).profiles.presence()
    assert _corr(P[0], P[1]) < -0.25


def test_coupling_beats_no_coupling_on_the_same_pair():
    """Phylogeny-controlled: on one tree, the coupled pair's correlation exceeds the same
    pair's correlation under J = 0 (differencing out the shared-ancestry baseline)."""
    tree = near_star_tree(8, age=6.0)
    coupled = simulate_coupled(
        tree, pathway_blocks([2, 1, 1], within=3.0, **_REGIME), seed=3).profiles.presence()
    null = simulate_coupled(
        tree, pathway_blocks([2, 1, 1], within=0.0, **_REGIME), seed=3).profiles.presence()
    assert _corr(coupled[0], coupled[1]) > _corr(null[0], null[1]) + 0.2


# ── 3. driver contract ───────────────────────────────────────────────────────
def test_profile_shape_and_panel_rows():
    tree = simulate_species_tree(BirthDeath(1.0, 0.2), n_tips=25, age=4.0, seed=5)
    spec = pathway_blocks([3, 2], within=2.0, **_REGIME)      # 5-family panel
    res = simulate_coupled(tree, spec, seed=1)
    # every panel family is a row (even any that went globally extinct), species = tips
    assert res.profiles.families == spec.panel_ids
    assert res.profiles.shape == (5, 25)
    assert res.profiles.shape[1] == len(tree.extant_leaves())


def test_reproducible_with_seed():
    tree = simulate_species_tree(BirthDeath(1.0, 0.2), n_tips=20, age=4.0, seed=2)
    spec = pathway_blocks([2, 2], within=2.0, **_REGIME)
    a = simulate_coupled(tree, spec, seed=7).profiles
    b = simulate_coupled(tree, spec, seed=7).profiles
    assert np.array_equal(a.matrix, b.matrix)


def test_initial_presence_mask():
    tree = simulate_species_tree(BirthDeath(1.0, 0.2), n_tips=15, age=1.0, seed=9)
    spec = pathway_blocks([2, 2], within=2.0, base_loss=0.0, transfer=0.0, h=0.0)
    # no loss, no transfer, tiny age → the initial presence pattern is preserved at the tips
    mask = np.array([1, 0, 1, 0])
    res = simulate_coupled(tree, spec, seed=1, initial_presence=mask)
    present_rows = set(np.unique(res.profiles.coo[0]))
    assert present_rows == {0, 2}                            # only the seeded families


def test_runs_on_realistic_tree():
    """Smoke test on an ordinary birth–death tree (with phylogenetic structure)."""
    tree = simulate_species_tree(BirthDeath(1.0, 0.3), n_tips=40, age=5.0, seed=4)
    spec = pathway_blocks([4, 4], within=2.5, between=-1.0, **_REGIME)
    res = simulate_coupled(tree, spec, seed=2)
    assert res.profiles.shape == (8, 40)
    assert len(res.event_log) > 0
