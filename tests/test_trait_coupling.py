"""Tests for trait-conditioned gene families (:mod:`zombi2.trait_coupling`).

Four layers:

1. **Unit — trajectory** — a discrete trait's exact stochastic map becomes exact per-branch
   segments and refresh points; a continuous trait is sub-segmented with interpolated values.
2. **Unit — coupling** — the ``build`` constructor (count / fraction / id list, signed weights)
   and the exact loss-rate formula ``base_loss·exp(-effect_loss·w_i·s)``.
3. **Ground-truth recovery** — paint a deterministic trait onto a near-star tree (two clades,
   one favoured by the trait, one not) and confirm the responsive families track the trait at
   the tips while inert families do not.
4. **Driver / CLI** — shape, reproducibility, and the ``coevolve --couple traits:genes`` edge.

The recovery test uses a near-star tree (all lineages split just below the root, then evolve
independently for the whole age) to isolate the injected trait→gene signal from the shared-
ancestry confounding a normal birth–death tree would add — the same device as
``tests/test_coupling.py``.
"""

import itertools
import math

import numpy as np
import pytest

from zombi2 import (
    BrownianMotion, Mk, MultivariateBrownian, TraitResult, simulate_species_tree,
    simulate_traits, BirthDeath,
)
from zombi2.events import EventType
from zombi2.genome import Gene, IdManager, UnorderedGenome
from zombi2.trait_coupling import (
    TraitGeneCoupling, TraitLinkedRates, TraitTrajectory, simulate_trait_linked_genomes,
)
from zombi2.tree import Tree, TreeNode


def near_star_tree(k: int, age: float, delta: float = 1e-3) -> Tree:
    """A balanced binary tree of ``2**k`` tips whose internal nodes all sit within
    ``[0, delta·k]`` of the root — every terminal branch is ~``age`` long, so tips are
    near-independent (minimal shared ancestry)."""
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


def _genome(families):
    ids = IdManager()
    g = UnorderedGenome(ids)
    for fam in families:
        g._add(Gene(ids.new_gene(), fam))
    return g


# ── 1. trajectory ────────────────────────────────────────────────────────────
def test_trajectory_discrete_is_exact_map():
    tree = simulate_species_tree(BirthDeath(1.0, 0.3), n_tips=25, age=4.0, seed=1)
    result = simulate_traits(tree, Mk.equal_rates(2, 0.7), seed=2)
    traj = TraitTrajectory.from_result(result)

    # every interior jump of the stochastic map is a refresh point
    interior_jumps = sum(max(0, len(segs) - 1)
                         for node, segs in result.history.items() if node.parent is not None)
    assert len(traj.refresh_times(0.0, tree.total_age)) == interior_jumps

    # value at a branch's start equals the first state of its history (== parent's end state)
    for node in tree.nodes_preorder():
        if node.parent is None or not result.history.get(node):
            continue
        first_state = result.history[node][0][0]
        assert traj.value(node.name, node.parent.time) == float(first_state)
        break


def test_trajectory_continuous_subsegments_and_interpolates():
    tree = simulate_species_tree(BirthDeath(1.0, 0.3), n_tips=12, age=4.0, seed=3)
    result = simulate_traits(tree, BrownianMotion(sigma2=0.6), seed=4)
    steps = 5
    traj = TraitTrajectory.from_result(result, steps=steps)

    node = next(n for n in tree.nodes_preorder() if n.parent is not None)
    p, c = float(result.node_values[node.parent]), float(result.node_values[node])
    b0, dt = node.parent.time, node.branch_length() / steps
    # midpoint interpolation on the first sub-segment
    mid_val = traj.value(node.name, b0 + dt * 0.5)
    assert mid_val == pytest.approx(p + (c - p) * (0.5 / steps))
    # each non-root branch contributes (steps - 1) interior refresh points
    n_branches = sum(1 for n in tree.nodes_preorder() if n.parent is not None)
    assert len(traj.refresh_times(0.0, tree.total_age)) == n_branches * (steps - 1)


def test_trajectory_rejects_multivariate():
    tree = simulate_species_tree(BirthDeath(1.0, 0.3), n_tips=8, age=3.0, seed=5)
    result = simulate_traits(tree, MultivariateBrownian(R=np.eye(2)), seed=6)
    with pytest.raises(ValueError, match="univariate"):
        TraitTrajectory.from_result(result)


def test_trajectory_state_values_remap():
    tree = simulate_species_tree(BirthDeath(1.0, 0.3), n_tips=10, age=3.0, seed=7)
    result = simulate_traits(tree, Mk.equal_rates(2, 0.5), seed=8)
    traj = TraitTrajectory.from_result(result, state_values=[-1.0, 1.0])
    # every value is one of the remapped numbers
    for node in tree.extant_leaves():
        assert traj.value(node.name, node.time - 1e-9) in (-1.0, 1.0)


# ── 2. coupling spec ─────────────────────────────────────────────────────────
def test_build_count_fraction_ids():
    c_int = TraitGeneCoupling.build(20, 5, weight=1.0, seed=0)
    assert c_int.n_responsive == 5 and set(c_int.weights) == {0.0, 1.0}

    c_frac = TraitGeneCoupling.build(20, 0.25, weight=1.0, seed=0)
    assert c_frac.n_responsive == 5

    c_ids = TraitGeneCoupling.build(20, ["F1", "F3", 7], weight=2.0, seed=0)
    assert sorted(c_ids.responsive_ids) == ["F1", "F3", "F7"]
    assert c_ids.weights[1] == 2.0 and c_ids.weights[3] == 2.0 and c_ids.weights[7] == 2.0


def test_build_signed_mixes_directions():
    c = TraitGeneCoupling.build(200, 1.0, weight=1.0, signed=True, seed=1)
    signs = np.sign(c.weights)
    assert (signs > 0).any() and (signs < 0).any()          # both directions present


def test_build_rejects_out_of_range_index():
    with pytest.raises(ValueError, match="out of range"):
        TraitGeneCoupling.build(5, ["F9"])


# ── 3. rate-model formula ────────────────────────────────────────────────────
def test_loss_modulation_formula():
    # panel F0 (responsive, w=+1) and F1 (inert); constant trait s via a default-only trajectory
    coupling = TraitGeneCoupling(n_families=2, weights=np.array([1.0, 0.0]),
                                 effect_loss=2.0, base_loss=1.5, transfer=0.4)
    s = 0.75
    traj = TraitTrajectory(starts={}, vals={}, boundaries=[], default=s)
    rates = TraitLinkedRates(coupling, traj)

    weights = {(ew.event, ew.family): ew.rate
               for ew in rates.event_weights(_genome(["F0", "F1"]), "b", 0.0)}
    assert weights[(EventType.LOSS, "F0")] == pytest.approx(1.5 * math.exp(-2.0 * 1.0 * s))
    assert weights[(EventType.LOSS, "F1")] == pytest.approx(1.5)                    # inert = base
    assert weights[(EventType.TRANSFER, None)] == pytest.approx(0.4 * 2)            # field-blind


def test_effect_gain_scales_transfer():
    coupling = TraitGeneCoupling(n_families=1, weights=np.array([1.0]),
                                 effect_gain=0.5, transfer=1.0, base_loss=1.0)
    traj = TraitTrajectory(starts={}, vals={}, boundaries=[], default=2.0)
    rates = TraitLinkedRates(coupling, traj)
    tr = next(ew.rate for ew in rates.event_weights(_genome(["F0"]), "b", 0.0)
              if ew.event is EventType.TRANSFER)
    assert tr == pytest.approx(1.0 * 1 * math.exp(0.5 * 2.0))


# ── 4. ground-truth recovery ─────────────────────────────────────────────────
def _painted_two_clade_trait(tree: Tree) -> TraitResult:
    """A deterministic continuous trait: one root subtree at +1 (favoured), the other at -1."""
    values = {tree.root: 1.0}

    def paint(node, val):
        values[node] = val
        for ch in node.children:
            paint(ch, val)

    c0, c1 = tree.root.children
    paint(c0, 1.0)
    paint(c1, -1.0)
    return TraitResult(tree=tree, model=None, node_values=values, history=None, kind="continuous")


def test_inject_recover_trait_tracks_responsive_families():
    tree = near_star_tree(6, age=6.0)                       # 64 near-independent tips
    trait = _painted_two_clade_trait(tree)
    # tip -> clade sign from the painted trait
    sign = {n.name: trait.node_values[n] for n in tree.extant_leaves()}

    coupling = TraitGeneCoupling.build(40, 0.5, weight=1.0, effect_loss=3.0,
                                       base_loss=1.0, transfer=0.4, seed=9)
    res = simulate_trait_linked_genomes(tree, trait, coupling, seed=13)

    pres = res.profiles.presence()                         # (40, n_species) 0/1
    order = res.profiles.species
    aer_cols = [j for j, s in enumerate(order) if sign[s] > 0]
    ana_cols = [j for j, s in enumerate(order) if sign[s] < 0]
    resp = [coupling.index[i] for i in coupling.responsive_ids]
    inert = [i for i in range(coupling.n_families) if i not in set(resp)]

    def mean(rows, cols):
        return pres[np.ix_(rows, cols)].mean()

    # responsive families are retained where the trait favours them, purged where it does not
    assert mean(resp, aer_cols) > 0.6
    assert mean(resp, ana_cols) < 0.2
    # inert families do not distinguish the two clades (no spurious trait signal)
    assert abs(mean(inert, aer_cols) - mean(inert, ana_cols)) < 0.2


# ── 5. driver / backward-compat ──────────────────────────────────────────────
def test_driver_shape_and_reproducible():
    tree = simulate_species_tree(BirthDeath(1.0, 0.3), n_tips=20, age=4.0, seed=2)
    mk = Mk.equal_rates(2, 0.5)
    coupling = TraitGeneCoupling.build(16, 0.5, effect_loss=2.0, seed=3)

    a = simulate_trait_linked_genomes(tree, mk, coupling, seed=100)
    b = simulate_trait_linked_genomes(tree, mk, coupling, seed=100)
    assert a.profiles.shape == (16, len(tree.extant_leaves()))
    assert np.array_equal(a.profiles.matrix, b.profiles.matrix)   # same seed → identical
    assert a.genomes().profiles.shape == a.profiles.shape          # promotes to Genomes


def test_zero_effect_is_uncoupled():
    # effect_loss=0 → every family evolves at the base rates (no trait dependence)
    tree = simulate_species_tree(BirthDeath(1.0, 0.3), n_tips=30, age=4.0, seed=4)
    trait = simulate_traits(tree, Mk.equal_rates(2, 0.6), seed=4)
    coupling = TraitGeneCoupling.build(30, 0.5, weight=1.0, effect_loss=0.0,
                                       base_loss=0.6, transfer=0.5, seed=5)
    res = simulate_trait_linked_genomes(tree, trait, coupling, seed=6)
    # responsive and inert families have indistinguishable prevalence (coupling is off)
    pres = res.profiles.presence()
    resp = [coupling.index[i] for i in coupling.responsive_ids]
    inert = [i for i in range(30) if i not in set(resp)]
    assert abs(pres[resp].mean() - pres[inert].mean()) < 0.25


def test_default_rate_model_has_no_refresh_times():
    # backward-compat: an ordinary rate model contributes no breakpoints
    from zombi2 import UniformRates
    assert UniformRates(0.1, 0.1, 0.1).refresh_times(0.0, 1.0) == []


# ── 6. CLI ───────────────────────────────────────────────────────────────────
def test_cli_coevolve_smoke(tmp_path):
    from zombi2.cli import main

    tree = simulate_species_tree(BirthDeath(1.0, 0.3), n_tips=20, age=4.0, seed=1)
    tpath = tmp_path / "sp.nwk"
    tpath.write_text(tree.to_newick() + "\n")
    out = tmp_path / "co"
    rc = main(["coevolve", "--couple", "traits:genes", "-t", str(tpath), "--trait-model", "mk",
               "--states", "2", "--trait-center", "--panel", "24", "--responsive", "0.5",
               "--effect-loss", "3", "--write", "profiles", "trees", "--seed", "7", "-o", str(out)])
    assert rc == 0
    for name in ("Profiles.tsv", "Presence.tsv", "traits.tsv", "trait_tree.nwk",
                 "coupling.tsv", "coevolve.log", "species_tree.nwk"):
        assert (out / name).exists(), name
    assert (out / "gene_trees").is_dir()
    # the manifest records one weight per panel family
    body = [ln for ln in (out / "coupling.tsv").read_text().splitlines()
            if ln and not ln.startswith("#") and not ln.startswith("family")]
    assert len(body) == 24
