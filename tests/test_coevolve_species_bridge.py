"""Tests for the into-species bridge (:mod:`zombi2.coevolve.species_bridge`).

The traits:species edge run through the grammar: a discrete trait drives species diversification.
Checked by (1) the MuSSE construction from grammar responses, (2) equivalence to a hand-built MuSSE,
(3) the SSE signal — the fast-speciation state is enriched among the extant tips, and (4) the grammar
classifies the edge as fuse / tree-growing.
"""

import numpy as np
import pytest

from zombi2.coevolve.grammar import CouplingGraph, Jump, Scalar, Table, couple
from zombi2.coevolve.species_bridge import (
    musse_from_responses, simulate_gene_driven_diversification,
    simulate_trait_driven_diversification,
)
from zombi2.coevolve.sse import MuSSE

_Q2 = [[0.0, 0.1], [0.1, 0.0]]                      # 2-state trait transition (MuSSE fills the diagonal)


def test_musse_built_from_table_responses():
    model = musse_from_responses([0, 1], _Q2,
                                 speciation=Table({0: 1.0, 1: 3.0}),
                                 extinction=Table({0: 0.2, 1: 0.2}))
    assert list(model.lambdas) == [1.0, 3.0]        # per-state birth = the Table's values
    assert list(model.mus) == [0.2, 0.2]
    assert model.k == 2


def test_matches_a_hand_built_musse():
    bridged = musse_from_responses([0, 1], _Q2, Table({0: 1.0, 1: 3.0}), Table({0: 0.2, 1: 0.2}))
    direct = MuSSE(birth=[1.0, 3.0], death=[0.2, 0.2], Q=_Q2, states=[0, 1])
    assert list(bridged.lambdas) == list(direct.lambdas)
    assert list(bridged.mus) == list(direct.mus)
    assert np.allclose(bridged.Q, direct.Q)


def test_fast_speciation_state_is_enriched_at_the_tips():
    # state 1 speciates far faster than state 0; starting in state 0, SSE should enrich state 1
    # among the survivors (the classic state-dependent-diversification signal).
    res = simulate_trait_driven_diversification(
        [0, 1], [[0.0, 0.15], [0.15, 0.0]],
        speciation=Table({0: 0.5, 1: 4.0}), extinction=Table({0: 0.2, 1: 0.2}),
        n_tips=300, root_state=0, seed=1)
    states = list(res.values.values())
    frac_fast = sum(s == 1 for s in states) / len(states)
    assert frac_fast > 0.5                          # the fast-speciation state dominates the tips


def test_reproducible_under_the_same_seed():
    kw = dict(speciation=Table({0: 1.0, 1: 3.0}), extinction=Table({0: 0.2, 1: 0.2}),
              n_tips=80, root_state=0)
    a = simulate_trait_driven_diversification([0, 1], _Q2, seed=7, **kw)
    b = simulate_trait_driven_diversification([0, 1], _Q2, seed=7, **kw)
    # res.values is keyed by TreeNode objects (distinct per run) — compare by stable node name
    assert {n.name: s for n, s in a.values.items()} == {n.name: s for n, s in b.values.items()}


def test_grammar_classifies_traits_species_as_fuse_and_tree_growing():
    # the bridge realizes an edge the grammar routes to the fuse / into-species path
    edge = couple("traits", "species", "speciation", 1.0)
    g = CouplingGraph([edge])
    assert g.grows_tree
    assert g.is_fused(edge)
    assert g.mode == "bidirectional"


# ── genes:species (key innovation) ────────────────────────────────────────────
def test_gene_driven_diversification_requires_a_scalar_response():
    # genes:species is an exp-link — a free per-state Table is not a valid response here
    with pytest.raises(TypeError, match="exp-link"):
        simulate_gene_driven_diversification(2, speciation=Table({0: 1.0, 1: 3.0}), n_tips=20)


def test_key_innovation_driver_is_enriched_at_the_tips():
    # driver 0 starts present (root_drivers=1) and boosts speciation (β=1.6) → it radiates and
    # dominates the tips (the key-innovation signal).
    res = simulate_gene_driven_diversification(
        2, speciation=Scalar(1.6), root_drivers=1, transfer=0.6, loss=0.1, n_tips=200, seed=1)
    prev = res.tip_prevalence()
    assert prev[0] > 0.8                                   # the boosted key-innovation driver


def test_gene_driven_diversification_reproducible():
    kw = dict(speciation=Scalar(1.2), root_drivers=1, n_tips=60)
    a = simulate_gene_driven_diversification(2, seed=3, **kw)
    b = simulate_gene_driven_diversification(2, seed=3, **kw)
    assert a.tip_prevalence() == b.tip_prevalence()       # deterministic given the seed


def test_grammar_classifies_genes_species_as_fuse_and_tree_growing():
    edge = couple("genomes", "species", "speciation", 1.0)   # genomes → species (into-species)
    g = CouplingGraph([edge])
    assert g.grows_tree and g.is_fused(edge)


# ── species:traits / species:genomes (cladogenetic overlays) ──────────────────
import zombi2 as z
from zombi2.coevolve.species_bridge import (
    simulate_cladogenetic_genomes, simulate_cladogenetic_trait,
)


def _species_tree(tips=60, seed=1):
    return z.simulate_species_tree(z.BirthDeath(1, 0.3), n_tips=tips, age=5, seed=seed)


def test_cladogenetic_trait_jump_adds_tip_variance():
    tree = _species_tree()
    no_jump = simulate_cladogenetic_trait(tree, z.BrownianMotion(sigma2=0.05), Jump(), seed=1)
    with_jump = simulate_cladogenetic_trait(tree, z.BrownianMotion(sigma2=0.05), Jump(scale=3.0), seed=1)
    var_no = np.var(list(no_jump.values.values()))
    var_yes = np.var(list(with_jump.values.values()))
    assert var_yes > var_no                               # cladogenetic jumps spread the tips


def test_cladogenetic_genomes_runs_and_is_reproducible():
    tree = _species_tree(tips=40)
    jump = Jump(probability=0.12, gain=2.0)
    a = simulate_cladogenetic_genomes(tree, jump, initial_families=30, seed=2)
    b = simulate_cladogenetic_genomes(tree, jump, initial_families=30, seed=2)
    assert a.profile_matrix().to_tsv() == b.profile_matrix().to_tsv()   # deterministic given seed


def test_grammar_classifies_species_traits_as_an_overlay():
    # the reverse direction (species → trait) does NOT grow the tree; it overlays a given one
    edge = couple("species", "traits", "value", 0.3, driver_kind="event")
    g = CouplingGraph([edge])
    assert not g.grows_tree
    assert not g.is_fused(edge)
    assert g.mode == "directional"


# ── joints: ClaSSE (traits↔species) and co-diversification (genes↔species) ─────
def test_classe_joint_jumps_the_trait_at_speciation():
    kw = dict(speciation=Table({0: 1.0, 1: 2.5}), extinction=Table({0: 0.2, 1: 0.2}),
              n_tips=120, root_state=0)
    plain = simulate_trait_driven_diversification([0, 1], _Q2, seed=4, **kw)
    classe = simulate_trait_driven_diversification([0, 1], _Q2, cladogenesis=Jump(probability=0.4),
                                                   seed=4, **kw)
    assert len(plain.values) == len(classe.values) == 120     # both grow the tree
    # ClaSSE also shifts the state at each split → the seeded runs diverge
    assert ({n.name: s for n, s in plain.values.items()}
            != {n.name: s for n, s in classe.values.items()})


def test_co_diversification_joint_runs_with_a_founder_burst():
    res = simulate_gene_driven_diversification(
        3, speciation=Scalar(1.2), root_drivers=1,
        cladogenesis=Jump(probability=0.1, gain=1.0), n_tips=120, seed=5)
    assert len(res.tree.extant_leaves()) == 120
    assert len(res.tip_prevalence()) == 3                     # the joint reshuffles the driver panel
