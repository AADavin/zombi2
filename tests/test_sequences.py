"""Sequences level, slice 1: a nucleotide substitution model evolving along the gene trees under a
strict clock, endpoint P-matrix sampling → SequencesResult(.alignments, .ancestral)."""

from __future__ import annotations

import json
import re

import numpy as np
import pytest

from zombi2 import species
from zombi2.genomes import FamilyGenomesResult, simulate_genomes_family
from zombi2.genomes.events import copy_label
from zombi2.genomes.gene_trees import GeneNode, GeneTree
from zombi2.rates import modifiers as mod
from zombi2.sequences import SequencesResult, simulate_sequences
from zombi2.sequences.substitution_models import (AMINO_ACIDS, SubstitutionModel, dayhoff, gtr,
                                                  hky85, jc69, jtt, k80, lg, poisson, reversible,
                                                  wag)


# --- hand-built gene trees: origination → speciation → two extant tips -----------------------------

def _run(gene_trees, *, t_split: float = 1.0, t_now: float = 2.0) -> FamilyGenomesResult:
    """The **genome run** the sequence level requires, wrapped around hand-built gene trees: a
    three-lineage species tree (root 0 splits at ``t_split`` into extant tips 1 and 2 at ``t_now``)
    carrying the given families. The gene trees are attached directly — these families are written by
    hand, so there is no event log for the run to derive them from."""
    tree = species.Tree({0: species.Node(0, None, 0.0, t_split, (1, 2), "speciation"),
                         1: species.Node(1, 0, t_split, t_now, None, "extant"),
                         2: species.Node(2, 0, t_split, t_now, None, "extant")}, 0)
    run = FamilyGenomesResult(complete_tree=tree, genomes={}, edges=[], seed=None)
    run.gene_trees = dict(gene_trees)      # a cached_property: the instance dict wins
    return run


def _pair_run(t_spec: float, t_tip: float) -> FamilyGenomesResult:
    """A minimal family (per-segment gene ids) in its genome run: the founding gene (id 0) on species
    0 ends by a ``speciation`` at ``t_spec``; its two daughters (ids 1, 2) reach ``extant`` tips
    (species 1, 2) at ``t_tip``. The root→tip branch has length ``t_tip - t_spec`` (so
    ``_pair_run(0.0, d)`` is one branch of length ``d`` from the root gene to each tip)."""
    root = GeneNode("speciation", 0, t_spec, 0)
    root.children = [GeneNode("extant", 1, t_tip, 1), GeneNode("extant", 2, t_tip, 2)]
    return _run({0: GeneTree(0, root, 0.0)}, t_split=t_spec, t_now=t_tip)


def _one_lineage(family: int, lineage: int, t_tip: float) -> GeneTree:
    """A family whose founding gene (id 0) on species ``lineage`` ends by a ``speciation`` at t=0,
    both daughters (ids 1, 2) staying on ``lineage`` to ``extant`` tips at ``t_tip`` — so both
    branches ride the *same* species lineage's clock. Used to test that the lineage clock is shared
    across families."""
    root = GeneNode("speciation", lineage, 0.0, 0)
    root.children = [GeneNode("extant", lineage, t_tip, 1), GeneNode("extant", lineage, t_tip, 2)]
    return GeneTree(family, root, 0.0)


def _iter_nodes(root):
    stack = [root]
    while stack:
        n = stack.pop()
        yield n
        stack.extend(n.children)


def _seqs(result: SequencesResult, fam: int = 0) -> list[str]:
    return list(result.alignments[fam].values()) + list(result.ancestral[fam].values())


# --- substitution models (the menu) ----------------------------------------------------------------

def test_p_matrix_is_a_valid_reversible_transition():
    models = [jc69(), k80(2.0), hky85(2.0, (0.1, 0.2, 0.3, 0.4)),
              gtr((1, 2, 1, 1, 2, 1), (0.15, 0.25, 0.3, 0.3))]
    for m in models:
        P = m.p_matrix(0.37)
        assert np.allclose(P.sum(1), 1.0)          # rows are distributions
        assert (P >= 0).all()
        assert np.allclose(m.stationary @ P, m.stationary)   # stationary is preserved
    assert np.allclose(jc69().p_matrix(0.0), np.eye(4))       # P(0) = identity


def test_models_are_normalised_to_one_substitution_per_unit_time():
    # -Σ π_i Q_ii == 1 (branch lengths are in substitutions/site)
    for m in (jc69(), k80(3.0), hky85(2.0, (0.2, 0.3, 0.3, 0.2)), gtr()):
        assert np.isclose(-(m.stationary * np.diag(m.Q)).sum(), 1.0)


# --- a matrix of your own: reversible() ------------------------------------------------------------

def _hky_exchangeabilities(kappa: float) -> np.ndarray:
    """HKY85's 4×4 exchangeability matrix written out by hand: transitions (A↔G, C↔T) at ``kappa``,
    transversions at 1. Built here rather than taken from the module so the test compares the public
    constructor against the menu, not the module against itself."""
    return np.array([[0.0, 1.0, kappa, 1.0],
                     [1.0, 0.0, 1.0, kappa],
                     [kappa, 1.0, 0.0, 1.0],
                     [1.0, kappa, 1.0, 0.0]])


def test_reversible_rebuilds_a_menu_model_exactly():
    """The public door is the same door the menu goes through: a hand-written HKY85 matrix and the
    same frequencies give the model `hky85()` gives, cell for cell — including the normalisation to
    one expected substitution per site per unit branch length."""
    frequencies = (0.3, 0.2, 0.2, 0.3)
    mine = reversible(_hky_exchangeabilities(2.0), frequencies, name="HKY85")
    theirs = hky85(kappa=2.0, frequencies=frequencies)
    assert np.allclose(mine.Q, theirs.Q)
    assert np.allclose(mine.stationary, theirs.stationary)
    assert np.allclose(mine.p_matrix(0.42), theirs.p_matrix(0.42))
    assert mine.alphabet == theirs.alphabet == "ACGT"


def test_a_custom_model_is_normalised_and_keeps_its_stationary():
    """The two properties the rest of the level leans on. Normalisation is what makes a branch
    length mean substitutions per site whatever matrix produced it — the reason one phylogram is
    exact for a partitioned run — and stationarity is what makes the founding draw a fixed point."""
    rng = np.random.default_rng(3)
    S = rng.random((4, 4))
    S = S + S.T
    np.fill_diagonal(S, 0.0)
    pi = rng.random(4) + 0.1
    pi = pi / pi.sum()
    m = reversible(S, pi, name="Random")
    assert np.isclose(-(m.stationary * np.diag(m.Q)).sum(), 1.0)
    assert np.allclose(m.stationary @ m.p_matrix(0.4), m.stationary)
    assert np.allclose(m.p_matrix(0.4).sum(1), 1.0)


def test_a_custom_model_can_use_the_amino_acid_alphabet():
    """`reversible()` is not a nucleotide constructor with a wider door — `poisson()` is literally
    this call, and a run under a hand-built 20-state matrix writes amino acids."""
    S = np.ones((20, 20)) - np.eye(20)
    mine = reversible(S, np.full(20, 0.05), name="Poisson", alphabet=AMINO_ACIDS)
    assert np.allclose(mine.Q, poisson().Q)
    r = simulate_sequences(_pair_run(1.0, 2.0), model=mine, length=60, seed=2)
    assert set("".join(_seqs(r))) <= set(AMINO_ACIDS)


#: bad exchangeability matrices, and the *detail* the message has to name — not merely "invalid".
_BAD_MATRICES = [
    (np.array([[0.0, 1.0, 1.0, 1.0], [2.0, 0.0, 1.0, 1.0],
               [1.0, 1.0, 0.0, 1.0], [1.0, 1.0, 1.0, 0.0]]), "asymmetric by up to 1"),
    (np.full((4, 4), -1.0) + np.eye(4), "minimum entry -1"),
    (np.ones((4, 4)), "zero diagonal"),
    (np.ones((3, 3)) - np.eye(3), r"K = 4"),                   # a 3×3 S against 4 frequencies
]


@pytest.mark.parametrize("S, match", _BAD_MATRICES)
def test_reversible_names_what_is_wrong_with_a_bad_matrix(S, match):
    with pytest.raises(ValueError, match=match):
        reversible(S, [0.25] * 4)


@pytest.mark.parametrize("frequencies, match", [
    ([0.4, 0.3, 0.2, 0.0], "strictly positive"),
    ([0.4, 0.3, 0.2, 0.2], "sum to 1"),
])
def test_reversible_names_what_is_wrong_with_the_frequencies(frequencies, match):
    with pytest.raises(ValueError, match=match):
        reversible(np.ones((4, 4)) - np.eye(4), frequencies)


@pytest.mark.parametrize("alphabet, match", [
    ("ACG", "characters but"),                                 # too few for a 4×4
    ("ACGA", "distinct"),                                      # a repeat: two states, one letter
    ("ACGé", "distinct"),                                      # non-ASCII cannot go in a FASTA
])
def test_a_custom_alphabet_must_name_every_state_exactly_once(alphabet, match):
    with pytest.raises(ValueError, match=match):
        reversible(np.ones((4, 4)) - np.eye(4), [0.25] * 4, alphabet=alphabet)


def test_a_non_reversible_rate_matrix_is_refused_rather_than_silently_wrong():
    """The load-bearing guard. `p_matrix()` eigendecomposes diag(√π)·Q·diag(1/√π), which is similar
    to Q only under detailed balance, so a cyclic A→C→G→T→A matrix — a perfectly good generator,
    just not a reversible one — would come back as a well-formed transition matrix for a *different*
    model. Nothing downstream could notice: the sequences would look like sequences."""
    Q = np.array([[-1.0, 1.0, 0.0, 0.0],
                  [0.0, -1.0, 1.0, 0.0],
                  [0.0, 0.0, -1.0, 1.0],
                  [1.0, 0.0, 0.0, -1.0]])
    assert np.allclose(Q.sum(1), 0.0)                          # it *is* a rate matrix
    with pytest.raises(ValueError, match="detailed balance") as e:
        SubstitutionModel("Cyclic", Q, np.full(4, 0.25))
    assert "not implemented" in str(e.value)                   # a fact about the code, not the model
    assert "reversible(" in str(e.value)                       # and where the thing does belong


def test_a_rate_matrix_whose_rows_do_not_sum_to_zero_is_refused():
    Q = np.array([[-1.0, 1.0, 0.0, 0.5],
                  [1.0, -1.0, 0.0, 0.0],
                  [0.0, 0.0, 0.0, 0.0],
                  [0.5, 0.0, 0.0, -0.5]])
    with pytest.raises(ValueError, match="sum to 0"):
        SubstitutionModel("Leaky", Q, np.full(4, 0.25))


@pytest.mark.parametrize("build", [jc69, lambda: k80(3.0), lambda: hky85(2.0, (0.2, 0.3, 0.3, 0.2)),
                                   gtr, poisson, jtt, dayhoff, wag, lg])
def test_every_menu_model_passes_the_reversibility_guard(build):
    """A tolerance regression guard. The published 20-state matrices are symmetric only to
    round-off, so this pins that the check that refuses a non-reversible model does not also refuse
    WAG."""
    m = build()
    flux = m.stationary[:, None] * m.Q
    assert np.allclose(flux, flux.T)


def test_across_sites_survives_the_reversibility_guard():
    """`across_sites()` rebuilds the model through `dataclasses.replace`, which re-runs
    ``__post_init__`` — so the guard sees every decorated model too, and must pass them."""
    m = hky85(kappa=2.0).across_sites(gamma_shape=0.5, invariant=0.1)
    assert m.name == "HKY85+I+G4"


# --- the engine: determinism, the strict clock, structure ------------------------------------------

def test_determinism_same_seed_identical_different_seed_differs():
    run = _pair_run(1.0, 2.0)
    a = simulate_sequences(run, model=jc69(), length=200, seed=7)
    b = simulate_sequences(run, model=jc69(), length=200, seed=7)
    assert a.alignments == b.alignments and a.ancestral == b.ancestral
    c = simulate_sequences(run, model=jc69(), length=200, seed=8)
    assert a.alignments != c.alignments


def test_zero_rate_leaves_every_node_equal_to_the_root():
    # substitution = 0 → every branch length is 0 → no site ever changes
    r = simulate_sequences(_pair_run(1.0, 2.0), model=hky85(kappa=3.0), length=150,
                           substitution=0.0, seed=1)
    assert len(set(_seqs(r))) == 1


def test_zero_length_branch_copies_its_parent():
    # the root gene (id 0) ends by a duplication at t=1; its daughter id 1 is an extant tip AT t=1 — a
    # zero-length branch — so it must copy the root's sequence
    root = GeneNode("duplication", 0, 1.0, 0)
    root.children = [GeneNode("extant", 1, 1.0, 1), GeneNode("extant", 2, 2.0, 2)]
    r = simulate_sequences(_run({0: GeneTree(0, root, 0.0)}), model=jc69(), length=120, seed=3)
    assert r.alignments[0]["n1_g1"] == r.ancestral[0]["n0_g0"]


def test_every_sequence_has_the_requested_length_over_the_alphabet():
    r = simulate_sequences(_pair_run(1.0, 2.0), model=gtr(), length=123, seed=1)
    for seq in _seqs(r):
        assert len(seq) == 123
        assert set(seq) <= set("ACGT")


def test_alignment_labels_are_exactly_the_extant_tips():
    r = simulate_sequences(_pair_run(1.0, 2.0), model=k80(2.0), length=10, seed=1)
    assert set(r.alignments[0]) == {"n1_g1", "n2_g2"}  # the two extant daughter genes: species 1, 2


def test_jc69_holds_uniform_base_composition():
    # root drawn from the uniform stationary; JC69 keeps it uniform, so tips stay ≈ 25% each
    r = simulate_sequences(_pair_run(1.0, 5.0), model=jc69(), length=20000, seed=42)
    seq = r.alignments[0]["n1_g1"]
    for base in "ACGT":
        assert abs(seq.count(base) / len(seq) - 0.25) < 0.03


def test_jc69_pdistance_matches_theory_and_rate_scales_it():
    # the endpoint distribution must match the JC69 process: p(d) = 3/4 (1 - exp(-4d/3)), with the
    # branch length d = substitution · Δt. Guards the model normalisation and the time→subs/site
    # conversion — a structural test would miss a rescaled or mis-normalised rate matrix.
    def pdist(a, b):
        return np.mean(np.frombuffer(a.encode(), np.uint8) != np.frombuffer(b.encode(), np.uint8))

    def root_tip_pdistance(*, t_tip, substitution):
        r = simulate_sequences(_pair_run(0.0, t_tip), model=jc69(), length=40000,
                               substitution=substitution, seed=1)
        # root gene n0_g0 → daughter tip n1_g1
        return pdist(r.ancestral[0]["n0_g0"], r.alignments[0]["n1_g1"])

    theory = lambda d: 0.75 * (1 - np.exp(-4 * d / 3))          # noqa: E731
    assert abs(root_tip_pdistance(t_tip=1.0, substitution=1.0) - theory(1.0)) < 0.01
    # substitution = 0.5 halves the effective distance: Δt=2 behaves like d=1.0
    assert abs(root_tip_pdistance(t_tip=2.0, substitution=0.5) - theory(1.0)) < 0.01


def test_hky85_transition_bias_makes_diverged_tips_still_reflect_frequencies():
    # a strongly skewed base composition is reproduced at the tips (endpoint stays near stationary)
    frequencies = (0.4, 0.1, 0.1, 0.4)
    r = simulate_sequences(_pair_run(1.0, 4.0), model=hky85(4.0, frequencies), length=20000, seed=5)
    seq = r.alignments[0]["n1_g1"]
    comp = [seq.count(b) / len(seq) for b in "ACGT"]
    assert comp[0] > comp[1] and comp[3] > comp[2]   # A,T (0.4) exceed C,G (0.1)


# --- a family with no surviving copy ---------------------------------------------------------------

def test_family_with_no_extant_copy_has_empty_alignment_but_full_ancestral():
    # the root gene (id 0) speciates, but both daughters are lost → nothing observable, yet all three
    # are real nodes of the tree with a reconstructed sequence apiece, the two dead tips included
    root = GeneNode("speciation", 0, 1.0, 0)
    root.children = [GeneNode("loss", 0, 2.0, 1), GeneNode("loss", 0, 2.0, 2)]
    r = simulate_sequences(_run({0: GeneTree(0, root, 0.0)}), model=jc69(), length=10, seed=1)
    assert r.alignments[0] == {}
    assert set(r.ancestral[0]) == {"n0_g0", "n0_g1", "n0_g2"}    # all three still on species 0


# --- integration: species → genomes → sequences ----------------------------------------------------

def test_a_real_genome_run_is_covered_node_for_node():
    sp = species.simulate_species_tree(birth=1.0, death=0.2, n_extant=8, seed=1)
    g = simulate_genomes_family(sp, duplication=0.2, loss=0.2, transfer=0.1,
                                initial_families=6, seed=2)
    r = simulate_sequences(g, model=hky85(kappa=2.5), length=300, seed=3)   # the genome run itself
    assert isinstance(r, SequencesResult)
    assert set(r.alignments) == set(g.gene_trees) == set(r.ancestral)       # one entry per family
    for fam, gt in g.gene_trees.items():
        nodes = list(_iter_nodes(gt.complete))
        n_extant = sum(1 for n in nodes if n.is_leaf and n.kind == "extant")
        assert len(r.alignments[fam]) == n_extant
        # everything that is not an extant tip: internal nodes and the tips where a copy or its
        # species died. Together they account for every node in the tree, exactly once.
        assert len(r.ancestral[fam]) == len(nodes) - n_extant
        for seq in list(r.alignments[fam].values()) + list(r.ancestral[fam].values()):
            assert len(seq) == 300


def test_integration_is_deterministic_given_the_seed():
    sp = species.simulate_species_tree(birth=1.0, death=0.3, n_extant=6, seed=11)
    g = simulate_genomes_family(sp, duplication=0.3, loss=0.2, initial_families=5, seed=12)
    a = simulate_sequences(g, model=gtr(), length=100, seed=13)
    b = simulate_sequences(g, model=gtr(), length=100, seed=13)
    assert a.alignments == b.alignments and a.ancestral == b.ancestral


# --- the lineage clock (ByLineage): the uncorrelated / relaxed clock --------------------------------

def _pdist(a: str, b: str) -> float:
    return sum(x != y for x, y in zip(a, b)) / len(a)


def test_bylineage_zero_spread_is_bit_identical_to_the_strict_clock():
    # spread=0 draws 1.0 without touching the rng, so the run matches the strict clock exactly
    run = _pair_run(1.0, 2.0)
    strict = simulate_sequences(run, model=jc69(), length=300, seed=5)
    clocked = simulate_sequences(run, model=jc69(), length=300,
                                 substitution=1.0 * mod.ByLineage(spread=0.0), seed=5)
    assert clocked.alignments == strict.alignments and clocked.ancestral == strict.ancestral


def test_bylineage_perturbs_the_output_and_stays_valid():
    run = _pair_run(1.0, 2.0)
    strict = simulate_sequences(run, model=jc69(), length=300, seed=5)
    clocked = simulate_sequences(run, model=jc69(), length=300,
                                 substitution=1.0 * mod.ByLineage(spread=0.5), seed=5)
    assert clocked.alignments != strict.alignments          # the clock rescales branch lengths
    for seq in _seqs(clocked):
        assert len(seq) == 300 and set(seq) <= set("ACGT")


def test_bylineage_is_deterministic():
    run = _pair_run(1.0, 2.0)
    spec = 1.0 * mod.ByLineage(spread=0.4)
    a = simulate_sequences(run, model=hky85(2.0), length=200, substitution=spec, seed=9)
    b = simulate_sequences(run, model=hky85(2.0), length=200, substitution=spec, seed=9)
    assert a.alignments == b.alignments and a.ancestral == b.ancestral


def test_bylineage_clock_is_shared_across_families_on_a_lineage():
    # 20 identical single-branch families, ALL on species lineage 0 → all feel the SAME clock[0].
    # A per-family clock would scatter their root→tip divergences by ~spread; a shared clock leaves
    # only sampling noise, so the across-family spread collapses.
    run = _run({f: _one_lineage(f, lineage=0, t_tip=1.0) for f in range(20)}, t_split=0.0, t_now=1.0)
    r = simulate_sequences(run, model=jc69(), length=5000,
                           substitution=1.0 * mod.ByLineage(spread=0.8), seed=4)
    ds = [_pdist(r.ancestral[f]["n0_g0"], r.alignments[f]["n0_g1"]) for f in range(20)]
    mean = sum(ds) / len(ds)
    std = (sum((d - mean) ** 2 for d in ds) / len(ds)) ** 0.5
    assert std < 0.02      # shared clock ⇒ ~0.006 sampling noise; a per-family clock would be far larger


def test_sequence_clock_rejects_multiple_or_unwired_modifiers(tmp_path):
    run = _pair_run(1.0, 2.0)
    with pytest.raises(ValueError, match="lineage clocks"):   # two clocks — a lineage has one
        simulate_sequences(run, model=jc69(), length=10,
                           substitution=1.0 * mod.FromParent(spread=0.3) * mod.ByLineage(spread=0.2))
    with pytest.raises(ValueError, match="lineage clocks"):   # two ByLineage
        simulate_sequences(run, model=jc69(), length=10,
                           substitution=1.0 * mod.ByLineage(spread=0.3) * mod.ByLineage(spread=0.2))
    with pytest.raises(ValueError, match="OnTime"):  # ByLineage × OnTime — a modifier this level
        simulate_sequences(run, model=jc69(), length=10,   # does not read
                           substitution=1.0 * mod.ByLineage(spread=0.3) * mod.OnTime({0: 1.0}))
    # A clock and a driver, though, are two different axes and DO compose (SPEC §5: modifiers
    # multiply) — so this one has to run rather than raise.
    driver = tmp_path / "d.tsv"
    driver.write_text("time\tkind\tlineage\tfrom\tto\n0.0\tinitial\tn0\t\ta\n", encoding="utf-8")
    simulate_sequences(run, model=jc69(), length=10, seed=1,
                       substitution=1.0 * mod.ByLineage(spread=0.3)
                       * mod.DrivenBy(str(driver), {"a": 2.0}))


def test_a_between_kernel_is_refused_on_the_substitution_rate(tmp_path):
    """A ``Between`` weights an ordered (donor, recipient) pair, and a rate is read on one lineage —
    there is no donor for the pair's first half to name. Meaningless rather than unimplemented, so the
    message has to name the slot it does belong in (SPEC §5)."""
    from zombi2.rates.mapping import Between

    run = _pair_run(1.0, 2.0)
    driver = tmp_path / "d.tsv"
    driver.write_text("time\tkind\tlineage\tfrom\tto\n0.0\tinitial\tn0\t\ta\n", encoding="utf-8")
    with pytest.raises(ValueError, match="transfer_to"):
        simulate_sequences(run, model=jc69(), length=10,
                           substitution=1.0 * mod.DrivenBy(str(driver),
                                                           Between({("a", "b"): 1.0})))


def test_a_live_level_source_is_refused_as_a_joint_run():
    """SPEC §3: Traits and Sequences can be conditioned and never joined. A live-level source is the
    joint spelling of DrivenBy, so it must come back as the modelling answer rather than as a missing
    file called 'trait'."""
    run = _pair_run(1.0, 2.0)
    with pytest.raises(ValueError, match="cannot be joined"):
        simulate_sequences(run, model=jc69(), length=10,
                           substitution=1.0 * mod.DrivenBy("trait", {"a": 2.0}))


def test_divergence_is_refused_alongside_a_driven_rate():
    """``divergence / height`` is the base only when the modifiers average to 1 along a root-to-tip
    path — which the two clocks are mean-corrected to do and a driver is not. Solving anyway would log
    a divergence the run did not realise, so it is refused with the base to write instead."""
    run = _pair_run(1.0, 2.0)
    with pytest.raises(ValueError, match="divergence"):
        simulate_sequences(run, model=jc69(), length=10, divergence=0.2,
                           substitution=mod.DrivenBy("trait_events.tsv", {"a": 2.0}))


# --- the lineage clock (FromParent): the autocorrelated clock ---------------------------------------

def test_fromparent_perturbs_the_output_and_stays_valid():
    run = _pair_run(1.0, 2.0)
    strict = simulate_sequences(run, model=jc69(), length=300, seed=5)
    clocked = simulate_sequences(run, model=jc69(), length=300,
                                 substitution=1.0 * mod.FromParent(spread=0.5), seed=5)
    assert clocked.alignments != strict.alignments          # the clock rescales branch lengths
    for seq in _seqs(clocked):
        assert len(seq) == 300 and set(seq) <= set("ACGT")


def test_fromparent_is_deterministic():
    run = _pair_run(1.0, 2.0)
    spec = 1.0 * mod.FromParent(spread=0.4)
    a = simulate_sequences(run, model=hky85(2.0), length=200, substitution=spec, seed=9)
    b = simulate_sequences(run, model=hky85(2.0), length=200, substitution=spec, seed=9)
    assert a.alignments == b.alignments and a.ancestral == b.ancestral


# --- phylograms: the gene / species trees in substitutions/site -------------------------------------

def _leaves(nwk: str) -> set[str]:
    # n<species>_g<copy> in leaf position (not after a ')')
    return set(re.findall(r"(?<![)\w])n\d+_g\d+", nwk))


def _total_bl(nwk: str) -> float:
    return sum(float(x) for x in re.findall(r":([0-9.eE+-]+)", nwk))


def _small_run(clock=1.0):
    sp = species.simulate_species_tree(birth=1.0, death=0.2, n_extant=8, seed=1)
    g = simulate_genomes_family(sp, duplication=0.3, loss=0.3, transfer=0.1, initial_families=8, seed=2)
    r = simulate_sequences(g, model=jc69(), length=10, substitution=clock, seed=3)
    return g, r


def test_result_carries_phylograms_and_species_phylogram():
    g, r = _small_run()
    assert set(r.phylograms) == set(g.gene_trees)
    for ph in r.phylograms.values():
        assert set(ph) == {"complete", "extant"} and isinstance(ph["complete"], str)
    assert set(r.species_phylogram) == {"complete", "extant"}


def test_strict_clock_phylogram_matches_the_chronogram_lengths():
    # base 1, strict clock -> subs/site == time, so the phylogram's branch lengths equal the
    # chronogram's (same topology; the phylogram labels every node by gene id, so internal labels
    # differ). The species phylogram keeps the species labels, so there it is byte-identical.
    # The root branch is included: the family's founding sequence evolves across the stem, so the
    # phylogram states it in subs/site just as the chronogram states it in time.
    def bls(nwk):
        return sorted(re.findall(r":([0-9.eE+-]+)", nwk))

    g, r = _small_run(clock=1.0)
    for fam, gt in g.gene_trees.items():
        assert re.search(r":[0-9.eE+-]+;$", r.phylograms[fam]["complete"])   # the root branch is there
        assert bls(r.phylograms[fam]["complete"]) == bls(gt.to_newick("complete"))
        if r.phylograms[fam]["extant"] is not None:
            assert bls(r.phylograms[fam]["extant"]) == bls(gt.to_newick("extant"))
    assert r.species_phylogram["complete"] == g.complete_tree.to_newick()


def test_the_founding_sequence_evolves_across_the_stem():
    # a family exists from its origination, so the sequence it started with is not the sequence its
    # root gene ended with — it evolved across the stem, and the phylogram's root branch is that.
    g, r = _small_run(clock=1.0)
    assert set(r.founding) == set(g.gene_trees)
    for fam, gt in g.gene_trees.items():
        assert len(r.founding[fam]) == len(next(iter(r.ancestral[fam].values()), r.founding[fam]))
        stem = gt.complete.time - gt.origination
        root_seq = r.ancestral[fam].get(copy_label(gt.complete.species, gt.complete.copy))
        if root_seq is not None and stem > 0:
            # not asserted equal or unequal site-by-site (a short stem may fix nothing), but the
            # branch the phylogram reports must be exactly the stem under a rate-1 strict clock
            written = float(r.phylograms[fam]["complete"].rsplit(":", 1)[1].rstrip(";"))
            assert written == pytest.approx(stem, rel=1e-5)


def test_every_phylogram_node_has_a_sequence():
    # the complete phylogram labels every node by its copy label (n<species>_g<copy>), and every
    # one of those names a sequence: the extant tips in the alignment, all the rest — internal
    # nodes and dead tips alike — in the ancestral set. No label points at nothing.
    g, r = _small_run()
    for fam in g.gene_trees:
        internal = set(re.findall(r"\)(n\d+_g\d+)", r.phylograms[fam]["complete"]))
        assert internal <= set(r.ancestral[fam])
        labelled = set(re.findall(r"(n\d+_g\d+)", r.phylograms[fam]["complete"]))
        assert labelled == set(r.ancestral[fam]) | set(r.alignments[fam])


def test_phylogram_extant_leaves_match_the_alignment():
    g, r = _small_run()
    for fam in g.gene_trees:
        e = r.phylograms[fam]["extant"]
        if e is not None:
            assert _leaves(e) == set(r.alignments[fam])


def test_base_rate_scales_the_phylogram_branch_lengths():
    _, r1 = _small_run(clock=1.0)
    _, r2 = _small_run(clock=2.0)
    assert _total_bl(r2.species_phylogram["complete"]) == pytest.approx(
        2 * _total_bl(r1.species_phylogram["complete"]), rel=1e-4)


def test_lineage_clock_reshapes_the_phylograms():
    g, strict = _small_run(clock=1.0)
    _, clocked = _small_run(clock=1.0 * mod.ByLineage(spread=0.7))
    assert any(clocked.phylograms[f]["complete"] != strict.phylograms[f]["complete"] for f in g.gene_trees)


def test_bare_gene_trees_are_rejected():
    # the level runs on a genome run, never on a loose {family: GeneTree}: without the species tree
    # the lineage clock has no branches to ride and the species phylogram cannot exist, so accepting
    # the mapping would hide that degradation instead of naming it
    g, _ = _small_run()
    with pytest.raises(TypeError, match="simulate_genomes_family"):
        simulate_sequences(g.gene_trees, model=jc69(), length=10, seed=3)


def test_write_emits_phylogram_newick(tmp_path):
    _, r = _small_run()
    r.write(tmp_path, outputs=("phylograms", "species_phylogram"))
    fam0 = tmp_path / "phylograms" / "phylogram_fam0_complete.nwk"
    assert fam0.exists() and fam0.read_text(encoding="utf-8").rstrip().endswith(";")
    assert (tmp_path / "clock_species_tree_complete.nwk").exists()


# --- validation ------------------------------------------------------------------------------------

def test_rejects_bad_arguments_and_unwired_rate_specs():
    from zombi2.rates.scope import PerLineage
    run = _pair_run(1.0, 2.0)
    with pytest.raises(TypeError):
        simulate_sequences(run, model="jc69", length=10)                  # not a SubstitutionModel
    with pytest.raises(ValueError):
        simulate_sequences(run, model=jc69(), length=0)                   # non-positive length
    with pytest.raises(ValueError):
        simulate_sequences(run, model=jc69(), length=10, substitution=1.0 * mod.OnTime({0: 1.0}))
    with pytest.raises(ValueError):
        simulate_sequences(run, model=jc69(), length=10, substitution=PerLineage(1.0))
    with pytest.raises(TypeError, match="genome run"):                    # a species run is not one
        simulate_sequences(species.simulate_species_tree(birth=1.0, n_extant=4, seed=1),
                           model=jc69(), length=10)


# --- writing ---------------------------------------------------------------------------------------

def test_write_emits_fasta_per_family(tmp_path):
    r = simulate_sequences(_pair_run(1.0, 2.0), model=jc69(), length=20, seed=1)
    r.write(tmp_path)
    aln = tmp_path / "alignments" / "fam0.fasta"
    assert aln.exists() and ">n1_g1" in aln.read_text(encoding="utf-8")
    r.write(tmp_path, outputs=("ancestral",))
    assert (tmp_path / "ancestral" / "sequences_ancestral_fam0.fasta").exists()
    with pytest.raises(ValueError):
        r.write(tmp_path, outputs=("bogus",))


# --- divergence: state the outcome, let the rate be solved for ---------------

def test_divergence_solves_for_the_rate_from_the_tree_height():
    """The rate is per unit time, so the same number means different things on trees of different
    heights. divergence states the outcome instead, and the base falls out of the height."""
    from zombi2.tree import rescale
    short = rescale(_tree_for_divergence(), height=2.0)
    tall = rescale(_tree_for_divergence(), height=20.0)
    out = []
    for ct in (short, tall):
        g = simulate_genomes_family(ct, initial_families=8, seed=1)
        r = simulate_sequences(g, model=jc69(), length=400, divergence=0.2, seed=7)
        out.append(_mean_identity(r))
    # a tenfold difference in height, and the alignments still come out alike
    assert abs(out[0] - out[1]) < 0.05, out
    assert 0.6 < out[0] < 0.95, out


def test_divergence_composes_with_a_clock_shape_but_refuses_a_base():
    from zombi2.tree import rescale
    ct = rescale(_tree_for_divergence(), height=10.0)
    g = simulate_genomes_family(ct, initial_families=6, seed=1)
    # the shape alone composes: divergence sets the scale of a relaxed clock
    r = simulate_sequences(g, model=jc69(), length=200,
                           divergence=0.2, substitution=mod.ByLineage(spread=0.3), seed=7)
    assert 0.5 < _mean_identity(r) < 0.98
    # a base alongside is refused rather than silently overridden
    for base in (1.0, 0.5, 1.0 * mod.ByLineage(spread=0.3)):
        with pytest.raises(ValueError, match="names a base"):
            simulate_sequences(g, model=jc69(), length=50,
                               divergence=0.2, substitution=base, seed=7)


def _tree_for_divergence():
    return species.simulate_species_tree(birth=1.0, death=0.0, n_extant=10, seed=3).complete_tree


def _mean_identity(result):
    import itertools
    import numpy as np
    tot = match = 0
    for aln in result.alignments.values():
        for a, b in itertools.combinations(list(aln.values()), 2):
            n = min(len(a), len(b))
            A = np.frombuffer(a[:n].encode(), dtype=np.uint8)
            B = np.frombuffer(b[:n].encode(), dtype=np.uint8)
            match += int((A == B).sum()); tot += n
    return match / tot if tot else 0.0


def test_the_reported_identity_is_not_rounded_to_a_whole_percent(tmp_path, capsys):
    # 32% and 32.4% are different alignments, and the whole point of reporting it is that the rate
    # says nothing about what came out
    from zombi2.cli.main import main

    run = tmp_path / "r"
    main(["species", str(run), "--birth", "1", "--n-extant", "6", "--seed", "1", "--quiet"])
    main(["genomes", str(run), "--duplication", "0.2", "--loss", "0.2", "--seed", "1", "--quiet"])
    capsys.readouterr()
    main(["sequences", str(run), "--model", "jc69", "--length", "200", "--substitution", "0.1",
          "--seed", "1", "--quiet"])
    out = capsys.readouterr().out
    import re
    m = re.search(r"mean identity (\d+\.\d)%", out)
    assert m, out


# --- streaming: the same run, written as it goes instead of assembled first -----------------------

def _stream_fixture(seed=42):
    sp = species.simulate_species_tree(birth=1.0, death=0.3, n_extant=15, seed=1)
    return simulate_genomes_family(sp, duplication=0.2, transfer=0.1, loss=0.25, origination=0.5,
                                   initial_families=25, seed=seed)


def _tree_of(directory):
    """Every file under `directory` as {relative path: bytes} — what a streamed run must reproduce."""
    return {str(p.relative_to(directory)): p.read_bytes()
            for p in sorted(directory.rglob("*")) if p.is_file()}


@pytest.mark.parametrize("parallel", [False, 2])
def test_a_streamed_run_is_the_same_dataset_as_an_in_memory_one(tmp_path, parallel):
    # the invariant the whole feature rests on: --stream is a memory choice, not a modelling one, so
    # the same seed must leave the same bytes on disk whichever way the run was assembled
    g = _stream_fixture()
    kw = dict(model=hky85(2.0), length=120, seed=7, parallel=parallel)
    wanted = ("alignments", "phylograms", "species_phylogram", "summary")

    simulate_sequences(g, **kw).write(tmp_path / "mem", outputs=wanted)
    handle = simulate_sequences(g, **kw, stream_to=tmp_path / "strm", outputs=wanted)

    # every DATA file byte-for-byte. The summary is held out of that comparison for one field and one
    # reason: a streamed run estimates mean identity from one sampled pair per family, because the
    # in-memory way of measuring it needs every alignment at once, which is exactly what streaming
    # does not keep. Same quantity, different sample — so it is compared on its own terms below.
    mem, strm = _tree_of(tmp_path / "mem"), _tree_of(tmp_path / "strm")
    summary = "sequences_summary.json"
    assert {k: v for k, v in mem.items() if k != summary} == \
           {k: v for k, v in strm.items() if k != summary}
    a, b = (json.loads(t[summary]) for t in (mem, strm))
    assert a.pop("mean_pairwise_identity") == pytest.approx(b.pop("mean_pairwise_identity"),
                                                           abs=0.05)
    assert a == b, "the summaries may differ only in how identity was sampled"
    assert handle.n_sequences > 0 and handle.n_families > 0
    assert handle.outputs == wanted


def test_the_streamed_handle_counts_what_an_in_memory_run_counts(tmp_path):
    # a family with no surviving copy writes no alignment, so it must not be counted — otherwise a
    # streamed run reports more families than the same run in memory
    g = _stream_fixture()
    kw = dict(model=jc69(), length=60, seed=3)
    mem = simulate_sequences(g, **kw)
    handle = simulate_sequences(g, **kw, stream_to=tmp_path / "s")
    assert handle.n_families == sum(1 for a in mem.alignments.values() if a)
    assert handle.n_sequences == sum(len(a) for a in mem.alignments.values())


def test_the_streamed_handle_estimates_the_identity_the_cli_reports(tmp_path):
    # the saturation warning is the most useful thing the CLI says about a sequence run, and the
    # in-memory way of measuring it needs every alignment at once. The sink samples a pair per
    # family instead: a different sample of the same quantity, so it must land in the same region.
    g = _stream_fixture()
    kw = dict(model=hky85(2.0), length=400, seed=11)
    from zombi2.sequences import mean_pairwise_identity

    mem = simulate_sequences(g, **kw)
    handle = simulate_sequences(g, **kw, stream_to=tmp_path / "s")
    # against the CLI's own measure, which is what it stands in for. Not the first two sequences of
    # each family: those are adjacent gene ids, so usually recent duplicates, and comparing them
    # measures something noticeably more similar than a random pair.
    assert handle.identity == pytest.approx(mean_pairwise_identity(mem.alignments), abs=0.05)


def test_write_narrows_a_streamed_run(tmp_path):
    g = _stream_fixture()
    handle = simulate_sequences(g, model=jc69(), length=60, seed=1, stream_to=tmp_path / "s",
                                outputs=("alignments",))
    assert (tmp_path / "s" / "alignments").is_dir()
    assert not (tmp_path / "s" / "phylograms").exists()
    assert not (tmp_path / "s" / "clock_species_tree_complete.nwk").exists()
    assert handle.outputs == ("alignments",)


def test_streaming_is_refused_where_it_cannot_work(tmp_path):
    g = _stream_fixture()
    with pytest.raises(ValueError, match="outputs applies to a streamed run"):
        simulate_sequences(g, model=jc69(), length=60, seed=1, outputs=("alignments",))
    with pytest.raises(ValueError, match="unknown stream outputs"):
        simulate_sequences(g, model=jc69(), length=60, seed=1, stream_to=tmp_path / "s",
                           outputs=("wibble",))


def test_the_cli_streams_and_reports_the_same_run(tmp_path, capsys):
    from zombi2.cli.main import main

    def run(where, *extra):
        assert main(["species", str(where), "--birth", "1", "--death", "0.3", "--n-extant", "12",
                     "--seed", "1", "--quiet"]) == 0
        assert main(["genomes", str(where), "--duplication", "0.2", "--transfer", "0.1",
                     "--initial-families", "15", "--seed", "1", "--quiet"]) == 0
        capsys.readouterr()
        assert main(["sequences", str(where), "--model", "hky85", "--length", "90", "--seed", "5",
                     "--quiet", *extra]) == 0
        return capsys.readouterr().out

    a, b = tmp_path / "a", tmp_path / "b"
    line_mem, line_str = run(a), run(b, "--stream")
    # the same counts and the same model, reported the same way. Identity is sampled differently and
    # the wall clock is the wall clock, so neither is part of the claim.
    def strip(s):
        return re.sub(r"(mean identity [\d.]+%|in [\d.e-]+ s)", "", s.split("(", 1)[1])
    assert strip(line_mem) == strip(line_str)
    # and the files themselves are the same run, .log aside (it records the flags and the timestamp)
    # .log records the flags and a timestamp; the summary's identity is sampled differently when
    # streaming (see the byte-identity test above), so neither is part of the data claim
    data = lambda d: {k: v for k, v in _tree_of(d).items()
                      if not k.endswith((".log", "_summary.json"))}
    assert data(a / "sequences") == data(b / "sequences")


# --- site profiles: one set of equilibrium frequencies per position ----------------------------

import numpy as _np  # noqa: E402


def _profile_run(profiles=None, *, length=40, seed=2, families=2):
    from zombi2.genomes import simulate_genomes_family
    from zombi2.species import simulate_species_tree
    sp = simulate_species_tree(birth=1.0, n_extant=20, seed=1)
    g = simulate_genomes_family(sp, initial_families=families, duplication=0.1, loss=0.1, seed=1)
    return g, simulate_sequences(g, model=lg(), length=length, divergence=0.6, seed=seed,
                                 profiles=profiles)


def test_a_profile_puts_the_residue_it_prefers_at_that_site():
    """The point of a profile: position i's composition follows row i, not the model's own
    frequencies. A site whose profile is nearly a point mass comes out nearly invariant."""
    L = 40
    pref = _np.full((L, 20), 1e-3)
    pref[:, 5] = 1.0                       # every site should be the sixth residue of the alphabet
    g, res = _profile_run({0: pref}, length=L)
    want = lg().alphabet[5]
    seqs = list(res.alignments[0].values())
    matching = sum(1 for s in seqs for c in s if c == want)
    assert matching / (len(seqs) * L) > 0.95


def test_a_flat_profile_is_the_model_it_was_built_from():
    """Every row equal to the model's own frequencies must give back the model — the property that
    says a profile *adds* something rather than replacing the run's model with a different one.
    Statistical, not byte-for-byte: L single-site models walk the random stream differently from one
    L-site model."""
    import collections
    base = lg()
    flat = _np.tile(base.stationary, (200, 1))
    counts = []
    for profiles in (None, {0: flat}):
        c: collections.Counter = collections.Counter()
        for seed in range(6):
            _, res = _profile_run(profiles, length=200, seed=100 + seed, families=1)
            for aln in res.alignments.values():
                for seq in aln.values():
                    c.update(seq)
        total = sum(c.values())
        counts.append(_np.array([c[a] / total for a in base.alphabet]))
    # the two agree with each other about as closely as each agrees with the model it came from
    assert _np.abs(counts[0] - counts[1]).max() < 0.03


def test_a_profile_on_one_family_leaves_the_others_alone():
    sharp = _np.random.default_rng(0).dirichlet(_np.full(20, 0.05), size=40) + 1e-4
    _, plain = _profile_run()
    _, one = _profile_run({0: sharp})
    assert plain.alignments[1] == one.alignments[1]
    assert plain.alignments[0] != one.alignments[0]


def test_a_profile_refuses_what_it_cannot_mean():
    from zombi2.sequences import _resolve_profiles
    good = _np.full((5, 20), 0.05)
    assert len(_resolve_profiles({1: good}, lg(), 5)[1]) == 5
    for bad, fragment in (
            (_np.full((5, 19), 0.05), "must be (L, 20)"),      # wrong width for the alphabet
            (_np.full((4, 20), 0.05), "but length=5"),          # disagrees with the run's length
            (_np.zeros((5, 20)), "no state it could be in"),    # a site with nowhere to be
            (_np.eye(5, 20), "exactly zero"),                   # a residue declared impossible
    ):
        with pytest.raises(ValueError, match=fragment.replace("(", r"\(").replace(")", r"\)")):
            _resolve_profiles({1: bad}, lg(), 5)
    with pytest.raises(TypeError, match="must be a dict"):
        _resolve_profiles([good], lg(), 5)


def test_profiles_are_refused_alongside_what_would_contradict_them():
    g, _ = _profile_run()
    flat = _np.tile(lg().stationary, (40, 1))
    with pytest.raises(ValueError, match="both decide which model each site"):
        simulate_sequences(g, partitions=((lg(), 40),), profiles={0: flat}, divergence=0.5, seed=1)
    with pytest.raises(ValueError, match="not implemented for the parallel engine"):
        simulate_sequences(g, model=lg(), length=40, profiles={0: flat}, parallel=2,
                           divergence=0.5, seed=1)


def test_a_profile_composes_with_gamma():
    """A profile says which residues belong at a site; `+G` says how fast sites change. Independent,
    so both apply — and the Gamma still slows the alignment down with a profile in place."""
    L = 60
    sharp = _np.random.default_rng(1).dirichlet(_np.full(20, 0.3), size=L) + 1e-4
    from zombi2.genomes import simulate_genomes_family
    from zombi2.species import simulate_species_tree
    sp = simulate_species_tree(birth=1.0, n_extant=20, seed=1)
    g = simulate_genomes_family(sp, initial_families=1, seed=1)

    def conserved(model):
        res = simulate_sequences(g, model=model, length=L, divergence=0.6, seed=2,
                                 profiles={0: sharp})
        seqs = list(res.alignments[0].values())
        return sum(1 for i in range(L) if len({s[i] for s in seqs}) == 1)

    assert conserved(lg().across_sites(gamma_shape=0.4)) > conserved(lg())


def test_a_profile_on_a_nucleotide_block_must_match_the_length_the_genome_fixed():
    """`length` is rejected on a nucleotide run, so the row-count check in `_resolve_profiles` has
    nothing to compare against and the block's own length is what must agree. Without this a short
    profile silently shortened the sequence and the alignment stopped matching the coordinates the
    genome run wrote."""
    from zombi2.genomes import simulate_genomes_nucleotide
    from zombi2.sequences.substitution_models import jc69
    from zombi2.species import simulate_species_tree
    sp = simulate_species_tree(birth=1.0, death=0.1, n_extant=5, seed=1)
    g = simulate_genomes_nucleotide(sp, root_length=900, genes=3, gene_length=100, seed=2)
    block = sorted(g.block_trees)[0]
    _, start, end = g.root_blocks[block]
    ok = simulate_sequences(g, model=jc69(), divergence=0.2, seed=3,
                            profiles={block: _np.full((end - start, 4), 0.25)})
    assert {len(s) for s in ok.alignments[block].values()} == {end - start}
    with pytest.raises(ValueError, match="rows but that block is"):
        simulate_sequences(g, model=jc69(), divergence=0.2, seed=3,
                           profiles={block: _np.full((7, 4), 0.25)})
