"""Sequences level, slice 1: a nucleotide substitution model evolving along the gene trees under a
strict clock, endpoint P-matrix sampling → SequencesResult(.alignments, .ancestral)."""

from __future__ import annotations

import numpy as np
import pytest

from zombi2 import species
from zombi2.genomes import simulate_genomes_unordered
from zombi2.genomes.gene_trees import GeneNode, GeneTree
from zombi2.rates import modifiers as mod
from zombi2.sequences import SequencesResult, simulate_sequences
from zombi2.sequences.substitution_models import gtr, hky85, jc69, k80


# --- hand-built gene trees: origination → speciation → two extant tips -----------------------------

def _pair_tree(t_spec: float, t_tip: float) -> GeneTree:
    """A minimal family: a founding ``origination`` at t=0, a ``speciation`` at ``t_spec``, then two
    ``extant`` tips (copy 0, species 1 and 2) at ``t_tip``. ``t_spec == 0`` gives a zero-length
    origination→speciation branch (the child then copies the root)."""
    root = GeneNode("origination", 0, 0.0, 0)
    spec = GeneNode("speciation", 0, t_spec, 0)
    spec.children = [GeneNode("extant", 1, t_tip, 0), GeneNode("extant", 2, t_tip, 0)]
    root.children = [spec]
    return GeneTree(0, root)


def _single_branch_family(species: int = 0, t_tip: float = 1.0) -> GeneTree:
    """One family = a founding ``origination`` on ``species`` at t=0 → one ``extant`` tip at ``t_tip``,
    on the same species lineage. Used to test that the lineage clock is shared across families."""
    root = GeneNode("origination", species, 0.0, 0)
    root.children = [GeneNode("extant", species, t_tip, 0)]
    return GeneTree(0, root)


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


# --- the engine: determinism, the strict clock, structure ------------------------------------------

def test_determinism_same_seed_identical_different_seed_differs():
    gts = {0: _pair_tree(1.0, 2.0)}
    a = simulate_sequences(gts, model=jc69(), length=200, seed=7)
    b = simulate_sequences(gts, model=jc69(), length=200, seed=7)
    assert a.alignments == b.alignments and a.ancestral == b.ancestral
    c = simulate_sequences(gts, model=jc69(), length=200, seed=8)
    assert a.alignments != c.alignments


def test_zero_rate_leaves_every_node_equal_to_the_root():
    # substitution = 0 → every branch length is 0 → no site ever changes
    r = simulate_sequences({0: _pair_tree(1.0, 2.0)}, model=hky85(kappa=3.0), length=150,
                           substitution=0.0, seed=1)
    assert len(set(_seqs(r))) == 1


def test_zero_length_branch_copies_its_parent():
    # t_spec = 0 → the origination→speciation branch has zero length, so they share a sequence
    r = simulate_sequences({0: _pair_tree(0.0, 1.0)}, model=jc69(), length=120, seed=3)
    anc = r.ancestral[0]
    orig = next(v for label, v in anc.items() if label.startswith("origination"))
    spec = next(v for label, v in anc.items() if label.startswith("speciation"))
    assert orig == spec


def test_every_sequence_has_the_requested_length_over_the_alphabet():
    r = simulate_sequences({0: _pair_tree(1.0, 2.0)}, model=gtr(), length=123, seed=1)
    for seq in _seqs(r):
        assert len(seq) == 123
        assert set(seq) <= set("ACGT")


def test_alignment_labels_are_exactly_the_extant_tips():
    r = simulate_sequences({0: _pair_tree(1.0, 2.0)}, model=k80(2.0), length=10, seed=1)
    assert set(r.alignments[0]) == {"g0_n1", "g0_n2"}


def test_jc69_holds_uniform_base_composition():
    # root drawn from the uniform stationary; JC69 keeps it uniform, so tips stay ≈ 25% each
    r = simulate_sequences({0: _pair_tree(1.0, 5.0)}, model=jc69(), length=20000, seed=42)
    seq = r.alignments[0]["g0_n1"]
    for base in "ACGT":
        assert abs(seq.count(base) / len(seq) - 0.25) < 0.03


def test_jc69_pdistance_matches_theory_and_rate_scales_it():
    # the endpoint distribution must match the JC69 process: p(d) = 3/4 (1 - exp(-4d/3)), with the
    # branch length d = substitution · Δt. Guards the model normalisation and the time→subs/site
    # conversion — a structural test would miss a rescaled or mis-normalised rate matrix.
    def pdist(a, b):
        return np.mean(np.frombuffer(a.encode(), np.uint8) != np.frombuffer(b.encode(), np.uint8))

    def root_tip_pdistance(*, t_tip, substitution):
        r = simulate_sequences({0: _pair_tree(0.0, t_tip)}, model=jc69(), length=40000,
                               substitution=substitution, seed=1)
        root = next(v for label, v in r.ancestral[0].items() if label.startswith("origination"))
        return pdist(root, r.alignments[0]["g0_n1"])

    theory = lambda d: 0.75 * (1 - np.exp(-4 * d / 3))          # noqa: E731
    assert abs(root_tip_pdistance(t_tip=1.0, substitution=1.0) - theory(1.0)) < 0.01
    # substitution = 0.5 halves the effective distance: Δt=2 behaves like d=1.0
    assert abs(root_tip_pdistance(t_tip=2.0, substitution=0.5) - theory(1.0)) < 0.01


def test_hky85_transition_bias_makes_diverged_tips_still_reflect_frequencies():
    # a strongly skewed base composition is reproduced at the tips (endpoint stays near stationary)
    freqs = (0.4, 0.1, 0.1, 0.4)
    r = simulate_sequences({0: _pair_tree(1.0, 4.0)}, model=hky85(4.0, freqs), length=20000, seed=5)
    seq = r.alignments[0]["g0_n1"]
    comp = [seq.count(b) / len(seq) for b in "ACGT"]
    assert comp[0] > comp[1] and comp[3] > comp[2]   # A,T (0.4) exceed C,G (0.1)


# --- a family with no surviving copy ---------------------------------------------------------------

def test_family_with_no_extant_copy_has_empty_alignment_but_full_ancestral():
    root = GeneNode("origination", 0, 0.0, 0)
    root.children = [GeneNode("loss", 0, 1.0, 0)]      # the only copy is lost mid-branch
    r = simulate_sequences({0: GeneTree(0, root)}, model=jc69(), length=10, seed=1)
    assert r.alignments[0] == {}
    assert len(r.ancestral[0]) == 1                    # the origination node still gets a sequence


# --- integration: species → genomes → sequences ----------------------------------------------------

def test_accepts_a_genomes_result_and_covers_every_node():
    sp = species.simulate_species_tree(birth=1.0, death=0.2, n_extant=8, seed=1)
    g = simulate_genomes_unordered(sp, duplication=0.2, loss=0.2, transfer=0.1,
                                   initial_families=6, seed=2)
    r = simulate_sequences(g, model=hky85(kappa=2.5), length=300, seed=3)   # GenomesResult input
    assert isinstance(r, SequencesResult)
    assert set(r.alignments) == set(g.gene_trees) == set(r.ancestral)       # one entry per family
    for fam, gt in g.gene_trees.items():
        nodes = list(_iter_nodes(gt.complete))
        n_internal = sum(1 for n in nodes if not n.is_leaf)
        n_extant = sum(1 for n in nodes if n.is_leaf and n.kind == "extant")
        assert len(r.ancestral[fam]) == n_internal
        assert len(r.alignments[fam]) == n_extant
        for seq in list(r.alignments[fam].values()) + list(r.ancestral[fam].values()):
            assert len(seq) == 300


def test_integration_is_deterministic_given_the_seed():
    sp = species.simulate_species_tree(birth=1.0, death=0.3, n_extant=6, seed=11)
    g = simulate_genomes_unordered(sp, duplication=0.3, loss=0.2, initial_families=5, seed=12)
    a = simulate_sequences(g, model=gtr(), length=100, seed=13)
    b = simulate_sequences(g, model=gtr(), length=100, seed=13)
    assert a.alignments == b.alignments and a.ancestral == b.ancestral


# --- the lineage clock (ByLineage): the uncorrelated / relaxed clock --------------------------------

def _pdist(a: str, b: str) -> float:
    return sum(x != y for x, y in zip(a, b)) / len(a)


def test_bylineage_zero_spread_is_bit_identical_to_the_strict_clock():
    # spread=0 draws 1.0 without touching the rng, so the run matches the strict clock exactly
    gts = {0: _pair_tree(1.0, 2.0)}
    strict = simulate_sequences(gts, model=jc69(), length=300, seed=5)
    clocked = simulate_sequences(gts, model=jc69(), length=300,
                                 substitution=1.0 * mod.ByLineage(spread=0.0), seed=5)
    assert clocked.alignments == strict.alignments and clocked.ancestral == strict.ancestral


def test_bylineage_perturbs_the_output_and_stays_valid():
    gts = {0: _pair_tree(1.0, 2.0)}
    strict = simulate_sequences(gts, model=jc69(), length=300, seed=5)
    clocked = simulate_sequences(gts, model=jc69(), length=300,
                                 substitution=1.0 * mod.ByLineage(spread=0.5), seed=5)
    assert clocked.alignments != strict.alignments          # the clock rescales branch lengths
    for seq in _seqs(clocked):
        assert len(seq) == 300 and set(seq) <= set("ACGT")


def test_bylineage_is_deterministic():
    gts = {0: _pair_tree(1.0, 2.0)}
    spec = 1.0 * mod.ByLineage(spread=0.4)
    a = simulate_sequences(gts, model=hky85(2.0), length=200, substitution=spec, seed=9)
    b = simulate_sequences(gts, model=hky85(2.0), length=200, substitution=spec, seed=9)
    assert a.alignments == b.alignments and a.ancestral == b.ancestral


def test_bylineage_clock_is_shared_across_families_on_a_lineage():
    # 20 identical single-branch families, ALL on species lineage 0 → all feel the SAME clock[0].
    # A per-family clock would scatter their root→tip divergences by ~spread; a shared clock leaves
    # only sampling noise, so the across-family spread collapses.
    gts = {f: _single_branch_family(species=0, t_tip=1.0) for f in range(20)}
    r = simulate_sequences(gts, model=jc69(), length=5000,
                           substitution=1.0 * mod.ByLineage(spread=0.8), seed=4)
    ds = [_pdist(next(v for lab, v in r.ancestral[f].items() if lab.startswith("origination")),
                 r.alignments[f]["g0_n0"]) for f in range(20)]
    mean = sum(ds) / len(ds)
    std = (sum((d - mean) ** 2 for d in ds) / len(ds)) ** 0.5
    assert std < 0.02      # shared clock ⇒ ~0.006 sampling noise; a per-family clock would be far larger


def test_bylineage_rejects_other_and_multiple_modifiers():
    gts = {0: _pair_tree(1.0, 2.0)}
    with pytest.raises(ValueError):                 # Inherited clock — a later slice
        simulate_sequences(gts, model=jc69(), length=10, substitution=1.0 * mod.Inherited(spread=0.3))
    with pytest.raises(ValueError):                 # two ByLineage — only a single clock is wired
        simulate_sequences(gts, model=jc69(), length=10,
                           substitution=1.0 * mod.ByLineage(spread=0.3) * mod.ByLineage(spread=0.2))
    with pytest.raises(ValueError):                 # ByLineage × Time — mixed modifiers
        simulate_sequences(gts, model=jc69(), length=10,
                           substitution=1.0 * mod.ByLineage(spread=0.3) * mod.Time({0: 1.0}))


# --- validation ------------------------------------------------------------------------------------

def test_rejects_bad_arguments_and_unwired_rate_specs():
    from zombi2.rates.scope import PerLineage
    gts = {0: _pair_tree(1.0, 2.0)}
    with pytest.raises(TypeError):
        simulate_sequences(gts, model="jc69", length=10)                  # not a SubstitutionModel
    with pytest.raises(ValueError):
        simulate_sequences(gts, model=jc69(), length=0)                   # non-positive length
    with pytest.raises(ValueError):
        simulate_sequences(gts, model=jc69(), length=10, substitution=1.0 * mod.Time({0: 1.0}))
    with pytest.raises(ValueError):
        simulate_sequences(gts, model=jc69(), length=10, substitution=PerLineage(1.0))


# --- writing ---------------------------------------------------------------------------------------

def test_write_emits_fasta_per_family(tmp_path):
    r = simulate_sequences({0: _pair_tree(1.0, 2.0)}, model=jc69(), length=20, seed=1)
    r.write(tmp_path)
    aln = tmp_path / "sequences_alignment_fam0.fasta"
    assert aln.exists() and ">g0_n1" in aln.read_text()
    r.write(tmp_path, outputs=("ancestral",))
    assert (tmp_path / "sequences_ancestral_fam0.fasta").exists()
    with pytest.raises(ValueError):
        r.write(tmp_path, outputs=("bogus",))
