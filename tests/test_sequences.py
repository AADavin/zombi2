"""Sequences level, slice 1: a nucleotide substitution model evolving along the gene trees under a
strict clock, endpoint P-matrix sampling → SequencesResult(.alignments, .ancestral)."""

from __future__ import annotations

import re

import numpy as np
import pytest

from zombi2 import species
from zombi2.genomes import GenomesResult, simulate_genomes_unordered
from zombi2.genomes.gene_trees import GeneNode, GeneTree
from zombi2.rates import modifiers as mod
from zombi2.sequences import SequencesResult, simulate_sequences
from zombi2.sequences.substitution_models import gtr, hky85, jc69, k80


# --- hand-built gene trees: origination → speciation → two extant tips -----------------------------

def _run(gene_trees, *, t_split: float = 1.0, t_now: float = 2.0) -> GenomesResult:
    """The **genome run** the sequence level requires, wrapped around hand-built gene trees: a
    three-lineage species tree (root 0 splits at ``t_split`` into extant tips 1 and 2 at ``t_now``)
    carrying the given families. The gene trees are attached directly — these families are written by
    hand, so there is no event log for the run to derive them from."""
    tree = species.Tree({0: species.Node(0, None, 0.0, t_split, (1, 2), "speciation"),
                         1: species.Node(1, 0, t_split, t_now, None, "extant"),
                         2: species.Node(2, 0, t_split, t_now, None, "extant")}, 0)
    run = GenomesResult(complete_tree=tree, genomes={}, events=[], seed=None)
    run.gene_trees = dict(gene_trees)      # a cached_property: the instance dict wins
    return run


def _pair_run(t_spec: float, t_tip: float) -> GenomesResult:
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
    assert r.alignments[0]["g1"] == r.ancestral[0]["g0"]


def test_every_sequence_has_the_requested_length_over_the_alphabet():
    r = simulate_sequences(_pair_run(1.0, 2.0), model=gtr(), length=123, seed=1)
    for seq in _seqs(r):
        assert len(seq) == 123
        assert set(seq) <= set("ACGT")


def test_alignment_labels_are_exactly_the_extant_tips():
    r = simulate_sequences(_pair_run(1.0, 2.0), model=k80(2.0), length=10, seed=1)
    assert set(r.alignments[0]) == {"g1", "g2"}             # the two extant daughter genes


def test_jc69_holds_uniform_base_composition():
    # root drawn from the uniform stationary; JC69 keeps it uniform, so tips stay ≈ 25% each
    r = simulate_sequences(_pair_run(1.0, 5.0), model=jc69(), length=20000, seed=42)
    seq = r.alignments[0]["g1"]
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
        return pdist(r.ancestral[0]["g0"], r.alignments[0]["g1"])   # root gene g0 → daughter tip g1

    theory = lambda d: 0.75 * (1 - np.exp(-4 * d / 3))          # noqa: E731
    assert abs(root_tip_pdistance(t_tip=1.0, substitution=1.0) - theory(1.0)) < 0.01
    # substitution = 0.5 halves the effective distance: Δt=2 behaves like d=1.0
    assert abs(root_tip_pdistance(t_tip=2.0, substitution=0.5) - theory(1.0)) < 0.01


def test_hky85_transition_bias_makes_diverged_tips_still_reflect_frequencies():
    # a strongly skewed base composition is reproduced at the tips (endpoint stays near stationary)
    freqs = (0.4, 0.1, 0.1, 0.4)
    r = simulate_sequences(_pair_run(1.0, 4.0), model=hky85(4.0, freqs), length=20000, seed=5)
    seq = r.alignments[0]["g1"]
    comp = [seq.count(b) / len(seq) for b in "ACGT"]
    assert comp[0] > comp[1] and comp[3] > comp[2]   # A,T (0.4) exceed C,G (0.1)


# --- a family with no surviving copy ---------------------------------------------------------------

def test_family_with_no_extant_copy_has_empty_alignment_but_full_ancestral():
    # the root gene (id 0) speciates, but both daughters are lost → no extant tip, yet the root is a
    # real internal node with a reconstructed sequence
    root = GeneNode("speciation", 0, 1.0, 0)
    root.children = [GeneNode("loss", 0, 2.0, 1), GeneNode("loss", 0, 2.0, 2)]
    r = simulate_sequences(_run({0: GeneTree(0, root, 0.0)}), model=jc69(), length=10, seed=1)
    assert r.alignments[0] == {}
    assert set(r.ancestral[0]) == {"g0"}               # the root gene still gets a sequence


# --- integration: species → genomes → sequences ----------------------------------------------------

def test_a_real_genome_run_is_covered_node_for_node():
    sp = species.simulate_species_tree(birth=1.0, death=0.2, n_extant=8, seed=1)
    g = simulate_genomes_unordered(sp, duplication=0.2, loss=0.2, transfer=0.1,
                                   initial_families=6, seed=2)
    r = simulate_sequences(g, model=hky85(kappa=2.5), length=300, seed=3)   # the genome run itself
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
    ds = [_pdist(r.ancestral[f]["g0"], r.alignments[f]["g1"]) for f in range(20)]
    mean = sum(ds) / len(ds)
    std = (sum((d - mean) ** 2 for d in ds) / len(ds)) ** 0.5
    assert std < 0.02      # shared clock ⇒ ~0.006 sampling noise; a per-family clock would be far larger


def test_bylineage_rejects_other_and_multiple_modifiers():
    run = _pair_run(1.0, 2.0)
    with pytest.raises(ValueError):                 # FromParent clock — a later slice
        simulate_sequences(run, model=jc69(), length=10, substitution=1.0 * mod.FromParent(spread=0.3))
    with pytest.raises(ValueError):                 # two ByLineage — only a single clock is wired
        simulate_sequences(run, model=jc69(), length=10,
                           substitution=1.0 * mod.ByLineage(spread=0.3) * mod.ByLineage(spread=0.2))
    with pytest.raises(ValueError):                 # ByLineage × OnTime — mixed modifiers
        simulate_sequences(run, model=jc69(), length=10,
                           substitution=1.0 * mod.ByLineage(spread=0.3) * mod.OnTime({0: 1.0}))


# --- phylograms: the gene / species trees in substitutions/site -------------------------------------

def _leaves(nwk: str) -> set[str]:
    return set(re.findall(r"(?<![)\w])g\d+", nwk))     # g<copy> in leaf position (not after a ')')


def _total_bl(nwk: str) -> float:
    return sum(float(x) for x in re.findall(r":([0-9.eE+-]+)", nwk))


def _small_run(clock=1.0):
    sp = species.simulate_species_tree(birth=1.0, death=0.2, n_extant=8, seed=1)
    g = simulate_genomes_unordered(sp, duplication=0.3, loss=0.3, transfer=0.1, initial_families=8, seed=2)
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
        root_seq = r.ancestral[fam].get(f"g{gt.complete.copy}")
        if root_seq is not None and stem > 0:
            # not asserted equal or unequal site-by-site (a short stem may fix nothing), but the
            # branch the phylogram reports must be exactly the stem under a rate-1 strict clock
            written = float(r.phylograms[fam]["complete"].rsplit(":", 1)[1].rstrip(";"))
            assert written == pytest.approx(stem, rel=1e-5)


def test_phylogram_internal_nodes_pair_with_ancestral():
    # every internal node in the complete phylogram is labelled by its gene id, matching an ancestral key
    g, r = _small_run()
    for fam in g.gene_trees:
        internal = set(re.findall(r"\)(g\d+)", r.phylograms[fam]["complete"]))
        assert internal == set(r.ancestral[fam])


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
    with pytest.raises(TypeError, match="simulate_genomes_unordered"):
        simulate_sequences(g.gene_trees, model=jc69(), length=10, seed=3)


def test_write_emits_phylogram_newick(tmp_path):
    _, r = _small_run()
    r.write(tmp_path, outputs=("phylograms", "species_phylogram"))
    fam0 = tmp_path / "sequences_phylogram_fam0_complete.nwk"
    assert fam0.exists() and fam0.read_text().rstrip().endswith(";")
    assert (tmp_path / "sequences_species_phylogram_complete.nwk").exists()


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
    aln = tmp_path / "sequences_alignment_fam0.fasta"
    assert aln.exists() and ">g1" in aln.read_text()
    r.write(tmp_path, outputs=("ancestral",))
    assert (tmp_path / "sequences_ancestral_fam0.fasta").exists()
    with pytest.raises(ValueError):
        r.write(tmp_path, outputs=("bogus",))
