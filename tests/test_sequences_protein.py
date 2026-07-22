"""Sequences level: the **protein** menu — the empirical 20-state models (poisson · jtt · dayhoff ·
wag · lg), their published matrices, and a protein sequence evolving along the gene trees.

The matrices are data, so the tests that matter are the ones a reader cannot do by eye: that the
20 rows of ``Q`` sit in the order :data:`AMINO_ACIDS` claims (a permuted matrix still *looks*
plausible), that the published numbers landed in the right cells, and that a simulated alignment
reproduces the model it was drawn from — composition, expected substitutions, and which residue
replaces which.
"""

from __future__ import annotations

import numpy as np
import pytest

from zombi2 import species
from zombi2.genomes import GenomesResult, simulate_genomes_unordered
from zombi2.genomes.gene_trees import GeneNode, GeneTree
from zombi2.rates import modifiers as mod
from zombi2.sequences import _aa_matrices, simulate_sequences
from zombi2.sequences.substitution_models import (
    AMINO_ACIDS, _lower_triangle, dayhoff, jc69, jtt, lg, poisson, wag,
)

PROTEIN_MODELS = (poisson, jtt, dayhoff, wag, lg)
EMPIRICAL = (("LG", lg, _aa_matrices._LG_EXCH, _aa_matrices._LG_PI),
             ("WAG", wag, _aa_matrices._WAG_EXCH, _aa_matrices._WAG_PI),
             ("JTT", jtt, _aa_matrices._JTT_EXCH, _aa_matrices._JTT_PI),
             ("Dayhoff", dayhoff, _aa_matrices._DAYHOFF_EXCH, _aa_matrices._DAYHOFF_PI))
IDX = {a: i for i, a in enumerate(AMINO_ACIDS)}


def _cell(tri, a: str, b: str) -> float:
    """The published exchangeability of residues ``a`` and ``b``, read straight out of the flat lower
    triangle (row ``i`` starts at ``i(i-1)/2``) — the file's own layout, not the assembled matrix."""
    i, j = max(IDX[a], IDX[b]), min(IDX[a], IDX[b])
    return tri[i * (i - 1) // 2 + j]


def _genome_run(gene_trees: dict[int, GeneTree], *, t_split: float = 1.0, t_now: float = 2.0):
    """The genome run the sequence level requires, around hand-built gene trees: a three-lineage
    species tree (root 0 splits at ``t_split`` into extant tips 1 and 2 at ``t_now``) carrying the
    given families. The gene trees are attached directly — these families are written by hand, so
    there is no event log to derive them from."""
    tree = species.Tree({0: species.Node(0, None, 0.0, t_split, (1, 2), "speciation"),
                         1: species.Node(1, 0, t_split, t_now, None, "extant"),
                         2: species.Node(2, 0, t_split, t_now, None, "extant")}, 0)
    run = GenomesResult(complete_tree=tree, genomes={}, events=[], seed=None)
    run.gene_trees = dict(gene_trees)      # a cached_property: the instance dict wins
    return run


def _one_branch(t_tip: float, families: int = 1):
    """``families`` copies of a minimal family: a root gene (id 0) at t=0 that speciates immediately
    into two tips (ids 1, 2) at ``t_tip`` — so each root→tip branch has length exactly ``t_tip``."""
    def tree(f):
        root = GeneNode("speciation", 0, 0.0, 0)
        root.children = [GeneNode("extant", 1, t_tip, 1), GeneNode("extant", 2, t_tip, 2)]
        return GeneTree(f, root, 0.0)
    return _genome_run({f: tree(f) for f in range(families)}, t_split=0.0, t_now=t_tip)


# --- the alphabet and the matrices: is every number where it claims to be? -------------------------

def test_amino_acid_alphabet_is_the_paml_column_order():
    # the order the published matrices are written in — get this wrong and every model is a silent
    # permutation of itself
    assert AMINO_ACIDS == "".join("A R N D C Q E G H I L K M F P S T W Y V".split())
    assert len(AMINO_ACIDS) == len(set(AMINO_ACIDS)) == 20


def test_lower_triangle_expands_into_the_cells_paml_names():
    for _, _, tri, _ in EMPIRICAL:
        S = _lower_triangle(tri, 20)
        assert S.shape == (20, 20)
        assert np.allclose(S, S.T) and (np.diag(S) == 0).all()
        assert S[IDX["R"], IDX["A"]] == tri[0]        # the first stored value is row R, column A
        assert S[IDX["V"], IDX["Y"]] == tri[-1]       # the last is row V, column Y
        assert (S[np.triu_indices(20, 1)] > 0).sum() > 150   # a mostly-filled triangle


def test_published_spot_values_sit_in_the_named_cells():
    # values read off the PAML .dat files (lg.dat, wag.dat, jones.dat, dayhoff.dat)
    published = {"LG": {("R", "A"): 0.425093, ("I", "V"): 10.649107, ("V", "Y"): 0.249313},
                 "WAG": {("R", "A"): 0.551571, ("I", "V"): 7.8213, ("V", "Y"): 0.31473},
                 "JTT": {("R", "A"): 58, ("I", "V"): 961, ("V", "Y"): 16},
                 "Dayhoff": {("R", "A"): 27, ("D", "E"): 1153, ("V", "Y"): 28}}
    for name, _, tri, _ in EMPIRICAL:
        for (a, b), value in published[name].items():
            assert _cell(tri, a, b) == value, f"{name} S({a},{b})"


def test_the_matrices_are_in_register_with_the_alphabet():
    """The ordering check that survives a permutation: chemistry. Every empirical model's most
    exchangeable pairs are the conservative substitutions (I↔V, D↔E, F↔Y, R↔K, S↔T), and its rarest
    residue is tryptophan — both true of the published matrices, neither true of a shuffled one."""
    conservative = {"IV", "VI", "DE", "ED", "FY", "YF", "RK", "KR", "ST", "TS", "LM", "ML",
                    "ND", "DN", "QE", "EQ", "QH", "HQ", "HY", "YH"}
    for name, build, tri, pi in EMPIRICAL:
        S = _lower_triangle(tri, 20)
        iu = np.triu_indices(20, 1)
        top = [AMINO_ACIDS[iu[0][k]] + AMINO_ACIDS[iu[1][k]]
               for k in np.argsort(S[iu])[::-1][:5]]
        assert set(top) <= conservative, f"{name} top exchangeabilities are {top}"
        assert "IV" in top, f"{name}: isoleucine↔valine should be among the most exchangeable"
        pi = np.asarray(pi)
        assert AMINO_ACIDS[pi.argmin()] == "W", f"{name}: the rarest residue should be tryptophan"
        assert AMINO_ACIDS[pi.argmax()] in "ALG", f"{name}: the commonest should be A / L / G"
        assert np.allclose(build().stationary, pi / pi.sum())   # π reaches the model unpermuted


def test_every_protein_model_is_a_valid_reversible_20_state_model():
    for build in PROTEIN_MODELS:
        m = build()
        assert m.k == 20 and m.alphabet == AMINO_ACIDS
        assert np.isclose(m.stationary.sum(), 1.0) and (m.stationary > 0).all()
        pi_Q = m.stationary[:, None] * m.Q
        assert np.allclose(pi_Q, pi_Q.T)                       # detailed balance: π_i Q_ij = π_j Q_ji
        P = m.p_matrix(0.42)
        assert np.allclose(P.sum(1), 1.0) and (P >= 0).all()
        assert np.allclose(m.stationary @ P, m.stationary)     # π is the stationary distribution
        assert np.allclose(m.p_matrix(0.0), np.eye(20))


def test_protein_models_are_normalised_to_one_substitution_per_unit_time():
    # -Σ π_i Q_ii == 1, exactly as the nucleotide models — branch lengths stay in substitutions/site
    for build in PROTEIN_MODELS:
        m = build()
        assert np.isclose(-(m.stationary * np.diag(m.Q)).sum(), 1.0)


def test_poisson_is_the_jc69_of_proteins():
    m = poisson()
    assert np.allclose(m.stationary, 1.0 / 20.0)
    off = m.Q[~np.eye(20, dtype=bool)]
    assert np.allclose(off, off[0])                            # every off-diagonal rate is equal
    assert np.isclose(off[0], 1.0 / 19.0)                      # ...and normalised: 19 · q = 1


# --- simulation: does an alignment reproduce the model it was drawn from? --------------------------

def test_sequences_use_all_twenty_residues_and_only_those():
    r = simulate_sequences(_one_branch(1.0), model=lg(), length=4000, seed=1)
    seqs = list(r.alignments[0].values()) + list(r.ancestral[0].values())
    assert len(seqs) == 3 and all(len(s) == 4000 for s in seqs)
    for s in seqs:
        assert set(s) == set(AMINO_ACIDS)                      # all 20 present, nothing else


def test_long_run_composition_matches_the_model_frequencies():
    # a long branch washes out the root draw, so the tip composition is the model's own π
    for m in (lg(), wag(), jtt(), dayhoff()):
        r = simulate_sequences(_one_branch(8.0), model=m, length=60000, seed=3)
        seq = r.alignments[0]["g1"]
        observed = np.array([seq.count(a) for a in AMINO_ACIDS]) / len(seq)
        assert np.abs(observed - m.stationary).max() < 0.008, m.name


def test_expected_substitutions_per_unit_branch_length_is_one():
    """The normalisation, measured on the sequences rather than on ``Q``: at branch length ``d`` the
    chance a site differs is ``1 - Σ π_i P(d)_ii`` — which is ``d`` itself as ``d → 0``, i.e. one
    expected substitution per unit branch length. Checked at a short branch and a long one."""
    m = lg()
    for d in (0.02, 1.5):
        theory = 1.0 - float((m.stationary * np.diag(m.p_matrix(d))).sum())
        r = simulate_sequences(_one_branch(d), model=m, length=60000, seed=7)
        root, tip = r.ancestral[0]["g0"], r.alignments[0]["g1"]
        observed = np.mean(np.frombuffer(root.encode(), np.uint8)
                           != np.frombuffer(tip.encode(), np.uint8))
        assert abs(observed - theory) < 0.008, d
        if d == 0.02:                      # short branch: the theory *is* the branch length
            assert abs(theory - d) < 0.02 * d
            assert abs(observed - d) < 0.002


def test_isoleucine_prefers_valine_in_the_simulated_alignment():
    """The register check, end to end: LG's largest exchangeability is I↔V, so among the sites that
    started as isoleucine and changed, valine must be the commonest destination — in the decoded
    sequences, not in the matrix. A permuted alphabet or a mis-decoded state array breaks this."""
    r = simulate_sequences(_one_branch(0.1), model=lg(), length=100000, seed=11)
    root, tip = r.ancestral[0]["g0"], r.alignments[0]["g1"]
    destinations: dict[str, int] = {}
    for a, b in zip(root, tip):
        if a == "I" and b != "I":
            destinations[b] = destinations.get(b, 0) + 1
    assert max(destinations, key=destinations.get) == "V"
    # I→V beats the runner-up I→L by the ratio of their fluxes, S·π: 10.649·0.0691 / 4.145·0.0991 ≈ 1.8
    assert destinations["V"] > 1.4 * destinations["L"]


def test_the_empirical_models_are_genuinely_different():
    # same tree, same seed, same length: LG, WAG, JTT and Dayhoff are different chemistry, so they
    # must not agree — and none of them is the nucleotide alphabet
    run = _one_branch(0.6)
    tips = {}
    for build in (lg, wag, jtt, dayhoff, poisson):
        r = simulate_sequences(run, model=build(), length=3000, seed=5)
        tips[build().name] = r.alignments[0]["g1"]
    assert len(set(tips.values())) == len(tips)              # five distinct alignments
    assert not any(set(s) <= set("ACGT") for s in tips.values())


def test_determinism_and_the_lineage_clock_carry_over_to_proteins():
    run = _one_branch(1.0)
    a = simulate_sequences(run, model=wag(), length=500, seed=2)
    b = simulate_sequences(run, model=wag(), length=500, seed=2)
    assert a.alignments == b.alignments and a.ancestral == b.ancestral
    clocked = simulate_sequences(run, model=wag(), length=500,
                                 substitution=1.0 * mod.ByLineage(spread=0.5), seed=2)
    assert clocked.alignments != a.alignments                # the clock rescales the branches
    assert set(clocked.alignments[0]["g1"]) <= set(AMINO_ACIDS)


def test_zero_rate_leaves_a_protein_alignment_untouched():
    r = simulate_sequences(_one_branch(1.0), model=lg(), length=200, substitution=0.0, seed=1)
    assert len(set(list(r.alignments[0].values()) + list(r.ancestral[0].values()))) == 1


# --- integration: species → genomes → protein sequences --------------------------------------------

def test_protein_run_over_a_real_genome_history(tmp_path):
    sp = species.simulate_species_tree(birth=1.0, death=0.2, n_extant=8, seed=1)
    g = simulate_genomes_unordered(sp, duplication=0.2, loss=0.2, transfer=0.1,
                                   initial_families=6, seed=2)
    r = simulate_sequences(g, model=lg(), length=120, seed=3)
    assert set(r.alignments) == set(g.gene_trees)
    for aln in r.alignments.values():
        for seq in aln.values():
            assert len(seq) == 120 and set(seq) <= set(AMINO_ACIDS)
    # the outputs are the same files as for a nucleotide run — only the residues differ
    r.write(tmp_path)
    fasta = next(p for p in tmp_path.iterdir() if p.name.startswith("fam"))
    assert fasta.read_text().startswith(">g")
    assert (tmp_path / "phylogram_fam0_complete.nwk").exists()


def test_a_nucleotide_and_a_protein_model_stay_in_their_own_alphabets():
    run = _one_branch(1.0)
    dna = simulate_sequences(run, model=jc69(), length=300, seed=4).alignments[0]["g1"]
    protein = simulate_sequences(run, model=lg(), length=300, seed=4).alignments[0]["g1"]
    assert set(dna) <= set("ACGT")
    assert not set(protein) <= set("ACGT")


def test_a_protein_model_takes_no_parameters():
    for build in PROTEIN_MODELS:
        with pytest.raises(TypeError):
            build(2.0)          # empirical: there is nothing to tune
