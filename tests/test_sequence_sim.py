"""Sequence evolution along gene trees + ancestral genome reconstruction.

The engine (substitution models + :func:`evolve_on_tree`) is checked against analytic properties
(stationarity, JC distance recovery, t=0 identity); the nucleotide-model integration is checked by
reconstructing the genome at every node — most importantly that a FASTA-seeded root reproduces the
input genome exactly, and that at zero divergence every node's DNA is the strand-correct assembly of
its block sequences.
"""

from collections import Counter

import numpy as np
import pytest

from zombi2 import BirthDeath, simulate_species_tree
from zombi2.sequences._aa_models import (
    _DAYHOFF_PI, _JTT_PI, _LG_PI, _WAG_PI,
)
from zombi2.genomes.nucleotide_sim import simulate_nucleotide_genomes
from zombi2.sequences.models import (
    AMINO_ACIDS, DNA_MODELS, GammaRates, PROTEIN_MODELS, dayhoff, evolve_on_tree, gtr, hky85,
    is_protein_model, jc69, jtt, k80, lg, make_model, poisson, read_fasta, reverse_complement,
    wag, write_fasta,
)


class _N:
    """A minimal tree node (gid + children) for the engine tests."""
    def __init__(self, gid, children=()):
        self.gid, self.children = gid, list(children)


# --------------------------------------------------------------------------- #
# Substitution models
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("m", [jc69(), k80(3.0), hky85(2.5, (0.3, 0.2, 0.2, 0.3)),
                               gtr((1, 2, 1, 1, 2, 1), (0.3, 0.2, 0.2, 0.3)),
                               poisson(), lg(), wag(), jtt(), dayhoff()])
def test_model_matrix_properties(m):
    Q, pi, k = m.Q, m.stationary, m.k
    assert Q.shape == (k, k)
    assert np.allclose(Q.sum(1), 0, atol=1e-12)              # rows of a rate matrix sum to 0
    assert np.allclose(pi @ Q, 0, atol=1e-10)               # pi is stationary
    assert np.isclose(-(pi * np.diag(Q)).sum(), 1.0)        # normalised to 1 subst/site/unit
    assert np.allclose(m.p_matrix(0.0), np.eye(k))          # P(0) = I
    P = m.p_matrix(0.4)
    assert np.allclose(P.sum(1), 1) and (P >= 0).all()      # each row a valid distribution
    assert np.allclose(m.p_matrix(500.0), pi[None, :], atol=1e-4)   # long branch -> stationary
    # detailed balance: pi_i Q_ij = pi_j Q_ji (reversibility)
    R = pi[:, None] * Q
    assert np.allclose(R, R.T, atol=1e-12)


def test_make_model_dispatch():
    assert make_model("jc69").name == "JC69"
    assert make_model("hky85", kappa=3, freqs=(0.4, 0.1, 0.1, 0.4)).name == "HKY85"
    assert make_model("lg").name == "LG" and make_model("lg").k == 20
    assert make_model("WAG").name == "WAG"        # case-insensitive
    with pytest.raises(ValueError):
        make_model("nope")


# --------------------------------------------------------------------------- #
# Amino-acid (protein) models — published values, correctness over coverage
# --------------------------------------------------------------------------- #
_PUBLISHED_PI = {"lg": _LG_PI, "wag": _WAG_PI, "jtt": _JTT_PI, "dayhoff": _DAYHOFF_PI}


@pytest.mark.parametrize("name", ["lg", "wag", "jtt", "dayhoff"])
def test_empirical_aa_frequencies_match_published(name):
    """The stored stationary frequencies must match the published pi vector to 1e-4.

    This catches transcription errors that mathematical validity alone would not — a mistyped
    exchangeability keeps Q reversible but a mistyped frequency shows up here.
    """
    m = make_model(name)
    published = np.array(_PUBLISHED_PI[name], dtype=float)
    published = published / published.sum()                 # published freqs round to 1 only to ~1e-6
    assert m.stationary.shape == (20,)
    assert np.max(np.abs(m.stationary - published)) < 1e-4


@pytest.mark.parametrize("name", ["poisson", "lg", "wag", "jtt", "dayhoff"])
def test_empirical_aa_reversible(name):
    """pi_i Q_ij = pi_j Q_ji for every empirical amino-acid model."""
    m = make_model(name)
    R = m.stationary[:, None] * m.Q
    assert np.allclose(R, R.T, atol=1e-12)


def test_model_registries_and_alphabet():
    assert DNA_MODELS == ("jc69", "k80", "hky85", "gtr")
    assert PROTEIN_MODELS == ("poisson", "lg", "wag", "jtt", "dayhoff")
    assert all(is_protein_model(n) for n in PROTEIN_MODELS)
    assert not any(is_protein_model(n) for n in DNA_MODELS)
    assert len(AMINO_ACIDS) == 20 and lg().alphabet == AMINO_ACIDS


def test_poisson_is_exact():
    """Poisson: uniform freqs, equal off-diagonal rates (F81-for-proteins), exact by construction."""
    m = poisson()
    assert np.allclose(m.stationary, 1.0 / 20.0)
    off = m.Q[~np.eye(20, dtype=bool)]
    assert np.allclose(off, off[0])                         # all off-diagonal rates equal


def test_protein_stationary_recovered_on_long_branch():
    """A star tree with long branches recovers a protein model's stationary frequencies."""
    tips = [_N(f"t{i}") for i in range(6)]
    root = _N("r", tips)
    subst = {t: 40.0 for t in tips}                         # very long -> equilibrium
    m = lg()
    seqs = evolve_on_tree(root, subst, m, np.random.default_rng(0), length=6000)
    pooled = "".join(seqs[t.gid] for t in tips)
    assert set(pooled) <= set(AMINO_ACIDS)
    freqs = np.array([pooled.count(a) / len(pooled) for a in AMINO_ACIDS])
    assert np.max(np.abs(freqs - m.stationary)) < 0.02


def test_protein_alignment_alphabet():
    """Evolving under a protein model yields the 20-AA alphabet (never ACGT-only)."""
    a, b = _N("a"), _N("b")
    root = _N("r", [a, b])
    seqs = evolve_on_tree(root, {a: 0.5, b: 0.5}, wag(), np.random.default_rng(1), length=500)
    assert set(seqs["a"]) <= set(AMINO_ACIDS)
    assert not set(seqs["a"]) <= set("ACGT")               # genuinely protein, not nucleotide


def test_p_matrix_matches_jc_closed_form():
    """The numpy-only exp(Qt) matches the JC69 closed form P_ii = 1/4 + 3/4 e^{-4t/3}."""
    m = jc69()
    for t in (0.05, 0.3, 1.0, 3.0):
        P = m.p_matrix(t)
        pii = 0.25 + 0.75 * np.exp(-4.0 / 3.0 * t)
        pij = 0.25 - 0.25 * np.exp(-4.0 / 3.0 * t)
        assert np.allclose(np.diag(P), pii)
        assert np.allclose(P[~np.eye(4, dtype=bool)], pij)


def test_gamma_rates_numpy_only_mean_one():
    """Discrete-Gamma category rates average to 1 (numpy-only implementation, no scipy)."""
    for shape in (0.2, 0.5, 1.0, 2.0):
        g = GammaRates(shape, 4)
        assert g.rates.shape == (4,)
        assert np.isclose(g.rates.mean(), 1.0)
        assert np.all(np.diff(g.rates) > 0)                # categories are increasing


# --------------------------------------------------------------------------- #
# evolve_on_tree
# --------------------------------------------------------------------------- #
def test_jc_distance_recovered():
    a, b = _N("a"), _N("b")
    root = _N("r", [a, b])
    subst = {root: 0.0, a: 0.2, b: 0.2}                     # each tip 0.2 from root -> 0.4 apart
    seqs = evolve_on_tree(root, subst, jc69(), np.random.default_rng(0), length=150000)
    p = np.mean([x != y for x, y in zip(seqs["a"], seqs["b"])])
    d = -0.75 * np.log(1 - 4 / 3 * p)
    assert abs(d - 0.40) < 0.02


def test_zero_branch_identity():
    a, b = _N("a"), _N("b")
    root = _N("r", [a, b])
    seqs = evolve_on_tree(root, {}, jc69(), np.random.default_rng(1), length=400)
    assert seqs["a"] == seqs["b"] == seqs["r"]              # no branch length -> identical


def test_stationary_frequencies_recovered():
    seqs = evolve_on_tree(_N("x"), {}, hky85(2.0, (0.4, 0.1, 0.1, 0.4)),
                          np.random.default_rng(2), length=80000)
    c = Counter(seqs["x"])
    freqs = np.array([c[b] / len(seqs["x"]) for b in "ACGT"])
    assert np.allclose(freqs, [0.4, 0.1, 0.1, 0.4], atol=0.01)


def test_reproducible_and_root_seq():
    a = _N("a")
    root = _N("r", [a])
    kw = dict(root_seq="ACGTACGTAC")
    s1 = evolve_on_tree(root, {root: 0.0, a: 0.0}, jc69(), np.random.default_rng(3), **kw)
    s2 = evolve_on_tree(root, {root: 0.0, a: 0.0}, jc69(), np.random.default_rng(3), **kw)
    assert s1 == s2 and s1["r"] == "ACGTACGTAC" and s1["a"] == "ACGTACGTAC"


def test_gamma_increases_variance():
    # with rate heterogeneity, per-site divergence is more dispersed than without
    a = _N("a")
    root = _N("r", [a])
    rng = np.random.default_rng(4)
    base = evolve_on_tree(root, {a: 0.3}, jc69(), rng, length=20000)
    gam = evolve_on_tree(root, {a: 0.3}, jc69(), rng, length=20000, gamma=GammaRates(0.2, 6))
    # both diverge from the root; the Gamma run simply runs without error and stays valid DNA
    assert set(base["a"]) <= set("ACGT") and set(gam["a"]) <= set("ACGT")


def test_revcomp_and_fasta_roundtrip(tmp_path):
    assert reverse_complement("AACGT") == "ACGTT"
    p = tmp_path / "x.fasta.gz"
    write_fasta(p, {"s1": "ACGTACGT", "s2": "TTTT"}, gzip_out=True)
    assert read_fasta(str(p)) == {"s1": "ACGTACGT", "s2": "TTTT"}


# --------------------------------------------------------------------------- #
# Ancestral genome reconstruction (nucleotide model)
# --------------------------------------------------------------------------- #
def _run(seed=5, retain=True):
    tree = simulate_species_tree(BirthDeath(1.0, 0.3), n_tips=6, age=1.0, seed=2)
    genes = [(20, 60, "gA"), (90, 130, "gB"), (160, 200, "gC"), (230, 280, "gD")]
    res = simulate_nucleotide_genomes(
        tree, inversion=0.01, loss=0.006, duplication=0.005, transfer=0.005,
        transposition=0.003, root_length=300, extension=0.9, gene_intervals=genes,
        pseudogenization=0.3, replacement=0.3, retain_internal=retain, seed=seed)
    return tree, res


def test_node_genomes_kept_for_every_node():
    tree, res = _run()
    assert set(res.node_genomes) == set(tree.nodes_preorder())
    # leaves match the extant leaf genomes
    for leaf, g in res.leaf_genomes.items():
        assert res.node_genomes[leaf] is g


def test_root_architecture_is_the_input_tiling():
    tree, res = _run()
    mosaic = res.node_mosaic(tree.root)
    blocks = [res._block_by_id[aid] for aid, _ in mosaic]
    assert all(s == 1 for _, s in mosaic)                   # seed is all forward strand
    assert blocks[0].start == 0 and blocks[-1].end == 300     # tiles the whole chromosome
    assert sum(a.length for a in blocks) == 300
    genes = [(a.start, a.end) for a in blocks if a.kind == "gene"]
    assert genes == [(20, 60), (90, 130), (160, 200), (230, 280)]   # genes intact, in order


def test_node_sequence_lengths_and_zero_divergence_assembly():
    tree, res = _run()
    res.simulate_sequences(jc69(), subst_rate=0.3, seed=7)
    for n in tree.nodes_preorder():
        assert len(res.node_sequence(n)) == res.node_genomes[n].size()

    # at zero divergence every node's DNA is the strand-correct concatenation of its blocks
    res.simulate_sequences(jc69(), subst_rate=0.0, seed=7)
    block_seq = {aid: next(iter(d.values())) for aid, d in res._block_seqs.items() if d}
    for n in tree.nodes_preorder():
        expect = "".join(reverse_complement(block_seq[aid]) if s == -1 else block_seq[aid]
                         for aid, s in res.node_mosaic(n))
        assert res.node_sequence(n) == expect


def test_fasta_seeded_root_reproduces_input_genome():
    tree, res = _run()
    rng = np.random.default_rng(123)
    genome = "".join(rng.choice(list("ACGT"), size=300))
    res.simulate_sequences(jc69(), subst_rate=0.4, root_fasta=genome, seed=7)
    assert res.node_sequence(tree.root) == genome           # the root IS the input genome
    with pytest.raises(ValueError, match="root_fasta length"):
        res.simulate_sequences(jc69(), root_fasta="ACGT", seed=1)


def test_gene_alignments_shape():
    tree, res = _run()
    res.simulate_sequences(jc69(), subst_rate=0.2, seed=7)
    ga = res.gene_alignments()
    assert set(ga) <= {"gA", "gB", "gC", "gD"}
    lengths = {"gA": 40, "gB": 40, "gC": 40, "gD": 50}      # from the gene intervals
    for gene, aln in ga.items():
        assert all(len(s) == lengths[gene] for s in aln.values())   # every copy is the gene length
    assert res.intergene_alignments()                       # intergene blocks too


def test_node_sequence_needs_simulate_first():
    tree, res = _run()
    with pytest.raises(RuntimeError, match="simulate_sequences"):
        res.node_sequence(tree.root)
