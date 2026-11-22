"""Codon models — ``gy94()``, a model whose states are codons and whose sequences are DNA.

A nucleotide model treats the three positions of a codon alike, so a gene evolved under one gains
stop codons and amino-acid changes far faster than a real gene does. A codon model has no stop
codon among its states, and ``omega`` says how often a change that alters the amino acid gets
through. These check the model itself, that everything a run reads and writes is still DNA, and
what is refused. Whether the *process* is the one Goldman and Yang wrote down is checked against
closed forms in ``test_validation.py``.
"""

import numpy as np
import pytest

from zombi2 import sequences
from zombi2.genomes import simulate_genomes_family, simulate_genomes_nucleotide
from zombi2.params import Clade
from zombi2.sequences import gy94, hky85, simulate_sequences
from zombi2.sequences.lineage_models import Models
from zombi2.sequences.substitution_models import GENETIC_CODES, letters_of, lg, states_of
from zombi2.species import simulate_species_tree

_TCAG = [a + b + c for a in "TCAG" for b in "TCAG" for c in "TCAG"]
AMINO = {code: dict(zip(_TCAG, table)) for code, table in GENETIC_CODES.items()}
REAL = "ATGGCGAAGCTGCGCGAGGCCATCGAGCGCCTGGAGCAGGCC"          # a gene, without its stop codon


def _genomes(*, seed=1, n_extant=10, families=12):
    sp = simulate_species_tree(birth=1.0, death=0.2, n_extant=n_extant, seed=seed)
    return simulate_genomes_family(sp, initial_families=families, duplication=0.2, loss=0.2,
                                   seed=seed)


def _codons(seq):
    return [seq[i:i + 3] for i in range(0, len(seq), 3)]


# --- the model ---------------------------------------------------------------------------------------

@pytest.mark.parametrize("code, states, tga", [(11, 61, "*"), (4, 62, "W")])
def test_the_states_are_the_codons_that_are_not_stops(code, states, tga):
    m = gy94(omega=0.1, code=code)
    assert m.k == len(m.codons) == states and m.alphabet == "ACGT"
    assert list(m.codons) == sorted(c for c in _TCAG if AMINO[code][c] != "*")
    assert AMINO[code]["TGA"] == tga and ("TGA" in m.codons) == (tga != "*")
    assert not hky85().codons and not lg().codons


def test_a_change_is_one_base_and_its_rate_says_what_kind():
    kappa, omega = 3.0, 0.2
    m = gy94(omega=omega, kappa=kappa)
    pi = m.stationary
    scale = m.Q[m.codons.index("AAA"), m.codons.index("AAC")] / pi[m.codons.index("AAC")]
    for i, x in enumerate(m.codons):
        for j, y in enumerate(m.codons):
            if i == j:
                continue
            changed = [(a, b) for a, b in zip(x, y) if a != b]
            rate = m.Q[i, j] / pi[j] / scale
            if len(changed) != 1:
                assert rate == 0.0                                  # never two bases at once
                continue
            transition = set(changed[0]) in ({"A", "G"}, {"C", "T"})
            same_protein = AMINO[11][x] == AMINO[11][y]
            # AAA→AAC is a transversion that changes the amino acid, so it is the unit: omega
            expected = (kappa if transition else 1.0) * (1.0 if same_protein else omega) / omega
            assert rate == pytest.approx(expected)


def test_a_branch_length_is_substitutions_per_base():
    # every change is one base, so three changes per codon is one per base — the unit of hky85()
    for m in (gy94(omega=0.1), gy94(omega=1.0, kappa=4.0, frequencies=(0.1, 0.4, 0.3, 0.2), code=4)):
        assert -(m.stationary * np.diag(m.Q)).sum() == pytest.approx(3.0)


def test_frequencies_are_four_bases_or_one_per_codon():
    assert gy94(omega=0.1).stationary == pytest.approx(np.full(61, 1 / 61))
    bases = dict(zip("ACGT", (0.1, 0.4, 0.3, 0.2)))
    m = gy94(omega=0.1, frequencies=tuple(bases.values()))
    product = np.array([bases[c[0]] * bases[c[1]] * bases[c[2]] for c in m.codons])
    assert m.stationary == pytest.approx(product / product.sum())
    own = np.arange(1, 62, dtype=float)
    assert gy94(omega=0.1, frequencies=own / own.sum()).stationary == pytest.approx(own / own.sum())


@pytest.mark.parametrize("kw, match", [
    (dict(omega=0), "omega must be a positive number"),
    (dict(omega=0.1, kappa=-1), "kappa must be a positive number"),
    (dict(omega=0.1, code=1), "NCBI translation table"),
    (dict(omega=0.1, frequencies=(0.5, 0.5)), "4 base frequencies"),
    (dict(omega=0.1, frequencies=(0.5, 0.5, 0.5, 0.5)), "sum to 1"),
])
def test_what_the_constructor_refuses(kw, match):
    with pytest.raises(ValueError, match=match):
        gy94(**kw)
    with pytest.raises(TypeError):
        gy94()                                   # omega has no default: nothing is assumed


def test_states_and_letters_go_both_ways():
    m = gy94(omega=0.1)
    states = states_of(m, REAL)
    assert len(states) == len(REAL) // 3 and letters_of(m, states) == REAL
    assert letters_of(m, np.stack([states, states])) == REAL + REAL
    for bad, match in (("ATGTAA", "stop codon"), ("ATGG", "whole number of codons"),
                       ("ATGNNN", "not three letters")):
        with pytest.raises(ValueError, match=match):
            states_of(m, bad)
    assert states_of(gy94(omega=0.1, code=4), "TGA").tolist() == [gy94(omega=0.1, code=4).codons.index("TGA")]


# --- a run -------------------------------------------------------------------------------------------

@pytest.mark.parametrize("code", [11, 4])
def test_a_run_writes_dna_with_no_stop_codon(code):
    g = _genomes()
    r = simulate_sequences(g, model=gy94(omega=0.2, code=code), length=300, divergence=0.3, seed=1)
    assert r.alphabet == "ACGT"
    every = [s for fam in r.alignments for s in (*r.alignments[fam].values(),
                                                 *r.ancestral[fam].values(), r.founding[fam])]
    assert len(every) > 100 and {len(s) for s in every} == {300}
    assert all(set(s) <= set("ACGT") for s in every)
    assert not any(AMINO[code][c] == "*" for s in every for c in _codons(s))
    r.gc()                                              # GC content is DNA's, and this is DNA
    assert r.summary()["sites"] == {"min": 300, "max": 300}


def test_omega_decides_how_often_the_protein_changes():
    g = _genomes(families=30)
    share = {}
    for omega in (0.05, 1.0):
        r = simulate_sequences(g, model=gy94(omega=omega), length=600, divergence=0.05, seed=1)
        same = changed = 0
        for fam, aln in r.alignments.items():
            start = _codons(r.founding[fam])
            for s in aln.values():
                for a, b in zip(start, _codons(s)):
                    if a != b:
                        changed += 1
                        same += AMINO[11][a] == AMINO[11][b]
        share[omega] = same / changed
    assert share[0.05] > 0.8 and share[1.0] < 0.35


def test_the_same_divergence_gives_the_identity_a_nucleotide_model_gives():
    # the reason a branch length is per base: `divergence` then means what it means for hky85()
    g = _genomes(families=30)
    codon = simulate_sequences(g, model=gy94(omega=1.0), length=600, divergence=0.05, seed=1)
    bases = simulate_sequences(g, model=hky85(), length=600, divergence=0.05, seed=1)
    a = sequences.mean_pairwise_identity(codon.alignments)
    b = sequences.mean_pairwise_identity(bases.alignments)
    assert abs(a - b) < 0.01 and 0.85 < a < 0.97


def test_a_length_counts_bases_and_is_a_whole_number_of_codons():
    g = _genomes()
    with pytest.raises(ValueError, match="multiple of 3"):
        simulate_sequences(g, model=gy94(omega=0.1), length=100, seed=1)
    r = simulate_sequences(g, partitions=[(gy94(omega=0.05), 90), (gy94(omega=1.0), 30)], seed=1)
    assert {len(s) for s in r.founding.values()} == {120}
    with pytest.raises(ValueError, match="not have the same states"):
        simulate_sequences(g, partitions=[(gy94(omega=0.1), 90), (hky85(), 30)], seed=1)
    with pytest.raises(ValueError, match="not have the same states"):
        simulate_sequences(g, partitions=[(gy94(omega=0.1), 90), (gy94(omega=0.1, code=4), 30)],
                           seed=1)


def test_a_gene_starts_from_its_real_sequence():
    g = _genomes()
    r = simulate_sequences(g, model=gy94(omega=0.1), founding={0: REAL, 1: REAL.lower()},
                           length=30, substitution=0.0, seed=1)
    assert r.founding[0] == r.founding[1] == REAL
    assert set(r.alignments[0].values()) == {REAL} and len(r.founding[2]) == 30
    for bad, match in ((REAL + "TAA", "stop codon"), (REAL + "A", "whole number of codons")):
        with pytest.raises(ValueError, match=match):
            simulate_sequences(g, model=gy94(omega=0.1), founding={0: bad}, length=30, seed=1)


def test_the_engines_and_the_decorations_take_a_codon_model(tmp_path):
    g = _genomes()
    model = gy94(omega=0.1).across_sites(gamma_shape=0.5, invariant=0.1)
    kw = dict(model=model, length=150, divergence=0.2, seed=3)
    one, two = simulate_sequences(g, parallel=1, **kw), simulate_sequences(g, parallel=2, **kw)
    assert one.alignments == two.alignments and one.founding == two.founding
    serial = simulate_sequences(g, **kw)
    assert {len(s) for a in serial.alignments.values() for s in a.values()} == {150}
    streamed = simulate_sequences(g, stream_to=tmp_path / "s", flat=True, parallel=1, **kw)
    assert streamed.sites == 150
    # a clade under its own selection: same codons, another omega
    sp = g.complete_tree
    clade = next(i for i, n in sp.nodes.items() if n.children and n.parent is not None)
    label = sp.labels()[clade]
    per_clade = Models().set_by(Clade({"relaxed": label}),
                                {"relaxed": gy94(omega=1.0), "rest": gy94(omega=0.05)})
    r = simulate_sequences(g, model=per_clade, length=150, divergence=0.2, seed=3)
    assert not any(AMINO[11][c] == "*" for a in r.alignments.values() for s in a.values()
                   for c in _codons(s))
    with pytest.raises(ValueError, match="same states"):
        Models().set_by(Clade({"relaxed": label}), {"relaxed": gy94(omega=1.0), "rest": hky85()})


# --- refusals ----------------------------------------------------------------------------------------

@pytest.mark.parametrize("kw, match", [
    (dict(length=300, insertion=0.1), "indels are not read with a codon model"),
    (dict(length=300, deletion=0.1), "indels are not read with a codon model"),
    (dict(length=300, record=True), "record=True is not read with a codon model"),
    (dict(length=300, profiles={0: np.full((100, 4), 0.25)}), "codon model"),
])
def test_what_a_codon_run_refuses(kw, match):
    with pytest.raises(ValueError, match=match):
        simulate_sequences(_genomes(), model=gy94(omega=0.1), seed=1, **kw)


def test_a_nucleotide_genome_run_does_not_read_a_codon_model_yet():
    sp = simulate_species_tree(birth=1.0, death=0.0, n_extant=4, seed=1)
    g = simulate_genomes_nucleotide(sp, root_length=600, genes=2, gene_length=99, seed=1)
    with pytest.raises(ValueError, match="not read on a nucleotide genome run"):
        simulate_sequences(g, model=gy94(omega=0.1), seed=1)
    with pytest.raises(ValueError, match="not read on a nucleotide genome run"):
        simulate_sequences(g, model=hky85(), intergene_model=gy94(omega=1.0), seed=1)


def test_a_joint_run_does_not_read_a_codon_model_yet():
    with pytest.raises(ValueError, match="joint run does not read"):
        sequences.gene(name="a", model=gy94(omega=0.1), length=300)
    with pytest.raises(ValueError, match="amino-acid profiles"):
        gy94(omega=0.1).profile_mixture("C10")


# --- the command line --------------------------------------------------------------------------------

def test_the_command_line_runs_a_codon_model(tmp_path, capsys):
    from zombi2.cli.main import main

    run = tmp_path / "r"
    assert main(["species", str(run), "--birth", "1", "--death", "0.2", "--n-extant", "6",
                 "--seed", "1", "--quiet"]) == 0
    assert main(["genomes", str(run), "--initial-families", "6", "--seed", "1", "--quiet"]) == 0
    assert main(["sequences", str(run), "--model", "gy94", "--omega", "0.1", "--kappa", "3",
                 "--code", "4", "--divergence", "0.1", "--write", "founding", "alignments",
                 "--seed", "1", "--quiet"]) == 0
    founding = (run / "sequences" / "sequences_founding.fasta").read_text(encoding="utf-8")
    seqs = ["".join(record.splitlines()[1:]) for record in founding.split(">")[1:]]
    assert seqs and {len(s) for s in seqs} == {999}              # the default, in whole codons
    assert not any(AMINO[4][c] == "*" for s in seqs for c in _codons(s))
    log = (run / "sequences" / "sequences.log").read_text(encoding="utf-8")
    assert "omega\t0.1" in log and "code\t4" in log and "kappa\t3.0" in log
    capsys.readouterr()
    for extra, match in ((["--model", "gy94"], "needs --omega"),
                         (["--model", "hky85", "--omega", "0.1"], "don't apply"),
                         (["--model", "jc69", "--code", "4"], "don't apply")):
        with pytest.raises(SystemExit):
            main(["sequences", str(run), "--force", "--seed", "1", *extra])
        assert match in capsys.readouterr().err
