"""Founding sequences and per-family rates that are **given** — ``founding=`` and
``family_multipliers=`` on `simulate_sequences`.

A run draws both by default: each family starts from a draw at the model's stationary frequencies,
and a rate that varies among families draws one multiplier apiece. Both are also written, as
``sequences_founding.fasta`` and ``family_multipliers.tsv``. These check the other direction: that a
run can take what a run writes, so the genes of a real genome can start from their real sequences
and run at the rates measured for them.
"""

import re

import pytest

from zombi2 import sequences
from zombi2.genomes import family, simulate_genomes_family, simulate_genomes_nucleotide
from zombi2.params import LogNormal, PerSite
from zombi2.sequences import simulate_sequences
from zombi2.sequences.substitution_models import AMINO_ACIDS, hky85, jc69, lg
from zombi2.species import simulate_species_tree

REAL = {0: "ATGGCGAAGCTGCGCGAGGCCATCGAGCGCCTGGAGCAGGCCTGA",
        1: "ATGAGCGAACTGACCGTGAAAGACCTGATCGCGCGCATCTGA",
        2: "ATGACCCGCATCGAGCTGCTGCAGAAGCTCGGCGTGTGA"}


def _genomes(*, seed=1, n_extant=8, **kw):
    sp = simulate_species_tree(birth=1.0, death=0.2, n_extant=n_extant, seed=seed)
    params = dict(initial_families=3, duplication=0.3, loss=0.2, transfer=0.3, seed=seed)
    params.update(kw)
    return simulate_genomes_family(sp, **params)


def _length(newick):
    return sum(float(x) for x in re.findall(r":([0-9.eE+-]+)", newick))


# --- founding ---------------------------------------------------------------------------------------

def test_a_family_starts_from_the_sequence_it_is_given():
    # at rate zero nothing changes, so every copy of every family is still the sequence handed in
    g = _genomes()
    r = simulate_sequences(g, model=hky85(), founding=REAL, substitution=0.0, seed=1)
    assert r.founding == REAL
    for fam, aln in r.alignments.items():
        assert aln and set(aln.values()) == {REAL[fam]}
        assert set(r.ancestral[fam].values()) == {REAL[fam]}


def test_the_given_sequence_is_where_the_family_begins_and_then_it_evolves():
    g = _genomes()
    r = simulate_sequences(g, model=hky85(), founding=REAL, divergence=0.3, seed=1)
    assert r.founding == REAL
    for fam, aln in r.alignments.items():
        assert {len(s) for s in aln.values()} == {len(REAL[fam])}     # its length is the family's
        assert set(aln.values()) != {REAL[fam]}


def test_families_not_given_a_sequence_still_draw_one():
    g = _genomes()
    r = simulate_sequences(g, model=hky85(), founding={1: REAL[1]}, length=30, seed=1)
    assert r.founding[1] == REAL[1]
    assert len(r.founding[0]) == len(r.founding[2]) == 30
    with pytest.raises(ValueError, match="length= is missing"):
        simulate_sequences(g, model=hky85(), founding={1: REAL[1]}, seed=1)


def test_a_family_is_named_by_number_by_file_spelling_or_by_declared_name():
    sp = simulate_species_tree(birth=1.0, death=0.2, n_extant=6, seed=2)
    g = simulate_genomes_family(sp, families=[family("rpoB"), family("dnaA"), family("gyrA")], seed=2)
    rpob, dnaa, gyra = (g.family_names[n] for n in ("rpoB", "dnaA", "gyrA"))
    # the run also holds the unnamed families of the initial genome, which draw at `length`
    r = simulate_sequences(g, model=hky85(), substitution=0.0, length=30, seed=1,
                           founding={"rpoB": REAL[0], f"fam{dnaa}": REAL[1], gyra: REAL[2]})
    assert (r.founding[rpob], r.founding[dnaa], r.founding[gyra]) == (REAL[0], REAL[1], REAL[2])
    with pytest.raises(ValueError, match="twice"):
        simulate_sequences(g, model=hky85(), founding={"rpoB": REAL[0], rpob: REAL[0]}, length=9)
    with pytest.raises(ValueError, match="declared name"):
        simulate_sequences(g, model=hky85(), founding={"recA": REAL[0]}, length=9)


def test_a_fasta_is_read_as_given_and_lower_case_is_taken(tmp_path):
    g = _genomes()
    path = tmp_path / "real.fasta"
    path.write_text(f">fam0 a real gene\n{REAL[0][:20]}\n{REAL[0][20:].lower()}\n>2\n{REAL[2]}\n",
                    encoding="utf-8")
    r = simulate_sequences(g, model=hky85(), founding=path, length=12, substitution=0.0, seed=1)
    assert r.founding[0] == REAL[0] and r.founding[2] == REAL[2] and len(r.founding[1]) == 12


def test_a_protein_run_takes_protein_sequences_under_a_mixture_too():
    # a profile mixture lays a family's sites out by the profile each drew, so a given sequence has
    # to follow its sites into that layout and come back out in order
    g = _genomes()
    protein = {f: "".join(AMINO_ACIDS[(7 * i + 3 * f) % 20] for i in range(60 + f)) for f in range(3)}
    for model in (lg(), lg().profile_mixture("C10")):
        r = simulate_sequences(g, model=model, founding=protein, substitution=0.0, seed=1)
        assert r.founding == protein
        assert all(set(aln.values()) == {protein[fam]} for fam, aln in r.alignments.items())
    r = simulate_sequences(g, model=lg().profile_mixture("C10"), founding=protein, divergence=0.4,
                           seed=1)
    assert r.founding == protein
    assert all(len(r.profile_assignments[fam]) == len(protein[fam]) for fam in protein)


def test_partitions_fix_the_length_so_a_given_sequence_has_to_match():
    g = _genomes()
    parts = [(hky85(), 20), (jc69(), 25)]
    r = simulate_sequences(g, partitions=parts, founding={0: REAL[0]}, substitution=0.0, seed=1)
    assert r.founding[0] == REAL[0] and len(r.founding[1]) == 45
    with pytest.raises(ValueError, match="add up to 45"):
        simulate_sequences(g, partitions=parts, founding={1: REAL[1]}, seed=1)


# --- family multipliers -----------------------------------------------------------------------------

def test_a_given_multiplier_is_the_speed_the_family_runs_at():
    g = _genomes()
    base = simulate_sequences(g, model=hky85(), length=40, seed=1)
    table = {0: 0.2, 1: 1.0, 2: 3.0}
    r = simulate_sequences(g, model=hky85(), length=40, family_multipliers=table, seed=1)
    assert {f: v["substitution"] for f, v in r.family_multipliers.items()} == table
    for fam, speed in table.items():
        ratio = (_length(r.phylograms[fam]["complete"])
                 / _length(base.phylograms[fam]["complete"]))
        assert ratio == pytest.approx(speed, rel=1e-5)
    assert not base.family_multipliers            # a run that is given none still writes none


def test_a_multiplier_of_one_for_every_family_is_the_run_without_any():
    g = _genomes()
    base = simulate_sequences(g, model=hky85(), length=40, seed=1)
    r = simulate_sequences(g, model=hky85(), length=40, seed=1,
                           family_multipliers={f: 1.0 for f in g.gene_trees})
    assert r.alignments == base.alignments and r.phylograms == base.phylograms


def test_a_listed_family_takes_its_number_and_the_others_keep_their_draw():
    g = _genomes()
    rate = PerSite(1.0).varying_among("families", LogNormal(0.0, 0.8))
    drawn = simulate_sequences(g, model=hky85(), length=40, substitution=rate, seed=1)
    r = simulate_sequences(g, model=hky85(), length=40, substitution=rate, seed=1,
                           family_multipliers={1: 0.25})
    assert r.family_multipliers[1]["substitution"] == 0.25
    for fam in (0, 2):
        assert r.family_multipliers[fam] == drawn.family_multipliers[fam]
    # and with nothing drawn, an unlisted family runs at the run's own rate
    alone = simulate_sequences(g, model=hky85(), length=40, seed=1, family_multipliers={1: 0.25})
    assert {f: v["substitution"] for f, v in alone.family_multipliers.items()} == \
           {0: 1.0, 1: 0.25, 2: 1.0}


def test_both_shapes_a_result_and_a_caller_write_are_read():
    g = _genomes()
    flat = simulate_sequences(g, model=hky85(), length=40, seed=1, family_multipliers={"fam0": 2})
    nested = simulate_sequences(g, model=hky85(), length=40, seed=1,
                                family_multipliers={0: {"substitution": 2.0}, 1: {"substitution": None}})
    assert flat.family_multipliers == nested.family_multipliers
    assert flat.alignments == nested.alignments


# --- what a run writes, a run can take --------------------------------------------------------------

def test_a_run_reads_back_the_two_files_another_run_wrote(tmp_path):
    g = _genomes()
    rate = PerSite(1.0).varying_among("families", LogNormal(0.0, 0.8))
    first = simulate_sequences(g, model=hky85(), length=40, substitution=rate, seed=1)
    first.write(tmp_path, outputs=("founding", "family_multipliers"))
    again = simulate_sequences(g, model=hky85(), seed=2,
                               founding=tmp_path / "sequences_founding.fasta",
                               family_multipliers=tmp_path / "family_multipliers.tsv")
    assert again.founding == first.founding
    assert again.family_multipliers == first.family_multipliers
    assert again.alignments != first.alignments                  # a new run from the same start


@pytest.mark.parametrize("model", [hky85(), None])
def test_every_engine_starts_from_the_same_given_sequences(tmp_path, model):
    g = _genomes()
    if model is None:                       # the per-family layout of a mixture, on each engine
        model = lg().profile_mixture("C10")
        given = {f: "".join(AMINO_ACIDS[(5 * i + f) % 20] for i in range(50 + f)) for f in range(3)}
    else:
        given = REAL
    kw = dict(model=model, founding=given, family_multipliers={0: 0.5, 2: 2.0}, divergence=0.2,
              seed=4)
    serial = simulate_sequences(g, **kw)
    one, two = simulate_sequences(g, parallel=1, **kw), simulate_sequences(g, parallel=2, **kw)
    assert serial.founding == one.founding == two.founding == given
    assert one.alignments == two.alignments and one.phylograms == two.phylograms
    assert serial.family_multipliers == one.family_multipliers == two.family_multipliers
    assert {len(s) for s in serial.alignments[1].values()} == {len(given[1])}
    streamed = simulate_sequences(g, stream_to=tmp_path / "s", flat=True,
                                  outputs=("founding", "family_multipliers", "alignments"),
                                  parallel=1, **kw)
    assert streamed.n_families == 3
    text = (tmp_path / "s" / "sequences_founding.fasta").read_text(encoding="utf-8")
    assert all(f">fam{f}\n{seq}\n" in text for f, seq in given.items())


# --- refusals ---------------------------------------------------------------------------------------

@pytest.mark.parametrize("kw, match", [
    (dict(founding={0: "ATGNNN"}, length=9), "not in the model's alphabet"),
    (dict(founding={7: "ATG"}, length=9), "does not have"),
    (dict(founding={}, length=9), "is empty"),
    (dict(founding={0: ""}, length=9), "no sequence"),
    (dict(founding=["ATG"], length=9), "takes"),
    (dict(founding=REAL, insertion=0.1), "indels"),
    (dict(family_multipliers={0: 0.0}, length=9), "positive number"),
    (dict(family_multipliers={0: float("nan")}, length=9), "positive number"),
    (dict(family_multipliers={0: "fast"}, length=9), "positive number"),
    (dict(family_multipliers={9: 2.0}, length=9), "does not have"),
    (dict(family_multipliers={}, length=9), "no family"),
])
def test_what_is_refused(kw, match):
    with pytest.raises((ValueError, TypeError), match=match):
        simulate_sequences(_genomes(), model=hky85(), seed=1, **kw)


def test_a_nucleotide_run_has_no_families_to_give_them_to():
    sp = simulate_species_tree(birth=1.0, death=0.0, n_extant=4, seed=1)
    g = simulate_genomes_nucleotide(sp, root_length=600, genes=2, gene_length=100, seed=1)
    with pytest.raises(ValueError, match="fasta="):
        simulate_sequences(g, model=jc69(), founding={0: "ACGT"}, seed=1)
    with pytest.raises(ValueError, match="blocks"):
        simulate_sequences(g, model=jc69(), family_multipliers={0: 2.0}, seed=1)


def test_a_joint_run_is_written_gene_by_gene():
    sp = simulate_species_tree(birth=1.0, death=0.0, n_extant=4, seed=1)
    g = simulate_genomes_family(sp, families=[family("a"), family("b")], seed=1)
    genes = [sequences.gene(name="a", model=lg(), length=30),
             sequences.gene(name="b", model=lg(), length=30)]
    with pytest.raises(ValueError, match="founding"):
        simulate_sequences(g, joint=True, genes=genes, founding={0: "ACGT"}, seed=1)


# --- the command line -------------------------------------------------------------------------------

def test_the_command_line_takes_both_files(tmp_path):
    from zombi2.cli.main import main

    run = tmp_path / "r"
    assert main(["species", str(run), "--birth", "1", "--death", "0.2", "--n-extant", "6",
                 "--seed", "1", "--quiet"]) == 0
    assert main(["genomes", str(run), "--initial-families", "3", "--duplication", "0.3",
                 "--loss", "0.2", "--seed", "1", "--quiet"]) == 0
    (tmp_path / "real.fasta").write_text(f">fam0\n{REAL[0]}\n>1\n{REAL[1]}\n", encoding="utf-8")
    (tmp_path / "rates.tsv").write_text("family\tsubstitution\n0\t0.2\n2\t3.0\n", encoding="utf-8")
    assert main(["sequences", str(run), "--model", "hky85", "--length", "30",
                 "--founding", str(tmp_path / "real.fasta"),
                 "--family-multipliers", str(tmp_path / "rates.tsv"),
                 "--write", "founding", "family_multipliers", "--seed", "1", "--quiet"]) == 0
    out = run / "sequences"
    founding = (out / "sequences_founding.fasta").read_text(encoding="utf-8")
    assert f">fam0\n{REAL[0]}\n" in founding and f">fam1\n{REAL[1]}\n" in founding
    assert (out / "family_multipliers.tsv").read_text(encoding="utf-8") == \
        "family\tsubstitution\n0\t0.2\n1\t1.0\n2\t3.0\n"
    log = (out / "sequences.log").read_text(encoding="utf-8")
    assert "real.fasta" in log and "rates.tsv" in log
