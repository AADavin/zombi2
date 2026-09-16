"""A substitution rate drawn among families (issue #443).

``substitution = PerSite(1.0).varying_among('families', law)`` draws one multiplier per gene
family, before any site evolves, and the family keeps that value for its whole life. The clock rides lineages
and this rides families, so the two are separate axes and multiply::

    branch = substitution · Δt · lineage clock · family multiplier

The run keeps the multipliers on its result and writes them in ``family_multipliers.tsv``, the one
genomes level writes its event-rate multipliers in. A run whose rate does not vary among families
writes the header alone, and takes no draw at all.
"""

import re

import numpy as np
import pytest

from zombi2 import sequences
from zombi2.genomes import family, simulate_genomes_family, simulate_genomes_nucleotide
from zombi2.genomes.multipliers import SEQUENCE_TARGETS, multipliers_from_tsv
from zombi2.params import Drift, LogNormal, PerSite
from zombi2.sequences import jc69
from zombi2.sequences.multipliers import family_multipliers as draw_multipliers
from zombi2.species import simulate_species_tree

HEADER = "family\tsubstitution\n"


@pytest.fixture(scope="module")
def genomes():
    tree = simulate_species_tree(birth=1.0, death=0.3, n_extant=12, seed=1).complete_tree
    return simulate_genomes_family(tree, duplication=0.15, transfer=0.05, loss=0.15,
                                   origination=0.5, initial_families=20, seed=2)


def _varying(base=0.05, sigma=0.8):
    return PerSite(base).varying_among("families", LogNormal(0.0, sigma))


def _lengths(newick):
    """Every branch length in a Newick string, in the order they are written."""
    return [float(x) for x in re.findall(r":([0-9.eE+-]+)", newick)]


def _ranks(values):
    values = np.asarray(values, dtype=float)
    ranks = np.empty(len(values))
    ranks[values.argsort(kind="stable")] = np.arange(len(values))
    return ranks


def test_a_run_whose_rate_does_not_vary_writes_the_header_alone(genomes, tmp_path):
    run = sequences.simulate_sequences(genomes, model=jc69(), length=40, substitution=0.05, seed=3)
    assert run.family_multipliers == {}
    run.write(tmp_path)
    assert (tmp_path / "family_multipliers.tsv").read_text(encoding="utf-8") == HEADER


def test_a_rate_that_does_not_vary_among_families_takes_no_draw():
    """The randomness a run consumes is what says whether a feature is really absent: a rate
    carrying no per-family modifier must leave the generator exactly where it found it."""
    untouched, used = np.random.default_rng(0), np.random.default_rng(0)
    assert draw_multipliers((), range(50), used) == {}
    assert untouched.random() == used.random()


def test_every_family_has_a_row(genomes):
    run = sequences.simulate_sequences(genomes, model=jc69(), length=40, seed=3,
                                       substitution=_varying())
    assert sorted(run.family_multipliers) == sorted(genomes.gene_trees)
    for row in run.family_multipliers.values():
        assert tuple(row) == SEQUENCE_TARGETS
        assert isinstance(row["substitution"], float) and row["substitution"] > 0
    drawn = {row["substitution"] for row in run.family_multipliers.values()}
    assert len(drawn) == len(run.family_multipliers)          # one independent draw apiece


def test_the_multiplier_scales_every_branch_of_that_family_alone(genomes):
    """The phylogram is the tree the sequences were drawn along, so it is where the multiplier has to
    show: each family's branches are the strict-clock run's, times the one number the table holds."""
    strict = sequences.simulate_sequences(genomes, model=jc69(), length=40, substitution=0.05, seed=3)
    varied = sequences.simulate_sequences(genomes, model=jc69(), length=40, seed=3,
                                          substitution=_varying())
    for fam, row in varied.family_multipliers.items():
        one = np.asarray(_lengths(strict.phylograms[fam]["complete"]))
        many = np.asarray(_lengths(varied.phylograms[fam]["complete"]))
        # a Newick branch length is written to about seven significant digits, so this tolerance is
        # the file format's; the arithmetic underneath is one multiplication
        assert np.allclose(many, one * row["substitution"], rtol=1e-5)


def test_the_clock_and_the_family_draw_compose(genomes):
    """SPEC §5: modifiers multiply. One says which lineages were dealt a fast tempo, the other which
    families were — and a run carrying both is the clock-only run scaled family by family."""
    clock = PerSite(0.05).varying_among("lineages", LogNormal(0.0, 0.4))
    both = PerSite(0.05).varying_among("lineages", LogNormal(0.0, 0.4)) \
                        .varying_among("families", LogNormal(0.0, 0.8))
    a = sequences.simulate_sequences(genomes, model=jc69(), length=40, substitution=clock, seed=4)
    b = sequences.simulate_sequences(genomes, model=jc69(), length=40, substitution=both, seed=4)
    assert b.family_multipliers and a.family_multipliers == {}
    for fam, row in b.family_multipliers.items():
        one = np.asarray(_lengths(a.phylograms[fam]["complete"]))
        many = np.asarray(_lengths(b.phylograms[fam]["complete"]))
        assert np.allclose(many, one * row["substitution"], rtol=1e-5)


def test_the_species_phylogram_stays_on_the_run_s_own_rate(genomes):
    """The species phylogram is the clock made visible, and the clock is what every family shares.
    A per-family multiplier belongs to one family, so it must not reach that tree."""
    strict = sequences.simulate_sequences(genomes, model=jc69(), length=40, substitution=0.05, seed=3)
    varied = sequences.simulate_sequences(genomes, model=jc69(), length=40, seed=3,
                                          substitution=_varying())
    assert varied.species_phylogram == strict.species_phylogram


def test_a_family_with_a_larger_multiplier_diverges_further(genomes):
    """The multiplier scales a rate, so it shows in the alignments too: a family drawn fast is one whose
    tips are less alike."""
    run = sequences.simulate_sequences(genomes, model=jc69(), length=400, seed=5,
                                       substitution=_varying(base=0.1, sigma=1.0))
    families = [f for f, aln in run.alignments.items() if len(aln) > 1]
    mult = [run.family_multipliers[f]["substitution"] for f in families]
    identity = [sequences.mean_pairwise_identity({f: run.alignments[f]}) for f in families]
    assert np.corrcoef(_ranks(mult), _ranks(identity))[0, 1] < -0.5


def test_a_restricted_run_draws_for_the_families_it_evolves():
    tree = simulate_species_tree(birth=1.0, death=0.2, n_extant=8, seed=6).complete_tree
    g = simulate_genomes_family(tree, duplication=0.1, transfer=0.05, loss=0.1, initial_families=6,
                                seed=7, families=[family("A"), family("B")])
    run = sequences.simulate_sequences(g, families=["A"], model=jc69(), length=40, seed=8,
                                       substitution=_varying())
    assert sorted(run.family_multipliers) == [g.family_names["A"]]


def test_the_table_reads_back(genomes, tmp_path):
    run = sequences.simulate_sequences(genomes, model=jc69(), length=40, seed=3,
                                       substitution=_varying())
    run.write(tmp_path)
    text = (tmp_path / "family_multipliers.tsv").read_text(encoding="utf-8")
    assert text.startswith(HEADER)
    assert multipliers_from_tsv(text) == run.family_multipliers


def test_family_multipliers_is_an_output_that_can_be_left_out(genomes, tmp_path):
    run = sequences.simulate_sequences(genomes, model=jc69(), length=40, seed=3,
                                       substitution=_varying())
    run.write(tmp_path, outputs=("phylograms",))
    assert not (tmp_path / "family_multipliers.tsv").exists()


def test_a_streamed_run_writes_the_table_the_run_in_memory_writes(genomes, tmp_path):
    kw = dict(model=jc69(), length=40, seed=3, substitution=_varying())
    sequences.simulate_sequences(genomes, **kw).write(tmp_path / "written")
    sequences.simulate_sequences(genomes, stream_to=tmp_path / "streamed", **kw)
    assert ((tmp_path / "streamed" / "family_multipliers.tsv").read_bytes()
            == (tmp_path / "written" / "family_multipliers.tsv").read_bytes())


def test_the_parallel_engine_gives_the_same_table_for_any_worker_count(genomes):
    kw = dict(model=jc69(), length=40, seed=3, substitution=_varying())
    a = sequences.simulate_sequences(genomes, parallel=2, **kw)
    b = sequences.simulate_sequences(genomes, parallel=3, **kw)
    assert a.family_multipliers == b.family_multipliers
    assert a.phylograms == b.phylograms


def test_a_nucleotide_run_refuses_a_draw_among_families():
    tree = simulate_species_tree(birth=1.0, death=0.2, n_extant=6, seed=9).complete_tree
    g = simulate_genomes_nucleotide(tree, duplication=0.1, duplication_extent=50.0, seed=10,
                                    root_length=900, genes=3, gene_length=100)
    with pytest.raises(ValueError, match="units are blocks"):
        sequences.simulate_sequences(g, model=jc69(), seed=11, substitution=_varying())


def test_an_inherited_value_among_families_is_still_refused(genomes):
    """A family is born by origination and has no parent family to inherit from, so there is no
    per-family Drift to carry — the level takes the draw and refuses the drift."""
    with pytest.raises(ValueError, match="does not read"):
        sequences.simulate_sequences(
            genomes, model=jc69(), length=40, seed=3,
            substitution=PerSite(0.05).varying_among("families", Drift(LogNormal(0.0, 0.5))))


def test_the_multiplier_reaches_the_indel_rates_too(genomes, monkeypatch):
    """Indel rates are relative to substitution, so a family's own multiplier reaches them through the
    same base: a family drawn fast substitutes fast and gains and loses sites fast. One speed for the
    family, not two. Read off the base each family's indel history was drawn at, because the counts
    themselves are Poisson and one family's draw says nothing."""
    seen = []
    real = sequences.draw_indel_history

    def record(*args, **kw):
        seen.append(kw["rate_base"])
        return real(*args, **kw)

    monkeypatch.setattr(sequences, "draw_indel_history", record)
    run = sequences.simulate_sequences(genomes, model=jc69(), length=100, seed=12,
                                       insertion=0.05, deletion=0.05, substitution=_varying(base=0.1))
    assert seen == pytest.approx([0.1 * run.family_multipliers[f]["substitution"]
                                  for f in sorted(run.family_multipliers)])
