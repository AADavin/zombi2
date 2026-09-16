"""``family_multipliers.tsv``: each family's drawn rate multipliers (issue #441).

A rate written with ``varying_among("families", law)`` draws a number for each family when the family is
born. The run keeps it on the result, writes it with its other files, and ``read_run`` reads it back.
A run whose rates do not vary among families writes the header alone.
"""

import numpy as np
import pytest

from zombi2.genomes import family, read_run, simulate_genomes_family, simulate_genomes_ordered
from zombi2.genomes.family import FamilyGenomesResult
from zombi2.genomes.multipliers import FAMILY_TARGETS, ORDERED_TARGETS, multipliers_from_tsv
from zombi2.genomes.ordered import OrderedGenomesResult
from zombi2.params import Gamma, LogNormal, PerCopy, Random
from zombi2.species import simulate_species_tree


@pytest.fixture(scope="module")
def tree():
    return simulate_species_tree(birth=1.0, death=0.3, n_extant=15, seed=1).complete_tree


def _varying():
    speed = Random("families", LogNormal(0.0, 0.5))    # one Random on duplication and loss: one draw
    return dict(duplication=PerCopy(0.2).varying_among(speed),
                transfer=PerCopy(0.3).varying_among("families", Gamma(shape=2.0, scale=0.5)),
                loss=PerCopy(0.25).varying_among(speed))


def _header(targets):
    return "\t".join(("family", *targets)) + "\n"


def _ranks(values):
    """Ranks from 0, ties sharing their mean rank."""
    values = np.asarray(values, dtype=float)
    ranks = np.empty(len(values))
    ranks[values.argsort(kind="stable")] = np.arange(len(values))
    for v in np.unique(values):
        at = values == v
        ranks[at] = ranks[at].mean()
    return ranks


def _extant_copies(run):
    """Copies of each family over the extant tips, from ``profiles.tsv``."""
    header, *rows = run.profiles.to_tsv().splitlines()
    return {int(row.split("\t")[0]): sum(map(int, row.split("\t")[1:])) for row in rows}


@pytest.mark.parametrize("simulate, targets", [(simulate_genomes_family, FAMILY_TARGETS),
                                               (simulate_genomes_ordered, ORDERED_TARGETS)],
                         ids=["family", "ordered"])
def test_a_run_whose_rates_do_not_vary_writes_the_header_alone(tree, tmp_path, simulate, targets):
    run = simulate(tree, duplication=0.2, transfer=0.3, loss=0.25, initial_families=5, seed=2)
    assert run.family_multipliers == {}
    run.write(tmp_path)
    assert (tmp_path / "family_multipliers.tsv").read_text(encoding="utf-8") == _header(targets)


def test_every_family_has_a_row(tree):
    run = simulate_genomes_family(tree, **_varying(), origination=2.0, initial_families=20, seed=3)
    born = {e.family for e in run.edges if e.kind == "origination"}
    assert sorted(run.family_multipliers) == sorted(born)
    assert len(born) == run.summary()["families"]["born"]
    for row in run.family_multipliers.values():
        assert tuple(row) == FAMILY_TARGETS
        assert all(isinstance(v, float) and v > 0 for v in row.values())
        assert row["duplication"] == row["loss"]              # one law object, one draw
    assert any(row["transfer"] != row["loss"] for row in run.family_multipliers.values())


def test_a_rate_that_does_not_vary_gives_every_family_one(tree):
    run = simulate_genomes_family(tree, duplication=PerCopy(0.2).varying_among("families", LogNormal(0.0, 0.5)),
                                  transfer=0.3, loss=0.25, initial_families=20, seed=3)
    assert {row["loss"] for row in run.family_multipliers.values()} == {1.0}
    assert {row["transfer"] for row in run.family_multipliers.values()} == {1.0}
    assert len({row["duplication"] for row in run.family_multipliers.values()}) > 1


@pytest.mark.parametrize("simulate", [simulate_genomes_family, simulate_genomes_ordered],
                         ids=["family", "ordered"])
def test_a_family_with_its_own_rate_has_an_empty_cell(tree, tmp_path, simulate):
    run = simulate(tree, **_varying(), initial_families=10, seed=4,
                   families=[family("A", loss=PerCopy(0.5)), family("B")])
    a, b = run.family_names["A"], run.family_names["B"]
    assert run.family_multipliers[a]["loss"] is None
    assert isinstance(run.family_multipliers[a]["duplication"], float)
    assert None not in run.family_multipliers[b].values()
    run.write(tmp_path)
    rows = {line.split("\t")[0]: line.split("\t")
            for line in (tmp_path / "family_multipliers.tsv").read_text(encoding="utf-8").splitlines()}
    assert rows[str(a)][1 + FAMILY_TARGETS.index("loss")] == ""


@pytest.mark.parametrize("engine", ["serial", "parallel", "ordered"])
def test_a_family_with_a_larger_loss_multiplier_keeps_fewer_copies(tree, engine):
    """The table holds the numbers the run used: with loss the only event, a family that drew a larger
    loss multiplier leaves fewer copies at the tips."""
    kw = dict(duplication=0.0, transfer=0.0, origination=0.0, initial_families=300, seed=5,
              loss=PerCopy(0.3).varying_among("families", LogNormal(0.0, 1.0)))
    if engine == "ordered":
        run = simulate_genomes_ordered(tree, **kw)
    else:
        run = simulate_genomes_family(tree, **kw, parallel=2 if engine == "parallel" else False)
    copies = _extant_copies(run)
    families = sorted(run.family_multipliers)
    loss = [run.family_multipliers[f]["loss"] for f in families]
    kept = [copies.get(f, 0) for f in families]
    assert np.corrcoef(_ranks(loss), _ranks(kept))[0, 1] < -0.3


def test_the_parallel_engine_gives_the_same_table_for_any_worker_count(tree):
    kw = dict(**_varying(), origination=2.0, initial_families=20, seed=6)
    assert (simulate_genomes_family(tree, **kw, parallel=2).family_multipliers
            == simulate_genomes_family(tree, **kw, parallel=3).family_multipliers)


def test_a_streamed_family_run_writes_the_table_the_run_in_memory_writes(tree, tmp_path):
    kw = dict(**_varying(), origination=2.0, initial_families=20, seed=7, parallel=2)
    simulate_genomes_family(tree, **kw).write(tmp_path / "written")
    handle = simulate_genomes_family(tree, **kw, stream_to=tmp_path / "streamed")
    written = (tmp_path / "written" / "family_multipliers.tsv").read_bytes()
    assert (tmp_path / "streamed" / "family_multipliers.tsv").read_bytes() == written
    assert handle.path("family_multipliers").endswith("family_multipliers.tsv")
    assert written.count(b"\n") == 1 + handle.n_families


def test_a_streamed_family_run_whose_rates_do_not_vary_writes_the_header_alone(tree, tmp_path):
    simulate_genomes_family(tree, duplication=0.2, loss=0.25, initial_families=5, seed=2, parallel=2,
                            stream_to=tmp_path)
    assert (tmp_path / "family_multipliers.tsv").read_text(encoding="utf-8") == _header(FAMILY_TARGETS)


def test_an_ordered_run_has_a_column_for_each_rearrangement(tree, tmp_path):
    kw = dict(duplication=0.2, transfer=0.3, loss=0.25, origination=2.0, initial_families=20, seed=8,
              inversion=PerCopy(0.1).varying_among("families", Gamma(shape=2.0, scale=0.5)))
    run = simulate_genomes_ordered(tree, **kw)
    assert all(tuple(row) == ORDERED_TARGETS for row in run.family_multipliers.values())
    assert len({row["inversion"] for row in run.family_multipliers.values()}) > 1
    assert {row["transposition"] for row in run.family_multipliers.values()} == {1.0}
    run.write(tmp_path / "written")
    simulate_genomes_ordered(tree, **kw, stream_to=tmp_path / "streamed")
    written = (tmp_path / "written" / "family_multipliers.tsv").read_bytes()
    # read as text for the header: Windows writes each line ending as \r\n
    assert (tmp_path / "written" / "family_multipliers.tsv").read_text(encoding="utf-8").startswith(
        _header(ORDERED_TARGETS))
    assert (tmp_path / "streamed" / "family_multipliers.tsv").read_bytes() == written


def test_read_run_gives_the_multipliers_back(tree, tmp_path):
    run = simulate_genomes_family(tree, **_varying(), initial_families=10, seed=9,
                                  families=[family("A", loss=PerCopy(0.5))])
    run.write(tmp_path)
    assert read_run(tmp_path).family_multipliers == run.family_multipliers


def test_family_multipliers_is_an_output_that_can_be_left_out(tree, tmp_path):
    assert "family_multipliers" in FamilyGenomesResult.OUTPUTS
    assert "family_multipliers" in OrderedGenomesResult.OUTPUTS
    run = simulate_genomes_family(tree, **_varying(), initial_families=5, seed=2)
    run.write(tmp_path, outputs=("events",))
    assert not (tmp_path / "family_multipliers.tsv").exists()


def test_a_table_with_another_header_is_refused():
    # the sequences level writes this table too, with its own column, so the message lists every
    # header the reader knows (issue #443)
    with pytest.raises(ValueError, match="must start with one of the headers"):
        multipliers_from_tsv("family\tloss\n0\t1.0\n")
