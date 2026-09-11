"""``links.tsv``: the links a family genome run read from its own gene content (issue #426).

Every run writes the table; a run that read no gene content writes its header alone. One row per
modifier on ``"genomes:…"``, with the family whose parameter it is (``run`` for the run's own), the
parameter, the driver, the modifier and the mapping as text.
"""

import pytest

from zombi2.genomes import family, read_run, simulate_genomes_family
from zombi2.genomes.family import FamilyGenomesResult
from zombi2.genomes.links import COLUMNS, Link, links_from_tsv
from zombi2.params import Between, PerCopy, Recipients, Table
from zombi2.species import simulate_species_tree

HEADER = "\t".join(COLUMNS) + "\n"


@pytest.fixture(scope="module")
def tree():
    return simulate_species_tree(birth=1.0, death=0.2, n_extant=15, seed=1).complete_tree


def _run(tree, **kw):
    base = dict(duplication=0.1, loss=0.3, transfer=0.2, initial_families=3, seed=2)
    return simulate_genomes_family(tree, **{**base, **kw})


def step(f):
    return 1.0 if f == 1.0 else 20.0


def by_size(n):
    return 1.0 + n / 10


def test_a_run_without_links_writes_the_header_alone(tree, tmp_path):
    g = _run(tree)
    assert g.links == ()
    g.write(tmp_path)
    assert (tmp_path / "links.tsv").read_text(encoding="utf-8") == HEADER


def test_each_kind_of_link_is_one_row(tree):
    g = _run(tree, joint=True,
             transfer_to=Recipients().weighted_by("genomes:count", by_size),
             families=[
                 family("A"),
                 family("B", loss=PerCopy(0.3).scaled_by("genomes:A", {"present": 1.0, "absent": 20.0}),
                        transfer_to=Recipients().weighted_by(
                            "genomes:A", Between({("present", "present"): 20.0}))),
                 family("C", module="pathway", loss=PerCopy().set_by("genomes:module:pathway", step)),
                 family("D", module="pathway")])
    assert g.links == (
        Link("run", "transfer_to", "genomes:count", "weighted_by", "Curve(by_size)"),
        Link("B", "loss", "genomes:A", "scaled_by", "present=1.0;absent=20.0"),
        Link("B", "transfer_to", "genomes:A", "weighted_by", "present>present=20.0"),
        Link("C", "loss", "genomes:module:pathway", "set_by", "0/2=20.0;1/2=20.0;2/2=1.0"),
    )


def test_a_run_rate_is_listed_under_run(tree):
    g = _run(tree, joint=True, families=[family("A")],
             loss=PerCopy(0.3).scaled_by("genomes:A", {"present": 0.5, "absent": 2.0}))
    assert g.links == (Link("run", "loss", "genomes:A", "scaled_by", "present=0.5;absent=2.0"),)


def test_two_modifiers_on_one_rate_are_two_rows_in_written_order(tree):
    g = _run(tree, joint=True, families=[
        family("A"), family("C"),
        family("B", loss=PerCopy(0.3).scaled_by("genomes:A", {"present": 1.0, "absent": 2.0})
                                     .scaled_by("genomes:C", {"present": 1.0, "absent": 3.0}))])
    assert [(x.family, x.driver) for x in g.links] == [("B", "genomes:A"), ("B", "genomes:C")]


def test_a_table_default_is_written(tree):
    g = _run(tree, joint=True, families=[
        family("A"), family("B", loss=PerCopy(0.3).scaled_by("genomes:A", Table({"present": 2.0}, default=0.5)))])
    assert g.links[0].mapping == "present=2.0;default=0.5"


def test_a_driver_read_from_a_trait_is_not_a_link(tree):
    from zombi2.traits import simulate_discrete

    habitat = simulate_discrete(tree, states=["soil", "water"], switch=0.3, seed=2)
    g = _run(tree, loss=PerCopy(0.3).scaled_by(habitat, {"soil": 1.0, "water": 2.0}),
             families=[family("B", transfer_to=Recipients().weighted_by(habitat, {"soil": 1.0, "water": 0.5}))])
    assert g.links == ()


def test_read_run_gives_the_links_back(tree, tmp_path):
    g = _run(tree, joint=True, families=[
        family("A"), family("B", loss=PerCopy(0.3).scaled_by("genomes:A", {"present": 1.0, "absent": 20.0}))])
    g.write(tmp_path)
    assert read_run(tmp_path).links == g.links


def test_links_is_an_output_that_can_be_left_out(tree, tmp_path):
    assert "links" in FamilyGenomesResult.OUTPUTS
    _run(tree).write(tmp_path, outputs=("profiles",))
    assert not (tmp_path / "links.tsv").exists()


def test_a_streamed_run_writes_the_header_alone(tree, tmp_path):
    simulate_genomes_family(tree, duplication=0.1, loss=0.3, initial_families=3, seed=2,
                            parallel=2, stream_to=tmp_path / "stream")
    assert (tmp_path / "stream" / "links.tsv").read_text(encoding="utf-8") == HEADER


def test_a_table_with_another_header_is_refused():
    with pytest.raises(ValueError, match="must start with the header"):
        links_from_tsv("family\ttarget\tdriver\tverb\tmapping\n")
