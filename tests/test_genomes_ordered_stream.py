"""A streamed ordered genome run (issue #436): ``stream_to=DIR`` writes the files as the run goes.

The run is the same run, so for one seed a streamed run's files hold what ``result.write(DIR)`` writes
for the run kept in memory, byte for byte. ``gene_order.tsv`` is the one file in a different order: a
stream writes a node's rows when its branch ends, so that file is compared with its rows sorted, and
its order is checked on its own.
"""

import pathlib

import pytest

from zombi2.genomes import StreamedRun, family, ordered, read_run, simulate_genomes_ordered
from zombi2.params import PerCopy, Recipients
from zombi2.species import simulate_species_tree


@pytest.fixture(scope="module")
def tree():
    # death > 0, so some branches die and the files name them e<id>
    return simulate_species_tree(birth=1.0, death=0.3, n_extant=25, seed=3).complete_tree


def _daughter(tree):
    return sorted(n for n in tree.nodes if tree.nodes[n].parent == tree.root)[0]


def _cases(tree):
    return {
        "plain": dict(duplication=0.2, loss=0.2, transfer=0.2, origination=1.0, initial_families=30),
        "segments, replacement and rearrangements": dict(
            duplication=0.3, loss=0.2, transfer=0.3, origination=1.0, initial_families=30,
            duplication_extent=3, loss_extent=2, transfer_extent=2, replacement=True,
            inversion=0.2, transposition=0.2),
        "chromosome events on linear chromosomes": dict(
            duplication=0.2, loss=0.2, origination=1.0, initial_families=30, chromosomes=3,
            topology="linear", fission=0.15, fusion=0.15, chromosome_origination=0.05,
            chromosome_loss=0.05, translocation=0.1),
        "joint, own rates, a driven transfer_to and a placed family": dict(
            duplication=0.2, loss=0.15, transfer=0.3, origination=1.0, initial_families=20, joint=True,
            transfer_to=Recipients().weighted_by("genomes:A", {"present": 5.0, "absent": 1.0}),
            families=[family("A", module="m"),
                      family("B", module="m", loss=PerCopy(0.4).scaled_by(
                          "genomes:module:m", lambda f: 0.2 if f == 1.0 else 1.0)),
                      family("C", origin=_daughter(tree), duplication=0.5)]),
        "a rate on a schedule": dict(duplication=PerCopy(0.2).changing_at({0.0: 1.0, 1.0: 3.0}),
                                     loss=0.2, origination=1.0, initial_families=20),
    }


def _files(directory):
    d = pathlib.Path(directory)
    return {p.relative_to(d).as_posix(): p.read_bytes() for p in sorted(d.rglob("*")) if p.is_file()}


def _assert_same_files(written, streamed):
    a, b = _files(written), _files(streamed)
    assert sorted(a) == sorted(b)
    for name in a:
        if name == "gene_order.tsv":
            head_a, *rows_a = a[name].decode().splitlines()
            head_b, *rows_b = b[name].decode().splitlines()
            assert head_a == head_b and sorted(rows_a) == sorted(rows_b), name
        else:
            assert a[name] == b[name], f"{name} differs"


CASE_NAMES = ["plain", "segments, replacement and rearrangements",
              "chromosome events on linear chromosomes",
              "joint, own rates, a driven transfer_to and a placed family", "a rate on a schedule"]


@pytest.mark.parametrize("case", CASE_NAMES)
def test_a_streamed_run_writes_what_the_run_in_memory_writes(tree, tmp_path, case):
    kw = _cases(tree)[case]
    kept = simulate_genomes_ordered(tree, seed=7, **kw)
    kept.write(tmp_path / "written")
    handle = simulate_genomes_ordered(tree, seed=7, stream_to=tmp_path / "streamed", **kw)
    _assert_same_files(tmp_path / "written", tmp_path / "streamed")
    assert isinstance(handle, StreamedRun)
    assert handle.n_events == len(kept.edges)
    assert handle.n_families == kept.summary()["families"]["born"]
    assert handle.outputs == ordered.OrderedGenomesResult.OUTPUTS
    assert pathlib.Path(handle.path("gene_order")).exists()


def test_gene_order_lists_the_nodes_as_their_branches_end(tree, tmp_path):
    simulate_genomes_ordered(tree, seed=2, stream_to=tmp_path, **_cases(tree)["plain"])
    names = tree.labels()
    rows = (tmp_path / "gene_order.tsv").read_text().splitlines()[1:]
    order = list(dict.fromkeys(row.split("\t", 1)[0] for row in rows))
    expected = [names[i] for _, i in sorted((tree.nodes[i].end_time, i) for i in tree.nodes)]
    assert order == [name for name in expected if name in set(order)]
    seen, previous = set(), None                      # each node's rows sit together
    for row in rows:
        name = row.split("\t", 1)[0]
        assert name == previous or name not in seen
        seen.add(name)
        previous = name


def test_outputs_picks_the_files(tree, tmp_path):
    simulate_genomes_ordered(tree, seed=1, stream_to=tmp_path, outputs=("events",),
                             **_cases(tree)["plain"])
    assert sorted(_files(tmp_path)) == ["genome_events.tsv", "rearrangement_events.tsv"]


def test_gene_trees_alone_leave_no_event_log_behind(tree, tmp_path):
    kw = _cases(tree)["segments, replacement and rearrangements"]
    kept = simulate_genomes_ordered(tree, seed=4, **kw)
    kept.write(tmp_path / "written", outputs=("gene_trees",))
    simulate_genomes_ordered(tree, seed=4, stream_to=tmp_path / "streamed", outputs=("gene_trees",),
                             **kw)
    _assert_same_files(tmp_path / "written", tmp_path / "streamed")


def test_gene_trees_built_in_many_groups_are_the_same_trees(tree, tmp_path, monkeypatch):
    """With a few rows per group, the log is split into many groups of families, each built alone."""
    monkeypatch.setattr(ordered, "_GENE_TREE_GROUP_ROWS", 40)
    kw = _cases(tree)["joint, own rates, a driven transfer_to and a placed family"]
    kept = simulate_genomes_ordered(tree, seed=5, **kw)
    kept.write(tmp_path / "written", outputs=("gene_trees", "events"))
    simulate_genomes_ordered(tree, seed=5, stream_to=tmp_path / "streamed",
                             outputs=("gene_trees", "events"), **kw)
    _assert_same_files(tmp_path / "written", tmp_path / "streamed")


def test_a_second_run_into_the_same_directory_leaves_nothing_of_the_first(tree, tmp_path):
    simulate_genomes_ordered(tree, seed=1, stream_to=tmp_path, initial_families=40, duplication=0.3,
                             origination=2.0)
    first = {p.name for p in (tmp_path / "gene_trees").iterdir()}
    simulate_genomes_ordered(tree, seed=1, stream_to=tmp_path, initial_families=3, origination=0.0)
    second = {p.name for p in (tmp_path / "gene_trees").iterdir()}
    assert len(second) < len(first) and not second - first


def test_a_streamed_run_reads_back(tree, tmp_path):
    kw = _cases(tree)["plain"]
    kept = simulate_genomes_ordered(tree, seed=3, **kw)
    back = read_run(simulate_genomes_ordered(tree, seed=3, stream_to=tmp_path, **kw))
    assert len(back.edges) == len(kept.edges)
    assert sorted(back.gene_trees) == sorted(kept.gene_trees)


def test_an_ordered_run_written_in_memory_reads_back(tree, tmp_path):
    """`read_run` read ``initial_genome.tsv`` as two columns, which only a family run writes, so it
    failed on every ordered run directory, streamed or not."""
    kept = simulate_genomes_ordered(tree, seed=3, **_cases(tree)["plain"])
    kept.write(tmp_path)
    back = read_run(tmp_path)
    assert len(back.edges) == len(kept.edges)
    assert [c.id for c in back.initial_genome] == [g.id for chrom in kept.initial_genome
                                                    for g in chrom.genes]


def test_outputs_without_stream_to_is_refused(tree):
    with pytest.raises(ValueError, match="outputs applies to a streamed run"):
        simulate_genomes_ordered(tree, seed=1, outputs=("events",))


def test_an_unknown_output_is_refused(tree, tmp_path):
    with pytest.raises(ValueError, match="unknown stream outputs"):
        simulate_genomes_ordered(tree, seed=1, stream_to=tmp_path, outputs=("genomes",))
