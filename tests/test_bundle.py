"""Bundled output (issue #463): ``bundle=True`` / ``--bundle`` writes each per-family output as one file.

The bundled layout holds exactly what the per-family layout holds. So for every writer, bundling a
run written without ``bundle`` gives the files that writer gives with it, byte for byte, and
unbundling gives the per-family files back. A streamed run bundles too, and gives what the run in
memory gives.
"""

import pathlib

import pytest

from zombi2.cli.main import main
from zombi2.genomes import (ordered, simulate_genomes_family, simulate_genomes_nucleotide,
                            simulate_genomes_ordered)
from zombi2.sequences import jc69, simulate_sequences
from zombi2.species import simulate_species_tree
from zombi2.tools.bundle import bundle_run, unbundle_run

SEQUENCE_OUTPUTS = ("alignments", "ancestral", "phylograms")


@pytest.fixture(scope="module")
def tree():
    # death > 0, so some families die out and write a complete tree with no extant one
    return simulate_species_tree(birth=1.0, death=0.3, n_extant=12, seed=3).complete_tree


@pytest.fixture(scope="module")
def genomes(tree):
    return simulate_genomes_family(tree, duplication=0.2, transfer=0.1, loss=0.6, origination=0.5,
                                   initial_families=10, seed=4)


def _files(directory):
    d = pathlib.Path(directory)
    return {p.relative_to(d).as_posix(): p.read_bytes() for p in sorted(d.rglob("*")) if p.is_file()}


def _assert_converts(per_family, bundled, tmp_path):
    """Bundling ``per_family`` gives ``bundled``, and unbundling that gives ``per_family`` back."""
    work = tmp_path / "work"
    work.mkdir()
    for name, data in _files(per_family).items():
        (work / name).parent.mkdir(parents=True, exist_ok=True)
        (work / name).write_bytes(data)
    assert bundle_run(work, remove=True)
    assert _files(work) == _files(bundled)
    assert unbundle_run(work, remove=True)
    assert _files(work) == _files(per_family)


def test_family_gene_trees_bundle_and_unbundle_exactly(genomes, tmp_path):
    genomes.write(tmp_path / "a")
    genomes.write(tmp_path / "b", bundle=True)
    assert "gene_trees.tsv" in _files(tmp_path / "b")
    assert not (tmp_path / "b" / "gene_trees").exists()
    # a family whose copies all died has a complete tree and an empty extant cell
    rows = (tmp_path / "b" / "gene_trees.tsv").read_text().splitlines()
    assert rows[0] == "family\tcomplete\textant"
    assert any(row.endswith("\t") for row in rows[1:])
    _assert_converts(tmp_path / "a", tmp_path / "b", tmp_path)


def test_ordered_and_nucleotide_gene_trees_bundle_exactly(tree, tmp_path):
    runs = {"ordered": simulate_genomes_ordered(tree, duplication=0.2, loss=0.2, transfer=0.2,
                                                origination=1.0, initial_families=15, seed=5),
            "nucleotide": simulate_genomes_nucleotide(tree, duplication=0.5, duplication_extent=40,
                                                      loss=0.5, loss_extent=40, root_length=600,
                                                      genes=3, gene_length=90, seed=6)}
    for name, run in runs.items():
        run.write(tmp_path / name / "a")
        run.write(tmp_path / name / "b", bundle=True)
        _assert_converts(tmp_path / name / "a", tmp_path / name / "b", tmp_path / name)


def test_sequences_bundle_and_unbundle_exactly(genomes, tmp_path):
    s = simulate_sequences(genomes, model=jc69(), length=40, seed=5)
    s.write(tmp_path / "a", outputs=SEQUENCE_OUTPUTS)
    s.write(tmp_path / "b", outputs=SEQUENCE_OUTPUTS, bundle=True)
    assert sorted(_files(tmp_path / "b")) == ["alignments.fasta", "ancestral.fasta",
                                               "phylograms.tsv"]
    # a header keeps the tip label as the sequence name, so it still matches the gene tree
    header = (tmp_path / "b" / "alignments.fasta").read_text().splitlines()[0]
    name, family = header[1:].split(" ")
    assert name.startswith("n") and "_g" in name and family.startswith("fam")
    _assert_converts(tmp_path / "a", tmp_path / "b", tmp_path)


def test_nucleotide_sequences_bundle_their_blocks(tree, tmp_path):
    g = simulate_genomes_nucleotide(tree, duplication=0.5, duplication_extent=40, loss=0.5,
                                    loss_extent=40, root_length=600, genes=3, gene_length=90, seed=6)
    s = simulate_sequences(g, model=jc69(), substitution=0.3, seed=7)
    s.write(tmp_path / "a", outputs=SEQUENCE_OUTPUTS)
    s.write(tmp_path / "b", outputs=SEQUENCE_OUTPUTS, bundle=True)
    assert (tmp_path / "b" / "phylograms.tsv").read_text().startswith("block\tcomplete\textant\n")
    _assert_converts(tmp_path / "a", tmp_path / "b", tmp_path)


def test_a_streamed_sequence_run_bundles_what_the_run_in_memory_bundles(genomes, tmp_path):
    kw = dict(model=jc69(), length=40, seed=5)
    simulate_sequences(genomes, **kw).write(tmp_path / "memory", outputs=SEQUENCE_OUTPUTS,
                                            bundle=True)
    simulate_sequences(genomes, **kw, stream_to=tmp_path / "streamed", outputs=SEQUENCE_OUTPUTS,
                       bundle=True)
    assert _files(tmp_path / "streamed") == _files(tmp_path / "memory")


@pytest.mark.parametrize("group_rows", [None, 40])
def test_a_streamed_ordered_run_bundles_what_the_run_in_memory_bundles(tree, tmp_path, monkeypatch,
                                                                     group_rows):
    # with 40 rows a group, the trees are built in many groups and their rows merged back in order
    if group_rows is not None:
        monkeypatch.setattr(ordered, "_GENE_TREE_GROUP_ROWS", group_rows)
        monkeypatch.setattr(ordered, "_GENE_TREE_OPEN_FILES", 3)
    kw = dict(duplication=0.2, loss=0.2, transfer=0.2, origination=1.0, initial_families=30, seed=5)
    simulate_genomes_ordered(tree, **kw).write(tmp_path / "memory", outputs=("gene_trees",),
                                               bundle=True)
    simulate_genomes_ordered(tree, **kw, stream_to=tmp_path / "streamed", outputs=("gene_trees",),
                             bundle=True)
    assert _files(tmp_path / "streamed") == _files(tmp_path / "memory")


@pytest.mark.parametrize("parallel", [False, 2])
def test_a_streamed_family_run_bundles_what_it_writes_unbundled(tree, tmp_path, parallel):
    kw = dict(duplication=0.2, transfer=0.1, loss=0.3, origination=0.5, initial_families=300,
              seed=4, parallel=parallel, outputs=("gene_trees",))
    simulate_genomes_family(tree, **kw, stream_to=tmp_path / "a")
    simulate_genomes_family(tree, **kw, stream_to=tmp_path / "b", bundle=True)
    assert sorted(_files(tmp_path / "b")) == ["gene_trees.tsv"]
    _assert_converts(tmp_path / "a", tmp_path / "b", tmp_path)


def test_bundle_without_stream_to_is_refused(tree, genomes):
    kw = dict(duplication=0.2, loss=0.2, origination=0.5, initial_families=5, seed=1, bundle=True)
    for simulate in (simulate_genomes_family, simulate_genomes_ordered):
        with pytest.raises(ValueError, match="bundle applies to a streamed run"):
            simulate(tree, **kw)
    with pytest.raises(ValueError, match="bundle applies to a streamed run"):
        simulate_sequences(genomes, model=jc69(), length=20, seed=1, bundle=True)


def test_a_write_in_one_layout_clears_the_other(genomes, tmp_path):
    # a run's directory describes that run: the previous run's layout must not survive beside it
    genomes.write(tmp_path, outputs=("gene_trees",))
    genomes.write(tmp_path, outputs=("gene_trees",), bundle=True)
    assert sorted(_files(tmp_path)) == ["gene_trees.tsv"]
    genomes.write(tmp_path, outputs=("gene_trees",))
    assert "gene_trees.tsv" not in _files(tmp_path)


def test_two_layouts_that_agree_are_accepted_and_two_that_differ_are_refused(genomes, tmp_path):
    genomes.write(tmp_path, outputs=("gene_trees",))
    bundle_run(tmp_path)                            # both layouts now, and they agree
    assert bundle_run(tmp_path, remove=True)
    assert sorted(_files(tmp_path)) == ["gene_trees.tsv"]
    unbundle_run(tmp_path)
    first = sorted((tmp_path / "gene_trees").iterdir())[0]
    first.write_text("(a,b);\n")
    with pytest.raises(ValueError, match="they differ"):
        bundle_run(tmp_path)
    with pytest.raises(ValueError, match="they differ"):
        unbundle_run(tmp_path)


def test_a_flat_run_converts_in_place(genomes, tmp_path):
    genomes.write(tmp_path / "a", outputs=("gene_trees",), flat=True)
    genomes.write(tmp_path / "b", outputs=("gene_trees",), bundle=True)
    bundle_run(tmp_path / "a", remove=True)
    assert _files(tmp_path / "a") == _files(tmp_path / "b")
    unbundle_run(tmp_path / "a", flat=True, remove=True)
    genomes.write(tmp_path / "c", outputs=("gene_trees",), flat=True)
    assert _files(tmp_path / "a") == _files(tmp_path / "c")


def test_the_cli_writes_and_converts_a_bundled_run(tmp_path, capsys):
    run = tmp_path / "run"
    main(["species", str(run), "--birth", "1.0", "--death", "0.3", "--n-extant", "10", "--seed", "1",
          "--quiet"])
    main(["genomes", str(run), "--duplication", "0.2", "--loss", "0.3", "--origination", "0.5",
          "--seed", "2", "--bundle", "--quiet"])
    main(["sequences", str(run), "--length", "30", "--seed", "3", "--bundle", "--quiet",
          "--write", "alignments", "phylograms"])
    assert (run / "genomes" / "gene_trees.tsv").is_file()
    assert not (run / "genomes" / "gene_trees").exists()
    assert sorted(p.name for p in (run / "sequences").iterdir()
                  if p.suffix in (".fasta", ".tsv")) == ["alignments.fasta", "phylograms.tsv"]
    before = _files(run)
    capsys.readouterr()
    assert main(["tools", "unbundle", str(run), "--remove"]) == 0
    assert "gene_trees" in capsys.readouterr().out
    assert (run / "genomes" / "gene_trees").is_dir()
    assert main(["tools", "bundle", str(run), "--remove"]) == 0
    assert _files(run) == before
