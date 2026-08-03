"""The second panel's findings, each pinned by the test that would have caught it.

Five newcomers — someone installing from a package index, someone who refuses to read a manual, a
lecturer building a class practical, a returning ZOMBI v1 user, and someone handed a folder of output
by a departed postdoc — tested v0.28.0. Where the first panel found silent paths in the engine, this
one found them in the **bookkeeping**: a directory that describes a run that no longer exists, a
summary that counts files it did not write, a default nobody was told about.
"""

from __future__ import annotations

import json
import re

import pytest

from zombi2 import genomes, sequences, species
from zombi2.cli.main import main
from zombi2.genomes.events import event_counts
from zombi2.sequences.substitution_models import jc69


def _tree(**kw):
    return species.simulate_species_tree(birth=1.0, death=0.3, n_extant=10, seed=1, **kw)


# --- 1. a run's directory describes that run ------------------------------------------------------

def test_a_second_run_does_not_leave_the_first_run_behind(tmp_path):
    """The defect two reviewers hit from two directions. A lecturer re-ran the genome level at a
    higher loss rate, got zero surviving families — announced clearly — and ``tools treedist`` still
    answered ``rf 0`` from a leftover gene tree. A student comparing transfer rates in one directory
    counted the leftovers and concluded gene trees were being written for families that had died out.

    Both saw a plausible number rather than an error, which is what makes it worth a test."""
    run = tmp_path / "run"
    main(["species", str(run), "--birth", "1", "--death", "0.3", "--n-extant", "20",
          "--seed", "7", "--quiet"])
    main(["genomes", str(run), "--duplication", "0.2", "--transfer", "1.0", "--loss", "0.2",
          "--origination", "0.5", "--initial-families", "60", "--seed", "3", "--quiet"])
    trees = run / "genomes" / "gene_trees"
    assert len(list(trees.glob("*_extant.nwk"))) > 10, "the first run should leave real gene trees"

    # the same directory, a run so lossy that nothing survives
    main(["genomes", str(run), "--duplication", "0.05", "--transfer", "0.0", "--loss", "3.0",
          "--origination", "0.05", "--initial-families", "5", "--seed", "4", "--quiet", "--force"])
    profiles = (run / "genomes" / "profiles.tsv").read_text(encoding="utf-8").splitlines()
    assert len(profiles) - 1 == 0, "this run is meant to have no surviving families"
    assert list(trees.glob("*_extant.nwk")) == [], (
        "gene trees from the previous run are still on disk; anything reading this directory will "
        "answer from a run that no longer exists")


def test_the_directory_matches_the_run_family_by_family(tmp_path):
    """The general property, not just the empty case: what is on disk is this run's families and no
    others."""
    tree = _tree()
    big = genomes.simulate_genomes_family(tree, duplication=0.3, loss=0.1, initial_families=40, seed=1)
    big.write(tmp_path / "out")
    small = genomes.simulate_genomes_family(tree, duplication=0.1, loss=0.4, initial_families=5, seed=2)
    small.write(tmp_path / "out")

    on_disk = {int(p.stem.split("fam")[1].split("_")[0])
               for p in (tmp_path / "out" / "gene_trees").glob("*_complete.nwk")}
    assert on_disk == set(small.gene_trees), "the directory holds families the second run never had"


def test_flat_is_left_alone(tmp_path):
    """Under ``--flat`` every output and every level share one directory, so nothing in it can safely
    be called this output's to remove. The guard must not reach there."""
    tree = _tree()
    g = genomes.simulate_genomes_family(tree, duplication=0.2, loss=0.2, initial_families=6, seed=1)
    g.write(tmp_path / "flat", flat=True)
    (tmp_path / "flat" / "notes.txt").write_text("mine", encoding="utf-8")
    g.write(tmp_path / "flat", flat=True)
    assert (tmp_path / "flat" / "notes.txt").exists(), "a flat write must not clear the directory"


def test_a_streamed_run_does_not_clear_itself_family_by_family(tmp_path):
    """A streamed run resolves its output directory once per family, so a clear on every resolution
    would leave only the last family. It clears once, when the sink opens."""
    tree = _tree()
    g = genomes.simulate_genomes_family(tree, duplication=0.2, loss=0.2, initial_families=8, seed=1)
    streamed = sequences.simulate_sequences(g, model=jc69(), length=20, seed=1,
                                            stream_to=str(tmp_path / "s"))
    files = list((tmp_path / "s" / "alignments").glob("*.fasta"))
    assert len(files) == streamed.n_families > 1


# --- 2. genome_summary.json at every resolution ---------------------------------------------------

@pytest.mark.parametrize("resolution", ["family", "ordered", "nucleotide"])
def test_every_resolution_writes_a_genome_summary(tmp_path, resolution):
    """The migration guide names the replacing-transfer loss undercount as the change most likely to
    hand a returning user a plausible wrong number, and points them at ``genome_summary.json``. That
    file was written only at the family resolution — so the remedy was missing at the two resolutions
    where the gap is larger (64% at ordered, measured)."""
    run = tmp_path / resolution
    main(["species", str(run), "--birth", "1", "--death", "0.3", "--n-extant", "8",
          "--seed", "1", "--quiet"])
    argv = ["genomes", str(run), "--resolution", resolution, "--seed", "2", "--quiet"]
    if resolution == "nucleotide":
        argv += ["--root-length", "2000", "--genes", "6", "--gene-length", "100"]
    else:
        argv += ["--initial-families", "20"]
    assert main(argv) == 0

    path = run / "genomes" / "genome_summary.json"
    assert path.exists(), f"no genome_summary.json at --resolution {resolution}"
    s = json.loads(path.read_text(encoding="utf-8"))
    assert s["level"] == "genomes" and s["resolution"] == resolution
    assert set(s["events"]) >= {"initial", "origination", "duplication", "transfer", "loss",
                                "speciation"}


def test_the_summary_corrects_the_loss_undercount_the_raw_log_has():
    """The number the guide tells a returning user to trust. A copy displaced by an arriving transfer
    has no ``loss`` row of its own — it is the second parent of the ``transfer_replacing`` row — so
    the log's loss rows are short by exactly the number of replacing transfers."""
    tree = _tree()
    g = genomes.simulate_genomes_ordered(tree, duplication=0.3, transfer=0.6, loss=0.2,
                                         replacement=True, initial_families=30, chromosomes=1, seed=3)
    raw = sum(1 for e in g.events if e.kind == "loss")
    replacing = sum(1 for e in g.events if e.kind == "transfer_replacing")
    assert replacing > 0, "this configuration is meant to exercise replacement"
    assert g.summary()["events"]["loss"] == raw + replacing


def test_the_three_resolutions_count_events_the_same_way():
    """One implementation, so they cannot drift — which they would, being three engines that each
    have their own idea of what an event is."""
    tree = _tree()
    fam = genomes.simulate_genomes_family(tree, duplication=0.2, loss=0.2, initial_families=12, seed=4)
    t0 = tree.complete_tree.nodes[tree.complete_tree.root].birth_time
    assert fam.summary()["events"] == event_counts(fam.edges, t0)


# --- 3. the summary describes the directory it sits in --------------------------------------------

def test_the_written_summary_does_not_count_ancestral_sequences_it_did_not_write(tmp_path):
    """Ancestral sequences are reconstructed either way but only written when asked for, so a default
    run reported a count beside a directory that had none. Whoever inherits the folder cannot tell
    "never written" from "lost in transfer"."""
    tree = _tree()
    g = genomes.simulate_genomes_family(tree, duplication=0.2, loss=0.2, initial_families=8, seed=1)
    r = sequences.simulate_sequences(g, model=jc69(), length=30, seed=1)

    r.write(tmp_path / "default")
    s = json.loads((tmp_path / "default" / "sequences_summary.json").read_text(encoding="utf-8"))
    assert not list((tmp_path / "default" / "ancestral").glob("*"))
    assert "ancestral_sequences" not in s, "the summary counts files that are not in the directory"

    r.write(tmp_path / "asked", outputs=("alignments", "ancestral", "summary"))
    s = json.loads((tmp_path / "asked" / "sequences_summary.json").read_text(encoding="utf-8"))
    assert s["ancestral_sequences"] > 0
    assert list((tmp_path / "asked" / "ancestral").glob("*")), "…and now they really are there"

    # `summary()` itself still describes the run, not the write — it is the Python API's answer
    assert r.summary()["ancestral_sequences"] > 0


# --- 4. a run nobody seeded still says which seed it drew -----------------------------------------

def test_an_unseeded_cli_run_announces_the_drawn_seed(tmp_path, capsys):
    """It was always in the report — but the report is a file the user has not opened. A class
    following a worksheet got 17, 19 and 15 families where the sheet said 20, and the only clue was
    that the answers were wrong."""
    main(["species", str(tmp_path / "a"), "--birth", "1", "--death", "0.3", "--n-extant", "10",
          "--quiet"])
    err = capsys.readouterr().err
    assert "no --seed given" in err and "drew" in err
    drawn = int(err.split("drew")[1].split()[0])

    # and it is the seed that was used: replaying it gives the same tree
    main(["species", str(tmp_path / "b"), "--birth", "1", "--death", "0.3", "--n-extant", "10",
          "--seed", str(drawn), "--quiet"])
    assert (tmp_path / "a" / "species" / "species_complete.nwk").read_text(encoding="utf-8") == \
           (tmp_path / "b" / "species" / "species_complete.nwk").read_text(encoding="utf-8")


def test_a_seeded_run_says_nothing(tmp_path, capsys):
    """A healthy run prints nothing here, so the line means what it says when it appears."""
    main(["species", str(tmp_path / "s"), "--birth", "1", "--death", "0.3", "--n-extant", "10",
          "--seed", "3", "--quiet"])
    assert "drew" not in capsys.readouterr().err


# --- 5. the defaults are stated where they are used -----------------------------------------------

def test_the_help_states_the_four_rate_defaults(capsys):
    """A pasted command that lost its tail still ran, at rates ``--help`` did not name, and looked
    like a successful run of the command the user meant."""
    with pytest.raises(SystemExit):
        main(["genomes", "--help"])
    help_text = capsys.readouterr().out
    for flag, default in (("--duplication", "0.2"), ("--transfer", "0.1"),
                          ("--loss", "0.25"), ("--origination", "0.5")):
        where = help_text.index(flag)
        assert f"default {default}" in help_text[where:where + 320], (
            f"{flag}'s default ({default}) is not stated in its help")


def test_the_stated_defaults_are_the_ones_actually_used(tmp_path):
    """The help and the code have to agree, and only a run can say whether they do."""
    run = tmp_path / "r"
    main(["species", str(run), "--birth", "1", "--death", "0.3", "--n-extant", "8",
          "--seed", "1", "--quiet"])
    main(["genomes", str(run), "--seed", "1", "--quiet"])
    log = (run / "genomes" / "genomes.log").read_text(encoding="utf-8")
    for flag, default in (("duplication", "0.2"), ("transfer", "0.1"),
                          ("loss", "0.25"), ("origination", "0.5")):
        assert f"{flag}\t{default}" in log, f"{flag} did not run at the documented default"


# --- 6. scoring a gene tree against the species tree it grew in -----------------------------------

def _run_with_a_partial_family(tmp_path):
    """A run holding a family that occupies only part of the tree — the ordinary case."""
    run = tmp_path / "td"
    main(["species", str(run), "--birth", "1", "--death", "0.3", "--n-extant", "12",
          "--seed", "42", "--quiet"])
    main(["genomes", str(run), "--duplication", "0.05", "--transfer", "0.10", "--loss", "0.08",
          "--origination", "0.2", "--initial-families", "20", "--seed", "3", "--quiet"])
    species_tree = run / "species" / "species_extant.nwk"
    n_species = species_tree.read_text(encoding="utf-8").count(",") + 1
    for path in sorted((run / "genomes" / "gene_trees").glob("*_extant.nwk")):
        tips = len(set(re.findall(r"n\d+_g\d+", path.read_text(encoding="utf-8"))))
        if 3 < tips < n_species:
            return species_tree, path
    pytest.skip("no partially-present family in this run")


def test_treedist_restricts_to_the_shared_taxa(tmp_path, capsys):
    """The comparison this tool is most often reached for is the one it refused. Only a universal,
    single-copy family occupies every genome — a lecturer found not one of 22 families qualified, so
    the step in her practical was impossible until she retuned the rates."""
    species_tree, gene_tree = _run_with_a_partial_family(tmp_path)

    with pytest.raises(SystemExit):        # argparse's error() exits rather than returning
        main(["tools", "treedist", str(species_tree), str(gene_tree), "--metric", "rf"])
    assert "--restrict" in capsys.readouterr().err, "the refusal must name the way out"

    assert main(["tools", "treedist", str(species_tree), str(gene_tree),
                 "--metric", "all", "--restrict"]) == 0
    out = capsys.readouterr()
    assert "scored on the" in out.err
    assert {ln.split("\t")[0] for ln in out.out.strip().splitlines()} == {
        "rf", "rf-normalized", "branch-score"}


def test_restrict_changes_nothing_when_the_taxa_already_match(tmp_path, capsys):
    """A no-op where it should be one: the flag must not quietly score a different question."""
    run = tmp_path / "u"
    main(["species", str(run), "--birth", "1", "--death", "0.3", "--n-extant", "10",
          "--seed", "5", "--quiet"])
    tree = str(run / "species" / "species_extant.nwk")
    capsys.readouterr()                    # drop the species run's own "wrote …" line
    main(["tools", "treedist", tree, tree, "--metric", "all"])
    plain = capsys.readouterr().out
    main(["tools", "treedist", tree, tree, "--metric", "all", "--restrict"])
    assert capsys.readouterr().out == plain


def test_restricting_to_too_few_taxa_is_refused(tmp_path, capsys):
    """Two trees sharing two taxa have no topology to disagree about; a number there would be noise."""
    a, b = tmp_path / "a.nwk", tmp_path / "b.nwk"
    a.write_text("((A:1,B:1):1,(C:1,D:1):1);\n", encoding="utf-8")
    b.write_text("((A:1,B:1):1,(X:1,Y:1):1);\n", encoding="utf-8")
    with pytest.raises(SystemExit):
        main(["tools", "treedist", str(a), str(b), "--metric", "rf", "--restrict"])
    assert "share only 2 taxa" in capsys.readouterr().err


def test_prune_to_a_named_tip_set_keeps_a_real_dated_tree():
    """`--restrict` prunes rather than intersecting clade sets, so a length-aware metric still means
    something. The pruned tree must be bifurcating, carry exactly the kept tips, and merge the branch
    lengths across the nodes it suppressed."""
    from zombi2.tree import prune

    sp = species.simulate_species_tree(birth=1.0, death=0.3, n_extant=12, seed=3)
    complete = sp.complete_tree
    keep = set(sorted(n.id for n in complete.extant_leaves())[:5])
    sub = prune(complete, tips=keep)

    assert {n.id for n in sub.leaves()} == keep
    assert all(n.children is None or len(n.children) == 2 for n in sub.nodes.values())
    # depth is preserved: a kept tip sits where it always did
    for i in keep:
        assert sub.nodes[i].end_time == pytest.approx(complete.nodes[i].end_time)
    # and the extant-tree behaviour is untouched
    assert {n.id for n in prune(complete).leaves()} == {n.id for n in complete.extant_leaves()}
