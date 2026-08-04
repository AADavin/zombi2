"""Tests for the ``zombi2`` command line (species + genomes) and the ``read_newick`` reader it
loads trees with. The CLI is a thin shell over the level ``simulate_*`` functions, so these tests
exercise argument wiring, output files, ``--params``, and the clean error paths — not the engines
(which their own suites cover)."""

import collections
import os
import re

import pytest

from zombi2.cli.main import main
from zombi2.genomes import simulate_genomes_family
from zombi2.sequences.substitution_models import AMINO_ACIDS
from zombi2.species import simulate_species_tree
from zombi2.tree import read_newick


# ── read_newick ─────────────────────────────────────────────────────────────────────

def test_read_newick_round_trips_a_complete_tree():
    # a run with survivors round-trips: ids, node count, and the extant/extinct split are preserved
    r = simulate_species_tree(birth=1.0, death=0.4, n_extant=40, seed=3)
    t = r.complete_tree
    back, names = read_newick(t.to_newick())
    assert names == {}                                           # a ZOMBI tree: labels ARE the ids
    assert set(back.nodes) == set(t.nodes)
    assert back.root == t.root
    assert len(back.extant_leaves()) == len(t.extant_leaves())
    assert len(back.extinct_leaves()) == len(t.extinct_leaves())
    # every branch length (duration) survives to Newick's 6 significant figures — the root's
    # included, so the stem is not silently dropped and the round-tripped tree keeps its full height
    for i, n in t.nodes.items():
        assert back.nodes[i].end_time - back.nodes[i].birth_time == pytest.approx(
            n.end_time - n.birth_time, rel=1e-5, abs=1e-9)
    assert back.nodes[back.root].end_time - back.nodes[back.root].birth_time > 0.0
    assert max(n.end_time for n in back.nodes.values()) == pytest.approx(
        max(n.end_time for n in t.nodes.values()), rel=1e-5)


def test_read_newick_zombi_tree_honours_an_authoritative_fate_table():
    # a ZOMBI tree records only branch lengths, so depth cannot tell an unsampled survivor (it sits at
    # the present) from an extant one. When the run's species_fates.tsv is passed, it is authoritative:
    # here a present-day survivor is declared unsampled and must come back unsampled, not extant.
    r = simulate_species_tree(birth=1.0, death=0.3, n_extant=12, seed=4)
    survivors = [n.id for n in r.complete_tree.extant_leaves()]
    names = r.complete_tree.labels()          # the table keys on the same label the tree carries
    fates = {names[i]: "extant" for i in survivors}
    fates[names[survivors[0]]] = "unsampled"                      # override one survivor
    for n in r.complete_tree.extinct_leaves():
        fates[names[n.id]] = "extinct"
    back, _ = read_newick(r.complete_tree.to_newick(), tip_fates=fates)
    assert back.nodes[survivors[0]].fate == "unsampled"          # the table won, not the depth-guess
    assert len(back.extant_leaves()) == len(survivors) - 1              # the unsampled one is no longer extant
    assert len(back.unsampled_leaves()) == 1


def test_read_newick_ultrametric_external_tree_is_all_extant_with_a_name_map():
    # ultrametric (every tip at depth 2) → every tip extant, and the user's labels come back mapped
    t, names = read_newick("((human:1,chimp:1):1,(mouse:0.8,rat:0.8):1.2);")
    assert len(t.extant_leaves()) == 4 and not t.extinct_leaves()
    assert all(n.fate == "speciation" for n in t.nodes.values() if n.children is not None)
    assert sorted(names.values()) == ["chimp", "human", "mouse", "rat"]
    assert t.nodes[t.root].birth_time == 0.0


def test_read_newick_nonultrametric_external_tree_refuses_to_guess():
    with pytest.raises(ValueError, match="not ultrametric"):
        read_newick("((a:1,b:1):1,c:1.5);")                      # c ends early: extinct or sampled?


def test_read_newick_nonultrametric_tree_uses_supplied_fates():
    t, names = read_newick("((a:1,b:1):1,c:1.5);", tip_fates={"a": "extant", "b": "extant",
                                                              "c": "extinct"})
    fate = {names[n.id]: n.fate for n in t.nodes.values() if n.children is None}
    assert fate == {"a": "extant", "b": "extant", "c": "extinct"}


@pytest.mark.parametrize("fates, msg", [
    ({"a": "extant", "b": "extant"}, "missing a fate"),          # c not covered
    ({"a": "extant", "b": "extant", "c": "extant", "z": "extinct"}, "not in the tree"),
    ({"a": "extant", "b": "extant", "c": "dead"}, "extant.*extinct"),
])
def test_read_newick_tip_fates_are_validated(fates, msg):
    with pytest.raises(ValueError, match=msg):
        read_newick("((a:1,b:1):1,c:1.5);", tip_fates=fates)


def test_species_fates_file_is_a_valid_tip_fates_input(tmp_path):
    # a species run's own species_fates.tsv is written in the --tip-fates format, so it feeds straight
    # back in: the 'lineage<TAB>fate' header is skipped and the 'unsampled' value is kept
    from zombi2.cli.genomes import _read_tip_fates
    r = simulate_species_tree(birth=1.0, death=0.4, n_extant=12, sampling=0.5, seed=5)
    r.write(tmp_path)
    parsed = _read_tip_fates(str(tmp_path / "species_fates.tsv"))
    assert "lineage" not in parsed                               # the header row did not leak in as a tip
    assert "unsampled" in set(parsed.values())                   # sampling<1 produced unsampled tips, kept
    names = r.complete_tree.labels()
    assert parsed == {names[n.id]: n.fate for n in r.complete_tree.leaves()}


def test_genomes_reads_the_runs_fate_table_so_unsampled_tips_are_not_extant(tmp_path, capsys):
    # the sampling handoff bug: without the fate table, genomes reads the tree by depth and counts every
    # survivor extant (extant + unsampled); consuming species_fates.tsv, it builds only the sampled ones
    run = tmp_path / "run"
    main(["species", str(run), "--birth", "1.0", "--death", "0.4", "--n-extant", "20",
          "--sampling", "0.5", "--seed", "7"])
    fates = (run / "species" / "species_fates.tsv").read_text(encoding="utf-8").splitlines()[1:]
    n_extant = sum(ln.endswith("\textant") for ln in fates)
    assert 0 < n_extant < 20                                     # sampling really thinned the survivors
    capsys.readouterr()
    main(["genomes", str(run), "--duplication", "0.1", "--loss", "0.1", "--origination", "0.5",
          "--seed", "1"])
    out = capsys.readouterr().out
    assert f"{n_extant} extant genomes" in out                   # the sampled survivors, not all 20

    # and with the table removed, it falls back to the depth-guess and over-counts (the old behaviour)
    (run / "species" / "species_fates.tsv").unlink()
    import shutil
    shutil.rmtree(run / "genomes")
    main(["genomes", str(run), "--duplication", "0.1", "--loss", "0.1", "--origination", "0.5",
          "--seed", "1"])
    assert "20 extant genomes" in capsys.readouterr().out


def test_unsampled_is_an_accepted_external_tip_fate():
    # 'unsampled' (a survivor not observed) is a legal fate for an external tree too, not just ZOMBI's own
    t, names = read_newick("((a:1,b:1):1,c:1.5);",
                           tip_fates={"a": "extant", "b": "unsampled", "c": "extinct"})
    fate = {names[n.id]: n.fate for n in t.nodes.values() if n.children is None}
    assert fate == {"a": "extant", "b": "unsampled", "c": "extinct"}


def test_read_newick_output_feeds_the_genomes_engine():
    r = simulate_species_tree(birth=1.0, death=0.3, n_extant=25, seed=5)
    back, _ = read_newick(r.complete_tree.to_newick())
    g = simulate_genomes_family(back, duplication=0.2, loss=0.2, origination=0.5, seed=9)
    assert g.profiles.shape[1] == len(back.extant_leaves())            # one column per extant tip


@pytest.mark.parametrize("bad, msg", [
    ("", "empty"),
    ("((a,b)", "unbalanced"),
    ("(a,b,c);", "bifurcating"),
    ("(a:x,b:1);", "branch length"),
])
def test_read_newick_rejects_malformed(bad, msg):
    with pytest.raises(ValueError, match=msg):
        read_newick(bad)


# ── zombi2 species ──────────────────────────────────────────────────────────────────

def test_species_run_writes_the_expected_files(tmp_path, capsys):
    rc = main(["species", str(tmp_path), "--birth", "1", "--death", "0.3", "--n-extant", "20", "--seed", "1", "--flat"])
    assert rc == 0
    written = {p.name for p in tmp_path.iterdir()}
    assert {"species_complete.nwk", "species_extant.nwk", "species_events.tsv",
            "species.log"} <= written
    assert "extant" in capsys.readouterr().out
    # the reproducibility log records the resolved parameters
    log = (tmp_path / "species.log").read_text(encoding="utf-8")
    assert "birth\t1.0" in log and "n_extant\t20" in log


def test_species_write_is_selective(tmp_path):
    main(["species", str(tmp_path), "--birth", "1", "--total-time", "3", "--seed", "1", "--write", "complete", "--flat"])
    nwk = {p.name for p in tmp_path.iterdir() if p.suffix == ".nwk"}
    assert nwk == {"species_complete.nwk"}                       # no extant/events file


def test_species_is_deterministic_given_the_seed(tmp_path):
    a, b = tmp_path / "a", tmp_path / "b"
    for out in (a, b):
        main(["species", str(out), "--birth", "1", "--death", "0.3", "--n-extant", "30", "--seed", "13", "--flat"])
    assert (a / "species_complete.nwk").read_text(encoding="utf-8") == (b / "species_complete.nwk").read_text(encoding="utf-8")


@pytest.mark.parametrize("argv", [
    ["species", "--n-extant", "5"],                              # no --birth
    ["species", "--birth", "1"],                                 # no stop condition
    ["species", "--birth", "1", "--n-extant", "5", "--total-time", "4"],  # both stops
])
def test_species_argument_errors_exit_2(argv, tmp_path):
    with pytest.raises(SystemExit) as e:
        main(argv + ["-o", str(tmp_path)])
    assert e.value.code == 2


def test_species_engine_error_is_reported_cleanly(tmp_path, capsys):
    # mass extinctions need a fixed end (--total-time); with --n-extant the engine raises, and the
    # CLI must report it as a one-line error (rc 1), never a traceback
    rc = main(["species", str(tmp_path), "--birth", "1", "--n-extant", "10", "--mass-extinction", "3", "0.5", "--flat"])
    assert rc == 1
    assert "zombi2: error:" in capsys.readouterr().err


def test_engine_error_is_respelt_in_flags(tmp_path, capsys, tree_file):
    # the engines raise for a Python caller ("pass trim_overlaps=True"), which is a dead end at a
    # shell prompt. The CLI must hand back the flag that reaches the same keyword.
    gff = tmp_path / "overlapping.gff"
    gff.write_text("##sequence-region c1 1 900\n"
                   "c1\t.\tgene\t100\t400\t.\t+\t.\tID=a\n"
                   "c1\t.\tgene\t300\t600\t.\t+\t.\tID=b\n", encoding="utf-8")
    rc = main(["genomes", str(tmp_path / "g"), "--from", str(tree_file),
               "--resolution", "nucleotide", "--gff", str(gff), "--flat"])
    assert rc == 1
    err = capsys.readouterr().err
    assert "--trim-overlaps" in err and "trim_overlaps=True" not in err


def test_an_equals_sign_in_the_data_is_left_alone(tmp_path, capsys, tree_file):
    # only names the command actually has are respelt, so a '=' inside a filename (or a mapping in
    # a rate expression) survives the trip verbatim
    rc = main(["genomes", str(tmp_path / "g"), "--from", str(tree_file),
               "--resolution", "nucleotide", "--gff", str(tmp_path / "no=such=file.gff"), "--flat"])
    assert rc == 1
    assert "no=such=file.gff" in capsys.readouterr().err


# ── zombi2 genomes ──────────────────────────────────────────────────────────────────

@pytest.fixture
def tree_file(tmp_path):
    """A species tree written to disk, for the genomes command to read."""
    main(["species", str(tmp_path), "--birth", "1", "--death", "0.3", "--n-extant", "25", "--seed", "1", "--flat"])
    return tmp_path / "species_complete.nwk"


def test_genomes_unordered_writes_events_and_profiles(tmp_path, tree_file):
    out = tmp_path / "g"
    rc = main(["genomes", str(out), "--from", str(tree_file), "--duplication", "0.2", "--transfer", "0.1", "--loss", "0.25", "--origination", "0.5", "--seed", "42", "--flat"])
    assert rc == 0
    # a genomes run written to its own directory carries the tree it evolved along, so it stays
    # replayable without the species run beside it
    written = {p.name for p in out.iterdir()}
    assert {"genome_events.tsv", "profiles.tsv", "genomes.tsv",
            "species_complete.nwk", "genomes.log"} <= written
    assert any(n.startswith("gene_tree_fam") for n in written)   # --flat, so not in a subdirectory


def test_genomes_ordered_writes_structured_outputs(tmp_path, tree_file):
    out = tmp_path / "g"
    rc = main(["genomes", str(out), "--from", str(tree_file), "--resolution", "ordered", "--duplication", "0.2", "--loss", "0.2", "--origination", "0.5", "--inversion", "0.3", "--chromosomes", "3", "--seed", "42", "--write", "gene_order", "events", "--flat"])
    assert rc == 0
    assert {p.name for p in out.iterdir()} == {"gene_order.tsv", "genome_events.tsv",
                                               "rearrangement_events.tsv", "species_complete.nwk",
                                               "species_fates.tsv", "genomes.log"}


def test_genomes_ordered_events_carry_where_each_one_happened(tmp_path, tree_file):
    # the genealogy and its positions are one table; the rearrangements are their own
    out = tmp_path / "g"
    rc = main(["genomes", str(out), "--from", str(tree_file), "--resolution", "ordered", "--duplication", "0.3", "--loss", "0.2", "--origination", "0.6", "--inversion", "0.4", "--seed", "42", "--write", "events", "--flat"])
    assert rc == 0
    assert not list(out.glob("rearrangements.tsv")) and not list(out.glob("genome_event_*.tsv"))
    rows = (out / "genome_events.tsv").read_text(encoding="utf-8").splitlines()
    cols = rows[0].split("\t")
    assert cols == ["time", "kind", "family", "parents", "children",
                    "chromosome", "start", "length", "dest_chromosome", "dest_position"]
    at = {c: i for i, c in enumerate(cols)}               # by name: a new column cannot shift these
    body = [r.split("\t") for r in rows[1:]]
    assert {r[at["kind"]] for r in body} >= {"origination", "duplication", "loss"}
    assert "inversion" not in {r[at["kind"]] for r in body}   # it ends no lineage: its own file
    assert [r for r in body if r[at["start"]]], "no event carries the place it happened"
    # a duplication is ONE row that ends one copy and starts two
    dups = [r for r in body if r[at["kind"]] == "duplication"]
    assert dups and all(len(r[at["parents"]].split(";")) == 1
                        and len(r[at["children"]].split(";")) == 2 for r in dups)
    rear = (out / "rearrangement_events.tsv").read_text(encoding="utf-8").splitlines()
    assert rear[0].split("\t")[:3] == ["time", "kind", "lineage"]
    assert {r.split("\t")[1] for r in rear[1:]} == {"inversion"}


def test_genomes_rejects_ordered_only_flag_under_unordered(tmp_path, tree_file):
    with pytest.raises(SystemExit) as e:
        main(["genomes", str(tmp_path / "g"), "--from", str(tree_file), "--inversion", "0.3", "--flat"])
    assert e.value.code == 2


def test_genomes_rejects_write_output_foreign_to_resolution(tmp_path, tree_file):
    with pytest.raises(SystemExit) as e:                         # gene_order is ordered-only
        main(["genomes", str(tmp_path / "g"), "--from", str(tree_file), "--write", "gene_order", "--flat"])
    assert e.value.code == 2


def test_genomes_nucleotide_runs_and_writes_its_own_outputs(tmp_path, tree_file):
    out = tmp_path / "g"
    rc = main(["genomes", str(out), "--from", str(tree_file), "--resolution", "nucleotide", "--root-length", "600", "--genes", "3", "--gene-length", "100", "--inversion", "1.0", "--duplication", "0.5", "--loss", "0.4", "--seed", "1", "--flat"])
    assert rc == 0
    # the nucleotide default is events + genes; blocks is opt-in
    written = {p.name for p in out.iterdir()}
    assert {"genome_events.tsv", "genes.tsv", "species_complete.nwk", "genomes.log"} <= written
    assert len((out / "genes.tsv").read_text(encoding="utf-8").splitlines()) > 1


def test_genomes_nucleotide_write_selects_blocks(tmp_path, tree_file):
    out = tmp_path / "g"
    rc = main(["genomes", str(out), "--from", str(tree_file), "--resolution", "nucleotide", "--root-length", "400", "--genes", "2", "--gene-length", "100", "--inversion", "1.0", "--seed", "1", "--write", "blocks", "events", "--flat"])
    assert rc == 0
    head = (out / "blocks.tsv").read_text(encoding="utf-8").splitlines()[0].split("\t")
    assert head == ["lineage", "chromosome", "position", "source", "start", "end", "strand",
                    "copy", "gene"]


def test_genomes_nucleotide_seeds_from_a_gff(tmp_path, tree_file):
    gff = tmp_path / "seed.gff"
    gff.write_text("##gff-version 3\n"
                   "##sequence-region chrom1 1 900\n"
                   "chrom1\tZOMBI2\tgene\t101\t200\t.\t+\t.\tID=dnaA\n"
                   "chrom1\tZOMBI2\tgene\t401\t500\t.\t-\t.\tID=recA\n", encoding="utf-8")
    out = tmp_path / "g"
    rc = main(["genomes", str(out), "--from", str(tree_file), "--resolution", "nucleotide", "--gff", str(gff), "--inversion", "1.0", "--seed", "1", "--flat"])
    assert rc == 0
    rows = [r.split("\t") for r in (out / "genes.tsv").read_text(encoding="utf-8").splitlines()[1:]]
    assert [r[1] for r in rows] == ["dnaA", "recA"]          # names survive to the output
    assert [r[5] for r in rows] == ["1", "-1"]               # ...and so does the coding strand


def test_genomes_is_deterministic_across_resolutions(tmp_path, tree_file):
    def run(tag):
        out = tmp_path / tag
        main(["genomes", str(out), "--from", str(tree_file), "--resolution", "nucleotide", "--root-length", "500", "--genes", "2", "--gene-length", "100", "--inversion", "1.0", "--duplication", "0.4", "--seed", "7", "--write", "events", "blocks", "--flat"])
        return {p.name: p.read_text(encoding="utf-8") for p in out.iterdir() if p.suffix == ".tsv"}
    assert run("a") == run("b")


@pytest.mark.parametrize("argv, why", [
    (["--initial-families", "5"], "nucleotide has no initial-families"),
    (["--replacement"], "nucleotide transfers are additive"),
    (["--loss", "0.2 * ByFamily(spread=0.5)"], "nucleotide wires OnTime and DrivenBy, not ByFamily"),
    (["--gff", "x.gff", "--genes", "3"], "gff and genes are mutually exclusive"),
    (["--write", "gene_order"], "gene_order is an ordered output"),
    (["--write", "profiles"], "the nucleotide resolution has no profiles"),
])
def test_genomes_nucleotide_rejects_foreign_options(tmp_path, tree_file, argv, why):
    with pytest.raises(SystemExit) as e:
        main(["genomes", str(tmp_path / "g"), "--from", str(tree_file), "--resolution", "nucleotide", *argv, "--flat"])
    assert e.value.code == 2, why


def test_genomes_nucleotide_accepts_a_skyline(tmp_path, tree_file):
    """The nucleotide resolution speaks the rate grammar now, so a skyline is written the same way
    there as anywhere else — one notation across the levels, not a per-resolution dialect."""
    rc = main(["genomes", str(tmp_path / "g"), "--from", str(tree_file), "--resolution", "nucleotide",
               "--root-length", "800", "--genes", "2", "--gene-length", "100",
               "--inversion", "3.0 * OnTime({0: 1.0, 1.0: 0.0})", "--seed", "1", "--flat", "--quiet"])
    assert rc == 0


@pytest.mark.parametrize("resolution", ["family", "ordered"])
def test_genomes_rejects_nucleotide_only_flags_elsewhere(tmp_path, tree_file, resolution):
    with pytest.raises(SystemExit) as e:                     # bp knobs need a nucleotide genome
        main(["genomes", str(tmp_path / "g"), "--from", str(tree_file), "--resolution", resolution, "--root-length", "500", "--flat"])
    assert e.value.code == 2


def test_sequences_reads_a_nucleotide_handoff_through_from(tmp_path, tree_file):
    # --from a nucleotide run works like any other: the handoff says which resolution wrote it (only
    # that one writes blocks.tsv), so nothing has to be repeated from the genomes command
    out = tmp_path / "g"
    main(["genomes", str(out), "--from", str(tree_file), "--resolution", "nucleotide", "--root-length", "400", "--genes", "2", "--gene-length", "100", "--seed", "1", "--flat"])
    s = tmp_path / "s"
    rc = main(["sequences", str(s), "--from", str(out), "--model", "jc69", "--seed", "1", "--flat"])
    assert rc == 0
    assert list(s.glob("genome_n*.fasta"))            # assembled genomes, which only this level has
    assert list(s.glob("block*.fasta"))               # blocks, not families: the files say so


def test_an_emptied_genome_run_says_so_on_stderr(tmp_path, capsys):
    """There is no genome floor at the family resolution, so a high --loss can strip every lineage of
    every gene. That is a real outcome, but nothing in the outputs shows it — profiles.tsv is simply
    empty — so the command says so. It goes to stderr and survives --quiet, because a scripted batch is
    exactly the caller who needs to hear that its data is degenerate."""
    run = tmp_path / "run"
    main(["species", str(run), "--birth", "1.0", "--death", "0.3", "--n-extant", "8", "--seed", "1",
          "--quiet"])
    capsys.readouterr()
    rc = main(["genomes", str(run), "--loss", "5.0", "--duplication", "0", "--transfer", "0",
               "--origination", "0", "--initial-families", "10", "--seed", "3", "--quiet"])
    err = capsys.readouterr().err
    assert rc == 0                                        # a degenerate run is not a failed run
    assert "empty" in err and "--loss" in err             # what happened, and the knob to turn


def test_sequences_survives_a_run_with_no_families(tmp_path):
    """A run with no gene families has no alignments, so mean pairwise identity is undefined. The run
    report rendered it as a percentage regardless and raised `TypeError`, after every output file had
    already been written — reachable from an emptied genome and from a plain --initial-families 0."""
    run = tmp_path / "run"
    main(["species", str(run), "--birth", "1.0", "--death", "0.3", "--n-extant", "8", "--seed", "1",
          "--quiet"])
    main(["genomes", str(run), "--initial-families", "0", "--loss", "0.1", "--duplication", "0",
          "--transfer", "0", "--origination", "0", "--seed", "1", "--quiet"])
    rc = main(["sequences", str(run), "--length", "60", "--seed", "1", "--quiet"])
    assert rc == 0
    assert (run / "run.zombi2").exists()                   # the report was written, not crashed on


def test_genomes_missing_tree_is_reported_cleanly(tmp_path, capsys):
    rc = main(["genomes", str(tmp_path / "g"), "--from", str(tmp_path / "nope.nwk"), "--duplication", "0.1", "--flat"])
    assert rc == 1
    assert "tree file not found" in capsys.readouterr().err


def test_genomes_on_ultrametric_external_tree_writes_a_name_map(tmp_path):
    (tmp_path / "ext.nwk").write_text("((human:1,chimp:1):1,(mouse:0.8,rat:0.8):1.2);\n")
    out = tmp_path / "g"
    rc = main(["genomes", str(out), "--from", str(tmp_path / "ext.nwk"), "--duplication", "0.3", "--origination", "1.0", "--seed", "1", "--flat"])
    assert rc == 0
    # all four tips are observed, so the profile matrix has four columns
    assert len((out / "profiles.tsv").read_text(encoding="utf-8").splitlines()[0].split("\t")) == 1 + 4
    # names.tsv maps ZOMBI's n<id> back to the user's labels
    mapped = dict(row.split("\t") for row in (out / "names.tsv").read_text(encoding="utf-8").splitlines()[1:])
    assert sorted(mapped.values()) == ["chimp", "human", "mouse", "rat"]


def test_genomes_on_nonultrametric_tree_needs_tip_fates(tmp_path, capsys):
    (tmp_path / "ext.nwk").write_text("((a:1,b:1):1,c:1.5);\n")     # c ends early
    rc = main(["genomes", str(tmp_path / "g"), "--from", str(tmp_path / "ext.nwk"), "--duplication", "0.3", "--seed", "1", "--flat"])
    assert rc == 1
    assert "not ultrametric" in capsys.readouterr().err


def test_genomes_nonultrametric_tree_runs_with_tip_fates_file(tmp_path):
    (tmp_path / "ext.nwk").write_text("((a:1,b:1):1,c:1.5);\n")
    (tmp_path / "fates.tsv").write_text("a\textant\nb\textant\nc\textinct\n", encoding="utf-8")
    out = tmp_path / "g"
    rc = main(["genomes", str(out), "--from", str(tmp_path / "ext.nwk"), "--tip-fates", str(tmp_path / "fates.tsv"), "--duplication", "0.3", "--origination", "1.0", "--seed", "1", "--flat"])
    assert rc == 0
    # c is extinct, so only a and b are observed → two profile columns
    assert len((out / "profiles.tsv").read_text(encoding="utf-8").splitlines()[0].split("\t")) == 1 + 2


# ── zombi2 sequences ────────────────────────────────────────────────────────────────

@pytest.fixture
def genomes_dir(tmp_path):
    """A completed species→genomes run on disk, for the sequences command to replay."""
    main(["species", str(tmp_path), "--birth", "1", "--death", "0.3", "--n-extant", "25", "--seed", "1", "--flat"])
    gdir = tmp_path / "g"
    main(["genomes", str(gdir), "--from", str(tmp_path / "species_complete.nwk"), "--duplication", "0.2", "--transfer", "0.1", "--loss", "0.25", "--origination", "0.6", "--seed", "42", "--flat"])
    return gdir


def test_sequences_writes_alignments_and_phylograms_by_default(tmp_path, genomes_dir):
    out = tmp_path / "s"
    rc = main(["sequences", str(out), "--from", str(genomes_dir), "--model", "hky85", "--kappa", "2", "--length", "300", "--seed", "1", "--flat"])
    assert rc == 0
    names = {p.name for p in out.iterdir()}
    assert "sequences.log" in names
    # the default write set is alignments + phylograms; ancestral / species_phylogram are opt-in
    assert any(n.startswith("fam") for n in names)
    assert any(n.startswith("phylogram_fam") for n in names)
    assert not any(n.startswith("sequences_ancestral_fam") for n in names)
    assert not any(n.startswith("sequences_species_phylogram") for n in names)


def test_sequences_write_selects_ancestral_and_species_phylogram(tmp_path, genomes_dir):
    out = tmp_path / "s"
    rc = main(["sequences", str(out), "--from", str(genomes_dir), "--model", "jc69", "--length", "200", "--seed", "1", "--write", "ancestral", "species_phylogram", "--flat"])
    assert rc == 0
    names = {p.name for p in out.iterdir()}
    # the species phylogram is produced only because the CLI hands the engine the species tree
    assert "clock_species_tree_complete.nwk" in names
    assert any(n.startswith("sequences_ancestral_fam") for n in names)
    assert not any(n.startswith("fam") for n in names)   # not requested


def test_sequences_relaxed_clock_runs_and_is_logged(tmp_path, genomes_dir):
    # the relaxed clock is not its own flag: it is a ByLineage modifier on the substitution rate
    out = tmp_path / "s"
    rc = main(["sequences", str(out), "--from", str(genomes_dir), "--model", "gtr", "--frequencies", "0.3", "0.2", "0.2", "0.3", "--substitution", "1.0 * ByLineage(spread=0.4, dist='gamma')", "--seed", "1", "--flat"])
    assert rc == 0
    log = (out / "sequences.log").read_text(encoding="utf-8")
    assert "gamma lineage clock, spread 0.4" in log
    # the rate is logged in its written form, so the log line pastes back into --substitution
    assert "substitution\t1.0 * ByLineage(spread=0.4, dist='gamma')" in log


def test_sequences_rejects_a_model_foreign_parameter(tmp_path, genomes_dir):
    with pytest.raises(SystemExit) as e:                         # --kappa is meaningless for jc69
        main(["sequences", str(tmp_path / "s"), "--from", str(genomes_dir), "--model", "jc69", "--kappa", "2", "--flat"])
    assert e.value.code == 2


@pytest.mark.parametrize("model", ["poisson", "jtt", "dayhoff", "wag", "lg"])
def test_sequences_protein_models_write_amino_acid_alignments(tmp_path, genomes_dir, model):
    out = tmp_path / model
    rc = main(["sequences", str(out), "--from", str(genomes_dir), "--model", model, "--length", "60", "--seed", "1", "--flat"])
    assert rc == 0
    fasta = next(p for p in out.iterdir() if p.name.startswith("fam"))
    residues = set("".join(ln for ln in fasta.read_text(encoding="utf-8").splitlines() if not ln.startswith(">")))
    assert residues <= set(AMINO_ACIDS) and not residues <= set("ACGT")


@pytest.mark.parametrize("flag,value", [("--kappa", ["2"]), ("--frequencies", ["0.25"] * 4),
                                        ("--exchangeabilities", ["1"] * 6)])
def test_sequences_protein_models_take_no_parameters(tmp_path, genomes_dir, flag, value):
    # an empirical matrix has nothing to tune, so a nucleotide knob must be rejected, not ignored
    with pytest.raises(SystemExit) as e:
        main(["sequences", str(tmp_path / "s"), "--from", str(genomes_dir), "--model", "lg", flag, *value, "--flat"])
    assert e.value.code == 2


def test_sequences_is_deterministic_given_the_seed(tmp_path, genomes_dir):
    a, b = tmp_path / "a", tmp_path / "b"
    for out in (a, b):
        main(["sequences", str(out), "--from", str(genomes_dir), "--model", "hky85", "--length", "250", "--seed", "7", "--flat"])
    assert _dir_seq_text(a) == _dir_seq_text(b)


def _dir_seq_text(d):
    return {p.name: p.read_text(encoding="utf-8") for p in d.iterdir() if p.suffix in (".fasta", ".nwk")}


def test_sequences_missing_genomes_dir_is_reported_cleanly(tmp_path, capsys):
    rc = main(["sequences", str(tmp_path / "s"), "--from", str(tmp_path / "nope"), "--model", "jc69", "--flat"])
    assert rc == 1
    assert "holds no genomes run" in capsys.readouterr().err


def test_sequences_needs_the_genome_event_log(tmp_path, capsys):
    # a genomes run written with --write profiles has no event log to replay
    main(["species", str(tmp_path), "--birth", "1", "--death", "0.3", "--n-extant", "15", "--seed", "1", "--flat"])
    gdir = tmp_path / "g"
    main(["genomes", str(gdir), "--from", str(tmp_path / "species_complete.nwk"), "--duplication", "0.2", "--seed", "1", "--write", "profiles", "--flat"])
    rc = main(["sequences", str(tmp_path / "s"), "--from", str(gdir), "--model", "jc69", "--flat"])
    assert rc == 1
    assert "holds no genomes run" in capsys.readouterr().err


# ── traits ──────────────────────────────────────────────────────────────────────────

def test_traits_continuous_writes_values_and_tree(tmp_path, tree_file):
    out = tmp_path / "t"
    rc = main(["traits", "--kind", "continuous", str(out), "--from", str(tree_file), "--rate", "1.0", "--seed", "1", "--flat"])
    assert rc == 0
    assert {p.name for p in out.iterdir()} == {"trait_values.tsv", "trait_tree.nwk", "traits.log",
                                               "trait_summary.json"}
    header, first = (out / "trait_values.tsv").read_text(encoding="utf-8").splitlines()[:2]
    assert header == "node\tkind\ttrait"
    assert first.split("\t")[0].startswith("n")          # n<id>, matching the Newick
    float(first.split("\t")[2])                          # a continuous trait is a number (3rd column)


def test_traits_ou_and_threshold_run(tmp_path, tree_file):
    # OU: the same diffusion pulled to an optimum
    assert main(["traits", "--kind", "continuous", str(tmp_path / "ou"), "--from", str(tree_file), "--rate", "1.0", "--reverts-to", "2", "--pull", "0.5", "--seed", "1", "--flat"]) == 0
    # the threshold model: a discrete state read off a continuous liability
    out = tmp_path / "th"
    assert main(["traits", str(out), "--from", str(tree_file), "--kind", "discrete", "--states", "absent,present", "--liability", "1.0", "--threshold", "0.0", "--seed", "1", "--flat"]) == 0
    states = {ln.split("\t")[2] for ln in (out / "trait_values.tsv").read_text(encoding="utf-8").splitlines()[1:]}
    assert states <= {"absent", "present"}


def test_traits_continuous_ou_with_a_modified_rate(tmp_path, tree_file):
    # the OU knobs and a σ² modifier in one run — the combination the engine used to refuse. Nothing
    # about the command changed: --rate takes the same written form it always did, so this is the
    # end-to-end check that the unlocked combination reaches the engine from the command line.
    out = tmp_path / "ou_eb"
    rc = main(["traits", "--kind", "continuous", str(out), "--from", str(tree_file),
               "--rate", "1.0 * OnTime({0: 4.0, 1: 1.0})", "--reverts-to", "2.0", "--pull", "0.5",
               "--seed", "1", "--flat"])
    assert rc == 0
    assert (out / "trait_values.tsv").exists()
    log = (out / "traits.log").read_text(encoding="utf-8")
    assert "rate\t1.0 * OnTime({0: 4, 1: 1})" in log            # the rate, as the run resolved it
    assert "reverts_to\t2.0" in log and "pull\t0.5" in log       # alongside the OU knobs, not instead


def test_traits_discrete_writes_the_event_log(tmp_path, tree_file):
    out = tmp_path / "t"
    rc = main(["traits", str(out), "--from", str(tree_file), "--kind", "discrete", "--states", "marine,terrestrial", "--switch", "0.3", "--start", "marine", "--seed", "1", "--write", "values", "events", "tree", "--flat"])
    assert rc == 0
    assert {p.name for p in out.iterdir()} == {"trait_values.tsv", "trait_events.tsv",
                                               "trait_tree.nwk", "traits.log"}
    lines = (out / "trait_events.tsv").read_text(encoding="utf-8").splitlines()
    assert lines[0] == "time\tkind\tlineage\tfrom\tto"
    assert lines[1].split("\t")[1] == "initial"       # the t=0 row, the conditioning file's anchor


def test_traits_at_speciation_logs_on_speciation_changes(tmp_path, tree_file):
    out = tmp_path / "t"
    main(["traits", "--kind", "continuous", str(out), "--from", str(tree_file), "--rate", "1.0", "--at-speciation", "0.5", "--seed", "1", "--write", "events", "--flat"])
    rows = (out / "trait_events.tsv").read_text(encoding="utf-8").splitlines()[1:]
    kinds = [r.split("\t")[1] for r in rows]
    assert kinds[0] == "initial"                                              # t=0, then the jumps
    assert kinds[1:] and all(k == "on_speciation" for k in kinds[1:])         # a diffusion has no
    #                                                    along-branch events, only the split jumps


@pytest.mark.parametrize("argv, msg", [
    (["--kind", "discrete", "--states", "a,b", "--switch", "0.1", "--rate", "2"], "--kind continuous"),
    (["--switch", "0.1"], "--kind discrete"),                       # discrete knob, continuous run
    (["--kind", "discrete", "--switch", "0.1"], "--states"),        # discrete without a state space
    (["--kind", "discrete", "--states", "a,b"], "--switch"),        # discrete without a model
    (["--write", "bogus"], "invalid choice"),                       # not a write token
    (["--start", "marine"], "must be a number"),                    # a label on a continuous trait
])
def test_traits_argument_errors_exit_2(argv, msg, tmp_path, tree_file, capsys):
    with pytest.raises(SystemExit) as e:
        main(["traits", "--kind", "continuous", str(tmp_path / "t"), "--from", str(tree_file), *argv, "--flat"])
    assert e.value.code == 2
    assert msg in capsys.readouterr().err


def test_traits_is_deterministic_given_the_seed(tmp_path, tree_file):
    written = []
    for name in ("a", "b"):
        out = tmp_path / name
        main(["traits", str(out), "--from", str(tree_file), "--kind", "discrete", "--states", "a,b", "--switch", "0.4", "--seed", "99", "--flat"])
        written.append((out / "trait_values.tsv").read_text(encoding="utf-8"))
    assert written[0] == written[1]


def test_traits_missing_tree_is_reported_cleanly(tmp_path, capsys):
    rc = main(["traits", "--kind", "continuous", str(tmp_path / "t"), "--from", str(tmp_path / "nope.nwk"), "--rate", "1", "--flat"])
    assert rc == 1                                        # a clean one-line error, not a traceback
    assert "tree file not found" in capsys.readouterr().err


def test_traits_on_external_tree_writes_a_name_map(tmp_path):
    (tmp_path / "ext.nwk").write_text("((human:1,chimp:1):1,(mouse:0.8,rat:0.8):1.2);\n")
    out = tmp_path / "t"
    rc = main(["traits", "--kind", "continuous", str(out), "--from", str(tmp_path / "ext.nwk"), "--rate", "1.0", "--seed", "1", "--flat"])
    assert rc == 0
    names = dict(ln.split("\t") for ln in (out / "names.tsv").read_text(encoding="utf-8").splitlines()[1:])
    assert sorted(names.values()) == ["chimp", "human", "mouse", "rat"]
    # trait_values now carries every node; the name map covers the tips, so each named tip is among them
    nodes = {ln.split("\t")[0] for ln in (out / "trait_values.tsv").read_text(encoding="utf-8").splitlines()[1:]}
    assert set(names) <= nodes


# ── --switch takes the written form of a rate ───────────────────────────────────────
#
# It used to be `type=float`, so argparse killed a rate expression before ZOMBI2 saw it —
# "invalid float value: '0.2 * DrivenBy(...)'" — while the engine had supported a driven switch rate
# all along. These pin the four shapes it now reads, and that a driven one is really driven.

def _switch_run(tmp_path, name, tree_file, switch, seed="3"):
    out = tmp_path / name
    rc = main(["traits", str(out), "--from", str(tree_file), "--kind", "discrete",
               "--states", "x,y", "--switch", switch, "--seed", seed, "--flat"])
    assert rc == 0, switch
    return out


def _driver_file(tmp_path, tree_file):
    """A discrete trait grown first, written out — the file a conditioned switch rate reads."""
    out = _switch_run(tmp_path, "habitat", tree_file, "0.5", seed="1")
    return out / "trait_events.tsv"


@pytest.mark.parametrize("switch", [
    "0.2",                                              # a bare number: the symmetric rate
    "{'x->y': 0.4, 'y->x': 0.1}",                       # only the named transitions
    "[[0, 0.4], [0.1, 0]]",                             # the k x k matrix
])
def test_traits_switch_takes_every_shape_the_keyword_does(switch, tmp_path, tree_file):
    out = _switch_run(tmp_path, re.sub(r"\W", "", switch)[:12] or "n", tree_file, switch)
    states = {ln.split("\t")[2] for ln in
              (out / "trait_values.tsv").read_text(encoding="utf-8").splitlines()[1:]}
    assert states <= {"x", "y"}


def test_traits_switch_takes_a_rate_expression_and_the_driver_reaches_the_engine(tmp_path, tree_file):
    from zombi2.rates.parse import parse_rate, written_form

    driver = _driver_file(tmp_path, tree_file)
    expression = f"0.2 * DrivenBy({str(driver)!r}, {{'x': 12.0, 'y': 1.0}})"
    plain = _switch_run(tmp_path, "plain", tree_file, "0.2")
    driven = _switch_run(tmp_path, "driven", tree_file, expression)
    # the driver's own states are x/y here (a second trait over the same labels), so every lineage in
    # state x switches twelve times as fast — the two runs cannot come out the same
    events = lambda d: (d / "trait_events.tsv").read_text(encoding="utf-8")
    assert events(plain) != events(driven)
    # the run records what drove it, by content as well as by name, so the pair can be reproduced.
    # Compare against the written form rather than a hand-built string: a Windows path is escaped in it
    log = (driven / "traits.log").read_text(encoding="utf-8")
    assert f"switch\t{written_form(parse_rate(expression))}" in log
    assert "input\t" in log and str(driver) in log      # the digest names the file as it is on disk


def test_traits_a_driven_switch_inside_a_per_transition_dict_is_seen_by_the_log(tmp_path, tree_file):
    """A rate hidden in a ``{'a->b': rate}`` dict is still a driver: the log has to render it in the
    written form and the input digest has to name the file, or a conditioned run reads as an
    unconditioned one and its driver can be re-run out from under it."""
    from zombi2.rates.parse import parse_rate, written_form

    driver = _driver_file(tmp_path, tree_file)
    inner = f"0.3 * DrivenBy({str(driver)!r}, {{'x': 5.0}})"
    out = _switch_run(tmp_path, "dict_driven", tree_file, f"{{'x->y': {inner}, 'y->x': 0.1}}")
    log = (out / "traits.log").read_text(encoding="utf-8")
    # the entry is rendered as a rate, not as a Rate repr, which no flag would take back
    assert f"switch\t{{'x->y': {written_form(parse_rate(inner))}, 'y->x': 0.1}}" in log
    assert "input\t" in log and str(driver) in log


@pytest.mark.parametrize("switch, msg", [
    ("0.2 * OnDiversity(3)", "unknown name"),               # a rate-grammar error, reported as one
    ("nonsense", "unknown name"),
])
def test_traits_a_bad_switch_expression_is_a_rate_error(switch, msg, tmp_path, tree_file, capsys):
    with pytest.raises(SystemExit) as e:
        main(["traits", str(tmp_path / "t"), "--from", str(tree_file), "--kind", "discrete",
              "--states", "x,y", "--switch", switch, "--seed", "1", "--flat"])
    assert e.value.code == 2
    assert msg in capsys.readouterr().err


def test_traits_liability_reads_the_grammar_and_the_engine_refuses_what_it_cannot_honour(tmp_path, tree_file, capsys):
    """The other bare-float rate flag. A modified liability variance-rate is not implemented, and the
    point of parsing it here is that the *engine* gets to say so — argparse's "invalid float value"
    named neither the model nor the way out."""
    base = ["traits", str(tmp_path / "l"), "--from", str(tree_file), "--kind", "discrete",
            "--states", "x,y", "--threshold", "0", "--seed", "1", "--flat"]
    assert main([*base, "--liability", "1.0"]) == 0
    assert main([*base, "--liability", "PerLineage(1.0)"]) == 0
    assert main([*base, "--liability", "1.0 * OnTime({0: 1.0, 2: 0.5})"]) == 1
    assert "modified liability variance-rate is not implemented" in capsys.readouterr().err
    # several liabilities are the correlated-traits form, which has no flag for the correlation
    with pytest.raises(SystemExit) as e:
        main([*base, "--liability", "{'size': 1.0, 'limb': 2.0}"])
    assert e.value.code == 2
    assert "--liability is one variance-rate" in capsys.readouterr().err


def test_traits_params_file_can_carry_a_switch_expression(tmp_path, tree_file):
    """A TOML value is a *default*, and argparse runs `type` over a string default — which is what
    makes the SPEC's "one written form, everywhere" true of --switch as well as --birth."""
    driver = _driver_file(tmp_path, tree_file)
    # a TOML *literal* string ('''), so the driver path is taken as written — a basic "..." string
    # processes backslash escapes, which eats a Windows path before ZOMBI2 sees it
    (tmp_path / "p.toml").write_text(
        'kind = "discrete"\nstates = "x,y"\nseed = 4\n'
        "switch = \'\'\'0.2 * DrivenBy(\'%s\', {\'x\': 9.0, \'y\': 1.0})\'\'\'\n" % driver,
        encoding="utf-8")
    out = tmp_path / "fromparams"
    assert main(["traits", str(out), "--from", str(tree_file), "--params", str(tmp_path / "p.toml"),
                 "--flat"]) == 0
    assert "switch\t0.2 * DrivenBy(" in (out / "traits.log").read_text(encoding="utf-8")


def test_traits_params_file_drives_the_run_and_cli_overrides(tmp_path, tree_file):
    # a [traits] table scopes one file to this command (so one file can drive a whole pipeline)
    (tmp_path / "p.toml").write_text('[traits]\nkind = "discrete"\n'
                                     'states = "marine,terrestrial"\nswitch = 0.15\n'
                                     'write = ["values", "events"]\nseed = 7\n', encoding="utf-8")
    argv = ["traits", "--params", str(tmp_path / "p.toml"), "--from", str(tree_file)]
    out = tmp_path / "a"
    # --kind is required, but the file supplies it — that is why it is validated in run() rather
    # than marked argparse-`required`, which no default could satisfy
    assert main(["traits", str(out), *argv[1:], "--flat"]) == 0
    assert {p.name for p in out.iterdir()} == {"trait_values.tsv", "trait_events.tsv", "traits.log"}
    states = {ln.split("\t")[2] for ln in (out / "trait_values.tsv").read_text(encoding="utf-8").splitlines()[1:]}
    assert states <= {"marine", "terrestrial"}          # the file's states reached the engine

    # a flag given on the command line still wins over the file
    other = tmp_path / "b"
    assert main([*argv, str(other), "--seed", "8", "--flat"]) == 0
    assert (other / "trait_values.tsv").read_text(encoding="utf-8") != (out / "trait_values.tsv").read_text(encoding="utf-8")


# ── --params ────────────────────────────────────────────────────────────────────────

def test_the_log_records_the_input_tree_by_content_not_only_by_name(tmp_path):
    # two different trees under the same filename give logs whose *parameters* are identical: same
    # command line, same seed, same rates. The digest is the only thing that tells the two runs
    # apart. (`result` is dropped with it — a different tree of course yields a different number of
    # families, and it is the parameter block this test is about.)
    logs = []
    for seed in (1, 2):
        sp, gn = tmp_path / f"sp{seed}", tmp_path / f"gn{seed}"
        main(["species", str(sp), "--birth", "1", "--death", "0.3", "--n-extant", "10",
              "--seed", str(seed), "--flat", "--quiet"])
        tree = tmp_path / "tree.nwk"                  # the SAME name, twice, holding two trees
        tree.write_text((sp / "species_complete.nwk").read_text(encoding="utf-8"))
        main(["genomes", str(gn), "--from", str(tree), "--duplication", "0.2", "--loss", "0.2",
              "--seed", "7", "--flat", "--quiet"])
        logs.append([ln for ln in (gn / "genomes.log").read_text(encoding="utf-8").splitlines()
                     if not ln.startswith(("timestamp", "command_line", "run\t", "source\t",
                                           "result\t"))])
    digests = [[ln for ln in log if ln.startswith("input\t")] for log in logs]
    assert len(digests[0]) == 1 and digests[0] != digests[1]     # the tree, and it is not the same tree
    assert [ln for ln in logs[0] if not ln.startswith("input\t")] == \
           [ln for ln in logs[1] if not ln.startswith("input\t")]   # nothing else differs at all


def test_the_log_records_a_driver_file_as_an_input(tmp_path):
    # a conditioned run reads the driver as surely as it reads the tree, so it is logged the same way
    run = tmp_path / "r"
    main(["species", str(run), "--birth", "1", "--death", "0.2", "--n-extant", "10", "--seed", "1",
          "--quiet"])
    main(["traits", str(run), "--kind", "discrete", "--states", "cave,surface", "--switch", "0.4",
          "--seed", "1", "--quiet"])
    driver = str(run / "traits" / "trait_events.tsv")
    main(["genomes", str(run), "--duplication", "0.2", "--seed", "1", "--quiet",
          "--loss", f"0.25 * DrivenBy({driver!r}, {{'cave': 4.0, 'surface': 1.0}})"])
    inputs = [ln.split("\t")[2] for ln in (run / "genomes" / "genomes.log").read_text(encoding="utf-8").splitlines()
              if ln.startswith("input\t")]
    assert driver in inputs


def test_params_file_supplies_defaults_and_cli_overrides(tmp_path):
    (tmp_path / "p.toml").write_text("birth = 2.0\ndeath = 0.3\nn-extant = 12\n", encoding="utf-8")
    out = tmp_path / "o"
    # birth comes from the file; the command line still overrides it
    main(["species", str(out), "--params", str(tmp_path / "p.toml"), "--birth", "1.0", "--seed", "1", "--flat"])
    log = (out / "species.log").read_text(encoding="utf-8")
    assert "birth\t1.0" in log and "n_extant\t12" in log


def test_params_file_scopes_by_command_table(tmp_path):
    (tmp_path / "pipeline.toml").write_text(
        "[species]\nbirth = 1.0\nn-extant = 15\n\n[genomes]\nduplication = 0.2\nwrite = "
        '["profiles"]\n', encoding="utf-8")
    sp, gn = tmp_path / "sp", tmp_path / "gn"
    main(["species", str(sp), "--params", str(tmp_path / "pipeline.toml"), "--seed", "1", "--flat"])
    main(["genomes", str(gn), "--params", str(tmp_path / "pipeline.toml"), "--from", str(sp / "species_complete.nwk"), "--seed", "1", "--flat"])
    assert {p.name for p in gn.iterdir()} == {"profiles.tsv",
                                              "species_complete.nwk", "species_fates.tsv", "genomes.log"}


def test_params_unknown_key_errors(tmp_path):
    (tmp_path / "bad.toml").write_text("birth = 1.0\nbogus = 3\n", encoding="utf-8")
    with pytest.raises(SystemExit) as e:
        main(["species", str(tmp_path / "o"), "--params", str(tmp_path / "bad.toml"), "--n-extant", "5", "--flat"])
    assert e.value.code == 2


def test_params_mistyped_section_errors(tmp_path, capsys):
    # a typo'd [table] must not be silently dropped (which would run every rate at its default)
    (tmp_path / "typo.toml").write_text("[speces]\nbirth = 1.0\nn-extant = 10\n", encoding="utf-8")
    with pytest.raises(SystemExit):
        main(["species", str(tmp_path / "o"), "--params", str(tmp_path / "typo.toml"), "--seed", "1", "--flat"])
    err = capsys.readouterr().err
    assert "unknown section" in err and "[speces]" in err


def test_params_top_level_key_broadcasts_under_a_table(tmp_path):
    # a top-level scalar is a shared base for every command; a [command] table overrides on conflict
    (tmp_path / "p.toml").write_text("seed = 99\n[species]\nbirth = 1.0\nn-extant = 8\n", encoding="utf-8")
    out = tmp_path / "o"
    main(["species", str(out), "--params", str(tmp_path / "p.toml"), "--flat"])
    assert "seed\t99" in (out / "species.log").read_text(encoding="utf-8")      # the shared seed was applied, not dropped


def test_params_command_table_overrides_the_shared_base(tmp_path):
    # on conflict the [command] table wins over a top-level key of the same name
    (tmp_path / "p.toml").write_text("seed = 1\nn-extant = 5\n[species]\nbirth = 1.0\nn-extant = 12\n", encoding="utf-8")
    out = tmp_path / "o"
    main(["species", str(out), "--params", str(tmp_path / "p.toml"), "--flat"])
    log = (out / "species.log").read_text(encoding="utf-8")
    assert "n_extant\t12" in log and "seed\t1" in log


def test_params_append_option_is_overridden_by_the_command_line(tmp_path):
    # --mass-extinction is an 'append' action: a params default plus a CLI flag must NOT concatenate
    (tmp_path / "p.toml").write_text(
        "[species]\nbirth = 1.0\ntotal-time = 5.0\nmass-extinction = [[2.0, 0.9]]\n", encoding="utf-8")
    out = tmp_path / "o"
    main(["species", str(out), "--params", str(tmp_path / "p.toml"),
          "--mass-extinction", "3.0", "0.1", "--seed", "3", "--flat"])
    log = (out / "species.log").read_text(encoding="utf-8")
    assert "mass_extinction\t[[3.0, 0.1]]" in log              # only the command line's pulse, not both


def test_params_last_of_two_files_wins(tmp_path):
    # two --params: the last file's values are used (and it is the one the log names)
    (tmp_path / "a.toml").write_text("[species]\nbirth = 1.0\nn-extant = 5\n", encoding="utf-8")
    (tmp_path / "b.toml").write_text("[species]\nbirth = 1.0\nn-extant = 40\n", encoding="utf-8")
    out = tmp_path / "o"
    main(["species", str(out), "--params", str(tmp_path / "a.toml"),
          "--params", str(tmp_path / "b.toml"), "--seed", "1", "--flat"])
    assert "n_extant\t40" in (out / "species.log").read_text(encoding="utf-8")


def test_params_invalid_choice_errors_cleanly(tmp_path, genomes_dir):
    # a bad choices= value in --params must be a clean error, not a KeyError deep in the command
    (tmp_path / "bad.toml").write_text('[genomes]\nresolution = "unordred"\nduplication = 0.1\n', encoding="utf-8")
    with pytest.raises(SystemExit):
        main(["genomes", str(tmp_path / "o"), "--from", str(genomes_dir / "species_complete.nwk"),
              "--params", str(tmp_path / "bad.toml"), "--seed", "1", "--flat"])


def test_params_can_supply_the_sequences_model(tmp_path, genomes_dir):
    # --model is not argparse-required, so a --params file can supply it (like --birth on species)
    (tmp_path / "seq.toml").write_text('[sequences]\nmodel = "hky85"\nlength = 200\n', encoding="utf-8")
    main(["sequences", str(genomes_dir), "--params", str(tmp_path / "seq.toml"), "--seed", "1"])
    assert (genomes_dir / "sequences" / "sequences.log").exists()


# ── rates in their written form (SPEC §5) ───────────────────────────────────────────
#
# Every rate flag takes the same expression the Python API takes, so there is one notation for a
# rate across Python, the command line and a --params file. These tests hold that line: the
# expression reaches the engine, it changes the run, an unwired modifier is refused, and the
# parameters log records something you can paste back.

def test_species_takes_a_rate_expression_and_it_bends_the_tree(tmp_path):
    # a skyline that collapses speciation at t=2 must give a smaller tree than the flat rate,
    # i.e. the modifier reached the engine rather than being parsed and dropped
    # seed 2 survives both flat and skyline to the present (the collapsed skyline rate can otherwise
    # let a run go fully extinct, which is now refused)
    flat, skyline = tmp_path / "flat", tmp_path / "sky"
    main(["species", str(flat), "--birth", "1.0", "--death", "0.2", "--total-time", "6", "--seed", "2", "--flat"])
    main(["species", str(skyline), "--birth", "1.0 * OnTime({0: 1.0, 2: 0.05})", "--death", "0.2", "--total-time", "6", "--seed", "2", "--flat"])
    n = {d: len(read_newick((d / "species_complete.nwk").read_text(encoding="utf-8"))[0].nodes)
         for d in (flat, skyline)}
    assert n[skyline] < n[flat]


def test_species_takes_a_scope_wrapper(tmp_path):
    # Global(base) is one budget for the whole tree: linear growth, so far fewer lineages
    out = tmp_path / "g"
    rc = main(["species", str(out), "--birth", "Global(2.0)", "--total-time", "5", "--seed", "1", "--flat"])
    assert rc == 0
    assert "birth\tGlobal(2.0)" in (out / "species.log").read_text(encoding="utf-8")


def test_species_records_the_rate_in_its_written_form(tmp_path):
    # the log line is the flag value again — a reproducibility record you can paste back
    out = tmp_path / "o"
    main(["species", str(out), "--birth", "1.0 * OnTime({0: 1.0, 3: 0.3})", "--total-time", "4", "--seed", "1", "--flat"])
    assert "birth\t1.0 * OnTime({0: 1, 3: 0.3})" in (out / "species.log").read_text(encoding="utf-8")


def test_species_refuses_a_modifier_it_does_not_wire(tmp_path, capsys):
    # DrivenBy would return a factor of 1.0 at this level — a run quietly not the model asked for.
    # (Driving speciation is a *joint* model, so it is `zombi2 joint`'s business, not this command's.)
    rc = main(["species", str(tmp_path / "o"), "--birth",
               "1.0 * DrivenBy('habitat.tsv', Table({'a': 2.0}))", "--total-time", "3",
               "--seed", "1", "--flat"])
    assert rc == 1
    assert "does not support" in capsys.readouterr().err


def test_a_typo_in_a_modifier_is_caught_at_the_flag(tmp_path, capsys):
    with pytest.raises(SystemExit) as e:
        main(["species", str(tmp_path / "o"), "--birth", "1.0 * OnDiversity(cap=10)", "--total-time", "3", "--flat"])
    assert e.value.code == 2
    assert "did you mean 'OnTotalDiversity'" in capsys.readouterr().err


def test_a_rate_expression_is_never_executed(tmp_path, capsys):
    with pytest.raises(SystemExit) as e:
        main(["species", str(tmp_path / "o"), "--birth", "__import__('os').system('true')", "--total-time", "3", "--flat"])
    assert e.value.code == 2
    assert "only call a scope or a modifier" in capsys.readouterr().err


def test_genomes_takes_a_rate_expression(tmp_path, tree_file):
    out = tmp_path / "g"
    rc = main(["genomes", str(out), "--from", str(tree_file), "--duplication", "0.2", "--loss", "0.25 * OnTime({0: 1.0, 2: 3.0})", "--origination", "0.5", "--seed", "42", "--flat"])
    assert rc == 0
    assert "loss\t0.25 * OnTime({0: 1, 2: 3})" in (out / "genomes.log").read_text(encoding="utf-8")


@pytest.fixture
def driver_file(tmp_path, tree_file):
    """A discrete habitat trait grown on ``tree_file`` and written as its event log — the file a
    conditioned ``DrivenBy`` names as its source and replays against the shared tree."""
    main(["traits", str(tmp_path), "--kind", "discrete", "--from", str(tree_file), "--states", "competent,normal", "--switch", "0.4", "--seed", "1", "--write", "events", "--flat"])
    return tmp_path / "trait_events.tsv"


def test_genomes_transfer_can_be_driven_from_the_cli(tmp_path, driver_file, tree_file):
    # the DONOR side: a rate, so it takes the ordinary written form and changes how much HGT happens
    out = tmp_path / "g"
    rc = main(["genomes", str(out), "--from", str(tree_file), "--initial-families", "5", "--transfer", f"0.2 * DrivenBy('{driver_file}', {{'competent': 4.0, 'normal': 1.0}})", "--seed", "2", "--flat"])
    assert rc == 0
    # the log records the rate in its written form, which is what pastes back into the flag — so
    # compare against that rather than a hand-built string (a Windows path is escaped in it)
    from zombi2.rates.parse import parse_rate, written_form
    written = written_form(parse_rate(f"0.2 * DrivenBy('{driver_file}', "
                                      f"{{'competent': 4.0, 'normal': 1.0}})"))
    assert f"transfer\t{written}" in (out / "genomes.log").read_text(encoding="utf-8")


def test_genomes_transfer_to_takes_a_driven_recipient_weight(tmp_path, driver_file, tree_file):
    # the RECIPIENT side: the choice slot, so the modifier is written on its own, with no base
    out = tmp_path / "g"
    rc = main(["genomes", str(out), "--from", str(tree_file), "--initial-families", "5", "--transfer", "0.5", "--transfer-to", f"DrivenBy('{driver_file}', {{'competent': 2.0, 'normal': 1.0}})", "--seed", "2", "--flat"])
    assert rc == 0
    from zombi2.cli.genomes import _transfer_to
    written = repr(_transfer_to(f"DrivenBy('{driver_file}', "
                                f"{{'competent': 2.0, 'normal': 1.0}})"))
    assert f"transfer_to\t{written}" in (out / "genomes.log").read_text(encoding="utf-8")


@pytest.mark.parametrize("resolution", ["ordered", "nucleotide"])
def test_genomes_transfer_to_takes_a_driven_weight_at_every_resolution(
        tmp_path, driver_file, tree_file, resolution):
    """The choice slot is not the family resolution's own: the same flag, the same written form and
    the same conditioning record at all three. ``--transfer-to`` was already parsed for every
    resolution and already passed to every engine — the engines were what refused it."""
    out = tmp_path / resolution
    extra = ["--root-length", "2000"] if resolution == "nucleotide" else ["--initial-families", "5"]
    rc = main(["genomes", str(out), "--from", str(tree_file), "--resolution", resolution,
               "--transfer", "0.5", "--transfer-to",
               f"DrivenBy('{driver_file}', {{'competent': 2.0, 'normal': 1.0}})",
               "--seed", "2", "--flat", *extra])
    assert rc == 0
    from zombi2.cli.genomes import _transfer_to
    written = repr(_transfer_to(f"DrivenBy('{driver_file}', "
                                f"{{'competent': 2.0, 'normal': 1.0}})"))
    assert f"transfer_to\t{written}" in (out / "genomes.log").read_text(encoding="utf-8")


def test_genomes_transfer_to_rejects_a_rate_expression(tmp_path, tree_file, capsys):
    with pytest.raises(SystemExit) as e:
        main(["genomes", str(tmp_path / "g"), "--from", str(tree_file), "--transfer-to", "1.0 * DrivenBy('d.tsv', {'a': 2})", "--flat"])
    assert e.value.code == 2
    assert "written on its own" in capsys.readouterr().err


def test_genomes_transfer_to_names_its_rules_for_a_misspelt_one(tmp_path, tree_file, capsys):
    # the flag lost argparse's `choices`, so it owes the reader the list itself
    with pytest.raises(SystemExit) as e:
        main(["genomes", str(tmp_path / "g"), "--from", str(tree_file), "--transfer-to", "uniforn", "--flat"])
    assert e.value.code == 2
    assert "'uniform', 'distance', or a DrivenBy" in capsys.readouterr().err


def test_genomes_params_file_carries_a_driven_transfer_to(tmp_path, driver_file, tree_file):
    # the written form is the same text in the file as on the flag — one notation (SPEC §5)
    # a TOML *literal* string ('''), so the driver path is taken as written — a basic
    # "..." string processes backslash escapes, which eats a Windows path before ZOMBI2 sees it
    (tmp_path / "p.toml").write_text(
        "transfer = 0.5\ninitial-families = 5\n"
        "transfer-to = \'\'\'DrivenBy(\'%s\', {\'competent\': 2.0, \'normal\': 1.0})\'\'\'\n" % driver_file,
        encoding="utf-8")
    out = tmp_path / "g"
    rc = main(["genomes", str(out), "--params", str(tmp_path / "p.toml"), "--from", str(tree_file), "--seed", "2", "--flat"])
    assert rc == 0
    assert "transfer_to\tDrivenBy(" in (out / "genomes.log").read_text(encoding="utf-8")


def test_traits_takes_a_rate_expression(tmp_path, tree_file):
    out = tmp_path / "t"
    rc = main(["traits", "--kind", "continuous", str(out), "--from", str(tree_file), "--rate", "1.0 * FromParent(spread=0.2)", "--seed", "1", "--flat"])
    assert rc == 0
    assert "rate\t1.0 * FromParent(spread=0.2)" in (out / "traits.log").read_text(encoding="utf-8")


def test_params_file_takes_a_rate_expression(tmp_path):
    # the same text as the flag, quoted as a TOML string — no second notation for a rate
    (tmp_path / "p.toml").write_text(
        'birth = "1.0 * OnTime({0: 1.0, 3: 0.3})"\ndeath = 0.3\ntotal-time = 5\n', encoding="utf-8")
    out = tmp_path / "o"
    rc = main(["species", str(out), "--params", str(tmp_path / "p.toml"), "--seed", "2", "--flat"])
    assert rc == 0
    assert "birth\t1.0 * OnTime({0: 1, 3: 0.3})" in (out / "species.log").read_text(encoding="utf-8")


def test_params_file_rate_expression_matches_the_flag(tmp_path):
    (tmp_path / "p.toml").write_text('birth = "1.0 * OnTotalDiversity(cap=20)"\n', encoding="utf-8")
    viafile, viaflag = tmp_path / "f", tmp_path / "c"
    main(["species", str(viafile), "--params", str(tmp_path / "p.toml"), "--total-time", "5", "--seed", "3", "--flat"])
    main(["species", str(viaflag), "--birth", "1.0 * OnTotalDiversity(cap=20)", "--total-time", "5", "--seed", "3", "--flat"])
    assert (viafile / "species_complete.nwk").read_text(encoding="utf-8") == \
        (viaflag / "species_complete.nwk").read_text(encoding="utf-8")


def test_the_rates_help_lists_only_what_the_level_wires(capsys):
    # the help is built from each level's IMPLEMENTED_MODIFIERS, so it cannot advertise the unwired
    for command, present, absent in [
            ("species", ["FromParent", "ByLineage"], ["ByFamily"]),   # both per-lineage forms wired
            # both clocks and the trait driver are wired; the diversity covariate is not
            ("sequences", ["ByLineage", "FromParent", "DrivenBy"], ["OnTotalDiversity"])]:
        with pytest.raises(SystemExit):
            main([command, "--help", "--flat"])
        out = capsys.readouterr().out
        block = out[out.index("RATES"):]
        assert all(p in block for p in present) and all(a not in block for a in absent)


# ── top-level dispatch ──────────────────────────────────────────────────────────────

def test_version_and_help_do_not_crash(capsys):
    for flag in ("--version", "--help"):
        with pytest.raises(SystemExit) as e:
            main([flag, "--flat"])
        assert e.value.code == 0
    assert "ZOMBI2" in capsys.readouterr().out


# ── output layout: grouped by default, one directory under --flat ───────────────────────────

def _tree_for(root):
    main(["species", str(root), "--birth", "1", "--death", "0.3", "--n-extant", "8", "--seed", "1"])
    return str(root / "species" / "species_complete.nwk")


def test_each_level_writes_into_its_own_directory(tmp_path):
    tree = _tree_for(tmp_path)
    main(["genomes", str(tmp_path), "--from", tree, "--initial-families", "4", "--duplication", "0.3", "--seed", "2", "--write", "events", "profiles", "gene_trees"])
    main(["sequences", str(tmp_path), "--from", str(tmp_path), "--model", "jc69", "--length", "20", "--seed", "1"])
    main(["traits", "--kind", "continuous", str(tmp_path), "--from", tree, "--rate", "1.0", "--seed", "1"])

    assert (tmp_path / "species" / "species_complete.nwk").exists()
    assert (tmp_path / "genomes" / "genome_events.tsv").exists()
    assert (tmp_path / "traits" / "trait_values.tsv").exists()
    # per-family outputs nest again: a hundred families would otherwise bury the rest
    assert list((tmp_path / "genomes" / "gene_trees").glob("gene_tree_fam*.nwk"))
    assert list((tmp_path / "sequences" / "alignments").glob("*.fasta"))
    assert list((tmp_path / "sequences" / "phylograms").glob("*.nwk"))
    # each run log sits with the level that wrote it
    for level, log in (("species", "species.log"), ("genomes", "genomes.log"),
                       ("sequences", "sequences.log"), ("traits", "traits.log")):
        assert (tmp_path / level / log).exists()
    # the one file at the top is the run report — it spans every level, so it belongs to the run, not
    # to any one level's directory
    assert [p.name for p in tmp_path.iterdir() if p.is_file()] == ["run.zombi2"]


def test_flat_puts_everything_in_one_directory(tmp_path):
    main(["species", str(tmp_path), "--birth", "1", "--death", "0.3", "--n-extant", "8", "--seed", "1", "--flat"])
    tree = str(tmp_path / "species_complete.nwk")
    assert (tmp_path / "species.log").exists()
    main(["genomes", str(tmp_path), "--from", tree, "--initial-families", "4", "--duplication", "0.3", "--seed", "2", "--flat", "--write", "events", "gene_trees"])
    assert (tmp_path / "genome_events.tsv").exists()
    assert list(tmp_path.glob("gene_tree_fam*.nwk"))            # not in a subdirectory
    assert not (tmp_path / "genomes").exists()
    assert (tmp_path / "genomes.log").exists()          # the log is loose too under --flat


def test_sequences_handoff_takes_either_layout(tmp_path):
    # grouped: point --genomes at the run directory, whose genomes/ holds the handoff files
    grouped = tmp_path / "grouped"
    tree = _tree_for(grouped)
    main(["genomes", str(grouped), "--from", tree, "--initial-families", "3", "--duplication", "0.3", "--seed", "2"])
    assert main(["sequences", str(grouped), "--from", str(grouped), "--model", "jc69", "--length", "20", "--seed", "1"]) == 0

    # flat: point it at the directory holding them directly
    flat = tmp_path / "flat"
    main(["species", str(flat), "--birth", "1", "--death", "0.3", "--n-extant", "8", "--seed", "1", "--flat"])
    main(["genomes", str(flat), "--from", str(flat / "species_complete.nwk"), "--initial-families", "3", "--duplication", "0.3", "--seed", "2", "--flat"])
    assert main(["sequences", str(flat), "--from", str(flat), "--model", "jc69", "--length", "20", "--seed", "1", "--flat"]) == 0
    assert list(flat.glob("fam*.fasta"))


def test_the_two_layouts_write_the_same_files(tmp_path):
    outs = {}
    for name, extra in (("grouped", []), ("flat", ["--flat"])):
        root = tmp_path / name
        main(["species", str(root), "--birth", "1", "--death", "0.3", "--n-extant", "8", "--seed", "1", *extra])
        tree = root / ("species_complete.nwk" if extra else "species/species_complete.nwk")
        main(["genomes", str(root), "--from", str(tree), "--initial-families", "4", "--duplication", "0.3", "--seed", "2", "--write", "events", "profiles", "gene_trees", *extra])
        outs[name] = {p.name for p in root.rglob("*") if p.is_file()}
    # the data files match — only the directories differ. The one extra is run.zombi2: the grouped
    # layout writes a report over its per-level records, which --flat (one commingled directory, no
    # per-level records to read) has no way to build.
    assert outs["grouped"] - {"run.zombi2"} == outs["flat"]
    assert "run.zombi2" in outs["grouped"] and "run.zombi2" not in outs["flat"]


def test_a_run_directory_is_read_as_well_as_written(tmp_path):
    # the whole point of the positional: name the directory once and each level finds the last
    main(["species", str(tmp_path), "--birth", "1", "--death", "0.3", "--n-extant", "8", "--seed", "1"])
    assert main(["genomes", str(tmp_path), "--initial-families", "3", "--seed", "2"]) == 0
    assert main(["traits", "--kind", "continuous", str(tmp_path), "--rate", "1.0", "--seed", "1"]) == 0
    assert main(["sequences", str(tmp_path), "--model", "jc69", "--length", "20", "--seed", "1"]) == 0
    assert (tmp_path / "genomes" / "genome_events.tsv").exists()
    assert (tmp_path / "traits" / "trait_values.tsv").exists()
    assert list((tmp_path / "sequences" / "alignments").glob("*.fasta"))


def test_a_flat_run_directory_is_read_as_well_as_written(tmp_path):
    main(["species", str(tmp_path), "--birth", "1", "--death", "0.3", "--n-extant", "8",
          "--seed", "1", "--flat"])
    assert main(["genomes", str(tmp_path), "--initial-families", "3", "--seed", "2", "--flat"]) == 0
    assert main(["sequences", str(tmp_path), "--model", "jc69", "--length", "20", "--seed", "1",
                 "--flat"]) == 0
    assert (tmp_path / "genome_events.tsv").exists()


def test_from_takes_a_newick_file(tmp_path):
    (tmp_path / "ext.nwk").write_text("((a:1,b:1):1,(c:1,d:1):1);\n")
    run = tmp_path / "r"
    assert main(["genomes", str(run), "--from", str(tmp_path / "ext.nwk"),
                 "--initial-families", "3", "--seed", "2"]) == 0
    assert (run / "genomes" / "genome_events.tsv").exists()


def test_from_reads_one_run_and_writes_another(tmp_path):
    # replicates off a single species tree: read out/, write rep1/ and rep2/, leaving out/ alone
    src = tmp_path / "src"
    main(["species", str(src), "--birth", "1", "--death", "0.3", "--n-extant", "10", "--seed", "1"])
    for name, seed in (("rep1", "1"), ("rep2", "2")):
        assert main(["genomes", str(tmp_path / name), "--from", str(src),
                     "--initial-families", "4", "--duplication", "0.3", "--seed", seed]) == 0
    assert not (src / "genomes").exists()                    # the source run is untouched
    a = (tmp_path / "rep1" / "genomes" / "genome_events.tsv").read_text(encoding="utf-8")
    b = (tmp_path / "rep2" / "genomes" / "genome_events.tsv").read_text(encoding="utf-8")
    assert a != b                                            # different seeds, same tree

    # and sequences the same way, reading a genomes run and writing elsewhere
    assert main(["sequences", str(tmp_path / "seqs"), "--from", str(tmp_path / "rep1"),
                 "--model", "jc69", "--length", "20", "--seed", "1"]) == 0
    assert list((tmp_path / "seqs" / "sequences" / "alignments").glob("*.fasta"))


def test_a_directory_without_a_tree_says_so(tmp_path, capsys):
    (tmp_path / "empty").mkdir()
    rc = main(["genomes", str(tmp_path / "empty"), "--seed", "1"])
    assert rc == 1
    err = capsys.readouterr().err
    assert "holds no species tree" in err and "species_complete.nwk" in err


def test_a_directory_without_a_genomes_run_says_so(tmp_path, capsys):
    (tmp_path / "empty").mkdir()
    rc = main(["sequences", str(tmp_path / "empty"), "--model", "jc69"])
    assert rc == 1
    assert "holds no genomes run" in capsys.readouterr().err


def test_the_directory_and_the_file_give_the_same_run(tmp_path):
    main(["species", str(tmp_path), "--birth", "1", "--death", "0.3", "--n-extant", "10", "--seed", "1"])
    a, b = tmp_path / "a", tmp_path / "b"
    main(["genomes", str(a), "--from", str(tmp_path), "--initial-families", "4",
          "--duplication", "0.3", "--seed", "9"])
    main(["genomes", str(b), "--from", str(tmp_path / "species" / "species_complete.nwk"),
          "--initial-families", "4", "--duplication", "0.3", "--seed", "9"])
    assert (a / "genomes" / "genome_events.tsv").read_text(encoding="utf-8") == \
           (b / "genomes" / "genome_events.tsv").read_text(encoding="utf-8")


# ── zombi2 joint ────────────────────────────────────────────────────────────────────────────

def test_joint_trait_writes_both_levels(tmp_path):
    # BiSSE: the trait drives which lineages speciate, so neither level can be grown first
    rc = main(["joint", str(tmp_path),
               "--birth", "1.0 * DrivenBy('trait', {'small': 1.0, 'large': 3.0})", "--death", "0.2",
               "--states", "small,large", "--switch", "0.3", "--n-extant", "30", "--seed", "1"])
    assert rc == 0
    assert (tmp_path / "species" / "species_complete.nwk").exists()
    # the trait is written the way `zombi2 traits` writes it, not TraitsResult.write's bare default
    assert {p.name for p in (tmp_path / "traits").iterdir()} == {
        "trait_values.tsv", "trait_events.tsv", "trait_tree.nwk", "trait_summary.json"}
    states = {ln.split("\t")[2] for ln in
              (tmp_path / "traits" / "trait_values.tsv").read_text(encoding="utf-8").splitlines()[1:]}
    assert states <= {"small", "large"}


def test_joint_genome_driver_nests_its_gene_trees(tmp_path):
    rc = main(["joint", str(tmp_path),
               "--birth", "1.0 * DrivenBy('genomes:toxin', {'present': 3.0, 'absent': 1.0})",
               "--origination", "0.2", "--loss", "0.1", "--family-names", "toxin",
               "--n-extant", "20", "--seed", "1"])
    assert rc == 0
    written = {p.name for p in (tmp_path / "genomes").iterdir()}
    assert {"genome_events.tsv", "profiles.tsv", "genomes.tsv"} <= written
    # gene trees get their own directory here too, as they do under `zombi2 genomes`
    assert list((tmp_path / "genomes" / "gene_trees").glob("gene_tree_fam*.nwk"))


def test_joint_is_deterministic_given_the_seed(tmp_path):
    argv = ["--birth", "1.0 * DrivenBy('trait', {'a': 1.0, 'b': 2.0})",
            "--states", "a,b", "--switch", "0.3", "--n-extant", "20", "--seed", "5"]
    for name in ("x", "y"):
        main(["joint", str(tmp_path / name), *argv])
    assert (tmp_path / "x" / "species" / "species_complete.nwk").read_text(encoding="utf-8") == \
           (tmp_path / "y" / "species" / "species_complete.nwk").read_text(encoding="utf-8")


@pytest.mark.parametrize("argv, msg", [
    (["--birth", "1.0", "--n-extant", "10"], "needs a driver"),
    (["--birth", "1.0", "--states", "a,b", "--origination", "0.2", "--n-extant", "10"],
     "give one driver"),
    (["--birth", "1.0", "--origination", "0.2", "--switch", "0.3", "--n-extant", "10"],
     "need --states"),
    (["--birth", "1.0", "--states", "a", "--n-extant", "10"], "at least two"),
    (["--states", "a,b", "--n-extant", "10"], "--birth is required"),
    (["--birth", "1.0", "--states", "a,b"], "exactly one stop condition"),
])
def test_joint_argument_errors_exit_2(argv, msg, tmp_path, capsys):
    with pytest.raises(SystemExit) as e:
        main(["joint", str(tmp_path), *argv])
    assert e.value.code == 2
    assert msg in capsys.readouterr().err


def test_traits_requires_an_explicit_kind(tmp_path, capsys, tree_file):
    # the state space decides which other flags apply, so there is no default to fall back on
    with pytest.raises(SystemExit) as e:
        main(["traits", str(tmp_path / "t"), "--from", str(tree_file), "--rate", "1.0"])
    assert e.value.code == 2
    assert "--kind is required" in capsys.readouterr().err


# ── sequences on a nucleotide genome run ────────────────────────────────────────────

def _nucleotide_run(tmp_path, *, extra=()):
    """species → nucleotide genomes → sequences, all three from the command line."""
    run = str(tmp_path / "run")
    main(["species", run, "--birth", "1.0", "--death", "0.3", "--n-extant", "5",
          "--seed", "3", "--quiet"])
    main(["genomes", run, "--resolution", "nucleotide", "--root-length", "2000", "--genes", "5",
          "--gene-length", "150", "--inversion", "3.0", "--inversion-extent", "250",
          "--duplication", "1.0", "--loss", "1.0", "--seed", "3", "--quiet"])
    main(["sequences", run, "--model", "hky85", "--kappa", "3.0", "--substitution", "0.05",
          "--seed", "3", "--quiet", *extra])
    return tmp_path / "run" / "sequences"


def test_sequences_runs_on_a_nucleotide_handoff_and_assembles_the_genomes(tmp_path):
    # the handoff says what it is — blocks.tsv is the nucleotide resolution's and no other's — so
    # nothing has to be repeated from the genomes command
    out = _nucleotide_run(tmp_path)
    aln = {p.name for p in (out / "alignments").iterdir()}
    assert aln and all(n.startswith("block") for n in aln)    # per-block, in its own subfolder
    gdir = out / "genomes"                                    # the assembled genomes, grouped
    genomes = sorted(p.name for p in gdir.glob("genome_n*.fasta"))
    assert len(genomes) == 9                                  # one FASTA per node of the tree
    assert not [p for p in out.iterdir() if p.name.startswith("genome_n")]   # none left flat
    for name in genomes:
        text = (gdir / name).read_text(encoding="utf-8")
        assert text.startswith(">") and set("".join(text.splitlines()[1:])) <= set("ACGT")


def test_sequences_writes_a_genome_for_every_node_and_the_initial_one(tmp_path):
    # all of them by default, named by whose they are: no node is a special case, and the one that
    # belongs to no node is called "initial"
    out = _nucleotide_run(tmp_path)
    tree, _ = read_newick((tmp_path / "run" / "species" / "species_complete.nwk").read_text(encoding="utf-8"))
    gdir = out / "genomes"
    names = {p.name for p in gdir.iterdir()}
    assert {f"genome_n{i}.fasta" for i in tree.nodes} <= names
    assert "genome_initial.fasta" in names                    # the initial genome joins them here
    assert not [n for n in names if "ancestral" in n]
    initial = "".join((gdir / "genome_initial.fasta").read_text(encoding="utf-8").splitlines()[1:])
    assert len(initial) == 2000                               # the genome the run was seeded with


def test_sequences_matches_the_python_api_on_a_nucleotide_run(tmp_path):
    # the CLI is a shell: rebuilt from disk, it must give exactly what the objects give in memory
    from zombi2.genomes import simulate_genomes_nucleotide
    from zombi2.sequences import simulate_sequences
    from zombi2.sequences.substitution_models import hky85

    out = _nucleotide_run(tmp_path)
    sp = simulate_species_tree(birth=1.0, death=0.3, n_extant=5, seed=3)
    g = simulate_genomes_nucleotide(sp, root_length=2000, genes=5, gene_length=150, inversion=3.0,
                                    inversion_extent=250, duplication=1.0, loss=1.0, seed=3)
    r = simulate_sequences(g, model=hky85(kappa=3.0), substitution=0.05, seed=3)
    for lineage, chroms in r.genomes.items():
        text = (out / "genomes" / f"genome_{lineage}.fasta").read_text(encoding="utf-8")
        assert "".join(text.splitlines()[1:]) == "".join(chroms.values())


def test_sequences_rejects_length_and_an_intergene_model_where_they_do_not_apply(tmp_path):
    with pytest.raises(SystemExit):
        _nucleotide_run(tmp_path, extra=("--length", "300"))   # the genome sets the lengths
    run = str(tmp_path / "family")
    main(["species", run, "--birth", "1.0", "--n-extant", "4", "--seed", "1", "--quiet"])
    main(["genomes", run, "--initial-families", "4", "--duplication", "0.2", "--loss", "0.2",
          "--seed", "1", "--quiet"])
    with pytest.raises(SystemExit):                            # no spacer in a gene-family run
        main(["sequences", run, "--model", "jc69", "--intergene-model", "jc69", "--quiet"])


# ── zombi2 tools format (homology tables) ────────────────────────────────────────────

def _dtl_run(tmp_path, *, name="run", extra=()):
    """A species + genomes run with all of duplication/transfer/loss, so its gene trees carry every
    kind of internal node. Returns the run directory."""
    run = tmp_path / name
    main(["species", str(run), "--birth", "1", "--death", "0.3", "--n-extant", "12", "--seed", "3",
          "--quiet"])
    main(["genomes", str(run), "--duplication", "0.4", "--transfer", "0.3", "--loss", "0.25",
          "--origination", "0.6", "--seed", "7", "--quiet", *extra])
    return run


def _check_table(text):
    """A homology TSV is a square, symmetric O/P/X grid with a dashed diagonal and matching headers."""
    lines = text.rstrip("\n").split("\n")
    headers = lines[0].split("\t")
    assert headers[0] == ""                                     # the blank corner cell
    labels = headers[1:]
    rows = [ln.split("\t") for ln in lines[1:]]
    assert [r[0] for r in rows] == labels                      # row labels equal the column labels
    grid = [r[1:] for r in rows]
    n = len(labels)
    assert all(len(r) == n for r in grid)                      # square
    assert all(grid[i][i] == "-" for i in range(n))            # dashed diagonal
    assert all(grid[i][j] == grid[j][i] for i in range(n) for j in range(n))   # symmetric
    off = {grid[i][j] for i in range(n) for j in range(n) if i != j}
    # the divergence event, plus x for a transfer SINCE — including Tx, a second transfer below a first
    assert off <= {"S", "D", "T", "Sx", "Dx", "Tx"}
    return labels, grid


def test_tools_format_writes_homology_tables(tmp_path, capsys):
    run = _dtl_run(tmp_path)
    assert main(["tools", "format", str(run), "--format", "homology"]) == 0
    hdir = run / "genomes" / "homology"
    tables = sorted(hdir.glob("homology_fam*.tsv"))
    assert tables, "no homology tables written"
    # one table per family that has an extant gene tree (exactly the families with an _extant.nwk)
    table_fams = {p.stem.removeprefix("homology_fam") for p in tables}
    extant_fams = {p.name.split("_")[2].removeprefix("fam")
                   for p in (run / "genomes" / "gene_trees").glob("*_extant.nwk")}
    assert table_fams == extant_fams
    letters = set()
    for t in tables:
        _, grid = _check_table(t.read_text(encoding="utf-8"))
        letters |= {c for row in grid for c in row}
    # a DTL run exercises every combination: each divergence event, with and without a transfer since
    assert {"S", "D", "T", "Sx", "Dx"} <= letters
    assert "wrote" in capsys.readouterr().out


def test_tools_format_flat_writes_into_the_run_directory(tmp_path):
    run = _dtl_run(tmp_path)
    assert main(["tools", "format", str(run), "--format", "homology", "--flat"]) == 0
    assert (run / "homology_fam0.tsv").exists()                # straight in, no homology/ subdir
    assert not (run / "genomes" / "homology").exists()


def test_tools_format_reads_another_run_with_from(tmp_path):
    src = _dtl_run(tmp_path, name="src")
    dest = tmp_path / "analysis"
    assert main(["tools", "format", str(dest), "--from", str(src), "--format", "homology"]) == 0
    assert sorted((dest / "genomes" / "homology").glob("homology_fam*.tsv"))
    assert not (src / "genomes" / "homology").exists()         # the source run is left untouched


def test_tools_format_matches_the_python_api(tmp_path):
    # the CLI is a shell over the library: each table it writes must equal the classifier applied to
    # the gene trees rebuilt from the run's own recorded history (its species tree + event log).
    from zombi2.genomes.events import edges_from_tsv
    from zombi2.genomes.gene_trees import gene_trees_from_edges
    from zombi2.tools.homology import homology_tsv

    run = _dtl_run(tmp_path)
    main(["tools", "format", str(run), "--format", "homology"])
    tree, _ = read_newick((run / "species" / "species_complete.nwk").read_text(encoding="utf-8"))
    events = edges_from_tsv((run / "genomes" / "genome_events.tsv").read_text(encoding="utf-8"))
    trees = gene_trees_from_edges(events, tree)
    written = 0
    for fam, gt in trees.items():
        if gt.extant is None:
            continue
        on_disk = (run / "genomes" / "homology" / f"homology_fam{fam}.tsv").read_text(encoding="utf-8")
        assert on_disk == homology_tsv(gt.complete)
        written += 1
    assert written == len(list((run / "genomes" / "homology").glob("*.tsv")))


def test_tools_format_on_a_directory_without_a_run_errors_cleanly(tmp_path, capsys):
    (tmp_path / "empty").mkdir()
    assert main(["tools", "format", str(tmp_path / "empty"), "--format", "homology"]) == 1
    assert "holds no genomes run" in capsys.readouterr().err


def test_tools_format_rejects_an_unknown_format(tmp_path):
    run = _dtl_run(tmp_path)
    with pytest.raises(SystemExit):
        main(["tools", "format", str(run), "--format", "bogus"])


def test_tools_format_requires_a_format_and_names_the_choices(tmp_path, capsys):
    # no default: the three answer different questions, so picking one silently would answer a
    # question nobody asked. Bare, it must say so and name what is on offer.
    run = _dtl_run(tmp_path)
    with pytest.raises(SystemExit):
        main(["tools", "format", str(run)])
    err = capsys.readouterr().err
    assert "required: --format" in err
    assert all(name in err for name in ("homology", "markers", "recphylo"))
    assert not (run / "genomes" / "homology").exists()          # and nothing was written


def test_tools_format_reports_one_complete_path_per_output(tmp_path, capsys):
    # the line has to be pasteable as it stands: a whole path, not a base plus a fragment to compose
    run = _dtl_run(tmp_path)
    capsys.readouterr()                                         # drop the setup chatter
    assert main(["tools", "format", str(run), "--format", "homology", "markers"]) == 0
    lines = capsys.readouterr().out.strip().split("\n")
    n = len(list((run / "genomes" / "homology").glob("*.tsv")))
    # built with os.path.join/os.sep, not literal slashes: the point of the line is that it pastes,
    # so it has to be the platform's own separator all the way through (this caught a "/" on Windows)
    homology = os.path.join(str(run), "genomes", "homology") + os.sep
    markers = os.path.join(str(run), "genomes", "markers.tsv")
    assert lines == [f"wrote {n} table(s) in {homology}", f"wrote 1 table in {markers}"]
    assert (run / "genomes" / "markers.tsv").exists()           # one file needs no folder of its own


def _nucleotide_gene_run(tmp_path, *, name="nucrun"):
    """A nucleotide genomes run WITH declared genes, its duplication/transfer segments long enough to
    span a whole gene (so some families gain paralogs/xenologs). Returns the run directory."""
    run = tmp_path / name
    main(["species", str(run), "--birth", "1", "--death", "0.2", "--n-extant", "6", "--seed", "1",
          "--quiet"])
    main(["genomes", str(run), "--resolution", "nucleotide", "--root-length", "4000", "--genes", "6",
          "--gene-length", "200", "--duplication", "1.0", "--transfer", "0.8", "--loss", "0.3",
          "--duplication-extent", "260", "--transfer-extent", "260", "--seed", "1", "--quiet"])
    return run


def test_tools_format_on_a_nucleotide_run_writes_one_table_per_declared_gene(tmp_path):
    from zombi2.genomes.nucleotide import read_nucleotide_genomes
    from zombi2.tools.homology import homology_tsv

    run = _nucleotide_gene_run(tmp_path)
    assert main(["tools", "format", str(run), "--format", "homology"]) == 0
    tables = sorted((run / "genomes" / "homology").glob("homology_fam*.tsv"))
    tree, _ = read_newick((run / "species" / "species_complete.nwk").read_text(encoding="utf-8"))
    genome = read_nucleotide_genomes(run / "genomes", tree)
    # one table per DECLARED gene some node still carries — never the intergenic spacer
    declared = {f for f, gt in genome.gene_trees.items() if gt.extant is not None}
    assert {int(p.stem.removeprefix("homology_fam")) for p in tables} == declared
    assert declared <= set(genome.gene_spans)                  # every keyed family is a declared gene
    for p in tables:
        fam = int(p.stem.removeprefix("homology_fam"))
        _check_table(p.read_text(encoding="utf-8"))
        assert p.read_text(encoding="utf-8") == homology_tsv(genome.gene_trees[fam].complete)


def test_tools_format_refuses_a_nucleotide_run_with_no_declared_genes(tmp_path, capsys):
    run = tmp_path / "spacer"
    main(["species", str(run), "--birth", "1", "--death", "0.2", "--n-extant", "5", "--seed", "1",
          "--quiet"])
    main(["genomes", str(run), "--resolution", "nucleotide", "--root-length", "2000", "--genes", "0",
          "--duplication", "0.5", "--loss", "0.4", "--seed", "1", "--quiet"])   # --genes 0: all-intergenic
    assert main(["tools", "format", str(run), "--format", "homology"]) == 1
    assert "declared no genes" in capsys.readouterr().err


# ── run-directory echo and --quiet ──────────────────────────────────────────────────

def test_completion_message_has_no_double_slash_on_a_trailing_slash_dir(tmp_path, capsys):
    # a dir written with a trailing slash (as the quickstart's `out/` does) must echo as `out/`, not
    # `out//` — the run directory is normalised so the completion line reads cleanly
    rc = main(["species", str(tmp_path) + "/", "--birth", "1", "--death", "0.3",
               "--n-extant", "10", "--seed", "1", "--flat"])
    assert rc == 0
    out = capsys.readouterr().out
    assert "//" not in out
    assert f"wrote {tmp_path}/" in out


def test_tools_format_still_says_what_it_wrote_under_quiet(tmp_path, capsys):
    # it used to print nothing at all under --quiet, alone among the commands, so a scripted run had
    # to go looking on disk to find out whether anything had been written. --quiet takes away the
    # progress bar; the one line naming what landed is not progress.
    run = tmp_path / "run"
    main(["species", str(run), "--birth", "1", "--death", "0.3", "--n-extant", "8", "--seed", "1"])
    main(["genomes", str(run), "--duplication", "0.2", "--loss", "0.2", "--origination", "0.5",
          "--initial-families", "5", "--seed", "3"])
    capsys.readouterr()                                   # drop the setup chatter
    assert main(["tools", "format", str(run), "--format", "homology", "--quiet"]) == 0
    quiet = capsys.readouterr().out
    assert main(["tools", "format", str(run), "--format", "homology"]) == 0
    assert quiet.startswith("wrote ") and quiet == capsys.readouterr().out   # same line either way


# ── the cross-level staleness guard ─────────────────────────────────────────────────

def _pipeline(run, seed=1):
    main(["species", str(run), "--birth", "1", "--death", "0.3", "--n-extant", "8",
          "--seed", str(seed), "--quiet"])
    main(["genomes", str(run), "--duplication", "0.2", "--loss", "0.2", "--origination", "0.4",
          "--initial-families", "5", "--seed", str(seed), "--quiet"])
    main(["sequences", str(run), "--model", "jc69", "--length", "20", "--seed", str(seed), "--quiet"])


def test_rerunning_an_upstream_level_refuses_when_a_downstream_exists(tmp_path, capsys):
    run = tmp_path / "run"
    _pipeline(run)
    capsys.readouterr()
    # re-running genomes would leave the sequences built from it mismatched → refuse, delete nothing
    rc = main(["genomes", str(run), "--duplication", "0.3", "--loss", "0.1", "--origination", "0.4",
               "--initial-families", "5", "--seed", "9", "--quiet"])
    assert rc == 1
    err = capsys.readouterr().err
    assert "zombi2: error:" in err and "sequences" in err and "--force" in err
    assert (run / "sequences").exists()                  # a refusal is non-destructive


def test_force_reruns_the_level_and_removes_the_stale_downstream(tmp_path, capsys):
    run = tmp_path / "run"
    _pipeline(run)
    capsys.readouterr()
    rc = main(["genomes", str(run), "--duplication", "0.3", "--loss", "0.1", "--origination", "0.4",
               "--initial-families", "5", "--seed", "9", "--force", "--quiet"])
    assert rc == 0
    assert not (run / "sequences").exists()              # the now-stale downstream is cleared
    assert (run / "genomes").exists()                    # the level itself was rebuilt


def test_forward_pipeline_and_last_level_rerun_are_not_blocked(tmp_path):
    run = tmp_path / "run"
    _pipeline(run)                                        # species → genomes → sequences, all fine
    # re-running the LAST level (nothing downstream) is safe — e.g. a different substitution model
    assert main(["sequences", str(run), "--model", "k80", "--length", "20", "--seed", "2",
                 "--quiet"]) == 0


def test_flat_layout_is_left_to_the_user(tmp_path):
    run = tmp_path / "run"
    for cmd in (["species", str(run), "--birth", "1", "--death", "0.3", "--n-extant", "8", "--seed",
                 "1", "--quiet", "--flat"],
                ["genomes", str(run), "--duplication", "0.2", "--loss", "0.2", "--origination", "0.4",
                 "--initial-families", "5", "--seed", "1", "--quiet", "--flat"]):
        main(cmd)
    # --flat commingles the levels in one directory, so there is no per-level folder to guard on:
    # re-running is not blocked (documented limitation)
    assert main(["species", str(run), "--birth", "1", "--death", "0.3", "--n-extant", "8", "--seed",
                 "2", "--quiet", "--flat"]) == 0


# ── the staleness guard: DrivenBy conditioning edges ────────────────────────────────

def _conditioned_pipeline(run):
    """species -> discrete trait -> genomes with a loss *driven by* that trait -> sequences."""
    main(["species", str(run), "--birth", "1", "--death", "0.3", "--n-extant", "8", "--seed", "1",
          "--quiet"])
    main(["traits", str(run), "--kind", "discrete", "--states", "cave,surface", "--switch", "0.4",
          "--seed", "5", "--quiet"])
    main(["genomes", str(run), "--duplication", "0.2", "--origination", "0.5", "--seed", "3", "--quiet",
          "--loss", f"0.25 * DrivenBy('{run}/traits/trait_events.tsv', {{'cave': 4.0}})"])
    main(["sequences", str(run), "--model", "jc69", "--length", "20", "--seed", "1", "--quiet"])


def test_conditioning_records_which_level_drove_a_rate(tmp_path):
    run = tmp_path / "run"
    _conditioned_pipeline(run)
    assert (run / "genomes" / "conditioned_on").read_text(encoding="utf-8").split() == ["traits"]


def test_genomes_ordered_takes_a_driven_rate_from_the_command_line(tmp_path):
    """The whole conditioned ordered story on the command line: grow the driver, then drive
    ``--inversion`` by it in its written form.

    The bookkeeping is the point as much as the run. ``conditioned_on`` and the input digests used to
    be read off ``--duplication``/``--transfer``/``--loss``/``--origination``/``--transfer-to`` only,
    so a run driven through any other rate recorded nothing: re-running the driver level did not know
    it had orphaned this one, and the log did not pin the driver file it read."""
    run = tmp_path / "run"
    main(["species", str(run), "--birth", "1", "--death", "0.3", "--n-extant", "8", "--seed", "1",
          "--quiet"])
    main(["traits", str(run), "--kind", "discrete", "--states", "cave,surface", "--switch", "0.4",
          "--seed", "5", "--quiet"])
    driver = f"{run}/traits/trait_events.tsv"
    rc = main(["genomes", str(run), "--resolution", "ordered", "--initial-families", "12",
               "--inversion", f"0.3 * DrivenBy('{driver}', {{'cave': 4.0, 'surface': 1.0}})",
               "--seed", "3", "--quiet"])
    assert rc == 0
    assert (run / "genomes" / "rearrangement_events.tsv").exists()
    assert (run / "genomes" / "conditioned_on").read_text(encoding="utf-8").split() == ["traits"]
    log = (run / "genomes" / "genomes.log").read_text(encoding="utf-8")
    assert "trait_events.tsv" in log and re.search(r"\b[0-9a-f]{64}\b", log)


def test_rerunning_a_trait_a_genome_was_conditioned_on_refuses(tmp_path, capsys):
    run = tmp_path / "run"
    _conditioned_pipeline(run)
    capsys.readouterr()
    rc = main(["traits", str(run), "--kind", "discrete", "--states", "cave,surface", "--switch", "0.9",
               "--seed", "2", "--quiet"])
    assert rc == 1
    err = capsys.readouterr().err
    # the cascade: the conditioned genomes AND the sequences beneath it are both named, nothing deleted
    assert "zombi2: error:" in err and "genomes" in err and "sequences" in err
    assert (run / "genomes").exists() and (run / "sequences").exists()


def test_force_reruns_the_trait_and_clears_its_conditioned_downstream(tmp_path, capsys):
    run = tmp_path / "run"
    _conditioned_pipeline(run)
    capsys.readouterr()
    rc = main(["traits", str(run), "--kind", "discrete", "--states", "cave,surface", "--switch", "0.9",
               "--seed", "2", "--force", "--quiet"])
    assert rc == 0
    assert not (run / "genomes").exists() and not (run / "sequences").exists()
    assert (run / "traits").exists()


def test_a_trait_drives_the_substitution_rate_from_the_command_line(tmp_path, capsys):
    """The whole conditioned pipeline for Traits → Sequences, in the written form of the rate: grow a
    habitat trait, then a sequences run whose substitution rate reads it. The run has to record the
    dependency, hash the driver file it read, and refuse a re-run of the trait beneath it."""
    run = tmp_path / "run"
    main(["species", str(run), "--birth", "1", "--death", "0.3", "--n-extant", "8", "--seed", "1",
          "--quiet"])
    main(["traits", str(run), "--kind", "discrete", "--states", "cave,surface", "--switch", "0.4",
          "--seed", "5", "--quiet", "--write", "events"])
    main(["genomes", str(run), "--duplication", "0.2", "--loss", "0.2", "--origination", "0.5",
          "--seed", "3", "--quiet"])
    driver = f"{run}/traits/trait_events.tsv"
    rate = f"0.05 * DrivenBy('{driver}', {{'cave': 0.2, 'surface': 1.0}})"
    capsys.readouterr()
    assert main(["sequences", str(run), "--model", "jc69", "--length", "200", "--seed", "1",
                 "--quiet", "--substitution", rate]) == 0
    assert "driven by trait_events.tsv" in capsys.readouterr().out

    assert (run / "sequences" / "conditioned_on").read_text(encoding="utf-8").split() == ["traits"]
    log = (run / "sequences" / "sequences.log").read_text(encoding="utf-8")
    assert "DrivenBy(" in log and "trait_events.tsv" in log       # the rate in its written form
    assert any(line.startswith("input\t") and driver in line for line in log.splitlines())

    # and the trait is now a driver of this run: re-running it without --force would orphan the run
    capsys.readouterr()
    assert main(["traits", str(run), "--kind", "discrete", "--states", "cave,surface", "--switch",
                 "0.9", "--seed", "2", "--quiet"]) == 1
    assert "sequences" in capsys.readouterr().err


def test_an_unconditioned_genome_leaves_the_trait_free_to_rerun(tmp_path):
    run = tmp_path / "run"
    main(["species", str(run), "--birth", "1", "--death", "0.3", "--n-extant", "8", "--seed", "1",
          "--quiet"])
    main(["traits", str(run), "--kind", "discrete", "--states", "cave,surface", "--switch", "0.4",
          "--seed", "5", "--quiet"])
    main(["genomes", str(run), "--duplication", "0.2", "--loss", "0.2", "--origination", "0.5",
          "--seed", "3", "--quiet"])                      # NOT conditioned on the trait
    assert not (run / "genomes" / "conditioned_on").exists()
    # trait and genomes are independent here, so re-running the trait is allowed
    assert main(["traits", str(run), "--kind", "discrete", "--states", "cave,surface", "--switch",
                 "0.9", "--seed", "2", "--quiet"]) == 0


# ── onboarding nudges + no empty genomes/ dir ───────────────────────────────────────

def test_completion_nudge_points_to_the_output_files(tmp_path, capsys):
    main(["species", str(tmp_path / "run"), "--birth", "1", "--death", "0.3", "--n-extant", "6",
          "--seed", "1"])
    out = capsys.readouterr().out
    assert "species_complete.nwk" in out and "species_extant.nwk" in out   # names both trees
    assert "next:" not in out                      # never prescribes what to do next


def test_quiet_suppresses_the_nudge_but_keeps_the_wrote_line(tmp_path, capsys):
    main(["species", str(tmp_path / "run"), "--birth", "1", "--death", "0.3", "--n-extant", "6",
          "--seed", "1", "--quiet"])
    out = capsys.readouterr().out
    assert "wrote" in out
    assert "→" not in out                          # no nudge lines under --quiet


def test_non_nucleotide_sequences_run_leaves_no_empty_genomes_dir(tmp_path):
    run = tmp_path / "run"
    main(["species", str(run), "--birth", "1", "--death", "0.3", "--n-extant", "6", "--seed", "1",
          "--quiet"])
    main(["genomes", str(run), "--duplication", "0.2", "--loss", "0.2", "--origination", "0.5",
          "--initial-families", "5", "--seed", "3", "--quiet"])
    main(["sequences", str(run), "--model", "jc69", "--length", "20", "--seed", "1", "--quiet"])
    assert not (run / "sequences" / "genomes").exists()   # nucleotide-only output, empty here


def test_sequences_log_records_effective_model_params(tmp_path):
    run = tmp_path / "run"
    main(["species", str(run), "--birth", "1", "--death", "0.3", "--n-extant", "5", "--seed", "1",
          "--quiet"])
    main(["genomes", str(run), "--duplication", "0.2", "--initial-families", "3", "--seed", "3",
          "--quiet"])
    main(["sequences", str(run), "--model", "hky85", "--length", "20", "--seed", "1", "--quiet"])
    log = (run / "sequences" / "sequences.log").read_text(encoding="utf-8")
    # the resolved values the run used, not the bare `None` that was on the command line
    assert "kappa\t2.0" in log
    assert "frequencies\t[0.25, 0.25, 0.25, 0.25]" in log
    assert "exchangeabilities\tNone" in log    # a knob this model does not have stays None


# ── the family-tier knobs on the command line ───────────────────────────────────────────────────

def test_max_family_size_bounds_a_run(tmp_path, tree_file):
    """The cap is reachable from the command line, as a plain count of copies in one genome."""
    import collections
    out = tmp_path / "capped"
    main(["genomes", str(out), "--from", str(tree_file), "--duplication", "1.2", "--loss", "0.1",
          "--origination", "0.3", "--initial-families", "6", "--max-family-size", "3",
          "--seed", "1", "--flat", "--quiet"])
    header, *rows = (out / "genomes.tsv").read_text(encoding="utf-8").splitlines()
    cols = header.split("\t")
    lineage, family = cols.index("lineage"), cols.index("family")
    per_genome = collections.defaultdict(collections.Counter)
    for r in rows:
        f = r.split("\t")
        per_genome[f[lineage]][f[family]] += 1
    worst = max(max(c.values()) for c in per_genome.values())
    assert worst <= 3, f"a family reached {worst} copies in one genome, over the cap of 3"


def test_max_family_size_none_removes_the_cap(tmp_path, tree_file):
    """`none` is how you ask for unbounded growth on purpose — the flag exists because the cap is
    on by default and was previously unreachable from the CLI."""
    rc = main(["genomes", str(tmp_path / "u"), "--from", str(tree_file), "--duplication", "0.5",
               "--loss", "0.2", "--initial-families", "4", "--max-family-size", "none",
               "--seed", "1", "--flat", "--quiet"])
    assert rc == 0
    assert "max_family_size\tNone" in (tmp_path / "u" / "genomes.log").read_text(encoding="utf-8")


def test_max_family_size_rejects_a_non_number(tmp_path, tree_file):
    with pytest.raises(SystemExit) as e:
        main(["genomes", str(tmp_path / "g"), "--from", str(tree_file),
              "--max-family-size", "lots", "--flat"])
    assert e.value.code == 2


def test_family_speed_takes_a_byfamily_draw(tmp_path, tree_file):
    rc = main(["genomes", str(tmp_path / "fs"), "--from", str(tree_file), "--duplication", "0.4",
               "--loss", "0.2", "--origination", "0.3", "--initial-families", "8",
               "--family-speed", "ByFamily(spread=0.6)", "--seed", "1", "--flat", "--quiet"])
    assert rc == 0


def test_family_speed_rejects_a_bare_number(tmp_path, tree_file):
    """It is a *draw*, not a factor — a plain number would silently mean something else."""
    with pytest.raises(SystemExit) as e:
        main(["genomes", str(tmp_path / "g"), "--from", str(tree_file),
              "--family-speed", "2.0", "--flat"])
    assert e.value.code == 2


@pytest.mark.parametrize("flag, value", [("--max-family-size", "5"),
                                         ("--family-speed", "ByFamily(spread=0.5)")])
def test_the_family_knobs_are_refused_under_nucleotide(tmp_path, tree_file, flag, value):
    """Each is absent there for its own reason, and the error says which rather than giving one
    blanket explanation that fits some of the group and not the rest."""
    with pytest.raises(SystemExit) as e:
        main(["genomes", str(tmp_path / "g"), "--from", str(tree_file), "--resolution", "nucleotide",
              "--root-length", "2000", flag, value, "--flat"])
    assert e.value.code == 2


def test_an_unseeded_run_records_the_seed_it_drew(tmp_path):
    """Leaving --seed out never meant "unseeded" — it meant the seed came from the OS and was thrown
    away, so the log said `seed None` and the dataset could not be regenerated by anyone, its author
    included, with nothing to warn them. The seed is now drawn here and written down, so the log
    stays a complete record: re-running its command line reproduces the run."""
    out = tmp_path / "unseeded"
    assert main(["species", str(out), "--birth", "1", "--death", "0.3", "--n-extant", "12",
                 "--quiet"]) == 0
    logged = dict(line.split("\t", 1) for line in
                  (out / "species" / "species.log").read_text(encoding="utf-8").splitlines() if "\t" in line)
    assert logged["seed"] not in ("None", ""), "an unseeded run still records no seed"

    # and replaying that seed reproduces the run exactly
    again = tmp_path / "replay"
    assert main(["species", str(again), "--birth", "1", "--death", "0.3", "--n-extant", "12",
                 "--seed", logged["seed"], "--quiet"]) == 0
    assert ((out / "species" / "species_complete.nwk").read_text(encoding="utf-8")
            == (again / "species" / "species_complete.nwk").read_text(encoding="utf-8"))

    # a second unseeded run still draws fresh randomness — the point is that it is recorded, not fixed
    third = tmp_path / "other"
    main(["species", str(third), "--birth", "1", "--death", "0.3", "--n-extant", "12", "--quiet"])
    assert ((third / "species" / "species_complete.nwk").read_text(encoding="utf-8")
            != (out / "species" / "species_complete.nwk").read_text(encoding="utf-8"))


def test_treedist_refuses_a_repeated_tip_label(tmp_path):
    """The leaf-set check compares *sets*, so a repeated label collapsed into it and passed; the
    relabelling then mapped both copies to one id and a distance came back — a plausible number
    computed on something that is not a tree. From a scoring tool that is worse than a refusal."""
    a = tmp_path / "a.nwk"; a.write_text("((A:1,B:1):1,(C:1,C:1):1);")
    b = tmp_path / "b.nwk"; b.write_text("((A:1,C:1):1,(B:1,C:1):1);")
    with pytest.raises(SystemExit) as e:
        main(["tools", "treedist", str(a), str(b), "--metric", "rf"])
    assert e.value.code == 2

    ok_a = tmp_path / "ok_a.nwk"; ok_a.write_text("((A:1,B:1):1,(C:1,D:1):1);")
    ok_b = tmp_path / "ok_b.nwk"; ok_b.write_text("((A:1,C:1):1,(B:1,D:1):1);")
    assert main(["tools", "treedist", str(ok_a), str(ok_b), "--metric", "rf"]) == 0


# --- the small papercuts a second round of testers reported ----------------

def test_a_missing_params_file_says_so(tmp_path, capsys):
    with pytest.raises(SystemExit):
        main(["species", str(tmp_path / "o"), "--params", str(tmp_path / "nope.toml"),
              "--n-extant", "5"])
    err = capsys.readouterr().err
    assert "params file not found" in err and "Errno" not in err


def test_a_broken_params_file_names_the_file_and_the_path_trap(tmp_path, capsys):
    # TOML reads a backslash in a "..." string as an escape, so a Windows path is rejected with a
    # message about hex values that says nothing about paths
    bad = tmp_path / "p.toml"
    bad.write_text('birth = 1.0\ntransfer-to = "DrivenBy(\'C:\\Users\\x\', {})"\n', encoding="utf-8")
    with pytest.raises(SystemExit):
        main(["genomes", str(tmp_path / "o"), "--params", str(bad)])
    err = capsys.readouterr().err
    assert "p.toml is not valid TOML" in err and "literal string" in err


def test_an_unknown_flag_suggests_the_real_one(tmp_path, capsys):
    with pytest.raises(SystemExit):
        main(["genomes", str(tmp_path / "o"), "--duplication", "0.2", "--transfers", "0.1"])
    err = capsys.readouterr().err
    assert "unrecognized arguments: --transfers 0.1" in err
    assert "did you mean --transfer?" in err
    assert "usage: zombi2 genomes" in err          # the COMMAND's usage, not the top level's


def test_a_flag_unlike_anything_gets_no_invented_suggestion(tmp_path, capsys):
    with pytest.raises(SystemExit):
        main(["species", str(tmp_path / "o"), "--birth", "1", "--n-extant", "5", "--wibble"])
    assert "did you mean" not in capsys.readouterr().err


def test_the_log_records_the_options_this_resolution_has(tmp_path, tree_file):
    # a family run has no --root-length and no --inversion; logging them at their defaults reads as
    # though it had them and chose those values
    fam, nuc = tmp_path / "f", tmp_path / "n"
    main(["genomes", str(fam), "--from", str(tree_file), "--duplication", "0.2", "--seed", "1",
          "--flat", "--quiet"])
    main(["genomes", str(nuc), "--from", str(tree_file), "--resolution", "nucleotide",
          "--root-length", "200", "--seed", "1", "--flat", "--quiet"])
    fam_log = (fam / "genomes.log").read_text(encoding="utf-8")
    nuc_log = (nuc / "genomes.log").read_text(encoding="utf-8")
    keys = lambda log: {ln.split("\t")[0] for ln in log.splitlines() if "\t" in ln}

    assert not keys(fam_log) & {"root_length", "gene_length", "inversion", "chromosomes", "gff"}
    assert not keys(nuc_log) & {"initial_families", "replacement", "family_speed"}
    assert "root_length" in keys(nuc_log)          # ...and it does record what it HAS
    assert {"duplication", "seed", "resolution"} <= keys(fam_log)


def test_stream_writes_the_run_and_logs_what_it_wrote(tmp_path, tree_file, capsys):
    # --stream had no CLI test at all, and shipped broken: the log's `write` line was read from a
    # variable only the in-memory branch set, so every streamed run printed its success line and then
    # died with UnboundLocalError, exit 1, having written no genomes.log.
    out = tmp_path / "s"
    assert main(["genomes", str(out), "--from", str(tree_file), "--duplication", "0.2",
                 "--transfer", "0.1", "--loss", "0.25", "--origination", "0.5",
                 "--initial-families", "10", "--seed", "1", "--stream", "--flat", "--quiet"]) == 0
    assert "streamed to disk" in capsys.readouterr().out
    log = (out / "genomes.log").read_text(encoding="utf-8")          # the log is written at all...
    write_line = next(ln for ln in log.splitlines() if ln.startswith("write\t"))
    for name in ("events", "profiles", "genomes", "gene_trees"):     # ...and says what landed
        assert name in write_line
        assert (out / "genome_events.tsv").exists()
    assert (out / "gene_trees").is_dir()


def test_stream_logs_only_the_outputs_asked_for(tmp_path, tree_file):
    # --write narrows a streamed run exactly as it narrows an in-memory one, and the log must record
    # the narrowed set rather than the default it did not use
    out = tmp_path / "s"
    assert main(["genomes", str(out), "--from", str(tree_file), "--duplication", "0.2",
                 "--initial-families", "10", "--seed", "1", "--stream", "--write", "events",
                 "profiles", "--flat", "--quiet"]) == 0
    log = (out / "genomes.log").read_text(encoding="utf-8")
    write_line = next(ln for ln in log.splitlines() if ln.startswith("write\t"))
    assert "events" in write_line and "profiles" in write_line
    assert "gene_trees" not in write_line and not (out / "gene_trees").exists()


def test_a_run_that_samples_none_of_its_survivors_is_refused(tmp_path, capsys):
    # sampling can observe none of the survivors, and the run then has no present at all. It used to
    # return normally: `zombi2 species` wrote no species_extant.nwk and exited 0, and `zombi2 genomes`
    # then produced an empty dataset, also at exit 0 — a silent empty replicate in a batch. It is the
    # same dead end the extinction guard refuses, so it is refused the same way.
    rc = main(["species", str(tmp_path / "e"), "--birth", "1", "--death", "0.3", "--n-extant", "20",
               "--sampling", "0.2", "--seed", "16", "--quiet"])
    assert rc == 1
    err = capsys.readouterr().err
    assert "zombi2: error:" in err and "sampling" in err
    assert "survivors" in err and "seed" in err            # says what happened and what to do
    assert not (tmp_path / "e").exists()                   # and leaves nothing behind


def test_the_summary_accounts_for_every_tip_it_counts(tmp_path, capsys):
    # "0 extant + 8 extinct (28 tips)" is arithmetic that does not work: under --sampling the
    # unsampled survivors are tips too, and leaving them out of the parts reads as a bug
    assert main(["species", str(tmp_path / "s"), "--birth", "1", "--death", "0.3", "--n-extant", "20",
                 "--sampling", "0.5", "--seed", "3", "--quiet"]) == 0
    line = capsys.readouterr().out
    counts = {word: int(n) for n, word in
              re.findall(r"(\d+) (extant|extinct|unsampled|fossils)", line)}
    n_tips = int(re.search(r"\((\d+) tips", line).group(1))
    assert counts.get("unsampled", 0) > 0                  # the seed leaves some unobserved...
    assert sum(counts.values()) == n_tips                  # ...and now the parts add up


def test_the_old_scope_spelling_says_what_the_number_used_to_mean(tmp_path, tree_file, capsys):
    # PerLineage(10) meant 10 x the node count of the species tree, so a stale command line ran at a
    # cap it did not read. It fails now, and the message does the arithmetic rather than just saying
    # the spelling changed.
    with pytest.raises(SystemExit) as e:
        main(["genomes", str(tmp_path / "g"), "--from", str(tree_file), "--duplication", "0.2",
              "--max-family-size", "PerLineage(10)", "--flat"])
    assert e.value.code == 2
    err = " ".join(capsys.readouterr().err.split())
    assert "--max-family-size 10" in err and "nodes in the species tree" in err

    with pytest.raises(SystemExit):                     # Global(N) was exactly N — say so
        main(["genomes", str(tmp_path / "h"), "--from", str(tree_file), "--duplication", "0.2",
              "--max-family-size", "Global(10)", "--flat"])
    assert "--max-family-size 10" in " ".join(capsys.readouterr().err.split())


def test_the_cap_can_be_set_from_a_params_file(tmp_path, tree_file):
    # a TOML integer reaches the engine directly (params seed argparse defaults, bypassing `type`),
    # so the plain-number form has to be what the engine accepts — under the scope spelling this
    # route could not set the cap at all
    params = tmp_path / "p.toml"
    params.write_text("[genomes]\nmax-family-size = 3\n", encoding="utf-8")
    out = tmp_path / "g"
    assert main(["genomes", str(out), "--from", str(tree_file), "--duplication", "1.2",
                 "--loss", "0.1", "--initial-families", "6", "--params", str(params),
                 "--seed", "1", "--flat", "--quiet"]) == 0
    assert "max_family_size\t3" in (out / "genomes.log").read_text(encoding="utf-8")
    counts = collections.Counter()
    header, *rows = (out / "genomes.tsv").read_text(encoding="utf-8").splitlines()
    cols = header.split("\t")
    for row in rows:
        cell = dict(zip(cols, row.split("\t")))
        counts[(cell["lineage"], cell["family"])] += 1
    assert counts and max(counts.values()) == 3




def test_a_run_written_from_python_carries_its_species_tree(tmp_path):
    # the README's Python route used to write gene trees and no species tree, so a beginner got a
    # dataset with no ground truth in it and nothing said so
    sp = simulate_species_tree(birth=1.0, death=0.3, n_extant=12, seed=1)
    g = simulate_genomes_family(sp, duplication=0.2, transfer=0.1, initial_families=5, seed=2)
    g.write(tmp_path / "py")
    written = tmp_path / "py" / "species_complete.nwk"
    assert written.exists()
    assert written.read_text(encoding="utf-8") == sp.complete_tree.to_newick() + "\n"


def test_the_cli_keeps_one_species_tree_not_one_per_level(tmp_path):
    # the Result writes it so a Python directory stands alone; a run already has a canonical copy
    # under species/, and two files for one fact is how they drift apart
    run = tmp_path / "r"
    assert main(["species", str(run), "--birth", "1", "--n-extant", "8", "--seed", "1",
                 "--quiet"]) == 0
    assert main(["genomes", str(run), "--duplication", "0.2", "--initial-families", "5",
                 "--seed", "1", "--quiet"]) == 0
    assert (run / "species" / "species_complete.nwk").exists()
    assert not (run / "genomes" / "species_complete.nwk").exists()
    log = (run / "genomes" / "genomes.log").read_text(encoding="utf-8")
    # the `write` line, not the whole log: the log also records the command line that invoked this
    write_line = next(ln for ln in log.splitlines() if ln.startswith("write\t"))
    assert "species_tree" not in write_line        # and it says what it actually wrote
    assert "gene_trees" in write_line


def test_force_says_what_it_deleted_even_under_quiet(tmp_path, capsys):
    # --force is mandatory in any pipeline that re-runs a level, so its notice was permanently
    # suppressed by the --quiet those pipelines also pass: 294 files could vanish silently
    run = tmp_path / "r"
    for cmd in (["species", str(run), "--birth", "1", "--n-extant", "8", "--seed", "1"],
                ["genomes", str(run), "--duplication", "0.2", "--initial-families", "5",
                 "--seed", "1"]):
        assert main(cmd + ["--quiet"]) == 0
    assert (run / "genomes").exists()
    capsys.readouterr()
    assert main(["species", str(run), "--birth", "1", "--n-extant", "8", "--seed", "2",
                 "--quiet", "--force"]) == 0
    captured = capsys.readouterr()
    assert "genomes" in captured.err and "--force removed" in captured.err
    assert "--force removed" not in captured.out           # a diagnostic, so stderr not stdout
    assert not (run / "genomes").exists()                  # and it really did remove it


def test_a_bare_run_announces_the_demo_values_it_invented(tmp_path, capsys):
    # a bare run fills the rates it was not given with round demonstration values. It no longer warns
    # on stderr (a paragraph on every bare run was noise); the invented values are recorded instead — in
    # the log and run.zombi2 — visible where they are useful rather than shouted where they are not.
    assert main(["species", str(tmp_path / "bare"), "--seed", "1", "--quiet"]) == 0
    assert capsys.readouterr().err == ""                     # silent now
    log = (tmp_path / "bare" / "species" / "species.log").read_text(encoding="utf-8")
    assert "birth\t1.0" in log and "n_extant\t20" in log     # but the invented values are recorded

    # supplying them yourself is likewise silent
    assert main(["species", str(tmp_path / "given"), "--birth", "1", "--n-extant", "8",
                 "--seed", "1", "--quiet"]) == 0
    assert capsys.readouterr().err == ""


def test_stream_is_refused_on_a_nucleotide_sequences_run(tmp_path, tree_file, capsys):
    # a nucleotide run reassembles every node's genome, which needs every block's sequence at once —
    # the opposite of keeping nothing. Refused up front rather than half-run.
    run = tmp_path / "n"
    assert main(["genomes", str(run), "--from", str(tree_file), "--resolution", "nucleotide",
                 "--root-length", "400", "--genes", "3", "--gene-length", "100", "--seed", "1", "--quiet"]) == 0
    with pytest.raises(SystemExit) as e:
        main(["sequences", str(run), "--model", "jc69", "--seed", "1", "--quiet", "--stream"])
    assert e.value.code == 2
    err = capsys.readouterr().err
    assert "--stream" in err and "nucleotide" in err and "--resolution family" in err

def test_write_offers_every_output_the_result_can_write():
    """The CLI's --write choices are the result's own write vocabulary, not a copy of it.

    They were a copy, with a comment saying so, and they drifted: `species_tree` and
    `initial_sequence` were writable from Python and unnameable on the command line, and `summary`
    was written by default at three levels and listed at none."""
    from zombi2.cli.genomes import _OUTPUTS
    from zombi2.genomes.family import FamilyGenomesResult
    from zombi2.genomes.nucleotide import NucleotideGenomesResult
    from zombi2.genomes.ordered import OrderedGenomesResult

    for name, result in (("family", FamilyGenomesResult), ("ordered", OrderedGenomesResult),
                         ("nucleotide", NucleotideGenomesResult)):
        assert tuple(_OUTPUTS[name]) == tuple(result.OUTPUTS), name
        # and the declared vocabulary is really what write() honours: every default is in it
        import inspect
        default = inspect.signature(result.write).parameters["outputs"].default
        assert set(default) <= set(result.OUTPUTS), f"{name}: a default output is not declared"

def test_named_traits_sit_side_by_side_and_one_can_drive_the_other(tmp_path):
    """A run directory holds one slot per level, which is wrong for traits: a tree can carry several,
    and one can drive another. --name gives each its own, so the driver's files survive the target."""
    run = tmp_path / "run"
    assert main(["species", str(run), "--birth", "1", "--death", "0.2", "--n-extant", "12",
                 "--seed", "1", "--quiet"]) == 0
    assert main(["traits", str(run), "--name", "habitat", "--kind", "discrete", "--states",
                 "slow,fast", "--switch", "0.5", "--seed", "1", "--quiet"]) == 0
    driver = run / "traits" / "habitat" / "trait_events.tsv"
    assert driver.is_file()
    before = driver.read_text(encoding="utf-8")

    assert main(["traits", str(run), "--name", "size", "--kind", "continuous", "--rate",
                 f"1.0 * DrivenBy('{driver}', {{'fast': 6.0}})", "--seed", "2", "--quiet"]) == 0
    # the driver is untouched, and the target is beside it rather than on top of it
    assert driver.read_text(encoding="utf-8") == before
    assert (run / "traits" / "size" / "trait_values.tsv").is_file()
    # the run report names each one, so a reader can tell them apart
    report = (run / "run.zombi2").read_text(encoding="utf-8")
    assert "TRAITS (habitat)" in report and "TRAITS (size)" in report


def test_name_and_flat_ask_for_opposite_things(tmp_path):
    run = tmp_path / "run"
    main(["species", str(run), "--birth", "1", "--n-extant", "6", "--seed", "1", "--quiet"])
    with pytest.raises(SystemExit):
        main(["traits", str(run), "--name", "a", "--flat", "--kind", "continuous", "--rate", "1.0",
              "--seed", "1", "--quiet"])


def test_traits_records_which_level_drove_its_rate(tmp_path):
    """`zombi2 traits` wrote no `conditioned_on` marker, so a trait driving another trait recorded no
    dependency at all — `genomes` and `sequences` both record theirs.

    The within-level edge is deliberately not followed by the staleness graph: both traits live under
    `traits/`, so the edge is traits -> traits, and a graph whose nodes are levels would mark the
    level stale against itself — `--force` would then delete every trait in the run, the one just
    written included. The marker is recorded (it is what the log and the run report show); automatic
    invalidation waits for the graph's nodes to be `traits/<name>`."""
    run = tmp_path / "r"
    assert main(["species", str(run), "--birth", "1.0", "--n-extant", "8", "--seed", "1", "--quiet"]) == 0
    assert main(["traits", str(run), "--name", "habitat", "--kind", "discrete", "--states",
                 "cave,surface", "--switch", "0.5", "--seed", "1", "--quiet"]) == 0
    driver = run / "traits" / "habitat" / "trait_events.tsv"
    assert main(["traits", str(run), "--name", "size", "--kind", "continuous", "--rate",
                 f"1.0 * DrivenBy('{driver}', {{'cave': 0.2, 'surface': 1.0}})",
                 "--seed", "1", "--quiet"]) == 0

    assert (run / "traits" / "size" / "conditioned_on").read_text().split() == ["traits"]
    assert not (run / "traits" / "habitat" / "conditioned_on").exists()   # the driver drives nothing

    # re-running the driver with --force must not take the trait beside it with it
    assert main(["traits", str(run), "--name", "habitat", "--kind", "discrete", "--states",
                 "cave,surface", "--switch", "0.9", "--seed", "7", "--quiet", "--force"]) == 0
    assert (run / "traits" / "size" / "trait_values.tsv").exists()


def test_fasta_is_refused_on_a_family_run(tmp_path):
    """Regression. ``--fasta`` was missing from the nucleotide-only list, so a family run accepted
    it, ignored it, and wrote it into genomes.log as though it had applied — the one silent-ignore
    the resolution guard exists to prevent, and its own sibling ``--gff`` was already caught."""
    main(["species", str(tmp_path), "--birth", "1", "--death", "0.3", "--n-extant", "12",
          "--seed", "1"])
    fasta = tmp_path / "t.fasta"
    fasta.write_text(">chr\nACGTACGTAC\n", encoding="utf-8")
    with pytest.raises(SystemExit) as e:
        main(["genomes", str(tmp_path), "--duplication", "0.2", "--loss", "0.25",
              "--origination", "0.5", "--fasta", str(fasta), "--seed", "1"])
    assert e.value.code == 2
