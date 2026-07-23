"""Tests for the ``zombi2`` command line (species + genomes) and the ``read_newick`` reader it
loads trees with. The CLI is a thin shell over the level ``simulate_*`` functions, so these tests
exercise argument wiring, output files, ``--params``, and the clean error paths — not the engines
(which their own suites cover)."""

import pytest

from zombi2.cli.main import main
from zombi2.genomes import simulate_genomes_unordered
from zombi2.sequences.substitution_models import AMINO_ACIDS
from zombi2.species import read_newick, simulate_species_tree


# ── read_newick ─────────────────────────────────────────────────────────────────────

def test_read_newick_round_trips_a_complete_tree():
    # a run with survivors round-trips: ids, node count, and the extant/extinct split are preserved
    r = simulate_species_tree(birth=1.0, death=0.4, n_extant=40, seed=3)
    t = r.complete_tree
    back, names = read_newick(t.to_newick())
    assert names == {}                                           # a ZOMBI tree: labels ARE the ids
    assert set(back.nodes) == set(t.nodes)
    assert back.root == t.root
    assert len(back.extant()) == len(t.extant())
    assert len(back.extinct()) == len(t.extinct())
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
    survivors = [n.id for n in r.complete_tree.extant()]
    fates = {f"n{i}": "extant" for i in survivors}
    fates[f"n{survivors[0]}"] = "unsampled"                       # override one survivor
    for n in r.complete_tree.extinct():
        fates[f"n{n.id}"] = "extinct"
    back, _ = read_newick(r.complete_tree.to_newick(), tip_fates=fates)
    assert back.nodes[survivors[0]].fate == "unsampled"          # the table won, not the depth-guess
    assert len(back.extant()) == len(survivors) - 1              # the unsampled one is no longer extant
    assert len(back.unsampled()) == 1


def test_read_newick_ultrametric_external_tree_is_all_extant_with_a_name_map():
    # ultrametric (every tip at depth 2) → every tip extant, and the user's labels come back mapped
    t, names = read_newick("((human:1,chimp:1):1,(mouse:0.8,rat:0.8):1.2);")
    assert len(t.extant()) == 4 and not t.extinct()
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
    assert parsed == {f"n{n.id}": n.fate for n in r.complete_tree.leaves()}


def test_genomes_reads_the_runs_fate_table_so_unsampled_tips_are_not_extant(tmp_path, capsys):
    # the sampling handoff bug: without the fate table, genomes reads the tree by depth and counts every
    # survivor extant (extant + unsampled); consuming species_fates.tsv, it builds only the sampled ones
    run = tmp_path / "run"
    main(["species", str(run), "--birth", "1.0", "--death", "0.4", "--n-extant", "20",
          "--sampling", "0.5", "--seed", "7"])
    fates = (run / "species" / "species_fates.tsv").read_text().splitlines()[1:]
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
    g = simulate_genomes_unordered(back, duplication=0.2, loss=0.2, origination=0.5, seed=9)
    assert g.profiles.shape[1] == len(back.extant())            # one column per extant tip


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
    log = (tmp_path / "species.log").read_text()
    assert "birth\t1.0" in log and "n_extant\t20" in log


def test_species_write_is_selective(tmp_path):
    main(["species", str(tmp_path), "--birth", "1", "--total-time", "3", "--seed", "1", "--write", "complete", "--flat"])
    nwk = {p.name for p in tmp_path.iterdir() if p.suffix == ".nwk"}
    assert nwk == {"species_complete.nwk"}                       # no extant/events file


def test_species_is_deterministic_given_the_seed(tmp_path):
    a, b = tmp_path / "a", tmp_path / "b"
    for out in (a, b):
        main(["species", str(out), "--birth", "1", "--death", "0.3", "--n-extant", "30", "--seed", "13", "--flat"])
    assert (a / "species_complete.nwk").read_text() == (b / "species_complete.nwk").read_text()


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
                   "c1\t.\tgene\t300\t600\t.\t+\t.\tID=b\n")
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
                                               "species_complete.nwk", "species_fates.tsv", "genomes.log"}


def test_genomes_ordered_events_carry_where_each_one_happened(tmp_path, tree_file):
    # the genealogy, the positions and the rearrangements are one table now, not three
    out = tmp_path / "g"
    rc = main(["genomes", str(out), "--from", str(tree_file), "--resolution", "ordered", "--duplication", "0.3", "--loss", "0.2", "--origination", "0.6", "--inversion", "0.4", "--seed", "42", "--write", "events", "--flat"])
    assert rc == 0
    assert not list(out.glob("rearrangements.tsv")) and not list(out.glob("genome_event_*.tsv"))
    rows = (out / "genome_events.tsv").read_text().splitlines()
    cols = rows[0].split("\t")
    assert cols == ["time", "kind", "lineage", "family", "copy", "parent", "recipient",
                    "donor", "dest_lineage", "chromosome", "position", "length",
                    "dest_chromosome", "dest_position", "flipped"]
    body = [r.split("\t") for r in rows[1:]]
    assert {r[1] for r in body} >= {"origination", "duplication", "loss", "inversion"}
    assert [r for r in body if r[10]], "no event carries the place it happened"


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
    rc = main(["genomes", str(out), "--from", str(tree_file), "--resolution", "nucleotide", "--root-length", "600", "--genes", "3", "--inversion", "1.0", "--duplication", "0.5", "--loss", "0.4", "--seed", "1", "--flat"])
    assert rc == 0
    # the nucleotide default is events + genes; blocks is opt-in
    written = {p.name for p in out.iterdir()}
    assert {"genome_events.tsv", "genes.tsv", "species_complete.nwk", "genomes.log"} <= written
    assert len((out / "genes.tsv").read_text().splitlines()) > 1


def test_genomes_nucleotide_write_selects_blocks(tmp_path, tree_file):
    out = tmp_path / "g"
    rc = main(["genomes", str(out), "--from", str(tree_file), "--resolution", "nucleotide", "--root-length", "400", "--genes", "2", "--inversion", "1.0", "--seed", "1", "--write", "blocks", "events", "--flat"])
    assert rc == 0
    head = (out / "blocks.tsv").read_text().splitlines()[0].split("\t")
    assert head == ["lineage", "chromosome", "position", "source", "start", "end", "strand",
                    "copy", "gene"]


def test_genomes_nucleotide_seeds_from_a_gff(tmp_path, tree_file):
    gff = tmp_path / "seed.gff"
    gff.write_text("##gff-version 3\n"
                   "##sequence-region chrom1 1 900\n"
                   "chrom1\tZOMBI2\tgene\t101\t200\t.\t+\t.\tID=dnaA\n"
                   "chrom1\tZOMBI2\tgene\t401\t500\t.\t-\t.\tID=recA\n")
    out = tmp_path / "g"
    rc = main(["genomes", str(out), "--from", str(tree_file), "--resolution", "nucleotide", "--gff", str(gff), "--inversion", "1.0", "--seed", "1", "--flat"])
    assert rc == 0
    rows = [r.split("\t") for r in (out / "genes.tsv").read_text().splitlines()[1:]]
    assert [r[1] for r in rows] == ["dnaA", "recA"]          # names survive to the output
    assert [r[5] for r in rows] == ["1", "-1"]               # ...and so does the coding strand


def test_genomes_is_deterministic_across_resolutions(tmp_path, tree_file):
    def run(tag):
        out = tmp_path / tag
        main(["genomes", str(out), "--from", str(tree_file), "--resolution", "nucleotide", "--root-length", "500", "--genes", "2", "--inversion", "1.0", "--duplication", "0.4", "--seed", "7", "--write", "events", "blocks", "--flat"])
        return {p.name: p.read_text() for p in out.iterdir() if p.suffix == ".tsv"}
    assert run("a") == run("b")


@pytest.mark.parametrize("argv, why", [
    (["--initial-families", "5"], "nucleotide has no initial-families"),
    (["--replacement"], "nucleotide transfers are additive"),
    (["--loss", "0.2 * OnTime({0: 1.0, 3: 2.0})"], "nucleotide takes constant rates only"),
    (["--gff", "x.gff", "--genes", "3"], "gff and genes are mutually exclusive"),
    (["--write", "gene_order"], "gene_order is an ordered output"),
    (["--write", "profiles"], "the nucleotide resolution has no profiles"),
])
def test_genomes_nucleotide_rejects_foreign_options(tmp_path, tree_file, argv, why):
    with pytest.raises(SystemExit) as e:
        main(["genomes", str(tmp_path / "g"), "--from", str(tree_file), "--resolution", "nucleotide", *argv, "--flat"])
    assert e.value.code == 2, why


@pytest.mark.parametrize("resolution", ["unordered", "ordered"])
def test_genomes_rejects_nucleotide_only_flags_elsewhere(tmp_path, tree_file, resolution):
    with pytest.raises(SystemExit) as e:                     # bp knobs need a nucleotide genome
        main(["genomes", str(tmp_path / "g"), "--from", str(tree_file), "--resolution", resolution, "--root-length", "500", "--flat"])
    assert e.value.code == 2


def test_sequences_reads_a_nucleotide_handoff_through_from(tmp_path, tree_file):
    # --from a nucleotide run works like any other: the handoff says which resolution wrote it (only
    # that one writes blocks.tsv), so nothing has to be repeated from the genomes command
    out = tmp_path / "g"
    main(["genomes", str(out), "--from", str(tree_file), "--resolution", "nucleotide", "--root-length", "400", "--genes", "2", "--seed", "1", "--flat"])
    s = tmp_path / "s"
    rc = main(["sequences", str(s), "--from", str(out), "--model", "jc69", "--seed", "1", "--flat"])
    assert rc == 0
    assert list(s.glob("genome_n*.fasta"))            # assembled genomes, which only this level has
    assert list(s.glob("block*.fasta"))               # blocks, not families: the files say so


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
    assert len((out / "profiles.tsv").read_text().splitlines()[0].split("\t")) == 1 + 4
    # names.tsv maps ZOMBI's n<id> back to the user's labels
    mapped = dict(row.split("\t") for row in (out / "names.tsv").read_text().splitlines()[1:])
    assert sorted(mapped.values()) == ["chimp", "human", "mouse", "rat"]


def test_genomes_on_nonultrametric_tree_needs_tip_fates(tmp_path, capsys):
    (tmp_path / "ext.nwk").write_text("((a:1,b:1):1,c:1.5);\n")     # c ends early
    rc = main(["genomes", str(tmp_path / "g"), "--from", str(tmp_path / "ext.nwk"), "--duplication", "0.3", "--seed", "1", "--flat"])
    assert rc == 1
    assert "not ultrametric" in capsys.readouterr().err


def test_genomes_nonultrametric_tree_runs_with_tip_fates_file(tmp_path):
    (tmp_path / "ext.nwk").write_text("((a:1,b:1):1,c:1.5);\n")
    (tmp_path / "fates.tsv").write_text("a\textant\nb\textant\nc\textinct\n")
    out = tmp_path / "g"
    rc = main(["genomes", str(out), "--from", str(tmp_path / "ext.nwk"), "--tip-fates", str(tmp_path / "fates.tsv"), "--duplication", "0.3", "--origination", "1.0", "--seed", "1", "--flat"])
    assert rc == 0
    # c is extinct, so only a and b are observed → two profile columns
    assert len((out / "profiles.tsv").read_text().splitlines()[0].split("\t")) == 1 + 2


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
    log = (out / "sequences.log").read_text()
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
    residues = set("".join(ln for ln in fasta.read_text().splitlines() if not ln.startswith(">")))
    assert residues <= set(AMINO_ACIDS) and not residues <= set("ACGT")


@pytest.mark.parametrize("flag,value", [("--kappa", ["2"]), ("--frequencies", ["0.25"] * 4),
                                        ("--gtr-rates", ["1"] * 6)])
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
    return {p.name: p.read_text() for p in d.iterdir() if p.suffix in (".fasta", ".nwk")}


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
    assert {p.name for p in out.iterdir()} == {"trait_values.tsv", "trait_tree.nwk", "traits.log"}
    header, first = (out / "trait_values.tsv").read_text().splitlines()[:2]
    assert header == "node\ttrait"
    assert first.split("\t")[0].startswith("n")          # n<id>, matching the Newick
    float(first.split("\t")[1])                          # a continuous trait is a number


def test_traits_ou_and_threshold_run(tmp_path, tree_file):
    # OU: the same diffusion pulled to an optimum
    assert main(["traits", "--kind", "continuous", str(tmp_path / "ou"), "--from", str(tree_file), "--rate", "1.0", "--reverts-to", "2", "--pull", "0.5", "--seed", "1", "--flat"]) == 0
    # the threshold model: a discrete state read off a continuous liability
    out = tmp_path / "th"
    assert main(["traits", str(out), "--from", str(tree_file), "--kind", "discrete", "--states", "absent,present", "--liability", "1.0", "--threshold", "0.0", "--seed", "1", "--flat"]) == 0
    states = {ln.split("\t")[1] for ln in (out / "trait_values.tsv").read_text().splitlines()[1:]}
    assert states <= {"absent", "present"}


def test_traits_discrete_writes_the_event_log(tmp_path, tree_file):
    out = tmp_path / "t"
    rc = main(["traits", str(out), "--from", str(tree_file), "--kind", "discrete", "--states", "marine,terrestrial", "--switch", "0.3", "--start", "marine", "--seed", "1", "--write", "values", "events", "tree", "--flat"])
    assert rc == 0
    assert {p.name for p in out.iterdir()} == {"trait_values.tsv", "trait_events.tsv",
                                               "trait_tree.nwk", "traits.log"}
    lines = (out / "trait_events.tsv").read_text().splitlines()
    assert lines[0] == "time\tkind\tlineage\tfrom\tto"
    assert lines[1].split("\t")[1] == "root"          # the origin row, the conditioning file's anchor


def test_traits_at_speciation_logs_on_speciation_changes(tmp_path, tree_file):
    out = tmp_path / "t"
    main(["traits", "--kind", "continuous", str(out), "--from", str(tree_file), "--rate", "1.0", "--at-speciation", "0.5", "--seed", "1", "--write", "events", "--flat"])
    rows = (out / "trait_events.tsv").read_text().splitlines()[1:]
    kinds = [r.split("\t")[1] for r in rows]
    assert kinds[0] == "root"                                                 # the origin, then jumps
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
        written.append((out / "trait_values.tsv").read_text())
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
    names = dict(ln.split("\t") for ln in (out / "names.tsv").read_text().splitlines()[1:])
    assert sorted(names.values()) == ["chimp", "human", "mouse", "rat"]
    # the name map joins the values table on its node column
    nodes = {ln.split("\t")[0] for ln in (out / "trait_values.tsv").read_text().splitlines()[1:]}
    assert nodes <= set(names)


def test_traits_params_file_drives_the_run_and_cli_overrides(tmp_path, tree_file):
    # a [traits] table scopes one file to this command (so one file can drive a whole pipeline)
    (tmp_path / "p.toml").write_text('[traits]\nkind = "discrete"\n'
                                     'states = "marine,terrestrial"\nswitch = 0.15\n'
                                     'write = ["values", "events"]\nseed = 7\n')
    argv = ["traits", "--params", str(tmp_path / "p.toml"), "--from", str(tree_file)]
    out = tmp_path / "a"
    # --kind is required, but the file supplies it — that is why it is validated in run() rather
    # than marked argparse-`required`, which no default could satisfy
    assert main(["traits", str(out), *argv[1:], "--flat"]) == 0
    assert {p.name for p in out.iterdir()} == {"trait_values.tsv", "trait_events.tsv", "traits.log"}
    states = {ln.split("\t")[1] for ln in (out / "trait_values.tsv").read_text().splitlines()[1:]}
    assert states <= {"marine", "terrestrial"}          # the file's states reached the engine

    # a flag given on the command line still wins over the file
    other = tmp_path / "b"
    assert main([*argv, str(other), "--seed", "8", "--flat"]) == 0
    assert (other / "trait_values.tsv").read_text() != (out / "trait_values.tsv").read_text()


# ── --params ────────────────────────────────────────────────────────────────────────

def test_params_file_supplies_defaults_and_cli_overrides(tmp_path):
    (tmp_path / "p.toml").write_text("birth = 2.0\ndeath = 0.3\nn-extant = 12\n")
    out = tmp_path / "o"
    # birth comes from the file; the command line still overrides it
    main(["species", str(out), "--params", str(tmp_path / "p.toml"), "--birth", "1.0", "--seed", "1", "--flat"])
    log = (out / "species.log").read_text()
    assert "birth\t1.0" in log and "n_extant\t12" in log


def test_params_file_scopes_by_command_table(tmp_path):
    (tmp_path / "pipeline.toml").write_text(
        "[species]\nbirth = 1.0\nn-extant = 15\n\n[genomes]\nduplication = 0.2\nwrite = "
        '["profiles"]\n')
    sp, gn = tmp_path / "sp", tmp_path / "gn"
    main(["species", str(sp), "--params", str(tmp_path / "pipeline.toml"), "--seed", "1", "--flat"])
    main(["genomes", str(gn), "--params", str(tmp_path / "pipeline.toml"), "--from", str(sp / "species_complete.nwk"), "--seed", "1", "--flat"])
    assert {p.name for p in gn.iterdir()} == {"profiles.tsv",
                                              "species_complete.nwk", "species_fates.tsv", "genomes.log"}


def test_params_unknown_key_errors(tmp_path):
    (tmp_path / "bad.toml").write_text("birth = 1.0\nbogus = 3\n")
    with pytest.raises(SystemExit) as e:
        main(["species", str(tmp_path / "o"), "--params", str(tmp_path / "bad.toml"), "--n-extant", "5", "--flat"])
    assert e.value.code == 2


def test_params_mistyped_section_errors(tmp_path, capsys):
    # a typo'd [table] must not be silently dropped (which would run every rate at its default)
    (tmp_path / "typo.toml").write_text("[speces]\nbirth = 1.0\nn-extant = 10\n")
    with pytest.raises(SystemExit):
        main(["species", str(tmp_path / "o"), "--params", str(tmp_path / "typo.toml"), "--seed", "1", "--flat"])
    err = capsys.readouterr().err
    assert "unknown section" in err and "[speces]" in err


def test_params_top_level_key_broadcasts_under_a_table(tmp_path):
    # a top-level scalar is a shared base for every command; a [command] table overrides on conflict
    (tmp_path / "p.toml").write_text("seed = 99\n[species]\nbirth = 1.0\nn-extant = 8\n")
    out = tmp_path / "o"
    main(["species", str(out), "--params", str(tmp_path / "p.toml"), "--flat"])
    assert "seed\t99" in (out / "species.log").read_text()      # the shared seed was applied, not dropped


def test_params_command_table_overrides_the_shared_base(tmp_path):
    # on conflict the [command] table wins over a top-level key of the same name
    (tmp_path / "p.toml").write_text("seed = 1\nn-extant = 5\n[species]\nbirth = 1.0\nn-extant = 12\n")
    out = tmp_path / "o"
    main(["species", str(out), "--params", str(tmp_path / "p.toml"), "--flat"])
    log = (out / "species.log").read_text()
    assert "n_extant\t12" in log and "seed\t1" in log


def test_params_append_option_is_overridden_by_the_command_line(tmp_path):
    # --mass-extinction is an 'append' action: a params default plus a CLI flag must NOT concatenate
    (tmp_path / "p.toml").write_text(
        "[species]\nbirth = 1.0\ntotal-time = 5.0\nmass-extinction = [[2.0, 0.9]]\n")
    out = tmp_path / "o"
    main(["species", str(out), "--params", str(tmp_path / "p.toml"),
          "--mass-extinction", "3.0", "0.1", "--seed", "3", "--flat"])
    log = (out / "species.log").read_text()
    assert "mass_extinction\t[[3.0, 0.1]]" in log              # only the command line's pulse, not both


def test_params_last_of_two_files_wins(tmp_path):
    # two --params: the last file's values are used (and it is the one the log names)
    (tmp_path / "a.toml").write_text("[species]\nbirth = 1.0\nn-extant = 5\n")
    (tmp_path / "b.toml").write_text("[species]\nbirth = 1.0\nn-extant = 40\n")
    out = tmp_path / "o"
    main(["species", str(out), "--params", str(tmp_path / "a.toml"),
          "--params", str(tmp_path / "b.toml"), "--seed", "1", "--flat"])
    assert "n_extant\t40" in (out / "species.log").read_text()


def test_params_invalid_choice_errors_cleanly(tmp_path, genomes_dir):
    # a bad choices= value in --params must be a clean error, not a KeyError deep in the command
    (tmp_path / "bad.toml").write_text('[genomes]\nresolution = "unordred"\nduplication = 0.1\n')
    with pytest.raises(SystemExit):
        main(["genomes", str(tmp_path / "o"), "--from", str(genomes_dir / "species_complete.nwk"),
              "--params", str(tmp_path / "bad.toml"), "--seed", "1", "--flat"])


def test_params_can_supply_the_sequences_model(tmp_path, genomes_dir):
    # --model is not argparse-required, so a --params file can supply it (like --birth on species)
    (tmp_path / "seq.toml").write_text('[sequences]\nmodel = "hky85"\nlength = 200\n')
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
    n = {d: len(read_newick((d / "species_complete.nwk").read_text())[0].nodes)
         for d in (flat, skyline)}
    assert n[skyline] < n[flat]


def test_species_takes_a_scope_wrapper(tmp_path):
    # Global(base) is one budget for the whole tree: linear growth, so far fewer lineages
    out = tmp_path / "g"
    rc = main(["species", str(out), "--birth", "Global(2.0)", "--total-time", "5", "--seed", "1", "--flat"])
    assert rc == 0
    assert "birth\tGlobal(2.0)" in (out / "species.log").read_text()


def test_species_records_the_rate_in_its_written_form(tmp_path):
    # the log line is the flag value again — a reproducibility record you can paste back
    out = tmp_path / "o"
    main(["species", str(out), "--birth", "1.0 * OnTime({0: 1.0, 3: 0.3})", "--total-time", "4", "--seed", "1", "--flat"])
    assert "birth\t1.0 * OnTime({0: 1, 3: 0.3})" in (out / "species.log").read_text()


def test_species_refuses_a_modifier_it_does_not_wire(tmp_path, capsys):
    # ByLineage would return a factor of 1.0 at this level — a run quietly not the model asked for
    rc = main(["species", str(tmp_path / "o"), "--birth", "1.0 * ByLineage(spread=0.3)", "--total-time", "3", "--seed", "1", "--flat"])
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
    assert "loss\t0.25 * OnTime({0: 1, 2: 3})" in (out / "genomes.log").read_text()


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
    assert f"transfer\t0.2 * DrivenBy('{driver_file}', Table({{'competent': 4, 'normal': 1}}))" \
        in (out / "genomes.log").read_text()


def test_genomes_transfer_to_takes_a_driven_recipient_weight(tmp_path, driver_file, tree_file):
    # the RECIPIENT side: the choice slot, so the modifier is written on its own, with no base
    out = tmp_path / "g"
    rc = main(["genomes", str(out), "--from", str(tree_file), "--initial-families", "5", "--transfer", "0.5", "--transfer-to", f"DrivenBy('{driver_file}', {{'competent': 2.0, 'normal': 1.0}})", "--seed", "2", "--flat"])
    assert rc == 0
    assert f"transfer_to\tDrivenBy('{driver_file}', Table({{'competent': 2, 'normal': 1}}))" \
        in (out / "genomes.log").read_text()


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
    (tmp_path / "p.toml").write_text(
        f'transfer = 0.5\ninitial-families = 5\n'
        f'transfer-to = "DrivenBy(\'{driver_file}\', {{\'competent\': 2.0, \'normal\': 1.0}})"\n')
    out = tmp_path / "g"
    rc = main(["genomes", str(out), "--params", str(tmp_path / "p.toml"), "--from", str(tree_file), "--seed", "2", "--flat"])
    assert rc == 0
    assert "transfer_to\tDrivenBy(" in (out / "genomes.log").read_text()


def test_traits_takes_a_rate_expression(tmp_path, tree_file):
    out = tmp_path / "t"
    rc = main(["traits", "--kind", "continuous", str(out), "--from", str(tree_file), "--rate", "1.0 * FromParent(spread=0.2)", "--seed", "1", "--flat"])
    assert rc == 0
    assert "rate\t1.0 * FromParent(spread=0.2)" in (out / "traits.log").read_text()


def test_params_file_takes_a_rate_expression(tmp_path):
    # the same text as the flag, quoted as a TOML string — no second notation for a rate
    (tmp_path / "p.toml").write_text(
        'birth = "1.0 * OnTime({0: 1.0, 3: 0.3})"\ndeath = 0.3\ntotal-time = 5\n')
    out = tmp_path / "o"
    rc = main(["species", str(out), "--params", str(tmp_path / "p.toml"), "--seed", "2", "--flat"])
    assert rc == 0
    assert "birth\t1.0 * OnTime({0: 1, 3: 0.3})" in (out / "species.log").read_text()


def test_params_file_rate_expression_matches_the_flag(tmp_path):
    (tmp_path / "p.toml").write_text('birth = "1.0 * OnTotalDiversity(cap=20)"\n')
    viafile, viaflag = tmp_path / "f", tmp_path / "c"
    main(["species", str(viafile), "--params", str(tmp_path / "p.toml"), "--total-time", "5", "--seed", "3", "--flat"])
    main(["species", str(viaflag), "--birth", "1.0 * OnTotalDiversity(cap=20)", "--total-time", "5", "--seed", "3", "--flat"])
    assert (viafile / "species_complete.nwk").read_text() == \
        (viaflag / "species_complete.nwk").read_text()


def test_the_rates_help_lists_only_what_the_level_wires(capsys):
    # the help is built from each level's WIRED_MODIFIERS, so it cannot advertise the unwired
    for command, present, absent in [("species", "FromParent", "ByLineage"),
                                     ("sequences", "ByLineage", "FromParent")]:
        with pytest.raises(SystemExit):
            main([command, "--help", "--flat"])
        out = capsys.readouterr().out
        block = out[out.index("RATES"):]
        assert present in block and absent not in block


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
    # and nothing loose at the top
    assert not [p for p in tmp_path.iterdir() if p.is_file()]


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
    assert outs["grouped"] == outs["flat"]      # same files, only the directories differ


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
    a = (tmp_path / "rep1" / "genomes" / "genome_events.tsv").read_text()
    b = (tmp_path / "rep2" / "genomes" / "genome_events.tsv").read_text()
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
    assert (a / "genomes" / "genome_events.tsv").read_text() == \
           (b / "genomes" / "genome_events.tsv").read_text()


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
        "trait_values.tsv", "trait_events.tsv", "trait_tree.nwk"}
    states = {ln.split("\t")[1] for ln in
              (tmp_path / "traits" / "trait_values.tsv").read_text().splitlines()[1:]}
    assert states <= {"small", "large"}


def test_joint_genome_driver_nests_its_gene_trees(tmp_path):
    rc = main(["joint", str(tmp_path),
               "--birth", "1.0 * DrivenBy('genomes:toxin', {'present': 3.0, 'absent': 1.0})",
               "--origination", "0.2", "--loss", "0.1", "--families", "toxin",
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
    assert (tmp_path / "x" / "species" / "species_complete.nwk").read_text() == \
           (tmp_path / "y" / "species" / "species_complete.nwk").read_text()


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
          "--gene-length", "150", "--inversion", "3.0", "--inversion-length", "250",
          "--duplication", "1.0", "--loss", "1.0", "--seed", "3", "--quiet"])
    main(["sequences", run, "--model", "hky85", "--kappa", "3.0", "--substitution", "0.05",
          "--seed", "3", "--quiet", *extra])
    return tmp_path / "run" / "sequences"


def test_sequences_runs_on_a_nucleotide_handoff_and_assembles_the_genomes(tmp_path):
    # the handoff says what it is — blocks.tsv is the nucleotide resolution's and no other's — so
    # nothing has to be repeated from the genomes command
    out = _nucleotide_run(tmp_path)
    names = {p.name for p in out.iterdir()}
    assert {p.name for p in (out / "alignments").iterdir()} and \
           all(n.startswith("block") for n in {p.name for p in (out / "alignments").iterdir()})
    genomes = sorted(n for n in names if n.startswith("genome_n"))
    assert len(genomes) == 9                                  # one FASTA per node of the tree
    for name in genomes:
        text = (out / name).read_text()
        assert text.startswith(">") and set("".join(text.splitlines()[1:])) <= set("ACGT")


def test_sequences_writes_a_genome_for_every_node_and_the_initial_one(tmp_path):
    # all of them by default, named by whose they are: no node is a special case, and the one that
    # belongs to no node is called "initial"
    out = _nucleotide_run(tmp_path)
    tree, _ = read_newick((tmp_path / "run" / "species" / "species_complete.nwk").read_text())
    names = {p.name for p in out.iterdir()}
    assert {f"genome_n{i}.fasta" for i in tree.nodes} <= names
    assert "genome_initial.fasta" in names
    assert not [n for n in names if "ancestral" in n]
    initial = "".join((out / "genome_initial.fasta").read_text().splitlines()[1:])
    assert len(initial) == 2000                               # the genome the run was seeded with


def test_sequences_matches_the_python_api_on_a_nucleotide_run(tmp_path):
    # the CLI is a shell: rebuilt from disk, it must give exactly what the objects give in memory
    from zombi2.genomes import simulate_genomes_nucleotide
    from zombi2.sequences import simulate_sequences
    from zombi2.sequences.substitution_models import hky85

    out = _nucleotide_run(tmp_path)
    sp = simulate_species_tree(birth=1.0, death=0.3, n_extant=5, seed=3)
    g = simulate_genomes_nucleotide(sp, root_length=2000, genes=5, gene_length=150, inversion=3.0,
                                    inversion_length=250, duplication=1.0, loss=1.0, seed=3)
    r = simulate_sequences(g, model=hky85(kappa=3.0), substitution=0.05, seed=3)
    for lineage, chroms in r.genomes.items():
        text = (out / f"genome_{lineage}.fasta").read_text()
        assert "".join(text.splitlines()[1:]) == "".join(chroms.values())


def test_sequences_rejects_length_and_an_intergene_model_where_they_do_not_apply(tmp_path):
    with pytest.raises(SystemExit):
        _nucleotide_run(tmp_path, extra=("--length", "300"))   # the genome sets the lengths
    run = str(tmp_path / "unordered")
    main(["species", run, "--birth", "1.0", "--n-extant", "4", "--seed", "1", "--quiet"])
    main(["genomes", run, "--initial-families", "4", "--duplication", "0.2", "--loss", "0.2",
          "--seed", "1", "--quiet"])
    with pytest.raises(SystemExit):                            # no spacer in a gene-family run
        main(["sequences", run, "--model", "jc69", "--intergene-model", "jc69", "--quiet"])
