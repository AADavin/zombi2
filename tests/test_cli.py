"""Tests for the ``zombi2`` command line (species + genomes) and the ``read_newick`` reader it
loads trees with. The CLI is a thin shell over the level ``simulate_*`` functions, so these tests
exercise argument wiring, output files, ``--params``, and the clean error paths — not the engines
(which their own suites cover)."""

import pytest

from zombi2.cli.main import main
from zombi2.genomes import simulate_genomes_unordered
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
    # every non-root branch length (duration) survives to Newick's 6 significant figures
    for i, n in t.nodes.items():
        if i == t.root:
            continue
        assert back.nodes[i].end_time - back.nodes[i].birth_time == pytest.approx(
            n.end_time - n.birth_time, rel=1e-5, abs=1e-9)
    # the crown-rooted convention drops the root's own branch length (SPEC §8): it reads back as a
    # zero-duration crown, so the tree starts at the first split
    assert back.nodes[back.root].end_time - back.nodes[back.root].birth_time == 0.0


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
    rc = main(["species", "--birth", "1", "--death", "0.3", "--n-extant", "20",
               "--seed", "1", "-o", str(tmp_path)])
    assert rc == 0
    written = {p.name for p in tmp_path.iterdir()}
    assert {"species_complete.nwk", "species_extant.nwk", "species_events.tsv",
            "species.log"} <= written
    assert "extant" in capsys.readouterr().out
    # the reproducibility log records the resolved parameters
    log = (tmp_path / "species.log").read_text()
    assert "birth\t1.0" in log and "n_extant\t20" in log


def test_species_write_is_selective(tmp_path):
    main(["species", "--birth", "1", "--total-time", "3", "--seed", "1",
          "--write", "complete", "-o", str(tmp_path)])
    nwk = {p.name for p in tmp_path.iterdir() if p.suffix == ".nwk"}
    assert nwk == {"species_complete.nwk"}                       # no extant/events file


def test_species_is_deterministic_given_the_seed(tmp_path):
    a, b = tmp_path / "a", tmp_path / "b"
    for out in (a, b):
        main(["species", "--birth", "1", "--death", "0.3", "--n-extant", "30",
              "--seed", "13", "-o", str(out)])
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
    rc = main(["species", "--birth", "1", "--n-extant", "10", "--mass-extinction", "3", "0.5",
               "-o", str(tmp_path)])
    assert rc == 1
    assert "zombi2: error:" in capsys.readouterr().err


# ── zombi2 genomes ──────────────────────────────────────────────────────────────────

@pytest.fixture
def tree_file(tmp_path):
    """A species tree written to disk, for the genomes command to read."""
    main(["species", "--birth", "1", "--death", "0.3", "--n-extant", "25", "--seed", "1",
          "-o", str(tmp_path)])
    return tmp_path / "species_complete.nwk"


def test_genomes_unordered_writes_events_and_profiles(tmp_path, tree_file):
    out = tmp_path / "g"
    rc = main(["genomes", "-t", str(tree_file), "--duplication", "0.2", "--transfer", "0.1",
               "--loss", "0.25", "--origination", "0.5", "--seed", "42", "-o", str(out)])
    assert rc == 0
    # genome_species_tree.nwk is the always-written handoff tree for `zombi2 sequences --genomes`
    assert {p.name for p in out.iterdir()} == {"genome_events.tsv", "profiles.tsv",
                                               "genome_species_tree.nwk", "genomes.log"}


def test_genomes_ordered_writes_structured_outputs(tmp_path, tree_file):
    out = tmp_path / "g"
    rc = main(["genomes", "-t", str(tree_file), "--resolution", "ordered", "--duplication", "0.2",
               "--loss", "0.2", "--origination", "0.5", "--inversion", "0.3", "--chromosomes", "3",
               "--seed", "42", "-o", str(out), "--write", "gene_order", "rearrangements"])
    assert rc == 0
    assert {p.name for p in out.iterdir()} == {"gene_order.tsv", "rearrangements.tsv",
                                               "genome_species_tree.nwk", "genomes.log"}


def test_genomes_rejects_ordered_only_flag_under_unordered(tmp_path, tree_file):
    with pytest.raises(SystemExit) as e:
        main(["genomes", "-t", str(tree_file), "--inversion", "0.3", "-o", str(tmp_path / "g")])
    assert e.value.code == 2


def test_genomes_rejects_write_output_foreign_to_resolution(tmp_path, tree_file):
    with pytest.raises(SystemExit) as e:                         # gene_order is ordered-only
        main(["genomes", "-t", str(tree_file), "--write", "gene_order", "-o", str(tmp_path / "g")])
    assert e.value.code == 2


def test_genomes_missing_tree_is_reported_cleanly(tmp_path, capsys):
    rc = main(["genomes", "-t", str(tmp_path / "nope.nwk"), "--duplication", "0.1",
               "-o", str(tmp_path / "g")])
    assert rc == 1
    assert "tree file not found" in capsys.readouterr().err


def test_genomes_on_ultrametric_external_tree_writes_a_name_map(tmp_path):
    (tmp_path / "ext.nwk").write_text("((human:1,chimp:1):1,(mouse:0.8,rat:0.8):1.2);\n")
    out = tmp_path / "g"
    rc = main(["genomes", "-t", str(tmp_path / "ext.nwk"), "--duplication", "0.3",
               "--origination", "1.0", "--seed", "1", "-o", str(out)])
    assert rc == 0
    # all four tips are observed, so the profile matrix has four columns
    assert len((out / "profiles.tsv").read_text().splitlines()[0].split("\t")) == 1 + 4
    # names.tsv maps ZOMBI's n<id> back to the user's labels
    mapped = dict(row.split("\t") for row in (out / "names.tsv").read_text().splitlines()[1:])
    assert sorted(mapped.values()) == ["chimp", "human", "mouse", "rat"]


def test_genomes_on_nonultrametric_tree_needs_tip_fates(tmp_path, capsys):
    (tmp_path / "ext.nwk").write_text("((a:1,b:1):1,c:1.5);\n")     # c ends early
    rc = main(["genomes", "-t", str(tmp_path / "ext.nwk"), "--duplication", "0.3",
               "--seed", "1", "-o", str(tmp_path / "g")])
    assert rc == 1
    assert "not ultrametric" in capsys.readouterr().err


def test_genomes_nonultrametric_tree_runs_with_tip_fates_file(tmp_path):
    (tmp_path / "ext.nwk").write_text("((a:1,b:1):1,c:1.5);\n")
    (tmp_path / "fates.tsv").write_text("a\textant\nb\textant\nc\textinct\n")
    out = tmp_path / "g"
    rc = main(["genomes", "-t", str(tmp_path / "ext.nwk"), "--tip-fates", str(tmp_path / "fates.tsv"),
               "--duplication", "0.3", "--origination", "1.0", "--seed", "1", "-o", str(out)])
    assert rc == 0
    # c is extinct, so only a and b are observed → two profile columns
    assert len((out / "profiles.tsv").read_text().splitlines()[0].split("\t")) == 1 + 2


# ── zombi2 sequences ────────────────────────────────────────────────────────────────

@pytest.fixture
def genomes_dir(tmp_path):
    """A completed species→genomes run on disk, for the sequences command to replay."""
    main(["species", "--birth", "1", "--death", "0.3", "--n-extant", "25", "--seed", "1",
          "-o", str(tmp_path)])
    gdir = tmp_path / "g"
    main(["genomes", "-t", str(tmp_path / "species_complete.nwk"), "--duplication", "0.2",
          "--transfer", "0.1", "--loss", "0.25", "--origination", "0.6", "--seed", "42",
          "-o", str(gdir)])
    return gdir


def test_sequences_writes_alignments_and_phylograms_by_default(tmp_path, genomes_dir):
    out = tmp_path / "s"
    rc = main(["sequences", "--genomes", str(genomes_dir), "--model", "hky85", "--kappa", "2",
               "--length", "300", "--seed", "1", "-o", str(out)])
    assert rc == 0
    names = {p.name for p in out.iterdir()}
    assert "sequences.log" in names
    # the default write set is alignments + phylograms; ancestral / species_phylogram are opt-in
    assert any(n.startswith("sequences_alignment_fam") for n in names)
    assert any(n.startswith("sequences_phylogram_fam") for n in names)
    assert not any(n.startswith("sequences_ancestral_fam") for n in names)
    assert not any(n.startswith("sequences_species_phylogram") for n in names)


def test_sequences_write_selects_ancestral_and_species_phylogram(tmp_path, genomes_dir):
    out = tmp_path / "s"
    rc = main(["sequences", "--genomes", str(genomes_dir), "--model", "jc69", "--length", "200",
               "--seed", "1", "-o", str(out), "--write", "ancestral", "species_phylogram"])
    assert rc == 0
    names = {p.name for p in out.iterdir()}
    # the species phylogram is produced only because the CLI hands the engine the species tree
    assert "sequences_species_phylogram_complete.nwk" in names
    assert any(n.startswith("sequences_ancestral_fam") for n in names)
    assert not any(n.startswith("sequences_alignment_fam") for n in names)   # not requested


def test_sequences_relaxed_clock_runs_and_is_logged(tmp_path, genomes_dir):
    # the relaxed clock is not its own flag: it is a ByLineage modifier on the substitution rate
    out = tmp_path / "s"
    rc = main(["sequences", "--genomes", str(genomes_dir), "--model", "gtr",
               "--frequencies", "0.3", "0.2", "0.2", "0.3",
               "--substitution", "1.0 * ByLineage(spread=0.4, dist='gamma')",
               "--seed", "1", "-o", str(out)])
    assert rc == 0
    log = (out / "sequences.log").read_text()
    assert "gamma lineage clock, spread 0.4" in log
    # the rate is logged in its written form, so the log line pastes back into --substitution
    assert "substitution\t1.0 * ByLineage(spread=0.4, dist='gamma')" in log


def test_sequences_rejects_a_model_foreign_parameter(tmp_path, genomes_dir):
    with pytest.raises(SystemExit) as e:                         # --kappa is meaningless for jc69
        main(["sequences", "--genomes", str(genomes_dir), "--model", "jc69", "--kappa", "2",
              "-o", str(tmp_path / "s")])
    assert e.value.code == 2


def test_sequences_is_deterministic_given_the_seed(tmp_path, genomes_dir):
    a, b = tmp_path / "a", tmp_path / "b"
    for out in (a, b):
        main(["sequences", "--genomes", str(genomes_dir), "--model", "hky85", "--length", "250",
              "--seed", "7", "-o", str(out)])
    assert _dir_seq_text(a) == _dir_seq_text(b)


def _dir_seq_text(d):
    return {p.name: p.read_text() for p in d.iterdir() if p.suffix in (".fasta", ".nwk")}


def test_sequences_missing_genomes_dir_is_reported_cleanly(tmp_path, capsys):
    rc = main(["sequences", "--genomes", str(tmp_path / "nope"), "--model", "jc69",
               "-o", str(tmp_path / "s")])
    assert rc == 1
    assert "genome_species_tree.nwk not found" in capsys.readouterr().err


def test_sequences_needs_the_genome_event_log(tmp_path, capsys):
    # a genomes run written with --write profiles has no event log to replay
    main(["species", "--birth", "1", "--death", "0.3", "--n-extant", "15", "--seed", "1",
          "-o", str(tmp_path)])
    gdir = tmp_path / "g"
    main(["genomes", "-t", str(tmp_path / "species_complete.nwk"), "--duplication", "0.2",
          "--seed", "1", "-o", str(gdir), "--write", "profiles"])
    rc = main(["sequences", "--genomes", str(gdir), "--model", "jc69", "-o", str(tmp_path / "s")])
    assert rc == 1
    assert "genome_events.tsv not found" in capsys.readouterr().err


# ── traits ──────────────────────────────────────────────────────────────────────────

def test_traits_continuous_writes_values_and_tree(tmp_path, tree_file):
    out = tmp_path / "t"
    rc = main(["traits", "-t", str(tree_file), "--rate", "1.0", "--seed", "1", "-o", str(out)])
    assert rc == 0
    assert {p.name for p in out.iterdir()} == {"trait_values.tsv", "trait_tree.nwk", "traits.log"}
    header, first = (out / "trait_values.tsv").read_text().splitlines()[:2]
    assert header == "node\ttrait"
    assert first.split("\t")[0].startswith("n")          # n<id>, matching the Newick
    float(first.split("\t")[1])                          # a continuous trait is a number


def test_traits_ou_and_threshold_run(tmp_path, tree_file):
    # OU: the same diffusion pulled to an optimum
    assert main(["traits", "-t", str(tree_file), "--rate", "1.0", "--reverts-to", "2",
                 "--pull", "0.5", "--seed", "1", "-o", str(tmp_path / "ou")]) == 0
    # the threshold model: a discrete state read off a continuous liability
    out = tmp_path / "th"
    assert main(["traits", "-t", str(tree_file), "--kind", "discrete", "--states", "absent,present",
                 "--liability", "1.0", "--threshold", "0.0", "--seed", "1", "-o", str(out)]) == 0
    states = {ln.split("\t")[1] for ln in (out / "trait_values.tsv").read_text().splitlines()[1:]}
    assert states <= {"absent", "present"}


def test_traits_discrete_writes_the_event_log_and_driver(tmp_path, tree_file):
    out = tmp_path / "t"
    rc = main(["traits", "-t", str(tree_file), "--kind", "discrete",
               "--states", "marine,terrestrial", "--switch", "0.3", "--start", "marine",
               "--seed", "1", "--write", "values", "changes", "tree", "driver", "-o", str(out)])
    assert rc == 0
    assert {p.name for p in out.iterdir()} == {"trait_values.tsv", "trait_changes.tsv",
                                               "trait_tree.nwk", "trait_driver.tsv", "traits.log"}
    assert (out / "trait_changes.tsv").read_text().splitlines()[0] == "time\tkind\tlineage\tfrom\tto"


def test_traits_at_speciation_logs_on_speciation_changes(tmp_path, tree_file):
    out = tmp_path / "t"
    main(["traits", "-t", str(tree_file), "--rate", "1.0", "--at-speciation", "0.5", "--seed", "1",
          "--write", "changes", "-o", str(out)])
    rows = (out / "trait_changes.tsv").read_text().splitlines()[1:]
    assert rows and all(r.split("\t")[1] == "on_speciation" for r in rows)   # a diffusion has no
    #                                                    along-branch events, only the split jumps


@pytest.mark.parametrize("argv, msg", [
    (["--kind", "discrete", "--states", "a,b", "--switch", "0.1", "--rate", "2"], "--kind continuous"),
    (["--switch", "0.1"], "--kind discrete"),                       # discrete knob, continuous run
    (["--kind", "discrete", "--switch", "0.1"], "--states"),        # discrete without a state space
    (["--kind", "discrete", "--states", "a,b"], "--switch"),        # discrete without a model
    (["--write", "driver"], "not available for --kind continuous"),  # driver is discrete-only
    (["--start", "marine"], "must be a number"),                    # a label on a continuous trait
])
def test_traits_argument_errors_exit_2(argv, msg, tmp_path, tree_file, capsys):
    with pytest.raises(SystemExit) as e:
        main(["traits", "-t", str(tree_file), *argv, "-o", str(tmp_path / "t")])
    assert e.value.code == 2
    assert msg in capsys.readouterr().err


def test_traits_is_deterministic_given_the_seed(tmp_path, tree_file):
    written = []
    for name in ("a", "b"):
        out = tmp_path / name
        main(["traits", "-t", str(tree_file), "--kind", "discrete", "--states", "a,b",
              "--switch", "0.4", "--seed", "99", "-o", str(out)])
        written.append((out / "trait_values.tsv").read_text())
    assert written[0] == written[1]


def test_traits_missing_tree_is_reported_cleanly(tmp_path, capsys):
    rc = main(["traits", "-t", str(tmp_path / "nope.nwk"), "--rate", "1", "-o", str(tmp_path / "t")])
    assert rc == 1                                        # a clean one-line error, not a traceback
    assert "tree file not found" in capsys.readouterr().err


def test_traits_on_external_tree_writes_a_name_map(tmp_path):
    (tmp_path / "ext.nwk").write_text("((human:1,chimp:1):1,(mouse:0.8,rat:0.8):1.2);\n")
    out = tmp_path / "t"
    rc = main(["traits", "-t", str(tmp_path / "ext.nwk"), "--rate", "1.0", "--seed", "1",
               "-o", str(out)])
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
                                     'write = ["values", "changes"]\nseed = 7\n')
    argv = ["traits", "--params", str(tmp_path / "p.toml"), "-t", str(tree_file)]
    out = tmp_path / "a"
    assert main([*argv, "-o", str(out)]) == 0
    assert {p.name for p in out.iterdir()} == {"trait_values.tsv", "trait_changes.tsv", "traits.log"}
    states = {ln.split("\t")[1] for ln in (out / "trait_values.tsv").read_text().splitlines()[1:]}
    assert states <= {"marine", "terrestrial"}          # the file's states reached the engine

    # a flag given on the command line still wins over the file
    other = tmp_path / "b"
    assert main([*argv, "--seed", "8", "-o", str(other)]) == 0
    assert (other / "trait_values.tsv").read_text() != (out / "trait_values.tsv").read_text()


# ── --params ────────────────────────────────────────────────────────────────────────

def test_params_file_supplies_defaults_and_cli_overrides(tmp_path):
    (tmp_path / "p.toml").write_text("birth = 2.0\ndeath = 0.3\nn-extant = 12\n")
    out = tmp_path / "o"
    # birth comes from the file; the command line still overrides it
    main(["species", "--params", str(tmp_path / "p.toml"), "--birth", "1.0", "--seed", "1",
          "-o", str(out)])
    log = (out / "species.log").read_text()
    assert "birth\t1.0" in log and "n_extant\t12" in log


def test_params_file_scopes_by_command_table(tmp_path):
    (tmp_path / "pipeline.toml").write_text(
        "[species]\nbirth = 1.0\nn-extant = 15\n\n[genomes]\nduplication = 0.2\nwrite = "
        '["profiles"]\n')
    sp, gn = tmp_path / "sp", tmp_path / "gn"
    main(["species", "--params", str(tmp_path / "pipeline.toml"), "--seed", "1", "-o", str(sp)])
    main(["genomes", "--params", str(tmp_path / "pipeline.toml"), "-t",
          str(sp / "species_complete.nwk"), "--seed", "1", "-o", str(gn)])
    assert {p.name for p in gn.iterdir()} == {"profiles.tsv", "genome_species_tree.nwk",
                                              "genomes.log"}   # write=["profiles"] (+ handoff tree)


def test_params_unknown_key_errors(tmp_path):
    (tmp_path / "bad.toml").write_text("birth = 1.0\nbogus = 3\n")
    with pytest.raises(SystemExit) as e:
        main(["species", "--params", str(tmp_path / "bad.toml"), "--n-extant", "5",
              "-o", str(tmp_path / "o")])
    assert e.value.code == 2


# ── rates in their written form (SPEC §5) ───────────────────────────────────────────
#
# Every rate flag takes the same expression the Python API takes, so there is one notation for a
# rate across Python, the command line and a --params file. These tests hold that line: the
# expression reaches the engine, it changes the run, an unwired modifier is refused, and the
# parameters log records something you can paste back.

def test_species_takes_a_rate_expression_and_it_bends_the_tree(tmp_path):
    # a skyline that collapses speciation at t=2 must give a smaller tree than the flat rate,
    # i.e. the modifier reached the engine rather than being parsed and dropped
    flat, skyline = tmp_path / "flat", tmp_path / "sky"
    main(["species", "--birth", "1.0", "--death", "0.2", "--total-time", "6",
          "--seed", "4", "-o", str(flat)])
    main(["species", "--birth", "1.0 * OnTime({0: 1.0, 2: 0.05})", "--death", "0.2",
          "--total-time", "6", "--seed", "4", "-o", str(skyline)])
    n = {d: len(read_newick((d / "species_complete.nwk").read_text())[0].nodes)
         for d in (flat, skyline)}
    assert n[skyline] < n[flat]


def test_species_takes_a_scope_wrapper(tmp_path):
    # Global(base) is one budget for the whole tree: linear growth, so far fewer lineages
    out = tmp_path / "g"
    rc = main(["species", "--birth", "Global(2.0)", "--total-time", "5", "--seed", "1",
               "-o", str(out)])
    assert rc == 0
    assert "birth\tGlobal(2.0)" in (out / "species.log").read_text()


def test_species_records_the_rate_in_its_written_form(tmp_path):
    # the log line is the flag value again — a reproducibility record you can paste back
    out = tmp_path / "o"
    main(["species", "--birth", "1.0 * OnTime({0: 1.0, 3: 0.3})", "--total-time", "4",
          "--seed", "1", "-o", str(out)])
    assert "birth\t1.0 * OnTime({0: 1, 3: 0.3})" in (out / "species.log").read_text()


def test_species_refuses_a_modifier_it_does_not_wire(tmp_path, capsys):
    # ByLineage would return a factor of 1.0 at this level — a run quietly not the model asked for
    rc = main(["species", "--birth", "1.0 * ByLineage(spread=0.3)", "--total-time", "3",
               "--seed", "1", "-o", str(tmp_path / "o")])
    assert rc == 1
    assert "does not support" in capsys.readouterr().err


def test_a_typo_in_a_modifier_is_caught_at_the_flag(tmp_path, capsys):
    with pytest.raises(SystemExit) as e:
        main(["species", "--birth", "1.0 * OnDiversity(cap=10)", "--total-time", "3",
              "-o", str(tmp_path / "o")])
    assert e.value.code == 2
    assert "did you mean 'OnTotalDiversity'" in capsys.readouterr().err


def test_a_rate_expression_is_never_executed(tmp_path, capsys):
    with pytest.raises(SystemExit) as e:
        main(["species", "--birth", "__import__('os').system('true')", "--total-time", "3",
              "-o", str(tmp_path / "o")])
    assert e.value.code == 2
    assert "only call a scope or a modifier" in capsys.readouterr().err


def test_genomes_takes_a_rate_expression(tmp_path, tree_file):
    out = tmp_path / "g"
    rc = main(["genomes", "-t", str(tree_file), "--duplication", "0.2",
               "--loss", "0.25 * OnTime({0: 1.0, 2: 3.0})", "--origination", "0.5",
               "--seed", "42", "-o", str(out)])
    assert rc == 0
    assert "loss\t0.25 * OnTime({0: 1, 2: 3})" in (out / "genomes.log").read_text()


def test_traits_takes_a_rate_expression(tmp_path, tree_file):
    out = tmp_path / "t"
    rc = main(["traits", "-t", str(tree_file), "--rate", "1.0 * FromParent(spread=0.2)",
               "--seed", "1", "-o", str(out)])
    assert rc == 0
    assert "rate\t1.0 * FromParent(spread=0.2)" in (out / "traits.log").read_text()


def test_params_file_takes_a_rate_expression(tmp_path):
    # the same text as the flag, quoted as a TOML string — no second notation for a rate
    (tmp_path / "p.toml").write_text(
        'birth = "1.0 * OnTime({0: 1.0, 3: 0.3})"\ndeath = 0.3\ntotal-time = 5\n')
    out = tmp_path / "o"
    rc = main(["species", "--params", str(tmp_path / "p.toml"), "--seed", "2", "-o", str(out)])
    assert rc == 0
    assert "birth\t1.0 * OnTime({0: 1, 3: 0.3})" in (out / "species.log").read_text()


def test_params_file_rate_expression_matches_the_flag(tmp_path):
    (tmp_path / "p.toml").write_text('birth = "1.0 * OnTotalDiversity(cap=20)"\n')
    viafile, viaflag = tmp_path / "f", tmp_path / "c"
    main(["species", "--params", str(tmp_path / "p.toml"), "--total-time", "5", "--seed", "3",
          "-o", str(viafile)])
    main(["species", "--birth", "1.0 * OnTotalDiversity(cap=20)", "--total-time", "5",
          "--seed", "3", "-o", str(viaflag)])
    assert (viafile / "species_complete.nwk").read_text() == \
        (viaflag / "species_complete.nwk").read_text()


def test_the_rates_help_lists_only_what_the_level_wires(capsys):
    # the help is built from each level's WIRED_MODIFIERS, so it cannot advertise the unwired
    for command, present, absent in [("species", "FromParent", "ByLineage"),
                                     ("sequences", "ByLineage", "FromParent")]:
        with pytest.raises(SystemExit):
            main([command, "--help"])
        out = capsys.readouterr().out
        block = out[out.index("RATES"):]
        assert present in block and absent not in block


# ── top-level dispatch ──────────────────────────────────────────────────────────────

def test_version_and_help_do_not_crash(capsys):
    for flag in ("--version", "--help"):
        with pytest.raises(SystemExit) as e:
            main([flag])
        assert e.value.code == 0
    assert "ZOMBI2" in capsys.readouterr().out
