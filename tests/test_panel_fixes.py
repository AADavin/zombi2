"""The panel review's findings, each pinned by the test that would have caught it.

A five-reviewer review of 0.27.0 found ten defects that shared one shape: the software fails loudly
and helpfully almost everywhere, which is exactly what makes its few *silent* paths dangerous — users
learn to trust the noise. Each test here is one of those paths, written so that the silence comes
back as a failure rather than as a wrong dataset.
"""

from __future__ import annotations

import subprocess
import sys
import textwrap

import numpy as np
import pytest

from zombi2 import genomes, sequences, species, traits
from zombi2.cli.main import main
from zombi2.params import PerCopy, PerLineage, PerSite
from zombi2.params.evaluate import Modifier
from zombi2.rng import seed_sequence, stream
from zombi2.sequences.substitution_models import hky85, jc69
from zombi2.tree import read_newick

PY = sys.executable

#: one small tree, shared by the engine-coverage cases below
_TREE = species.simulate_species_tree(birth=1.0, death=0.3, n_extant=8, seed=1)


# --- 1. one seed, one stream per level ------------------------------------------------------------

def test_two_levels_at_the_same_seed_draw_different_numbers():
    """The heart of it. Every level used to open ``default_rng(seed)``, so two levels handed the same
    integer replayed the *same* PCG64 stream — and a tree's height and its genome's copy count came
    out correlated at Spearman −0.79 over 6000 seeds, with nothing in the output to say so. SPEC §2
    calls two levels that do not read each other independent; that is a claim about two streams."""
    for a, b in (("species", "genomes"), ("species", "traits"), ("genomes", "sequences"),
                 ("sequences", "traits"), ("traits", "joint")):
        first = stream(a, 42)[0].random(64)
        second = stream(b, 42)[0].random(64)
        assert not np.array_equal(first, second), f"{a} and {b} share a stream at seed 42"


def test_the_levels_are_uncorrelated_at_equal_seeds():
    """The property the streams exist for, measured the way the reviewer measured it: a species tree's
    height and a genome's copy count, both at seed ``s``, over many seeds. Correlated at −0.79 before;
    the bound here is loose enough never to flake and tight enough to catch a shared stream."""
    stem = species.simulate_species_tree(birth=0.0, death=0.0, total_time=2.0, seed=1).complete_tree
    heights, copies = [], []
    for s in range(400):
        heights.append(species.simulate_species_tree(
            birth=1.0, death=0.0, n_extant=8, seed=s).summary()["tree"]["height"])
        copies.append(len(genomes.simulate_genomes_family(
            stem, duplication=0.5, initial_families=1, max_family_size=None, seed=s).genomes[0]))
    r = float(np.corrcoef(np.argsort(np.argsort(heights)), np.argsort(np.argsort(copies)))[0, 1])
    assert abs(r) < 0.2, f"levels at equal seeds are correlated (Spearman {r:+.3f})"


def test_the_species_level_keeps_the_stream_it_always_had():
    """Independence needs the levels' streams to *differ*, not to all be new — so one level can stay
    put, and species is the one worth keeping: every tree ever published under a seed still
    reproduces. Its sequence is the unkeyed root, which is what ``default_rng(seed)`` builds."""
    keyed, _ = seed_sequence("species", 7)
    assert np.array_equal(np.random.default_rng(keyed).random(32),
                          np.random.default_rng(7).random(32))


# --- 2. a seed the caller did not give is still written down --------------------------------------

@pytest.mark.parametrize("run", [
    lambda seed: species.simulate_species_tree(birth=1.0, death=0.3, n_extant=8, seed=seed),
    lambda seed: genomes.simulate_genomes_family(
        species.simulate_species_tree(birth=1.0, death=0.3, n_extant=6, seed=1),
        duplication=0.2, loss=0.2, initial_families=4, seed=seed),
    lambda seed: traits.simulate_continuous(
        species.simulate_species_tree(birth=1.0, death=0.3, n_extant=6, seed=1),
        rate=1.0, seed=seed),
])
def test_an_unseeded_python_run_records_the_seed_it_drew(run):
    """``docs/reproducibility.md`` says a seed you did not give is still written down, and that there
    is no such thing as an unrepeatable ZOMBI2 run. That was true of the CLI and false of the API:
    ``result.seed`` stayed ``None`` and the run was gone — which is precisely what happens while
    exploring interactively, the way anyone finds an interesting realisation."""
    first = run(None)
    assert isinstance(first.seed, int), "an unseeded run must still record the seed it drew"
    assert repr(first.seed) in repr(first)
    again = run(first.seed)
    assert again.seed == first.seed


# --- 3. Newick round-trips exactly, so the two front doors are one run ----------------------------

def test_a_written_tree_reads_back_as_the_same_tree():
    """The CLI hands a tree between levels through this file, so a rounded length is a different
    tree. At the old fixed 12 digits every branch shifted by ~2e−12 on the round trip, which moved
    every downstream waiting time: the same tree and the same seed then gave *different* histories
    through Python and through the CLI."""
    tree = species.simulate_species_tree(birth=1.0, death=0.3, n_extant=25, seed=5).complete_tree
    back, _ = read_newick(tree.to_newick())
    for i, node in tree.nodes.items():
        assert back.nodes[i].end_time == node.end_time, f"node {i} moved on the round trip"
    assert back.to_newick() == tree.to_newick()


def test_python_and_the_cli_are_the_same_run(tmp_path):
    """The whole point of the exact round trip, end to end."""
    run = tmp_path / "cli"
    main(["species", str(run), "--birth", "1.0", "--death", "0.3", "--n-extant", "20",
          "--seed", "5", "--quiet"])
    main(["genomes", str(run), "--duplication", "0.3", "--loss", "0.3", "--initial-families", "10",
          "--seed", "7", "--quiet"])
    sp = species.simulate_species_tree(birth=1.0, death=0.3, n_extant=20, seed=5)
    g = genomes.simulate_genomes_family(sp, duplication=0.3, loss=0.3, initial_families=10, seed=7)
    g.write(tmp_path / "py")
    assert (tmp_path / "py" / "genome_events.tsv").read_text(encoding="utf-8") == \
           (run / "genomes" / "genome_events.tsv").read_text(encoding="utf-8")


# --- 4. an interrupted species run is not a species run -------------------------------------------

def test_a_species_run_missing_its_fates_is_refused(tmp_path):
    """What SIGINT during the write phase leaves, and equally a partial rsync or a tidied directory.
    It used to be consumed with exit 0: under ``--sampling`` every unsampled tip read back as
    sampled, so the same seed and rates gave 30 extant genomes where the complete run gave 14."""
    run = tmp_path / "part"
    main(["species", str(run), "--birth", "1", "--death", "0.3", "--n-extant", "30",
          "--sampling", "0.5", "--seed", "8", "--quiet"])
    for path in run.rglob("*"):
        if path.is_file() and path.name != "species_complete.nwk":
            path.unlink()
    assert main(["genomes", str(run), "--duplication", "0.2", "--loss", "0.2",
                 "--seed", "1", "--quiet"]) != 0


def test_an_external_tree_file_still_goes_straight_in(tmp_path):
    """The other half: a published phylogeny has no fates file and never will, so the guard must not
    catch it. It is a *file*, which is how the two are told apart."""
    ext = tmp_path / "ext.nwk"
    ext.write_text("((A:1,B:1):1,C:2);\n", encoding="utf-8")
    main(["genomes", str(tmp_path / "out"), "--from", str(ext), "--duplication", "0.2",
          "--loss", "0.2", "--origination", "0.5", "--seed", "1", "--quiet"])
    assert (tmp_path / "out" / "genomes" / "genome_events.tsv").exists()


# --- 5. the family-size cap says so, on stderr ----------------------------------------------------

def test_a_cap_that_bound_is_reported_on_stderr(tmp_path, capsys):
    """It was in ``run.zombi2`` and in Chapter 4's prose, and nowhere a job log would keep it — while
    the far less damaging all-genomes-empty case got a loud stderr warning. At ``--duplication 0.8``
    the default cap has been measured pulling the realised rate to 0.32."""
    run = tmp_path / "r"
    main(["species", str(run), "--birth", "1", "--death", "0.3", "--n-extant", "20",
          "--seed", "4242", "--quiet"])
    main(["genomes", str(run), "--duplication", "0.8", "--loss", "0.1", "--initial-families", "30",
          "--seed", "3", "--quiet"])
    err = capsys.readouterr().err
    assert "max-family-size" in err and "below the ones you declared" in err


def test_a_cap_that_did_not_bind_says_nothing(tmp_path, capsys):
    """A healthy run prints nothing, so an empty stderr still means "nothing to report"."""
    run = tmp_path / "r"
    main(["species", str(run), "--birth", "1", "--death", "0.3", "--n-extant", "10",
          "--seed", "1", "--quiet"])
    main(["genomes", str(run), "--duplication", "0.05", "--loss", "0.05",
          "--initial-families", "5", "--seed", "1", "--quiet"])
    assert "max-family-size" not in capsys.readouterr().err


# --- 6. a mistyped output token is not a silent no-op ---------------------------------------------

def test_a_bogus_write_token_raises_at_every_level(tmp_path):
    """Species, sequences and traits always raised; the two genome results wrote nothing and exited
    clean — silent data loss you discover three pipeline steps later, when the next tool has no
    input."""
    sp = species.simulate_species_tree(birth=1.0, death=0.3, n_extant=6, seed=1)
    fam = genomes.simulate_genomes_family(sp, duplication=0.2, loss=0.2, initial_families=4, seed=1)
    ordered = genomes.simulate_genomes_ordered(sp, duplication=0.2, loss=0.2, initial_families=4,
                                               chromosomes=1, seed=1)
    for result in (sp, fam, ordered,
                   sequences.simulate_sequences(fam, model=jc69(), length=20, seed=1),
                   traits.simulate_continuous(sp, rate=1.0, seed=1)):
        with pytest.raises(ValueError, match="unknown write outputs"):
            result.write(tmp_path / "nope", outputs=("gene_treez",))


# --- 7. duplicate tip labels ----------------------------------------------------------------------

def test_a_duplicate_tip_label_is_refused():
    """The simulation itself would be fine, which is what makes this worth refusing rather than
    warning about: nothing downstream would ever look wrong. ``names.tsv`` — the documented join back
    to the caller's taxa — would map two node ids to one name, and every merge by taxon name would
    silently duplicate or drop rows."""
    with pytest.raises(ValueError, match="duplicate tip label"):
        read_newick("((A:1,B:1):1,A:2);")
    tree, names = read_newick("((A:1,B:1):1,C:2);")
    assert sorted(names.values()) == ["A", "B", "C"]


def test_a_zombi_tree_is_unaffected():
    """A complete tree's labels are ids, so it cannot reach the new check — and must not."""
    tree = species.simulate_species_tree(birth=1.0, death=0.3, n_extant=8, seed=1).complete_tree
    assert read_newick(tree.to_newick())[0].to_newick() == tree.to_newick()


# --- 8. an unguarded parallel script says what is wrong -------------------------------------------

def _unguarded_parallel_script(tmp_path, start_method: str) -> subprocess.CompletedProcess:
    """A `.py` whose simulate call sits at the top level, run under a chosen start method."""
    script = tmp_path / f"unguarded_{start_method}.py"
    script.write_text(textwrap.dedent(f"""
        import multiprocessing
        multiprocessing.set_start_method({start_method!r}, force=True)
        from zombi2 import species, genomes
        sp = species.simulate_species_tree(birth=1.0, death=0.3, n_extant=12, seed=2)
        genomes.simulate_genomes_family(sp, duplication=0.3, loss=0.3, origination=0.5,
                                        initial_families=20, seed=4, parallel=2)
        print("finished")
    """), encoding="utf-8")
    return subprocess.run([PY, str(script)], capture_output=True, text=True, timeout=300)


@pytest.mark.skipif("spawn" not in __import__("multiprocessing").get_all_start_methods(),
                    reason="no spawn start method on this platform")
def test_an_unguarded_parallel_script_names_the_main_guard(tmp_path):
    """Under ``spawn`` a worker re-imports the caller's script and runs its top level again. Left
    alone that is worse than a crash: every child silently repeated the whole simulation and the run
    "succeeded", having done N× the work and printed N copies of its output. The parent used to see
    only ``BrokenProcessPool``, which names nothing."""
    done = _unguarded_parallel_script(tmp_path, "spawn")
    assert done.returncode != 0
    assert "__main__" in done.stderr
    assert done.stdout.count("finished") == 0, "no worker may have re-run the whole program"


@pytest.mark.skipif("fork" not in __import__("multiprocessing").get_all_start_methods(),
                    reason="no fork start method on this platform")
def test_the_guard_is_not_demanded_where_fork_makes_it_pointless(tmp_path):
    """The other half, and the reason the check is on the start method rather than on the script:
    ``fork`` (the Linux default) inherits the parent instead of re-importing it, so there is nothing
    to guard against and an unguarded script is perfectly correct. Refusing it there would be us
    inventing a requirement the platform does not have."""
    done = _unguarded_parallel_script(tmp_path, "fork")
    assert done.returncode == 0, done.stderr
    assert done.stdout.count("finished") == 1


# --- 9. a written run reads back ------------------------------------------------------------------

def test_a_written_run_reads_back_and_feeds_the_next_level(tmp_path):
    """``zombi2 sequences --from DIR`` always reopened a genomes run; from Python the same handoff
    was a dead end. It mattered most for ``stream_to=``, whose whole point is that the run does not
    fit in memory — and whose handle could then not be passed on."""
    sp = species.simulate_species_tree(birth=1.0, death=0.3, n_extant=10, seed=2)
    grown = genomes.simulate_genomes_family(sp, duplication=0.2, loss=0.2, initial_families=5, seed=1)
    grown.write(tmp_path / "run")

    back = genomes.read_run(tmp_path / "run")
    assert len(back.gene_trees) == len(grown.gene_trees)
    assert len(back.events) == len(grown.events)
    assert back.profiles.shape == grown.profiles.shape          # not a KeyError
    assert sequences.simulate_sequences(tmp_path / "run", model=hky85(), length=20, seed=1).seed == 1


def test_a_streamed_run_feeds_the_next_level(tmp_path):
    streamed = genomes.simulate_genomes_family(
        species.simulate_species_tree(birth=1.0, death=0.3, n_extant=10, seed=2),
        duplication=0.2, loss=0.2, initial_families=5, seed=1, stream_to=str(tmp_path / "s"))
    assert (tmp_path / "s" / "species_complete.nwk").exists(), "a streamed run must be self-contained"
    from_handle = sequences.simulate_sequences(streamed, model=hky85(), length=20, seed=1)
    from_path = sequences.simulate_sequences(str(tmp_path / "s"), model=hky85(), length=20, seed=1)
    assert sorted(from_handle.alignments) == sorted(from_path.alignments)


def test_a_run_read_back_without_its_gene_content_says_so(tmp_path):
    """``genomes.tsv`` is an optional output. Without it there are no profiles to derive — which used
    to surface as a bare ``KeyError`` from an accessor."""
    sp = species.simulate_species_tree(birth=1.0, death=0.3, n_extant=8, seed=2)
    g = genomes.simulate_genomes_family(sp, duplication=0.2, loss=0.2, initial_families=4, seed=1)
    g.write(tmp_path / "thin", outputs=("events", "species_tree"))
    back = genomes.read_run(tmp_path / "thin")
    assert back.gene_trees                                   # the genealogy is all there
    assert "gene content not loaded" in repr(back)           # and the repr does not claim 0 nodes
    with pytest.raises(ValueError, match="no per-node gene content"):
        back.profiles


# --- 10. a modifier of your own ------------------------------------------------------------------

class OnLogTime(Modifier):
    """A third-party modifier: it declares the engine it is wired for, and takes ``**_`` because the
    engine supplies whatever context it happens to have."""

    implemented_for = ("species",)

    def factor(self, *, time: float = 0.0, **_) -> float:
        return 1.0 / (1.0 + time)


def test_a_declared_third_party_modifier_runs():
    """`Modifier` is public with a clean ``factor(**context)`` contract, and a subclass composed into
    a `Rate` correctly — and was then refused by every level, with no registry and no entry point, so
    extending the grammar meant forking the package."""
    result = species.simulate_species_tree(birth=PerLineage(2.0)._and(OnLogTime()),
                                           n_extant=20, seed=1)
    assert len(result.complete_tree.extant_leaves()) == 20


def test_the_gate_still_holds_for_everything_undeclared():
    """The gate itself is right — a modifier the engine never reads returns its default 1.0 and gives
    a run that is quietly not the model asked for (SPEC §5). Only an explicit declaration opens it,
    and only for the engine named."""
    class Undeclared(Modifier):
        def factor(self, **_) -> float:
            return 2.0

    with pytest.raises(ValueError, match="does not support"):
        species.simulate_species_tree(birth=PerLineage(1.0)._and(Undeclared()),
                                      n_extant=8, seed=1)

    tree = species.simulate_species_tree(birth=1.0, death=0.3, n_extant=6, seed=1)
    with pytest.raises(ValueError, match="does not support"):   # wired for species, not for traits
        traits.simulate_continuous(tree, rate=PerLineage(1.0)._and(OnLogTime()), seed=1)


@pytest.mark.parametrize("engine, run", [
    ("species", lambda m: species.simulate_species_tree(
        birth=PerLineage(1.0)._and(m), death=0.2, n_extant=8, seed=1)),
    ("genomes.family", lambda m: genomes.simulate_genomes_family(
        _TREE, duplication=PerCopy(0.2)._and(m), loss=0.2, initial_families=4, seed=1)),
    ("genomes.ordered", lambda m: genomes.simulate_genomes_ordered(
        _TREE, duplication=PerCopy(0.2)._and(m), loss=0.2, initial_families=4, chromosomes=1,
        seed=1)),
    ("genomes.nucleotide", lambda m: genomes.simulate_genomes_nucleotide(
        _TREE, loss=PerLineage(0.5)._and(m), root_length=400, seed=1)),
    ("traits.continuous", lambda m: traits.simulate_continuous(
        _TREE, rate=PerLineage(1.0)._and(m), seed=1)),
    ("traits.discrete", lambda m: traits.simulate_discrete(
        _TREE, states=["a", "b"], switch=PerLineage(0.5)._and(m), seed=1)),
])
def test_every_advertised_engine_actually_calls_a_third_party_modifier(engine, run):
    """The gate opening is only half of it: the engine has to *call* the thing it let through.

    `traits.discrete` did not. It chose between a constant generator and a rebuilt-per-stretch one by
    asking "are there drivers?", so a modifier that was not a `Driven` fell between the two and
    crashed the driver resolver looking for a `.key`. This asserts the whole advertised list, because
    the failure mode is per engine and invisible from the outside."""
    calls: list[tuple[str, ...]] = []

    class Spy(Modifier):
        implemented_for = (engine,)

        def factor(self, **context) -> float:
            calls.append(tuple(sorted(context)))
            return 1.0

    run(Spy())
    assert calls, f"the {engine} engine accepted the modifier and never called factor()"


def test_the_sequence_level_refuses_rather_than_ignoring():
    """The one engine that does not offer the hatch, and must say so. It reads its modifiers itself —
    the clock is drawn per lineage before any site evolves — so one it did not ship could be accepted
    and then never called, which is the silence the whole mechanism exists to prevent."""
    class Spy(Modifier):
        implemented_for = ("sequences",)

        def factor(self, **_) -> float:
            return 1.0

    g = genomes.simulate_genomes_family(_TREE, duplication=0.2, loss=0.2, initial_families=4, seed=1)
    with pytest.raises(ValueError, match="does not read"):
        sequences.simulate_sequences(g, model=jc69(), length=20,
                                     substitution=PerSite(1.0)._and(Spy()), seed=1)


def test_the_built_in_modifiers_are_unaffected():
    tree = species.simulate_species_tree(birth=PerLineage(1.0).changing_at({0: 1.0, 3: 0.3}), death=0.2,
                                         n_extant=15, seed=1)
    assert len(tree.complete_tree.extant_leaves()) == 15
