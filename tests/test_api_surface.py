"""Tests for what a Python caller meets first: the type checks, the reprs, the typing marker.

None of this is about the models. It is about the level's *edge* — what happens when a caller hands
it the wrong thing, and what a result says about itself when you type its name in a session.
"""

import pathlib

import pytest

from zombi2 import genomes, sequences, species, traits
from zombi2.tree import as_tree


@pytest.fixture(scope="module")
def run():
    sp = species.simulate_species_tree(birth=1.0, death=0.3, n_extant=12, seed=1)
    g = genomes.simulate_genomes_family(sp, duplication=0.3, transfer=0.1, loss=0.2,
                                        origination=0.5, initial_families=6, seed=1)
    return sp, g


# --- the wrong argument, named where the caller can see it -----------------

@pytest.mark.parametrize("call, level", [
    (lambda t: genomes.simulate_genomes_family(t, duplication=0.2), "genomes"),
    (lambda t: genomes.simulate_genomes_ordered(t, duplication=0.2), "genomes"),
    (lambda t: traits.simulate_continuous(t, rate=1.0), "traits"),
    (lambda t: traits.simulate_discrete(t, states=["a", "b"], switch=0.1), "traits"),
])
@pytest.mark.parametrize("wrong", [None, 42, ["a", "b"]])
def test_a_level_refuses_a_non_tree_by_name(call, level, wrong):
    with pytest.raises(TypeError, match=f"the {level} level runs on a species tree"):
        call(wrong)


def test_a_newick_string_is_told_to_parse_itself_first():
    # the mistake worth naming: it IS a tree to the person holding it, just not one yet
    with pytest.raises(TypeError, match="Newick string: parse it first"):
        genomes.simulate_genomes_family("((a:1,b:1):1);", duplication=0.2)


def test_a_path_is_told_to_read_the_file_first():
    with pytest.raises(TypeError, match="looks like a path: read the file first"):
        genomes.simulate_genomes_family("tree.nwk", duplication=0.2)


def test_a_tree_and_a_species_result_are_both_accepted(run):
    sp, _ = run
    assert as_tree(sp, level="genomes") is sp.complete_tree
    assert as_tree(sp.complete_tree, level="genomes") is sp.complete_tree


def test_the_sequence_level_still_names_what_it_wants(run):
    sp, _ = run
    with pytest.raises(TypeError, match="runs on a genome run"):
        sequences.simulate_sequences(sp, model=sequences.jc69(), length=10)


# --- a result says what it is, not what it holds ---------------------------

def test_a_result_reprs_as_a_summary_not_as_its_contents(run):
    sp, g = run
    s = sequences.simulate_sequences(g, model=sequences.jc69(), length=20, seed=1)
    t = traits.simulate_continuous(sp, rate=1.0, seed=1)
    for obj, opening in [(sp, "SpeciesResult("), (sp.complete_tree, "Tree("),
                         (g, "FamilyGenomesResult("), (s, "SequencesResult("),
                         (t, "TraitsResult(")]:
        text = repr(obj)
        assert text.startswith(opening)
        # the point: the dataclass repr of these runs to megabytes, and an interactive session gets
        # it just for typing the name
        assert len(text) < 200, f"{opening} repr is {len(text)} chars"
        assert "seed=" in text or obj is sp.complete_tree   # a tree has no seed of its own


def test_the_repr_reports_the_run_it_came_from(run):
    sp, g = run
    assert f"{sp.n_extant} extant tips" in repr(sp)
    assert f"{len(g.events)} events" in repr(g)


# --- the typing marker ------------------------------------------------------

def test_the_package_ships_a_py_typed_marker():
    # without it a type checker in someone else's project ignores every annotation in ZOMBI2 and
    # treats the whole package as Any (PEP 561) — the annotations are written, so they should count
    assert (pathlib.Path(genomes.__file__).parent.parent / "py.typed").exists()


# --- docstrings render on the docs site -------------------------------------

def test_no_docstring_carries_an_unrendered_sphinx_role():
    # `docs/reference/api.md` is generated from these docstrings by mkdocstrings, which renders
    # markdown and knows nothing about Sphinx's :func:`x` roles — so one left in a docstring
    # publishes verbatim, colons and backticks included. Backticks are what both readers want.
    import re

    root = pathlib.Path(genomes.__file__).parent.parent
    role = re.compile(r":(func|meth|class|attr|mod|data|obj|exc):`")
    offenders = [f"{p.relative_to(root)}:{i}" for p in sorted(root.rglob("*.py"))
                 for i, line in enumerate(p.read_text(encoding="utf-8").splitlines(), 1) if role.search(line)]
    assert not offenders, f"Sphinx roles left in {len(offenders)} place(s): {offenders[:5]}"


# --- one place decides how a node is written --------------------------------

def test_every_node_label_is_minted_through_the_one_helper():
    # `n<id>` is the most load-bearing token in the project — the species Newick, every event log,
    # the profile headers, the gene-tree labels, the FASTA records. It was minted by an inline
    # f-string in eighteen places, so changing how a node is written meant finding all eighteen.
    # It now goes through zombi2.tree.node_label, and this is what keeps it that way.
    import re

    root = pathlib.Path(genomes.__file__).parent.parent
    # `f"n{i}"`, and `f"...\tn{e.node}\t..."` too — but not the `\n{` of an escaped newline
    inline = re.compile(r'(?<![\\\w])n\{')
    offenders = [f"{p.relative_to(root)}:{i}" for p in sorted(root.rglob("*.py"))
                 for i, line in enumerate(p.read_text(encoding="utf-8").splitlines(), 1)
                 if inline.search(line) and p.name != "tree.py"]
    assert not offenders, f"a node label minted by hand in {offenders}"


def test_the_node_label_pair_round_trips():
    from zombi2.tree import node_from_label, node_label

    for i in (0, 7, 12345):
        assert node_from_label(node_label(i)) == i
    assert node_label(None) == ""
    assert node_from_label("12") == 12          # a log written before the columns carried their n


# --- the CLI never does anything Python cannot ----------------------------

#: Flags that say how a run is *driven or written*, not what is simulated: the CLI's own business.
#: A Python caller passes objects and calls `.write()`, so these have nothing to correspond to.
_PLUMBING = {"help", "params", "seed", "from", "write", "flat", "quiet", "force", "parallel",
             "stream", "version", "resolution", "kind", "name", "out", "tol", "o"}

#: Flags whose Python route is a **constructor** rather than a keyword on the entry point. A
#: substitution model's physical parameters belong to the model (`hky85(kappa=…)`), and a joint run's
#: driver is built before the call (`traits.discrete(states=…)`), so the capability is there, one
#: object further in. Every flag here spells its parameter *identically*; where it did not — `freqs`
#: for `--frequencies`, `rates` for `--gtr-rates`, `start` for `--trait-start` — the names were
#: brought together rather than mapped, because two names for one thing is the drift this test is
#: meant to catch and not a thing to encode in the test that catches it.
_VIA_CONSTRUCTOR = {
    "kappa": "zombi2.sequences.substitution_models:hky85",
    "frequencies": "zombi2.sequences.substitution_models:hky85",
    "exchangeabilities": "zombi2.sequences.substitution_models:gtr",
    "gamma-shape": "zombi2.sequences.substitution_models:SubstitutionModel.across_sites",
    "invariant": "zombi2.sequences.substitution_models:SubstitutionModel.across_sites",
    "rate-categories": "zombi2.sequences.substitution_models:SubstitutionModel.across_sites",
    "tip-fates": "zombi2.tree:read_newick",
    "duplication": "zombi2.genomes:family",
    "loss": "zombi2.genomes:family",
    "origination": "zombi2.genomes:family",
    "initial-families": "zombi2.genomes:family",
    "family-names": "zombi2.genomes:family",
    "states": "zombi2.traits:discrete",
    "switch": "zombi2.traits:discrete",
    "start": "zombi2.traits:discrete",
    "at-speciation": "zombi2.traits:discrete",
}

#: The one flag whose name deliberately differs, and why: it is `action="append"`, so each use names
#: ONE pulse while the parameter holds the list of them. Pluralising the flag would misname the
#: single use it is written for, and singularising the parameter would misname a list.
_PLURAL = {"mass-extinction": ("zombi2.species:simulate_species_tree", "mass_extinctions")}


def _params_of(path: str) -> set:
    import importlib
    import inspect

    module, _, attr = path.partition(":")
    obj = importlib.import_module(module)
    for part in attr.split("."):
        obj = getattr(obj, part)
    return set(inspect.signature(obj).parameters)


def test_every_cli_flag_is_reachable_from_python():
    """The command line must never be able to simulate something the Python API cannot.

    The two are one library with two front doors, and the CLI's flags are the API's keywords by
    design — so a flag with no Python route is a capability that exists only behind a subprocess.
    That is exactly what happened with the ordered resolution: `zombi2 genomes --resolution ordered`
    into `zombi2 sequences` worked, while the same two calls in Python raised, because the CLI's
    directory handoff quietly rebuilt a different result class on the way in. The documentation said
    the Python route worked, in four places.

    The reverse is not required and is not checked: plenty of Python is not worth a flag (a matrix,
    a callable, a correlated trait), and SPEC §9 says so.
    """
    import contextlib
    import inspect
    import io
    import re

    from zombi2 import genomes, joint, sequences, species, traits
    from zombi2.cli import main as cli_main

    entry_points = {
        "species": [species.simulate_species_tree],
        "genomes": [genomes.simulate_genomes_family, genomes.simulate_genomes_ordered,
                    genomes.simulate_genomes_nucleotide],
        "sequences": [sequences.simulate_sequences],
        "traits": [traits.simulate_continuous, traits.simulate_discrete],
        "joint": [joint.simulate_joint],
    }
    unreachable = []
    for command, fns in entry_points.items():
        buf = io.StringIO()
        with contextlib.suppress(SystemExit), contextlib.redirect_stdout(buf), \
                contextlib.redirect_stderr(buf):
            cli_main.main([command, "--help"])
        flags = {f[2:] for f in re.findall(r"(--[a-z][a-z0-9-]+)", buf.getvalue())}
        assert flags, f"the help scrape found no flags for {command}; this test would pass vacuously"
        direct = set()
        for fn in fns:
            direct |= set(inspect.signature(fn).parameters)
        for flag in sorted(flags - _PLUMBING):
            keyword = flag.replace("-", "_")
            if keyword in direct:
                continue
            provider = _VIA_CONSTRUCTOR.get(flag)
            if provider and keyword in _params_of(provider):
                continue
            plural = _PLURAL.get(flag)
            if plural and plural[1] in _params_of(plural[0]):
                continue
            unreachable.append(f"{command}: --{flag}")
    assert not unreachable, ("these CLI flags have no Python route — the command line can simulate "
                             "something the API cannot:\n  " + "\n  ".join(unreachable))


def test_the_sequence_level_takes_an_ordered_run(run):
    """`simulate_sequences` refused an `OrderedGenomesResult` while its own docstring, Chapter 7 and
    Appendix B all said "a family or ordered run". The gate tested the class rather than what the
    level needs, which is `gene_trees` and `complete_tree` — both of which an ordered result has.

    Only the CLI was unaffected, and only by accident: its directory handoff rebuilds a
    `FamilyGenomesResult` from `genome_events.tsv`, so the two front doors disagreed."""
    from zombi2.sequences.substitution_models import jc69

    sp, _ = run
    ordered = genomes.simulate_genomes_ordered(sp, duplication=0.2, loss=0.1, initial_families=5,
                                               inversion=0.3, inversion_extent=2, seed=1)
    result = sequences.simulate_sequences(ordered, model=jc69(), length=40, seed=1)
    assert result.alignments and result.unit == "family"
    assert set(result.alignments) == set(ordered.gene_trees)
