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
                 for i, line in enumerate(p.read_text().splitlines(), 1) if role.search(line)]
    assert not offenders, f"Sphinx roles left in {len(offenders)} place(s): {offenders[:5]}"
