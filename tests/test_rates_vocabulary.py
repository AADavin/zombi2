"""The grid spelling — the drivers, the laws and the verbs.

A value is two independent facts: what it varies among, and how its number is made. A verb says what
a parameter does with that number. Keeping the two apart is what lets one class cover a whole row,
so the *next* cell needs no new class and no invented name.
"""

from __future__ import annotations

import re

import pytest

from zombi2 import genomes
from zombi2.params import (Drift, LogNormal, PerCopy, PerLineage, Random, Recipients, Time,
)
from zombi2.params import driver as drv
from zombi2.params import law as law
from zombi2.params import connection as conn
from zombi2.species import simulate_species_tree


def test_what_you_can_import_and_what_you_can_write_are_one_surface():
    """The package's own claim, checked. Five names were writable in a ``--birth`` flag and absent
    from the package, so a user following the manual's `Between` example had to discover
    ``zombi2.params.mapping`` — the two-surface split the comment above ``__all__`` says does not
    exist.

    Two entries are deliberately exported and not writable, and both for a reason the comment
    states: `Curve` takes a function, which no flag can carry, and ``UNITS`` is data rather than a
    name you call."""
    import zombi2.params
    from zombi2.params.parse import _NAMES

    assert not set(_NAMES) - set(zombi2.params.__all__), "writable in a flag, missing from the package"
    assert set(zombi2.params.__all__) - set(_NAMES) == {"Curve", "UNITS"}


class TestTheGridBuildsWhatTheEnginesRun:
    """A value is a unit plus a law, and the verbs build the objects the engines already dispatch
    on — so the grammar is a spelling, not a second machine."""

    def test_a_drawn_cell_reports_its_kind_and_unit(self):
        for unit in ("families", "lineages", "chromosomes", "sites"):
            assert Random(unit, LogNormal(0.0, 0.5)).reads == ("drawn", unit)

    def test_an_inherited_cell_reports_its_kind_and_unit(self):
        assert Random("lineages", Drift(LogNormal(0.0, 0.2))).reads == ("inherited", "lineages")

    @pytest.mark.parametrize("written, expected", [
        (PerLineage(0.5).changing_at({0: 1.0, 3: 0.3}), drv.OnTime),
        (PerLineage().set_by(Time(), {0: 0.5, 3: 0.15}), drv.OnTime),
        (PerCopy(0.25).scaled_by("habitat.tsv", {"cave": 4.0}), conn.Driven),
        (Recipients().weighted_by("habitat.tsv", {"cave": 4.0}), conn.Driven),
        (PerCopy().set_by("habitat.tsv", {"cave": 1.0}), conn.SetBy),
    ])
    def test_a_verb_builds_the_modifier_the_engines_dispatch_on(self, written, expected):
        built = (written.weights if hasattr(written, "weights") else written.modifiers)[0]
        assert type(built) is expected

    def test_the_law_is_what_makes_two_cells_the_same(self):
        """A cell is its unit and its law — nothing about how either was typed."""
        assert Random("families", LogNormal(0.0, 0.5)) == Random("families", LogNormal(0.0, 0.5))
        assert Random("families", LogNormal(0.0, 0.5)) != Random("families", LogNormal(0.0, 0.6))
        assert Random("families", LogNormal(0.0, 0.5)) != Random("lineages", LogNormal(0.0, 0.5))

    def test_a_cell_without_a_law_is_refused(self):
        """There is no default shape: an unstated distribution would be a model nobody wrote."""
        with pytest.raises(ValueError, match="needs the law its value follows"):
            Random("families")
        with pytest.raises(ValueError, match="needs the distribution of its per-split step"):
            Random("lineages", Drift())


class TestACellNobodyBuiltRefusesByName:
    """The grid makes the empty cells visible, so the refusal can say *which* cell and why, rather
    than reading as an unknown name."""

    def test_a_cell_nobody_built_constructs_and_the_engine_refuses_it(self):
        """The grid is real, so any unit builds a value — no new class, no invented name. What
        decides whether a run happens is the *level*, which declares the cells it can carry a number
        for, and refuses the rest naming which cell it was rather than which class."""
        cell = Random("chromosomes", LogNormal(0.0, 0.5))
        assert cell.reads == ("drawn", "chromosomes")

        with pytest.raises(ValueError, match=re.escape("varying_among('chromosomes', ...)")):
            simulate_species_tree(birth=PerLineage(1.0).varying_among(cell), death=0.2,
                                  n_extant=5, seed=1)

    def test_a_cell_the_wrong_level_carries_is_named_by_its_cell(self):
        """One class covers a whole row, so "carries a Random" would be true and useless — the
        question is always *among what*. The cell is named by the expression that writes it, so the
        refusal and the list of what the level does take are in one vocabulary."""
        with pytest.raises(ValueError, match=re.escape("varying_among('lineages', ...)")):
            genomes.simulate_genomes_family(
                simulate_species_tree(birth=1.0, death=0.2, n_extant=6, seed=1).complete_tree,
                loss=PerCopy(0.2).varying_among("lineages", LogNormal(0.0, 0.5)),
                initial_families=3, seed=1)

    def test_a_unit_that_is_not_a_unit_lists_the_units(self):
        with pytest.raises(ValueError, match="unknown unit 'banana'"):
            Random("banana", LogNormal(0.0, 0.5))

    def test_a_smooth_curve_of_time_is_refused_rather_than_approximated(self):
        """It is a real model — it makes the rate vary continuously, which needs the engine to
        integrate rather than sample. Refusing beats a quiet approximation."""
        with pytest.raises(ValueError, match="smooth function of|takes a schedule"):
            PerLineage(1.0).changing_at(lambda t: 2.0 ** -t)
        with pytest.raises(ValueError, match="smooth function of|takes a schedule"):
            PerLineage().set_by(Time(), lambda t: 2.0 ** -t)

    def test_setby_is_now_implemented(self):
        """It refused, with its reason, until the engines could take a replaced base. They can, so
        it builds — see `test_rates_setby.py` for what it does."""
        built = PerCopy().set_by("habitat.tsv", {"cave": 1.0}).modifiers[0]
        assert built.reads == ("driven", "lineages")


class TestTheVerbsAreForValuesNotFactors:

    def test_wrapping_a_draw_in_a_verb_is_refused(self):
        """A draw is already a factor; a verb is for a value a mapping has to turn into one. Catching
        it here is worth a line, because ``scaled_by(Random('families', ...))`` reads plausibly and
        would otherwise be threaded as a driver that never resolves."""
        with pytest.raises(TypeError, match="already a factor"):
            PerCopy(0.25).scaled_by(Random("families", LogNormal(0.0, 0.5)), {"x": 1.0})

    def test_a_verb_without_a_mapping_says_what_a_mapping_is(self):
        with pytest.raises(ValueError, match="needs a mapping"):
            PerCopy(0.25).scaled_by("habitat.tsv")
        with pytest.raises(ValueError, match="needs a mapping"):
            Recipients().weighted_by("habitat.tsv")


class TestTheTwoShortcutsAreTheOnlySpellingForTheirDriver:
    """`Random` and `Time` are the only built-in drivers, and each has one verb of its own.
    `scaled_by` refuses both and names it, so there is exactly one way to write each (SPEC §4.1)."""

    def test_scaled_by_refuses_the_clock_and_names_changing_at(self):
        with pytest.raises(ValueError, match=r"changing_at\(\{0: 1.0, 3: 0.3\}\)"):
            PerLineage(0.5).scaled_by(Time(), {0: 1.0, 3: 0.3})

    def test_the_clock_refusal_also_names_the_other_reading(self):
        with pytest.raises(ValueError, match=r"set_by\(Time\(\)"):
            PerLineage(0.5).scaled_by(Time(), {0: 1.0, 3: 0.3})

    def test_scaled_by_refuses_a_random_and_names_varying_among(self):
        with pytest.raises(TypeError, match="varying_among"):
            PerCopy(0.25).scaled_by(Random("families", LogNormal(0.0, 0.5)), {"x": 1.0})

    def test_weighted_by_refuses_the_clock_too(self):
        with pytest.raises(ValueError, match="changing_at"):
            Recipients().weighted_by(Time(), {0: 1.0, 3: 0.3})


class TestARetiredSpellingIsAnsweredInPythonToo:
    """One table, two surfaces. The text form has named the replacement for every retired spelling
    since the grammar changed; Python said only that the name was absent, which is the one thing the
    reader porting a script already knew."""

    @pytest.mark.parametrize("name, replacement", [
        ("ScaledBy", "the verbs are methods on the parameter"),
        ("SetBy", "written from the bare scope"),
        ("Weights", r"Recipients\(\).weighted_by"),
        ("Drawn", r"varying_among\('families'"),
        ("Inherited", r"Drift\(LogNormal"),
        ("OnTime", "changing_at"),
        ("OnTotalDiversity", r"scaled_by\(TotalDiversity\(cap=100\)\)"),
        ("DrivenBy", "write the verb that says what the number does"),
        ("ByFamily", r"varying_among\('families'"),
        ("ByLineage", r"varying_among\('lineages'"),
        ("FromParent", r"Drift\(LogNormal"),
    ])
    def test_an_import_of_a_retired_name_names_its_replacement(self, name, replacement):
        """An `ImportError` rather than an `AttributeError`, because ``from … import`` discards the
        second and substitutes its own generic "cannot import name" — losing the sentence exactly
        where a port hits it first."""
        with pytest.raises(ImportError, match=replacement):
            __import__("zombi2.params", fromlist=[name]).__getattr__(name)

    def test_a_name_that_was_never_in_the_grammar_is_still_an_attribute_error(self):
        with pytest.raises(AttributeError, match="has no attribute 'Wobble'"):
            import zombi2.params
            zombi2.params.Wobble

    def test_the_retired_names_are_the_same_list_the_text_form_reads(self):
        from zombi2.params import RETIRED
        from zombi2.params.parse import RETIRED as PARSED
        assert RETIRED is PARSED


def test_one_object_read_twice_is_still_one_draw_in_the_new_spelling():
    """The sharing rule is a property of the object, so it survives the renaming untouched."""
    tree = simulate_species_tree(birth=1.0, death=0.2, n_extant=16, seed=4).complete_tree
    speed = Random("families", LogNormal(0.0, 0.6))

    shared = genomes.simulate_genomes_family(
        tree, duplication=PerCopy(0.3).varying_among(speed),
        loss=PerCopy(0.3).varying_among(speed), initial_families=20, seed=2)
    apart = genomes.simulate_genomes_family(
        tree, duplication=PerCopy(0.3).varying_among("families", LogNormal(0.0, 0.6)),
        loss=PerCopy(0.3).varying_among("families", LogNormal(0.0, 0.6)),
        initial_families=20, seed=2)

    assert [e.time for e in shared.events] != [e.time for e in apart.events]


def test_a_named_random_carries_its_law_so_a_second_argument_is_refused():
    """Sharing one object is what makes two rates share one draw — passing a law beside it would
    read as building a second."""
    speed = Random("families", LogNormal(0.0, 0.6))
    with pytest.raises(TypeError, match="already carries its law"):
        PerCopy(0.3).varying_among(speed, LogNormal(0.0, 0.6))
