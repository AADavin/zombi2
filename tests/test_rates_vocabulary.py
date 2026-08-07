"""The grid spelling — `values` and `verbs`.

A value is two independent facts: what it is attached to, and how its number is made. A verb says
what a parameter does with that number. Keeping the two apart is what lets one class cover a whole
row, so the *next* cell needs no new class and no invented name.
"""

from __future__ import annotations

import pytest

from zombi2 import genomes
from zombi2.rates import Drawn, Inherited, ScaledBy, SetBy, Time, Weights
from zombi2.rates import modifiers as mod
from zombi2.species import simulate_species_tree


class TestTheGridBuildsWhatTheEnginesRun:
    """A value is a unit plus a way of making its number, and the verbs build the objects the engines
    already dispatch on — so the grammar is a spelling, not a second machine."""

    def test_a_drawn_cell_reports_its_kind_and_unit(self):
        for unit in ("family", "lineage", "chromosome", "site"):
            assert Drawn(per=unit, spread=0.5).reads == ("drawn", unit)

    def test_an_inherited_cell_reports_its_kind_and_unit(self):
        assert Inherited(per="lineage", spread=0.2).reads == ("inherited", "lineage")

    @pytest.mark.parametrize("built, expected", [
        (ScaledBy(Time(), {0: 1.0, 3: 0.3}), mod.OnTime),
        (ScaledBy("habitat.tsv", {"cave": 4.0}), mod.DrivenBy),
        (Weights("habitat.tsv", {"cave": 4.0}), mod.DrivenBy),
        (SetBy("habitat.tsv", {"cave": 1.0}), mod.SetBy),
    ])
    def test_a_verb_builds_the_modifier_the_engines_dispatch_on(self, built, expected):
        assert type(built) is expected

    def test_spread_is_a_lognormal_and_dist_is_anything_else(self):
        """``spread=σ`` *is* ``dist=LogNormal(0.0, σ)``; the two spellings are one modifier."""
        from zombi2.rates.distributions import LogNormal
        assert Drawn(per="family", spread=0.5) == Drawn(per="family", dist=LogNormal(0.0, 0.5))

    def test_giving_both_or_neither_is_refused(self):
        from zombi2.rates.distributions import Gamma
        for kwargs in ({}, {"spread": 0.5, "dist": Gamma(2.0, 0.5)}):
            with pytest.raises(ValueError, match="a spread or a dist"):
                Drawn(per="family", **kwargs)


class TestACellNobodyBuiltRefusesByName:
    """The grid makes the empty cells visible, so the refusal can say *which* cell and why, rather
    than reading as an unknown name."""

    def test_a_cell_nobody_built_constructs_and_the_engine_refuses_it(self):
        """The grid is real, so any unit builds a value — no new class, no invented name. What
        decides whether a run happens is the *level*, which declares the cells it can carry a number
        for, and refuses the rest naming which cell it was rather than which class."""
        cell = Drawn(per="chromosome", spread=0.5)
        assert cell.reads == ("drawn", "chromosome")

        with pytest.raises(ValueError, match="drawn per chromosome"):
            simulate_species_tree(birth=1.0 * cell, death=0.2, n_extant=5, seed=1)

    def test_a_cell_the_wrong_level_carries_is_named_by_its_cell(self):
        """`Drawn` covers a whole row, so "carries Drawn" would be true and useless — the question is
        always *per what*."""
        with pytest.raises(ValueError, match="drawn per lineage"):
            genomes.simulate_genomes_family(
                simulate_species_tree(birth=1.0, death=0.2, n_extant=6, seed=1).complete_tree,
                loss=0.2 * Drawn(per="lineage", spread=0.5), initial_families=3, seed=1)

    def test_a_unit_that_is_not_a_unit_lists_the_units(self):
        with pytest.raises(ValueError, match="unknown unit 'banana'"):
            Drawn(per="banana", spread=0.5)

    def test_a_smooth_curve_of_time_is_refused_rather_than_approximated(self):
        """It is a real model — it makes the rate vary continuously, which needs the engine to
        integrate rather than sample. Refusing beats a quiet approximation."""
        with pytest.raises(ValueError, match="smooth function of|takes a schedule"):
            ScaledBy(Time(), lambda t: 2.0 ** -t)

    def test_setby_is_now_implemented(self):
        """It refused, with its reason, until the engines could take a replaced base. They can, so
        it builds — see `test_rates_setby.py` for what it does."""
        assert SetBy("habitat.tsv", {"cave": 1.0}).reads == ("driven", "lineage")


class TestTheVerbsAreForValuesNotFactors:

    def test_wrapping_a_draw_in_a_verb_is_refused(self):
        """A draw is already a factor; a verb is for a value a mapping has to turn into one. Catching
        it here is worth a line, because `ScaledBy(Drawn(per='family', ...))` reads plausibly and would
        otherwise be threaded as a driver that never resolves."""
        with pytest.raises(TypeError, match="already a factor"):
            ScaledBy(Drawn(per="family", spread=0.5))

    def test_a_verb_without_a_mapping_says_what_a_mapping_is(self):
        with pytest.raises(ValueError, match="needs a mapping"):
            ScaledBy("habitat.tsv")
        with pytest.raises(ValueError, match="needs a mapping"):
            Weights("habitat.tsv")


def test_one_object_read_twice_is_still_one_draw_in_the_new_spelling():
    """The sharing rule is a property of the object, so it survives the renaming untouched."""
    tree = simulate_species_tree(birth=1.0, death=0.2, n_extant=16, seed=4).complete_tree
    speed = Drawn(per="family", spread=0.6)

    shared = genomes.simulate_genomes_family(tree, duplication=0.3 * speed, loss=0.3 * speed,
                                             initial_families=20, seed=2)
    apart = genomes.simulate_genomes_family(
        tree, duplication=0.3 * Drawn(per="family", spread=0.6),
        loss=0.3 * Drawn(per="family", spread=0.6), initial_families=20, seed=2)

    assert [e.time for e in shared.events] != [e.time for e in apart.events]
