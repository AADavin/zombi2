"""The grid spelling — `values` and `verbs`.

A value is two independent facts, what it is attached to and how its number is made, and a verb says
what a parameter does with that number. The six modifiers are six cells of that grid rather than six
mechanisms, so both spellings build the same object and a run is unchanged. What the grid buys is
that the *next* cell needs no new class and no invented name.
"""

from __future__ import annotations

import pytest

from zombi2 import genomes
from zombi2.rates import Drawn, Inherited, ScaledBy, SetBy, Time, Weights
from zombi2.rates import modifiers as mod
from zombi2.species import simulate_species_tree


class TestTheGridBuildsTodaysObjects:
    """Each cell is a name for something that already exists, so nothing downstream changes."""

    @pytest.mark.parametrize("grid, listed", [
        (Drawn(per="family", spread=0.5), mod.ByFamily(spread=0.5)),
        (Drawn(per="lineage", spread=0.3), mod.ByLineage(spread=0.3)),
        (Drawn(per="family", spread=0.5, dist="gamma"), mod.ByFamily(spread=0.5, dist="gamma")),
        (Inherited(per="lineage", spread=0.2), mod.FromParent(spread=0.2)),
        (Inherited(per="lineage", spread=0.2, bins=4), mod.FromParent(spread=0.2, bins=4)),
        (ScaledBy(Time(), {0: 1.0, 3: 0.3}), mod.OnTime({0: 1.0, 3: 0.3})),
        (ScaledBy("habitat.tsv", {"cave": 4.0}), mod.DrivenBy("habitat.tsv", {"cave": 4.0})),
        (Weights("habitat.tsv", {"cave": 4.0}), mod.DrivenBy("habitat.tsv", {"cave": 4.0})),
    ])
    def test_a_cell_is_the_modifier_it_names(self, grid, listed):
        assert grid == listed
        assert type(grid) is type(listed)


class TestARunIsUnchanged:
    """The point of a renaming that is not a rewrite: seed for seed, the same events."""

    def test_inherited_matches_fromparent(self):
        def run(m):
            return [(e.time, e.kind, e.node) for e in
                    simulate_species_tree(birth=1.0 * m, death=0.2, n_extant=40, seed=9).events]

        assert run(Inherited(per="lineage", spread=0.2)) == run(mod.FromParent(spread=0.2))

    def test_scaledby_time_matches_ontime(self):
        def run(m):
            return [(e.time, e.kind) for e in
                    simulate_species_tree(birth=1.0 * m, death=0.2, n_extant=40, seed=9).events]

        assert run(ScaledBy(Time(), {0: 1.0, 1: 0.3})) == run(mod.OnTime({0: 1.0, 1: 0.3}))

    def test_drawn_per_family_matches_byfamily(self):
        tree = simulate_species_tree(birth=1.0, death=0.2, n_extant=20, seed=4).complete_tree

        def run(m):
            g = genomes.simulate_genomes_family(tree, duplication=0.3, loss=0.3 * m,
                                                initial_families=30, seed=2)
            return [(e.time, e.kind, e.family) for e in g.events]

        assert run(Drawn(per="family", spread=0.6)) == run(mod.ByFamily(spread=0.6))


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
        with pytest.raises(ValueError, match="ByLineage"):
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
        it here is worth a line, because `ScaledBy(ByFamily(...))` reads plausibly and would
        otherwise be threaded as a driver that never resolves."""
        with pytest.raises(TypeError, match="already a factor"):
            ScaledBy(mod.ByFamily(spread=0.5))

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
