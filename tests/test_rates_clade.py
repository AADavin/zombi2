"""`Clade` — a value read off the tree itself, so a rate can depend on where a lineage sits.

Every other driver was grown by another level and handed over. A clade is a fact about the tree the
run is already walking, so there is nothing to grow first — which makes it the cheapest driver there
is, and the one that needed no engine change: it answers the same `as_driver_trajectory` protocol a
genome's `presence()` does.
"""

from __future__ import annotations

import collections
import math

import pytest

from zombi2 import species
from zombi2 import genomes, traits
from zombi2.genomes._transfer import resolve_groups
from zombi2.params import Clade, PerCopy, PerLineage, Table
from zombi2.params.parse import parse_rate
from zombi2.params.conditioned import resolve_driver
from zombi2.genomes import simulate_genomes_family
from zombi2.species import simulate_species_tree


@pytest.fixture(scope="module")
def tree():
    return simulate_species_tree(birth=1.0, death=0.2, n_extant=16, seed=4).complete_tree


@pytest.fixture(scope="module")
def halves(tree):
    """The root's two daughters, so the clades are disjoint and cover everything but the root."""
    a, b = tree.nodes[tree.root].children
    return Clade({"fast": a, "slow": b})


class TestTheTrajectory:

    def test_it_agrees_with_the_painting_on_every_lineage(self, tree, halves):
        traj = resolve_driver(halves, tree)
        painted = resolve_groups(tree, halves.groups)
        for i, node in tree.nodes.items():
            assert traj.value(i, node.birth_time) == painted[i]

    def test_a_lineage_in_no_named_clade_is_rest(self, tree, halves):
        """The root sits above both daughters, so it belongs to neither. ``rest`` is a state like any
        other and a mapping may name it."""
        assert resolve_driver(halves, tree).value(tree.root, 0.0) == "rest"
        assert "rest" in resolve_driver(halves, tree).states()

    def test_membership_never_switches_along_a_branch(self, tree, halves):
        """Which is why this driver costs the Gillespie nothing: a driver that switches mid-branch
        forces the race to stop at each switch, and a clade never switches."""
        traj = resolve_driver(halves, tree)
        for i, node in tree.nodes.items():
            assert traj.next_change(i, node.birth_time) == float("inf")


class TestItDrivesARate:

    def _losses(self, tree, halves, loss):
        painted = resolve_groups(tree, halves.groups)
        run = genomes.simulate_genomes_family(tree, loss=loss, initial_families=60, seed=3)
        return collections.Counter(painted[e.lineage] for e in run.edges if e.kind == "loss")

    def test_a_zero_weight_clade_loses_nothing(self, tree, halves):
        """The decisive test, because a ratio is confounded by copy count: set one clade's factor to
        zero and it must record no losses at all, while the other is untouched."""
        flat = self._losses(tree, halves, 0.3)
        assert flat["fast"] > 0 and flat["slow"] > 0

        no_slow = self._losses(tree, halves, PerCopy(0.3).scaled_by(halves, {"fast": 1.0, "slow": 0.0,
                                                       "rest": 1.0}))
        assert no_slow["slow"] == 0
        assert no_slow["fast"] > 0

        no_fast = self._losses(tree, halves, PerCopy(0.3).scaled_by(halves, {"fast": 0.0, "slow": 1.0,
                                                       "rest": 1.0}))
        assert no_fast["fast"] == 0
        assert no_fast["slow"] > 0

    def test_it_drives_a_trait_too(self, tree, halves):
        """Nothing about it is genome-specific — any level that reads a driver reads this one.

        Asserting that the run produced events would pass with the clade doing nothing at all, so
        the mapping switches one group off: no lineage of the "slow" clade may switch, and some
        lineage outside it must."""
        painted = resolve_groups(tree, halves.groups)
        slow = {i for i, label in painted.items() if label == "slow"}
        run = traits.simulate_discrete(
            tree, states=["a", "b"],
            switch=PerLineage(0.2).scaled_by(halves, {"fast": 4.0, "slow": 0.0, "rest": 1.0}), seed=1)

        switches = [c for c in run.events if c.kind != "initial"]
        assert switches, "nothing switched anywhere, so the zero proves nothing"
        assert not [c for c in switches if c.lineage in slow], "a clade weighted 0.0 still switched"


class TestWhatItRefuses:

    def test_a_growing_tree_has_no_clades(self, tree, halves):
        """A joint run grows the tree as it goes, so a clade is not defined while the run is
        happening — the refusal is about the model, not the implementation."""
        with pytest.raises(TypeError, match="live level"):
            from zombi2.joint import simulate
            simulate(species.birth_death(birth=PerLineage(1.0).scaled_by(halves, {"fast": 2.0}), death=0.2, n_extant=10), traits.discrete(states=["x", "y"], switch=0.1), seed=1)

    def test_overlapping_clades_are_refused(self, tree):
        with pytest.raises(ValueError, match="disjoint"):
            resolve_driver(Clade({"outer": tree.root, "inner": tree.nodes[tree.root].children[0]}),
                           tree)

    def test_it_needs_a_non_empty_mapping_of_labels(self):
        with pytest.raises(ValueError, match="non-empty"):
            Clade({})
        with pytest.raises(ValueError, match="labels must be non-empty strings"):
            Clade({"": 3})


# --- a clade AND a time window: a scheduled Table entry ------------------------------------------
#
# Chaining two verbs cannot express this: scaled_by(clade, {...}).changing_at({...}) multiplies two
# factors that each apply to every lineage, so the time window lands on the whole tree. A Table entry
# written as a schedule scopes the factor to the state AND to the time.

def test_a_scheduled_entry_reads_its_state_and_the_clock():
    t = Table({"endo": {0: 1.0, 6.0: 20.0}, "rest": 1.0})
    assert [t.multiplier("endo", time=x) for x in (0.0, 5.9, 6.0, 9.0)] == [1.0, 1.0, 20.0, 20.0]
    assert [t.multiplier("rest", time=x) for x in (0.0, 9.0)] == [1.0, 1.0]
    # the breakpoint has to reach the engine's horizon or the Gillespie steps straight past it
    assert t.next_change(0.0) == 6.0 and t.next_change(6.0) == math.inf
    assert Table({"a": 1.0}).next_change(0.0) == math.inf          # an ordinary table never moves


def test_a_scheduled_entry_survives_the_written_form():
    # the promise the whole grammar rests on: one notation across Python, the flag and the TOML file
    written = "PerCopy(0.02).scaled_by(Clade({'e': ['n1']}), {'e': {0: 1.0, 6.0: 20.0}, 'rest': 1.0})"
    rate = parse_rate(written)
    assert repr(rate) == ("PerCopy(0.02).scaled_by(Clade({'e': ['n1']}), "
                          "Table({'e': {0.0: 1.0, 6.0: 20.0}, 'rest': 1.0}))")
    assert parse_rate(repr(rate)) == rate                          # and it parses back to itself


def test_the_clade_and_the_window_both_bite():
    # the model itself: x20 inside this clade, and only after t=6. Both halves have to hold — a
    # schedule that leaked outside the clade, or a clade factor that ignored the clock, would show
    # up as losses in the wrong box.
    sp = simulate_species_tree(birth=0.55, death=0.15, total_time=8.0, seed=85)
    clade = Clade({"endo": ["n76", "n112"]})
    inside = set(clade.resolve(sp)["endo"])
    g = simulate_genomes_family(
        sp, initial_families=300,
        loss=PerCopy(0.02).scaled_by(clade, {"endo": {0: 1.0, 6.0: 20.0}, "rest": 1.0}), seed=4)
    box = collections.Counter((e.lineage in inside, e.time >= 6.0)
                              for e in g.edges if e.kind == "loss")
    # inside the clade the rate is x20 after t=6, so the late losses must dwarf the early ones ...
    assert box[(True, True)] > 20 * box[(True, False)]
    # ... while outside it nothing changed at t=6: 6 time units before against 2 after, on a growing
    # tree, so the two boxes stay within a factor of a few of each other rather than 20-fold apart
    assert box[(False, False)] > box[(False, True)] > 0.2 * box[(False, False)]
