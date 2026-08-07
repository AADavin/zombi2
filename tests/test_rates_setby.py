"""`SetBy` — a driver that **replaces** a parameter's base rather than multiplying it.

The literature usually states a driven rate absolutely: *the loss rate is 1.0 in caves*, not *four
times a background nobody wrote down*. `ScaledBy` can only say the second, so saying the first meant
inventing a background and dividing by it. `SetBy` says it directly.

It is a `DrivenBy`, so every engine that resolves drivers resolves this one — the trajectory, the
mid-branch switches and the mapping checks are all the same machinery. One line in `Rate.effective`
asks a `SetBy` for the base and everything else for a factor.
"""

from __future__ import annotations

import pytest

from zombi2 import genomes, traits
from zombi2.rates import ScaledBy, SetBy
from zombi2.rates import scope
from zombi2.rates.rate import as_rate
from zombi2.species import simulate_species_tree


@pytest.fixture(scope="module")
def tree():
    return simulate_species_tree(birth=1.0, death=0.2, n_extant=16, seed=4).complete_tree


@pytest.fixture(scope="module")
def habitat(tree):
    return traits.simulate_discrete(tree, states=["cave", "surface"], switch=0.4, seed=2)


def _rate(spec, default=scope.PerCopy):
    return as_rate(spec, default_scope=default)


class TestItReplacesTheBase:

    def test_the_driver_supplies_the_number_itself(self):
        r = _rate(SetBy("h.tsv", {"cave": 1.0, "surface": 0.25}))
        assert r.effective(copies=1, drivers={"h.tsv": "cave"}) == pytest.approx(1.0)
        assert r.effective(copies=1, drivers={"h.tsv": "surface"}) == pytest.approx(0.25)

    def test_the_scope_still_applies(self):
        """`SetBy` replaces the base, not the *per what?*. A per-copy rate set to 1.0 is 1.0 per
        copy, so four copies make four — otherwise a replaced base would quietly become a total."""
        r = _rate(SetBy("h.tsv", {"cave": 1.0}))
        assert r.effective(copies=4, drivers={"h.tsv": "cave"}) == pytest.approx(4.0)

    def test_a_replaced_base_can_still_be_scaled(self):
        r = _rate(SetBy("h.tsv", {"cave": 1.0}) * ScaledBy("s.tsv", {"big": 3.0}))
        got = r.effective(copies=2, drivers={"h.tsv": "cave", "s.tsv": "big"})
        assert got == pytest.approx(6.0)          # 1.0 set, tripled, times two copies

    def test_a_written_base_is_ignored_and_so_is_refused(self):
        """Silently discarding a number someone wrote is the failure mode this grammar exists to
        avoid, so `0.25 * SetBy(...)` raises rather than dropping the 0.25."""
        with pytest.raises(TypeError, match="no base to write in front of it"):
            0.25 * SetBy("h.tsv", {"cave": 1.0})


class TestInARealRun:

    def test_an_absolute_statement_matches_the_same_model_stated_as_a_multiple(self, tree, habitat):
        """The correctness check that matters: `SetBy(driver, {...})` and `1.0 * ScaledBy(driver,
        {...})` are the same model written two ways, so they must be the same run, event for event."""
        def run(loss):
            g = genomes.simulate_genomes_family(tree, loss=loss, initial_families=40, seed=3)
            return [(e.time, e.kind, e.family) for e in g.events]

        assert run(SetBy(habitat, {"cave": 1.0, "surface": 0.05})) == \
            run(1.0 * ScaledBy(habitat, {"cave": 1.0, "surface": 0.05}))

    def test_it_drives_a_trait_rate_too(self, tree, habitat):
        run = traits.simulate_discrete(tree, states=["a", "b"],
                                       switch=SetBy(habitat, {"cave": 2.0, "surface": 0.1}), seed=1)
        assert run.events


class TestWhatItRefuses:

    def test_two_bases_on_one_rate(self):
        """Each claims to *be* the number, and no order of application is more right than the other,
        so it raises rather than letting whichever was written last win."""
        r = _rate(SetBy("h.tsv", {"cave": 1.0}) * SetBy("h.tsv", {"cave": 2.0}))
        with pytest.raises(ValueError, match="a base can only be replaced once"):
            r.check_one_base("loss")

    def test_two_bases_are_refused_by_the_engine(self, tree, habitat):
        with pytest.raises(ValueError, match="a base can only be replaced once"):
            genomes.simulate_genomes_family(
                tree, loss=SetBy(habitat, {"cave": 1.0}) * SetBy(habitat, {"cave": 2.0}),
                initial_families=5, seed=1)

    def test_one_base_and_any_number_of_factors_is_fine(self):
        r = _rate(SetBy("h.tsv", {"cave": 1.0}) * ScaledBy("a", {"x": 2.0}) * ScaledBy("b", {"y": 3.0}))
        r.check_one_base("loss")          # does not raise


def test_it_is_written_the_way_it_is_read():
    """The repr is what a run's log records and a reader pastes back into a flag, so it has to be an
    expression that reproduces the rate."""
    assert repr(SetBy("h.tsv", {"cave": 1.0})).startswith("SetBy('h.tsv'")

    from zombi2.rates.parse import parse_rate
    # a bare modifier parses to itself; only a product becomes a Rate
    assert isinstance(parse_rate('SetBy("h.tsv", {"cave": 1.0})'), SetBy)
    assert isinstance(parse_rate('SetBy("h.tsv", {"cave": 1.0}) * ScaledBy("s", {"x": 2.0})')
                      .modifiers[0], SetBy)
