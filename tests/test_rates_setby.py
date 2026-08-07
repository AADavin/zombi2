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

    def test_it_drives_a_continuous_trait_rate_too(self, tree, habitat):
        run = traits.simulate_continuous(
            tree, rate=SetBy(habitat, {"cave": 2.0, "surface": 0.1}), seed=1)
        assert run.node_values

    def test_a_level_that_cannot_replace_a_base_refuses_it(self, tree, habitat):
        """A `SetBy` is a `DrivenBy`, so a gate listing DrivenBy would let it in anywhere a driver
        goes — and four levels admitted it that way and could not honour it. A level has to name
        `SetBy` to accept one."""
        with pytest.raises(ValueError, match="does not support|does not read"):
            traits.simulate_discrete(tree, states=["a", "b"],
                                     switch=SetBy(habitat, {"cave": 2.0, "surface": 0.1}), seed=1)


class TestWhatItRefuses:

    def test_two_bases_on_one_rate(self):
        """Each claims to *be* the number, and no order of application is more right than the other,
        so it raises rather than letting whichever was written last win."""
        with pytest.raises(TypeError, match="a rate carries one SetBy"):
            SetBy("h.tsv", {"cave": 1.0}) * SetBy("h.tsv", {"cave": 2.0})

    def test_two_bases_are_refused_however_they_are_assembled(self):
        """The product refuses them, and `as_rate` refuses them again for anything that reaches it
        another way — the guard lives at the choke point every level already calls, because a rule
        enforced by whoever remembers to call it is a rule three levels did not have."""
        from zombi2.rates.rate import Rate
        from zombi2.rates.scope import PerCopy

        smuggled = Rate(1.0, PerCopy(1.0),
                        (SetBy("h", {"c": 1.0}), SetBy("h", {"c": 2.0})))
        with pytest.raises(ValueError, match="a base can only be replaced once"):
            as_rate(smuggled, default_scope=PerCopy, label="loss")

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


class TestTheHolesAnAdversarialReviewFound:
    """Every one of these shipped broken for a few hours. `SetBy` is a `DrivenBy`, so it passed every
    gate that listed `DrivenBy` — including four levels that could not honour a replaced base — and
    the "no base in front of it" guard only ever saw the operand immediately to its left."""

    @pytest.mark.parametrize("build", [
        pytest.param(lambda: 0.25 * SetBy("h", {"c": 1.0}), id="number * SetBy"),
        pytest.param(lambda: 0.25 * ScaledBy("s", {"b": 1.0}) * SetBy("h", {"c": 1.0}),
                     id="number * ScaledBy * SetBy"),
        pytest.param(lambda: scope.PerCopy(0.25) * SetBy("h", {"c": 1.0}), id="scope * SetBy"),
        pytest.param(lambda: SetBy("h", {"c": 1.0}) * SetBy("h", {"c": 2.0}), id="SetBy * SetBy"),
    ])
    def test_a_base_can_never_be_written_where_it_would_be_discarded(self, build):
        """Only the first of these was caught. The others each built a Rate whose written base was
        then silently overwritten by the driver's number."""
        with pytest.raises(TypeError):
            build()

    @pytest.mark.parametrize("level", ["discrete-trait", "nucleotide", "sequences"])
    def test_a_level_that_cannot_replace_a_base_refuses_rather_than_mangling_it(self, level, tree):
        """Three of these took two `SetBy` and ran the last one written; the sequence level, which
        reads its modifiers itself, multiplied them together instead. Both are the silence the whole
        declaration mechanism exists to prevent."""
        one = SetBy("h.tsv", {"cave": 1.0, "surface": 0.5})
        if level == "discrete-trait":
            with pytest.raises(ValueError, match="does not support"):
                traits.simulate_discrete(tree, states=["a", "b"], switch=one, seed=1)
        elif level == "nucleotide":
            with pytest.raises(ValueError, match="does not support"):
                genomes.simulate_genomes_nucleotide(tree, loss=one, root_length=1000, seed=1)
        else:
            from zombi2.sequences import jc69, simulate_sequences
            g = genomes.simulate_genomes_family(tree, duplication=0.1, loss=0.1,
                                                initial_families=3, seed=5)
            with pytest.raises(ValueError, match="does not read"):
                simulate_sequences(g, model=jc69(), length=40, substitution=one, seed=1)

    def test_an_extent_has_no_base_to_replace(self):
        """An extent is already an absolute size drawn from a distribution. A `SetBy` there was
        admitted and applied as a multiplier, which is a different model wearing the same words."""
        from zombi2.rates.extent import as_extent
        with pytest.raises(ValueError, match="an extent cannot be SetBy"):
            as_extent(SetBy("h", {"c": 5.0}))

    def test_a_clade_driven_rate_is_written_so_it_can_be_pasted_back(self):
        """It recorded `DrivenBy('<Clade>', ...)` — which parses, as a *filename*, so the log looked
        reproducible and was not. A driver that can write itself now does; one that cannot records an
        unquoted placeholder that fails loudly."""
        from zombi2.rates import Clade
        from zombi2.rates.parse import parse_rate

        written = repr(ScaledBy(Clade({"fast": ["n1", "n2"]}), {"fast": 3.0}))
        assert "<Clade>" not in written
        assert parse_rate(f"0.2 * {written}").modifiers[0].driver == Clade({"fast": ["n1", "n2"]})

    def test_a_driver_that_cannot_be_written_says_so_rather_than_looking_like_a_file(self, habitat):
        assert repr(ScaledBy(habitat, {"cave": 2.0})).startswith("DrivenBy(<TraitsResult>")
