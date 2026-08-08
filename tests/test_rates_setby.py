"""`SetBy` — a driver that **replaces** a parameter's base rather than multiplying it.

The literature usually states a driven rate absolutely: *the loss rate is 1.0 in caves*, not *four
times a background nobody wrote down*. `ScaledBy` can only say the second, so saying the first meant
inventing a background and dividing by it. `SetBy` says it directly.

It is a `Driven`, so every engine that resolves drivers resolves this one — the trajectory, the
mid-branch switches and the mapping checks are all the same machinery. One line in `Rate.effective`
asks a `SetBy` for the base and everything else for a factor.
"""

from __future__ import annotations

import pytest

from zombi2 import genomes, traits
from zombi2.rates import ScaledBy, SetBy, Weights
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
        """`node_values` is non-empty whatever the rate is, so the assertion has to be about the
        number `SetBy` supplied: the same model written as a multiple of 1.0 must give the same
        run, and a different set of numbers a different one."""
        driven = traits.simulate_continuous(
            tree, rate=SetBy(habitat, {"cave": 2.0, "surface": 0.1}), seed=1)
        same = traits.simulate_continuous(
            tree, rate=1.0 * ScaledBy(habitat, {"cave": 2.0, "surface": 0.1}), seed=1)
        other = traits.simulate_continuous(
            tree, rate=SetBy(habitat, {"cave": 0.1, "surface": 2.0}), seed=1)

        assert driven.node_values == pytest.approx(same.node_values)
        assert driven.node_values != pytest.approx(other.node_values)

    def test_a_level_that_cannot_replace_a_base_refuses_it(self, tree, habitat):
        """A `SetBy` is a `Driven`, so a gate listing Driven would let it in anywhere a driver
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
    """Every one of these shipped broken for a few hours. `SetBy` is a `Driven`, so it passed every
    gate that listed `Driven` — including four levels that could not honour a replaced base — and
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
        """It recorded `ScaledBy('<Clade>', ...)` — which parses, as a *filename*, so the log looked
        reproducible and was not. A driver that can write itself now does; one that cannot records an
        unquoted placeholder that fails loudly."""
        from zombi2.rates import Clade
        from zombi2.rates.parse import parse_rate

        written = repr(ScaledBy(Clade({"fast": ["n1", "n2"]}), {"fast": 3.0}))
        assert "<Clade>" not in written
        assert parse_rate(f"0.2 * {written}").modifiers[0].driver == Clade({"fast": ["n1", "n2"]})

    def test_a_driver_that_cannot_be_written_says_so_rather_than_looking_like_a_file(self, habitat):
        assert repr(ScaledBy(habitat, {"cave": 2.0})).startswith("ScaledBy(<TraitsResult>")

    def test_a_replaced_base_is_written_so_it_can_be_pasted_back(self):
        """A run records its rates in the written form, and a `SetBy` rate was recorded with a base
        in front of it — `1.0 * SetBy(...)`, or `PerCopy(1.0) * SetBy(...)` — which is exactly the
        spelling `SetBy` refuses. So every log naming a replaced base named a rate that would not
        parse, and the record of the model was not one you could run again.

        The scope goes with the base, and loses nothing: a scope cannot be written round a `SetBy`
        (there is no base for it to wrap), so a `SetBy` rate always carries its level's default, and
        reading the text back at that level restores it."""
        from zombi2.rates.parse import parse_rate, written_form

        alone = SetBy("h.tsv", {"cave": 1.0})
        assert written_form(alone) == repr(alone)
        assert parse_rate(written_form(alone)) == alone

        both = _rate(SetBy("h.tsv", {"cave": 1.0}) * ScaledBy("g.tsv", {"x": 2.0}))
        text = written_form(both)
        assert text.startswith("SetBy(")                    # written first: nothing may precede it
        assert parse_rate(text).modifiers == both.modifiers

    def test_a_choice_has_no_base_to_replace_either(self, tree, habitat):
        """The same hole as the extent above, on the other kind of target. `transfer_to` weights the
        candidate recipients against each other, so there is no base for `SetBy` to replace and the
        word means nothing there — but `SetBy` is a `Driven`, so the check that admits a driven
        `transfer_to` admitted it, and the run went ahead treating it as an ordinary weighting."""
        with pytest.raises(ValueError, match="transfer_to cannot be SetBy"):
            genomes.simulate_genomes_family(
                tree, transfer=0.5, initial_families=6, seed=2,
                transfer_to=SetBy(habitat, {"cave": 3.0, "surface": 1.0}))

        # the same numbers, spelled as what they are, still run
        genomes.simulate_genomes_family(
            tree, transfer=0.5, initial_families=6, seed=2,
            transfer_to=Weights(habitat, {"cave": 3.0, "surface": 1.0}))


def test_a_modifier_of_your_own_cannot_vouch_for_a_replaced_base():
    """`implemented_for` lets a third-party modifier declare which engines read it, and it may
    promise that for a factor it *computes*. Replacing a base is not that: it is a capability three
    levels have and four do not, so a `SetBy` subclass admitted through the hatch would be honoured
    nowhere. The same reason a carried value cannot go through it."""
    from zombi2.rates.modifiers import is_implemented, matches_declared

    class Mine(SetBy):
        implemented_for = ("species", "genomes.family", "traits.discrete")

    m = Mine("h.tsv", {"cave": 1.0})
    for engine in Mine.implemented_for:
        assert not is_implemented(m, (), engine), f"the hatch vouched for a SetBy at {engine}"

    # a level that names SetBy still takes it — the hatch is what closed, not the gate
    assert matches_declared(SetBy("h.tsv", {"cave": 1.0}), (SetBy,))


def test_the_verb_has_to_match_what_it_is_attached_to(tree, habitat):
    """Retiring `DrivenBy` made this checkable. `ScaledBy` and `Weights` build the same object, so
    before the verbs there was nothing to test — a weight written on a rate simply behaved as a
    factor, and a factor written on `transfer_to` as a weight, both in silence. The verb says which
    the number is, so each side now refuses the other and names the one that fits."""
    from zombi2 import genomes
    from zombi2.rates import Weights

    with pytest.raises(ValueError, match="carries Weights"):
        genomes.simulate_genomes_family(
            tree, loss=0.25 * Weights(habitat, {"cave": 2.0}), initial_families=5, seed=1)

    with pytest.raises(ValueError, match="transfer_to takes Weights"):
        genomes.simulate_genomes_family(
            tree, transfer=0.4, initial_families=5, seed=1,
            transfer_to=ScaledBy(habitat, {"cave": 2.0}))

    # each with the verb that fits: both run
    genomes.simulate_genomes_family(
        tree, loss=0.25 * ScaledBy(habitat, {"cave": 2.0}), initial_families=5, seed=1)
    genomes.simulate_genomes_family(
        tree, transfer=0.4, initial_families=5, seed=1,
        transfer_to=Weights(habitat, {"cave": 2.0}))


def test_all_three_verbs_name_their_first_argument_the_same():
    """SPEC §7 and every error message call the thing a driven parameter reads a **driver**, and the
    manual writes `ScaledBy(driver, mapping)` in the sentence that introduces it. Two of the three
    verbs called it `value`, so `SetBy(driver=…)` worked and `ScaledBy(driver=…)` raised — one
    keyword that lands on one verb and not on its siblings, which is the drift this project's
    one-concept-one-word rule exists to stop."""
    from zombi2.rates import Weights

    for verb in (ScaledBy, Weights, SetBy):
        made = verb(driver="h.tsv", mapping={"cave": 2.0})
        assert made.driver == "h.tsv", verb
        assert repr(made).startswith(f"{verb.__name__}("), repr(made)
