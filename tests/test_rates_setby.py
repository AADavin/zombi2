"""`set_by` — a driver that **replaces** a parameter's base rather than multiplying it.

The literature usually states a driven rate absolutely: *the loss rate is 1.0 in caves*, not *four
times a background nobody wrote down*. `scaled_by` can only say the second, so saying the first meant
inventing a background and dividing by it. `set_by` says it directly.

It builds a `Driven`, so every engine that resolves drivers resolves this one — the trajectory, the
mid-branch switches and the mapping checks are all the same machinery. One line in `Rate.effective`
asks a `SetBy` for the base and everything else for a factor.
"""

from __future__ import annotations

import pytest

from zombi2 import genomes, traits
from zombi2.params import Extent, PerCopy, PerLineage, Recipients
from zombi2.params import scope
from zombi2.params.connection import SetBy
from zombi2.params.parameter import as_rate
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
        r = _rate(PerCopy().set_by("h.tsv", {"cave": 1.0, "surface": 0.25}))
        assert r.effective(copies=1, drivers={"h.tsv": "cave"}) == pytest.approx(1.0)
        assert r.effective(copies=1, drivers={"h.tsv": "surface"}) == pytest.approx(0.25)

    def test_the_scope_still_applies(self):
        """`set_by` replaces the base, not the *per what?*. A per-copy rate set to 1.0 is 1.0 per
        copy, so four copies make four — otherwise a replaced base would quietly become a total."""
        r = _rate(PerCopy().set_by("h.tsv", {"cave": 1.0}))
        assert r.effective(copies=4, drivers={"h.tsv": "cave"}) == pytest.approx(4.0)

    def test_a_replaced_base_can_still_be_scaled(self):
        r = _rate(PerCopy().set_by("h.tsv", {"cave": 1.0}).scaled_by("s.tsv", {"big": 3.0}))
        got = r.effective(copies=2, drivers={"h.tsv": "cave", "s.tsv": "big"})
        assert got == pytest.approx(6.0)          # 1.0 set, tripled, times two copies

    def test_a_written_base_is_ignored_and_so_is_refused(self):
        """Silently discarding a number someone wrote is the failure mode this grammar exists to
        avoid, so a base in front of `set_by` raises rather than dropping the 0.25."""
        with pytest.raises(TypeError, match="silently discard"):
            PerCopy(0.25).set_by("h.tsv", {"cave": 1.0})


class TestInARealRun:

    def test_an_absolute_statement_matches_the_same_model_stated_as_a_multiple(self, tree, habitat):
        """The correctness check that matters: ``PerCopy().set_by(driver, {...})`` and
        ``PerCopy(1.0).scaled_by(driver, {...})`` are the same model written two ways, so they must
        be the same run, event for event."""
        def run(loss):
            g = genomes.simulate_genomes_family(tree, loss=loss, initial_families=40, seed=3)
            return [(e.time, e.kind, e.family) for e in g.events]

        assert run(PerCopy().set_by(habitat, {"cave": 1.0, "surface": 0.05})) == \
            run(PerCopy(1.0).scaled_by(habitat, {"cave": 1.0, "surface": 0.05}))

    def test_it_drives_a_continuous_trait_rate_too(self, tree, habitat):
        """`node_values` is non-empty whatever the rate is, so the assertion has to be about the
        number `set_by` supplied: the same model written as a multiple of 1.0 must give the same
        run, and a different set of numbers a different one."""
        driven = traits.simulate_continuous(
            tree, rate=PerLineage().set_by(habitat, {"cave": 2.0, "surface": 0.1}), seed=1)
        same = traits.simulate_continuous(
            tree, rate=PerLineage(1.0).scaled_by(habitat, {"cave": 2.0, "surface": 0.1}), seed=1)
        other = traits.simulate_continuous(
            tree, rate=PerLineage().set_by(habitat, {"cave": 0.1, "surface": 2.0}), seed=1)

        assert driven.node_values == pytest.approx(same.node_values)
        assert driven.node_values != pytest.approx(other.node_values)

    def test_a_level_that_cannot_replace_a_base_refuses_it(self, tree, habitat):
        """A `SetBy` is a `Driven`, so a gate listing Driven would let it in anywhere a driver
        goes — and four levels admitted it that way and could not honour it. A level has to name
        `SetBy` to accept one."""
        with pytest.raises(ValueError, match="does not support|does not read"):
            traits.simulate_discrete(
                tree, states=["a", "b"],
                switch=PerLineage().set_by(habitat, {"cave": 2.0, "surface": 0.1}), seed=1)


class TestWhatItRefuses:

    def test_two_bases_on_one_rate(self):
        """Each claims to *be* the number, and no order of application is more right than the other,
        so it raises rather than letting whichever was written last win."""
        with pytest.raises(TypeError, match="silently discard"):
            PerCopy().set_by("h.tsv", {"cave": 1.0}).set_by("h.tsv", {"cave": 2.0})

    def test_two_bases_are_refused_however_they_are_assembled(self):
        """The verb refuses them, and `as_rate` refuses them again for anything that reaches it
        another way — the guard lives at the choke point every level already calls, because a rule
        enforced by whoever remembers to call it is a rule three levels did not have."""
        from zombi2.params.parameter import Rate
        from zombi2.params.scope import PerCopy as _PerCopy

        smuggled = Rate(1.0, _PerCopy,
                        (SetBy("h", {"c": 1.0}), SetBy("h", {"c": 2.0})))
        with pytest.raises(ValueError, match="a base can only be replaced once"):
            as_rate(smuggled, default_scope=_PerCopy, label="loss")

    def test_one_base_and_any_number_of_factors_is_fine(self):
        r = _rate(PerCopy().set_by("h.tsv", {"cave": 1.0})
                  .scaled_by("a", {"x": 2.0}).scaled_by("b", {"y": 3.0}))
        r.check_one_base("loss")          # does not raise


def test_it_is_written_the_way_it_is_read():
    """The repr is what a run's log records and a reader pastes back into a flag, so it has to be an
    expression that reproduces the rate."""
    from zombi2.params.parse import parse_rate

    written = repr(PerCopy().set_by("h.tsv", {"cave": 1.0}))
    assert written == "PerCopy().set_by('h.tsv', Table({'cave': 1.0}))"
    assert parse_rate(written) == PerCopy().set_by("h.tsv", {"cave": 1.0})

    both = PerCopy().set_by("h.tsv", {"cave": 1.0}).scaled_by("s", {"x": 2.0})
    assert isinstance(parse_rate(repr(both)).modifiers[0], SetBy)


class TestTheHolesAnAdversarialReviewFound:
    """Every one of these shipped broken for a few hours. `SetBy` is a `Driven`, so it passed every
    gate that listed `Driven` — including four levels that could not honour a replaced base — and
    the "no base in front of it" guard only ever saw the operand immediately to its left.

    The verbs made the rule sayable once: a `set_by` is refused whenever anything is already to its
    left, which is one check on one object rather than a check per operand position."""

    @pytest.mark.parametrize("build", [
        pytest.param(lambda: PerCopy(0.25).set_by("h", {"c": 1.0}), id="a base, then set_by"),
        pytest.param(lambda: PerCopy(0.25).scaled_by("s", {"b": 1.0}).set_by("h", {"c": 1.0}),
                     id="a base and a factor, then set_by"),
        pytest.param(lambda: PerCopy().changing_at({0: 1.0}).set_by("h", {"c": 1.0}),
                     id="a schedule, then set_by"),
        pytest.param(lambda: PerCopy().set_by("h", {"c": 1.0}).set_by("h", {"c": 2.0}),
                     id="set_by twice"),
    ])
    def test_a_base_can_never_be_written_where_it_would_be_discarded(self, build):
        """Only the first of these was caught. The others each built a Rate whose written base was
        then silently overwritten by the driver's number."""
        with pytest.raises(TypeError):
            build()

    @pytest.mark.parametrize("level", ["discrete-trait", "nucleotide", "sequences"])
    def test_a_level_that_cannot_replace_a_base_refuses_rather_than_mangling_it(self, level, tree):
        """Three of these took two `set_by` and ran the last one written; the sequence level, which
        reads its modifiers itself, multiplied them together instead. Both are the silence the whole
        declaration mechanism exists to prevent."""
        if level == "discrete-trait":
            one = PerLineage().set_by("h.tsv", {"cave": 1.0, "surface": 0.5})
            with pytest.raises(ValueError, match="does not support"):
                traits.simulate_discrete(tree, states=["a", "b"], switch=one, seed=1)
        elif level == "nucleotide":
            one = PerLineage().set_by("h.tsv", {"cave": 1.0, "surface": 0.5})
            with pytest.raises(ValueError, match="does not support"):
                genomes.simulate_genomes_nucleotide(tree, loss=one, root_length=1000, seed=1)
        else:
            from zombi2.params import PerSite
            from zombi2.sequences import jc69, simulate_sequences
            one = PerSite().set_by("h.tsv", {"cave": 1.0, "surface": 0.5})
            g = genomes.simulate_genomes_family(tree, duplication=0.1, loss=0.1,
                                                initial_families=3, seed=5)
            with pytest.raises(ValueError, match="does not read"):
                simulate_sequences(g, model=jc69(), length=40, substitution=one, seed=1)

    def test_an_extent_has_no_base_to_replace(self):
        """An extent is already an absolute size drawn from a distribution. A replaced base there was
        admitted and applied as a multiplier, which is a different model wearing the same words."""
        from zombi2.params.parameter import as_extent
        with pytest.raises(ValueError, match="an extent cannot be set_by"):
            Extent(500).set_by("h", {"c": 5.0})
        with pytest.raises(ValueError, match="an extent cannot be set_by"):
            as_extent(PerCopy().set_by("h", {"c": 5.0}))

    def test_a_clade_driven_rate_is_written_so_it_can_be_pasted_back(self):
        """It recorded a quoted ``'<Clade>'`` — which parses, as a *filename*, so the log looked
        reproducible and was not. A driver that can write itself now does; one that cannot records an
        unquoted placeholder that fails loudly."""
        from zombi2.params import Clade
        from zombi2.params.parse import parse_rate

        written = repr(PerCopy(0.2).scaled_by(Clade({"fast": ["n1", "n2"]}), {"fast": 3.0}))
        assert "<Clade>" not in written
        assert parse_rate(written).modifiers[0].driver == Clade({"fast": ["n1", "n2"]})

    def test_a_driver_that_cannot_be_written_says_so_rather_than_looking_like_a_file(self, habitat):
        assert repr(PerCopy(0.25).scaled_by(habitat, {"cave": 2.0})) \
            == "PerCopy(0.25).scaled_by(<TraitsResult>, Table({'cave': 2.0}))"

    def test_a_replaced_base_is_written_so_it_can_be_pasted_back(self):
        """A run records its rates in the written form, and a replaced base was recorded with a
        number in front of it — the one spelling `set_by` refuses. So every log naming a replaced
        base named a rate that would not parse, and the record of the model was not one you could
        run again. It is written first now, on a scope with no number, which is what a replaced base
        means."""
        from zombi2.params.parse import parse_rate, written_form

        alone = PerCopy().set_by("h.tsv", {"cave": 1.0})
        assert written_form(alone) == repr(alone)
        assert parse_rate(written_form(alone)) == alone

        both = _rate(PerCopy().set_by("h.tsv", {"cave": 1.0}).scaled_by("g.tsv", {"x": 2.0}))
        text = written_form(both)
        assert text.startswith("PerCopy().set_by(")     # written first: nothing may precede it
        assert parse_rate(text).modifiers == both.modifiers

    def test_a_choice_has_no_base_to_replace_either(self, tree, habitat):
        """The same hole as the extent above, on the other kind of target. `transfer_to` weights the
        candidate recipients against each other, so there is no base for `set_by` to replace and the
        word means nothing there."""
        with pytest.raises(ValueError, match="written from Recipients"):
            genomes.simulate_genomes_family(
                tree, transfer=0.5, initial_families=6, seed=2,
                transfer_to=PerCopy().set_by(habitat, {"cave": 3.0, "surface": 1.0}))

        # a choice refuses the verb where it is written, too
        with pytest.raises(TypeError, match="no base to replace"):
            Recipients().set_by(habitat, {"cave": 3.0})

        # the same numbers, spelled as what they are, still run
        genomes.simulate_genomes_family(
            tree, transfer=0.5, initial_families=6, seed=2,
            transfer_to=Recipients().weighted_by(habitat, {"cave": 3.0, "surface": 1.0}))


def test_a_modifier_of_your_own_cannot_vouch_for_a_replaced_base():
    """`implemented_for` lets a third-party modifier declare which engines read it, and it may
    promise that for a factor it *computes*. Replacing a base is not that: it is a capability three
    levels have and four do not, so a `SetBy` subclass admitted through the hatch would be honoured
    nowhere. The same reason a carried value cannot go through it."""
    from zombi2.params.evaluate import is_implemented, matches_declared

    class Mine(SetBy):
        implemented_for = ("species", "genomes.family", "traits.discrete")

    m = Mine("h.tsv", {"cave": 1.0})
    for engine in Mine.implemented_for:
        assert not is_implemented(m, (), engine), f"the hatch vouched for a SetBy at {engine}"

    # a level that names SetBy still takes it — the hatch is what closed, not the gate
    assert matches_declared(SetBy("h.tsv", {"cave": 1.0}), (SetBy,))


def test_the_verb_has_to_match_what_it_is_attached_to(tree, habitat):
    """`scaled_by` and `weighted_by` build the same object, so before the verbs there was nothing to
    test — a weight written on a rate simply behaved as a factor, and a factor written on
    `transfer_to` as a weight, both in silence. The verb is recorded on the object, so each side now
    refuses the other and names the one that fits."""
    from zombi2.params import connection as verbs

    with pytest.raises(ValueError, match="carries weighted_by"):
        genomes.simulate_genomes_family(
            tree, loss=PerCopy(0.25)._and(verbs.weighted_by(habitat, {"cave": 2.0})),
            initial_families=5, seed=1)

    with pytest.raises(ValueError, match="not scaled_by"):
        genomes.simulate_genomes_family(
            tree, transfer=0.4, initial_families=5, seed=1,
            transfer_to=verbs.scaled_by(habitat, {"cave": 2.0}))

    # each with the verb that fits: both run
    genomes.simulate_genomes_family(
        tree, loss=PerCopy(0.25).scaled_by(habitat, {"cave": 2.0}), initial_families=5, seed=1)
    genomes.simulate_genomes_family(
        tree, transfer=0.4, initial_families=5, seed=1,
        transfer_to=Recipients().weighted_by(habitat, {"cave": 2.0}))


def test_the_verb_is_refused_where_it_is_written_too():
    """Before the engine ever sees it: a rate has no candidates to compare, and a choice has no base
    to scale, so each names the verb that fits."""
    with pytest.raises(TypeError, match="the verb is scaled_by"):
        PerCopy(0.25).weighted_by("h.tsv", {"cave": 2.0})
    with pytest.raises(TypeError, match="The verb is weighted_by"):
        Recipients().scaled_by("h.tsv", {"cave": 2.0})


def test_all_three_verbs_name_their_first_argument_the_same():
    """SPEC §7 and every error message call the thing a driven parameter reads a **driver**, and the
    manual writes ``.scaled_by(driver, mapping)`` in the sentence that introduces it. Two of the
    three verbs called it `value`, so one keyword landed on one verb and not on its siblings — the
    drift this project's one-concept-one-word rule exists to stop."""
    assert PerCopy(0.25).scaled_by(driver="h.tsv", mapping={"cave": 2.0}) \
        == PerCopy(0.25).scaled_by("h.tsv", {"cave": 2.0})
    assert PerCopy().set_by(driver="h.tsv", mapping={"cave": 2.0}) \
        == PerCopy().set_by("h.tsv", {"cave": 2.0})
    assert Recipients().weighted_by(driver="h.tsv", mapping={"cave": 2.0}) \
        == Recipients().weighted_by("h.tsv", {"cave": 2.0})
