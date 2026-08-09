"""Every rate the grammar can build must survive being written down and read back.

`tests/test_rate_parse.py` checks a handful of hand-written expressions. That is the right check on
the wrong sample: the written form is a **reproducibility record**, so a run's log has to reproduce
the model whatever the model was, and a list someone remembered to extend covers only the shapes
someone remembered.

So this file enumerates the space instead. Every scope is crossed with every verb call and with every
ordered *pair* of them; whatever the grammar refuses to build is skipped, and whatever it builds is
written, reparsed, and compared four ways:

1. **the written form is a fixed point** — rendering the reparsed rate gives the same text back;
2. **the objects are equal** — the cheap pre-filter;
3. **the rates agree numerically**, over a grid of times, counts and driver states;
4. **the carried values agree**, drawn from a seeded generator.

Numbers 3 and 4 are the ones that matter, and 2 alone would not do. `Rate` is a frozen dataclass, so
``==`` compares fields — and some fields are deliberately outside ``__eq__`` (a `Driven`'s verb, an
`OnTime`'s) — so a rendering can compare equal while behaving differently. Number 4 exists because a
carried modifier's factor never passes through `Rate.effective` at all: the engine draws it and hands
it back, so a mangled law is invisible to number 3.

The expressions are built through the **verbs**, not by assembling modifiers, so the sample is rates
a user could have typed rather than objects only this file knows how to make.
"""

from __future__ import annotations

import itertools

import numpy as np
import pytest

from zombi2.params import Clade, Extent, Recipients, Time
from zombi2.params import modifiers as mod
from zombi2.params import scope
from zombi2.params import values as values_mod
from zombi2.params.distributions import Exponential, Fixed, Gamma, Geometric, LogNormal, Uniform
from zombi2.params.mapping import Between, Scalar, Table
from zombi2.params.modifiers import Driven, Drift, values_at_birth
from zombi2.params.parse import parse_rate, written_form
from zombi2.params.rate import Rate, as_rate

#: (name, chain one verb onto a rate, the values this verb's driver can take, does it replace the base)
#:
#: A fresh object every time, never a shared instance: two rates reading one object share a draw by
#: identity (SPEC §6), and a fixture that handed out the same `Random` twice would be testing that.
VERB_CALLS: list[tuple[str, object, tuple, bool]] = [
    ("skyline", lambda r: r.changing_at({0: 1.0, 1.5: 0.3, 4: 2.25}), (), False),
    ("skyline-flat", lambda r: r.changing_at({0: 0.5}), (), False),
    ("diversity", lambda r: r.scaled_by(mod.TotalDiversity(cap=100)), (), False),
    ("drawn-lineage", lambda r: r.varying_among("lineages", LogNormal(0.0, 0.35)), (), False),
    ("drawn-family", lambda r: r.varying_among("families", LogNormal(0.0, 0.6)), (), False),
    ("drawn-gamma", lambda r: r.varying_among("families", Gamma(shape=4.0, scale=0.25)), (), False),
    ("drawn-lognormal", lambda r: r.varying_among("lineages", LogNormal(0.0, 0.4)), (), False),
    ("drawn-exponential", lambda r: r.varying_among("families", Exponential(1.0)), (), False),
    ("drawn-uniform", lambda r: r.varying_among("lineages", Uniform(0.5, 1.5)), (), False),
    ("drawn-geometric", lambda r: r.varying_among("families", Geometric(3.0)), (), False),
    ("drawn-fixed", lambda r: r.varying_among("lineages", Fixed(2.0)), (), False),
    ("inherited", lambda r: r.varying_among("lineages", Drift(LogNormal(0.0, 0.2))), (), False),
    ("inherited-binned",
     lambda r: r.varying_among("lineages", Drift(LogNormal(0.0, 0.2), bins=8)), (), False),
    ("scaled-table", lambda r: r.scaled_by("h.tsv", {"a": 2.0, "b": 0.5}), ("a", "b"), False),
    ("scaled-table-default",
     lambda r: r.scaled_by("h.tsv", Table({"a": 2.0}, default=0.75)), ("a", "z"), False),
    ("scaled-scalar", lambda r: r.scaled_by("x.tsv", Scalar(0.7)), (0.0, 1.25, -2.0), False),
    ("scaled-clade",
     lambda r: r.scaled_by(Clade({"A": ["n1", "n2"], "B": 3}), {"A": 2.0}), ("A", "B", "rest"),
     False),
    ("scaled-step",
     lambda r: r.scaled_by("h.tsv", {"a": 2.0, "b": 0.5}, step=0.05), ("a", "b"), False),
    ("set-by", lambda r: r.set_by("h.tsv", {"a": 1.0, "b": 0.25}), ("a", "b"), True),
    ("set-by-step",
     lambda r: r.set_by("h.tsv", {"a": 1.0, "b": 0.25}, step=0.05), ("a", "b"), True),
    ("set-by-time", lambda r: r.set_by(Time(), {0: 1.0, 2: 0.4}), (), True),
]

SCOPES = (scope.Global, scope.PerLineage, scope.PerCopy, scope.PerSite, scope.PerChromosome)

#: counts and clock readings to evaluate at. The times straddle every breakpoint above, because a
#: skyline that reparsed with a shifted breakpoint agrees everywhere else.
COUNTS = [
    {"copies": 1, "lineages": 1, "sites": 1, "chromosomes": 1, "diversity": 1},
    {"copies": 37, "lineages": 5, "sites": 900, "chromosomes": 3, "diversity": 42},
    {"copies": 0, "lineages": 2, "sites": 12, "chromosomes": 1, "diversity": 250},
]
TIMES = (0.0, 0.75, 1.5, 1.6, 2.0, 3.999, 4.0, 9.0)


def _build(scope_cls, calls):
    """``scope(base)`` with those verbs chained on, or ``None`` where the grammar refuses that
    combination.

    Refusals are the point rather than an inconvenience: one memory structure per axis, one base per
    rate, a `set_by` with nothing in front of it. Whatever survives is a rate a user could have
    written, so it is a rate the written form has to reproduce."""
    try:
        # a replaced base is written first, from the bare scope: everything to its left is a base it
        # would discard, so a rate that starts with `set_by` starts with no number
        rate = scope_cls() if calls[0][1] else scope_cls(0.37)
        for apply, _replaces in calls:
            rate = apply(rate)
        rate = as_rate(rate, default_scope=scope_cls)
        rate.check_one_base("this rate")
        mod.check_one_memory(tuple(rate.modifiers), label="this rate", unit="lineages")
        mod.check_one_memory(tuple(rate.modifiers), label="this rate", unit="families")
    except (ValueError, TypeError):
        return None
    return rate


def _cases():
    """Every scope crossed with every verb call and every ordered pair of them."""
    singles = [((name,), ((f, sets),), probes) for name, f, probes, sets in VERB_CALLS]
    pairs = [((a[0], b[0]), ((a[1], a[3]), (b[1], b[3])), a[2] + b[2])
             for a, b in itertools.permutations(VERB_CALLS, 2)]
    for scope_cls in SCOPES:
        for names, calls, probes in singles + pairs:
            rate = _build(scope_cls, calls)
            if rate is not None:
                yield pytest.param(rate, probes,
                                   id=f"{scope_cls.__name__}-{'+'.join(names)}")


def _drivers(rate: Rate, value) -> dict:
    """Every driven modifier on ``rate`` pointed at ``value``.

    Built from the rate's own modifiers rather than from a fixed dict, because a driver's key is not
    stable across a round trip — a `Clade` keys on `id()`, so the reparsed rate's key is a different
    number for the same clade. Keying per rate is what makes the two comparable at all."""
    return {m.key: value for m in rate.modifiers if isinstance(m, Driven)}


def _values(rate: Rate, probes) -> list[float]:
    """The rate evaluated across the grid — the numbers a run would actually use."""
    out = []
    for counts, t in itertools.product(COUNTS, TIMES):
        for probe in (probes or (None,)):
            try:
                out.append(rate.effective(**counts, time=t, drivers=_drivers(rate, probe)))
            except (KeyError, TypeError, ValueError) as e:
                out.append(f"raised {type(e).__name__}")
    return out


def _carried(rate: Rate, seed: int = 20260809) -> list:
    """What the engine would draw and keep for one unit, and how it descends.

    A carried modifier's factor never reaches `Rate.effective` — the engine draws it once per unit
    and hands it back — so `_values` above cannot see it at all. Same seed, same sequence: a law or a
    distribution that came back changed shows up here and nowhere else."""
    out: list = []
    for m, _unit in rate.carried_modifiers():
        rng = np.random.default_rng(seed)
        out.append(tuple(round(v, 12) for v in values_at_birth((m,), rng)))
        if m.reads[0] == mod.INHERITED:
            v = m.initial()
            walk = [v]
            for _ in range(6):
                v = m.descend(v, rng)
                walk.append(round(v, 12))
            out.append(tuple(walk))
    return out


@pytest.mark.parametrize("rate,probes", list(_cases()))
def test_a_rate_survives_being_written_down_and_read_back(rate, probes):
    text = written_form(rate)
    # Resolved the way a level resolves the text it is handed: every written rate now carries its
    # scope, including a replaced base — `PerCopy().set_by(...)` — so the default only ever fills in
    # a bare number, which is exactly what a level does.
    back = as_rate(parse_rate(text), default_scope=rate.scope)

    assert written_form(back) == text, "the written form is not a fixed point"
    assert back == rate, f"reparsed to a different object: {text}"
    assert _values(back, probes) == _values(rate, probes), f"same object, different numbers: {text}"
    assert _carried(back) == _carried(rate), f"a carried value came back changed: {text}"
    assert [back.next_change(t) for t in TIMES] == [rate.next_change(t) for t in TIMES], text


#: Spellings a user may type that render back in the canonical short form. There is one written form
#: per rate, so these are not second spellings of the *record* — they are the longer things a person
#: writes, which have to read as the same rate.
ALTERNATE_SPELLINGS = [
    # a named Random, which is how two rates share one draw (SPEC §6); on one rate it is the same
    # model as the unit-and-law form, and renders as it
    ("PerCopy(0.25).varying_among(Random('families', LogNormal(0.0, 0.5)))",
     "PerCopy(0.25).varying_among('families', LogNormal(0.0, 0.5))"),
    ("PerLineage(1.0).varying_among(Random('lineages', Drift(LogNormal(0.0, 0.2))))",
     "PerLineage(1.0).varying_among('lineages', Drift(LogNormal(0.0, 0.2)))"),
    # the Python qualifier a manual snippet carries
    ("scope.PerCopy(0.25).changing_at({0: 1.0, 3: 0.3})",
     "PerCopy(0.25).changing_at({0.0: 1.0, 3.0: 0.3})"),
]


@pytest.mark.parametrize("typed,canonical", ALTERNATE_SPELLINGS,
                         ids=[t.split("(")[0] + "-" + str(i)
                              for i, (t, _) in enumerate(ALTERNATE_SPELLINGS)])
def test_a_longer_spelling_reads_as_the_rate_it_names(typed, canonical):
    assert written_form(parse_rate(typed)) == canonical
    assert parse_rate(typed) == parse_rate(canonical)


def test_the_enumeration_is_actually_covering_the_grammar():
    """A generative test that silently stopped generating would pass forever.

    So assert the shape of the sample itself: every verb call reaches at least one case, and the case
    count stays in the hundreds rather than collapsing to the handful this file replaced."""
    cases = list(_cases())
    assert len(cases) > 1000, f"only {len(cases)} rates generated"

    # every verb call in the fixture reaches at least one built rate (an id is "Scope-a+b")
    covered = {part for case in cases
               for part in case.id.split("-", 1)[1].split("+")}
    assert {name for name, _, _, _ in VERB_CALLS} <= covered

    # every verb a rate can take is written by some case, so a verb added without a fixture entry
    # fails here rather than going quietly untested. `weighted_by` is the one a rate refuses — it
    # belongs to a choice — and it is covered by the choice round trips below.
    from zombi2.params import verbs
    rendered = " ".join(written_form(case.values[0]) for case in cases)
    on_a_rate = [v for v in verbs.VERBS if v != verbs.WEIGHTED_BY]
    assert not [v for v in on_a_rate if v not in rendered], rendered[:200]

    # and so is every name the grammar says is writable. `Random` is written through its verb, so it
    # appears in the alternate spellings rather than in a rendering — one canonical form per rate is
    # the point, not an omission.
    corpus = rendered + " " + " ".join(t for t, _ in ALTERNATE_SPELLINGS)
    missing = [name for name in mod.WRITABLE + values_mod.WRITABLE if name not in corpus]
    assert not missing, f"writable but never written: {missing}"


def test_a_written_form_that_cannot_be_reparsed_is_a_known_hole_not_a_surprise():
    """A user function has no written form, and that is the one thing this file cannot promise.

    Kept as a test so the limit is stated where the guarantee is, rather than discovered by someone
    whose run log turned out not to reproduce their model."""
    rate = as_rate(scope.PerLineage(0.1).scaled_by("gc.tsv", lambda x: 1 + 2 * x),
                   default_scope=scope.PerLineage)
    text = written_form(rate)
    assert "<lambda>" in text
    with pytest.raises(Exception):
        parse_rate(text)


# --- the other two parameter kinds ----------------------------------------
#
# The `set_by`-drops-`step` hole above was found by checking rates. An extent and a choice are the
# other two things you can attach a driver to, and neither was checked anywhere.

EXTENTS = [
    ("bare-number", 500),
    ("fixed", Fixed(3.0)),
    ("geometric", Geometric(250.0)),
    ("gamma", Gamma(shape=2.0, scale=250.0)),
    ("scaled", Extent(800).scaled_by("h.tsv", {"host": 3.0, "free": 1.0})),
    ("scaled-step", Extent(800).scaled_by("h.tsv", {"host": 3.0}, step=0.25)),
    ("skyline", Extent(800).changing_at({0: 1.0, 3: 0.5})),
    ("both", Extent(800).changing_at({0: 1.0, 3: 0.5}).scaled_by("h.tsv", {"host": 3.0})),
]


@pytest.mark.parametrize("name,spec", EXTENTS, ids=[n for n, _ in EXTENTS])
def test_an_extent_survives_being_written_down_and_read_back(name, spec):
    """An extent is ``base × modifiers`` with no scope (SPEC §6), so it renders through the same
    writer a rate does and has to come back the same size."""
    from zombi2.params.extent import as_extent

    text = written_form(spec)
    back = parse_rate(text)
    assert written_form(back) == text, "the written form is not a fixed point"

    a, b = as_extent(spec), as_extent(back)
    assert a == b, f"reparsed to a different extent: {text}"
    ctx = {"copies": 12, "lineages": 3, "chromosomes": 1, "sites": 90, "time": 1.0,
           "drivers": {m.key: "host" for m in a.modifiers if isinstance(m, Driven)}}
    # same seed, same sizes: an extent's number is drawn, so equality of objects is not enough
    for t in (0.0, 2.0, 4.0):
        ra, rb = np.random.default_rng(11), np.random.default_rng(11)
        assert [a.sample(ra, **{**ctx, "time": t}) for _ in range(5)] == \
               [b.sample(rb, **{**ctx, "time": t}) for _ in range(5)], text


def test_a_transfer_choice_is_written_as_the_form_that_takes_it_back():
    """A choice is not a rate, and its written form differs in ways that bit twice.

    A named rule is written bare, because `--transfer-to` takes `uniform` and not `'uniform'`. And a
    weighting is written from `Recipients()`, without the base a rate carries — a choice has none and
    the flag refuses one in front, so the rate writer's output was an expression the CLI rejects."""
    from zombi2.genomes._transfer import resolve_transfer_to
    from zombi2.params.choice import Clades, Distance
    from zombi2.params.parse import written_choice

    assert written_choice("uniform") == "uniform"
    assert written_choice("distance") == "distance"
    assert written_choice(Distance(decay=2.0)) == "Distance(decay=2.0)"
    assert written_choice(Clades({"A": ["n1"]}, Between({("A", "A"): 1.0}, default=0.0))) == \
        "Clades({'A': ['n1']}, Between({('A', 'A'): 1.0}, default=0.0))"

    w = Recipients().weighted_by("gc.tsv", {"a": 2.0})
    assert written_choice(w) == "Recipients().weighted_by('gc.tsv', Table({'a': 2.0}))"
    assert not written_choice(w).startswith("1.0"), "a choice has no base to write in front of it"
    # and what it writes is what the engine takes back
    assert resolve_transfer_to(parse_rate(written_choice(w))) == resolve_transfer_to(w)
    for rule in ("uniform", "distance"):
        resolve_transfer_to(written_choice(rule))


def test_every_choice_rule_reads_back_from_its_own_written_form():
    """One notation across Python, the CLI and a --params file — for choices too.

    `Distance` and `Clades` used to live at the genome level, where the parser could not see them, so
    a non-default `Distance(decay=…)` was Python-only and `Clades(...)` could not be typed at all.
    They are grammar objects — things a user writes — so they now live beside the rest of it."""
    from zombi2.params.choice import Clades, Distance
    from zombi2.params.parse import written_choice

    for spec in (Distance(decay=3.0),
                 Clades({"A": ["n1", "n2"], "B": 40},
                        Between({("A", "B"): 1.0, ("B", "A"): 1.0}, default=0.0)),
                 Recipients(),
                 Recipients().weighted_by("gc.tsv", Table({"a": 2.0}, default=0.5)),
                 Recipients().weighted_by("gc.tsv", {"a": 2.0}).weighted_by("h.tsv", {"b": 3.0})):
        text = written_choice(spec)
        assert parse_rate(text) == spec, text
        assert written_choice(parse_rate(text)) == text, "not a fixed point"


def test_a_clades_repr_is_the_call_that_reads_back_not_the_dataclass_one():
    """A dataclass gives `Clades(groups=…, between=…)`. That is valid Python but not what the parser
    or the docs use, and a run's log is a record you paste back."""
    from zombi2.params.choice import Clades

    c = Clades({"A": ["n1"]}, Between({("A", "A"): 1.0}, default=0.0))
    assert repr(c).startswith("Clades({'A': ['n1']}, Between(")
    assert "groups=" not in repr(c)
