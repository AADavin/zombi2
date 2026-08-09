"""Every rate the grammar can build must survive being written down and read back.

`tests/test_rate_parse.py` checks a handful of hand-written expressions. That is the right check on
the wrong sample: the written form is a **reproducibility record**, so a run's log has to reproduce
the model whatever the model was, and a list someone remembered to extend covers only the shapes
someone remembered.

So this file enumerates the space instead. Every scope is crossed with every modifier and with every
ordered *pair* of modifiers; whatever the grammar refuses to build is skipped, and whatever it builds
is written, reparsed, and compared four ways:

1. **the written form is a fixed point** — rendering the reparsed rate gives the same text back;
2. **the objects are equal** — the cheap pre-filter;
3. **the rates agree numerically**, over a grid of times, counts and driver states;
4. **the carried values agree**, drawn from a seeded generator.

Numbers 3 and 4 are the ones that matter, and 2 alone would not do. `Rate` is a frozen dataclass, so
``==`` compares fields — and some fields are deliberately outside ``__eq__`` (a `Driven`'s verb, a
`Drawn`'s spread), so a rendering can compare equal while behaving differently. Number 4 exists
because a carried modifier's factor never passes through `Rate.effective` at all: the engine draws it
and hands it back, so a mangled ``spread`` is invisible to number 3.

Written against the grammar as it stands, deliberately. A test that guards a migration has to be
known-good *before* the migration, not written in the same change it is meant to check.
"""

from __future__ import annotations

import itertools

import numpy as np
import pytest

from zombi2.rates import Clade, ScaledBy, SetBy, Time
from zombi2.rates import modifiers as mod
from zombi2.rates import scope
from zombi2.rates import values as values_mod
from zombi2.rates.distributions import Exponential, Fixed, Gamma, Geometric, LogNormal, Uniform
from zombi2.rates.mapping import Scalar, Table
from zombi2.rates.modifiers import CARRIED_KINDS, Driven, values_at_birth
from zombi2.rates.parse import parse_rate, written_form
from zombi2.rates.rate import Rate, as_rate

#: (name, build a fresh one, the values this modifier's driver can take)
#:
#: Fresh objects each time, never a shared instance: two rates reading one object share a draw by
#: identity (SPEC §5), and a fixture that handed out the same `Drawn` twice would be testing that.
MODIFIERS: list[tuple[str, object, tuple]] = [
    ("skyline", lambda: mod.OnTime({0: 1.0, 1.5: 0.3, 4: 2.25}), ()),
    ("skyline-flat", lambda: mod.OnTime({0: 0.5}), ()),
    ("diversity", lambda: mod.OnTotalDiversity(cap=100), ()),
    ("drawn-lineage", lambda: mod.Drawn(per="lineage", dist=LogNormal(0.0, 0.35)), ()),
    ("drawn-family", lambda: mod.Drawn(per="family", dist=LogNormal(0.0, 0.6)), ()),
    ("drawn-gamma", lambda: mod.Drawn(per="family", dist=Gamma(shape=4.0, scale=0.25)), ()),
    ("drawn-lognormal", lambda: mod.Drawn(per="lineage", dist=LogNormal(0.0, 0.4)), ()),
    ("drawn-exponential", lambda: mod.Drawn(per="family", dist=Exponential(1.0)), ()),
    ("drawn-uniform", lambda: mod.Drawn(per="lineage", dist=Uniform(0.5, 1.5)), ()),
    ("drawn-geometric", lambda: mod.Drawn(per="family", dist=Geometric(3.0)), ()),
    ("drawn-fixed", lambda: mod.Drawn(per="lineage", dist=Fixed(2.0)), ()),
    ("inherited", lambda: mod.Inherited(per="lineage", dist=LogNormal(0.0, 0.2)), ()),
    ("inherited-binned", lambda: mod.Inherited(per="lineage", dist=LogNormal(0.0, 0.2), bins=8), ()),
    ("driven-table", lambda: ScaledBy("h.tsv", {"a": 2.0, "b": 0.5}), ("a", "b")),
    ("driven-table-default", lambda: ScaledBy("h.tsv", Table({"a": 2.0}, default=0.75)), ("a", "z")),
    ("driven-scalar", lambda: ScaledBy("x.tsv", Scalar(0.7)), (0.0, 1.25, -2.0)),
    ("driven-clade", lambda: ScaledBy(Clade({"A": ["n1", "n2"], "B": 3}), {"A": 2.0}),
     ("A", "B", "rest")),
    ("driven-time", lambda: ScaledBy(Time(), {0: 1.0, 2: 0.4}), ()),
    ("driven-step", lambda: ScaledBy("h.tsv", {"a": 2.0, "b": 0.5}, step=0.05), ("a", "b")),
    ("set-by", lambda: SetBy("h.tsv", {"a": 1.0, "b": 0.25}), ("a", "b")),
    ("set-by-step", lambda: SetBy("h.tsv", {"a": 1.0, "b": 0.25}, step=0.05), ("a", "b")),
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


def _build(scope_cls, factories):
    """``scope(base) × modifiers``, or ``None`` where the grammar refuses that combination.

    Refusals are the point rather than an inconvenience: one memory structure per axis, one base per
    rate, a `SetBy` with nothing in front of it. Whatever survives is a rate a user could have
    written, so it is a rate the written form has to reproduce."""
    mods = [f() for f in factories]
    try:
        if isinstance(mods[0], mod.SetBy):     # a replaced base is written first, with no number
            spec: object = mods[0]
            for m in mods[1:]:
                spec = spec * m
        else:
            spec = scope_cls(0.37)
            for m in mods:
                spec = spec * m
        rate = as_rate(spec, default_scope=scope_cls)
        rate.check_one_base("this rate")
        mod.check_one_memory(tuple(rate.modifiers), label="this rate", unit="lineage")
        mod.check_one_memory(tuple(rate.modifiers), label="this rate", unit="family")
    except (ValueError, TypeError):
        return None
    return rate


def _cases():
    """Every scope crossed with every modifier and every ordered pair of them."""
    singles = [((name,), (f,), probes) for name, f, probes in MODIFIERS]
    pairs = [((a[0], b[0]), (a[1], b[1]), a[2] + b[2])
             for a, b in itertools.permutations(MODIFIERS, 2)]
    for scope_cls in SCOPES:
        for names, factories, probes in singles + pairs:
            rate = _build(scope_cls, factories)
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
    and hands it back — so `_values` above cannot see it at all. Same seed, same sequence: a spread
    or a distribution that came back changed shows up here and nowhere else."""
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
    # Resolved the way a level resolves the text it is handed. A bare number comes back a number, a
    # lone `SetBy` comes back a modifier, and `SetBy(...) * OnTime(...)` comes back a Rate with no
    # scope at all — because the written form of a **replaced base carries no scope**, there being no
    # way to write one in front of `SetBy` today. Applying the same default the original was resolved
    # with is therefore the honest comparison, and it is exactly what a level does.
    back = as_rate(parse_rate(text), default_scope=type(rate.scope))

    assert written_form(back) == text, "the written form is not a fixed point"
    assert back == rate, f"reparsed to a different object: {text}"
    assert _values(back, probes) == _values(rate, probes), f"same object, different numbers: {text}"
    assert _carried(back) == _carried(rate), f"a carried value came back changed: {text}"
    assert [back.next_change(t) for t in TIMES] == [rate.next_change(t) for t in TIMES], text


def test_the_enumeration_is_actually_covering_the_grammar():
    """A generative test that silently stopped generating would pass forever.

    So assert the shape of the sample itself: every writable modifier reaches at least one case, and
    the case count stays in the hundreds rather than collapsing to the handful this file replaced."""
    cases = list(_cases())
    assert len(cases) > 1000, f"only {len(cases)} rates generated"

    # every modifier in the fixture reaches at least one built rate (an id is "Scope-a+b")
    covered = {part for case in cases
               for part in case.id.split("-", 1)[1].split("+")}
    assert {name for name, _, _ in MODIFIERS} <= covered

    # and every name the grammar says is writable appears in some rendered expression, so a modifier
    # added to `WRITABLE` without a fixture entry fails here rather than going quietly untested
    rendered = " ".join(written_form(case.values[0]) for case in cases)
    missing = [name for name in mod.WRITABLE + values_mod.WRITABLE if name not in rendered]
    assert not missing, f"writable but never rendered: {missing}"


def test_a_written_form_that_cannot_be_reparsed_is_a_known_hole_not_a_surprise():
    """A user function has no written form, and that is the one thing this file cannot promise.

    Kept as a test so the limit is stated where the guarantee is, rather than discovered by someone
    whose run log turned out not to reproduce their model."""
    rate = as_rate(0.1 * ScaledBy("gc.tsv", lambda x: 1 + 2 * x), default_scope=scope.PerLineage)
    text = written_form(rate)
    assert "<lambda>" in text
    with pytest.raises(Exception):
        parse_rate(text)


# --- the other two parameter kinds ----------------------------------------
#
# The `SetBy`-drops-`step` hole above was found by checking rates. An extent and a choice are the
# other two things you can attach a driver to, and neither was checked anywhere.

EXTENTS = [
    ("bare-number", 500),
    ("fixed", Fixed(3.0)),
    ("geometric", Geometric(250.0)),
    ("gamma", Gamma(shape=2.0, scale=250.0)),
    ("driven", 800 * ScaledBy("h.tsv", {"host": 3.0, "free": 1.0})),
    ("skyline", 800 * mod.OnTime({0: 1.0, 3: 0.5})),
    ("drawn", 800 * mod.Drawn(per="family", dist=LogNormal(0.0, 0.4))),
]


@pytest.mark.parametrize("name,spec", EXTENTS, ids=[n for n, _ in EXTENTS])
def test_an_extent_survives_being_written_down_and_read_back(name, spec):
    """An extent is ``base × modifiers`` with no scope (SPEC §6), so it renders through the same
    writer a rate does and has to come back the same size."""
    from zombi2.rates.extent import as_extent

    text = written_form(spec)
    back = parse_rate(text)
    assert written_form(back) == text, "the written form is not a fixed point"

    a, b = as_extent(spec), as_extent(back)
    assert a == b, f"reparsed to a different extent: {text}"
    if any(m.reads and m.reads[0] in CARRIED_KINDS for m in a.modifiers):
        return   # a carried factor never reaches Extent.sample: the engine weights the segment by
                 # what it covers (SPEC §6), so there is no number here to compare
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
    `Weights` is written on its own, without the ``1.0 *`` a rate carries — a choice has no base and
    the flag refuses one in front, so the rate writer's output was an expression the CLI rejects."""
    from zombi2.genomes._transfer import Clades, Distance, resolve_transfer_to
    from zombi2.rates import Weights
    from zombi2.rates.mapping import Between
    from zombi2.rates.parse import written_choice

    assert written_choice("uniform") == "uniform"
    assert written_choice("distance") == "distance"
    assert written_choice(Distance(decay=2.0)) == "Distance(decay=2.0)"
    assert written_choice(Clades({"A": ["n1"]}, Between({("A", "A"): 1.0}, default=0.0))) == \
        "Clades({'A': ['n1']}, Between({('A', 'A'): 1.0}, default=0.0))"

    w = Weights("gc.tsv", {"a": 2.0})
    assert written_choice(w) == "Weights('gc.tsv', Table({'a': 2.0}))"
    assert not written_choice(w).startswith("1.0 *"), "a choice has no base to write in front of it"
    # and what it writes is what the engine takes back
    assert resolve_transfer_to(parse_rate(written_choice(w))) == w
    for rule in ("uniform", "distance"):
        resolve_transfer_to(written_choice(rule))


def test_the_cli_cannot_yet_write_two_choice_rules_that_python_can():
    """`Distance(decay=)` and `Clades(...)` are not in the parser's whitelist, so they are Python-only.

    Recorded as a test rather than left to be discovered: the project's rule is one notation across
    Python, the CLI and a --params file, and these two are the exception. The day they are added,
    this test fails and says so."""
    for text in ("Distance(decay=3.0)", "Clades({'A': ['n1']}, Between({('A', 'A'): 1.0}))"):
        with pytest.raises(ValueError, match="unknown name"):
            parse_rate(text)
