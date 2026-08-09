"""The verbs — what reading something does to a parameter (SPEC §4).

A parameter reads a **driver** through a **mapping**, and a verb says what the resulting number
does. There are three, because the engine can only do three things with a number::

    loss        = PerCopy(0.25).scaled_by(habitat, {"cave": 4.0, "surface": 1.0})   # multiply
    loss        = PerCopy().set_by(habitat, {"cave": 1.0, "surface": 0.25})         # replace
    transfer_to = Recipients().weighted_by(competence, {"competent": 3.0})          # compare

`scaled_by` multiplies, so its number is a dimensionless factor and the parameter keeps its base.
`set_by` replaces, so its number carries the parameter's own units and there is no base to write.
`weighted_by` is for a **choice** — an argument that decides *who*, not how fast — which has no base
at all, because only the ratios between candidates are read.

Two drivers are built in, take no naming, and are by far the most written, so each has **one verb of
its own** and that verb is the only way to write it::

    substitution = PerSite(0.01).varying_among('lineages', LogNormal(0.0, 0.3))   # Random
    birth        = PerLineage(0.5).changing_at({0: 1.0, 3: 0.3})                  # Time

`scaled_by` refuses both and names the shortcut, so there is exactly one spelling for each. There
are two shortcuts and there will not be a third: a new driver uses `scaled_by` and gets no verb.

**The verbs are methods on the parameters** (`Rate`, `Extent`, `Choice`), not names you call,
because a driver, a law or a mapping is added by dropping in a class while a verb cannot be added
without changing how parameters compose. What lives here is the one implementation each method
calls, so the same verb on a rate and on an extent cannot drift apart, and the list of names the
written form whitelists.
"""

from __future__ import annotations

from .modifiers import (CHANGING_AT, SCALED_BY, SET_BY, VARYING_AMONG, WEIGHTED_BY, Drawn, Driven,
                        Drift, Inherited, Modifier, OnTime, OnTotalDiversity, Random, SetBy,
                        TotalDiversity, describe)
from .retired import check_no_retired_keywords
from .values import Measured, Time

#: The verb names, and so the attribute names `zombi2.params.parse` will follow in an expression.
#: The written form is a call on an attribute, so this is the parser's whitelist for the attribute
#: half — the other half being `_NAMES`, the things it may call by name.
#:
#: Built from the constants defined in `zombi2.params.modifiers`, beside the ``verb`` field that
#: records which verb was written, and re-exported here so a reader outside this package compares
#: against ``verbs.SCALED_BY`` rather than against a string it spelled out itself.
VERBS = (SCALED_BY, SET_BY, WEIGHTED_BY, VARYING_AMONG, CHANGING_AT)

__all__ = ["VERBS", "SCALED_BY", "SET_BY", "WEIGHTED_BY", "VARYING_AMONG", "CHANGING_AT",
           "written_with", "scaled_by", "set_by", "weighted_by", "varying_among", "changing_at"]


def written_with(m: object, verb: str) -> bool:
    """Whether ``m`` records ``verb`` as the verb that wrote it.

    The same object serves several verbs — a `Driven` is what `scaled_by` and `weighted_by` both
    build — so which one was typed is a fact only the object remembers, and it is what tells a
    mismatched verb from a right one. Asking through here rather than by reaching for ``.verb``
    keeps the reading in one place and works for a modifier that records none.
    """
    return getattr(m, "verb", None) == verb


def _refuse_a_factor(driver: object, verb: str) -> None:
    """Something that is already a dimensionless factor, handed to a verb as though it were a
    driver. A `Random` has its own verb; anything else here is what a verb *builds*, so in neither
    case is there a value for a mapping to turn into a number."""
    if isinstance(driver, (Drawn, Inherited, Drift)):
        raise TypeError(
            f"a Random is already a factor and has its own verb: write "
            f".varying_among('families', LogNormal(0.0, 0.5)) rather than {verb}(...). Verbs are "
            f"for drivers a mapping has to turn into a number — a trait, a clock, a count.")
    if isinstance(driver, Modifier):
        raise TypeError(
            f"{describe(driver)} is what a verb builds, not something a verb reads: it is already a "
            f"dimensionless factor, so {verb}(...) has nothing to map. Write the verb that produces "
            f"it on the parameter itself.")


def _refuse_time(driver: object, verb: str) -> None:
    """`Time` handed to a general verb. Its two readings are `changing_at` (factors) and
    ``set_by(Time(), ...)`` (the rates themselves), and those are the only two spellings."""
    if isinstance(driver, Time):
        raise ValueError(
            f"the run's clock has its own verb: write .changing_at({{0: 1.0, 3: 0.3}}) rather than "
            f"{verb}(Time(), ...). If you meant the schedule to hold the rates themselves rather "
            f"than multiples of a base, that is PerLineage().set_by(Time(), {{0: 0.5, 3: 0.15}}).")


def _refuse_a_whole_rule(driver: object) -> None:
    """`Distance` or `Clades` handed to `weighted_by` as though it were a driver.

    Both are complete ``transfer_to`` rules — a rule read off the tree, mapping and all — rather
    than a value a mapping turns into a weight, and each is written on its own. Passing one as a
    driver builds a `Driven` whose lookup key no engine ever fills, and which renders as
    ``weighted_by(<Distance>, ...)``: a run's log recording a rule nobody can paste back, and one
    the round-trip check cannot catch because the text never parses. Splitting `Distance` into a
    driver and a mapping of its own is designed and not built, so this refuses rather than accepting
    something that would quietly weight nothing."""
    from .choice import Clades, Distance

    if isinstance(driver, Distance):
        raise ValueError(
            "Distance() is a whole transfer_to rule, not a driver a mapping reads: write "
            "transfer_to=Distance(decay=1.0), which weights a recipient by exp(-decay × d / depth). "
            "Choosing the shape of that fall yourself — weighted_by(Distance(), lambda d: 1 / (1 + "
            "d ** 2)) — is not implemented; today the decay is the one knob.")
    if isinstance(driver, Clades):
        raise ValueError(
            "Clades({...}, Between({...})) is a whole transfer_to rule, not a driver: it already "
            "carries the kernel that weights the (donor, recipient) clade pair, and it is written "
            "on its own — transfer_to=Clades({'A': ['n12'], 'B': ['n40']}, "
            "Between({('A', 'B'): 10.0})). "
            "weighted_by reads a value another level evolved, and clade membership is a fact about "
            "the tree.")


def _schedule(mapping: object, verb: str) -> dict:
    """A time driver's mapping. A schedule — ``{0: 1.0, 3: 0.3}``, the value from each breakpoint on
    — is the only shape an engine reads, because a rate that is piecewise-constant in time can be
    stepped to exactly. A smooth curve of time is a real model and is refused rather than
    approximated: it makes the rate vary continuously, which needs the engine to integrate its
    hazard rather than sample it."""
    if isinstance(mapping, dict):
        return mapping
    raise ValueError(
        f"{verb} takes a schedule — {{0: 1.0, 3: 0.3}}, the value from each breakpoint on — because "
        f"a rate that changes in steps can be stepped to exactly. A smooth function of time is not "
        f"implemented: it makes the rate vary continuously between events, which needs the engine "
        f"to integrate the rate rather than read it at a point.")


def scaled_by(driver: object, mapping: object = None, *, step: float | None = None) -> Modifier:
    """Multiply the parameter's base by a factor read from ``driver``.

    The factor is dimensionless, and almost every parameter takes one::

        loss  = PerCopy(0.25).scaled_by(habitat, {"cave": 4.0, "surface": 1.0})   # a grown trait
        birth = PerLineage(1.0).scaled_by(TotalDiversity(cap=100))                # the standing LTT

    ``mapping`` turns the driver's value into that factor, and its shape follows the driver's
    **type**: a categorical driver takes a table (a dict), a numerical one a curve (a callable) or a
    ``Scalar`` log-link.

    ``step`` is the resolution a **continuous** driver is read at, in the tree's own time units. A
    categorical driver switches at moments the engine can step to exactly and ignores it.
    """
    _refuse_time(driver, "scaled_by")
    _refuse_a_factor(driver, "scaled_by")
    if isinstance(driver, TotalDiversity):
        if mapping is not None:
            raise ValueError(
                "TotalDiversity carries its own shape — the linear fall to a cap — so there is no "
                "mapping to write beside it: scaled_by(TotalDiversity(cap=100)). A general curve "
                "of standing diversity is not implemented (SPEC §10).")
        assert driver.cap is not None            # TotalDiversity refuses a driver without its cap
        return OnTotalDiversity(driver.cap)
    if isinstance(driver, Measured):
        raise ValueError(
            f"scaled_by({type(driver).__name__}(), ...) is not implemented — that driver exists in "
            f"the grammar but no engine supplies it yet.")
    if mapping is None:
        raise ValueError(
            "scaled_by(driver, mapping) needs a mapping: a dict for a categorical driver, a "
            "callable for a numerical one.")
    return Driven(driver, mapping, step, verb=SCALED_BY)


def set_by(driver: object, mapping: object = None, *, step: float | None = None) -> Modifier:
    """Replace the parameter's base with a number read from ``driver``, in the parameter's own
    units::

        loss  = PerCopy().set_by("habitat.tsv", {"aquatic": 1.0, "terrestrial": 0.25})
        birth = PerLineage().set_by(Time(), {0: 0.5, 3: 0.15})

    Written with no base in front, because the driver supplies the whole number; the scope still
    stands, because replacing *how fast* says nothing about *per what*.

    The clock is the one driver that can replace as well as scale, and it needs no new machinery to:
    a schedule on a base of 1.0 *is* the rate, so this builds the same `OnTime` `changing_at` does
    and records which verb wrote it (`OnTime.verb`), so a run's log says back what was typed.
    """
    if isinstance(driver, Time):
        return OnTime(_schedule(mapping, "set_by(Time(), ...)"), verb=SET_BY)
    _refuse_a_factor(driver, "set_by")
    if isinstance(driver, (TotalDiversity, Measured)):
        raise ValueError(
            f"set_by({type(driver).__name__}(), ...) is not implemented — no engine can be handed "
            f"a base from that driver. Scale a base you write yourself instead: "
            f"scaled_by({type(driver).__name__}(...)).")
    if mapping is None:
        raise ValueError(
            "set_by(driver, mapping) needs a mapping: a dict for a categorical driver, a callable "
            "for a numerical one. Its numbers are the rate itself, not factors.")
    return SetBy(driver, mapping, step, verb=SET_BY)


def weighted_by(driver: object, mapping: object = None, *, step: float | None = None) -> Driven:
    """Weight the candidates of a **choice** — an argument that decides *who*, not how fast.

    ``transfer_to``, the recipient of a horizontal transfer, is the only choice today. A choice has
    no base, because only the ratios between candidates are read, and a weight of zero means that
    candidate cannot be chosen::

        transfer_to = Recipients().weighted_by(competence, {"competent": 3.0, "normal": 1.0})

    A weight may read **both ends** — the donor's group and the recipient's — through a ``Between``
    kernel, which is the mapping for a driver that sits on a pair rather than on one lineage.
    """
    _refuse_time(driver, "weighted_by")
    _refuse_a_factor(driver, "weighted_by")
    _refuse_a_whole_rule(driver)
    if isinstance(driver, (TotalDiversity, Measured)):
        raise ValueError(
            f"weighted_by({type(driver).__name__}(), ...) is not implemented — a choice weights "
            f"each candidate by something that candidate has, and that driver is a property of the "
            f"run rather than of a lineage, so every candidate would weigh the same.")
    if mapping is None:
        raise ValueError(
            "weighted_by(driver, mapping) needs a mapping: a dict of per-candidate weights, a "
            "callable, or a Between kernel to weight the (donor, recipient) pair.")
    return Driven(driver, mapping, step, verb=WEIGHTED_BY)


def varying_among(among: object = None, law: object = None, **retired: object) -> Modifier:
    """Let the parameter vary at random among the units of one kind (SPEC §6)::

        loss = PerCopy(0.25).varying_among('families', LogNormal(0.0, 0.5))
        rate = PerLineage(1.0).varying_among('lineages', Drift(LogNormal(0.0, 0.2)))

    ``among`` is the plural unit name and ``law`` says what happens to the drawn value afterwards —
    a bare distribution for a value drawn and held, a `Drift` for one carried down the tree and
    perturbed at each split.

    It also takes a **named** `Random`, with no second argument, which is how two rates share one
    draw: the engine caches a unit's value by object identity, so one object read twice is one
    number and two objects are two.

    ``**retired`` catches ``per=`` and ``spread=``, the keywords this verb replaced, so Python
    answers them with the same sentence a flag does. Every parameter's ``varying_among`` passes its
    own through to here, so the answer cannot be good at one level and absent at another. The unit
    has a default for that reason alone: ``varying_among(per='families', …)`` writes the unit into a
    keyword, and a required positional would make Python complain about the missing argument before
    anything here could say what ``per=`` became.
    """
    check_no_retired_keywords(retired, where="varying_among")
    if among is None:
        raise TypeError(
            "varying_among takes the plural unit to vary among and the law it follows — "
            "varying_among('families', LogNormal(0.0, 0.5)) — or a Random built by name, "
            "varying_among(family_speed), which is how two rates share one draw.")
    if isinstance(among, (Drawn, Inherited)):
        if law is not None:
            raise TypeError(
                "a named Random already carries its law, so varying_among takes it alone: "
                "varying_among(family_speed). Sharing one object is what makes two rates share one "
                "draw — build a second one and there is nothing to share.")
        return among
    if isinstance(among, Modifier):
        raise TypeError(
            f"varying_among takes a Random — a unit and a law — and {describe(among)} is not one.")
    if not isinstance(among, str):
        raise TypeError(
            f"varying_among takes the plural unit to vary among and the law it follows — "
            f"varying_among('families', LogNormal(0.0, 0.5)) — or a Random built by name. "
            f"Got {among!r}.")
    return Random(among, law)


def changing_at(schedule: object) -> OnTime:
    """Let the parameter change in time — a skyline, the run's clock read as a schedule of factors::

        birth = PerLineage(0.5).changing_at({0: 1.0, 3: 0.3})   # 1.0, then 30% of it from time 3

    The numbers are multiples of the base. For the other reading — the schedule holding the rates
    themselves — write ``set_by(Time(), ...)``, which builds the same thing on a base of 1.0.
    """
    return OnTime(_schedule(schedule, "changing_at"))
