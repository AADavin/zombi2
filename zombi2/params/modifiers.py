"""Rate modifiers — the connections a parameter makes to something else (SPEC §5).

A rate is ``scope(base)`` with verbs chained onto it, and each verb records a **modifier**: a
reading of some other quantity, turned into a dimensionless factor that multiplies the base. An
**extent** takes the same modifiers (SPEC §6).

**Nobody writes these classes.** They are what the verbs on `Rate` build::

    birth = PerLineage(1.0).changing_at({0: 1.0, 3: 0.3})              # -> OnTime
    birth = PerLineage(1.0).scaled_by(TotalDiversity(cap=100))         # -> OnTotalDiversity
    loss  = PerCopy(0.25).varying_among('families', LogNormal(0.0, 0.5))   # -> Drawn
    loss  = PerCopy(0.25).scaled_by('habitat.tsv', {'aquatic': 3.0})   # -> Driven

The **deterministic** modifiers (`OnTime`, `OnTotalDiversity`) compute their factor as a pure
function of the context. The **carried** ones hold a value the engine draws per unit and hands
back: `Inherited` drifts parent→child (via ``initial``/``descend``) and `Drawn` takes an
independent draw for each unit (via ``draw``). The unit is an argument, not a class — a draw among
families and a draw among lineages are one model at two attachments, which is why `Random` writes
both.

A modifier knows how to produce its own factor, or its own draw, and how it is **written**; how
parameters compose is `Rate`'s business, because a verb cannot be added without changing that.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any, ClassVar, Mapping

from .distributions import Distribution, LogNormal, as_distribution
from .retired import check_no_retired_keywords

#: The kinds of value a modifier can read, as the second half of `Modifier.reads`.
#:
#: A modifier is a *reading* of a **value**, and a value is made in one of a few ways. The kind says
#: which, and it decides who is responsible for the number: ``measured`` values the engine already
#: knows (the clock, a count) and a modifier computes its factor from the context; ``drawn`` and
#: ``inherited`` values have to be produced once per unit, remembered, and handed back, which only
#: the engine can do because it owns the generator and knows when units are born.
MEASURED = "measured"      # computed from the run's state whenever it is read
DRAWN = "drawn"            # drawn at random when the unit is created, then fixed
INHERITED = "inherited"    # the parent's, perturbed at each split
DRIVEN = "driven"          # another level's value: recorded beforehand, or growing alongside

#: The kinds the **engine** must draw and carry per unit, rather than a modifier computing them.
#: `Rate.carried_modifiers` is the query that finds them.
CARRIED_KINDS = (DRAWN, INHERITED)

#: The units a value can be attached to, in the order they nest — the other half of `Modifier.reads`.
#: A value's unit decides what may read it: a parameter may read a value when the parameter's own
#: units include the value's (SPEC §5), so a trait on a lineage can drive gene loss, and a family's
#: tempo cannot drive speciation. A unit here that no engine carries is a cell nobody has built, not
#: a different kind of model.
#:
#: **Plural**, because a value varies *among* families rather than being counted *per* family: "per"
#: is the scope word and nothing else (SPEC §5). The change is also a small safety win — ``'family'``
#: and ``'families'`` are different strings, so a rate written in the old vocabulary fails loudly
#: instead of quietly meaning something new.
UNITS = ("run", "lineages", "chromosomes", "families", "copies", "sites")

#: The name of each verb, as a constant. A verb is not only called: it is **recorded** on what it
#: builds (``Driven.verb``, ``OnTime.verb``), because the same object serves more than one verb and
#: only the writer knows which was typed. That makes the name a contract read outside this module —
#: the transfer engine refuses a `Driven` written with `SCALED_BY` where a choice's weights belong —
#: and a literal spelled out at each end is a contract that renames in silence: when these strings
#: last changed, that refusal stopped firing and nothing said so. They live here rather than in
#: `zombi2.params.connection` only because that module imports this one; it re-exports them, so
#: ``verbs.SCALED_BY`` is the name to compare against.
SCALED_BY = "scaled_by"
SET_BY = "set_by"
WEIGHTED_BY = "weighted_by"
VARYING_AMONG = "varying_among"
CHANGING_AT = "changing_at"


def values_at_birth(mods: "tuple[Modifier, ...]", rng,
                    shared: "dict[int, float] | None" = None) -> tuple[float, ...]:
    """The value a newly created unit carries, one per modifier, in written order.

    An `INHERITED` value starts from its own beginning (`Inherited.initial`); a `DRAWN` one is drawn.
    The dispatch reads `Modifier.reads`, not the class, so a carried modifier an engine has never
    heard of is drawn like the ones it has. Drawing in written order is what keeps a run
    reproducible, and drawing from **every** modifier is the point: taking only the first was how a
    second one silently left the model.

    ``shared`` makes one value shared between the rates of a **single unit**. It is a cache keyed by
    modifier identity: pass the same dict while producing each of that unit's rates, and a modifier
    written on two of them is drawn once and both rates get the same number. That is how "a family
    that loses fast also duplicates fast" is said — one object, read twice — against "fast at losing
    only", which is two objects. Two modifiers that merely compare equal are still two values,
    because the question is whether you wrote one thing or two. Omit the cache and each draws for
    itself.

    Callers wanting the combined factor take ``math.prod`` of the result; a unit that never splits
    (a gene family) needs only that, while one that does keeps the values apart, because an
    inherited value has to perturb its parent's own number rather than a product."""
    out = []
    for m in mods:
        key = id(m)
        if shared is None or key not in shared:
            value = m.initial() if m.reads and m.reads[0] == INHERITED else m.draw(rng)
            if shared is None:
                out.append(value)
                continue
            shared[key] = value
        out.append(shared[key])
    return tuple(out)


def values_at_split(mods: "tuple[Modifier, ...]", parent_values: tuple[float, ...], rng,
                    shared: "dict[int, float] | None" = None) -> tuple[float, ...]:
    """A daughter's carried values: its parent's, perturbed (`INHERITED`), or a fresh independent
    draw that ignores the parent (`DRAWN`). That one line is the whole autocorrelated / uncorrelated
    split (SPEC §5). ``shared`` works as in `values_at_birth`."""
    out = []
    for i, m in enumerate(mods):
        key = id(m)
        if shared is None or key not in shared:
            value = (m.descend(parent_values[i], rng)
                     if m.reads and m.reads[0] == INHERITED else m.draw(rng))
            if shared is None:
                out.append(value)
                continue
            shared[key] = value
        out.append(shared[key])
    return tuple(out)


def check_one_memory(mods: "tuple[Modifier, ...]", *, label: str, unit: str) -> None:
    """SPEC §5's **one memory structure per axis**: a value on one unit is either drawn afresh each
    time (no memory) or inherited and perturbed (continuous memory), and those are two accounts of
    the same thing rather than a composition.

    So mixing the two kinds on one unit raises. Several of the **same** kind do not: two drawn
    factors multiply to one drawn factor, which is an ordinary composition and is what modifiers do.
    Every level calls this rather than writing its own count, so the rule cannot be strict in one
    place and lax in another — it used to be three different rules in three engines."""
    kinds = {m.reads[0] for m in mods if m.reads}
    if DRAWN in kinds and INHERITED in kinds:
        names = ", ".join(sorted(describe(m) for m in mods))
        raise ValueError(
            f"{label} carries both a drawn and an inherited value among {unit} ({names}), which are "
            f"the two answers to the same question — where that unit's factor comes from. An "
            f"inherited one starts from its parent's and is perturbed (autocorrelated); a drawn one "
            f"starts afresh with no memory of the parent (uncorrelated). Pick one — a law is either "
            f"a bare distribution or a Drift, never both. Several of the same kind are fine and "
            f"multiply.")


#: How each modifier class is **written**, for the messages that list what a level accepts. A
#: declaration is a promise about what you may write, and a class name is not something anyone
#: writes any more: the verbs are methods on the parameter, so the promise is spelled as one.
_WRITTEN_AS = {}      # filled in below, once the classes exist


def cell_name(entry) -> str:
    """What to call one entry of a level's ``IMPLEMENTED_MODIFIERS`` in a message — how that class is
    written, or how a ``(kind, unit)`` cell is written. Shared so an error and the CLI's help cannot
    describe the same declaration two different ways.

    Every entry is named by the **expression that writes it**, with ``...`` where the argument the
    user chooses goes: ``varying_among('families', ...)``, ``scaled_by(TotalDiversity(...))``. A
    declaration is a promise about what you may write, and this list is what an engine's refusal and
    ``zombi2 <command> -h`` both print, so a name nobody can type sends the reader to a syntax error.
    The earlier wording for a cell, ``drawn among families``, did exactly that.

    One verb writes both cells and the **law** is what differs: a bare distribution is drawn afresh
    for each unit, a `Drift` starts from the parent's value (SPEC §5). Only those two kinds reach
    here, because a cell is the grain for exactly the carried ones (`CARRIED_KINDS`) and everything
    else is declared by class."""
    if not isinstance(entry, tuple):
        return _WRITTEN_AS.get(entry, entry.__name__)
    kind, unit = entry
    return f"{VARYING_AMONG}({unit!r}, {'Drift(...)' if kind == INHERITED else '...'})"


def _driver_form(driver: object) -> str:
    """How a driver is written in a run's log — which has to be either the expression that
    reproduces it, or something that plainly is not one.

    A string driver is a filename or a level name and writes itself. A driver that knows its own
    written form (a `Clade`, which is built from literals) gives it. Anything else — a grown
    ``TraitsResult``, a genome's ``presence(...)`` — is an object from an earlier run and cannot be
    written at all, so it is recorded as ``<TraitsResult>``: a placeholder that fails loudly if
    pasted back, rather than a quoted ``'<TraitsResult>'`` that would be read as a *filename* and
    look like a run someone could reproduce."""
    if isinstance(driver, str):
        return repr(driver)
    written = getattr(driver, "written_form", None)
    if callable(written):
        return written()
    return f"<{type(driver).__name__}>"


def describe(m: "Modifier") -> str:
    """What to call one modifier **instance** in a message.

    A `Drawn` or an `Inherited` covers a whole row of the grid and is named by its cell —
    ``varying_among('families', ...)``, ``varying_among('lineages', Drift(...))`` — because "carries
    a Random" would be true and useless when the whole question is *among what*, and the law is what
    separates the two. Anything else is named by the **verb** that built it. Either way the name is
    the spelling that writes it, so a refusal and the list of what the level does take are in one
    vocabulary."""
    if isinstance(m, (Drawn, Inherited)):
        return cell_name(m.reads)
    verb = getattr(m, "verb", None)
    if verb is not None:
        return verb
    return _WRITTEN_AS.get(type(m), type(m).__name__)


def is_implemented(m: "Modifier", engines: tuple, engine: str) -> bool:
    """Whether ``engine`` may run modifier ``m``: it matches one entry of that level's
    ``IMPLEMENTED_MODIFIERS``, or it names that engine in its own `Modifier.implemented_for`. Every
    engine gate goes through here, so the escape hatch cannot be honoured in one level and forgotten
    in another.

    An entry is **a class** or **a cell**. A class is the right grain for `OnTime` against
    `OnTotalDiversity`: both read a measured value on the run, yet an engine can thread a schedule's
    breakpoints without threading the standing diversity, so the two are separately declarable. A
    cell — ``(DRAWN, "families")`` — is the right grain for `Drawn` and `Inherited`, which cover
    every unit, where what an engine supports is the *unit* it can carry a number for."""
    if matches_declared(m, engines):
        return True
    if engine not in getattr(m, "implemented_for", ()):
        return False
    # The hatch lets a modifier of your own vouch for itself, and it can — for a factor it *computes*
    # from the context, which is a promise only the modifier has to keep. It cannot vouch for a
    # **carried** value: that number has to be drawn when a unit is born, kept, and handed back, and
    # only the engine can do those, for the units it declares. Accepting one on a unit the level does
    # not carry would draw nothing and skip its factor, so the rate would run undriven in silence —
    # the exact failure this whole gate exists to prevent, so the hatch stops here.
    # A `SetBy` is refused here for the same reason wearing a different hat: replacing a base is a
    # capability an engine has or has not, and only three declare it. A subclass of `SetBy` vouching
    # for itself would be admitted at the four that cannot honour one.
    reads = getattr(m, "reads", None)
    return not (getattr(m, "replaces_base", False)
                or (reads is not None and reads[0] in CARRIED_KINDS))


def matches_declared(m: "Modifier", entries: tuple) -> bool:
    """Whether ``m`` is one of the entries a level declares — **without** the third-party escape
    hatch of `Modifier.implemented_for`.

    The sequences level needs this rather than `is_implemented`, and the reason is worth keeping: it
    is the one engine that reads its modifiers itself instead of evaluating them through
    `Rate.effective`, because its clock is drawn per lineage before any site
    evolves. A modifier of someone else's could therefore be accepted by the hatch and then never
    called, which is exactly the silence the whole declaration mechanism exists to prevent."""
    for entry in entries:
        if isinstance(entry, tuple):
            if m.reads == entry:
                return True
        elif getattr(m, "replaces_base", False) or getattr(entry, "replaces_base", False):
            # A `SetBy` is a `Driven`, so a plain isinstance would let it in wherever a driver is
            # allowed — and replacing a base is a capability an engine has or has not, which
            # Driven's declaration says nothing about. Four levels admitted it that way and then
            # could not honour it: three overwrote the base in a loop so the last one written won,
            # and the sequence level multiplied them together. Match it by exact type instead, so a
            # level has to name `SetBy` to accept one.
            if type(m) is entry:
                return True
        elif isinstance(m, entry):
            return True
    return False


def _no_longer_multiplied(left: object, right: object) -> Exception:
    """The refusal for ``a * b`` anywhere in the grammar, written once.

    ``*`` used to be how a rate was composed, so every retired expression in a paper, a notebook or
    an old ``--params`` file starts by multiplying. Left to CPython that fails with "unsupported
    operand type(s)", which names two classes nobody wrote and says nothing about what to do; this
    names the verbs instead."""
    from .rate import RateCompositionError

    return RateCompositionError(
        f"'*' no longer composes a rate — the verbs do, and each returns a new rate so they chain: "
        f"PerCopy(0.25).scaled_by(driver, mapping), .set_by(driver, mapping), "
        f".varying_among('families', LogNormal(0.0, 0.5)), .changing_at({{0: 1.0, 3: 0.3}}). "
        f"Got {left!r} * {right!r}.")


class Modifier:
    """Base for rate modifiers.

    A modifier reads the context keys it cares about (``time``, ``diversity``,
    ``branch``, ``family``, …) and returns a dimensionless, non-negative multiplier;
    it ignores the rest. Abstract — use a subclass, and write a verb rather than a subclass.
    """

    #: What this modifier reads, as ``(kind, unit)`` — the value's kind (one of `MEASURED`,
    #: `DRAWN`, `INHERITED`, `DRIVEN`) and the unit it lives on (``"run"``, ``"lineages"``,
    #: ``"families"``, …). It is the one thing an engine dispatches on, so a modifier a level has
    #: never heard of is threaded like the ones it has.
    #:
    #: The split it records is the useful one. A **measured** value is one the engine already has,
    #: so the modifier computes its own factor from the context and the engine does nothing. A
    #: **drawn** or **inherited** value has to be produced once per unit, remembered for that
    #: unit's life, and handed back at every evaluation — which only the engine can do, because it
    #: owns the generator and knows when a unit is born. Those are the kinds in `CARRIED_KINDS`,
    #: and `Rate.carried_modifiers` is how an engine asks for them without knowing which classes exist.
    reads: ClassVar[tuple[str, str] | None] = None

    #: Whether this modifier **replaces** the base rather than multiplying it. Read by
    #: `is_implemented` and `matches_declared` instead of an isinstance check, so a gate can ask what
    #: a modifier does without importing the class that does it — which is what let `SetBy` move to
    #: `connection` beside the verb that writes it.
    replaces_base: ClassVar[bool] = False

    #: The engines a **third-party** modifier declares itself implemented for. Each engine ships an
    #: ``IMPLEMENTED_MODIFIERS`` tuple and refuses anything outside it, because a modifier it never
    #: reads would return its default 1.0 and give a run that is quietly not the model you asked for
    #: (SPEC §5). That gate is right, but it was also a closed door: a `Modifier` subclass of your own
    #: composed into a `Rate` correctly and was then refused by every level, with no registry and no
    #: entry point — so extending the grammar meant forking the package. Naming an engine here is the
    #: opt-in::
    #:
    #:     class OnLogTime(Modifier):
    #:         implemented_for = ("species",)
    #:         def factor(self, *, time: float = 0.0, **_): return 1.0 / (1.0 + time)
    #:
    #: The engine names, and the context each one supplies to `factor`:
    #:
    #: =====================  =================================================================
    #: ``species``            ``time``, ``lineages``, ``diversity``
    #: ``genomes.family``     ``time``, ``lineages``, ``copies``, ``drivers``
    #: ``genomes.ordered``    ``time``, ``lineages``, ``copies``, ``chromosomes``, ``drivers``
    #: ``genomes.nucleotide`` ``time``, ``lineages``, ``copies``, ``chromosomes``, ``drivers``
    #: ``traits.continuous``  ``time``, ``lineages``, ``diversity``, ``drivers``
    #: ``traits.discrete``    ``time``, ``lineages``, ``drivers``
    #: ``joint``              ``time``, ``lineages``, ``diversity``, ``drivers``
    #: =====================  =================================================================
    #:
    #: The genome engines thread ``drivers`` only when some rate or extent in the run is driven, so a
    #: modifier that reads it must default its key. At the nucleotide resolution ``copies`` is always
    #: 0 — gene events there are counted per lineage, not per copy, so there is no copy count to
    #: pass. An **extent** is read in the same context as its level's rates, so one list serves both.
    #:
    #: Two things this cannot vouch for, whatever it names. A **carried** value — drawn or inherited
    #: — is produced by the engine when a unit is born and handed back, which only the engine can do,
    #: for the units it declares. And a **`SetBy`**, which replaces a base rather than scaling one,
    #: is a capability three levels have and four do not. Either is admitted by a level naming it and
    #: by nothing else.
    #:
    #: **``sequences`` is not on that list, deliberately.** Every engine above evaluates its rate
    #: through `Rate.effective`, which multiplies in whatever `factor` returns. The sequence level
    #: reads its two kinds of modifier itself — the clock is *drawn among lineages* before any site
    #: evolves, not evaluated at an event — so a modifier declaring itself implemented there would be
    #: accepted and then never called, which is the silence this whole mechanism exists to prevent.
    #: It refuses instead, and says why.
    #:
    #: **The hatch cannot vouch for a carried value.** It works for a modifier that computes its own
    #: factor from the context — ``reads`` unset, or `MEASURED` / `DRIVEN`. A modifier declaring
    #: `DRAWN` or `INHERITED` needs the *engine* to draw its number for each unit and hand it back,
    #: which an engine can only do for the units it declares, so such a modifier is admitted by a
    #: level naming its cell and by nothing else.
    #:
    #: Declaring an engine is a claim you are making: it calls `factor` with the context above and
    #: nothing more, so take ``**_`` and default every key you read. Built-in modifiers leave this
    #: empty; the engine lists them by type. The rate *text* grammar (a `--birth` flag, a ``--params``
    #: file) knows only the built-in names, so a modifier of your own is Python-only — as an object
    #: you constructed has to be. Worked examples: Chapter 2, "Writing your own".
    implemented_for: tuple[str, ...] = ()

    def factor(self, **context: Any) -> float:
        raise NotImplementedError

    def draw(self, rng) -> float:
        """One value for a newly created unit — what a modifier reading a `DRAWN` value provides.

        A modifier reading an `INHERITED` value implements `initial` and `descend` instead, because a
        daughter's number starts from its parent's rather than from nothing. Everything else needs
        neither, so the default says so rather than returning a plausible 1.0."""
        raise NotImplementedError(
            f"{type(self).__name__} does not draw a value per unit; it reads {self.reads!r}")

    def initial(self) -> float:
        """The value a **root** unit starts with, for a modifier reading an `INHERITED` value —
        where the walk down the tree begins. A `DRAWN` one has `draw` instead."""
        raise NotImplementedError(
            f"{type(self).__name__} does not inherit a value per unit; it reads {self.reads!r}")

    def descend(self, parent_value: float, rng) -> float:
        """A daughter's value from its parent's, for a modifier reading an `INHERITED` value. This is
        the whole autocorrelated / uncorrelated split: an inherited value starts here, a drawn one
        ignores its parent entirely."""
        raise NotImplementedError(
            f"{type(self).__name__} does not inherit a value per unit; it reads {self.reads!r}")

    def next_change(self, time: float) -> float:
        """The next time strictly after ``time`` at which this modifier's factor changes on
        its own — a skyline breakpoint. ``inf`` if it never changes with time (the default;
        most modifiers change only at events, not autonomously)."""
        return math.inf

    def written_call(self) -> str:
        """How this modifier is written as a **verb call** on a parameter — ``scaled_by(...)``,
        ``varying_among(...)``, ``changing_at(...)``. `Rate.__repr__` joins these onto the scope to
        render the whole expression, so this is the one place each connection says how it is typed.

        A modifier of someone else's cannot be written at all — the text grammar whitelists names
        and knows only the built-in ones — so the default is a placeholder that fails loudly if
        pasted back, rather than an expression that looks reproducible and is not.

        The placeholder is built from the **class name**, and never from ``repr(self)``: `__repr__`
        below calls this, so a subclass overriding neither — which is exactly the third-party
        modifier `implemented_for` invites — sent the two into each other, and every log line, every
        ``--params`` record and every error message that named the rate died of a RecursionError
        while the run itself carried on fine. It is the same shape `_driver_form` uses for a driver
        that cannot be written, for the same reason."""
        return f"<{type(self).__name__}>"

    def __repr__(self) -> str:
        """The verb call, because for most modifiers that is the only way to write one. `Drawn` and
        `Inherited` override this with their standalone ``Random(...)`` form, which is a thing you
        can name and share (SPEC §6); nothing else has a name of its own."""
        return self.written_call()

    def __mul__(self, other: object):
        raise _no_longer_multiplied(self, other)

    def __rmul__(self, other: object):
        # `0.25 * Drawn(...)` and friends. `*` composed a rate until the verbs replaced it, so it
        # raises here with the sentence that says what to write instead: a bare TypeError from
        # CPython would name two types and nothing a reader of a rate can act on.
        raise _no_longer_multiplied(other, self)


class OnTime(Modifier):
    """The rate changes in time — a skyline / episodic schedule. Written ``changing_at``::

        birth = PerLineage(0.5).changing_at({0: 1.0, 3: 0.3})    # 1.0 on [0, 3), then 0.3 on

    ``schedule`` maps each interval's start time to a relative factor, dimensionless: on a base of
    ``2.0`` the schedule scales it. Before the earliest breakpoint the earliest factor applies
    (define the schedule from time 0 to avoid surprise).

    The **same object** carries the other half of `Time`, ``set_by(Time(), ...)``, where the
    schedule holds the rates themselves rather than multiples of a base::

        birth = PerLineage(0.5).changing_at({0: 1.0, 3: 0.3})    # 30% of what it was
        birth = PerLineage().set_by(Time(), {0: 0.5, 3: 0.15})   # the rates themselves

    Those two are the same model, because a schedule on a base of 1.0 *is* the rate — which is why
    `Rate.set_by` builds this rather than a `SetBy` and there is nothing for an engine to learn.
    What differs is the sentence the reader typed, so `verb` records which, exactly as a `Driven`
    records its own: a run's log has to say back what was written, not an equivalent.
    """

    reads: ClassVar[tuple[str, str] | None] = (MEASURED, "run")

    #: `CHANGING_AT`, or `SET_BY` when the schedule holds the rates rather than factors.
    verb: str = CHANGING_AT

    def __init__(self, schedule: Mapping[float, float], *, verb: str | None = None) -> None:
        steps = tuple(sorted((float(t), float(f)) for t, f in schedule.items()))
        if not steps:
            raise ValueError(
                "a time schedule cannot be empty, e.g. .changing_at({0: 1.0, 3: 0.3})")
        for t, f in steps:
            if not math.isfinite(t):
                raise ValueError(f"schedule times must be finite, got {t!r}")
            if not math.isfinite(f) or f < 0:
                raise ValueError(f"schedule values must be finite and non-negative, got {f!r}")
        self._steps = steps
        if verb is not None:
            self.verb = verb

    # `factor` narrows the base signature: this modifier cannot answer without its key, and
    # giving it a default would make a level that forgot to thread it return a plausible
    # wrong number in silence. `IMPLEMENTED_MODIFIERS` is what guarantees the key arrives — a level
    # that does not thread it rejects the modifier outright rather than reaching here.
    def factor(self, *, time: float, **_: Any) -> float:  # type: ignore[override]
        f = self._steps[0][1]  # before the first breakpoint, the earliest factor applies
        for t, fac in self._steps:
            if t <= time:
                f = fac
            else:
                break
        return f

    def next_change(self, time: float) -> float:
        for t, _ in self._steps:  # steps are sorted; the first breakpoint strictly after `time`
            if t > time:
                return t
        return math.inf

    def written_call(self) -> str:
        # `repr(float)` rather than `:g`, which rounds to six significant figures: this text is what
        # a run's log records and what a reader pastes back into a flag, so a rate that prints
        # rounded is a log naming a model that is not the one that ran.
        inner = ", ".join(f"{float(t)!r}: {float(f)!r}" for t, f in self._steps)
        if self.verb == SET_BY:
            return f"{SET_BY}(Time(), {{{inner}}})"
        return f"{CHANGING_AT}({{{inner}}})"

    def __eq__(self, other: object) -> bool:
        """By the schedule alone. `verb` is outside this on purpose: the two spellings run the same
        model, and two rates that behave identically should compare equal."""
        return isinstance(other, OnTime) and other._steps == self._steps

    def __hash__(self) -> int:
        return hash((OnTime, self._steps))


@dataclass(frozen=True)
class TotalDiversity:
    """The lineages standing right now, as a **driver** a rate can be scaled by::

        birth = PerLineage(1.0).scaled_by(TotalDiversity(cap=100))

    The carrying capacity is written on the driver rather than in a mapping beside it, and that is
    a limit rather than a design: the factor an engine reads is the linear fall to a cap, and a
    general curve of diversity would have to be integrated rather than read at a point, exactly as
    a smooth function of time would (SPEC §10). So there is one shape, and it takes its one number
    here. `scaled_by` refuses a mapping alongside, rather than accepting one and ignoring it.
    """

    cap: float | None = None

    def __post_init__(self) -> None:
        if self.cap is None:
            raise ValueError(
                "TotalDiversity needs its carrying capacity — TotalDiversity(cap=100), the "
                "diversity at which the rate reaches zero. The linear fall to a cap is the one "
                "shape an engine reads, so there is no mapping to write instead.")
        if isinstance(self.cap, bool) or not isinstance(self.cap, (int, float)):
            raise TypeError(f"TotalDiversity cap must be a real number, got {self.cap!r}")
        if not math.isfinite(self.cap) or self.cap <= 0:
            raise ValueError(f"TotalDiversity cap must be finite and positive, got {self.cap!r}")

    def __repr__(self) -> str:
        assert self.cap is not None            # __post_init__ refuses a driver without its cap
        return f"TotalDiversity(cap={float(self.cap)!r})"


@dataclass(frozen=True, repr=False)            # repr=False: the written form, from `Modifier`
class OnTotalDiversity(Modifier):
    """The rate slows as standing diversity grows — diversity-dependence. Built by
    ``scaled_by(TotalDiversity(cap=100))``.

    The factor falls linearly from 1 toward 0 as diversity rises to ``cap`` (a carrying
    capacity), and stays 0 beyond it: a cap of 100 halves the rate at 50 lineages and stops it
    at 100.
    """

    reads: ClassVar[tuple[str, str] | None] = (MEASURED, "run")

    cap: float

    def __post_init__(self) -> None:
        if isinstance(self.cap, bool) or not isinstance(self.cap, (int, float)):
            raise TypeError(f"TotalDiversity cap must be a real number, got {self.cap!r}")
        if not math.isfinite(self.cap) or self.cap <= 0:
            raise ValueError(f"TotalDiversity cap must be finite and positive, got {self.cap!r}")

    # `factor` narrows the base signature: this modifier cannot answer without its key, and
    # giving it a default would make a level that forgot to thread it return a plausible
    # wrong number in silence. `IMPLEMENTED_MODIFIERS` is what guarantees the key arrives — a level
    # that does not thread it rejects the modifier outright rather than reaching here.
    def factor(self, *, diversity: float, **_: Any) -> float:  # type: ignore[override]
        return max(0.0, 1.0 - diversity / self.cap)

    def written_call(self) -> str:
        return f"{SCALED_BY}(TotalDiversity(cap={float(self.cap)!r}))"


@dataclass(frozen=True)
class Drift:
    """A **law** for `Random`: the parent's value times a mean-corrected draw at each split —
    continuous memory down a genealogy::

        birth = PerLineage(1.0).varying_among('lineages', Drift(LogNormal(0.0, 0.2)))  # ClaDS

    ``dist`` is the **per-split step**, not the value. Each law owns and documents its own
    argument — a bare distribution is the value itself, ``Drift``'s is the step — so nobody has to
    infer the role from the slot, which is what the retired ``spread=`` made impossible.

    ``bins`` discretises the drift onto a ladder: the rate takes one of ``bins`` values and a
    daughter moves to a **neighbouring rung**, or stays — the rate-category clock. A knob rather
    than a law of its own, because the model is unchanged: a daughter starts from its parent and is
    perturbed.
    """

    #: `Any` rather than `Distribution`: a law may be written as anything `as_distribution` accepts,
    #: and `Inherited` is what coerces it — checking here would only move the same check earlier and
    #: leave two places to keep in step.
    dist: Any = None
    bins: int | None = None

    def __repr__(self) -> str:
        extra = f", bins={self.bins}" if self.bins is not None else ""
        return f"Drift({self.dist!r}{extra})"


@dataclass(frozen=True)
class Inherited(Modifier):
    """What ``varying_among(unit, Drift(dist))`` builds — a value inherited and perturbed.

    Geometric Brownian motion on the rate: the step is **mean-corrected** so ``E[factor] = 1``.
    Without the correction the rate inflates down the tree (``E[rate] ≈ e^{σ²/2}``) — a real
    historical bug. Any distribution that states a mean works, and is corrected by dividing by it,
    exactly as `Drawn` does. The engine drives `initial` / `descend`, threading each unit's current
    factor back through `Rate.effective`'s ``carried_factor``.

    ``bins`` takes a `LogNormal` step and nothing else, because the ladder's spacing *is* that
    step's σ — with any other shape there is no rung spacing to read, and inventing one would be a
    model nobody wrote. ``None`` is the continuous form.
    """

    unit: str = ""
    # `Distribution` after `__post_init__` coerces it; a caller may hand in anything
    # `as_distribution` accepts, exactly as `Drawn` does
    dist: "Distribution" = None    # type: ignore[assignment]
    bins: int | None = None

    @property
    def reads(self) -> tuple[str, str]:      # type: ignore[override]
        return (INHERITED, self.unit)

    def __post_init__(self) -> None:
        if self.unit not in UNITS:
            raise ValueError(
                f"unknown unit {self.unit!r}; a value varies among one of {list(UNITS)}")
        if self.dist is None:
            raise ValueError(
                f"Drift needs the distribution of its per-split step — "
                f"varying_among({self.unit!r}, Drift(LogNormal(0.0, 0.2))).")
        object.__setattr__(self, "dist", as_distribution(self.dist))
        mean = self.dist.mean()       # raises for a distribution that cannot state one
        if not math.isfinite(mean) or mean <= 0:
            raise ValueError(
                f"an inherited step is mean-corrected by its distribution's mean, which must be "
                f"finite and positive; {self.dist!r} has a mean of {mean!r}")
        if self.bins is not None:
            if isinstance(self.bins, bool) or not isinstance(self.bins, int):
                raise TypeError(f"Drift bins must be a whole number, got {self.bins!r}")
            if self.bins < 2:
                raise ValueError(f"Drift bins must be at least 2, got {self.bins}")
            if not isinstance(self.dist, LogNormal):
                raise ValueError(
                    f"Drift(..., bins=) takes a LogNormal step and got {self.dist!r}. The ladder's "
                    f"rungs are spaced by that step's sigma, so another shape leaves no spacing to "
                    f"read — write Drift(LogNormal(0.0, sigma), bins=n), or drop bins for the "
                    f"continuous drift, which takes any distribution.")

    def written_call(self) -> str:
        return f"{VARYING_AMONG}({self.unit!r}, {Drift(self.dist, self.bins)!r})"

    def __repr__(self) -> str:
        """The standalone form — the `Random` a user can name and share between two rates. ``bins``
        is omitted when it is unset, because the parser refuses ``bins=None``: a written form that
        will not parse is not a written form."""
        return f"Random({self.unit!r}, {Drift(self.dist, self.bins)!r})"

    def __eq__(self, other: object) -> bool:
        return (isinstance(other, Inherited) and other.unit == self.unit
                and other.dist == self.dist and other.bins == self.bins)

    def __hash__(self) -> int:
        return hash((Inherited, self.unit, self.dist, self.bins))

    def factor(self, **_: Any) -> float:
        """A carried modifier has no factor to compute. Its number is drawn when the unit is born and
        kept by the engine, which hands it back through `Rate.effective`'s ``carried_factor`` — so asking
        for one here, without that stored value, could only ever return a plausible 1.0 for a rate
        that should have varied."""
        raise NotImplementedError(
            f"a Random is carried, not computed: the engine draws its value for each of the "
            f"{self.unit} and passes it to Rate.effective(carried_factor=...). It has no factor of "
            f"its own.")

    def _ladder(self) -> list[float]:
        """The ``bins`` rungs, geometric in the step's sigma and scaled so their mean is 1.

        The mean is taken over the rungs because a nearest-neighbour walk with reflecting ends is
        uniform at stationarity, so an even ladder gives ``E[factor] = 1`` — the same promise the
        continuous form keeps with its lognormal correction."""
        n = self.bins or 0
        assert isinstance(self.dist, LogNormal)   # __post_init__ refuses bins with anything else
        sigma = self.dist.sigma
        raw = [math.exp(sigma * (k - (n - 1) / 2)) for k in range(n)]
        mean = sum(raw) / n
        return [r / mean for r in raw]

    def initial(self) -> float:
        """The root's factor: 1.0 for the continuous form; the middle rung when binned."""
        if self.bins is None:
            return 1.0
        rungs = self._ladder()
        return rungs[len(rungs) // 2]

    def descend(self, parent_value: float, rng) -> float:
        """A daughter's factor: the parent's, times one mean-corrected step — or, when binned, its
        parent's rung or one either side, with the ends reflecting.

        The parent's rung is recovered by nearest match rather than by inverting the ladder, so no
        float drift down a deep tree can put a lineage between two rungs.

        A lognormal step keeps its closed form, ``exp(N(-σ²/2, σ))``, rather than going through
        ``sample()/mean()``. The two are the same distribution and consume the same one draw, but they
        differ in the last bit — and a one-ULP change to a step is enough to move a Gillespie
        comparison, so every autocorrelated run in the repository would have shifted for nothing.
        ``mu`` cancels in the correction, which is why only ``sigma`` appears."""
        if self.bins is None:
            if isinstance(self.dist, LogNormal):
                sigma = self.dist.sigma
                return parent_value * math.exp(rng.normal(-0.5 * sigma * sigma, sigma))
            return parent_value * (self.dist.sample(rng) / self.dist.mean())
        rungs = self._ladder()
        here = min(range(len(rungs)), key=lambda i: abs(math.log(rungs[i] / parent_value)))
        step = int(rng.integers(-1, 2))                      # -1, 0 or +1
        return rungs[min(max(here + step, 0), len(rungs) - 1)]


class Drawn(Modifier):
    """What ``varying_among(unit, dist)`` builds when the law is a **bare distribution** — a value
    drawn once for each unit and then fixed for that unit's life, i.i.d. heterogeneity with no
    memory. The unit is an argument, so one class covers every cell::

        loss         = PerCopy(0.25).varying_among('families', LogNormal(0.0, 0.5))
        substitution = PerSite(1.0).varying_among('lineages', LogNormal(0.0, 0.3))
        loss         = PerCopy(0.25).varying_among('families', Gamma(shape=4.0, scale=0.25))

    The law is the distribution the **value** is drawn from — any built-in `Distribution`:
    ``Fixed``, ``Exponential``, ``Gamma``, ``LogNormal``, ``Uniform``, ``Geometric``. A callable or a
    scipy frozen distribution is refused here, because neither states a mean to normalise by; an
    extent takes both, being a size rather than a multiplier.

    **Every draw is normalised to mean 1**, by dividing by the distribution's own mean. A drawn value
    is a *multiplier*, and one that does not average to 1 changes what the base means — a base of 0.25
    would stop being the average rate. So widening the spread separates units and leaves the average
    one where you put it, and a distribution's **location is normalised away**: what it contributes is
    its shape. ``Exponential(1.0)`` and ``Exponential(7.0)`` are therefore the same law. A
    distribution that cannot state its mean — a bare callable, a scipy frozen distribution — is
    refused rather than normalised by a guess. (A number that *is* the rate rather than a factor is
    `set_by`, where nothing is normalised.)

    **Writing one object or two decides what varies together.** One `Random` read by several rates is
    one draw for that unit, so a family that loses fast also duplicates fast; two separately built
    ones are independent even with identical arguments (SPEC §5).
    """

    def __init__(self, unit: str, dist: object = None) -> None:
        if unit not in UNITS:
            raise ValueError(f"unknown unit {unit!r}; a value varies among one of {list(UNITS)}")
        if dist is None:
            raise ValueError(
                f"a Random needs the law its value follows — "
                f"varying_among({unit!r}, LogNormal(0.0, 0.5)).")
        self.unit = unit
        self.dist = as_distribution(dist)
        mean = self.dist.mean()       # raises for a distribution that cannot state one
        if not math.isfinite(mean) or mean <= 0:
            raise ValueError(
                f"a drawn multiplier is normalised by its distribution's mean, which must be finite "
                f"and positive; {self.dist!r} has a mean of {mean!r}")

    @property
    def reads(self) -> tuple[str, str]:      # type: ignore[override]
        return (DRAWN, self.unit)

    def factor(self, **_: Any) -> float:
        """A carried modifier has no factor to compute — see `Inherited.factor`."""
        raise NotImplementedError(
            f"a Random is carried, not computed: the engine draws its value for each of the "
            f"{self.unit} and passes it to Rate.effective(carried_factor=...). It has no factor of "
            f"its own.")

    def draw(self, rng) -> float:
        """One independent multiplier for a unit, with mean 1.

        A degenerate distribution gives exactly 1.0 and consumes no randomness, so switching
        heterogeneity off leaves every other draw in the run where it was — which is what makes a
        zero-variation control comparable to the run it controls for."""
        if isinstance(self.dist, LogNormal) and self.dist.sigma == 0.0:
            return 1.0
        return self.dist.sample(rng) / self.dist.mean()

    def written_call(self) -> str:
        return f"{VARYING_AMONG}({self.unit!r}, {self.dist!r})"

    def __repr__(self) -> str:
        """The standalone form — the `Random` a user can name and share between two rates."""
        return f"Random({self.unit!r}, {self.dist!r})"

    def __eq__(self, other: object) -> bool:
        """Two draws are the same modifier when they are on the same unit with the same
        distribution."""
        return isinstance(other, Drawn) and other.unit == self.unit and other.dist == self.dist

    def __hash__(self) -> int:
        return hash((Drawn, self.unit, self.dist))


def Random(unit: str | None = None, law: object = None, **retired: Any) -> Modifier:
    """A value drawn for each unit of that kind — the one driver that is not measured anywhere.

    ``unit`` is plural (``'lineages'``, ``'families'``, ``'copies'``, ``'sites'``,
    ``'chromosomes'``). The **law** says what happens to the value afterwards, which is a separate
    question from what it starts as::

        Random('families', LogNormal(0.0, 0.5))          # drawn once, held for that family's life
        Random('lineages', Drift(LogNormal(0.0, 0.3)))   # the parent's, perturbed at each split
        Random('lineages', Drift(LogNormal(0.0, 0.3), bins=8))    # the rate-category clock

    A bare distribution is deliberate, not an oversight: it follows the convention the grammar
    already uses everywhere — a bare dict is a table, a bare function a curve — where the plain case
    is written plainly and anything else is named.

    Usually written through the verb, ``rate.varying_among('families', law)``, which builds one of
    these and attaches it. Building it **by name** is how two rates share one draw::

        family_speed = Random('families', LogNormal(0.0, 0.5))
        duplication  = PerCopy(0.20).varying_among(family_speed)
        loss         = PerCopy(0.10).varying_among(family_speed)   # exactly half, in every family

    because the engine caches a unit's draw by object identity. Two separately built ``Random``
    objects are two draws even with identical arguments: the question is whether you wrote one
    thing or two.

    ``**retired`` catches ``spread=`` and ``per=``, the two keywords this replaced, so Python
    answers them with the same sentence a flag does rather than with "unexpected keyword argument".
    The unit has a default for that reason alone: ``Random(per='family', …)`` writes it into a
    keyword, and a required positional would make Python complain about the missing argument before
    anything here could say what ``per=`` became.
    """
    check_no_retired_keywords(retired, where="Random")
    if unit is None:
        raise TypeError(
            f"a Random needs the plural unit its value varies among — Random('families', "
            f"LogNormal(0.0, 0.5)); one of {list(UNITS)}.")
    if isinstance(law, Drift):
        return Inherited(unit, law.dist, law.bins)
    return Drawn(unit, law)


WRITABLE = ("Random", "Drift", "TotalDiversity")

_WRITTEN_AS.update({
    OnTime: CHANGING_AT,
    OnTotalDiversity: f"{SCALED_BY}(TotalDiversity(...))",
    Drawn: VARYING_AMONG,
    Inherited: VARYING_AMONG,
})   # `connection` adds its own two — it imports this module, not the other way round


__all__ = ["Modifier", "OnTime", "OnTotalDiversity", "TotalDiversity", "Drawn", "Inherited",
           "Drift", "Random",
           "MEASURED", "DRAWN", "INHERITED", "DRIVEN", "CARRIED_KINDS", "UNITS",
           "SCALED_BY", "SET_BY", "WEIGHTED_BY", "VARYING_AMONG", "CHANGING_AT",
           "values_at_birth", "values_at_split", "check_one_memory",
           "cell_name", "describe", "is_implemented", "matches_declared", "WRITABLE"]
