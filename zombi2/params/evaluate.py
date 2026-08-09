"""What an **engine** calls, and the base every built modifier shares.

A user writes a parameter and chains verbs onto it (`parameter`, `connection`). What those verbs
build is a `Modifier`, and what an engine does with one is here: draw and carry its value per unit,
ask what it reads, ask whether this level supports it, and name it in a message.

The split is by audience rather than by subject. Nothing here is a thing anyone writes; everything
here is something an engine calls — which is why this module can be read without knowing the grammar,
and the grammar can be read without knowing this module.
"""


from __future__ import annotations

import math
from typing import Any, ClassVar


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

#: The units a *value* varies among — ``'run'`` is a measured driver's attachment, not one of these.
VARYING_UNITS = tuple(u for u in UNITS if u != "run")

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
_WRITTEN_AS: dict[type, str] = {}      # filled in below, once the classes exist


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

    A **carried** value covers a whole row of the grid and is named by its cell —
    ``varying_among('families', ...)``, ``varying_among('lineages', Drift(...))`` — because "carries
    a Random" would be true and useless when the whole question is *among what*, and the law is what
    separates the two. Anything else is named by the **verb** that built it. Either way the name is
    the spelling that writes it, so a refusal and the list of what the level does take are in one
    vocabulary.

    Carried-ness is read off `Modifier.reads` rather than off the two classes that have
    it, so this module needs no import from `law` — which is what lets `law` import
    this one."""
    reads = getattr(m, "reads", None)
    if reads is not None and reads[0] in CARRIED_KINDS:
        return cell_name(reads)
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
    from .parameter import RateCompositionError

    return RateCompositionError(
        f"'*' no longer composes a rate — the verbs do, and each returns a new rate so they chain: "
        f"PerCopy(0.25).scaled_by(driver, mapping), .set_by(driver, mapping), "
        f".varying_among('families', LogNormal(0.0, 0.5)), .changing_at({{0: 1.0, 3: 0.3}}). "
        f"Got {left!r} * {right!r}.")


class Modifier:
    """Base for rate modifiers.

    A modifier reads the context keys it cares about (``time``, ``lineages``, ``diversity``,
    ``copies``, ``chromosomes``, ``drivers``, …) and returns a dimensionless, non-negative
    multiplier; it ignores the rest. Abstract — use a subclass, and write a verb rather than a
    subclass.
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
    #: is a capability three of the engines above have and four do not. Either is admitted by a
    #: level naming it and by nothing else.
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
    #: you constructed has to be. Worked example: Appendix A, "Writing your own".
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
        can name and share (SPEC §5); nothing else has a name of its own."""
        return self.written_call()

    def __mul__(self, other: object):
        raise _no_longer_multiplied(self, other)

    def __rmul__(self, other: object):
        # `0.25 * Drawn(...)` and friends. `*` composed a rate until the verbs replaced it, so it
        # raises here with the sentence that says what to write instead: a bare TypeError from
        # CPython would name two types and nothing a reader of a rate can act on.
        raise _no_longer_multiplied(other, self)


