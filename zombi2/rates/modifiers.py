"""Rate modifiers — the context multipliers of a rate (SPEC §5).

Every rate in ZOMBI2 is ``scope(base) × modifiers``. A *modifier* multiplies a rate
by a dimensionless factor that depends on context — the current time, the standing
diversity, the branch, the family, a driver's value. An **extent** takes the same
modifiers (SPEC §6). Modifiers **multiply** (that is
the whole difference from scope wrappers, which *wrap*), and the word *"per"* is
reserved for scope, so a modifier never starts with "per".

You reach them through ``mod``::

    birth = 1.0 * mod.OnTime({0: 1.0, 3: 0.3})   # a skyline: 1.0, then 0.3 from time 3 on
    birth = 1.0 * mod.OnTotalDiversity(cap=100)       # slows to 0 as diversity approaches 100

The **deterministic** modifiers (``OnTime``, ``OnTotalDiversity``) compute their factor as a pure
function of the context. The **carried** ones hold a value the engine draws per unit and hands back:
``Inherited(per=…)`` drifts parent→child (via ``initial``/``descend``) and ``Drawn(per=…)`` takes an
independent draw for each unit (via ``draw``). The unit is an argument, not a class — a draw per
family and a draw per lineage are one model at two attachments.

Composition (``*``) belongs to the Rate module; a modifier here knows only how to produce its own
factor, or its own draw.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any, ClassVar, Mapping

from .distributions import LogNormal, as_distribution

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
UNITS = ("run", "lineage", "chromosome", "family", "copy", "site")



def values_at_birth(mods: "tuple[Modifier, ...]", rng,
                    shared: "dict[int, float] | None" = None) -> tuple[float, ...]:
    """The value a newly created unit carries, one per modifier, in written order.

    An `INHERITED` value starts from its own beginning (`Inherited.initial`); a `DRAWN` one is drawn.
    The dispatch reads `Modifier.reads`, not the class, so a per-unit modifier an engine has never
    heard of is handled like the ones it has. Drawing in written order is what keeps a run
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
            f"{label} carries both a drawn and an inherited value per {unit} ({names}), which are "
            f"the two answers to the same question — where that {unit}'s factor comes from. An "
            f"inherited one starts from its parent's and is perturbed (autocorrelated); a drawn one "
            f"starts afresh with no memory of the parent (uncorrelated). Pick one. Several of the "
            f"same kind are fine and multiply.")


def cell_name(entry) -> str:
    """What to call one entry of a level's ``IMPLEMENTED_MODIFIERS`` in a message — a class's name,
    or the named cell for a ``(kind, unit)`` pair. Shared so an error and the CLI's help cannot
    describe the same declaration two different ways."""
    if not isinstance(entry, tuple):
        return entry.__name__
    return f"{entry[0]} per {entry[1]}"


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

    Its class name, except for the two classes that cover a whole row of the grid: a `Drawn` or an
    `Inherited` is named by its cell — ``drawn per family``, ``inherited per lineage`` — because
    saying "carries Drawn" would be true and useless when the whole question is *per what*."""
    if isinstance(m, (Drawn, Inherited)):
        return cell_name(m.reads)
    return type(m).__name__


def is_implemented(m: "Modifier", engines: tuple, engine: str) -> bool:
    """Whether ``engine`` may run modifier ``m``: it matches one entry of that level's
    ``IMPLEMENTED_MODIFIERS``, or it names that engine in its own `Modifier.implemented_for`. Every
    engine gate goes through here, so the escape hatch cannot be honoured in one level and forgotten
    in another.

    An entry is **a class** or **a cell**. A class is the right grain for `OnTime` against
    `OnTotalDiversity`: both read a measured value on the run, yet an engine can thread a schedule's
    breakpoints without threading the standing diversity, so the two are separately declarable. A
    cell — ``(DRAWN, "family")`` — is the right grain for `Drawn` and `Inherited`, which are one class
    covering every unit, where what an engine supports is the *unit* it can carry a number for."""
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
    return not (isinstance(m, SetBy) or (reads is not None and reads[0] in CARRIED_KINDS))


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
        elif isinstance(m, SetBy) or isinstance(entry, type) and issubclass(entry, SetBy):
            # A `SetBy` is a `DrivenBy`, so a plain isinstance would let it in wherever a driver is
            # allowed — and replacing a base is a capability an engine has or has not, which
            # DrivenBy's declaration says nothing about. Four levels admitted it that way and then
            # could not honour it: three overwrote the base in a loop so the last one written won,
            # and the sequence level multiplied them together. Match it by exact type instead, so a
            # level has to name `SetBy` to accept one.
            if type(m) is entry:
                return True
        elif isinstance(m, entry):
            return True
    return False


class Modifier:
    """Base for rate modifiers.

    A modifier reads the context keys it cares about (``time``, ``diversity``,
    ``branch``, ``family``, …) and returns a dimensionless, non-negative multiplier;
    it ignores the rest. Abstract — use a subclass.
    """

    #: What this modifier reads, as ``(kind, unit)`` — the value's kind (one of `MEASURED`,
    #: `DRAWN`, `INHERITED`, `DRIVEN`) and the unit it lives on (``"run"``, ``"lineage"``,
    #: ``"family"``, …). It is SPEC §5's preposition table, written where the code can read it:
    #: ``On`` is a measured value, ``By`` a drawn one, ``From`` an inherited one, and ``DrivenBy``
    #: another level's.
    #:
    #: The split it records is the useful one. A **measured** value is one the engine already has,
    #: so the modifier computes its own factor from the context and the engine does nothing. A
    #: **drawn** or **inherited** value has to be produced once per unit, remembered for that
    #: unit's life, and handed back at every evaluation — which only the engine can do, because it
    #: owns the generator and knows when a unit is born. Those are the kinds in `CARRIED_KINDS`,
    #: and `Rate.carried_modifiers` is how an engine asks for them without knowing which classes exist.
    reads: ClassVar[tuple[str, str] | None] = None

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
    #: reads its two kinds of modifier itself — the clock is *drawn per lineage* before any site
    #: evolves, not evaluated at an event — so a modifier declaring itself implemented there would be
    #: accepted and then never called, which is the silence this whole mechanism exists to prevent.
    #: It refuses instead, and says why.
    #:
    #: **The hatch cannot vouch for a carried value.** It works for a modifier that computes its own
    #: factor from the context — ``reads`` unset, or `MEASURED` / `DRIVEN`. A modifier declaring
    #: `DRAWN` or `INHERITED` needs the *engine* to draw its number per unit and hand it back, which
    #: an engine can only do for the units it declares, so such a modifier is admitted by a level
    #: naming its cell and by nothing else.
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

    def __rmul__(self, other: object):
        # `number * mod`, `scope * mod`, `mod * mod`, `Rate * mod` all build a Rate (see zombi2.rate)
        from .rate import Rate
        from .scope import Scope

        if isinstance(other, bool):
            return NotImplemented
        if isinstance(other, (int, float)):
            return Rate(float(other), None, (self,))
        if isinstance(other, Scope):
            return Rate(other.base, other, (self,))
        if isinstance(other, Modifier):
            return Rate(1.0, None, (other, self))
        if isinstance(other, Rate):
            return Rate(other.base, other.scope, other.modifiers + (self,))
        return NotImplemented

    def __mul__(self, other: object):
        from .rate import Rate

        if isinstance(other, Modifier):
            return Rate(1.0, None, (self, other))
        return self.__rmul__(other)


class OnTime(Modifier):
    """The rate changes in time — a skyline / episodic schedule.

    ``schedule`` maps each interval's start time to a relative factor::

        OnTime({0: 1.0, 3: 0.3})   # factor 1.0 on [0, 3), then 0.3 from time 3 on

    Factors are relative (dimensionless): on a base of ``2.0`` the schedule scales it.
    Before the earliest breakpoint the earliest factor applies (define the schedule
    from time 0 to avoid surprise).
    """

    reads: ClassVar[tuple[str, str] | None] = (MEASURED, "run")

    def __init__(self, schedule: Mapping[float, float]) -> None:
        steps = tuple(sorted((float(t), float(f)) for t, f in schedule.items()))
        if not steps:
            raise ValueError("OnTime needs a non-empty schedule, e.g. OnTime({0: 1.0, 3: 0.3})")
        for t, f in steps:
            if not math.isfinite(t):
                raise ValueError(f"OnTime schedule times must be finite, got {t!r}")
            if not math.isfinite(f) or f < 0:
                raise ValueError(f"OnTime factors must be finite and non-negative, got {f!r}")
        self._steps = steps

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

    def __repr__(self) -> str:
        # `repr(float)` rather than `:g`, which rounds to six significant figures: this text is what
        # a run's log records and what a reader pastes back into a flag, so a rate that prints
        # rounded is a log naming a model that is not the one that ran.
        inner = ", ".join(f"{float(t)!r}: {float(f)!r}" for t, f in self._steps)
        return f"OnTime({{{inner}}})"

    def __eq__(self, other: object) -> bool:
        return isinstance(other, OnTime) and other._steps == self._steps

    def __hash__(self) -> int:
        return hash((OnTime, self._steps))


@dataclass(frozen=True)
class OnTotalDiversity(Modifier):
    """The rate slows as standing diversity grows — diversity-dependence.

    The factor falls linearly from 1 toward 0 as diversity rises to ``cap`` (a carrying
    capacity), and stays 0 beyond it: ``OnTotalDiversity(cap=100)`` halves the rate at 50
    lineages and stops it at 100.
    """

    reads: ClassVar[tuple[str, str] | None] = (MEASURED, "run")

    cap: float

    def __post_init__(self) -> None:
        if isinstance(self.cap, bool) or not isinstance(self.cap, (int, float)):
            raise TypeError(f"OnTotalDiversity cap must be a real number, got {self.cap!r}")
        if not math.isfinite(self.cap) or self.cap <= 0:
            raise ValueError(f"OnTotalDiversity cap must be finite and positive, got {self.cap!r}")

    # `factor` narrows the base signature: this modifier cannot answer without its key, and
    # giving it a default would make a level that forgot to thread it return a plausible
    # wrong number in silence. `IMPLEMENTED_MODIFIERS` is what guarantees the key arrives — a level
    # that does not thread it rejects the modifier outright rather than reaching here.
    def factor(self, *, diversity: float, **_: Any) -> float:  # type: ignore[override]
        return max(0.0, 1.0 - diversity / self.cap)


@dataclass(frozen=True)
class Inherited(Modifier):
    """A value **inherited from a parent and perturbed** at each split — continuous memory down a
    genealogy. ``per`` names the unit it is carried on, so the class is one model and the unit is
    data::

        birth = 1.0 * Inherited(per="lineage", spread=0.2)     # rate drift down a clade (ClaDS)

    Geometric Brownian motion on the rate: the per-split factor is lognormal and **mean-corrected**
    so ``E[factor] = 1``. Without the correction the rate inflates down the tree
    (``E[rate] ≈ e^{σ²/2}``) — a real historical bug. The engine drives `initial` / `descend`,
    threading each unit's current factor back through `Rate.effective`'s ``carried_factor``.

    ``bins`` discretises the drift: the rate takes one of ``bins`` values on a geometric ladder and a
    daughter moves to a **neighbouring rung**, or stays — the rate-category clock. It is a knob rather
    than a model of its own because the model is unchanged: a daughter starts from its parent and is
    perturbed. ``None`` is the continuous form.

    """

    per: str
    spread: float
    bins: int | None = None

    @property
    def reads(self) -> tuple[str, str]:      # type: ignore[override]
        return (INHERITED, self.per)

    def __post_init__(self) -> None:
        if self.per not in UNITS:
            raise ValueError(
                f"unknown unit {self.per!r}; a value is attached to one of {list(UNITS)}")
        if isinstance(self.spread, bool) or not isinstance(self.spread, (int, float)):
            raise TypeError(f"Inherited spread must be a real number, got {self.spread!r}")
        if not math.isfinite(self.spread) or self.spread < 0:
            raise ValueError(f"Inherited spread must be finite and non-negative, got {self.spread!r}")
        if self.bins is not None:
            if isinstance(self.bins, bool) or not isinstance(self.bins, int):
                raise TypeError(f"Inherited bins must be a whole number, got {self.bins!r}")
            if self.bins < 2:
                raise ValueError(f"Inherited bins must be at least 2, got {self.bins}")

    def __repr__(self) -> str:
        """The written form a run's log records and a reader pastes back into a flag, so it has to be
        an expression that reproduces the rate. ``bins`` is omitted when it is unset, because the
        parser refuses ``bins=None``: a written form that will not parse is not a written form."""
        extra = f", bins={self.bins}" if self.bins is not None else ""
        return f"Inherited(per={self.per!r}, spread={self.spread!r}{extra})"

    def factor(self, **_: Any) -> float:
        """A carried modifier has no factor to compute. Its number is drawn when the unit is born and
        kept by the engine, which hands it back through `Rate.effective`'s ``carried_factor`` — so asking
        for one here, without that stored value, could only ever return a plausible 1.0 for a rate
        that should have varied."""
        raise NotImplementedError(
            f"{type(self).__name__} is carried, not computed: the engine draws its value per "
            f"{self.per} and passes it to Rate.effective(carried_factor=...). It has no factor of its own.")

    def _ladder(self) -> list[float]:
        """The ``bins`` rungs, geometric in ``spread`` and scaled so their mean is 1.

        The mean is taken over the rungs because a nearest-neighbour walk with reflecting ends is
        uniform at stationarity, so an even ladder gives ``E[factor] = 1`` — the same promise the
        continuous form keeps with its lognormal correction."""
        n = self.bins or 0
        raw = [math.exp(self.spread * (k - (n - 1) / 2)) for k in range(n)]
        mean = sum(raw) / n
        return [r / mean for r in raw]

    def initial(self) -> float:
        """The root's factor: 1.0 for the continuous form; the middle rung when binned."""
        if self.bins is None:
            return 1.0
        rungs = self._ladder()
        return rungs[len(rungs) // 2]

    def descend(self, parent_value: float, rng) -> float:
        """A daughter's factor: the parent's, times one mean-corrected lognormal step — or, when
        binned, its parent's rung or one either side, with the ends reflecting.

        The parent's rung is recovered by nearest match rather than by inverting the ladder, so no
        float drift down a deep tree can put a lineage between two rungs."""
        sigma = self.spread
        if self.bins is None:
            return parent_value * math.exp(rng.normal(-0.5 * sigma * sigma, sigma))
        rungs = self._ladder()
        here = min(range(len(rungs)), key=lambda i: abs(math.log(rungs[i] / parent_value)))
        step = int(rng.integers(-1, 2))                      # -1, 0 or +1
        return rungs[min(max(here + step, 0), len(rungs) - 1)]


class Drawn(Modifier):
    """A value **drawn once per unit** and then fixed for that unit's life — i.i.d. heterogeneity
    with no memory. ``per`` names the unit, so one class covers every cell::

        loss = 0.25 * Drawn(per="family", spread=0.5)          # family rate heterogeneity
        substitution = 1.0 * Drawn(per="lineage", spread=0.3)  # the uncorrelated relaxed clock

    ``spread`` is the common case and means a **lognormal** of that log-scale σ. For any other shape,
    pass a distribution instead — any built-in `Distribution`:
    ``Fixed``, ``Exponential``, ``Gamma``, ``LogNormal``, ``Uniform``, ``Geometric``. A callable or a
    scipy frozen distribution is refused here, because neither states a mean to normalise by; an
    extent takes both, being a size rather than a multiplier::

        loss = 0.25 * Drawn(per="family", dist=Gamma(shape=4.0, scale=0.25))

    Give one or the other, never both: ``spread=σ`` *is* ``dist=LogNormal(0.0, σ)``.

    **Every draw is normalised to mean 1**, by dividing by the distribution's own mean. A drawn value
    is a *multiplier*, and one that does not average to 1 changes what the base means — a base of 0.25
    would stop being the average rate. So widening the spread separates units and leaves the average
    one where you put it, and a distribution's **location is normalised away**: what it contributes is
    its shape. ``Exponential(1.0)`` and ``Exponential(7.0)`` are therefore the same modifier. A
    distribution that cannot state its mean — a bare callable, a scipy frozen distribution — is
    refused rather than normalised by a guess. (A number that *is* the rate rather than a factor is
    `SetBy`, where nothing is normalised.)

    **Writing one object or two decides what varies together.** One `Drawn` read by several rates is
    one draw for that unit, so a family that loses fast also duplicates fast; two separately built
    ones are independent even with identical arguments (SPEC §5).
    """

    def __init__(self, *, per: str, spread: float | None = None, dist: object = None) -> None:
        if per not in UNITS:
            raise ValueError(f"unknown unit {per!r}; a value is attached to one of {list(UNITS)}")
        if (spread is None) == (dist is None):
            raise ValueError(
                "Drawn takes a spread or a dist, not both and not neither: spread=0.5 is a lognormal "
                "of that log-scale, and dist=Gamma(...) is any distribution, normalised to mean 1.")
        if spread is not None:
            if isinstance(spread, bool) or not isinstance(spread, (int, float)) \
                    or not math.isfinite(spread) or spread < 0:
                raise ValueError(
                    f"Drawn spread must be a finite non-negative number, got {spread!r}")
            spread = float(spread)
            dist = LogNormal(0.0, spread)
        self.per = per
        self.spread = spread          # kept only so the written form can use the short spelling
        self.dist = as_distribution(dist)
        mean = self.dist.mean()       # raises for a distribution that cannot state one
        if not math.isfinite(mean) or mean <= 0:
            raise ValueError(
                f"a drawn multiplier is normalised by its distribution's mean, which must be finite "
                f"and positive; {self.dist!r} has a mean of {mean!r}")

    @property
    def reads(self) -> tuple[str, str]:      # type: ignore[override]
        return (DRAWN, self.per)

    def factor(self, **_: Any) -> float:
        """A carried modifier has no factor to compute — see `Inherited.factor`."""
        raise NotImplementedError(
            f"Drawn is carried, not computed: the engine draws its value per {self.per} and passes "
            f"it to Rate.effective(carried_factor=...). It has no factor of its own.")

    def draw(self, rng) -> float:
        """One independent multiplier for a unit, with mean 1. A zero spread gives exactly 1.0 and
        consumes no randomness, so switching heterogeneity off leaves every other draw where it was."""
        if self.spread == 0.0:
            return 1.0
        return self.dist.sample(rng) / self.dist.mean()

    def __repr__(self) -> str:
        """The written form a run's log records and a reader pastes back into a flag — in the short
        ``spread=`` spelling where that is how it was written."""
        inner = f"spread={self.spread!r}" if self.spread is not None else f"dist={self.dist!r}"
        return f"Drawn(per={self.per!r}, {inner})"

    def __eq__(self, other: object) -> bool:
        """Two draws are the same modifier when they are on the same unit with the same
        distribution. ``spread`` is not compared: it records only *how it was written*, and
        ``spread=0.5`` and ``dist=LogNormal(0.0, 0.5)`` are one modifier spelled two ways."""
        return isinstance(other, Drawn) and other.per == self.per and other.dist == self.dist

    def __hash__(self) -> int:
        return hash((Drawn, self.per, self.dist))


class DrivenBy(Modifier):
    """The factor is read from **another evolved value** — the one mechanism behind both conditioning and
    joining (SPEC §2).

    It is Ch2's definition made literal: *a rate that reads a value which varies from lineage to
    lineage, rather than a fixed number*. ``DrivenBy`` reads the driver's value on each lineage and
    the mapping turns it into a number — a multiplier on a rate or an extent, a weight on
    ``transfer_to``::

        loss = 0.25 * mod.DrivenBy("habitat.tsv", {"aquatic": 3.0, "terrestrial": 1.0})
        birth = 1.0 * mod.DrivenBy("trait", {"small": 1.0, "large": 2.0})   # a joint model

    ``driver`` says where the driven value comes from, and that single choice splits *conditioned*
    from *joint* — the chapter's spine, *can the driver be grown first?*:

    - a **filename** (``"habitat.tsv"``), a **grown driver result** (a ``TraitsResult``, discrete or
      continuous), or a genome result's ``presence(...)`` / ``completion(...)`` —
      the driver was grown first and handed over (**conditioned**): two ordinary runs. The result
      object is the file's in-memory shortcut — same conditioning, no ``write``/read step;
    - a **level name** (``"trait"``, ``"genomes:count"``) — the driver co-evolves in one run
      (**joint**): neither level can be grown first.

    ``mapping`` says how the driver's value becomes the factor — a `Table`
    (a dict, for a discrete driver), a `Curve` (a callable, continuous),
    a `Scalar` (a log-link coefficient), or a `Between` (a weight per donor/recipient pair, which
    only ``transfer_to`` takes); a raw dict / callable / number is coerced (`as_mapping()`).

    What a ``DrivenBy`` can be attached to comes in **three kinds**, and only the first is a rate:
    *how often* an event fires (a rate, e.g. ``loss``), *how much* it takes (an extent, e.g.
    ``loss_extent``, at the ordered and nucleotide resolutions), and a **choice** of who receives it
    (``transfer_to``, a weight per candidate rather than a multiplier). It always maps a value to a
    number; it never drives a *value*, such as an OU optimum.

    Like a carried modifier, ``DrivenBy`` reads
    a value the **engine** threads per lineage — here a ``drivers`` mapping ``{key: value}`` — and
    is otherwise dumb: it just maps the value to a factor. The engine owns *where* the value comes from
    (a file it loaded, or the live level growing beside the tree) and *when* it changes (a discrete
    driver switches mid-branch, so the engine steps its Gillespie at each switch); a rate reaching an
    engine that has not threaded its ``driver`` gets a factor of 1.0 (inert).
    """

    #: A driven value is *not* carried: the engine resolves it per lineage into ``drivers`` and this
    #: modifier maps it, whether it was recorded beforehand (conditioned) or is growing alongside
    #: (joint). Which of those it is depends on the ``driver`` argument, not on the class, so the
    #: kind here is the pair's shared name rather than one of them.
    reads: ClassVar[tuple[str, str] | None] = (DRIVEN, "lineage")

    def __init__(self, driver: object, mapping: object, step: float | None = None) -> None:
        from .mapping import as_mapping

        if isinstance(driver, str):
            if not driver.strip():
                raise ValueError("DrivenBy driver must be a non-empty string (a filename or level name)")
            base: object = driver                    # a string driver is its own context key
        else:
            base = id(driver)                        # an in-memory driver result (conditioning): key by identity
        if step is not None:
            step = float(step)
            if not (step > 0.0) or step == float("inf"):
                raise ValueError(
                    f"DrivenBy step is the resolution a CONTINUOUS driver is read at, in the tree's own "
                    f"time units, so it must be finite and positive; got {step!r}.")
        # the step is part of the key: the same driver read at two resolutions is two trajectories, and
        # keying on the driver alone would silently resolve it once and share the first one
        self.key: object = base if step is None else (base, step)
        self.driver = driver
        self.step = step
        self.mapping = as_mapping(mapping)

    def factor(self, *, drivers: Mapping | None = None, **_: Any) -> float:
        """The mapped multiplier for this lineage's driver value — the engine threads the value under
        ``drivers[key]`` (``key`` is the driver string, or the identity of an in-memory driver). No
        ``drivers`` (or this driver absent) ⇒ 1.0, so an unthreaded rate is inert (the engine is
        responsible for supplying the value where a driven rate is supported)."""
        if drivers is None:
            return 1.0
        value = drivers.get(self.key)
        if value is None:
            return 1.0
        return self.mapping.multiplier(value)

    def __repr__(self) -> str:
        return f"DrivenBy({_driver_form(self.driver)}, {self.mapping!r})"

    def __eq__(self, other: object) -> bool:
        return (isinstance(other, DrivenBy) and other.key == self.key
                and other.mapping == self.mapping)

    def __hash__(self) -> int:
        # by key only (a mapping — a dict or callable — need not be hashable); equal DrivenBy share a
        # key, so this stays consistent with __eq__ and keeps a Rate carrying it hashable.
        return hash((DrivenBy, self.key))


#: The names a rate may be **written** with — what `zombi2.rates.parse` whitelists, and the only
#: things in this module a user ever calls. Kept explicit rather than derived from ``__all__``, which
#: also carries the helpers an engine uses: a whitelist that grows whenever a helper is exported is
#: a whitelist that stops meaning anything.
WRITABLE = ("OnTime", "OnTotalDiversity", "Drawn", "Inherited", "DrivenBy", "SetBy")

class SetBy(DrivenBy):
    """**Replace** the parameter's base with a value read from a driver, rather than multiplying it::

        loss = SetBy(habitat, {"cave": 1.0, "surface": 0.25})   # the rate itself, per state

    Written with **no base in front**, because there is none to write: the driver supplies the whole
    number, in the parameter's own units rather than as a dimensionless factor. That is what the
    literature usually means — "the loss rate is 1.0 in caves", not "four times a background nobody
    stated" — and spelling an absolute statement as a multiple of an invented background is the kind
    of quiet mismatch this grammar exists to avoid.

    **The scope still applies.** ``SetBy`` replaces the base, not the *per what?*: a per-copy rate set
    to 1.0 is still 1.0 per copy, so it is multiplied by the copies present exactly as a written base
    would be. Only the number changes.

    It is a `DrivenBy`, so every engine that resolves drivers resolves this one too — the trajectory,
    the mid-branch switches, the mapping checks are all the same machinery. What differs is one line
    in `Rate.effective`, which asks a ``SetBy`` for the base and every
    other modifier for a factor. The two compose: a replaced base may still be scaled.

        loss = SetBy(habitat, {...}) * ScaledBy(size, Scalar(0.5))

    A rate may carry **one** ``SetBy``. Two would be two answers to the same question, and neither
    order of application is more right than the other, so it raises rather than picking.
    """

    def __rmul__(self, other: object):
        if isinstance(other, (int, float)) and not isinstance(other, bool):
            from .rate import RateCompositionError
            raise RateCompositionError(
                f"SetBy replaces the base, so there is no base to write in front of it: use "
                f"`SetBy(driver, mapping)` on its own rather than `{other!r} * SetBy(...)`. If you "
                f"meant to scale a base you state yourself, that is ScaledBy.")
        return super().__rmul__(other)

    def __mul__(self, other: object):
        if isinstance(other, SetBy):
            from .rate import RateCompositionError
            raise RateCompositionError(
                "a rate carries one SetBy: each of two would claim to be the whole number, and no "
                "order of application is more right than the other. Keep one; if you meant to scale "
                "the result, that is ScaledBy, which multiplies and composes freely.")
        return super().__mul__(other)

    def __repr__(self) -> str:
        return f"SetBy({_driver_form(self.driver)}, {self.mapping!r})"


__all__ = ["Modifier", "OnTime", "OnTotalDiversity", "Drawn", "Inherited", "DrivenBy", "SetBy",
           "MEASURED", "DRAWN", "INHERITED", "DRIVEN", "CARRIED_KINDS", "UNITS",
           "values_at_birth", "values_at_split", "check_one_memory",
           "cell_name", "describe", "is_implemented", "matches_declared", "WRITABLE"]
