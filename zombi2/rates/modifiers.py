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
function of the context. The **stochastic** ones also carry a draw the engine drives with a random
generator: ``FromParent`` (the rate drifts parent→child, via ``initial``/``descend``), ``ByLineage``
(one i.i.d. draw per lineage) and ``ByFamily`` (one per family), the last two via ``draw``.

Composition (``*``) belongs to the Rate module; a modifier here knows only how to produce its own
factor, or its own draw.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any, ClassVar, Mapping

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
#: `Rate.carried` is the query that finds them.
CARRIED_KINDS = (DRAWN, INHERITED)

#: The units a value can be attached to, in the order they nest — the other half of `Modifier.reads`.
#: A value's unit decides what may read it: a parameter may read a value when the parameter's own
#: units include the value's (SPEC §5), so a trait on a lineage can drive gene loss, and a family's
#: tempo cannot drive speciation. A unit here that no engine carries is a cell nobody has built, not
#: a different kind of model.
UNITS = ("run", "lineage", "chromosome", "family", "copy", "site")

#: The cells that have a name of their own, because each names a model people cite. A cell without an
#: entry here is written the general way — ``Drawn(per="chromosome", …)`` — which is the point: a new
#: unit needs no name invented for it.
_CELL_NAMES = {(DRAWN, "family"): "ByFamily",
               (DRAWN, "lineage"): "ByLineage",
               (INHERITED, "lineage"): "FromParent"}


def draw_product(mods: "tuple[Modifier, ...]", rng,
                 drawn: "dict[int, float] | None" = None) -> float:
    """One draw from each modifier, multiplied — the factor a newly created unit carries.

    For a unit that never splits, the carried value is fixed for its whole life and only the product
    is ever needed, so the engine keeps one number per unit (a gene family). Where a unit *does*
    split, the engine keeps the factors separately instead, because `FromParent` has to
    nudge its parent's own number rather than a product.

    Drawing in written order is what keeps a run reproducible, and drawing from **every** modifier
    is the point: taking only the first was how a second one silently left the model.

    ``drawn`` makes one value shared. It is a cache for a **single unit**, keyed by modifier
    identity: pass the same dict while drawing each of that unit's rates, and a modifier written on
    two of them is drawn once and both rates get the same number. That is how "a family that loses
    fast also duplicates fast" is said — one object, read twice — against "fast at losing only",
    which is two objects. Two modifiers that merely compare equal are still two draws, because the
    question is whether you wrote one thing or two, not whether their spreads match. Omit the cache
    and every modifier draws for itself.
    """
    out = 1.0
    for m in mods:
        if drawn is None:
            out *= m.draw(rng)
            continue
        key = id(m)
        if key not in drawn:
            drawn[key] = m.draw(rng)
        out *= drawn[key]
    return out


def product(values) -> float:
    """The carried factors multiplied out — what `Rate.effective` takes as
    ``carried``. Empty (a rate carrying nothing per unit) is 1.0."""
    out = 1.0
    for v in values:
        out *= v
    return out


def carried_at_birth(mods: "tuple[Modifier, ...]", rng,
                     drawn: "dict[int, float] | None" = None) -> tuple[float, ...]:
    """The factors a newly created unit carries, one per modifier, in written order.

    An `INHERITED` value starts from its own beginning (`FromParent.initial`); a `DRAWN` one is drawn.
    The dispatch reads `Modifier.reads`, not the class, so a per-unit modifier this
    engine has never heard of is handled like the ones it has. ``drawn`` shares one number between
    the rates of a single unit — see `draw_product`."""
    out = []
    for m in mods:
        key = id(m)
        if drawn is None or key not in drawn:
            value = m.initial() if m.reads and m.reads[0] == INHERITED else m.draw(rng)
            if drawn is None:
                out.append(value)
                continue
            drawn[key] = value
        out.append(drawn[key])
    return tuple(out)


def carried_at_split(mods: "tuple[Modifier, ...]", parent_values: tuple[float, ...], rng,
                     drawn: "dict[int, float] | None" = None) -> tuple[float, ...]:
    """A daughter's carried factors: its parent's, perturbed (`INHERITED`), or a fresh independent
    draw that ignores the parent (`DRAWN`). That one line is the whole autocorrelated / uncorrelated
    split (SPEC §5)."""
    out = []
    for i, m in enumerate(mods):
        key = id(m)
        if drawn is None or key not in drawn:
            value = (m.descend(parent_values[i], rng)
                     if m.reads and m.reads[0] == INHERITED else m.draw(rng))
            if drawn is None:
                out.append(value)
                continue
            drawn[key] = value
        out.append(drawn[key])
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
    return _CELL_NAMES.get(entry, f"{entry[0]} per {entry[1]}")


def describe(m: "Modifier") -> str:
    """What to call one modifier **instance** in a message.

    Its class name, except for the two classes that cover a whole row of the grid: a `Drawn` or an
    `Inherited` is named by its cell, so a refusal says ``ByFamily`` where that cell has a name and
    ``drawn per chromosome`` where it does not. Saying "carries Drawn" would be true and useless,
    since the whole question is *per what*."""
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
    return matches_declared(m, engines) or engine in getattr(m, "implemented_for", ())


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
    #: and `Rate.carried` is how an engine asks for them without knowing which classes exist.
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
    #: ``traits.continuous``  ``time``, ``lineages``, ``diversity``, ``inherited``, ``drivers``
    #: ``traits.discrete``    ``time``, ``lineages``, ``drivers``
    #: ``joint``              ``time``, ``lineages``, ``diversity``, ``drivers``
    #: =====================  =================================================================
    #:
    #: The genome engines thread ``drivers`` only when some rate or extent in the run is driven, so a
    #: modifier that reads it must default its key.
    #:
    #: **``sequences`` is not on that list, deliberately.** Every engine above evaluates its rate
    #: through `Rate.effective`, which multiplies in whatever `factor` returns. The sequence level
    #: reads its two kinds of modifier itself — the clock is *drawn per lineage* before any site
    #: evolves, not evaluated at an event — so a modifier declaring itself implemented there would be
    #: accepted and then never called, which is the silence this whole mechanism exists to prevent.
    #: It refuses instead, and says why.
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
        inner = ", ".join(f"{t:g}: {f:g}" for t, f in self._steps)
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
    threading each unit's current factor back through `Rate.effective`'s ``carried``.

    ``bins`` discretises the drift: the rate takes one of ``bins`` values on a geometric ladder and a
    daughter moves to a **neighbouring rung**, or stays — the rate-category clock. It is a knob rather
    than a model of its own because the model is unchanged: a daughter starts from its parent and is
    perturbed. ``None`` is the continuous form.

    `FromParent` is this with ``per="lineage"``, kept as the name of that cell.
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
        an expression that reproduces the rate. The named cell is used where one exists, because that
        is the spelling the parser and the manual have always shown."""
        extra = f", bins={self.bins}" if self.bins is not None else ""
        if self.per == "lineage":
            return f"FromParent(spread={self.spread!r}{extra})"
        return f"Inherited(per={self.per!r}, spread={self.spread!r}{extra})"

    def factor(self, **_: Any) -> float:
        """A carried modifier has no factor to compute. Its number is drawn when the unit is born and
        kept by the engine, which hands it back through `Rate.effective`'s ``carried`` — so asking
        for one here, without that stored value, could only ever return a plausible 1.0 for a rate
        that should have varied."""
        raise NotImplementedError(
            f"{type(self).__name__} is carried, not computed: the engine draws its value per "
            f"{self.per} and passes it to Rate.effective(carried=...). It has no factor of its own.")

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


@dataclass(frozen=True)
class Drawn(Modifier):
    """A value **drawn once per unit** and then fixed for that unit's life — i.i.d. heterogeneity
    with no memory. ``per`` names the unit, so one class covers every cell::

        loss = 0.25 * Drawn(per="family", spread=0.5)          # family rate heterogeneity
        substitution = 1.0 * Drawn(per="lineage", spread=0.3)  # the uncorrelated relaxed clock

    Each draw is **mean-corrected** so ``E[factor] = 1`` — widening ``spread`` spreads the units out
    without moving the average one off the base. ``dist`` is ``"lognormal"`` (σ = the log-scale) or
    ``"gamma"`` (σ = the coefficient of variation); the two agree to first order in σ.

    **Writing one object or two decides what varies together.** One `Drawn` read by several rates is
    one draw for that unit, so a family that loses fast also duplicates fast; two separately built
    ones are independent even with identical arguments (SPEC §5).

    `ByFamily` and `ByLineage` are this with ``per`` fixed, kept as the names of those cells. A unit
    no engine carries is refused at that level's gate, by name, rather than ignored.
    """

    per: str
    spread: float
    dist: str = "lognormal"

    @property
    def reads(self) -> tuple[str, str]:      # type: ignore[override]
        return (DRAWN, self.per)

    def __post_init__(self) -> None:
        if self.per not in UNITS:
            raise ValueError(
                f"unknown unit {self.per!r}; a value is attached to one of {list(UNITS)}")
        if isinstance(self.spread, bool) or not isinstance(self.spread, (int, float)) \
                or not math.isfinite(self.spread) or self.spread < 0:
            raise ValueError(f"Drawn spread must be a finite non-negative number, got {self.spread!r}")
        if self.dist not in ("lognormal", "gamma"):
            raise ValueError(f"Drawn dist must be 'lognormal' or 'gamma', got {self.dist!r}")

    def __repr__(self) -> str:
        """The written form — the named cell where one exists, so a run's log and a ``--params`` file
        keep the spelling the parser and the manual have always shown."""
        named = {"family": "ByFamily", "lineage": "ByLineage"}.get(self.per)
        if named is not None:
            return f"{named}(spread={self.spread!r}, dist={self.dist!r})"
        return f"Drawn(per={self.per!r}, spread={self.spread!r}, dist={self.dist!r})"

    def factor(self, **_: Any) -> float:
        """A carried modifier has no factor to compute. Its number is drawn when the unit is born and
        kept by the engine, which hands it back through `Rate.effective`'s ``carried`` — so asking
        for one here, without that stored value, could only ever return a plausible 1.0 for a rate
        that should have varied."""
        raise NotImplementedError(
            f"{type(self).__name__} is carried, not computed: the engine draws its value per "
            f"{self.per} and passes it to Rate.effective(carried=...). It has no factor of its own.")

    def draw(self, rng) -> float:
        """One independent, mean-1 multiplier for a unit. ``spread = 0`` gives 1.0 (no variation)."""
        s = self.spread
        if s == 0.0:
            return 1.0
        if self.dist == "lognormal":
            return math.exp(rng.normal(-0.5 * s * s, s))     # mean-corrected lognormal
        return float(rng.gamma(1.0 / (s * s), s * s))        # mean-1 gamma, coefficient of variation = s


def FromParent(spread: float, bins: int | None = None) -> Inherited:
    """The inherited cell on a lineage — ``Inherited(per="lineage", …)``.

    Kept as its own name because it names a model people cite: clade drift at the species level
    (ClaDS), the autocorrelated clock at the sequence level, variable-rates BM at the trait level."""
    return Inherited(per="lineage", spread=spread, bins=bins)


def ByLineage(spread: float, dist: str = "lognormal") -> Drawn:
    """The drawn cell on a lineage — ``Drawn(per="lineage", …)``, the uncorrelated ("relaxed") clock.

    At the sequence level the engine draws one value per **species lineage** and shares it across
    every gene family passing through that lineage."""
    return Drawn(per="lineage", spread=spread, dist=dist)


def ByFamily(spread: float, dist: str = "lognormal") -> Drawn:
    """The drawn cell on a gene family — ``Drawn(per="family", …)``.

    Not accepted on ``origination``, which is the rate at which families are *created*: at the moment
    it is read there is no family to have drawn a factor for. The engine rejects it rather than
    quietly ignoring it."""
    return Drawn(per="family", spread=spread, dist=dist)


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

    Like `FromParent` (``inherited``) and `ByLineage` (``bylineage``), ``DrivenBy`` reads
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
        drv = self.driver if isinstance(self.driver, str) else f"<{type(self.driver).__name__}>"
        return f"DrivenBy({drv!r}, {self.mapping!r})"

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
WRITABLE = ("OnTime", "OnTotalDiversity", "Drawn", "Inherited",
            "FromParent", "ByLineage", "ByFamily", "DrivenBy", "SetBy")

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
            raise TypeError(
                f"SetBy replaces the base, so there is no base to write in front of it: use "
                f"`SetBy(driver, mapping)` on its own rather than `{other!r} * SetBy(...)`. If you "
                f"meant to scale a base you state yourself, that is ScaledBy.")
        return super().__rmul__(other)

    def __repr__(self) -> str:
        drv = self.driver if isinstance(self.driver, str) else f"<{type(self.driver).__name__}>"
        return f"SetBy({drv!r}, {self.mapping!r})"


__all__ = ["Modifier", "OnTime", "OnTotalDiversity", "Drawn", "Inherited",
           "FromParent", "ByLineage", "ByFamily", "DrivenBy", "SetBy",
           "MEASURED", "DRAWN", "INHERITED", "DRIVEN", "CARRIED_KINDS", "UNITS",
           "product", "draw_product", "carried_at_birth", "carried_at_split", "check_one_memory",
           "cell_name", "describe", "is_implemented", "matches_declared", "WRITABLE"]
