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


def is_implemented(m: "Modifier", engines: tuple[type, ...], engine: str) -> bool:
    """Whether ``engine`` may run modifier ``m``: it is one of the types that engine threads
    (``engines``, the level's ``IMPLEMENTED_MODIFIERS``), or it names that engine in its own
    `Modifier.implemented_for`. Every engine gate goes through here, so the escape hatch cannot be
    honoured in one level and forgotten in another."""
    return isinstance(m, engines) or engine in getattr(m, "implemented_for", ())


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
class FromParent(Modifier):
    """The rate drifts along the tree — each lineage inherits its parent's rate times a random
    factor drawn at the split (geometric Brownian motion on the rate: clade drift at the species
    level, the autocorrelated clock at the sequence level). ``spread`` (σ) sets the drift width.

    The per-split factor is lognormal, **mean-corrected** so ``E[factor] = 1``. Without the
    correction the rate inflates down the tree (``E[rate] ≈ e^{σ²/2}`` instead of 1) — a real
    historical bug. The draw logic (`initial()` / `descend()`) is driven by the engine,
    which threads each lineage's current factor and passes it back to `factor()` as
    ``inherited``.

    ``bins`` discretises the drift: instead of a continuous step, the rate takes one of ``bins``
    values on a geometric ladder and a daughter moves to a **neighbouring rung**, or stays. That is
    the discrete-bin (rate-category) clock — the same inherit-and-perturb idea, in steps rather than
    continuously, which is what a lineage-category model assumes. It is a knob rather than a modifier
    of its own because the model is `FromParent`'s: a daughter starts from its parent and is
    perturbed. ``None`` (the default) is the continuous form, so a run written before this existed
    draws exactly as it did.
    """

    reads: ClassVar[tuple[str, str] | None] = (INHERITED, "lineage")

    spread: float
    bins: int | None = None

    def __post_init__(self) -> None:
        if isinstance(self.spread, bool) or not isinstance(self.spread, (int, float)):
            raise TypeError(f"FromParent spread must be a real number, got {self.spread!r}")
        if not math.isfinite(self.spread) or self.spread < 0:
            raise ValueError(f"FromParent spread must be finite and non-negative, got {self.spread!r}")
        if self.bins is not None:
            if isinstance(self.bins, bool) or not isinstance(self.bins, int):
                raise TypeError(f"FromParent bins must be a whole number, got {self.bins!r}")
            if self.bins < 2:
                raise ValueError(f"FromParent bins must be at least 2, got {self.bins}")

    def __repr__(self) -> str:
        """Omit ``bins`` when it is unset. The repr is what `written_form()` records in a run's log
        and what a reader pastes back into a flag, so it has to be the expression that reproduces the
        rate — and the dataclass default would render ``bins=None``, which the rate grammar rejects."""
        extra = f", bins={self.bins}" if self.bins is not None else ""
        return f"FromParent(spread={self.spread!r}{extra})"

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

    def factor(self, *, inherited: float = 1.0, **_: Any) -> float:
        """The lineage's current factor — the engine threads it and passes it back as ``inherited``."""
        return inherited


@dataclass(frozen=True)
class ByLineage(Modifier):
    """The rate varies independently from lineage to lineage — an *uncorrelated* ("relaxed") clock.

    Each lineage draws **one** i.i.d. multiplier with **no memory** of its parent (contrast
    `FromParent`, whose rate drifts parent→child). The draw is **mean-corrected** so
    ``E[factor] = 1`` — without it the mean rate inflates down the tree (the historical lognormal-clock
    bug). ``spread`` (σ) sets the width; ``dist`` is ``"lognormal"`` (default; σ = the log-scale) or
    ``"gamma"`` (σ = the coefficient of variation) — the two agree to first order in σ.

    At the sequence level this is the lineage clock: the engine draws one value per **species lineage**
    (via `draw()`) and shares it across every gene family passing through that lineage.
    It is the lineage-twin of the genome level's ``ByFamily``.
    """

    reads: ClassVar[tuple[str, str] | None] = (DRAWN, "lineage")

    spread: float
    dist: str = "lognormal"

    def __post_init__(self) -> None:
        if isinstance(self.spread, bool) or not isinstance(self.spread, (int, float)) \
                or not math.isfinite(self.spread) or self.spread < 0:
            raise ValueError(f"ByLineage spread must be a finite non-negative number, got {self.spread!r}")
        if self.dist not in ("lognormal", "gamma"):
            raise ValueError(f"ByLineage dist must be 'lognormal' or 'gamma', got {self.dist!r}")

    def draw(self, rng) -> float:
        """One independent, mean-1 multiplier for a lineage. ``spread = 0`` gives 1.0 (a strict clock)."""
        s = self.spread
        if s == 0.0:
            return 1.0
        if self.dist == "lognormal":
            return math.exp(rng.normal(-0.5 * s * s, s))     # mean-corrected lognormal
        return float(rng.gamma(1.0 / (s * s), s * s))        # mean-1 gamma, coefficient of variation = s

    def factor(self, *, bylineage: float = 1.0, **_: Any) -> float:
        """The lineage's drawn factor — the engine threads it and passes it back as ``bylineage``."""
        return bylineage


@dataclass(frozen=True)
class ByFamily(Modifier):
    """The rate varies independently from gene family to gene family.

    The family-twin of `ByLineage`, and the same i.i.d.-heterogeneity idea: each **family**
    draws one multiplier with no memory, mean-corrected so ``E[factor] = 1`` — so widening ``spread``
    spreads the families out without moving the average one off the base rate. ``dist`` is
    ``"lognormal"`` (default; σ = the log-scale) or ``"gamma"`` (σ = the coefficient of variation).

    **Where you put it decides what varies together**. On a single rate, that rate
    varies by family on its own::

        loss = 0.25 * mod.ByFamily(spread=0.5)      # a family that loses fast is not thereby
        duplication = 0.2 * mod.ByFamily(spread=0.5)   # duplicating fast — independent draws

    As ``family_speed=``, one draw scales **every** rate that family has, so a
    fast family is fast at everything::

        simulate_genomes_family(tree, duplication=0.2, loss=0.25,
                                family_speed=mod.ByFamily(spread=0.5))

    The two compose: a family-wide tempo, plus extra variation on one rate.

    Not accepted on ``origination``, which is the rate at which families are *created* — at the moment
    it is read there is no family to have drawn a factor for. The engine rejects it rather than
    quietly ignoring it.
    """

    reads: ClassVar[tuple[str, str] | None] = (DRAWN, "family")

    spread: float
    dist: str = "lognormal"

    def __post_init__(self) -> None:
        if isinstance(self.spread, bool) or not isinstance(self.spread, (int, float)) \
                or not math.isfinite(self.spread) or self.spread < 0:
            raise ValueError(f"ByFamily spread must be a finite non-negative number, got {self.spread!r}")
        if self.dist not in ("lognormal", "gamma"):
            raise ValueError(f"ByFamily dist must be 'lognormal' or 'gamma', got {self.dist!r}")

    def draw(self, rng) -> float:
        """One independent, mean-1 multiplier for a family. ``spread = 0`` gives 1.0 (no variation)."""
        s = self.spread
        if s == 0.0:
            return 1.0
        if self.dist == "lognormal":
            return math.exp(rng.normal(-0.5 * s * s, s))     # mean-corrected lognormal
        return float(rng.gamma(1.0 / (s * s), s * s))        # mean-1 gamma, coefficient of variation = s

    def factor(self, *, byfamily: float = 1.0, **_: Any) -> float:
        """The family's drawn factor — the engine threads it and passes it back as ``byfamily``."""
        return byfamily


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


__all__ = ["Modifier", "OnTime", "OnTotalDiversity", "FromParent", "ByLineage",
           "ByFamily", "DrivenBy",
           "MEASURED", "DRAWN", "INHERITED", "DRIVEN", "CARRIED_KINDS"]
