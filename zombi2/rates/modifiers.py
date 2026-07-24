"""Rate modifiers — the context multipliers of a rate (SPEC §5).

Every rate in ZOMBI2 is ``scope(base) × modifiers``. A *modifier* multiplies a rate
by a dimensionless factor that depends on context — the current time, the standing
diversity, the branch, the family, a driver level. Modifiers **multiply** (that is
the whole difference from scope wrappers, which *wrap*), and the word *"per"* is
reserved for scope, so a modifier never starts with "per".

You reach them through ``mod``::

    birth = 1.0 * mod.OnTime({0: 1.0, 3: 0.3})   # a skyline: 1.0, then 0.3 from time 3 on
    birth = 1.0 * mod.OnTotalDiversity(cap=100)       # slows to 0 as diversity approaches 100

The **deterministic** modifiers (``OnTime``, ``OnTotalDiversity``) have a factor that is a pure function of the
context. The **stochastic** ones additionally carry a draw method the engine drives with a random
generator: ``FromParent`` (the rate drifts parent→child along the tree, via ``initial``/``descend``)
and ``ByLineage`` (one i.i.d. draw per lineage, via ``draw``). Still to come: ``ByFamily`` (i.i.d. per
family) and ``Markov`` (a chain of rate categories). Composition (``*``), which turns
``scope(base) × modifiers`` into a Rate, is the Rate module; here each modifier only knows how to
produce its own factor (or, for the stochastic ones, its own draw).
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Mapping


class Modifier:
    """Base for rate modifiers.

    A modifier reads the context keys it cares about (``time``, ``diversity``,
    ``branch``, ``family``, …) and returns a dimensionless, non-negative multiplier;
    it ignores the rest. Abstract — use a subclass.
    """

    def factor(self, **context: float) -> float:
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

    def factor(self, *, time: float, **_: float) -> float:
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

    cap: float

    def __post_init__(self) -> None:
        if isinstance(self.cap, bool) or not isinstance(self.cap, (int, float)):
            raise TypeError(f"OnTotalDiversity cap must be a real number, got {self.cap!r}")
        if not math.isfinite(self.cap) or self.cap <= 0:
            raise ValueError(f"OnTotalDiversity cap must be finite and positive, got {self.cap!r}")

    def factor(self, *, diversity: float, **_: float) -> float:
        return max(0.0, 1.0 - diversity / self.cap)


@dataclass(frozen=True)
class FromParent(Modifier):
    """The rate drifts along the tree — each lineage inherits its parent's rate times a random
    factor drawn at the split (geometric Brownian motion on the rate: clade drift at the species
    level, the autocorrelated clock at the sequence level). ``spread`` (σ) sets the drift width.

    The per-split factor is lognormal, **mean-corrected** so ``E[factor] = 1``. Without the
    correction the rate inflates down the tree (``E[rate] ≈ e^{σ²/2}`` instead of 1) — a real
    historical bug. The draw logic (:meth:`initial` / :meth:`descend`) is driven by the engine,
    which threads each lineage's current factor and passes it back to :meth:`factor` as
    ``inherited``.
    """

    spread: float

    def __post_init__(self) -> None:
        if isinstance(self.spread, bool) or not isinstance(self.spread, (int, float)):
            raise TypeError(f"FromParent spread must be a real number, got {self.spread!r}")
        if not math.isfinite(self.spread) or self.spread < 0:
            raise ValueError(f"FromParent spread must be finite and non-negative, got {self.spread!r}")

    def initial(self) -> float:
        """The root's factor: 1.0 — the rate starts at its base."""
        return 1.0

    def descend(self, parent_value: float, rng) -> float:
        """A daughter's factor: the parent's, times one mean-corrected lognormal step."""
        sigma = self.spread
        return parent_value * math.exp(rng.normal(-0.5 * sigma * sigma, sigma))

    def factor(self, *, inherited: float = 1.0, **_: float) -> float:
        """The lineage's current factor — the engine threads it and passes it back as ``inherited``."""
        return inherited


@dataclass(frozen=True)
class ByLineage(Modifier):
    """The rate varies independently from lineage to lineage — an *uncorrelated* ("relaxed") clock.

    Each lineage draws **one** i.i.d. multiplier with **no memory** of its parent (contrast
    :class:`FromParent`, whose rate drifts parent→child). The draw is **mean-corrected** so
    ``E[factor] = 1`` — without it the mean rate inflates down the tree (the historical lognormal-clock
    bug). ``spread`` (σ) sets the width; ``dist`` is ``"lognormal"`` (default; σ = the log-scale) or
    ``"gamma"`` (σ = the coefficient of variation) — the two agree to first order in σ.

    At the sequence level this is the lineage clock: the engine draws one value per **species lineage**
    (via :meth:`draw`) and shares it across every gene family passing through that lineage.
    It is the lineage-twin of the genome level's ``ByFamily`` — the same
    i.i.d.-heterogeneity idea, by lineage instead of by family. (A fully per-gene-tree-branch clock is
    the deferred ``ByBranch``.)
    """

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

    def factor(self, *, bylineage: float = 1.0, **_: float) -> float:
        """The lineage's drawn factor — the engine threads it and passes it back as ``bylineage``."""
        return bylineage


@dataclass(frozen=True)
class ByFamily(Modifier):
    """The rate varies independently from gene family to gene family.

    The family-twin of :class:`ByLineage`, and the same i.i.d.-heterogeneity idea: each **family**
    draws one multiplier with no memory, mean-corrected so ``E[factor] = 1`` — so widening ``spread``
    spreads the families out without moving the average one off the base rate. ``dist`` is
    ``"lognormal"`` (default; σ = the log-scale) or ``"gamma"`` (σ = the coefficient of variation).

    **Where you put it decides what varies together**. On a single rate, that rate
    varies by family on its own::

        loss = 0.25 * mod.ByFamily(spread=0.5)      # a family that loses fast is not thereby
        duplication = 0.2 * mod.ByFamily(spread=0.5)   # duplicating fast — independent draws

    In the family-wide ``family_speed=`` slot, one draw scales **every** rate that family has, so a
    fast family is fast at everything::

        simulate_genomes_unordered(tree, duplication=0.2, loss=0.25,
                                   family_speed=mod.ByFamily(spread=0.5))

    The two compose: a family-wide tempo, plus extra variation on one rate.

    Not accepted on ``origination``, which is the rate at which families are *created* — at the moment
    it is read there is no family to have drawn a factor for. The engine rejects it rather than
    quietly ignoring it.
    """

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

    def factor(self, *, byfamily: float = 1.0, **_: float) -> float:
        """The family's drawn factor — the engine threads it and passes it back as ``byfamily``."""
        return byfamily


class DrivenBy(Modifier):
    """The rate is **driven by another level** — the one coupling mechanism (SPEC §2).

    A coupling is Ch2's definition made literal: *a parameter that reads its value from another
    level instead of a number you type*. ``DrivenBy`` reads the driver's value on each lineage and
    multiplies the base rate by the mapped factor::

        loss = 0.25 * mod.DrivenBy("habitat.tsv", {"aquatic": 3.0, "terrestrial": 1.0})
        birth = 1.0 * mod.DrivenBy("trait", {"small": 1.0, "large": 2.0})   # a joint model

    ``source`` says where the driver comes from, and that single choice splits *conditioned* from
    *joint* — the chapter's spine, *can the driver be grown first?*:

    - a **filename** (``"habitat.tsv"``) or a **grown driver result** (a discrete ``TraitsResult``) —
      the driver was grown first and handed over (**conditioned**): two ordinary runs. The result
      object is the file's in-memory shortcut — same conditioning, no ``write``/read step;
    - a **level name** (``"trait"``, ``"genomes:count"``) — the driver co-evolves in one run
      (**joint**): neither level can be grown first.

    ``mapping`` says how the driver's value becomes the factor — a :class:`~zombi2.rates.mapping.Table`
    (a dict, for a discrete driver), a :class:`~zombi2.rates.mapping.Curve` (a callable, continuous),
    or a :class:`~zombi2.rates.mapping.Scalar` (a log-link coefficient); a raw dict / callable / number
    is coerced (:func:`~zombi2.rates.mapping.as_mapping`).

    Like :class:`FromParent` (``inherited``) and :class:`ByLineage` (``bylineage``), ``DrivenBy`` reads
    a value the **engine** threads per lineage — here a ``drivers`` mapping ``{source: value}`` — and
    is otherwise dumb: it just maps the value to a factor. The engine owns *where* the value comes from
    (a file it loaded, or the live level growing beside the tree) and *when* it changes (a discrete
    driver switches mid-branch, so the engine steps its Gillespie at each switch); a rate reaching an
    engine that has not threaded its ``source`` gets a factor of 1.0 (inert). ``DrivenBy`` targets a
    **rate** (a "how often") and **multiplies**; driving a *value* (an OU optimum) is a different verb,
    deferred to experimental for v1.
    """

    def __init__(self, source: object, mapping: object) -> None:
        from .mapping import as_mapping

        if isinstance(source, str):
            if not source.strip():
                raise ValueError("DrivenBy source must be a non-empty string (a filename or level name)")
            self.key: object = source                # a string source is its own context key
        else:
            self.key = id(source)                    # an in-memory driver result (conditioning): key by identity
        self.source = source
        self.mapping = as_mapping(mapping)

    def factor(self, *, drivers: Mapping | None = None, **_: float) -> float:
        """The mapped multiplier for this lineage's driver value — the engine threads the value under
        ``drivers[key]`` (``key`` is the source string, or the identity of an in-memory driver). No
        ``drivers`` (or this source absent) ⇒ 1.0, so an unthreaded rate is inert (the engine is
        responsible for supplying the value where the coupling is supported)."""
        if drivers is None:
            return 1.0
        value = drivers.get(self.key)
        if value is None:
            return 1.0
        return self.mapping.multiplier(value)

    def __repr__(self) -> str:
        src = self.source if isinstance(self.source, str) else f"<{type(self.source).__name__}>"
        return f"DrivenBy({src!r}, {self.mapping!r})"

    def __eq__(self, other: object) -> bool:
        return (isinstance(other, DrivenBy) and other.key == self.key
                and other.mapping == self.mapping)

    def __hash__(self) -> int:
        # by key only (a mapping — a dict or callable — need not be hashable); equal DrivenBy share a
        # key, so this stays consistent with __eq__ and keeps a Rate carrying it hashable.
        return hash((DrivenBy, self.key))


__all__ = ["Modifier", "OnTime", "OnTotalDiversity", "FromParent", "ByLineage",
           "ByFamily", "DrivenBy"]
