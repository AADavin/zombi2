"""Rate modifiers — the context multipliers of a rate (SPEC §5).

Every rate in ZOMBI2 is ``scope(base) × modifiers``. A *modifier* multiplies a rate
by a dimensionless factor that depends on context — the current time, the standing
diversity, the branch, the family, a driver level. Modifiers **multiply** (that is
the whole difference from scope wrappers, which *wrap*), and the word *"per"* is
reserved for scope, so a modifier never starts with "per".

You reach them through ``mod``::

    birth = 1.0 * mod.Time({0: 1.0, 3: 0.3})   # a skyline: 1.0, then 0.3 from time 3 on
    birth = 1.0 * mod.Diversity(cap=100)       # slows to 0 as diversity approaches 100

This module holds the **deterministic** modifiers, whose factor is a pure function of
the context. The stochastic ones — ``Inherited`` (the rate drifts along the tree),
``ByLineage`` / ``ByFamily`` (i.i.d. draws), ``Markov`` (a Markov chain of
rate categories) — need a random generator and the tree, so they arrive in the next
module. Composition (``*``), which turns ``scope(base) × modifiers`` into a Rate, is
the Rate module; here each modifier only knows how to produce its own factor.
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


class Time(Modifier):
    """The rate changes in time — a skyline / episodic schedule.

    ``schedule`` maps each interval's start time to a relative factor::

        Time({0: 1.0, 3: 0.3})   # factor 1.0 on [0, 3), then 0.3 from time 3 on

    Factors are relative (dimensionless): on a base of ``2.0`` the schedule scales it.
    Before the earliest breakpoint the earliest factor applies (define the schedule
    from time 0 to avoid surprise).
    """

    def __init__(self, schedule: Mapping[float, float]) -> None:
        steps = tuple(sorted((float(t), float(f)) for t, f in schedule.items()))
        if not steps:
            raise ValueError("Time needs a non-empty schedule, e.g. Time({0: 1.0, 3: 0.3})")
        for t, f in steps:
            if not math.isfinite(t):
                raise ValueError(f"Time schedule times must be finite, got {t!r}")
            if not math.isfinite(f) or f < 0:
                raise ValueError(f"Time factors must be finite and non-negative, got {f!r}")
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
        return f"Time({{{inner}}})"

    def __eq__(self, other: object) -> bool:
        return isinstance(other, Time) and other._steps == self._steps

    def __hash__(self) -> int:
        return hash((Time, self._steps))


@dataclass(frozen=True)
class Diversity(Modifier):
    """The rate slows as standing diversity grows — diversity-dependence.

    The factor falls linearly from 1 toward 0 as diversity rises to ``cap`` (a carrying
    capacity), and stays 0 beyond it: ``Diversity(cap=100)`` halves the rate at 50
    lineages and stops it at 100.
    """

    cap: float

    def __post_init__(self) -> None:
        if isinstance(self.cap, bool) or not isinstance(self.cap, (int, float)):
            raise TypeError(f"Diversity cap must be a real number, got {self.cap!r}")
        if not math.isfinite(self.cap) or self.cap <= 0:
            raise ValueError(f"Diversity cap must be finite and positive, got {self.cap!r}")

    def factor(self, *, diversity: float, **_: float) -> float:
        return max(0.0, 1.0 - diversity / self.cap)


@dataclass(frozen=True)
class Inherited(Modifier):
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
            raise TypeError(f"Inherited spread must be a real number, got {self.spread!r}")
        if not math.isfinite(self.spread) or self.spread < 0:
            raise ValueError(f"Inherited spread must be finite and non-negative, got {self.spread!r}")

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


__all__ = ["Modifier", "Time", "Diversity", "Inherited"]
