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
``ByBranch`` / ``ByFamily`` / ``Speed`` (i.i.d. draws), ``Markov`` (a Markov chain of
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


__all__ = ["Modifier", "Time", "Diversity"]
