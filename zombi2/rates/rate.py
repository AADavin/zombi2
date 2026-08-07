"""Internal plumbing: a rate is ``scope(base) × modifiers``, evaluated at run time.

**This is not part of the public API.** Users write a number (``1.0``), a scope wrapper
(``scope.PerLineage(0.25)``), or a product (``1.0 * mod.OnTime({...}) * mod.OnTotalDiversity(cap=100)``).
The ``*`` produces a `Rate` — the glue that *defers* ``base × scope × modifiers`` until
the engine knows the current moment (time, diversity, the branch, the counts) and can multiply
it out. There is no user-facing "Rate" concept; it is the thing a rate expression evaluates to.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any

from .modifiers import CARRIED_KINDS, Modifier
from .scope import Scope


@dataclass(frozen=True)
class Rate:
    """``base × scope × modifiers``, not yet evaluated. Internal — never built by users directly."""

    base: float
    scope: Scope | None = None
    modifiers: tuple[Modifier, ...] = ()

    def with_default_scope(self, default: type[Scope]) -> "Rate":
        """Fill in the level's default scope (per lineage, per copy, …) when none was set explicitly."""
        if self.scope is not None:
            return self
        return Rate(self.base, default(self.base), self.modifiers)

    def effective(self, *, carried: float = 1.0, **context: Any) -> float:
        """The rate *right now*: the scope-applied base times the product of the modifier factors.

        ``context`` carries the current state (``time``, ``diversity``, the counts ``lineages`` /
        ``copies`` / …); the scope reads the count it needs and each modifier the keys it needs.
        Requires a scope — resolve a bare-number rate with `with_default_scope()` first.

        ``carried`` is the product of this rate's **carried** factors for the unit being evaluated —
        the values the engine drew and kept per lineage or per family (`carried`). Those modifiers
        are skipped in the loop below, because their number does not come from the context: the
        engine already holds it and hands it in here, multiplied out. One float rather than a value
        per modifier, so a rate carrying several costs no more to evaluate than one carrying one.
        """
        if self.scope is None:
            raise ValueError("this rate has no scope yet; resolve it with with_default_scope(...)")
        value = self.scope.total(**context)
        for m in self.modifiers:
            reads = getattr(m, "reads", None)
            if reads is not None and reads[0] in CARRIED_KINDS:
                continue  # its factor arrives through `carried`, drawn and kept by the engine
            value *= m.factor(**context)
        return value * carried

    def carried(self, unit: str | None = None) -> tuple[tuple[Modifier, str], ...]:
        """Every modifier on this rate that reads a value the **engine** has to draw and carry,
        paired with the unit it is carried per (see `Modifier.reads`).

        A modifier reading a *measured* value computes its own factor from the context and needs
        nothing from the engine. A *drawn* or *inherited* one does: its number is produced once
        when a unit is born, kept for that unit's life, and handed back at every evaluation, and
        only the engine can do that. This is the one query for finding them, so a level does not
        have to know which modifier classes exist to thread them.

        ``unit`` narrows the answer to one kind of unit (``"lineage"``, ``"family"``). The result
        keeps the order the modifiers were written in, and it keeps **all** of them — a rate
        carrying two drawn values answers with two.
        """
        found = []
        for m in self.modifiers:
            reads = getattr(m, "reads", None)
            if reads is None or reads[0] not in CARRIED_KINDS:
                continue
            if unit is None or reads[1] == unit:
                found.append((m, reads[1]))
        return tuple(found)

    def next_change(self, time: float) -> float:
        """The next time a component of this rate changes on its own — the earliest skyline
        breakpoint across its modifiers. ``inf`` if the rate never changes with time."""
        nc = math.inf
        for m in self.modifiers:
            nc = min(nc, m.next_change(time))
        return nc

    def __mul__(self, other: object):
        if isinstance(other, Modifier):
            return Rate(self.base, self.scope, self.modifiers + (other,))
        return NotImplemented

    __rmul__ = __mul__  # a number/scope on the left is handled there; only Modifier*Rate reaches here


def as_rate(spec: object, *, default_scope: type[Scope]) -> Rate:
    """Coerce a user rate spec into a resolved `Rate`, filling the level's default scope.

    Accepts a number, a scope wrapper, a modifier (product), or an already-built ``Rate``.
    """
    if isinstance(spec, Rate):
        return spec.with_default_scope(default_scope)
    if isinstance(spec, Scope):
        return Rate(spec.base, spec, ())
    if isinstance(spec, Modifier):
        return Rate(1.0, None, (spec,)).with_default_scope(default_scope)
    if isinstance(spec, bool) or not isinstance(spec, (int, float)):
        raise TypeError(
            f"a rate must be a number, a scope wrapper, or a modifier product, got {spec!r}"
        )
    return Rate(float(spec)).with_default_scope(default_scope)


__all__ = ["Rate", "as_rate"]
