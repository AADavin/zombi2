"""Internal plumbing: a rate is ``scope(base) × modifiers``, evaluated at run time.

**This is not part of the public API.** Users write a number (``1.0``), a scope wrapper
(``scope.PerLineage(0.25)``), or a product (``1.0 * mod.OnTime({...}) * mod.OnTotalDiversity(cap=100)``).
The ``*`` produces a :class:`Rate` — the glue that *defers* ``base × scope × modifiers`` until
the engine knows the current moment (time, diversity, the branch, the counts) and can multiply
it out. There is no user-facing "Rate" concept; it is the thing a rate expression evaluates to.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from .modifiers import Modifier
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

    def effective(self, **context: float) -> float:
        """The rate *right now*: the scope-applied base times the product of the modifier factors.

        ``context`` carries the current state (``time``, ``diversity``, the counts ``lineages`` /
        ``copies`` / …); the scope reads the count it needs and each modifier the keys it needs.
        Requires a scope — resolve a bare-number rate with :meth:`with_default_scope` first.
        """
        if self.scope is None:
            raise ValueError("this rate has no scope yet; resolve it with with_default_scope(...)")
        value = self.scope.total(**context)
        for m in self.modifiers:
            value *= m.factor(**context)
        return value

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
    """Coerce a user rate spec into a resolved :class:`Rate`, filling the level's default scope.

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
