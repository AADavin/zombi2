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

from .modifiers import CARRIED_KINDS, Modifier, SetBy
from .scope import Scope


class RateCompositionError(TypeError):
    """A `*` the grammar itself refuses, with a message written for that mistake.

    A `TypeError` so nothing that catches one stops catching it, and its own class so the written
    form can tell it from CPython's ``unsupported operand type(s)``: the parser re-raises ours
    verbatim, because the sentence says what to write instead, and answers CPython's generically,
    because "Rate and float" is about types rather than about the rate.
    """


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

    def effective(self, *, carried_factor: float = 1.0, **context: Any) -> float:
        """The rate *right now*: the scope-applied base times the product of the modifier factors.

        ``context`` carries the current state (``time``, ``diversity``, the counts ``lineages`` /
        ``copies`` / …); the scope reads the count it needs and each modifier the keys it needs.
        Requires a scope — resolve a bare-number rate with `with_default_scope()` first.

        ``carried_factor`` is the product of the values the engine drew and kept for the unit being
        evaluated — per lineage, per family (`carried_modifiers`). Those modifiers
        are skipped in the loop below, because their number does not come from the context: the
        engine already holds it and hands it in here, multiplied out. One float rather than a value
        per modifier, so a rate carrying several costs no more to evaluate than one carrying one.
        """
        if self.scope is None:
            raise ValueError("this rate has no scope yet; resolve it with with_default_scope(...)")
        base = self.base
        for m in self.modifiers:
            if isinstance(m, SetBy):
                base = m.factor(**context)   # the driver supplies the number itself, not a factor
        value = self.scope.total_of(base, **context)
        for m in self.modifiers:
            if isinstance(m, SetBy):
                continue  # already used, as the base
            reads = getattr(m, "reads", None)
            if reads is not None and reads[0] in CARRIED_KINDS:
                continue  # its factor arrives through `carried_factor`, drawn and kept by the engine
            value *= m.factor(**context)
        return value * carried_factor

    def check_one_base(self, label: str = "this rate") -> None:
        """A rate may carry **one** `SetBy`. Two would each claim to *be* the
        base, and no order of application is more right than another, so this raises rather than
        letting the last one written win in silence.

        And none of its modifiers may be a `Weights`: that verb says the numbers are compared
        against each other and normalised, which only a **choice** does — on a rate they would
        silently behave as factors, which is the other model. The class is the same either way, so
        before the verbs there was nothing here to check."""
        set_by = [m for m in self.modifiers if isinstance(m, SetBy)]
        if len(set_by) > 1:
            raise ValueError(
                f"{label} carries {len(set_by)} SetBy modifiers, and a base can only be replaced "
                f"once — each of them claims to be the whole number. Keep one; if you meant to scale "
                f"the result, that is ScaledBy, which multiplies and composes freely.")
        if any(getattr(m, "verb", None) == "Weights" for m in self.modifiers):
            raise ValueError(
                f"{label} carries Weights, which weights the candidates of a choice against each "
                f"other — transfer_to is the only one. On a rate the number multiplies a base, so "
                f"the verb is ScaledBy: the same driver and the same mapping, read as a factor.")

    def carried_modifiers(self, unit: str | None = None) -> tuple[tuple[Modifier, str], ...]:
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
        if isinstance(other, SetBy):
            raise RateCompositionError(
                "SetBy replaces the base, so it cannot follow one: everything to its left — the "
                "number, the scope, any factors — is a base it would silently discard. Write it "
                "first instead: SetBy(driver, mapping) * ScaledBy(...). `0.25 * SetBy(...)` is "
                "refused for the same reason, and this is that hole one operand further along.")
        if isinstance(other, Modifier):
            return Rate(self.base, self.scope, self.modifiers + (other,))
        return NotImplemented

    __rmul__ = __mul__  # a number/scope on the left is handled there; only Modifier*Rate reaches here


def as_rate(spec: object, *, default_scope: type[Scope], label: str = "this rate") -> Rate:
    """Coerce a user rate spec into a resolved `Rate`, filling the level's default scope.

    Accepts a number, a scope wrapper, a modifier (product), or an already-built ``Rate``.

    Every level coerces its rates through here, which is why the one-base rule is checked here
    rather than in each level's own validation: a rule enforced by whoever remembers to call it is a
    rule three levels did not have.
    """
    if isinstance(spec, Rate):
        spec.check_one_base(label)
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
