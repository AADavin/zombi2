"""The mapping of a :class:`~zombi2.rates.modifiers.DrivenBy` coupling — how a driver
level's value becomes a rate multiplier (SPEC §2, ``coupling-api.md``).

A coupling is *a rate that reads its value from another level*. The ``DrivenBy`` modifier
reads the driver's value on a lineage; the **mapping** turns that value into the
dimensionless factor the rate is multiplied by. It is the "response" of the old coevolve
grammar, distilled to the three shapes that **multiply a rate**:

- :class:`Table`  — a **discrete** driver → a dict of factors: ``{"aquatic": 3.0, "terrestrial": 1.0}``.
- :class:`Curve`  — a **continuous** driver → a function: ``lambda x: math.exp(0.5 * x)``.
- :class:`Scalar` — a single log-link coefficient: ``multiplier = exp(strength · value)``.

You rarely name these — pass a raw ``dict`` / callable / number as ``mapping=`` and
:func:`as_mapping` coerces it (a dict → ``Table``, a callable → ``Curve``, a number →
``Scalar``), exactly as :func:`~zombi2.rates.rate.as_rate` coerces a rate spec.

The fourth grammar shape, **Jump** (a burst fired *at an event*, e.g. a pulse of gene
change at each split), is **not** a rate multiplier — it changes a state at a moment, not a
"how often" — so it does not live here and is not reachable through ``DrivenBy`` (SPEC §4).
"""

from __future__ import annotations

import math

_MAX_EXPONENT = 40.0  # clamp the log-link argument so a large driver value cannot overflow exp()


class Mapping:
    """Base for a driver-value → rate-multiplier response. Abstract — use :class:`Table`,
    :class:`Curve`, or :class:`Scalar` (or pass a raw dict / callable / number, which
    :func:`as_mapping` coerces). A mapping returns a **dimensionless, non-negative** factor."""

    def multiplier(self, value: object) -> float:
        raise NotImplementedError


class Table(Mapping):
    """A **discrete** driver → a lookup of factors, one per driver state::

        Table({"aquatic": 3.0, "terrestrial": 1.0})   # 3× the rate in aquatic lineages

    ``default`` (1.0) is the factor for any state not named — so an unlisted state leaves the
    rate unchanged. This is the primary ``DrivenBy`` mapping (MuSSE-style per-state rates).

    States are matched by their **string form** — ``Table({0: 3.0, 1: 1.0})`` and ``Table({"0":
    3.0, "1": 1.0})`` behave identically, and both match a driver whose value is ``0`` or ``"0"``.
    A conditioned driver arrives from a text file (always a string), and a live joint driver arrives
    as its native label; string-matching makes the two agree, so an ``int``-labelled trait does not
    silently miss its mapping."""

    def __init__(self, per_state, default: float = 1.0) -> None:
        if not isinstance(per_state, dict) or not per_state:
            raise ValueError(f"Table needs a non-empty {{state: factor}} dict, got {per_state!r}")
        table = {}
        for state, factor in per_state.items():
            key = str(state)  # states matched by string form (a driver file is text); see the class docstring
            if key in table:
                raise ValueError(
                    f"Table states collide as strings: {state!r} and an earlier key both map to {key!r}"
                )
            table[key] = _check_factor(factor, f"Table factor for {state!r}")
        self.per_state = table
        self.default = _check_factor(default, "Table default")

    def multiplier(self, value: object) -> float:
        return self.per_state.get(str(value), self.default)

    def __repr__(self) -> str:
        inner = ", ".join(f"{s!r}: {f:g}" for s, f in self.per_state.items())
        tail = "" if self.default == 1.0 else f", default={self.default:g}"
        return f"Table({{{inner}}}{tail})"

    def __eq__(self, other: object) -> bool:
        return (isinstance(other, Table) and other.per_state == self.per_state
                and other.default == self.default)


class Curve(Mapping):
    """A **continuous** driver → an arbitrary function of the value, optionally capped::

        Curve(lambda x: math.exp(0.5 * x))          # exponential response
        Curve(lambda x: 1.0 + x, bound=5.0)          # linear, capped at 5×

    ``bound`` (a ceiling on the factor) is what an exact Gillespie thinner needs when the
    driver is unbounded; omit it for a naturally-bounded ``fn``. The function must return a
    finite, non-negative number for every driver value it sees (a rate cannot go negative)."""

    def __init__(self, fn, bound: float | None = None) -> None:
        if not callable(fn):
            raise TypeError(f"Curve needs a callable value→factor function, got {fn!r}")
        if bound is not None:
            if isinstance(bound, bool) or not isinstance(bound, (int, float)) \
                    or not math.isfinite(bound) or bound < 0:
                raise ValueError(f"Curve bound must be a finite non-negative number, got {bound!r}")
            bound = float(bound)
        self.fn = fn
        self.bound = bound

    def multiplier(self, value: object) -> float:
        f = self.fn(_numeric(value, "Curve"))
        if isinstance(f, bool) or not isinstance(f, (int, float)) or not math.isfinite(f) or f < 0:
            raise ValueError(
                f"the Curve returned {f!r} for driver value {value!r}; a rate multiplier must be a "
                f"finite non-negative number"
            )
        f = float(f)
        return min(f, self.bound) if self.bound is not None else f

    def __repr__(self) -> str:
        tail = "" if self.bound is None else f", bound={self.bound:g}"
        return f"Curve({getattr(self.fn, '__name__', 'fn')}{tail})"


class Scalar(Mapping):
    """A single log-link coefficient — ``multiplier = exp(strength · value)``::

        Scalar(0.0)    # null: factor 1 for every value
        Scalar(0.7)    # a binary 0/1 driver gives factor 1 (off) or exp(0.7) ≈ 2.0 (on)

    The natural response when the driver is already a 0/1 indicator or a single continuous
    covariate: one knob, ``strength`` (0 ⇒ no coupling). The exponent is clamped so a large
    value cannot overflow."""

    def __init__(self, strength: float) -> None:
        if isinstance(strength, bool) or not isinstance(strength, (int, float)) \
                or not math.isfinite(strength):
            raise ValueError(f"Scalar strength must be a finite number, got {strength!r}")
        self.strength = float(strength)

    def multiplier(self, value: object) -> float:
        x = self.strength * _numeric(value, "Scalar")
        x = max(-_MAX_EXPONENT, min(_MAX_EXPONENT, x))  # guard exp() against overflow
        return math.exp(x)

    def __repr__(self) -> str:
        return f"Scalar(strength={self.strength:g})"

    def __eq__(self, other: object) -> bool:
        return isinstance(other, Scalar) and other.strength == self.strength


def _check_factor(x: object, where: str) -> float:
    """Coerce ``x`` to a finite, non-negative float (a rate multiplier) or raise naming ``where``."""
    if isinstance(x, bool) or not isinstance(x, (int, float)) or not math.isfinite(x) or x < 0:
        raise ValueError(f"{where} must be a finite non-negative number, got {x!r}")
    return float(x)


def _numeric(value: object, cls: str) -> float:
    """A driver value as a float, for a **continuous** response (``Curve`` / ``Scalar``). Raises a
    clear error if the driver is a discrete label — the usual sign a discrete driver was given a
    continuous mapping (use a ``Table`` / dict for a discrete driver)."""
    try:
        return float(value)
    except (TypeError, ValueError):
        raise ValueError(
            f"{cls} is a continuous-driver response but got the discrete driver value {value!r}; use a "
            f"Table (a dict mapping) for a discrete driver such as a habitat state."
        ) from None


def as_mapping(spec: object) -> Mapping:
    """Coerce a ``DrivenBy`` mapping spec into a :class:`Mapping`.

    Accepts an already-built mapping (returned unchanged), a ``dict`` (→ :class:`Table`), a
    callable (→ :class:`Curve`), or a number (→ :class:`Scalar`). Mirrors
    :func:`~zombi2.rates.rate.as_rate` / :func:`~zombi2.rates.distributions.as_distribution`.
    """
    if isinstance(spec, Mapping):
        return spec
    if isinstance(spec, dict):
        return Table(spec)
    if isinstance(spec, bool):
        raise TypeError(f"a DrivenBy mapping cannot be a bool, got {spec!r}")
    if isinstance(spec, (int, float)):
        return Scalar(float(spec))
    if callable(spec):
        return Curve(spec)
    raise TypeError(
        f"a DrivenBy mapping must be a dict (Table), a callable (Curve), a number (Scalar), or a "
        f"Table/Curve/Scalar, got {spec!r}"
    )


__all__ = ["Mapping", "Table", "Curve", "Scalar", "as_mapping"]
