"""Verbs — the *what the number does* half of a driven parameter (SPEC §5).

A parameter reads a **value** (`values`) through a **mapping** (`mapping`), and a
verb says what the resulting number does to the parameter. There are three, and which one you use is
decided by what you are attaching to rather than by taste::

    loss        = 0.25 * ScaledBy(habitat, {"cave": 4.0, "surface": 1.0})   # multiply the base
    loss        =        SetBy(habitat,    {"cave": 1.0, "surface": 0.25})  # replace the base
    transfer_to =        Weights(competence, {"competent": 3.0})            # compare candidates

`ScaledBy` multiplies, so its number is a dimensionless factor and the parameter keeps its base.
`SetBy` replaces, so its number carries the parameter's own units and there is no base to write.
`Weights` is for a **choice** — an argument that decides *who*, not how fast — which has no base at
all, because only the ratios between candidates are read.

**These build the objects the engines already run.** ``ScaledBy(Time(), {...})`` *is*
``OnTime({...})``; ``ScaledBy(habitat, {...})`` *is* ``DrivenBy(habitat, {...})``. The verb chooses
which by looking at the value, so nothing downstream changes and a run is identical either way.
"""

from __future__ import annotations

from .modifiers import DrivenBy, Modifier, OnTime, SetBy, describe
from .values import Measured, Time


def _mapping_for_time(mapping: object) -> dict:
    """A time value's mapping. A schedule — ``{0: 1.0, 3: 0.3}``, factor from each breakpoint on — is
    the only shape an engine reads today, because a rate that is piecewise-constant in time can be
    stepped to exactly. A smooth curve of time is a real model and is refused rather than
    approximated: it makes the rate vary continuously, which needs the engine to integrate its hazard
    rather than sample it."""
    if isinstance(mapping, dict):
        return mapping
    raise ValueError(
        "ScaledBy(Time(), ...) takes a schedule — {0: 1.0, 3: 0.3}, the factor from each breakpoint "
        "on — because a rate that changes in steps can be stepped to exactly. A smooth function of "
        "time is not implemented: it makes the rate vary continuously between events, which needs "
        "the engine to integrate the rate rather than read it at a point.")


def ScaledBy(value: object, mapping: object = None) -> Modifier:
    """Multiply the parameter's base by a factor read from ``value``.

    The factor is dimensionless, and almost every parameter takes one::

        loss  = 0.25 * ScaledBy(habitat, {"cave": 4.0, "surface": 1.0})   # a grown trait
        birth = 1.0  * ScaledBy(Time(), {0: 1.0, 3: 0.3})                 # the run's clock

    ``mapping`` turns the value into that factor, and its shape follows the value's **type**: a
    categorical value takes a table (a dict), a numerical one takes a curve (a callable) or a
    ``Scalar`` log-link. A `Drawn` or `Inherited`
    value is already a factor, so it is written on its own without a verb.
    """
    if isinstance(value, Time):
        return OnTime(_mapping_for_time(mapping))
    if isinstance(value, Measured):
        raise ValueError(
            f"ScaledBy({type(value).__name__}(), ...) is not implemented — that value exists in the "
            f"grammar but no engine supplies it yet.")
    if isinstance(value, Modifier):
        name = describe(value)
        raise TypeError(
            f"{name} is already a factor, so it needs no verb: write `rate = base * {name}(...)` "
            f"rather than wrapping it in ScaledBy. Verbs are for values a mapping has to turn into a "
            f"number — a trait, a clock, a count.")
    if mapping is None:
        raise ValueError(
            "ScaledBy(value, mapping) needs a mapping: a dict for a categorical value, a callable "
            "for a numerical one.")
    return DrivenBy(value, mapping)


def Weights(value: object, mapping: object = None) -> Modifier:
    """Weight the candidates of a **choice** — an argument that decides *who*, not how fast.

    ``transfer_to``, the recipient of a horizontal transfer, is the only choice today. A choice has no
    base, because only the ratios between candidates are read, and a weight of zero means that
    candidate cannot be chosen::

        transfer_to = Weights(competence, {"competent": 3.0, "normal": 1.0})

    A weight may read **both ends** — the donor's group and the recipient's — through a ``Between``
    kernel, which is the mapping for a value that sits on a pair rather than on one lineage.
    """
    if mapping is None:
        raise ValueError(
            "Weights(value, mapping) needs a mapping: a dict of per-candidate weights, a callable, or "
            "a Between kernel to weight the (donor, recipient) pair.")
    return DrivenBy(value, mapping)


#: The verbs a rate may be **written** with — see `zombi2.rates.modifiers.WRITABLE`.
WRITABLE = ("ScaledBy", "Weights")

__all__ = ["ScaledBy", "SetBy", "Weights", "WRITABLE"]
