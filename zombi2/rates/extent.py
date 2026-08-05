"""The extent — how much a segmental event takes once it has started (SPEC §6).

An extent is the second axis of a segmental event. The **rate** says how often one starts; the extent
says how much it takes, and the two multiply. It is written the way a rate is, **minus the scope**,
because it is already an absolute quantity and has no "per what?" to answer::

    extent  =  base × modifiers

``base`` is a size — a number (the mean) or a `Distribution` — and
the modifiers are the same dimensionless multipliers a rate takes, scaling the size::

    loss_extent = 800 * mod.DrivenBy(habitat, {"host": 3.0, "free": 1.0})

reads as "host-restricted lineages delete stretches three times as long", which is a different
statement from raising the loss *rate*: one deletes more often, the other deletes in bigger chunks.

**An extent carries no Gillespie breakpoints.** It is read when an event fires, not while the race
runs, so a modifier on an extent never changes a total rate and never has to be stepped to — unlike
the same modifier on a rate. That is why nothing here contributes to a horizon.
"""

from __future__ import annotations

from dataclasses import dataclass

from .distributions import Distribution, Geometric, as_distribution
from .modifiers import Modifier

__all__ = ["Extent", "as_extent"]


@dataclass(frozen=True)
class Extent:
    """``base × modifiers``, not yet drawn. Internal — built by `as_extent()`, never by users."""

    base: Distribution
    modifiers: tuple[Modifier, ...] = ()

    def _factor(self, **context: object) -> float:
        f = 1.0
        for m in self.modifiers:
            f *= m.factor(**context)
        return f

    def sample(self, rng, **context: object) -> float:
        """One drawn size, scaled by the modifiers in this context.

        Scaling the **draw** rather than the distribution's parameter is what lets any base work —
        ``Fixed``, a scipy frozen distribution, a bare callable — while still meaning what it says:
        the expected size is the base's mean times the factor."""
        return float(self.base.sample(rng)) * self._factor(**context)

    def mean(self, **context: object) -> float:
        """The mean size in this context, for an engine parameterised by the mean rather than by a
        drawn value (the nucleotide one samples an arc's far end from a geometric of this mean).

        Requires a `Geometric` base, which is the only shape that
        engine supports; scaling its mean is the same statement in expectation as scaling a draw."""
        if not isinstance(self.base, Geometric):
            raise ValueError(
                f"this extent's base is {type(self.base).__name__}, which has no mean to scale — an "
                f"engine parameterised by the mean takes a geometric extent only.")
        return self.base.mean * self._factor(**context)

    @property
    def is_driven(self) -> bool:
        """Whether anything about this extent varies with context. ``False`` is the common case, and
        lets an engine skip building a context it would not read."""
        return bool(self.modifiers)


def as_extent(spec) -> Extent:
    """Coerce an extent spec (SPEC §6) — a number, a distribution, or ``base × modifiers``.

    A bare number is the **mean**, not an exact size: ``3`` is ``Geometric(mean=3)``, so runs vary
    around three. Write ``Fixed(3)`` for exactly three every time. ``None`` is ``Geometric(mean=1)``,
    a single unit, the default wherever an extent is optional.

    This is where an extent parts company with `as_distribution()`,
    where a bare number is a *fixed* value. The readings differ because the quantities do: a sampled
    per-family rate given as ``0.1`` means that rate, whereas an extent given as ``500`` means runs of
    about 500 — nobody wants every inversion to be exactly the same size.

    A **scope** is refused. ``PerLineage(500) * …`` asks "per what?", and an extent has no answer: it
    is already an absolute size.
    """
    from .rate import Rate                       # a `500 * modifier` product arrives as a Rate

    if isinstance(spec, Extent):
        return spec
    if isinstance(spec, Rate):
        if spec.scope is not None:
            raise ValueError(
                "an extent takes no scope — it is already an absolute size, and there is no "
                "'per what?' to answer (SPEC §6). Write the size and its modifiers alone, "
                "e.g. 500 * DrivenBy(...).")
        return Extent(_base(spec.base), tuple(spec.modifiers))
    return Extent(_base(spec))


def _base(spec) -> Distribution:
    """The size itself: a bare number is a geometric of that mean, ``None`` a single unit."""
    if spec is None:
        return Geometric(mean=1)
    if isinstance(spec, Distribution):
        return spec
    if isinstance(spec, (int, float)) and not isinstance(spec, bool):
        return Geometric(mean=spec)
    return as_distribution(spec)
