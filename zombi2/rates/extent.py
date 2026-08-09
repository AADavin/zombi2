"""The extent — how much a segmental event takes once it has started (SPEC §6).

An extent is the second axis of a segmental event. The **rate** says how often one starts; the extent
says how much it takes, and the two multiply. It is written the way a rate is, **minus the scope**,
because it is already an absolute quantity and has no "per what?" to answer::

    loss_extent = Gamma(2.0, 250.0)
    loss_extent = Extent(Gamma(2.0, 250.0)).scaled_by(habitat, {"host": 3.0, "free": 1.0})

The base is a size — a number (the mean) or a `Distribution` — and the verbs read the same drivers a
rate's do, scaling the size. The second line reads as "host-restricted lineages delete stretches
three times as long", which is a different statement from raising the loss *rate*: one deletes more
often, the other deletes in bigger chunks.

``Extent(...)`` is written only when something is chained onto it. A bare size needs no wrapper,
which is why the plain case above is the distribution alone.

**An extent carries no Gillespie breakpoints.** It is read when an event fires, not while the race
runs, so a modifier on an extent never changes a total rate and never has to be stepped to — unlike
the same modifier on a rate. That is why nothing here contributes to a horizon.
"""

from __future__ import annotations

from dataclasses import dataclass

from . import verbs
from .distributions import Distribution, Geometric, as_distribution
from .modifiers import Modifier

__all__ = ["Extent", "as_extent"]


@dataclass(frozen=True, repr=False)     # repr=False: `__repr__` below is the written form
class Extent:
    """``base × modifiers``, not yet drawn. ``Extent(...)`` is the entry point a verb chains onto;
    a bare number or distribution is coerced by `as_extent()`."""

    # `Distribution` after `__post_init__` coerces it; a caller may hand in anything `_base` accepts
    base: Distribution = None    # type: ignore[assignment]
    modifiers: tuple[Modifier, ...] = ()

    def __post_init__(self) -> None:
        object.__setattr__(self, "base", _base(self.base))

    # --- the verbs (SPEC §4) ---------------------------------------------------------------------

    def scaled_by(self, driver: object, mapping: object = None, *,
                  step: float | None = None) -> "Extent":
        """Multiply the size by a factor read from ``driver`` — see `verbs.scaled_by`."""
        return self._and(verbs.scaled_by(driver, mapping, step=step))

    def changing_at(self, schedule: object) -> "Extent":
        """Let the size change in time, on a schedule of factors — see `verbs.changing_at`."""
        return self._and(verbs.changing_at(schedule))

    def varying_among(self, among: object = None, law: object = None,
                      **retired: object) -> "Extent":
        """Refused. An extent carries no drawn or inherited value.

        Such a value is not computed from the context: it is drawn by the engine when a unit is
        born, kept for that unit's life, and handed back at every reading — and no level does that
        for an extent. There is nowhere for it to arrive, either, since `sample` and `mean` read
        their modifiers through ``factor()``, which a carried one deliberately has none of. Every
        level's gate refuses one, so this only ever built an object whose every read path raised
        mid-run; the refusal belongs where it is written. It takes the retired keywords too, so
        ``per=`` reaches this sentence rather than a complaint about an unexpected argument — an
        extent that cannot vary at all is the more useful thing to be told."""
        raise ValueError(
            "an extent cannot vary at random among units: such a value is drawn by the engine when "
            "a unit is born and handed back at each reading, and no level carries one for an "
            "extent — an extent takes scaled_by and changing_at. Two things say what this usually "
            "means — vary the RATE, so events start more often in some units: "
            "PerCopy(0.25).varying_among('families', LogNormal(0.0, 0.5)); or scale the extent by "
            "a driver, so events take more when they do start: "
            "Extent(500).scaled_by(driver, {...}).")

    def set_by(self, driver: object, mapping: object = None, *,
               step: float | None = None) -> "Extent":
        """Refused. An extent's base is a *distribution* over sizes, so one scalar cannot replace
        it — the replacement would fix every event to the same size, which is not what any of the
        sizes here mean."""
        raise _cannot_be_set_by()

    def _and(self, m: Modifier) -> "Extent":
        return Extent(self.base, self.modifiers + (m,))

    # --- what the engines call --------------------------------------------------------------

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
        return self.base.mean() * self._factor(**context)

    @property
    def has_modifiers(self) -> bool:
        """Whether anything about this extent varies with context. ``False`` is the common case, and
        lets an engine skip building a context it would not read.

        Not called ``is_driven``: *driven* is one of the four modifier kinds (SPEC §5), and this is
        true of a schedule too, which is measured."""
        return bool(self.modifiers)

    # --- the written form ---------------------------------------------------------------------

    def __repr__(self) -> str:
        """The expression that produced this extent. A plain size renders as the distribution alone,
        because that is how it is written — the wrapper appears only when a verb is chained on."""
        if not self.modifiers:
            return repr(self.base)
        return f"Extent({self.base!r})" + "".join(f".{m.written_call()}" for m in self.modifiers)

    def __mul__(self, other: object):
        from .modifiers import _no_longer_multiplied
        raise _no_longer_multiplied(self, other)

    def __rmul__(self, other: object):
        from .modifiers import _no_longer_multiplied
        raise _no_longer_multiplied(other, self)


def as_extent(spec) -> Extent:
    """Coerce an extent spec (SPEC §6) — a number, a distribution, or an ``Extent`` with verbs.

    A bare number is the **mean**, not an exact size: ``3`` is ``Geometric(mean=3)``, so runs vary
    around three. Write ``Fixed(3)`` for exactly three every time. ``None`` is ``Geometric(mean=1)``,
    a single unit, the default wherever an extent is optional.

    This is where an extent parts company with `as_distribution()`,
    where a bare number is a *fixed* value. The readings differ because the quantities do: a sampled
    per-family rate given as ``0.1`` means that rate, whereas an extent given as ``500`` means runs of
    about 500 — nobody wants every inversion to be exactly the same size.

    A **rate** is refused. ``PerLineage(500)`` asks "per what?", and an extent has no answer: it is
    already an absolute size.
    """
    from .modifiers import SetBy
    from .rate import Rate

    if isinstance(spec, Extent):
        return spec
    if isinstance(spec, SetBy):
        raise _cannot_be_set_by()
    if isinstance(spec, Rate):
        if any(isinstance(m, SetBy) for m in spec.modifiers):
            raise _cannot_be_set_by()
        raise ValueError(
            "an extent takes no scope — it is already an absolute size, and there is no 'per what?' "
            "to answer (SPEC §6). Write the size alone, or Extent(500).scaled_by(...) when a verb "
            "is chained onto it.")
    return Extent(spec)


def _cannot_be_set_by() -> ValueError:
    """One sentence for the refusal, wherever a replaced extent is written — on the verb, or on a
    rate that carries one and was handed to an extent."""
    return ValueError(
        "an extent cannot be set_by. An extent is already an absolute size drawn from a "
        "distribution, so there is no base for a driver to replace — and a replaced one would fix "
        "every event to the same size, which is not what any of the sizes here mean. Scale the "
        "distribution instead: Extent(500).scaled_by(driver, {...}) (SPEC §6).")


def _base(spec) -> Distribution:
    """The size itself: a bare number is a geometric of that mean, ``None`` a single unit."""
    if spec is None:
        return Geometric(mean=1)
    if isinstance(spec, Distribution):
        return spec
    if isinstance(spec, (int, float)) and not isinstance(spec, bool):
        return Geometric(mean=spec)
    return as_distribution(spec)
