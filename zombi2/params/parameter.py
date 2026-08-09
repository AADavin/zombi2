"""The **parameters** — the three things you set, and the only things a verb attaches to.

    rate    how often an event fires  — written from its scope: ``PerCopy(0.25)``
    extent  how much it takes         — ``Extent(Gamma(2.0, 250.0))``: no scope, already absolute
    choice  which candidate receives  — ``Recipients()``: no base, only ratios are read

Each opens with its own constructor and takes verbs from there. What differs between them is which
verbs are legal, and that is a fact about the **parameter** rather than about the driver: a rate can
be scaled or replaced, an extent only scaled, a choice only weighted. `Choice` is in `choice`, with
the two topological rules that are written the same way.
"""

from __future__ import annotations

from . import connection as verbs
from .distributions import Distribution, Geometric, as_distribution
from .evaluate import CARRIED_KINDS, Modifier
from .scope import Scope
from dataclasses import dataclass
from typing import Any
import math







class RateCompositionError(TypeError):
    """A composition the grammar itself refuses, with a message written for that mistake.

    A `TypeError` so nothing that catches one stops catching it, and its own class so the written
    form can tell it from CPython's ``unsupported operand type(s)``: the parser re-raises ours
    verbatim, because the sentence says what to write instead, and answers CPython's generically,
    because "Rate and float" is about types rather than about the rate.
    """


@dataclass(frozen=True, repr=False)     # repr=False: `__repr__` below is the written form
class Rate:
    """``base × scope × modifiers``, not yet evaluated.

    ``base`` is ``None`` for a rate whose number comes from a driver — ``PerCopy().set_by(...)`` —
    where there is no base to write and the scope still stands.
    """

    base: float | None = None
    scope: type[Scope] | None = None
    modifiers: tuple[Modifier, ...] = ()

    def __post_init__(self) -> None:
        """The base is checked here rather than at each scope, because every rate passes through
        this one constructor and a rule enforced per call site is a rule some level forgets."""
        if self.base is None:
            return
        if isinstance(self.base, bool) or not isinstance(self.base, (int, float)):
            raise TypeError(f"a rate base must be a real number, got {self.base!r}")
        if not math.isfinite(self.base) or self.base < 0:
            raise ValueError(f"a rate base must be finite and non-negative, got {self.base!r}")
        object.__setattr__(self, "base", float(self.base))

    # --- the verbs (SPEC §5): each returns a NEW rate, so they chain and nothing is mutated -------

    def scaled_by(self, driver: object, mapping: object = None, *,
                  step: float | None = None) -> "Rate":
        """Multiply this rate by a factor read from ``driver`` — see `verbs.scaled_by`."""
        return self._and(verbs.scaled_by(driver, mapping, step=step))

    def set_by(self, driver: object, mapping: object = None, *,
               step: float | None = None) -> "Rate":
        """Replace this rate's base with a number read from ``driver`` — see `verbs.set_by`.

        Written first, on a scope with no number in front of it, because everything to its left is
        a base it would silently discard. That is the one rule, stated here rather than at the two
        operand positions the retired ``*`` needed it at.
        """
        if self.base is not None or self.modifiers:
            raise RateCompositionError(
                f"set_by replaces the base, so it cannot follow one: everything to its left — the "
                f"number, any factors — is a base it would silently discard. Write it first, on the "
                f"bare scope: {self._scope_name()}().set_by(driver, mapping).scaled_by(...). "
                f"Got {self!r}.")
        m = verbs.set_by(driver, mapping, step=step)
        # `set_by(Time(), ...)` builds an `OnTime`, whose schedule holds the rates themselves; that
        # is a base of 1.0 times those factors, which is the same run and needs no new machinery.
        base = None if getattr(m, "replaces_base", False) else 1.0
        return Rate(base, self.scope, (m,))

    def varying_among(self, among: object = None, law: object = None, **retired: object) -> "Rate":
        """Let this rate vary at random among the units of one kind — see `verbs.varying_among`.

        ``**retired`` is passed straight through, so the keywords this verb replaced (``per=``,
        ``spread=``) are answered by the one table rather than by "unexpected keyword argument"."""
        return self._and(verbs.varying_among(among, law, **retired))

    def changing_at(self, schedule: object) -> "Rate":
        """Let this rate change in time, on a schedule of factors — see `verbs.changing_at`."""
        return self._and(verbs.changing_at(schedule))

    def weighted_by(self, driver: object, mapping: object = None, *,
                    step: float | None = None) -> "Rate":
        """Refused. Weights are compared against each other and normalised across candidates, which
        only a **choice** does — ``transfer_to`` is the only one."""
        raise RateCompositionError(
            "weighted_by weights the candidates of a choice against each other — transfer_to is "
            "the only one, written from Recipients(). On a rate the number multiplies a base, so "
            "the verb is scaled_by: the same driver and the same mapping, read as a factor.")

    def _and(self, m: Modifier) -> "Rate":
        return Rate(self.base, self.scope, self.modifiers + (m,))

    def _scope_name(self) -> str:
        return self.scope.__name__ if self.scope is not None else "PerLineage"

    # --- what the engines call --------------------------------------------------------------

    def with_default_scope(self, default: type[Scope]) -> "Rate":
        """Fill in the level's default scope (per lineage, per copy, …) when none was written."""
        if self.scope is not None:
            return self
        return Rate(self.base, default, self.modifiers)

    def effective(self, *, carried_factor: float = 1.0, **context: Any) -> float:
        """The rate *right now*: the scope-applied base times the product of the modifier factors.

        ``context`` carries the current state (``time``, ``diversity``, the counts ``lineages`` /
        ``copies`` / …); the scope reads the count it needs and each modifier the keys it needs.
        Requires a scope — resolve a bare-number rate with `with_default_scope()` first.

        ``carried_factor`` is the product of the values the engine drew and kept for the unit being
        evaluated — among lineages, among families (`carried_modifiers`). Those modifiers
        are skipped in the loop below, because their number does not come from the context: the
        engine already holds it and hands it in here, multiplied out. One float rather than a value
        per modifier, so a rate carrying several costs no more to evaluate than one carrying one.
        """
        if self.scope is None:
            raise ValueError("this rate has no scope yet; resolve it with with_default_scope(...)")
        base = self.base
        for m in self.modifiers:
            if getattr(m, "replaces_base", False):
                base = m.factor(**context)   # the driver supplies the number itself, not a factor
        if base is None:
            raise ValueError(
                "this rate has no number: a scope written on its own is only half a rate, and the "
                "set_by that would supply the rest is missing.")
        value = self.scope.total_of(base, **context)
        for m in self.modifiers:
            if getattr(m, "replaces_base", False):
                continue  # already used, as the base
            reads = getattr(m, "reads", None)
            if reads is not None and reads[0] in CARRIED_KINDS:
                continue  # its factor arrives through `carried_factor`, drawn and kept by the engine
            value *= m.factor(**context)
        return value * carried_factor

    def check_one_base(self, label: str = "this rate") -> None:
        """A rate has **one** base: a number written in front of the scope, or one `SetBy` supplying
        it. Two `SetBy`s would each claim to *be* the base, and no order of application is more
        right than another, so this raises rather than letting the last one written win in silence;
        none at all leaves a scope that says per what but not how fast.

        Every level coerces through `as_rate`, which calls this, so the rule cannot be strict in one
        place and lax in another."""
        set_by = [m for m in self.modifiers if getattr(m, "replaces_base", False)]
        if len(set_by) > 1:
            raise ValueError(
                f"{label} carries {len(set_by)} set_by verbs, and a base can only be replaced "
                f"once — each of them claims to be the whole number. Keep one; if you meant to scale "
                f"the result, that is scaled_by, which multiplies and composes freely.")
        if self.base is None and not set_by:
            raise ValueError(
                f"{label} is a scope with no number: {self._scope_name()}() says per what but not "
                f"how fast. Write the number — {self._scope_name()}(0.25) — or the driver that "
                f"supplies it, {self._scope_name()}().set_by(driver, mapping).")
        if any(verbs.written_with(m, verbs.WEIGHTED_BY) for m in self.modifiers):
            raise ValueError(
                f"{label} carries weighted_by, which weights the candidates of a choice against "
                f"each other — transfer_to is the only one. On a rate the number multiplies a base, "
                f"so the verb is scaled_by: the same driver and the same mapping, read as a factor.")

    def carried_modifiers(self, unit: str | None = None) -> tuple[tuple[Modifier, str], ...]:
        """Every modifier on this rate that reads a value the **engine** has to draw and carry,
        paired with the unit it is carried among (see `Modifier.reads`).

        A modifier reading a *measured* value computes its own factor from the context and needs
        nothing from the engine. A *drawn* or *inherited* one does: its number is produced once
        when a unit is born, kept for that unit's life, and handed back at every evaluation, and
        only the engine can do that. This is the one query for finding them, so a level does not
        have to know which modifier classes exist to thread them.

        ``unit`` narrows the answer to one kind of unit (``"lineages"``, ``"families"``). The result
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

    # --- the written form ---------------------------------------------------------------------

    def __repr__(self) -> str:
        """The expression that produced this rate. It is the ``__repr__`` and not a separate
        renderer so that the two cannot drift apart: a log, an error and a ``--params`` record all
        say the same thing, and it is the thing that was typed.

        A `SetBy` is written first, because anything to its left is a base it would discard, and the
        head then carries the scope with no number — ``PerCopy().set_by(...)`` — which is exactly
        what a replaced base means.
        """
        mods = self._ordered()
        bare = bool(mods) and verbs.written_with(mods[0], verbs.SET_BY)
        if self.scope is None:
            # A rate built directly rather than written: every written rate starts from its scope.
            # With a number there is still something true to print; with neither, a placeholder that
            # fails loudly if pasted back beats a plausible default that would name a different rate.
            head = repr(float(self.base)) if self.base is not None else "<no scope, no base>"
        elif bare or self.base is None:
            head = f"{self.scope.__name__}()"
        else:
            head = f"{self.scope.__name__}({float(self.base)!r})"
        return head + "".join(f".{m.written_call()}" for m in mods)

    def _ordered(self) -> tuple[Modifier, ...]:
        """The modifiers as written, except that a `SetBy` moves to the front (SPEC §5). ``sorted``
        is stable, so everything else keeps the order it was written in."""
        return tuple(sorted(self.modifiers, key=lambda m: not getattr(m, "replaces_base", False)))

    # --- the retired spelling -------------------------------------------------------------------

    def __mul__(self, other: object):
        from .evaluate import _no_longer_multiplied
        raise _no_longer_multiplied(self, other)

    def __rmul__(self, other: object):
        from .evaluate import _no_longer_multiplied
        raise _no_longer_multiplied(other, self)


def as_rate(spec: object, *, default_scope: type[Scope], label: str = "this rate") -> Rate:
    """Coerce a user rate spec into a resolved `Rate`, filling the level's default scope.

    Two cases, and there is no third: a bare number, which gets the level's default scope, and an
    already-built ``Rate`` — because a scope constructor returns one of those and so does every
    verb.

    Every level coerces its rates through here, which is why the one-base rule is checked here
    rather than in each level's own validation: a rule enforced by whoever remembers to call it is a
    rule three levels did not have.
    """
    if isinstance(spec, Rate):
        spec.check_one_base(label)
        return spec.with_default_scope(default_scope)
    if isinstance(spec, bool) or not isinstance(spec, (int, float)):
        raise TypeError(
            f"a rate is a number, or a scope with verbs chained onto it — {default_scope.__name__}"
            f"(0.25), PerLineage(0.5).changing_at({{0: 1.0, 3: 0.3}}) — got {spec!r}")
    return Rate(float(spec)).with_default_scope(default_scope)


__all__ = ["Rate", "RateCompositionError", "as_rate", "Extent", "as_extent"]


@dataclass(frozen=True, repr=False)     # repr=False: `__repr__` below is the written form
class Extent:
    """``base × modifiers``, not yet drawn. ``Extent(...)`` is the entry point a verb chains onto;
    a bare number or distribution is coerced by `as_extent()`."""

    # `Distribution` after `__post_init__` coerces it; a caller may hand in anything `_base` accepts
    base: Distribution = None    # type: ignore[assignment]
    modifiers: tuple[Modifier, ...] = ()

    def __post_init__(self) -> None:
        object.__setattr__(self, "base", _base(self.base))

    # --- the verbs (SPEC §5) ---------------------------------------------------------------------

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
        from .evaluate import _no_longer_multiplied
        raise _no_longer_multiplied(self, other)

    def __rmul__(self, other: object):
        from .evaluate import _no_longer_multiplied
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
    from .parameter import Rate

    if isinstance(spec, Extent):
        return spec
    if getattr(spec, "replaces_base", False):
        raise _cannot_be_set_by()
    if isinstance(spec, Rate):
        if any(getattr(m, "replaces_base", False) for m in spec.modifiers):
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
