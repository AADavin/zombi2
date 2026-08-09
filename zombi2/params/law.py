"""The **laws** — what a value drawn per unit does from there.

``Random(unit, law)`` gives every unit of one kind its own multiplier. The law says what happens to
it afterwards, which is a separate question from what it starts as:

- a bare distribution — drawn once and held for that unit's life, with no memory of its parent;
- `Drift(dist)` — the parent's value times a draw from ``dist`` at each split, so the value wanders
  down the genealogy and close relatives keep similar ones.

`Drawn` and `Inherited` are what those two build. Neither is written by anyone: a law is, and the
law is what decides which of them a `Random` becomes. They are here rather than beside the base
class because reading a built object without the thing that writes it is what made the old
``modifiers`` hard.
"""

from __future__ import annotations

from dataclasses import dataclass

import math
from typing import Any

from .distributions import Distribution, LogNormal, as_distribution
from .evaluate import (DRAWN, INHERITED, UNITS, VARYING_AMONG, Modifier, _WRITTEN_AS)


@dataclass(frozen=True)
class Drift:
    """A **law** for `Random`: the parent's value times a mean-corrected draw at each split —
    continuous memory down a genealogy::

        birth = PerLineage(1.0).varying_among('lineages', Drift(LogNormal(0.0, 0.2)))  # ClaDS

    ``dist`` is the **per-split step**, not the value. Each law owns and documents its own
    argument — a bare distribution is the value itself, ``Drift``'s is the step — so nobody has to
    infer the role from the slot, which is what the retired ``spread=`` made impossible.

    ``bins`` discretises the drift onto a ladder: the rate takes one of ``bins`` values and a
    daughter moves to a **neighbouring rung**, or stays — the rate-category clock. A knob rather
    than a law of its own, because the model is unchanged: a daughter starts from its parent and is
    perturbed.
    """

    #: `Any` rather than `Distribution`: a law may be written as anything `as_distribution` accepts,
    #: and `Inherited` is what coerces it — checking here would only move the same check earlier and
    #: leave two places to keep in step.
    dist: Any = None
    bins: int | None = None

    def __repr__(self) -> str:
        extra = f", bins={self.bins}" if self.bins is not None else ""
        return f"Drift({self.dist!r}{extra})"


@dataclass(frozen=True)
class Inherited(Modifier):
    """What ``varying_among(unit, Drift(dist))`` builds — a value inherited and perturbed.

    Geometric Brownian motion on the rate: the step is **mean-corrected** so ``E[factor] = 1``.
    Without the correction the rate inflates down the tree (``E[rate] ≈ e^{σ²/2}``) — a real
    historical bug. Any distribution that states a mean works, and is corrected by dividing by it,
    exactly as `Drawn` does. The engine drives `initial` / `descend`, threading each unit's current
    factor back through `Rate.effective`'s ``carried_factor``.

    ``bins`` takes a `LogNormal` step and nothing else, because the ladder's spacing *is* that
    step's σ — with any other shape there is no rung spacing to read, and inventing one would be a
    model nobody wrote. ``None`` is the continuous form.
    """

    unit: str = ""
    # `Distribution` after `__post_init__` coerces it; a caller may hand in anything
    # `as_distribution` accepts, exactly as `Drawn` does
    dist: "Distribution" = None    # type: ignore[assignment]
    bins: int | None = None

    @property
    def reads(self) -> tuple[str, str]:      # type: ignore[override]
        return (INHERITED, self.unit)

    def __post_init__(self) -> None:
        if self.unit not in UNITS:
            raise ValueError(
                f"unknown unit {self.unit!r}; a value varies among one of {list(UNITS)}")
        if self.dist is None:
            raise ValueError(
                f"Drift needs the distribution of its per-split step — "
                f"varying_among({self.unit!r}, Drift(LogNormal(0.0, 0.2))).")
        object.__setattr__(self, "dist", as_distribution(self.dist))
        mean = self.dist.mean()       # raises for a distribution that cannot state one
        if not math.isfinite(mean) or mean <= 0:
            raise ValueError(
                f"an inherited step is mean-corrected by its distribution's mean, which must be "
                f"finite and positive; {self.dist!r} has a mean of {mean!r}")
        if self.bins is not None:
            if isinstance(self.bins, bool) or not isinstance(self.bins, int):
                raise TypeError(f"Drift bins must be a whole number, got {self.bins!r}")
            if self.bins < 2:
                raise ValueError(f"Drift bins must be at least 2, got {self.bins}")
            if not isinstance(self.dist, LogNormal):
                raise ValueError(
                    f"Drift(..., bins=) takes a LogNormal step and got {self.dist!r}. The ladder's "
                    f"rungs are spaced by that step's sigma, so another shape leaves no spacing to "
                    f"read — write Drift(LogNormal(0.0, sigma), bins=n), or drop bins for the "
                    f"continuous drift, which takes any distribution.")

    def written_call(self) -> str:
        return f"{VARYING_AMONG}({self.unit!r}, {Drift(self.dist, self.bins)!r})"

    def __repr__(self) -> str:
        """The standalone form — the `Random` a user can name and share between two rates. ``bins``
        is omitted when it is unset, because the parser refuses ``bins=None``: a written form that
        will not parse is not a written form."""
        return f"Random({self.unit!r}, {Drift(self.dist, self.bins)!r})"

    def __eq__(self, other: object) -> bool:
        return (isinstance(other, Inherited) and other.unit == self.unit
                and other.dist == self.dist and other.bins == self.bins)

    def __hash__(self) -> int:
        return hash((Inherited, self.unit, self.dist, self.bins))

    def factor(self, **_: Any) -> float:
        """A carried modifier has no factor to compute. Its number is drawn when the unit is born and
        kept by the engine, which hands it back through `Rate.effective`'s ``carried_factor`` — so asking
        for one here, without that stored value, could only ever return a plausible 1.0 for a rate
        that should have varied."""
        raise NotImplementedError(
            f"a Random is carried, not computed: the engine draws its value for each of the "
            f"{self.unit} and passes it to Rate.effective(carried_factor=...). It has no factor of "
            f"its own.")

    def _ladder(self) -> list[float]:
        """The ``bins`` rungs, geometric in the step's sigma and scaled so their mean is 1.

        The mean is taken over the rungs because a nearest-neighbour walk with reflecting ends is
        uniform at stationarity, so an even ladder gives ``E[factor] = 1`` — the same promise the
        continuous form keeps with its lognormal correction."""
        n = self.bins or 0
        assert isinstance(self.dist, LogNormal)   # __post_init__ refuses bins with anything else
        sigma = self.dist.sigma
        raw = [math.exp(sigma * (k - (n - 1) / 2)) for k in range(n)]
        mean = sum(raw) / n
        return [r / mean for r in raw]

    def initial(self) -> float:
        """The root's factor: 1.0 for the continuous form; the middle rung when binned."""
        if self.bins is None:
            return 1.0
        rungs = self._ladder()
        return rungs[len(rungs) // 2]

    def descend(self, parent_value: float, rng) -> float:
        """A daughter's factor: the parent's, times one mean-corrected step — or, when binned, its
        parent's rung or one either side, with the ends reflecting.

        The parent's rung is recovered by nearest match rather than by inverting the ladder, so no
        float drift down a deep tree can put a lineage between two rungs.

        A lognormal step keeps its closed form, ``exp(N(-σ²/2, σ))``, rather than going through
        ``sample()/mean()``. The two are the same distribution and consume the same one draw, but they
        differ in the last bit — and a one-ULP change to a step is enough to move a Gillespie
        comparison, so every autocorrelated run in the repository would have shifted for nothing.
        ``mu`` cancels in the correction, which is why only ``sigma`` appears."""
        if self.bins is None:
            if isinstance(self.dist, LogNormal):
                sigma = self.dist.sigma
                return parent_value * math.exp(rng.normal(-0.5 * sigma * sigma, sigma))
            return parent_value * (self.dist.sample(rng) / self.dist.mean())
        rungs = self._ladder()
        here = min(range(len(rungs)), key=lambda i: abs(math.log(rungs[i] / parent_value)))
        step = int(rng.integers(-1, 2))                      # -1, 0 or +1
        return rungs[min(max(here + step, 0), len(rungs) - 1)]


class Drawn(Modifier):
    """What ``varying_among(unit, dist)`` builds when the law is a **bare distribution** — a value
    drawn once for each unit and then fixed for that unit's life, i.i.d. heterogeneity with no
    memory. The unit is an argument, so one class covers every cell::

        loss         = PerCopy(0.25).varying_among('families', LogNormal(0.0, 0.5))
        substitution = PerSite(1.0).varying_among('lineages', LogNormal(0.0, 0.3))
        loss         = PerCopy(0.25).varying_among('families', Gamma(shape=4.0, scale=0.25))

    The law is the distribution the **value** is drawn from — any built-in `Distribution`:
    ``Fixed``, ``Exponential``, ``Gamma``, ``LogNormal``, ``Uniform``, ``Geometric``. A callable or a
    scipy frozen distribution is refused here, because neither states a mean to normalise by; an
    extent takes both, being a size rather than a multiplier.

    **Every draw is normalised to mean 1**, by dividing by the distribution's own mean. A drawn value
    is a *multiplier*, and one that does not average to 1 changes what the base means — a base of 0.25
    would stop being the average rate. So widening the spread separates units and leaves the average
    one where you put it, and a distribution's **location is normalised away**: what it contributes is
    its shape. ``Exponential(1.0)`` and ``Exponential(7.0)`` are therefore the same law. A
    distribution that cannot state its mean — a bare callable, a scipy frozen distribution — is
    refused rather than normalised by a guess. (A number that *is* the rate rather than a factor is
    `set_by`, where nothing is normalised.)

    **Writing one object or two decides what varies together.** One `Random` read by several rates is
    one draw for that unit, so a family that loses fast also duplicates fast; two separately built
    ones are independent even with identical arguments (SPEC §5).
    """

    def __init__(self, unit: str, dist: object = None) -> None:
        if unit not in UNITS:
            raise ValueError(f"unknown unit {unit!r}; a value varies among one of {list(UNITS)}")
        if dist is None:
            raise ValueError(
                f"a Random needs the law its value follows — "
                f"varying_among({unit!r}, LogNormal(0.0, 0.5)).")
        self.unit = unit
        self.dist = as_distribution(dist)
        mean = self.dist.mean()       # raises for a distribution that cannot state one
        if not math.isfinite(mean) or mean <= 0:
            raise ValueError(
                f"a drawn multiplier is normalised by its distribution's mean, which must be finite "
                f"and positive; {self.dist!r} has a mean of {mean!r}")

    @property
    def reads(self) -> tuple[str, str]:      # type: ignore[override]
        return (DRAWN, self.unit)

    def factor(self, **_: Any) -> float:
        """A carried modifier has no factor to compute — see `Inherited.factor`."""
        raise NotImplementedError(
            f"a Random is carried, not computed: the engine draws its value for each of the "
            f"{self.unit} and passes it to Rate.effective(carried_factor=...). It has no factor of "
            f"its own.")

    def draw(self, rng) -> float:
        """One independent multiplier for a unit, with mean 1.

        A degenerate distribution gives exactly 1.0 and consumes no randomness, so switching
        heterogeneity off leaves every other draw in the run where it was — which is what makes a
        zero-variation control comparable to the run it controls for."""
        if isinstance(self.dist, LogNormal) and self.dist.sigma == 0.0:
            return 1.0
        return self.dist.sample(rng) / self.dist.mean()

    def written_call(self) -> str:
        return f"{VARYING_AMONG}({self.unit!r}, {self.dist!r})"

    def __repr__(self) -> str:
        """The standalone form — the `Random` a user can name and share between two rates."""
        return f"Random({self.unit!r}, {self.dist!r})"

    def __eq__(self, other: object) -> bool:
        """Two draws are the same modifier when they are on the same unit with the same
        distribution."""
        return isinstance(other, Drawn) and other.unit == self.unit and other.dist == self.dist

    def __hash__(self) -> int:
        return hash((Drawn, self.unit, self.dist))


_WRITTEN_AS.update({Drawn: VARYING_AMONG, Inherited: VARYING_AMONG})

#: The laws a `Random` may be **written** with. A bare distribution is the other one and needs no
#: name — it is the law of holding what you drew.
WRITABLE = ("Drift",)

__all__ = ["Drift", "Drawn", "Inherited", "WRITABLE"]
