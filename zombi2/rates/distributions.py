"""The shapes a run draws its numbers from: a `Random` value's law, and a segmental event's extent.

    from zombi2.rates import Gamma, PerCopy

    loss        = PerCopy(0.25).varying_among("families", Gamma(shape=4.0, scale=0.25))
    loss_extent = Gamma(shape=2.0, scale=250.0)

A handful of distributions ship here, and `as_distribution()` also takes any scipy.stats frozen
distribution (anything with an ``.rvs``) or a plain callable ``rng -> float``. There is no hard
scipy dependency.

Where the two uses differ is the mean. An **extent** is an absolute size, so its distribution is
used as written. A **drawn multiplier** is normalised to mean 1, so only its shape counts — which
is why one that cannot state its ``mean()`` is refused there rather than normalised by a guess.
"""

from __future__ import annotations

import math
from abc import ABC, abstractmethod

#: The distributions a rate or an extent may be **written** with — every one is built from literals,
#: so it round-trips through a run's log and a ``--params`` file. ``Distribution`` itself is absent:
#: it is abstract, and not a thing a user writes.
WRITABLE = ("Fixed", "Exponential", "Gamma", "LogNormal", "Uniform", "Geometric")

__all__ = [
    "Distribution", "Fixed", "Exponential", "Gamma", "LogNormal", "Uniform", "Geometric",
    "as_distribution", "WRITABLE",
]


class Distribution(ABC):
    """Something that yields a float given a numpy Generator."""

    #: What to call this in a message. A user writes ``scipy.stats.expon(...)`` or a lambda and
    #: never sees the wrapper's class name, so the wrappers below override it; everything else is
    #: named by its own class, which is exactly what was written.
    written: str = ""

    @abstractmethod
    def sample(self, rng) -> float:
        ...

    def __init_subclass__(cls, **kw) -> None:
        super().__init_subclass__(**kw)
        if not cls.__dict__.get("written"):
            cls.written = cls.__name__

    def mean(self) -> float:
        """This distribution's mean, where it is known in closed form.

        A `Drawn` value divides by it, because a drawn value is a
        *multiplier* and a multiplier whose average is not 1 changes what the base means. A
        distribution whose mean cannot be computed — a bare callable, a scipy frozen distribution —
        raises here rather than being normalised by a guess."""
        raise NotImplementedError(
            f"a {self.written} does not state its own mean, so it cannot be normalised to a "
            f"multiplier. Use one of the built-in distributions, or scale it yourself and read it "
            f"with SetBy, where the number is the value rather than a factor.")


class Fixed(Distribution):
    """A degenerate distribution — always the same value (also what a bare float becomes)."""

    def __init__(self, value: float):
        self.value = float(value)

    def sample(self, rng) -> float:
        return self.value

    def mean(self) -> float:
        return self.value

    def __repr__(self) -> str:
        return f'Fixed({self.value!r})'

    def __eq__(self, other: object) -> bool:
        return type(other) is type(self) and other.__dict__ == self.__dict__

    def __hash__(self) -> int:
        return hash((type(self), tuple(sorted(self.__dict__.items()))))


class Exponential(Distribution):
    """Exponential with the given mean."""

    def __init__(self, mean: float):
        if mean <= 0:
            raise ValueError(f"Exponential mean must be > 0, got {mean}")
        self.mean_ = float(mean)

    def sample(self, rng) -> float:
        return float(rng.exponential(self.mean_))

    def mean(self) -> float:
        return self.mean_

    def __repr__(self) -> str:
        return f'Exponential({self.mean_!r})'

    def __eq__(self, other: object) -> bool:
        return type(other) is type(self) and other.__dict__ == self.__dict__

    def __hash__(self) -> int:
        return hash((type(self), tuple(sorted(self.__dict__.items()))))


class Gamma(Distribution):
    """Gamma with shape ``k`` and scale ``theta`` (mean = k*theta)."""

    def __init__(self, shape: float, scale: float):
        if shape <= 0 or scale <= 0:
            raise ValueError("Gamma shape and scale must be > 0")
        self.shape = float(shape)
        self.scale = float(scale)

    def sample(self, rng) -> float:
        return float(rng.gamma(self.shape, self.scale))

    def mean(self) -> float:
        return self.shape * self.scale

    def __repr__(self) -> str:
        return f'Gamma(shape={self.shape!r}, scale={self.scale!r})'

    def __eq__(self, other: object) -> bool:
        return type(other) is type(self) and other.__dict__ == self.__dict__

    def __hash__(self) -> int:
        return hash((type(self), tuple(sorted(self.__dict__.items()))))


class LogNormal(Distribution):
    """Log-normal parameterised by the underlying normal's ``mu`` and ``sigma``."""

    def __init__(self, mu: float, sigma: float):
        if sigma < 0:
            raise ValueError("LogNormal sigma must be >= 0")
        self.mu = float(mu)
        self.sigma = float(sigma)

    def sample(self, rng) -> float:
        return float(rng.lognormal(self.mu, self.sigma))

    def mean(self) -> float:
        return math.exp(self.mu + 0.5 * self.sigma * self.sigma)

    def __repr__(self) -> str:
        return f'LogNormal({self.mu!r}, {self.sigma!r})'

    def __eq__(self, other: object) -> bool:
        return type(other) is type(self) and other.__dict__ == self.__dict__

    def __hash__(self) -> int:
        return hash((type(self), tuple(sorted(self.__dict__.items()))))


class Uniform(Distribution):
    """Uniform on [low, high]."""

    def __init__(self, low: float, high: float):
        if high < low:
            raise ValueError("Uniform requires high >= low")
        self.low = float(low)
        self.high = float(high)

    def sample(self, rng) -> float:
        return float(rng.uniform(self.low, self.high))

    def mean(self) -> float:
        return 0.5 * (self.low + self.high)

    def __repr__(self) -> str:
        return f'Uniform({self.low!r}, {self.high!r})'

    def __eq__(self, other: object) -> bool:
        return type(other) is type(self) and other.__dict__ == self.__dict__

    def __hash__(self) -> int:
        return hash((type(self), tuple(sorted(self.__dict__.items()))))


class Geometric(Distribution):
    """Geometric on ``{1, 2, 3, …}`` with the given ``mean`` (≥ 1) — a positive integer count, e.g. a
    segment's extent in genes. ``Geometric(mean=1)`` is degenerate at 1 (single-gene events)."""

    def __init__(self, mean: float):
        if mean < 1:
            raise ValueError(f"Geometric mean must be >= 1, got {mean}")
        self.mean_ = float(mean)

    def sample(self, rng) -> float:
        return float(rng.geometric(1.0 / self.mean_))

    def mean(self) -> float:
        return self.mean_

    def __repr__(self) -> str:
        return f'Geometric({self.mean_!r})'

    def __eq__(self, other: object) -> bool:
        return type(other) is type(self) and other.__dict__ == self.__dict__

    def __hash__(self) -> int:
        return hash((type(self), tuple(sorted(self.__dict__.items()))))


class _ScipyDist(Distribution):
    written = "scipy frozen distribution"

    def __init__(self, frozen):
        self._frozen = frozen

    def sample(self, rng) -> float:
        return float(self._frozen.rvs(random_state=rng))


class _CallableDist(Distribution):
    written = "callable"

    def __init__(self, fn):
        self._fn = fn

    def sample(self, rng) -> float:
        return float(self._fn(rng))


def as_distribution(spec) -> Distribution:
    """Coerce ``spec`` (Distribution | float | scipy frozen dist | callable) to a Distribution."""
    if isinstance(spec, Distribution):
        return spec
    if isinstance(spec, (int, float)):
        return Fixed(spec)
    if hasattr(spec, "rvs"):  # scipy.stats frozen distribution
        return _ScipyDist(spec)
    if callable(spec):
        return _CallableDist(spec)
    raise TypeError(f"cannot interpret {spec!r} as a distribution")

