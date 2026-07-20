"""Small distribution helpers for per-family sampled rates.

Ships a handful of built-in distributions and also accepts, via
:func:`as_distribution`, any scipy.stats frozen distribution (anything with an
``.rvs`` method) or a plain callable ``rng -> float``. No hard scipy dependency.

    z.FamilySampledRates(duplication=z.Gamma(2, 0.1),           # built-in
                         transfer=scipy.stats.expon(scale=0.1),  # scipy frozen dist
                         loss=lambda rng: rng.gamma(2, 0.05))    # callable
"""

from __future__ import annotations

from abc import ABC, abstractmethod

__all__ = [
    "Distribution", "Fixed", "Exponential", "Gamma", "LogNormal", "Uniform", "Geometric",
    "as_distribution",
]


class Distribution(ABC):
    """Something that yields a float given a numpy Generator."""

    @abstractmethod
    def sample(self, rng) -> float:
        ...


class Fixed(Distribution):
    """A degenerate distribution — always the same value (also what a bare float becomes)."""

    def __init__(self, value: float):
        self.value = float(value)

    def sample(self, rng) -> float:
        return self.value


class Exponential(Distribution):
    """Exponential with the given mean."""

    def __init__(self, mean: float):
        if mean <= 0:
            raise ValueError(f"Exponential mean must be > 0, got {mean}")
        self.mean = float(mean)

    def sample(self, rng) -> float:
        return float(rng.exponential(self.mean))


class Gamma(Distribution):
    """Gamma with shape ``k`` and scale ``theta`` (mean = k*theta)."""

    def __init__(self, shape: float, scale: float):
        if shape <= 0 or scale <= 0:
            raise ValueError("Gamma shape and scale must be > 0")
        self.shape = float(shape)
        self.scale = float(scale)

    def sample(self, rng) -> float:
        return float(rng.gamma(self.shape, self.scale))


class LogNormal(Distribution):
    """Log-normal parameterised by the underlying normal's ``mu`` and ``sigma``."""

    def __init__(self, mu: float, sigma: float):
        if sigma < 0:
            raise ValueError("LogNormal sigma must be >= 0")
        self.mu = float(mu)
        self.sigma = float(sigma)

    def sample(self, rng) -> float:
        return float(rng.lognormal(self.mu, self.sigma))


class Uniform(Distribution):
    """Uniform on [low, high]."""

    def __init__(self, low: float, high: float):
        if high < low:
            raise ValueError("Uniform requires high >= low")
        self.low = float(low)
        self.high = float(high)

    def sample(self, rng) -> float:
        return float(rng.uniform(self.low, self.high))


class Geometric(Distribution):
    """Geometric on ``{1, 2, 3, …}`` with the given ``mean`` (≥ 1) — a positive integer count, e.g. a
    segment/extension length in genes. ``Geometric(mean=1)`` is degenerate at 1 (single-gene events)."""

    def __init__(self, mean: float):
        if mean < 1:
            raise ValueError(f"Geometric mean must be >= 1, got {mean}")
        self.mean = float(mean)

    def sample(self, rng) -> float:
        return float(rng.geometric(1.0 / self.mean))


class _ScipyDist(Distribution):
    def __init__(self, frozen):
        self._frozen = frozen

    def sample(self, rng) -> float:
        return float(self._frozen.rvs(random_state=rng))


class _CallableDist(Distribution):
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
