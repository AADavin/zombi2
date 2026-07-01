"""The numeric hot loop, isolated behind a narrow protocol.

Everything numerically hot in the forward simulation — the exponential waiting time and
the weighted choice of (branch, event) — goes through :class:`EventSampler`. It knows
nothing about genomes or event types (it takes a total rate and a weight vector), which
keeps it a clean boundary for a future Rust/PyO3 drop-in replacement.
"""

from __future__ import annotations

import math
from abc import ABC, abstractmethod

import numpy as np


class EventSampler(ABC):
    """Numeric sampling primitives for the Gillespie loop."""

    @abstractmethod
    def next_waiting_time(self, total_rate: float, rng: np.random.Generator) -> float:
        """Time to the next event ~ Exp(total_rate)."""

    @abstractmethod
    def choose_index(self, weights, rng: np.random.Generator) -> int:
        """Index chosen with probability proportional to ``weights`` (non-negative)."""


class NumpyEventSampler(EventSampler):
    """Default numpy implementation."""

    def next_waiting_time(self, total_rate: float, rng: np.random.Generator) -> float:
        if total_rate <= 0.0:
            return math.inf
        return float(rng.exponential(1.0 / total_rate))

    def choose_index(self, weights, rng: np.random.Generator) -> int:
        w = np.asarray(weights, dtype=float)
        cumulative = np.cumsum(w)
        total = cumulative[-1]
        if total <= 0.0:
            raise ValueError("choose_index called with non-positive total weight")
        r = rng.random() * total
        return int(np.searchsorted(cumulative, r, side="right"))
