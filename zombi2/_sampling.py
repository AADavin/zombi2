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


class Fenwick:
    """A Fenwick (binary indexed) tree for weighted sampling over a fixed index set.

    Supports O(log n) point-assignment (``set``) and O(log n) sampling (``find`` — the
    smallest index whose prefix sum reaches a value), with the running ``total`` kept in
    O(1). Used by the simulator to pick the next event's branch in O(log branches) instead
    of scanning every branch; slots for non-alive branches simply hold 0.
    """

    def __init__(self, n: int):
        self.n = n
        self.tree = [0.0] * (n + 1)  # 1-indexed internal storage
        self.vals = [0.0] * n        # current value at each 0-indexed slot
        self.total = 0.0

    def set(self, i: int, value: float) -> None:
        delta = value - self.vals[i]
        if delta == 0.0:
            return
        self.vals[i] = value
        self.total += delta
        j = i + 1
        while j <= self.n:
            self.tree[j] += delta
            j += j & (-j)

    def find(self, value: float) -> int:
        """Smallest 0-indexed slot whose inclusive prefix sum is >= ``value`` (0 < value <= total)."""
        pos = 0
        for k in range(self.n.bit_length(), -1, -1):
            nxt = pos + (1 << k)
            if nxt <= self.n and self.tree[nxt] < value:
                pos = nxt
                value -= self.tree[nxt]
        return pos  # 0-indexed position of the found element
