"""Timing primitives for the ZOMBI2 performance suite.

The whole suite rests on one idea: a *benchmark* is a function you can call
repeatedly, timed with a monotonic clock, after an untimed warm-up. Everything
else (scaling curves, engine comparisons) is just calling :func:`measure` over a
grid of problem sizes and stashing the raw per-repeat times.

We keep every repeat (not just the mean) so the plotting layer can draw honest
error bands and so a re-analysis never needs to re-run the sims.
"""

from __future__ import annotations

import time
from dataclasses import dataclass, field, asdict
from statistics import median as _median, pstdev
from typing import Callable


# --- one measured point --------------------------------------------------

@dataclass
class Point:
    """One benchmarked configuration: a series label, an x value, and the raw times.

    ``times`` holds every timed repeat in seconds. ``work`` carries problem-size
    metadata (number of tips, families, events, ...) so figures can report
    throughput, and so a stale result is self-describing.
    """

    series: str                     # which curve this belongs to, e.g. "Rust · profiles"
    x: float                        # independent variable, e.g. number of tips
    times: list[float]              # every timed repeat, seconds
    work: dict = field(default_factory=dict)

    # summary statistics ---------------------------------------------------
    @property
    def best(self) -> float:
        """Fastest repeat — the cleanest estimate of intrinsic cost (least noise)."""
        return min(self.times)

    @property
    def median(self) -> float:
        return _median(self.times)

    @property
    def worst(self) -> float:
        return max(self.times)

    @property
    def spread(self) -> float:
        """Population standard deviation across repeats (0 for a single repeat)."""
        return pstdev(self.times) if len(self.times) > 1 else 0.0

    def to_dict(self) -> dict:
        return asdict(self)

    @classmethod
    def from_dict(cls, d: dict) -> "Point":
        return cls(series=d["series"], x=d["x"], times=list(d["times"]),
                   work=dict(d.get("work", {})))


# --- the measurement loop ------------------------------------------------

def measure(
    fn: Callable[[], object],
    *,
    repeat: int = 5,
    warmup: int = 1,
    max_seconds: float | None = None,
) -> tuple[list[float], object]:
    """Time ``fn`` ``repeat`` times after ``warmup`` untimed calls.

    ``fn`` takes no arguments and does the work to be timed (close over inputs).
    Returns ``(times, last_return_value)`` — the return value lets a caller pull
    out work metrics (e.g. how many events the sim produced).

    ``max_seconds`` is a soft budget: once the elapsed timed total exceeds it we
    stop early (but always do at least one timed repeat), so a single very large
    point can't blow up the whole run.
    """
    for _ in range(max(0, warmup)):
        fn()

    times: list[float] = []
    result = None
    elapsed = 0.0
    for i in range(max(1, repeat)):
        t0 = time.perf_counter()
        result = fn()
        dt = time.perf_counter() - t0
        times.append(dt)
        elapsed += dt
        if max_seconds is not None and elapsed >= max_seconds and i + 1 >= 1:
            break
    return times, result
