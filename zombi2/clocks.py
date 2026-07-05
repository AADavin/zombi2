"""Relaxed-molecular-clock namespace (scikit-learn-style).

Re-exports the clock family so users can write::

    from zombi2.clocks import StrictClock, UncorrelatedLogNormalClock, RateVariation

Every name here is the *same object* as the corresponding top-level ``zombi2`` attribute --
this module is a thin, additive namespace over the implementation module
(:mod:`zombi2.rate_variation`); it does not redefine anything. The clocks all share the
:class:`~zombi2.rate_variation.Clock` interface and turn a chronogram (timetree) into a
phylogram (branch lengths in expected substitutions per site).
"""

from __future__ import annotations

from .rate_variation import (
    Clock,
    RateScaledTree,
    StrictClock,
    UncorrelatedLogNormalClock,
    UncorrelatedGammaClock,
    WhiteNoiseClock,
    AutocorrelatedLogNormalClock,
    CIRClock,
    RateVariation,
)

__all__ = [
    "Clock",
    "RateScaledTree",
    "StrictClock",
    "UncorrelatedLogNormalClock",
    "UncorrelatedGammaClock",
    "WhiteNoiseClock",
    "AutocorrelatedLogNormalClock",
    "CIRClock",
    "RateVariation",
]
