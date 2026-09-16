"""The sequences level's per-family rate variation — one substitution factor per gene family.

The `clock` rides **lineages**: every gene passing through a species branch evolves at that branch's
rate, so a fast species is fast in all of its genes. This module is the other axis. A rate written
``PerSite(1.0).varying_among('families', LogNormal(0.0, 0.8))`` draws one factor per **family**,
before any site evolves, and the family keeps it for the whole of its life. A ribosomal protein and a
phage tail gene then run at different speeds down the same branch.

The two axes are separate and they multiply, which is what SPEC §5 says modifiers do::

    branch = substitution · Δt · lineage clock · family factor

That is the point of having both. Under the shared clock alone a family's depth and its breadth are
locked together: a family that spans many species is always deep. Letting each family carry its own
speed takes them apart, so a slow family can be wide and shallow and a fast one narrow and deep.

The factor multiplies the rate's **base**, so it reaches everything the base reaches — the alignment,
the phylogram, and the indel rates, which are relative to substitution. There is one entry point,
`family_factors`, and the engine folds its number into the per-family base it already computes, so
the alignment and the phylogram cannot disagree about the tree the sequences were drawn along.

The factors are **written**: ``family_multipliers.tsv``, one row per family, in the format the
genomes level writes its event-rate multipliers in (`zombi2.genomes.multipliers`). A method that
fits a rate per family can then be compared against the rates each family was simulated with.
"""

from __future__ import annotations

import math

from ..params.evaluate import values_at_birth


def family_factors(mods, families, rng) -> "dict[int, float]":
    """One factor per family, drawn once and then fixed: ``{family: factor}``.

    ``mods`` is what the substitution rate carries among families, from `Rate.carried_modifiers`, and
    **every** one of them is drawn and multiplied in — taking only the first is how a second one
    silently leaves the model. ``{}`` for a rate carrying none, and then no randomness is consumed, so
    a run without a family draw is bit-identical to one from before this existed.

    ``families`` are drawn in sorted order, the same order the engine evolves them in, so the run is
    reproducible from its seed. Factored out so the serial loop and the parallel engine draw the same
    numbers — each hands in its own ``rng``, and the parallel engine draws here in the parent, which
    is what keeps it worker-count invariant.
    """
    if not mods:
        return {}
    return {f: math.prod(values_at_birth(mods, rng)) for f in sorted(families)}


__all__ = ["family_factors"]
