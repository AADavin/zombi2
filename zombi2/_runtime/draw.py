"""Drawing one of several things in proportion to its weight — the one sampler every engine needs.

Whenever a rate differs across the living lineages — because a driver varies it, or because a
per-family draw weights it — the total is summed *with* those weights, so whatever the event lands
on has to be drawn with them too. Summing one way and drawing another would make the rate say one
thing and the run do another.

It lived three times over, once in `zombi2.species`, once in `zombi2.genomes._live` and once in
`zombi2.joint`, and the genome copy was imported across packages from the species one. Here it is
plumbing rather than domain code, which is what this package is for, and it can be imported from
anywhere without a cycle: `zombi2.genomes` imports `zombi2.species`, so the reverse could not.
"""

from __future__ import annotations

__all__ = ["weighted_index"]


def weighted_index(rng, weights, total: float) -> int:
    """The index of one of ``weights``, drawn in proportion to it. ``total`` is their sum, passed in
    because the caller has already computed it to race the events."""
    r = float(rng.random()) * total
    acc = 0.0
    for i, w in enumerate(weights):
        acc += w
        if r < acc:
            return i
    return len(weights) - 1     # floating-point guard: r == total lands on the last one
