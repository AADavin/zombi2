"""The live lineage set — parallel arrays every genome engine grows over, the ``species._grow`` shape.

Both resolutions evolve **all lineages alive at once** along one global clock (a transfer couples two
contemporaneous lineages), so both keep the same three parallel structures: ``alive`` (species node
ids), ``gen`` (each lineage's working genome, at whatever resolution), and ``pos`` (node id → its
index, so a lineage can be retired in O(1)). These two helpers manage that set and are
genome-shape-agnostic — ``gen`` may hold multisets or lists of chromosomes — so they live here, one
home, shared by the unordered core and the ordered engine.
"""

from __future__ import annotations

import functools
import gc


def without_cyclic_gc(fn):
    """Run ``fn`` with Python's *cyclic* garbage collector paused, restoring the caller's setting
    afterwards (even on error).

    A genome run builds millions of ``Event`` and ``GeneCopy`` objects that live until the run ends
    and hold only ints, floats, strings and ``None`` — so they can never form a reference cycle.
    The cyclic collector cannot know that, so as the pile grows it keeps re-scanning all of it looking
    for cycles that by construction do not exist; that scanning measures at roughly a third of a run.
    Pausing it removes the waste. Reference counting still frees every transient object the instant it
    falls out of use, so memory is unaffected, and the collector touches neither the logic nor the
    RNG, so the output is byte-for-byte identical — this is purely a speed change.

    Nested engine calls compose correctly: an inner call sees GC already disabled, so it leaves it to
    the outer call to re-enable."""
    @functools.wraps(fn)
    def wrapper(*args, **kwargs):
        was_enabled = gc.isenabled()
        gc.disable()
        try:
            return fn(*args, **kwargs)
        finally:
            if was_enabled:
                gc.enable()
    return wrapper


def enter(alive, gen, pos, node_id, genome) -> None:
    """Enter a lineage into the alive set with its (mutable) working genome."""
    pos[node_id] = len(alive)
    alive.append(node_id)
    gen.append(genome)


def retire(alive, gen, pos, k) -> None:
    """Retire the lineage at index ``k`` (swap the last into its slot), keeping ``pos`` in sync."""
    removed = alive[k]
    alive[k] = alive[-1]
    gen[k] = gen[-1]
    pos[alive[k]] = k          # the moved lineage now lives at k (a self-assign if k was last)
    alive.pop()
    gen.pop()
    del pos[removed]
