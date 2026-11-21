"""The live lineage set — parallel arrays every genome engine grows over, the ``species._grow`` shape.

All three resolutions evolve **all lineages alive at once** along one global clock (a transfer couples
two contemporaneous lineages), so all keep the same three parallel structures: ``alive`` (species node
ids), ``gen`` (each lineage's working genome, at whatever resolution), and ``pos`` (node id → its
index, so a lineage can be retired in O(1)). These two helpers manage that set and are
genome-shape-agnostic — ``gen`` may hold multisets or lists of chromosomes — so they live here, one
home, shared by the family core, the ordered engine and the nucleotide one.
"""

from __future__ import annotations

# re-exported: the three genome engines import it from here, and it is plumbing, not genome code
from .._runtime.draw import SummedWeights, WeightedIndex, weighted_index  # noqa: F401

import functools
import gc
import heapq
import math


def without_cyclic_gc(fn):
    """Run ``fn`` with Python's *cyclic* garbage collector paused, restoring the caller's setting
    afterwards (even on error).

    A genome run builds millions of ``GeneEdge`` and ``GeneCopy`` objects that live until the run ends
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


class DriverSwitches:
    """The next moment a finished driver switches on a living lineage, and which lineages switch
    then, without asking every lineage at every event.

    A driven rate is constant on a lineage until its driver's next switch, so a driven Gillespie
    loop has to stop there. Asking each lineage for its next switch at each event costs the whole
    live set per event; the answer changes only when a lineage enters or one of them switches. So
    each living lineage holds one entry per driver in a heap, keyed by its next switch: the top is
    the horizon, and popping what has come due names the lineages whose rate has just changed.

    A lineage that left is dropped when its entry surfaces (``pos`` says who is living), and a
    lineage that comes back, as a family's footprint allows, replaces its entry: ``_armed`` holds
    the one time that counts for each lineage and driver, so an older entry is recognised and skipped.
    """

    __slots__ = ("_trajs", "_heap", "_armed")

    def __init__(self, trajs) -> None:
        self._trajs = trajs
        self._heap: list[tuple[float, int, str]] = []
        self._armed: dict[tuple[int, str], float] = {}

    def _arm(self, node_id: int, key: str, time: float) -> None:
        nxt = self._trajs[key].next_change(node_id, time)
        if nxt < math.inf:
            self._armed[(node_id, key)] = nxt
            heapq.heappush(self._heap, (nxt, node_id, key))
        else:
            self._armed.pop((node_id, key), None)

    def entered(self, node_id: int, time: float) -> None:
        """A lineage joins the live set at ``time``: note each driver's next switch on it."""
        for key in self._trajs:
            self._arm(node_id, key, time)

    def due(self, time: float, pos) -> list[int]:
        """The living lineages on which a driver has switched by ``time``, each one's next switch
        noted in its place. Called before the rates are read at ``time``, so a switch falling on the
        very instant the loop arrived at is taken up there rather than one step late."""
        heap, armed = self._heap, self._armed
        out = []
        while heap and heap[0][0] <= time:
            at, node_id, key = heapq.heappop(heap)
            if armed.get((node_id, key)) != at:
                continue                                    # an older entry of a lineage that came back
            if node_id not in pos:
                del armed[(node_id, key)]                   # the lineage left: nothing to switch
                continue
            out.append(node_id)
            self._arm(node_id, key, time)
        return out

    def next(self, pos) -> float:
        """The next switch on a living lineage, or ``inf``. Call after `due`."""
        heap, armed = self._heap, self._armed
        while heap:
            at, node_id, key = heap[0]
            if armed.get((node_id, key)) == at and node_id in pos:
                return at
            heapq.heappop(heap)
            if armed.get((node_id, key)) == at:
                del armed[(node_id, key)]
        return math.inf

