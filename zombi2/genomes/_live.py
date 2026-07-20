"""The live lineage set — parallel arrays every genome engine grows over, the ``species._grow`` shape.

Both resolutions evolve **all lineages alive at once** along one global clock (a transfer couples two
contemporaneous lineages), so both keep the same three parallel structures: ``alive`` (species node
ids), ``gen`` (each lineage's working genome, at whatever resolution), and ``pos`` (node id → its
index, so a lineage can be retired in O(1)). These two helpers manage that set and are
genome-shape-agnostic — ``gen`` may hold multisets or lists of chromosomes — so they live here, one
home, shared by the unordered core and the ordered engine.
"""

from __future__ import annotations


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
