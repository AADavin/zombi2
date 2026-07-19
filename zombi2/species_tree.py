"""Species trees — the forward birth-death engine.

Constant-per-lineage birth and death grow a tree forward in time and record every
speciation and extinction. ``birth`` and ``death`` are full **rate specs** — a number,
a scope wrapper, or a product — so ``birth = scope.Global(1.0)`` gives one shared
tree-wide budget (linear growth) and ``birth = 1.0 * mod.Diversity(cap=100)`` slows the
tree as it fills up. The rate is re-evaluated after every event, which is exact for
``scope`` and ``Diversity`` (they only change at events); ``Time`` needs the
interval-aware sampler and is the next slice, so the engine does not pass ``time`` yet.

Still to come: the extant tree, Newick, sampling, fossils, and the move to
``zombi2.species``. This lives here for now so the old package is untouched.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from .rate import as_rate
from .scope import PerLineage


@dataclass
class Node:
    """One lineage segment: born at ``birth_time``, ended at ``end_time`` by a split, a
    death, or reaching the present. A split has two ``children``; a leaf has none."""

    id: int
    parent: int | None
    birth_time: float
    end_time: float | None = None
    children: tuple[int, int] | None = None
    fate: str = "alive"  # alive → "extant" | "extinct"; internal splits are "speciation"


@dataclass(frozen=True)
class Event:
    """A recorded event in the true history: a speciation (with its two children) or an extinction."""

    time: float
    kind: str  # "speciation" | "extinction"
    node: int
    children: tuple[int, int] | None = None


@dataclass
class Tree:
    """The complete tree: every lineage that ever lived, keyed by id, rooted at ``root``."""

    nodes: dict[int, Node]
    root: int

    def leaves(self) -> list[Node]:
        """Every lineage with no descendants — extant **and** extinct."""
        return [n for n in self.nodes.values() if n.children is None]

    def extant(self) -> list[Node]:
        """The lineages alive at the present."""
        return [n for n in self.nodes.values() if n.fate == "extant"]

    def extinct(self) -> list[Node]:
        """The lineages that died before the present."""
        return [n for n in self.nodes.values() if n.fate == "extinct"]


@dataclass
class SpeciesResult:
    """Minimal result: the complete tree, the event log, the seed. (The extant tree,
    Newick, sampling and the full <Level>Result spine are later slices.)"""

    complete_tree: Tree
    events: list[Event]
    seed: int | None

    @property
    def n_extant(self) -> int:
        return len(self.complete_tree.extant())


_MAX_ATTEMPTS = 1000  # survival-conditioned retries before giving up on n_extant


def _grow(rng, birth_rate, death_rate, n_extant: int | None, age: float | None) -> tuple[Tree, list[Event]]:
    """Grow one forward birth-death tree until it reaches ``n_extant`` living lineages,
    reaches ``age``, or dies out. Returns the complete tree and the event log."""
    nodes: dict[int, Node] = {}
    counter = 0

    def new_node(parent: int | None, t: float) -> int:
        nonlocal counter
        i = counter
        counter += 1
        nodes[i] = Node(i, parent, t)
        return i

    root = new_node(None, 0.0)
    alive = [root]  # a list so picks are reproducible given the seed
    t = 0.0
    events: list[Event] = []

    while alive:
        n = len(alive)
        if n_extant is not None and n >= n_extant:
            break
        # standing diversity = the living lineages; the scope reads `lineages`, Diversity `diversity`
        ctx = {"lineages": n, "diversity": n, "time": t}
        total_birth = birth_rate.effective(**ctx)
        total_death = death_rate.effective(**ctx)
        total = total_birth + total_death
        # the total rate is constant until the next skyline breakpoint (or the age limit)
        next_change = min(birth_rate.next_change(t), death_rate.next_change(t))
        horizon = next_change if age is None else min(age, next_change)

        if total > 0.0:
            t_event = t + float(rng.exponential(1.0 / total))
            if t_event < horizon:  # an event fires before the rate changes
                t = t_event
                i = int(rng.integers(n))
                node = alive[i]
                alive[i] = alive[-1]  # swap-remove keeps picks O(1) and reproducible
                alive.pop()
                if rng.random() < total_birth / total:
                    nodes[node].end_time = t
                    nodes[node].fate = "speciation"
                    c1, c2 = new_node(node, t), new_node(node, t)
                    nodes[node].children = (c1, c2)
                    alive.extend((c1, c2))
                    events.append(Event(t, "speciation", node, (c1, c2)))
                else:
                    nodes[node].end_time = t
                    nodes[node].fate = "extinct"
                    events.append(Event(t, "extinction", node))
                continue

        # no event fired before the horizon
        if math.isinf(horizon):
            break  # no age limit and the rate never changes again → nothing more can happen
        if age is not None and horizon == age:
            t = age
            break
        t = horizon  # a skyline breakpoint: advance and re-evaluate the (now changed) rate

    for i in alive:  # whoever is still alive reached the present
        nodes[i].end_time = t
        nodes[i].fate = "extant"

    return Tree(nodes, root), events


def simulate_species_tree(birth, death=0.0, *, n_extant=None, age=None, seed=None) -> SpeciesResult:
    """Grow a forward birth-death tree.

    ``birth`` and ``death`` are rate specs (a number, a ``scope`` wrapper, or a product
    with modifiers); the default scope is **per lineage** (each lineage speciates/dies at
    the base rate, so the tree grows exponentially). Yule = ``death=0``.

    Stop at exactly ``n_extant`` living lineages, **or** at time ``age`` — give exactly
    one. ``n_extant`` is **conditioned on survival**: a birth-death tree can die out, so we
    restart (advancing the same generator) until one reaches ``n_extant``. ``age`` is not
    conditioned. Deterministic given ``seed``.
    """
    birth_rate = as_rate(birth, default_scope=PerLineage)
    death_rate = as_rate(death, default_scope=PerLineage)
    if (n_extant is None) == (age is None):
        raise ValueError("give exactly one of n_extant or age")
    if n_extant is not None and (isinstance(n_extant, bool) or not isinstance(n_extant, int) or n_extant < 1):
        raise ValueError(f"n_extant must be a positive integer, got {n_extant!r}")
    if age is not None and (not isinstance(age, (int, float)) or not math.isfinite(age) or age <= 0):
        raise ValueError(f"age must be a positive finite number, got {age!r}")

    rng = np.random.default_rng(seed)

    if age is not None:
        tree, events = _grow(rng, birth_rate, death_rate, None, age)
        return SpeciesResult(tree, events, seed)

    for _ in range(_MAX_ATTEMPTS):
        tree, events = _grow(rng, birth_rate, death_rate, n_extant, None)
        if sum(1 for nd in tree.nodes.values() if nd.fate == "extant") == n_extant:
            return SpeciesResult(tree, events, seed)
    raise RuntimeError(
        f"could not grow a tree to {n_extant} extant lineages in {_MAX_ATTEMPTS} attempts; "
        "birth must comfortably exceed death for large n_extant"
    )


__all__ = ["simulate_species_tree", "SpeciesResult", "Tree", "Node", "Event"]
