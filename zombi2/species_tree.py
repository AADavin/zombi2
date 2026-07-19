"""Species trees — the forward birth-death engine (species slice 1: the engine only).

Minimal on purpose: constant **per-lineage** birth and death grow a tree forward in time
and record every speciation and extinction. Rates are plain numbers here; the scope
wrappers and modifiers (``Time`` / ``Diversity`` / ``Inherited``) wire in as the next
slice, and the reconstructed tree, Newick, sampling and fossils follow after that. The
final home is ``zombi2.species``; it lives here for now so it cannot collide with the
old package while the rewrite is in flight.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np


@dataclass
class Node:
    """One lineage segment: born at ``birth_time``, ended at ``end_time`` by a split, a
    death, or reaching the present. A split has two ``children``; a tip has none."""

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

    def tips(self) -> list[Node]:
        return [n for n in self.nodes.values() if n.children is None]

    def extant(self) -> list[Node]:
        return [n for n in self.nodes.values() if n.fate == "extant"]

    def extinct(self) -> list[Node]:
        return [n for n in self.nodes.values() if n.fate == "extinct"]


@dataclass
class SpeciesResult:
    """Minimal result: the complete tree, the event log, the seed. (The reconstructed tree,
    Newick, sampling and the full <Level>Result spine are later slices.)"""

    complete: Tree
    events: list[Event]
    seed: int | None

    @property
    def n_extant(self) -> int:
        return len(self.complete.extant())


def _check_rate(name: str, value: float) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise TypeError(f"{name} must be a real number, got {value!r}")
    if not math.isfinite(value) or value < 0:
        raise ValueError(f"{name} must be finite and non-negative, got {value!r}")
    return float(value)


_MAX_ATTEMPTS = 1000  # survival-conditioned retries before giving up on n_tips


def _grow(rng, birth: float, death: float, n_tips: int | None, age: float | None) -> tuple[Tree, list[Event]]:
    """Grow one forward birth-death tree until it hits ``n_tips`` living lineages, reaches
    ``age``, or dies out. Returns the complete tree and the event log."""
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
        if n_tips is not None and n >= n_tips:
            break
        total = (birth + death) * n
        if total <= 0.0:  # nothing can happen (e.g. birth=death=0)
            break
        t_next = t + float(rng.exponential(1.0 / total))
        if age is not None and t_next >= age:
            t = age
            break
        t = t_next

        i = int(rng.integers(n))
        node = alive[i]
        alive[i] = alive[-1]  # swap-remove keeps picks O(1) and reproducible
        alive.pop()

        if rng.random() < birth / (birth + death):
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

    for i in alive:  # whoever is still alive reached the present
        nodes[i].end_time = t
        nodes[i].fate = "extant"

    return Tree(nodes, root), events


def simulate_species_tree(birth, death=0.0, *, n_tips=None, age=None, seed=None) -> SpeciesResult:
    """Grow a forward birth-death tree at constant per-lineage rates.

    Each living lineage speciates at rate ``birth`` and dies at rate ``death`` (Yule =
    ``death=0``). Stop at exactly ``n_tips`` living lineages, **or** at time ``age`` —
    give exactly one. Deterministic given ``seed``.

    ``n_tips`` is **conditioned on survival**: a birth-death tree can die out, so we
    restart (advancing the same generator) until one reaches ``n_tips``. ``age`` is not
    conditioned — the tree grows for the given time and may die out.
    """
    birth = _check_rate("birth", birth)
    death = _check_rate("death", death)
    if (n_tips is None) == (age is None):
        raise ValueError("give exactly one of n_tips or age")
    if n_tips is not None and (isinstance(n_tips, bool) or not isinstance(n_tips, int) or n_tips < 1):
        raise ValueError(f"n_tips must be a positive integer, got {n_tips!r}")
    if age is not None and (not isinstance(age, (int, float)) or not math.isfinite(age) or age <= 0):
        raise ValueError(f"age must be a positive finite number, got {age!r}")

    rng = np.random.default_rng(seed)

    if age is not None:
        tree, events = _grow(rng, birth, death, None, age)
        return SpeciesResult(tree, events, seed)

    for _ in range(_MAX_ATTEMPTS):
        tree, events = _grow(rng, birth, death, n_tips, None)
        if sum(1 for nd in tree.nodes.values() if nd.fate == "extant") == n_tips:
            return SpeciesResult(tree, events, seed)
    raise RuntimeError(
        f"could not grow a tree to {n_tips} tips in {_MAX_ATTEMPTS} attempts; "
        "birth must comfortably exceed death for large n_tips"
    )


__all__ = ["simulate_species_tree", "SpeciesResult", "Tree", "Node", "Event"]
