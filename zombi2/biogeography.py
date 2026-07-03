"""Historical biogeography — the DEC model of geographic-range evolution.

A species' trait here is its **geographic range**: a non-empty subset of a fixed set of
discrete *areas*. The Dispersal–Extinction–Cladogenesis model (Ree & Smith 2008) evolves it
with two processes:

* **anagenetic** change along branches — *dispersal* (gain an area) and local *extinction*
  (lose an area), a continuous-time Markov chain over the range states, exactly like
  :class:`~zombi2.Mk`;
* **cladogenetic** change at speciation nodes — the ancestral range is inherited/divided
  between the two daughters (narrow sympatry, subset sympatry, or vicariance).

Because of the node process, DEC needs its own driver, :func:`simulate_biogeography`, rather
than :func:`~zombi2.simulate_traits`. The result is a discrete
:class:`~zombi2.traits.TraitResult` whose states are area-label tuples (the ranges).

    import zombi2 as z
    tree = z.simulate_species_tree(z.BirthDeath(1.0, 0.3), n_tips=20, age=5.0, seed=1)
    dec = z.DEC(areas=["A", "B", "C"], dispersal=0.05, extinction=0.1)
    res = z.simulate_biogeography(tree, dec, root_state={"A"}, seed=1)
    res.labeled_values()          # {extant leaf: ('A', 'B') ...} — the observed ranges
"""

from __future__ import annotations

import numpy as np

from .tree import Tree
from .traits import Mk, TraitResult


def _popcount(mask: int) -> int:
    return bin(mask).count("1")


class DEC(Mk):
    """The Dispersal–Extinction–Cladogenesis model (Ree & Smith 2008).

    Ranges are the non-empty subsets of ``areas`` with at most ``max_range_size`` areas. Along a
    branch, a range in area set ``R`` gains an absent area ``a`` (dispersal) at rate
    ``Σ_{b∈R} dispersal[b, a]`` and, if ``|R| ≥ 2``, loses an area ``a`` (local extinction) at
    rate ``extinction[a]`` — a range never becomes empty. This anagenetic process is a
    :class:`~zombi2.Mk` chain over the enumerated ranges. At a speciation the range is split by
    :meth:`cladogenesis`.

    Parameters
    ----------
    areas:
        The number of areas, or a list of area labels.
    dispersal:
        Dispersal rate: a scalar applied to every ordered area pair, or an ``n x n`` matrix
        ``dispersal[i, j]`` (rate of dispersal from area ``i`` to area ``j``).
    extinction:
        Local (per-area) extinction rate: a scalar, or a length-``n`` vector.
    max_range_size:
        Maximum number of areas a range may span (default: all areas).
    root:
        Root policy as in :class:`~zombi2.Mk` (over the enumerated ranges), or a ``set`` /
        ``frozenset`` of area labels giving an explicit root range.
    """

    def __init__(self, areas, dispersal, extinction, max_range_size=None, root="uniform"):
        if isinstance(areas, int):
            n = areas
            labels = list(range(n))
        else:
            labels = list(areas)
            n = len(labels)
        if n < 2:
            raise ValueError("DEC needs at least 2 areas")
        maxr = n if max_range_size is None else int(max_range_size)
        if not (1 <= maxr <= n):
            raise ValueError(f"max_range_size must be in [1, {n}]")

        D = (np.full((n, n), float(dispersal)) if np.isscalar(dispersal)
             else np.asarray(dispersal, dtype=float))
        if D.shape != (n, n):
            raise ValueError("dispersal must be a scalar or an n x n matrix")
        np.fill_diagonal(D, 0.0)
        E = (np.full(n, float(extinction)) if np.isscalar(extinction)
             else np.asarray(extinction, dtype=float))
        if E.shape != (n,):
            raise ValueError("extinction must be a scalar or a length-n vector")
        if np.any(D < 0) or np.any(E < 0):
            raise ValueError("dispersal and extinction rates must be >= 0")

        masks = sorted((m for m in range(1, 1 << n) if _popcount(m) <= maxr),
                       key=lambda m: (_popcount(m), m))
        index_of = {m: i for i, m in enumerate(masks)}
        k = len(masks)
        Q = np.zeros((k, k))
        for i, R in enumerate(masks):
            inside = [a for a in range(n) if R >> a & 1]
            if len(inside) < maxr:                       # dispersal: gain an absent area
                for a in range(n):
                    if not (R >> a & 1):
                        rate = float(sum(D[b, a] for b in inside))
                        if rate > 0:
                            Q[i, index_of[R | (1 << a)]] += rate
            if len(inside) >= 2:                         # local extinction: lose an area
                for a in inside:
                    if E[a] > 0:
                        Q[i, index_of[R & ~(1 << a)]] += E[a]

        self._n = n
        self._labels = labels
        self._maxr = maxr
        self._masks = masks
        self._index_of = index_of
        states = [tuple(labels[a] for a in range(n) if m >> a & 1) for m in masks]

        if isinstance(root, (set, frozenset)):
            root = self.encode(root)                       # clean error for an invalid root range
        super().__init__(Q, states=states, root=root)

    # --- range <-> index ---------------------------------------------------
    def _encode(self, range_labels) -> int:
        mask = 0
        for a in range_labels:
            idx = self._labels.index(a) if a in self._labels else int(a)
            mask |= (1 << idx)
        return mask

    def encode(self, range_labels) -> int:
        """Index of a range given as an iterable of area labels."""
        mask = self._encode(range_labels)
        if mask not in self._index_of:
            raise ValueError("range is empty or exceeds max_range_size")
        return self._index_of[mask]

    # --- cladogenesis ------------------------------------------------------
    def cladogenesis(self, index, rng):
        """Split the range at a speciation into the two daughters' ranges (as state indices).

        A single-area range is inherited whole by both daughters (narrow sympatry). A widespread
        range yields, with equal probability, one daughter with a single area ``a`` and the other
        with either the *full* ancestral range (subset sympatry) or the *complement* ``R\\{a}``
        (vicariance); the daughter order is randomized.
        """
        R = self._masks[index]
        inside = [a for a in range(self._n) if R >> a & 1]
        if len(inside) == 1:
            return index, index
        outcomes = []
        for a in inside:
            single = 1 << a
            outcomes.append((single, R))              # subset sympatry: {a} | full range
            outcomes.append((single, R & ~single))    # vicariance: {a} | complement
        m1, m2 = outcomes[int(rng.integers(len(outcomes)))]
        i1, i2 = self._index_of[m1], self._index_of[m2]
        if rng.random() < 0.5:
            i1, i2 = i2, i1
        return i1, i2

    def __repr__(self) -> str:
        return f"DEC(areas={self._n}, max_range_size={self._maxr})"


def simulate_biogeography(tree: Tree, model: DEC, *, root_state=None,
                         seed: int | None = None, rng=None) -> TraitResult:
    """Evolve a geographic range down ``tree`` under a :class:`DEC` model.

    Ranges evolve by dispersal/extinction along branches (the anagenetic :class:`~zombi2.Mk`
    chain) and are split by :meth:`DEC.cladogenesis` at every speciation, so the two daughters
    can inherit different ranges. Returns a discrete :class:`~zombi2.traits.TraitResult` whose
    ``.labeled_values()`` are the tips' ranges (tuples of area labels), ``.ancestral_states()``
    the ancestral ranges, and ``.history`` the anagenetic (dispersal/extinction) map along each
    branch.

    ``root_state`` sets the root range — an ``Mk`` state index, or an iterable of area labels
    (e.g. ``{"A"}``); if ``None`` it follows the model's root policy.
    """
    if rng is None:
        rng = np.random.default_rng(seed)

    if root_state is None:
        root_idx = model.root_value(rng)
    elif isinstance(root_state, (int, np.integer)):
        root_idx = int(root_state)
    else:
        root_idx = model.encode(root_state)

    node_values: dict = {}
    history: dict = {}
    stack = [(tree.root, root_idx)]
    while stack:
        node, start = stack.pop()
        if node.parent is None:
            end, segs = start, []
        else:
            end, segs = model.evolve(start, node.branch_length(), node.parent.time, rng)
        node_values[node] = end
        history[node] = segs
        children = node.children
        if len(children) == 2:
            r1, r2 = model.cladogenesis(end, rng)
            stack.append((children[0], r1))
            stack.append((children[1], r2))
        elif len(children) == 1:                     # degree-two node: pass the range through
            stack.append((children[0], end))
        elif len(children) > 2:                      # cladogenesis is defined only for bifurcations
            raise ValueError(
                f"DEC requires a binary tree; node {node.name!r} has {len(children)} children "
                "(resolve polytomies before simulating biogeography)"
            )

    return TraitResult(tree=tree, model=model, node_values=node_values,
                       history=history, kind="discrete")
