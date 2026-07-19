"""Traits — a value riding the species tree (level 4).

A trait is not a genealogy of events like the other three levels; it is a **value that rides the
tree** — a body size, a habitat, a presence/absence — and you observe the value itself, not an event
count (``docs/design/trait-api.md``). So the trait level has no "rate of events": its compact
**source of truth is the value at every node** (``node_values``), not an event log, and the rich
views — ``values`` at the extant tips, the discrete stochastic ``history``, the realized ``events`` —
are derived from it. That is a real seam, named rather than papered over. What keeps traits inside the
one framework is that the *ways* a value evolves reuse the same ``scope(base) × modifiers`` rate
grammar (SPEC §5).

This slice is the **continuous** core: Brownian motion, the native process. A continuous trait
diffuses along each branch — over a branch of duration ``dt`` the value moves by ``Normal(0, σ²·dt)``
— so simulated node-by-node in preorder it reproduces the exact tip law (Felsenstein 1985): the
extant tips are multivariate-normal with variance ``σ² ×`` (root-to-tip depth) and covariance
``σ² ×`` (shared path length). ``rate`` is the variance-rate σ², a rate spec that is *per lineage*:
each lineage carries its own independent diffusion, never pooled across the tree. That non-pooling is
the trait seam in the rate grammar — the engine evaluates the rate one lineage at a time
(``lineages=1``), where the event levels sum a per-unit rate over everything alive at once.

Still to come, each its own slice: the OU (``reverts_to`` / ``pull``) and early-burst (``Time``)
knobs on this same function; the discrete twin ``simulate_discrete`` (Mk / threshold) with its
stochastic-map ``history``; the ``correlation=`` overlay for traits that drift together; and the
named-and-deferred cases (``at_speciation`` jumps, ``regimes``, hidden states, DEC → experimental).
SSE (BiSSE/MuSSE/QuaSSE) is **not** a trait model — it is trait↔species *joint*, Part III.
"""

from __future__ import annotations

import math
import pathlib
from dataclasses import dataclass

import numpy as np

from ..rates.rate import as_rate
from ..rates.scope import PerLineage
from ..species import SpeciesResult, Tree

_WRITE_OUTPUTS = ("values",)  # the write vocabulary this slice supports


@dataclass
class TraitsResult:
    """What ``simulate_continuous`` returns: the ``complete_tree`` it ran on, ``node_values`` at
    **every** node (the trait's compact source of truth — extant, extinct, and internal alike), the
    ``seed``, and (for a discrete trait, a later slice) the stochastic-map ``history``. The observed
    trait dataset is the extant tips, ``.values``; ``.write`` materialises the chosen outputs.

    The trait seam: unlike the event-log levels, ``node_values`` *is* the source of truth here — a
    continuous value has no instantaneous events to log — so ``.events`` is a derived view (the
    realized discrete state-changes) and is **empty for a continuous trait**.
    """

    complete_tree: Tree
    node_values: dict[int, float]
    seed: int | None
    kind: str = "continuous"
    history: dict[int, list] | None = None

    @property
    def values(self) -> dict[int, float]:
        """The observed trait dataset — the value at each **extant** tip (the comparative-data
        vector). Internal and extinct nodes keep their exact ancestral / lineage values in
        ``node_values``."""
        return {n.id: self.node_values[n.id] for n in self.complete_tree.extant()}

    @property
    def events(self) -> list:
        """The realized discrete state-changes along the tree — a derived view, **empty for a
        continuous trait** (which diffuses with no instantaneous events). It is populated when the
        discrete twin (Mk / threshold) lands in a later slice."""
        return []

    def write(self, directory, outputs=_WRITE_OUTPUTS) -> None:
        """Write chosen ``outputs`` to ``directory`` (created if needed): ``"values"`` →
        ``trait_values.tsv``, the ``node<TAB>trait`` table over the extant tips."""
        unknown = [o for o in outputs if o not in _WRITE_OUTPUTS]
        if unknown:
            raise ValueError(f"unknown write outputs {unknown}; choose from {list(_WRITE_OUTPUTS)}")
        d = pathlib.Path(directory)
        d.mkdir(parents=True, exist_ok=True)
        if "values" in outputs:
            (d / "trait_values.tsv").write_text(_values_tsv(self.values))


def _values_tsv(values: dict[int, float]) -> str:
    """The extant-tip values as a two-column ``node<TAB>trait`` table, one row per tip in id order.
    Tips are named ``n<id>`` to match the tree's Newick leaf labels."""
    rows = ["node\ttrait"]
    for i in sorted(values):
        rows.append(f"n{i}\t{values[i]:.6g}")
    return "\n".join(rows) + "\n"


def _preorder(tree: Tree) -> list[int]:
    """Node ids in an order that visits every node **after its parent** (a valid preorder). The
    forward engine always gives a child a higher id than its parent, so ascending id order suffices
    — the same monotonic-id fact ``genomes.prune`` relies on in reverse. No recursion needed."""
    return sorted(tree.nodes)


def simulate_continuous(tree, *, start=0.0, rate=1.0, seed=None) -> TraitsResult:
    """Evolve a continuous trait down a tree by Brownian motion and return a :class:`TraitsResult`.

    ``tree`` is the **complete** species tree (a :class:`~zombi2.species.Tree`, or a
    :class:`~zombi2.species.SpeciesResult` whose ``complete_tree`` is used). The trait evolves on
    **every** lineage, extant and extinct alike, so the ancestral states are exact and complete;
    the observed dataset is the extant tips, ``result.values``.

    ``start`` is the value at ``t = 0`` (the crown origin, ``root.birth_time``): the root lineage
    diffuses over its own branch ``[0, first split]`` like any other, so a trait and a genome evolve
    over the **same** branch set (the genome engine fires D/T/L/O on that branch too), and each
    node's stored value is the trait at that node's ``end_time`` — ``node_values[root]`` is the value
    at the first split, not ``start``. ``rate`` is the variance-rate σ² (a ``scope(base) ×
    modifiers`` rate spec — a bare number this slice; the Time / early-burst and OU knobs are later
    slices). It is *per lineage*: each lineage diffuses independently at σ², never pooled across the
    tree — so over a branch of duration ``dt`` the value moves by ``Normal(0, σ²·dt)``.
    Deterministic given ``seed``.
    """
    tree = tree.complete_tree if isinstance(tree, SpeciesResult) else tree
    if isinstance(start, bool) or not isinstance(start, (int, float)) or not math.isfinite(start):
        raise ValueError(f"start must be a finite number, got {start!r}")
    # this slice wires only the default scope (per lineage) and no modifiers; a non-default scope or
    # any modifier (Time = early burst, OU knobs) is a later slice, so reject it loudly rather than
    # silently ignore it — the same discipline the genome engine uses per slice.
    r = as_rate(rate, default_scope=PerLineage)
    if not isinstance(r.scope, PerLineage):
        raise ValueError(
            f"rate has a {type(r.scope).__name__} scope, but a continuous trait's variance-rate is "
            f"per lineage — drop the scope wrapper (per lineage is the default)."
        )
    if r.modifiers:
        carried = ", ".join(type(m).__name__ for m in r.modifiers)
        raise ValueError(
            f"rate carries {carried}, which this slice does not wire — only a bare variance-rate σ² "
            f"is supported (Time / early burst and the OU knobs are later slices)."
        )
    sigma2 = r.effective(lineages=1)  # the per-lineage variance-rate; constant with no modifiers

    rng = np.random.default_rng(seed)
    node_values: dict[int, float] = {}
    for i in _preorder(tree):
        node = tree.nodes[i]
        # the root starts from `start` at t=0 and diffuses over its own branch; every other node
        # starts from its parent's end value (parent < i, so it is already set). One uniform rule:
        # node_values[i] is the trait at node i's end_time.
        x = float(start) if node.parent is None else node_values[node.parent]
        dt = node.end_time - node.birth_time
        std = math.sqrt(sigma2 * dt) if sigma2 > 0.0 and dt > 0.0 else 0.0
        node_values[i] = x + (float(rng.normal(0.0, std)) if std > 0.0 else 0.0)

    return TraitsResult(tree, node_values, seed)


__all__ = ["simulate_continuous", "TraitsResult"]
