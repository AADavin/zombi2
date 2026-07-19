"""Traits — a value riding the species tree (level 4).

A trait is not a genealogy of events like the other three levels; it is a **value that rides the
tree** — a body size, a habitat, a presence/absence — and you observe the value itself, not an event
count (``docs/design/trait-api.md``). So the trait level has no "rate of events": its compact
**source of truth is the value at every node** (``node_values``), not an event log, and the rich
views — ``values`` at the extant tips, the discrete stochastic ``history``, the realized ``events`` —
are derived from it. That is a real seam, named rather than papered over. What keeps traits inside the
one framework is that the *ways* a value evolves reuse the same ``scope(base) × modifiers`` rate
grammar (SPEC §5).

This is the **continuous** trait level — ``simulate_continuous`` — and its three variants are the
same diffusion wearing different knobs, not three classes (SPEC §4):

- **Brownian motion**, the native process: over a branch the value moves by ``Normal(0, σ²·dt)``, so
  node-by-node in preorder it reproduces the exact tip law (Felsenstein 1985): the extant tips are
  multivariate-normal with variance ``σ² ×`` (root-to-tip depth) and covariance ``σ² ×`` (shared
  path length). ``rate`` is the variance-rate σ².
- **Ornstein–Uhlenbeck**: add ``reverts_to`` (the optimum θ) and ``pull`` (the strength α) and the
  diffusion is pulled toward θ — stabilizing selection. The exact per-branch transition is normal
  with mean ``θ + (x−θ)·e^{−α·dt}`` and variance ``σ²/(2α)·(1−e^{−2α·dt})``. (These are the same two
  knobs the CIR clock grows one level over — a shared vocabulary, not shared code.)
- **Early burst / ACDC**: give ``rate`` a ``Time`` skyline (``rate = σ² * mod.Time({0: 1, 5: 0.2})``)
  and σ² changes through time — the *same* ``Time`` modifier that gives the species tree its skyline.
  The per-branch variance is then the exact integral ``∫ σ²(t) dt`` over the branch.

``rate`` is *per lineage*: each lineage carries its own independent diffusion, never pooled across the
tree. That non-pooling is the trait seam in the rate grammar — the engine evaluates the rate one
lineage at a time (``lineages=1``), where the event levels sum a per-unit rate over everything alive
at once. (OU with a time-varying σ² — the two knob-sets at once — is deferred; use one or the other.)

Still to come, each its own slice: the discrete twin ``simulate_discrete`` (Mk / threshold) with its
stochastic-map ``history``; the ``correlation=`` overlay for traits that drift together; and the
named-and-deferred cases (``at_speciation`` jumps, ``regimes``, hidden states, DEC → experimental).
SSE (BiSSE/MuSSE/QuaSSE) is **not** a trait model — it is trait↔species *joint*, Part III.
"""

from __future__ import annotations

import math
import pathlib
from dataclasses import dataclass

import numpy as np

from ..rates.modifiers import Time
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


def _accrued_variance(rate, t0: float, t1: float) -> float:
    """The variance a diffusing trait accrues over a branch spanning ``[t0, t1]`` — the integral
    ``∫ σ²(t) dt`` of the variance-rate. For a bare σ² this is ``σ²·(t1−t0)`` (Brownian motion); for a
    ``Time`` skyline (early burst) it sums σ² over each interval the branch crosses, stepping at the
    schedule's breakpoints. The same breakpoint walk the species/genome engines use — integrated over
    the branch rather than sampled at a point (σ² is piecewise-constant, so the integral is exact)."""
    total = 0.0
    t = t0
    while t < t1:
        nxt = min(rate.next_change(t), t1)  # constant rate → inf → one step of length (t1−t0)
        total += rate.effective(lineages=1, time=t) * (nxt - t)
        t = nxt
    return total


def simulate_continuous(tree, *, start=0.0, rate=1.0, reverts_to=None, pull=None,
                        seed=None) -> TraitsResult:
    """Evolve a continuous trait down a tree and return a :class:`TraitsResult`. One process, three
    variants selected by knobs (SPEC §4): **Brownian motion** (bare ``rate``), **Ornstein–Uhlenbeck**
    (add ``reverts_to`` + ``pull``), **early burst** (give ``rate`` a ``Time`` skyline).

    ``tree`` is the **complete** species tree (a :class:`~zombi2.species.Tree`, or a
    :class:`~zombi2.species.SpeciesResult` whose ``complete_tree`` is used). The trait evolves on
    **every** lineage, extant and extinct alike, so the ancestral states are exact and complete; the
    observed dataset is the extant tips, ``result.values``.

    ``start`` is the value at ``t = 0`` (the crown origin, ``root.birth_time``): the root lineage
    diffuses over its own branch ``[0, first split]`` like any other, so a trait and a genome evolve
    over the **same** branch set, and each node's stored value is the trait at that node's
    ``end_time`` (``node_values[root]`` is the value at the first split, not ``start``).

    ``rate`` is the variance-rate σ² (a ``scope(base) × modifiers`` rate spec), *per lineage*: each
    lineage diffuses independently at σ², never pooled across the tree. A bare number is Brownian
    motion (``Normal(0, σ²·dt)`` over a branch); a ``Time`` modifier makes σ² change through time —
    early burst / ACDC — with the per-branch variance the exact integral ``∫ σ²(t) dt``.

    ``reverts_to`` (the optimum θ) and ``pull`` (the strength α > 0) turn the diffusion into
    Ornstein–Uhlenbeck — the value is pulled toward θ while it diffuses, the exact per-branch
    transition being ``Normal(θ + (x−θ)·e^{−α·dt}, σ²/(2α)·(1−e^{−2α·dt}))``. Give **both** or
    neither. OU with a time-varying σ² (both knob-sets at once) is not wired yet — use one or the
    other. Deterministic given ``seed``.
    """
    tree = tree.complete_tree if isinstance(tree, SpeciesResult) else tree
    if isinstance(start, bool) or not isinstance(start, (int, float)) or not math.isfinite(start):
        raise ValueError(f"start must be a finite number, got {start!r}")
    r = as_rate(rate, default_scope=PerLineage)
    if not isinstance(r.scope, PerLineage):
        raise ValueError(
            f"rate has a {type(r.scope).__name__} scope, but a continuous trait's variance-rate is "
            f"per lineage — drop the scope wrapper (per lineage is the default)."
        )
    # Time (early burst) is wired; every other modifier (Diversity, Inherited / clade drift) is a
    # later slice, so reject it loudly rather than silently ignore it — the genome engine's discipline.
    for m in r.modifiers:
        if not isinstance(m, Time):
            raise ValueError(
                f"rate carries {type(m).__name__}, which the continuous trait engine does not support "
                f"yet — only Time (early burst) is wired (Diversity / Inherited are later slices)."
            )
    has_time = bool(r.modifiers)  # all remaining modifiers are Time skylines

    # OU: reverts_to (θ) + pull (α) turn the diffusion into mean-reversion — both or neither.
    is_ou = reverts_to is not None or pull is not None
    if is_ou:
        if reverts_to is None or pull is None:
            raise ValueError(
                "Ornstein–Uhlenbeck needs both reverts_to (the optimum) and pull (the strength); "
                "give both, or neither for Brownian motion."
            )
        if isinstance(reverts_to, bool) or not isinstance(reverts_to, (int, float)) \
                or not math.isfinite(reverts_to):
            raise ValueError(f"reverts_to must be a finite number, got {reverts_to!r}")
        if isinstance(pull, bool) or not isinstance(pull, (int, float)) \
                or not math.isfinite(pull) or pull <= 0:
            raise ValueError(
                f"pull must be a finite positive number (omit it for Brownian motion), got {pull!r}"
            )
        if has_time:
            raise ValueError(
                "a time-varying variance-rate (early burst, via Time) combined with OU "
                "(reverts_to / pull) is not wired yet — use one or the other this slice."
            )
        theta, alpha = float(reverts_to), float(pull)
        sigma2 = r.effective(lineages=1)  # σ² is constant under OU (Time is rejected above)

    rng = np.random.default_rng(seed)
    node_values: dict[int, float] = {}
    for i in _preorder(tree):
        node = tree.nodes[i]
        # the root starts from `start` at t=0; every other node from its parent's end value (parent
        # < i, already set). One uniform rule: node_values[i] is the trait at node i's end_time.
        x = float(start) if node.parent is None else node_values[node.parent]
        t0, t1 = node.birth_time, node.end_time
        if is_ou:
            e = math.exp(-alpha * (t1 - t0))       # mean-reversion toward θ over the branch
            mean = theta + (x - theta) * e
            var = sigma2 / (2.0 * alpha) * (1.0 - e * e)
        else:
            mean = x                                # pure diffusion (BM / early burst)
            var = _accrued_variance(r, t0, t1)
        std = math.sqrt(var) if var > 0.0 else 0.0
        node_values[i] = mean + (float(rng.normal(0.0, std)) if std > 0.0 else 0.0)

    return TraitsResult(tree, node_values, seed)


__all__ = ["simulate_continuous", "TraitsResult"]
