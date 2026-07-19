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
- **Variable-rates BM** ("ClaDS for traits"): give ``rate`` an ``Inherited`` modifier
  (``rate = σ² * mod.Inherited(spread=0.3)``) and σ² drifts branch-to-branch — each lineage inherits
  its parent's σ² times a lognormal kick at the split — the *same* ``Inherited`` modifier that drifts
  the species rate (ClaDS) and the autocorrelated clock, one level over. (``reverts_to`` / ``pull`` are
  OU function arguments that revert the trait *value*, **not** a modifier — a rate modifier reverts a
  *rate*, which is the sequences level's CIR clock, a different mechanism.)

``rate`` is *per lineage*: each lineage carries its own independent diffusion, never pooled across the
tree. That non-pooling is the trait seam in the rate grammar — the engine evaluates the rate one
lineage at a time (``lineages=1``), where the event levels sum a per-unit rate over everything alive
at once. (OU with a time-varying σ² — the two knob-sets at once — is deferred; use one or the other.)

The **discrete** twin is ``simulate_discrete`` — a state switching along the tree (the Mk model). Its
jumps are simulated *exactly* by the Gillespie algorithm along every branch, so each node's
``(state, duration)`` segments *are* the realized history (a stochastic character map, ``.history``)
and ``.events`` reads off the transitions — the trait level's first genuine event log. ``switch``
gives the rates (symmetric shortcut, ``{"a->b": rate}`` dict, or a ``k×k`` matrix).

Still to come, each its own slice: the ``correlation=`` overlay for traits that drift together, and
the threshold model (a discrete state read off a continuous liability) that underpins its discrete
case; then the named-and-deferred cases (``at_speciation`` jumps, ``regimes``, hidden states, DEC →
experimental). SSE (BiSSE/MuSSE/QuaSSE) is **not** a trait model — it is trait↔species *joint*, Part III.
"""

from __future__ import annotations

import math
import pathlib
from dataclasses import dataclass

import numpy as np

from ..rates.modifiers import Inherited, Time
from ..rates.rate import as_rate
from ..rates.scope import PerLineage
from ..species import SpeciesResult, Tree

_WRITE_OUTPUTS = ("values", "changes")  # the write vocabulary; "changes" is the discrete event log


@dataclass(frozen=True)
class Change:
    """A realized discrete-trait transition — an event of the stochastic character map: on lineage
    ``lineage`` at ``time`` (crown-forward, the species-tree clock), the state changed from
    ``from_state`` to ``to_state``. Derived from :attr:`TraitsResult.history`; a continuous trait has
    none."""

    time: float
    lineage: int
    from_state: object
    to_state: object


@dataclass
class TraitsResult:
    """What ``simulate_continuous`` / ``simulate_discrete`` returns: the ``complete_tree`` it ran on,
    ``node_values`` at **every** node (the trait's compact source of truth — extant, extinct, and
    internal alike; a float for a continuous trait, a state label for a discrete one), the ``seed``,
    the ``kind`` (``"continuous"`` / ``"discrete"``), and — for a discrete trait — the stochastic-map
    ``history`` (each node's branch as ``(state, duration)`` segments). The observed trait dataset is
    the extant tips, ``.values``; ``.write`` materialises the chosen outputs.

    The trait seam: unlike the event-log levels, ``node_values`` (continuous) / ``history``
    (discrete) *is* the source of truth here, and ``.events`` — the realized discrete state-changes —
    is a **derived view**, empty for a continuous trait (which diffuses with no instantaneous events).
    """

    complete_tree: Tree
    node_values: dict[int, object]
    seed: int | None
    kind: str = "continuous"
    history: dict[int, list] | None = None

    @property
    def values(self) -> dict[int, object]:
        """The observed trait dataset — the value at each **extant** tip (the comparative-data
        vector): a float for a continuous trait, a state label for a discrete one. Internal and
        extinct nodes keep their exact ancestral / lineage values in ``node_values``."""
        return {n.id: self.node_values[n.id] for n in self.complete_tree.extant()}

    @property
    def events(self) -> list[Change]:
        """The realized discrete state-changes across the whole tree, in time order — the events of
        the stochastic character map, derived from ``history``. **Empty for a continuous trait**
        (which diffuses with no instantaneous events)."""
        if self.history is None:
            return []
        out: list[Change] = []
        for i in self.history:
            segs = self.history[i]
            t = self.complete_tree.nodes[i].birth_time
            for (s_from, dur), (s_to, _d) in zip(segs, segs[1:]):
                t += dur
                out.append(Change(t, i, s_from, s_to))
        out.sort(key=lambda c: c.time)
        return out

    def write(self, directory, outputs=("values",)) -> None:
        """Write chosen ``outputs`` to ``directory`` (created if needed): ``"values"`` →
        ``trait_values.tsv`` (the ``node<TAB>trait`` table over the extant tips); ``"changes"`` →
        ``trait_changes.tsv`` (the realized discrete transitions — header-only for a continuous
        trait)."""
        unknown = [o for o in outputs if o not in _WRITE_OUTPUTS]
        if unknown:
            raise ValueError(f"unknown write outputs {unknown}; choose from {list(_WRITE_OUTPUTS)}")
        d = pathlib.Path(directory)
        d.mkdir(parents=True, exist_ok=True)
        if "values" in outputs:
            (d / "trait_values.tsv").write_text(_values_tsv(self.values))
        if "changes" in outputs:
            (d / "trait_changes.tsv").write_text(_changes_tsv(self.events))


def _fmt(v) -> str:
    """A trait value for TSV: a continuous float compactly, a discrete state label as-is."""
    return f"{v:.6g}" if isinstance(v, float) else str(v)


def _values_tsv(values: dict[int, object]) -> str:
    """The extant-tip values as a two-column ``node<TAB>trait`` table, one row per tip in id order.
    Tips are named ``n<id>`` to match the tree's Newick leaf labels."""
    rows = ["node\ttrait"]
    for i in sorted(values):
        rows.append(f"n{i}\t{_fmt(values[i])}")
    return "\n".join(rows) + "\n"


def _changes_tsv(changes: list[Change]) -> str:
    """The realized discrete transitions as ``time<TAB>lineage<TAB>from<TAB>to``, in time order."""
    rows = ["time\tlineage\tfrom\tto"]
    for c in changes:
        rows.append(f"{c.time:.6g}\tn{c.lineage}\t{c.from_state}\t{c.to_state}")
    return "\n".join(rows) + "\n"


def _preorder(tree: Tree) -> list[int]:
    """Node ids in an order that visits every node **after its parent** (a valid preorder). The
    forward engine always gives a child a higher id than its parent, so ascending id order suffices
    — the same monotonic-id fact ``genomes.prune`` relies on in reverse. No recursion needed."""
    return sorted(tree.nodes)


def _accrued_variance(rate, t0: float, t1: float, inherited: float = 1.0) -> float:
    """The variance a diffusing trait accrues over a branch spanning ``[t0, t1]`` — the integral
    ``∫ σ²(t) dt`` of the variance-rate. For a bare σ² this is ``σ²·(t1−t0)`` (Brownian motion); for a
    ``Time`` skyline (early burst) it sums σ² over each interval the branch crosses, stepping at the
    schedule's breakpoints. The same breakpoint walk the species/genome engines use — integrated over
    the branch rather than sampled at a point (σ² is piecewise-constant, so the integral is exact).

    ``inherited`` is the lineage's :class:`~zombi2.rates.modifiers.Inherited` factor (variable-rates
    BM), constant along the branch, threaded in by the caller and passed through to the rate; it
    factors straight out of the integral. A rate with no ``Inherited`` modifier ignores it."""
    total = 0.0
    t = t0
    while t < t1:
        nxt = min(rate.next_change(t), t1)  # constant rate → inf → one step of length (t1−t0)
        total += rate.effective(lineages=1, time=t, inherited=inherited) * (nxt - t)
        t = nxt
    return total


def simulate_continuous(tree, *, start=0.0, rate=1.0, reverts_to=None, pull=None,
                        seed=None) -> TraitsResult:
    """Evolve a continuous trait down a tree and return a :class:`TraitsResult`. One process, its
    variants selected by knobs (SPEC §4): **Brownian motion** (bare ``rate``), **Ornstein–Uhlenbeck**
    (add ``reverts_to`` + ``pull``), **early burst** (a ``Time`` skyline on ``rate``), and
    **variable-rates BM** (an ``Inherited`` modifier on ``rate``).

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
    early burst / ACDC — with the per-branch variance the exact integral ``∫ σ²(t) dt``; an
    ``Inherited(spread=…)`` modifier makes σ² **drift branch-to-branch** — variable-rates BM ("ClaDS
    for traits") — each lineage inheriting its parent's σ² times a lognormal kick drawn at the split.

    ``reverts_to`` (the optimum θ) and ``pull`` (the strength α > 0) turn the diffusion into
    Ornstein–Uhlenbeck — the value is pulled toward θ while it diffuses, the exact per-branch
    transition being ``Normal(θ + (x−θ)·e^{−α·dt}, σ²/(2α)·(1−e^{−2α·dt}))``. Give **both** or
    neither. OU with a *modified* σ² (early burst or variable rates on ``rate``) is not wired yet —
    use one or the other. Deterministic given ``seed``.
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
    # Time (early burst) and Inherited (variable-rates BM) are wired; other modifiers (Diversity) are
    # a later slice, so reject them loudly rather than silently ignore them — the genome engine's
    # discipline.
    for m in r.modifiers:
        if not isinstance(m, (Time, Inherited)):
            raise ValueError(
                f"rate carries {type(m).__name__}, which the continuous trait engine does not support "
                f"yet — Time (early burst) and Inherited (variable-rates BM) are wired; Diversity is a "
                f"later slice."
            )
    drifts = [m for m in r.modifiers if isinstance(m, Inherited)]
    if len(drifts) > 1:
        raise ValueError("rate carries more than one Inherited modifier; a variance-rate drifts one way")
    drift = drifts[0] if drifts else None  # the per-lineage σ² drift (variable-rates BM), or None

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
        if r.modifiers:
            raise ValueError(
                "a modified variance-rate (early burst via Time, or variable rates via Inherited) "
                "combined with OU (reverts_to / pull) is not wired yet — use one or the other."
            )
        theta, alpha = float(reverts_to), float(pull)
        sigma2 = r.effective(lineages=1)  # σ² is constant under OU (modifiers are rejected above)

    rng = np.random.default_rng(seed)
    node_values: dict[int, float] = {}
    inh: dict[int, float] = {}  # each lineage's σ² drift factor (variable-rates BM), constant per branch
    for i in _preorder(tree):
        node = tree.nodes[i]
        # the root starts from `start` at t=0; every other node from its parent's end value (parent
        # < i, already set). One uniform rule: node_values[i] is the trait at node i's end_time.
        x = float(start) if node.parent is None else node_values[node.parent]
        # thread the Inherited factor: the root's is 1.0, each daughter's is its parent's times a
        # lognormal kick drawn at the split (so σ² is autocorrelated down the tree). None ⇒ 1.0, no draw.
        if node.parent is None:
            inh[i] = drift.initial() if drift else 1.0
        else:
            inh[i] = drift.descend(inh[node.parent], rng) if drift else 1.0
        t0, t1 = node.birth_time, node.end_time
        if is_ou:
            e = math.exp(-alpha * (t1 - t0))       # mean-reversion toward θ over the branch
            mean = theta + (x - theta) * e
            var = sigma2 / (2.0 * alpha) * (1.0 - e * e)
        else:
            mean = x                                # pure diffusion (BM / early burst / variable-rates)
            var = _accrued_variance(r, t0, t1, inherited=inh[i])
        std = math.sqrt(var) if var > 0.0 else 0.0
        node_values[i] = mean + (float(rng.normal(0.0, std)) if std > 0.0 else 0.0)

    return TraitsResult(tree, node_values, seed)


# --- discrete traits: a state switching along the tree (Mk) ------------------------------------

def _q_matrix(states, switch) -> np.ndarray:
    """Build the ``k×k`` transition-rate matrix ``Q`` from ``switch`` — the CTMC generator whose
    off-diagonal ``Q[i, j] ≥ 0`` is the rate ``state i → state j``. ``switch`` is one of:

    - a **number** — the symmetric equal-rates shortcut: every ``i → j`` (``i ≠ j``) at that rate;
    - a ``{"from->to": rate}`` **dict** — only the named transitions, others zero (asymmetric);
    - a ``k×k`` **matrix** — the off-diagonal rates directly (the diagonal is ignored).

    Each diagonal is then set to minus its row sum, so rows sum to zero (a proper generator)."""
    k = len(states)
    idx = {s: i for i, s in enumerate(states)}
    Q = np.zeros((k, k))
    if isinstance(switch, bool):
        raise ValueError(f"switch must be a rate, not a bool, got {switch!r}")
    if isinstance(switch, (int, float)):
        if not math.isfinite(switch) or switch < 0:
            raise ValueError(f"switch rate must be finite and non-negative, got {switch!r}")
        Q[:] = float(switch)
        np.fill_diagonal(Q, 0.0)
    elif isinstance(switch, dict):
        for key, rate in switch.items():
            parts = [p.strip() for p in str(key).split("->")]
            if len(parts) != 2:
                raise ValueError(f"switch keys must read 'from->to', got {key!r}")
            frm, to = parts
            if frm not in idx or to not in idx:
                raise ValueError(f"switch key {key!r} names a state not in states={list(states)}")
            if frm == to:
                raise ValueError(f"switch key {key!r} is a self-transition; only i→j (i≠j) is a rate")
            if isinstance(rate, bool) or not isinstance(rate, (int, float)) \
                    or not math.isfinite(rate) or rate < 0:
                raise ValueError(f"switch rate for {key!r} must be finite and non-negative, got {rate!r}")
            Q[idx[frm], idx[to]] = float(rate)
    elif isinstance(switch, (list, tuple, np.ndarray)):
        arr = np.asarray(switch, dtype=float)
        if arr.shape != (k, k):
            raise ValueError(f"switch matrix must be {k}×{k} for {k} states, got shape {arr.shape}")
        Q = arr.copy()
        np.fill_diagonal(Q, 0.0)
        if np.any(Q < 0) or not np.all(np.isfinite(Q)):
            raise ValueError("switch matrix off-diagonals must be finite and non-negative")
    else:
        raise ValueError(
            "switch must be a number (symmetric rate), a {'from->to': rate} dict, or a k×k matrix"
        )
    np.fill_diagonal(Q, -Q.sum(axis=1))  # rows sum to zero
    return Q


def _gillespie(state: int, dt: float, Q: np.ndarray, rng) -> tuple[int, list]:
    """Exact CTMC simulation along a branch of duration ``dt`` from integer ``state`` (Gillespie).
    Returns ``(end_state, segments)`` where ``segments`` is a list of ``(state, duration)`` pieces
    summing to ``dt`` — the realized character history on this branch (a stochastic character map)."""
    k = Q.shape[0]
    segments: list[tuple[int, float]] = []
    elapsed = 0.0
    current = state
    while True:
        rate_out = -Q[current, current]
        if rate_out <= 0.0:  # absorbing state: no further jumps
            segments.append((current, dt - elapsed))
            return current, segments
        wait = float(rng.exponential(1.0 / rate_out))
        if elapsed + wait >= dt:  # the next jump falls past the branch end
            segments.append((current, dt - elapsed))
            return current, segments
        segments.append((current, wait))
        elapsed += wait
        probs = Q[current].copy()
        probs[current] = 0.0
        probs /= rate_out  # the embedded jump chain: where to, given a jump happened
        current = int(rng.choice(k, p=probs))


def simulate_discrete(tree, *, states, switch=None, start=None, liability=None, threshold=None,
                      seed=None) -> TraitsResult:
    """Evolve a discrete-state trait down a tree as a continuous-time Markov chain (the Mk model) and
    return a :class:`TraitsResult`. The jumps are simulated **exactly** by the Gillespie algorithm
    along every branch, so each node's ``(state, duration)`` segments *are* the realized history — a
    stochastic character map — and ``.events`` reads off the transitions.

    ``tree`` is the **complete** species tree (a :class:`~zombi2.species.Tree` or a
    :class:`~zombi2.species.SpeciesResult`); the trait evolves on every lineage, and ``.values`` reads
    the extant tips. ``states`` is the list of state labels (≥ 2, unique). ``switch`` gives the
    transition rates — a symmetric rate (``switch=0.1``), a ``{"marine->terrestrial": 0.1}`` dict of
    asymmetric rates, or a ``k×k`` matrix; see :func:`_q_matrix`. ``start`` is the root state (a label
    in ``states``); ``None`` draws one uniformly at random. As under convention B for continuous
    traits, the root evolves over its own branch, so ``node_values[root]`` is the state at the first
    split. Deterministic given ``seed``.

    ``liability`` / ``threshold`` (a discrete state read off an underlying continuous liability — the
    Wright–Felsenstein threshold model) are reserved but not wired yet; they arrive with the
    ``correlation=`` overlay, whose discrete case they underpin.
    """
    tree = tree.complete_tree if isinstance(tree, SpeciesResult) else tree
    if liability is not None or threshold is not None:
        raise ValueError(
            "threshold traits (liability= / threshold=) are a later slice (they arrive with the "
            "correlation overlay); give switch= for an Mk discrete-state trait."
        )
    states = list(states)
    if len(states) < 2:
        raise ValueError(f"a discrete trait needs at least 2 states, got {states!r}")
    if len(set(states)) != len(states):
        raise ValueError(f"states must be unique, got {states!r}")
    if switch is None:
        raise ValueError("give switch= — the transition rate(s) between the discrete states.")
    Q = _q_matrix(states, switch)

    rng = np.random.default_rng(seed)
    idx = {s: i for i, s in enumerate(states)}
    if start is None:
        start_i = int(rng.integers(len(states)))
    elif start in idx:
        start_i = idx[start]
    else:
        raise ValueError(
            f"start must be one of states={states} (or None for a uniform draw), got {start!r}"
        )

    node_values: dict[int, object] = {}
    history: dict[int, list] = {}
    for i in _preorder(tree):
        node = tree.nodes[i]
        # the root starts from `start` at t=0 and evolves over its own branch; every other node from
        # its parent's end state (parent < i, already set) — the same convention-B walk as continuous.
        cur = start_i if node.parent is None else idx[node_values[node.parent]]
        end_i, segs = _gillespie(cur, node.end_time - node.birth_time, Q, rng)
        node_values[i] = states[end_i]
        history[i] = [(states[s], d) for s, d in segs]

    return TraitsResult(tree, node_values, seed, kind="discrete", history=history)


__all__ = ["simulate_continuous", "simulate_discrete", "TraitsResult", "Change"]
