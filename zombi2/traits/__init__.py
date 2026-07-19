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
- **Early burst / ACDC**: give ``rate`` a ``OnTime`` skyline (``rate = σ² * mod.OnTime({0: 1, 5: 0.2})``)
  and σ² changes through time — the *same* ``OnTime`` modifier that gives the species tree its skyline.
  The per-branch variance is then the exact integral ``∫ σ²(t) dt`` over the branch.
- **Variable-rates BM** ("ClaDS for traits"): give ``rate`` an ``FromParent`` modifier
  (``rate = σ² * mod.FromParent(spread=0.3)``) and σ² drifts branch-to-branch — each lineage inherits
  its parent's σ² times a lognormal kick at the split — the *same* ``FromParent`` modifier that drifts
  the species rate (ClaDS) and the autocorrelated clock, one level over. (``reverts_to`` / ``pull`` are
  OU function arguments that revert the trait *value*, **not** a modifier — a rate modifier reverts a
  *rate*, which is the sequences level's CIR clock, a different mechanism.)
- **Diversity-dependent** (ecological limits): give ``rate`` a ``OnTotalDiversity`` modifier
  (``rate = σ² * mod.OnTotalDiversity(cap=100)``) and σ² slows as the clade fills — scaled by
  ``(1 − standing_diversity/cap)`` as the tree's lineages-through-time grows — the *same* ``OnTotalDiversity``
  modifier that slows species diversification, read here off the fixed tree (one-way, tree → trait).

``rate`` thus takes the whole modifier vocabulary — ``OnTime``, ``FromParent``, ``OnTotalDiversity`` — like any
other rate, and they compose (``σ² * OnTime({…}) * FromParent(spread=…)``).

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

import bisect
import math
import pathlib
from dataclasses import dataclass

import numpy as np

from ..rates.modifiers import OnTotalDiversity, FromParent, OnTime
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
    """The extant-tip values as a ``node<TAB>…`` table, one row per tip in id order (tips named
    ``n<id>`` to match the Newick). A single trait gives a ``node<TAB>trait`` table; correlated traits
    (per-node ``{trait: value}`` dicts) give one column per trait."""
    if values and isinstance(next(iter(values.values())), dict):  # correlated / multi-trait
        cols = list(next(iter(values.values())))
        rows = ["node\t" + "\t".join(str(c) for c in cols)]
        for i in sorted(values):
            rows.append(f"n{i}\t" + "\t".join(_fmt(values[i][c]) for c in cols))
        return "\n".join(rows) + "\n"
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


class _LTT:
    """The tree's lineages-through-time step function — how many lineages are alive at time ``t``
    (``birth ≤ t < end``), the *standing diversity* a :class:`~zombi2.rates.modifiers.OnTotalDiversity`
    modifier reads. Built once per run and used to integrate a diversity-dependent σ² over each
    branch, stepping at the tree's own speciation / extinction times (where the diversity changes)."""

    def __init__(self, tree: Tree) -> None:
        deltas: dict[float, int] = {}
        for n in tree.nodes.values():
            deltas[n.birth_time] = deltas.get(n.birth_time, 0) + 1  # a lineage starts
            deltas[n.end_time] = deltas.get(n.end_time, 0) - 1      # a lineage ends
        self._times = sorted(deltas)
        self._div: list[int] = []
        running = 0
        for t in self._times:
            running += deltas[t]
            self._div.append(running)  # standing diversity on the interval [times[k], times[k+1])

    def count(self, t: float) -> int:
        """The standing diversity at time ``t``."""
        k = bisect.bisect_right(self._times, t) - 1
        return self._div[k] if k >= 0 else 0

    def next_change(self, t: float) -> float:
        """The next time strictly after ``t`` at which the standing diversity changes, else ``inf``."""
        k = bisect.bisect_right(self._times, t)
        return self._times[k] if k < len(self._times) else math.inf


def _accrued_variance(rate, t0: float, t1: float, inherited: float = 1.0, ltt: "_LTT | None" = None) -> float:
    """The variance a diffusing trait accrues over a branch spanning ``[t0, t1]`` — the integral
    ``∫ σ²(t) dt`` of the variance-rate. For a bare σ² this is ``σ²·(t1−t0)`` (Brownian motion); for a
    ``OnTime`` skyline (early burst) it sums σ² over each interval the branch crosses, stepping at the
    schedule's breakpoints. The same breakpoint walk the species/genome engines use — integrated over
    the branch rather than sampled at a point (σ² is piecewise-constant, so the integral is exact).

    ``inherited`` is the lineage's :class:`~zombi2.rates.modifiers.FromParent` factor (variable-rates
    BM), constant along the branch, threaded in by the caller and passed through to the rate; it
    factors straight out of the integral. A rate with no ``FromParent`` modifier ignores it.

    ``ltt`` is the tree's lineages-through-time function when the rate carries a ``OnTotalDiversity`` modifier
    (diversity-dependent σ²): the integral then also steps at the tree's speciation / extinction times,
    reading the standing diversity on each sub-interval. ``None`` when σ² does not depend on diversity.
    (Stepping is O(events the branch crosses); fine for the trait level's one value per branch.)"""
    total = 0.0
    t = t0
    while t < t1:
        nxt = min(rate.next_change(t), t1)  # constant rate → inf → one step of length (t1−t0)
        div = 1.0
        if ltt is not None:                 # diversity-dependent σ²: also step where the LTT changes
            div = ltt.count(t)
            nxt = min(nxt, ltt.next_change(t))
        total += rate.effective(lineages=1, time=t, inherited=inherited, diversity=div) * (nxt - t)
        t = nxt
    return total


def _symmetric_sqrt(matrix: np.ndarray) -> np.ndarray:
    """A symmetric square root ``L`` of a symmetric PSD matrix (``L @ L == matrix``), via its
    eigendecomposition; tiny negative eigenvalues from round-off are clipped to zero."""
    w, V = np.linalg.eigh((matrix + matrix.T) / 2.0)
    return (V * np.sqrt(np.clip(w, 0.0, None))) @ V.T


def _correlation_matrix(traits: list, correlation) -> np.ndarray:
    """The ``k×k`` correlation matrix from pairwise ``{(a, b): ρ}`` — 1 on the diagonal, ρ off it
    (symmetric), 0 for unspecified pairs. Validates each ρ ∈ [−1, 1] and that the matrix is
    positive-semidefinite (the ρ values must be jointly consistent)."""
    idx = {t: i for i, t in enumerate(traits)}
    R = np.eye(len(traits))
    for pair, rho in (correlation or {}).items():
        if not (isinstance(pair, tuple) and len(pair) == 2):
            raise ValueError(f"correlation keys must be (trait_a, trait_b) pairs, got {pair!r}")
        a, b = pair
        if a not in idx or b not in idx:
            raise ValueError(f"correlation key {pair!r} names a trait not in {traits}")
        if a == b:
            raise ValueError(f"correlation key {pair!r} is a self-correlation")
        if isinstance(rho, bool) or not isinstance(rho, (int, float)) or not -1.0 <= rho <= 1.0:
            raise ValueError(f"correlation for {pair!r} must be a number in [−1, 1], got {rho!r}")
        R[idx[a], idx[b]] = R[idx[b], idx[a]] = float(rho)
    if float(np.linalg.eigvalsh(R).min()) < -1e-9:
        raise ValueError(
            "the correlation matrix is not positive-semidefinite — the given ρ values are jointly "
            "inconsistent (e.g. three traits cannot all be strongly negatively correlated)."
        )
    return R


def _simulate_correlated(tree, start, rate, reverts_to, pull, correlation, seed) -> TraitsResult:
    """Correlated continuous traits in **one call** (the joint rule inside a level). ``start`` and
    ``rate`` are dicts over the same trait names; the branch increment is drawn from ``MVN(0, Σ·dt)``
    with ``Σ = D R D`` (``D = diag(σ_i)``, ``R`` from ``correlation``), so at a tip the correlation
    between two traits equals their ρ. Correlated Brownian motion this slice — per-trait modifiers and
    multivariate OU are later slices. Convention B holds: the root diffuses over its own branch."""
    if not isinstance(start, dict) or not isinstance(rate, dict):
        raise ValueError(
            "for correlated traits give both start and rate as dicts keyed by trait name, e.g. "
            "start={'size': 0.0, 'limb': 0.0}, rate={'size': 1.0, 'limb': 0.8}."
        )
    traits = list(start)
    if set(rate) != set(traits):
        raise ValueError(
            f"start and rate must name the same traits; got {sorted(start)} vs {sorted(rate)}"
        )
    if len(traits) < 2:
        raise ValueError("correlated traits need ≥ 2 traits; one trait is a plain simulate_continuous call")
    if reverts_to is not None or pull is not None:
        raise ValueError("multivariate OU (reverts_to / pull with correlated traits) is not wired yet.")

    sigma2 = np.empty(len(traits))
    for i, name in enumerate(traits):
        if isinstance(start[name], bool) or not isinstance(start[name], (int, float)) \
                or not math.isfinite(start[name]):
            raise ValueError(f"start[{name!r}] must be a finite number, got {start[name]!r}")
        r = as_rate(rate[name], default_scope=PerLineage)
        if not isinstance(r.scope, PerLineage):
            raise ValueError(f"rate[{name!r}] must be per lineage — drop the scope wrapper.")
        if r.modifiers:
            raise ValueError(
                f"rate[{name!r}] carries a modifier; per-trait modifiers combined with correlation are "
                f"a later slice — use bare per-trait rates here."
            )
        sigma2[i] = r.effective(lineages=1)

    R = _correlation_matrix(traits, correlation)
    sd = np.sqrt(sigma2)
    sigma = _symmetric_sqrt((sd[:, None] * R) * sd[None, :])  # sqrt of Σ = D R D (Σ_ij = σ_i σ_j ρ_ij)
    start_vec = np.array([float(start[t]) for t in traits])
    k = len(traits)

    rng = np.random.default_rng(seed)
    node_values: dict[int, dict] = {}
    for i in _preorder(tree):
        node = tree.nodes[i]
        x = start_vec if node.parent is None else np.array([node_values[node.parent][t] for t in traits])
        dt = node.end_time - node.birth_time
        vec = x + (math.sqrt(dt) * (sigma @ rng.standard_normal(k)) if dt > 0.0 else 0.0)
        node_values[i] = {t: float(vec[j]) for j, t in enumerate(traits)}

    return TraitsResult(tree, node_values, seed)


def simulate_continuous(tree, *, start=0.0, rate=1.0, reverts_to=None, pull=None,
                        correlation=None, seed=None) -> TraitsResult:
    """Evolve a continuous trait down a tree and return a :class:`TraitsResult`. One process, its
    variants selected by knobs (SPEC §4): **Brownian motion** (bare ``rate``), **Ornstein–Uhlenbeck**
    (add ``reverts_to`` + ``pull``), **early burst** (a ``OnTime`` skyline on ``rate``), and
    **variable-rates BM** (an ``FromParent`` modifier on ``rate``).

    **Correlated traits** ride together in **one call** (the joint rule inside a level): pass
    ``start`` and ``rate`` as dicts keyed by trait name and a ``correlation={(a, b): ρ}`` overlay
    (each ρ ∈ [−1, 1]). The traits then diffuse jointly — the branch increment is drawn from
    ``MVN(0, Σ·dt)`` with ``Σ = D R D`` (``D = diag(σ_i)``, ``R`` the correlation matrix), so at a tip
    the correlation between two traits is exactly their ρ. This slice wires correlated **Brownian
    motion** (bare per-trait rates); per-trait modifiers and multivariate OU are later slices.

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
    motion (``Normal(0, σ²·dt)`` over a branch); a ``OnTime`` modifier makes σ² change through time —
    early burst / ACDC — with the per-branch variance the exact integral ``∫ σ²(t) dt``; an
    ``FromParent(spread=…)`` modifier makes σ² **drift branch-to-branch** — variable-rates BM ("ClaDS
    for traits") — each lineage inheriting its parent's σ² times a lognormal kick drawn at the split;
    a ``OnTotalDiversity(cap=…)`` modifier makes σ² **slow as the clade fills up** — diversity-dependent /
    ecological-limits trait evolution — σ² scaled by ``(1 − standing_diversity/cap)`` as the tree's
    lineages-through-time grows (the tree is a fixed input the trait reads).

    ``reverts_to`` (the optimum θ) and ``pull`` (the strength α > 0) turn the diffusion into
    Ornstein–Uhlenbeck — the value is pulled toward θ while it diffuses, the exact per-branch
    transition being ``Normal(θ + (x−θ)·e^{−α·dt}, σ²/(2α)·(1−e^{−2α·dt}))``. Give **both** or
    neither. OU with a *modified* σ² (early burst or variable rates on ``rate``) is not wired yet —
    use one or the other. Deterministic given ``seed``.
    """
    tree = tree.complete_tree if isinstance(tree, SpeciesResult) else tree
    if isinstance(start, dict) or isinstance(rate, dict) or correlation is not None:
        return _simulate_correlated(tree, start, rate, reverts_to, pull, correlation, seed)
    if isinstance(start, bool) or not isinstance(start, (int, float)) or not math.isfinite(start):
        raise ValueError(f"start must be a finite number, got {start!r}")
    r = as_rate(rate, default_scope=PerLineage)
    if not isinstance(r.scope, PerLineage):
        raise ValueError(
            f"rate has a {type(r.scope).__name__} scope, but a continuous trait's variance-rate is "
            f"per lineage — drop the scope wrapper (per lineage is the default)."
        )
    # OnTime (early burst), FromParent (variable-rates BM), and OnTotalDiversity (diversity-dependent) are the
    # wired σ² modifiers; anything else is rejected loudly — the genome engine's discipline.
    for m in r.modifiers:
        if not isinstance(m, (OnTime, FromParent, OnTotalDiversity)):
            raise ValueError(
                f"rate carries {type(m).__name__}, which the continuous trait engine does not support "
                f"— OnTime (early burst), FromParent (variable-rates BM), and OnTotalDiversity "
                f"(diversity-dependent) are wired."
            )
    drifts = [m for m in r.modifiers if isinstance(m, FromParent)]
    if len(drifts) > 1:
        raise ValueError("rate carries more than one FromParent modifier; a variance-rate drifts one way")
    drift = drifts[0] if drifts else None  # the per-lineage σ² drift (variable-rates BM), or None
    has_diversity = any(isinstance(m, OnTotalDiversity) for m in r.modifiers)  # σ² reads the standing LTT

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
                "a modified variance-rate (early burst via OnTime, variable rates via FromParent, or "
                "diversity-dependence via OnTotalDiversity) combined with OU (reverts_to / pull) is not "
                "wired yet — use one or the other."
            )
        theta, alpha = float(reverts_to), float(pull)
        sigma2 = r.effective(lineages=1)  # σ² is constant under OU (modifiers are rejected above)

    rng = np.random.default_rng(seed)
    ltt = _LTT(tree) if has_diversity else None  # the standing-diversity curve, when σ² reads it
    node_values: dict[int, float] = {}
    inh: dict[int, float] = {}  # each lineage's σ² drift factor (variable-rates BM), constant per branch
    for i in _preorder(tree):
        node = tree.nodes[i]
        # the root starts from `start` at t=0; every other node from its parent's end value (parent
        # < i, already set). One uniform rule: node_values[i] is the trait at node i's end_time.
        x = float(start) if node.parent is None else node_values[node.parent]
        # thread the inherited factor: the root's is 1.0, each daughter's is its parent's times a
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
            var = _accrued_variance(r, t0, t1, inherited=inh[i], ltt=ltt)
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


def _finite(x, name: str) -> float:
    """Coerce ``x`` to a finite float or raise a clear ``ValueError`` naming ``name``."""
    if isinstance(x, bool) or not isinstance(x, (int, float)) or not math.isfinite(x):
        raise ValueError(f"{name} must be a finite number, got {x!r}")
    return float(x)


def _liability_sigma2(spec, name) -> float:
    """The liability's variance-rate σ² — a bare per-lineage rate this slice (no modifiers)."""
    r = as_rate(spec, default_scope=PerLineage)
    where = "liability" if name is None else f"liability[{name!r}]"
    if not isinstance(r.scope, PerLineage) or r.modifiers:
        raise ValueError(f"{where} must be a bare variance-rate (per lineage, no modifiers) this slice.")
    return r.effective(lineages=1)


def _threshold_start_vec(start, traits: list) -> np.ndarray:
    """The starting liability vector: ``None`` → zeros, a scalar → the same for every trait, a dict →
    per-trait (missing traits default to 0.0)."""
    if start is None:
        return np.zeros(len(traits))
    if isinstance(start, dict):
        extra = set(start) - set(traits)
        if extra:
            raise ValueError(f"start names traits not in the liabilities: {sorted(extra)}")
        return np.array([_finite(start.get(t, 0.0), f"start[{t!r}]") for t in traits])
    return np.full(len(traits), _finite(start, "start"))


def _threshold_cuts(states: list, threshold) -> np.ndarray:
    """The ``k−1`` strictly-increasing cut points for ``k`` states (a scalar is allowed for ``k=2``)."""
    if threshold is None:
        raise ValueError("a threshold trait needs threshold= — a cut point (k−1 increasing cuts for k states)")
    thr = np.atleast_1d(np.asarray(threshold, dtype=float))
    if thr.ndim != 1 or len(thr) != len(states) - 1:
        raise ValueError(
            f"threshold must give {len(states) - 1} cut point(s) for {len(states)} states, got {threshold!r}"
        )
    if len(thr) > 1 and bool(np.any(np.diff(thr) <= 0)):
        raise ValueError(f"threshold cut points must be strictly increasing, got {threshold!r}")
    return thr


def _simulate_threshold(tree, states, liability, threshold, start, correlation, seed) -> TraitsResult:
    """The Wright–Felsenstein **threshold** model — a discrete state read off an underlying continuous
    Brownian liability: the liability diffuses (convention B), and the observed state is which
    ``threshold``-cut interval it lands in. With ``liability`` a dict + ``correlation``, several traits'
    liabilities diffuse **jointly** (correlated discrete traits, one call), each cut by the shared
    thresholds. There is no Gillespie stochastic map — the discreteness is a cut through a continuous
    path — so ``.history`` is ``None`` and ``.events`` is empty."""
    thr = _threshold_cuts(states, threshold)

    def label(x):
        return states[int(np.searchsorted(thr, x))]

    rng = np.random.default_rng(seed)
    node_values: dict[int, object] = {}

    if isinstance(liability, dict):  # correlated discrete traits — joint liabilities, shared thresholds
        traits = list(liability)
        if len(traits) < 2:
            raise ValueError("a correlated threshold trait needs ≥ 2 liabilities; one is a plain threshold trait")
        sigma2 = np.array([_liability_sigma2(liability[t], t) for t in traits])
        R = _correlation_matrix(traits, correlation)
        sd = np.sqrt(sigma2)
        chol = _symmetric_sqrt((sd[:, None] * R) * sd[None, :])
        start_vec = _threshold_start_vec(start, traits)
        k = len(traits)
        liab: dict[int, np.ndarray] = {}
        for i in _preorder(tree):
            node = tree.nodes[i]
            x = start_vec if node.parent is None else liab[node.parent]
            dt = node.end_time - node.birth_time
            liab[i] = x + (math.sqrt(dt) * (chol @ rng.standard_normal(k)) if dt > 0.0 else 0.0)
            node_values[i] = {t: label(liab[i][j]) for j, t in enumerate(traits)}
    else:  # a single threshold trait
        if correlation is not None:
            raise ValueError("correlation= needs a dict liability= (one liability per trait)")
        sig2 = _liability_sigma2(liability, None)
        start_liab = 0.0 if start is None else _finite(start, "start")
        liab_s: dict[int, float] = {}
        for i in _preorder(tree):
            node = tree.nodes[i]
            x = start_liab if node.parent is None else liab_s[node.parent]
            dt = node.end_time - node.birth_time
            std = math.sqrt(sig2 * dt) if sig2 > 0.0 and dt > 0.0 else 0.0
            liab_s[i] = x + (float(rng.normal(0.0, std)) if std > 0.0 else 0.0)
            node_values[i] = label(liab_s[i])

    return TraitsResult(tree, node_values, seed, kind="discrete", history=None)


def simulate_discrete(tree, *, states, switch=None, start=None, liability=None, threshold=None,
                      correlation=None, seed=None) -> TraitsResult:
    """Evolve a discrete-state trait down a tree and return a :class:`TraitsResult`. Two mechanisms:

    - **Mk** (``switch=``) — a continuous-time Markov chain over the ``states``, simulated **exactly**
      by Gillespie along every branch, so each node's ``(state, duration)`` segments *are* the realized
      history (``.history``) and ``.events`` reads off the transitions. ``switch`` is a symmetric rate
      (``0.1``), a ``{"marine->terrestrial": 0.1}`` dict, or a ``k×k`` matrix (see :func:`_q_matrix`).
      ``start`` is the root state (a label in ``states``; ``None`` draws one uniformly).
    - **Threshold** (``liability=`` + ``threshold=``) — the Wright–Felsenstein model: a discrete state
      read off an underlying continuous Brownian **liability** (variance-rate ``liability``), cut into
      ``states`` by the ``threshold`` cut point(s) (``k−1`` increasing cuts for ``k`` states). ``start``
      is the initial *liability* (a number, default 0.0). Give ``liability`` as a dict + a
      ``correlation={(a, b): ρ}`` overlay to evolve **correlated** discrete traits jointly — their
      liabilities diffuse together (``Σ = D R D``) and each is cut by the shared thresholds. A threshold
      trait has no Gillespie map, so ``.history`` is ``None`` and ``.events`` empty.

    ``tree`` is the **complete** species tree (a :class:`~zombi2.species.Tree` or
    :class:`~zombi2.species.SpeciesResult`); the trait evolves on every lineage (convention B: the root
    diffuses over its own branch), and ``.values`` reads the extant tips. Deterministic given ``seed``.
    """
    tree = tree.complete_tree if isinstance(tree, SpeciesResult) else tree
    states = list(states)
    if len(states) < 2:
        raise ValueError(f"a discrete trait needs at least 2 states, got {states!r}")
    if len(set(states)) != len(states):
        raise ValueError(f"states must be unique, got {states!r}")
    if liability is not None or threshold is not None:
        if switch is not None:
            raise ValueError("give switch= (an Mk trait) OR liability=/threshold= (a threshold trait), not both")
        return _simulate_threshold(tree, states, liability, threshold, start, correlation, seed)
    if correlation is not None:
        raise ValueError("correlation= on a discrete trait needs the threshold model — give liability= and threshold=")
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
