"""Traits — the continuous (diffusion) engine: Brownian motion and its OU / early-burst / variable-rate / diversity-dependent / driven variants (simulate_continuous)."""

from __future__ import annotations

import bisect
import math

import numpy as np

from ..rates.modifiers import ByFamily, DrivenBy, FromParent, OnTime, OnTotalDiversity
from ..rates.rate import as_rate
from ..rates.scope import PerLineage
from ..tree import Tree, as_tree

from ._shared import _correlation_matrix, _driven_mods, _preorder, _resolve_drivers, _symmetric_sqrt
from .result import Change, TraitsResult

WIRED_MODIFIERS = (OnTime, FromParent, OnTotalDiversity, DrivenBy)  #: the modifiers a continuous rate takes

class _LTT:
    """The tree's lineages-through-time step function — how many lineages are alive at time ``t``
    (``birth ≤ t < end``), the *standing diversity* a `OnTotalDiversity`
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



def _accrued_variance(rate, t0: float, t1: float, inherited: float = 1.0, ltt: "_LTT | None" = None,
                      trajs: dict | None = None, node_id: int | None = None) -> float:
    """The variance a diffusing trait accrues over a branch spanning ``[t0, t1]`` — the integral
    ``∫ σ²(t) dt`` of the variance-rate. For a bare σ² this is ``σ²·(t1−t0)`` (Brownian motion); for a
    ``OnTime`` skyline (early burst) it sums σ² over each interval the branch crosses, stepping at the
    schedule's breakpoints. The same breakpoint walk the species/genome engines use — integrated over
    the branch rather than sampled at a point (σ² is piecewise-constant, so the integral is exact).

    ``inherited`` is the lineage's `FromParent` factor (variable-rates
    BM), constant along the branch, threaded in by the caller and passed through to the rate; it
    factors straight out of the integral. A rate with no ``FromParent`` modifier ignores it.

    ``ltt`` is the tree's lineages-through-time function when the rate carries a ``OnTotalDiversity`` modifier
    (diversity-dependent σ²): the integral then also steps at the tree's speciation / extinction times,
    reading the standing diversity on each sub-interval. ``None`` when σ² does not depend on diversity.

    ``trajs`` (with ``node_id``, the lineage this branch is) are the driver trajectories when the rate
    carries a `DrivenBy` — σ² read off another trait. The driver's
    value on this lineage is threaded in as ``drivers``, and the integral **steps where the driver
    switches** (``next_change``) exactly as it steps at a skyline breakpoint: a discrete driver
    switches *mid-branch*, so a single sample per branch would credit the whole branch to whichever
    state happened to be read. ``None`` when σ² is not driven.
    (Stepping is O(events the branch crosses); fine for the trait level's one value per branch.)"""
    total = 0.0
    t = t0
    while t < t1:
        nxt = min(rate.next_change(t), t1)  # constant rate → inf → one step of length (t1−t0)
        div = 1.0
        if ltt is not None:                 # diversity-dependent σ²: also step where the LTT changes
            div = ltt.count(t)
            nxt = min(nxt, ltt.next_change(t))
        drivers = None
        if trajs:                           # driven σ²: also step where this lineage's driver switches
            drivers = {key: traj.value(node_id, t) for key, traj in trajs.items()}
            for traj in trajs.values():
                nxt = min(nxt, traj.next_change(node_id, t))
        total += rate.effective(lineages=1, time=t, inherited=inherited, diversity=div,
                                drivers=drivers) * (nxt - t)
        t = nxt
    return total



def _at_speciation_jump_sd(at_speciation) -> float:
    """The on-speciation jump width (√variance) from ``at_speciation`` — ``0.0`` if not requested. A
    jump of ``Normal(0, at_speciation)`` is added to each daughter's value at every speciation (the
    punctuational / speciational mode); it *reads* the tree it rides on without feeding back into it,
    so it is a trait-level option, not a joint model (SPEC §4)."""
    if at_speciation is None:
        return 0.0
    if isinstance(at_speciation, bool) or not isinstance(at_speciation, (int, float)) \
            or not math.isfinite(at_speciation) or at_speciation < 0:
        raise ValueError(
            f"at_speciation must be a non-negative number (the jump variance), got {at_speciation!r}"
        )
    return math.sqrt(at_speciation)



def _simulate_regimes(tree, start, rate, reverts_to, pull, regimes, at_speciation, seed,
                      progress=False) -> TraitsResult:
    """Multi-optimum OU — the optimum shifts by **regime**, a discrete stochastic map painted on the
    *same* tree (typically a `simulate_discrete()` run). Along each branch the value follows OU
    toward the current regime's optimum, integrated **exactly** across the regime's ``(state,
    duration)`` segments (a regime may switch part-way along a branch); convention B paints the root
    branch too. ``reverts_to`` is a dict ``{regime_state: θ}``; ``pull`` (α > 0) and ``rate`` (σ²) are
    shared across regimes (the OUM variant)."""
    if getattr(regimes, "kind", None) != "discrete" or regimes.history is None:
        raise ValueError(
            "regimes must be a discrete TraitsResult carrying a stochastic map — paint them with "
            "simulate_discrete(...) on this same tree."
        )
    if set(regimes.history) != set(tree.nodes):
        raise ValueError("regimes must be painted on the SAME tree this trait rides (node ids differ).")
    if at_speciation is not None:
        raise ValueError("at_speciation combined with regimes is not implemented yet.")
    if isinstance(start, dict) or isinstance(rate, dict):
        raise ValueError("regimes is a single-trait model; give scalar start and rate.")
    if not isinstance(reverts_to, dict):
        raise ValueError("with regimes, reverts_to is a dict {regime_state: optimum θ}, one per regime.")
    if isinstance(pull, bool) or not isinstance(pull, (int, float)) or not math.isfinite(pull) or pull is None \
            or pull <= 0:
        raise ValueError(f"regimes needs pull (the OU strength α > 0), got {pull!r}")
    if isinstance(start, bool) or not isinstance(start, (int, float)) or not math.isfinite(start):
        raise ValueError(f"start must be a finite number, got {start!r}")
    r = as_rate(rate, default_scope=PerLineage)
    if not isinstance(r.scope, PerLineage) or r.modifiers:
        raise ValueError("rate must be a bare variance-rate (per lineage, no modifiers) with regimes this slice.")
    sigma2, alpha = r.effective(lineages=1), float(pull)
    used = {seg[0] for segs in regimes.history.values() for seg in segs}
    missing = used - set(reverts_to)
    if missing:
        raise ValueError(f"reverts_to is missing an optimum for regime state(s) {sorted(missing)}")

    rng = np.random.default_rng(seed)
    node_values: dict[int, float] = {}
    for i in _preorder(tree, progress):
        node = tree.nodes[i]
        x = float(start) if node.parent is None else node_values[node.parent]
        for regime, dt in regimes.history[i]:  # integrate OU piece-by-piece across the regime segments
            theta = float(reverts_to[regime])
            e = math.exp(-alpha * dt)
            var = sigma2 / (2.0 * alpha) * (1.0 - e * e)
            x = theta + (x - theta) * e + (float(rng.normal(0.0, math.sqrt(var))) if var > 0.0 else 0.0)
        node_values[i] = x
    return TraitsResult(tree, node_values, [], seed)  # correlated / regimes: no along-branch event log



def _simulate_correlated(tree, start, rate, reverts_to, pull, correlation, at_speciation, seed,
                         progress=False) -> TraitsResult:
    """Correlated continuous traits in **one call** (the joint rule inside a level). ``start`` and
    ``rate`` are dicts over the same trait names; the branch increment is drawn from ``MVN(0, Σ·dt)``
    with ``Σ = D R D`` (``D = diag(σ_i)``, ``R`` from ``correlation``), so at a tip the correlation
    between two traits equals their ρ. Correlated Brownian motion only. Convention B holds: the root diffuses over its own branch."""
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
        raise ValueError("multivariate OU (reverts_to / pull with correlated traits) is not implemented yet.")
    if at_speciation is not None:
        raise ValueError("at_speciation (on-speciation jumps) with correlated traits is not implemented yet.")

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
    for i in _preorder(tree, progress):
        node = tree.nodes[i]
        x = start_vec if node.parent is None else np.array([node_values[node.parent][t] for t in traits])
        dt = node.end_time - node.birth_time
        vec = x + (math.sqrt(dt) * (sigma @ rng.standard_normal(k)) if dt > 0.0 else 0.0)
        node_values[i] = {t: float(vec[j]) for j, t in enumerate(traits)}

    return TraitsResult(tree, node_values, [], seed)  # correlated / regimes: no along-branch event log



def simulate_continuous(tree, *, start=0.0, rate=1.0, reverts_to=None, pull=None,
                        correlation=None, at_speciation=None, regimes=None, seed=None,
                        progress=False) -> TraitsResult:
    """Evolve a continuous trait down a tree and return a `TraitsResult`. One process, its
    variants selected by knobs (SPEC §4): **Brownian motion** (bare ``rate``), **Ornstein–Uhlenbeck**
    (add ``reverts_to`` + ``pull``), **early burst** (a ``OnTime`` skyline on ``rate``), and
    **variable-rates BM** (an ``FromParent`` modifier on ``rate``).

    **Correlated traits** ride together in **one call** (the joint rule inside a level): pass
    ``start`` and ``rate`` as dicts keyed by trait name and a ``correlation={(a, b): ρ}`` overlay
    (each ρ ∈ [−1, 1]). The traits then diffuse jointly — the branch increment is drawn from
    ``MVN(0, Σ·dt)`` with ``Σ = D R D`` (``D = diag(σ_i)``, ``R`` the correlation matrix), so at a tip
    the correlation between two traits is exactly their ρ. Correlated **Brownian motion** only, with bare per-trait rates.

    ``tree`` is the **complete** species tree (a `Tree`, or a
    `SpeciesResult` whose ``complete_tree`` is used). The trait evolves on
    **every** lineage, extant and extinct alike, so the ancestral states are exact and complete; the
    observed dataset is the extant tips, ``result.values``.

    ``start`` is the value at ``t = 0`` (the origin, ``root.birth_time``): the root lineage
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
    lineages-through-time grows (the tree is a fixed input the trait reads); a
    ``DrivenBy(driver, {…})`` modifier makes σ² **read another trait** — the driver grown first on this
    same tree and handed over as its result object or its written ``trait_events.tsv``, so a lineage
    diffuses faster while the driver is in one state than another. A discrete driver switches
    *mid-branch*, and the per-branch variance is the integral across those pieces, so a branch that
    spends half its length in the fast state accrues exactly half the fast variance.

    ``reverts_to`` (the optimum θ) and ``pull`` (the strength α > 0) turn the diffusion into
    Ornstein–Uhlenbeck — the value is pulled toward θ while it diffuses, the exact per-branch
    transition being ``Normal(θ + (x−θ)·e^{−α·dt}, σ²/(2α)·(1−e^{−2α·dt}))``. Give **both** or
    neither. That transition law is exact only for a σ² that is constant along the branch, so OU with a
    *modified* σ² — a ``OnTime``, ``FromParent``, ``OnTotalDiversity`` or ``DrivenBy`` on ``rate`` — is
    not implemented yet; use one or the other.

    ``at_speciation`` adds an **on-speciation** jump — ``Normal(0, at_speciation)`` on each daughter at
    every speciation (the punctuational mode), layered on top of the along-branch anagenesis.
    ``regimes`` gives **multi-optimum OU**: pass a discrete `TraitsResult` (a stochastic map
    painted by `simulate_discrete()` on this same tree) and a per-regime ``reverts_to={regime: θ}``,
    and the value follows OU toward whichever regime's optimum a branch is in. Deterministic given
    ``seed``.
    """
    tree = as_tree(tree, level="traits")
    if regimes is not None:
        return _simulate_regimes(tree, start, rate, reverts_to, pull, regimes, at_speciation, seed,
                                 progress)
    if isinstance(start, dict) or isinstance(rate, dict) or correlation is not None:
        return _simulate_correlated(tree, start, rate, reverts_to, pull, correlation,
                                    at_speciation, seed, progress)
    if isinstance(start, bool) or not isinstance(start, (int, float)) or not math.isfinite(start):
        raise ValueError(f"start must be a finite number, got {start!r}")
    r = as_rate(rate, default_scope=PerLineage)
    if not isinstance(r.scope, PerLineage):
        raise ValueError(
            f"rate has a {type(r.scope).__name__} scope, but a continuous trait's variance-rate is "
            f"per lineage — drop the scope wrapper (per lineage is the default)."
        )
    # OnTime (early burst), FromParent (variable-rates BM), OnTotalDiversity (diversity-dependent)
    # and DrivenBy (σ² driven by another trait) are the σ² modifiers this engine supports; anything
    # else is rejected loudly — the genome engine's discipline.
    for m in r.modifiers:
        if isinstance(m, ByFamily):
            # not a missing feature: there is nothing here for it to mean
            raise ValueError(
                "rate carries ByFamily, but a trait has no gene families — ByFamily belongs on a "
                "genomes rate. For per-lineage heterogeneity here use FromParent (variable-rates BM)."
            )
        if not isinstance(m, WIRED_MODIFIERS):
            raise ValueError(
                f"rate carries {type(m).__name__}, which the continuous trait engine does not "
                f"support. It takes OnTime (early burst), FromParent (variable-rates BM), "
                f"OnTotalDiversity (diversity-dependent), and DrivenBy (driven by another trait)."
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
            # a driven σ² is a modified σ², so it falls under this guard like the other three — and
            # the message names the one actually given rather than making the reader match it against
            # a list. The OU transition law is exact only for a σ² that is constant along the branch.
            carried = ", ".join(dict.fromkeys(type(m).__name__ for m in r.modifiers))
            raise ValueError(
                f"rate carries {carried}, and a modified variance-rate combined with OU (reverts_to / "
                f"pull) is not implemented yet — use one or the other."
            )
        theta, alpha = float(reverts_to), float(pull)
        sigma2 = r.effective(lineages=1)  # σ² is constant under OU (modifiers are rejected above)

    # conditioning: a σ² carrying DrivenBy reads another trait, grown first on this same tree. Resolve
    # each driver once into a trajectory (value + next-switch, keyed by the shared node id), from a
    # written trait log or a grown result handed over in memory. Undriven ⇒ empty, and the walk below
    # is exactly the walk it was — no driver, no lookup, no change to the draw order.
    trajs = _resolve_drivers(_driven_mods(r), tree)

    jump_sd = _at_speciation_jump_sd(at_speciation)  # on-speciation jump width (0 if not requested)

    rng = np.random.default_rng(seed)
    ltt = _LTT(tree) if has_diversity else None  # the standing-diversity curve, when σ² reads it
    node_values: dict[int, float] = {}
    root = tree.nodes[tree.root]
    # the initial value at t=0 — the origin the log reconstructs from (SPEC §2). A diffusion cannot be
    # rebuilt from events, but the row keeps the file's shape uniform across trait kinds.
    events: list[Change] = [Change(root.birth_time, "initial", tree.root, None, float(start))]
    inh: dict[int, float] = {}  # each lineage's σ² drift factor (variable-rates BM), constant per branch
    for i in _preorder(tree, progress):
        node = tree.nodes[i]
        # the root starts from `start` at t=0; every other node from its parent's end value (parent
        # < i, already set). One uniform rule: node_values[i] is the trait at node i's end_time.
        x = float(start) if node.parent is None else node_values[node.parent]
        if node.parent is not None and jump_sd > 0.0:
            jumped = x + float(rng.normal(0.0, jump_sd))  # on speciation: a jump at the split…
            events.append(Change(node.birth_time, "on_speciation", i, x, jumped))
            x = jumped                                    # …then anagenesis along the branch
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
            var = _accrued_variance(r, t0, t1, inherited=inh[i], ltt=ltt, trajs=trajs, node_id=i)
        std = math.sqrt(var) if var > 0.0 else 0.0
        node_values[i] = mean + (float(rng.normal(0.0, std)) if std > 0.0 else 0.0)

    return TraitsResult(tree, node_values, events, seed)


# --- discrete traits: a state switching along the tree (Mk) ------------------------------------

