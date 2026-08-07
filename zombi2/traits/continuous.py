"""Traits — the continuous (diffusion) engine: Brownian motion and its OU / early-burst / variable-rate / diversity-dependent / driven variants (simulate_continuous)."""

from __future__ import annotations

import bisect
import math
from typing import cast

import numpy as np

from ..rates.mapping import check_not_a_kernel
from ..rng import stream
from ..rates.modifiers import (describe, DRAWN, INHERITED, DrivenBy, OnTime, OnTotalDiversity,
                               carried_at_birth, carried_at_split, check_one_memory,
                               is_implemented, product)
from ..rates.rate import as_rate
from ..rates.scope import PerLineage
from ..tree import Tree, as_tree

from ._shared import _correlation_matrix, _driven_mods, _preorder, _resolve_drivers, _symmetric_sqrt
from .result import Change, TraitsResult

IMPLEMENTED_MODIFIERS = (OnTime, (INHERITED, "lineage"), OnTotalDiversity, DrivenBy)  #: the cells a continuous rate takes

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
                      trajs: dict | None = None, node_id: int | None = None, *,
                      pull: float = 0.0) -> float:
    """The variance a diffusing trait accrues over a branch spanning ``[t0, t1]`` — the integral
    ``∫ σ²(t) dt`` of the variance-rate. For a bare σ² this is ``σ²·(t1−t0)`` (Brownian motion); for a
    ``OnTime`` skyline (early burst) it sums σ² over each interval the branch crosses, stepping at the
    schedule's breakpoints. The same breakpoint walk the species/genome engines use — integrated over
    the branch rather than sampled at a point (σ² is piecewise-constant, so the integral is exact).

    ``pull`` is the OU strength α, and it changes *which* integral this is. Under
    ``dX = −α(X−θ)dt + σ(t)dW`` the noise injected at time ``s`` has already decayed by
    ``e^{−α(t1−s)}`` when the branch ends, so it enters the end-of-branch variance squared:

    ====================  ==================================================
    ``pull == 0`` (BM)    ``∫_{t0}^{t1} σ²(s) ds``
    ``pull == α`` (OU)    ``∫_{t0}^{t1} e^{−2α(t1−s)} σ²(s) ds``
    ====================  ==================================================

    These are **not** the same number and one is not a small correction on the other: on a branch of
    length 2 with α = 1.3 the weighted integral can be an order of magnitude smaller, because
    everything the trait accrued early has been pulled back toward θ by the time the branch ends.
    Reusing the Brownian integral under OU would give a run that still looks like OU — mean-reverting,
    normally distributed — while being the wrong process, which is why the weight lives here rather
    than in the caller. On the interval ``[a, b)`` where σ² is constant the weighted contribution is
    exactly ``σ²·e^{−2α(t1−b)}·(1−e^{−2α(b−a)})/(2α)``, written below with ``expm1`` because for a
    small α the difference of two exponentials both close to 1 loses most of its significant digits.

    A modifier that is constant along the branch — ``FromParent``, the variable-rates σ² drift —
    factors straight out of either integral, so it composes with OU for free.

    ``inherited`` is the lineage's `FromParent` factor (variable-rates
    BM), constant along the branch, threaded in by the caller and passed through to the rate; it
    factors straight out of the integral. A rate with no ``FromParent`` modifier ignores it.

    ``ltt`` is the tree's lineages-through-time function when the rate carries a ``OnTotalDiversity`` modifier
    (diversity-dependent σ²): the integral then also steps at the tree's speciation / extinction times,
    reading the standing diversity on each sub-interval. ``None`` when σ² does not depend on diversity.

    ``trajs`` (with ``node_id``, the lineage this branch is) are the driver trajectories when the rate
    carries a `DrivenBy` — σ² read off another level. The driver's
    value on this lineage is threaded in as ``drivers``, and the integral **steps where the driver
    switches** (``next_change``) exactly as it steps at a skyline breakpoint: a discrete driver
    switches *mid-branch*, so a single sample per branch would credit the whole branch to whichever
    state happened to be read. ``None`` when σ² is not driven.
    (Stepping is O(events the branch crosses); fine for the trait level's one value per branch.)"""
    if pull > 0.0 and not rate.modifiers:
        # plain OU: σ² is constant along the branch, so the weighted integral collapses to the closed
        # form the OU transition law is written with. Kept as its own expression, in that exact
        # spelling, so a run written before the weighted integral existed draws exactly the numbers
        # it drew then — the expm1 form below is the same quantity but not the same last bits.
        e = math.exp(-pull * (t1 - t0))
        return rate.effective(lineages=1) / (2.0 * pull) * (1.0 - e * e)
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
        sigma2 = rate.effective(lineages=1, time=t, carried=inherited, diversity=div,
                                drivers=drivers)
        if pull > 0.0:
            # the OU weight for this piece: what survives of the noise injected in [t, nxt) once the
            # pull has acted on it for the rest of the branch. Both factors matter — the first is the
            # decay after the piece ends, the second the decay accumulated inside it.
            total += sigma2 * math.exp(-2.0 * pull * (t1 - nxt)) \
                * (-math.expm1(-2.0 * pull * (nxt - t))) / (2.0 * pull)
        else:
            total += sigma2 * (nxt - t)
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



def _per_trait(value, traits: list, name: str, *, positive: bool = False,
               nonnegative: bool = False) -> np.ndarray:
    """One number per trait, from either spelling of a per-trait argument: a **bare number** means
    the same value for every trait, a **dict** ``{trait: value}`` gives them one at a time.

    This is the shape ``start`` and ``rate`` already use for a correlated trait set, and SPEC §7's
    one-concept-one-word rule says the OU optimum, the OU strength and the on-speciation jump
    variance take it too rather than inventing parallel spellings. The bare number is the same
    shortcut a bare number is in the rate grammar (SPEC §5): the common case says it once."""
    if isinstance(value, dict):
        if set(value) != set(traits):
            raise ValueError(
                f"{name} must name the same traits as start; got {sorted(value)} vs {sorted(traits)}"
            )
        given = [(f"{name}[{t!r}]", value[t]) for t in traits]
    elif isinstance(value, (list, tuple, np.ndarray)):
        # a bare sequence carries no trait names, so which entry belongs to which trait would depend
        # on the order a dict happened to be written in. The dict says it.
        raise ValueError(
            f"{name} must be one number (shared across the traits) or a dict {{trait: value}} over "
            f"{sorted(traits)}; a bare sequence names no traits, got {value!r}"
        )
    else:
        given = [(name, value) for _ in traits]
    out = np.empty(len(traits))
    for i, (where, v) in enumerate(given):
        if isinstance(v, bool) or not isinstance(v, (int, float)) or not math.isfinite(v):
            raise ValueError(f"{where} must be a finite number, got {v!r}")
        if positive and v <= 0:
            raise ValueError(f"{where} must be a finite positive number, got {v!r}")
        if nonnegative and v < 0:
            raise ValueError(f"{where} must be a non-negative number, got {v!r}")
        out[i] = float(v)
    return out



def _simulate_regimes(tree, start, rate, reverts_to, pull, regimes, at_speciation, seed,
                      progress=False) -> TraitsResult:
    """Multi-optimum OU — the optimum shifts by **regime**, a discrete stochastic map painted on the
    *same* tree (typically a `simulate_discrete()` run). Along each branch the value follows OU
    toward the current regime's optimum, integrated **exactly** across the regime's ``(state,
    duration)`` segments (a regime may switch part-way along a branch); convention B paints the root
    branch too. ``reverts_to`` is a dict ``{regime_state: θ}``; ``pull`` (α > 0) and ``rate`` (σ²) are
    shared across regimes (the OUM variant).

    ``at_speciation`` layers the on-speciation jump on top, exactly as it does for a single-optimum
    run: one jump variance, drawn at each split before the branch's own OU integration. It is shared
    across regimes — a per-regime jump size would be a second table to keep consistent with
    ``reverts_to``, and nothing yet asks for one."""
    if getattr(regimes, "kind", None) != "discrete" or regimes.history is None:
        raise ValueError(
            "regimes must be a discrete TraitsResult carrying a stochastic map — paint them with "
            "simulate_discrete(...) on this same tree."
        )
    if set(regimes.history) != set(tree.nodes):
        raise ValueError("regimes must be painted on the SAME tree this trait rides (node ids differ).")
    if isinstance(at_speciation, dict):
        raise ValueError(
            "at_speciation is one jump variance, shared across regimes; a jump size per regime "
            "is not implemented yet."
        )
    if isinstance(start, dict) or isinstance(rate, dict):
        raise ValueError(
            "multi-optimum OU over several correlated traits is not implemented yet — it would need "
            "an optimum per regime per trait, reverts_to={regime: {trait: θ}}. Give a scalar start "
            "and rate for one trait."
        )
    if not isinstance(reverts_to, dict):
        raise ValueError("with regimes, reverts_to is a dict {regime_state: optimum θ}, one per regime.")
    if isinstance(pull, bool) or not isinstance(pull, (int, float)) or not math.isfinite(pull) or pull is None \
            or pull <= 0:
        raise ValueError(f"regimes needs pull (the OU strength α > 0), got {pull!r}")
    if isinstance(start, bool) or not isinstance(start, (int, float)) or not math.isfinite(start):
        raise ValueError(f"start must be a finite number, got {start!r}")
    r = as_rate(rate, default_scope=PerLineage)
    if not isinstance(r.scope, PerLineage):
        raise ValueError("rate must be per lineage with regimes — drop the scope wrapper.")
    if r.modifiers:
        # this path walks the regime map's (state, duration) segments, not absolute time, so it has
        # none of the plumbing a modified σ² needs (the schedule's breakpoints, the standing
        # diversity, a driver's switch times). The single-optimum OU path does take them.
        carried = ", ".join(dict.fromkeys(describe(m) for m in r.modifiers))
        raise ValueError(
            f"rate carries {carried}, and a modified variance-rate combined with regimes is not "
            f"implemented yet — give a bare variance-rate here, or drop regimes and use "
            f"reverts_to / pull, which do take a modified σ²."
        )
    sigma2, alpha = r.effective(lineages=1), float(pull)
    used = {seg[0] for segs in regimes.history.values() for seg in segs}
    missing = used - set(reverts_to)
    if missing:
        raise ValueError(f"reverts_to is missing an optimum for regime state(s) {sorted(missing)}")

    jump_sd = _at_speciation_jump_sd(at_speciation)  # on-speciation jump width (0 if not requested)

    rng, seed = stream("traits", seed)      # own stream, and a drawn seed if none was given
    node_values: dict[int, float] = {}
    root = tree.nodes[tree.root]
    # the same log a single-optimum run carries (SPEC §2): the value at t=0, then each jump at a
    # split. A diffusion has no along-branch events to record, whichever optimum it is heading for.
    events: list[Change] = [Change(root.birth_time, "initial", tree.root, None, float(start))]
    for i in _preorder(tree, progress):
        node = tree.nodes[i]
        x = float(start) if node.parent is None else node_values[node.parent]
        if node.parent is not None and jump_sd > 0.0:
            jumped = x + float(rng.normal(0.0, jump_sd))  # on speciation: a jump at the split…
            events.append(Change(node.birth_time, "on_speciation", i, x, jumped))
            x = jumped                                    # …then the regime-by-regime anagenesis
        for regime, dt in regimes.history[i]:  # integrate OU piece-by-piece across the regime segments
            theta = float(reverts_to[regime])
            e = math.exp(-alpha * dt)
            var = sigma2 / (2.0 * alpha) * (1.0 - e * e)
            x = theta + (x - theta) * e + (float(rng.normal(0.0, math.sqrt(var))) if var > 0.0 else 0.0)
        node_values[i] = x
    return TraitsResult(tree, cast("dict[int, object]", node_values), events, seed)



def _simulate_correlated(tree, start, rate, reverts_to, pull, correlation, at_speciation, seed,
                         progress=False) -> TraitsResult:
    """Correlated continuous traits in **one call** (the joint rule inside a level). ``start`` and
    ``rate`` are dicts over the same trait names; the branch increment is drawn from ``MVN(0, Σ·dt)``
    with ``Σ = D R D`` (``D = diag(σ_i)``, ``R`` from ``correlation``), so at a tip the correlation
    between two traits equals their ρ. Convention B holds: the root diffuses over its own branch.

    ``reverts_to`` and ``pull`` turn it into **multivariate Ornstein–Uhlenbeck in its diagonal-drift
    restriction**: each trait reverts to its own optimum at its own strength, and the correlation
    stays where it was, in the diffusion. The branch transition is then exactly

        ``mean_i = θ_i + (x_i − θ_i)·e^{−α_i·dt}``
        ``Cov_ij = Σ_ij·(1 − e^{−(α_i+α_j)·dt}) / (α_i + α_j)``

    which is closed-form, needs no matrix exponential, and reduces to the univariate OU law on the
    diagonal. What it **excludes** is a full drift matrix ``A`` — one trait's deviation pulling
    another, the off-diagonal terms of mvMORPH's general OU — because that is a different model, not
    a parameterisation of this one: it needs ``e^{−A·dt}`` and a Lyapunov solve for the covariance,
    and it raises questions (a non-symmetric or defective ``A``) that a diagonal drift never poses.
    A ``pull`` given as a matrix is refused by name rather than quietly read as its diagonal.

    ``at_speciation`` adds the jump at each split, and the jump is drawn under the **same**
    ``correlation`` overlay the diffusion uses (``MVN(0, D_v R D_v)`` with ``v`` the per-trait jump
    variance). Traits declared correlated whose cladogenetic component was independent would be a
    silent mismatch between what the argument says and what the run does."""
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
    is_ou = reverts_to is not None or pull is not None
    if is_ou:
        if isinstance(pull, (list, tuple, np.ndarray)) and len(pull) > 0 \
                and isinstance(pull[0], (list, tuple, np.ndarray)):
            # the honest rejection: not "you spelled it wrong" but "this is a different model".
            raise ValueError(
                "pull is a full drift matrix (one trait's deviation pulling another), which is not "
                "implemented yet — give pull as one number (the strength shared across traits) or a "
                "dict {trait: strength}, which is the diagonal-drift multivariate OU this engine has."
            )
        if reverts_to is None or pull is None:
            raise ValueError(
                "Ornstein–Uhlenbeck needs both reverts_to (the optimum) and pull (the strength); "
                "give both, or neither for Brownian motion."
            )
        theta_vec = _per_trait(reverts_to, traits, "reverts_to")
        alpha_vec = _per_trait(pull, traits, "pull", positive=True)
    jump_var = None if at_speciation is None else \
        _per_trait(at_speciation, traits, "at_speciation", nonnegative=True)

    sigma2 = np.empty(len(traits))
    for i, name in enumerate(traits):
        if isinstance(start[name], bool) or not isinstance(start[name], (int, float)) \
                or not math.isfinite(start[name]):
            raise ValueError(f"start[{name!r}] must be a finite number, got {start[name]!r}")
        r = as_rate(rate[name], default_scope=PerLineage)
        if not isinstance(r.scope, PerLineage):
            raise ValueError(f"rate[{name!r}] must be per lineage — drop the scope wrapper.")
        if r.modifiers:
            carried = ", ".join(dict.fromkeys(describe(m) for m in r.modifiers))
            raise ValueError(
                f"rate[{name!r}] carries {carried}; per-trait modifiers combined with correlation "
                f"are not implemented yet — with a σ² that moves, the branch covariance is "
                f"∫ ρ_ij·σ_i(t)·σ_j(t) dt, which is not each trait's own integral recombined. Use "
                f"bare per-trait rates here."
            )
        sigma2[i] = r.effective(lineages=1)

    R = _correlation_matrix(traits, correlation)
    sd = np.sqrt(sigma2)
    Sigma = (sd[:, None] * R) * sd[None, :]                    # Σ = D R D (Σ_ij = σ_i σ_j ρ_ij)
    sigma = _symmetric_sqrt(Sigma)
    # the jump rides the same overlay: its covariance is D_v R D_v, so two traits declared correlated
    # jump together too. All-zero variances mean no jump at all — and, deliberately, no draw, so
    # at_speciation=0.0 leaves a run byte-identical to one that never asked for jumps.
    jump_sqrt = None
    if jump_var is not None and float(jump_var.max()) > 0.0:
        jsd = np.sqrt(jump_var)
        jump_sqrt = _symmetric_sqrt((jsd[:, None] * R) * jsd[None, :])
    # a shared α scales Σ by one number, so its square root is the scaled Σ-root and the one eigh
    # above serves every branch; a per-trait α mixes the traits differently in each entry of the
    # covariance, so that case builds its own root per branch length.
    shared_alpha = float(alpha_vec[0]) if is_ou and bool(np.all(alpha_vec == alpha_vec[0])) else None
    if is_ou:
        alpha_sum = np.add.outer(alpha_vec, alpha_vec)         # (α_i + α_j), the OU covariance divisor
    start_vec = np.array([float(start[t]) for t in traits])
    k = len(traits)

    rng, seed = stream("traits", seed)      # own stream, and a drawn seed if none was given
    node_values: dict[int, dict[str, float]] = {}
    # The log every other continuous run carries: the value the run started from, and each jump at a
    # split. A diffusion cannot be rebuilt from events — that is as true here as for one trait — so
    # this is a record of the run's discrete moments, not a driver file. The rows hold the whole
    # **vector**, because a correlated jump moves every trait at once and splitting it into one row per
    # trait would suggest they were separate events.
    def _vec(v) -> dict:
        return {t: float(v[j]) for j, t in enumerate(traits)}

    events: list[Change] = [Change(tree.nodes[tree.root].birth_time, "initial", tree.root, None,
                                   _vec(start_vec))]
    for i in _preorder(tree, progress):
        node = tree.nodes[i]
        x = start_vec if node.parent is None else np.array([node_values[node.parent][t] for t in traits])
        if node.parent is not None and jump_sqrt is not None:
            before = x
            x = x + jump_sqrt @ rng.standard_normal(k)   # on speciation: a jump at the split…
            events.append(Change(node.birth_time, "on_speciation", i, _vec(before), _vec(x)))
        dt = node.end_time - node.birth_time              # …then anagenesis along the branch
        if not is_ou:
            vec = x + (math.sqrt(dt) * (sigma @ rng.standard_normal(k)) if dt > 0.0 else 0.0)
        elif dt <= 0.0:
            vec = x
        else:
            mean = theta_vec + (x - theta_vec) * np.exp(-alpha_vec * dt)
            if shared_alpha is not None:
                root = math.sqrt(-math.expm1(-2.0 * shared_alpha * dt) / (2.0 * shared_alpha)) * sigma
            else:
                root = _symmetric_sqrt(Sigma * (-np.expm1(-alpha_sum * dt)) / alpha_sum)
            vec = mean + root @ rng.standard_normal(k)
        node_values[i] = {t: float(vec[j]) for j, t in enumerate(traits)}

    return TraitsResult(tree, cast("dict[int, object]", node_values), events, seed)



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
    the correlation between two traits is exactly their ρ. Add ``reverts_to`` and ``pull`` — one
    value shared, or a dict of one per trait — and it is **multivariate Ornstein–Uhlenbeck in its
    diagonal-drift restriction**: each trait reverts to its own optimum at its own strength, the
    correlation stays in the diffusion, and the branch covariance is
    ``Σ_ij·(1 − e^{−(α_i+α_j)·dt})/(α_i + α_j)``. One trait's deviation pulling *another* — a full
    drift matrix — is a different model and is refused by name, not read as a diagonal. A correlated
    run takes bare per-trait rates. Its log is **widened** rather than absent: a value is a per-trait
    vector, so ``trait_events.tsv`` gets one ``from``/``to`` column pair per trait.

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
    ``DrivenBy(driver, {…})`` modifier makes σ² **read another level** — the driver grown first on this
    same tree and handed over as its result object or its written ``trait_events.tsv``, so a lineage
    diffuses faster while the driver is in one state than another. A discrete driver switches
    *mid-branch*, and the per-branch variance is the integral across those pieces, so a branch that
    spends half its length in the fast state accrues exactly half the fast variance.

    ``reverts_to`` (the optimum θ) and ``pull`` (the strength α > 0) turn the diffusion into
    Ornstein–Uhlenbeck — the value is pulled toward θ while it diffuses, the exact per-branch
    transition being ``Normal(θ + (x−θ)·e^{−α·dt}, σ²/(2α)·(1−e^{−2α·dt}))``. Give **both** or
    neither. The optimum and the pull **compose with the σ² modifiers**: a trait that bursts early
    and reverts to an optimum is one rate with one modifier and two arguments. A σ² that moves along
    the branch leaves the mean untouched (it never read σ²) and makes the variance the exact
    pull-weighted integral ``∫ e^{−2α(t₁−s)}·σ²(s) ds``, stepping where the schedule, the standing
    diversity or the driver steps. That weight is the whole difference from Brownian motion's
    ``∫ σ²(s) ds``: under OU, variance accrued early has been pulled back toward θ by the time the
    branch ends, so the two integrals differ by an order of magnitude on a typical branch.

    ``at_speciation`` adds an **on-speciation** jump — ``Normal(0, at_speciation)`` on each daughter at
    every speciation (the punctuational mode), layered on top of the along-branch anagenesis. Under
    ``correlation=`` it takes one variance per trait (or one shared) and the jump is drawn under the
    same overlay the diffusion uses.
    ``regimes`` gives **multi-optimum OU**: pass a discrete `TraitsResult` (a stochastic map
    painted by `simulate_discrete()` on this same tree) and a per-regime ``reverts_to={regime: θ}``,
    and the value follows OU toward whichever regime's optimum a branch is in; it takes
    ``at_speciation`` too, one jump variance shared across regimes, and it takes a bare σ² — a
    modified variance-rate with ``regimes`` is not implemented yet. Deterministic given ``seed``.
    """
    tree = as_tree(tree, level="traits")
    if regimes is not None:
        if correlation is not None:
            # `regimes` dispatches before the correlated engine and threads no correlation, so a
            # correlation passed here would be read by nothing and the run would silently be the
            # uncorrelated model (SPEC §5: refuse, never ignore). The obstacle is the same one that
            # keeps a modified σ² out of `regimes`: the branch covariance under multi-optimum OU is
            # ∫ρ_ij·σ_i·σ_j weighted by each regime's own pull, which this engine does not integrate.
            raise ValueError(
                "correlation= with regimes= is not implemented yet: multi-optimum OU evolves one "
                "trait, so there is no second trait for a correlation to be with. Use correlation= "
                "on its own (correlated BM/OU), or regimes= on its own (multi-optimum OU).")
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
    # and DrivenBy (σ² driven by another level) are the σ² modifiers this engine supports; anything
    # else is rejected loudly — the genome engine's discipline.
    for m in r.modifiers:
        if m.reads == (DRAWN, "family"):
            # not a missing feature: there is nothing here for it to mean
            raise ValueError(
                "rate carries ByFamily, but a trait has no gene families — ByFamily belongs on a "
                "genomes rate. For per-lineage heterogeneity here use FromParent (variable-rates BM)."
            )
        if not is_implemented(m, IMPLEMENTED_MODIFIERS, "traits.continuous"):
            raise ValueError(
                f"rate carries {describe(m)}, which the continuous trait engine does not "
                f"support. It takes OnTime (early burst), FromParent (variable-rates BM), "
                f"OnTotalDiversity (diversity-dependent), and DrivenBy (driven by another level)."
            )
        if isinstance(m, DrivenBy):
            check_not_a_kernel(m.mapping, label="rate")
    # the per-lineage modifiers σ² carries (variable-rates BM), asked for the same way every level
    # asks — and every one of them is kept, so two compose rather than the second going quietly.
    drift = tuple(m for m, _ in r.carried(unit="lineage"))
    check_one_memory(drift, label="rate", unit="lineage")
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
        theta, alpha = float(reverts_to), float(pull)

    # conditioning: a σ² carrying DrivenBy reads another level, grown first on this same tree. Resolve
    # each driver once into a trajectory (value + next-switch, keyed by the shared node id), from a
    # written trait log or a grown result handed over in memory. Undriven ⇒ empty, and the walk below
    # is exactly the walk it was — no driver, no lookup, no change to the draw order.
    trajs = _resolve_drivers(_driven_mods(r), tree, "traits.continuous")

    jump_sd = _at_speciation_jump_sd(at_speciation)  # on-speciation jump width (0 if not requested)

    rng, seed = stream("traits", seed)      # own stream, and a drawn seed if none was given
    ltt = _LTT(tree) if has_diversity else None  # the standing-diversity curve, when σ² reads it
    node_values: dict[int, float] = {}
    root = tree.nodes[tree.root]
    # the initial value at t=0 — the origin the log reconstructs from (SPEC §2). A diffusion cannot be
    # rebuilt from events, but the row keeps the file's shape uniform across trait kinds.
    events: list[Change] = [Change(root.birth_time, "initial", tree.root, None, float(start))]
    inh: dict[int, tuple[float, ...]] = {}  # each lineage's σ² drift factors (variable-rates BM),
                                           # one per carried modifier, constant along its branch
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
            inh[i] = carried_at_birth(drift, rng)
        else:
            inh[i] = carried_at_split(drift, inh[node.parent], rng)
        t0, t1 = node.birth_time, node.end_time
        if is_ou:
            e = math.exp(-alpha * (t1 - t0))       # mean-reversion toward θ over the branch
            mean = theta + (x - theta) * e         # the mean does not read σ², so a modified σ² leaves it
            # …and the variance is the pull-weighted integral, which for a bare σ² is exactly the
            # closed form σ²/(2α)·(1−e^{−2α·dt}) this branch used before the weight existed.
            var = _accrued_variance(r, t0, t1, inherited=product(inh[i]), ltt=ltt, trajs=trajs, node_id=i,
                                    pull=alpha)
        else:
            mean = x                                # pure diffusion (BM / early burst / variable-rates)
            var = _accrued_variance(r, t0, t1, inherited=product(inh[i]), ltt=ltt, trajs=trajs, node_id=i)
        std = math.sqrt(var) if var > 0.0 else 0.0
        node_values[i] = mean + (float(rng.normal(0.0, std)) if std > 0.0 else 0.0)

    return TraitsResult(tree, cast("dict[int, object]", node_values), events, seed)


# --- discrete traits: a state switching along the tree (Mk) ------------------------------------

