"""Traits — the discrete engine: the Mk state-switching model (Gillespie, its switch rate optionally driven by another trait) and the threshold model (simulate_discrete), plus the DiscreteTrait process spec for joint runs."""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from ..rates.mapping import check_not_a_kernel
from ..rng import stream
from ..rates.modifiers import DrivenBy, Modifier, is_wired
from ..rates.rate import Rate, as_rate
from ..rates.scope import PerLineage
from ..tree import as_tree

from ._shared import _correlation_matrix, _preorder, _resolve_drivers, _symmetric_sqrt
from .result import Change, TraitsResult

#: the modifiers a discrete trait's ``switch`` rate takes — declared, as every level declares its
#: own, so Appendix A and the CLI help are built from the engine rather than kept by hand. It is a
#: shorter list than the continuous rate's: a switch rate reads a driver and nothing else.
WIRED_MODIFIERS = (DrivenBy,)


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


# --- a driven switch rate: the generator reads another trait -------------------------------------

def _switch_specs(switch) -> list:
    """The per-transition rate specs a ``switch`` carries — the values of a ``{'from->to': rate}``
    dict, or the single symmetric spec. (A ``k×k`` matrix is numbers by construction, so it has
    none: a matrix entry cannot carry a modifier.)"""
    if isinstance(switch, dict):
        return list(switch.values())
    return [switch]


def _switch_drivers(switch) -> list:
    """The `DrivenBy` modifiers the switch rates carry — a switch rate written as a rate expression
    (``0.4 * mod.DrivenBy(habitat, {"aquatic": 3.0})``) rather than as a bare number, so the trait
    switches faster on the lineages where the driver is in one state than another.

    ``DrivenBy`` is the only modifier a switch rate takes; anything else on it would be read by no
    part of this engine, so it is refused by name rather than silently ignored."""
    mods = []
    for spec in _switch_specs(switch):
        if not isinstance(spec, (Rate, Modifier)):
            continue                       # a bare number (or a matrix): nothing to drive
        r = as_rate(spec, default_scope=PerLineage)
        if not isinstance(r.scope, PerLineage):
            raise ValueError(
                f"a switch rate has a {type(r.scope).__name__} scope, but a discrete trait switches "
                f"per lineage — drop the scope wrapper (per lineage is the default).")
        for m in r.modifiers:
            if not is_wired(m, WIRED_MODIFIERS, "traits.discrete"):
                raise ValueError(
                    f"a switch rate carries {type(m).__name__}, which the discrete trait engine does "
                    f"not support. It takes DrivenBy (the switch rate driven by another trait).")
            if isinstance(m, DrivenBy):   # a modifier wired in from outside carries no mapping
                check_not_a_kernel(m.mapping, label="a switch rate")
            mods.append(m)
    return mods


def _driven_entries(states, switch) -> list:
    """The off-diagonal switch rates as ``[(i, j, Rate)]`` — the driven twin of `_q_matrix()`, whose
    entries are rate *specs* rather than settled numbers, so the generator can be rebuilt wherever
    the driver switches. The same two shapes a driven rate can be written in: one symmetric rate
    (every ``i → j`` alike) or a ``{'from->to': rate}`` dict (only the named transitions, the rest
    zero). Validation mirrors `_q_matrix()`'s, the rate itself being checked by its scope wrapper."""
    k = len(states)
    idx = {s: i for i, s in enumerate(states)}
    if isinstance(switch, dict):
        entries = []
        for key, spec in switch.items():
            parts = [p.strip() for p in str(key).split("->")]
            if len(parts) != 2:
                raise ValueError(f"switch keys must read 'from->to', got {key!r}")
            frm, to = parts
            if frm not in idx or to not in idx:
                raise ValueError(f"switch key {key!r} names a state not in states={list(states)}")
            if frm == to:
                raise ValueError(f"switch key {key!r} is a self-transition; only i→j (i≠j) is a rate")
            entries.append((idx[frm], idx[to], as_rate(spec, default_scope=PerLineage)))
        return entries
    r = as_rate(switch, default_scope=PerLineage)
    return [(i, j, r) for i in range(k) for j in range(k) if i != j]


def _driven_q(entries, k: int, drivers: dict, time: float) -> np.ndarray:
    """The generator at one instant — each off-diagonal entry evaluated against the driver values on
    this lineage right now, the diagonal then set to minus its row sum (rows sum to zero)."""
    Q = np.zeros((k, k))
    for i, j, r in entries:
        Q[i, j] = r.effective(lineages=1, time=time, drivers=drivers)
    np.fill_diagonal(Q, -Q.sum(axis=1))
    return Q


def _merge_runs(segments: list) -> list:
    """Adjacent equal-state stretches summed into one. Cutting a branch at the driver's switches can
    leave two pieces of the same state either side of a cut the trait never crossed; merging them
    keeps ``.history`` and the event log reading exactly as the undriven walk writes them (only real
    transitions recorded)."""
    merged: list[tuple[int, float]] = []
    for state, dur in segments:
        if merged and merged[-1][0] == state:
            merged[-1] = (state, merged[-1][1] + dur)
        else:
            merged.append((state, dur))
    return merged


def _gillespie_driven(state: int, node, node_id: int, entries, k: int, trajs: dict,
                      rng) -> tuple[int, list]:
    """Exact CTMC along one branch whose switch rates are **driven** by another trait.

    The driver is piecewise-constant along the branch, so the generator is too: the branch is cut at
    the driver's own switches (`~zombi2.rates.driver.DriverTrajectory.next_change`) and the plain
    `_gillespie()` runs over each stretch under that stretch's ``Q``, carrying the state across.
    Memorylessness makes the concatenation exactly the CTMC with the time-varying generator. Reading
    the driver once per branch would be wrong — a driver that flips half-way would have the whole
    branch run at whichever state happened to be read.

    Returns the same ``(end_state, segments)`` `_gillespie()` returns."""
    t0, t1 = node.birth_time, node.end_time
    if t1 <= t0:
        return state, [(state, 0.0)]
    segments: list[tuple[int, float]] = []
    t, cur = t0, state
    while t < t1:
        nxt = t1
        for traj in trajs.values():
            nxt = min(nxt, traj.next_change(node_id, t))
        drivers = {key: traj.value(node_id, t) for key, traj in trajs.items()}
        cur, segs = _gillespie(cur, nxt - t, _driven_q(entries, k, drivers, t), rng)
        segments.extend(segs)
        t = nxt
    return cur, _merge_runs(segments)



def _finite(x, name: str) -> float:
    """Coerce ``x`` to a finite float or raise a clear ``ValueError`` naming ``name``."""
    if isinstance(x, bool) or not isinstance(x, (int, float)) or not math.isfinite(x):
        raise ValueError(f"{name} must be a finite number, got {x!r}")
    return float(x)



def _liability_sigma2(spec, name) -> float:
    """The liability's variance-rate σ² — a bare per-lineage rate (no modifiers)."""
    r = as_rate(spec, default_scope=PerLineage)
    where = "liability" if name is None else f"liability[{name!r}]"
    if not isinstance(r.scope, PerLineage) or r.modifiers:
        raise ValueError(
            f"{where} must be a bare variance-rate, per lineage; a modified liability variance-rate "
            f"is not implemented yet."
        )
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



def _simulate_threshold(tree, states, liability, threshold, start, correlation, seed,
                        progress=False) -> TraitsResult:
    """The Wright–Felsenstein **threshold** model — a discrete state read off an underlying continuous
    Brownian liability: the liability diffuses (convention B), and the observed state is which
    ``threshold``-cut interval it lands in. With ``liability`` a dict + ``correlation``, several traits'
    liabilities diffuse **jointly** (correlated discrete traits, one call), each cut by the shared
    thresholds. There is no Gillespie stochastic map — the discreteness is a cut through a continuous
    path — so ``.history`` is ``None`` and ``.events`` is empty."""
    thr = _threshold_cuts(states, threshold)

    def label(x):
        return states[int(np.searchsorted(thr, x))]

    rng, seed = stream("traits", seed)      # own stream, and a drawn seed if none was given
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
        for i in _preorder(tree, progress):
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
        for i in _preorder(tree, progress):
            node = tree.nodes[i]
            x_s = start_liab if node.parent is None else liab_s[node.parent]
            dt = node.end_time - node.birth_time
            std = math.sqrt(sig2 * dt) if sig2 > 0.0 and dt > 0.0 else 0.0
            liab_s[i] = x_s + (float(rng.normal(0.0, std)) if std > 0.0 else 0.0)
            node_values[i] = label(liab_s[i])

    return TraitsResult(tree, node_values, [], seed, kind="threshold")



def simulate_discrete(tree, *, states, switch=None, start=None, liability=None, threshold=None,
                      correlation=None, at_speciation=None, seed=None,
                      progress=False) -> TraitsResult:
    """Evolve a discrete-state trait down a tree and return a `TraitsResult`. Two mechanisms:

    - **Mk** (``switch=``) — a continuous-time Markov chain over the ``states``, simulated **exactly**
      by Gillespie along every branch, so each node's ``(state, duration)`` segments *are* the realized
      history (``.history``) and ``.events`` reads off the transitions. ``switch`` is a symmetric rate
      (``0.1``), a ``{"marine->terrestrial": 0.1}`` dict, or a ``k×k`` matrix (see `_q_matrix()`).
      ``start`` is the root state (a label in ``states``; ``None`` draws one uniformly). A switch rate
      may be **driven by another trait** grown first on this same tree — write it as a rate
      expression, ``switch=0.4 * mod.DrivenBy(habitat, {"aquatic": 3.0})`` or per transition,
      ``switch={"a->b": 0.2 * mod.DrivenBy(habitat, {"aquatic": 3.0}), "b->a": 0.2}``. The driver
      switches mid-branch, so the generator is rebuilt at each of its switches and the branch is
      simulated piece by piece — the exact CTMC with a time-varying generator, not one sample per
      branch.
    - **Threshold** (``liability=`` + ``threshold=``) — the Wright–Felsenstein model: a discrete state
      read off an underlying continuous Brownian **liability** (variance-rate ``liability``), cut into
      ``states`` by the ``threshold`` cut point(s) (``k−1`` increasing cuts for ``k`` states). ``start``
      is the initial *liability* (a number, default 0.0). Give ``liability`` as a dict + a
      ``correlation={(a, b): ρ}`` overlay to evolve **correlated** discrete traits jointly — their
      liabilities diffuse together (``Σ = D R D``) and each is cut by the shared thresholds. A threshold
      trait has no Gillespie map, so ``.history`` is ``None`` and ``.events`` empty.

    ``tree`` is the **complete** species tree (a `Tree` or
    `SpeciesResult`); the trait evolves on every lineage (convention B: the root
    diffuses over its own branch), and ``.values`` reads the extant tips. On an Mk trait,
    ``at_speciation`` (a probability in ``[0, 1]``) adds an **on-speciation** shift — each daughter hops
    to a uniformly-chosen other state with that chance at every speciation. Deterministic given ``seed``.
    """
    tree = as_tree(tree, level="traits")
    states = list(states)
    if len(states) < 2:
        raise ValueError(f"a discrete trait needs at least 2 states, got {states!r}")
    if len(set(states)) != len(states):
        raise ValueError(f"states must be unique, got {states!r}")
    if liability is not None or threshold is not None:
        if switch is not None:
            raise ValueError("give switch= (an Mk trait) OR liability=/threshold= (a threshold trait), not both")
        if at_speciation is not None:
            raise ValueError("at_speciation is not implemented for threshold traits yet — it applies to Mk (switch=) traits")
        return _simulate_threshold(tree, states, liability, threshold, start, correlation, seed,
                                   progress)
    if correlation is not None:
        raise ValueError("correlation= on a discrete trait needs the threshold model — give liability= and threshold=")
    if switch is None:
        raise ValueError("give switch= — the transition rate(s) between the discrete states.")
    # conditioning: a switch rate carrying DrivenBy reads another trait, grown first on this same
    # tree. The generator is then a function of the driver, so it is built per stretch rather than
    # once. No DrivenBy ⇒ `trajs` is empty and the walk below is exactly the walk it was.
    sw_mods = _switch_drivers(switch)
    entries = _driven_entries(states, switch) if sw_mods else None
    Q = None if sw_mods else _q_matrix(states, switch)
    trajs = _resolve_drivers(sw_mods, tree)
    if at_speciation is not None and (isinstance(at_speciation, bool)
            or not isinstance(at_speciation, (int, float)) or not 0.0 <= at_speciation <= 1.0):
        raise ValueError(f"at_speciation must be a probability in [0, 1] (the shift chance), got {at_speciation!r}")
    shift = 0.0 if at_speciation is None else float(at_speciation)

    rng, seed = stream("traits", seed)      # own stream, and a drawn seed if none was given
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
    root = tree.nodes[tree.root]
    # the initial state at t=0 — the origin the log reconstructs from: tree + this + the switches give
    # the driver on every lineage, so the event log is the conditioning file (no separate driver).
    events: list[Change] = [Change(root.birth_time, "initial", tree.root, None, states[start_i])]
    for i in _preorder(tree, progress):
        node = tree.nodes[i]
        # the root starts from `start` at t=0 and evolves over its own branch; every other node from
        # its parent's end state (parent < i, already set) — the same convention-B walk as continuous.
        cur = start_i if node.parent is None else idx[node_values[node.parent]]
        if node.parent is not None and shift > 0.0 and float(rng.random()) < shift:
            j = int(rng.integers(len(states) - 1))  # on speciation: hop to a uniform *other* state
            new = j if j < cur else j + 1
            events.append(Change(node.birth_time, "on_speciation", i, states[cur], states[new]))
            cur = new
        if Q is not None:
            end_i, segs = _gillespie(cur, node.end_time - node.birth_time, Q, rng)
        else:  # a driven switch rate: no constant generator — cut the branch where the driver switches
            end_i, segs = _gillespie_driven(cur, node, i, entries, len(states), trajs, rng)
        t = node.birth_time  # the transitions between the Gillespie segments are the on-branch events
        for (s1, d1), (s2, _d) in zip(segs, segs[1:]):
            t += d1
            events.append(Change(t, "on_branch", i, states[s1], states[s2]))
        node_values[i] = states[end_i]
    events.sort(key=lambda c: c.time)
    return TraitsResult(tree, node_values, events, seed, kind="discrete")


# --- process specs: a trait bundled but UNEXECUTED, for a joint model to grow with the tree --------


@dataclass(frozen=True)
class DiscreteTrait:
    """A discrete (Mk) trait **process** — its parameters bundled but not yet run (SPEC §4).
    ``simulate_discrete(tree, ...)`` is the runner that grows this on a *fixed*
    tree; a **joint** model instead takes this spec and grows the trait *with* the tree it drives
    (``joint.simulate_joint(trait=traits.discrete(...))``), so neither can be simulated first. Same
    parameters as `simulate_discrete()` (the Mk half): ``states``, ``switch`` (the rate spec),
    ``start`` (the root state, ``None`` = uniform), ``at_speciation`` (the on-speciation shift
    probability)."""

    states: tuple
    switch: object
    start: object = None
    at_speciation: object = None

    def _resolve(self, rng):
        """Build the concrete CTMC the engine grows: ``(states_list, Q, start_index, shift_prob)`` —
        the same setup `simulate_discrete()` does, so a joint run and a fixed-tree run share one
        trait model. ``rng`` draws the root state when ``start`` is ``None``."""
        states = list(self.states)
        Q = _q_matrix(states, self.switch)
        idx = {s: i for i, s in enumerate(states)}
        if self.start is None:
            start_i = int(rng.integers(len(states)))
        elif self.start in idx:
            start_i = idx[self.start]
        else:
            raise ValueError(f"start must be one of states={states} (or None for a uniform draw), got {self.start!r}")
        if self.at_speciation is not None and (isinstance(self.at_speciation, bool)
                or not isinstance(self.at_speciation, (int, float)) or not 0.0 <= self.at_speciation <= 1.0):
            raise ValueError(f"at_speciation must be a probability in [0, 1] (the shift chance), got {self.at_speciation!r}")
        shift = 0.0 if self.at_speciation is None else float(self.at_speciation)
        return states, Q, start_i, shift



def discrete(*, states, switch=None, start=None, at_speciation=None) -> DiscreteTrait:
    """A discrete-trait **process spec** — `DiscreteTrait`, unexecuted — for a joint model to
    grow with the tree it drives (``joint.simulate_joint(trait=traits.discrete(states=[...], switch=...))``).
    A thin bundle of `simulate_discrete()`'s Mk parameters; validated when the joint run resolves
    it. (Threshold traits are not a driving process; there is no ``discrete`` spec for them.)"""
    states = list(states)
    if len(states) < 2:
        raise ValueError(f"a discrete trait needs at least 2 states, got {states!r}")
    if len(set(states)) != len(states):
        raise ValueError(f"states must be unique, got {states!r}")
    if switch is None:
        raise ValueError("give switch= — the transition rate(s) between the discrete states.")
    return DiscreteTrait(tuple(states), switch, start, at_speciation)


