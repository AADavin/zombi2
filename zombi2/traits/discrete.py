"""Traits — the discrete engine: the Mk state-switching model (Gillespie, its switch rate optionally driven by another level) and the threshold model (simulate_discrete), plus the DiscreteTrait process spec for joint runs."""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from ..params.mapping import check_not_a_kernel
from ..rng import stream
from ..params.evaluate import Modifier, describe, is_implemented
from ..params.connection import Driven
from ..params.parameter import Rate, as_rate
from ..params.scope import PerLineage
from ..tree import as_tree

from ._shared import _correlation_matrix, _preorder, _resolve_drivers, _symmetric_sqrt
from .result import Change, TraitsResult

#: the modifiers a discrete trait's ``switch`` rate takes — declared, as every level declares its
#: own, so Appendix A and the CLI help are built from the engine rather than kept by hand. It is a
#: shorter list than the continuous rate's: a switch rate reads a driver and nothing else.
IMPLEMENTED_MODIFIERS = (Driven,)


def _switch_rate(spec, where: str) -> Rate:
    """``spec`` as a resolved `Rate`, refusing the scopes a switch rate cannot have.

    A switch rate is written from its scope like every other rate in the library —
    ``PerLineage(0.4)``, or a bare ``0.4``, which means the same thing. `PerLineage` is the only
    scope that means anything here: a lineage's state changes on its own clock, and there is no
    larger pool of chances for a wider scope to count over. Both runners coerce through this one
    function, so `simulate_discrete()` on a fixed tree and `DiscreteTrait` inside a joint run cannot
    answer the same spec differently.

    ``where`` names the rate in any refusal — "a switch rate", or the transition it belongs to."""
    r = as_rate(spec, default_scope=PerLineage, label=where)
    assert r.scope is not None          # `as_rate` filled the level's default just above
    if r.scope is not PerLineage:
        raise ValueError(
            f"{where} has a {r.scope.__name__} scope, but a discrete trait switches per lineage: "
            f"each lineage's state changes on its own clock, so there is no wider pool of chances "
            f"to count over. Write PerLineage(...), or a bare number, which is per lineage here.")
    return r


def _switch_number(spec, where: str) -> float:
    """The number ``spec`` settles to, for a generator that is one constant matrix.

    A rate written from its scope is accepted and reduced to the number it stands for, because
    ``PerLineage(0.4)`` and ``0.4`` are the same rate and the first is how the grammar teaches every
    other rate to be written. A rate carrying a verb settles to no single number, so it is refused
    here and handled by `_driven_entries()` instead."""
    r = _switch_rate(spec, where)
    if r.modifiers:
        carried = ", ".join(dict.fromkeys(describe(m) for m in r.modifiers))
        raise ValueError(
            f"{where} carries {carried}, and the generator built here is one constant matrix, so "
            f"there is nothing for it to read: a driven switch rate needs the driver's own "
            f"trajectory along each branch. Grow the driver first and pass the finished tree to "
            f"simulate_discrete(tree, ...), which rebuilds the generator wherever the driver moves.")
    return r.effective(lineages=1)


def _q_matrix(states, switch) -> np.ndarray:
    """Build the ``k×k`` transition-rate matrix ``Q`` from ``switch`` — the CTMC generator whose
    off-diagonal ``Q[i, j] ≥ 0`` is the rate ``state i → state j``. ``switch`` is one of:

    - a **rate** — the symmetric equal-rates shortcut: every ``i → j`` (``i ≠ j``) at that rate.
      ``0.4`` or ``PerLineage(0.4)``, which are the same rate written two ways;
    - a ``{"from->to": rate}`` **dict** — only the named transitions, others zero (asymmetric),
      each value a rate in those same two spellings;
    - a ``k×k`` **matrix** — the off-diagonal rates directly, as plain numbers (the diagonal is
      ignored). Every entry is counted per lineage, so a cell has no scope of its own to write.

    Each diagonal is then set to minus its row sum, so rows sum to zero (a proper generator)."""
    k = len(states)
    idx = {s: i for i, s in enumerate(states)}
    Q = np.zeros((k, k))
    if isinstance(switch, Rate):
        switch = _switch_number(switch, "a switch rate")
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
            if isinstance(rate, Rate):
                rate = _switch_number(rate, f"the switch rate for {key!r}")
            if isinstance(rate, bool) or not isinstance(rate, (int, float)) \
                    or not math.isfinite(rate) or rate < 0:
                raise ValueError(f"switch rate for {key!r} must be finite and non-negative, got {rate!r}")
            Q[idx[frm], idx[to]] = float(rate)
    elif isinstance(switch, (list, tuple, np.ndarray)):
        try:
            arr = np.asarray(switch, dtype=float)
        except TypeError:  # a rate object in a cell: numpy's own message is about float(), not switch
            raise ValueError(
                "a switch matrix holds plain numbers — every entry is counted per lineage, so a "
                "cell has no scope of its own to write. Give the {'from->to': rate} dict instead "
                "when an entry is a rate."
            ) from None
        if arr.shape != (k, k):
            raise ValueError(f"switch matrix must be {k}×{k} for {k} states, got shape {arr.shape}")
        Q = arr.copy()
        np.fill_diagonal(Q, 0.0)
        if np.any(Q < 0) or not np.all(np.isfinite(Q)):
            raise ValueError("switch matrix off-diagonals must be finite and non-negative")
    else:
        raise ValueError(
            "switch must be a rate — 0.4, or PerLineage(0.4) written out — for the symmetric case, "
            "a {'from->to': rate} dict for named transitions, or a k×k matrix of numbers."
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


def _switch_modifiers(switch) -> list:
    """Every modifier the switch rates carry — a switch rate written as a rate expression
    (``PerLineage(0.4).scaled_by(habitat, {"aquatic": 3.0})``) rather than as a bare number, so the
    trait switches faster on the lineages where the driver is in one state than another.

    ``scaled_by`` is the only verb this engine ships for a switch rate; anything else built in
    would be read by no part of it, so it is refused by name rather than silently ignored. A
    third-party modifier that named this engine in its `Modifier.implemented_for` is returned here
    too — carrying *any* modifier is what puts the run on the rebuild-per-stretch path below, where
    the generator is a function of the context instead of a constant."""
    mods = []
    for spec in _switch_specs(switch):
        if not isinstance(spec, (Rate, Modifier)):
            continue                       # a bare number (or a matrix): nothing to drive
        r = _switch_rate(spec, "a switch rate")
        for m in r.modifiers:
            if not is_implemented(m, IMPLEMENTED_MODIFIERS, "traits.discrete"):
                raise ValueError(
                    f"a switch rate carries {describe(m)}, which the discrete trait engine does "
                    f"not support. It takes scaled_by (the switch rate driven by another level).")
            if isinstance(m, Driven):     # only a driver has a mapping to check
                check_not_a_kernel(m.mapping, label="a switch rate")
            mods.append(m)
    return mods


def _driven_entries(states, switch) -> list:
    """The off-diagonal switch rates as ``[(i, j, Rate)]`` — the driven twin of `_q_matrix()`, whose
    entries are rate *specs* rather than settled numbers, so the generator can be rebuilt wherever
    the driver switches. The same two shapes a driven rate can be written in: one symmetric rate
    (every ``i → j`` alike) or a ``{'from->to': rate}`` dict (only the named transitions, the rest
    zero). Validation mirrors `_q_matrix()`'s, the rate itself being checked by `_switch_rate()`."""
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
            entries.append((idx[frm], idx[to], _switch_rate(spec, f"the switch rate for {key!r}")))
        return entries
    r = _switch_rate(switch, "a switch rate")
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
    """Exact CTMC along one branch whose switch rates are **driven** by another level.

    The driver is piecewise-constant along the branch, so the generator is too: the branch is cut at
    the driver's own switches (`~zombi2.params.conditioned.DriverTrajectory.next_change`) and the plain
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
        # The stretch ends at whichever comes first: a driver switching, or a rate changing on its
        # own clock. Only the drivers were asked, so a modifier that varies with time — a skyline, or
        # one of your own admitted by `implemented_for` — was read once at the start of the stretch
        # and then held for the rest of the branch, which every other engine steps to.
        nxt = t1
        for _i, _j, r in entries:
            nxt = min(nxt, r.next_change(t))
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
    if r.scope is not PerLineage or r.modifiers:
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
      (``0.1``, or ``PerLineage(0.1)`` written out), a ``{"marine->terrestrial": 0.1}`` dict, or a
      ``k×k`` matrix of numbers (see `_q_matrix()`).
      ``start`` is the root state (a label in ``states``; ``None`` draws one uniformly). A switch rate
      may be **driven by another level** grown first on this same tree — write it as a rate
      expression, ``switch=PerLineage(0.4).scaled_by(habitat, {"aquatic": 3.0})`` or per transition,
      ``switch={"a->b": PerLineage(0.2).scaled_by(habitat, {"aquatic": 3.0}), "b->a": 0.2}``. The driver
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
    # conditioning: a switch rate written with scaled_by reads another level, grown first on this same
    # tree. The generator is then a function of the driver, so it is built per stretch rather than
    # once. No modifier at all ⇒ one constant Q and the walk below is exactly the walk it was.
    #
    # The two questions are separate, and conflating them was a bug: *any* modifier puts the run on
    # the rebuild-per-stretch path (that is what makes the generator a function of the context),
    # while only a `Driven` names a driver to resolve a trajectory from. A third-party
    # modifier used to pass the gate, land in the driver list, and crash the resolver looking for a
    # `.key` it does not have.
    sw_mods = _switch_modifiers(switch)
    entries = _driven_entries(states, switch) if sw_mods else None
    Q = None if sw_mods else _q_matrix(states, switch)
    trajs = _resolve_drivers([m for m in sw_mods if isinstance(m, Driven)], tree, "traits.discrete")
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
    # the driver on every lineage, so the event log is the driver file (no separate driver).
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
    (``joint.simulate(species.birth_death(...), traits.discrete(...))``), so neither can be simulated
    first. Same
    parameters as `simulate_discrete()` (the Mk half): ``states``, ``switch`` (the rate spec),
    ``start`` (the root state, ``None`` = uniform), ``at_speciation`` (the on-speciation shift
    probability)."""

    states: tuple
    switch: object
    start: object = None
    at_speciation: object = None
    #: what a driver calls this trait — ``scaled_by("traits:size", …)``. A run with one trait in it
    #: needs no name and takes ``"trait"``; naming is what lets a run hold two.
    name: "str | None" = None

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



def discrete(*, states, switch=None, start=None, at_speciation=None, name=None) -> DiscreteTrait:
    """A discrete-trait **process spec** — `DiscreteTrait`, unexecuted — for a joint model to
    simulate with the tree it drives (``joint.simulate(species.birth_death(...),
    traits.discrete(states=[...], switch=...))``).
    A thin bundle of `simulate_discrete()`'s Mk parameters; validated when the joint run resolves
    it. (Threshold traits are not a driving process; there is no ``discrete`` spec for them.)"""
    states = list(states)
    if len(states) < 2:
        raise ValueError(f"a discrete trait needs at least 2 states, got {states!r}")
    if len(set(states)) != len(states):
        raise ValueError(f"states must be unique, got {states!r}")
    if switch is None:
        raise ValueError("give switch= — the transition rate(s) between the discrete states.")
    if name is not None and (not isinstance(name, str) or not name.strip()):
        raise ValueError(f"a trait's name is what a driver calls it, so it must be a non-empty "
                         f"string; got {name!r}")
    return DiscreteTrait(tuple(states), switch, start, at_speciation, name)




# --- two traits, each reading the other: the trait level joined to itself --------------------------

def _product_generator(specs, resolved):
    """The generator of the **pair**, over every combination of the two traits' states.

    Two traits that read each other are one Markov chain on the product of their state spaces, and
    that is not an approximation of the pair — it *is* the pair. From ``(i, j)`` the only moves are
    to ``(i', j)`` at trait A's rate, read with B sitting in ``j``, and to ``(i, j')`` at B's rate
    read with A sitting in ``i``. Nothing moves both at once, because two switches never coincide.

    Because each rate depends only on the *other* trait's current state, the whole matrix can be
    built once and handed to `_gillespie()` — the same exact branch walk a single trait takes. What
    makes the pair joint is that neither column of it can be filled in without the other."""
    (a_states, a_entries), (b_states, b_entries) = resolved
    ka, kb = len(a_states), len(b_states)
    n = ka * kb
    Q = np.zeros((n, n))
    at = lambda i, j: i * kb + j
    for i in range(ka):
        for j in range(kb):
            row = at(i, j)
            # trait A moves, reading B's state right now
            drivers_b = {k: b_states[j] for k in _driver_keys(specs[1])}
            for x, y, r in a_entries:
                if x == i:
                    Q[row, at(y, j)] += r.effective(lineages=1, drivers=drivers_b)
            # trait B moves, reading A's
            drivers_a = {k: a_states[i] for k in _driver_keys(specs[0])}
            for x, y, r in b_entries:
                if x == j:
                    Q[row, at(i, y)] += r.effective(lineages=1, drivers=drivers_a)
    np.fill_diagonal(Q, 0.0)
    np.fill_diagonal(Q, -Q.sum(axis=1))
    return Q


def _driver_keys(spec):
    """What a rate may call this trait: its name, and the bare word when a run holds one of it."""
    return ("trait",) if spec.name is None else ("trait", f"traits:{spec.name}")


def simulate_traits(tree, traits, *, joint=False, seed=None, progress=False):
    """Evolve **several traits at once** along a fixed tree, each one able to read the others.

    One trait is `simulate_discrete()` / `simulate_continuous()`, and a trait grown first and then
    read is conditioning — two ordinary runs. This is for the case with no order: trait A's switch
    rate reads trait B while B's reads A, so neither can be finished before the other starts. That is
    the trait level joined to itself (SPEC §3), and being one level with one kind of result it stays
    here rather than going to `zombi2.joint.simulate`.

    ``traits`` is a list of `discrete()` specs, each with a ``name``, and a rate reads another by
    ``scaled_by("traits:<name>", …)``. ``joint=True`` says the run is what it is, and is checked both
    ways: asking for it when no trait reads another is an error, and reading another without it is an
    error too.

    ``at_speciation`` works here as it does in a run of its own: each trait carrying one hops on its
    own at the split, and the pair lands wherever the two hops leave it.

    Returns ``{name: TraitsResult}`` — one complete result per trait, exactly what the single-trait
    runners return, so every reader of one works on these unchanged. Deterministic given ``seed``.
    """
    from ..params.connection import Driven

    tree = as_tree(tree, level="traits")
    specs = list(traits)
    if len(specs) < 2:
        raise ValueError(
            f"simulate_traits evolves several traits at once, and got {len(specs)}. One trait is "
            f"simulate_discrete(tree, ...); two that read each other are what this is for.")
    if len(specs) > 2:
        raise NotImplementedError(
            f"two traits reading each other is what this evolves today, and got {len(specs)}. Three "
            f"is the same idea over a bigger product of states, and is not built.")
    for spec in specs:
        if not isinstance(spec, DiscreteTrait):
            raise TypeError(
                f"simulate_traits takes discrete trait specs — traits.discrete(name='size', "
                f"states=[...], switch=...) — and got {spec!r}. A continuous trait in a cycle needs "
                f"its diffusion held still over short stretches, which is not built.")
        if not spec.name:
            raise ValueError(
                "each trait needs a name here, because a rate reads the other one by it: "
                "traits.discrete(name='habitat', ...) and scaled_by('traits:habitat', ...).")
    if len({s.name for s in specs}) != len(specs):
        raise ValueError(f"trait names must be unique, got {[s.name for s in specs]}")

    # each trait's own alphabet and its rate specs, left unsettled: a switch rate reading the other
    # trait is not one number, which is exactly what `_q_matrix` would demand
    resolved = [(list(s.states), _driven_entries(list(s.states), s.switch)) for s in specs]
    names = {k for s in specs for k in _driver_keys(s)}
    reads = 0
    for spec, other in ((specs[0], specs[1]), (specs[1], specs[0])):
        for entry in _driven_entries(list(spec.states), spec.switch):
            for m in entry[2].modifiers:
                if isinstance(m, Driven):
                    if m.driver not in names:
                        raise ValueError(
                            f"trait {spec.name!r} reads {m.driver!r}, which is not a trait in this "
                            f"run. The traits here are {sorted(s.name for s in specs)}, read as "
                            f'scaled_by("traits:<name>", ...).')
                    reads += 1
    if reads and not joint:
        raise ValueError(
            "a trait here reads another trait in this same run, so the two are joint — neither can "
            "be finished before the other starts. Say so with joint=True. To read a trait grown "
            "EARLIER, pass that run's result rather than a name, which is conditioning.")
    if joint and not reads:
        raise ValueError(
            "joint=True says the traits drive each other, but none reads another. Give a switch "
            'rate a scaled_by("traits:<name>", ...), or evolve them as separate runs.')

    rng, seed = stream("traits", seed)
    Q = _product_generator(specs, resolved)
    (a_states, _a), (b_states, _b) = resolved
    kb = len(b_states)
    starts = []
    for spec, (states, _e) in zip(specs, resolved):
        idx = {s: i for i, s in enumerate(states)}
        if spec.start is None:
            starts.append(int(rng.integers(len(states))))
        elif spec.start in idx:
            starts.append(idx[spec.start])
        else:
            raise ValueError(f"start must be one of states={states} (or None for a uniform draw), "
                             f"got {spec.start!r}")
    start = starts[0] * kb + starts[1]

    root = tree.nodes[tree.root]
    node_pairs: dict[int, int] = {}
    per_trait: list[list] = [
        [Change(root.birth_time, "initial", tree.root, None, a_states[starts[0]])],
        [Change(root.birth_time, "initial", tree.root, None, b_states[starts[1]])]]
    shifts = [0.0 if s.at_speciation is None else float(s.at_speciation) for s in specs]
    for i in _preorder(tree, progress):
        node = tree.nodes[i]
        cur = start if node.parent is None else node_pairs[node.parent]
        if node.parent is not None and any(shifts):
            # each trait hops on its own at the split, exactly as it would in a run of its own; the
            # pair simply lands wherever the two hops leave it
            parts = [cur // kb, cur % kb]
            for k, (shift, (states, _e)) in enumerate(zip(shifts, resolved)):
                if shift > 0.0 and float(rng.random()) < shift:
                    j = int(rng.integers(len(states) - 1))   # to a uniform *other* state
                    new = j if j < parts[k] else j + 1
                    per_trait[k].append(Change(node.birth_time, "on_speciation", i,
                                               states[parts[k]], states[new]))
                    parts[k] = new
            cur = parts[0] * kb + parts[1]
        end, segs = _gillespie(cur, node.end_time - node.birth_time, Q, rng)
        # one product move is one trait switching, so unpacking the segments splits the pair's
        # history back into the two the reader asked for, with no state left ambiguous
        t = node.birth_time
        for (s1, d1), (s2, _d) in zip(segs, segs[1:]):
            t += d1
            for k, (states, changes) in enumerate(((a_states, per_trait[0]),
                                                   (b_states, per_trait[1]))):
                was, now = (s1 // kb, s2 // kb) if k == 0 else (s1 % kb, s2 % kb)
                if was != now:
                    changes.append(Change(t, "on_branch", i, states[was], states[now]))
        node_pairs[i] = end
    out = {}
    for k, (spec, (states, _e)) in enumerate(zip(specs, resolved)):
        values = {i: states[(p // kb) if k == 0 else (p % kb)] for i, p in node_pairs.items()}
        per_trait[k].sort(key=lambda c: c.time)
        out[spec.name] = TraitsResult(tree, values, per_trait[k], seed, kind="discrete")
    return out
