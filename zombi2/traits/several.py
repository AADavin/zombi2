"""Traits — several traits at once on one tree, each able to read the others (simulate_traits).

One trait grown first and then read is **conditioning**: two ordinary runs in order. This module is
for the case with no order — trait A reads trait B while B reads A — which is the trait level joined
to itself (SPEC §3). Being one level with one kind of result, it stays on the level's own function
rather than going to `zombi2.joint.simulate`.

Two engines sit behind the one door, and which one runs follows from the specs:

- **all discrete** — the pair is one Markov chain over the product of the two state spaces, walked
  exactly by the same Gillespie the single-trait engine uses (`_product_generator`). No step, no
  error.
- **anything else** — a continuous trait is one of the two, the feedback has no exact solution, and
  `zombi2.traits.stepped` holds the traits still over short steps instead.
"""

from __future__ import annotations

import numpy as np

from ..params.connection import Driven
from ..params.parameter import as_rate
from ..params.scope import PerLineage
from ..rng import stream
from ..tree import as_tree

from ._shared import _preorder
from .continuous import ContinuousTrait
from .discrete import DiscreteTrait, _driven_entries, _gillespie
from .result import Change, TraitsResult


def _driver_keys(spec):
    """What a rate may call this trait: its name, and the bare word when a run holds one of it."""
    return ("trait",) if spec.name is None else ("trait", f"traits:{spec.name}")


def _readings(spec) -> list[tuple[object, str]]:
    """Every live-driver reading a spec carries, as ``(what is read, where it was written)``.

    ``what is read`` is the driver's name. ``where`` names the parameter it was written on, so a
    refusal can say which line to look at. A discrete trait reads through its switch rates; a
    continuous one through σ², through an optimum written with ``set_by``, and through ``regimes=``,
    which names a discrete trait rather than carrying a modifier."""
    out: list[tuple[object, str]] = []
    if isinstance(spec, ContinuousTrait):
        for m in as_rate(spec.rate, default_scope=PerLineage).modifiers:
            if isinstance(m, Driven) and isinstance(m.driver, str):
                out.append((m.driver, "rate"))
        if isinstance(spec.reverts_to, Driven) and isinstance(spec.reverts_to.driver, str):
            out.append((spec.reverts_to.driver, "reverts_to"))
        if isinstance(spec.regimes, str):
            out.append((spec.regimes, "regimes"))
    else:
        for entry in _driven_entries(list(spec.states), spec.switch):
            for m in entry[2].modifiers:
                if isinstance(m, Driven) and isinstance(m.driver, str):
                    out.append((m.driver, "switch"))
    return out


def simulate_traits(tree, traits, *, joint=False, seed=None, progress=False):
    """Evolve **several traits at once** along a fixed tree, each one able to read the others.

    One trait is `~zombi2.traits.discrete.simulate_discrete` /
    `~zombi2.traits.continuous.simulate_continuous`, and a trait grown first and then read is
    conditioning — two ordinary runs. This is for the case with no order: trait A's rate reads trait
    B while B's reads A, so neither can be finished before the other starts. That is the trait level
    joined to itself (SPEC §3), and being one level with one kind of result it stays here rather
    than going to `zombi2.joint.simulate`.

    ``traits`` is a list of process specs — `~zombi2.traits.discrete.discrete` or
    `~zombi2.traits.continuous.continuous` — each with a ``name``, and a rate reads another by
    ``"traits:<name>"``. ``joint=True`` says the run is what it is, and is checked both ways: asking
    for it when no trait reads another is an error, and reading another without it is an error too.

    **Two discrete traits** are one Markov chain over the product of their state spaces, walked
    exactly::

        traits.simulate_traits(tree, [
            traits.discrete(name="habitat", states=["surface", "cave"],
                            switch={"surface->cave": PerLineage(0.05).scaled_by(
                                        "traits:size", {"small": 1.0, "large": 8.0}),
                                    "cave->surface": 0.1}),
            traits.discrete(name="size", states=["small", "large"], switch=...)], joint=True, seed=1)

    **A continuous trait in the pair** has no exact solution, so the run holds the traits still over
    short steps (`zombi2.traits.stepped`, which holds the method and the error). ``step`` rides on
    the connection, never on this call. Two models run there. A **curved optimum**, each trait's
    optimum a curved function of the other::

        traits.simulate_traits(tree, [
            traits.continuous(name="brain", start=0.0, rate=0.1, pull=1.0,
                              reverts_to=set_by("traits:group", lambda v: 2.0 * math.tanh(v),
                                                step=0.01)),
            traits.continuous(name="group", start=0.0, rate=0.1, pull=0.5,
                              reverts_to=set_by("traits:brain", lambda v: 1.5 * math.tanh(0.7 * v),
                                                step=0.01))], joint=True, seed=1)

    and **a continuous trait with a discrete one**, where the discrete trait paints the continuous
    trait's optimum and the continuous trait's value drives the discrete trait's switch rate::

        traits.simulate_traits(tree, [
            traits.discrete(name="habitat", states=["surface", "cave"],
                            switch=PerLineage(0.2).scaled_by("traits:size", Scalar(0.4),
                                                             step=0.02)),
            traits.continuous(name="size", start=0.0, rate=0.4, pull=1.0,
                              regimes="traits:habitat",
                              reverts_to={"surface": 0.0, "cave": 2.0})], joint=True, seed=1)

    Straight-line optima are **refused** rather than run here: they are one multivariate
    Ornstein–Uhlenbeck with a full drift matrix, which `simulate_continuous` solves exactly, and the
    refusal names the call to write. Two traits setting each other's σ² are refused too; that one is
    not built.

    ``at_speciation`` works here as it does in a run of its own: each trait jumps on its own at the
    split, and the pair lands wherever the jumps leave it.

    Returns ``{name: TraitsResult}`` — one complete result per trait, exactly what the single-trait
    runners return, so every reader of one works on these unchanged. Deterministic given ``seed``.
    Where a continuous trait is in the pair, halving ``step`` **changes the numbers** for a given
    seed as well as the error, so convergence is checked across seeds; `zombi2.traits.stepped` holds
    the recipe.
    """
    from . import stepped

    tree = as_tree(tree, level="traits")
    specs = list(traits)
    if len(specs) < 2:
        raise ValueError(
            f"simulate_traits evolves several traits at once, and got {len(specs)}. One trait is "
            f"simulate_discrete(tree, ...) or simulate_continuous(tree, ...); two that read each "
            f"other are what this is for.")
    if len(specs) > 2:
        raise NotImplementedError(
            f"two traits reading each other is what this evolves today, and got {len(specs)}. "
            f"Three is a scope choice now rather than a cost: the product generator two discrete "
            f"traits use grows exponentially in the number of traits, but the stepped engine is "
            f"linear in it, so nothing about three traits is expensive. It is simply not built.")
    for spec in specs:
        if not isinstance(spec, (DiscreteTrait, ContinuousTrait)):
            raise TypeError(
                f"simulate_traits takes trait process specs — traits.discrete(name='habitat', "
                f"states=[...], switch=...) or traits.continuous(name='size', rate=...) — and got "
                f"{spec!r}. A finished result is a driver you already have, which is conditioning: "
                f"pass it to the driven trait's own run instead.")
        if not spec.name:
            raise ValueError(
                "each trait needs a name here, because a rate reads the other one by it: "
                "traits.discrete(name='habitat', ...) and scaled_by('traits:habitat', ...).")
    if len({s.name for s in specs}) != len(specs):
        raise ValueError(f"trait names must be unique, got {[s.name for s in specs]}")

    names = {k for s in specs for k in _driver_keys(s)}
    all_discrete = all(isinstance(s, DiscreteTrait) for s in specs)
    reads = 0
    for spec in specs:
        for driver, where in _readings(spec):
            if driver not in names:
                raise ValueError(
                    f"trait {spec.name!r} reads {driver!r} on its {where}, which is not a trait in "
                    f"this run. The traits here are {sorted(s.name for s in specs)}, read as "
                    f'"traits:<name>".')
            if driver == "trait" and not all_discrete:
                raise ValueError(
                    f"trait {spec.name!r} reads {driver!r} on its {where}, and this run holds two "
                    f"traits, so the bare word does not say which. Name it: "
                    f'"traits:<name>".')
            reads += 1
    if reads and not joint:
        raise ValueError(
            "a trait here reads another trait in this same run, so the two are joint — neither can "
            "be finished before the other starts. Say so with joint=True. To read a trait grown "
            "EARLIER, pass that run's result rather than a name, which is conditioning.")
    if joint and not reads:
        raise ValueError(
            "joint=True says the traits drive each other, but none reads another. Give a rate a "
            'scaled_by("traits:<name>", ...) or an optimum a set_by("traits:<name>", ...), or '
            "evolve them as separate runs.")

    if not all_discrete:
        return stepped.simulate(tree, specs, seed=seed, progress=progress)
    return _product_run(tree, specs, seed=seed, progress=progress)


# --- two discrete traits: one chain over the product of their states -------------------------------

def _product_generator(specs, resolved):
    """The generator of the **pair**, over every combination of the two traits' states.

    Two traits that read each other are one Markov chain on the product of their state spaces, and
    that is not an approximation of the pair — it *is* the pair. From ``(i, j)`` the only moves are
    to ``(i', j)`` at trait A's rate, read with B sitting in ``j``, and to ``(i, j')`` at B's rate
    read with A sitting in ``i``. Nothing moves both at once, because two switches never coincide.

    Because each rate depends only on the *other* trait's current state, the whole matrix can be
    built once and handed to `~zombi2.traits.discrete._gillespie` — the same exact branch walk a
    single trait takes. What makes the pair joint is that neither column of it can be filled in
    without the other."""
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


def _product_run(tree, specs, *, seed, progress):
    """Two discrete traits, walked exactly as one chain over the product of their states."""
    # each trait's own alphabet and its rate specs, left unsettled: a switch rate reading the other
    # trait is not one number, which is exactly what `_q_matrix` would demand
    resolved = [(list(s.states), _driven_entries(list(s.states), s.switch)) for s in specs]

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
