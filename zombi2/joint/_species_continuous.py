"""Species with a **diffusing** trait, grown together — QuaSSE.

The one joint model that does not race exactly. A diffusion moves at every instant, so this
run slices: each lineage's value is held fixed across a step and released at the boundary,
where the exact transition law applies."""

from __future__ import annotations

import math


from .._runtime.draw import weighted_index as _weighted_index
from .._runtime.grown import Grown
from ..params.connection import Driven
from ..params.mapping import Table
from ..species import Event as SpeciesEvent
from ..tree import Node, Tree
from ..traits import Change
from ._runaway import ONE_TRAIT as _ONE_TRAIT
from ._runaway import runaway as _runaway

#: How many slices a sliced run will walk before it decides nothing is going to happen. A sliced run
#: cannot use the other engines' "nothing scheduled, so stop" test: a rate of zero now says nothing
#: about the rate one slice later, because the driver is still moving. So the walk needs an end of
#: its own, and this is it. It is large enough that no run reaching ``n_extant`` meets it and small
#: enough to answer in seconds when the rates really are dead.
_MAX_SLICES = 2_000_000


def grow(rng, birth_rate, death_rate, trait, step: float, n_extant, total_time,
                           max_lineages=None, keys=(_ONE_TRAIT,)):
    """Grow a forward birth-death tree whose birth/death read a **continuously diffusing** trait —
    QuaSSE. Returns a `~zombi2._runtime.grown.Grown`.

    Every other joint model here races exactly, because its driver only changes at events and an
    event ends the step. A diffusion changes at every instant, so there is no interval over which
    the birth rate holds still and nothing for a Gillespie step to be drawn against.

    **This one slices.** Time is cut into steps of ``step``; inside a step every lineage's value is
    held where it was, so the rates *are* constant and the race inside the slice is the ordinary
    one. At the boundary each lineage's value is advanced by the exact transition law of its own
    diffusion — ``Normal(0, ∫σ²)`` for Brownian motion, the pull-weighted form under
    Ornstein–Uhlenbeck — over the time since that lineage last moved, which is the slice for an old
    lineage and the remainder of it for one born mid-slice. So the **trait** is exact and only its
    coupling to speciation is approximated: a lineage speciates at the rate its value had at the top
    of the slice rather than the rate it has at that instant.

    The error is first-order in ``step`` and one-sided in a way worth knowing: the value is carried
    forward, never interpolated, because a growing tree has no future to read. The conditioned
    version of the same model — a trait grown first and handed over — interpolates between two known
    ends, so for the same ``step`` it is the more accurate of the two.
    """
    from ..traits.continuous import _accrued_variance

    rate, start, theta, alpha, jump_sd = trait._resolve()
    is_ou = alpha is not None
    # What to thread the value under. A `Driven` looks its value up by `m.key`, which is the driver
    # name alone only while there is no step — with one it is `(name, step)`, because the same driver
    # read at two resolutions is two trajectories. Threading the bare name here left every factor at
    # 1.0 and the run silently undriven, which is the failure `Driven.factor`'s inert default makes
    # quiet: the tree came out identical whatever curve was written.
    lookup = set()
    for label, r in (("birth", birth_rate), ("death", death_rate)):
        for m in r.modifiers:
            if not isinstance(m, Driven) or m.driver not in keys:
                continue
            lookup.add(m.key)
            # a Curve or a Scalar is defined on the whole line, so there is no alphabet to check it
            # against; a table is the one mapping a diffusing value can never match a key of
            if isinstance(m.mapping, Table):
                raise ValueError(
                    f"{label} reads a continuous trait through a {{state: factor}} table, and a "
                    f"diffusing value is a number that never equals a state name. Map it with a "
                    f"Curve (value → factor) or a Scalar (a log-link): "
                    f'scaled_by("trait", Curve(lambda x: math.exp(0.5 * x)), step=0.05).')

    nodes: dict[int, Node] = {}
    counter = 0

    def new_node(parent, t):
        nonlocal counter
        i = counter
        counter += 1
        nodes[i] = Node(i, parent, t)
        return i

    root = new_node(None, 0.0)
    alive = [root]              # living lineage ids
    xs = [float(start)]         # each lineage's trait value, kept in lock-step with `alive`
    since = [0.0]               # …and when it last diffused, so a mid-slice birth is not over-diffused
    t = 0.0
    slice_end = step
    species_events: list[SpeciesEvent] = []
    # a diffusion cannot be rebuilt from its events, so the log carries what the trait level's own log
    # carries: the value at t=0 and each jump at a split. The values themselves ride in node_values.
    trait_events: list[Change] = [Change(0.0, "initial", root, None, float(start))]
    end_value: dict[int, float] = {}

    ceiling = None if max_lineages is None else max(max_lineages, n_extant or 0)
    slices = 0
    while alive:
        n = len(alive)
        if ceiling is not None and n > ceiling:
            raise _runaway(ceiling, t)
        ctx = {"diversity": n, "time": t}
        vals = [{key: xs[k] for key in lookup} for k in range(n)]
        wb = [birth_rate.effective(lineages=1, drivers=vals[k], **ctx) for k in range(n)]
        wd = [death_rate.effective(lineages=1, drivers=vals[k], **ctx) for k in range(n)]
        total_b, total_d = sum(wb), sum(wd)
        total = total_b + total_d

        # the slice boundary is a horizon like any other: the rates change there, so the step stops
        next_change = min(birth_rate.next_change(t), death_rate.next_change(t), slice_end)
        horizon = next_change if total_time is None else min(next_change, total_time)

        if total > 0.0:
            t_ev = t + float(rng.exponential(1.0 / total))
            if t_ev < horizon:
                t = t_ev
                if n == n_extant:
                    break
                if float(rng.random()) * total < total_b:      # speciation
                    i = _weighted_index(rng, wb, total_b)
                    node_id, x = alive[i], xs[i]
                    alive[i] = alive[-1]; alive.pop()
                    xs[i] = xs[-1]; xs.pop()
                    since[i] = since[-1]; since.pop()
                    node = nodes[node_id]
                    node.end_time = t
                    node.fate = "speciation"
                    end_value[node_id] = x
                    c1, c2 = new_node(node_id, t), new_node(node_id, t)
                    node.children = (c1, c2)
                    for c in (c1, c2):
                        d = x
                        if jump_sd > 0.0:      # the punctuational jump, drawn per daughter
                            d = x + float(rng.normal(0.0, jump_sd))
                            trait_events.append(Change(t, "on_speciation", c, x, d))
                        alive.append(c); xs.append(d); since.append(t)
                    species_events.append(SpeciesEvent(t, "speciation", node_id, (c1, c2)))
                else:                                          # extinction
                    i = _weighted_index(rng, wd, total_d)
                    node_id, x = alive[i], xs[i]
                    alive[i] = alive[-1]; alive.pop()
                    xs[i] = xs[-1]; xs.pop()
                    since[i] = since[-1]; since.pop()
                    node = nodes[node_id]
                    node.end_time = t
                    node.fate = "extinct"
                    end_value[node_id] = x
                    species_events.append(SpeciesEvent(t, "extinction", node_id))
                continue

        if total_time is not None and horizon == total_time:
            t = total_time
            break
        if horizon == slice_end:
            # the slice ends: let every living lineage diffuse over the time since it last moved
            for k in range(len(alive)):
                xs[k] = _advance(rng, rate, xs[k], since[k], slice_end, theta, alpha, is_ou,
                                 _accrued_variance)
                since[k] = slice_end
            t = slice_end
            slice_end += step
            slices += 1
            if slices > _MAX_SLICES:
                raise RuntimeError(
                    f"the run walked {_MAX_SLICES} slices of step={step:g} — to time {t:.6g} — "
                    f"without reaching n_extant={n_extant}. Either the driven birth rate is "
                    f"effectively zero over the values the trait reaches, or step is far finer than "
                    f"the timescale the tree grows on.")
            continue
        t = horizon      # a skyline breakpoint on birth/death: advance and re-read them

    for k, node_id in enumerate(alive):
        nodes[node_id].end_time = t
        # whoever is still alive reached the present, and did so with the last part-slice unrun
        end_value[node_id] = _advance(rng, rate, xs[k], since[k], t, theta, alpha, is_ou,
                                      _accrued_variance)
        nodes[node_id].fate = "extant"

    return Grown(Tree(nodes, root), species_events, trait_values=dict(end_value),
                 trait_events=trait_events)


def _advance(rng, rate, x, t0, t1, theta, alpha, is_ou, accrued):
    """One lineage's diffusion from ``t0`` to ``t1`` — the same transition law `simulate_continuous`
    walks a branch with, over a slice instead of a branch."""
    if t1 <= t0:
        return x
    if is_ou:
        e = math.exp(-alpha * (t1 - t0))
        mean = theta + (x - theta) * e
        var = accrued(rate, t0, t1, pull=alpha)
    else:
        mean, var = x, accrued(rate, t0, t1)
    return mean + (float(rng.normal(0.0, math.sqrt(var))) if var > 0.0 else 0.0)


