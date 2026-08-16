"""Species with a discrete trait, grown together — BiSSE and MuSSE.

The tree comes out of the run. Birth and death read the trait state on each lineage, while
the trait switches by its own Mk process on the tree it is shaping. The race is exact: the
trait changes only at events, and an event ends the step."""

from __future__ import annotations

import math


from .._runtime.draw import weighted_index as _weighted_index
from .._runtime.grown import Grown
from ..params.connection import Driven
from ..species import Event as SpeciesEvent
from ..tree import Node, Tree
from ..traits import Change, DiscreteTrait
from ._runaway import ONE_TRAIT as _ONE_TRAIT
from ._runaway import runaway as _runaway

def grow(rng, birth_rate, death_rate, trait: DiscreteTrait, n_extant, total_time,
                max_lineages=None, keys=(_ONE_TRAIT,)):
    """Grow a forward birth-death tree whose birth/death read a discrete trait that evolves on it.
    Returns a `~zombi2._runtime.grown.Grown`: the complete tree, the speciation/extinction log, the
    trait state at every node, and the trait's switch log."""
    states, Q, start_i, shift = trait._resolve(rng)
    k_states = len(states)
    out_rate = [float(-Q[s, s]) for s in range(k_states)]  # the trait's total switch-out rate per state

    # birth/death are driven by the trait; the trait's declared states are known up front, so the check
    # is exhaustive — every mapping key must be one of them. A key outside the alphabet is a state that
    # can never occur (a typo whose factor would silently never apply), and a mapping matching none of
    # them would be a silently undriven run; both are refused here rather than run as if driven.
    from ..params.conditioned import check_mapping_fires
    for label, rate in (("birth", birth_rate), ("death", death_rate)):
        for m in rate.modifiers:
            if isinstance(m, Driven):
                check_mapping_fires(m.mapping, states, driver_label=f"{label} (trait)", exhaustive=True)

    nodes: dict[int, Node] = {}
    counter = 0

    def new_node(parent, t):
        nonlocal counter
        i = counter
        counter += 1
        nodes[i] = Node(i, parent, t)
        return i

    root = new_node(None, 0.0)
    alive = [root]      # living lineage ids
    st = [start_i]      # each lineage's trait state index, kept in lock-step with `alive`
    t = 0.0
    species_events: list[SpeciesEvent] = []
    # the initial state at t=0, exactly as the standalone traits engine seeds its log: tree + this row
    # + the switches give the trait on every lineage, which is what lets the log be read back as a
    # driver. Without it a joint run writes a trait_events.tsv no conditioned run can read.
    trait_events: list[Change] = [Change(0.0, "initial", root, None, states[start_i])]
    end_state: dict[int, int] = {}  # node id → its trait state index when it ended (→ node_values)

    ceiling = None if max_lineages is None else max(max_lineages, n_extant or 0)
    while alive:
        n = len(alive)
        if ceiling is not None and n > ceiling:
            raise _runaway(ceiling, t)
        ctx = {"diversity": n, "time": t}
        # per-lineage rates: birth/death read the lineage's trait state (scaled_by("trait", …)); the
        # trait switch rate is the CTMC out-rate for that state (the trait's own dynamics, undriven).
        vals = [{key: states[st[k]] for key in keys} for k in range(n)]
        wb = [birth_rate.effective(lineages=1, drivers=vals[k], **ctx) for k in range(n)]
        wd = [death_rate.effective(lineages=1, drivers=vals[k], **ctx) for k in range(n)]
        ws = [out_rate[st[k]] for k in range(n)]
        total_b, total_d, total_s = sum(wb), sum(wd), sum(ws)
        total = total_b + total_d + total_s

        # the trait switch rate is constant between events; only a skyline (changing_at) on birth/death or
        # the total_time limit advances the clock on its own.
        next_change = min(birth_rate.next_change(t), death_rate.next_change(t))
        horizon = next_change if total_time is None else min(next_change, total_time)

        if total > 0.0:
            t_ev = t + float(rng.exponential(1.0 / total))
            if t_ev < horizon:
                t = t_ev
                if n == n_extant:  # already at the target; stop at this next event's time, unapplied
                    break
                r = float(rng.random()) * total
                if r < total_b:  # speciation
                    i = _weighted_index(rng, wb, total_b)
                    node_id, cur = alive[i], st[i]
                    alive[i] = alive[-1]; alive.pop()          # swap-remove keeps the state array in step
                    st[i] = st[-1]; st.pop()
                    node = nodes[node_id]
                    node.end_time = t
                    node.fate = "speciation"
                    end_state[node_id] = cur
                    c1, c2 = new_node(node_id, t), new_node(node_id, t)
                    node.children = (c1, c2)
                    for c in (c1, c2):  # each daughter inherits the parent's state (+ optional split shift)
                        d = cur
                        if shift > 0.0 and float(rng.random()) < shift:
                            j = int(rng.integers(k_states - 1))  # hop to a uniform *other* state
                            d = j if j < cur else j + 1
                            trait_events.append(Change(t, "on_speciation", c, states[cur], states[d]))
                        alive.append(c); st.append(d)
                    species_events.append(SpeciesEvent(t, "speciation", node_id, (c1, c2)))
                elif r < total_b + total_d:  # extinction
                    i = _weighted_index(rng, wd, total_d)
                    node_id, cur = alive[i], st[i]
                    alive[i] = alive[-1]; alive.pop()
                    st[i] = st[-1]; st.pop()
                    node = nodes[node_id]
                    node.end_time = t
                    node.fate = "extinct"
                    end_state[node_id] = cur
                    species_events.append(SpeciesEvent(t, "extinction", node_id))
                else:  # trait switch — change one lineage's state, no topology change
                    i = _weighted_index(rng, ws, total_s)
                    node_id, cur = alive[i], st[i]
                    probs = Q[cur].copy()
                    probs[cur] = 0.0
                    probs /= out_rate[cur]          # the embedded jump chain: where to, given a jump
                    new = int(rng.choice(k_states, p=probs))
                    st[i] = new
                    trait_events.append(Change(t, "on_branch", node_id, states[cur], states[new]))
                continue

        if math.isinf(horizon):
            break  # nothing scheduled and no skyline change → nothing more can happen
        if total_time is not None and horizon == total_time:
            t = total_time
            break
        t = horizon  # a skyline breakpoint: advance and re-evaluate the (now changed) birth/death

    for k, node_id in enumerate(alive):  # whoever is still alive reached the present
        nodes[node_id].end_time = t
        nodes[node_id].fate = "extant"
        end_state[node_id] = st[k]

    node_values = {i: states[end_state[i]] for i in nodes}
    return Grown(Tree(nodes, root), species_events, trait_values=node_values,
                 trait_events=trait_events)


