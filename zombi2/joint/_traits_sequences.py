"""A trait and a gene's sequence, each driving the other, on a tree the run is handed (design note
§7, the new arrow).

Both directions already run as conditioning. A trait drives a gene's substitution rate, and that
gene's composition drives how fast the trait switches. Write both at once and there is a cycle: to
evolve the sequence you would need the trait on every stretch of every branch, and the trait's own
rate is what the sequence decides.

**The walk is by time.** Species time is cut into slices of ``step``. In each slice:

1. every gene copy's composition is read, per species lineage, and held for the slice;
2. the trait is walked across the slice with its generator built from that — an ordinary Gillespie,
   so within the slice the trait's own dynamics are exact, switches mid-slice included;
3. every live gene copy is advanced across the slice at a rate reading the trait, **integrated**
   across the switches step 2 just produced.

Step 3 is exact, and not by approximation: the substitution model does not change when the trait
switches, only the rate does, and ``exp(Q·b₁)·exp(Q·b₂) = exp(Q·(b₁+b₂))``. So the branch length is
summed across the trait's segments and drawn once, which is the same distribution as drawing segment
by segment.

The one approximation is step 1 — the composition a lineage's trait rate reads belongs to the top of
the slice rather than to each instant. That is the same slicing a diffusing trait driving speciation
makes, and the same check applies: halve ``step``, rerun the same seed, and see whether the answer
moves.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from ..genomes.gene_trees import GeneNode, GeneTree
from ..sequences.evolution import _cdf_for, _sample
from ..traits.discrete import _driven_q, _gillespie
from ..traits.result import Change


@dataclass
class _Copy:
    """One gene copy in flight — the node its branch ends at, and where the walk has carried it."""

    node: GeneNode
    at: float
    states: np.ndarray
    accrued: float = 0.0


@dataclass
class _Out:
    """What the walk keeps."""

    states: dict[int, np.ndarray] = field(default_factory=dict)
    length_of: dict[int, float] = field(default_factory=dict)
    values: dict[int, object] = field(default_factory=dict)      # species node -> trait state at its end
    changes: list[Change] = field(default_factory=list)


def grow(rng, tree, *, gene_tree: GeneTree, model, length: int, founder, letters: str,
         absent: float, gene_keys, base_rate: float, gene_factors,
         trait_states, trait_entries, trait_start: int, trait_shift: float, trait_keys,
         step: float):
    """Walk the trait and the gene's sequence together, slice by slice.

    ``gene_factors`` are the ``(lookup key, mapping)`` pairs the substitution rate reads off the
    trait; ``gene_keys`` / ``trait_keys`` are the lookup keys each level is threaded under, which is
    ``(name, step)`` rather than the bare name wherever a connection declared a step.

    Returns an `_Out`: the sequence at every gene-tree node, each gene branch's length in
    substitutions per site, the trait's state at every species node, and the trait's event log.
    """
    out = _Out()
    k = len(trait_states)
    idx = tuple(model.alphabet.index(c) for c in letters)
    cache: dict[tuple[int, float], np.ndarray] = {}

    root = tree.nodes[tree.root]
    out.changes.append(Change(root.birth_time, "initial", tree.root, None, trait_states[trait_start]))
    state: dict[int, int] = {tree.root: trait_start}
    live: list[_Copy] = []
    pending_gene: GeneTree | None = gene_tree
    t = root.birth_time
    end = max(n.end_time for n in tree.nodes.values())
    order = sorted(tree.nodes)              # a parent's id precedes its children's, so preorder

    while t < end - 1e-12:
        nxt = min(t + step, end)
        if pending_gene is not None and pending_gene.origination < nxt:
            live.append(_Copy(pending_gene.complete, pending_gene.origination, founder))
            pending_gene = None

        # 1. the driver the TRAIT reads, frozen for this slice
        shares = _shares(live, idx)

        # 2. the trait, across the slice — its own Gillespie, so switches mid-slice are exact
        traj: dict[int, list[tuple[int, float, float]]] = {}
        for i in order:
            node = tree.nodes[i]
            lo, hi = max(node.birth_time, t), min(node.end_time, nxt)
            if hi <= lo:
                continue
            if i not in state:              # born inside this slice: inherit, then the optional jump
                cur = state[node.parent]
                if trait_shift > 0.0 and float(rng.random()) < trait_shift:
                    j = int(rng.integers(k - 1))
                    new = j if j < cur else j + 1
                    out.changes.append(Change(node.birth_time, "on_speciation", i,
                                              trait_states[cur], trait_states[new]))
                    cur = new
                state[i] = cur
            drivers = {key: shares.get(i, absent) for key in gene_keys}
            Q = _driven_q(trait_entries, k, drivers, lo)
            cur, segs = _gillespie(state[i], hi - lo, Q, rng)
            at, pieces = lo, []
            for st, dur in segs:
                pieces.append((st, at, at + dur))
                at += dur
            for (s1, _a, b), (s2, _c, _d) in zip(pieces, pieces[1:]):
                out.changes.append(Change(b, "on_branch", i, trait_states[s1], trait_states[s2]))
            traj[i] = pieces
            state[i] = cur

        # 3. the gene, across the slice — its rate reads the trait the step above just walked
        live = _advance(live, out, model, base_rate, gene_factors, trait_keys, trait_states, traj,
                        nxt, rng, cache)
        t = nxt

    for copy in live:                       # whatever is still in flight reached the present
        out.states[id(copy.node)] = copy.states
        out.length_of[id(copy.node)] = copy.accrued
    out.values = {i: trait_states[s] for i, s in state.items()}
    out.changes.sort(key=lambda c: c.time)
    return out


def _shares(live: list[_Copy], idx) -> dict[int, float]:
    """Each species lineage's share of the counted letters, pooled over the copies it carries."""
    hits: dict[int, int] = {}
    total: dict[int, int] = {}
    for copy in live:
        s = copy.node.species
        hits[s] = hits.get(s, 0) + int(np.isin(copy.states, idx).sum())
        total[s] = total.get(s, 0) + copy.states.size
    return {s: hits[s] / total[s] for s in total if total[s]}


def _advance(live, out: _Out, model, base: float, factors, trait_keys, trait_states, traj,
             until: float, rng, cache):
    """Carry the gene's copies to ``until``, splitting at every gene-tree node they pass.

    A copy's branch length over the slice is the trait's factor **integrated** across the segments
    the trait spent in each state — summed, then drawn once, because the matrix is the same
    throughout and only the length differs."""
    todo, out_live = list(live), []
    while todo:
        copy = todo.pop()
        stop = min(copy.node.time, until)
        if stop > copy.at:
            bl = 0.0
            for st, lo, hi in traj.get(copy.node.species, ()):
                lo, hi = max(lo, copy.at), min(hi, stop)
                if hi <= lo:
                    continue
                f = 1.0
                for key, mapping in factors:
                    f *= mapping.multiplier(trait_states[st]) if key in trait_keys else 1.0
                bl += base * f * (hi - lo)
            if bl > 0.0:
                copy.states = _sample(copy.states, _cdf_for(cache, model, bl), rng)
            copy.accrued += bl
            copy.at = stop
        if copy.at >= copy.node.time - 1e-12:
            out.states[id(copy.node)] = copy.states
            out.length_of[id(copy.node)] = copy.accrued
            for child in copy.node.children:
                todo.append(_Copy(child, copy.node.time, copy.states))
        else:
            out_live.append(copy)
    return out_live
