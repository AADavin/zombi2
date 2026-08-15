"""The sequence level joined to **itself** — two genes, each one's substitution rate reading the
other's composition, in one run (design note §6, §8, §9).

Every other run at this level evolves one gene tree at a time: a family is walked from its
origination to its tips, then the next family starts. That order is impossible here. If the
chaperone's rate reads the client's composition and the client's reads the chaperone's, then neither
family can be finished before the other starts, and a family walked to its tips would have to read a
driver that does not exist yet.

So the walk is **by time rather than by family**. Species time is cut into slices of ``step``, and
inside a slice every living gene copy of every family advances by the same stretch, at a rate read
off the *other* family's composition as it stood at the top of the slice. At the boundary the
compositions are recomputed and the next slice begins.

That is the approximation, and it is the same one a diffusing trait driving speciation makes
(`zombi2.joint`): a composition moves with every substitution, so there is no interval over which
either rate holds still and nothing exact to draw against. Inside a slice the rate is constant, so
the transition matrix is the ordinary one; what is approximated is that the rate a copy evolves at
belongs to the top of its slice rather than to each instant. Halve ``step``, rerun the same seed, and
see whether the answer moves.

A composition is read **per species lineage**, pooled over whatever copies of that family the
lineage carries — the same statistic `zombi2.sequences._composition.Composition` computes off a
finished run, computed here as the run goes. A lineage carrying none of the family reads the
``absent`` the gene declared.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np

from ..genomes.gene_trees import GeneNode, GeneTree
from .evolution import _cdf_for, _sample


@dataclass
class _Live:
    """One gene copy in flight: the node its branch ends at, where that branch started, and the
    states it holds at ``at`` — the time the walk has carried it to so far."""

    node: GeneNode
    at: float
    states: np.ndarray
    accrued: float = 0.0        # substitutions/site down this branch, summed over the slices it spans


@dataclass
class _Family:
    """One family's walk: its tree, its model, and everything the run keeps from it."""

    gene: object                                  # the GeneSpec
    tree: GeneTree
    live: list[_Live] = field(default_factory=list)
    states: dict[int, np.ndarray] = field(default_factory=dict)   # id(node) -> states at that node
    length_of: dict[int, float] = field(default_factory=dict)     # id(node) -> its branch, in subs/site
    founding: np.ndarray | None = None


def _composition_by_lineage(fam: _Family, letters: tuple[int, ...], absent: float) -> dict[int, float]:
    """Each species lineage's share of ``letters``, pooled over the copies of this family it is
    carrying right now. A lineage with no copy in flight is simply not a key; the caller falls back
    to ``absent``, which is the declared answer rather than a guess."""
    hits: dict[int, int] = {}
    total: dict[int, int] = {}
    for live in fam.live:
        s = live.node.species
        hits[s] = hits.get(s, 0) + int(np.isin(live.states, letters).sum())
        total[s] = total.get(s, 0) + live.states.size
    return {s: hits[s] / total[s] for s in total if total[s]}


def grow(genes, trees, models, lengths, rates, offers, step, rng, *, starts=None,
         progress=False):
    """Walk every family forward together, slice by slice.

    ``genes`` is the list of gene specs in the order they were written; the others are dicts keyed by
    the gene's **name**: its `GeneTree`, its model, its site count, its ``(base, factors)`` rate — a
    number and the list of ``(driver name, mapping)`` it reads — and the ``(letters, absent)`` it
    offers, or ``None`` where it offers nothing.

    Returns ``{name: (states_by_id, founding_states, length_by_id)}``: the arrays every node of that
    family ended up with, the sequence it began with, and each branch's length in substitutions per
    site — accumulated slice by slice, because a driven rate is not constant down a branch and a
    phylogram scaled by one sample of it would not be the tree its own alignment was drawn along.
    """
    fams = {g.name: _Family(g, trees[g.name]) for g in genes}
    letters = {g.name: (tuple(models[g.name].alphabet.index(c) for c in offers[g.name][0]),
                        offers[g.name][1]) if offers[g.name] else None
               for g in genes}
    cache: dict[tuple[int, float], np.ndarray] = {}

    # every family founds from a set of stationary frequencies, in the order written, so the run is
    # reproducible from its seed and the order is the one the reader can see. `starts` gives a gene a
    # DIFFERENT model's frequencies to found from: it then arrives with a foreign composition and
    # ameliorates toward its own, which is the one way a composition has anywhere to go.
    for g in genes:
        m = models[g.name]
        founder = m if not starts or starts.get(g.name) is None else starts[g.name]
        fams[g.name].founding = rng.choice(
            m.k, size=lengths[g.name], p=founder.stationary).astype(np.int8)

    # A family exists from its origination, so nothing of it is in flight before then. The walk
    # therefore starts at the earliest origination and admits each family when its own comes round.
    pending = {g.name: fams[g.name].tree for g in genes}
    t = min(gt.origination for gt in pending.values())
    end = max(_deepest(gt.complete) for gt in pending.values())
    slices = 0
    while t < end - 1e-12:
        nxt = min(t + step, end)
        for name, gt in list(pending.items()):          # a family joins the walk at its origination
            if gt.origination < nxt:
                fam = fams[name]
                fam.live.append(_Live(gt.complete, gt.origination, fam.founding))
                del pending[name]
        # the drivers, frozen for the whole slice: this is the approximation, and the only one
        shares = {name: _composition_by_lineage(fams[name], *letters[name])
                  for name in fams if letters[name]}
        for g in genes:                                 # gene order, so the draws are reproducible
            _advance(fams[g.name], models[g.name], rates[g.name], shares, offers, nxt, rng, cache)
        t = nxt
        slices += 1
    for fam in fams.values():                           # whatever is still in flight ends here
        for live in fam.live:
            fam.states[id(live.node)] = live.states
            fam.length_of[id(live.node)] = live.accrued
    return {name: (fam.states, fam.founding, fam.length_of) for name, fam in fams.items()}


def _deepest(node: GeneNode) -> float:
    """When the last thing in this tree happens — the walk's horizon."""
    best = node.time
    stack = list(node.children)
    while stack:
        n = stack.pop()
        best = max(best, n.time)
        stack.extend(n.children)
    return best


def _advance(fam: _Family, model, rate, shares, offers, until: float, rng, cache) -> None:
    """Carry one family's live copies forward to ``until``, splitting at every node they pass.

    A copy that reaches the end of its branch inside this slice does not wait: its node is recorded,
    its children start from its states, and they carry on within the same slice. The *rate* they
    carry on at is still the slice's, which is what makes the slice the unit of approximation rather
    than the branch."""
    base, factors = rate
    out: list[_Live] = []
    todo = fam.live
    while todo:
        live = todo.pop()
        stop = min(live.node.time, until)
        dt = stop - live.at
        if dt > 0.0:
            f = 1.0
            for driver, mapping in factors:
                # the driver's share on THIS lineage, or the value it declared for a lineage that
                # carries none of it
                per_lineage = shares[driver]
                f *= mapping.multiplier(per_lineage.get(live.node.species, offers[driver][1]))
            bl = base * f * dt
            if bl > 0.0:
                live.states = _sample(live.states, _cdf_for(cache, model, bl), rng)
            live.accrued += bl
            live.at = stop
        if live.at >= live.node.time - 1e-12:           # the branch ended: record it and split
            fam.states[id(live.node)] = live.states
            fam.length_of[id(live.node)] = live.accrued
            for child in live.node.children:
                todo.append(_Live(child, live.node.time, live.states))
        else:
            out.append(live)
    fam.live = out


def scaled_tree(gt: GeneTree, length_of: dict[int, float]) -> GeneTree:
    """The family's tree with node ``time`` holding cumulative **substitutions/site** — a phylogram,
    built from the lengths the walk actually accrued rather than from one sample of a rate that moved
    down every branch."""
    root = gt.complete
    scaled_root = GeneNode(root.kind, root.species, length_of[id(root)], root.copy)
    stack = [(root, scaled_root)]
    while stack:
        onode, snode = stack.pop()
        for ochild in onode.children:
            schild = GeneNode(ochild.kind, ochild.species,
                              snode.time + length_of[id(ochild)], ochild.copy)
            snode.children.append(schild)
            stack.append((ochild, schild))
    return GeneTree(gt.family, scaled_root, 0.0)


def check_step(step, tallest: float) -> float:
    """``step`` is the stretch of species time a composition is held fixed across."""
    if isinstance(step, bool) or not isinstance(step, (int, float)) or not math.isfinite(step) \
            or step <= 0:
        raise ValueError(
            f"step is the stretch of time each gene's composition is held fixed across, in the "
            f"tree's own units, so it must be finite and positive; got {step!r}.")
    if step >= tallest:
        raise ValueError(
            f"step={step:g} is not shorter than the tree itself ({tallest:.3g}), so every gene "
            f"would read the other's composition once, at the start, and the loop would be one "
            f"conditioned run in each direction. Pick a step the composition moves little within.")
    return float(step)
