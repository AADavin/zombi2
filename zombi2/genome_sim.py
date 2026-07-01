"""Algorithm 2 — forward gene families along a fixed species tree.

One global continuous-time (Gillespie) process runs over all currently-alive branches,
interleaved with species-tree node times. This is the *single* simulator loop; it talks
only to the :class:`~zombi2.genome.Genome`, :class:`~zombi2.rates.RateModel` and
:class:`~zombi2.events.EventSampler` interfaces, so it never changes when a new genome
representation, rate model or event type is added.

Speciation is *implicit*: at a species-tree node the parent branch's genome is cloned
into each child branch. No speciation event is written to the log — a gene lineage's
splits are recovered later from the species tree itself (this is what keeps v1 minimal
while leaving gene-tree reconstruction possible in v1.1).
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from ._sampling import EventSampler, NumpyEventSampler
from .events import EventLog, EventRecord, EventType, GeneOp
from .genome import Genome, IdManager, UnorderedGenome
from .rates import RateModel
from .tree import Tree, TreeNode


@dataclass
class GenomeResult:
    """Output of a gene-family simulation."""

    event_log: EventLog
    leaf_genomes: dict[TreeNode, Genome]  # extant leaf -> its final genome
    ids: IdManager


class GenomeSimulator:
    """Forward D/T/L/O simulation along a fixed species tree."""

    def __init__(self, sampler: EventSampler | None = None, *, max_events_per_interval: int = 1_000_000):
        self.sampler = sampler or NumpyEventSampler()
        self.max_events_per_interval = max_events_per_interval

    def simulate(
        self,
        tree: Tree,
        rate_model: RateModel,
        rng: np.random.Generator,
        *,
        initial_size: int = 20,
        genome_factory=UnorderedGenome,
    ) -> GenomeResult:
        """Simulate gene families on ``tree``.

        ``genome_factory(ids)`` builds the root genome; every other genome is produced
        by :meth:`Genome.clone`, so a different representation is a one-argument swap
        with no change to this loop.
        """
        ids = IdManager()
        log = EventLog()
        root = tree.root
        rate_model.bind_rng(rng)  # lets stateful rate models seed / reset per run

        # --- seed the root genome ------------------------------------------
        root_genome = genome_factory(ids)
        for _ in range(initial_size):
            params = rate_model.target_params(EventType.ORIGINATION, root_genome, root.name, root.time)
            ops = root_genome.originate(rng, params)
            log.add(EventRecord(EventType.ORIGINATION, root.name, root.time, ops))

        # --- root speciation: seed the child branches ----------------------
        alive: dict[TreeNode, Genome] = {}
        self._speciate(root_genome, root, alive, log)

        leaf_genomes: dict[TreeNode, Genome] = {}

        # --- walk species-tree node events in time order -------------------
        node_events = sorted(
            (n for n in tree.nodes_preorder() if n.parent is not None),
            key=lambda n: n.time,
        )
        t = root.time
        for node in node_events:
            self._evolve_interval(alive, t, node.time, rate_model, log, rng)
            t = node.time
            genome = alive.pop(node)
            if node.is_leaf():
                if node.is_extant:
                    leaf_genomes[node] = genome
            else:  # speciation: re-mint lineage ids into each child and log it
                self._speciate(genome, node, alive, log)

        return GenomeResult(event_log=log, leaf_genomes=leaf_genomes, ids=ids)

    # --- speciation: re-mint lineage ids into both children, log it --------
    @staticmethod
    def _speciate(genome, node, alive, log):
        child1, child2 = node.children
        g1, map1 = genome.clone_reminting()
        g2, map2 = genome.clone_reminting()
        for (old, new1, fam), (_o2, new2, _f2) in zip(map1, map2):
            log.add(EventRecord(
                EventType.SPECIATION, node.name, node.time,
                [GeneOp(old, fam, "parent"), GeneOp(new1, fam, "child"),
                 GeneOp(new2, fam, "child")],
            ))
        alive[child1] = g1
        alive[child2] = g2

    # --- Gillespie over a constant-membership interval ---------------------
    def _evolve_interval(
        self,
        alive: dict[TreeNode, Genome],
        t0: float,
        t1: float,
        rate_model: RateModel,
        log: EventLog,
        rng: np.random.Generator,
    ) -> None:
        t = t0
        branches = list(alive.keys())  # membership is constant across (t0, t1)
        for _ in range(self.max_events_per_interval):
            entries = []  # (branch, EventWeight)
            weights: list[float] = []
            total = 0.0
            for b in branches:
                genome = alive[b]
                supported = genome.supported_events()
                for ew in rate_model.event_weights(genome, b.name, t):
                    if ew.rate > 0.0 and ew.event in supported:
                        entries.append((b, ew))
                        weights.append(ew.rate)
                        total += ew.rate
            if total <= 0.0:
                return
            dt = self.sampler.next_waiting_time(total, rng)
            if not math.isfinite(dt) or t + dt >= t1:
                return
            t += dt
            b, ew = entries[self.sampler.choose_index(weights, rng)]
            self._fire(ew, b, alive, t, rate_model, log, rng)
        raise RuntimeError(
            f"exceeded max_events_per_interval={self.max_events_per_interval}; "
            "gene families are likely growing without bound (duplication > loss). "
            "Set carrying_capacity= or max_copies= on the rate model to regulate growth."
        )

    # --- apply a single event ---------------------------------------------
    def _fire(
        self,
        ew,  # EventWeight(event, family, rate)
        branch: TreeNode,
        alive: dict[TreeNode, Genome],
        t: float,
        rate_model: RateModel,
        log: EventLog,
        rng: np.random.Generator,
    ) -> None:
        genome = alive[branch]
        event, family = ew.event, ew.family
        params = rate_model.target_params(event, genome, branch.name, t)

        if event is EventType.ORIGINATION:
            ops = genome.originate(rng, params)
            log.add(EventRecord(EventType.ORIGINATION, branch.name, t, ops))
            return

        if event is EventType.TRANSFER:
            recipients = [x for x in alive if x is not branch]
            if not recipients:  # no co-existing lineage (should not happen for N>=2)
                return
            selection = genome.draw_target(EventType.TRANSFER, rng, params, family=family)
            segment = genome.extract_segment(selection, rng)  # re-mints donor + copy
            recipient = recipients[int(rng.integers(len(recipients)))]
            at = alive[recipient].choose_insertion_point(segment, rng)
            alive[recipient].insert_segment(segment, at, rng)
            fam = segment.family
            for old, cont, g in zip(segment.donor_old_gids, segment.donor_cont_gids, segment.genes):
                log.add(EventRecord(
                    EventType.TRANSFER, branch.name, t,
                    [GeneOp(old, fam, "parent"), GeneOp(cont, fam, "donor_copy"),
                     GeneOp(g.gid, fam, "transfer_copy")],
                    donor=branch.name, recipient=recipient.name,
                ))
            return

        # duplication or loss
        selection = genome.draw_target(event, rng, params, family=family)
        ops = genome.apply(event, selection, rng, params)
        log.add(EventRecord(event, branch.name, t, ops))
