"""Algorithm 2 — forward gene families along a fixed species tree.

One global continuous-time (Gillespie) process runs over all currently-alive branches,
interleaved with species-tree node times. This is the *single* simulator loop; it talks
only to the :class:`~zombi2.genome.Genome`, :class:`~zombi2.rates.RateModel` and
:class:`~zombi2.events.EventSampler` interfaces, so it never changes when a new genome
representation, rate model or event type is added.

Transfer mechanics (recipient choice by phylogenetic distance, replacement vs additive,
self-transfer) live in a :class:`~zombi2.transfers.TransferModel`; a hard family-size cap
(``max_family_size``) bounds growth across duplication and transfer alike.
"""

from __future__ import annotations

import math
from collections import Counter
from dataclasses import dataclass

import numpy as np

from ._sampling import EventSampler, NumpyEventSampler
from .events import EventLog, EventRecord, EventType, GeneOp, Selection
from .genome import Genome, IdManager, UnorderedGenome
from .rates import RateModel
from .transfers import TransferModel
from .tree import Tree, TreeNode


@dataclass
class GenomeResult:
    """Output of a gene-family simulation."""

    event_log: EventLog
    leaf_genomes: dict[TreeNode, Genome]  # extant leaf -> its final genome
    ids: IdManager


def resolve_max_family_size(max_family_size, n_species: int) -> int | None:
    """Resolve the cap: ``int`` -> absolute; ``float`` -> multiple of the species count."""
    if max_family_size is None:
        return None
    if isinstance(max_family_size, bool):  # guard: bool is an int subclass
        raise TypeError("max_family_size must be an int or float, not bool")
    if isinstance(max_family_size, float):
        return max(1, round(max_family_size * n_species))
    return max(1, int(max_family_size))


class GenomeSimulator:
    """Forward gene-family simulation along a fixed species tree."""

    def __init__(self, sampler: EventSampler | None = None, *, max_events_per_interval: int = 1_000_000):
        self.sampler = sampler or NumpyEventSampler()
        self.max_events_per_interval = max_events_per_interval
        self._transfers = TransferModel()
        self._cap: int | None = None

    def simulate(
        self,
        tree: Tree,
        rate_model: RateModel,
        rng: np.random.Generator,
        *,
        initial_size: int = 20,
        transfers: TransferModel | None = None,
        max_family_size=None,
        genome_factory=UnorderedGenome,
    ) -> GenomeResult:
        """Simulate gene families on ``tree``.

        ``genome_factory(ids)`` builds the root genome; every other genome is produced by
        :meth:`Genome.clone_reminting`. ``transfers`` sets transfer mechanics;
        ``max_family_size`` (int absolute, or float as a multiple of the species count)
        bounds family growth.
        """
        self._transfers = transfers or TransferModel()
        n_species = len(tree.extant_leaves())
        self._cap = resolve_max_family_size(max_family_size, n_species)

        ids = IdManager()
        log = EventLog()
        root = tree.root
        rate_model.bind(rng, max_family_size=self._cap, tree=tree)

        # --- seed the root genome ------------------------------------------
        root_genome = genome_factory(ids)
        for _ in range(initial_size):
            params = rate_model.target_params(EventType.ORIGINATION, root_genome, root.name, root.time)
            ops = root_genome.originate(rng, params)
            log.add(EventRecord(EventType.ORIGINATION, root.name, root.time, ops))

        # --- root speciation: seed the child branches ----------------------
        alive: dict[TreeNode, Genome] = {}
        # persistent per-branch rate cache: node -> (candidate events, subtotal). It is
        # updated only where the genome changes (events) or membership changes
        # (speciations), so it is never rebuilt wholesale — the key to O(N) scaling.
        cache: dict[TreeNode, tuple] = {}
        self._speciate(root_genome, root, alive, log)
        for child in root.children:
            cache[child] = self._branch_weights(alive[child], child, rate_model, root.time)

        leaf_genomes: dict[TreeNode, Genome] = {}

        # --- walk species-tree node events in time order -------------------
        node_events = sorted(
            (n for n in tree.nodes_preorder() if n.parent is not None),
            key=lambda n: n.time,
        )
        t = root.time
        for node in node_events:
            self._evolve_interval(alive, cache, t, node.time, rate_model, log, rng)
            t = node.time
            genome = alive.pop(node)
            cache.pop(node, None)
            if node.is_leaf():
                if node.is_extant:
                    leaf_genomes[node] = genome
            else:  # speciation: re-mint lineage ids into each child and log it
                self._speciate(genome, node, alive, log)
                for child in node.children:
                    cache[child] = self._branch_weights(alive[child], child, rate_model, t)

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
    def _evolve_interval(self, alive, cache, t0, t1, rate_model, log, rng):
        """Gillespie loop over one inter-speciation interval, using the persistent cache.

        Branch membership is constant across ``(t0, t1)``, and for the shipped rate models
        a branch's weights depend only on its genome (not on ``t``). We assemble a subtotal
        array from the cache (O(branches) cheap reads — no rate recomputation), then after
        each event recompute only the branch(es) whose genome changed. A rate model that
        varies continuously with time can set ``time_dependent = True`` to force a full
        refresh each step (correctness over speed).
        """
        refresh_all = getattr(rate_model, "time_dependent", False)
        branches = list(alive.keys())
        index = {b: i for i, b in enumerate(branches)}
        t = t0
        if refresh_all:
            for b in branches:
                cache[b] = self._branch_weights(alive[b], b, rate_model, t)
        sub = np.array([cache[b][1] for b in branches], dtype=float)
        total = float(sub.sum())

        for _ in range(self.max_events_per_interval):
            if total <= 0.0:
                return
            dt = self.sampler.next_waiting_time(total, rng)
            if not math.isfinite(dt) or t + dt >= t1:
                return
            t += dt

            bi = int(np.searchsorted(np.cumsum(sub), rng.random() * total, side="right"))
            bi = min(bi, len(branches) - 1)
            branch = branches[bi]
            ew = self._pick_entry(cache[branch][0], float(sub[bi]), rng)
            changed = self._fire(ew, branch, alive, t, rate_model, log, rng)

            to_refresh = branches if refresh_all else changed
            for cb in to_refresh:
                i = index.get(cb)
                if i is None:
                    continue
                cache[cb] = self._branch_weights(alive[cb], cb, rate_model, t)
                total += cache[cb][1] - sub[i]
                sub[i] = cache[cb][1]
        raise RuntimeError(
            f"exceeded max_events_per_interval={self.max_events_per_interval}; "
            "gene families are likely growing without bound (duplication/transfer > loss). "
            "Set max_family_size= (or carrying_capacity= on the rate model) to regulate growth."
        )

    @staticmethod
    def _branch_weights(genome, branch, rate_model, t):
        supported = genome.supported_events()
        kept = [ew for ew in rate_model.event_weights(genome, branch.name, t)
                if ew.rate > 0.0 and ew.event in supported]
        return kept, math.fsum(ew.rate for ew in kept)

    def _pick_entry(self, kept, subtotal, rng):
        r = rng.random() * subtotal
        acc = 0.0
        for ew in kept:
            acc += ew.rate
            if acc >= r:
                return ew
        return kept[-1]

    # --- recipient choice (uniform, or phylogenetic-distance weighted) -----
    def _choose_recipient(self, donor, alive, t, rng):
        candidates = list(alive) if self._transfers.allow_self else [x for x in alive if x is not donor]
        if not candidates:
            return None
        decay = self._transfers.distance_decay
        if decay is None:
            return candidates[int(rng.integers(len(candidates)))]

        # patristic distance at time t between two co-existing lineages is 2*(t - t_MRCA);
        # we only need t_MRCA. Mark the donor's ancestor chain once, then walk each
        # candidate up to the first marked node. O(alive * depth) — the optimisable hotspot.
        ancestors = set()
        node = donor
        while node is not None:
            ancestors.add(node)
            node = node.parent

        def mrca_time(r):
            if r is donor:
                return t  # self-transfer: distance 0
            node = r
            while node not in ancestors:
                node = node.parent
            return node.time

        distances = [2.0 * (t - mrca_time(r)) for r in candidates]
        dmin = min(distances)
        weights = [math.exp(-decay * (d - dmin)) for d in distances]  # softmax-stable
        return candidates[self.sampler.choose_index(weights, rng)]

    # --- apply a single event; return the set of branches whose genome changed ----
    def _fire(self, ew, branch, alive, t, rate_model, log, rng):
        genome = alive[branch]
        event, family = ew.event, ew.family
        params = rate_model.target_params(event, genome, branch.name, t)

        if event is EventType.ORIGINATION:
            ops = genome.originate(rng, params)
            log.add(EventRecord(EventType.ORIGINATION, branch.name, t, ops))
            return (branch,)

        if event is EventType.TRANSFER:
            recipient = self._choose_recipient(branch, alive, t, rng)
            if recipient is None:
                return ()
            selection = genome.draw_target(EventType.TRANSFER, rng, params, family=family)
            segment = genome.extract_segment(selection, rng)  # re-mints donor + copy
            rec_genome = alive[recipient]
            at = rec_genome.choose_insertion_point(segment, rng)
            rec_genome.insert_segment(segment, at, rng)
            for old, cont, g in zip(segment.donor_old_gids, segment.donor_cont_gids, segment.genes):
                log.add(EventRecord(
                    EventType.TRANSFER, branch.name, t,
                    [GeneOp(old, g.family, "parent"), GeneOp(cont, g.family, "donor_copy"),
                     GeneOp(g.gid, g.family, "transfer_copy")],
                    donor=branch.name, recipient=recipient.name,
                ))
            self._reconcile_recipient(rec_genome, segment, recipient, t, params, log, rng)
            return (branch, recipient)

        # duplication / loss / inversion / transposition (one log record per group)
        selection = genome.draw_target(event, rng, params, family=family)
        if (event is EventType.DUPLICATION and self._cap is not None
                and genome.copy_number(selection.genes[0].family) >= self._cap):
            return ()  # family already at the cap — skip this duplication
        for group in genome.apply(event, selection, rng, params):
            log.add(EventRecord(event, branch.name, t, group))
        return (branch,)

    def _reconcile_recipient(self, genome, segment, branch, t, params, log, rng):
        """After a transfer: apply replacement (by chance) and enforce the family cap.

        Works per family in the transferred segment, so it is correct for single-gene and
        multi-gene (ordered) transfers alike.
        """
        protected = {g.gid for g in segment.genes}
        for family, n_inserted in Counter(g.family for g in segment.genes).items():
            total = genome.copy_number(family)
            pre_existing = total - n_inserted
            removals = 0
            if self._cap is not None and total > self._cap:  # forced down to the cap
                removals = total - self._cap
            elif pre_existing >= 1 and self._transfers.replacement > 0:
                removals = sum(rng.random() < self._transfers.replacement for _ in range(n_inserted))
            for _ in range(removals):
                self._remove_one(genome, family, protected, branch, t, params, log, rng)

    def _remove_one(self, genome, family, protected, branch, t, params, log, rng):
        """Remove one copy of ``family``, preferring a pre-existing (non-transferred) one."""
        candidates = [g for g in genome.genes() if g.family == family]
        preferred = [g for g in candidates if g.gid not in protected] or candidates
        if not preferred:
            return
        victim = preferred[int(rng.integers(len(preferred)))]
        for group in genome.apply(EventType.LOSS, Selection(genes=(victim,)), rng, params):
            log.add(EventRecord(EventType.LOSS, branch.name, t, group))
