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
from dataclasses import dataclass, field

import numpy as np

from zombi2._sampling import EventSampler, Fenwick, NumpyEventSampler
from zombi2.genomes.events import EventLog, EventRecord, EventType, GeneOp, Selection
from zombi2.genomes.genome import Gene, Genome, IdManager, UnorderedGenome
from zombi2.genomes.rates import RateModel
from zombi2.genomes.transfers import TransferModel
from zombi2.tree import Tree, TreeNode


@dataclass
class GenomeResult:
    """Output of a gene-family simulation."""

    event_log: EventLog
    leaf_genomes: dict[TreeNode, Genome]  # extant leaf -> its final genome
    ids: IdManager
    node_genomes: dict = field(default_factory=dict)  # every node -> its genome (retain_internal)


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
        self._conversions = None
        self._cap: int | None = None

    def simulate(
        self,
        tree: Tree,
        rate_model: RateModel,
        rng: np.random.Generator,
        *,
        initial_size: int = 20,
        transfers: TransferModel | None = None,
        conversions=None,
        max_family_size=None,
        genome_factory=UnorderedGenome,
        retain_internal: bool = False,
        log_seed_originations: bool = False,
        seed_originations: dict[str, list[str]] | None = None,
    ) -> GenomeResult:
        """Simulate gene families on ``tree``.

        ``genome_factory(ids)`` builds the root genome; every other genome is produced by
        :meth:`Genome.clone_reminting`. ``transfers`` sets transfer mechanics; ``conversions`` sets
        intra-genome gene-conversion directionality (any object with a ``bias`` attribute — the
        experimental :class:`zombi2.experimental.ConversionModel`; ``None`` = unbiased);
        ``max_family_size`` (int absolute, or float as a multiple of the species count) bounds family
        growth.

        ``log_seed_originations``: when the factory *pre-seeds* families into the root genome
        (e.g. the coupling panel) rather than originating them via ``initial_size``, set this to
        log an ORIGINATION for each at the root — otherwise those families have no birth record
        and gene-tree reconstruction (which starts from the origination) yields nothing. The
        default ``False`` leaves every existing caller (empty factory + ``initial_size`` loop,
        or the nucleotide chromosome factory) byte-identical.

        ``seed_originations`` maps a node *name* to a list of family ids to originate *at that
        node* — one fresh copy each, logged as an ORIGINATION, injected into the node's genome
        the moment the walk reaches it (before it speciates), so the family is absent above that
        node and evolves normally within its subtree. This lets a family be born at an internal
        clade instead of only the root (the coupling panel seeds each family at its inferred
        origin). No RNG is drawn, so the default ``None`` leaves every existing caller
        byte-identical. Requires a multiset genome (one exposing ``_add``).
        """
        self._transfers = transfers or TransferModel()
        self._conversions = conversions
        n_species = len(tree.extant_leaves())
        self._cap = resolve_max_family_size(max_family_size, n_species)

        ids = IdManager()
        log = EventLog()
        root = tree.root
        rate_model.bind(rng, max_family_size=self._cap, tree=tree)

        # --- seed the root genome ------------------------------------------
        root_genome = genome_factory(ids)
        if log_seed_originations:
            # families the factory pre-seeded (e.g. the coupling panel) get a root origination
            # so the genealogy carries their birth record (else gene trees can't be built).
            for gene in list(root_genome.genes()):
                log.add(EventRecord(EventType.ORIGINATION, root.name, root.time,
                                    [GeneOp(gene.gid, gene.family, "origin")]))
        # families a caller asks to be born *at the root* (internal-node births are injected
        # during the walk below, when each node is reached).
        self._seed_at_node(root_genome, root, seed_originations, ids, log)
        for _ in range(initial_size):
            params = rate_model.target_params(EventType.ORIGINATION, root_genome, root.name, root.time)
            ops = root_genome.originate(rng, params)
            log.add(EventRecord(EventType.ORIGINATION, root.name, root.time, ops))

        # --- root speciation: seed the child branches ----------------------
        alive: dict[TreeNode, Genome] = {}
        # persistent per-branch rate cache: node -> (candidate events, subtotal), updated
        # only where a genome changes (events) or membership changes (speciations).
        cache: dict[TreeNode, tuple] = {}
        # a persistent Fenwick tree over every node's subtotal (0 when not alive) gives
        # O(log branches) event-branch selection and O(log) updates — no per-event scan.
        nodes_by_index = list(tree.nodes_preorder())
        index = {node: i for i, node in enumerate(nodes_by_index)}
        name_to_node = {node.name: node for node in nodes_by_index}
        if seed_originations:
            unknown = set(seed_originations) - set(name_to_node)
            if unknown:
                raise ValueError(
                    f"seed_originations references unknown node name(s): {sorted(unknown)}")
        fenwick = Fenwick(len(nodes_by_index))

        def activate(branch):
            cache[branch] = self._branch_weights(alive[branch], branch, rate_model, t)
            fenwick.set(index[branch], cache[branch][1])

        t = root.time
        node_genomes: dict = {}
        if retain_internal:
            node_genomes[root] = root_genome        # the seed genome == the user's input genome
        self._speciate(root_genome, root, alive, log)
        for child in root.children:
            activate(child)

        leaf_genomes: dict[TreeNode, Genome] = {}

        # --- walk species-tree node events in time order -------------------
        node_events = sorted(
            (n for n in tree.nodes_preorder() if n.parent is not None),
            key=lambda n: n.time,
        )
        for node in node_events:
            self._evolve_interval(alive, cache, fenwick, index, nodes_by_index, name_to_node,
                                  t, node.time, rate_model, log, rng)
            t = node.time
            genome = alive.pop(node)
            cache.pop(node, None)
            fenwick.set(index[node], 0.0)  # this branch ends here
            # families born at this node join the genome before it speciates onward
            self._seed_at_node(genome, node, seed_originations, ids, log)
            if retain_internal:
                node_genomes[node] = genome         # the ancestral genome at this node
            if node.is_leaf():
                if node.is_extant:
                    leaf_genomes[node] = genome
            elif len(node.children) == 1:
                # a degree-two species node (e.g. an FBD sampled ancestor): the lineage — and
                # its genome — simply continues, so pass the genome straight to the child
                child = node.children[0]
                alive[child] = genome
                activate(child)
            else:  # speciation: re-mint lineage ids into each child and log it
                self._speciate(genome, node, alive, log)
                for child in node.children:
                    activate(child)

        return GenomeResult(event_log=log, leaf_genomes=leaf_genomes, ids=ids,
                            node_genomes=node_genomes)

    # --- per-node origin seeding: born-at-this-node families join the genome ----
    @staticmethod
    def _seed_at_node(genome, node, seed_originations, ids, log):
        """Originate each family scheduled to be *born at* ``node`` (per ``seed_originations``)
        into ``genome`` as one fresh copy, logging an ORIGINATION so the genealogy carries its
        birth. No RNG is drawn (ids are deterministic), so ``seed_originations=None`` is a no-op
        for every existing caller. Requires a multiset genome exposing ``_add(Gene)``."""
        if not seed_originations:
            return
        families = seed_originations.get(node.name)
        if not families:
            return
        add = getattr(genome, "_add", None)
        if add is None:
            raise TypeError(
                f"seed_originations requires a multiset genome supporting _add(Gene); "
                f"{type(genome).__name__} does not")
        for fam in families:
            gene = Gene(ids.new_gene(), fam)
            add(gene)
            log.add(EventRecord(EventType.ORIGINATION, node.name, node.time,
                                [GeneOp(gene.gid, fam, "origin")]))

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
    def _evolve_interval(self, alive, cache, fenwick, index, nodes_by_index, name_to_node,
                         t0, t1, rate_model, log, rng):
        """Gillespie loop over one inter-speciation interval.

        Membership is constant across ``(t0, t1)``; the persistent Fenwick already holds
        every alive branch's subtotal, so selecting the next event's branch is O(log) and
        needs no per-interval setup. After each event only the changed branch(es) are
        refreshed (cache + Fenwick). A rate model whose weights vary continuously with time
        can set ``time_dependent = True`` to force a full refresh each step.

        A rate model may also expose ``refresh_times(t0, t1)`` — a sorted list of
        ``(time, branch_name)`` at which a branch's weights change on their own (a trait
        drifting/jumping along the branch, say). We interleave those into the loop: the
        integration advances to each breakpoint before it would overshoot it, refreshes just
        that branch, and continues — exact for piecewise-constant schedules, and free when the
        list is empty (every existing rate model).
        """
        refresh_all = getattr(rate_model, "time_dependent", False)
        t = t0
        if refresh_all:
            for b in list(alive):
                cache[b] = self._branch_weights(alive[b], b, rate_model, t)
                fenwick.set(index[b], cache[b][1])

        breaks = rate_model.refresh_times(t0, t1)
        bi, nb = 0, len(breaks)

        for _ in range(self.max_events_per_interval):
            total = fenwick.total
            next_bt = breaks[bi][0] if bi < nb else math.inf
            if total <= 0.0:
                # nothing can fire until a scheduled change alters the rates; skip to it
                if next_bt < t1:
                    t = next_bt
                    bi = self._apply_breaks(breaks, bi, t, name_to_node, alive,
                                            cache, fenwick, index, rate_model)
                    continue
                return
            dt = self.sampler.next_waiting_time(total, rng)
            if not math.isfinite(dt):
                return
            if t + dt >= min(t1, next_bt):
                if next_bt < t1:  # a scheduled refresh falls before the drawn event
                    t = next_bt
                    bi = self._apply_breaks(breaks, bi, t, name_to_node, alive,
                                            cache, fenwick, index, rate_model)
                    continue
                return  # crossed t1 with no earlier breakpoint — interval done
            t += dt

            # (1 - U) keeps the draw in (0, total], so we always land on a live branch
            branch = nodes_by_index[fenwick.find((1.0 - rng.random()) * total)]
            ew = self._pick_entry(cache[branch][0], cache[branch][1], rng)
            changed = self._fire(ew, branch, alive, t, rate_model, log, rng)

            for cb in (list(alive) if refresh_all else changed):
                cache[cb] = self._branch_weights(alive[cb], cb, rate_model, t)
                fenwick.set(index[cb], cache[cb][1])
        raise RuntimeError(
            f"exceeded max_events_per_interval={self.max_events_per_interval}; "
            "gene families are likely growing without bound (duplication/transfer > loss). "
            "Set max_family_size= (or carrying_capacity= on the rate model) to regulate growth."
        )

    def _apply_breaks(self, breaks, bi, t, name_to_node, alive, cache, fenwick, index, rate_model):
        """Refresh every branch scheduled to change at time ``t`` (there may be several), and
        return the advanced breakpoint cursor. Branches no longer alive are skipped."""
        while bi < len(breaks) and breaks[bi][0] == t:
            node = name_to_node.get(breaks[bi][1])
            if node is not None and node in alive:
                cache[node] = self._branch_weights(alive[node], node, rate_model, t)
                fenwick.set(index[node], cache[node][1])
            bi += 1
        return bi

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

    # --- recipient choice (uniform, distance-weighted, and/or per-branch receptivity) -----
    def _choose_recipient(self, donor, alive, t, rng):
        candidates = list(alive) if self._transfers.allow_self else [x for x in alive if x is not donor]
        if not candidates:
            return None
        decay = self._transfers.distance_decay
        recept = self._transfers.receptivity
        if decay is None and recept is None:
            return candidates[int(rng.integers(len(candidates)))]  # uniform fast path (unchanged)

        if decay is None:
            weights = [1.0] * len(candidates)
        else:
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

        if recept is not None:  # per-branch absorption weight (branches not listed default to 1)
            weights = [w * recept.get(c.name, 1.0) for w, c in zip(weights, candidates)]
        if sum(weights) <= 0.0:  # every eligible recipient has zero receptivity — no transfer
            return None
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
            rec_genome = alive[recipient]
            # Field-biased establishment (coupled models only): the copy comes from a real
            # donor, then establishes with a probability set by the recipient's local field —
            # putting part of the coupling on the gain channel while staying donor-limited.
            # Default models return 1.0 and the draw is skipped, so their RNG stream is intact.
            p = rate_model.establishment_probability(selection, rec_genome, t)
            if p < 1.0 and rng.random() >= p:
                return ()  # transfer fired but the copy did not establish — no-op
            segment = genome.extract_segment(selection, rng)  # re-mints donor + copy
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

        if event is EventType.CONVERSION:
            # Intra-genome gene conversion (activated by zombi2.experimental.GeneConversionRates,
            # which emits the CONVERSION weights; dormant otherwise). The donor lineage bifurcates
            # and the recipient copy of the same family is overwritten, so copy number is unchanged
            # and the converted copy coalesces with the donor at time t (concerted evolution). The
            # rate model only fires this on a family with >= 2 copies, so `family` is always set. It
            # is logged as a CONVERSION record (the donor bifurcation; donor == recipient == this
            # branch, an intra-genome event) paired with a LOSS record (the overwritten copy) — a
            # shape the existing gene-tree reconstruction already understands. ``conversions`` is any
            # object with a ``bias`` attribute (or None for unbiased).
            bias = self._conversions.bias if self._conversions is not None else 0.0
            donor_group, loss_group = genome.convert(family, rng, bias)
            log.add(EventRecord(EventType.CONVERSION, branch.name, t, donor_group,
                                donor=branch.name, recipient=branch.name))
            log.add(EventRecord(EventType.LOSS, branch.name, t, loss_group))
            return (branch,)

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

        A genome that performs its own *homologous* replacement (the nucleotide genic model)
        exposes ``pop_replaced_segments``: a list (possibly empty) means it already swapped the
        recipient's syntenic copy, and those removed segments are logged here as recipient losses,
        bypassing the random-removal path below. ``None`` (the default for every other genome, and
        for the plain nucleotide model) leaves the random-removal replacement unchanged.
        """
        pop = getattr(genome, "pop_replaced_segments", None)
        if pop is not None:
            removed = pop()
            if removed is not None:
                for seg in removed:
                    log.add(EventRecord(EventType.LOSS, branch.name, t,
                                        [GeneOp(seg.seg_id, seg.family, "lost")]))
                return

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
