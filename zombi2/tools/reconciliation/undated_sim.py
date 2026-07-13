"""Forward simulator for the ALEml_undated / GeneRax ``UndatedDTL`` model — the generative
twin of :mod:`zombi2.tools.reconciliation.undated`.

:mod:`~zombi2.tools.reconciliation.undated` computes ``P(gene_tree | species_tree, rates)`` under
the undated ALE model; this module samples in the other direction — it draws gene-family
histories *from* that exact model. Same rates, both directions, so a family simulated here and
then scored there round-trips, and the pair is a self-contained benchmark generator for
ALE/ALErax and a validation harness for ALElite.

Why this is a **tool** and not another entry to the forward ``genomes`` engine: the undated model
is a fundamentally different process from ZOMBI2's native, dated, contemporaneous-transfer
Gillespie simulator. The rates here are dimensionless **per-branch odds relative to a speciation**
(as in :class:`~zombi2.tools.reconciliation.undated.UndatedDTL`), *not* per-unit-time rates, and
the species tree carries no meaningful dates. Rather than launder those odds into fake times, we
simulate the undated process natively and keep its own semantics.

The model (identical to the likelihood's). A single gene copy sitting at the top of species
branch ``e`` consumes exactly one of four slots, with ``denom = 1 + d + t + l``::

    pD = d/denom   duplication  — two copies re-enter e (they stack: this is why the expected
                                  number of duplications on a branch is not pD)
    pT = t/denom   transfer     — one copy re-enters e, one lands on another branch
                                  (any branch for plain undated; a time-overlapping branch for
                                  reldated, exactly as :func:`undated._transfer_neighbors` picks)
    pL = l/denom   loss         — the copy dies
    pS = 1/denom   speciate/    — at an internal branch, one copy enters each daughter; at a
                   sample         leaf branch, the copy is a sampled survivor

Because D and T re-enter the same branch, a lineage keeps resolving fresh slots on ``e`` until it
is lost or speciates/samples — the geometric within-branch stacking the DP resums analytically.

Output. The sampler emits, per family, the same :class:`~zombi2.genomes.events.EventRecord` log the
native ``genomes`` engine produces, then hands it to
:func:`zombi2.genomes.reconciliation.reconcile`. So the ground truth is a real ZOMBI2
:class:`~zombi2.genomes.reconciliation.Reconciliation` (complete + extant annotated Newicks and an
S/D/T/L event table) that drops straight into ``recon-accuracy``, ``zombi2 sequence`` and the
reconparser — no bespoke format. Fully extinct families (no surviving copy) are counted and feed
the ``n_extinct`` term of :func:`undated.undated_joint_loglik`.

Two exact limits pin the sampler against the DP (see ``tests``): with ``d = t = 0`` a family that
survives in all ``k`` tips of a matched subtree is drawn with probability ``pS**(2k-1)`` (one slot
per speciation and per tip), and the ``((A,B))``-present-only-in-A history with probability
``pS**2 * pL`` — the same closed forms the likelihood module checks.
"""

from __future__ import annotations

import random
from collections import Counter
from dataclasses import dataclass, field

from zombi2.genomes.events import EventRecord, EventType, GeneOp
from zombi2.genomes.reconciliation import Reconciliation, reconcile

from .genetree import GeneTree
from .species import SpeciesTree
from .undated import UndatedDTL, _transfer_neighbors


@dataclass(slots=True)
class UndatedSimResult:
    """The outcome of :func:`simulate_undated`.

    ``reconciliations`` has one entry per family that originated (in simulation order); a fully
    extinct family has ``.extant is None``. ``records`` maps each family id to its raw event log
    (the same shape ``genomes`` writes, so the existing writers/replayers consume it unchanged).
    ``n_extinct`` is the number of families that left no surviving copy — pass it straight to
    :func:`undated.undated_joint_loglik`'s ``n_extinct`` when scoring the survivors.
    """

    reconciliations: list[Reconciliation]
    records: dict[str, list[EventRecord]]
    gid2species: dict[str, str]
    species_tree: SpeciesTree
    model: UndatedDTL
    origination: str
    transfers: str
    total_age: float
    event_counts: dict[str, int] = field(default_factory=dict)
    leaf_names: list[str] = field(default_factory=list)
    family_copies: list[dict[str, int]] = field(default_factory=list)

    @property
    def n_families(self) -> int:
        return len(self.reconciliations)

    @property
    def n_extinct(self) -> int:
        return sum(1 for r in self.reconciliations if r.extant is None)

    @property
    def n_surviving(self) -> int:
        return self.n_families - self.n_extinct

    def gene_trees(self) -> list[GeneTree]:
        """The extant (survivors-only) gene trees, ready for :func:`undated.undated_loglik`.

        Extinct families are skipped (they have no observable gene tree). The tip labels are the
        reconciled ``<species>|<gid>``, so :meth:`GeneTree.from_newick` reads the species with the
        default ``sep="|"``.
        """
        return [GeneTree.from_newick(r.extant, sep="|")
                for r in self.reconciliations if r.extant is not None]

    def profile_rows(self):
        """Yield ``(family_id, counts)`` for each SURVIVING family, where ``counts`` is the
        per-species copy number aligned with :attr:`leaf_names` — the gene-family profile, i.e. the
        phyletic pattern that is the classic undated observable. Fully extinct families are omitted."""
        for i, copies in enumerate(self.family_copies, 1):
            if copies:
                yield str(i), [copies.get(name, 0) for name in self.leaf_names]


def _ensure_dated(tree, transfers: str) -> None:
    """Guarantee the species tree carries usable node times.

    The forward ``genomes`` engine silently produces *nothing* on a length-less tree (every branch
    spans zero time). Undated simulation does not need real dates for plain (``"global"``)
    transfers — branch lengths are cosmetic there — so if the tree has none we lay down unit-length
    branches (node time = depth) rather than fail quietly. Reldated (``"dated"``) transfers *do*
    need the dates to decide time overlap, so a length-less tree is a hard error there.
    """
    nodes = list(tree.nodes_preorder())
    if max((n.time for n in nodes), default=0.0) > 0.0:
        return  # already dated
    if transfers == "dated":
        raise ValueError(
            "transfers='dated' (reldated) needs a dated species tree; the given tree has no "
            "branch lengths. Provide dates, or use transfers='global' (plain undated)."
        )
    # root first (preorder), then children = parent depth + 1
    for node in nodes:
        node.time = 0.0 if node.parent is None else node.parent.time + 1.0


def _interior(rng: random.Random, lo: float, hi: float) -> float:
    """A time strictly inside ``(lo, hi]`` for an intra-branch event (cosmetic — the undated model
    has no real event times; only the ordering parent-before-child matters, which this preserves)."""
    return lo + rng.random() * (hi - lo) if hi > lo else hi


def _simulate_family(family: str, origin_node, sp: SpeciesTree, nb, node_by_index, all_nodes,
                     index_by_name, model: UndatedDTL, rng: random.Random, max_events: int):
    """Draw one family's complete event log under the undated model. Returns
    ``(records, gid2species)`` — the raw :class:`EventRecord` list and the surviving-lineage
    ``{gid: species}`` map, both in exactly the shape :func:`reconcile` expects."""
    pD, pT, pL, pS = model.probs()
    records: list[EventRecord] = []
    gid2species: dict[str, str] = {}
    counter = 0

    def fresh() -> str:
        nonlocal counter
        counter += 1
        return f"{family}.{counter}"

    def branch_top(node) -> float:
        return node.parent.time if node.parent is not None else 0.0

    root_gid = fresh()
    origin_t = branch_top(origin_node)
    records.append(EventRecord(EventType.ORIGINATION, origin_node.name, origin_t,
                               [GeneOp(root_gid, family, "origin")]))

    # An explicit stack (not recursion): duplication stacking can be deep for high odds, and a
    # species tree can be tall. Each item is a lineage entering the top of a branch; popping it
    # resolves exactly one slot, so D/T re-push continuing copies onto the same branch.
    stack: list[tuple[str, object, float]] = [(root_gid, origin_node, origin_t)]
    n_events = 0
    while stack:
        gid, node, t_enter = stack.pop()
        n_events += 1
        if n_events > max_events:
            raise RuntimeError(
                f"family {family!r} exceeded max_events={max_events}: the undated model is "
                f"supercritical at these odds (d+t large relative to l). Lower the duplication/"
                f"transfer odds or raise max_events."
            )
        end = node.time
        u = rng.random()
        if u < pL:                                   # ---- loss ----
            records.append(EventRecord(EventType.LOSS, node.name, _interior(rng, t_enter, end),
                                       [GeneOp(gid, family, "lost")]))
        elif u < pL + pD:                            # ---- duplication (both copies stay on e) ----
            t = _interior(rng, t_enter, end)
            a, b = fresh(), fresh()
            records.append(EventRecord(EventType.DUPLICATION, node.name, t,
                [GeneOp(gid, family, "parent"), GeneOp(a, family, "left"),
                 GeneOp(b, family, "right")]))
            stack.append((a, node, t))
            stack.append((b, node, t))
        elif u < pL + pD + pT:                       # ---- transfer (donor stays, copy leaves) ----
            t = _interior(rng, t_enter, end)
            recipient = _choose_recipient(node, nb, node_by_index, all_nodes, index_by_name, rng)
            cont = fresh()
            genes = [GeneOp(gid, family, "parent"), GeneOp(cont, family, "donor_copy")]
            if recipient is not None:
                tc = fresh()
                genes.append(GeneOp(tc, family, "transfer_copy"))
                records.append(EventRecord(EventType.TRANSFER, node.name, t, genes,
                                           donor=node.name, recipient=recipient.name))
                stack.append((tc, recipient, branch_top(recipient)))
            else:
                # No time-overlapping recipient (e.g. a root stem with no contemporary): the
                # transferred copy has nowhere to land and leaves no descendant. The donor copy
                # continues; the unary transfer node is suppressed when the extant tree is pruned.
                records.append(EventRecord(EventType.TRANSFER, node.name, t, genes,
                                           donor=node.name, recipient=None))
            stack.append((cont, node, t))
        else:                                        # ---- speciation / sample (pS) ----
            if node.children:
                child_gids = [fresh() for _ in node.children]
                records.append(EventRecord(EventType.SPECIATION, node.name, end,
                    [GeneOp(gid, family, "parent"),
                     *(GeneOp(cg, family, "child") for cg in child_gids)]))
                for cg, child in zip(child_gids, node.children):
                    stack.append((cg, child, end))
            else:
                gid2species[gid] = node.name         # sampled survivor in this leaf species

    return records, gid2species


def _choose_recipient(donor, nb, node_by_index, all_nodes, index_by_name, rng: random.Random):
    """Pick a transfer recipient branch (returns a tree node, or ``None``). ``nb is None`` (plain
    undated) → uniform over every other branch; otherwise (reldated) → uniform over the donor's
    time-overlapping branches, or ``None`` if it has none. Matches the recipient set
    :func:`undated._transfer_neighbors` averages over."""
    if nb is None:
        if len(all_nodes) <= 1:
            return None
        while True:
            r = rng.choice(all_nodes)
            if r is not donor:
                return r
    neighbors = nb[index_by_name[donor.name]]
    if not neighbors:
        return None
    return node_by_index[rng.choice(neighbors)]


def simulate_undated(tree, model: UndatedDTL, *, n_families: int = 100,
                     origination: str = "root", transfers: str = "global",
                     seed: int | None = None, max_events: int = 1_000_000) -> UndatedSimResult:
    """Sample ``n_families`` gene families under the undated ALE model on ``tree``.

    ``tree`` is a :class:`zombi2.tree.Tree`. ``model`` carries the per-branch odds
    (:class:`~zombi2.tools.reconciliation.undated.UndatedDTL`). ``origination`` is ``"root"`` (every
    family enters on the root branch — the ZOMBI2 root-seeded convention) or ``"uniform"`` (each
    family enters on a uniformly chosen branch). ``transfers`` is ``"global"`` (plain undated: a
    transfer may land on any branch) or ``"dated"`` (reldated: only on a time-overlapping branch —
    needs a dated ``tree``). ``seed`` makes the draw reproducible.

    Returns an :class:`UndatedSimResult`; each family's ground-truth reconciliation is a real
    :class:`~zombi2.genomes.reconciliation.Reconciliation`.
    """
    if origination not in ("root", "uniform"):
        raise ValueError(f"origination must be 'root' or 'uniform', got {origination!r}")
    if transfers not in ("global", "dated"):
        raise ValueError(f"transfers must be 'global' or 'dated', got {transfers!r}")
    if n_families < 0:
        raise ValueError("n_families must be >= 0")

    _ensure_dated(tree, transfers)
    sp = SpeciesTree.from_tree(tree)
    nb = _transfer_neighbors(sp, transfers)

    rng = random.Random(seed)
    all_nodes = list(tree.nodes_preorder())
    nodes_by_name = {n.name: n for n in all_nodes}
    node_by_index = [nodes_by_name[b.name] for b in sp.branches]
    index_by_name = {b.name: i for i, b in enumerate(sp.branches)}
    total_age = max((n.time for n in all_nodes if not n.children), default=0.0)

    reconciliations: list[Reconciliation] = []
    records_by_family: dict[str, list[EventRecord]] = {}
    all_gid2species: dict[str, str] = {}
    family_copies: list[dict[str, int]] = []
    counts = {"O": 0, "D": 0, "T": 0, "L": 0, "S": 0}

    for k in range(n_families):
        family = str(k + 1)
        origin_node = tree.root if origination == "root" else rng.choice(all_nodes)
        records, gid2species = _simulate_family(
            family, origin_node, sp, nb, node_by_index, all_nodes, index_by_name,
            model, rng, max_events)
        records_by_family[family] = records
        all_gid2species.update(gid2species)
        family_copies.append(dict(Counter(gid2species.values())))
        recon = reconcile(records, gid2species, total_age)
        reconciliations.append(recon)
        counts["O"] += 1
        for ev in recon.events:
            counts[ev.event] = counts.get(ev.event, 0) + 1

    return UndatedSimResult(
        reconciliations=reconciliations,
        records=records_by_family,
        gid2species=all_gid2species,
        species_tree=sp,
        model=model,
        origination=origination,
        transfers=transfers,
        total_age=total_age,
        event_counts=counts,
        leaf_names=sorted(sp.leaf_index),
        family_copies=family_copies,
    )
