"""Genomes I — the unordered D/T/L/O gene-family core.

A genome is a multiset of gene families that evolves along the species tree by four events:
**origination** (a new family arises in a lineage — per lineage), **duplication** (a gene copy
duplicates, +1 copy in its family — per copy), **loss** (a gene copy is lost — per copy), and
**transfer** (a copy is donated to a *contemporaneous* lineage — per copy). Rates are the
cross-level ``scope(base) × modifiers`` grammar (``SPEC §5``); the defaults are the natural
"per what?" for each event.

This reads as the genome twin of :mod:`zombi2.species_tree`: one forward Gillespie over the
**complete** tree, plain/frozen dataclasses, an event log as the source of truth (per-family gene
trees are derived from it later), ``as_rate``/``.effective`` for every rate. Because a transfer at
time ``t`` couples two lineages alive at ``t``, the engine evolves **all lineages alive at once**
along one global clock — exactly like ``species_tree._grow`` over its ``alive`` list, except the
species tree is a fixed input (its ``end_time``s form a schedule that decides who is alive), so
there is no birth-death race, no survival conditioning. Speciations and extinctions from that
schedule enter/retire lineages; between them one Gillespie fires D/T/L/O. ``transfer=0`` is the
special case where the lineages are independent — same law as evolving each segment alone.

Still to come: per-family heterogeneity (``ByFamily`` + ``Speed``), the sparse profiles and lazy
gene-tree views behind the ``record=`` memory dial, and the Rust core. This lives here for now so
the legacy ``zombi2/genomes`` package is untouched.
"""

from __future__ import annotations

import collections
import math
import pathlib
from dataclasses import dataclass, field
from functools import cached_property

import numpy as np

from ..rates.modifiers import DrivenBy, OnTime
from ..rates.rate import as_rate
from ..rates.scope import PerCopy, PerLineage
from ..species import SpeciesResult, Tree
from ._live import enter, retire
from ._transfer import Distance, mean_root_to_tip, recipient_index

from .events import Event, events_tsv
from .gene_trees import GeneNode, GeneTree, gene_trees_from_events
from .chromosomes import ChromosomeEvent
from .nucleotide import NucleotideGenome, NucleotideGenomesResult, simulate_genomes_nucleotide
from .ordered import (
    Chromosome,
    EventPosition,
    Gene,
    Inversion,
    OrderedGenomesResult,
    Translocation,
    Transposition,
    simulate_genomes_ordered,
)
from .profiles import Profiles, profiles_from_genomes

#: The rate grammar this level wires (SPEC §5) — read by the engine gates below and by the CLI's
#: help, so a modifier is never advertised without being implemented. Each rate keeps its natural
#: scope this slice, ``DrivenBy`` is wired for the single-lineage events, and the ordered engine
#: wires ``OnTime`` only; the gates below say so per rate.
WIRED_MODIFIERS = (OnTime, DrivenBy)


@dataclass(frozen=True)
class GeneCopy:
    """One gene copy: a member of family ``family``, identified by a globally-unique ``id``. Its
    birth/death times and parentage live in the event log (the source of truth); the copy carries
    only what a genome snapshot needs to be self-describing — who it is and which family it is in. A
    genome may hold several copies sharing a ``family`` (that family's copy count)."""

    id: int
    family: int


@dataclass
class GenomesResult:
    """What ``simulate_genomes_unordered`` returns: the ``complete_tree`` it ran on, the final
    ``genomes`` at **every** node (extant and extinct), the ``events`` log (the compact source of
    truth), and the ``seed``. The observed genomes are the extant tips —
    ``{n.id: genomes[n.id] for n in complete_tree.extant()}``. The phyletic ``profiles`` are derived
    from those tips on access, and ``write`` materialises the chosen outputs to disk. (Lazy gene
    trees and the ``record=`` scale dial are later slices.)"""

    complete_tree: Tree
    genomes: dict[int, tuple[GeneCopy, ...]]
    events: list[Event]
    seed: int | None
    #: ``{name: family id}`` for families seeded by ``families=[…]`` — the handle to a *named* family
    #: (a toxin, an operon) that you can look up in the genome; empty when only anonymous families were used.
    family_names: dict[str, int] = field(default_factory=dict)

    def family_counts(self, node_id: int) -> collections.Counter:
        """A multiset view of one node's genome: ``family id → copy count``."""
        return collections.Counter(c.family for c in self.genomes[node_id])

    def has_family(self, node_id: int, name: str) -> bool:
        """Whether the named family ``name`` (seeded via ``families=``) is present — has ≥ 1 copy — in
        the genome at ``node_id``. The presence signal a joint ``DrivenBy("genomes:<name>", …)`` reads."""
        if name not in self.family_names:
            raise KeyError(f"no named family {name!r}; seeded families are {sorted(self.family_names)}")
        fid = self.family_names[name]
        return any(c.family == fid for c in self.genomes[node_id])

    @cached_property
    def profiles(self) -> Profiles:
        """The phyletic profiles — each gene family's copy count in each extant species — derived
        from the observed genomes (the classic comparative-genomics matrix). See :mod:`.profiles`."""
        extant = [n.id for n in self.complete_tree.extant()]
        return profiles_from_genomes(self.genomes, extant)

    @cached_property
    def gene_trees(self) -> dict[int, GeneTree]:
        """``{family id: GeneTree}`` — each family's true genealogy inside the complete tree,
        derived from the event log. Each ``GeneTree`` exposes ``.complete`` and ``.extant``. See
        :mod:`.gene_trees`."""
        return gene_trees_from_events(self.events, self.complete_tree)

    def write(self, directory, outputs=("events", "profiles")) -> None:
        """Materialise chosen ``outputs`` to ``directory`` (created if needed):

        - ``"events"`` → ``genome_events.tsv``, the event log (the source of truth).
        - ``"profiles"`` → ``profiles.tsv``, the family × extant-species copy-count matrix.
        """
        d = pathlib.Path(directory)
        d.mkdir(parents=True, exist_ok=True)
        if "events" in outputs:
            (d / "genome_events.tsv").write_text(events_tsv(self.events))
        if "profiles" in outputs:
            (d / "profiles.tsv").write_text(self.profiles.to_tsv())


# --- the live genomes: parallel arrays under swap-remove, the ``species_tree._grow`` shape --------

def _pick_copy(rng, gen, total_copies) -> tuple[int, int]:
    """A uniform global copy pick → ``(lineage index k, copy index j in gen[k])``. Realises
    per-copy scope across the whole pool: every copy, in any lineage, is equally likely."""
    j = int(rng.integers(total_copies))
    for k, g in enumerate(gen):
        if j < len(g):
            return k, j
        j -= len(g)
    raise AssertionError("total_copies out of sync with the genomes")  # unreachable


def _weighted_index(rng, weights: list[float], total: float) -> int:
    """Pick a lineage index in proportion to ``weights`` (which sum to ``total``) — the per-lineage
    pick a driven rate needs, the twin of ``species_tree._weighted_index``. When a rate is driven by
    a trait, lineages differ in their loss/duplication/origination rate, so the affected lineage is
    drawn weighted by its own effective rate rather than uniformly from the pool."""
    r = rng.random() * total
    acc = 0.0
    for i, w in enumerate(weights):
        acc += w
        if r < acc:
            return i
    return len(weights) - 1  # floating-point guard: r == total lands on the last lineage


def _driven_mods(rate) -> list:
    """The :class:`~zombi2.rates.modifiers.DrivenBy` modifiers a rate carries, or ``[]`` when it is a
    plain number/scope/OnTime. A non-empty list means the rate is *per-lineage*: each lineage's factor
    depends on the driver value on that branch, so the engine evaluates the rate lineage-by-lineage and
    picks the affected lineage weighted (the ``species_tree._grow`` shape). Each modifier's ``key``
    identifies its driver in the threaded ``drivers`` dict; its ``source`` resolves to a trajectory."""
    return [m for m in rate.modifiers if isinstance(m, DrivenBy)]


# --- the D/T/L/O mutators (each records to the event log; ids from the minters) -------------------

def _originate(genome, node, t, events, new_copy, new_family) -> None:
    """A new gene family arises: mint a founding copy in a fresh family and record it."""
    c = new_copy(new_family())
    genome.append(c)
    events.append(Event(t, "origination", node.id, c.family, c.id))


def _duplicate(genome, j, node, t, events, new_copy) -> None:
    """The gene at index ``j`` duplicates. In the ZOMBI1 per-segment model every event re-ids: the
    gene *ends* and **two** fresh copies descend from it, so both carry new ids (and the id in any
    node is that node's own)."""
    old = genome[j]
    cont, dup = new_copy(old.family), new_copy(old.family)
    genome[j] = cont                                   # the continuing lineage (a fresh id)
    genome.append(dup)                                 # the new copy (a fresh id)
    events.append(Event(t, "duplication", node.id, old.family, cont.id, parent=old.id))
    events.append(Event(t, "duplication", node.id, old.family, dup.id, parent=old.id))


def _lose_at(genome, j, node, t, events) -> None:
    """The copy at index ``j`` is lost (swap-remove — the genome is an order-agnostic multiset)."""
    lost = genome[j]
    genome[j] = genome[-1]
    genome.pop()
    events.append(Event(t, "loss", node.id, lost.family, lost.id))


def _do_transfer(rng, tree, alive, gen, total_copies, t, events, new_copy,
                 transfer_to, replacement, self_transfer, depth) -> int:
    """A gene transfers from a donor copy to a contemporaneous recipient lineage. Returns the change
    in total copy count: +1 additive, 0 replacement (the arriving copy displaces a resident)."""
    kd, jd = _pick_copy(rng, gen, total_copies)
    donor = alive[kd]
    src = gen[kd][jd]
    fam = src.family
    cand = [k for k in range(len(alive)) if self_transfer or k != kd]
    kr = recipient_index(rng, tree, alive, cand, donor, t, transfer_to, depth)
    recipient = alive[kr]
    rg = gen[kr]
    # the donor gene ends; two fresh copies descend from it (ZOMBI1 re-id): the continuation on the
    # donor branch and the transferred copy on the recipient branch — a horizontal edge in the gene tree.
    cont, xfer = new_copy(fam), new_copy(fam)
    gen[kd][jd] = cont
    delta = 1
    if replacement:
        residents = [p for p, c in enumerate(rg) if c.family == fam and c.id != cont.id]
        if residents:  # homologous overwrite; empty ⇒ additive fallback (the gene still arrives)
            p = residents[int(rng.integers(len(residents)))]
            victim = rg[p]
            rg[p] = rg[-1]
            rg.pop()
            events.append(Event(t, "loss", recipient, fam, victim.id))
            delta = 0
    rg.append(xfer)
    events.append(Event(t, "transfer", donor, fam, cont.id, parent=src.id))
    events.append(Event(t, "transfer", recipient, fam, xfer.id, parent=src.id, recipient=recipient))
    return delta


def simulate_genomes_unordered(tree, *, duplication=0.0, transfer=0.0, loss=0.0, origination=0.0,
                               transfer_to="uniform", replacement=False, self_transfer=False,
                               initial_families=0, families=None, seed=None) -> GenomesResult:
    """Evolve a multiset of gene families along a species tree by duplication, transfer, loss, and
    origination.

    ``tree`` is the **complete** species tree (a :class:`~zombi2.species_tree.Tree`, or a
    :class:`~zombi2.species_tree.SpeciesResult` whose ``complete_tree`` is used). Genomes evolve on
    **every** lineage, extant and extinct alike, so the true gene-tree history is complete and a
    transfer can arrive "from the dead"; the observed genomes are the extant tips.

    Rates (each a ``scope(base) × modifiers`` spec): ``duplication``/``transfer``/``loss`` default
    **per copy**, ``origination`` **per lineage**. When a transfer fires it moves a copy from a
    uniformly-chosen donor copy to a recipient lineage alive at that instant, chosen by
    ``transfer_to`` — ``"uniform"`` (any other contemporaneous lineage) or ``"distance"`` /
    ``Distance(decay=)`` (closer relatives likelier). ``replacement=True`` overwrites a homologous
    copy in the recipient (additive fallback if it has none); ``self_transfer=True`` lets a lineage
    donate to itself. The root starts with ``initial_families`` families of one copy each, recorded
    as originations at the crown. ``families=["toxin", …]`` additionally seeds **named** families —
    each gets a normal (integer) family id, but its name is remembered in ``result.family_names`` so
    you can track a specific family (``result.has_family(node, "toxin")``); this is the handle a joint
    ``DrivenBy("genomes:toxin", …)`` reads. Deterministic given ``seed``.

    **Conditioning (a trait drives a rate).** ``loss``, ``duplication``, or ``origination`` may be
    *driven by another level* — ``loss = 0.25 * mod.DrivenBy("habitat.tsv", {"aquatic": 3.0,
    "terrestrial": 1.0})`` scales each lineage's loss by the habitat on that branch, read from a
    driver file grown first (``traits.simulate_discrete(...).write(dir, outputs=("driver",))``). A
    driven rate is then *per-lineage*: it is summed over the living lineages (each with its own copy
    count and driver value), the affected lineage is drawn weighted by its rate, and the Gillespie
    steps at every mid-branch switch of the driver (SPEC §2, ``coupling-api.md``). Driven transfer
    (which couples two lineages) is a later slice.
    """
    tree = tree.complete_tree if isinstance(tree, SpeciesResult) else tree
    dup = as_rate(duplication, default_scope=PerCopy)
    tra = as_rate(transfer, default_scope=PerCopy)
    los = as_rate(loss, default_scope=PerCopy)
    org = as_rate(origination, default_scope=PerLineage)
    # this slice wires only the default scope (D/T/L per copy, origination per lineage). The wired
    # modifiers are OnTime (skyline) and DrivenBy (a conditioned/joint driver). A non-default scope
    # would set the *total* rate one way while the engine still picks the affected copy/lineage the
    # default way — a silent mismatch (e.g. a PerCopy origination is base×0 copies, a no-op) — so
    # reject it. DrivenBy is a per-lineage driver; conditioned driving is wired for the single-lineage
    # events (loss/duplication/origination) but not yet for transfer (two lineages, donor and recipient).
    _drivable = {"duplication", "loss", "origination"}
    for label, rate, want in (("duplication", dup, PerCopy), ("transfer", tra, PerCopy),
                              ("loss", los, PerCopy), ("origination", org, PerLineage)):
        if not isinstance(rate.scope, want):
            raise ValueError(
                f"{label} has a {type(rate.scope).__name__} scope, but the unordered genome engine "
                f"wires only {want.__name__} for {label} this slice — scope overrides are a later slice."
            )
        for m in rate.modifiers:
            if isinstance(m, OnTime):
                continue
            if isinstance(m, DrivenBy):
                if label not in _drivable:
                    raise ValueError(
                        f"{label} carries DrivenBy, but conditioned driving of {label} is a later slice "
                        f"— loss, duplication, and origination are wired (transfer couples two lineages)."
                    )
                continue
            raise ValueError(
                f"{label} carries {type(m).__name__}, which the unordered genome engine does not "
                f"support yet — only OnTime (skyline) and DrivenBy (a conditioned/joint driver) are "
                f"wired. Per-family heterogeneity (ByFamily, Speed) and clade drift are later slices."
            )
    if transfer_to == "distance":
        transfer_to = Distance()
    if transfer_to != "uniform" and not isinstance(transfer_to, Distance):
        raise ValueError(f"transfer_to must be 'uniform', 'distance', or Distance(decay=), got {transfer_to!r}")
    if isinstance(initial_families, bool) or not isinstance(initial_families, int) or initial_families < 0:
        raise ValueError(f"initial_families must be a non-negative integer, got {initial_families!r}")
    families = list(families) if families is not None else []
    for name in families:
        if not isinstance(name, str) or not name.strip():
            raise ValueError(f"families must be a list of non-empty family names (strings), got {name!r}")
    if len(set(families)) != len(families):
        raise ValueError(f"family names must be unique, got {families}")

    rng = np.random.default_rng(seed)
    copy_counter = 0
    family_counter = 0

    def new_copy(family: int) -> GeneCopy:
        nonlocal copy_counter
        c = GeneCopy(copy_counter, family)
        copy_counter += 1
        return c

    def new_family() -> int:
        nonlocal family_counter
        f = family_counter
        family_counter += 1
        return f

    depth = mean_root_to_tip(tree)  # timescale for Distance weighting (unused by "uniform")
    schedule = sorted((tree.nodes[i].end_time, i) for i in tree.nodes)  # (end_time, node_id)

    root = tree.nodes[tree.root]
    t = root.birth_time
    alive: list[int] = []
    gen: list[list[GeneCopy]] = []
    pos: dict[int, int] = {}
    genomes: dict[int, tuple[GeneCopy, ...]] = {}
    events: list[Event] = []
    enter(alive, gen, pos, root.id, [])
    for _ in range(initial_families):  # seed the crown as originations at t = root.birth_time
        _originate(gen[0], root, t, events, new_copy, new_family)
    family_names: dict[str, int] = {}  # named families: a minted id per name (so GeneCopy.family stays int)
    for name in families:
        fid = new_family()
        family_names[name] = fid
        c = new_copy(fid)
        gen[0].append(c)
        events.append(Event(t, "origination", root.id, fid, c.id))
    total_copies = len(gen[0])

    # conditioning: a rate carrying DrivenBy reads a driver per lineage. Resolve each driver once into
    # a DriverTrajectory (value + next-switch lookups, keyed by the shared species node id) — from a
    # file (str source) or an in-memory trait result (object source). No driven rate ⇒ this is empty
    # and the loop stays byte-identical to an uncoupled run.
    dup_mods, los_mods, org_mods = _driven_mods(dup), _driven_mods(los), _driven_mods(org)
    by_key = {}  # driver key → its source (deduped, so a driver shared across rates resolves once)
    for m in (*dup_mods, *los_mods, *org_mods):
        by_key.setdefault(m.key, m.source)
    trajs = {}
    if by_key:
        from ..rates.driver import resolve_driver
        trajs = {key: resolve_driver(src) for key, src in by_key.items()}
    any_driven = bool(by_key)

    si = 0
    while si < len(schedule):
        n = total_copies
        k_alive = len(alive)
        ctx = {"copies": n, "lineages": k_alive, "time": t}
        # a driven rate is per-lineage: sum its effective rate over the living lineages (each read with
        # its own copy count and its branch's driver value), keeping the weights for the affected-lineage
        # pick — the species_tree._grow shape. An undriven rate stays pooled (one .effective, uniform
        # pick), so a run with no coupling is byte-identical to before.
        w_dup = w_los = w_org = None
        if any_driven:
            drivers = [{key: trajs[key].value(alive[k], t) for key in trajs} for k in range(k_alive)]
            if dup_mods:
                w_dup = [dup.effective(copies=len(gen[k]), lineages=1, time=t, drivers=drivers[k])
                         for k in range(k_alive)]
            if los_mods:
                w_los = [los.effective(copies=len(gen[k]), lineages=1, time=t, drivers=drivers[k])
                         for k in range(k_alive)]
            if org_mods:
                w_org = [org.effective(copies=len(gen[k]), lineages=1, time=t, drivers=drivers[k])
                         for k in range(k_alive)]
        r_dup = sum(w_dup) if w_dup is not None else (dup.effective(**ctx) if n else 0.0)
        r_los = sum(w_los) if w_los is not None else (los.effective(**ctx) if n else 0.0)
        r_org = sum(w_org) if w_org is not None else org.effective(**ctx)
        can_xfer = n > 0 and (k_alive >= 2 or self_transfer)  # a recipient must be able to exist
        r_tra = tra.effective(**ctx) if can_xfer else 0.0    # per copy (transfer is not driven this slice)
        total = r_dup + r_los + r_org + r_tra

        next_species = schedule[si][0]  # the tree's own next event: who is alive changes only here
        horizon = min(next_species, dup.next_change(t), los.next_change(t),
                      org.next_change(t), tra.next_change(t))
        if any_driven:  # a driven rate also changes when the driver switches mid-branch — step there
            driver_next = min((trajs[key].next_change(alive[k], t) for key in trajs
                               for k in range(k_alive)), default=math.inf)
            horizon = min(horizon, driver_next)

        if total > 0.0:
            t_ev = t + float(rng.exponential(1.0 / total))
            if t_ev < horizon:  # a genome event fires before the alive set or the rate changes
                t = t_ev
                r = float(rng.random()) * total
                if r < r_dup:
                    if w_dup is not None:  # driven: weighted lineage, then a uniform copy within it
                        k = _weighted_index(rng, w_dup, r_dup)
                        j = int(rng.integers(len(gen[k])))
                    else:
                        k, j = _pick_copy(rng, gen, n)
                    _duplicate(gen[k], j, tree.nodes[alive[k]], t, events, new_copy)
                    total_copies += 1
                elif r < r_dup + r_los:
                    if w_los is not None:
                        k = _weighted_index(rng, w_los, r_los)
                        j = int(rng.integers(len(gen[k])))
                    else:
                        k, j = _pick_copy(rng, gen, n)
                    _lose_at(gen[k], j, tree.nodes[alive[k]], t, events)
                    total_copies -= 1
                elif r < r_dup + r_los + r_org:
                    k = (_weighted_index(rng, w_org, r_org) if w_org is not None
                         else int(rng.integers(k_alive)))  # origination is per lineage
                    _originate(gen[k], tree.nodes[alive[k]], t, events, new_copy, new_family)
                    total_copies += 1
                else:
                    total_copies += _do_transfer(rng, tree, alive, gen, n, t, events, new_copy,
                                                 transfer_to, replacement, self_transfer, depth)
                continue

        if horizon == next_species:  # advance to the tree's next event(s); process the whole tie-batch
            t = next_species
            while si < len(schedule) and schedule[si][0] == t:
                i = schedule[si][1]
                g = gen[pos[i]]
                genomes[i] = tuple(g)  # finalise this lineage (extant, extinct, or unsampled)
                total_copies -= len(g)
                retire(alive, gen, pos, pos[i])
                node = tree.nodes[i]
                if node.children is not None:  # a speciation: each gene re-ids into each daughter
                    for c in node.children:
                        child_genome = []
                        for old in g:  # ZOMBI1: the gene ends here and continues under a fresh id
                            nc = new_copy(old.family)
                            child_genome.append(nc)
                            events.append(Event(t, "speciation", c, old.family, nc.id, parent=old.id))
                        enter(alive, gen, pos, c, child_genome)
                        total_copies += len(child_genome)
                si += 1
        else:
            t = horizon  # a skyline breakpoint: advance and re-evaluate the (now changed) rate

    return GenomesResult(tree, genomes, events, seed, family_names)


# --- process spec: a genome bundled but UNEXECUTED, for a joint model to grow with the tree --------

@dataclass(frozen=True)
class UnorderedGenome:
    """An unordered-genome **process** — its D/T/L/O parameters bundled but not yet run (the genome
    twin of :class:`~zombi2.traits.DiscreteTrait`). ``simulate_genomes_unordered(tree, ...)`` runs
    this on a *fixed* tree; a **joint** model (``joint.simulate_joint(genome=genomes.unordered(...))``)
    grows the genome *with* the tree whose speciation its gene content drives. Duplication, loss, and
    origination (each a ``scope(base) × modifiers`` rate, ``OnTime`` allowed) plus ``initial_families``
    and named ``families`` (the handle a ``DrivenBy("genomes:<name>", …)`` reads). Transfer is deferred
    for joint runs (a growing tree's contemporaneous set is still forming as events fire)."""

    duplication: object
    loss: object
    origination: object
    initial_families: int
    families: tuple

    def _resolve(self):
        """Coerce and validate the three rates for the joint engine — ``(duplication, loss,
        origination)`` as resolved :class:`~zombi2.rates.rate.Rate`s. The genome is the **driver**
        here, not the target, so its own rates carry no coupling (``OnTime`` is the only modifier)."""
        dup = as_rate(self.duplication, default_scope=PerCopy)
        los = as_rate(self.loss, default_scope=PerCopy)
        org = as_rate(self.origination, default_scope=PerLineage)
        for label, rate, want in (("duplication", dup, PerCopy), ("loss", los, PerCopy),
                                  ("origination", org, PerLineage)):
            if not isinstance(rate.scope, want):
                raise ValueError(
                    f"{label} has a {type(rate.scope).__name__} scope, but a joint genome wires only "
                    f"{want.__name__} for {label} — drop the scope wrapper."
                )
            for m in rate.modifiers:
                if not isinstance(m, OnTime):
                    raise ValueError(
                        f"{label} carries {type(m).__name__}; a joint genome's own rates take only "
                        f"OnTime — the genome is the DRIVER of speciation here, not a driven target."
                    )
        return dup, los, org


def unordered(*, duplication=0.0, loss=0.0, origination=0.0, initial_families=0,
              families=None) -> UnorderedGenome:
    """An unordered-genome **process spec** — :class:`UnorderedGenome`, unexecuted — for a joint model
    to grow with the tree its gene content drives (``joint.simulate_joint(genome=genomes.unordered(
    origination=0.2, loss=0.1, families=["toxin"]))``). Duplication / loss / origination and named
    ``families``; transfer is a later slice for joint runs."""
    fams = tuple(families) if families is not None else ()
    if isinstance(initial_families, bool) or not isinstance(initial_families, int) or initial_families < 0:
        raise ValueError(f"initial_families must be a non-negative integer, got {initial_families!r}")
    for name in fams:
        if not isinstance(name, str) or not name.strip():
            raise ValueError(f"families must be a list of non-empty family names (strings), got {name!r}")
    if len(set(fams)) != len(fams):
        raise ValueError(f"family names must be unique, got {list(fams)}")
    return UnorderedGenome(duplication, loss, origination, initial_families, fams)


__all__ = ["simulate_genomes_unordered", "GenomesResult", "Event", "GeneCopy", "Distance",
           "Profiles", "GeneTree", "GeneNode", "UnorderedGenome", "unordered",
           "simulate_genomes_ordered", "OrderedGenomesResult", "Gene", "Chromosome",
           "ChromosomeEvent", "Inversion", "Transposition", "Translocation", "EventPosition",
           "simulate_genomes_nucleotide", "NucleotideGenomesResult", "NucleotideGenome"]
