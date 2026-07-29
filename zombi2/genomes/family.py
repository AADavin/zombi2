"""Genomes I — the D/T/L/O gene-family core.

A genome is a multiset of gene families that evolves along the species tree by four events:
**origination** (a new family arises in a lineage — per lineage), **duplication** (a gene copy
duplicates, +1 copy in its family — per copy), **loss** (a gene copy is lost — per copy), and
**transfer** (a copy is donated to a *contemporaneous* lineage — per copy). Rates are the
cross-level ``scope(base) × modifiers`` grammar (``SPEC §5``); the defaults are the natural
"per what?" for each event.

This reads as the genome twin of `zombi2.species`: one forward Gillespie over the
**complete** tree, plain/frozen dataclasses, an event log as the source of truth (per-family gene
trees are derived from it later), ``as_rate``/``.effective`` for every rate. Because a transfer at
time ``t`` couples two lineages alive at ``t``, the engine evolves **all lineages alive at once**
along one global clock — exactly like ``species_tree._grow`` over its ``alive`` list, except the
species tree is a fixed input (its ``end_time``s form a schedule that decides who is alive), so
there is no birth-death race, no survival conditioning. Speciations and extinctions from that
schedule enter/retire lineages; between them one Gillespie fires D/T/L/O. ``transfer=0`` is the
special case where the lineages are independent — same law as evolving each segment alone.

"""

from __future__ import annotations

import collections
import math
import pathlib
from dataclasses import dataclass, field
from functools import cached_property
from typing import TYPE_CHECKING

import numpy as np

from ..rates.mapping import Between, check_kernel_fires
from ..rates.modifiers import ByFamily, DrivenBy, OnTime
from ..rates.rate import Rate, as_rate
from ..rates.scope import PerCopy, PerLineage, Scope
from ..tree import Tree, as_tree
from ._live import enter, retire, weighted_index, without_cyclic_gc
from ._transfer import Clades, Distance, mean_root_to_tip, recipient_index, resolve_groups

from .._runtime.outputs import grouped_dir
from .._runtime.progress import progress_bar
from .events import Event, events_tsv, gene_label
from .gene_trees import GeneTree, gene_trees_from_events, write_gene_trees
from .profiles import Profiles, profiles_from_genomes

if TYPE_CHECKING:  # a streamed run returns a StreamedRun (built by the per-family engine); type-only
    from ._perfamily import StreamedRun

#: The rate grammar this level wires (SPEC §5) — read by the engine gates below and by the CLI's
#: help, so a modifier is never advertised without being implemented. Each rate keeps its natural
#: scope this slice, ``DrivenBy`` is wired for the single-lineage events; the ordered engine wires
#: ``OnTime`` and ``ByFamily``, the nucleotide one ``OnTime`` and ``DrivenBy``. The gates say so per rate.
WIRED_MODIFIERS = (OnTime, DrivenBy, ByFamily)


@dataclass(frozen=True)
class GeneCopy:
    """One gene copy: a member of family ``family``, identified by a globally-unique ``id``. Its
    birth/death times and parentage live in the event log (the source of truth); the copy carries
    only what a genome snapshot needs to be self-describing — who it is and which family it is in. A
    genome may hold several copies sharing a ``family`` (that family's copy count)."""

    id: int
    family: int


@dataclass
class FamilyGenomesResult:
    """What ``simulate_genomes_family`` returns: the ``complete_tree`` it ran on, the final
    ``genomes`` at **every** node (extant and extinct), the ``events`` log (the compact source of
    truth), and the ``seed``. The observed genomes are the extant tips —
    ``{n.id: genomes[n.id] for n in complete_tree.extant()}``. The phyletic ``profiles`` are derived
    from those tips on access, and ``write`` materialises the chosen outputs to disk."""

    complete_tree: Tree
    genomes: dict[int, tuple[GeneCopy, ...]]
    events: list[Event]
    seed: int | None
    #: ``{name: family id}`` for families declared by ``family_names=[…]`` — the handle to a *named* family
    #: (a toxin, an operon) that you can look up in the genome; empty when only anonymous families were used.
    family_names: dict[str, int] = field(default_factory=dict)
    #: The genome the run **started** with, at the root lineage's origination — before any event.
    #: It is not in `genomes`, which holds a genome per *node*, and a node sits at the **end**
    #: of its branch: the root branch is real simulated time, so ``genomes[root]`` is this genome plus
    #: whatever happened along the stem. The same reason ``GeneTree.origination`` is its own field.
    initial_genome: tuple[GeneCopy, ...] = ()

    def __repr__(self) -> str:
        return (f"FamilyGenomesResult({len(self.complete_tree.extant())} extant genomes, "
                f"{len(self.genomes)} nodes, {len(self.events)} events, seed={self.seed})")

    def family_counts(self, node_id: int) -> collections.Counter:
        """A multiset view of one node's genome: ``family id → copy count``."""
        return collections.Counter(c.family for c in self.genomes[node_id])

    def has_family(self, node_id: int, name: str) -> bool:
        """Whether the named family ``name`` (declared via ``family_names=``) is present — has ≥ 1 copy — in
        the genome at ``node_id``. The presence signal a joint ``DrivenBy("genomes:<name>", …)`` reads."""
        if name not in self.family_names:
            raise KeyError(f"no named family {name!r}; declared families are {sorted(self.family_names)}")
        fid = self.family_names[name]
        return any(c.family == fid for c in self.genomes[node_id])

    @cached_property
    def profiles(self) -> Profiles:
        """The phyletic profiles — each gene family's copy count in each extant species — derived
        from the observed genomes (the classic comparative-genomics matrix). See `profiles`."""
        extant = [n.id for n in self.complete_tree.extant()]
        return profiles_from_genomes(self.genomes, extant)

    @cached_property
    def gene_trees(self) -> dict[int, GeneTree]:
        """``{family id: GeneTree}`` — each family's true genealogy inside the complete tree,
        derived from the event log. Each ``GeneTree`` exposes ``.complete`` and ``.extant``. See
        `gene_trees`."""
        return gene_trees_from_events(self.events, self.complete_tree)

    def write(self, directory, outputs=("events", "profiles", "genomes", "initial_genome",
                                        "gene_trees", "species_tree"), *, flat: bool = False) -> None:
        """Materialise chosen ``outputs`` to ``directory`` (created if needed):

        - ``"events"`` → ``genome_events.tsv``, the event log (the source of truth).
        - ``"profiles"`` → ``profiles.tsv``, the family × extant-species copy-count matrix.
        - ``"genomes"`` → ``genomes.tsv``, every node's gene content, one row per gene copy —
          **ancestors included**, where ``profiles.tsv`` counts only the extant tips.
        - ``"initial_genome"`` → ``initial_genome.tsv``, the genome the run started with. Its own
          file, not a row in ``genomes.tsv``, because it belongs to no node: it sits at the start of
          the root branch, and every ``lineage`` in that table is a node at the end of one.
        - ``"gene_trees"`` → ``gene_tree_fam<family>_{complete,extant}.nwk`` under ``gene_trees/``,
          each family's true genealogy. A family with no surviving copy writes no ``_extant`` file.

        - ``"species_tree"`` → ``species_complete.nwk``, the tree the run evolved along. Written
          because a directory of gene trees with no species tree is not a dataset anyone can use:
          every one of these outputs is *indexed by* that tree's node labels, and the truth a gene
          tree is compared against is the species tree it grew inside. A run written from Python
          used to leave it out entirely, so the quickstart handed back gene trees with nothing to
          compare them to and said nothing about it.
        The gene trees are two files per family, so they get a subdirectory rather than burying the
        tables above; ``flat=True`` writes everything into ``directory`` instead.
        """
        d = pathlib.Path(directory)
        d.mkdir(parents=True, exist_ok=True)
        names = self.complete_tree.labels()      # e<id> for a lineage that died; n<id> for the rest
        if "events" in outputs:
            (d / "genome_events.tsv").write_text(events_tsv(self.events, names), encoding="utf-8")
        if "profiles" in outputs:
            (d / "profiles.tsv").write_text(self.profiles.to_tsv(), encoding="utf-8")
        if "genomes" in outputs:
            (d / "genomes.tsv").write_text(self._genomes_tsv(), encoding="utf-8")
        if "initial_genome" in outputs:
            (d / "initial_genome.tsv").write_text(self._initial_genome_tsv(), encoding="utf-8")
        if "gene_trees" in outputs:
            write_gene_trees(self.gene_trees, grouped_dir(d, "gene_trees", flat), names)
        if "species_tree" in outputs:
            (d / "species_complete.nwk").write_text(self.complete_tree.to_newick() + "\n",
                                                    encoding="utf-8")

    def _genomes_tsv(self) -> str:
        """Every node's gene content, one row per copy, in the order the genome holds them. The
        family counterpart of the ordered resolution's ``gene_order.tsv`` — without a chromosome
        or a position, because at this resolution a genome is a set, not a sequence."""
        cols = ("lineage", "family", "copy")
        names = self.complete_tree.labels()
        rows = [f"{names[s]}\t{c.family}\t{gene_label(c.id)}"
                for s in sorted(self.genomes)
                for c in sorted(self.genomes[s], key=lambda c: (c.family, c.id))]
        return "\n".join(["\t".join(cols), *rows]) + "\n"

    def _initial_genome_tsv(self) -> str:
        """The genome the run started with — ``genomes.tsv``'s columns without ``lineage``, which is
        the whole point: it belongs to the start of the root branch, not to a node."""
        cols = ("family", "copy")
        rows = [f"{c.family}\t{gene_label(c.id)}"
                for c in sorted(self.initial_genome, key=lambda c: (c.family, c.id))]
        return "\n".join(["\t".join(cols), *rows]) + "\n"


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


def _pick_copy_by_family(rng, genome, mult: dict[int, float]) -> int:
    """A copy index within one lineage, drawn in proportion to each copy's family multiplier.

    The within-lineage twin of `_weighted_index()`. Needed whenever families carry different
    rates: the totals are summed with those multipliers, so the copy has to be drawn with them too,
    or the rate would say one thing and the picking another."""
    total = sum(mult[c.family] for c in genome)
    r = float(rng.random()) * total
    acc = 0.0
    for j, c in enumerate(genome):
        acc += mult[c.family]
        if r < acc:
            return j
    return len(genome) - 1                    # float guard: r == total lands on the last copy


def _sum_mult(mult: dict[int, float], genome) -> float:
    """A lineage's per-family multipliers, summed over its live copies."""
    return sum(mult[c.family] for c in genome)


class _FamilyWeights:
    """Each lineage's summed per-family multipliers, held as arrays parallel to ``gen``.

    With per-family multipliers a per-copy rate is the unit rate times the sum of those multipliers
    over a lineage's live copies, so the Gillespie loop needs that sum for **every** living lineage
    on **every** event. Recomputing all of them each time costs the whole live gene pool per event,
    which is quadratic in genome size — and it is nearly all waste, because one event changes one
    lineage by one copy and leaves the rest untouched. So the sums are kept here across events and
    only the lineage an event actually touched is rebuilt. Rates that share a multiplier table
    (`simulate_genomes_family()` hands the same dict to each rate carrying no ``ByFamily`` of its
    own) share the array too, and so are summed once between them rather than once each.

    A rebuilt sum is the same expression over the same list in the same order, so it is the same
    float to the last bit: a run is byte-for-byte what recomputing everything gave."""

    def __init__(self, mult: dict[str, dict[int, float]], gen) -> None:
        by_table: dict[int, tuple[dict[int, float], list[str]]] = {}
        for key, m in mult.items():
            by_table.setdefault(id(m), (m, []))[1].append(key)
        #: (multipliers, the rates using them, their shared per-lineage sums)
        self._groups = [(m, keys, [_sum_mult(m, g) for g in gen]) for m, keys in by_table.values()]
        #: rate → its sums; the arrays are only ever mutated in place, so this stays valid
        self._view = {key: arr for _m, keys, arr in self._groups for key in keys}
        self._dirty: set[int] = set()

    def current(self, gen) -> dict[str, list[float]]:
        """The sums, with any lineage marked by `touched()` rebuilt first."""
        for k in self._dirty:
            for m, _keys, arr in self._groups:
                arr[k] = _sum_mult(m, gen[k])
        self._dirty.clear()
        return self._view

    def touched(self, k: int) -> None:
        """Lineage ``k``'s genome changed: rebuild its sums on the next `current()`."""
        self._dirty.add(k)

    # enter/retire arrive only from the schedule, and the event branch always loops back through
    # current() before reaching it, so these two never see a pending mark.
    def entered(self, genome) -> None:
        for m, _keys, arr in self._groups:
            arr.append(_sum_mult(m, genome))

    def retired(self, k: int) -> None:
        for _m, _keys, arr in self._groups:
            arr[k] = arr[-1]            # mirror _live.retire's swap-remove, or the arrays desync
            arr.pop()


def _driven_mods(rate) -> list:
    """The `DrivenBy` modifiers a rate carries, or ``[]`` when it is a
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


def resolve_max_family_size(max_family_size) -> int | None:
    """Validate the per-genome family cap — **a plain count of copies in one genome**, or ``None``
    for no cap.

    There is no scope wrapper, because there is no "per what?" left to answer: the cap is compared
    against one genome's own copy count, so it is per genome by construction. It used to take one,
    and that was a mistake worth recording. ``PerLineage(n)`` multiplied the number by the size of
    the *species tree* — so the shipped default of ``PerLineage(10)`` was 1470 copies in one genome
    on a 147-node tree, and more on a bigger one. A per-genome bound that moves when you add species
    is not a bound the person who set it can predict, and one genome does not grow because the tree
    did.

    Migration: ``Global(n)`` is exactly ``n``; ``PerLineage(n)`` was ``n`` × the node count of the
    complete tree, so match the old number by multiplying it out, or (better) pick the cap you
    actually want.
    """
    if max_family_size is None:
        return None
    if isinstance(max_family_size, Scope):
        was = (f"{type(max_family_size).__name__}(n) was n × the size of the species tree"
               if isinstance(max_family_size, PerLineage) else
               f"{type(max_family_size).__name__}(n) was exactly n")
        raise ValueError(
            f"max_family_size is a plain count of copies in one genome — write "
            f"max_family_size={max_family_size.base:g}, not {max_family_size!r}. A cap on copies in "
            f"one genome has only one 'per what?', so the scope said nothing ({was}). None removes "
            f"the cap.")
    # a float is refused rather than rounded: 10 against 10.0 used to mean two different caps, and
    # that is exactly the ambiguity this parameter is done with
    if isinstance(max_family_size, bool) or not isinstance(max_family_size, int):
        raise ValueError(f"max_family_size must be a whole number of copies (or None for no cap), "
                         f"got {max_family_size!r}")
    if max_family_size < 1:
        raise ValueError(f"max_family_size must be at least 1 copy, got {max_family_size!r}")
    return max_family_size


def _at_cap(genome, family: int, cap: int | None) -> bool:
    """Whether ``family`` already fills its quota in this genome — the condition that zeroes the
    duplication rate, and a transfer's arrival, for that family.

    A plain scan, which is right where a genome is one family's copies (the per-family engine) and
    wrong where it is the whole genome — see `_FamilyCounts`, which answers the same question
    for the global loop without walking anything."""
    if cap is None:
        return False
    n = 0
    for c in genome:
        if c.family == family:
            n += 1
            if n >= cap:
                return True
    return False


class _FamilyCounts:
    """How many copies of each family every living lineage holds, kept beside ``gen``.

    ``max_family_size`` asks this on **every** duplication and every arriving transfer, and it used
    to be answered by scanning the lineage's whole genome. That put an O(genome) step in the inner
    loop of a run whose genome is what grows, and it dominated: 70% of a 4000-family run, and the
    whole of the level's superlinear scaling — with the cap removed the same run was linear.

    A counter per lineage answers it by lookup. The count is exact (integers in, integers out), so
    the cap binds where it bound before and a run is byte-for-byte what the scan gave."""

    def __init__(self, gen) -> None:
        self._counts = [collections.Counter(c.family for c in g) for g in gen]

    def at_cap(self, k: int, family: int, cap: int | None) -> bool:
        return cap is not None and self._counts[k][family] >= cap

    def added(self, k: int, family: int) -> None:
        self._counts[k][family] += 1

    def removed(self, k: int, family: int) -> None:
        counts = self._counts[k]
        counts[family] -= 1
        if not counts[family]:
            del counts[family]                  # or the counter grows a key per family ever held

    def entered(self, genome) -> None:
        self._counts.append(collections.Counter(c.family for c in genome))

    def entered_like(self, counts) -> None:
        """Enter a lineage holding the same families in the same numbers as ``counts`` — a daughter
        at a speciation, whose genome is its parent's re-identified. A copy rather than another walk."""
        self._counts.append(counts.copy())

    def retired(self, k: int):
        """Retire lineage ``k`` and give back what it was holding, which is what its daughters
        inherit (the parent is retired before they enter, so the index is gone by then)."""
        counts = self._counts[k]
        self._counts[k] = self._counts[-1]      # mirror _live.retire's swap-remove
        self._counts.pop()
        return counts


def _do_transfer(rng, tree, alive, gen, counts, kd, jd, t, events, new_copy,
                 transfer_to, replacement, self_transfer, depth, to_traj=None, cap=None,
                 groups=None) -> tuple[int, int | None]:
    """The copy ``jd`` of the donor lineage ``kd`` transfers to a contemporaneous recipient lineage.
    The donor is picked by the caller (uniformly across the copy pool, or weighted by lineage when
    the transfer rate is driven), the recipient by ``transfer_to``. Returns the change in total copy
    count (+1 additive, 0 replacement — the arriving copy displaces a resident) and **which** lineage
    received it, or ``None`` when nothing happened.

    **No eligible recipient ⇒ nothing happens.** Under a driven ``transfer_to`` a candidate mapped to
    weight 0 cannot receive, and at some instants that is every candidate. The event is then dropped
    before anything is minted or logged — which is not an approximation: rejecting an event on a
    condition that depends only on the current state is Poisson thinning, so the kept transfers are
    exactly the process whose transfer rate is zero while no recipient is eligible."""
    donor = alive[kd]
    cand = [k for k in range(len(alive)) if self_transfer or k != kd]
    kr = recipient_index(rng, tree, alive, cand, donor, t, transfer_to, depth, to_traj, groups)
    if kr is None:                                     # every candidate weighs 0 — no-op (see above)
        return 0, None
    src = gen[kd][jd]
    fam = src.family
    recipient = alive[kr]
    rg = gen[kr]
    if counts.at_cap(kr, fam, cap):  # the recipient is full of this family: same thinning as above
        return 0, None
    # the donor gene ends; two fresh copies descend from it (ZOMBI1 re-id): the continuation on the
    # donor branch and the transferred copy on the recipient branch — a horizontal edge in the gene tree.
    # The donor swaps one copy for another of the same family, so only the recipient's genome changes
    # composition — which is what the caller has to mark for the per-family weights.
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
            counts.removed(kr, fam)
            delta = 0
    rg.append(xfer)
    counts.added(kr, fam)
    events.append(Event(t, "transfer", donor, fam, cont.id, parent=src.id, donor=donor))
    events.append(Event(t, "transfer", recipient, fam, xfer.id, parent=src.id, recipient=recipient,
                        donor=donor))
    return delta, kr


@without_cyclic_gc
def simulate_genomes_family(tree, *, duplication=0.0, transfer=0.0, loss=0.0, origination=0.0,
                            transfer_to="uniform", replacement=False, self_transfer=False,
                            initial_families=100, family_names=None, family_speed=None,
                            max_family_size=10, seed=None, parallel=False, stream_to=None,
                            outputs=None, progress=False) -> "FamilyGenomesResult | StreamedRun":
    """Evolve a multiset of gene families along a species tree by duplication, transfer, loss, and
    origination.

    ``tree`` is the **complete** species tree (a `Tree`, or a
    `SpeciesResult` whose ``complete_tree`` is used). Genomes evolve on
    **every** lineage, extant and extinct alike, so the true gene-tree history is complete and a
    transfer can arrive "from the dead"; the observed genomes are the extant tips.

    Rates (each a ``scope(base) × modifiers`` spec): ``duplication``/``transfer``/``loss`` default
    **per copy**, ``origination`` **per lineage**. When a transfer fires it moves a copy from a
    uniformly-chosen donor copy to a recipient lineage alive at that instant, chosen by
    ``transfer_to`` — ``"uniform"`` (any other contemporaneous lineage), ``"distance"`` /
    ``Distance(decay=)`` (closer relatives likelier), ``Clades({...}, Between({...}))`` (weighted by
    the donor's and recipient's **named clade**, so transfer can run *between* two clades — see below),
    or ``mod.DrivenBy(source, mapping)`` (weighted by another level; see below). ``replacement=True``
    overwrites a homologous
    copy in the recipient (additive fallback if it has none); ``self_transfer=True`` lets a lineage
    donate to itself. The root starts with ``initial_families`` families of one copy each, recorded
    as originations at the origin. ``family_names=["toxin", …]`` additionally declares **named** families —
    each gets a normal (integer) family id, but its name is remembered in ``result.family_names`` so
    you can track a specific family (``result.has_family(node, "toxin")``); this is the handle a joint
    ``DrivenBy("genomes:toxin", …)`` reads. Deterministic given ``seed``.

    **Conditioning (a trait drives a rate).** Any of the four rates may be *driven by another level* —
    ``loss = 0.25 * mod.DrivenBy("habitat.tsv", {"aquatic": 3.0, "terrestrial": 1.0})`` scales each
    lineage's loss by the habitat on that branch, read from a driver file grown first
    (``traits.simulate_discrete(...).write(dir, outputs=("driver",))``). A driven rate is then
    *per-lineage*: it is summed over the living lineages (each with its own copy count and driver
    value), the affected lineage is drawn weighted by its rate, and the Gillespie steps at every
    mid-branch switch of the driver (SPEC §2). For ``transfer`` the affected
    lineage is the **donor**, so a driven ``transfer`` says how often a lineage *donates*.

    **Conditioning (a trait drives who receives).** ``transfer_to = mod.DrivenBy(source, mapping)`` is
    the other half, and a different model: the mapping's numbers are per-candidate **weights**, not
    rate multipliers, so they leave the total amount of transfer alone and only redistribute it
    (SPEC §5, the choice slot). Candidate lineage ``k`` gets weight ``mapping(driver value on k now)``
    and receives with probability ``w_k / Σw`` — five candidates at weight 1 and five at weight 2 send
    two thirds of transfers to the weight-2 group. Weight 0 means "cannot receive"; when every
    candidate weighs 0 the transfer does not happen (see `_do_transfer()`). The two couplings are
    independent and may be used together or apart.

    ``parallel`` opts into a **separate** engine that evolves the families concurrently, one per worker
    process — the families are independent (a transfer roams a copy across lineages but never mixes two
    families), so it enumerates every family's origination first and then evolves each on its own. It is
    worker-count invariant (each family draws from its own stream spawned off ``seed``) but gives a
    different-though-equally-valid draw than the serial default for a given seed. ``False`` (default)
    runs the serial loop; ``True`` uses every core; an ``int`` sets the worker count. A **driven** rate
    or ``transfer_to`` (``DrivenBy``) has no per-family engine yet — the run announces this and falls
    back to the serial loop. The gain is real but modest (a merge over the whole event log stays
    serial), so a handful of workers is the sweet spot; unlike the sequences level it does not scale far.
    Because it spawns processes, a calling script must guard its entry with ``if __name__ ==
    "__main__":`` (the ``zombi2`` CLI already does).

    ``stream_to=DIR`` takes the same engine to the many-families regime: each family is written straight
    to disk as it finishes — no whole-run merge, no run held in memory (a run that fills gigabytes in
    memory streams in tens of megabytes) — and a light `StreamedRun` handle comes
    back instead of a ``FamilyGenomesResult``. ``outputs=`` picks which files, exactly as
    `FamilyGenomesResult.write()` takes them (default: all of them). It is the per-family engine, and
    ``outputs`` without ``stream_to`` is an error.
    """
    tree = as_tree(tree, level="genomes")
    dup = as_rate(duplication, default_scope=PerCopy)
    tra = as_rate(transfer, default_scope=PerCopy)
    los = as_rate(loss, default_scope=PerCopy)
    org = as_rate(origination, default_scope=PerLineage)
    # this slice wires only the default scope (D/T/L per copy, origination per lineage). The wired
    # modifiers are OnTime (skyline) and DrivenBy (a conditioned/joint driver). A non-default scope
    # would set the *total* rate one way while the engine still picks the affected copy/lineage the
    # default way — a silent mismatch (e.g. a PerCopy origination is base×0 copies, a no-op) — so
    # reject it. DrivenBy is a per-lineage driver, wired for all four events; on transfer the driven
    # lineage is the DONOR (who receives is the separate transfer_to slot, below).
    for label, rate, want in (("duplication", dup, PerCopy), ("transfer", tra, PerCopy),
                              ("loss", los, PerCopy), ("origination", org, PerLineage)):
        if not isinstance(rate.scope, want):
            raise ValueError(
                f"{label} has a {type(rate.scope).__name__} scope, but the family genome engine "
                f"wires only {want.__name__} for {label} this slice — scope overrides are a later slice."
            )
        for m in rate.modifiers:
            if isinstance(m, ByFamily) and label == "origination":
                raise ValueError(
                    "origination carries ByFamily, but origination is the rate at which families are "
                    "CREATED — when it is read there is no family yet to have drawn a factor for. "
                    "Put ByFamily on duplication, transfer or loss, or use family_speed= for a "
                    "family-wide tempo.")
            if isinstance(m, DrivenBy) and isinstance(m.mapping, Between):
                raise ValueError(
                    f"{label} carries DrivenBy(…, Between(…)); a Between kernel is donor-conditioned — "
                    f"it weights a recipient by the (donor, recipient) group pair — so it belongs only "
                    f"in transfer_to (who RECEIVES), never in a rate (which has no donor to condition "
                    f"on). Drive a rate with a Table (a plain dict).")
            if isinstance(m, (OnTime, DrivenBy, ByFamily)):
                continue
            raise ValueError(
                f"{label} carries {type(m).__name__}, which the family genome engine does not "
                f"support yet — only OnTime (skyline), DrivenBy (a conditioned/joint driver) and "
                f"ByFamily (per-family heterogeneity) are wired. Clade drift is a later slice."
            )
    if any(isinstance(m, ByFamily) for rate in (dup, tra, los) for m in rate.modifiers) and \
            any(isinstance(m, DrivenBy) for rate in (dup, tra, los, org) for m in rate.modifiers):
        raise ValueError(
            "ByFamily and DrivenBy on the same run is a later slice: one weights lineages by a "
            "driver and the other weights copies by their family, and combining them means "
            "weighting by the product. Use one or the other for now.")
    if family_speed is not None and not isinstance(family_speed, ByFamily):
        raise ValueError(
            f"family_speed takes a ByFamily modifier — family_speed=mod.ByFamily(spread=0.5) — "
            f"got {family_speed!r}. It is the family-wide slot: one draw per family scaling every "
            f"rate that family has.")
    if transfer_to == "distance":
        transfer_to = Distance()
    if isinstance(transfer_to, Rate):
        # `1.0 * mod.DrivenBy(...)` — the rate spelling, in a slot that is not a rate. Say so rather
        # than let a Rate fall through to the generic "must be …" message.
        raise ValueError(
            "transfer_to takes the DrivenBy modifier on its own, not a rate — write "
            "transfer_to=mod.DrivenBy(source, {...}) with no base number. In this slot the mapping's "
            "numbers are relative WEIGHTS over the candidate recipients (normalised), not a rate "
            "multiplier: they change who receives, never how much transfer happens."
        )
    if isinstance(transfer_to, (list, tuple)):
        raise ValueError(
            "transfer_to takes one recipient rule, not several — combining Distance (relatedness) "
            "with a DrivenBy weighting is a later slice. Give 'uniform', 'distance' / "
            "Distance(decay=), or mod.DrivenBy(source, {...})."
        )
    if transfer_to != "uniform" and not isinstance(transfer_to, (Distance, DrivenBy, Clades)):
        raise ValueError(
            f"transfer_to must be 'uniform', 'distance' / Distance(decay=), "
            f"Clades({{...}}, Between({{...}})) (weight by named clade), or "
            f"mod.DrivenBy(source, {{...}}) (a recipient weight driven by another level), "
            f"got {transfer_to!r}")
    if isinstance(initial_families, bool) or not isinstance(initial_families, int) or initial_families < 0:
        raise ValueError(f"initial_families must be a non-negative integer, got {initial_families!r}")
    family_names = list(family_names) if family_names is not None else []
    for name in family_names:
        if not isinstance(name, str) or not name.strip():
            raise ValueError(f"family_names must be a list of non-empty family names (strings), got {name!r}")
    if len(set(family_names)) != len(family_names):
        raise ValueError(f"family names must be unique, got {family_names}")

    # A family's copies in one genome are capped. Growth compounds — a duplication rate above the
    # loss rate multiplies without bound — so a run needs a ceiling somewhere. An int is that number
    # of copies; a float (the default, 10.0) is that multiple of the lineages in the complete tree,
    # so the bound travels with the size of the run. Refusing an event on a condition that depends
    # only on the current state is Poisson thinning, so what is kept is exactly the process whose
    # duplication rate is zero for a family already at its quota — a declared ceiling, not a
    # truncated run. ``None`` removes it.
    cap = resolve_max_family_size(max_family_size)

    # A Clades rule paints every lineage with its clade once (membership is a fact about the tree, not
    # a driver, so it is constant along a branch and adds no Gillespie breakpoints). A kernel naming
    # only absent groups weights every candidate at its default — secretly uniform — so refuse it here.
    group_of = resolve_groups(tree, transfer_to.groups) if isinstance(transfer_to, Clades) else None
    if group_of is not None:
        check_kernel_fires(transfer_to.between, set(group_of.values()), source_label="clades")

    # conditioning: a rate carrying DrivenBy reads a driver per lineage. Resolve each driver once into
    # a DriverTrajectory (value + next-switch lookups, keyed by the shared species node id) — from a
    # file (str source) or an in-memory trait result (object source). No driven rate ⇒ this is empty
    # and the loop stays byte-identical to an uncoupled run.
    dup_mods, los_mods = _driven_mods(dup), _driven_mods(los)
    org_mods, tra_mods = _driven_mods(org), _driven_mods(tra)
    by_key = {}  # driver key → its source (deduped, so a driver shared across rates resolves once)
    for m in (*dup_mods, *los_mods, *org_mods, *tra_mods):
        by_key.setdefault(m.key, m.source)
    rate_keys = list(by_key)     # the drivers that move a RATE: they alone set the Gillespie horizon
    to_mod = transfer_to if isinstance(transfer_to, DrivenBy) else None
    if to_mod is not None:       # the transfer_to driver is read only at the instant a transfer fires
        by_key.setdefault(to_mod.key, to_mod.source)
    resolved = {}
    if by_key:
        from ..rates.driver import check_mapping_fires, resolve_driver
        resolved = {key: resolve_driver(src, tree) for key, src in by_key.items()}
        # a mapping whose states never occur in the driver leaves every lineage at the default factor,
        # so the rate is never driven and the run is secretly the uncoupled model — refuse it here,
        # naming the driver, rather than let it pass as a coupled run
        for m in (*dup_mods, *los_mods, *org_mods, *tra_mods, *( (to_mod,) if to_mod else () )):
            label = m.source if isinstance(m.source, str) else f"<{type(m.source).__name__}>"
            check_mapping_fires(m.mapping, resolved[m.key].states(), source_label=label)
    # a driven transfer_to changes no rate — the weights are evaluated when a transfer fires, so the
    # recipient driver deliberately stays OUT of trajs (no per-lineage rate weights, no extra horizon
    # breakpoints for it). Only the rate drivers make the loop per-lineage.
    trajs = {key: resolved[key] for key in rate_keys}
    to_traj = resolved[to_mod.key] if to_mod is not None else None
    if to_mod is not None and isinstance(to_mod.mapping, Between):  # a donor-conditioned trait kernel:
        label = to_mod.source if isinstance(to_mod.source, str) else f"<{type(to_mod.source).__name__}>"
        check_kernel_fires(to_mod.mapping, to_traj.states(), source_label=label)

    # Parallel is a *separate* engine (opt-in): families are independent, so it evolves them one per
    # process (SPEC-style — serial by default). `stream_to` takes the same engine one step further —
    # each family is written straight to disk and a light StreamedRun handle comes back, so a run of a
    # million families never has to fit in memory (`outputs` picks which files, as `.write` does). The
    # drivers and clade painting above are resolved once, before the split, and handed to whichever
    # engine runs — they are shared validation and shared input, not one engine's business. A
    # configuration neither engine covers still returns None there, so this serial reference loop runs
    # unchanged (decision A); a streamed run raises instead, being unable to fall back without pulling
    # the whole thing into memory.
    if outputs is not None and stream_to is None:
        raise ValueError(
            "outputs applies to a streamed run (stream_to=DIR), which writes the files itself; for an "
            "in-memory run choose them when you call result.write(outputs=...).")
    if parallel or stream_to is not None:
        from ._perfamily import run_parallel_family
        result = run_parallel_family(
            tree, dup=dup, tra=tra, los=los, org=org, transfer_to=transfer_to,
            replacement=replacement, self_transfer=self_transfer, initial_families=initial_families,
            family_names=family_names, family_speed=family_speed, cap=cap, seed=seed, parallel=parallel,
            progress=progress, stream_to=stream_to, outputs=outputs,
            trajs=trajs, to_traj=to_traj, group_of=group_of,
            driven={"duplication": bool(dup_mods), "transfer": bool(tra_mods),
                    "loss": bool(los_mods), "origination": bool(org_mods)})
        if result is not None:
            return result

    rng = np.random.default_rng(seed)
    copy_counter = 0
    family_counter = 0

    def new_copy(family: int) -> GeneCopy:
        nonlocal copy_counter
        c = GeneCopy(copy_counter, family)
        copy_counter += 1
        return c

    # Per-family multipliers, drawn once when a family is minted and then fixed for its whole life:
    # family_speed scales every rate that family has (one draw), and a ByFamily on a single rate
    # varies that rate on its own (a separate draw). Placement is what decides whether a family's
    # rates move together. Empty unless one of them is used, and then the engine takes its weighted
    # path; a run without either draws nothing here and is byte-identical to before.
    fam_by = {"duplication": next((m for m in dup.modifiers if isinstance(m, ByFamily)), None),
              "transfer": next((m for m in tra.modifiers if isinstance(m, ByFamily)), None),
              "loss": next((m for m in los.modifiers if isinstance(m, ByFamily)), None)}
    any_family = family_speed is not None or any(fam_by.values())
    # A rate with no ByFamily of its own carries family_speed and nothing else, so every such rate
    # holds the same multiplier for the same family: one table shared between them, which is what
    # lets _FamilyWeights sum them once rather than once per rate.
    speed_only: dict[int, float] = {}
    fam_mult: dict[str, dict[int, float]] = {
        key: ({} if m is not None else speed_only) for key, m in fam_by.items()}

    def new_family() -> int:
        nonlocal family_counter
        f = family_counter
        family_counter += 1
        if any_family:
            speed = family_speed.draw(rng) if family_speed is not None else 1.0
            for key, m in fam_by.items():
                fam_mult[key][f] = speed * (m.draw(rng) if m is not None else 1.0)
        return f

    depth = mean_root_to_tip(tree)  # timescale for Distance weighting (unused by "uniform")
    # a Clades rule paints every lineage with its clade once (membership is a fact about the tree, not
    # a driver, so it is constant along a branch and adds no Gillespie breakpoints). A kernel naming
    # only absent groups weights every candidate at its default — secretly uniform — so refuse it here.
    schedule = sorted((tree.nodes[i].end_time, i) for i in tree.nodes)  # (end_time, node_id)

    root = tree.nodes[tree.root]
    t = root.birth_time
    alive: list[int] = []
    gen: list[list[GeneCopy]] = []
    pos: dict[int, int] = {}
    genomes: dict[int, tuple[GeneCopy, ...]] = {}
    events: list[Event] = []
    enter(alive, gen, pos, root.id, [])
    for _ in range(initial_families):  # lay down the origin's genome as originations at t = root.birth_time
        _originate(gen[0], root, t, events, new_copy, new_family)
    named: dict[str, int] = {}  # a minted id per declared name (so GeneCopy.family stays an int)
    for name in family_names:
        fid = new_family()
        named[name] = fid
        c = new_copy(fid)
        gen[0].append(c)
        events.append(Event(t, "origination", root.id, fid, c.id))
    total_copies = len(gen[0])
    initial_genome = tuple(gen[0])   # the run's starting genome: a snapshot before the stem runs

    any_driven = bool(rate_keys)
    # the per-family weight sums, carried across events rather than rebuilt each time (see the class)
    weights = _FamilyWeights(fam_mult, gen) if any_family else None
    counts = _FamilyCounts(gen)      # the family cap's question, answered without walking a genome

    # the species tree's schedule is the run's spine: one entry per speciation/extinction, so how
    # far through it we are is how far through the tree the genomes have got
    bar = progress_bar(len(schedule), "genomes", unit="branch", enabled=progress)
    si = 0
    while si < len(schedule):
        bar.to(si)
        n = total_copies
        k_alive = len(alive)
        ctx = {"copies": n, "lineages": k_alive, "time": t}
        can_xfer = n > 0 and (k_alive >= 2 or self_transfer)  # a recipient must be able to exist
        # a driven rate is per-lineage: sum its effective rate over the living lineages (each read with
        # its own copy count and its branch's driver value), keeping the weights for the affected-lineage
        # pick — the species_tree._grow shape. An undriven rate stays pooled (one .effective, uniform
        # pick), so a run with no coupling is byte-identical to before. For transfer the affected
        # lineage is the donor, so a driven transfer weights who donates.
        w_dup = w_los = w_org = w_tra = None
        if any_family:
            # A per-copy rate pools over copies, so with per-family multipliers the total is the
            # unit rate times the sum of those multipliers over the live copies — and the copy has
            # to be drawn with the same weights, or the rates would say one thing and the picking
            # another. Summed per lineage, so the existing weighted-lineage pick can be reused.
            fw = weights.current(gen)
            unit = {"duplication": dup.effective(copies=1, lineages=1, time=t),
                    "loss": los.effective(copies=1, lineages=1, time=t),
                    "transfer": tra.effective(copies=1, lineages=1, time=t) if can_xfer else 0.0}
            w_dup = [unit["duplication"] * s for s in fw["duplication"]]
            w_los = [unit["loss"] * s for s in fw["loss"]]
            if can_xfer:
                w_tra = [unit["transfer"] * s for s in fw["transfer"]]
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
            if tra_mods and can_xfer:
                w_tra = [tra.effective(copies=len(gen[k]), lineages=1, time=t, drivers=drivers[k])
                         for k in range(k_alive)]
        r_dup = sum(w_dup) if w_dup is not None else (dup.effective(**ctx) if n else 0.0)
        r_los = sum(w_los) if w_los is not None else (los.effective(**ctx) if n else 0.0)
        r_org = sum(w_org) if w_org is not None else org.effective(**ctx)
        r_tra = sum(w_tra) if w_tra is not None else (tra.effective(**ctx) if can_xfer else 0.0)
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
                    if w_dup is not None:  # weighted lineage, then a copy within it
                        k = weighted_index(rng, w_dup, r_dup)
                        j = (_pick_copy_by_family(rng, gen[k], fam_mult["duplication"])
                             if any_family else int(rng.integers(len(gen[k]))))
                    else:
                        k, j = _pick_copy(rng, gen, n)
                    fam = gen[k][j].family
                    if not counts.at_cap(k, fam, cap):
                        _duplicate(gen[k], j, tree.nodes[alive[k]], t, events, new_copy)
                        counts.added(k, fam)
                        total_copies += 1
                        if weights is not None:
                            weights.touched(k)
                elif r < r_dup + r_los:
                    if w_los is not None:
                        k = weighted_index(rng, w_los, r_los)
                        j = (_pick_copy_by_family(rng, gen[k], fam_mult["loss"])
                             if any_family else int(rng.integers(len(gen[k]))))
                    else:
                        k, j = _pick_copy(rng, gen, n)
                    counts.removed(k, gen[k][j].family)      # before the copy leaves the genome
                    _lose_at(gen[k], j, tree.nodes[alive[k]], t, events)
                    total_copies -= 1
                    if weights is not None:
                        weights.touched(k)
                elif r < r_dup + r_los + r_org:
                    k = (weighted_index(rng, w_org, r_org) if w_org is not None
                         else int(rng.integers(k_alive)))  # origination is per lineage
                    _originate(gen[k], tree.nodes[alive[k]], t, events, new_copy, new_family)
                    counts.added(k, gen[k][-1].family)       # the copy _originate just appended
                    total_copies += 1
                    if weights is not None:
                        weights.touched(k)
                else:
                    if w_tra is not None:  # driven: weighted DONOR lineage, then a uniform copy in it
                        kd = weighted_index(rng, w_tra, r_tra)
                        if not gen[kd]:    # only via _weighted_index's r == total float guard: a
                            # zero-weight lineage has no copies to donate, so take the heaviest instead
                            kd = max(range(k_alive), key=lambda k: w_tra[k])
                        jd = (_pick_copy_by_family(rng, gen[kd], fam_mult["transfer"])
                              if any_family else int(rng.integers(len(gen[kd]))))
                    else:
                        kd, jd = _pick_copy(rng, gen, n)
                    delta, kr = _do_transfer(rng, tree, alive, gen, counts, kd, jd, t, events,
                                             new_copy, transfer_to, replacement, self_transfer,
                                             depth, to_traj, cap, group_of)
                    total_copies += delta
                    if weights is not None and kr is not None:
                        weights.touched(kr)   # only the recipient's composition changed (see there)
                continue

        if horizon == next_species:  # advance to the tree's next event(s); process the whole tie-batch
            t = next_species
            while si < len(schedule) and schedule[si][0] == t:
                i = schedule[si][1]
                k_out = pos[i]
                g = gen[k_out]
                genomes[i] = tuple(g)  # finalise this lineage (extant, extinct, or unsampled)
                total_copies -= len(g)
                retire(alive, gen, pos, k_out)
                inherited = counts.retired(k_out)   # what the daughters below inherit, if any
                if weights is not None:
                    weights.retired(k_out)
                node = tree.nodes[i]
                if node.children is not None:  # a speciation: each gene re-ids into each daughter
                    for c in node.children:
                        child_genome = []
                        for old in g:  # ZOMBI1: the gene ends here and continues under a fresh id
                            nc = new_copy(old.family)
                            child_genome.append(nc)
                            events.append(Event(t, "speciation", c, old.family, nc.id, parent=old.id))
                        enter(alive, gen, pos, c, child_genome)
                        counts.entered_like(inherited)   # a re-id of the parent: same families
                        if weights is not None:
                            weights.entered(child_genome)
                        total_copies += len(child_genome)
                si += 1
        else:
            t = horizon  # a skyline breakpoint: advance and re-evaluate the (now changed) rate

    bar.close()
    return FamilyGenomesResult(tree, genomes, events, seed, named, initial_genome)


# --- process spec: a genome bundled but UNEXECUTED, for a joint model to grow with the tree --------

@dataclass(frozen=True)
class FamilyGenome:
    """A family-genome **process** — its D/T/L/O parameters bundled but not yet run (the genome
    twin of `DiscreteTrait`). ``simulate_genomes_family(tree, ...)`` runs
    this on a *fixed* tree; a **joint** model (``joint.simulate_joint(genome=genomes.family(...))``)
    grows the genome *with* the tree whose speciation its gene content drives. Duplication, loss, and
    origination (each a ``scope(base) × modifiers`` rate, ``OnTime`` allowed) plus ``initial_families``
    and named ``family_names`` (the handle a ``DrivenBy("genomes:<name>", …)`` reads). Transfer is not
    available in a joint run: a growing tree's contemporaneous set is still forming as events fire."""

    duplication: object
    loss: object
    origination: object
    initial_families: int
    family_names: tuple

    def _resolve(self):
        """Coerce and validate the three rates for the joint engine — ``(duplication, loss,
        origination)`` as resolved `Rate`s. The genome is the **driver**
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


def family(*, duplication=0.0, loss=0.0, origination=0.0, initial_families=100,
              family_names=None) -> FamilyGenome:
    """A family-genome **process spec** — `FamilyGenome`, unexecuted — for a joint model
    to grow with the tree its gene content drives (``joint.simulate_joint(genome=genomes.family(
    origination=0.2, loss=0.1, family_names=["toxin"]))``). Duplication / loss / origination and named
    ``family_names``; a joint run takes no transfer."""
    fams = tuple(family_names) if family_names is not None else ()
    if isinstance(initial_families, bool) or not isinstance(initial_families, int) or initial_families < 0:
        raise ValueError(f"initial_families must be a non-negative integer, got {initial_families!r}")
    for name in fams:
        if not isinstance(name, str) or not name.strip():
            raise ValueError(f"family_names must be a list of non-empty family names (strings), got {name!r}")
    if len(set(fams)) != len(fams):
        raise ValueError(f"family names must be unique, got {list(fams)}")
    return FamilyGenome(duplication, loss, origination, initial_families, fams)


__all__ = ["simulate_genomes_family", "FamilyGenomesResult", "GeneCopy", "FamilyGenome", "family"]
