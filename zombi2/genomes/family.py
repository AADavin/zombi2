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
from typing import ClassVar, TYPE_CHECKING


from ..params.conditioned import check_mapping_fires, names_a_live_level, resolve_driver
from ..params.mapping import check_not_a_kernel
from ..rng import resolve_seed, stream
from ..params.driver import OnTime
from ..params.evaluate import DRAWN, describe, is_implemented, matches_declared, values_at_birth
from ..params.connection import Driven, SetBy
from ..params.parameter import Rate, as_rate
from ..params.retired import check_no_retired_keywords
from ..params.scope import PerCopy, PerLineage
from ..tree import Tree, as_tree, node_from_label
from ._live import enter, retire, weighted_index, without_cyclic_gc
from ._transfer import (mean_root_to_tip, prepare_transfer_to, recipient_index,
                        resolve_transfer_to)

from .._runtime.outputs import fresh_dirs, grouped_dir
from .._runtime.progress import progress_bar
from .._runtime.summary import _stats, write_summary
from .events import Event, GeneEdge, event_counts, events_from_edges, events_tsv, gene_label
from .gene_trees import GeneTree, gene_trees_from_edges, write_gene_trees
from .profiles import Profiles, profiles_from_genomes

if TYPE_CHECKING:  # a streamed run returns a StreamedRun (built by the per-family engine); type-only
    from ._perfamily import StreamedRun

#: The rate grammar this level supports (SPEC §5) — read by the engine gates below and by the CLI's
#: help, so a modifier is never advertised without being implemented. ``duplication`` / ``transfer``
#: / ``loss`` also take a ``PerLineage`` scope override; ``changing_at``, ``scaled_by`` and
#: ``set_by`` are implemented on all four rates (on ``transfer`` the driven lineage is the donor),
#: and a per-family draw on all but ``origination``. The ordered engine declares the same four
#: modifiers as this one, the nucleotide engine only ``changing_at`` and ``scaled_by``. The gates
#: say so per rate.
IMPLEMENTED_MODIFIERS = (OnTime, Driven, SetBy, (DRAWN, "families"))


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
    ``node_genomes`` at **every** node (extant and extinct), the ``events`` log (the compact source of
    truth), and the ``seed``. The observed dataset is the extant tips, ``genomes``. The phyletic
    ``profiles`` are derived from those tips on access, and ``write`` materialises the chosen outputs
    to disk."""

    complete_tree: Tree
    node_genomes: dict[int, tuple[GeneCopy, ...]]
    edges: list[GeneEdge]
    seed: int | None
    #: ``{name: family id}`` for families declared by ``families=[family(…)]`` — the handle to a *named* family
    #: (a toxin, an operon) that you can look up in the genome; empty when only anonymous families were used.
    family_names: dict[str, int] = field(default_factory=dict)
    #: ``{module name: (family name, …)}`` for groups declared by ``modules=`` — a pathway or a
    #: complex, whose *completion* in a lineage (`completion`) is a driver. Empty when none were
    #: declared; a module changes nothing about how the genome evolves.
    modules: dict[str, tuple[str, ...]] = field(default_factory=dict)
    #: The genome the run **started** with, at the root lineage's origination — before any event.
    #: It is not in `node_genomes`, which holds a genome per *node*, and a node sits at the **end**
    #: of its branch: the root branch is real simulated time, so ``node_genomes[root]`` is this genome
    #: plus whatever happened along the stem. The same reason ``GeneTree.origination`` is its own field.
    initial_genome: tuple[GeneCopy, ...] = ()
    #: The per-genome family cap this run actually ran under, resolved (``None`` for no cap). Kept
    #: because the cap is otherwise invisible in the output: when it binds it discards duplications
    #: and arriving transfers, so realised rates fall below the declared ones, and `summary()` is
    #: where a reader finds out that happened to them.
    max_family_size: int | None = None

    def __repr__(self) -> str:
        # "0 nodes" beside a real tip count reads as a broken run; a reopened run that did not write
        # genomes.tsv has the genealogy and not the gene content, so it says which.
        content = f"{len(self.node_genomes)} nodes" if self.node_genomes else "gene content not loaded"
        return (f"FamilyGenomesResult({len(self.complete_tree.extant_leaves())} extant genomes, "
                f"{content}, {len(self.events)} events, seed={self.seed})")

    @property
    def genomes(self) -> dict[str, tuple[GeneCopy, ...]]:
        """The observed dataset — the genome at each **extant** tip, keyed by the **tip name** the
        tree writes: ``n5``.

        Keyed by name because the only thing anyone does with this is join it to the tree, or to
        another level grown on that tree, and both name their tips. ``genomes.tsv`` on disk is
        keyed by name too, so what you get in Python and what you get from a file are the same
        dataset — they used to share no keys at all, and nothing said so. This is the trait level's
        `TraitsResult.values` for gene content.

        `node_genomes` is the run's own record: every node, extant and extinct and internal alike,
        keyed by node id. Use that one to join against ``complete_tree.nodes`` or the event log."""
        extant = list(self.complete_tree.extant_leaves())
        if extant and not self.node_genomes:
            # a run reopened by `read_run` from a directory whose 'genomes' output was not written:
            # the genealogy is all there, the gene content is not. Say that, rather than hand back
            # an empty dict that reads as a run in which nothing survived.
            raise ValueError(
                "this run has no per-node gene content, so there are no genomes to hand back — it "
                "was read back from a directory whose genomes.tsv was not written. Re-run the "
                "genomes level with 'genomes' among its outputs. The gene trees and the event log "
                "are unaffected.")
        name = self.complete_tree.labels()
        return {name[i]: self.node_genomes[i] for i in extant}

    def family_counts(self, node_id: int) -> collections.Counter:
        """A multiset view of one node's genome: ``family id → copy count``."""
        return collections.Counter(c.family for c in self.node_genomes[node_id])

    def completion(self, name: str):
        """A module's completion as a **conditioning driver** — `ModuleCompletion`, a number in
        ``[0, 1]``: the fraction of the module's families a lineage carries.

        Read it with a `Curve`, the way any continuous driver is read; a threshold goes there rather
        than here (``lambda f: 8.0 if f > 0.8 else 1.0``)."""
        from .presence import ModuleCompletion
        if name not in self.modules:
            raise KeyError(f"no module {name!r}; declared modules are {sorted(self.modules)}")
        return ModuleCompletion(self, name)

    def presence(self, name: str):
        """The named family's presence as a **conditioning driver** — `GenePresence`.

        ``has_family`` answers for one node; this answers for every lineage at every instant, which
        is what a driven rate needs::

            switch=PerLineage(0.1).scaled_by(g.presence("tox"), {"present": 5.0, "absent": 1.0})
        """
        from .presence import GenePresence
        if name not in self.family_names:
            raise KeyError(f"no named family {name!r}; declared families are "
                           f"{sorted(self.family_names)}")
        return GenePresence(self, name)

    def has_family(self, node_id: int, name: str) -> bool:
        """Whether the named family ``name`` (declared via ``families=``) is present — has ≥ 1 copy — in
        the genome at ``node_id``. The presence a joint ``scaled_by("genomes:<name>", …)`` reads as its driver."""
        if name not in self.family_names:
            raise KeyError(f"no named family {name!r}; declared families are {sorted(self.family_names)}")
        fid = self.family_names[name]
        return any(c.family == fid for c in self.node_genomes[node_id])

    @cached_property
    def profiles(self) -> Profiles:
        """The phyletic profiles — each gene family's copy count in each extant species — derived
        from the observed genomes (the classic comparative-genomics matrix). See `profiles`."""
        extant = list(self.complete_tree.extant_leaves())
        if extant and not self.node_genomes:
            # a run reopened by `read_run` from a directory whose 'genomes' output was not written:
            # the genealogy is all there, the gene content is not. Say that, rather than KeyError.
            raise ValueError(
                "this run has no per-node gene content, so there are no profiles to derive — it was "
                "read back from a directory whose genomes.tsv was not written. Re-run the genomes "
                "level with 'genomes' among its outputs, or read profiles.tsv if that one is there. "
                "The gene trees and the event log are unaffected.")
        return profiles_from_genomes(self.node_genomes, extant)

    @cached_property
    def events(self) -> list[Event]:
        """The genome events — **one per row of ``genome_events.tsv``**, the same objects the writer
        formats.

        `edges` is the finer record this is grouped from: one entry per gene-tree *edge*, so a
        duplication is two of them and a transfer likewise. That is the shape a gene tree is built
        out of, and it used to be what this attribute returned — which meant counting duplications in
        Python gave twice the file's number, and a filter on ``kind == "transfer"`` matched everything
        here and nothing there. One word, one meaning: an event is what the log has a row for.
        """
        return events_from_edges(self.edges)

    @cached_property
    def gene_trees(self) -> dict[int, GeneTree]:
        """``{family id: GeneTree}`` — each family's true genealogy inside the complete tree,
        derived from the event log. Each ``GeneTree`` exposes ``.complete`` and ``.extant``. See
        `gene_trees`."""
        return gene_trees_from_edges(self.edges, self.complete_tree)

    def summary(self) -> dict:
        """What this run produced, as a plain dict — the payload of ``genome_summary.json``.

        **Every count here is one per event**, which is also what ``genome_events.tsv`` is now one row
        of. They are counted from the `GeneEdge` objects, and those are one per gene-tree *edge*: a
        duplication, a transfer and a speciation each end one gene and start two, so counting edges
        inflates them exactly 2×. A duplication's two edges share the ``parent`` gene they descend
        from, and a gene ends at exactly one event, so distinct parents *are* the events.

        ``loss`` counts every gene that died, which under ``replacement`` is more than the log's
        ``loss`` rows: a copy displaced by an arriving transfer has no row of its own — it is the
        second parent of that ``transfer_replacing`` row, because its death and the transfer are one
        event. It is a loss of the gene tree all the same, and that is what this counts.

        The family counts are the other thing nobody could reconcile: ``gene_trees/`` holds a file pair
        per family that ever existed, while the run's summary line counts the ones that survived, so
        "96 gene families" sat next to 213 files with nothing to explain the gap. Both numbers are
        here, named.

        ``origination`` counts only the families that arose *during* the run: the initial genome is
        logged as origination at the root's own start time, so a bare count of that kind is de-novo
        arrivals plus ``initial_families``, which is a number nobody asked for. They are separate here.

        And the cap, which was invisible. When ``max_family_size`` binds it discards duplications and
        arriving transfers, so realised rates fall below the declared ones — so this reports which
        families are sitting at it, because that is the signal a reader can act on.

        ``empty_genomes`` is the other end of the same story. There is **no floor** at this
        resolution: loss is counted per gene copy and the last copy is a copy like any other, so a
        high loss rate can strip a lineage of every gene it has. That is a real outcome of the model,
        not a failure — but it is invisible in the outputs, because a genome with no genes writes no
        row in ``profiles.tsv`` and leaves no gene tree for a sequence to run down. This is the
        number that says it happened, before a reader wonders why the matrix is short. (The ordered
        and nucleotide resolutions do have a floor, but it is a statement about what a chromosome is
        — a loss never takes a chromosome below its last gene — not a bound on genome size.)"""
        t0 = self.complete_tree.nodes[self.complete_tree.root].birth_time
        counted = event_counts(self.edges, t0)    # shared with the ordered and nucleotide summaries

        extant = list(self.complete_tree.extant_leaves())
        born = {e.family for e in self.edges}
        surviving = {c.family for i in extant for c in self.node_genomes.get(i, ())}
        genes_per_genome = [len(self.node_genomes.get(i, ())) for i in extant]
        cells = [collections.Counter(c.family for c in self.node_genomes.get(i, ())) for i in extant]
        copies = [n for cell in cells for n in cell.values()]

        cap = self.max_family_size
        at_cap = sorted({fam for cell in cells for fam, n in cell.items()
                         if cap is not None and n >= cap})
        return {
            "level": "genomes",
            "seed": self.seed,
            "resolution": "family",
            # one number per EVENT — the same thing a row of genome_events.tsv is
            "events": counted,
            "families": {"born": len(born), "surviving": len(surviving),
                         "died_out": len(born) - len(surviving),
                         "named": len(self.family_names)},
            "extant_genomes": len(extant),
            "empty_genomes": sum(1 for i in extant if not self.node_genomes.get(i, ())),
            "genes_per_genome": _stats(genes_per_genome),
            "copies_per_family_per_genome": _stats(copies),
            # the cap made visible. `families_at_cap` is what to look at: a family sitting at the
            # ceiling had events discarded, so its realised rates are below the ones you declared.
            "family_size_cap": {
                "cap": cap,
                "families_at_cap": len(at_cap),
                "cells_at_cap": sum(1 for cell in cells for n in cell.values()
                                    if cap is not None and n >= cap),
                "family_ids_at_cap": at_cap},
        }

    #: Every token ``write()`` honours — the write vocabulary, declared rather than left
    #: implicit in the method body. The CLI builds ``--write``'s choices from this, so the two
    #: cannot drift: they did, and `initial_sequence` and `species_tree` were writable from
    #: Python and unnameable on the command line.
    OUTPUTS = ("events", "profiles", "genomes", "initial_genome", "gene_trees",
               "species_tree", "summary")

    def write(self, directory, outputs=("events", "profiles", "genomes", "initial_genome",
                                        "gene_trees", "species_tree", "summary"), *,
              flat: bool = False) -> None:
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
        # An unknown token used to write nothing and exit clean — silent data loss you discover
        # three pipeline steps later, when the next tool has no input. The other levels have always
        # raised; these two did not.
        if unknown := [o for o in outputs if o not in self.OUTPUTS]:
            raise ValueError(f"unknown write outputs {unknown}; choose from {list(self.OUTPUTS)}")
        d = pathlib.Path(directory)
        d.mkdir(parents=True, exist_ok=True)
        # a run's directory describes that run: clear the per-unit directories this write is
        # about to fill, so nothing from a previous run survives inside them (see fresh_dirs)
        fresh_dirs(d, ("gene_trees",), flat)
        names = self.complete_tree.labels()      # e<id> for a lineage that died; n<id> for the rest
        if "events" in outputs:
            (d / "genome_events.tsv").write_text(events_tsv(self.edges, names), encoding="utf-8")
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
        if "summary" in outputs:
            write_summary(d / "genome_summary.json", self.summary())

    def _genomes_tsv(self) -> str:
        """Every node's gene content, one row per copy, in the order the genome holds them. The
        family counterpart of the ordered resolution's ``gene_order.tsv`` — without a chromosome
        or a position, because at this resolution a genome is a set, not a sequence."""
        cols = ("lineage", "family", "copy")
        names = self.complete_tree.labels()
        rows = [f"{names[s]}\t{c.family}\t{gene_label(c.id)}"
                for s in sorted(self.node_genomes)
                for c in sorted(self.node_genomes[s], key=lambda c: (c.family, c.id))]
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


def _pick_host(rng, gen, n_hosts: int) -> int:
    """A uniform pick among the lineages holding at least one copy — the per-lineage twin of
    `_pick_copy`.

    Under a per-lineage scope every occupied genome has the same chance whatever its size, which is
    the whole difference from per-copy. Empty genomes are skipped rather than picked and re-drawn:
    they contribute no share of the total either (see ``n_hosts``), and a lineage that contributed
    nothing must not be able to be chosen."""
    i = int(rng.integers(n_hosts))
    for k, g in enumerate(gen):
        if g:
            if i == 0:
                return k
            i -= 1
    raise AssertionError("n_hosts out of sync with the genomes")  # unreachable


def _pick_copy_by_family(rng, genome, mult: dict[int, float], fixed=None, unit: float = 1.0) -> int:
    """A copy index within one lineage, drawn in proportion to each copy's family weight.

    The within-lineage twin of `weighted_index()`. Needed whenever families carry different
    rates: the totals are summed with those weights, so the copy has to be drawn with them too,
    or the rate would say one thing and the picking another.

    ``fixed`` is the table of families that set a rate of **their own** (`GeneFamily`), and it is
    ``None`` unless some family in the run does. Then a copy's weight is whichever of the two applies
    to its family, which `_family_weights()` explains — and the arithmetic here is the same
    expression, one lineage at a time, so the pick and the total cannot disagree."""
    if fixed is None:
        weights = [mult[c.family] for c in genome]
    else:
        weights = [unit * mult[c.family] + fixed[c.family] for c in genome]
    total = sum(weights)
    r = float(rng.random()) * total
    acc = 0.0
    for j, w in enumerate(weights):
        acc += w
        if r < acc:
            return j
    return len(genome) - 1                    # float guard: r == total lands on the last copy


#: Key suffix for a `_FamilyWeights` table holding the families that set their own rate. The class is
#: generic over its keys and groups by table identity, so the two kinds of table ride in one dict.
_FIXED = "$fixed"


def _family_weights(unit: float, sums, fixed) -> list[float]:
    """One event's rate on every lineage: the run's unit rate times the multipliers of the families
    using it, **plus** the families that set a rate of their own.

    Two sums rather than one because the two are different quantities. A family with no rate of its
    own runs at the run's rate, times whatever it drew — a dimensionless multiplier, so the run's
    rate is factored out of the sum and applied once. A family that writes its own rate *is* that
    rate, in the rate's own units, so there is no run rate to factor out and nothing to multiply. That
    is why a written rate is never stored as a multiple of the run's: the encoding would need a
    division, and it would break the moment the run's rate is 0 or carries a verb of its own.

    ``fixed`` is ``None`` when no family in the run writes a rate, and then this is the expression it
    always was — the same floats in the same order, so such a run is unchanged to the last bit."""
    if fixed is None:
        return [unit * s for s in sums]
    return [unit * s + f for s, f in zip(sums, fixed)]


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
    (`simulate_genomes_family()` hands the same dict to each rate carrying no per-family draw of its
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


def _driven_weights(rate, gen, k_alive: int, t: float, drivers, fam_sums, own=None) -> list[float]:
    """One rate's per-lineage weights under a driver — and, when the run also carries a per-family
    draw, the family multipliers folded in.

    ``lineages=1`` reads one lineage's own share, so a PerLineage rate contributes its base and a
    PerCopy one its base times that lineage's copies. An empty genome is zeroed explicitly: under
    PerCopy `effective` already returns 0 for it, so this changes nothing there, and under PerLineage
    it is what stops a lineage with no victim taking a share.

    With ``fam_sums`` — that rate's per-lineage sums of the family multipliers over the live copies —
    the two weights **multiply**: the driver's factor belongs to the lineage and the multipliers to
    its contents, so the lineage's total is the unit rate read on this lineage times those sums.
    ``copies=1`` because the sums already count the copies (a family carrying no draw weighs 1).

    ``own`` is the per-lineage sum over the families that wrote a rate of their own, and it is
    **added** rather than multiplied. The driver scales the run's rate; a family that states its own
    rate has no run rate in it to scale, which is the same reading `_family_weights()` gives."""
    if fam_sums is None:
        return [rate.effective(copies=len(gen[k]), lineages=1, time=t, drivers=drivers[k])
                if gen[k] else 0.0 for k in range(k_alive)]
    return [rate.effective(copies=1, lineages=1, time=t, drivers=drivers[k]) * fam_sums[k]
            + (own[k] if own is not None else 0.0)
            if gen[k] else 0.0 for k in range(k_alive)]


def _driven_mods(rate) -> list:
    """The driven modifiers a rate carries — what ``scaled_by`` and ``set_by`` build — or ``[]`` when
    it is a plain number, a bare scope or a schedule. A non-empty list means the rate is
    *per-lineage*: each lineage's factor
    depends on the driver value on that branch, so the engine evaluates the rate lineage-by-lineage and
    picks the affected lineage weighted (the ``species_tree._grow`` shape). Each modifier's ``key``
    identifies it in the threaded ``drivers`` dict; its ``driver`` resolves to a trajectory."""
    return [m for m in rate.modifiers if isinstance(m, Driven)]


# --- the D/T/L/O mutators (each records to the event log; ids from the minters) -------------------

def _originate(genome, node, t, events, new_copy, new_family) -> None:
    """A new gene family arises: mint a founding copy in a fresh family and record it."""
    c = new_copy(new_family())
    genome.append(c)
    events.append(GeneEdge(t, "origination", node.id, c.family, c.id))


def _duplicate(genome, j, node, t, events, new_copy) -> None:
    """The gene at index ``j`` duplicates. In the ZOMBI1 per-segment model every event re-ids: the
    gene *ends* and **two** fresh copies descend from it, so both carry new ids (and the id in any
    node is that node's own)."""
    old = genome[j]
    cont, dup = new_copy(old.family), new_copy(old.family)
    genome[j] = cont                                   # the continuing lineage (a fresh id)
    genome.append(dup)                                 # the new copy (a fresh id)
    events.append(GeneEdge(t, "duplication", node.id, old.family, cont.id, parent=old.id))
    events.append(GeneEdge(t, "duplication", node.id, old.family, dup.id, parent=old.id))


def _lose_at(genome, j, node, t, events) -> None:
    """The copy at index ``j`` is lost (swap-remove — the genome is an order-agnostic multiset)."""
    lost = genome[j]
    genome[j] = genome[-1]
    genome.pop()
    events.append(GeneEdge(t, "loss", node.id, lost.family, lost.id))


def resolve_modules(modules, family_names) -> dict:
    """Validate ``modules={"flagellum": ["flgA", …]}`` — a **named group of named families**.

    A module changes nothing about how the genome evolves: it is a way of asking one question of the
    result, *how much of this group does a lineage carry* (`ModuleCompletion`). It is declared with
    the run rather than at read time so the grouping is part of the record — the report and the
    summary can name it, and two runs of the same command mean the same thing by "flagellum".

    Members must be declared families: an anonymous family has an integer id
    that is an artefact of the order events fired in, so a module built on one would not survive a
    change of seed.
    """
    if modules is None:
        return {}
    if not isinstance(modules, dict):
        raise TypeError(
            f"modules must be a dict of {{name: [family names]}}, got {type(modules).__name__}.")
    declared = set(family_names)
    out: dict[str, tuple[str, ...]] = {}
    for name, members in modules.items():
        if not isinstance(name, str) or not name.strip():
            raise ValueError(f"a module needs a non-empty name, got {name!r}")
        members = tuple(members)
        if not members:
            raise ValueError(f"module {name!r} has no families in it, so its completion would be "
                             f"undefined (0 out of 0).")
        if len(set(members)) != len(members):
            raise ValueError(f"module {name!r} names a family twice: {members}")
        missing = [m for m in members if m not in declared]
        if missing:
            raise ValueError(
                f"module {name!r} names {missing}, which are not declared families. Every member has "
                f"to be a declared family — an anonymous family's id comes from the order events "
                f"fired in, so a module built on one would mean something else at another seed. "
                f"Declared: {sorted(declared)}.")
        out[name] = members
    return out


def resolve_families(families, tree):
    """Everything a run declares about **named** families, from the one list that declares them.

    ``families=[family("IS1", loss=0.02), family("toxin", origin=("n5", 0.4)), …]`` says a family's
    name, its rates, where it starts and the group it belongs to. Those were four things once, said
    by ``family_names=``, ``origins=`` and ``modules=``; one declaration says all of them, and only
    it can carry rates.

    Returns ``(declared, module_map, planted)`` — the `GeneFamily` list in the order their ids are
    minted, the ``{module: (family name, …)}`` grouping `resolve_modules` builds, and the
    ``{index into declared: (time, lineage)}`` of the families given an ``origin``.
    """
    declared: list[GeneFamily] = []
    for spec in (list(families) if families is not None else []):
        if not isinstance(spec, GeneFamily):
            raise TypeError(
                f"families takes gene-family declarations — families=[family('IS1', loss=0.02), …] — "
                f"got {spec!r}. A bare name is family('IS1').")
        declared.append(spec)
    names = [f.name for f in declared]
    if len(set(names)) != len(names):
        raise ValueError(f"family names must be unique, got {names}")

    groups: dict[str, list[str]] = {}
    for f in declared:
        if f.module:
            groups.setdefault(f.module, []).append(f.name)
    module_map = resolve_modules(groups or None, names)

    # a declared family may be planted at a chosen point instead of starting at the origin
    with_origin = [(i, f) for i, f in enumerate(declared) if f.origin is not None]
    resolved = resolve_origins([f.origin for _i, f in with_origin], tree) if with_origin else []
    return declared, module_map, {i: pair for (i, _f), pair in zip(with_origin, resolved)}


def resolve_family_rates(declared, run_rates):
    """The per-copy rate each declared family sets for **itself**, by event and family index.

    ``run_rates`` is ``{event: the run's resolved Rate}``. A family that writes nothing for an event
    is absent here and runs at the run's rate; a family that writes one runs at exactly that, in the
    rate's own units, with no multiplier from the run applied (SPEC §5's argument for ``set_by``,
    read one level down).

    Two things are refused for now, and both come back at the joint step. A family's rate must be
    `PerCopy`, because the run's own rate for these three events is what it is summed beside. And it
    must carry no verb: the value here is read once, before the run starts, while a verb makes it a
    function of the context.
    """
    out: dict[str, dict[int, float]] = {}
    for i, f in enumerate(declared):
        for key, spec in f.written().items():
            rate = as_rate(spec, default_scope=PerCopy, label=f"{f.name}'s {key}")
            if rate.scope is not PerCopy:
                raise ValueError(
                    f"family {f.name!r} writes a {rate.scope.__name__} {key} rate, but a family's own "
                    f"rate is per copy — it is summed beside the run's rate for the same event, and "
                    f"that is counted per copy. Write PerCopy(...), or a bare number.")
            if rate.modifiers:
                raise ValueError(
                    f"family {f.name!r}'s {key} carries {describe(rate.modifiers[0])}, which a "
                    f"family's own rate does not take yet: it is read once before the run starts, and "
                    f"a verb makes it a function of the context. Write a plain number here, and put "
                    f"the verb on the run's {key} to move every family together.")
            rate.check_one_base(f"family {f.name!r}'s {key}")
            out.setdefault(key, {})[i] = rate.effective(copies=1, lineages=1, time=0.0)
    return out


#: the live gene-content driver reading a lineage's whole gene count, as `zombi2.joint` spells it
LIVE_COUNT = "genomes:count"

#: A joint genome run's ceiling on live copies. A rate that reads the genome's own content can feed
#: itself — more copies raise the rate, which makes more copies — so a run that looks calm on paper
#: can have no realistic end. It RAISES rather than stopping early, for the reason the species engine
#: gives: a run cut off at a size is no longer a sample from the process asked for.
MAX_LIVE_COPIES = 2_000_000


def resolve_live_drivers(mods, declared_names, *, joint: bool) -> list[str]:
    """Validate the **live** drivers a run's rates read, and return their keys.

    A live driver names gene content growing in this same run — ``"genomes:count"`` for a lineage's
    whole gene count, ``"genomes:<name>"`` for whether a declared family is there. That makes the run
    joint (SPEC §2): the driver cannot be finished first, because it is what the run is producing.

    ``joint`` is the run's own declaration, and it is checked both ways. Asking for a joint run with
    nothing reading a live driver is as much a mistake as reading one without saying so.
    """
    keys = []
    for m in mods:
        src = m.driver
        if src == LIVE_COUNT:
            check_mapping_fires(m.mapping, {0}, driver_label=f"the driver {src!r}")
        else:
            name = src.split(":", 1)[1]
            if name not in declared_names:
                raise ValueError(
                    f'scaled_by("{src}", ...) reads family {name!r}, which this run does not declare '
                    f"— add families=[…, family({name!r})]. Declared: {sorted(declared_names)}.")
            check_mapping_fires(m.mapping, {"present", "absent"},
                                driver_label=f"the driver {src!r}", exhaustive=True)
        keys.append(src)
    if keys and not joint:
        raise ValueError(
            f"a rate reads {sorted(set(keys))}, which is gene content this same run is producing, so "
            f"the run is joint — neither the driver nor what it drives can be finished first. Say so "
            f"with joint=True. To read gene content grown by an EARLIER run instead, pass that run's "
            f"presence(...) rather than a name, which is conditioning.")
    if joint and not keys:
        raise ValueError(
            "joint=True says two things in this run drive each other, but no rate reads live gene "
            'content. Give a rate a scaled_by("genomes:<family>", …) or scaled_by("genomes:count", …), '
            "or drop joint=True.")
    return keys


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
    if isinstance(max_family_size, Rate):
        # A scope constructor returns a `Rate` — `PerLineage(10)` is one, not a `Scope` instance —
        # so this is where a written scope arrives, and the scope itself is the class on it.
        name = max_family_size.scope.__name__ if max_family_size.scope is not None else "a scope"
        was = (f"{name}(n) was n × the size of the species tree"
               if max_family_size.scope is PerLineage else f"{name}(n) was exactly n")
        number = f"{max_family_size.base:g}" if max_family_size.base is not None else "n"
        raise ValueError(
            f"max_family_size is a plain count of copies in one genome — write "
            f"max_family_size={number}, not {max_family_size!r}. A cap on copies in "
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


def resolve_origins(origins, tree) -> list[tuple[float, int]]:
    """Validate ``origins=[("n5", 0.4), …]`` — families that originate at a **chosen** point of the
    tree — into ``[(time, lineage), …]``, **in the order they were written** (which is the order
    their family ids are then minted in, so the third origin is the third placed family whether or
    not it is the third to happen).

    Origination is otherwise drawn: the ``origination`` rate says how often a family arises and the
    engine picks where. This says *here*, on this branch at this time, which is the question a run
    asks when the family is the subject rather than a sample — one gene family born on the branch
    leading to a clade, and what duplication, transfer and loss then do to it.

    It **adds to** the run rather than replacing part of it: whatever ``initial_families`` and
    ``origination`` produce is produced as well, so ``origins`` alone (with both of those at 0) is a
    tree carrying exactly the families you placed, and ``origins`` beside them is an ordinary genome
    with one family planted where you want it. Nothing is switched off behind your back.

    A lineage is written as the tree writes it — ``n5``, or ``e5`` for one that went extinct — or as
    the bare node id. The time is the run's own clock (the origin at 0, as ``changing_at`` and
    `Time` read it), and must fall **inside** that lineage's life: the branch's end is where it is
    already over and its daughters have started, so a family originating there belongs to one of
    them. ``None`` puts it at the branch's start.
    """
    if origins is None:
        return []
    # the single-pair slip, caught here rather than as "too many values to unpack" three frames down
    if (isinstance(origins, tuple) and len(origins) == 2
            and not isinstance(origins[0], (list, tuple))):
        raise ValueError(f"origins takes a list of (lineage, time) pairs — write "
                         f"origins=[{origins!r}] to plant the one family.")
    names = tree.labels()
    out: list[tuple[float, int]] = []
    for item in origins:
        try:
            lineage, time = item
        except (TypeError, ValueError):
            raise ValueError(f"each origin is a (lineage, time) pair, got {item!r}") from None
        if isinstance(lineage, str):
            try:
                node_id = node_from_label(lineage)
            except ValueError:
                raise ValueError(
                    f"origin lineage {lineage!r} is not a lineage label; write it as the tree writes "
                    f"it — n5, or e5 for a lineage that went extinct — or as the bare node id.") from None
        elif isinstance(lineage, bool) or not isinstance(lineage, int):
            raise ValueError(f"origin lineage must be a label like 'n5' or a node id, got {lineage!r}")
        else:
            node_id = lineage
        if node_id not in tree.nodes:
            # named rather than given as a range: a pruned tree keeps the ids it had, so its
            # lineages need not be a contiguous block and "0–22" would then be a wrong answer
            raise ValueError(
                f"origin lineage {lineage!r} is not one of this tree's {len(tree.nodes)} lineages, "
                f"which are named as the tree writes them "
                f"({', '.join(names[i] for i in sorted(tree.nodes)[:4])}, …). A genome run walks the "
                f"COMPLETE tree, so the lineages that went extinct — the e-prefixed ones — are in it "
                f"too.")
        node = tree.nodes[node_id]
        if time is None:
            time = node.birth_time
        elif isinstance(time, bool) or not isinstance(time, (int, float)):
            raise ValueError(f"origin time must be a number on the run's clock (or None for the "
                             f"start of the branch), got {time!r}")
        if not node.birth_time <= time < node.end_time:
            raise ValueError(
                f"origin time {float(time):g} is outside lineage {names[node_id]}, which runs from "
                f"{node.birth_time:g} to {node.end_time:g}. The end is excluded: the branch is over "
                f"there and its daughters have begun, so a family originating at that instant "
                f"belongs to one of them.")
        out.append((float(time), node_id))
    return out


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

    def holds(self, k: int, family: int) -> bool:
        """Whether lineage ``k`` carries this family right now — what a live ``"genomes:<name>"``
        driver reads. Exact rather than a scan, because `removed` drops a family's key when its last
        copy goes, so presence is the key being there."""
        return family in self._counts[k]

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
    if transfer_to == "uniform":
        # O(1) uniform recipient: the same single draw recipient_index makes as
        # ``cand[rng.integers(len(cand))]``, where ``cand`` is every alive lineage except the donor
        # (unless self_transfer). Skipping the donor's slot is a +1 index shift, so this returns the
        # identical recipient without allocating the candidate list every transfer — a small, byte-
        # identical simplification of the undriven hot path.
        m = len(alive) if self_transfer else len(alive) - 1
        i = int(rng.integers(m))
        kr = i if (self_transfer or i < kd) else i + 1
    else:  # weighted rules (Distance / Clades / Driven) weigh every candidate — inherently O(alive)
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
    replaced = None
    if replacement:
        residents = [p for p, c in enumerate(rg) if c.family == fam and c.id != cont.id]
        if residents:  # homologous overwrite; empty ⇒ additive fallback (the gene still arrives)
            p = residents[int(rng.integers(len(residents)))]
            victim = rg[p]
            rg[p] = rg[-1]
            rg.pop()
            # the displaced copy dies, so it is a loss of the gene tree like any other. It is *also*
            # part of this transfer, and naming it on the transfer is what keeps the two linked: the
            # log writes one `transfer_replacing` row carrying both parents, and nobody downstream has
            # to recognise the pair by their shared timestamp.
            replaced = victim.id
            events.append(GeneEdge(t, "loss", recipient, fam, replaced))
            counts.removed(kr, fam)
            delta = 0
    rg.append(xfer)
    counts.added(kr, fam)
    events.append(GeneEdge(t, "transfer", donor, fam, cont.id, parent=src.id, donor=donor,
                        replaced=replaced))
    events.append(GeneEdge(t, "transfer", recipient, fam, xfer.id, parent=src.id, recipient=recipient,
                        donor=donor, replaced=replaced))
    return delta, kr


@without_cyclic_gc
def simulate_genomes_family(tree, *, duplication=0.0, transfer=0.0, loss=0.0, origination=0.0,
                            transfer_to="uniform", replacement=False, self_transfer=False,
                            initial_families=100, families=None, max_family_size=10, joint=False,
                            seed=None, parallel=False, stream_to=None, outputs=None,
                            progress=False, **retired) -> "FamilyGenomesResult | StreamedRun":
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
    or ``Recipients().weighted_by(driver, mapping)`` (weighted by an evolved value; see below). ``replacement=True``
    overwrites a homologous
    copy in the recipient (additive fallback if it has none); ``self_transfer=True`` lets a lineage
    donate to itself. The root starts with ``initial_families`` families of one copy each, recorded
    as originations at the origin. ``families=[family("toxin")]`` additionally declares **named** families —
    each gets a normal (integer) family id, but its name is remembered in ``result.family_names`` so
    you can track a specific family (``result.has_family(node, "toxin")``); this is the handle a joint
    ``scaled_by("genomes:toxin", …)`` reads. Deterministic given ``seed``.

    **A family with rates of its own.** ``families=[family("IS1", transfer=PerCopy(1.5), loss=0.02)]``
    declares a named family and gives it its own rates, so one family can behave nothing like the rest
    of the genome — a mobile element that transfers constantly, a core gene that is almost never lost.
    Whatever a family leaves out falls back to the run's rate for that event, so
    ``families=[family("toxin")]`` is a family declared for its name alone and is exactly
    ``families=[family("toxin")]``. `family()` also takes ``origin=`` and ``module=``, which are
    ``origins=`` and ``modules=`` said on the family itself. Origination takes no per-family value:
    when it is read the family does not exist yet to have one.

    A written rate **is** that family's rate, in the rate's own units. The run's rate does not reach
    it, and neither does a ``varying_among('families', …)`` draw meant to vary the run's rate among
    families (SPEC §5's argument for ``set_by``, one level down). Two things are refused for now: a
    family's rate carrying a verb of its own, and a family's rate beside a *driven* run rate; both
    arrive with the joint step. So does the per-family engine — ``parallel=`` and ``stream_to=``
    build one set of rates for the whole run, so they refuse a family that writes its own.

    **A family placed by hand.** ``origins=[("n5", 0.4)]`` originates a family on lineage ``n5`` at
    time ``0.4`` — the same event the ``origination`` rate produces, at a point you choose rather
    than one that is drawn. It **adds to** the run: whatever ``initial_families`` and ``origination``
    give you is still there, so ``origins`` with both of those at 0 is a tree carrying exactly the
    families you placed, and ``origins`` beside them is an ordinary genome with one family planted
    where you want it. The time is the run's own clock and must fall inside that lineage's life
    (``None`` puts it at the branch's start); a placed family is an ordinary family from that instant
    on, so it duplicates, transfers and is lost like any other, and gets its gene tree the same way.
    The ids come straight after the initial and named ones, in the order you wrote the origins, on
    either engine — so the family you placed is ``initial_families + len(family_names)``. See
    `resolve_origins`.

    **Conditioning (a trait drives a rate).** Any of the four rates may be *driven by another level* —
    ``loss = PerCopy(0.25).scaled_by("trait_events.tsv", {"aquatic": 3.0, "terrestrial": 1.0})`` scales each
    lineage's loss by the habitat on that branch, read from a driver file grown first
    (``traits.simulate_discrete(...).write(dir, outputs=("events",))``, which writes
    ``trait_events.tsv``). A driven rate is then *per-lineage*: it is summed over the living lineages (each with its own copy count and driver
    value), the affected lineage is drawn weighted by its rate, and the Gillespie steps at every
    mid-branch switch of the driver (SPEC §2). For ``transfer`` the affected
    lineage is the **donor**, so a driven ``transfer`` says how often a lineage *donates*.

    **Conditioning (a trait drives who receives).** ``transfer_to =
    Recipients().weighted_by(driver, mapping)`` is
    the other half, and a different model: the mapping's numbers are per-candidate **weights**, not
    rate multipliers, so they leave the total amount of transfer alone and only redistribute it
    (SPEC §5, a weight, not a rate). Candidate lineage ``k`` gets weight ``mapping(driver value on k now)``
    and receives with probability ``w_k / Σw`` — five candidates at weight 1 and five at weight 2 send
    two thirds of transfers to the weight-2 group. Weight 0 means "cannot receive"; when every
    candidate weighs 0 the transfer does not happen (see `_do_transfer()`). The two driven arguments are
    independent and may be used together or apart.

    ``parallel`` opts into a **separate** engine that evolves the families concurrently, one per worker
    process — the families are independent (a transfer roams a copy across lineages but never mixes two
    families), so it enumerates every family's origination first and then evolves each on its own. It is
    worker-count invariant (each family draws from its own stream spawned off ``seed``) but gives a
    different-though-equally-valid draw than the serial default for a given seed. ``False`` (default)
    runs the serial loop; ``True`` uses every core; an ``int`` sets the worker count. A **driven** rate
    or a weighted ``transfer_to`` runs on the per-family engine too — conditioning does not couple
    families, so nothing here forces a fallback. The gain is real but modest (a merge over the whole event log stays
    serial), so a handful of workers is the sweet spot; unlike the sequences level it does not scale far.
    Because it spawns processes, a calling script must guard its entry with ``if __name__ ==
    "__main__":`` (the ``zombi2`` CLI already does).

    ``stream_to=DIR`` takes the same engine to the many-families regime: each family is written straight
    to disk as it finishes — no whole-run merge, no run held in memory (a run that fills gigabytes in
    memory streams in tens of megabytes) — and a light `StreamedRun` handle comes
    back instead of a ``FamilyGenomesResult``. ``outputs=`` picks which files, as
    `FamilyGenomesResult.write()` takes them minus ``summary`` (a streamed run writes no
    ``genome_summary.json``); the default is all six. It is the per-family engine, and
    ``outputs`` without ``stream_to`` is an error.
    """
    tree = as_tree(tree, level="genomes")
    dup = as_rate(duplication, default_scope=PerCopy)
    tra = as_rate(transfer, default_scope=PerCopy)
    los = as_rate(loss, default_scope=PerCopy)
    org = as_rate(origination, default_scope=PerLineage)
    # A scope says *per what* a rate is counted, and the three copy-consuming events implement two
    # answers. **Per copy** (the default) puts every copy independently at risk, so a bigger genome
    # turns over faster. **Per lineage** is a fixed budget: the lineage duplicates or loses at its
    # rate whatever its genome size — and then the total and the pick both have to change, or the
    # rate would say one thing and the picking another. Origination is per lineage alone, because it
    # is the rate at which families are CREATED: per copy it would be base × 0 in an empty genome, a
    # silent no-op. A driver is read per lineage on all four events; on transfer the driven
    # lineage is the DONOR (who receives is the separate transfer_to choice, below).
    for label, rate, legal in (("duplication", dup, (PerCopy, PerLineage)),
                               ("transfer", tra, (PerCopy, PerLineage)),
                               ("loss", los, (PerCopy, PerLineage)),
                               ("origination", org, (PerLineage,))):
        # `rate.scope` holds the scope CLASS, not an instance — a scope constructor returns the rate
        # itself — so this is an identity test against the legal set rather than an isinstance one.
        assert rate.scope is not None            # as_rate fills the level's default where none was written
        if rate.scope not in legal:
            raise ValueError(
                f"{label} has a {rate.scope.__name__} scope, but the family genome engine "
                f"takes {' or '.join(s.__name__ for s in legal)} for {label}."
            )
        for m in rate.modifiers:
            if m.reads == (DRAWN, "families") and label == "origination":
                raise ValueError(
                    "origination carries a per-family draw, but origination is the rate at which families are "
                    "CREATED — when it is read there is no family yet to have drawn a factor for. "
                    "Write varying_among('families', …) on duplication, transfer or loss; writing one "
                    "such object on several of them gives a family-wide tempo, since one object is "
                    "one draw.")
            if isinstance(m, Driven):
                check_not_a_kernel(m.mapping, label=label)
            if is_implemented(m, IMPLEMENTED_MODIFIERS, "genomes.family"):
                continue
            raise ValueError(
                f"{label} carries {describe(m)}, which the family genome engine does not "
                f"support. It takes changing_at (skyline), scaled_by (a conditioned or joint driver), "
                f"set_by (a driver that replaces the base) and varying_among('families', …) "
                f"(per-family heterogeneity). Clade drift is not implemented yet."
            )
    for label, rate in (("duplication", dup), ("transfer", tra), ("loss", los),
                        ("origination", org)):
        rate.check_one_base(label)
    # A per-family draw under a per-lineage scope is refused because what it should MEAN is a
    # modelling decision nobody has taken, not because it is hard. Under PerCopy a family's
    # multiplier scales each copy's own rate, so a genome full of fast families turns over faster —
    # the multiplier moves the total. Under PerLineage the total is fixed by definition, so the same
    # multiplier could only decide *which* copy the event takes, normalised within the lineage. Those
    # are two different models, and running one while the user wrote the other is exactly the silent
    # mismatch this engine refuses everywhere else.
    # The check is over the whole RUN, not per rate, and getting that wrong was a real bug: a
    # per-family draw anywhere makes the engine take its per-family path for **every** gene rate,
    # summing each one over the live copies — so a `PerLineage` rate elsewhere in the same run was
    # silently counted per copy while nothing on the page said so. One draw on duplication was
    # enough to turn a `PerLineage` loss back into a `PerCopy` one.
    per_lineage_here = [lbl for lbl, r in (("duplication", dup), ("transfer", tra), ("loss", los))
                        if r.scope is PerLineage]
    drawn_here = [lbl for lbl, r in (("duplication", dup), ("transfer", tra), ("loss", los))
                  if any(m.reads == (DRAWN, "families") for m in r.modifiers)]
    if per_lineage_here and drawn_here:
        raise ValueError(
            f"{', '.join(per_lineage_here)} is PerLineage while {', '.join(drawn_here)} carries a "
            f"per-family draw, and the two cannot share a run. Under PerCopy a family's multiplier "
            f"scales each copy's rate, so it changes the lineage's total; under PerLineage the total "
            f"is fixed whatever the genome holds, so the multiplier could only choose which copy is "
            f"taken. Those are different models and the choice is not made yet — write PerCopy "
            f"throughout for the first, or drop the per-family draw for the second.")
    # the choice (SPEC §5), validated in the one place all three resolutions share: the mapping's
    # numbers are weights over the candidate recipients, never a rate multiplier
    transfer_to = resolve_transfer_to(transfer_to)
    if isinstance(initial_families, bool) or not isinstance(initial_families, int) or initial_families < 0:
        raise ValueError(f"initial_families must be a non-negative integer, got {initial_families!r}")
    check_no_retired_keywords(retired, where="simulate_genomes_family")
    # every named family, from the one list that declares them
    declared, module_map, planted_named = resolve_families(families, tree)
    family_names = [f.name for f in declared]
    # what each family sets for itself; empty unless some family writes a rate, and then the engine
    # takes the path that sums those beside the run's (see `_family_weights`)
    fam_own = resolve_family_rates(declared, {"duplication": dup, "transfer": tra, "loss": los})
    any_written = bool(fam_own)
    if any_written:
        for key, rate in (("duplication", dup), ("transfer", tra), ("loss", los)):
            if key not in fam_own:
                continue
            if rate.scope is not PerCopy:
                raise ValueError(
                    f"a family writes its own {key}, but the run's {key} is {rate.scope.__name__}. "
                    f"The two are summed over the same copies, so both are counted per copy — write "
                    f"PerCopy for the run's {key}, or drop the family's.")
            # a driven run rate composes: the driver scales the run's rate, and a family that states
            # its own has no run rate in it to scale (see `_driven_weights`)

    # A family's copies in one genome are capped. Growth compounds — a duplication rate above the
    # loss rate multiplies without bound — so a run needs a ceiling somewhere. An int is that number
    # of copies; a float (the default, 10.0) is that multiple of the lineages in the complete tree,
    # so the bound travels with the size of the run. Refusing an event on a condition that depends
    # only on the current state is Poisson thinning, so what is kept is exactly the process whose
    # duplication rate is zero for a family already at its quota — a declared ceiling, not a
    # truncated run. ``None`` removes it.
    cap = resolve_max_family_size(max_family_size)

    # conditioning: a rate written with scaled_by reads a driver per lineage. Resolve each driver once into
    # a DriverTrajectory (value + next-switch lookups, keyed by the shared species node id) — from a
    # file (a str driver) or an object handed over in memory (a trait result, or a genome's presence /
    # completion). No driven rate ⇒ this is empty
    # and the loop stays byte-identical to an undriven run.
    dup_mods, los_mods = _driven_mods(dup), _driven_mods(los)
    org_mods, tra_mods = _driven_mods(org), _driven_mods(tra)
    all_mods = (*dup_mods, *los_mods, *org_mods, *tra_mods)
    # Two kinds of driver, told apart by what the driver *is* (SPEC §5). A finished one — a file, a
    # grown result — was produced before this run started, and the run is conditioned. A **live**
    # name is gene content this run is itself producing, and the run is joint: there is no order to
    # simulate the two in, because they are the same run. Only the first can be resolved to a
    # trajectory; the second is read off the live genome as the loop goes.
    live_mods = [m for m in all_mods if names_a_live_level(m.driver)]
    file_mods = [m for m in all_mods if not names_a_live_level(m.driver)]
    live_keys = resolve_live_drivers(live_mods, set(family_names), joint=joint)
    # driver key → its Driven (deduped, so a driver shared across rates resolves once);
    # the modifier rather than the driver itself, because the driver's step rides on the modifier
    by_key: dict[object, "Driven"] = {}
    for m in file_mods:
        by_key.setdefault(m.key, m)
    resolved = {}
    if by_key:
        resolved = {key: resolve_driver(m.driver, tree, step=m.step, level="genomes.family")
                    for key, m in by_key.items()}
        # a mapping whose states never occur in the driver leaves every lineage at the default factor,
        # so the rate is never driven and the run is secretly the undriven model — refuse it here,
        # naming the driver, rather than let it pass as a driven run
        for m in file_mods:
            label = m.driver if isinstance(m.driver, str) else f"<{type(m.driver).__name__}>"
            check_mapping_fires(m.mapping, resolved[m.key].states(), driver_label=label)
    # `trajs` is the drivers that move a RATE: they alone make the loop per-lineage and set the
    # Gillespie horizon. It is built BEFORE transfer_to is prepared, and that order is
    # load-bearing — a driven transfer_to changes no rate, so its trajectory must not end up here
    # adding horizon breakpoints (see prepare_transfer_to). `resolved` is passed along as the driver
    # cache, so a driver shared between a rate and transfer_to is loaded once.
    trajs = dict(resolved)
    group_of, to_traj = prepare_transfer_to(tree, transfer_to, resolved, level="genomes.family")

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
    seed = resolve_seed(seed)     # drawn if none was given, so either engine below records it
    if (parallel or stream_to is not None) and live_keys:
        # The one thing that engine's whole design rests on: a family's history depends on no other
        # family, so each can be evolved alone. A rate reading the genome's own content is exactly
        # that dependence, so this is a refusal rather than a fallback.
        raise ValueError(
            "a joint genome run cannot use the per-family engine (parallel= / stream_to=), which "
            "evolves each family in its own process because families do not affect each other. A "
            "rate reading live gene content is that effect. Drop parallel / stream_to.")
    if (parallel or stream_to is not None) and planted_named:
        # Pass 1 of that engine enumerates every family's origination up front, seeding the declared
        # ones at the root; a family that arrives partway down is not in that enumeration, and
        # threading it through would renumber the families the serial engine mints. Stated rather
        # than silently dropped.
        raise ValueError(
            "family(origin=...) does not run on the per-family engine (parallel= / stream_to=), "
            "which enumerates every family's origination before it starts. Drop parallel / "
            "stream_to, or let the origination rate place the family.")
    if (parallel or stream_to is not None) and any_written:
        # not a fallback: that engine evolves one family per process against a context built once for
        # the whole run, so a per-family rate has nowhere to live in it. Saying so beats silently
        # running every family at the run's rate.
        raise ValueError(
            "a family writing its own rate cannot run on the per-family engine (parallel= / "
            "stream_to=), which builds one set of rates for the whole run and evolves each family "
            "against it. Drop parallel / stream_to, or give every family the run's rate.")
    if parallel or stream_to is not None:
        from ._perfamily import run_parallel_family
        result = run_parallel_family(
            tree, dup=dup, tra=tra, los=los, org=org, transfer_to=transfer_to,
            replacement=replacement, self_transfer=self_transfer, initial_families=initial_families,
            family_names=family_names, placed=[], modules=module_map, cap=cap,
            seed=seed, parallel=parallel,
            progress=progress, stream_to=stream_to, outputs=outputs,
            trajs=trajs, to_traj=to_traj, group_of=group_of,
            driven={"duplication": bool(dup_mods), "transfer": bool(tra_mods),
                    "loss": bool(los_mods), "origination": bool(org_mods)})
        if result is not None:
            return result

    rng, seed = stream("genomes", seed)     # the genomes level's own stream
    copy_counter = 0
    family_counter = 0

    def new_copy(family: int) -> GeneCopy:
        nonlocal copy_counter
        c = GeneCopy(copy_counter, family)
        copy_counter += 1
        return c

    # Per-family multipliers, drawn once when a family is created and then fixed for its whole life.
    # Whether a family's rates move together is decided by what was written: one `Random` object read
    # by two rates is one draw for both, two objects are two draws. Empty unless some rate carries
    # one, and then the engine takes its weighted path; a run carrying none draws nothing here.
    fam_by = {"duplication": tuple(m for m, _ in dup.carried_modifiers(unit="families")),
              "transfer": tuple(m for m, _ in tra.carried_modifiers(unit="families")),
              "loss": tuple(m for m, _ in los.carried_modifiers(unit="families"))}
    any_family = any(fam_by.values()) or any_written
    # A rate carrying nothing per family holds 1.0 for every family, so all such rates share one
    # empty table rather than each filling its own — which is what lets _FamilyWeights sum them once.
    # Sharing is off once some family writes a rate: a family may write its loss and not its
    # duplication, so the two tables then hold different numbers for the same family.
    no_variation: dict[int, float] = {}
    fam_mult: dict[str, dict[int, float]] = {
        key: ({} if (mods or any_written) else no_variation) for key, mods in fam_by.items()}
    #: per event, the families that set a rate of their own — their own per-copy rate, and 0.0 for
    #: every other family, which contributes through `fam_mult` instead. `None` when nobody does,
    #: and then every expression below is the one it always was (see `_family_weights`).
    fam_fixed: "dict[str, dict[int, float]] | None" = (
        {key: {} for key in fam_by} if any_written else None)
    #: that table for one event, or ``None`` — what a copy pick reads, against the per-lineage *sums*
    #: of the same table that `own_sums` reads. Two different shapes of the same information: a rate
    #: per family here, and per lineage the total over its live copies there.
    own_rates = (lambda key: fam_fixed[key]) if fam_fixed is not None else (lambda key: None)

    def new_family(declared_at: "int | None" = None) -> int:
        """Mint a family id. ``declared_at`` is its index in ``declared`` for a named family, which is
        how it finds the rates it wrote; an anonymous family passes nothing and runs at the run's."""
        nonlocal family_counter
        f = family_counter
        family_counter += 1
        if any_family:
            # one draw per distinct modifier *object* for this family, shared across its rates: the
            # same `Random` object written on duplication and on loss means one number, so a fast
            # family is fast at both. Two separately built ones are two draws even with the same law.
            shared: dict[int, float] = {}
            for key, mods in fam_by.items():
                own = None if declared_at is None else fam_own.get(key, {}).get(declared_at)
                if own is not None:
                    # this family's rate IS the number it wrote, so the run's rate does not reach it
                    # and neither does a draw meant to vary the run's rate among families
                    fam_mult[key][f] = 0.0
                    fam_fixed[key][f] = own
                    continue
                fam_mult[key][f] = math.prod(values_at_birth(mods, rng, shared))
                if fam_fixed is not None:
                    fam_fixed[key][f] = 0.0
        return f

    # Which of the three copy-consuming rates is a fixed per-lineage budget rather than a per-copy
    # risk. Read once: it decides both how the total is counted and how the affected lineage is
    # picked, and those two must never disagree.
    dup_per_lineage = dup.scope is PerLineage
    los_per_lineage = los.scope is PerLineage
    tra_per_lineage = tra.scope is PerLineage
    any_per_lineage = dup_per_lineage or los_per_lineage or tra_per_lineage

    depth = mean_root_to_tip(tree)  # timescale for Distance weighting (unused by "uniform")
    schedule = sorted((tree.nodes[i].end_time, i) for i in tree.nodes)  # (end_time, node_id)

    root = tree.nodes[tree.root]
    t = root.birth_time
    alive: list[int] = []
    gen: list[list[GeneCopy]] = []
    pos: dict[int, int] = {}
    genomes: dict[int, tuple[GeneCopy, ...]] = {}
    events: list[GeneEdge] = []
    enter(alive, gen, pos, root.id, [])
    for _ in range(initial_families):  # lay down the origin's genome as originations at t = root.birth_time
        _originate(gen[0], root, t, events, new_copy, new_family)
    named: dict[str, int] = {}  # a minted id per declared name (so GeneCopy.family stays an int)
    named_plants: list[tuple[float, int, int]] = []
    for i, spec in enumerate(declared):
        fid = new_family(i)
        named[spec.name] = fid
        if i in planted_named:
            # a declared family given an `origin` is planted there rather than seeded at the origin,
            # which is what `origins=` does for an anonymous one — the same event, with a name on it
            t_p, lineage = planted_named[i]
            named_plants.append((t_p, lineage, fid))
            continue
        c = new_copy(fid)
        gen[0].append(c)
        events.append(GeneEdge(t, "origination", root.id, fid, c.id))
    # A planted family is not seeded into a genome: it arrives at its own time, in the loop below.
    plants = sorted(named_plants)
    plant_i = 0                                                              # walked in time order
    total_copies = len(gen[0])
    initial_genome = tuple(gen[0])   # the run's starting genome: a snapshot before the stem runs

    any_driven = bool(trajs) or bool(live_keys)

    def live_value(src: str, k: int):
        """What a live driver reads on lineage ``k`` **right now** — the joint half of the driver
        mechanism (SPEC §2). A finished driver answers from a trajectory built before the run; this
        one answers from the genome the run is building.

        It needs no horizon breakpoint, and that is what makes the race exact rather than thinned:
        gene content changes only when a genome event fires, and an event ends the current step, so
        every rate is already constant between two events."""
        if src == LIVE_COUNT:
            return len(gen[k])
        return "present" if counts.holds(k, named[src.split(":", 1)[1]]) else "absent"
    # the per-family weight sums, carried across events rather than rebuilt each time (see the class).
    # The families that wrote their own rate ride in the same structure under a suffixed key, because
    # summing a table over a lineage's copies is the same work whichever kind of number is in it.
    _tables = dict(fam_mult)
    if fam_fixed is not None:
        _tables.update({key + _FIXED: m for key, m in fam_fixed.items()})
    weights = _FamilyWeights(_tables, gen) if any_family else None
    counts = _FamilyCounts(gen)      # the family cap's question, answered without walking a genome

    # the species tree's schedule is the run's spine: one entry per speciation/extinction, so how
    # far through it we are is how far through the tree the genomes have got
    bar = progress_bar(len(schedule), "genomes", unit="branch", enabled=progress)
    si = 0
    while si < len(schedule):
        bar.to(si)
        n = total_copies
        if live_keys and n > MAX_LIVE_COPIES:
            raise RuntimeError(
                f"the run passed {MAX_LIVE_COPIES} live gene copies at time {t:.3g} and is still "
                f"growing — a rate reading the genome's own content is feeding itself. Lower the "
                f"rates, flatten the mapping the driver is read through, or set a max_family_size.")
        k_alive = len(alive)
        ctx = {"copies": n, "lineages": k_alive, "time": t}
        # A copy-consuming event counted *per lineage* is counted per lineage that HOLDS a copy: an
        # empty genome offers nothing to duplicate or lose, so it must not contribute its share of
        # the total and then be picked with no victim inside it. Origination keeps `ctx` — an empty
        # genome can still gain a family. Computed only when some rate needs it, so the per-copy
        # path does exactly the work it did before.
        if any_per_lineage:
            n_hosts = sum(1 for g in gen if g)
            host_ctx = {"copies": n, "lineages": n_hosts, "time": t}
        else:
            n_hosts, host_ctx = 0, ctx
        can_xfer = n > 0 and (k_alive >= 2 or self_transfer)  # a recipient must be able to exist
        # a driven rate is per-lineage: sum its effective rate over the living lineages (each read with
        # its own copy count and its branch's driver value), keeping the weights for the affected-lineage
        # pick — the species_tree._grow shape. An undriven rate stays pooled (one .effective, uniform
        # pick), so a run with no driver is byte-identical to before. For transfer the affected
        # lineage is the donor, so a driven transfer weights who donates. A run carrying BOTH a driver
        # and a per-family draw multiplies the two — the driver's factor is the lineage's, the
        # multipliers are its contents' — which is what `_driven_weights` does with `fam_sums`.
        w_dup = w_los = w_org = w_tra = None
        fw = None
        if any_driven:  # each lineage's driver values, read before the weights that multiply them in
            drivers = [{**{key: trajs[key].value(alive[k], t) for key in trajs},
                        **{src: live_value(src, k) for src in live_keys}} for k in range(k_alive)]
        if any_family:
            # A per-copy rate pools over copies, so with per-family multipliers the total is the
            # unit rate times the sum of those multipliers over the live copies — and the copy has
            # to be drawn with the same weights, or the rates would say one thing and the picking
            # another. Summed per lineage, so the existing weighted-lineage pick can be reused.
            assert weights is not None       # `any_family` is exactly when it was built
            fw = weights.current(gen)
            unit = {"duplication": dup.effective(copies=1, lineages=1, time=t),
                    "loss": los.effective(copies=1, lineages=1, time=t),
                    "transfer": tra.effective(copies=1, lineages=1, time=t) if can_xfer else 0.0}
            own_sums = (lambda key: fw[key + _FIXED]) if fam_fixed is not None else (lambda key: None)

            def unit_at(key, k, _rates={"duplication": dup, "loss": los, "transfer": tra}):
                """The run's unit rate as lineage ``k`` reads it. Identical to ``unit[key]`` unless
                the rate is driven, and then it is the number `_driven_weights` used for that
                lineage — which the copy pick has to use too, or the totals and the pick disagree
                about how a written family rate compares with the run's."""
                if not any_driven:
                    return unit[key]
                return _rates[key].effective(copies=1, lineages=1, time=t, drivers=drivers[k])
            w_dup = _family_weights(unit["duplication"], fw["duplication"], own_sums("duplication"))
            w_los = _family_weights(unit["loss"], fw["loss"], own_sums("loss"))
            if can_xfer:
                w_tra = _family_weights(unit["transfer"], fw["transfer"], own_sums("transfer"))
        if any_driven:
            if dup_mods:
                w_dup = _driven_weights(dup, gen, k_alive, t, drivers,
                                        fw["duplication"] if fw is not None else None,
                                        own_sums("duplication") if fw is not None else None)
            if los_mods:
                w_los = _driven_weights(los, gen, k_alive, t, drivers,
                                        fw["loss"] if fw is not None else None,
                                        own_sums("loss") if fw is not None else None)
            if org_mods:
                # origination can never carry a per-family draw (refused above: when it is read there
                # is no family yet), so it needs no fam_sums branch
                w_org = [org.effective(copies=len(gen[k]), lineages=1, time=t, drivers=drivers[k])
                         for k in range(k_alive)]
            if tra_mods and can_xfer:
                w_tra = _driven_weights(tra, gen, k_alive, t, drivers,
                                        fw["transfer"] if fw is not None else None,
                                        own_sums("transfer") if fw is not None else None)
        r_dup = sum(w_dup) if w_dup is not None else (
            dup.effective(**(host_ctx if dup_per_lineage else ctx)) if n else 0.0)
        r_los = sum(w_los) if w_los is not None else (
            los.effective(**(host_ctx if los_per_lineage else ctx)) if n else 0.0)
        r_org = sum(w_org) if w_org is not None else org.effective(**ctx)
        r_tra = sum(w_tra) if w_tra is not None else (
            tra.effective(**(host_ctx if tra_per_lineage else ctx)) if can_xfer else 0.0)
        total = r_dup + r_los + r_org + r_tra

        next_species = schedule[si][0]  # the tree's own next event: who is alive changes only here
        # a family placed by `origins=` originates at a fixed instant, so it joins the horizon like
        # any other breakpoint: the waiting time can never step over it
        next_plant = plants[plant_i][0] if plant_i < len(plants) else math.inf
        horizon = min(next_species, next_plant, dup.next_change(t), los.next_change(t),
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
                        j = (_pick_copy_by_family(rng, gen[k], fam_mult["duplication"],
                                                  own_rates("duplication"), unit_at("duplication", k))
                             if any_family else int(rng.integers(len(gen[k]))))
                    elif dup_per_lineage:  # every occupied genome equally likely, then a copy in it
                        k = _pick_host(rng, gen, n_hosts)
                        j = int(rng.integers(len(gen[k])))
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
                        j = (_pick_copy_by_family(rng, gen[k], fam_mult["loss"],
                                                  own_rates("loss"), unit_at("loss", k))
                             if any_family else int(rng.integers(len(gen[k]))))
                    elif los_per_lineage:
                        k = _pick_host(rng, gen, n_hosts)
                        j = int(rng.integers(len(gen[k])))
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
                        if not gen[kd]:    # only via weighted_index's r == total float guard: a
                            # zero-weight lineage has no copies to donate, so take the heaviest instead
                            kd = max(range(k_alive), key=lambda k: w_tra[k])
                        jd = (_pick_copy_by_family(rng, gen[kd], fam_mult["transfer"],
                                                   own_rates("transfer"), unit_at("transfer", kd))
                              if any_family else int(rng.integers(len(gen[kd]))))
                    elif tra_per_lineage:  # every occupied genome donates equally often
                        kd = _pick_host(rng, gen, n_hosts)
                        jd = int(rng.integers(len(gen[kd])))
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
                if node.children:  # a speciation: each gene re-ids into each daughter
                    per_daughter = []
                    for c in node.children:
                        child_genome, rows = [], []
                        for old in g:  # ZOMBI1: the gene ends here and continues under a fresh id
                            nc = new_copy(old.family)
                            child_genome.append(nc)
                            rows.append(GeneEdge(t, "speciation", c, old.family, nc.id, parent=old.id))
                        per_daughter.append(rows)
                        enter(alive, gen, pos, c, child_genome)
                        counts.entered_like(inherited)   # a re-id of the parent: same families
                        if weights is not None:
                            weights.entered(child_genome)
                        total_copies += len(child_genome)
                    # the ids are minted daughter by daughter (which is what fixes them), but a gene's
                    # two rows are recorded together: one gene ending is one event, and the log writes
                    # it as one row with both daughters in it.
                    for pair in zip(*per_daughter):
                        events.extend(pair)
                si += 1
        elif plant_i < len(plants) and horizon == next_plant:
            # a placed family arrives. The lineage is live by construction (its time was checked
            # against that branch's own life), and a tie with the tree's schedule falls to the
            # branch above: the daughters have entered by the time this runs.
            t = horizon
            while plant_i < len(plants) and plants[plant_i][0] == t:
                _, lineage, fam = plants[plant_i]
                k = pos[lineage]
                c = new_copy(fam)
                gen[k].append(c)
                events.append(GeneEdge(t, "origination", lineage, fam, c.id))
                counts.added(k, fam)
                total_copies += 1
                if weights is not None:
                    weights.touched(k)
                plant_i += 1
        else:
            t = horizon  # a skyline breakpoint: advance and re-evaluate the (now changed) rate

    bar.close()
    return FamilyGenomesResult(tree, genomes, events, seed, named, module_map, initial_genome,
                               cap)


# --- process spec: a genome bundled but UNEXECUTED, for a joint model to grow with the tree --------

@dataclass(frozen=True)
class FamilyGenome:
    """A family-genome **process** — its D/T/L/O parameters bundled but not yet run (the genome
    twin of `DiscreteTrait`). ``simulate_genomes_family(tree, ...)`` runs
    this on a *fixed* tree; a **joint** model (``joint.simulate(species.birth_death(...),
    genomes.genome(...))``)
    grows the genome *with* the tree whose speciation its gene content drives. Duplication, loss, and
    origination (each a ``scope(base) × modifiers`` rate, ``changing_at`` allowed) plus ``initial_families``
    and named ``family_names`` (the handle a ``scaled_by("genomes:<name>", …)`` reads). Transfer is not
    available in a joint run: a growing tree's contemporaneous set is still forming as events fire."""

    duplication: object
    loss: object
    origination: object
    initial_families: int
    family_names: tuple
    #: Transfer is available only where the run is handed its tree: a growing tree's contemporaneous
    #: set is still forming as events fire, so there is no "who else is alive now" to draw from.
    transfer: object = 0.0
    max_family_size: "int | None" = 10

    def _resolve(self):
        """Coerce and validate the three rates for the joint engine — ``(duplication, loss,
        origination)`` as resolved `Rate`s. The genome is the **driver**
        here, not the target, so its own rates carry no driver (``changing_at`` is the only verb)."""
        dup = as_rate(self.duplication, default_scope=PerCopy)
        los = as_rate(self.loss, default_scope=PerCopy)
        org = as_rate(self.origination, default_scope=PerLineage)
        for label, rate, want in (("duplication", dup, PerCopy), ("loss", los, PerCopy),
                                  ("origination", org, PerLineage)):
            if rate.scope is not want:
                raise ValueError(
                    f"{label} has a {rate.scope.__name__} scope, but a joint genome takes only "
                    f"{want.__name__} for {label} — drop the scope wrapper."
                )
            for m in rate.modifiers:
                if not matches_declared(m, JOINT_IMPLEMENTED_MODIFIERS):
                    raise ValueError(
                        f"{label} carries {describe(m)}; a joint genome's own rates take only "
                        f"changing_at — the genome is the DRIVER of speciation here, not a driven "
                        f"target."
                    )
        return dup, los, org


#: What a **joint** genome's own rates take (SPEC §5). Declared, like every other level's gate,
#: rather than tested by hand: the genome is the driver of speciation in a joint run, not a driven
#: target, so it takes a schedule and nothing else. `matches_declared` rather than `is_implemented`,
#: because a third-party modifier vouching for itself would still not be threaded by this loop.
JOINT_IMPLEMENTED_MODIFIERS = (OnTime,)


def genome(*, duplication=0.0, loss=0.0, origination=0.0, transfer=0.0,
           initial_families=100, families=None, max_family_size=10,
           **retired) -> FamilyGenome:
    """A **whole-genome** process spec — `FamilyGenome`, unexecuted — for a joint run to simulate
    alongside the tree its gene content drives::

        joint.simulate(species.birth_death(birth=faster_with_toxin, n_extant=100),
                       genomes.genome(origination=0.2, loss=0.1, families=[family("toxin")]))

    Duplication / loss / origination, and ``families=[family("toxin")]`` for the declarations a
    ``scaled_by("genomes:toxin", …)`` reads. A joint run takes no transfer, and a family's own rates
    are not read here — this level's rates apply to the whole genome.

    It is ``genome`` and not ``family`` because it describes a genome. `family()` describes **one gene
    family**, which is what the word means everywhere else in ZOMBI2."""
    check_no_retired_keywords(retired, where="genomes.genome")
    declared = list(families) if families is not None else []
    for spec in declared:
        if not isinstance(spec, GeneFamily):
            raise TypeError(
                f"families takes gene-family declarations — families=[family('toxin')] — got {spec!r}.")
        if spec.written() or spec.origin is not None:
            raise ValueError(
                f"family {spec.name!r} sets rates or an origin, which a joint genome does not read: "
                f"the tree is being simulated with it, so every family runs at this spec's rates. "
                f"Declare it by name alone — family({spec.name!r}).")
    fams = tuple(f.name for f in declared)
    if isinstance(initial_families, bool) or not isinstance(initial_families, int) or initial_families < 0:
        raise ValueError(f"initial_families must be a non-negative integer, got {initial_families!r}")
    if len(set(fams)) != len(fams):
        raise ValueError(f"family names must be unique, got {list(fams)}")
    return FamilyGenome(duplication, loss, origination, initial_families, fams,
                        transfer, resolve_max_family_size(max_family_size))


# --- one named gene family, and the rates it runs at ----------------------------------------------

@dataclass(frozen=True)
class GeneFamily:
    """**One** named gene family: what it is called, where it starts, and the rates it runs at.

    What `family()` builds. A run's rates apply to every family in it, and this is how one family is
    given rates of its own — a mobile element that transfers far more often than the genes around it,
    a core gene that is almost never lost. Anything left ``None`` falls back to the run's rate for
    that event, so a family declared for its name alone behaves exactly as the rest of the genome."""

    name: str
    duplication: object = None
    transfer: object = None
    loss: object = None
    #: ``(lineage, time)`` — where and when this family is planted, instead of at the origin
    origin: object = None
    #: the named group this family belongs to, read back by ``result.completion(...)``
    module: "str | None" = None

    #: the three events a family may set for itself. Origination is not one of them: it is the rate at
    #: which families are *created*, so when it is read this family does not exist to have a rate.
    KEYS: "ClassVar[tuple[str, ...]]" = ("duplication", "transfer", "loss")

    def written(self) -> dict:
        """The rates this family sets for itself, by event name — empty when it sets none."""
        return {k: getattr(self, k) for k in self.KEYS if getattr(self, k) is not None}


def family(name=None, *, duplication=None, transfer=None, loss=None, origin=None,
           module=None) -> GeneFamily:
    """Declare **one gene family** by name, optionally with rates of its own (`GeneFamily`)::

        simulate_genomes_family(tree, initial_families=100, duplication=0.2, loss=0.25, seed=1,
            families=[family("IS1", transfer=PerCopy(1.5), loss=0.02),
                      family("toxin", loss=0.4, origin=("n5", 0.4)),
                      family("nuoA", module="aerobic")])

    The name is the handle everything else reads it by — ``result.presence("IS1")``,
    ``result.has_family(node, "IS1")``. ``duplication`` / ``transfer`` / ``loss`` are that family's
    own rates, and what is left out falls back to the run's. ``origin=(lineage, time)`` plants the
    family there instead of at the origin, and ``module=`` puts it in a named group.

    Origination takes no per-family value: it is the rate at which families are *created*, so when it
    is read this family does not exist yet to have one.
    """
    if name is None or not isinstance(name, str) or not name.strip():
        raise ValueError(
            "family() declares one gene family and needs its name: family('IS1'). For the "
            "whole-genome process spec a joint run takes, the name is genomes.genome(...).")
    if module is not None and (not isinstance(module, str) or not module.strip()):
        raise ValueError(f"module must be a non-empty group name, got {module!r}")
    # a bare lineage means the start of that branch, which is what `origin=(lineage, None)` says;
    # normalising here keeps `resolve_origins` the one place a lineage and a time are checked
    if origin is not None and not isinstance(origin, (tuple, list)):
        origin = (origin, None)
    return GeneFamily(name, duplication, transfer, loss, origin, module)


__all__ = ["simulate_genomes_family", "FamilyGenomesResult", "GeneCopy", "FamilyGenome", "genome",
           "GeneFamily", "family"]
