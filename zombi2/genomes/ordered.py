"""Genomes II — ordered: genes carry a position and an orientation, on chromosomes.

The ordered resolution layers **position** over the family D/T/L/O core (Chapter 4). A genome is
no longer a multiset of gene copies but a list of **chromosomes**, each an ordered run of oriented
`Gene`\\ s.

**Every gene-level event acts on an extent** — a run of consecutive genes (the ZOMBI1 model), its
extent drawn per event from a distribution (default ``Geometric(mean=1)`` — a single gene). The run
starts at a drawn gene and goes rightwards, and where it stops is set by the chromosome's
**topology**: a **circular** chromosome has no ends, so a run that reaches the last gene continues
from the first; a **linear** one has ends, so a run stops at the last gene. Over that
segment: **duplication** copies it in tandem, **loss** removes it, **transfer** sends it to a
contemporaneous recipient as a block, **inversion** reverses it (flipping strands), **transposition**
relocates it elsewhere on the same chromosome, and **translocation** moves it to a different
chromosome; a moved block lands inverted with probability ``inversion_probability``.  (Origination is
the exception — a family is born once, a single gene.)  Inversion/transposition/translocation never
re-mint gene ids: they reshape order and cross genes between chromosome lineages without ending them,
so they live in the ``rearrangements`` log, not the gene genealogy.

Chromosomes carry a genuine **identity** — a chromosome id re-minted at every event that reshapes it —
so ``chromosome_events`` is the true reticulating **chromosome network**: fission (a bifurcation),
fusion (the reticulation), chromosome origination (a de-novo replicon) and chromosome loss, rooted at
the initial and de-novo originations, recorded as an edge list (its ground truth — a network is a graph,
not eNewick).

It is the genome twin of the family core and shares its spine: one forward Gillespie over the
**complete** species tree, the same ``scope(base) × modifiers`` rate grammar, the same gene-genealogy
`GeneEdge` log (position-blind, so ``gene_trees`` and ``profiles`` are
derived from it unchanged), and the same live-lineage bookkeeping. What differs is the state (a list
of chromosomes) and the segmental, position-aware mutators, plus the ``rearrangements`` and
``chromosome_events`` logs. The nucleotide resolution (genes/intergenes, indels) is
`simulate_genomes_nucleotide()`.

**A rate here may be driven by another level** (``scaled_by``, SPEC §2 and §5) — and so may an extent
(SPEC §6). A driven rate is *per lineage*: it is summed over the living lineages, each read with its
own driver value, and the lineage an event lands on is then drawn with those same weights. Because
the gene-level rates are **per copy**, each lineage's weight carries its own gene count, so the pick
is two-stage — a lineage by its weight, then a gene inside it. A driven extent is read at the instant
an event fires, so it changes how much that event takes and never how often one starts.

``transfer_to`` — who receives — is the third place a driver can sit, and the one that is **not** a
rate: there the mapping's numbers are weights normalised across the candidate recipients, so they
redistribute the same transfers rather than change how many happen (SPEC §5, a weight). Its
four rules and the kernel that reads them are the family core's, shared through ``_transfer``. A
declared family can carry its own rule, and a transferred segment goes only where every rule among
its genes allows (SPEC §6).
"""

from __future__ import annotations

import array
import collections
import math
import pathlib
from typing import Any, Sequence, cast
from dataclasses import dataclass, field
from functools import cached_property


from ..params.conditioned import check_mapping_fires, names_a_live_level, resolve_driver
from ..rng import stream
from ..params.parameter import Extent, as_extent
from ..params.mapping import check_not_a_kernel
from ..params.choice import Distance
from ..params.driver import OnTime
from ..params.evaluate import (DRAWN, INHERITED, cell_name, check_one_memory, describe,
                               is_implemented, values_at_birth)
from ..params.connection import Driven, SetBy
from ..params.parameter import Rate, as_rate
from ..params.scope import PerChromosome, PerCopy, PerLineage
from ..tree import Tree, as_tree
from .chromosomes import (CHROMOSOME_EVENTS_HEADER, REARRANGEMENT_HEADER, ChromosomeEvent,
                          chromosome_child_branches, chromosome_event_rows, chromosome_events_tsv,
                          rearrangement_event_rows, rearrangement_events_tsv)
from ..params.retired import check_no_retired_keywords
from .family import (_FamilyCounts, _LiveGeneContent, live_target, resolve_families,
                     resolve_family_rates, resolve_family_transfer_to, resolve_live_drivers,
                     resolve_max_family_size)
from ._live import WeightedIndex, enter, retire, weighted_index, without_cyclic_gc
from ._transfer import (mean_root_to_tip, prepare_transfer_to, recipient_index,
                        recipient_index_all, resolve_transfer_to)
from .._runtime.outputs import fresh_dirs, grouped_dir
from .._runtime.summary import _stats, write_summary
from .._runtime.progress import progress_bar
from .events import (_COLS, Event, EventTally, GeneEdge, _branches, _name, edges_from_tsv,
                     event_counts, event_rows, events_from_edges, gene_label)
from .gene_trees import GeneTree, gene_trees_from_edges, write_gene_trees
from .links import Link, links_of, links_tsv
from .multipliers import (ORDERED_LINEAGE_TARGETS, ORDERED_TARGETS, draw_lineage_multipliers,
                          lineage_multipliers_of, lineage_multipliers_tsv, multipliers_of,
                          multipliers_tsv)
from .profiles import Profiles, profiles_from_genomes, profiles_header, profiles_row
from ._perfamily import StreamedRun

#: The rate grammar this engine supports (SPEC §5) — read by the gate below and by the CLI's help, so
#: a modifier is never advertised without being implemented. The same set the family core takes,
#: because the two are the same model at two resolutions: ``changing_at`` (a skyline in time),
#: ``scaled_by`` (a conditioned or joint driver), ``set_by`` (a driver that replaces the base rather
#: than scaling it), a per-family draw (per-family heterogeneity, weighted on the segment an event
#: covers rather than on the gene it started from — SPEC §6) and a per-lineage draw, drawn or
#: inherited (one multiplier per species branch, on the lineage as a driver's factor is). One
#: combination is refused: see the gate.
IMPLEMENTED_MODIFIERS = (OnTime, Driven, SetBy, (DRAWN, "families"), (DRAWN, "lineages"),
                         (INHERITED, "lineages"))

#: What an **extent** takes here (SPEC §6). An extent takes the modifiers a rate does, minus three.
#: A per-family draw attaches to the *contents*, and an extent is drawn before the run's genes are
#: known — a run covers several families, so there is no one family to draw a factor for; ``set_by``
#: replaces a base, and an extent has none. A per-lineage draw has no such reason and is simply not
#: built here: an extent is sampled on the acting lineage, so the branch is known when it is read,
#: but carrying a drawn factor would mean threading one through `Extent.sample`, and the size an
#: event covers would then want a column of its own in ``lineage_multipliers.tsv``. The two lists
#: are declared separately rather than hidden in an ``if``, because two of the three differences are
#: modelling facts and the third is a piece nobody has built, which is worth reading as such.
IMPLEMENTED_EXTENT_MODIFIERS = (OnTime, Driven)


@dataclass(frozen=True)
class Gene:
    """One gene copy with an **orientation**: a member of family ``family``, identified by a
    globally-unique ``id`` (per segment, the ZOMBI1 model), lying on its chromosome on the ``strand``
    ``+1`` or ``-1``. It is the family `GeneCopy` with the one thing that
    only makes sense once genes are ordered — which way it points. Its position is implicit: the index
    of the gene in its chromosome's ordered list. Birth/death and parentage live in the event log."""

    id: int
    family: int
    strand: int  # +1 / -1


@dataclass
class Chromosome:
    """One chromosome: an ordered run of `Gene`\\ s, identified by ``id`` (re-minted at every
    speciation, so it names a chromosome *lineage*), with a ``topology`` — ``"circular"`` or
    ``"linear"``.

    Topology decides where a segmental event's run stops. A **circular** chromosome has no ends, so a
    run that reaches the last gene continues from the first — it wraps position 0 — and is limited
    only by the whole chromosome. A **linear** one has ends, so a run stops at the last gene. Position
    0 is therefore a real boundary on a linear chromosome and pure bookkeeping on a circular one,
    where it may be re-anchored freely (see `_anchor()`).

    Topology also decides which chromosome events are legal. A **fusion** joins two chromosomes
    of the same topology only, because a ring and a molecule with two ends cannot become one
    molecule (`_fusion()`). A **fission** is legal on either, and gives both halves the parent's
    topology."""

    id: int
    topology: str
    #: a `Sequence`, not a ``list``, because a chromosome has two states: the engine's live one, which
    #: it mutates in place, and the frozen snapshot a finished run hands back, whose genes are a
    #: tuple so a result cannot be edited under the reader. `_live()` narrows back where the engine
    #: does the mutating.
    genes: Sequence[Gene]


# Every rearrangement record names its run the same way: ``start`` is the run's first position in the
# chromosome's frame just *before* the event, and ``length`` is how many genes it covers. The run is
# those positions counted rightwards **modulo the chromosome's gene count**, so ``start + length``
# greater than that count means the run wrapped position 0 — possible only on a circular chromosome.
# Destination fields (``dest``, ``dest_position``) are insertion indices in the frame that exists at
# the moment of insertion, i.e. after the run has been excised.

@dataclass(frozen=True)
class Inversion:
    """A recorded inversion: on species branch ``lineage`` at ``time``, the run of ``length`` genes
    starting at position ``start`` of chromosome ``chromosome`` was reversed and its strands flipped.
    On a circular chromosome the run may wrap position 0 (``start + length`` exceeds the chromosome's
    gene count). Gene ids are untouched — an inversion reshapes order, it does not end lineages — so
    it is logged here, separate from the gene-genealogy `GeneEdge`
    stream."""

    time: float
    lineage: int
    chromosome: int
    start: int
    length: int


@dataclass(frozen=True)
class Transposition:
    """A recorded transposition: on branch ``lineage`` at ``time``, the ``length`` genes starting at
    ``start`` on chromosome ``chromosome`` were excised and reinserted at position ``dest`` on the
    **same** chromosome, ``flipped`` (reversed + strands) or not. The run may wrap position 0 on a
    circular chromosome; ``dest`` indexes what was left after the excision, so it can never fall
    inside the run. Gene ids are untouched."""

    time: float
    lineage: int
    chromosome: int
    start: int
    length: int
    dest: int
    flipped: bool


@dataclass(frozen=True)
class Translocation:
    """A recorded translocation: on branch ``lineage`` at ``time``, the ``length`` genes starting at
    ``start`` on chromosome ``source`` were moved to position ``dest_position`` on chromosome ``dest``
    (a **different** chromosome of the same genome), ``flipped`` or not. The run may wrap position 0
    on a circular ``source``. Gene ids are untouched — a gene lineage crosses to another chromosome
    lineage, which is *not* a chromosome-network edge."""

    time: float
    lineage: int
    source: int
    dest: int
    start: int
    length: int
    dest_position: int
    flipped: bool


@dataclass(frozen=True)
class EventPosition:
    """**Where** one gene-genealogy event happened — the positional companion to an
    `GeneEdge`.

    The event log is position-blind on purpose (it records identity and descent, which is the same
    whatever the resolution), so the ordered engine records position here instead.

    **Every row belongs to exactly one branch.** ``lineage`` names it, and ``chromosome`` /
    ``start`` / ``length`` are coordinates in *that* branch's genome, as it stood just before the
    event. So a reader can take the rows for one branch and know everything that happened to it,
    without holding the rest of the run. One row covers a whole event, even when it acted on a run
    of ``length`` genes — it is not per gene.

    - ``"origination"`` — one new gene of family ``family`` inserted at ``start`` (``length`` 1).
      The only kind that carries a ``family``, because it is the only one whose material does not
      come from a genome the reader already holds.
    - ``"duplication"`` — the run at ``[start, start+length)`` copied in tandem, the copy block
      landing at ``dest_position`` (always ``start+length``; stated so the file needs no outside
      knowledge).
    - ``"loss"`` — the run at ``[start, start+length)`` removed.
    - ``"transfer_donor"`` — the run at ``[start, start+length)`` was copied **out** of this branch.
      The branch itself is unchanged; the row says what left and where it went.
    - ``"transfer_recipient"`` — a block of ``length`` genes arrived **at** ``start`` of this branch.

    A transfer spans two branches, so it is **two records** — one on each — and both name the whole
    edge in ``donor`` and ``recipient``. (This follows Krister Swenson's fork, which splits a transfer
    into a leaving and an arriving event, except that the branches are named outright rather than
    matched by timestamp.) In the *written* table the two are one row, the departing record filling
    ``chromosome`` / ``start`` / ``length`` and the arriving one ``dest_chromosome`` /
    ``dest_position``, so a reader never pairs anything.

    Together with the genomes (``gene_order``) and the rearrangement log this is **sufficient to
    replay a run**. The one record that does not reach the file is the ``loss`` of a copy displaced by
    a replacing transfer: that death is part of the transfer, so it is written as that row's second
    parent and named by copy id rather than by position — a replay tracking ids removes it without
    being told where it sat."""

    time: float
    kind: str  # origination | duplication | loss | transfer_donor | transfer_recipient
    lineage: int  # the species branch these coordinates are in
    chromosome: int
    start: int
    length: int
    family: int | None = None  # origination only: the family the new gene founds
    donor: int | None = None  # both transfer rows: the branch the block was copied out of
    recipient: int | None = None  # both transfer rows: the branch it arrived on
    dest_position: int | None = None  # duplication only: where the tandem copy block lands


@dataclass
class OrderedGenomesResult:
    """What `simulate_genomes_ordered()` returns: the ``complete_tree`` it ran on, the final
    ``genomes`` at **every** node as tuples of `Chromosome`\\ s, the shared gene-genealogy
    ``events`` log, the ``rearrangements`` (inversions, transpositions and translocations) and
    ``chromosome_events`` (the chromosome genealogy) logs, and the ``seed``. The observed genomes are
    the extant tips; ``profiles`` and ``gene_trees`` are derived from the (position-blind) genealogy
    exactly as for the family core; ``gene_order`` reads a node's layout, and ``write`` materialises
    the chosen outputs."""

    complete_tree: Tree
    node_genomes: dict[int, tuple[Chromosome, ...]]
    edges: list[GeneEdge]
    rearrangements: list[Inversion | Transposition | Translocation]
    chromosome_events: list[ChromosomeEvent]
    seed: int | None
    #: ``{name: family id}`` for families declared by ``families=[family(…)]`` — the handle to a *named* family.
    family_names: dict[str, int] = field(default_factory=dict)
    #: ``{module name: (family name, …)}`` for groups declared by ``modules=`` — a pathway or a
    #: complex, whose *completion* in a lineage (`completion`) is a driver. Empty when none were
    #: declared; a module changes nothing about how the genome evolves.
    modules: dict[str, tuple[str, ...]] = field(default_factory=dict)
    #: where each gene-genealogy `GeneEdge` happened — the positional
    #: companion to `events`, which is position-blind. See `EventPosition`.
    event_positions: list[EventPosition] = field(default_factory=list)
    #: The genome the run **started** with, at the root lineage's origination — before any event.
    #: It is not in `genomes`, which holds a genome per *node*, and a node sits at the **end**
    #: of its branch: the root branch is real simulated time, so ``genomes[root]`` is this genome plus
    #: whatever happened along the stem.
    initial_genome: tuple[Chromosome, ...] = ()
    #: The links this run read from its own gene content, one `Link` per modifier reading
    #: ``"genomes:…"`` on a rate, an extent or a family's own rate, which ``write`` puts in
    #: ``links.tsv``. Empty when the run read none.
    links: "tuple[Link, ...]" = ()
    #: Each family's drawn rate multipliers, ``{family: {target: multiplier}}``, as at the family
    #: resolution, with a column for each rearrangement too; ``write`` puts them in
    #: ``family_multipliers.tsv``. Empty when no rate varies among families.
    family_multipliers: "dict[int, dict[str, float | None]]" = field(default_factory=dict)
    #: Each species branch's drawn rate multipliers, ``{node: {target: multiplier}}``, as at the
    #: family resolution, with a column for each rearrangement and each chromosome event too;
    #: ``write`` puts them in ``lineage_multipliers.tsv``. Keyed by node id, as ``node_genomes`` is.
    #: Empty when no rate varies among lineages.
    lineage_multipliers: "dict[int, dict[str, float]]" = field(default_factory=dict)

    def __repr__(self) -> str:
        return (f"OrderedGenomesResult({len(self.complete_tree.extant_leaves())} extant genomes, "
                f"{len(self.node_genomes)} nodes, {len(self.edges)} events, "
                f"{len(self.rearrangements)} rearrangements, seed={self.seed})")

    @property
    def genomes(self) -> dict[str, tuple[Chromosome, ...]]:
        """The observed dataset — the genome at each **extant** tip, keyed by the **tip name** the
        tree writes: ``n5``.

        Keyed by name because the only thing anyone does with this is join it to the tree, or to
        another level grown on that tree, and both name their tips. The written output is keyed by
        name too, so what you get in Python and what you get from a file are the same dataset.
        This is the trait level's `TraitsResult.values` for gene content.

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
        """A multiset view of one node's genome: ``family id → copy count`` (across all chromosomes)."""
        return collections.Counter(g.family for chrom in self.node_genomes[node_id] for g in chrom.genes)

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
        """Whether the named family ``name`` (declared via ``family_names=``) has ≥ 1 copy in the genome at
        ``node_id`` (across all chromosomes)."""
        if name not in self.family_names:
            raise KeyError(f"no named family {name!r}; declared families are {sorted(self.family_names)}")
        fid = self.family_names[name]
        return any(g.family == fid for chrom in self.node_genomes[node_id] for g in chrom.genes)

    def gene_order(self, node_id: int) -> list[tuple[int, int, int, int, int]]:
        """One node's layout as ``(chromosome, position, strand, family, gene id)`` rows, chromosome
        by chromosome and left to right within each — the ordered analogue of ``family_counts``."""
        return [(chrom.id, pos, g.strand, g.family, g.id)
                for chrom in self.node_genomes[node_id] for pos, g in enumerate(chrom.genes)]

    @cached_property
    def _extant_genes(self) -> dict[int, tuple[Gene, ...]]:
        """The observed genomes flattened to gene multisets (chromosomes dropped) — the view the
        genealogy-derived, position-blind outputs read."""
        extant = list(self.complete_tree.extant_leaves())
        return {s: tuple(g for chrom in self.node_genomes[s] for g in chrom.genes) for s in extant}

    @cached_property
    def profiles(self) -> Profiles:
        """The phyletic profiles — each gene family's copy count in each extant species — derived
        from the observed genomes, flattening across chromosomes (position does not enter). See
        `profiles`."""
        return profiles_from_genomes(self._extant_genes, self._extant_genes.keys())

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
        """``{family id: GeneTree}`` — each family's true genealogy inside the complete tree, derived
        from the (position-blind) event log exactly as for the family core. See `gene_trees`."""
        return gene_trees_from_edges(self.edges, self.complete_tree)

    #: Every token ``write()`` honours — the write vocabulary, declared rather than left
    #: implicit in the method body. The CLI builds ``--write``'s choices from this, so the two
    #: cannot drift: they did, and `initial_sequence` and `species_tree` were writable from
    #: Python and unnameable on the command line.
    OUTPUTS = ("events", "profiles", "gene_order", "initial_genome",
               "chromosome_events", "gene_trees", "species_tree", "summary", "links",
               "family_multipliers", "lineage_multipliers")

    def write(self, directory, outputs=("events", "profiles", "gene_order", "initial_genome",
                                        "gene_trees", "chromosome_events", "species_tree",
                                        "summary", "links", "family_multipliers",
                                        "lineage_multipliers"), *,
              flat: bool = False) -> None:
        """Materialise chosen ``outputs`` to ``directory`` (created if needed):

        - ``"events"`` → **two** tables, because a run does two different things to a genome.
          ``genome_events.tsv`` is the gene genealogy — one row per event, in the format every
          resolution writes — with **where** each event happened beside it.
          ``rearrangement_events.tsv`` is the ancestry-neutral rearrangements: an inversion, a
          transposition or a translocation begins and ends no gene lineage, so it has no parents and
          no children and nothing to say in those columns. The two used to be one table, which meant
          nine columns empty on every rearrangement row and six on every genealogy row. Together with
          ``gene_order`` they are enough to replay the run.
        - ``"profiles"`` → ``profiles.tsv``, the family × extant-species copy-count matrix.
        - ``"gene_order"`` → ``gene_order.tsv``, every node's layout (one row per gene), ancestors
          included — so a branch's rearrangements can be replayed from its parent's genome.
        - ``"initial_genome"`` → ``initial_genome.tsv``, the layout the run started with. Its own
          file, not a row in ``gene_order.tsv``, because it belongs to no node: it sits at the start
          of the root branch, and every ``lineage`` in that table is a node at the end of one.
        - ``"chromosome_events"`` → ``chromosome_events.tsv``, the chromosome genealogy edges. The
          one log kept apart: it is a network over chromosome **ids**, with list-valued parents and
          children, joined on a different key from everything above.
        - ``"gene_trees"`` → ``gene_tree_fam<family>_{complete,extant}.nwk`` under ``gene_trees/``,
          each family's true genealogy — unchanged from the family resolution, position being
          orthogonal to it.

        The gene trees are two files per family, so they get a subdirectory rather than burying the
        tables above; ``flat=True`` writes everything into ``directory`` instead.
        - ``"links"`` → ``links.tsv``, the links the run read from its own gene content, one row per
          link (see `zombi2.genomes.links`); the header alone when there are none.
        - ``"family_multipliers"`` → ``family_multipliers.tsv``, each family's drawn rate
          multipliers, one row per family (see `zombi2.genomes.multipliers`); the header alone when no
          rate varies among families.
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
        names = self.complete_tree.labels()   # e<id> for a lineage that died; n<id> for the rest
        if "events" in outputs:
            (d / "genome_events.tsv").write_text(
                _events_tsv(self.edges, self.event_positions, names), encoding="utf-8")
            (d / "rearrangement_events.tsv").write_text(
                rearrangement_events_tsv(self.rearrangements, names), encoding="utf-8")
        if "profiles" in outputs:
            (d / "profiles.tsv").write_text(self.profiles.to_tsv(), encoding="utf-8")
        if "gene_order" in outputs:
            (d / "gene_order.tsv").write_text(self._gene_order_tsv(names), encoding="utf-8")
        if "initial_genome" in outputs:
            (d / "initial_genome.tsv").write_text(self._initial_genome_tsv(), encoding="utf-8")
        if "chromosome_events" in outputs:
            (d / "chromosome_events.tsv").write_text(
                chromosome_events_tsv(self.chromosome_events, self.complete_tree, names),
                encoding="utf-8")
        if "gene_trees" in outputs:
            write_gene_trees(self.gene_trees, grouped_dir(d, "gene_trees", flat),
                             self.complete_tree.labels())
        if "species_tree" in outputs:            # the tree everything here is indexed by: without
            (d / "species_complete.nwk").write_text(   # it a directory of gene trees is not a dataset
                self.complete_tree.to_newick() + "\n", encoding="utf-8")
        if "summary" in outputs:
            write_summary(d / "genome_summary.json", self.summary())
        if "links" in outputs:
            (d / "links.tsv").write_text(links_tsv(self.links), encoding="utf-8")
        if "family_multipliers" in outputs:
            (d / "family_multipliers.tsv").write_text(
                multipliers_tsv(self.family_multipliers, ORDERED_TARGETS), encoding="utf-8")
        if "lineage_multipliers" in outputs:
            (d / "lineage_multipliers.tsv").write_text(
                lineage_multipliers_tsv(self.lineage_multipliers, self.complete_tree.labels(),
                                        ORDERED_LINEAGE_TARGETS), encoding="utf-8")

    def summary(self) -> dict:
        """What this run produced, as a plain dict — the payload of ``genome_summary.json``.

        The **corrected event counts** are the reason this exists. ``genome_events.tsv``'s ``loss``
        rows undercount real losses whenever ``replacement`` is on, because a copy displaced by an
        arriving transfer has no row of its own; the migration guide names that as the change most
        likely to hand a returning user a plausible wrong number, and points them here. This file used
        to be written only at the family resolution — so the advice was sound and the remedy was
        absent at the two resolutions where the gap is *larger* (64% at ordered, measured).

        `event_counts` is shared with the other two resolutions, so the three cannot drift. The rest
        is what this resolution has and the family core does not: where the genes sit, and what moved
        them."""
        t0 = self.complete_tree.nodes[self.complete_tree.root].birth_time
        extant = list(self.complete_tree.extant_leaves())
        return _ordered_summary(
            seed=self.seed, events=event_counts(self.edges, t0),
            born=len({e.family for e in self.edges}),
            surviving=len({g.family for i in extant for c in self.node_genomes.get(i, ())
                           for g in c.genes}),
            named=len(self.family_names), extant=len(extant),
            empty=sum(1 for i in extant if not any(c.genes for c in self.node_genomes.get(i, ()))),
            genes_per_genome=[sum(len(c.genes) for c in self.node_genomes.get(i, ())) for i in extant],
            chromosomes_per_genome=[len(self.node_genomes.get(i, ())) for i in extant],
            rearrangements=collections.Counter(type(r).__name__.lower() for r in self.rearrangements),
            chromosome_events=collections.Counter(e.kind for e in self.chromosome_events))

    def _initial_genome_tsv(self) -> str:
        """The layout the run started with (`_initial_genome_tsv`)."""
        return _initial_genome_tsv(self.initial_genome)

    def _gene_order_tsv(self, names=None) -> str:
        """Every node's gene arrangement, with each chromosome's **topology** beside its id.

        Topology is written here because it is load-bearing and was recoverable from nothing else the
        run wrote: it decides where a segmental event stops (a run wraps position 0 on a ring and
        stops at the end of a linear molecule) and which chromosomes may fuse. It also decides how
        the arrangement is read *out* — the standard rearrangement formats need a per-chromosome
        terminator that differs between a ring and a linear molecule — so a reader handed the output
        directory alone could not produce one. Repeating it on every gene's row is redundant, and the
        alternative was a file of its own for one column.

        A chromosome carrying no genes has no rows here and so no topology, as it has no position or
        strand either: this is the gene arrangement, and an empty replicon has none."""
        rows: list[str] = []
        for s in sorted(self.node_genomes):
            rows.extend(_gene_order_rows(_name(names, s), self.node_genomes[s]))
        return "\n".join(["\t".join(_GENE_ORDER_COLS), *rows]) + "\n"


#: the columns of ``gene_order.tsv``
_GENE_ORDER_COLS = ("lineage", "chromosome", "topology", "position", "strand", "family", "copy")


def _gene_order_rows(name: str, genome) -> list[str]:
    """One node's rows of ``gene_order.tsv``: ``name`` is the node as the files write it, ``genome``
    its chromosomes, read chromosome by chromosome and left to right within each."""
    return [f"{name}\t{chrom.id}\t{chrom.topology}\t{p}\t{g.strand}\t{g.family}\t{gene_label(g.id)}"
            for chrom in genome for p, g in enumerate(chrom.genes)]


def _initial_genome_tsv(initial_genome) -> str:
    """The layout the run started with — ``gene_order.tsv``'s columns without ``lineage``, which is
    the whole point: it belongs to the start of the root branch, not to a node."""
    cols = ("chromosome", "topology", "position", "strand", "family", "copy")
    rows = [f"{chrom.id}\t{chrom.topology}\t{pos}\t{g.strand}\t{g.family}\t{gene_label(g.id)}"
            for chrom in initial_genome for pos, g in enumerate(chrom.genes)]
    return "\n".join(["\t".join(cols), *rows]) + "\n"


def _ordered_summary(*, seed, events, born, surviving, named, extant, empty, genes_per_genome,
                     chromosomes_per_genome, rearrangements, chromosome_events) -> dict:
    """The payload of ``genome_summary.json``, from the numbers it reports. One function, so a run
    kept in memory (`OrderedGenomesResult.summary`) and a streamed run (`_OrderedStream`) cannot report
    them differently. ``rearrangements`` and ``chromosome_events`` count records by kind."""
    return {
        "level": "genomes",
        "seed": seed,
        "resolution": "ordered",
        "events": events,
        "families": {"born": born, "surviving": surviving, "died_out": born - surviving,
                     "named": named},
        "extant_genomes": extant,
        "empty_genomes": empty,
        "genes_per_genome": _stats(genes_per_genome),
        "chromosomes_per_genome": _stats(chromosomes_per_genome),
        # this resolution's own two records: what moved genes without changing their ancestry, and
        # what happened to the replicons carrying them
        "rearrangements": {k: rearrangements.get(k, 0)
                           for k in ("inversion", "transposition", "translocation")},
        "chromosome_events": dict(sorted(chromosome_events.items())),
    }


#: ``genome_events.tsv`` here: the shared genealogy columns (`_COLS` — one row per event, its
#: participants written ``n<species>_g<copy>``) with **where** the event happened beside them. The
#: coordinates are the one thing this resolution has that the family core does not, so they are the
#: one thing it adds; the genealogy half is written by `event_rows()`, not repeated here, because
#: `edges_from_tsv()` reads this table by requiring `_COLS` as a literal **prefix** of the header
#: and spelling them out twice let the two drift.
#:
#: ``chromosome`` / ``start`` / ``length`` are coordinates in the branch's own genome just before the
#: event, as ``gene_order`` numbers it. ``dest_position`` is where the material *landed*: the tandem
#: copy block for a duplication, the arriving block for a transfer — and for a transfer
#: ``dest_chromosome`` names the recipient's chromosome, the recipient branch itself already being
#: inside the arriving copy's token. So one row carries a whole transfer, both ends of it, where the
#: old table spent a row on each side and left a reader to pair them.
#:
#: A segmental event acts on a run of genes of several families, and this table has one ``family``
#: column, so it writes one row per gene lineage — each carrying the **same arc**, the one the event
#: acted on. Rows of one event therefore repeat their coordinates; that is the price of a row being
#: about one gene, and it buys a row that stands alone. (The rearrangements that used to be
#: interleaved here, with the nine genealogy columns blank, are now ``rearrangement_events.tsv``.)
_EVENT_COLS = _COLS + ("chromosome", "start", "length", "dest_chromosome", "dest_position")


def _position_key(kind, lineage, family):
    """What pairs an event with the `EventPosition` recorded for it: its kind, the branch whose
    coordinates those are, and — for an origination alone — the family. The initial originations all
    fire at t=0 on the root branch, and nothing else separates them."""
    return (kind, lineage, family if kind == "origination" else None)


def _coordinates(events, event_positions, branches=None) -> list[str]:
    """The coordinate cells of every written row, in `event_rows()`'s order — the two are zipped into
    one table, so this walks the same `events_from_edges()` the genealogy writer does.

    A transfer needs *both* of its `EventPosition`\\ s: the departing one for the donor's arc and the
    arriving one for where the block landed. Which branch each is on comes from the copies themselves
    (`_branches()`) — the donor's continuation leads ``children``, the arriving copy follows — so a
    self-transfer, where the two branches are the same, still resolves by kind.

    A copy displaced by a replacing transfer has no row of its own (it is that transfer's second
    parent), so its position is not written: it is named by id, and a replay that tracks ids removes
    it without needing to be told where it sat.

    ``branches`` is as in `event_rows`: read off ``events`` when it is not given."""
    where: dict = {}
    for p in event_positions:
        where.setdefault((p.time, *_position_key(p.kind, p.lineage, p.family)), p)
    branch = _branches(events) if branches is None else branches
    out = []
    for time, kind, family, parents, children in (
            (e.time, e.kind, e.family, e.parents, e.children) for e in events_from_edges(events)):
        cells: tuple = ("", "", "", "", "")
        if kind.startswith("transfer"):
            left = where.get((time, "transfer_donor", branch[children[0]], None))
            landed = where.get((time, "transfer_recipient", branch[children[1]], None))
            if left is not None and landed is not None:
                cells = (left.chromosome, left.start, left.length, landed.chromosome, landed.start)
        elif kind != "speciation":               # a speciation copies a genome whole: no arc
            lineage = branch[parents[0] if kind == "loss" else children[0]]
            p = where.get((time, *_position_key(kind, lineage, family)))
            if p is not None:
                cells = (p.chromosome, p.start, p.length, "",
                         "" if p.dest_position is None else p.dest_position)
        out.append("\t".join(str(c) for c in cells))
    return out


def _events_tsv(events, event_positions, names=None) -> str:
    """The genealogy with the place each event happened (see `_EVENT_COLS`)."""
    rows = event_rows(events, names)
    return "\n".join(["\t".join(_EVENT_COLS),
                      *[f"{r}\t{c}" for r, c in zip(rows, _coordinates(events, event_positions))]]
                     ) + "\n"


# --- a run written to disk as it goes (stream_to=) ------------------------------------------------

#: How many records a streamed run holds before it writes them. Writing each step's few records on
#: their own spent more time building rows than the run spent simulating. Any batch of whole steps
#: gives the same rows, the whole log included, and this one is a few megabytes.
_FLUSH_RECORDS = 50_000

#: How many rows of the event log one group of families may hold when a streamed run builds its gene
#: trees at the end: a group is parsed whole, so this is what that step holds in memory. And how many
#: group files are open at once, within a default open-file limit; past that the log is read again
#: for the next files.
_GENE_TREE_GROUP_ROWS = 100_000
_GENE_TREE_OPEN_FILES = 128


class _OrderedStream:
    """Where a streamed ordered run records what it does: its files, written as the run goes.

    The engine appends to the same four lists an in-memory run keeps — `edges`, `positions`,
    `rearrangements` and `chromosome_events` — and `flush` turns what they hold into rows and empties
    them, whenever they hold `_FLUSH_RECORDS` records (`flush_if_full`) and at the end.

    A row names each copy with the branch it lived on, and a copy a row ends was often born many steps
    before. So the stream keeps the branch of every living copy and chromosome, and forgets one when it
    ends: a copy when an event ends it, a chromosome when an edge ends it, and every copy and
    chromosome of a tip when the tip's branch ends. A tip's are forgotten at the next `flush`, after
    its rows are written: some may have been born in the records not yet written, and forgetting them
    before those were written would leave them kept for good.

    `branch_ended` writes a node's ``gene_order.tsv`` rows when its branch ends, so that file lists the
    nodes in the order their branches end. For an extant tip it also keeps what ``profiles.tsv`` and
    the summary need: which families the tip holds, and how many copies of each. `close` writes the
    remaining files and builds the gene trees from the event log (`_write_gene_trees_from_log`)."""

    def __init__(self, directory, outputs, tree) -> None:
        self.directory = pathlib.Path(directory)
        self.directory.mkdir(parents=True, exist_ok=True)
        # the gene trees are one file pair per family, so a previous run's must not survive beside these
        fresh_dirs(self.directory, ("gene_trees",), flat=False)
        self.outputs = tuple(outputs)
        self.tree, self.names = tree, tree.labels()
        self._extant = set(tree.extant_leaves())
        # the buffers the engine appends to, emptied by `flush`
        self.edges: list[GeneEdge] = []
        self.positions: list[EventPosition] = []
        self.rearrangements: list = []
        self.chromosome_events: list[ChromosomeEvent] = []
        self._copy_branch: dict[int, int] = {}
        self._chromosome_branch: dict[int, int] = {}
        self._ended_copies: list[int] = []        # a tip's, forgotten at the next flush
        self._ended_chromosomes: list[int] = []
        self._tally = EventTally(tree.nodes[tree.root].birth_time)
        self._rearranged: collections.Counter = collections.Counter()
        self._chromosome_kinds: collections.Counter = collections.Counter()
        self.n_edges = 0
        self._n_rows = 0
        self._born = 0
        self._tips: dict[int, tuple] = {}       # extant tip -> (families, copies), for profiles.tsv
        self._surviving: set[int] = set()
        self._genes_per_genome: list[int] = []
        self._chromosomes_per_genome: list[int] = []
        self._empty = 0
        want = set(self.outputs)
        # the gene trees are built from the event log, so it is written when either is asked for, and
        # removed at the end when only the gene trees were
        self._events_path = self.directory / ("genome_events.tsv" if "events" in want
                                              else "_genome_events_for_gene_trees.tsv")
        self._events = (self._open(self._events_path, "\t".join(_EVENT_COLS))
                        if want & {"events", "gene_trees"} else None)
        self._rearrangement_file = (self._open(self.directory / "rearrangement_events.tsv",
                                               REARRANGEMENT_HEADER) if "events" in want else None)
        self._chromosome_file = (self._open(self.directory / "chromosome_events.tsv",
                                            CHROMOSOME_EVENTS_HEADER)
                                 if "chromosome_events" in want else None)
        self._gene_order = (self._open(self.directory / "gene_order.tsv", "\t".join(_GENE_ORDER_COLS))
                            if "gene_order" in want else None)

    @staticmethod
    def _open(path, header: str):
        f = open(path, "w", encoding="utf-8")
        f.write(header + "\n")
        return f

    def flush_if_full(self) -> None:
        """`flush`, once the buffers hold `_FLUSH_RECORDS` records."""
        if len(self.edges) + len(self.rearrangements) + len(self.chromosome_events) >= _FLUSH_RECORDS:
            self.flush()

    def flush(self) -> None:
        """Write what the buffers hold, and empty them."""
        edges = self.edges
        if edges:
            births = {e.copy: e.lineage for e in edges}
            branch = collections.ChainMap(births, self._copy_branch)
            if self._events is not None:
                rows = event_rows(edges, self.names, branch)
                cells = _coordinates(edges, self.positions, branch)
                self._events.writelines(f"{row}\t{cell}\n" for row, cell in zip(rows, cells))
                self._n_rows += len(rows)
            self._tally.add(edges)
            self.n_edges += len(edges)
            self._copy_branch.update(births)
            for e in edges:
                if e.kind == "origination":
                    self._born += 1             # every family begins with exactly one origination
                ended = e.parent if e.parent is not None else (e.copy if e.kind == "loss" else None)
                if ended is not None:
                    self._copy_branch.pop(ended, None)
            edges.clear()
        self.positions.clear()
        for copy in self._ended_copies:
            self._copy_branch.pop(copy, None)
        self._ended_copies.clear()
        if self.rearrangements:
            if self._rearrangement_file is not None:
                self._rearrangement_file.writelines(
                    row + "\n" for row in rearrangement_event_rows(self.rearrangements, self.names))
            self._rearranged.update(type(r).__name__.lower() for r in self.rearrangements)
            self.rearrangements.clear()
        if self.chromosome_events:
            born = {child: lineage for chromosome_event in self.chromosome_events
                    for child, lineage in chromosome_child_branches(chromosome_event, self.tree)}
            where = collections.ChainMap(born, self._chromosome_branch)
            if self._chromosome_file is not None:
                self._chromosome_file.writelines(
                    row + "\n" for row in chromosome_event_rows(self.chromosome_events, where,
                                                                 self.names))
            self._chromosome_kinds.update(chromosome_event.kind
                                          for chromosome_event in self.chromosome_events)
            self._chromosome_branch.update(born)
            for chromosome_event in self.chromosome_events:
                for parent in chromosome_event.parents:
                    self._chromosome_branch.pop(parent, None)
            self.chromosome_events.clear()
        for chromosome in self._ended_chromosomes:
            self._chromosome_branch.pop(chromosome, None)
        self._ended_chromosomes.clear()

    def branch_ended(self, node_id: int, genome) -> None:
        """Node ``node_id``'s branch has ended with ``genome``: write its gene order, and for a tip keep
        what the profiles and the summary need and forget its copies and chromosomes, which end here
        without a row. A speciation's own rows end those of an internal node."""
        if self._gene_order is not None:
            self._gene_order.writelines(
                row + "\n" for row in _gene_order_rows(_name(self.names, node_id), genome))
        if self.tree.nodes[node_id].children:
            return
        self._ended_chromosomes.extend(chrom.id for chrom in genome)
        self._ended_copies.extend(g.id for chrom in genome for g in chrom.genes)
        if node_id not in self._extant:
            return
        held = collections.Counter(g.family for chrom in genome for g in chrom.genes)
        self._surviving.update(held)
        self._genes_per_genome.append(sum(held.values()))
        self._chromosomes_per_genome.append(len(genome))
        self._empty += not any(chrom.genes for chrom in genome)
        if "profiles" in self.outputs:
            families = sorted(held)
            self._tips[node_id] = (array.array("q", families),
                                   array.array("q", (held[f] for f in families)))

    def close(self, *, seed, links, initial_genome, named: int, family_multipliers,
              lineage_multipliers) -> StreamedRun:
        """Write what is left, build the gene trees, and hand back the run's `StreamedRun`."""
        self.flush()
        for f in (self._events, self._rearrangement_file, self._chromosome_file, self._gene_order):
            if f is not None:
                f.close()
        d, want = self.directory, set(self.outputs)
        if "gene_trees" in want:
            _write_gene_trees_from_log(self._events_path, self._n_rows, d / "gene_trees", d,
                                       self.tree, self.names)
            if "events" not in want:
                self._events_path.unlink()
        if "profiles" in want:
            self._write_profiles(d / "profiles.tsv")
        if "initial_genome" in want:
            (d / "initial_genome.tsv").write_text(_initial_genome_tsv(initial_genome), encoding="utf-8")
        if "species_tree" in want:
            (d / "species_complete.nwk").write_text(self.tree.to_newick() + "\n", encoding="utf-8")
        if "summary" in want:
            write_summary(d / "genome_summary.json", _ordered_summary(
                seed=seed, events=self._tally.counts(), born=self._born,
                surviving=len(self._surviving), named=named, extant=len(self._extant),
                empty=self._empty, genes_per_genome=self._genes_per_genome,
                chromosomes_per_genome=self._chromosomes_per_genome,
                rearrangements=self._rearranged, chromosome_events=self._chromosome_kinds))
        if "links" in want:
            (d / "links.tsv").write_text(links_tsv(links), encoding="utf-8")
        if "family_multipliers" in want:
            (d / "family_multipliers.tsv").write_text(
                multipliers_tsv(family_multipliers, ORDERED_TARGETS), encoding="utf-8")
        if "lineage_multipliers" in want:
            (d / "lineage_multipliers.tsv").write_text(
                lineage_multipliers_tsv(lineage_multipliers, self.tree.labels(),
                                        ORDERED_LINEAGE_TARGETS), encoding="utf-8")
        return StreamedRun(str(d), seed, self._born, self.n_edges, self.outputs)

    def _write_profiles(self, path) -> None:
        """``profiles.tsv`` from the extant tips' family counts: one row per family present at a tip, in
        family order, and one column per extant species, in node order."""
        species = sorted(self._extant)
        by_family: dict[int, list] = collections.defaultdict(list)
        for column, s in enumerate(species):
            families, copies = self._tips.pop(s, ((), ()))
            for family, n in zip(families, copies):
                by_family[family].append((column, n))
        with open(path, "w", encoding="utf-8") as out:
            out.write(profiles_header(species) + "\n")
            for family in sorted(by_family):
                values = [0] * len(species)
                for column, n in by_family.pop(family):
                    values[column] = n
                out.write(profiles_row(family, values) + "\n")


def _write_gene_trees_from_log(events_path, n_rows: int, directory, scratch, tree, names) -> None:
    """Build every family's gene tree from a written event log and write it into ``directory``,
    holding one group of families in memory at a time.

    Every row of a family goes to the group ``family % groups``, so each group is a whole log for the
    families in it, and building a group's trees gives each family the tree the whole log would. The
    log is read once for every `_GENE_TREE_OPEN_FILES` groups, which are written under ``scratch``,
    built, and removed before the next are written. A group holds about `_GENE_TREE_GROUP_ROWS` rows,
    and never less than its largest family."""
    groups = max(1, -(-n_rows // _GENE_TREE_GROUP_ROWS))
    if groups == 1:
        edges = edges_from_tsv(pathlib.Path(events_path).read_text(encoding="utf-8"))
        write_gene_trees(gene_trees_from_edges(edges, tree), directory, names)
        return
    parts = pathlib.Path(scratch) / "_gene_tree_groups"
    parts.mkdir(exist_ok=True)
    for first in range(0, groups, _GENE_TREE_OPEN_FILES):
        these = range(first, min(first + _GENE_TREE_OPEN_FILES, groups))
        with open(events_path, encoding="utf-8") as log:
            header = log.readline()
            files = {i: open(parts / f"{i}.tsv", "w", encoding="utf-8") for i in these}
            try:
                for f in files.values():
                    f.write(header)
                for line in log:
                    group_file = files.get(int(line.split("\t", 3)[2]) % groups)   # the family column
                    if group_file is not None:
                        group_file.write(line)
            finally:
                for f in files.values():
                    f.close()
        for i in these:
            part = parts / f"{i}.tsv"
            edges = edges_from_tsv(part.read_text(encoding="utf-8"))
            write_gene_trees(gene_trees_from_edges(edges, tree), directory, names)
            part.unlink()
    parts.rmdir()


# --- picking, over the chromosome-nested state ----------------------------------------------------

def _gene_in(genome, m: int) -> tuple[int, int]:
    """The ``m``-th gene of one genome → ``(chromosome index ci, position j)``, counting chromosome
    by chromosome and left to right within each.

    The inner half of a per-copy pick, on its own because two callers need it at different scopes:
    `_pick_gene()` draws ``m`` over the whole live pool, while a **driven** rate has already drawn the
    lineage (weighted by its own rate) and draws ``m`` only inside *that* genome."""
    for ci, chrom in enumerate(genome):
        if m < len(chrom.genes):
            return ci, m
        m -= len(chrom.genes)
    raise AssertionError("the gene index is past the end of the genome")  # unreachable


def _genome_size(genome) -> int:
    """How many genes one genome holds, across all its chromosomes — the ``copies`` a per-copy rate
    is counted per when it is read on a single lineage."""
    return sum(len(c.genes) for c in genome)


def _live_value(target, genome, counts):
    """What a live gene-content driver reads on one lineage right now: its gene count for
    ``"genomes:count"`` (``target`` None), ``"present"`` / ``"absent"`` for a family's id, or the share
    of a module's families (a tuple of ids) that ``counts``, the lineage's copies per family, holds."""
    if target is None:
        return _genome_size(genome)
    if isinstance(target, tuple):
        return sum(1 for f in target if counts[f]) / len(target)
    return "present" if counts[target] else "absent"


#: Tests set this to check the running counts against a fresh count of every genome before each step
#: of the loop and once after it. Off in a real run, where it would cost the walk the counts replace.
_CHECK_COUNTS = False


def _families(genome) -> collections.Counter:
    """How many genes of each family a genome holds, read gene by gene across its chromosomes."""
    return collections.Counter(g.family for chrom in genome for g in chrom.genes)


class _GeneCounts(_FamilyCounts):
    """The family resolution's running counts (`_FamilyCounts`) for ordered genomes: how many genes of
    each family every living lineage holds, kept beside ``gen``.

    The engine changes a row wherever it adds or removes genes: a duplication, a loss, an origination,
    a transfer's arrival and the residents it replaces, and a chromosome loss. At a speciation the
    daughters copy their parent's row. Rearrangements, fissions, fusions and a new empty chromosome
    move no gene in or out, so they leave the rows alone. The family-size cap and every rule that
    reads gene content then look a family up instead of walking a genome.

    With ``recording`` on it also keeps, for each lineage, the families whose copy number changed and
    by how much, until `take_changed` hands them over. Every count change passes through `added` and
    `removed`, so that is where it is written down, and it is how a lineage's row learns what an event
    did to it without reading every family it carries (`_LineageRows`)."""

    def __init__(self, gen) -> None:
        self._counts = [_families(genome) for genome in gen]
        self._changed: list[dict[int, int]] = [{} for _ in gen]
        self.recording = False

    def entered(self, genome) -> None:
        self._counts.append(_families(genome))
        self._changed.append({})

    def entered_like(self, counts) -> None:
        super().entered_like(counts)
        self._changed.append({})

    def retired(self, k: int):
        held = super().retired(k)
        self._changed[k] = self._changed[-1]    # mirror the swap-remove
        self._changed.pop()
        return held

    def of(self, k: int) -> collections.Counter:
        """Lineage ``k``'s genes per family, to read. A family with no gene has no key."""
        return self._counts[k]

    def added(self, k: int, family: int) -> None:
        super().added(k, family)
        if self.recording:
            changed = self._changed[k]
            changed[family] = changed.get(family, 0) + 1

    def removed(self, k: int, family: int) -> None:
        super().removed(k, family)
        if self.recording:
            changed = self._changed[k]
            changed[family] = changed.get(family, 0) - 1

    def added_all(self, k: int, families) -> None:
        for family in families:
            self.added(k, family)

    def removed_all(self, k: int, families) -> None:
        for family in families:
            self.removed(k, family)

    def take_changed(self, k: int) -> dict[int, int]:
        """The families whose copy number in lineage ``k`` changed since this was last asked, each with
        its net change, in the order they first changed — and forget them. A family that went up and
        back down is listed with 0, which costs a look and changes nothing."""
        changed, self._changed[k] = self._changed[k], {}
        return changed


def _check_counts(gen, counts, rows=None) -> None:
    """Raise when a running count differs from a fresh count of its genome (see `_CHECK_COUNTS`) —
    the counts per family, and the gene and chromosome counts a uniform pick reads off `rows`."""
    fresh = [_families(genome) for genome in gen]
    if fresh != counts._counts:
        k = next(k for k in range(max(len(fresh), len(counts._counts)))
                 if k >= len(fresh) or k >= len(counts._counts) or fresh[k] != counts._counts[k])
        raise AssertionError(f"running counts differ from the genomes at lineage index {k}")
    if rows is None:
        return
    for k, genome in enumerate(gen):
        for name, kept, now in (("gene", rows.genes[k], sum(len(c.genes) for c in genome)),
                                ("chromosome", rows.chromosomes[k], len(genome))):
            if kept != now:
                raise AssertionError(
                    f"the kept {name} count differs from the genome at lineage index {k}: "
                    f"{kept!r} against {now!r}")


#: Tests set this to build every lineage's row whole after each step and check it against the row the
#: engine built — whole, or brought up to date from what changed (see `_LineageRows`). Off in a real
#: run, where it would cost the whole builds the rows are there to avoid.
_CHECK_ROWS = False


def _same(kept, fresh) -> bool:
    """Whether a kept row entry agrees with a freshly built one: numbers to a relative 1e-9, a table of
    them entry by entry, anything else exactly."""
    if isinstance(kept, float) and isinstance(fresh, float):
        return kept == fresh or abs(kept - fresh) <= 1e-9 * max(1.0, abs(kept), abs(fresh))
    if isinstance(kept, dict) and isinstance(fresh, dict):
        return kept.keys() == fresh.keys() and all(_same(kept[x], fresh[x]) for x in kept)
    return kept == fresh


def _bounds(kept, fresh) -> bool:
    """Whether a kept largest rate still bounds the fresh one. It may be larger — a rate that left the
    lineage is not taken back out of it — but never smaller."""
    if kept is None or fresh is None:
        return kept is fresh
    return kept >= fresh - 1e-9 * max(1.0, abs(fresh))


class _LineageRows:
    """What a step reads **per living lineage**, kept beside ``gen`` and built one lineage at a time.

    A run whose rates are driven reads every living lineage on its own: its driver values, its rate
    for each event class, its per-family draws summed over its genes, and its summed own rates. Every
    step used to build all of that for every lineage, so a step cost the living lineages times the
    declared families, and a run's cost grew with the square of the tree.

    An event changes the gene content of one lineage, two when a transfer arrives, so no other row
    can differ from the step before: the engine marks the ones it changed with `touched`. A lineage
    that has just entered is built whole. One built before is brought up to date from what changed in
    it (`_GeneCounts.take_changed`): the families whose copy number moved, and the families whose own
    rate reads a driver that moved — not every family it carries. Time is the other thing a row
    reads, and only through schedules, whose every breakpoint is already an instant the loop stops
    at, so crossing one builds every row whole (`touched_all`).

    A row brought up to date holds the numbers a whole build gives, except that its sums were taken by
    subtracting and adding what changed rather than in one pass, so they can differ in the last
    digit; and its largest rate only ever grows, which keeps it a bound on every rate in the row."""

    def __init__(self, w_labels, fam_keys, own_keys) -> None:
        self.drivers: list = []                      # each lineage's driver values, by driver name
        self.own_units: dict[str, list] = {key: [] for key in own_keys}   # the run's rate on its genes
        self.own_owned: dict[str, list] = {key: [] for key in own_keys}   # the own rates it carries
        self.own_largest: dict[str, list] = {key: [] for key in own_keys} # a bound on every gene's rate
        self.own_written: dict[str, list] = {key: [] for key in own_keys} # Σ copies × own rate
        self.own_covered: dict[str, list] = {key: [] for key in own_keys} # the weight those copies take
        # A weight a step totals and an event draws by is kept in a tree rather than a list, so
        # neither costs the living lineages (`WeightedIndex`).
        self.w = {label: WeightedIndex() for label in w_labels}           # its driven rate, per class
        self.fw = {key: WeightedIndex() for key in fam_keys}              # its per-family draws, summed
        self.own_sums = {key: WeightedIndex() for key in own_keys}        # its summed own rates
        # Its gene and chromosome counts. A uniform pick over the whole pool — the plain path's, and
        # every run's for an event class carrying no weight — used to walk the lineages to find the
        # one the draw landed in, which is the same cost per event as adding weights up was per step.
        self.genes = WeightedIndex()
        self.chromosomes = WeightedIndex()
        # every row, with how `_CHECK_ROWS` compares it against a whole build
        self._named = ([("drivers", self.drivers, _same)]
                       + [(f"own unit:{key}", row, _same) for key, row in self.own_units.items()]
                       + [(f"own rates:{key}", row, _same) for key, row in self.own_owned.items()]
                       + [(f"own largest:{key}", row, _bounds) for key, row in self.own_largest.items()]
                       + [(f"own written:{key}", row, _same) for key, row in self.own_written.items()]
                       + [(f"own covered:{key}", row, _same) for key, row in self.own_covered.items()]
                       + [(f"w:{key}", row, _same) for key, row in self.w.items()]
                       + [(f"fw:{key}", row, _same) for key, row in self.fw.items()]
                       + [(f"own sum:{key}", row, _same) for key, row in self.own_sums.items()])
        self._rows = [self.drivers, *self.own_units.values(), *self.own_owned.values(),
                      *self.own_largest.values(), *self.own_written.values(),
                      *self.own_covered.values()]
        self._weights = [*self.w.values(), *self.fw.values(), *self.own_sums.values(),
                         self.genes, self.chromosomes]
        self._stale: set[int] = set()
        self._fresh: set[int] = set()     # entered and not yet built: built whole
        self._all = False

    def __len__(self) -> int:
        return len(self.drivers)

    def entered(self, genome) -> None:
        """A lineage enters the alive set with this genome: a slot on every row, to be built whole at
        the next step, and its counts entered."""
        for row in self._rows:
            row.append(None)
        for weights in self._weights:
            weights.append(0.0)
        k = len(self.drivers) - 1
        self._stale.add(k)
        self._fresh.add(k)
        self._count(k, genome)

    def retired(self, k: int) -> None:
        """Retire lineage ``k``, mirroring `_live.retire`'s swap-remove: the last lineage moves into
        slot ``k``, and brings its marks with it."""
        last = len(self.drivers) - 1
        for row in self._rows:
            row[k] = row[last]
            row.pop()
        for weights in self._weights:
            weights.remove(k)
        for marks in (self._stale, self._fresh):
            moved = last in marks
            marks.discard(last)
            marks.discard(k)
            if moved and k != last:
                marks.add(k)

    def touched(self, k: int, genome) -> None:
        """Lineage ``k``'s genome changed: its row has to be built again, and its counts move with it.

        Every event that adds or removes a gene or a chromosome says so here, and nothing else has
        to, which is what keeps the counts the uniform picks read equal to the genomes themselves."""
        self._stale.add(k)
        self._count(k, genome)

    def _count(self, k: int, genome) -> None:
        self.genes.set(k, float(sum(len(c.genes) for c in genome)))
        self.chromosomes.set(k, float(len(genome)))

    def touched_all(self) -> None:
        """Time moved past a breakpoint, so every row has to be built whole. A flag rather than every
        index, because a run stops at a breakpoint far more often than it reads the rows."""
        self._all = True

    def take_stale(self) -> list[tuple[int, bool]]:
        """The rows to build now, each with whether to build it whole — a lineage that has just
        entered, or every lineage once time has passed a breakpoint — or to bring it up to date from
        what changed in it. They count as current the moment they are handed over, so a caller that
        takes them must build every one."""
        if self._all:
            stale = [(k, True) for k in range(len(self.drivers))]
        else:
            stale = [(k, k in self._fresh) for k in sorted(self._stale)]
        self._all, self._stale, self._fresh = False, set(), set()
        return stale

    def snapshot(self) -> list[list]:
        """Every row as it stands, for `_CHECK_ROWS` to compare a whole build against."""
        return [list(row) for _, row, _ in self._named]

    def check_against(self, before: list[list]) -> None:
        """Raise when a row the engine built differs from the one a whole build just produced — the
        claim the rows rest on: what the engine kept and brought up to date is what building every
        row from scratch gives. Numbers are compared to a relative 1e-9, because a sum brought up to
        date subtracts and adds where a whole build adds once; a largest rate only has to bound the
        fresh one."""
        for (name, row, agrees), kept in zip(self._named, before):
            for k, (fresh, then) in enumerate(zip(row, kept)):
                if not agrees(then, fresh):
                    raise AssertionError(
                        f"the kept {name} row differs from a fresh one at lineage index {k}: "
                        f"{then!r} was kept, {fresh!r} is what the lineage now reads")


class _LiveOrderedContent(_LiveGeneContent):
    """The family resolution's `_LiveGeneContent` for ordered genomes, where a genome is a list of
    chromosomes. A recipient rule reads it when a transfer fires, from the running counts as they are
    then. ``gen``, ``pos`` and ``counts`` are the engine's own, changed in place as the run goes."""

    def value(self, node_id: int, time: float) -> object:
        k = self._pos[node_id]
        return _live_value(self._family, self._gen[k], self._counts.of(k))


#: Tests set this to check a lineage's summed own rates against adding them up family by family.
#: Off in a real run, where it would cost the walk over every family the sum is there to avoid.
_CHECK_OWN_SUMS = False


#: How many sets of driver values one rate remembers before it stops adding more. A family's own rate
#: read on a module meets at most one value more than the module has families, and one read on a
#: presence meets two, so a run stays far below this. The limit is for a driver with many values — a
#: gene count, a continuous trait — where what is remembered would otherwise only grow.
_REMEMBERED = 4096


class _GeneRate:
    """A rate as one gene of a lineage carries it — `Rate.effective` with one copy, one lineage and one
    chromosome — remembered by the driver values it reads.

    A lineage's row holds the rate of every declared family it carries, and every event that changed
    the lineage computed all of them again: with 200 declared families, half of a run. Such a rate
    depends only on the values of the drivers its modifiers name, and an event changes one of those
    at most — one module's completion, one family's presence — so nearly every rate a row asks for
    was computed before from the same values. This returns that number: the result of the same call
    with the same arguments, not an estimate of it.

    A rate is remembered only where that holds: every modifier on it reads a driver, and none of them
    changes with time. Any other rate is computed at every call, as it was."""

    __slots__ = ("rate", "reads", "_remember", "_known")

    def __init__(self, rate: Rate) -> None:
        self.rate = rate
        # the drivers its modifiers name: between two breakpoints, the only things that move its number
        self.reads = tuple(m.key for m in rate.modifiers if isinstance(m, Driven))
        self._remember = (bool(rate.modifiers) and len(self.reads) == len(rate.modifiers)
                          and not math.isfinite(rate.next_change(-math.inf)))
        self._known: dict = {}

    def value(self, time: float, context: dict) -> float:
        """The rate one gene carries now. ``context`` is what the row passes to `Rate.effective`:
        ``{"drivers": the lineage's driver values}``, or empty when the run reads no driver."""
        if not self._remember:
            return self.rate.effective(copies=1, lineages=1, chromosomes=1, time=time, **context)
        drivers = context.get("drivers") or {}
        seen = tuple(drivers.get(key) for key in self.reads)
        try:
            return self._known[seen]
        except KeyError:
            pass
        except TypeError:           # a driver value that cannot be a key: compute, remember nothing
            return self.rate.effective(copies=1, lineages=1, chromosomes=1, time=time, **context)
        value = self.rate.effective(copies=1, lineages=1, chromosomes=1, time=time, **context)
        if len(self._known) < _REMEMBERED:
            self._known[seen] = value
        return value


class _FamilyWeights(dict):
    """Each family's drawn weight for one event class, keyed by family, and the largest drawn so far
    (``largest``) — the bound `_gene_by_rate` keeps or turns down a gene against. Families only
    arrive, so the largest only grows, and it is kept as they are written rather than looked for."""

    def __init__(self) -> None:
        super().__init__()
        self.largest = 0.0

    def __setitem__(self, family: int, weight: float) -> None:
        super().__setitem__(family, weight)
        if weight > self.largest:
            self.largest = weight


class _OwnRates(dict):
    """One lineage's rate per gene, keyed by family: a family's own rate where it writes one, the
    run's rate on every other gene. The run's families are filled in as the genes are read, so a pick
    touches only the families that lineage carries. ``largest`` bounds every rate in it, which is what
    lets a gene be drawn by rejection (`_gene_by_rate`); 0 says nothing is known and the genome is
    scored instead."""

    def __init__(self, own, unit: float, mult, largest: float = 0.0) -> None:
        super().__init__(own)
        self._unit, self._mult = unit, mult
        self.largest = largest

    def __missing__(self, family: int) -> float:
        value = self._unit * (self._mult[family] if self._mult is not None else 1.0)
        self[family] = value
        return value


def _own_table(units, owned, mult, largest):
    """``table(k)`` → lineage ``k``'s rate per gene (an `_OwnRates`), built when an event picks that
    lineage rather than for every lineage at every step."""
    def table(k: int) -> _OwnRates:
        return _OwnRates(owned[k], units[k], mult, largest[k])
    return table


def _check_own_sums(sums, held, table) -> None:
    """Raise when a lineage's summed own rates differ from adding them up family by family, which is
    what the engine does in the long way (see `_CHECK_OWN_SUMS`). The two are the same number
    algebraically, so they are compared to a tolerance rather than exactly."""
    for k, row in enumerate(held):
        rates = table(k)
        direct = sum(row[family] * rates[family] for family in row)
        if abs(direct - sums[k]) > 1e-9 * max(1.0, abs(direct)):
            raise AssertionError(
                f"the own-rate sum differs from the family-by-family sum at lineage index {k}: "
                f"{sums[k]} against {direct}")


def _pick_gene(rng, gen, total_copies, counted) -> tuple[int, int, int]:
    """A uniform global gene pick → ``(lineage k, chromosome index ci in gen[k], position j)``.
    Realises per-copy scope across the whole pool: every gene, in any chromosome of any lineage, is
    equally likely.

    The draw lands on the ``m``-th gene of the pool, and ``counted`` — the living lineages' gene
    counts, as `_LineageRows` keeps them — says which lineage holds it without reading any other
    lineage's chromosomes. Only the chosen genome is then walked, by `_gene_in`. Both halves answer
    what a left-to-right walk of every lineage answered, and the counts are whole numbers, so the
    gene is the gene that walk reached and a run is byte-identical."""
    m = int(rng.integers(total_copies))
    k, within = counted.find(float(m))
    ci, j = _gene_in(gen[k], int(within))
    return k, ci, j


def _pick_chromosome(rng, gen, total_chromosomes, counted, w=None) -> tuple[int, int] | None:
    """A chromosome pick → ``(lineage k, chromosome index ci in gen[k])``, or ``None`` when there is
    nothing to draw.

    Uniform over the whole pool when ``w`` is ``None``, which realises per-chromosome scope: every
    chromosome, in any lineage, is equally likely — found through ``counted``, the living lineages'
    chromosome counts, rather than by walking them (see `_pick_gene`). With ``w`` — the per-lineage totals of a **driven**
    per-chromosome rate — the lineage is drawn by its own weight and the chromosome uniformly inside
    it. That is the same two-stage shape, because a driven lineage's weight already carries its
    chromosome count: ``base × chromosomes_k × factor_k``. Drawing the lineage uniformly instead
    would say one thing in the total and another in the pick."""
    if w is not None:
        if w.total <= 0.0:
            return None                     # every living lineage weighs 0: the event cannot happen
        k = w.pick(rng)
        return (k, int(rng.integers(len(gen[k])))) if gen[k] else None
    k, within = counted.find(float(rng.integers(total_chromosomes)))
    return k, int(within)


# --- extent: every gene-level event acts on a run of consecutive genes (the ZOMBI1 model) ------------

def _extent(rng, ext, chrom, start, ctx=None) -> int:
    """A segment size in genes: sample the event's extent distribution, then clamp it to what
    the chromosome can carry from ``start``.

    A **linear** chromosome has ends, so the run stops at the last gene: ``1 <= m <= n - start``. A
    **circular** one has none, so the run wraps past position 0 and only the whole chromosome bounds
    it: ``1 <= m <= n``. That difference is the point of ``topology``. Clamping a circular run at the
    end of the gene array — as if the array boundary were a real end — would truncate every run that
    started near it, pull the realised mean extent below the nominal one, and leave the genes
    around position 0 covered less often than the rest."""
    m = max(1, int(ext.sample(rng, **(ctx or {}))))
    n = len(chrom.genes)
    return min(m, n) if chrom.topology == "circular" else min(m, n - start)


def _run_means(chrom, mult, m) -> list[float]:
    """For each start, the **mean** family weight of the run of ``m`` genes it opens (SPEC §6).

    Prefix-summed, so this is one pass over the chromosome rather than one per candidate run. A
    circular run wraps position 0; a linear one is clamped by its start, exactly as `_extent()`
    clamps it."""
    w = [mult[g.family] for g in chrom.genes]
    n = len(w)
    pre = [0.0] * (n + 1)
    for i, x in enumerate(w):
        pre[i + 1] = pre[i] + x
    circular = chrom.topology == "circular"
    out = []
    for s in range(n):
        if circular:
            tot = pre[n] - pre[s] + pre[s + m - n] if s + m > n else pre[s + m] - pre[s]
            ln = m
        else:
            ln = min(m, n - s)
            tot = pre[s + ln] - pre[s]
        out.append(tot / ln)
    return out


#: How many genes a rejection draw may try before scoring the genome is the cheaper way: one try per
#: ``_GENES_PER_TRY`` genes, and at least ``_MIN_TRIES``. One try reads one gene and costs about what
#: scoring 5 to 18 genes costs (measured on genomes of 500 and 3000 genes, on one and three
#: chromosomes), so the limit stops before rejection could cost more than the scan it stands in for.
_GENES_PER_TRY = 16
_MIN_TRIES = 8


def _pick_run_by_family(rng, genome, weights, ext, ctx=None) -> tuple[int, int, int] | None:
    """A run drawn with the per-family weight on the **segment**, not on its starting gene (SPEC §6).

    Returns ``(chromosome index, start, run size)``, or ``None`` when the genome has no genes to act
    on. The chance of a start is the **mean** weight of the run it opens: a run of heavily-weighted
    genes is favoured, a mixed one sits in between, an ordinary one is unweighted. Weighting the
    *starting* gene instead would apply a family's own rate to its **neighbours**, and the
    neighbourhood is reshuffled by every rearrangement, so the parameter would not even name a fixed
    thing over a run.

    It is drawn **without scoring every start**. A start's chance is the summed weight of the genes
    its run covers, and a gene is covered by the ``m`` runs of size ``m`` that start on it or up to
    ``m - 1`` genes before it. So drawing one gene in proportion to its weight, and then one of the
    runs covering it uniformly, gives every start exactly its share on a circular chromosome, where
    every gene has ``m`` such runs. On a linear chromosome the genes within ``m`` of an end have fewer,
    and the two draws differ because of them: not at all at the default size of one gene, and by 0.65%
    of the probability for runs of 10 genes on a chromosome of 3000.

    The gene is drawn by rejection (`_gene_by_rate`), which reads a few genes rather than all of them.
    Where the weights are so uneven that it runs out of tries, `_pick_run_by_scan` scores the genome
    instead. Both draws are exact, so which one an event took changes nothing about the process —
    only which random numbers it used.
    """
    size = _genome_size(genome)
    if not size:
        return None
    largest = getattr(weights, "largest", 0.0)
    gene = (_gene_by_rate(rng, genome, size, weights, largest,
                          max(_MIN_TRIES, size // _GENES_PER_TRY)) if largest > 0.0 else None)
    if gene is None:
        return _pick_run_by_scan(rng, genome, weights, ext, ctx)
    ci, j = gene
    chrom = genome[ci]
    n = len(chrom.genes)
    m = min(max(1, int(ext.sample(rng, **(ctx or {})))), n)
    if chrom.topology == "circular":
        return ci, (j - int(rng.integers(m))) % n, m
    first = max(0, j - m + 1)                 # the earliest start whose run still reaches gene j
    start = first + int(rng.integers(j - first + 1))
    return ci, start, min(m, n - start)


def _gene_by_rate(rng, genome, size, weights, largest, tries) -> tuple[int, int] | None:
    """One gene of ``genome`` drawn in proportion to its weight, as ``(chromosome index, position)``,
    or ``None`` when all ``tries`` draws were turned down.

    Rejection: a gene is drawn uniformly and kept with probability its weight over ``largest``, the
    largest weight a gene of this genome can carry; otherwise another is drawn. A kept gene has
    exactly the chance its weight gives it, however many were turned down before it. On average it
    takes ``largest`` over the genome's mean weight draws: one when every gene weighs the same, seven
    when a tenth of the genes weigh 20 times the rest."""
    for _ in range(tries):
        ci, j = _gene_in(genome, int(rng.integers(size)))
        if float(rng.random()) * largest < weights[genome[ci].genes[j].family]:
            return ci, j
    return None


def _pick_run_by_scan(rng, genome, weights, ext, ctx=None) -> tuple[int, int, int] | None:
    """The draw `_pick_run_by_family` makes, made by scoring every start: the chromosome by its summed
    weight, then the size, then the start by the mean weight of its run (`_run_means`). It reads
    every gene, so it is taken only where rejection ran out of tries, or where a table says nothing
    about its largest weight."""
    sums = [sum(weights[g.family] for g in c.genes) for c in genome]
    total = sum(sums)
    if total <= 0.0:
        return None
    ci = weighted_index(rng, sums, total)
    chrom = genome[ci]
    n = len(chrom.genes)
    m = min(max(1, int(ext.sample(rng, **(ctx or {})))), n)
    means = _run_means(chrom, weights, m)
    s = weighted_index(rng, means, sum(means))
    return ci, s, (m if chrom.topology == "circular" else min(m, n - s))


def _run_families(chrom, start, m) -> list[int]:
    """The family of each gene in the run ``[start, start+m)``, in order. The run may wrap position 0
    on a circular chromosome. Read before the event changes the chromosome, so the running counts
    can be moved by exactly what the event took or copied."""
    n = len(chrom.genes)
    return [chrom.genes[(start + i) % n].family for i in range(m)]


def _run_over_cap(held, families, cap) -> bool:
    """Whether copying a run that carries ``families`` would take any of them past ``cap``, where
    ``held`` is the lineage's genes per family (its `_GeneCounts` row).

    The segmental answer to the per-genome family quota. At the family resolution the unit is one
    gene, so the question is simply "is this family already full?"; here a run may carry several
    families, and several copies of one, so the test is *current + carried > cap* for each of them.
    That reduces to exactly the family-resolution condition when the run is a single gene.

    **The whole run is refused, never part of it.** Clipping the run to the genes still under quota
    would be a different process — it would quietly reshape the extent distribution, making runs
    shorter precisely where the genome is crowded. Refusing outright is Poisson thinning on a
    condition that reads only the current state, so what is kept is a clean process; a clipped run
    would not be.

    ``held`` is a count kept as the run goes rather than a walk over the genome, which this used to
    be, on every duplication and every arriving transfer. The count is exact, so the cap binds
    where it bound before."""
    if cap is None:
        return False
    return any(held[f] + k > cap for f, k in collections.Counter(families).items())


def _pick_event_run(rng, gen, n, counted, fw, fam_mult, key, ext, ext_ctx, w=None, hosts=None,
                    own=None):
    """``(lineage, chromosome index, start, run size)`` for one gene-level event, or ``None`` when
    there is nothing to act on.

    Three ways of drawing the same thing, and which one applies is fixed by what the rate carries.
    ``ext_ctx(k)`` builds the context the extent is sampled in, per lineage — it cannot be built
    before the lineage is known, because a driven extent is read on the lineage the event lands on.

    - **driven** (``w`` given, the per-lineage totals of a `Driven` rate) — the lineage is drawn by
      its own rate, then a gene uniformly inside it. Two stages, because a per-copy rate's per-lineage
      total is ``base × copies_k × factor_k``: the gene count is already in the weight, so within the
      lineage every gene is equally likely.
    - **per-family** (``fw`` given) — the lineage by its summed family weights, then the run by
      `_pick_run_by_family()`, so the weight reaches the segment rather than its starting gene.
    - **plain** — one uniform draw over the whole live gene pool.

    A per-family draw and a **driven** rate are refused in one run, because combining them would
    weight by the product of a lineage factor and a segment factor, which is a model neither of them
    is on its own. A per-family draw and a **per-lineage draw** do combine, and then ``w`` and ``fw``
    arrive together: the lineage is drawn by ``w``, which already holds its summed family draws times
    its branch's multiplier, and the segment inside it by those draws.

    - **a family's own rate** (``own`` given, ``(per-lineage sums, a table per lineage)``) — the
      lineage by the summed rates of its genes, then the run by `_pick_run_by_family()` over that
      lineage's rate per gene, so the weight reaches the segment, as a draw's does. The table is asked
      for the drawn lineage alone (see `_own_table`)."""
    if own is not None:
        sums, table = own
        if sums.total <= 0.0:
            return None
        k = sums.pick(rng)
        picked = _pick_run_by_family(rng, gen[k], table(k), ext, ext_ctx(k))
        if picked is None:
            return None
        ci, j, m = picked
        return k, ci, j, m
    if w is not None:
        if w.total <= 0.0:
            return None                     # every living lineage weighs 0: the event cannot happen
        k = w.pick(rng)
        if fw is not None:
            # the lineage's weight already holds its summed family draws, so the segment inside it
            # still has to be drawn by those draws rather than uniformly
            picked = _pick_run_by_family(rng, gen[k], fam_mult[key], ext, ext_ctx(k))
            if picked is None:
                return None
            ci, j, m = picked
            return k, ci, j, m
        size = _genome_size(gen[k])
        if not size:  # only via weighted_index's r == total float guard — a zero-weight lineage has
            return None                     # no gene to act on, so the event is declined (thinning)
        ci, j = _gene_in(gen[k], int(rng.integers(size)))
        return k, ci, j, _extent(rng, ext, gen[k][ci], j, ext_ctx(k))
    if hosts is not None:
        # per-lineage scope: the total counted one share per occupied genome, so the pick is one
        # occupied genome uniformly and then a gene uniformly inside it. A big genome is no more
        # likely to be chosen than a small one — which is the whole difference from per copy.
        if not hosts:
            return None
        k = hosts[int(rng.integers(len(hosts)))]
        ci, j = _gene_in(gen[k], int(rng.integers(_genome_size(gen[k]))))
        return k, ci, j, _extent(rng, ext, gen[k][ci], j, ext_ctx(k))
    if fw is None:
        k, ci, j = _pick_gene(rng, gen, n, counted)
        return k, ci, j, _extent(rng, ext, gen[k][ci], j, ext_ctx(k))
    lw = fw[key]
    if lw.total <= 0.0:
        return None
    k = lw.pick(rng)
    picked = _pick_run_by_family(rng, gen[k], fam_mult[key], ext, ext_ctx(k))
    if picked is None:
        return None
    ci, j, m = picked
    return k, ci, j, m


def _live(chrom: Chromosome) -> list[Gene]:
    """A live chromosome's genes, as the list they are. Only the engine calls this, and only while
    building: a snapshot's genes are a frozen tuple and must not be reached through here."""
    return cast("list[Gene]", chrom.genes)


def _anchor(chrom, start, m) -> int:
    """Make the run ``[start, start+m)`` one contiguous slice, and return the index it now begins at.

    A run that wraps position 0 — only possible on a circular chromosome — is brought to the front by
    rotating the gene list, so it becomes ``[0, m)``; every mutator can then work on a plain slice
    instead of two. Rotating a ring changes nothing biological: on a circular chromosome position 0
    is an index, not a feature of the molecule, so it is free to move (an event that changes the
    run's length, like a segmental duplication or loss, has to move it anyway). A run that does not
    wrap is left where it is and ``start`` comes back unchanged."""
    if start + m <= len(chrom.genes):
        return start
    chrom.genes[:] = chrom.genes[start:] + chrom.genes[:start]
    return 0


def _oriented(segment, flip):
    """The segment as inserted: reversed with each strand flipped if ``flip`` (a moved block that
    landed inverted), else unchanged. Ids are always preserved."""
    return [Gene(g.id, g.family, -g.strand) for g in reversed(segment)] if flip else list(segment)


# --- the mutators (position-, chromosome-, and extent-aware; each records to its log) ----------------

def _originate(genome, node, t, events, positions, new_gene, new_family, rng, family=None) -> int:
    """A new gene family arises: mint a single founding gene (a family is born once — no extent) on
    a uniformly-chosen chromosome at a uniformly-chosen position (strand ``+1``), and record it.
    Returns the family's id, so the caller can count its gene.

    ``family`` names a family whose id was minted earlier — the one case being a family ``origins=``
    places at a chosen point, whose id is fixed before the run walks the tree. The event is the same
    one either way; only where the id came from differs."""
    chrom = genome[int(rng.integers(len(genome)))]
    fam = new_family() if family is None else family
    g = new_gene(fam, +1)
    at = int(rng.integers(len(chrom.genes) + 1))
    _live(chrom).insert(at, g)
    events.append(GeneEdge(t, "origination", node.id, fam, g.id))
    positions.append(EventPosition(t, "origination", node.id, chrom.id, at, 1, family=fam))
    return fam


def _duplicate(chrom, j, m, node, t, events, positions, new_gene) -> int:
    """The ``m`` genes at ``[j, j+m)`` duplicate **in tandem**: each ends and two fresh copies (same
    strand) descend — the continuation in place, the copy block inserted immediately after the
    segment (order preserved). The run may wrap position 0 on a circular chromosome. Returns the
    ``m`` copies added."""
    j = _anchor(chrom, j, m)
    segment = chrom.genes[j:j + m]
    conts = [new_gene(g.family, g.strand) for g in segment]
    copies = [new_gene(g.family, g.strand) for g in segment]
    chrom.genes[j:j + m] = conts + copies              # [.. conts .., .. copies .., ...]
    for old, cont, cp in zip(segment, conts, copies):
        events.append(GeneEdge(t, "duplication", node.id, old.family, cont.id, parent=old.id))
        events.append(GeneEdge(t, "duplication", node.id, old.family, cp.id, parent=old.id))
    positions.append(EventPosition(t, "duplication", node.id, chrom.id, j, m, dest_position=j + m))
    return m


def _lose_at(chrom, j, m, node, t, events, positions) -> int:
    """The ``m`` genes at ``[j, j+m)`` are lost together, removed in place; the run may wrap position
    0 on a circular chromosome. Returns the number removed, which is ``0`` when the loss does not
    happen.

    **A loss never takes a chromosome below its last gene.** A run covering everything still on the
    chromosome does not fire — the same floor the nucleotide resolution enforces in
    `Chromosome.delete()`, so the two resolutions agree on what a chromosome is. Emptying the
    karyotype is `_chromosome_lose()`'s job, an event of its own, counted per chromosome.

    The refusal happens before `_anchor()`, which rotates a wrapping run to the front: a declined
    event must leave the gene order exactly as it found it. It happens after the draw, so the random
    stream is untouched and a run in which no loss is ever declined is byte-identical."""
    if m >= len(chrom.genes):
        return 0
    j = _anchor(chrom, j, m)
    for g in chrom.genes[j:j + m]:
        events.append(GeneEdge(t, "loss", node.id, g.family, g.id))
    del chrom.genes[j:j + m]
    positions.append(EventPosition(t, "loss", node.id, chrom.id, j, m))
    return m


def _invert(chrom, i, m, node, t, rearrangements) -> None:
    """Invert the segment ``[i, i+m)``: reverse the run and flip each gene's strand. On a circular
    chromosome the run may wrap position 0 — reversal on a ring is well defined, and an inversion
    spanning the origin is a real event; a run covering the whole chromosome reverses the
    entire ring, which is the same molecule read the other way round. Ids untouched — identity
    persists through an inversion — so only the rearrangement log is written, and it records the run
    in the frame it had **before** the event."""
    a = _anchor(chrom, i, m)
    chrom.genes[a:a + m] = [Gene(g.id, g.family, -g.strand) for g in reversed(chrom.genes[a:a + m])]
    rearrangements.append(Inversion(t, node.id, chrom.id, i, m))


def _transpose(chrom, i, m, node, t, rearrangements, rng, inversion_probability) -> None:
    """Excise the segment ``[i, i+m)`` and reinsert it elsewhere on the **same** chromosome, flipped
    (reversed + strands) with probability ``inversion_probability``. The run may wrap position 0 on a
    circular chromosome. The destination is drawn *after* the excision, over what is left, so it can
    never land inside the run itself; a run covering the whole chromosome leaves nothing behind, so
    the block goes straight back and only its orientation can change. Ids untouched."""
    a = _anchor(chrom, i, m)
    segment = chrom.genes[a:a + m]
    del chrom.genes[a:a + m]
    flipped = bool(rng.random() < inversion_probability)
    dest = int(rng.integers(len(chrom.genes) + 1))
    chrom.genes[dest:dest] = _oriented(segment, flipped)
    rearrangements.append(Transposition(t, node.id, chrom.id, i, m, dest, flipped))


def _translocate(genome, ci, i, m, node, t, rearrangements, rng, inversion_probability) -> None:
    """Move the segment ``[i, i+m)`` from chromosome ``ci`` to a **different** chromosome of the same
    genome, flipped with probability ``inversion_probability``. No-op if the genome has one
    chromosome. The run may wrap position 0 on a circular source; the destination is on another
    chromosome, so it never falls inside the run, and a run covering the whole source chromosome
    empties it (it survives as an empty replicon). Ids untouched — a gene lineage crosses to another
    chromosome lineage."""
    if len(genome) < 2:
        return
    source = genome[ci]
    a = _anchor(source, i, m)
    segment = source.genes[a:a + m]
    del source.genes[a:a + m]
    flipped = bool(rng.random() < inversion_probability)
    dj = int(rng.integers(len(genome) - 1))
    if dj >= ci:
        dj += 1                                        # a chromosome index distinct from ci
    dest = genome[dj]
    pos = int(rng.integers(len(dest.genes) + 1))
    dest.genes[pos:pos] = _oriented(segment, flipped)
    rearrangements.append(Translocation(t, node.id, source.id, dest.id, i, m, pos, flipped))


def _do_transfer(rng, tree, alive, gen, counts, kd, cdi, jd, m, t, events, positions, new_gene,
                 transfer_to, replacement, self_transfer, depth, cap=None,
                 to_traj=None, groups=None, fam_choice=None) -> tuple[int, int | None]:
    """The segment ``[jd, jd+m)`` on the donor's chromosome ``cdi`` transfers to a contemporaneous
    recipient chosen by ``transfer_to``: each gene ends → a continuation on the donor branch and a
    transferred copy on the recipient (a horizontal gene-tree edge). The run may wrap position 0 on a
    circular donor chromosome. The transferred copies arrive as a block at a random position on a
    uniformly-chosen recipient chromosome (strands travel with them).

    Returns the change in total gene count — ``+m`` additive, minus one per homologous copy displaced
    under ``replacement`` — and the **recipient's index**, which is the one lineage whose gene content
    this changed, so the caller can mark its row (`_LineageRows`). The donor's is not changed: the
    segment it sends is replaced in place by continuations of the same families. A transfer that
    no-ops returns ``(0, None)``.

    **No eligible recipient ⇒ nothing happens.** Under a `Clades` kernel or a driven ``transfer_to``
    a candidate at weight 0 cannot receive, and at some instants that is every candidate. The event is
    then dropped before anything is minted, moved or logged — which is not an approximation: rejecting
    an event on a condition that depends only on the current state is Poisson thinning, so the kept
    transfers are exactly the process whose transfer rate is zero while no recipient is eligible.

    That argument is why the recipient is drawn **first**, above `_anchor()`. Anchoring rotates the
    donor's gene list in place. The rotation is free on a ring biologically, but it renumbers every
    position the run writes out, so a drop after it would leave the donor changed by an event that did
    not happen — the one thing the thinning argument says cannot occur. The pick consumes the rng and
    the anchoring does not, so drawing it first leaves the draw order, and every existing run,
    untouched.

    **A family's own rule.** ``fam_choice`` is ``{family id: (rule, groups, trajectory)}`` for the
    families that carry their own ``transfer_to``. Every gene the segment covers carries its family's
    rule or the run's, and the segment arrives whole, so a lineage can receive it only where every one
    of those rules allows: their weights multiply, each rule once (see `recipient_index_all`). A
    segment whose genes all carry one rule is picked by that rule alone, as before."""
    donor = alive[kd]
    rules = None
    if fam_choice:
        genes, size = gen[kd][cdi].genes, len(gen[kd][cdi].genes)
        run_rule = (transfer_to, groups, to_traj)
        # read before `_anchor`, so a segment that wraps a ring is read across position 0
        covered = (fam_choice.get(genes[(jd + i) % size].family, run_rule) for i in range(m))
        rules = list({id(rule[0]): rule for rule in covered}.values())   # each rule once
        if len(rules) == 1:
            (transfer_to, groups, to_traj), rules = rules[0], None
    if rules is not None:
        cand = [k for k in range(len(alive)) if self_transfer or k != kd]
        if not cand:
            return 0, None
        kr = recipient_index_all(rng, tree, alive, cand, donor, t, rules, depth)
        if kr is None:                                 # no lineage every rule allows — no-op (see above)
            return 0, None
    elif transfer_to == "uniform":
        # O(1) uniform recipient — the same single draw as recipient_index's
        # cand[rng.integers(len(cand))] over every alive lineage but the donor; the donor-skip is a
        # +1 index shift, so no O(alive) candidate list is built per transfer (see family._do_transfer).
        npool = len(alive) if self_transfer else len(alive) - 1
        if npool <= 0:
            return 0, None
        i = int(rng.integers(npool))
        kr = i if (self_transfer or i < kd) else i + 1
    else:  # the weighted rules (Distance / Clades / Driven) weigh every candidate — O(alive)
        cand = [k for k in range(len(alive)) if self_transfer or k != kd]
        if not cand:                                   # the uniform branch's npool guard, restated
            return 0, None
        kr = recipient_index(rng, tree, alive, cand, donor, t, transfer_to, depth, to_traj, groups)
        if kr is None:                                 # every candidate weighs 0 — no-op (see above)
            return 0, None
    recipient = alive[kr]
    rgenome = gen[kr]
    jd = _anchor(gen[kd][cdi], jd, m)
    segment = gen[kd][cdi].genes[jd:jd + m]
    carried = [g.family for g in segment]
    if _run_over_cap(counts.of(kr), carried, cap):      # the recipient is full: same thinning
        return 0, None
    conts = [new_gene(g.family, g.strand) for g in segment]
    xfers = [new_gene(g.family, g.strand) for g in segment]
    gen[kd][cdi].genes[jd:jd + m] = conts               # continuations replace the segment on the donor
    # the donor's row first, then any displacements it causes, then the arrival: within one timestamp
    # the rows are written in the order a replayer must apply them
    positions.append(EventPosition(t, "transfer_donor", donor, gen[kd][cdi].id, jd, m,
                                   donor=donor, recipient=recipient))
    delta = m
    displaced: dict[int, int] = {}                      # arriving copy -> the resident it overwrote
    if replacement:
        cont_ids = {c.id for c in conts}                # self-transfer: never overwrite our own conts
        for x in xfers:                                 # each arriving copy may displace a homolog
            residents = [(ci, p) for ci, ch in enumerate(rgenome) for p, c in enumerate(ch.genes)
                         if c.family == x.family and c.id not in cont_ids]
            if residents:
                ci, p = residents[int(rng.integers(len(residents)))]
                victim = rgenome[ci].genes[p]
                del rgenome[ci].genes[p]
                counts.removed(kr, victim.family)
                displaced[x.id] = victim.id
                positions.append(EventPosition(t, "loss", recipient, rgenome[ci].id, p, 1))
                delta -= 1
    rchrom = rgenome[int(rng.integers(len(rgenome)))]   # arrive as a block on a random recipient chromosome
    pos = int(rng.integers(len(rchrom.genes) + 1))
    rchrom.genes[pos:pos] = xfers
    counts.added_all(kr, carried)                       # the donor keeps its families: continuations
    positions.append(EventPosition(t, "transfer_recipient", recipient, rchrom.id, pos, m,
                                   donor=donor, recipient=recipient))
    # A gene's three edges are recorded together — the resident it displaced dies *of* this transfer,
    # so its `loss` sits with the two transfer edges rather than in a block of its own ahead of them,
    # and `replaced` names it on both. That is what lets the log fold the two into one
    # `transfer_replacing` row and read it back unchanged. The displacements are drawn above, before
    # the arrival is placed, so the random stream is untouched by recording them here.
    for old, cont, xf in zip(segment, conts, xfers):
        replaced = displaced.get(xf.id)
        if replaced is not None:
            events.append(GeneEdge(t, "loss", recipient, old.family, replaced))
        events.append(GeneEdge(t, "transfer", donor, old.family, cont.id, parent=old.id, donor=donor,
                            replaced=replaced))
        events.append(GeneEdge(t, "transfer", recipient, old.family, xf.id, parent=old.id,
                            recipient=recipient, donor=donor, replaced=replaced))
    return delta, kr


# --- the chromosome events: they change chromosome number (the network dynamics) ------------------
# Each re-mints every chromosome id it touches (so no id spans an event) and records one
# ``ChromosomeEvent`` — the edge that makes the genealogy a network. Genes keep their ids: a
# rearrangement moves genes between chromosome lineages, it does not end gene lineages. Each returns
# ``(Δchromosomes, Δgenes)`` for the caller's running totals; ``(0, 0)`` means the event no-op'd.

def _fission(genome, ci, node, t, chromosome_events, new_chromosome, rng) -> tuple[int, int]:
    """Chromosome ``ci`` splits into two at a random cut, both re-minted — a **bifurcation** (one
    parent, two children). No-op on a chromosome of fewer than two genes (nothing to split)."""
    src = genome[ci]
    if len(src.genes) < 2:
        return (0, 0)
    cut = int(rng.integers(1, len(src.genes)))         # 1..len-1: both daughters non-empty
    a = Chromosome(new_chromosome(), src.topology, src.genes[:cut])
    b = Chromosome(new_chromosome(), src.topology, src.genes[cut:])
    genome[ci] = a
    genome.insert(ci + 1, b)
    chromosome_events.append(ChromosomeEvent(t, "fission", node.id, (src.id,), (a.id, b.id)))
    return (1, 0)


def _fusion(genome, ci, node, t, chromosome_events, new_chromosome, rng) -> tuple[int, int]:
    """Chromosome ``ci`` merges with another chromosome **of the same topology** — the
    **reticulation** (two parents, one child): the fused child re-mints, both parents end.

    **A fusion joins two chromosomes of the same topology only.** A ring and a molecule with two
    ends cannot become one molecule, so the partner is drawn from the same-topology chromosomes
    alone — the same rule the nucleotide resolution enforces in `_do_fusion()`, so the two
    resolutions agree on what a chromosome is. Drawing the partner uniformly over the whole
    karyotype and handing the child ``a.topology``, as this did, made a circular chromosome and a
    linear one into one chromosome whose topology was whichever of the two the chromosome pick
    happened to land on first.

    **No same-topology partner ⇒ nothing happens.** The event is dropped before anything is minted
    or logged, which is not an approximation: refusing an event on a condition that reads only the
    current state is Poisson thinning, so what is kept is exactly the process whose fusion rate is
    zero while the chosen chromosome has no legal partner. A genome of one circular and one linear
    chromosome therefore never fuses, however high ``fusion`` is set. (A genome of one chromosome is
    the same case — there is no other chromosome at all — so it needs no separate test.)

    In a genome of a single topology ``partners`` is every other chromosome in index order, so the
    one ``rng.integers`` draw is the one the old arithmetic made and maps to the same chromosome:
    such a run is byte-identical to the run before the rule existed."""
    a = genome[ci]
    partners = [k for k in range(len(genome)) if k != ci and genome[k].topology == a.topology]
    if not partners:
        return (0, 0)
    cj = partners[int(rng.integers(len(partners)))]
    b = genome[cj]
    fused = Chromosome(new_chromosome(), a.topology, a.genes + b.genes)
    genome[:] = [c for idx, c in enumerate(genome) if idx not in (ci, cj)] + [fused]
    chromosome_events.append(ChromosomeEvent(t, "fusion", node.id, (a.id, b.id), (fused.id,)))
    return (-1, 0)


def _chromosome_originate(genome, node, t, chromosome_events, new_chromosome) -> tuple[int, int]:
    """A de-novo replicon (a plasmid) appears: a fresh empty circular chromosome — a **root** of the
    chromosome network (no parent)."""
    new = Chromosome(new_chromosome(), "circular", [])
    genome.append(new)
    chromosome_events.append(ChromosomeEvent(t, "origination", node.id, (), (new.id,)))
    return (1, 0)


def _chromosome_lose(genome, ci, node, t, events, positions, chromosome_events) -> tuple[int, int]:
    """A whole chromosome and its genes die — a **leaf** of the chromosome network (no child); each
    gene on it ends as a gene ``loss``.

    **No-op if it would leave the genome with no genes.** Two cases, and the second is the one that
    is easy to miss: it is refused when it is the genome's last chromosome, and also when it is the
    last chromosome that *has* genes on it. A lineage can be carrying empty replicons — a de-novo
    plasmid from `_chromosome_originate()` starts empty, and `_translocate()` can empty one — so
    "not the last chromosome" is not enough on its own to keep a genome alive. Without the second
    check a lineage holding one gene-bearing chromosome beside an empty one loses everything.

    This is the same floor `_lose_at()` enforces on genes, for the same reason: a chromosome
    without a gene is still a replicon, but a genome without a gene has nothing left for any level
    below to read. Refusing on a condition that reads only the current state is Poisson thinning, so
    what is kept is exactly the process whose ``chromosome_loss`` is zero while the genome is down to
    its last genes; and the refusal happens after the draw, so the random stream is untouched and a
    run that never reaches that state is byte-identical."""
    if len(genome) < 2:
        return (0, 0)
    lost = genome[ci]
    if lost.genes and not any(c.genes for i, c in enumerate(genome) if i != ci):
        return (0, 0)
    for g in lost.genes:
        events.append(GeneEdge(t, "loss", node.id, g.family, g.id))
    if lost.genes:  # the whole chromosome goes, so its genes are one run starting at 0
        positions.append(EventPosition(t, "loss", node.id, lost.id, 0, len(lost.genes)))
    del genome[ci]
    chromosome_events.append(ChromosomeEvent(t, "loss", node.id, (lost.id,), ()))
    return (-1, -len(lost.genes))


# --- initial genome + validation -------------------------------------------------------------------------

def _topologies(chromosomes, topology) -> list[str]:
    """Resolve the ``topology`` argument to one label per initial chromosome."""
    if isinstance(chromosomes, bool) or not isinstance(chromosomes, int) or chromosomes < 1:
        raise ValueError(f"chromosomes must be a positive integer, got {chromosomes!r}")
    if isinstance(topology, str):
        labels = [topology] * chromosomes
    else:
        labels = list(topology)
        if len(labels) != chromosomes:
            raise ValueError(
                f"topology has {len(labels)} entries but chromosomes={chromosomes}; give one label "
                f"per chromosome or a single string for all"
            )
    for label in labels:
        if label not in ("circular", "linear"):
            raise ValueError(f"topology must be 'circular' or 'linear', got {label!r}")
    return labels


# --- the engine -----------------------------------------------------------------------------------

@without_cyclic_gc
def simulate_genomes_ordered(tree, *, duplication=0.0, transfer=0.0, loss=0.0, origination=0.0,
                             inversion=0.0, transposition=0.0, translocation=0.0,
                             chromosomes=1, topology="circular",
                             fission=0.0, fusion=0.0, chromosome_origination=0.0, chromosome_loss=0.0,
                             duplication_extent=None, loss_extent=None, transfer_extent=None,
                             inversion_extent=None, transposition_extent=None,
                             translocation_extent=None, inversion_probability=0.0,
                             transfer_to="uniform", replacement=False, self_transfer=False,
                             initial_families=100, families=None, joint=False,
                             max_family_size=10, seed=None, stream_to=None, outputs=None,
                             progress=False, **retired) -> "OrderedGenomesResult | StreamedRun":
    """Evolve ordered genomes — genes with a position and an orientation, on chromosomes — along a
    species tree, by the D/T/L/O core plus segmental rearrangements and the chromosome events.

    **Every gene-level event acts on an *extent*** — a run of consecutive genes (the ZOMBI1 model):
    ``duplication`` copies the run in tandem, ``loss`` removes it, ``transfer`` sends it to a
    contemporaneous recipient as a block, ``inversion`` reverses it (flipping strands), ``transposition``
    relocates it elsewhere on the same chromosome, and ``translocation`` moves it to a different
    chromosome. The run's **extent** is drawn per event from ``<event>_extent`` (a distribution,
    default ``Geometric(mean=1)`` — usually a single gene; dial the mean up for larger blocks).
    ``origination`` is the exception: a family is born once, a single gene, no extent.
    ``transposition`` and ``translocation`` land the moved block inverted with probability
    ``inversion_probability``.

    **Where a run stops is set by the chromosome's ``topology``.** A run goes rightwards from the gene
    it starts at. On a ``"circular"`` chromosome there are no ends, so a run that reaches the last
    gene continues from the first, and only the whole chromosome bounds it; on a ``"linear"`` one the
    run stops at the last gene. So on a circular chromosome every gene is covered by segmental events
    at the same rate, and the nominal mean extent is the realised one.

    Scopes follow the cross-level grammar, which counts an event per the thing it acts on: the
    gene-level events — ``duplication``/``transfer``/``loss`` and the rearrangements
    ``inversion``/``transposition``/``translocation`` — are **per copy**, since each acts on a run of
    genes that starts at one of them; the chromosome rates ``fission``/``fusion``/``chromosome_loss``
    are **per chromosome**; and the two events that make something from nothing,
    ``origination``/``chromosome_origination``, are **per lineage**. The
    run starts with ``chromosomes`` chromosomes of the given ``topology``, across which the
    ``initial_families`` founding genes are dealt **round-robin**; ``families=[family("toxin")]`` additionally
    declares **named** families (remembered in ``result.family_names`` for ``result.has_family(node,
    "toxin")``), as in the family core; ``replacement`` / ``self_transfer`` behave as in the family
    core. So does ``transfer_to``, which **chooses who receives** — ``"uniform"``,
    ``"distance"`` / ``Distance(decay=)`` (closer relatives likelier), ``Clades({...}, Between({...}))``
    (weight by the donor's and recipient's named clade) or ``Recipients().weighted_by(driver, mapping)`` (weight by
    another level; see below). What moves is a block of genes rather than a single copy, and the block
    arrives whole, so the rule chooses the recipient lineage exactly as it does at the family
    resolution.

    The **chromosome events** change chromosome *number*: ``fission`` (split), ``fusion`` (merge,
    between two chromosomes of the **same topology** — the reticulation; a ring and a molecule with
    two ends cannot become one molecule, so a genome of one of each never fuses),
    ``chromosome_origination`` (a de-novo replicon), ``chromosome_loss`` (a whole
    chromosome and its genes die; never the genome's last). Chromosomes carry identity — re-minted at
    every event that reshapes them — so ``chromosome_events`` is the true reticulating chromosome
    genealogy, rooted at the initial and de-novo originations. Deterministic given ``seed``.

    **A family placed by hand.** ``origins=[("n5", 0.4)]`` originates a family on lineage ``n5`` at
    time ``0.4`` — the ordinary origination event, at a point you choose rather than one that is
    drawn, and here too the founding gene lands on a uniformly-chosen chromosome and position. It
    adds to whatever ``initial_families`` and ``origination`` already give you. See
    `resolve_origins`.

    **Conditioning (a trait drives a rate).** Any rate here may be *driven by another level* —
    ``inversion = PerCopy(0.3).scaled_by(habitat, {"host": 4.0, "free": 1.0})`` scales each lineage's
    inversion rate by the habitat on that branch, read from a trait grown first (the finished
    ``TraitsResult``, or the ``trait_events.tsv`` it wrote). A driven rate is then *per lineage*: it is
    summed over the living lineages, each read with its own gene count, chromosome count and driver
    value; the lineage an event lands on is drawn with those same weights, and the gene inside it
    uniformly, because the gene count is already in the weight. The Gillespie steps at **every**
    mid-branch switch of the driver rather than averaging over a branch (SPEC §2). For ``transfer``
    the driven lineage is the **donor**, so a driven ``transfer`` says how often a lineage *donates*.

    **Conditioning (a trait drives who receives).** ``transfer_to =
    Recipients().weighted_by(driver, mapping)`` is
    the other half, and a different model: the mapping's numbers are per-candidate **weights**, not
    rate multipliers, so they leave the total amount of transfer alone and only redistribute it
    (SPEC §5, a weight, not a rate). Candidate lineage ``k`` gets weight ``mapping(driver value on k now)``
    and receives with probability ``w_k / Σw``. Weight 0 means "cannot receive"; when every candidate
    weighs 0 the transfer does not happen at all, and the donor's chromosome is left untouched. A
    ``Between({...})`` mapping reads the **donor's** value too, so transfer can be steered between
    guilds; ``Clades({...}, Between({...}))`` is the same steering by named clade, read off the tree
    instead of a driver. Because a weight is not a rate, a driven ``transfer_to`` adds no Gillespie
    breakpoint and composes freely with a driven ``transfer`` rate.

    **Conditioning (a trait drives an extent).** An extent takes the same modifiers a rate does
    (SPEC §6) — ``inversion_extent = Extent(4).scaled_by(habitat, {"host": 3.0, "free": 1.0})`` makes a
    host-restricted lineage invert *longer runs of genes*, which is a different statement from raising
    its inversion rate. An extent's modifier is read at the instant an event fires, so it changes how
    much a run takes and never how often one starts, and it adds no Gillespie breakpoint.

    **Joint runs (gene content drives a rate).** With ``joint=True`` a rate, an extent or
    ``transfer_to`` can read the gene content the run is building: ``"genomes:<family>"`` (whether a
    declared family is there), ``"genomes:count"`` (how many genes the lineage has) or
    ``"genomes:module:<group>"`` (the share of a module's families it carries), read on each lineage
    as the run goes. On a rate or an extent the factor belongs to the lineage, so it composes with
    extents as a trait's does (SPEC §6). On ``transfer_to`` each candidate is weighted by its own
    genome when the transfer fires.

    **A family's own rate.** ``families=[family("B", loss=0.8)]`` gives one family its own
    ``duplication``, ``transfer`` or ``loss``, fixed or read from a driver, as at the family
    resolution. It applies to the segment, as a per-family draw does (SPEC §6): every gene carries its
    family's rate or the run's, an event's total is their sum, and a run of genes is chosen by the mean
    rate of the genes it covers. With extents longer than one gene, a family's own rate decides where
    events start, not which genes they take: a family whose rate is 0 can still be removed by a
    segment that starts on a neighbour.

    **A family's own transfer_to.** ``family("B", transfer_to=Recipients().weighted_by("genomes:A",
    {"present": 20.0, "absent": 1.0}))`` gives one family its own recipient rule, in any form the
    run's ``transfer_to`` takes, as at the family resolution. A transferred segment arrives whole, and
    every gene it covers carries its family's rule or the run's. A lineage can receive the segment
    only where every one of those rules allows: their weights multiply, each rule once, before they
    are normalised. With the run's rule uniform, a segment carrying ``B`` keeps ``B``'s preference
    exactly, unless it also carries a family with a rule of its own. A weight of 0 in any rule
    excludes that lineage, and when no lineage is left the transfer does not happen.

    a per-family draw and a driven rate cannot be set in the same run: one weights lineages by a driver and
    the other weights the segment by what it covers.

    **Streaming to disk.** ``stream_to=DIR`` writes the run's files as the run goes and keeps none of
    its record in memory, for a run whose record would not fit. The run is the same run: the same seed
    gives the same events, and the files hold what ``result.write(DIR)`` writes for it. Two things
    differ. ``gene_order.tsv`` lists each node's rows when its branch ends rather than in node order,
    and a `StreamedRun` comes back instead of an `OrderedGenomesResult`. The gene trees are built at the
    end from the event log, one group of families at a time. ``outputs=`` picks the files, as
    ``write()`` takes them, and without ``stream_to`` it is an error. The genomes still alive are kept,
    as the run needs them, so at the end every extant genome is in memory at once.
    """
    tree = as_tree(tree, level="genomes")
    if outputs is not None and stream_to is None:
        raise ValueError(
            "outputs applies to a streamed run (stream_to=DIR), which writes the files itself; for an "
            "in-memory run choose them when you call result.write(outputs=...).")
    if stream_to is not None:
        outputs = tuple(OrderedGenomesResult.OUTPUTS if outputs is None else outputs)
        if unknown := [o for o in outputs if o not in OrderedGenomesResult.OUTPUTS]:
            raise ValueError(f"unknown stream outputs {unknown}; choose from "
                             f"{list(OrderedGenomesResult.OUTPUTS)}")
    labels = _topologies(chromosomes, topology)
    n_initial_chrom = chromosomes
    # this slice implements each event's default scope and the cells IMPLEMENTED_MODIFIERS
    # declares: changing_at (skyline), scaled_by (a conditioned/joint driver, per lineage), set_by (a
    # driver that replaces the base), a per-family draw and a draw among lineages —
    # the per-family one with the weight on the SEGMENT rather than on its starting gene (SPEC §6,
    # and _pick_run_by_family), the per-lineage one on the lineage, as a driver's factor is. The
    # per-lineage draw reaches every rate here, the chromosome events included: it scales whatever
    # the event acts on, where a per-family weight has to reach the genes a segment covers.
    _rates: dict[str, Rate] = {}
    for label, spec, want in (("duplication", duplication, PerCopy), ("transfer", transfer, PerCopy),
                              ("loss", loss, PerCopy), ("origination", origination, PerLineage),
                              ("inversion", inversion, PerCopy),
                              ("transposition", transposition, PerCopy),
                              ("translocation", translocation, PerCopy),
                              ("fission", fission, PerChromosome), ("fusion", fusion, PerChromosome),
                              ("chromosome_origination", chromosome_origination, PerLineage),
                              ("chromosome_loss", chromosome_loss, PerChromosome)):
        rate = as_rate(spec, default_scope=want)
        # An event that acts on **genes** takes either answer to *per what?*: per copy (the default —
        # each gene independently at risk, so a bigger genome turns over faster) or per lineage (a
        # fixed budget, the same however much the genome holds). The chromosome rates and origination
        # keep one scope each: origination creates families, so per copy it would be base × 0 in an
        # empty genome, and the chromosome rates are not implemented per lineage.
        legal = (want, PerLineage) if want is PerCopy else (want,)
        # `rate.scope` holds the scope CLASS, not an instance — a scope constructor returns the rate
        # itself — so this is an identity test against the legal set rather than an isinstance one.
        assert rate.scope is not None            # as_rate fills the level's default where none was written
        if rate.scope not in legal:
            raise ValueError(
                f"{label} has a {rate.scope.__name__} scope, but the ordered genome engine "
                f"takes {' or '.join(s.__name__ for s in legal)} for {label}."
            )
        for m in rate.modifiers:
            if m.reads == (DRAWN, "families") and label == "origination":
                raise ValueError(
                    "origination carries a per-family draw, but origination is the rate at which families are "
                    "CREATED — when it is read there is no family yet to have drawn a factor for. "
                    "Write varying_among('families', …) on duplication, transfer, loss, inversion, "
                    "transposition or translocation; writing one such object on several of them gives "
                    "a family-wide tempo, since one object is one draw.")
            if m.reads == (DRAWN, "families") and rate.scope is not PerCopy:
                raise ValueError(
                    f"{label} carries a per-family draw on a {rate.scope.__name__} scope. A per-family "
                    f"weight has to reach the genes an event covers, so it applies to the per-copy "
                    f"gene events only — not to the chromosome events, which act on whole replicons.")
            if isinstance(m, Driven):
                check_not_a_kernel(m.mapping, label=label)
            if not is_implemented(m, IMPLEMENTED_MODIFIERS, "genomes.ordered"):
                raise ValueError(
                    f"{label} carries {describe(m)}, which the ordered genome engine does not "
                    f"support. It takes changing_at (skyline), scaled_by (a conditioned or joint "
                    f"driver), set_by (a driver that replaces the base), varying_among('families', "
                    f"…) (per-family heterogeneity, weighted on the segment an event covers) and "
                    f"varying_among('lineages', …) (one multiplier per species branch, drawn or "
                    f"inherited)."
                )
        _rates[label] = rate
    # the eleven rates keep short names in the Gillespie loop below; the dict is what the driver
    # resolution and the per-lineage weights walk, so neither has to name all eleven again
    dup, tra, los, org = (_rates["duplication"], _rates["transfer"], _rates["loss"],
                          _rates["origination"])
    inv, trp, trl = _rates["inversion"], _rates["transposition"], _rates["translocation"]
    fis, fus = _rates["fission"], _rates["fusion"]
    cor, clo = _rates["chromosome_origination"], _rates["chromosome_loss"]
    for label, r in _rates.items():
        r.check_one_base(label)
        # SPEC §5: one memory structure per axis. A bare distribution has no memory and a Drift has
        # a continuous one, so a rate carrying both asks for a branch's multiplier to be independent
        # of its parent's and inherited from it at once — there is no model there to implement.
        check_one_memory(tuple(m for m, _ in r.carried_modifiers(unit="lineages")),
                         label=label, unit="lineages")
    # Over the whole RUN, not per rate, and getting that wrong was a real bug: a per-family draw
    # anywhere makes the engine take its per-family path for **every** gene rate, summing each one
    # over the live genes — so a `PerLineage` rate elsewhere in the same run had its total counted
    # per copy while its acting lineage was still drawn uniformly among occupied genomes. The total
    # and the pick then said different things, which is the one failure this engine must not have.
    _GENE_EVENTS = ("duplication", "transfer", "loss", "inversion", "transposition", "translocation")
    per_lineage_here = [lbl for lbl in _GENE_EVENTS if _rates[lbl].scope is PerLineage]
    drawn_here = [lbl for lbl in _GENE_EVENTS
                  if any(m.reads == (DRAWN, "families") for m in _rates[lbl].modifiers)]
    if per_lineage_here and drawn_here:
        raise ValueError(
            f"{', '.join(per_lineage_here)} is PerLineage while {', '.join(drawn_here)} carries a "
            f"per-family draw, and the two cannot share a run. Under PerCopy a family's multiplier "
            f"scales each gene's rate, so it changes the lineage's total; under PerLineage the total "
            f"is fixed whatever the genome holds, so the multiplier could only choose which segment "
            f"is taken. Those are different models and the choice is not made yet — write PerCopy "
            f"throughout for the first, or drop the per-family draw for the second.")
    if any(m.reads == (DRAWN, "families") for r in _rates.values() for m in r.modifiers) and \
            any(isinstance(m, Driven) for r in _rates.values() for m in r.modifiers):
        raise ValueError(
            "a per-family draw and a driver on the same run is not wired at the ordered resolution: "
            "a driver weights the lineage, and here a per-family draw has to weight the SEGMENT an "
            "event covers rather than the gene it started from, so the two are not one "
            "multiplication. The family resolution runs the pair — there a family's multiplier is "
            "the copy's, and the weight is simply the product. Use it, or use one of the two here.")
    # per-event extent distributions (segment size in genes); a bare number is the mean, None a single gene
    def _ext_spec(spec, label):
        """One event's extent (SPEC §6): ``base × modifiers``, no scope, in **genes** here. An extent
        takes the modifiers a rate takes at this resolution minus a per-family draw (see
        `IMPLEMENTED_EXTENT_MODIFIERS`), and they scale the size drawn — ``changing_at`` in time,
        ``scaled_by`` on the lineage the event lands on."""
        e = as_extent(spec)
        rate_slot = label.removesuffix("_extent")
        for m in e.modifiers:
            if isinstance(m, Driven):
                check_not_a_kernel(m.mapping, label=label)
            if m.reads == (DRAWN, "families"):
                raise ValueError(
                    f"{label} carries a per-family draw, which an extent cannot mean: the size is drawn before "
                    f"the run's genes are known, and a run covers several families, so there is no "
                    f"one family to draw a factor for. Put it on {rate_slot}, where it weights "
                    f"the segment by what it covers.")
            if not is_implemented(m, IMPLEMENTED_EXTENT_MODIFIERS, "genomes.ordered"):
                raise ValueError(
                    f"{label} carries {describe(m)}, which the ordered genome engine does not "
                    f"support on an extent — it takes "
                    f"{', '.join(cell_name(w) for w in IMPLEMENTED_EXTENT_MODIFIERS)}.")
        return e

    dup_ext, los_ext, tra_ext = (_ext_spec(duplication_extent, "duplication_extent"),
                                 _ext_spec(loss_extent, "loss_extent"),
                                 _ext_spec(transfer_extent, "transfer_extent"))
    inv_ext, trp_ext, trl_ext = (_ext_spec(inversion_extent, "inversion_extent"),
                                 _ext_spec(transposition_extent, "transposition_extent"),
                                 _ext_spec(translocation_extent, "translocation_extent"))
    _extents: dict[str, Extent] = {
        "duplication_extent": dup_ext, "loss_extent": los_ext, "transfer_extent": tra_ext,
        "inversion_extent": inv_ext, "transposition_extent": trp_ext,
        "translocation_extent": trl_ext}
    if not 0.0 <= inversion_probability <= 1.0:
        raise ValueError(f"inversion_probability must be in [0, 1], got {inversion_probability!r}")
    # the choice (SPEC §5), validated in the one place all three resolutions share: the mapping's
    # numbers are weights over the candidate recipients, never a rate multiplier
    transfer_to = resolve_transfer_to(transfer_to)
    if isinstance(initial_families, bool) or not isinstance(initial_families, int) or initial_families < 0:
        raise ValueError(f"initial_families must be a non-negative integer, got {initial_families!r}")
    check_no_retired_keywords(retired, where="simulate_genomes_ordered")
    # the family resolution's own resolver, so the two engines cannot disagree about what a
    # declaration means
    declared, module_map, planted_named = resolve_families(families, tree)
    family_names = [f.name for f in declared]
    # A family's own duplication, transfer or loss, read by the family resolution's resolver. Here it
    # applies to the segment, as a per-family draw does (SPEC §6): every gene carries its family's rate
    # or the run's, and the extent stays the run's.
    fam_own, fam_driven_rates = resolve_family_rates(
        declared, {"duplication": dup, "transfer": tra, "loss": los})
    own_keys = sorted(set(fam_own) | set(fam_driven_rates))
    any_written = bool(own_keys)
    for key in own_keys:
        run_scope = _rates[key].scope
        assert run_scope is not None             # as_rate fills the level's default where none was written
        if run_scope is not PerCopy:
            raise ValueError(
                f"a family writes its own {key}, but the run's {key} is "
                f"{run_scope.__name__}. The two are summed over the same genes, so both are "
                f"counted per copy: write PerCopy for the run's {key}, or drop the family's.")
    # a family's own recipient rule, checked by the same resolver as the run's transfer_to
    fam_transfer_to = resolve_family_transfer_to(declared)

    # The growth guard, as at the family resolution: duplication compounds, so a run whose rate sits
    # above its loss rate — or a family that drew a high a per-family draw factor — multiplies without bound
    # unless something stops it. A segment may carry several families, and several copies of one, so
    # the run is refused when it would take *any* of them past the quota (see _run_over_cap).
    cap = resolve_max_family_size(max_family_size)

    # Conditioning: a rate written with scaled_by reads a driver **per lineage**, so its rate stops being
    # one number for the whole live set and becomes one per lineage. Same machinery as the other two
    # resolutions — each driver resolves once into a DriverTrajectory keyed by the shared species node
    # id, from a file or an in-memory trait result. With no driven rate and no driven extent this is
    # empty and the loop stays exactly the pooled one, so an undriven run is untouched.
    driven = {label: [m for m in r.modifiers if isinstance(m, Driven)]
              for label, r in _rates.items()}
    ext_driven = {label: [m for m in e.modifiers if isinstance(m, Driven)]
                  for label, e in _extents.items()}
    # A **live** driver names gene content this run is building ("genomes:<family>", "genomes:count",
    # "genomes:module:<group>"), which makes the run joint. It is read off the genomes as the loop goes
    # rather than resolved into a trajectory, and its names are checked by the family resolution's own
    # resolver, so the two engines accept and refuse the same spellings.
    # a family's own rate reads its drivers as a run rate does, live or finished
    fam_rate_driven = [m for table in fam_driven_rates.values() for rate in table.values()
                       for m in rate.modifiers if isinstance(m, Driven)]
    live_rate_mods = [m for mods in (*driven.values(), fam_rate_driven) for m in mods
                      if names_a_live_level(m.driver)]
    live_ext_mods = [m for mods in ext_driven.values() for m in mods if names_a_live_level(m.driver)]
    # A recipient rule reading live gene content is read when a transfer fires. It is checked with the
    # rates' live drivers and makes the run joint too, but it moves no rate (see `prepare_transfer_to`).
    live_choices = [r for r in (transfer_to, *fam_transfer_to.values())
                    if isinstance(r, Driven) and names_a_live_level(r.driver)]
    resolve_live_drivers([*live_rate_mods, *live_ext_mods], set(family_names), joint=joint,
                         choice_mods=live_choices,
                         modules={name: len(members) for name, members in (module_map or {}).items()})
    by_key: dict = {}                   # driver key → its Driven (deduped: one driver resolves once)
    for mods in (*driven.values(), *ext_driven.values(), fam_rate_driven):
        for m in mods:
            if not names_a_live_level(m.driver):
                by_key.setdefault(m.key, m)
    resolved: dict = {}
    if by_key:
        resolved = {key: resolve_driver(m.driver, tree, step=m.step, level="genomes.ordered")
                    for key, m in by_key.items()}
        # a mapping whose states never occur leaves every lineage on the default factor, so the run
        # would secretly be the undriven model — refuse it here, naming the driver
        for mods in (*driven.values(), *ext_driven.values(), fam_rate_driven):
            for m in mods:
                if names_a_live_level(m.driver):
                    continue                    # checked above, and read live rather than resolved
                src = m.driver if isinstance(m.driver, str) else f"<{type(m.driver).__name__}>"
                check_mapping_fires(m.mapping, resolved[m.key].states(), driver_label=src)
    # Only a driver on a **rate** makes the loop per-lineage and adds a Gillespie breakpoint. A driver
    # on an **extent** is read at the instant an event fires — it changes how much that event takes,
    # never how often one happens — so it deliberately stays out of `trajs`: no per-lineage rate
    # weights, no extra horizon steps. (SPEC §6.)
    _rate_keys = {m.key for mods in (*driven.values(), fam_rate_driven) for m in mods}
    trajs = {key: traj for key, traj in resolved.items() if key in _rate_keys}
    # a live driver is a name by construction (names_a_live_level only admits strings)
    live_rate_keys = list(dict.fromkeys(cast(str, m.driver) for m in live_rate_mods))
    live_ext_keys = list(dict.fromkeys(cast(str, m.driver) for m in live_ext_mods))
    any_driven = bool(trajs) or bool(live_rate_keys)
    any_ext_driven = any(ext_driven.values())
    # The transfer_to slot is prepared **after** `trajs` is fixed, for the same reason: a driven
    # transfer_to is a weight, not a rate, so its trajectory must not join `trajs` and start adding
    # horizon breakpoints. `resolved` doubles as the driver cache, so a trait that drives both a rate
    # and who receives is loaded once and read from one trajectory.
    def prepare_choice(rule):
        """``(groups, trajectory)`` for a recipient rule. A rule reading live gene content has nothing
        to prepare here: its reader is built once the family ids and the genomes exist (below)."""
        if isinstance(rule, Driven) and names_a_live_level(rule.driver):
            return None, None
        return prepare_transfer_to(tree, rule, resolved, level="genomes.ordered")

    group_of, to_traj = prepare_choice(transfer_to)
    # A family's rule that equals the run's rule, or an earlier family's, is replaced by that rule, so
    # a segment carrying both counts it once (see `_do_transfer`). A name or a `Distance` is the same
    # rule when it is equal; any other rule only when it is the same object, as elsewhere in the
    # grammar one object read twice is one reading.
    choices: list[tuple] = [(transfer_to, group_of, to_traj)]
    fam_choice_prepared: dict[int, tuple] = {}
    for i, rule in fam_transfer_to.items():
        same = next((c for c in choices if c[0] is rule
                     or (isinstance(rule, (str, Distance)) and type(c[0]) is type(rule) and c[0] == rule)),
                    None)
        if same is None:
            same = (rule, *prepare_choice(rule))
            choices.append(same)
        fam_choice_prepared[i] = same

    rng, seed = stream("genomes", seed)     # own stream, and a drawn seed if none was given

    # Per-lineage multipliers: one per species branch, shared by every family that passes through it
    # (`zombi2.genomes.multipliers`). Drawn before any family is minted, because they come from the
    # tree, which exists before the run does; a run carrying none draws nothing, so it is
    # bit-identical to one from before this existed. On transfer the branch is the DONOR's, as a
    # driven transfer's is.
    lin_by = {label: tuple(m for m, _ in r.carried_modifiers(unit="lineages"))
              for label, r in _rates.items()}
    lin_mult = draw_lineage_multipliers(lin_by, tree, rng, targets=ORDERED_LINEAGE_TARGETS)
    any_lineage = bool(lin_mult)
    varying = {label: bool(mods) for label, mods in lin_by.items()}

    copy_counter = 0
    family_counter = 0
    chrom_counter = 0

    def new_gene(family: int, strand: int) -> Gene:
        nonlocal copy_counter
        g = Gene(copy_counter, family, strand)
        copy_counter += 1
        return g

    # Per-family multipliers, drawn once when a family is created and fixed for its whole life,
    # exactly as at the family resolution: one `Random` object read by two rates is one draw for
    # both, two objects are two draws. What differs here is where the weight lands — on the run an
    # event covers, not on the gene it started from (SPEC §6). Empty unless some rate carries one.
    fam_by = {"duplication": tuple(m for m, _ in dup.carried_modifiers(unit="families")),
              "transfer": tuple(m for m, _ in tra.carried_modifiers(unit="families")),
              "loss": tuple(m for m, _ in los.carried_modifiers(unit="families")),
              "inversion": tuple(m for m, _ in inv.carried_modifiers(unit="families")),
              "transposition": tuple(m for m, _ in trp.carried_modifiers(unit="families")),
              "translocation": tuple(m for m, _ in trl.carried_modifiers(unit="families"))}
    any_family = any(fam_by.values())
    fam_mult: dict[str, _FamilyWeights] = {key: _FamilyWeights() for key in fam_by}

    # Which gene-level rates are a fixed per-lineage budget rather than a per-gene risk. Read once:
    # it decides both how the total is counted and how the acting lineage is picked, and those two
    # must never disagree.
    per_lineage = {label: _rates[label].scope is PerLineage
                   for label in ("duplication", "transfer", "loss",
                                 "inversion", "transposition", "translocation")}
    any_per_lineage = any(per_lineage.values())

    def new_family() -> int:
        nonlocal family_counter
        f = family_counter
        family_counter += 1
        if any_family:
            # one draw per distinct modifier object for this family, shared across its rates (see
            # `values_at_birth`): one object written on two rates is one number.
            shared: dict[int, float] = {}
            for key, mods in fam_by.items():
                fam_mult[key][f] = math.prod(values_at_birth(mods, rng, shared))
        return f

    def new_chromosome() -> int:
        nonlocal chrom_counter
        cid = chrom_counter
        chrom_counter += 1
        return cid

    depth = mean_root_to_tip(tree)  # timescale for Distance weighting (unused by "uniform")
    schedule = sorted((tree.nodes[i].end_time, i) for i in tree.nodes)  # (end_time, node_id)

    root = tree.nodes[tree.root]
    t = root.birth_time
    alive: list[int] = []
    gen: list[list[Chromosome]] = []
    pos: dict[int, int] = {}
    genomes: dict[int, tuple[Chromosome, ...]] = {}
    # A streamed run writes these four as it goes and keeps none of them (`_OrderedStream`): the engine
    # appends to the stream's own lists, and the stream writes and empties them in batches.
    to_disk = _OrderedStream(stream_to, outputs, tree) if stream_to is not None else None
    events: list[GeneEdge] = to_disk.edges if to_disk is not None else []
    event_positions: list[EventPosition] = to_disk.positions if to_disk is not None else []
    rearrangements: list[Inversion | Transposition | Translocation] = (
        to_disk.rearrangements if to_disk is not None else [])
    chromosome_events: list[ChromosomeEvent] = (
        to_disk.chromosome_events if to_disk is not None else [])

    initial_chroms = []
    for label in labels:  # lay down the initial karyotype; each initial chromosome is a network root
        cid = new_chromosome()
        initial_chroms.append(Chromosome(cid, label, []))
        # `initial`, not `origination`: a replicon the run *starts* with is not something it did, so
        # counting `origination` in the log gives the de-novo replicons alone
        chromosome_events.append(ChromosomeEvent(t, "initial", root.id, (), (cid,)))
    # the origin's initial genome is logged like any other origination — each founding gene appended in turn —
    # so the position table is total over gene-content events and a replay of the root branch can
    # start from an empty karyotype (every other branch starts from its parent's gene_order rows)
    for i in range(initial_families):  # deal the founding genes round-robin across the chromosomes
        fam = new_family()
        chrom = initial_chroms[i % n_initial_chrom]
        _live(chrom).append(new_gene(fam, +1))
        events.append(GeneEdge(t, "origination", root.id, fam, chrom.genes[-1].id))
        event_positions.append(EventPosition(t, "origination", root.id, chrom.id,
                                             len(chrom.genes) - 1, 1, family=fam))
    # a family's own rate by the id the family is given: fixed ones as numbers, driven ones as the Rate
    fam_fixed_by_id: dict[str, dict[int, float]] = {key: {} for key in own_keys}
    fam_driven_by_id: dict[str, dict[int, _GeneRate]] = {key: {} for key in own_keys}
    named: dict[str, int] = {}  # a minted id per declared name, dealt round-robin after the anonymous ones
    named_plants: list[tuple[float, int, int]] = []
    for j, name in enumerate(family_names):
        fam = new_family()
        named[name] = fam
        for key in own_keys:
            if j in fam_own.get(key, {}):
                fam_fixed_by_id[key][fam] = fam_own[key][j]
            elif j in fam_driven_rates.get(key, {}):
                fam_driven_by_id[key][fam] = _GeneRate(fam_driven_rates[key][j])
        if j in planted_named:
            # given an `origin`, so it arrives there rather than at the tree's origin — the same
            # event, at a point chosen instead of drawn
            t_p, lineage = planted_named[j]
            named_plants.append((t_p, lineage, fam))
            continue
        chrom = initial_chroms[(initial_families + j) % n_initial_chrom]
        _live(chrom).append(new_gene(fam, +1))
        events.append(GeneEdge(t, "origination", root.id, fam, chrom.genes[-1].id))
        event_positions.append(EventPosition(t, "origination", root.id, chrom.id,
                                             len(chrom.genes) - 1, 1, family=fam))
    # the run's starting genome: a deep snapshot, so the live genome's events never reach it
    # the ids of the families `origins=` places: minted here, straight after the initial and named
    # ones and in the order they were written, so the same origins name the same families at either
    # resolution. Each is planted at its own time, in the loop below.
    plants = sorted(named_plants)
    plant_i = 0
    initial_genome = tuple(Chromosome(c.id, c.topology, list(c.genes)) for c in initial_chroms)
    enter(alive, gen, pos, root.id, initial_chroms)
    counts = _GeneCounts(gen)       # genes per family on every living lineage, changed with the genomes
    # What a step reads per lineage, kept across steps so a step rebuilds only the lineages an event
    # changed (`_LineageRows`). A run with no driven rate, no draw of either kind and no own rate
    # reads nothing per lineage, and then the rows are there only to count the living lineages.
    # A rate is read lineage by lineage when something makes it differ between lineages, and two
    # things do: a driver, whose value is the one on that branch, and a per-lineage draw, whose
    # multiplier was drawn for that branch. Either way the rate gets a `w` row, which the step totals
    # and an event draws the acting lineage by.
    lineage_rates = {label: rate for label, rate in _rates.items()
                     if driven[label] or varying[label]}
    # Where each `w` row comes from. A rate carrying a per-family draw as well keeps its weight in
    # `fw` — the family multipliers summed over its genes, which the segment pick needs unscaled —
    # so its row is that sum times the branch's multiplier, built beside `fw` rather than from the
    # rate alone. A driver and a per-family draw are refused in one run, so these two never overlap.
    w_from_family = {label for label in lineage_rates if any_family and fam_by.get(label)}
    w_from_rate = {label for label in lineage_rates if label not in w_from_family}
    any_rows = any_driven or any_family or any_written or any_lineage
    rows = _LineageRows(lineage_rates if (any_driven or any_lineage) else (),
                        fam_mult if any_family else (), own_keys if any_written else ())
    for genome in gen:
        rows.entered(genome)
    # Whether anything in this run changes with time on its own. A row reads time only through a
    # schedule — on a rate, on a family's own rate, or inside a driver's mapping — and through a
    # driver grown before the run. With none of those, the tree's own schedule moves `t` without
    # moving a single number a row holds, so a speciation marks the daughters that enter and leaves
    # every other lineage's row alone. That is the difference between a step at a speciation costing
    # the living lineages and costing two.
    # The rates that change on their own at some instant: a schedule on the rate, or on an entry of a
    # driver's mapping. `next_change` names the first change strictly after an instant, so a rate with
    # none after -inf has none at all. The loop sets its horizon from these alone, instead of asking
    # every rate at every step — with 200 declared families that was 211 questions a step, each
    # answered "never".
    timed_rates = [r for r in (dup, los, org, tra, inv, trp, trl, fis, fus, cor, clo,
                               *(fr for table in fam_driven_rates.values() for fr in table.values()))
                   if math.isfinite(r.next_change(-math.inf))]
    time_varying = bool(trajs) or bool(timed_rates)
    # one entry per event class whose families write their own rate, reading the rows as they stand:
    # the lineage's summed rates to draw it by, and `table(k)` for the rate of each of its genes
    own_pick = {key: (rows.own_sums[key],
                      _own_table(rows.own_units[key], rows.own_owned[key],
                                 fam_mult[key] if any_family else None, rows.own_largest[key]))
                for key in (own_keys if any_written else ())}
    # the run's own rate for each of those event classes, as the genes without an own rate carry it
    run_gene_rates = {key: _GeneRate(_rates[key]) for key in (own_keys if any_written else ())}
    # a family given an `origin` is not in the root genome — it arrives later, in the loop
    total_copies = initial_families + len(family_names) - len(named_plants)
    total_chromosomes = n_initial_chrom
    # each live driver paired with what it reads, resolved once: the names and the family ids are fixed
    live_rate_reads = [(src, live_target(src, named, module_map or {})) for src in live_rate_keys]
    live_ext_reads = [(src, live_target(src, named, module_map or {})) for src in live_ext_keys]

    # Recipient rules that read live gene content get their reader now, when the family ids and the
    # genomes exist. A family's own rule is keyed by the id the family was given, which is what each
    # of its genes carries.
    def with_reader(choice):
        rule, groups_c, traj_c = choice
        if isinstance(rule, Driven) and names_a_live_level(rule.driver):
            traj_c = _LiveOrderedContent(live_target(rule.driver, named, module_map or {}), gen, pos,
                                         counts)
        return rule, groups_c, traj_c

    transfer_to, group_of, to_traj = with_reader((transfer_to, group_of, to_traj))
    fam_choice = {named[declared[i].name]: with_reader(c) for i, c in fam_choice_prepared.items()}

    # eleven bare numbers on their default scopes — no modifier on any rate, none on any extent,
    # no per-lineage budget — is the common run, and it needs none of the loop's context machinery:
    # each total is scope(base) exactly, resolved here once. A rate whose base is None carries a
    # set_by, which is a modifier, so `plain` is False and its 0.0 is never read.
    plain = (not any_written and not any(r.modifiers for r in _rates.values()) and not any_per_lineage
             and not any(e.has_modifiers for e in (dup_ext, los_ext, tra_ext, inv_ext, trp_ext,
                                                   trl_ext)))
    dup_base, los_base, tra_base = dup.base or 0.0, los.base or 0.0, tra.base or 0.0
    org_base, inv_base, trp_base = org.base or 0.0, inv.base or 0.0, trp.base or 0.0
    trl_base, fis_base, fus_base = trl.base or 0.0, fis.base or 0.0, fus.base or 0.0
    cor_base, clo_base = cor.base or 0.0, clo.base or 0.0
    no_weights: dict = {}    # what `w` is when nothing is driven: read by .get, never written
    if plain:
        def _ext_ctx(k):
            # nothing in this run carries a modifier (`plain` pinned every extent to
            # has_modifiers=False), so an extent's sample() reads none of this: an empty context
            # is the same multiplication by 1.0, without rebuilding the loop's context per event.
            return {}

    # --- building a lineage's row ------------------------------------------------------------------
    # For bringing a row up to date from what changed in it: which live drivers read each family (its
    # presence, its module's completion), which read the gene count, and, per event class, which
    # declared families' own rates read each driver. All fixed for the run, so worked out once, in
    # lists, so the order a row visits them in is the run's and not a set's.
    drivers_of_family: dict[int, list] = collections.defaultdict(list)
    count_drivers: list = []
    for src, target in live_rate_reads:
        if target is None:
            count_drivers.append((src, target))
        else:
            for fam in (target if isinstance(target, tuple) else (target,)):
                drivers_of_family[fam].append((src, target))
    readers_of: dict[str, dict] = {key: collections.defaultdict(list) for key in own_keys}
    if any_written:
        for key in own_keys:
            for fam, gene_rate in fam_driven_by_id[key].items():
                for driver in gene_rate.reads:
                    readers_of[key][driver].append(fam)
    counts.recording = any_rows       # a row is brought up to date from the families that changed

    def run_unit_on(key, k, dk) -> float:
        """The run's rate per gene as lineage ``k`` reads it, for the families that did not write one
        of their own. The branch's multiplier belongs here and nowhere else in the own-rate
        arithmetic: it was written on the run's rate, so it scales the genes carrying that rate and
        not the genes whose family replaced it. Read through one function, so `own_sums` and the
        per-gene table `_own_table` hands an event cannot disagree about it."""
        unit = run_gene_rates[key].value(t, dk)
        return unit * lin_mult[key][alive[k]] if varying[key] else unit

    def w_from_rate_row(label, rate, k, size, n_chrom, values) -> float:
        """Lineage ``k``'s weight for one event class, read off the rate itself.

        Three shapes rather than one context built per call: a driven rate reads the branch's driver
        value, a rate varying among lineages reads the branch's drawn multiplier, and a rate doing
        both reads both. The driven-only shape is the expression it always was, so a conditioned run
        is unchanged to the last bit."""
        if not varying[label]:
            return rate.effective(copies=size, lineages=1, chromosomes=n_chrom, time=t,
                                  drivers=values)
        factor = lin_mult[label][alive[k]]
        if not driven[label]:
            return rate.effective(copies=size, lineages=1, chromosomes=n_chrom, time=t,
                                  carried_factor=factor)
        return rate.effective(copies=size, lineages=1, chromosomes=n_chrom, time=t,
                              drivers=values, carried_factor=factor)

    def w_from_family_row(label, k) -> float:
        """The same weight for a rate carrying a per-family draw **and** a per-lineage one: the
        family multipliers summed over the lineage's genes, times the branch's multiplier, times the
        run's unit rate. `fw` is left unscaled, because the segment pick reads it and a family that
        writes its own rate is measured against it (`set_own`)."""
        unit = _rates[label].effective(copies=1, lineages=1, chromosomes=1, time=t)
        return unit * lin_mult[label][alive[k]] * rows.fw[label][k]

    def family_total(label, rate, one, live=True) -> float:
        """One event class's total when the run draws per family: the unit rate times the weights
        summed over the live genes. A rate that also varies among lineages keeps that product per
        lineage, in its `w` row, because the branch's multiplier cannot be factored out of a sum over
        branches. Defined here rather than in the loop, which runs once per event."""
        if not live:
            return 0.0
        if label in rows.w:
            return rows.w[label].total
        return rate.effective(**one) * rows.fw[label].total

    def set_own(k: int, key: str, unit: float, owned: dict, written: float, covered: float,
                whole: float, largest: float) -> None:
        """Write lineage ``k``'s own-rate row for one event class. Every gene carries its family's own
        rate or the run's: the families that write one take ``covered`` of the lineage's ``whole``
        weight and give ``written``, and every other gene carries the run's rate on the rest."""
        if not owned:
            written = covered = 0.0       # nothing is left to carry a rounding error forward
        rows.own_sums[key].set(k, max(0.0, unit * (whole - covered) + written))
        rows.own_units[key][k] = unit
        rows.own_owned[key][k] = owned
        rows.own_written[key][k] = written
        rows.own_covered[key][k] = covered
        rows.own_largest[key][k] = largest

    def build_whole(k: int) -> None:
        """Build lineage ``k``'s row from its genome and counts alone."""
        counts.take_changed(k)            # a row built whole already holds whatever changed
        genome = gen[k]
        if any_driven:
            # A driven rate differs from lineage to lineage, so it is summed **over the living
            # lineages**, each read with its own driver value, its own gene count and its own
            # chromosome count — and the weights are kept, because the affected lineage must then be
            # drawn with them too. The gene count sits inside the weight, which is what makes a driven
            # per-copy rate a two-stage pick (a lineage, then a gene in it) rather than the one-stage
            # lineage draw a per-lineage rate takes. A per-family draw and a Driven cannot both be
            # set, so a `w` row built from a driver never sits beside an `fw` one; a `w` row built
            # from a per-lineage draw can, and `w_from_family_row` is where it is.
            values = {**{key: trajs[key].value(alive[k], t) for key in trajs},
                      **{src: _live_value(target, genome, counts.of(k))
                         for src, target in live_rate_reads}}
            rows.drivers[k] = values
        if any_driven or any_lineage:
            # the driver values this lineage reads, or None for a run with no driver at all: a rate
            # varying among lineages needs a row whether or not anything here is driven
            reads = rows.drivers[k] if any_driven else None
            size, n_chrom = _genome_size(genome), len(genome)
            for label in w_from_rate:
                rows.w[label].set(k, w_from_rate_row(label, _rates[label], k, size, n_chrom, reads))
        if any_family:
            # A per-copy rate pools over genes, so with per-family weights the total is the unit rate
            # times those weights summed over the live genes — and the run must then be drawn with the
            # same weights, or the rate would say one thing and the picking another. Summed per
            # lineage, so the lineage pick can reuse them. On a circular chromosome
            # ``Σ_s mean_w(s, m)`` is exactly this sum for every run size, which is why no per-size
            # term appears here (SPEC §6).
            for key, mult in fam_mult.items():
                rows.fw[key].set(k, sum(mult[g.family] for chrom in genome for g in chrom.genes))
            for label in w_from_family:
                rows.w[label].set(k, w_from_family_row(label, k))
        if any_written:
            # A family's own rate (SPEC §6): every gene carries a rate, its family's own when the
            # family writes one and the run's otherwise, and an event's total is their sum. The
            # acting lineage is drawn by its genes' summed rates and the segment by the mean rate of
            # the genes it covers, through the path the per-family draws take. The sum runs over the
            # families that write their own rate; every other gene carries the run's rate, so their
            # share is the lineage's whole weight less what the writing families take.
            held = counts.of(k)
            dk: dict[str, Any] = {"drivers": rows.drivers[k]} if any_driven else {}
            for key in own_keys:
                own_mult = fam_mult[key] if any_family else None
                unit = run_unit_on(key, k, dk)
                owned = {fam: value for fam, value in fam_fixed_by_id[key].items() if held[fam]}
                for fam, own_rate in fam_driven_by_id[key].items():
                    if held[fam]:
                        owned[fam] = own_rate.value(t, dk)
                written = covered = 0.0
                for fam, value in owned.items():
                    written += held[fam] * value
                    covered += held[fam] * (own_mult[fam] if own_mult is not None else 1.0)
                # the lineage's whole weight for this event: its summed per-family draws where the
                # run has them, its gene count otherwise
                whole = float(rows.fw[key][k] if any_family else _genome_size(genome))
                # every gene carries its family's own rate or the run's rate times its family's
                # draw, so the larger of those bounds all of them
                largest = max(unit * (own_mult.largest if own_mult is not None else 1.0),
                              max(owned.values(), default=0.0))
                set_own(k, key, unit, owned, written, covered, whole, largest)

    def bring_up_to_date(k: int) -> None:
        """Bring lineage ``k``'s row, built before, up to date from what changed in it since: the
        families whose copy number moved, the drivers those families move, and the families whose own
        rate reads a driver that moved. Nothing else in the row can have changed."""
        changed = counts.take_changed(k)
        genome = gen[k]
        held = counts.of(k)
        moved: list = []                  # the driver names whose value moved, in a fixed order
        if any_driven:
            values = rows.drivers[k]
            # the drivers that read a family that changed, and the gene count whenever the genome did
            asked = dict.fromkeys([pair for fam in changed for pair in drivers_of_family.get(fam, ())]
                                  + count_drivers)
            for src, target in asked:
                now = _live_value(target, genome, held)
                if now != values[src]:
                    if not moved:
                        values = dict(values)     # a new entry, so a kept snapshot is not edited
                    values[src] = now
                    moved.append(src)
            rows.drivers[k] = values
        if any_driven or any_lineage:
            # the driver values this lineage reads, or None for a run with no driver at all: a rate
            # varying among lineages needs a row whether or not anything here is driven
            reads = rows.drivers[k] if any_driven else None
            size, n_chrom = _genome_size(genome), len(genome)
            for label in w_from_rate:
                rows.w[label].set(k, w_from_rate_row(label, _rates[label], k, size, n_chrom, reads))
        if any_family and changed:
            for key, mult in fam_mult.items():
                weights = rows.fw[key]
                weights.set(k, weights[k] + sum(delta * mult[fam] for fam, delta in changed.items()))
            for label in w_from_family:
                rows.w[label].set(k, w_from_family_row(label, k))
        if any_written:
            dk: dict[str, Any] = {"drivers": rows.drivers[k]} if any_driven else {}
            size_now = float(_genome_size(genome))
            for key in own_keys:
                own_mult = fam_mult[key] if any_family else None
                fixed, driven_own, readers = fam_fixed_by_id[key], fam_driven_by_id[key], readers_of[key]
                unit = run_unit_on(key, k, dk)
                owned = rows.own_owned[key][k]
                written, covered = rows.own_written[key][k], rows.own_covered[key][k]
                largest = max(rows.own_largest[key][k],
                              unit * (own_mult.largest if own_mult is not None else 1.0))
                affected = [fam for fam in changed if fam in fixed or fam in driven_own]
                for src in moved:
                    affected.extend(readers.get(src, ()))
                if affected:
                    owned = dict(owned)
                    for fam in dict.fromkeys(affected):      # each family once, in a fixed order
                        weight = own_mult[fam] if own_mult is not None else 1.0
                        if fam in owned:                     # take out what it carried before
                            before = held[fam] - changed.get(fam, 0)
                            written -= before * owned.pop(fam)
                            covered -= before * weight
                        if held[fam]:                        # and put in what it carries now
                            value = fixed[fam] if fam in fixed else driven_own[fam].value(t, dk)
                            owned[fam] = value
                            written += held[fam] * value
                            covered += held[fam] * weight
                            largest = max(largest, value)
                whole = rows.fw[key][k] if any_family else size_now
                set_own(k, key, unit, owned, written, covered, whole, largest)

    bar = progress_bar(len(schedule), "genomes", unit="branch", enabled=progress)
    si = 0
    while si < len(schedule):
        if to_disk is not None:
            to_disk.flush_if_full()      # between steps, so a batch always holds whole events
        if _CHECK_COUNTS:
            _check_counts(gen, counts, rows)
        bar.to(si)
        n = total_copies
        k_alive = len(alive)
        c = total_chromosomes
        can_xfer = n > 0 and (k_alive >= 2 or self_transfer)
        next_species = schedule[si][0]
        # a family placed by `origins=` originates at a fixed instant, so it joins the horizon like
        # any other breakpoint: the waiting time can never step over it
        next_plant = plants[plant_i][0] if plant_i < len(plants) else math.inf
        if plain:
            # no modifier on any rate or extent: each total is scope(base) exactly — a gene rate
            # times the live genes, a chromosome rate times the standing chromosomes, the two
            # originations per living lineage — none of them ever changes on its own (next_change
            # is inf), and nothing below reads a context or a weight.
            w = no_weights
            fw = None
            own_pick = {}
            r_dup = dup_base * n
            r_los = los_base * n
            r_tra = tra_base * n if can_xfer else 0.0
            r_inv = inv_base * n
            r_trp = trp_base * n
            r_trl = trl_base * n
            r_org = org_base * k_alive
            r_fis = fis_base * c
            r_fus = fus_base * c
            r_cor = cor_base * k_alive
            r_clo = clo_base * c
            horizon = min(next_species, next_plant)
        else:
            ctx = {"copies": n, "lineages": k_alive, "chromosomes": total_chromosomes, "time": t}
            # A gene-level event counted PER LINEAGE is counted per lineage that HOLDS a gene: an
            # empty genome offers nothing to act on, so it must not take a share of the total and
            # then be picked with no victim inside it. Built only when some rate needs it, so the
            # per-copy path does exactly the work it did before.
            if any_per_lineage:
                gene_hosts = [k for k in range(k_alive) if _genome_size(gen[k])]
                gene_ctx = {**ctx, "lineages": len(gene_hosts)}
            else:
                gene_hosts, gene_ctx = None, ctx
            # Every number below reads one lineage on its own, so the step builds the rows of the
            # lineages an event changed and reads the rest as they stand (`_LineageRows`).
            if any_rows:
                for k, whole in rows.take_stale():
                    if whole:
                        build_whole(k)
                    else:
                        bring_up_to_date(k)
                if _CHECK_ROWS:
                    kept = rows.snapshot()
                    rows.touched_all()
                    for k, _whole in rows.take_stale():
                        build_whole(k)
                    rows.check_against(kept)
            w = rows.w
            fw = rows.fw if any_family else None

            def _r(label, pooled, live=True):
                """The total for one event class: summed per-lineage when driven, pooled when not."""
                if not live:
                    return 0.0
                return w[label].total if label in w else pooled

            if fw is not None:   # the run draws per family: the same test as `any_family`
                one = {"copies": 1, "lineages": 1, "chromosomes": 1, "time": t}
                r_dup = family_total("duplication", dup, one, live=bool(n))
                r_los = family_total("loss", los, one, live=bool(n))
                r_tra = family_total("transfer", tra, one, live=can_xfer)
                r_inv = family_total("inversion", inv, one, live=bool(n))
                r_trp = family_total("transposition", trp, one, live=bool(n))
                r_trl = family_total("translocation", trl, one, live=bool(n))
            else:
                # each gene-level rate is read in the context its own scope asks for: `gene_ctx` counts
                # only the occupied genomes, which is what a per-lineage budget is counted over
                def _gc(label):
                    return gene_ctx if per_lineage[label] else ctx

                r_dup = _r("duplication", dup.effective(**_gc("duplication")) if n else 0.0, live=bool(n))
                r_los = _r("loss", los.effective(**_gc("loss")) if n else 0.0, live=bool(n))
                r_tra = _r("transfer", tra.effective(**_gc("transfer")) if can_xfer else 0.0,
                           live=can_xfer)
                r_inv = _r("inversion", inv.effective(**_gc("inversion")) if n else 0.0, live=bool(n))
                r_trp = _r("transposition", trp.effective(**_gc("transposition")) if n else 0.0,
                           live=bool(n))
                r_trl = _r("translocation", trl.effective(**_gc("translocation")) if n else 0.0,
                           live=bool(n))
            if any_written:
                if _CHECK_OWN_SUMS:
                    held_now = [counts.of(k) for k in range(k_alive)]
                    for key in own_keys:
                        _check_own_sums(own_pick[key][0], held_now, own_pick[key][1])
                if "duplication" in own_pick:
                    r_dup = own_pick["duplication"][0].total
                if "loss" in own_pick:
                    r_los = own_pick["loss"][0].total
                if "transfer" in own_pick:
                    r_tra = own_pick["transfer"][0].total if can_xfer else 0.0
            r_org = _r("origination", org.effective(**ctx))                 # per lineage
            r_fis = _r("fission", fis.effective(**ctx) if c else 0.0, live=bool(c))  # per chromosome
            r_fus = _r("fusion", fus.effective(**ctx) if c else 0.0, live=bool(c))
            r_cor = _r("chromosome_origination", cor.effective(**ctx))      # per lineage (de-novo replicon)
            r_clo = _r("chromosome_loss", clo.effective(**ctx) if c else 0.0, live=bool(c))
            # the next instant a rate changes on its own, asked only of the rates that ever do — a
            # family's own schedule among them (`timed_rates`)
            horizon = min(next_species, next_plant)
            for timed in timed_rates:
                horizon = min(horizon, timed.next_change(t))
            if any_driven:  # a driven rate also changes when its driver switches mid-branch — step there
                horizon = min(horizon, min((trajs[key].next_change(alive[k], t) for key in trajs
                                            for k in range(k_alive)), default=math.inf))
            def _ext_ctx(k):
                """The context an extent is sampled in, on the lineage the event landed on.

                It cannot be built before the lineage is drawn, because a driven extent is read on
                the **acting** lineage at the instant the event fires — which is also why an extent
                adds no Gillespie breakpoint and never enters the horizon above (SPEC §6). With no
                driven extent this is the same context the rates were read in.

                The rest of `ctx` — the gene, lineage and chromosome counts — goes with it, because
                `Modifier.implemented_for` promises this engine supplies them and a modifier of
                your own is admitted onto an extent by the same gate that admits it onto a rate.
                Handing an extent a thinner context meant one gate certifying two different
                contracts: a modifier written the documented way read zeros, and one with a
                required keyword died mid-run."""
                # `ctx` was snapshotted at the top of the loop, before `t` advanced to the firing
                # instant, so `time` has to be taken fresh: an extent's own breakpoints are kept out
                # of the horizon, so a schedule's breakpoint routinely falls inside a stretch, and
                # reading the stale `t` would size the event on the wrong side of it.
                if not any_ext_driven:
                    return {**ctx, "time": t}
                return {**ctx, "time": t,
                        "drivers": {**{key: resolved[key].value(alive[k], t) for key in resolved},
                                    **{src: _live_value(target, gen[k], counts.of(k))
                                       for src, target in live_ext_reads}}}
        total = (r_dup + r_los + r_org + r_tra + r_inv + r_trp + r_trl
                 + r_fis + r_fus + r_cor + r_clo)

        if total > 0.0:
            t_ev = t + float(rng.exponential(1.0 / total))
            if t_ev < horizon:  # a genome event fires before the alive set or a rate changes
                t = t_ev
                r = float(rng.random()) * total
                b_los = r_dup + r_los                    # cumulative bounds, in the firing order below
                b_org = b_los + r_org
                b_tra = b_org + r_tra
                b_inv = b_tra + r_inv
                b_trp = b_inv + r_trp
                b_trl = b_trp + r_trl
                b_fis = b_trl + r_fis
                b_fus = b_fis + r_fus
                b_cor = b_fus + r_cor                    # ... and the remainder (to total) is clo
                if r < r_dup:                            # every gene-level event acts on an extent
                    picked = _pick_event_run(rng, gen, n, rows.genes, fw, fam_mult, "duplication", dup_ext,
                                             _ext_ctx, w.get("duplication"),
                                             gene_hosts if per_lineage["duplication"] else None,
                                             own=own_pick.get("duplication"))
                    if picked is not None:
                        k, ci, j, m = picked
                        copied = _run_families(gen[k][ci], j, m)
                        if not _run_over_cap(counts.of(k), copied, cap):
                            total_copies += _duplicate(gen[k][ci], j, m, tree.nodes[alive[k]], t,
                                                       events, event_positions, new_gene)
                            counts.added_all(k, copied)
                            rows.touched(k, gen[k])
                elif r < b_los:
                    picked = _pick_event_run(rng, gen, n, rows.genes, fw, fam_mult, "loss", los_ext,
                                             _ext_ctx, w.get("loss"),
                                             gene_hosts if per_lineage["loss"] else None,
                                             own=own_pick.get("loss"))
                    if picked is not None:
                        k, ci, j, m = picked
                        taken = _run_families(gen[k][ci], j, m)
                        if _lose_at(gen[k][ci], j, m, tree.nodes[alive[k]], t, events, event_positions):
                            total_copies -= m
                            counts.removed_all(k, taken)
                            rows.touched(k, gen[k])
                elif r < b_org:
                    # origination is per lineage: a uniform lineage, or one drawn by its own rate
                    # when that rate is driven (the same weights the total was summed with)
                    k = (w["origination"].pick(rng) if "origination" in w
                         else int(rng.integers(k_alive)))
                    counts.added(k, _originate(gen[k], tree.nodes[alive[k]], t, events, event_positions,
                                               new_gene, new_family, rng))
                    total_copies += 1
                    rows.touched(k, gen[k])
                elif r < b_tra:
                    picked = _pick_event_run(rng, gen, n, rows.genes, fw, fam_mult, "transfer", tra_ext,
                                             _ext_ctx, w.get("transfer"),
                                             gene_hosts if per_lineage["transfer"] else None,
                                             own=own_pick.get("transfer"))
                    if picked is not None:                # driven: the weighted lineage is the DONOR
                        kd, cdi, jd, m = picked
                        delta, kr = _do_transfer(rng, tree, alive, gen, counts, kd, cdi, jd, m, t,
                                                 events, event_positions, new_gene, transfer_to,
                                                 replacement, self_transfer, depth, cap,
                                                 to_traj, group_of, fam_choice)
                        total_copies += delta
                        if kr is not None:   # the recipient's gene content changed, the donor's did not
                            rows.touched(kr, gen[kr])
                elif r < b_inv:
                    picked = _pick_event_run(rng, gen, n, rows.genes, fw, fam_mult, "inversion", inv_ext,
                                             _ext_ctx, w.get("inversion"),
                                             gene_hosts if per_lineage["inversion"] else None)
                    if picked is not None:                # the run starts at a gene, so: per copy
                        k, ci, i0, m = picked
                        _invert(gen[k][ci], i0, m, tree.nodes[alive[k]], t, rearrangements)
                elif r < b_trp:
                    picked = _pick_event_run(rng, gen, n, rows.genes, fw, fam_mult, "transposition", trp_ext,
                                             _ext_ctx, w.get("transposition"),
                                             gene_hosts if per_lineage["transposition"] else None)
                    if picked is not None:
                        k, ci, i0, m = picked
                        _transpose(gen[k][ci], i0, m, tree.nodes[alive[k]], t, rearrangements, rng,
                                   inversion_probability)
                elif r < b_trl:
                    picked = _pick_event_run(rng, gen, n, rows.genes, fw, fam_mult, "translocation", trl_ext,
                                             _ext_ctx, w.get("translocation"),
                                             gene_hosts if per_lineage["translocation"] else None)
                    if picked is not None:
                        k, ci, j, m = picked
                        _translocate(gen[k], ci, j, m, tree.nodes[alive[k]], t, rearrangements, rng,
                                     inversion_probability)
                elif r < b_fis:
                    picked = _pick_chromosome(rng, gen, c, rows.chromosomes, w.get("fission"))
                    if picked is not None:
                        k, ci = picked
                        dc, dg = _fission(gen[k], ci, tree.nodes[alive[k]], t, chromosome_events,
                                          new_chromosome, rng)
                        total_chromosomes += dc
                        total_copies += dg
                        if dc or dg:
                            rows.touched(k, gen[k])
                elif r < b_fus:
                    picked = _pick_chromosome(rng, gen, c, rows.chromosomes, w.get("fusion"))
                    if picked is not None:
                        k, ci = picked
                        dc, dg = _fusion(gen[k], ci, tree.nodes[alive[k]], t, chromosome_events,
                                         new_chromosome, rng)
                        total_chromosomes += dc
                        total_copies += dg
                        if dc or dg:
                            rows.touched(k, gen[k])
                elif r < b_cor:
                    # chromosome origination is per lineage, uniform or driven, exactly as origination
                    k = (w["chromosome_origination"].pick(rng)
                         if "chromosome_origination" in w else int(rng.integers(k_alive)))
                    dc, dg = _chromosome_originate(gen[k], tree.nodes[alive[k]], t, chromosome_events,
                                                   new_chromosome)
                    total_chromosomes += dc
                    total_copies += dg
                    if dc or dg:
                        rows.touched(k, gen[k])
                else:
                    picked = _pick_chromosome(rng, gen, c, rows.chromosomes, w.get("chromosome_loss"))
                    if picked is not None:
                        k, ci = picked
                        taken = [g.family for g in gen[k][ci].genes]
                        dc, dg = _chromosome_lose(gen[k], ci, tree.nodes[alive[k]], t, events,
                                                  event_positions, chromosome_events)
                        if dc:
                            counts.removed_all(k, taken)
                        total_chromosomes += dc
                        total_copies += dg
                        if dc or dg:
                            rows.touched(k, gen[k])
                continue

        if horizon == next_species:  # advance to the tree's next event(s); process the whole tie-batch
            t = next_species
            if time_varying:             # time moved, and something here moves with it
                rows.touched_all()
            while si < len(schedule) and schedule[si][0] == t:
                i = schedule[si][1]
                g = gen[pos[i]]
                if to_disk is None:
                    genomes[i] = tuple(Chromosome(c.id, c.topology, tuple(c.genes)) for c in g)  # freeze
                else:
                    to_disk.branch_ended(i, g)     # written now, and not kept
                total_copies -= sum(len(c.genes) for c in g)
                total_chromosomes -= len(g)
                inherited = counts.retired(pos[i])  # what the daughters below inherit, if any
                rows.retired(pos[i])
                retire(alive, gen, pos, pos[i])
                node = tree.nodes[i]
                if node.children:  # a speciation: re-mint every chromosome and gene id
                    child_genomes: dict[int, list[Chromosome]] = {c: [] for c in node.children}
                    for pchrom in g:
                        dcids = []
                        per_daughter: list[list[GeneEdge]] = []
                        for c in node.children:
                            dcid = new_chromosome()
                            dcids.append(dcid)
                            dgenes, edges = [], []
                            for old in pchrom.genes:  # ZOMBI1: the gene ends and continues, fresh id
                                ng = new_gene(old.family, old.strand)
                                dgenes.append(ng)
                                edges.append(GeneEdge(t, "speciation", c, old.family, ng.id,
                                                   parent=old.id))
                            per_daughter.append(edges)
                            child_genomes[c].append(Chromosome(dcid, pchrom.topology, dgenes))
                        # the ids are minted daughter by daughter (which is what fixes them), but a
                        # gene's two edges are recorded together: one gene ending is one event, and
                        # the log writes it as one row naming both daughters
                        for gene_edges in zip(*per_daughter):
                            events.extend(gene_edges)
                        chromosome_events.append(
                            ChromosomeEvent(t, "speciation", node.id, (pchrom.id,), tuple(dcids)))
                    for c in node.children:
                        cg = child_genomes[c]
                        enter(alive, gen, pos, c, cg)
                        counts.entered_like(inherited)   # a re-id of the parent: same families
                        rows.entered(cg)
                        total_copies += sum(len(ch.genes) for ch in cg)
                        total_chromosomes += len(cg)
                si += 1
        elif plant_i < len(plants) and horizon == next_plant:
            # a placed family arrives — the ordinary origination event, at a time and on a lineage
            # that were chosen rather than drawn. The lineage is live by construction (its time was
            # checked against that branch's own life), and a tie with the tree's schedule falls to
            # the branch above, so the daughters have entered by the time this runs.
            t = horizon
            if time_varying:
                rows.touched_all()
            while plant_i < len(plants) and plants[plant_i][0] == t:
                _, lineage, fam = plants[plant_i]
                _originate(gen[pos[lineage]], tree.nodes[lineage], t, events, event_positions,
                           new_gene, new_family, rng, family=fam)
                counts.added(pos[lineage], fam)
                rows.touched(pos[lineage], gen[pos[lineage]])
                total_copies += 1
                plant_i += 1
        else:
            t = horizon  # a skyline breakpoint: advance and re-evaluate the (now changed) rate
            rows.touched_all()

    bar.close()
    if _CHECK_COUNTS:
        _check_counts(gen, counts, rows)
    links = links_of({**_rates, **_extents}, transfer_to, declared, fam_driven_rates, fam_transfer_to,
                     module_map or {})
    # a declared family's own rate, fixed or driven, replaces the run's, so no draw reaches it there
    own = {key: set(fam_fixed_by_id[key]) | set(fam_driven_by_id[key]) for key in own_keys}
    multipliers = multipliers_of(fam_mult, ORDERED_TARGETS, own) if any_family else {}
    per_lineage_multipliers = (lineage_multipliers_of(lin_mult, tree.labels(),
                                                     ORDERED_LINEAGE_TARGETS)
                               if any_lineage else {})
    if to_disk is not None:
        return to_disk.close(seed=seed, links=links, initial_genome=initial_genome, named=len(named),
                             family_multipliers=multipliers,
                             lineage_multipliers=per_lineage_multipliers)
    return OrderedGenomesResult(tree, genomes, events, rearrangements, chromosome_events, seed,
                                named, module_map, event_positions, initial_genome, links,
                                multipliers, per_lineage_multipliers)


__all__ = ["simulate_genomes_ordered", "OrderedGenomesResult", "Gene", "Chromosome",
           "ChromosomeEvent", "Inversion", "Transposition", "Translocation"]
