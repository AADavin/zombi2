"""Genomes II — ordered: genes carry a position and an orientation, on chromosomes.

The ordered resolution layers **position** over the family D/T/L/O core (Chapter 4). A genome is
no longer a multiset of gene copies but a list of **chromosomes**, each an ordered run of oriented
`Gene`\\ s.

**Every gene-level event acts on an extent** — a run of consecutive genes (the ZOMBI1 model), its
length drawn per event from a distribution (default ``Geometric(mean=1)`` — a single gene). The run
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

**A rate here may be driven by another level** (``ScaledBy``, SPEC §2 and §5) — and so may an extent
(SPEC §6). A driven rate is *per lineage*: it is summed over the living lineages, each read with its
own driver value, and the lineage an event lands on is then drawn with those same weights. Because
the gene-level rates are **per copy**, each lineage's weight carries its own gene count, so the pick
is two-stage — a lineage by its weight, then a gene inside it. A driven extent is read at the instant
an event fires, so it changes how much that event takes and never how often one starts.

``transfer_to`` — who receives — is the third place a driver can sit, and the one that is **not** a
rate: there the mapping's numbers are weights normalised across the candidate recipients, so they
redistribute the same transfers rather than change how many happen (SPEC §5, a weight). Its
four rules and the kernel that reads them are the family core's, shared through ``_transfer``.
"""

from __future__ import annotations

import collections
import math
import pathlib
from typing import Sequence, cast
from dataclasses import dataclass, field
from functools import cached_property


from ..rates.driver import check_mapping_fires, resolve_driver
from ..rng import stream
from ..rates.extent import Extent, as_extent
from ..rates.mapping import check_not_a_kernel
from ..rates.modifiers import (describe, DRAWN, Driven, OnTime, SetBy, cell_name, is_implemented,
                               values_at_birth)
from ..rates.rate import Rate, as_rate
from ..rates.scope import PerChromosome, PerCopy, PerLineage
from ..tree import Tree, as_tree
from .chromosomes import ChromosomeEvent, chromosome_events_tsv, rearrangement_events_tsv
from .family import resolve_modules, resolve_max_family_size
from ._live import enter, retire, weighted_index, without_cyclic_gc
from ._transfer import (mean_root_to_tip, prepare_transfer_to, recipient_index,
                        resolve_transfer_to)
from .._runtime.outputs import fresh_dirs, grouped_dir
from .._runtime.summary import _stats, write_summary
from .._runtime.progress import progress_bar
from .events import (_COLS, Event, GeneEdge, _branches, _name, event_counts, event_rows,
                     events_from_edges, gene_label)
from .gene_trees import GeneTree, gene_trees_from_edges, write_gene_trees
from .profiles import Profiles, profiles_from_genomes

#: The rate grammar this engine supports (SPEC §5) — read by the gate below and by the CLI's help, so
#: a modifier is never advertised without being implemented. The same four the family core takes,
#: because the two are the same model at two resolutions: ``OnTime`` (a skyline in time), ``ScaledBy``
#: (a conditioned or joint driver), ``SetBy`` (a driver that replaces the base rather than scaling
#: it) and a per-family draw (per-family heterogeneity, weighted on the segment an event covers
#: rather than on the gene it started from — SPEC §6). One combination is refused: see the gate.
IMPLEMENTED_MODIFIERS = (OnTime, Driven, SetBy, (DRAWN, "family"))

#: What an **extent** takes here (SPEC §6). An extent takes the modifiers a rate does, and at this
#: resolution that is one fewer: a per-family draw attaches to the *contents*, and an extent is drawn
#: before the run's genes are known — a run covers several families, so there is no one family to
#: draw a factor for. The two lists are declared separately rather than hidden in an ``if``, because
#: the difference is a modelling fact, not an implementation detail.
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

    Topology also decides which chromosome-tier events are legal. A **fusion** joins two chromosomes
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
    ``events`` log, the ``rearrangements`` (inversions) and ``chromosome_events`` (the chromosome
    genealogy) logs, and the ``seed``. The observed genomes are the extant tips; ``profiles`` and
    ``gene_trees`` are derived from the (position-blind) genealogy exactly as for the family core;
    ``gene_order`` reads a node's layout, and ``write`` materialises the chosen outputs."""

    complete_tree: Tree
    genomes: dict[int, tuple[Chromosome, ...]]
    edges: list[GeneEdge]
    rearrangements: list[Inversion | Transposition | Translocation]
    chromosome_events: list[ChromosomeEvent]
    seed: int | None
    #: ``{name: family id}`` for families declared by ``family_names=[…]`` — the handle to a *named* family.
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

    def __repr__(self) -> str:
        return (f"OrderedGenomesResult({len(self.complete_tree.extant_leaves())} extant genomes, "
                f"{len(self.genomes)} nodes, {len(self.edges)} events, "
                f"{len(self.rearrangements)} rearrangements, seed={self.seed})")

    def family_counts(self, node_id: int) -> collections.Counter:
        """A multiset view of one node's genome: ``family id → copy count`` (across all chromosomes)."""
        return collections.Counter(g.family for chrom in self.genomes[node_id] for g in chrom.genes)

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

            switch=0.1 * ScaledBy(g.presence("tox"), {"present": 5.0, "absent": 1.0})
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
        return any(g.family == fid for chrom in self.genomes[node_id] for g in chrom.genes)

    def gene_order(self, node_id: int) -> list[tuple[int, int, int, int, int]]:
        """One node's layout as ``(chromosome, position, strand, family, gene id)`` rows, chromosome
        by chromosome and left to right within each — the ordered analogue of ``family_counts``."""
        return [(chrom.id, pos, g.strand, g.family, g.id)
                for chrom in self.genomes[node_id] for pos, g in enumerate(chrom.genes)]

    @cached_property
    def _extant_genes(self) -> dict[int, tuple[Gene, ...]]:
        """The observed genomes flattened to gene multisets (chromosomes dropped) — the view the
        genealogy-derived, position-blind outputs read."""
        extant = [n.id for n in self.complete_tree.extant_leaves()]
        return {s: tuple(g for chrom in self.genomes[s] for g in chrom.genes) for s in extant}

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
               "chromosome_events", "gene_trees", "species_tree", "summary")

    def write(self, directory, outputs=("events", "profiles", "gene_order", "initial_genome",
                                        "gene_trees", "chromosome_events", "species_tree",
                                        "summary"), *,
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
        extant = [n.id for n in self.complete_tree.extant_leaves()]
        born = {e.family for e in self.edges}
        surviving = {g.family for i in extant for c in self.genomes.get(i, ()) for g in c.genes}
        genes_per_genome = [sum(len(c.genes) for c in self.genomes.get(i, ())) for i in extant]
        chrom_per_genome = [len(self.genomes.get(i, ())) for i in extant]
        rearrangements = collections.Counter(type(r).__name__.lower() for r in self.rearrangements)
        chromosome = collections.Counter(e.kind for e in self.chromosome_events)
        return {
            "level": "genomes",
            "seed": self.seed,
            "resolution": "ordered",
            "events": event_counts(self.edges, t0),
            "families": {"born": len(born), "surviving": len(surviving),
                         "died_out": len(born) - len(surviving),
                         "named": len(self.family_names)},
            "extant_genomes": len(extant),
            "empty_genomes": sum(1 for i in extant
                                 if not any(c.genes for c in self.genomes.get(i, ()))),
            "genes_per_genome": _stats(genes_per_genome),
            "chromosomes_per_genome": _stats(chrom_per_genome),
            # this resolution's own two records: what moved genes without changing their ancestry,
            # and what happened to the replicons carrying them
            "rearrangements": {k: rearrangements.get(k, 0)
                               for k in ("inversion", "transposition", "translocation")},
            "chromosome_events": dict(sorted(chromosome.items())),
        }

    def _initial_genome_tsv(self) -> str:
        """The layout the run started with — ``gene_order.tsv``'s columns without ``lineage``, which
        is the whole point: it belongs to the start of the root branch, not to a node."""
        cols = ("chromosome", "topology", "position", "strand", "family", "copy")
        rows = [f"{chrom.id}\t{chrom.topology}\t{pos}\t{g.strand}\t{g.family}\t{gene_label(g.id)}"
                for chrom in self.initial_genome for pos, g in enumerate(chrom.genes)]
        return "\n".join(["\t".join(cols), *rows]) + "\n"

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
        cols = ("lineage", "chromosome", "topology", "position", "strand", "family", "copy")
        rows: list[str] = []
        for s in sorted(self.genomes):
            topology = {c.id: c.topology for c in self.genomes[s]}
            rows.extend(f"{_name(names, s)}\t{ch}\t{topology.get(ch, '')}\t{p}\t{st}\t{fam}\t"
                        f"{gene_label(gid)}"
                        for (ch, p, st, fam, gid) in self.gene_order(s))
        return "\n".join(["\t".join(cols), *rows]) + "\n"


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


def _coordinates(events, event_positions) -> list[str]:
    """The coordinate cells of every written row, in `event_rows()`'s order — the two are zipped into
    one table, so this walks the same `events_from_edges()` the genealogy writer does.

    A transfer needs *both* of its `EventPosition`\\ s: the departing one for the donor's arc and the
    arriving one for where the block landed. Which branch each is on comes from the copies themselves
    (`_branches()`) — the donor's continuation leads ``children``, the arriving copy follows — so a
    self-transfer, where the two branches are the same, still resolves by kind.

    A copy displaced by a replacing transfer has no row of its own (it is that transfer's second
    parent), so its position is not written: it is named by id, and a replay that tracks ids removes
    it without needing to be told where it sat."""
    where: dict = {}
    for p in event_positions:
        where.setdefault((p.time, *_position_key(p.kind, p.lineage, p.family)), p)
    branch = _branches(events)
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


def _pick_gene(rng, gen, total_copies) -> tuple[int, int, int]:
    """A uniform global gene pick → ``(lineage k, chromosome index ci in gen[k], position j)``.
    Realises per-copy scope across the whole pool: every gene, in any chromosome of any lineage, is
    equally likely.

    One pass over the chromosomes, deliberately: the obvious spelling — ``_genome_size(genome)`` per
    lineage to decide whether the draw lands in it, then ``_gene_in()`` to walk the chosen one again —
    reads every chromosome of every skipped lineage and then re-reads the chosen lineage's, which is a
    per-event cost in the hot Gillespie loop. Counting chromosome by chromosome finds the same gene in
    a single walk. The draw is unchanged, so a run is byte-identical either way."""
    m = int(rng.integers(total_copies))
    for k, genome in enumerate(gen):
        for ci, chrom in enumerate(genome):
            n = len(chrom.genes)
            if m < n:
                return k, ci, m
            m -= n
    raise AssertionError("total_copies out of sync with the genomes")  # unreachable


def _pick_chromosome(rng, gen, total_chromosomes, w=None) -> tuple[int, int] | None:
    """A chromosome pick → ``(lineage k, chromosome index ci in gen[k])``, or ``None`` when there is
    nothing to draw.

    Uniform over the whole pool when ``w`` is ``None``, which realises per-chromosome scope: every
    chromosome, in any lineage, is equally likely. With ``w`` — the per-lineage totals of a **driven**
    tier rate — the lineage is drawn by its own weight and the chromosome uniformly inside it. That is
    the same two-stage shape, because a driven lineage's weight already carries its chromosome count:
    ``base × chromosomes_k × factor_k``. Drawing the lineage uniformly instead would say one thing in
    the total and another in the pick."""
    if w is not None:
        total = sum(w)
        if total <= 0.0:
            return None                     # every living lineage weighs 0: the event cannot happen
        k = weighted_index(rng, w, total)
        return (k, int(rng.integers(len(gen[k])))) if gen[k] else None
    m = int(rng.integers(total_chromosomes))
    for k, genome in enumerate(gen):
        if m < len(genome):
            return k, m
        m -= len(genome)
    raise AssertionError("total_chromosomes out of sync with the genomes")  # unreachable


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


def _pick_run_by_family(rng, genome, mult, ext, ctx=None) -> tuple[int, int, int] | None:
    """A run drawn with the per-family weight on the **segment**, not on its starting gene (SPEC §6).

    Returns ``(chromosome index, start, run size)``, or ``None`` when the genome has no genes to act
    on. The size is drawn first, then the start in proportion to the run's **mean** weight: a run of
    heavily-weighted genes is favoured, a mixed one sits in between, an ordinary one is unweighted.
    Weighting the *starting* gene instead — the obvious implementation — would apply a family's own
    rate to its **neighbours**, and the neighbourhood is reshuffled by every rearrangement, so the
    parameter would not even name a fixed thing over a run.

    Drawing the size before the start is **exact on a circular chromosome**: there ``Σ_s mean_w(s, m)``
    equals ``Σ_g w_g`` for every ``m``, so the total rate carries no per-size term and the two draws
    factorise cleanly. On a **linear** chromosome the run is clamped by its start, so that identity
    holds only approximately — the same edge effect, from the same cause, that clamping already gives
    a linear run.

    With uniform weights every mean is 1, the start pick is uniform and the size distribution is
    untouched, so a run that sets no per-family weight is byte-identical to one taking the plain path.
    """
    sums = [sum(mult[g.family] for g in c.genes) for c in genome]
    total = sum(sums)
    if total <= 0.0:
        return None
    ci = weighted_index(rng, sums, total)
    chrom = genome[ci]
    n = len(chrom.genes)
    m = min(max(1, int(ext.sample(rng, **(ctx or {})))), n)
    means = _run_means(chrom, mult, m)
    s = weighted_index(rng, means, sum(means))
    return ci, s, (m if chrom.topology == "circular" else min(m, n - s))


def _run_over_cap(genome, chrom, start, m, cap) -> bool:
    """Whether copying the run ``[start, start+m)`` would take any family it covers past ``cap``.

    The segmental answer to the per-genome family quota. At the family resolution the unit is one
    gene, so the question is simply "is this family already full?"; here a run may carry several
    families, and several copies of one, so the test is *current + carried > cap* for each of them.
    That reduces to exactly the family-resolution condition when the run is a single gene.

    **The whole run is refused, never part of it.** Clipping the run to the genes still under quota
    would be a different process — it would quietly reshape the extent distribution, making runs
    shorter precisely where the genome is crowded. Refusing outright is Poisson thinning on a
    condition that reads only the current state, so what is kept is a clean process; a clipped run
    would not be."""
    if cap is None:
        return False
    n = len(chrom.genes)
    carried: dict[int, int] = {}
    for i in range(m):
        f = chrom.genes[(start + i) % n].family
        carried[f] = carried.get(f, 0) + 1
    for f, k in carried.items():
        have = sum(1 for c in genome for g in c.genes if g.family == f)
        if have + k > cap:
            return True
    return False


def _pick_event_run(rng, gen, n, fw, fam_mult, key, ext, ext_ctx, w=None):
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

    The two weighted paths are mutually exclusive: the engine refuses a per-family draw and ``ScaledBy`` in
    one run, because combining them would weight by the product of a lineage factor and a segment
    factor, which is a model neither of them is on its own."""
    if w is not None:
        total = sum(w)
        if total <= 0.0:
            return None                     # every living lineage weighs 0: the event cannot happen
        k = weighted_index(rng, w, total)
        size = _genome_size(gen[k])
        if not size:  # only via weighted_index's r == total float guard — a zero-weight lineage has
            return None                     # no gene to act on, so the event is declined (thinning)
        ci, j = _gene_in(gen[k], int(rng.integers(size)))
        return k, ci, j, _extent(rng, ext, gen[k][ci], j, ext_ctx(k))
    if fw is None:
        k, ci, j = _pick_gene(rng, gen, n)
        return k, ci, j, _extent(rng, ext, gen[k][ci], j, ext_ctx(k))
    lw = fw[key]
    total = sum(lw)
    if total <= 0.0:
        return None
    k = weighted_index(rng, lw, total)
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

def _originate(genome, node, t, events, positions, new_gene, new_family, rng) -> None:
    """A new gene family arises: mint a single founding gene (a family is born once — no extent) on
    a uniformly-chosen chromosome at a uniformly-chosen position (strand ``+1``), and record it."""
    chrom = genome[int(rng.integers(len(genome)))]
    fam = new_family()
    g = new_gene(fam, +1)
    at = int(rng.integers(len(chrom.genes) + 1))
    _live(chrom).insert(at, g)
    events.append(GeneEdge(t, "origination", node.id, fam, g.id))
    positions.append(EventPosition(t, "origination", node.id, chrom.id, at, 1, family=fam))


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
    karyotype is `_chromosome_lose()`'s job, an event of its own at the chromosome tier.

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


def _do_transfer(rng, tree, alive, gen, kd, cdi, jd, m, t, events, positions, new_gene,
                 transfer_to, replacement, self_transfer, depth, cap=None,
                 to_traj=None, groups=None) -> int:
    """The segment ``[jd, jd+m)`` on the donor's chromosome ``cdi`` transfers to a contemporaneous
    recipient chosen by ``transfer_to``: each gene ends → a continuation on the donor branch and a
    transferred copy on the recipient (a horizontal gene-tree edge). The run may wrap position 0 on a
    circular donor chromosome. The transferred copies arrive as a block at a random position on a
    uniformly-chosen recipient chromosome (strands travel with them). Returns the change in total gene
    count: ``+m`` additive, minus one per homologous copy displaced under ``replacement``.

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
    untouched."""
    donor = alive[kd]
    if transfer_to == "uniform":
        # O(1) uniform recipient — the same single draw as recipient_index's
        # cand[rng.integers(len(cand))] over every alive lineage but the donor; the donor-skip is a
        # +1 index shift, so no O(alive) candidate list is built per transfer (see family._do_transfer).
        npool = len(alive) if self_transfer else len(alive) - 1
        if npool <= 0:
            return 0
        i = int(rng.integers(npool))
        kr = i if (self_transfer or i < kd) else i + 1
    else:  # the weighted rules (Distance / Clades / Driven) weigh every candidate — O(alive)
        cand = [k for k in range(len(alive)) if self_transfer or k != kd]
        if not cand:                                   # the uniform branch's npool guard, restated
            return 0
        kr = recipient_index(rng, tree, alive, cand, donor, t, transfer_to, depth, to_traj, groups)
        if kr is None:                                 # every candidate weighs 0 — no-op (see above)
            return 0
    recipient = alive[kr]
    rgenome = gen[kr]
    jd = _anchor(gen[kd][cdi], jd, m)
    segment = gen[kd][cdi].genes[jd:jd + m]
    if _run_over_cap(rgenome, gen[kd][cdi], jd, m, cap):   # the recipient is full: same thinning
        return 0
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
                displaced[x.id] = victim.id
                positions.append(EventPosition(t, "loss", recipient, rgenome[ci].id, p, 1))
                delta -= 1
    rchrom = rgenome[int(rng.integers(len(rgenome)))]   # arrive as a block on a random recipient chromosome
    pos = int(rng.integers(len(rchrom.genes) + 1))
    rchrom.genes[pos:pos] = xfers
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
    return delta


# --- the chromosome tier: events that change chromosome number (the network dynamics) -------------
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

    This is the same floor `_lose_at()` enforces one tier down, for the same reason: a chromosome
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
                             initial_families=100, family_names=None, modules=None,
                             max_family_size=10, seed=None,
                             progress=False) -> OrderedGenomesResult:
    """Evolve ordered genomes — genes with a position and an orientation, on chromosomes — along a
    species tree, by the D/T/L/O core plus segmental rearrangements and the chromosome tier.

    **Every gene-level event acts on an *extent*** — a run of consecutive genes (the ZOMBI1 model):
    ``duplication`` copies the run in tandem, ``loss`` removes it, ``transfer`` sends it to a
    contemporaneous recipient as a block, ``inversion`` reverses it (flipping strands), ``transposition``
    relocates it elsewhere on the same chromosome, and ``translocation`` moves it to a different
    chromosome. The run's length is drawn per event from ``<event>_extent`` (a distribution,
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
    genes that starts at one of them; the chromosome tier ``fission``/``fusion``/``chromosome_loss``
    is **per chromosome**; and the two events that make something from nothing,
    ``origination``/``chromosome_origination``, are **per lineage**. The
    run starts with ``chromosomes`` chromosomes of the given ``topology``, across which the
    ``initial_families`` founding genes are dealt **round-robin**; ``family_names=["toxin", …]`` additionally
    declares **named** families (remembered in ``result.family_names`` for ``result.has_family(node,
    "toxin")``), as in the family core; ``replacement`` / ``self_transfer`` behave as in the family
    core. So does ``transfer_to``, which **chooses who receives** — ``"uniform"``,
    ``"distance"`` / ``Distance(decay=)`` (closer relatives likelier), ``Clades({...}, Between({...}))``
    (weight by the donor's and recipient's named clade) or ``ScaledBy(driver, mapping)`` (weight by
    another level; see below). What moves is a block of genes rather than a single copy, and the block
    arrives whole, so the rule chooses the recipient lineage exactly as it does at the family
    resolution.

    The **chromosome tier** changes chromosome *number*: ``fission`` (split), ``fusion`` (merge,
    between two chromosomes of the **same topology** — the reticulation; a ring and a molecule with
    two ends cannot become one molecule, so a genome of one of each never fuses),
    ``chromosome_origination`` (a de-novo replicon), ``chromosome_loss`` (a whole
    chromosome and its genes die; never the genome's last). Chromosomes carry identity — re-minted at
    every event that reshapes them — so ``chromosome_events`` is the true reticulating chromosome
    genealogy, rooted at the initial and de-novo originations. Deterministic given ``seed``.

    **Conditioning (a trait drives a rate).** Any rate here may be *driven by another level* —
    ``inversion = 0.3 * ScaledBy(habitat, {"host": 4.0, "free": 1.0})`` scales each lineage's
    inversion rate by the habitat on that branch, read from a trait grown first (the finished
    ``TraitsResult``, or the ``trait_events.tsv`` it wrote). A driven rate is then *per lineage*: it is
    summed over the living lineages, each read with its own gene count, chromosome count and driver
    value; the lineage an event lands on is drawn with those same weights, and the gene inside it
    uniformly, because the gene count is already in the weight. The Gillespie steps at **every**
    mid-branch switch of the driver rather than averaging over a branch (SPEC §2). For ``transfer``
    the driven lineage is the **donor**, so a driven ``transfer`` says how often a lineage *donates*.

    **Conditioning (a trait drives who receives).** ``transfer_to = Weights(driver, mapping)`` is
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
    (SPEC §6) — ``inversion_extent = 4 * ScaledBy(habitat, {"host": 3.0, "free": 1.0})`` makes a
    host-restricted lineage invert *longer runs of genes*, which is a different statement from raising
    its inversion rate. An extent's modifier is read at the instant an event fires, so it changes how
    much a run takes and never how often one starts, and it adds no Gillespie breakpoint.

    a per-family draw and ``ScaledBy`` cannot be set in the same run: one weights lineages by a driver and
    the other weights the segment by what it covers.
    """
    tree = as_tree(tree, level="genomes")
    labels = _topologies(chromosomes, topology)
    n_initial_chrom = chromosomes
    # this slice implements each event's default scope and the three modifiers IMPLEMENTED_MODIFIERS
    # declares: OnTime (skyline), Driven (a conditioned/joint driver, per lineage) and a per-family draw —
    # the last with the weight on the SEGMENT rather than on its starting gene (SPEC §6, and
    # _pick_run_by_family). A scope override or a clade-drift modifier is a later slice, so reject it
    # rather than silently mis-scale (see the family engine for the reasoning).
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
        if not isinstance(rate.scope, want):
            raise ValueError(
                f"{label} has a {type(rate.scope).__name__} scope, but the ordered genome engine "
                f"takes only {want.__name__} for {label} this slice — scope overrides are a later slice."
            )
        for m in rate.modifiers:
            if m.reads == (DRAWN, "family") and label == "origination":
                raise ValueError(
                    "origination carries a per-family draw, but origination is the rate at which families are "
                    "CREATED — when it is read there is no family yet to have drawn a factor for. "
                    "Put Drawn(per='family') on duplication, transfer, loss, inversion, transposition "
                    "or translocation; writing one such object on several of them gives a "
                    "family-wide tempo, since one object is one draw.")
            if m.reads == (DRAWN, "family") and not isinstance(rate.scope, PerCopy):
                raise ValueError(
                    f"{label} carries a per-family draw on a {type(rate.scope).__name__} scope. A per-family "
                    f"weight has to reach the genes an event covers, so it applies to the per-copy "
                    f"gene events only — not to the chromosome tier, which acts on whole replicons.")
            if isinstance(m, Driven):
                check_not_a_kernel(m.mapping, label=label)
            if not is_implemented(m, IMPLEMENTED_MODIFIERS, "genomes.ordered"):
                raise ValueError(
                    f"{label} carries {describe(m)}, which the ordered genome engine does not "
                    f"support. It takes OnTime (skyline), ScaledBy (a conditioned or joint driver) and "
                    f"Drawn(per='family') (per-family heterogeneity, weighted on the segment an event covers). "
                    f"Clade drift is not implemented yet."
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
    if any(m.reads == (DRAWN, "family") for r in _rates.values() for m in r.modifiers) and \
            any(isinstance(m, Driven) for r in _rates.values() for m in r.modifiers):
        raise ValueError(
            "a per-family draw and a driver on the same run is a later slice: one weights lineages by a "
            "driver and the other weights the segment by what it covers, so combining them means "
            "weighting by the product. Use one or the other for now.")
    # per-event extent distributions (segment size in genes); a bare number is the mean, None a single gene
    def _ext_spec(spec, label):
        """One event's extent (SPEC §6): ``base × modifiers``, no scope, in **genes** here. An extent
        takes the modifiers a rate takes at this resolution minus a per-family draw (see
        `IMPLEMENTED_EXTENT_MODIFIERS`), and they scale the size drawn — ``OnTime`` in time, ``ScaledBy``
        on the lineage the event lands on."""
        e = as_extent(spec)
        rate_slot = label.removesuffix("_extent")
        for m in e.modifiers:
            if isinstance(m, Driven):
                check_not_a_kernel(m.mapping, label=label)
            if m.reads == (DRAWN, "family"):
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
    family_names = list(family_names) if family_names is not None else []
    for name in family_names:
        if not isinstance(name, str) or not name.strip():
            raise ValueError(f"family_names must be a list of non-empty family names (strings), got {name!r}")
    if len(set(family_names)) != len(family_names):
        raise ValueError(f"family names must be unique, got {family_names}")
    module_map = resolve_modules(modules, family_names)

    # The growth guard, as at the family resolution: duplication compounds, so a run whose rate sits
    # above its loss rate — or a family that drew a high a per-family draw factor — multiplies without bound
    # unless something stops it. A segment may carry several families, and several copies of one, so
    # the run is refused when it would take *any* of them past the quota (see _run_over_cap).
    cap = resolve_max_family_size(max_family_size)

    # Conditioning: a rate carrying ScaledBy reads a driver **per lineage**, so its rate stops being
    # one number for the whole live set and becomes one per lineage. Same machinery as the other two
    # resolutions — each driver resolves once into a DriverTrajectory keyed by the shared species node
    # id, from a file or an in-memory trait result. With no driven rate and no driven extent this is
    # empty and the loop stays exactly the pooled one, so an undriven run is untouched.
    driven = {label: [m for m in r.modifiers if isinstance(m, Driven)]
              for label, r in _rates.items()}
    ext_driven = {label: [m for m in e.modifiers if isinstance(m, Driven)]
                  for label, e in _extents.items()}
    by_key: dict = {}                   # driver key → its Driven (deduped: one driver resolves once)
    for mods in (*driven.values(), *ext_driven.values()):
        for m in mods:
            by_key.setdefault(m.key, m)
    resolved: dict = {}
    if by_key:
        resolved = {key: resolve_driver(m.driver, tree, step=m.step, level="genomes.ordered")
                    for key, m in by_key.items()}
        # a mapping whose states never occur leaves every lineage on the default factor, so the run
        # would secretly be the undriven model — refuse it here, naming the driver
        for mods in (*driven.values(), *ext_driven.values()):
            for m in mods:
                src = m.driver if isinstance(m.driver, str) else f"<{type(m.driver).__name__}>"
                check_mapping_fires(m.mapping, resolved[m.key].states(), driver_label=src)
    # Only a driver on a **rate** makes the loop per-lineage and adds a Gillespie breakpoint. A driver
    # on an **extent** is read at the instant an event fires — it changes how much that event takes,
    # never how often one happens — so it deliberately stays out of `trajs`: no per-lineage rate
    # weights, no extra horizon steps. (SPEC §6.)
    _rate_keys = {m.key for mods in driven.values() for m in mods}
    trajs = {key: traj for key, traj in resolved.items() if key in _rate_keys}
    any_driven = bool(trajs)
    any_ext_driven = any(ext_driven.values())
    # The transfer_to slot is prepared **after** `trajs` is fixed, for the same reason: a driven
    # transfer_to is a weight, not a rate, so its trajectory must not join `trajs` and start adding
    # horizon breakpoints. `resolved` doubles as the driver cache, so a trait that drives both a rate
    # and who receives is loaded once and read from one trajectory.
    group_of, to_traj = prepare_transfer_to(tree, transfer_to, resolved, level="genomes.ordered")

    rng, seed = stream("genomes", seed)     # own stream, and a drawn seed if none was given
    copy_counter = 0
    family_counter = 0
    chrom_counter = 0

    def new_gene(family: int, strand: int) -> Gene:
        nonlocal copy_counter
        g = Gene(copy_counter, family, strand)
        copy_counter += 1
        return g

    # Per-family multipliers, drawn once when a family is created and fixed for its whole life,
    # exactly as at the family resolution: one `Drawn` object read by two rates is one draw for
    # both, two objects are two draws. What differs here is where the weight lands — on the run an
    # event covers, not on the gene it started from (SPEC §6). Empty unless some rate carries one.
    fam_by = {"duplication": tuple(m for m, _ in dup.carried_modifiers(unit="family")),
              "transfer": tuple(m for m, _ in tra.carried_modifiers(unit="family")),
              "loss": tuple(m for m, _ in los.carried_modifiers(unit="family")),
              "inversion": tuple(m for m, _ in inv.carried_modifiers(unit="family")),
              "transposition": tuple(m for m, _ in trp.carried_modifiers(unit="family")),
              "translocation": tuple(m for m, _ in trl.carried_modifiers(unit="family"))}
    any_family = any(fam_by.values())
    fam_mult: dict[str, dict[int, float]] = {key: {} for key in fam_by}

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
    events: list[GeneEdge] = []
    event_positions: list[EventPosition] = []
    rearrangements: list[Inversion | Transposition | Translocation] = []
    chromosome_events: list[ChromosomeEvent] = []

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
    named: dict[str, int] = {}  # a minted id per declared name, dealt round-robin after the anonymous ones
    for j, name in enumerate(family_names):
        fam = new_family()
        named[name] = fam
        chrom = initial_chroms[(initial_families + j) % n_initial_chrom]
        _live(chrom).append(new_gene(fam, +1))
        events.append(GeneEdge(t, "origination", root.id, fam, chrom.genes[-1].id))
        event_positions.append(EventPosition(t, "origination", root.id, chrom.id,
                                             len(chrom.genes) - 1, 1, family=fam))
    # the run's starting genome: a deep snapshot, so the live genome's events never reach it
    initial_genome = tuple(Chromosome(c.id, c.topology, list(c.genes)) for c in initial_chroms)
    enter(alive, gen, pos, root.id, initial_chroms)
    total_copies = initial_families + len(family_names)
    total_chromosomes = n_initial_chrom

    bar = progress_bar(len(schedule), "genomes", unit="branch", enabled=progress)
    si = 0
    while si < len(schedule):
        bar.to(si)
        n = total_copies
        k_alive = len(alive)
        ctx = {"copies": n, "lineages": k_alive, "chromosomes": total_chromosomes, "time": t}
        c = total_chromosomes
        can_xfer = n > 0 and (k_alive >= 2 or self_transfer)
        # A driven rate differs from lineage to lineage, so it is summed **over the living lineages**,
        # each read with its own driver value, its own gene count and its own chromosome count — and
        # the weights are kept, because the affected lineage must then be drawn with them too. The
        # gene count sits inside the weight, which is what makes a driven per-copy rate a two-stage
        # pick (a lineage, then a gene in it) rather than the one-stage lineage draw a per-lineage
        # rate takes. A per-family draw and a Driven cannot both be set, so `w` and `fw` never
        # coexist.
        w: dict[str, list[float]] = {}
        if any_driven:
            drivers = [{key: trajs[key].value(alive[k], t) for key in trajs} for k in range(k_alive)]
            for label, rate in _rates.items():
                if driven[label]:
                    w[label] = [rate.effective(copies=_genome_size(gen[k]), lineages=1,
                                               chromosomes=len(gen[k]), time=t, drivers=drivers[k])
                                for k in range(k_alive)]

        def _r(label, pooled, live=True):
            """The total for one event class: summed per-lineage when driven, pooled when not."""
            if not live:
                return 0.0
            return sum(w[label]) if label in w else pooled

        # A per-copy rate pools over genes, so with per-family weights the total is the unit rate
        # times those weights summed over the live genes — and the run must then be drawn with the
        # same weights, or the rate would say one thing and the picking another. Summed per lineage,
        # so the lineage pick can reuse them. On a circular chromosome ``Σ_s mean_w(s, m)`` is exactly
        # this sum for every run size, which is why no per-size term appears here (SPEC §6).
        fw = None
        if any_family:
            fw = {key: [sum(mult[g.family] for chrom in gen[k] for g in chrom.genes)
                        for k in range(k_alive)]
                  for key, mult in fam_mult.items()}
            one = {"copies": 1, "lineages": 1, "chromosomes": 1, "time": t}
            r_dup = dup.effective(**one) * sum(fw["duplication"]) if n else 0.0
            r_los = los.effective(**one) * sum(fw["loss"]) if n else 0.0
            r_tra = tra.effective(**one) * sum(fw["transfer"]) if can_xfer else 0.0
            r_inv = inv.effective(**one) * sum(fw["inversion"]) if n else 0.0
            r_trp = trp.effective(**one) * sum(fw["transposition"]) if n else 0.0
            r_trl = trl.effective(**one) * sum(fw["translocation"]) if n else 0.0
        else:
            r_dup = _r("duplication", dup.effective(**ctx) if n else 0.0, live=bool(n))
            r_los = _r("loss", los.effective(**ctx) if n else 0.0, live=bool(n))
            r_tra = _r("transfer", tra.effective(**ctx) if can_xfer else 0.0, live=can_xfer)
            r_inv = _r("inversion", inv.effective(**ctx) if n else 0.0, live=bool(n))  # per copy
            r_trp = _r("transposition", trp.effective(**ctx) if n else 0.0, live=bool(n))  # per copy
            r_trl = _r("translocation", trl.effective(**ctx) if n else 0.0, live=bool(n))  # per copy
        r_org = _r("origination", org.effective(**ctx))                 # per lineage
        r_fis = _r("fission", fis.effective(**ctx) if c else 0.0, live=bool(c))  # per chromosome
        r_fus = _r("fusion", fus.effective(**ctx) if c else 0.0, live=bool(c))
        r_cor = _r("chromosome_origination", cor.effective(**ctx))      # per lineage (de-novo replicon)
        r_clo = _r("chromosome_loss", clo.effective(**ctx) if c else 0.0, live=bool(c))
        total = (r_dup + r_los + r_org + r_tra + r_inv + r_trp + r_trl
                 + r_fis + r_fus + r_cor + r_clo)

        next_species = schedule[si][0]
        horizon = min(next_species, dup.next_change(t), los.next_change(t), org.next_change(t),
                      tra.next_change(t), inv.next_change(t), trp.next_change(t), trl.next_change(t),
                      fis.next_change(t), fus.next_change(t), cor.next_change(t), clo.next_change(t))
        if any_driven:  # a driven rate also changes when its driver switches mid-branch — step there
            horizon = min(horizon, min((trajs[key].next_change(alive[k], t) for key in trajs
                                        for k in range(k_alive)), default=math.inf))

        def _ext_ctx(k):
            """The context an extent is sampled in, on the lineage the event landed on.

            It cannot be built before the lineage is drawn, because a driven extent is read on the
            **acting** lineage at the instant the event fires — which is also why an extent adds no
            Gillespie breakpoint and never enters the horizon above (SPEC §6). With no driven extent
            this is the same context the rates were read in.

            The rest of `ctx` — the gene, lineage and chromosome counts — goes with it, because
            `Modifier.implemented_for` promises this engine supplies them and a modifier of your own
            is admitted onto an extent by the same gate that admits it onto a rate. Handing an
            extent a thinner context meant one gate certifying two different contracts: a modifier
            written the documented way read zeros, and one with a required keyword died mid-run."""
            # `ctx` was snapshotted at the top of the loop, before `t` advanced to the firing
            # instant, so `time` has to be taken fresh: an extent's own breakpoints are kept out of
            # the horizon, so a schedule's breakpoint routinely falls inside a stretch, and reading
            # the stale `t` would size the event on the wrong side of it.
            if not any_ext_driven:
                return {**ctx, "time": t}
            return {**ctx, "time": t,
                    "drivers": {key: resolved[key].value(alive[k], t) for key in resolved}}

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
                    picked = _pick_event_run(rng, gen, n, fw, fam_mult, "duplication", dup_ext,
                                             _ext_ctx, w.get("duplication"))
                    if picked is not None:
                        k, ci, j, m = picked
                        if not _run_over_cap(gen[k], gen[k][ci], j, m, cap):
                            total_copies += _duplicate(gen[k][ci], j, m, tree.nodes[alive[k]], t,
                                                       events, event_positions, new_gene)
                elif r < b_los:
                    picked = _pick_event_run(rng, gen, n, fw, fam_mult, "loss", los_ext,
                                             _ext_ctx, w.get("loss"))
                    if picked is not None:
                        k, ci, j, m = picked
                        total_copies -= _lose_at(gen[k][ci], j, m, tree.nodes[alive[k]], t, events,
                                                 event_positions)
                elif r < b_org:
                    # origination is per lineage: a uniform lineage, or one drawn by its own rate
                    # when that rate is driven (the same weights the total was summed with)
                    k = (weighted_index(rng, w["origination"], r_org) if "origination" in w
                         else int(rng.integers(k_alive)))
                    _originate(gen[k], tree.nodes[alive[k]], t, events, event_positions, new_gene,
                               new_family, rng)
                    total_copies += 1
                elif r < b_tra:
                    picked = _pick_event_run(rng, gen, n, fw, fam_mult, "transfer", tra_ext,
                                             _ext_ctx, w.get("transfer"))
                    if picked is not None:                # driven: the weighted lineage is the DONOR
                        kd, cdi, jd, m = picked
                        total_copies += _do_transfer(rng, tree, alive, gen, kd, cdi, jd, m, t, events,
                                                     event_positions, new_gene, transfer_to,
                                                     replacement, self_transfer, depth, cap,
                                                     to_traj, group_of)
                elif r < b_inv:
                    picked = _pick_event_run(rng, gen, n, fw, fam_mult, "inversion", inv_ext,
                                             _ext_ctx, w.get("inversion"))
                    if picked is not None:                # the run starts at a gene, so: per copy
                        k, ci, i0, m = picked
                        _invert(gen[k][ci], i0, m, tree.nodes[alive[k]], t, rearrangements)
                elif r < b_trp:
                    picked = _pick_event_run(rng, gen, n, fw, fam_mult, "transposition", trp_ext,
                                             _ext_ctx, w.get("transposition"))
                    if picked is not None:
                        k, ci, i0, m = picked
                        _transpose(gen[k][ci], i0, m, tree.nodes[alive[k]], t, rearrangements, rng,
                                   inversion_probability)
                elif r < b_trl:
                    picked = _pick_event_run(rng, gen, n, fw, fam_mult, "translocation", trl_ext,
                                             _ext_ctx, w.get("translocation"))
                    if picked is not None:
                        k, ci, j, m = picked
                        _translocate(gen[k], ci, j, m, tree.nodes[alive[k]], t, rearrangements, rng,
                                     inversion_probability)
                elif r < b_fis:
                    picked = _pick_chromosome(rng, gen, c, w.get("fission"))
                    if picked is not None:
                        k, ci = picked
                        dc, dg = _fission(gen[k], ci, tree.nodes[alive[k]], t, chromosome_events,
                                          new_chromosome, rng)
                        total_chromosomes += dc
                        total_copies += dg
                elif r < b_fus:
                    picked = _pick_chromosome(rng, gen, c, w.get("fusion"))
                    if picked is not None:
                        k, ci = picked
                        dc, dg = _fusion(gen[k], ci, tree.nodes[alive[k]], t, chromosome_events,
                                         new_chromosome, rng)
                        total_chromosomes += dc
                        total_copies += dg
                elif r < b_cor:
                    # chromosome origination is per lineage, uniform or driven, exactly as origination
                    k = (weighted_index(rng, w["chromosome_origination"], r_cor)
                         if "chromosome_origination" in w else int(rng.integers(k_alive)))
                    dc, dg = _chromosome_originate(gen[k], tree.nodes[alive[k]], t, chromosome_events,
                                                   new_chromosome)
                    total_chromosomes += dc
                    total_copies += dg
                else:
                    picked = _pick_chromosome(rng, gen, c, w.get("chromosome_loss"))
                    if picked is not None:
                        k, ci = picked
                        dc, dg = _chromosome_lose(gen[k], ci, tree.nodes[alive[k]], t, events,
                                                  event_positions, chromosome_events)
                        total_chromosomes += dc
                        total_copies += dg
                continue

        if horizon == next_species:  # advance to the tree's next event(s); process the whole tie-batch
            t = next_species
            while si < len(schedule) and schedule[si][0] == t:
                i = schedule[si][1]
                g = gen[pos[i]]
                genomes[i] = tuple(Chromosome(c.id, c.topology, tuple(c.genes)) for c in g)  # freeze
                total_copies -= sum(len(c.genes) for c in g)
                total_chromosomes -= len(g)
                retire(alive, gen, pos, pos[i])
                node = tree.nodes[i]
                if node.children is not None:  # a speciation: re-mint every chromosome and gene id
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
                        total_copies += sum(len(ch.genes) for ch in cg)
                        total_chromosomes += len(cg)
                si += 1
        else:
            t = horizon  # a skyline breakpoint: advance and re-evaluate the (now changed) rate

    bar.close()
    return OrderedGenomesResult(tree, genomes, events, rearrangements, chromosome_events, seed,
                                named, module_map, event_positions, initial_genome)


__all__ = ["simulate_genomes_ordered", "OrderedGenomesResult", "Gene", "Chromosome",
           "ChromosomeEvent", "Inversion", "Transposition", "Translocation"]
