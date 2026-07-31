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
`Event` log (position-blind, so ``gene_trees`` and ``profiles`` are
derived from it unchanged), and the same live-lineage bookkeeping. What differs is the state (a list
of chromosomes) and the segmental, position-aware mutators, plus the ``rearrangements`` and
``chromosome_events`` logs. The nucleotide resolution (genes/intergenes, indels) is
`simulate_genomes_nucleotide()`.
"""

from __future__ import annotations

import collections
import pathlib
from dataclasses import dataclass, field
from functools import cached_property

import numpy as np

from ..rates.extent import as_extent
from ..rates.modifiers import ByFamily, OnTime
from ..rates.rate import as_rate
from ..rates.scope import PerChromosome, PerCopy, PerLineage
from ..tree import Tree, as_tree
from .chromosomes import ChromosomeEvent, chromosome_events_tsv, rearrangement_events_tsv
from .family import resolve_max_family_size
from ._live import enter, retire, weighted_index, without_cyclic_gc
from ._transfer import Distance, mean_root_to_tip, recipient_index
from .._runtime.outputs import grouped_dir
from .._runtime.progress import progress_bar
from .events import _COLS, Event, _branches, _grouped, _name, event_rows, gene_label
from .gene_trees import GeneTree, gene_trees_from_events, write_gene_trees
from .profiles import Profiles, profiles_from_genomes


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
    where it may be re-anchored freely (see `_anchor()`). Topology does not yet gate which
    fissions and fusions are legal."""

    id: int
    topology: str
    genes: list[Gene]


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
    it is logged here, separate from the gene-genealogy `Event`
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
    `Event`.

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
    events: list[Event]
    rearrangements: list[Inversion | Transposition | Translocation]
    chromosome_events: list[ChromosomeEvent]
    seed: int | None
    #: ``{name: family id}`` for families declared by ``family_names=[…]`` — the handle to a *named* family.
    family_names: dict[str, int] = field(default_factory=dict)
    #: where each gene-genealogy `Event` happened — the positional
    #: companion to `events`, which is position-blind. See `EventPosition`.
    event_positions: list[EventPosition] = field(default_factory=list)
    #: The genome the run **started** with, at the root lineage's origination — before any event.
    #: It is not in `genomes`, which holds a genome per *node*, and a node sits at the **end**
    #: of its branch: the root branch is real simulated time, so ``genomes[root]`` is this genome plus
    #: whatever happened along the stem.
    initial_genome: tuple[Chromosome, ...] = ()

    def __repr__(self) -> str:
        return (f"OrderedGenomesResult({len(self.complete_tree.extant_leaves())} extant genomes, "
                f"{len(self.genomes)} nodes, {len(self.events)} events, "
                f"{len(self.rearrangements)} rearrangements, seed={self.seed})")

    def family_counts(self, node_id: int) -> collections.Counter:
        """A multiset view of one node's genome: ``family id → copy count`` (across all chromosomes)."""
        return collections.Counter(g.family for chrom in self.genomes[node_id] for g in chrom.genes)

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
    def gene_trees(self) -> dict[int, GeneTree]:
        """``{family id: GeneTree}`` — each family's true genealogy inside the complete tree, derived
        from the (position-blind) event log exactly as for the family core. See `gene_trees`."""
        return gene_trees_from_events(self.events, self.complete_tree)

    #: Every token ``write()`` honours — the write vocabulary, declared rather than left
    #: implicit in the method body. The CLI builds ``--write``'s choices from this, so the two
    #: cannot drift: they did, and `initial_sequence` and `species_tree` were writable from
    #: Python and unnameable on the command line.
    OUTPUTS = ("events", "profiles", "gene_order", "initial_genome",
               "chromosome_events", "gene_trees", "species_tree")

    def write(self, directory, outputs=("events", "profiles", "gene_order", "initial_genome",
                                        "gene_trees", "chromosome_events", "species_tree"), *,
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
        d = pathlib.Path(directory)
        d.mkdir(parents=True, exist_ok=True)
        names = self.complete_tree.labels()   # e<id> for a lineage that died; n<id> for the rest
        if "events" in outputs:
            (d / "genome_events.tsv").write_text(
                _events_tsv(self.events, self.event_positions, names), encoding="utf-8")
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

    def _initial_genome_tsv(self) -> str:
        """The layout the run started with — ``gene_order.tsv``'s columns without ``lineage``, which
        is the whole point: it belongs to the start of the root branch, not to a node."""
        cols = ("chromosome", "position", "strand", "family", "copy")
        rows = [f"{chrom.id}\t{pos}\t{g.strand}\t{g.family}\t{gene_label(g.id)}"
                for chrom in self.initial_genome for pos, g in enumerate(chrom.genes)]
        return "\n".join(["\t".join(cols), *rows]) + "\n"

    def _gene_order_tsv(self, names=None) -> str:
        cols = ("lineage", "chromosome", "position", "strand", "family", "copy")
        rows = [f"{_name(names, s)}\t{ch}\t{p}\t{st}\t{fam}\t{gene_label(gid)}"
                for s in sorted(self.genomes)
                for (ch, p, st, fam, gid) in self.gene_order(s)]
        return "\n".join(["\t".join(cols), *rows]) + "\n"


#: ``genome_events.tsv`` here: the shared genealogy columns (`_COLS` — one row per event, its
#: participants written ``n<species>_g<copy>``) with **where** the event happened beside them. The
#: coordinates are the one thing this resolution has that the family core does not, so they are the
#: one thing it adds; the genealogy half is written by `event_rows()`, not repeated here, because
#: `events_from_tsv()` reads this table by requiring `_COLS` as a literal **prefix** of the header
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
    one table, so this walks the same `_grouped()` the genealogy writer does.

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
    for time, kind, family, parents, children in _grouped(events):
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

def _pick_gene(rng, gen, total_copies) -> tuple[int, int, int]:
    """A uniform global gene pick → ``(lineage k, chromosome index ci in gen[k], position j)``.
    Realises per-copy scope across the whole pool: every gene, in any chromosome of any lineage, is
    equally likely."""
    m = int(rng.integers(total_copies))
    for k, genome in enumerate(gen):
        for ci, chrom in enumerate(genome):
            if m < len(chrom.genes):
                return k, ci, m
            m -= len(chrom.genes)
    raise AssertionError("total_copies out of sync with the genomes")  # unreachable


def _pick_chromosome(rng, gen, total_chromosomes) -> tuple[int, int]:
    """A uniform global chromosome pick → ``(lineage k, chromosome index ci in gen[k])``. Realises
    per-chromosome scope across the whole pool."""
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


def _pick_event_run(rng, gen, n, fw, fam_mult, key, ext, ctx=None):
    """``(lineage, chromosome index, start, run size)`` for one gene-level event.

    Uniform over genes when no per-family weight is set — the plain path, untouched. With one set, the
    lineage is drawn by its summed weight and the run by `_pick_run_by_family()`, so the weight
    reaches the segment. ``None`` when the drawn lineage has nothing left to act on."""
    if fw is None:
        k, ci, j = _pick_gene(rng, gen, n)
        return k, ci, j, _extent(rng, ext, gen[k][ci], j, ctx)
    w = fw[key]
    total = sum(w)
    if total <= 0.0:
        return None
    k = weighted_index(rng, w, total)
    picked = _pick_run_by_family(rng, gen[k], fam_mult[key], ext, ctx)
    if picked is None:
        return None
    ci, j, m = picked
    return k, ci, j, m


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
    chrom.genes.insert(at, g)
    events.append(Event(t, "origination", node.id, fam, g.id))
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
        events.append(Event(t, "duplication", node.id, old.family, cont.id, parent=old.id))
        events.append(Event(t, "duplication", node.id, old.family, cp.id, parent=old.id))
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
        events.append(Event(t, "loss", node.id, g.family, g.id))
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
                 transfer_to, replacement, self_transfer, depth, cap=None) -> int:
    """The segment ``[jd, jd+m)`` on the donor's chromosome ``cdi`` transfers to a contemporaneous
    recipient: each gene ends → a continuation on the donor branch and a transferred copy on the
    recipient (a horizontal gene-tree edge). The run may wrap position 0 on a circular donor
    chromosome. The transferred copies arrive as a block at a random position on a uniformly-chosen
    recipient chromosome (strands travel with them). Returns the change in total gene count: ``+m``
    additive, minus one per homologous copy displaced under ``replacement``."""
    donor = alive[kd]
    jd = _anchor(gen[kd][cdi], jd, m)
    segment = gen[kd][cdi].genes[jd:jd + m]
    if transfer_to == "uniform":
        # O(1) uniform recipient — the same single draw as recipient_index's
        # cand[rng.integers(len(cand))] over every alive lineage but the donor; the donor-skip is a
        # +1 index shift, so no O(alive) candidate list is built per transfer (see family._do_transfer).
        npool = len(alive) if self_transfer else len(alive) - 1
        if npool <= 0:
            return 0
        i = int(rng.integers(npool))
        kr = i if (self_transfer or i < kd) else i + 1
    else:  # a Distance weighting must weigh every candidate — inherently O(alive)
        cand = [k for k in range(len(alive)) if self_transfer or k != kd]
        kr = recipient_index(rng, tree, alive, cand, donor, t, transfer_to, depth)
    recipient = alive[kr]
    rgenome = gen[kr]
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
            events.append(Event(t, "loss", recipient, old.family, replaced))
        events.append(Event(t, "transfer", donor, old.family, cont.id, parent=old.id, donor=donor,
                            replaced=replaced))
        events.append(Event(t, "transfer", recipient, old.family, xf.id, parent=old.id,
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
    """Chromosome ``ci`` merges with another chromosome of the same genome — the **reticulation**
    (two parents, one child): the fused child re-mints, both parents end. No-op if the genome has
    fewer than two chromosomes."""
    if len(genome) < 2:
        return (0, 0)
    cj = int(rng.integers(len(genome) - 1))
    if cj >= ci:
        cj += 1                                        # a uniform chromosome index distinct from ci
    a, b = genome[ci], genome[cj]
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
    gene on it ends as a gene ``loss``. No-op if it is the genome's last chromosome (a lineage never
    loses its entire genome this way)."""
    if len(genome) < 2:
        return (0, 0)
    lost = genome[ci]
    for g in lost.genes:
        events.append(Event(t, "loss", node.id, g.family, g.id))
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
                             initial_families=100, family_names=None, family_speed=None,
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
    "toxin")``), as in the family core; ``transfer_to`` / ``replacement`` / ``self_transfer`` behave
    as in the family core.

    The **chromosome tier** changes chromosome *number*: ``fission`` (split), ``fusion`` (merge — the
    reticulation), ``chromosome_origination`` (a de-novo replicon), ``chromosome_loss`` (a whole
    chromosome and its genes die; never the genome's last). Chromosomes carry identity — re-minted at
    every event that reshapes them — so ``chromosome_events`` is the true reticulating chromosome
    genealogy, rooted at the initial and de-novo originations. Deterministic given ``seed``.
    """
    tree = as_tree(tree, level="genomes")
    labels = _topologies(chromosomes, topology)
    n_initial_chrom = chromosomes
    dup = as_rate(duplication, default_scope=PerCopy)
    tra = as_rate(transfer, default_scope=PerCopy)
    los = as_rate(loss, default_scope=PerCopy)
    org = as_rate(origination, default_scope=PerLineage)
    inv = as_rate(inversion, default_scope=PerCopy)
    trp = as_rate(transposition, default_scope=PerCopy)
    trl = as_rate(translocation, default_scope=PerCopy)
    fis = as_rate(fission, default_scope=PerChromosome)
    fus = as_rate(fusion, default_scope=PerChromosome)
    cor = as_rate(chromosome_origination, default_scope=PerLineage)
    clo = as_rate(chromosome_loss, default_scope=PerChromosome)
    # this slice implements each event's default scope, OnTime (skyline), and ByFamily — the last with
    # the weight on the SEGMENT rather than on its starting gene (SPEC §6, and _pick_run_by_family). A
    # scope override or a clade/driven modifier is a later slice, so reject it rather than silently
    # mis-scale (see the family engine for the reasoning).
    for label, rate, want in (("duplication", dup, PerCopy), ("transfer", tra, PerCopy),
                              ("loss", los, PerCopy), ("origination", org, PerLineage),
                              ("inversion", inv, PerCopy), ("transposition", trp, PerCopy),
                              ("translocation", trl, PerCopy), ("fission", fis, PerChromosome),
                              ("fusion", fus, PerChromosome), ("chromosome_loss", clo, PerChromosome),
                              ("chromosome_origination", cor, PerLineage)):
        if not isinstance(rate.scope, want):
            raise ValueError(
                f"{label} has a {type(rate.scope).__name__} scope, but the ordered genome engine "
                f"takes only {want.__name__} for {label} this slice — scope overrides are a later slice."
            )
        for m in rate.modifiers:
            if isinstance(m, ByFamily) and label == "origination":
                raise ValueError(
                    "origination carries ByFamily, but origination is the rate at which families are "
                    "CREATED — when it is read there is no family yet to have drawn a factor for. "
                    "Put ByFamily on duplication, transfer, loss, inversion, transposition or "
                    "translocation, or use family_speed= for a family-wide tempo.")
            if isinstance(m, ByFamily) and not isinstance(rate.scope, PerCopy):
                raise ValueError(
                    f"{label} carries ByFamily on a {type(rate.scope).__name__} scope. A per-family "
                    f"weight has to reach the genes an event covers, so it applies to the per-copy "
                    f"gene events only — not to the chromosome tier, which acts on whole replicons.")
            if not isinstance(m, (OnTime, ByFamily)):
                raise ValueError(
                    f"{label} carries {type(m).__name__}, which the ordered genome engine does not "
                    f"support. It takes OnTime (skyline) and ByFamily. Clade drift and driven rates "
                    f"are not implemented yet at this resolution."
                )
    # per-event extent distributions (segment size in genes); a bare number is the mean, None a single gene
    def _ext_spec(spec, label):
        """One event's extent (SPEC §6): ``base × modifiers``, no scope, in **genes** here. An extent
        takes the modifiers this resolution supports on a rate — ``OnTime`` — and they scale the size
        drawn. A driver is not among them: this engine takes no ``DrivenBy`` on its rates either, and
        for trait-driven rearrangement the nucleotide resolution is the one that has it."""
        e = as_extent(spec)
        for m in e.modifiers:
            if not isinstance(m, OnTime):
                raise ValueError(
                    f"{label} carries {type(m).__name__}, which the ordered genome engine does not "
                    f"support on an extent — it takes only OnTime (a skyline in time). For a "
                    f"trait-driven extent use --resolution nucleotide.")
        return e

    dup_ext, los_ext, tra_ext = (_ext_spec(duplication_extent, "duplication_extent"),
                                 _ext_spec(loss_extent, "loss_extent"),
                                 _ext_spec(transfer_extent, "transfer_extent"))
    inv_ext, trp_ext, trl_ext = (_ext_spec(inversion_extent, "inversion_extent"),
                                 _ext_spec(transposition_extent, "transposition_extent"),
                                 _ext_spec(translocation_extent, "translocation_extent"))
    if not 0.0 <= inversion_probability <= 1.0:
        raise ValueError(f"inversion_probability must be in [0, 1], got {inversion_probability!r}")
    if transfer_to == "distance":
        transfer_to = Distance()
    if transfer_to != "uniform" and not isinstance(transfer_to, Distance):
        raise ValueError(f"transfer_to must be 'uniform', 'distance', or Distance(decay=), got {transfer_to!r}")
    if isinstance(initial_families, bool) or not isinstance(initial_families, int) or initial_families < 0:
        raise ValueError(f"initial_families must be a non-negative integer, got {initial_families!r}")
    family_names = list(family_names) if family_names is not None else []
    for name in family_names:
        if not isinstance(name, str) or not name.strip():
            raise ValueError(f"family_names must be a list of non-empty family names (strings), got {name!r}")
    if len(set(family_names)) != len(family_names):
        raise ValueError(f"family names must be unique, got {family_names}")

    # The growth guard, as at the family resolution: duplication compounds, so a run whose rate sits
    # above its loss rate — or a family that drew a high ByFamily factor — multiplies without bound
    # unless something stops it. A segment may carry several families, and several copies of one, so
    # the run is refused when it would take *any* of them past the quota (see _run_over_cap).
    cap = resolve_max_family_size(max_family_size)
    if family_speed is not None and not isinstance(family_speed, ByFamily):
        raise ValueError(
            f"family_speed must be a ByFamily(...) draw, got {type(family_speed).__name__}.")

    rng = np.random.default_rng(seed)
    copy_counter = 0
    family_counter = 0
    chrom_counter = 0

    def new_gene(family: int, strand: int) -> Gene:
        nonlocal copy_counter
        g = Gene(copy_counter, family, strand)
        copy_counter += 1
        return g

    # Per-family multipliers, drawn once when a family is minted and fixed for its whole life, exactly
    # as at the family resolution: family_speed scales every rate that family has (one draw), a
    # ByFamily on a single rate varies that rate alone (its own draw). What differs here is where the
    # weight lands — on the run an event covers, not on the gene it started from (SPEC §6). Empty
    # unless one of them is used, and a run without either is byte-identical to the plain path.
    fam_by = {"duplication": next((m for m in dup.modifiers if isinstance(m, ByFamily)), None),
              "transfer": next((m for m in tra.modifiers if isinstance(m, ByFamily)), None),
              "loss": next((m for m in los.modifiers if isinstance(m, ByFamily)), None),
              "inversion": next((m for m in inv.modifiers if isinstance(m, ByFamily)), None),
              "transposition": next((m for m in trp.modifiers if isinstance(m, ByFamily)), None),
              "translocation": next((m for m in trl.modifiers if isinstance(m, ByFamily)), None)}
    any_family = family_speed is not None or any(fam_by.values())
    fam_mult: dict[str, dict[int, float]] = {key: {} for key in fam_by}

    def new_family() -> int:
        nonlocal family_counter
        f = family_counter
        family_counter += 1
        if any_family:
            speed = family_speed.draw(rng) if family_speed is not None else 1.0
            for key, m in fam_by.items():
                fam_mult[key][f] = speed * (m.draw(rng) if m is not None else 1.0)
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
    events: list[Event] = []
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
        chrom.genes.append(new_gene(fam, +1))
        events.append(Event(t, "origination", root.id, fam, chrom.genes[-1].id))
        event_positions.append(EventPosition(t, "origination", root.id, chrom.id,
                                             len(chrom.genes) - 1, 1, family=fam))
    named: dict[str, int] = {}  # a minted id per declared name, dealt round-robin after the anonymous ones
    for j, name in enumerate(family_names):
        fam = new_family()
        named[name] = fam
        chrom = initial_chroms[(initial_families + j) % n_initial_chrom]
        chrom.genes.append(new_gene(fam, +1))
        events.append(Event(t, "origination", root.id, fam, chrom.genes[-1].id))
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
            r_dup = dup.effective(**ctx) if n else 0.0
            r_los = los.effective(**ctx) if n else 0.0
            r_tra = tra.effective(**ctx) if can_xfer else 0.0
            r_inv = inv.effective(**ctx) if n else 0.0                  # per copy (the run's start)
            r_trp = trp.effective(**ctx) if n else 0.0                  # per copy (the run's start)
            r_trl = trl.effective(**ctx) if n else 0.0                  # per copy; needs >=2 chromosomes
        r_org = org.effective(**ctx)                                    # per lineage
        r_fis = fis.effective(**ctx) if c else 0.0                      # per chromosome (the tier)
        r_fus = fus.effective(**ctx) if c else 0.0
        r_cor = cor.effective(**ctx)                                    # per lineage (de-novo replicon)
        r_clo = clo.effective(**ctx) if c else 0.0
        total = (r_dup + r_los + r_org + r_tra + r_inv + r_trp + r_trl
                 + r_fis + r_fus + r_cor + r_clo)

        next_species = schedule[si][0]
        horizon = min(next_species, dup.next_change(t), los.next_change(t), org.next_change(t),
                      tra.next_change(t), inv.next_change(t), trp.next_change(t), trl.next_change(t),
                      fis.next_change(t), fus.next_change(t), cor.next_change(t), clo.next_change(t))

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
                    picked = _pick_event_run(rng, gen, n, fw, fam_mult, "duplication", dup_ext, {"time": t})
                    if picked is not None:
                        k, ci, j, m = picked
                        if not _run_over_cap(gen[k], gen[k][ci], j, m, cap):
                            total_copies += _duplicate(gen[k][ci], j, m, tree.nodes[alive[k]], t,
                                                       events, event_positions, new_gene)
                elif r < b_los:
                    picked = _pick_event_run(rng, gen, n, fw, fam_mult, "loss", los_ext, {"time": t})
                    if picked is not None:
                        k, ci, j, m = picked
                        total_copies -= _lose_at(gen[k][ci], j, m, tree.nodes[alive[k]], t, events,
                                                 event_positions)
                elif r < b_org:
                    k = int(rng.integers(k_alive))  # origination is per lineage: a uniform lineage
                    _originate(gen[k], tree.nodes[alive[k]], t, events, event_positions, new_gene,
                               new_family, rng)
                    total_copies += 1
                elif r < b_tra:
                    picked = _pick_event_run(rng, gen, n, fw, fam_mult, "transfer", tra_ext, {"time": t})
                    if picked is not None:
                        kd, cdi, jd, m = picked
                        total_copies += _do_transfer(rng, tree, alive, gen, kd, cdi, jd, m, t, events,
                                                     event_positions, new_gene, transfer_to,
                                                     replacement, self_transfer, depth, cap)
                elif r < b_inv:
                    picked = _pick_event_run(rng, gen, n, fw, fam_mult, "inversion", inv_ext, {"time": t})
                    if picked is not None:                # the run starts at a gene, so: per copy
                        k, ci, i0, m = picked
                        _invert(gen[k][ci], i0, m, tree.nodes[alive[k]], t, rearrangements)
                elif r < b_trp:
                    picked = _pick_event_run(rng, gen, n, fw, fam_mult, "transposition", trp_ext, {"time": t})
                    if picked is not None:
                        k, ci, i0, m = picked
                        _transpose(gen[k][ci], i0, m, tree.nodes[alive[k]], t, rearrangements, rng,
                                   inversion_probability)
                elif r < b_trl:
                    picked = _pick_event_run(rng, gen, n, fw, fam_mult, "translocation", trl_ext, {"time": t})
                    if picked is not None:
                        k, ci, j, m = picked
                        _translocate(gen[k], ci, j, m, tree.nodes[alive[k]], t, rearrangements, rng,
                                     inversion_probability)
                elif r < b_fis:
                    k, ci = _pick_chromosome(rng, gen, c)
                    dc, dg = _fission(gen[k], ci, tree.nodes[alive[k]], t, chromosome_events,
                                      new_chromosome, rng)
                    total_chromosomes += dc
                    total_copies += dg
                elif r < b_fus:
                    k, ci = _pick_chromosome(rng, gen, c)
                    dc, dg = _fusion(gen[k], ci, tree.nodes[alive[k]], t, chromosome_events,
                                     new_chromosome, rng)
                    total_chromosomes += dc
                    total_copies += dg
                elif r < b_cor:
                    k = int(rng.integers(k_alive))  # chromosome origination is per lineage
                    dc, dg = _chromosome_originate(gen[k], tree.nodes[alive[k]], t, chromosome_events,
                                                   new_chromosome)
                    total_chromosomes += dc
                    total_copies += dg
                else:
                    k, ci = _pick_chromosome(rng, gen, c)
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
                    child_genomes = {c: [] for c in node.children}
                    for pchrom in g:
                        dcids = []
                        per_daughter: list[list[Event]] = []
                        for c in node.children:
                            dcid = new_chromosome()
                            dcids.append(dcid)
                            dgenes, edges = [], []
                            for old in pchrom.genes:  # ZOMBI1: the gene ends and continues, fresh id
                                ng = new_gene(old.family, old.strand)
                                dgenes.append(ng)
                                edges.append(Event(t, "speciation", c, old.family, ng.id,
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
                                named, event_positions, initial_genome)


__all__ = ["simulate_genomes_ordered", "OrderedGenomesResult", "Gene", "Chromosome",
           "ChromosomeEvent", "Inversion", "Transposition", "Translocation"]
