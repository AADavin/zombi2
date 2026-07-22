"""Genomes II — ordered: genes carry a position and an orientation, on chromosomes.

The ordered resolution layers **position** over the unordered D/T/L/O core (Chapter 6). A genome is
no longer a multiset of gene copies but a list of **chromosomes**, each an ordered run of oriented
:class:`Gene`\\ s.

**Every gene-level event acts on an extension** — a run of consecutive genes (the ZOMBI1 model), its
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
the seed and de-novo originations, recorded as an edge list (its ground truth — a network is a graph,
not eNewick; see ``chromosome-network.md``).

It is the genome twin of the unordered core and shares its spine: one forward Gillespie over the
**complete** species tree, the same ``scope(base) × modifiers`` rate grammar, the same gene-genealogy
:class:`~zombi2.genomes.events.Event` log (position-blind, so ``gene_trees`` and ``profiles`` are
derived from it unchanged), and the same live-lineage bookkeeping. What differs is the state (a list
of chromosomes) and the segmental, position-aware mutators, plus the ``rearrangements`` and
``chromosome_events`` logs. Still to come: the nucleotide resolution (genes/intergenes, indels).
"""

from __future__ import annotations

import collections
import pathlib
from dataclasses import dataclass, field
from functools import cached_property

import numpy as np

from ..rates.distributions import Geometric, as_distribution
from ..rates.modifiers import OnTime
from ..rates.rate import as_rate
from ..rates.scope import PerChromosome, PerCopy, PerLineage
from ..species import SpeciesResult, Tree
from .chromosomes import ChromosomeEvent, chromosome_events_tsv
from ._live import enter, retire
from ._transfer import Distance, mean_root_to_tip, recipient_index
from ..progress import progress_bar
from .events import Event, events_tsv, node_label
from .gene_trees import GeneTree, gene_trees_from_events, write_gene_trees
from .profiles import Profiles, profiles_from_genomes


@dataclass(frozen=True)
class Gene:
    """One gene copy with an **orientation**: a member of family ``family``, identified by a
    globally-unique ``id`` (per segment, the ZOMBI1 model), lying on its chromosome on the ``strand``
    ``+1`` or ``-1``. It is the unordered :class:`~zombi2.genomes.GeneCopy` with the one thing that
    only makes sense once genes are ordered — which way it points. Its position is implicit: the index
    of the gene in its chromosome's ordered list. Birth/death and parentage live in the event log."""

    id: int
    family: int
    strand: int  # +1 / -1


@dataclass
class Chromosome:
    """One chromosome: an ordered run of :class:`Gene`\\ s, identified by ``id`` (re-minted at every
    speciation, so it names a chromosome *lineage*), with a ``topology`` — ``"circular"`` or
    ``"linear"``.

    Topology decides where a segmental event's run stops. A **circular** chromosome has no ends, so a
    run that reaches the last gene continues from the first — it wraps position 0 — and is limited
    only by the whole chromosome. A **linear** one has ends, so a run stops at the last gene. Position
    0 is therefore a real boundary on a linear chromosome and pure bookkeeping on a circular one,
    where it may be re-anchored freely (see :func:`_anchor`). Topology does not yet gate which
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
    it is logged here, separate from the gene-genealogy :class:`~zombi2.genomes.events.Event`
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
    :class:`~zombi2.genomes.events.Event`.

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

    A transfer spans two branches, so it writes **two rows** — one on each — and both name the whole
    edge in ``donor`` and ``recipient``. Pair them on ``(time, donor, recipient)``; the donor row is
    written first. (This follows Krister Swenson's fork, which splits a transfer into a leaving and
    an arriving event, except that the branches are named outright rather than matched by timestamp.)

    Together with the genomes (``gene_order``) and the rearrangement log this is **sufficient to
    replay a run**: no join back to the genealogy is needed. A join is still possible — on
    ``(time, lineage, kind)`` — but is not one-to-one, since the crown's founding originations all
    fire at the same instant on the same lineage.

    Rows sharing a ``time`` apply **in the order written** — a replacing transfer displaces its
    homologs (each a ``loss`` on the recipient) before the arriving block is inserted."""

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
    """What :func:`simulate_genomes_ordered` returns: the ``complete_tree`` it ran on, the final
    ``genomes`` at **every** node as tuples of :class:`Chromosome`\\ s, the shared gene-genealogy
    ``events`` log, the ``rearrangements`` (inversions) and ``chromosome_events`` (the chromosome
    genealogy) logs, and the ``seed``. The observed genomes are the extant tips; ``profiles`` and
    ``gene_trees`` are derived from the (position-blind) genealogy exactly as for the unordered core;
    ``gene_order`` reads a node's layout, and ``write`` materialises the chosen outputs."""

    complete_tree: Tree
    genomes: dict[int, tuple[Chromosome, ...]]
    events: list[Event]
    rearrangements: list[Inversion | Transposition | Translocation]
    chromosome_events: list[ChromosomeEvent]
    seed: int | None
    #: ``{name: family id}`` for families seeded by ``families=[…]`` — the handle to a *named* family.
    family_names: dict[str, int] = field(default_factory=dict)
    #: where each gene-genealogy :class:`~zombi2.genomes.events.Event` happened — the positional
    #: companion to :attr:`events`, which is position-blind. See :class:`EventPosition`.
    event_positions: list[EventPosition] = field(default_factory=list)

    def family_counts(self, node_id: int) -> collections.Counter:
        """A multiset view of one node's genome: ``family id → copy count`` (across all chromosomes)."""
        return collections.Counter(g.family for chrom in self.genomes[node_id] for g in chrom.genes)

    def has_family(self, node_id: int, name: str) -> bool:
        """Whether the named family ``name`` (seeded via ``families=``) has ≥ 1 copy in the genome at
        ``node_id`` (across all chromosomes)."""
        if name not in self.family_names:
            raise KeyError(f"no named family {name!r}; seeded families are {sorted(self.family_names)}")
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
        extant = [n.id for n in self.complete_tree.extant()]
        return {s: tuple(g for chrom in self.genomes[s] for g in chrom.genes) for s in extant}

    @cached_property
    def profiles(self) -> Profiles:
        """The phyletic profiles — each gene family's copy count in each extant species — derived
        from the observed genomes, flattening across chromosomes (position does not enter). See
        :mod:`.profiles`."""
        return profiles_from_genomes(self._extant_genes, self._extant_genes.keys())

    @cached_property
    def gene_trees(self) -> dict[int, GeneTree]:
        """``{family id: GeneTree}`` — each family's true genealogy inside the complete tree, derived
        from the (position-blind) event log exactly as for the unordered core. See :mod:`.gene_trees`."""
        return gene_trees_from_events(self.events, self.complete_tree)

    def write(self, directory,
              outputs=("events", "profiles", "gene_order", "gene_trees", "rearrangements",
                       "chromosome_events", "event_positions")) -> None:
        """Materialise chosen ``outputs`` to ``directory`` (created if needed):

        - ``"events"`` → ``genome_events.tsv``, the gene-genealogy log (the source of truth).
        - ``"profiles"`` → ``profiles.tsv``, the family × extant-species copy-count matrix.
        - ``"gene_order"`` → ``gene_order.tsv``, every node's layout (one row per gene), ancestors
          included — so a branch's rearrangements can be replayed from its parent's genome.
        - ``"rearrangements"`` → ``rearrangements.tsv``, the inversion/transposition/translocation log.
        - ``"chromosome_events"`` → ``chromosome_events.tsv``, the chromosome genealogy edges.
        - ``"event_positions"`` → ``genome_event_positions.tsv``, where each gene-genealogy event
          happened. With ``gene_order`` and ``rearrangements`` this completes the replayable history:
          every event that changes a genome's layout now carries its coordinates.
        - ``"gene_trees"`` → ``gene_tree_fam<family>_{complete,extant}.nwk``, each family's true
          genealogy — unchanged from the unordered resolution, position being orthogonal to it.
        """
        d = pathlib.Path(directory)
        d.mkdir(parents=True, exist_ok=True)
        if "events" in outputs:
            (d / "genome_events.tsv").write_text(events_tsv(self.events))
        if "profiles" in outputs:
            (d / "profiles.tsv").write_text(self.profiles.to_tsv())
        if "gene_order" in outputs:
            (d / "gene_order.tsv").write_text(self._gene_order_tsv())
        if "rearrangements" in outputs:
            (d / "rearrangements.tsv").write_text(_rearrangements_tsv(self.rearrangements))
        if "chromosome_events" in outputs:
            (d / "chromosome_events.tsv").write_text(chromosome_events_tsv(self.chromosome_events))
        if "event_positions" in outputs:
            (d / "genome_event_positions.tsv").write_text(_event_positions_tsv(self.event_positions))
        if "gene_trees" in outputs:
            write_gene_trees(self.gene_trees, d)

    def _gene_order_tsv(self) -> str:
        cols = ("lineage", "chromosome", "position", "strand", "family", "copy")
        rows = [f"{node_label(s)}\t{ch}\t{p}\t{st}\t{fam}\t{gid}"
                for s in sorted(self.genomes)
                for (ch, p, st, fam, gid) in self.gene_order(s)]
        return "\n".join(["\t".join(cols), *rows]) + "\n"


def _rearrangements_tsv(rearrangements) -> str:
    # one table for all identity-preserving rearrangements; empty cells for fields a kind does not use
    cols = ("time", "kind", "lineage", "chromosome", "start", "length",
            "dest_chromosome", "dest_position", "flipped")

    def row(r):
        ln = node_label(r.lineage)
        if isinstance(r, Inversion):
            return (r.time, "inversion", ln, r.chromosome, r.start, r.length, "", "", "")
        if isinstance(r, Transposition):
            return (r.time, "transposition", ln, r.chromosome, r.start, r.length,
                    "", r.dest, int(r.flipped))
        return (r.time, "translocation", ln, r.source, r.start, r.length,
                r.dest, r.dest_position, int(r.flipped))
    rows = ["\t".join(str(v) for v in row(r)) for r in rearrangements]
    return "\n".join(["\t".join(cols), *rows]) + "\n"


_POSITION_COLS = ("time", "kind", "lineage", "chromosome", "start", "length", "family",
                  "donor", "recipient", "dest_position")


def _event_positions_tsv(event_positions: list[EventPosition]) -> str:
    # the positional companion to genome_events.tsv; empty cells for fields a kind does not use
    def cell(p, c):
        v = getattr(p, c)
        if v is None:
            return ""
        return node_label(v) if c in ("lineage", "donor", "recipient") else str(v)

    rows = ["\t".join(cell(p, c) for c in _POSITION_COLS) for p in event_positions]
    return "\n".join(["\t".join(_POSITION_COLS), *rows]) + "\n"


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


# --- extension: every gene-level event acts on a run of consecutive genes (the ZOMBI1 model) -------

def _extent(rng, dist, chrom, start) -> int:
    """A segment length in genes: sample the event's extension distribution, then clamp it to what
    the chromosome can carry from ``start``.

    A **linear** chromosome has ends, so the run stops at the last gene: ``1 <= m <= n - start``. A
    **circular** one has none, so the run wraps past position 0 and only the whole chromosome bounds
    it: ``1 <= m <= n``. That difference is the point of ``topology``. Clamping a circular run at the
    end of the gene array — as if the array boundary were a real end — would truncate every run that
    started near it, pull the realised mean extension below the nominal one, and leave the genes
    around position 0 covered less often than the rest."""
    m = max(1, int(dist.sample(rng)))
    n = len(chrom.genes)
    return min(m, n) if chrom.topology == "circular" else min(m, n - start)


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


# --- the mutators (position-, chromosome-, and extension-aware; each records to its log) -----------

def _originate(genome, node, t, events, positions, new_gene, new_family, rng) -> None:
    """A new gene family arises: mint a single founding gene (a family is born once — no extension) on
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
    0 on a circular chromosome. A run covering the whole chromosome empties it — the chromosome
    itself survives as an empty replicon, exactly as a de-novo one starts out; only
    :func:`_chromosome_lose` removes a chromosome from the karyotype. Returns the ``m`` removed."""
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
                 transfer_to, replacement, self_transfer, depth) -> int:
    """The segment ``[jd, jd+m)`` on the donor's chromosome ``cdi`` transfers to a contemporaneous
    recipient: each gene ends → a continuation on the donor branch and a transferred copy on the
    recipient (a horizontal gene-tree edge). The run may wrap position 0 on a circular donor
    chromosome. The transferred copies arrive as a block at a random position on a uniformly-chosen
    recipient chromosome (strands travel with them). Returns the change in total gene count: ``+m``
    additive, minus one per homologous copy displaced under ``replacement``."""
    donor = alive[kd]
    jd = _anchor(gen[kd][cdi], jd, m)
    segment = gen[kd][cdi].genes[jd:jd + m]
    cand = [k for k in range(len(alive)) if self_transfer or k != kd]
    kr = recipient_index(rng, tree, alive, cand, donor, t, transfer_to, depth)
    recipient = alive[kr]
    rgenome = gen[kr]
    conts = [new_gene(g.family, g.strand) for g in segment]
    xfers = [new_gene(g.family, g.strand) for g in segment]
    gen[kd][cdi].genes[jd:jd + m] = conts               # continuations replace the segment on the donor
    # the donor's row first, then any displacements it causes, then the arrival: within one timestamp
    # the rows are written in the order a replayer must apply them
    positions.append(EventPosition(t, "transfer_donor", donor, gen[kd][cdi].id, jd, m,
                                   donor=donor, recipient=recipient))
    delta = m
    if replacement:
        cont_ids = {c.id for c in conts}                # self-transfer: never overwrite our own conts
        for x in xfers:                                 # each arriving copy may displace a homolog
            residents = [(ci, p) for ci, ch in enumerate(rgenome) for p, c in enumerate(ch.genes)
                         if c.family == x.family and c.id not in cont_ids]
            if residents:
                ci, p = residents[int(rng.integers(len(residents)))]
                victim = rgenome[ci].genes[p]
                del rgenome[ci].genes[p]
                events.append(Event(t, "loss", recipient, victim.family, victim.id))
                positions.append(EventPosition(t, "loss", recipient, rgenome[ci].id, p, 1))
                delta -= 1
    rchrom = rgenome[int(rng.integers(len(rgenome)))]   # arrive as a block on a random recipient chromosome
    pos = int(rng.integers(len(rchrom.genes) + 1))
    rchrom.genes[pos:pos] = xfers
    positions.append(EventPosition(t, "transfer_recipient", recipient, rchrom.id, pos, m,
                                   donor=donor, recipient=recipient))
    for old, cont, xf in zip(segment, conts, xfers):
        events.append(Event(t, "transfer", donor, old.family, cont.id, parent=old.id))
        events.append(Event(t, "transfer", recipient, old.family, xf.id, parent=old.id, recipient=recipient))
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


# --- seeding + validation -------------------------------------------------------------------------

def _topologies(chromosomes, topology) -> list[str]:
    """Resolve the ``topology`` argument to one label per seeded chromosome."""
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

def simulate_genomes_ordered(tree, *, duplication=0.0, transfer=0.0, loss=0.0, origination=0.0,
                             inversion=0.0, transposition=0.0, translocation=0.0,
                             chromosomes=1, topology="circular",
                             fission=0.0, fusion=0.0, chromosome_origination=0.0, chromosome_loss=0.0,
                             duplication_extension=None, loss_extension=None, transfer_extension=None,
                             inversion_extension=None, transposition_extension=None,
                             translocation_extension=None, inversion_probability=0.0,
                             transfer_to="uniform", replacement=False, self_transfer=False,
                             initial_families=0, families=None, seed=None,
                             progress=False) -> OrderedGenomesResult:
    """Evolve ordered genomes — genes with a position and an orientation, on chromosomes — along a
    species tree, by the D/T/L/O core plus segmental rearrangements and the chromosome tier.

    **Every gene-level event acts on an *extension*** — a run of consecutive genes (the ZOMBI1 model):
    ``duplication`` copies the run in tandem, ``loss`` removes it, ``transfer`` sends it to a
    contemporaneous recipient as a block, ``inversion`` reverses it (flipping strands), ``transposition``
    relocates it elsewhere on the same chromosome, and ``translocation`` moves it to a different
    chromosome. The run's length is drawn per event from ``<event>_extension`` (a distribution,
    default ``Geometric(mean=1)`` — usually a single gene; dial the mean up for larger blocks).
    ``origination`` is the exception: a family is born once, a single gene, no extension.
    ``transposition`` and ``translocation`` land the moved block inverted with probability
    ``inversion_probability``.

    **Where a run stops is set by the chromosome's ``topology``.** A run goes rightwards from the gene
    it starts at. On a ``"circular"`` chromosome there are no ends, so a run that reaches the last
    gene continues from the first, and only the whole chromosome bounds it; on a ``"linear"`` one the
    run stops at the last gene. So on a circular chromosome every gene is covered by segmental events
    at the same rate, and the nominal mean extension is the realised one.

    Scopes follow the cross-level grammar, which counts an event per the thing it acts on: the
    gene-level events — ``duplication``/``transfer``/``loss`` and the rearrangements
    ``inversion``/``transposition``/``translocation`` — are **per copy**, since each acts on a run of
    genes that starts at one of them; the chromosome tier ``fission``/``fusion``/``chromosome_loss``
    is **per chromosome**; and the two events that make something from nothing,
    ``origination``/``chromosome_origination``, are **per lineage**. The
    root is seeded with ``chromosomes`` chromosomes of the given ``topology``, across which the
    ``initial_families`` founding genes are dealt **round-robin**; ``families=["toxin", …]`` additionally
    seeds **named** families (remembered in ``result.family_names`` for ``result.has_family(node,
    "toxin")``), as in the unordered core; ``transfer_to`` / ``replacement`` / ``self_transfer`` behave
    as in the unordered core.

    The **chromosome tier** changes chromosome *number*: ``fission`` (split), ``fusion`` (merge — the
    reticulation), ``chromosome_origination`` (a de-novo replicon), ``chromosome_loss`` (a whole
    chromosome and its genes die; never the genome's last). Chromosomes carry identity — re-minted at
    every event that reshapes them — so ``chromosome_events`` is the true reticulating chromosome
    genealogy, rooted at the seed and de-novo originations. Deterministic given ``seed``.
    """
    tree = tree.complete_tree if isinstance(tree, SpeciesResult) else tree
    labels = _topologies(chromosomes, topology)
    n_chrom_seed = chromosomes
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
    # like the unordered core, this slice wires only the default scope of each event and OnTime
    # (skyline) modifiers; a scope override or per-family/clade modifier is a later slice, so reject
    # them rather than silently mis-scale (see the unordered engine for the reasoning).
    for label, rate, want in (("duplication", dup, PerCopy), ("transfer", tra, PerCopy),
                              ("loss", los, PerCopy), ("origination", org, PerLineage),
                              ("inversion", inv, PerCopy), ("transposition", trp, PerCopy),
                              ("translocation", trl, PerCopy), ("fission", fis, PerChromosome),
                              ("fusion", fus, PerChromosome), ("chromosome_loss", clo, PerChromosome),
                              ("chromosome_origination", cor, PerLineage)):
        if not isinstance(rate.scope, want):
            raise ValueError(
                f"{label} has a {type(rate.scope).__name__} scope, but the ordered genome engine "
                f"wires only {want.__name__} for {label} this slice — scope overrides are a later slice."
            )
        for m in rate.modifiers:
            if not isinstance(m, OnTime):
                raise ValueError(
                    f"{label} carries {type(m).__name__}, which the ordered genome engine does not "
                    f"support yet — only OnTime (skyline) is wired. Per-family heterogeneity and clade "
                    f"drift are later slices."
                )
    # per-event extension distributions (segment length in genes); default to a single gene
    def _ext(spec):
        return as_distribution(spec) if spec is not None else Geometric(mean=1)
    dup_ext, los_ext, tra_ext = _ext(duplication_extension), _ext(loss_extension), _ext(transfer_extension)
    inv_ext, trp_ext, trl_ext = (_ext(inversion_extension), _ext(transposition_extension),
                                 _ext(translocation_extension))
    if not 0.0 <= inversion_probability <= 1.0:
        raise ValueError(f"inversion_probability must be in [0, 1], got {inversion_probability!r}")
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
    chrom_counter = 0

    def new_gene(family: int, strand: int) -> Gene:
        nonlocal copy_counter
        g = Gene(copy_counter, family, strand)
        copy_counter += 1
        return g

    def new_family() -> int:
        nonlocal family_counter
        f = family_counter
        family_counter += 1
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

    root_chroms = []
    for label in labels:  # seed the root karyotype; each seeded chromosome is a network root
        cid = new_chromosome()
        root_chroms.append(Chromosome(cid, label, []))
        chromosome_events.append(ChromosomeEvent(t, "origination", root.id, (), (cid,)))
    # the crown seeding is logged like any other origination — each founding gene appended in turn —
    # so the position table is total over gene-content events and a replay of the root branch can
    # start from an empty karyotype (every other branch starts from its parent's gene_order rows)
    for i in range(initial_families):  # deal the founding genes round-robin across the chromosomes
        fam = new_family()
        chrom = root_chroms[i % n_chrom_seed]
        chrom.genes.append(new_gene(fam, +1))
        events.append(Event(t, "origination", root.id, fam, chrom.genes[-1].id))
        event_positions.append(EventPosition(t, "origination", root.id, chrom.id,
                                             len(chrom.genes) - 1, 1, family=fam))
    family_names: dict[str, int] = {}  # named crown families, dealt round-robin after the anonymous ones
    for j, name in enumerate(families):
        fam = new_family()
        family_names[name] = fam
        chrom = root_chroms[(initial_families + j) % n_chrom_seed]
        chrom.genes.append(new_gene(fam, +1))
        events.append(Event(t, "origination", root.id, fam, chrom.genes[-1].id))
        event_positions.append(EventPosition(t, "origination", root.id, chrom.id,
                                             len(chrom.genes) - 1, 1, family=fam))
    enter(alive, gen, pos, root.id, root_chroms)
    total_copies = initial_families + len(families)
    total_chromosomes = n_chrom_seed

    bar = progress_bar(len(schedule), "genomes", unit="branch", enabled=progress)
    si = 0
    while si < len(schedule):
        bar.to(si)
        n = total_copies
        k_alive = len(alive)
        ctx = {"copies": n, "lineages": k_alive, "chromosomes": total_chromosomes, "time": t}
        c = total_chromosomes
        r_dup = dup.effective(**ctx) if n else 0.0
        r_los = los.effective(**ctx) if n else 0.0
        r_org = org.effective(**ctx)                                    # per lineage
        can_xfer = n > 0 and (k_alive >= 2 or self_transfer)
        r_tra = tra.effective(**ctx) if can_xfer else 0.0
        r_inv = inv.effective(**ctx) if n else 0.0                      # per copy (the run's start)
        r_trp = trp.effective(**ctx) if n else 0.0                      # per copy (the run's start)
        r_trl = trl.effective(**ctx) if n else 0.0                      # per copy; needs >=2 chromosomes
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
                if r < r_dup:                            # every gene-level event acts on an extension
                    k, ci, j = _pick_gene(rng, gen, n)
                    chrom = gen[k][ci]
                    m = _extent(rng, dup_ext, chrom, j)
                    total_copies += _duplicate(chrom, j, m, tree.nodes[alive[k]], t, events,
                                               event_positions, new_gene)
                elif r < b_los:
                    k, ci, j = _pick_gene(rng, gen, n)
                    chrom = gen[k][ci]
                    m = _extent(rng, los_ext, chrom, j)
                    total_copies -= _lose_at(chrom, j, m, tree.nodes[alive[k]], t, events,
                                             event_positions)
                elif r < b_org:
                    k = int(rng.integers(k_alive))  # origination is per lineage: a uniform lineage
                    _originate(gen[k], tree.nodes[alive[k]], t, events, event_positions, new_gene,
                               new_family, rng)
                    total_copies += 1
                elif r < b_tra:
                    kd, cdi, jd = _pick_gene(rng, gen, n)
                    m = _extent(rng, tra_ext, gen[kd][cdi], jd)
                    total_copies += _do_transfer(rng, tree, alive, gen, kd, cdi, jd, m, t, events,
                                                 event_positions, new_gene, transfer_to, replacement,
                                                 self_transfer, depth)
                elif r < b_inv:
                    k, ci, i0 = _pick_gene(rng, gen, n)   # the run starts at a gene, so: per copy
                    chrom = gen[k][ci]
                    _invert(chrom, i0, _extent(rng, inv_ext, chrom, i0),
                            tree.nodes[alive[k]], t, rearrangements)
                elif r < b_trp:
                    k, ci, i0 = _pick_gene(rng, gen, n)
                    chrom = gen[k][ci]
                    _transpose(chrom, i0, _extent(rng, trp_ext, chrom, i0),
                               tree.nodes[alive[k]], t, rearrangements, rng, inversion_probability)
                elif r < b_trl:
                    k, ci, j = _pick_gene(rng, gen, n)
                    m = _extent(rng, trl_ext, gen[k][ci], j)
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
                        for c in node.children:
                            dcid = new_chromosome()
                            dcids.append(dcid)
                            dgenes = []
                            for old in pchrom.genes:  # ZOMBI1: the gene ends and continues, fresh id
                                ng = new_gene(old.family, old.strand)
                                dgenes.append(ng)
                                events.append(Event(t, "speciation", c, old.family, ng.id, parent=old.id))
                            child_genomes[c].append(Chromosome(dcid, pchrom.topology, dgenes))
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
                                family_names, event_positions)


__all__ = ["simulate_genomes_ordered", "OrderedGenomesResult", "Gene", "Chromosome",
           "ChromosomeEvent", "Inversion", "Transposition", "Translocation"]
