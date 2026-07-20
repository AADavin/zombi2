"""Genomes II — ordered: genes carry a position and an orientation, on chromosomes.

The ordered resolution layers **position** over the unordered D/T/L/O core (Chapter 6). A genome is
no longer a multiset of gene copies but a list of **chromosomes**, each an ordered run of oriented
:class:`Gene`\\ s. This slice adds one rearrangement — the **inversion** (reverse a span of a
chromosome, flipping each gene's strand: the classic signed-permutation move) — and gives
chromosomes a genuine **identity**: a chromosome id is re-minted at every speciation and the
speciation edge is recorded, so the chromosome genealogy exists from the start (a tree for now; it
reticulates once fission/fusion arrive — the next slice).

It is the genome twin of the unordered core and shares its spine: one forward Gillespie over the
**complete** species tree, the same ``scope(base) × modifiers`` rate grammar, the same gene-genealogy
:class:`~zombi2.genomes.events.Event` log (position-blind, so ``gene_trees`` and ``profiles`` are
derived from it unchanged), and the same live-lineage bookkeeping. What differs is the state (a list
of chromosomes, not a flat multiset) and the mutators (position-aware: duplication is tandem, loss is
in-place, transfers arrive at a position) plus two extra logs — ``rearrangements`` (inversions) and
``chromosome_events`` (the chromosome genealogy).

Slice 2 adds the **chromosome tier** — the events that change chromosome *number*: fission (a
bifurcation), fusion (the reticulation), chromosome origination (a de-novo replicon) and chromosome
loss. The chromosome genealogy is now a genuine reticulating network, recorded as the
``chromosome_events`` edge list (its ground truth). Still to come: the **eNewick** serialisation of
that network and the gene→chromosome path (``chromosome-network.md``); the gene-movement
rearrangements — transposition and translocation — that shuffle genes between chromosome lineages
without ending them; then the nucleotide resolution.
"""

from __future__ import annotations

import collections
import pathlib
from dataclasses import dataclass
from functools import cached_property

import numpy as np

from ..rates.modifiers import Time
from ..rates.rate import as_rate
from ..rates.scope import PerChromosome, PerCopy, PerLineage
from ..species import SpeciesResult, Tree
from ._live import enter, retire
from ._transfer import Distance, mean_root_to_tip, recipient_index
from .events import Event, events_tsv
from .gene_trees import GeneTree, gene_trees_from_events
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
    ``"linear"``. In this slice topology is just a label (it will gate which fissions/fusions are
    legal later); it does not affect inversions."""

    id: int
    topology: str
    genes: list[Gene]


@dataclass(frozen=True)
class Inversion:
    """A recorded inversion: on species branch ``lineage`` at ``time``, the genes in positions
    ``start``..``end`` (inclusive) of chromosome ``chromosome`` were reversed and their strands
    flipped. Gene ids are untouched — an inversion reshapes order, it does not end lineages — so it is
    logged here, separate from the gene-genealogy :class:`~zombi2.genomes.events.Event` stream."""

    time: float
    lineage: int
    chromosome: int
    start: int
    end: int

    @property
    def length(self) -> int:
        """Number of genes inverted (``end - start + 1``)."""
        return self.end - self.start + 1


@dataclass(frozen=True)
class ChromosomeEvent:
    """One edge of the **chromosome genealogy** — a chromosome lineage's birth, split, merge, or
    death, fired on species branch ``lineage`` at ``time``. ``parents`` → ``children`` are chromosome
    ids and the arity names the event: ``"origination"`` (``()`` → one child: a seed or de-novo
    replicon, a **root**), ``"speciation"`` and ``"fission"`` (one parent → two children, a
    **bifurcation**), ``"fusion"`` (two parents → one child, the **reticulation** — in-degree 2, what
    makes this a network and not a tree), ``"loss"`` (one parent → ``()``, a **leaf**). The edge list
    is the network's ground truth; the eNewick serialisation is derived from it."""

    time: float
    kind: str  # "origination" | "speciation" | "fission" | "fusion" | "loss"
    lineage: int
    parents: tuple[int, ...]
    children: tuple[int, ...]


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
    rearrangements: list[Inversion]
    chromosome_events: list[ChromosomeEvent]
    seed: int | None

    def family_counts(self, node_id: int) -> collections.Counter:
        """A multiset view of one node's genome: ``family id → copy count`` (across all chromosomes)."""
        return collections.Counter(g.family for chrom in self.genomes[node_id] for g in chrom.genes)

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

    def write(self, directory, outputs=("events", "profiles", "gene_order")) -> None:
        """Materialise chosen ``outputs`` to ``directory`` (created if needed):

        - ``"events"`` → ``genome_events.tsv``, the gene-genealogy log (the source of truth).
        - ``"profiles"`` → ``profiles.tsv``, the family × extant-species copy-count matrix.
        - ``"gene_order"`` → ``gene_order.tsv``, the observed genomes' layout (one row per gene).
        - ``"rearrangements"`` → ``rearrangements.tsv``, the inversion log.
        - ``"chromosome_events"`` → ``chromosome_events.tsv``, the chromosome genealogy edges.
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
            (d / "chromosome_events.tsv").write_text(_chromosome_events_tsv(self.chromosome_events))

    def _gene_order_tsv(self) -> str:
        cols = ("species", "chromosome", "position", "strand", "family", "gene")
        rows = [f"{s}\t{ch}\t{p}\t{st}\t{fam}\t{gid}"
                for s in sorted(n.id for n in self.complete_tree.extant())
                for (ch, p, st, fam, gid) in self.gene_order(s)]
        return "\n".join(["\t".join(cols), *rows]) + "\n"


def _rearrangements_tsv(rearrangements: list[Inversion]) -> str:
    cols = ("time", "kind", "lineage", "chromosome", "start", "end", "length")
    rows = [f"{r.time}\tinversion\t{r.lineage}\t{r.chromosome}\t{r.start}\t{r.end}\t{r.length}"
            for r in rearrangements]
    return "\n".join(["\t".join(cols), *rows]) + "\n"


def _chromosome_events_tsv(chromosome_events: list[ChromosomeEvent]) -> str:
    cols = ("time", "kind", "lineage", "parents", "children")
    rows = [f"{e.time}\t{e.kind}\t{e.lineage}\t{';'.join(map(str, e.parents))}\t"
            f"{';'.join(map(str, e.children))}" for e in chromosome_events]
    return "\n".join(["\t".join(cols), *rows]) + "\n"


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


# --- the mutators (position- and chromosome-aware; each records to its log) ------------------------

def _originate(genome, node, t, events, new_gene, new_family, rng) -> None:
    """A new gene family arises in this lineage: mint a founding gene on a uniformly-chosen
    chromosome at a uniformly-chosen position (strand ``+1``), and record it."""
    chrom = genome[int(rng.integers(len(genome)))]
    fam = new_family()
    g = new_gene(fam, +1)
    chrom.genes.insert(int(rng.integers(len(chrom.genes) + 1)), g)
    events.append(Event(t, "origination", node.id, fam, g.id))


def _duplicate(chrom, j, node, t, events, new_gene) -> None:
    """The gene at position ``j`` duplicates **in tandem**: the gene ends and two fresh copies (same
    strand) descend from it — the continuation in place and the new copy immediately after."""
    old = chrom.genes[j]
    cont, dup = new_gene(old.family, old.strand), new_gene(old.family, old.strand)
    chrom.genes[j] = cont
    chrom.genes.insert(j + 1, dup)                      # tandem: adjacent, same orientation
    events.append(Event(t, "duplication", node.id, old.family, cont.id, parent=old.id))
    events.append(Event(t, "duplication", node.id, old.family, dup.id, parent=old.id))


def _lose_at(chrom, j, node, t, events) -> None:
    """The gene at position ``j`` is lost — removed **in place** (order is preserved, unlike the
    unordered core's swap-remove)."""
    lost = chrom.genes[j]
    del chrom.genes[j]
    events.append(Event(t, "loss", node.id, lost.family, lost.id))


def _invert(chrom, i, j, node, t, rearrangements) -> None:
    """Invert positions ``i``..``j`` (inclusive): reverse the run and flip each gene's strand. Gene
    ids are untouched — identity persists through an inversion — so nothing is written to the gene
    event log, only the rearrangement log."""
    chrom.genes[i:j + 1] = [Gene(g.id, g.family, -g.strand) for g in reversed(chrom.genes[i:j + 1])]
    rearrangements.append(Inversion(t, node.id, chrom.id, i, j))


def _do_transfer(rng, tree, alive, gen, total_copies, t, events, new_gene,
                 transfer_to, replacement, self_transfer, depth) -> int:
    """A gene transfers from a donor gene to a contemporaneous recipient lineage, arriving on a
    uniformly-chosen chromosome at a uniformly-chosen position (its strand travels with it). Returns
    the change in total gene count: +1 additive, 0 replacement (the arriving copy displaces a
    homologous resident)."""
    kd, cdi, jd = _pick_gene(rng, gen, total_copies)
    donor = alive[kd]
    src = gen[kd][cdi].genes[jd]
    fam = src.family
    cand = [k for k in range(len(alive)) if self_transfer or k != kd]
    kr = recipient_index(rng, tree, alive, cand, donor, t, transfer_to, depth)
    recipient = alive[kr]
    rgenome = gen[kr]
    # the donor gene ends; two fresh copies descend (ZOMBI1 re-id), same strand: the continuation on
    # the donor branch and the transferred copy on the recipient branch — a horizontal gene-tree edge.
    cont, xfer = new_gene(fam, src.strand), new_gene(fam, src.strand)
    gen[kd][cdi].genes[jd] = cont
    delta = 1
    if replacement:
        residents = [(ci, p) for ci, chrom in enumerate(rgenome)
                     for p, c in enumerate(chrom.genes) if c.family == fam and c.id != cont.id]
        if residents:  # homologous overwrite; empty ⇒ additive fallback (the gene still arrives)
            ci, p = residents[int(rng.integers(len(residents)))]
            victim = rgenome[ci].genes[p]
            del rgenome[ci].genes[p]
            events.append(Event(t, "loss", recipient, fam, victim.id))
            delta = 0
    rchrom = rgenome[int(rng.integers(len(rgenome)))]   # arrive on a uniformly-chosen chromosome
    rchrom.genes.insert(int(rng.integers(len(rchrom.genes) + 1)), xfer)
    events.append(Event(t, "transfer", donor, fam, cont.id, parent=src.id))
    events.append(Event(t, "transfer", recipient, fam, xfer.id, parent=src.id, recipient=recipient))
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


def _chromosome_lose(genome, ci, node, t, events, chromosome_events) -> tuple[int, int]:
    """A whole chromosome and its genes die — a **leaf** of the chromosome network (no child); each
    gene on it ends as a gene ``loss``. No-op if it is the genome's last chromosome (a lineage never
    loses its entire genome this way)."""
    if len(genome) < 2:
        return (0, 0)
    lost = genome[ci]
    for g in lost.genes:
        events.append(Event(t, "loss", node.id, g.family, g.id))
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
                             inversion=0.0, chromosomes=1, topology="circular",
                             fission=0.0, fusion=0.0, chromosome_origination=0.0, chromosome_loss=0.0,
                             transfer_to="uniform", replacement=False, self_transfer=False,
                             initial_families=0, seed=None) -> OrderedGenomesResult:
    """Evolve ordered genomes — genes with a position and an orientation, on chromosomes — along a
    species tree by duplication, transfer, loss, origination, and **inversion**.

    Everything the unordered core does, it does here (same ``tree`` input — a
    :class:`~zombi2.species.Tree` or :class:`~zombi2.species.SpeciesResult`; evolution on **every**
    lineage; the ``duplication``/``transfer``/``loss`` per-copy and ``origination`` per-lineage
    defaults; the ``transfer_to`` / ``replacement`` / ``self_transfer`` mechanics), with position
    added: duplication is **tandem**, loss removes in place, transfers and originations arrive at a
    random position. The root is seeded with ``chromosomes`` chromosomes of the given ``topology``,
    across which the ``initial_families`` founding genes are dealt **round-robin**.

    ``inversion`` (default **per chromosome**) reverses a random span of a chromosome and flips each
    gene's strand. The **chromosome tier** changes chromosome *number*: ``fission`` (a chromosome
    splits in two — per chromosome), ``fusion`` (two chromosomes in a genome merge — per chromosome —
    the reticulation), ``chromosome_origination`` (a de-novo replicon — per lineage), and
    ``chromosome_loss`` (a whole chromosome and its genes die — per chromosome; never the genome's
    last). Chromosomes carry identity — re-minted at every event that reshapes them — so
    ``chromosome_events`` is the true chromosome genealogy: a tree of speciations/fissions with fusion
    reticulations, rooted at the seed and de-novo originations. Deterministic given ``seed``.
    """
    tree = tree.complete_tree if isinstance(tree, SpeciesResult) else tree
    labels = _topologies(chromosomes, topology)
    n_chrom_seed = chromosomes
    dup = as_rate(duplication, default_scope=PerCopy)
    tra = as_rate(transfer, default_scope=PerCopy)
    los = as_rate(loss, default_scope=PerCopy)
    org = as_rate(origination, default_scope=PerLineage)
    inv = as_rate(inversion, default_scope=PerChromosome)
    fis = as_rate(fission, default_scope=PerChromosome)
    fus = as_rate(fusion, default_scope=PerChromosome)
    cor = as_rate(chromosome_origination, default_scope=PerLineage)
    clo = as_rate(chromosome_loss, default_scope=PerChromosome)
    # like the unordered core, this slice wires only the default scope of each event and Time
    # (skyline) modifiers; a scope override or per-family/clade modifier is a later slice, so reject
    # them rather than silently mis-scale (see the unordered engine for the reasoning).
    for label, rate, want in (("duplication", dup, PerCopy), ("transfer", tra, PerCopy),
                              ("loss", los, PerCopy), ("origination", org, PerLineage),
                              ("inversion", inv, PerChromosome), ("fission", fis, PerChromosome),
                              ("fusion", fus, PerChromosome), ("chromosome_loss", clo, PerChromosome),
                              ("chromosome_origination", cor, PerLineage)):
        if not isinstance(rate.scope, want):
            raise ValueError(
                f"{label} has a {type(rate.scope).__name__} scope, but the ordered genome engine "
                f"wires only {want.__name__} for {label} this slice — scope overrides are a later slice."
            )
        for m in rate.modifiers:
            if not isinstance(m, Time):
                raise ValueError(
                    f"{label} carries {type(m).__name__}, which the ordered genome engine does not "
                    f"support yet — only Time (skyline) is wired. Per-family heterogeneity and clade "
                    f"drift are later slices."
                )
    if transfer_to == "distance":
        transfer_to = Distance()
    if transfer_to != "uniform" and not isinstance(transfer_to, Distance):
        raise ValueError(f"transfer_to must be 'uniform', 'distance', or Distance(decay=), got {transfer_to!r}")
    if isinstance(initial_families, bool) or not isinstance(initial_families, int) or initial_families < 0:
        raise ValueError(f"initial_families must be a non-negative integer, got {initial_families!r}")

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
    rearrangements: list[Inversion] = []
    chromosome_events: list[ChromosomeEvent] = []

    root_chroms = []
    for label in labels:  # seed the root karyotype; each seeded chromosome is a network root
        cid = new_chromosome()
        root_chroms.append(Chromosome(cid, label, []))
        chromosome_events.append(ChromosomeEvent(t, "origination", root.id, (), (cid,)))
    for i in range(initial_families):  # deal the founding genes round-robin across the chromosomes
        fam = new_family()
        root_chroms[i % n_chrom_seed].genes.append(new_gene(fam, +1))
        events.append(Event(t, "origination", root.id, fam, root_chroms[i % n_chrom_seed].genes[-1].id))
    enter(alive, gen, pos, root.id, root_chroms)
    total_copies = initial_families
    total_chromosomes = n_chrom_seed

    si = 0
    while si < len(schedule):
        n = total_copies
        k_alive = len(alive)
        ctx = {"copies": n, "lineages": k_alive, "chromosomes": total_chromosomes, "time": t}
        c = total_chromosomes
        r_dup = dup.effective(**ctx) if n else 0.0
        r_los = los.effective(**ctx) if n else 0.0
        r_org = org.effective(**ctx)                                    # per lineage
        can_xfer = n > 0 and (k_alive >= 2 or self_transfer)
        r_tra = tra.effective(**ctx) if can_xfer else 0.0
        r_inv = inv.effective(**ctx) if n else 0.0                      # per chromosome; needs a gene
        r_fis = fis.effective(**ctx) if c else 0.0                      # per chromosome (the tier)
        r_fus = fus.effective(**ctx) if c else 0.0
        r_cor = cor.effective(**ctx)                                    # per lineage (de-novo replicon)
        r_clo = clo.effective(**ctx) if c else 0.0
        total = r_dup + r_los + r_org + r_tra + r_inv + r_fis + r_fus + r_cor + r_clo

        next_species = schedule[si][0]
        horizon = min(next_species, dup.next_change(t), los.next_change(t), org.next_change(t),
                      tra.next_change(t), inv.next_change(t), fis.next_change(t), fus.next_change(t),
                      cor.next_change(t), clo.next_change(t))

        if total > 0.0:
            t_ev = t + float(rng.exponential(1.0 / total))
            if t_ev < horizon:  # a genome event fires before the alive set or a rate changes
                t = t_ev
                r = float(rng.random()) * total
                b_los = r_dup + r_los                    # cumulative bounds, in the firing order below
                b_org = b_los + r_org
                b_tra = b_org + r_tra
                b_inv = b_tra + r_inv
                b_fis = b_inv + r_fis
                b_fus = b_fis + r_fus
                b_cor = b_fus + r_cor                    # ... and the remainder (to total) is clo
                if r < r_dup:
                    k, ci, j = _pick_gene(rng, gen, n)
                    _duplicate(gen[k][ci], j, tree.nodes[alive[k]], t, events, new_gene)
                    total_copies += 1
                elif r < b_los:
                    k, ci, j = _pick_gene(rng, gen, n)
                    _lose_at(gen[k][ci], j, tree.nodes[alive[k]], t, events)
                    total_copies -= 1
                elif r < b_org:
                    k = int(rng.integers(k_alive))  # origination is per lineage: a uniform lineage
                    _originate(gen[k], tree.nodes[alive[k]], t, events, new_gene, new_family, rng)
                    total_copies += 1
                elif r < b_tra:
                    total_copies += _do_transfer(rng, tree, alive, gen, n, t, events, new_gene,
                                                 transfer_to, replacement, self_transfer, depth)
                elif r < b_inv:
                    k, ci = _pick_chromosome(rng, gen, c)
                    chrom = gen[k][ci]
                    if chrom.genes:  # an empty chromosome has nothing to invert (a no-op; rare)
                        i0 = int(rng.integers(len(chrom.genes)))       # two uniform cut points
                        i1 = int(rng.integers(len(chrom.genes)))
                        _invert(chrom, min(i0, i1), max(i0, i1), tree.nodes[alive[k]], t, rearrangements)
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
                                              chromosome_events)
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

    return OrderedGenomesResult(tree, genomes, events, rearrangements, chromosome_events, seed)


__all__ = ["simulate_genomes_ordered", "OrderedGenomesResult", "Gene", "Chromosome",
           "ChromosomeEvent", "Inversion"]
