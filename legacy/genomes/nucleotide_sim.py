"""Forward simulation + trace-back post-processing for the nucleotide genome (M1).

Runs the *unchanged* :class:`~zombi2.genome_sim.GenomeSimulator` over
:class:`~zombi2.nucleotide_genome.NucleotideGenome` with a genome-level inversion rate,
then decomposes the result:

* **blocks** — the finest intervals of the ancestral genome that are never cut by a
  breakpoint in *any* extant leaf. Each block is a "segment with one shared history".
* **trace-back** — for each extant leaf, the ancestral origin of every nucleotide (this
  is just the leaf genome's :meth:`~zombi2.nucleotide_genome.NucleotideGenome.to_cells`).
* **mosaic** — each extant genome as an ordered, signed sequence of blocks.
* **block histories** — the inversions (branch, time) that touched each block.

For inversion-only M1 every block survives in every leaf exactly once (inversion neither
creates nor destroys sequence), so a block's genealogy is simply the species tree;
genuine per-block gene trees arrive with duplication/loss/transfer in later milestones.
"""

from __future__ import annotations

from bisect import bisect_left
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from zombi2._sampling import EventSampler
from zombi2.genomes.events import EventType, EventRecord, GeneOp
from zombi2.genomes.genome_sim import GenomeSimulator
from zombi2.genomes.nucleotide_genome import NucleotideGenome, SegmentRegistry
from zombi2.genomes.rates import Rates
from zombi2.genomes.reconciliation import build_gene_trees, reconcile
from zombi2.tree import Tree, TreeNode

# Runaway-growth ceiling for the length-scaled nucleotide engine. Structural-event rates here
# are proportional to genome length, so an additive gain rate (transfer/duplication) that outruns
# loss compounds without bound and the walk never terminates. A single genome that crosses this
# many segments is unambiguously in that regime: the largest genome any balanced (bounded) run
# reaches in the test suite is a few hundred segments, so this leaves a >~30x margin while turning
# a would-be infinite hang into a prompt, actionable error. Override via ``max_segments_per_genome``.
DEFAULT_MAX_SEGMENTS_PER_GENOME = 20_000


@dataclass(frozen=True)
class Block:
    """A maximal uncut interval ``[start, end)`` of one ancestral ``source``.

    ``kind`` is ``"gene"`` when the interval is (part of) a gene annotation, else
    ``"intergene"``; ``gene_id`` names the gene (``None`` for intergene). Classification is by
    *ancestral coordinate*, so it is unaffected by pseudogenization (an ancestrally-gene block
    stays a gene tree; a pseudogenization surfaces as a state-change edge inside that tree).
    """

    block_id: int
    source: str
    start: int
    end: int
    kind: str = "intergene"
    gene_id: str | None = None

    @property
    def length(self) -> int:
        return self.end - self.start


@dataclass
class NucleotideResult:
    """Output of :func:`simulate_nucleotide_genomes`."""

    species_tree: Tree
    leaf_genomes: dict  # extant leaf TreeNode -> NucleotideGenome
    event_log: object   # EventLog
    registry: SegmentRegistry  # segment provenance + split parent-links
    blocks: list         # list[Block], tiling every source's [0, len)
    root_length: int
    node_genomes: dict = field(default_factory=dict)  # every node -> genome (retain_internal)
    _by_source: dict = field(default_factory=dict, repr=False)

    def __post_init__(self):
        by_source: dict[str, list[Block]] = {}
        for a in self.blocks:
            by_source.setdefault(a.source, []).append(a)
        for blocks in by_source.values():
            blocks.sort(key=lambda a: a.start)
        self._by_source = by_source
        # parallel start arrays so any [src_start, src_end) segment finds its blocks by bisection
        self._starts = {src: [a.start for a in blocks] for src, blocks in by_source.items()}
        self._block_by_id = {a.block_id: a for a in self.blocks}

    def _covered(self, source: str, ss: int, se: int) -> list:
        """Blocks of ``source`` lying in ``[ss, se)`` (segment boundaries are block boundaries)."""
        blocks = self._by_source.get(source)
        if not blocks:
            return []
        st = self._starts[source]
        return blocks[bisect_left(st, ss):bisect_left(st, se)]

    # --- per-leaf views ----------------------------------------------------
    def trace_back(self, leaf: TreeNode) -> list[tuple[str, int, int]]:
        """Ancestral origin ``(source, src_pos, strand)`` of every nucleotide at ``leaf``."""
        return self.leaf_genomes[leaf].to_cells()

    def leaf_mosaic(self, leaf: TreeNode) -> list[tuple[int, int]]:
        """The leaf genome as an ordered, signed sequence of blocks: ``[(block_id, strand)]``."""
        out: list[tuple[int, int]] = []
        for seg in self.leaf_genomes[leaf]._iter_segments():
            covered = self._covered(seg.source, seg.src_start, seg.src_end)
            if seg.strand == -1:
                covered = reversed(covered)
            out.extend((a.block_id, seg.strand) for a in covered)
        return out

    # --- emergent phylogenetic profile (blocks as gene families) ------------
    def profile_matrix(self):
        """``(block_ids, species, matrix)`` — copy number of each block per extant leaf.

        Blocks are the emergent "gene families": a maximal block with one shared history.
        With loss only, entries are 0/1 (presence); duplication will lift them above 1.
        This is the phylogenetic-profile dataset the block decomposition produces for free.
        """
        leaves = sorted(self.leaf_genomes, key=lambda n: n.name)
        block_ids = [a.block_id for a in self.blocks]
        row = {a.block_id: i for i, a in enumerate(self.blocks)}
        matrix = np.zeros((len(self.blocks), len(leaves)), dtype=int)
        for j, leaf in enumerate(leaves):
            for seg in self.leaf_genomes[leaf]._iter_segments():
                for a in self._covered(seg.source, seg.src_start, seg.src_end):
                    matrix[row[a.block_id], j] += 1
        return block_ids, [n.name for n in leaves], matrix

    # --- per-block gene trees (steps 6-7: reconstruct the gene of each segment) ---
    def _top(self, sid: str, cache: dict) -> str:
        """The real-born ancestor of ``sid`` — climb split parent-links (degree-2 nodes)."""
        if sid in cache:
            return cache[sid]
        chain = []
        split_parent = self.registry.split_parent
        while sid in split_parent:
            chain.append(sid)
            sid = split_parent[sid]
        for c in chain:
            cache[c] = sid
        cache[sid] = sid
        return sid

    def _block_records(self):
        """Per-block ``(records, gid2species)`` for gene-tree / reconciliation reconstruction.

        Each block (a never-cut ancestral block) is an emergent gene family; its genealogy is
        the segment lineage tree restricted to segments covering the block, with breakpoint
        splits contracted (each split sends the block to exactly one child). The log is
        indexed **once**: every event's segment covers a contiguous run of blocks (source
        coordinates are ordered), found by bisection, so each record is appended straight to
        those blocks — no per-block rescan of the whole log. Records carry the species branch
        and (for transfers) the recipient, so both ``build_gene_trees`` and ``reconcile``
        can consume them.
        """
        prov = self.registry.provenance
        top_cache: dict[str, str] = {}
        records_by_block: dict[int, list] = {a.block_id: [] for a in self.blocks}
        species_by_block: dict[int, dict] = {a.block_id: {} for a in self.blocks}

        for r in self.event_log:
            ev = r.event
            if ev is EventType.ORIGINATION or ev is EventType.INSERTION:
                # ORIGINATION: the seed row per seed segment (gene / intergene tiling) or a novel
                # gene; INSERTION: a run of novel intergene nucleotides. Both root the blocks they
                # cover, so they enter the block records identically (an origin of new sequence).
                for op in r.genes:
                    e2 = prov.get(op.gid)
                    if e2 is None:
                        continue
                    s2, ss2, se2 = e2
                    rec = EventRecord(EventType.ORIGINATION, r.branch, r.time,
                                      [GeneOp(op.gid, s2, "origin")])
                    for a in self._covered(s2, ss2, se2):
                        records_by_block[a.block_id].append(rec)
                continue

            g0 = r.genes[0].gid
            entry = prov.get(g0)
            if entry is None:
                continue
            source, ss, se = entry
            blocks = self._covered(source, ss, se)
            if not blocks:
                continue
            if ev in (EventType.SPECIATION, EventType.DUPLICATION, EventType.TRANSFER):
                rep = self._top(g0, top_cache)
                rec = EventRecord(ev, r.branch, r.time,
                                  [GeneOp(rep, source, "parent"),
                                   *(GeneOp(op.gid, source, "child") for op in r.genes[1:])],
                                  donor=r.donor, recipient=r.recipient)
            elif ev is EventType.LOSS or ev is EventType.DELETION:
                # DELETION removes a run of intergene sequence — same as a LOSS for the block
                # genealogy (the covered intergene block copy terminates on this branch).
                if ev is EventType.LOSS and len(r.genes) == 2 and r.genes[1].role == "pseudogenized":
                    # pseudogenization was logged as a LOSS with a continuation row; rewrite it
                    # to a state-change edge (gene -> intergene) on the continuing lineage.
                    rec = EventRecord(EventType.PSEUDOGENIZATION, r.branch, r.time,
                                      [GeneOp(self._top(g0, top_cache), source, "parent"),
                                       GeneOp(r.genes[1].gid, source, "child")])
                else:
                    rec = EventRecord(EventType.LOSS, r.branch, r.time,
                                      [GeneOp(self._top(g0, top_cache), source, "lost")])
            else:
                continue  # inversion / transposition / translocation never re-mint a lineage
            for a in blocks:
                records_by_block[a.block_id].append(rec)

        for leaf, genome in self.leaf_genomes.items():
            name = leaf.name
            for seg in genome._iter_segments():
                blocks = self._covered(seg.source, seg.src_start, seg.src_end)
                if not blocks:
                    continue
                rep = self._top(seg.seg_id, top_cache)
                for a in blocks:
                    species_by_block[a.block_id][rep] = name

        return records_by_block, species_by_block

    def block_gene_trees(self) -> dict[int, tuple]:
        """``block_id -> (complete_newick, extant_newick)`` — one gene tree per block.

        Speciations, duplications and transfers are bifurcations, losses terminate; built
        with the shared :func:`~zombi2.reconciliation.build_gene_trees` so the Newick output
        matches the rest of ZOMBI2. ``extant_newick`` is ``None`` if nothing survives.
        """
        records_by_block, species_by_block = self._block_records()
        total_age = self.species_tree.total_age
        return {a.block_id: build_gene_trees(records_by_block[a.block_id],
                                            species_by_block[a.block_id], total_age)
                for a in self.blocks}

    # --- gene vs intergene partition (recover both tree sets) ---------------
    def gene_blocks(self) -> list:
        """The ancestrally-gene blocks (one per surviving gene copy lineage)."""
        return [a for a in self.blocks if a.kind == "gene"]

    def intergene_blocks(self) -> list:
        return [a for a in self.blocks if a.kind == "intergene"]

    def gene_trees(self) -> dict[int, tuple]:
        """``block_id -> (complete, extant)`` for the **gene** blocks only."""
        trees = self.block_gene_trees()
        return {aid: t for aid, t in trees.items() if self._block_by_id[aid].kind == "gene"}

    def intergene_trees(self) -> dict[int, tuple]:
        """``block_id -> (complete, extant)`` for the **intergene** blocks only."""
        trees = self.block_gene_trees()
        return {aid: t for aid, t in trees.items() if self._block_by_id[aid].kind == "intergene"}

    def pseudogenizations(self) -> list[tuple]:
        """``[(block_id, gene_id, species_branch, time, gene_lineage)]`` — each gene->intergene flip.

        Read off the per-block reconciliation events (the ``"G"`` rows), so it reports both where
        (species branch) and when a gene lost function while its sequence continued as intergene.
        """
        out: list[tuple] = []
        for block_id, rec in self.block_reconciliations().items():
            block = self._block_by_id[block_id]
            for e in rec.events:
                if e.event == "G":
                    out.append((block_id, block.gene_id, e.species, e.time, e.gene))
        return out

    def block_reconciliations(self) -> dict:
        """``block_id -> Reconciliation(complete, extant, events)`` — each block reconciled
        against the species tree.

        ``complete`` reconciles the complete gene tree (every event, including the real
        losses); ``extant`` reconciles the observable (pruned) gene tree — the cherries, no
        losses. See :func:`~zombi2.reconciliation.reconcile`.
        """
        records_by_block, species_by_block = self._block_records()
        total_age = self.species_tree.total_age
        return {a.block_id: reconcile(records_by_block[a.block_id],
                                     species_by_block[a.block_id], total_age)
                for a in self.blocks}

    def write_reconciliations(self, outdir) -> dict:
        """Write the reconciled trees + the events table to ``outdir``.

        ``reconciled_complete.nwk`` / ``reconciled_extant.nwk`` — one ``block_id<TAB>newick``
        line per block (the complete history with losses, and the observable cherries);
        ``reconciliation_events.tsv`` — the events (S/D/T/L) with their species location,
        transfer recipient, time and gene lineage. Returns a small summary dict.
        """
        out = Path(outdir)
        out.mkdir(parents=True, exist_ok=True)
        recon = self.block_reconciliations()
        complete_lines, extant_lines = [], []
        ev_lines, n_events = ["block\tevent\tspecies\trecipient\ttime\tgene"], 0
        for block_id, rec in recon.items():
            if rec.complete:
                complete_lines.append(f"{block_id}\t{rec.complete}")
            if rec.extant:
                extant_lines.append(f"{block_id}\t{rec.extant}")
            for e in rec.events:
                ev_lines.append(f"{block_id}\t{e.event}\t{e.species}\t{e.recipient or ''}\t"
                                f"{e.time:.10g}\t{e.gene or ''}")
                n_events += 1
        (out / "reconciled_complete.nwk").write_text("\n".join(complete_lines) + "\n")
        (out / "reconciled_extant.nwk").write_text("\n".join(extant_lines) + "\n")
        (out / "reconciliation_events.tsv").write_text("\n".join(ev_lines) + "\n")
        return {"path": str(out), "n_blocks": len(complete_lines), "n_events": n_events}

    # --- per-block history (step 7: the events that touched each segment) ---
    def block_histories(self) -> dict[int, list[tuple[str, float]]]:
        """``block_id -> [(branch, time), ...]`` inversions whose arc covered the block.

        Branch-tagged, so a consumer can restrict to a leaf's ancestral lineage.
        """
        out: dict[int, list[tuple[str, float]]] = {a.block_id: [] for a in self.blocks}
        for r in self.event_log:
            if r.event is not EventType.INVERSION:
                continue
            for op in r.genes:
                source, ss, se = self.registry.provenance[op.gid]
                for a in self._covered(source, ss, se):
                    out[a.block_id].append((r.branch, r.time))
        return out

    # --- ancestral genomes + sequences (needs retain_internal + simulate_sequences) -----
    def _chrom_mosaic(self, chrom) -> list[tuple[int, int]]:
        """One chromosome as an ordered, signed block sequence ``[(block_id, strand), ...]``."""
        out: list[tuple[int, int]] = []
        for seg in chrom.elements:
            covered = self._covered(seg.source, seg.src_start, seg.src_end)
            if seg.strand == -1:
                covered = list(reversed(covered))
            out.extend((a.block_id, seg.strand) for a in covered)
        return out

    def node_mosaics(self, node) -> dict:
        """Per-chromosome block mosaic at ``node``: ``{chrom_id: [(block_id, strand), ...]}``.

        Requires ``retain_internal`` (see :func:`simulate_nucleotide_genomes`).
        """
        return {cid: self._chrom_mosaic(chrom)
                for cid, chrom in self.node_genomes[node].chromosomes.items()}

    def node_mosaic(self, node) -> list[tuple[int, int]]:
        """The whole genome at ``node`` as one signed block sequence (every chromosome, in order).

        Generalises :meth:`leaf_mosaic` to any node; use :meth:`node_mosaics` to keep the
        chromosomes apart.
        """
        out: list[tuple[int, int]] = []
        for mosaic in self.node_mosaics(node).values():
            out.extend(mosaic)
        return out

    def _seed_source(self) -> str:
        """Source id of the seed (input) chromosome — the blocks that map to the real genome/FASTA."""
        g = self.node_genomes.get(self.species_tree.root)
        first = next(iter(g._iter_segments()), None) if g is not None else None
        if first is not None:
            return first.source
        by: dict = {}
        for a in self.blocks:
            by[a.source] = by.get(a.source, 0) + 1
        return max(by, key=by.get) if by else "1"

    def simulate_sequences(self, model, *, gamma=None, root_fasta=None, subst_rate: float = 1.0,
                           clock=None, rng=None, seed: int | None = None) -> dict:
        """Evolve a DNA sequence for every block lineage; cache + return ``{block_id: {gid: seq}}``.

        Each block's *complete* gene tree is scaled to substitutions/site (a strict clock by default,
        or a supplied :class:`~zombi2.SequenceEvolution` ``clock``; ``subst_rate`` scales the overall
        divergence) and a sequence is evolved down it under ``model`` (optionally with across-site
        :class:`~zombi2.sequence_sim.GammaRates`). A seed-chromosome block takes its root sequence from
        ``root_fasta`` (the real genome) when given, else a random root of the block's length.
        ``root_fasta`` is either a single string (seeds the main chromosome, keyed by ``root_length``)
        or a ``{source: dna}`` map — one entry per replicon — for a multi-chromosome genome.
        :meth:`node_sequence` then assembles these into the DNA at any node.
        """
        from zombi2.genomes.reconciliation import _node_tree
        from zombi2.sequences.evolution import SequenceEvolution, _annotate
        from zombi2.sequences.models import evolve_on_tree
        if rng is None:
            rng = np.random.default_rng(seed)
        zero = clock is None and subst_rate <= 0.0   # no substitutions: root propagates unchanged
        se = clock or (None if zero else SequenceEvolution(root_rate=subst_rate))
        segments = None if zero else se._lineage_segments(self.species_tree, rng)[0]
        total_age = self.species_tree.total_age
        records_by_block, species_by_block = self._block_records()
        # normalise root_fasta to a {source: dna} map: a bare string seeds the main chromosome.
        if isinstance(root_fasta, str):
            if len(root_fasta) != self.root_length:
                raise ValueError(f"root_fasta length {len(root_fasta)} != root_length {self.root_length}")
            root_seqs = {self._seed_source(): root_fasta}
        else:
            root_seqs = root_fasta or {}

        block_seqs: dict[int, dict] = {}
        for a in self.blocks:
            root_node = _node_tree(records_by_block[a.block_id], species_by_block[a.block_id], total_age)
            if root_node is None:
                block_seqs[a.block_id] = {}
                continue
            subst: dict = {}
            if not zero:
                _annotate(root_node, segments, max(0.0, se.family_speed.sample(rng)), subst)
            src_dna = root_seqs.get(a.source)
            root_seq = src_dna[a.start:a.end] if src_dna is not None else None
            block_seqs[a.block_id] = evolve_on_tree(root_node, subst, model, rng,
                                                  root_seq=root_seq, length=a.length, gamma=gamma)
        self._block_seqs = block_seqs
        self._seq_species = species_by_block
        return block_seqs

    def _chrom_sequence(self, chrom, block_seqs, top_cache, node) -> str:
        """Assemble one chromosome's DNA from its blocks' evolved sequences."""
        from zombi2.sequences.models import reverse_complement
        parts: list[str] = []
        for seg in chrom.elements:
            rep = self._top(seg.seg_id, top_cache)
            covered = self._covered(seg.source, seg.src_start, seg.src_end)
            if seg.strand == -1:
                covered = list(reversed(covered))
            for block in covered:
                s = block_seqs.get(block.block_id, {}).get(rep)
                if s is None:
                    raise KeyError(f"no sequence for block {block.block_id} lineage {rep!r} "
                                   f"at node {getattr(node, 'name', node)!r}")
                parts.append(reverse_complement(s) if seg.strand == -1 else s)
        return "".join(parts)

    def node_sequences(self, node) -> dict:
        """Per-chromosome assembled DNA at ``node``: ``{chrom_id: dna}`` (call
        :meth:`simulate_sequences` first). One entry per replicon — the natural output for a
        multi-chromosome genome (a chromosome plus its plasmids)."""
        block_seqs = getattr(self, "_block_seqs", None)
        if block_seqs is None:
            raise RuntimeError("call simulate_sequences(...) before node_sequences(...)")
        top_cache: dict = {}
        return {cid: self._chrom_sequence(chrom, block_seqs, top_cache, node)
                for cid, chrom in self.node_genomes[node].chromosomes.items()}

    def node_sequence(self, node) -> str:
        """Assemble the full DNA of the genome at ``node`` — every chromosome, concatenated in order
        (call :meth:`simulate_sequences` first). The root node reproduces the input genome; extant
        leaves give the observed genomes. Use :meth:`node_sequences` to keep the chromosomes apart.
        """
        return "".join(self.node_sequences(node).values())

    def gene_alignments(self) -> dict:
        """``{gene_id: {species_gid: seq}}`` extant alignments for gene blocks (needs sequences)."""
        return self._alignments("gene")

    def intergene_alignments(self) -> dict:
        """``{blockN: {species_gid: seq}}`` extant alignments for intergene blocks."""
        return self._alignments("intergene")

    def _alignments(self, kind: str) -> dict:
        block_seqs = getattr(self, "_block_seqs", None)
        if block_seqs is None:
            raise RuntimeError("call simulate_sequences(...) first")
        out: dict = {}
        for a in self.blocks:
            if a.kind != kind:
                continue
            g2s = self._seq_species.get(a.block_id, {})
            seqs = block_seqs.get(a.block_id, {})
            aln = {f"{g2s[gid]}_{gid}": seqs[gid] for gid in g2s if gid in seqs}
            if aln:
                out[a.gene_id if (kind == "gene" and a.gene_id) else f"block{a.block_id}"] = aln
        return out


def _merge(intervals: list[tuple[int, int]]) -> list[tuple[int, int]]:
    """Merge a list of ``[start, end)`` intervals into disjoint covered spans."""
    merged: list[list[int]] = []
    for lo, hi in sorted(intervals):
        if merged and lo <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], hi)
        else:
            merged.append([lo, hi])
    return [(lo, hi) for lo, hi in merged]


def _build_blocks(leaf_genomes: dict, root_length: int, registry=None) -> list[Block]:
    """Partition each source into blocks at the union of all extant-leaf breakpoints.

    With deletion, some ancestral positions survive in no extant leaf; those gaps carry
    no block. So an interval between consecutive breakpoints becomes a block only if some
    leaf still covers it. When a gene annotation is present (``registry``), gene boundaries
    are seeded as breakpoints — genes are their own segments, so this only makes classification
    robust — and each block is tagged ``"gene"``/``"intergene"`` by ancestral coordinate.
    """
    bounds: dict[str, set[int]] = {}
    spans: dict[str, list[tuple[int, int]]] = {}
    for genome in leaf_genomes.values():
        for seg in genome._iter_segments():
            bounds.setdefault(seg.source, {0}).add(seg.src_start)
            bounds[seg.source].add(seg.src_end)
            spans.setdefault(seg.source, []).append((seg.src_start, seg.src_end))
    if registry is not None:
        for source, gis in registry.genes.items():
            for gi in gis:
                bounds.setdefault(source, {0}).update((gi.start, gi.end))
    blocks: list[Block] = []
    aid = 0
    for source in sorted(bounds):
        covered = _merge(spans.get(source, []))
        cuts = sorted(bounds[source])
        for a, b in zip(cuts, cuts[1:]):
            if any(lo <= a and b <= hi for lo, hi in covered):  # surviving material only
                gid = registry.gene_id_at(source, a, b) if registry is not None else None
                blocks.append(Block(aid, source, a, b, "gene" if gid is not None else "intergene", gid))
                aid += 1
    return blocks


def _normalize_gene_intervals(gene_intervals, root_length: int):
    """Validate + normalise user gene intervals to sorted ``[(start, end, name|None)]``.

    Genes must be non-overlapping half-open intervals inside ``[0, root_length)``. Returns an
    empty list when no genes are supplied (the plain nucleotide model).
    """
    if not gene_intervals:
        return []
    out = []
    for gi in gene_intervals:
        if len(gi) == 2:
            a, b, name = gi[0], gi[1], None
        elif len(gi) == 3:
            a, b, name = gi
        else:
            raise ValueError(f"gene interval must be (start, end) or (start, end, name), got {gi!r}")
        a, b = int(a), int(b)
        if not (0 <= a < b <= root_length):
            raise ValueError(f"gene interval {gi!r} must satisfy 0 <= start < end <= "
                             f"root_length ({root_length})")
        out.append((a, b, str(name) if name is not None else None))
    out.sort()
    for (a1, b1, _), (a2, b2, _) in zip(out, out[1:]):
        if a2 < b1:
            raise ValueError(f"gene intervals must not overlap: [{a1},{b1}) and [{a2},{b2})")
    return out


def simulate_nucleotide_genomes(
    species_tree: Tree,
    *,
    inversion: float = 0.001,
    loss: float = 0.0,
    duplication: float = 0.0,
    transfer: float = 0.0,
    transposition: float = 0.0,
    translocation: float = 0.0,
    origination: float = 0.0,
    insertion: float = 0.0,
    deletion: float = 0.0,
    fission: float = 0.0,
    fusion: float = 0.0,
    chromosome_origination: float = 0.0,
    chromosome_loss: float = 0.0,
    indel_mean_length: float = 10.0,
    root_length: int = 1000,
    extension: float | None = 0.99,
    initial_chromosomes: int = 1,
    root_chromosomes: list[tuple] | None = None,
    circular: bool = True,
    transfers=None,
    gene_intervals=None,
    pseudogenization: float = 0.0,
    replacement: float = 0.0,
    retain_internal: bool = False,
    seed: int | None = None,
    rng: np.random.Generator | None = None,
    sampler: EventSampler | None = None,
    max_segments_per_genome: int | None = DEFAULT_MAX_SEGMENTS_PER_GENOME,
    output: str = "genomes",
) -> NucleotideResult:
    """Simulate variable-length structural events forward along ``species_tree``.

    ``inversion``, ``loss``, ``duplication``, ``transfer`` and ``transposition`` are
    **per-nucleotide** rates (total genome rate = ``rate * current_length``);
    ``origination`` is a **per-branch** rate that inserts a novel gene under a fresh source.
    ``translocation`` (per-nucleotide, needs >= 2 chromosomes) *moves* an arc to a different
    chromosome of the same genome — the intra-genome counterpart of transposition.
    ``extension`` sets the geometric event-length model (mean ``1/(1-extension)``
    nucleotides). ``initial_chromosomes`` is the number of root chromosomes seeded at the
    root of the tree (default 1); each is an independent full-length copy of the root chromosome
    under its own source namespace (in genic mode all copies share the same gene layout), so an
    ``N``-chromosome genome starts at ``N * root_length`` bp. Multi-chromosome genomes also arise
    dynamically through the chromosome-tier events below. Chromosomes are circular by default;
    ``circular=False`` makes every seeded chromosome **linear** (two ends, no origin wrap; requires
    ``output="genomes"``). For **heterogeneous** root chromosomes — a real chromosome-plus-plasmids
    genome, each replicon its own length, genes and topology — pass ``root_chromosomes`` instead: a
    list of ``(length, gene_intervals)`` or ``(length, gene_intervals, circular)`` (e.g. from
    :func:`~zombi2.read_gff_all` over a multi-sequence GFF, carrying each sequence's ``Is_circular``
    flag — so a *Borrelia*-style linear chromosome + circular plasmids seeds a mixed-topology genome).
    It is mutually exclusive with ``gene_intervals`` / ``initial_chromosomes`` and requires
    ``output="genomes"``.
    ``transfers`` is an optional :class:`~zombi2.TransferModel` (default:
    additive, uniform recipient, no self-transfer). Returns a :class:`NucleotideResult`
    carrying the extant leaf genomes, the event log, the segment registry, and the block
    partition (over the surviving ancestral material).

    ``insertion`` and ``deletion`` are **per-nucleotide** rates for intergenic indels: an
    insertion lays down a run of novel nucleotides (a fresh source, its own block) inside an
    intergene stretch, lengthening it; a deletion removes a run from *within a single* intergene
    stretch, clamped so it never reaches into a gene and never shrinks the chromosome below a
    small floor (:data:`~zombi2.nucleotide_genome.MIN_GENOME_LENGTH`). Both default to 0 (off).
    In genic mode indels act only in intergenes (a gene is never split, spanned or deleted); with
    no genes declared they act anywhere. ``indel_mean_length`` is the mean of the indel run's
    geometric length (a separate parameter from ``extension``; default 10).

    Genes & intergenes (genic mode): pass ``gene_intervals`` — a list of non-overlapping
    ``(start, end)`` or ``(start, end, name)`` intervals on the root chromosome. Event
    breakpoints then only fall in intergene positions, so genes are never split; each gene is
    exactly one block (``kind="gene"``) and its own genealogy, while intergene stretches
    decompose into ``kind="intergene"`` blocks. ``pseudogenization`` in ``[0, 1]`` is the
    probability that a loss hitting a gene *demotes* it to intergene (sequence retained, a
    ``G`` state-change edge in its tree) rather than deleting it. ``replacement`` in ``[0, 1]``
    is the probability that a transfer is a *homologous* replacement: the copy replaces the
    recipient's syntenic material (located via flanking genes of matching identity **and**
    orientation), falling back to additive insertion when the recipient has no such homolog.
    See :meth:`NucleotideResult.gene_trees` /
    :meth:`~NucleotideResult.intergene_trees` / :meth:`~NucleotideResult.pseudogenizations`.

    Chromosome-tier events (all **per-branch** rates, default 0 → a single chromosome end to
    end, byte-identical to the pre-tier engine): ``fission`` splits a chromosome into two by
    excising the arc between two breakpoints into a new circular replicon (breakpoints snap to
    segment boundaries, so genes stay whole and segment ids are preserved); ``fusion`` merges two
    chromosomes into one; ``chromosome_origination`` spawns a new empty circular replicon (a
    de-novo plasmid); ``chromosome_loss`` deletes a whole chromosome and every gene on it. These
    move genes between chromosomes without touching their identity, so per-block gene trees are
    unaffected — only chromosome loss ends lineages. The karyotype history is recorded in
    ``result.event_log.chromosome_records`` (one :class:`ChromosomeEvent` per fission / fusion /
    origination / loss). Requires ``output="genomes"`` (the Rust profiles path is single-chromosome).

    Duplication and additive transfer grow the genome with no cap, so keep them at or below
    ``loss`` over long ages to avoid runaway growth. Because a structural event's rate is
    proportional to genome length, an unbalanced gain rate compounds and the walk would never
    terminate; ``max_segments_per_genome`` (default
    :data:`DEFAULT_MAX_SEGMENTS_PER_GENOME`) is a safety ceiling that turns such a run into a
    prompt, actionable ``RuntimeError`` instead of an endless hang. It is far above any bounded
    (balanced) run, so it never fires in normal use and leaves seeded outputs unchanged; pass a
    larger value (or ``None`` to disable) only if a genome is genuinely meant to grow that big.

    ``output``:
        ``"genomes"`` (default) runs the pure-Python engine and returns the full result
        (event log, per-block gene trees and histories). ``"profiles"`` runs the compiled
        ``zombi2_core`` Rust engine over leaf segments only — much faster, and enough for
        ``profile_matrix()`` / ``leaf_mosaic()`` / ``trace_back()`` — but emits no event log,
        so ``block_gene_trees()`` / ``block_histories()`` are unavailable, it **requires**
        the extension, and it does **not** support the genic model (``gene_intervals`` /
        ``pseudogenization``).
    """
    if output not in ("genomes", "profiles"):
        raise ValueError(f"output must be 'genomes' or 'profiles', got {output!r}")
    if not (0.0 <= pseudogenization <= 1.0):
        raise ValueError(f"pseudogenization must be in [0, 1], got {pseudogenization}")
    if not (0.0 <= replacement <= 1.0):
        raise ValueError(f"replacement must be in [0, 1], got {replacement}")
    pending_genes = _normalize_gene_intervals(gene_intervals, root_length)
    if root_chromosomes is not None:
        # explicit heterogeneous replicons (e.g. a multi-sequence GFF): each is its own chromosome,
        # ``(length, gene_intervals)`` or ``(length, gene_intervals, circular)`` — so mixed
        # circular/linear genomes are possible. Mutually exclusive with the identical-copy knobs.
        if gene_intervals is not None:
            raise ValueError("give either gene_intervals or root_chromosomes, not both")
        if initial_chromosomes != 1:
            raise ValueError("root_chromosomes sets the chromosomes explicitly; "
                             "leave initial_chromosomes=1")
        if not root_chromosomes:
            raise ValueError("root_chromosomes must be a non-empty list of "
                             "(length, gene_intervals[, circular])")
        root_chromosomes = [(int(spec[0]), _normalize_gene_intervals(spec[1], int(spec[0])),
                             bool(spec[2]) if len(spec) > 2 else circular)
                            for spec in root_chromosomes]
    if output == "profiles":
        if pending_genes or root_chromosomes is not None:
            raise ValueError("gene intervals / explicit root_chromosomes require the Python engine "
                             "(output='genomes'); the Rust profiles path does not model "
                             "genes/intergenes or heterogeneous chromosomes")
        if pseudogenization:
            raise ValueError("pseudogenization requires output='genomes' (the Python engine)")
        if insertion or deletion:
            raise ValueError("intergenic indels (insertion/deletion) require output='genomes' "
                             "(the Python engine); the Rust profiles path does not model them")
        if fission or fusion or chromosome_origination or chromosome_loss or translocation:
            raise ValueError("chromosome-tier events (fission/fusion/chromosome_origination/"
                             "chromosome_loss) and translocation require output='genomes' (the "
                             "Python engine); the Rust profiles path is single-chromosome")
        if not circular:
            raise ValueError("linear chromosomes require output='genomes' (the Python engine); the "
                             "Rust profiles path models a circular chromosome only")
        if sampler is not None:
            raise ValueError("output='profiles' uses the Rust engine and ignores a custom sampler")
        if extension is None:
            raise ValueError("output='profiles' requires a numeric `extension` (the geometric "
                             "event-length parameter); the None mode is Python-engine only")
        if seed is None and rng is not None:
            seed = int(rng.integers(0, 2**63 - 1))
        from zombi2 import _rust
        return _rust.nucleotide(
            species_tree, inversion=inversion, loss=loss, duplication=duplication,
            transfer=transfer, transposition=transposition, origination=origination,
            root_length=root_length, extension=extension, initial_size=initial_chromosomes,
            transfers=transfers, seed=seed,
        )

    if rng is None:
        rng = np.random.default_rng(seed)
    rates = Rates(inversion=inversion, loss=loss, duplication=duplication,
                         transfer=transfer, transposition=transposition,
                         translocation=translocation,
                         insertion=insertion, deletion=deletion,
                         origination=origination, fission=fission, fusion=fusion,
                         chromosome_origination=chromosome_origination,
                         chromosome_loss=chromosome_loss)
    registry = SegmentRegistry(pending_genes=pending_genes)

    def factory(ids):
        return NucleotideGenome(ids, root_length=root_length, extension=extension,
                                registry=registry, pseudogenization=pseudogenization,
                                replacement=replacement, indel_mean_length=indel_mean_length,
                                initial_chromosomes=initial_chromosomes,
                                root_chromosomes=root_chromosomes, circular=circular)

    # One seed origination lays down all `initial_chromosomes` root chromosomes (the genome owns
    # the count); the walk then evolves them and fires any further per-branch originations.
    result = GenomeSimulator(
        sampler, max_segments_per_genome=max_segments_per_genome,
    ).simulate(
        species_tree, rates, rng, initial_size=1, transfers=transfers,
        genome_factory=factory, retain_internal=retain_internal,
    )
    # For ancestral reconstruction, blocks must tile every node's segments (not just the leaves'):
    # build from all node genomes so internal breakpoints and ancestral-only material are covered.
    block_genomes = result.node_genomes if retain_internal else result.leaf_genomes
    blocks = _build_blocks(block_genomes, root_length, registry)
    # with explicit heterogeneous replicons the root has no single length; report the genome total
    # (only used to validate a supplied root_fasta, which is a single-chromosome feature).
    effective_root_length = (sum(spec[0] for spec in root_chromosomes)
                             if root_chromosomes is not None else root_length)
    return NucleotideResult(
        species_tree=species_tree,
        leaf_genomes=result.leaf_genomes,
        event_log=result.event_log,
        registry=registry,
        blocks=blocks,
        root_length=effective_root_length,
        node_genomes=result.node_genomes,
    )
