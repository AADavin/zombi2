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

from ._sampling import EventSampler
from .events import EventType, EventRecord, GeneOp
from .genome_sim import GenomeSimulator
from .nucleotide_genome import NucleotideGenome, SegmentRegistry
from .rates import UniformRates
from .reconciliation import build_gene_trees, reconcile
from .tree import Tree, TreeNode


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
        for seg in self.leaf_genomes[leaf]._segments:
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
            for seg in self.leaf_genomes[leaf]._segments:
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
            if ev is EventType.ORIGINATION:
                # the seed origination has one row per seed segment (gene / intergene tiling),
                # each the root of the blocks it covers; a novel origination has a single row.
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
            elif ev is EventType.LOSS:
                if len(r.genes) == 2 and r.genes[1].role == "pseudogenized":
                    # pseudogenization was logged as a LOSS with a continuation row; rewrite it
                    # to a state-change edge (gene -> intergene) on the continuing lineage.
                    rec = EventRecord(EventType.PSEUDOGENIZATION, r.branch, r.time,
                                      [GeneOp(self._top(g0, top_cache), source, "parent"),
                                       GeneOp(r.genes[1].gid, source, "child")])
                else:
                    rec = EventRecord(EventType.LOSS, r.branch, r.time,
                                      [GeneOp(self._top(g0, top_cache), source, "lost")])
            else:
                continue  # inversion / transposition never re-mint a lineage
            for a in blocks:
                records_by_block[a.block_id].append(rec)

        for leaf, genome in self.leaf_genomes.items():
            name = leaf.name
            for seg in genome._segments:
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

        ``Reconciled_complete.nwk`` / ``Reconciled_extant.nwk`` — one ``block_id<TAB>newick``
        line per block (the complete history with losses, and the observable cherries);
        ``Reconciliation_events.tsv`` — the events (S/D/T/L) with their species location,
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
        (out / "Reconciled_complete.nwk").write_text("\n".join(complete_lines) + "\n")
        (out / "Reconciled_extant.nwk").write_text("\n".join(extant_lines) + "\n")
        (out / "Reconciliation_events.tsv").write_text("\n".join(ev_lines) + "\n")
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
    def node_mosaic(self, node) -> list[tuple[int, int]]:
        """The genome at ``node`` as an ordered, signed block sequence ``[(block_id, strand), ...]``.

        Generalises :meth:`leaf_mosaic` to any node; the root's mosaic is the input genome's
        gene/intergene tiling. Requires ``retain_internal`` (see :func:`simulate_nucleotide_genomes`).
        """
        out: list[tuple[int, int]] = []
        for seg in self.node_genomes[node]._segments:
            covered = self._covered(seg.source, seg.src_start, seg.src_end)
            if seg.strand == -1:
                covered = list(reversed(covered))
            out.extend((a.block_id, seg.strand) for a in covered)
        return out

    def _seed_source(self) -> str:
        """Source id of the seed (input) chromosome — the blocks that map to the real genome/FASTA."""
        g = self.node_genomes.get(self.species_tree.root)
        if g is not None and g._segments:
            return g._segments[0].source
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
        :meth:`node_sequence` then assembles these into the DNA at any node.
        """
        from .reconciliation import _node_tree
        from .sequence_evolution import SequenceEvolution, _annotate
        from .sequence_sim import evolve_on_tree
        if rng is None:
            rng = np.random.default_rng(seed)
        zero = clock is None and subst_rate <= 0.0   # no substitutions: root propagates unchanged
        se = clock or (None if zero else SequenceEvolution(root_rate=subst_rate))
        segments = None if zero else se._lineage_segments(self.species_tree, rng)[0]
        total_age = self.species_tree.total_age
        records_by_block, species_by_block = self._block_records()
        seed_source = self._seed_source()
        if root_fasta is not None and len(root_fasta) != self.root_length:
            raise ValueError(f"root_fasta length {len(root_fasta)} != root_length {self.root_length}")

        block_seqs: dict[int, dict] = {}
        for a in self.blocks:
            root_node = _node_tree(records_by_block[a.block_id], species_by_block[a.block_id], total_age)
            if root_node is None:
                block_seqs[a.block_id] = {}
                continue
            subst: dict = {}
            if not zero:
                _annotate(root_node, segments, max(0.0, se.family_speed.sample(rng)), subst)
            root_seq = (root_fasta[a.start:a.end]
                        if (root_fasta is not None and a.source == seed_source) else None)
            block_seqs[a.block_id] = evolve_on_tree(root_node, subst, model, rng,
                                                  root_seq=root_seq, length=a.length, gamma=gamma)
        self._block_seqs = block_seqs
        self._seq_species = species_by_block
        return block_seqs

    def node_sequence(self, node) -> str:
        """Assemble the full DNA of the genome at ``node`` (call :meth:`simulate_sequences` first).

        Concatenates, in genome order, each segment's block-lineage sequences (reverse-complemented
        on the − strand). The root node reproduces the input genome; extant leaves give the observed
        genomes.
        """
        from .sequence_sim import reverse_complement
        block_seqs = getattr(self, "_block_seqs", None)
        if block_seqs is None:
            raise RuntimeError("call simulate_sequences(...) before node_sequence(...)")
        top_cache: dict = {}
        parts: list[str] = []
        for seg in self.node_genomes[node]._segments:
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
        for seg in genome._segments:
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
    origination: float = 0.0,
    root_length: int = 1000,
    extension: float | None = 0.99,
    initial_chromosomes: int = 1,
    transfers=None,
    gene_intervals=None,
    pseudogenization: float = 0.0,
    replacement: float = 0.0,
    retain_internal: bool = False,
    seed: int | None = None,
    rng: np.random.Generator | None = None,
    sampler: EventSampler | None = None,
    output: str = "genomes",
) -> NucleotideResult:
    """Simulate variable-length structural events forward along ``species_tree``.

    ``inversion``, ``loss``, ``duplication``, ``transfer`` and ``transposition`` are
    **per-nucleotide** rates (total genome rate = ``rate * current_length``);
    ``origination`` is a **per-branch** rate that inserts a novel gene under a fresh source.
    ``extension`` sets the geometric event-length model (mean ``1/(1-extension)``
    nucleotides). ``initial_chromosomes`` is the number of root chromosomes seeded at the
    root of the tree (default 1); each is an independent copy of the root chromosome.
    ``transfers`` is an optional :class:`~zombi2.TransferModel` (default:
    additive, uniform recipient, no self-transfer). Returns a :class:`NucleotideResult`
    carrying the extant leaf genomes, the event log, the segment registry, and the block
    partition (over the surviving ancestral material).

    Genes & intergenes (genic mode): pass ``gene_intervals`` — a list of non-overlapping
    ``(start, end)`` or ``(start, end, name)`` intervals on the root chromosome. Event
    breakpoints then only fall in intergene positions, so genes are never split; each gene is
    exactly one block (``kind="gene"``) and its own genealogy, while intergene stretches
    decompose into ``kind="intergene"`` blocks. ``pseudogenization`` in ``[0, 1]`` is the
    probability that a loss hitting a gene *demotes* it to intergene (sequence retained, a
    ``G`` state-change edge in its tree) rather than deleting it. ``replacement`` in ``[0, 1]``
    is the probability that a transfer is a *homologous* replacement: the copy replaces the
    recipient's syntenic material (located via flanking genes), falling back to additive
    insertion when the recipient has no homolog. See :meth:`NucleotideResult.gene_trees` /
    :meth:`~NucleotideResult.intergene_trees` / :meth:`~NucleotideResult.pseudogenizations`.

    Duplication and additive transfer grow the genome with no cap, so keep them at or below
    ``loss`` over long ages to avoid runaway growth.

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
    if output == "profiles":
        if pending_genes:
            raise ValueError("gene intervals require the Python engine (output='genomes'); the "
                             "Rust profiles path does not model genes/intergenes")
        if pseudogenization:
            raise ValueError("pseudogenization requires output='genomes' (the Python engine)")
        if sampler is not None:
            raise ValueError("output='profiles' uses the Rust engine and ignores a custom sampler")
        if extension is None:
            raise ValueError("output='profiles' requires a numeric `extension` (the geometric "
                             "event-length parameter); the None mode is Python-engine only")
        if seed is None and rng is not None:
            seed = int(rng.integers(0, 2**63 - 1))
        from . import _rust
        return _rust.nucleotide(
            species_tree, inversion=inversion, loss=loss, duplication=duplication,
            transfer=transfer, transposition=transposition, origination=origination,
            root_length=root_length, extension=extension, initial_size=initial_chromosomes,
            transfers=transfers, seed=seed,
        )

    if rng is None:
        rng = np.random.default_rng(seed)
    rates = UniformRates(inversion=inversion, loss=loss, duplication=duplication,
                         transfer=transfer, transposition=transposition,
                         origination=origination)
    registry = SegmentRegistry(pending_genes=pending_genes)

    def factory(ids):
        return NucleotideGenome(ids, root_length=root_length, extension=extension,
                                registry=registry, pseudogenization=pseudogenization,
                                replacement=replacement)

    result = GenomeSimulator(sampler).simulate(
        species_tree, rates, rng, initial_size=initial_chromosomes, transfers=transfers,
        genome_factory=factory, retain_internal=retain_internal,
    )
    # For ancestral reconstruction, blocks must tile every node's segments (not just the leaves'):
    # build from all node genomes so internal breakpoints and ancestral-only material are covered.
    block_genomes = result.node_genomes if retain_internal else result.leaf_genomes
    blocks = _build_blocks(block_genomes, root_length, registry)
    return NucleotideResult(
        species_tree=species_tree,
        leaf_genomes=result.leaf_genomes,
        event_log=result.event_log,
        registry=registry,
        blocks=blocks,
        root_length=root_length,
        node_genomes=result.node_genomes,
    )
