"""Forward simulation + trace-back post-processing for the nucleotide genome (M1).

Runs the *unchanged* :class:`~zombi2.genome_sim.GenomeSimulator` over
:class:`~zombi2.nucleotide_genome.NucleotideGenome` with a genome-level inversion rate,
then decomposes the result:

* **atoms** — the finest intervals of the ancestral genome that are never cut by a
  breakpoint in *any* extant leaf. Each atom is a "segment with one shared history".
* **trace-back** — for each extant leaf, the ancestral origin of every nucleotide (this
  is just the leaf genome's :meth:`~zombi2.nucleotide_genome.NucleotideGenome.to_cells`).
* **mosaic** — each extant genome as an ordered, signed sequence of atoms.
* **atom histories** — the inversions (branch, time) that touched each atom.

For inversion-only M1 every atom survives in every leaf exactly once (inversion neither
creates nor destroys sequence), so an atom's genealogy is simply the species tree;
genuine per-atom gene trees arrive with duplication/loss/transfer in later milestones.
"""

from __future__ import annotations

from bisect import bisect_left
from dataclasses import dataclass, field

import numpy as np

from ._sampling import EventSampler
from .events import EventType, EventRecord, GeneOp
from .genome_sim import GenomeSimulator
from .nucleotide_genome import NucleotideGenome, SegmentRegistry
from .rates import UniformRates
from .reconciliation import build_gene_trees
from .tree import Tree, TreeNode


@dataclass(frozen=True)
class Atom:
    """A maximal uncut interval ``[start, end)`` of one ancestral ``source``."""

    atom_id: int
    source: str
    start: int
    end: int

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
    atoms: list         # list[Atom], tiling every source's [0, len)
    root_length: int
    _by_source: dict = field(default_factory=dict, repr=False)

    def __post_init__(self):
        by_source: dict[str, list[Atom]] = {}
        for a in self.atoms:
            by_source.setdefault(a.source, []).append(a)
        for atoms in by_source.values():
            atoms.sort(key=lambda a: a.start)
        self._by_source = by_source

    # --- per-leaf views ----------------------------------------------------
    def trace_back(self, leaf: TreeNode) -> list[tuple[str, int, int]]:
        """Ancestral origin ``(source, src_pos, strand)`` of every nucleotide at ``leaf``."""
        return self.leaf_genomes[leaf].to_cells()

    def leaf_mosaic(self, leaf: TreeNode) -> list[tuple[int, int]]:
        """The leaf genome as an ordered, signed sequence of atoms: ``[(atom_id, strand)]``."""
        out: list[tuple[int, int]] = []
        for seg in self.leaf_genomes[leaf]._segments:
            covered = [a for a in self._by_source[seg.source]
                       if a.start >= seg.src_start and a.end <= seg.src_end]
            if seg.strand == -1:
                covered = list(reversed(covered))
            out.extend((a.atom_id, seg.strand) for a in covered)
        return out

    # --- emergent phylogenetic profile (atoms as gene families) ------------
    def profile_matrix(self):
        """``(atom_ids, species, matrix)`` — copy number of each atom per extant leaf.

        Atoms are the emergent "gene families": a maximal block with one shared history.
        With loss only, entries are 0/1 (presence); duplication will lift them above 1.
        This is the phylogenetic-profile dataset the atom decomposition produces for free.
        """
        leaves = sorted(self.leaf_genomes, key=lambda n: n.name)
        atom_ids = [a.atom_id for a in self.atoms]
        row = {a.atom_id: i for i, a in enumerate(self.atoms)}
        matrix = np.zeros((len(self.atoms), len(leaves)), dtype=int)
        for j, leaf in enumerate(leaves):
            for seg in self.leaf_genomes[leaf]._segments:
                for a in self._by_source.get(seg.source, ()):
                    if a.start >= seg.src_start and a.end <= seg.src_end:
                        matrix[row[a.atom_id], j] += 1
        return atom_ids, [n.name for n in leaves], matrix

    # --- per-atom gene trees (steps 6-7: reconstruct the gene of each segment) ---
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

    def atom_gene_trees(self) -> dict[int, tuple]:
        """``atom_id -> (complete_newick, extant_newick)`` — one gene tree per atom.

        Each atom (a never-cut ancestral block) is an emergent gene family. Its genealogy
        is the segment lineage tree restricted to segments covering the atom: speciations,
        duplications and transfers are bifurcations, losses terminate, and breakpoint
        splits are contracted away (each split sends the atom to exactly one child). We
        collapse the split-chains into single lineages and hand the real events to the
        shared :func:`~zombi2.reconciliation.build_gene_trees`, so the reconstruction and
        Newick output match the rest of ZOMBI2. ``extant_newick`` is ``None`` if nothing
        survives.

        The log is indexed **once**: every event's segment covers a contiguous run of
        atoms (source coordinates are ordered), found by bisection, so each record is
        appended straight to those atoms — no per-atom rescans of the whole log.
        """
        prov = self.registry.provenance
        top_cache: dict[str, str] = {}
        total_age = self.species_tree.total_age

        # per-source atom lists sorted by start, with parallel start arrays for bisection
        starts = {src: [a.start for a in atoms] for src, atoms in self._by_source.items()}

        def covered(source: str, ss: int, se: int):
            atoms = self._by_source.get(source)
            if not atoms:
                return ()
            st = starts[source]
            return atoms[bisect_left(st, ss):bisect_left(st, se)]  # start in [ss, se)

        records_by_atom: dict[int, list] = {a.atom_id: [] for a in self.atoms}
        species_by_atom: dict[int, dict] = {a.atom_id: {} for a in self.atoms}

        for r in self.event_log:
            ev = r.event
            g0 = r.genes[0].gid
            entry = prov.get(g0)
            if entry is None:
                continue
            source, ss, se = entry
            atoms = covered(source, ss, se)
            if not atoms:
                continue
            if ev is EventType.ORIGINATION:
                rec = r  # gid == root == its own lineage rep
            elif ev in (EventType.SPECIATION, EventType.DUPLICATION, EventType.TRANSFER):
                rep = self._top(g0, top_cache)
                rec = EventRecord(ev, r.branch, r.time,
                                  [GeneOp(rep, source, "parent"),
                                   *(GeneOp(op.gid, source, "child") for op in r.genes[1:])])
            elif ev is EventType.LOSS:
                rec = EventRecord(ev, r.branch, r.time,
                                  [GeneOp(self._top(g0, top_cache), source, "lost")])
            else:
                continue  # inversion / transposition never re-mint a lineage
            for a in atoms:
                records_by_atom[a.atom_id].append(rec)

        for leaf, genome in self.leaf_genomes.items():
            name = leaf.name
            for seg in genome._segments:
                atoms = covered(seg.source, seg.src_start, seg.src_end)
                if not atoms:
                    continue
                rep = self._top(seg.seg_id, top_cache)
                for a in atoms:
                    species_by_atom[a.atom_id][rep] = name

        return {a.atom_id: build_gene_trees(records_by_atom[a.atom_id],
                                            species_by_atom[a.atom_id], total_age)
                for a in self.atoms}

    # --- per-atom history (step 7: the events that touched each segment) ---
    def atom_histories(self) -> dict[int, list[tuple[str, float]]]:
        """``atom_id -> [(branch, time), ...]`` inversions whose arc covered the atom.

        Branch-tagged, so a consumer can restrict to a leaf's ancestral lineage.
        """
        out: dict[int, list[tuple[str, float]]] = {a.atom_id: [] for a in self.atoms}
        for r in self.event_log:
            if r.event is not EventType.INVERSION:
                continue
            for op in r.genes:
                source, ss, se = self.registry.provenance[op.gid]
                for a in self._by_source.get(source, ()):
                    if a.start >= ss and a.end <= se:
                        out[a.atom_id].append((r.branch, r.time))
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


def _build_atoms(leaf_genomes: dict, root_length: int) -> list[Atom]:
    """Partition each source into atoms at the union of all extant-leaf breakpoints.

    With deletion, some ancestral positions survive in no extant leaf; those gaps carry
    no atom. So an interval between consecutive breakpoints becomes an atom only if some
    leaf still covers it.
    """
    bounds: dict[str, set[int]] = {}
    spans: dict[str, list[tuple[int, int]]] = {}
    for genome in leaf_genomes.values():
        for seg in genome._segments:
            bounds.setdefault(seg.source, {0}).add(seg.src_start)
            bounds[seg.source].add(seg.src_end)
            spans.setdefault(seg.source, []).append((seg.src_start, seg.src_end))
    atoms: list[Atom] = []
    aid = 0
    for source in sorted(bounds):
        covered = _merge(spans[source])
        cuts = sorted(bounds[source])
        for a, b in zip(cuts, cuts[1:]):
            if any(lo <= a and b <= hi for lo, hi in covered):  # surviving material only
                atoms.append(Atom(aid, source, a, b))
                aid += 1
    return atoms


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
    initial_size: int = 1,
    transfers=None,
    seed: int | None = None,
    rng: np.random.Generator | None = None,
    sampler: EventSampler | None = None,
) -> NucleotideResult:
    """Simulate variable-length structural events forward along ``species_tree``.

    ``inversion``, ``loss``, ``duplication``, ``transfer`` and ``transposition`` are
    **per-nucleotide** rates (total genome rate = ``rate * current_length``);
    ``origination`` is a **per-branch** rate that inserts a novel gene under a fresh source.
    ``extension`` sets the geometric event-length model (mean ``1/(1-extension)``
    nucleotides). ``transfers`` is an optional :class:`~zombi2.TransferModel` (default:
    additive, uniform recipient, no self-transfer). Returns a :class:`NucleotideResult`
    carrying the extant leaf genomes, the event log, the segment registry, and the atom
    partition (over the surviving ancestral material).

    Duplication and additive transfer grow the genome with no cap, so keep them at or below
    ``loss`` over long ages to avoid runaway growth.
    """
    if rng is None:
        rng = np.random.default_rng(seed)
    rates = UniformRates(inversion=inversion, loss=loss, duplication=duplication,
                         transfer=transfer, transposition=transposition,
                         origination=origination)
    registry = SegmentRegistry()

    def factory(ids):
        return NucleotideGenome(ids, root_length=root_length, extension=extension,
                                registry=registry)

    result = GenomeSimulator(sampler).simulate(
        species_tree, rates, rng, initial_size=initial_size, transfers=transfers,
        genome_factory=factory,
    )
    atoms = _build_atoms(result.leaf_genomes, root_length)
    return NucleotideResult(
        species_tree=species_tree,
        leaf_genomes=result.leaf_genomes,
        event_log=result.event_log,
        registry=registry,
        atoms=atoms,
        root_length=root_length,
    )
