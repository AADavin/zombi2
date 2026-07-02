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

from dataclasses import dataclass, field

import numpy as np

from ._sampling import EventSampler
from .events import EventType
from .genome_sim import GenomeSimulator
from .nucleotide_genome import NucleotideGenome
from .rates import UniformRates
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
    registry: dict      # seg_id -> (source, src_start, src_end)
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
                source, ss, se = self.registry[op.gid]
                for a in self._by_source.get(source, ()):
                    if a.start >= ss and a.end <= se:
                        out[a.atom_id].append((r.branch, r.time))
        return out


def _build_atoms(leaf_genomes: dict, root_length: int) -> list[Atom]:
    """Partition each source into atoms at the union of all extant-leaf breakpoints."""
    bounds: dict[str, set[int]] = {}
    for genome in leaf_genomes.values():
        for seg in genome._segments:
            b = bounds.setdefault(seg.source, {0})
            b.add(seg.src_start)
            b.add(seg.src_end)
    atoms: list[Atom] = []
    aid = 0
    for source in sorted(bounds):
        cuts = sorted(bounds[source])
        for a, b in zip(cuts, cuts[1:]):
            atoms.append(Atom(aid, source, a, b))
            aid += 1
    return atoms


def simulate_nucleotide_genomes(
    species_tree: Tree,
    *,
    inversion: float = 0.001,
    root_length: int = 1000,
    extension: float | None = 0.99,
    initial_size: int = 1,
    seed: int | None = None,
    rng: np.random.Generator | None = None,
    sampler: EventSampler | None = None,
) -> NucleotideResult:
    """Simulate variable-length inversions forward along ``species_tree`` (M1).

    ``inversion`` is a per-nucleotide rate: the total genome inversion rate is
    ``inversion * root_length``. ``extension`` sets the geometric inversion-length model
    (mean ``1/(1-extension)`` nucleotides). Returns a :class:`NucleotideResult` carrying
    the extant leaf genomes, the event log, the segment registry, and the atom partition.
    """
    if rng is None:
        rng = np.random.default_rng(seed)
    rates = UniformRates(inversion=inversion)
    registry: dict[str, tuple[str, int, int]] = {}

    def factory(ids):
        return NucleotideGenome(ids, root_length=root_length, extension=extension,
                                registry=registry)

    result = GenomeSimulator(sampler).simulate(
        species_tree, rates, rng, initial_size=initial_size, genome_factory=factory,
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
