"""Genomes — the gene-family core, at three resolutions (SPEC §7).

A genome evolves along the species tree by **duplication, transfer, loss and origination**. The three
resolutions share one spine (a forward Gillespie over the complete tree, the same rate grammar, the
same gene-genealogy event log) and differ only in the state they carry:

- `family` — a genome is a **multiset** of gene families (`simulate_genomes_family()`);
- `ordered` — genes gain a **position and orientation** on chromosomes
  (`simulate_genomes_ordered()`);
- `nucleotide` — genes and intergenes carry **sequence** (`simulate_genomes_nucleotide()`).

This package's ``__init__`` is only the public surface; each engine lives in its own module, with the
recipient rules in `_transfer`, the per-family parallel engine in `_perfamily`, and the
derived views (gene trees, profiles) in their own modules.
"""

from __future__ import annotations

from ..params.mapping import Between
from ._transfer import Clades, Distance
from .events import GeneEdge
from .gene_trees import GeneNode, GeneTree
from .profiles import Profiles
from .family import (GeneCopy, FamilyGenomesResult, FamilyGenome, GeneFamily,
                     simulate_genomes_family, family, genome)
# re-exported on the package path for the CLI / tests, but kept out of __all__ (not public API):
from .family import IMPLEMENTED_MODIFIERS, resolve_max_family_size  # noqa: F401
from .chromosomes import ChromosomeEvent
from .ordered import (
    Chromosome,
    EventPosition,
    Gene,
    Inversion,
    OrderedGenomesResult,
    Translocation,
    Transposition,
    simulate_genomes_ordered,
)
from .nucleotide import NucleotideGenome, NucleotideGenomesResult, simulate_genomes_nucleotide
from ._perfamily import StreamedRun
from .read import read_run

__all__ = ["simulate_genomes_family", "read_run", "FamilyGenomesResult", "GeneEdge", "GeneCopy", "Distance",
           "Clades", "Between",
           "Profiles", "GeneTree", "GeneNode", "FamilyGenome", "genome", "GeneFamily", "family",
           "simulate_genomes_ordered", "OrderedGenomesResult", "Gene", "Chromosome",
           "ChromosomeEvent", "Inversion", "Transposition", "Translocation", "EventPosition",
           "simulate_genomes_nucleotide", "NucleotideGenomesResult", "NucleotideGenome",
           "StreamedRun"]
