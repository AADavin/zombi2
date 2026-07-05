"""Gene-family / genome namespace (scikit-learn-style).

Re-exports the genome-simulation public API so users can write::

    from zombi2.genomes import simulate_genomes, FamilySampledRates

Every name here is the *same object* as the corresponding top-level
``zombi2`` attribute -- this module is a thin, additive namespace over the
implementation modules (:mod:`zombi2.genome`, :mod:`zombi2.nucleotide_genome`,
:mod:`zombi2.nucleotide_sim`, :mod:`zombi2.gff`, :mod:`zombi2.rates`,
:mod:`zombi2.transfers`, :mod:`zombi2.genome_sim`, :mod:`zombi2.profiles`,
:mod:`zombi2.reconciliation`, :mod:`zombi2.simulation`,
:mod:`zombi2.parallel`); it does not redefine anything.
"""

from __future__ import annotations

from .genome import Gene, Genome, UnorderedGenome, OrderedGene, OrderedGenome
from .nucleotide_genome import NucleotideGenome, Segment
from .nucleotide_sim import simulate_nucleotide_genomes, NucleotideResult, Atom
from .gff import read_gff, GffGenome
from .rates import (
    RateModel, UniformRates, GenomeWiseRates, FamilySampledRates, BranchRates,
    EventWeight,
)
from .transfers import TransferModel
from .genome_sim import GenomeSimulator, GenomeResult
from .profiles import ProfileMatrix
from .reconciliation import build_gene_trees
from .simulation import simulate_genomes, Genomes, GenomeTrace, read_events_trace
from .parallel import run_replicates

__all__ = [
    "Gene", "Genome", "UnorderedGenome", "OrderedGene", "OrderedGenome",
    "NucleotideGenome", "Segment", "simulate_nucleotide_genomes",
    "NucleotideResult", "Atom", "read_gff", "GffGenome",
    "RateModel", "UniformRates", "GenomeWiseRates", "FamilySampledRates",
    "BranchRates", "EventWeight", "TransferModel",
    "GenomeSimulator", "GenomeResult", "ProfileMatrix",
    "simulate_genomes", "Genomes", "GenomeTrace", "read_events_trace",
    "build_gene_trees", "run_replicates",
]
