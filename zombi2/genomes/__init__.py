"""Gene-family / genome namespace (scikit-learn-style).

Re-exports the genome-simulation public API so users can write::

    from zombi2.genomes import simulate_genomes, FamilySampledRates

Every name here is the *same object* as the corresponding top-level
``zombi2`` attribute -- this module is a thin, additive namespace over the
implementation modules (:mod:`zombi2.genomes.genome`,
:mod:`zombi2.genomes.nucleotide_genome`, :mod:`zombi2.genomes.nucleotide_sim`,
:mod:`zombi2.genomes.gff`, :mod:`zombi2.genomes.rates`,
:mod:`zombi2.genomes.transfers`, :mod:`zombi2.genomes.genome_sim`,
:mod:`zombi2.genomes.profiles`, :mod:`zombi2.genomes.reconciliation`,
:mod:`zombi2.genomes.simulation`, :mod:`zombi2.parallel`); it does not
redefine anything.
"""

from __future__ import annotations

from zombi2.genomes.genome import Gene, Genome, UnorderedGenome, OrderedGene, OrderedGenome
from zombi2.genomes.nucleotide_genome import NucleotideGenome, Segment
from zombi2.genomes.nucleotide_sim import simulate_nucleotide_genomes, NucleotideResult, Block
from zombi2.genomes.gff import read_gff, read_gff_all, GffGenome
from zombi2.genomes.rates import (
    RateModel, PerCopyRates, SharedRates, PerGenomeRates, FamilySampledRates, BranchRates,
    Modifier, ModifiedRates, BranchModifier, FamilyModifier, EventWeight,
)
from zombi2.genomes.transfers import TransferModel, PairModifier
from zombi2.genomes.conversion import ConversionModel
from zombi2.genomes.read_rates import read_family_rates, read_branch_rates
from zombi2.genomes.genome_sim import GenomeSimulator, GenomeResult
from zombi2.genomes.profiles import ProfileMatrix
from zombi2.genomes.reconciliation import build_gene_trees
from zombi2.genomes.simulation import simulate_genomes, Genomes, GenomeTrace, read_events_trace
from zombi2.parallel import run_replicates

__all__ = [
    "Gene", "Genome", "UnorderedGenome", "OrderedGene", "OrderedGenome",
    "NucleotideGenome", "Segment", "simulate_nucleotide_genomes",
    "NucleotideResult", "Block", "read_gff", "read_gff_all", "GffGenome",
    "RateModel", "PerCopyRates", "SharedRates", "PerGenomeRates", "FamilySampledRates",
    "BranchRates", "Modifier", "ModifiedRates", "BranchModifier", "FamilyModifier",
    "EventWeight", "TransferModel", "PairModifier", "ConversionModel",
    "read_family_rates", "read_branch_rates",
    "GenomeSimulator", "GenomeResult", "ProfileMatrix",
    "simulate_genomes", "Genomes", "GenomeTrace", "read_events_trace",
    "build_gene_trees", "run_replicates",
]
