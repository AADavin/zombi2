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

import warnings

from zombi2.genomes.genome import Gene, Genome, UnorderedGenome, OrderedGene, OrderedGenome
from zombi2.genomes.nucleotide_genome import NucleotideGenome, Segment
from zombi2.genomes.nucleotide_sim import simulate_nucleotide_genomes, NucleotideResult, Block
from zombi2.genomes.gff import read_gff, read_gff_all, GffGenome
from zombi2.genomes.rates import (
    RateModel, Rates, Per,
    PerCopyRates, PerLineageRates,  # deprecated presets for Rates(per=…); importable, not in __all__
    FamilySampledRates,
    LineageRates, Modifier, ModifiedRates, LineageModifier,
    FamilyModifier, EventWeight,
)
from zombi2.genomes.transfers import TransferModel, PairModifier
from zombi2.genomes.conversion import ConversionModel
from zombi2.genomes.read_rates import read_family_rates, read_lineage_rates
from zombi2.genomes.genome_sim import GenomeSimulator, GenomeResult
from zombi2.genomes.profiles import ProfileMatrix
from zombi2.genomes.reconciliation import build_gene_trees
from zombi2.genomes.simulation import simulate_genomes, Genomes, GenomeTrace, read_events_trace
from zombi2.parallel import run_replicates

__all__ = [
    "Gene", "Genome", "UnorderedGenome", "OrderedGene", "OrderedGenome",
    "NucleotideGenome", "Segment", "simulate_nucleotide_genomes",
    "NucleotideResult", "Block", "read_gff", "read_gff_all", "GffGenome",
    "RateModel", "Rates", "Per",
    "FamilySampledRates", "LineageRates", "Modifier", "ModifiedRates",
    "LineageModifier", "FamilyModifier",
    "EventWeight", "TransferModel", "PairModifier", "ConversionModel",
    "read_family_rates", "read_lineage_rates",
    "GenomeSimulator", "GenomeResult", "ProfileMatrix",
    "simulate_genomes", "Genomes", "GenomeTrace", "read_events_trace",
    "build_gene_trees", "run_replicates",
]

# --- deprecated aliases (PEP 562) -------------------------------------------
# Kept working for one minor version but marked: they resolve here with a
# DeprecationWarning and are absent from __all__/dir(). Removal: 0.4.0.
# See docs/design/naming-consolidation.md.
_DEPRECATED_ALIASES = {
    "SharedRates": "PerCopyRates",
    "PerGenomeRates": "PerLineageRates",
    "BranchRates": "LineageRates",
    "BranchModifier": "LineageModifier",
    "read_branch_rates": "read_lineage_rates",
}


def __getattr__(name):
    canonical = _DEPRECATED_ALIASES.get(name)
    if canonical is not None:
        warnings.warn(
            f"zombi2.genomes.{name} was renamed to zombi2.genomes.{canonical}; the old name "
            f"still works but is deprecated and will be removed in 0.4.0.",
            DeprecationWarning, stacklevel=2,
        )
        return globals()[canonical]
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
