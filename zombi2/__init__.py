"""ZOMBI2 — simulation of species trees (backward) and gene families (forward).

Public API::

    import zombi2 as z

    tree = z.simulate_species_tree(z.BirthDeath(birth=1.0, death=0.3),
                                   n_tips=20, age=5.0, seed=1)

    # every family the same rates:
    genomes = z.simulate_genomes(tree, duplication=0.2, transfer=0.1, loss=0.25,
                                 origination=0.5, seed=42)

    # every family its own rates, sampled from distributions:
    genomes = z.simulate_genomes(tree, z.FamilySampledRates(
        duplication=z.Gamma(2, 0.1), transfer=z.Exponential(0.1),
        loss=z.Gamma(2, 0.12), origination=0.5), seed=42)

    genomes.profiles.matrix   # families x extant-species copy numbers
    genomes.write("out/")
"""

from __future__ import annotations

__version__ = "0.2.0.dev0"

from .events import EventType, GeneOp, EventRecord, Selection, Region, TargetParams
from .tree import Tree, TreeNode, read_newick, prune
from .species_model import BirthDeath, Yule, EpisodicBirthDeath
from .species_sim import simulate_species_tree
from .ghosts import add_ghost_lineages
from .genome import Gene, Genome, UnorderedGenome, OrderedGene, OrderedGenome
from .nucleotide_genome import NucleotideGenome, Segment
from .nucleotide_sim import simulate_nucleotide_genomes, NucleotideResult, Atom
from .rates import (
    RateModel,
    UniformRates,
    GenomeWiseRates,
    FamilySampledRates,
    BranchRates,
    EventWeight,
)
from .transfers import TransferModel
from .distributions import (
    Distribution,
    Fixed,
    Exponential,
    Gamma,
    LogNormal,
    Uniform,
    as_distribution,
)
from .genome_sim import GenomeSimulator, GenomeResult
from .profiles import ProfileMatrix
from .reconciliation import build_gene_trees
from .rate_variation import RateVariation, RateScaledTree
from .traits import (
    BrownianMotion,
    OrnsteinUhlenbeck,
    MultivariateBrownian,
    MultivariateOU,
    MultiOptimumOU,
    ThresholdModel,
    EarlyBurst,
    Mk,
    CorrelatedBinary,
    HiddenStateMk,
    simulate_traits,
    replicate_traits,
    TraitResult,
    pagel_lambda,
    pagel_delta,
    pagel_kappa,
)
from .biogeography import DEC, simulate_biogeography
from .simulation import simulate_genomes, Genomes
from .matching import (
    match_profiles,
    match_profiles_smc,
    ABCFit,
    default_summary,
    default_gene_tree_summary,
    event_count_summary,
    frequency_spectrum,
    genome_sizes,
    copy_number_spectrum,
)
from .parallel import run_replicates
from ._rust import available as rust_available

__all__ = [
    "__version__",
    # events
    "EventType", "GeneOp", "EventRecord", "Selection", "Region", "TargetParams",
    # tree
    "Tree", "TreeNode", "read_newick", "prune",
    # species tree
    "BirthDeath", "Yule", "EpisodicBirthDeath",
    "simulate_species_tree", "add_ghost_lineages",
    # genome
    "Gene", "Genome", "UnorderedGenome", "OrderedGene", "OrderedGenome",
    # nucleotide genome (structural events at nucleotide resolution)
    "NucleotideGenome", "Segment", "simulate_nucleotide_genomes", "NucleotideResult", "Atom",
    # rates & transfers
    "RateModel", "UniformRates", "GenomeWiseRates", "FamilySampledRates",
    "BranchRates", "EventWeight", "TransferModel",
    # distributions
    "Distribution", "Fixed", "Exponential", "Gamma", "LogNormal", "Uniform", "as_distribution",
    # simulation
    "GenomeSimulator", "GenomeResult", "ProfileMatrix", "simulate_genomes", "Genomes",
    "build_gene_trees",
    # profile matching (rejection ABC + SMC)
    "match_profiles", "match_profiles_smc", "ABCFit", "default_summary", "default_gene_tree_summary",
    "event_count_summary", "frequency_spectrum", "genome_sizes", "copy_number_spectrum",
    # rate variation (relaxed clock)
    "RateVariation", "RateScaledTree",
    # trait evolution (phylogenetic comparative models)
    "BrownianMotion", "OrnsteinUhlenbeck", "MultivariateBrownian", "MultivariateOU",
    "MultiOptimumOU", "ThresholdModel", "EarlyBurst", "Mk", "CorrelatedBinary",
    "HiddenStateMk", "simulate_traits", "replicate_traits", "TraitResult",
    "pagel_lambda", "pagel_delta", "pagel_kappa",
    # historical biogeography (range evolution)
    "DEC", "simulate_biogeography",
    # parallelism
    "run_replicates",
    # Rust engine diagnostic (the built-in model requires the compiled extension)
    "rust_available",
]
