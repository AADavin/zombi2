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

# Low-level primitives (kept at top level only; no scikit-style namespace).
from .events import EventType, GeneOp, EventRecord, Selection, Region, TargetParams
from ._rust import available as rust_available

# scikit-learn-style namespaces are the single source of truth for the
# public API: importing from them here guarantees that, e.g.,
# ``zombi2.BirthDeath is zombi2.species.BirthDeath`` (same object). The
# namespace modules are thin re-exports over the implementation modules and do
# not redefine anything; see zombi2/species.py, zombi2/genomes.py, etc.
from .species import (
    Tree, TreeNode, read_newick, prune,
    BirthDeath, Yule, EpisodicBirthDeath, ClaDS, DiversityDependent,
    CladeShiftBirthDeath, simulate_species_tree, add_ghost_lineages,
)
from .genomes import (
    Gene, Genome, UnorderedGenome, OrderedGene, OrderedGenome,
    NucleotideGenome, Segment, simulate_nucleotide_genomes, NucleotideResult, Block,
    read_gff, GffGenome,
    RateModel, SharedRates, PerGenomeRates, FamilySampledRates, BranchRates,
    EventWeight, TransferModel,
    GenomeSimulator, GenomeResult, ProfileMatrix,
    simulate_genomes, Genomes, GenomeTrace, read_events_trace,
    build_gene_trees, run_replicates,
)
from .distributions import (
    Distribution, Fixed, Exponential, Gamma, LogNormal, Uniform, as_distribution,
)
from .clocks import (
    Clock, RateScaledTree, StrictClock, UncorrelatedLogNormalClock,
    UncorrelatedGammaClock, WhiteNoiseClock, AutocorrelatedLogNormalClock,
    CIRClock, RateVariation,
)
from .sequences import (
    SequenceEvolution, GenePhylograms,
    SubstitutionModel, GammaRates, jc69, k80, hky85, gtr,
    poisson, lg, wag, jtt, dayhoff, make_model, is_protein_model,
    DNA_MODELS, PROTEIN_MODELS, AMINO_ACIDS,
    evolve_on_tree, read_fasta, write_fasta,
)
from .traits import (
    BrownianMotion, OrnsteinUhlenbeck, MultivariateBrownian, MultivariateOU,
    MultiOptimumOU, ThresholdModel, EarlyBurst, Mk, CorrelatedBinary,
    HiddenStateMk, Cladogenesis, simulate_traits, replicate_traits, TraitResult,
    pagel_lambda, pagel_delta, pagel_kappa,
    DEC, simulate_biogeography,
)
from .coevolve import (
    CouplingSpec, PottsRates, pathway_blocks, simulate_coupled, CoupledResult,
    TraitGeneCoupling, TraitTrajectory, TraitLinkedRates, TraitLinkedResult,
    simulate_trait_linked_genomes,
    GeneDiversification, GeneDiversificationResult, simulate_gene_diversification,
    CladogeneticGenome, CladogeneticGenomeResult, simulate_cladogenetic_genome,
    GeneConditionedTrait, GeneConditionedTraitResult, simulate_gene_conditioned_trait,
    BiSSE, MuSSE, HiSSE, QuaSSE, simulate_sse,
)
# NOTE: ABC profile-matching inference (zombi2.matching / zombi2.abc) is withheld from the v1
# public surface — it is not yet documented/stabilised. The modules remain in-tree (import via
# ``from zombi2.matching import match_profiles``) and fully tested; re-add the export here, the
# 'abc' CLI command, and the docs nav entry to promote it back to the public API.

__all__ = [
    "__version__",
    # events
    "EventType", "GeneOp", "EventRecord", "Selection", "Region", "TargetParams",
    # tree
    "Tree", "TreeNode", "read_newick", "prune",
    # species tree
    "BirthDeath", "Yule", "EpisodicBirthDeath", "ClaDS", "DiversityDependent",
    "CladeShiftBirthDeath", "simulate_species_tree", "add_ghost_lineages",
    # genome
    "Gene", "Genome", "UnorderedGenome", "OrderedGene", "OrderedGenome",
    # nucleotide genome (structural events at nucleotide resolution)
    "NucleotideGenome", "Segment", "simulate_nucleotide_genomes", "NucleotideResult", "Block",
    "read_gff", "GffGenome",
    # rates & transfers
    "RateModel", "SharedRates", "PerGenomeRates", "FamilySampledRates",
    "BranchRates", "EventWeight", "TransferModel",
    # distributions
    "Distribution", "Fixed", "Exponential", "Gamma", "LogNormal", "Uniform", "as_distribution",
    # simulation
    "GenomeSimulator", "GenomeResult", "ProfileMatrix", "simulate_genomes", "Genomes", "GenomeTrace",
    "read_events_trace", "build_gene_trees",
    # gene-family coupling (Potts/Ising non-independence)
    "CouplingSpec", "PottsRates", "pathway_blocks", "simulate_coupled", "CoupledResult",
    # trait-conditioned gene families (trait <-> gene-family coupling)
    "TraitGeneCoupling", "TraitTrajectory", "TraitLinkedRates", "TraitLinkedResult",
    "simulate_trait_linked_genomes",
    # profile matching / ABC inference — withheld from the v1 public surface (see note above)
    # relaxed molecular clocks (chronogram -> phylogram; the shared lineage clock family)
    "Clock", "RateScaledTree", "StrictClock", "UncorrelatedLogNormalClock",
    "UncorrelatedGammaClock", "WhiteNoiseClock", "AutocorrelatedLogNormalClock",
    "CIRClock", "RateVariation",
    # family sequence evolution (gene x lineage substitution clock)
    "SequenceEvolution", "GenePhylograms",
    # sequence simulation (evolve DNA or protein along a gene tree; ancestral genome reconstruction)
    "SubstitutionModel", "GammaRates", "jc69", "k80", "hky85", "gtr",
    "poisson", "lg", "wag", "jtt", "dayhoff", "make_model", "is_protein_model",
    "DNA_MODELS", "PROTEIN_MODELS", "AMINO_ACIDS",
    "evolve_on_tree", "read_fasta", "write_fasta",
    # trait evolution (phylogenetic comparative models)
    "BrownianMotion", "OrnsteinUhlenbeck", "MultivariateBrownian", "MultivariateOU",
    "MultiOptimumOU", "ThresholdModel", "EarlyBurst", "Mk", "CorrelatedBinary",
    "HiddenStateMk", "Cladogenesis", "simulate_traits", "replicate_traits", "TraitResult",
    "pagel_lambda", "pagel_delta", "pagel_kappa",
    # historical biogeography (range evolution)
    "DEC", "simulate_biogeography",
    # state-dependent diversification (trait drives the tree — coevolve traits:species)
    "BiSSE", "MuSSE", "HiSSE", "QuaSSE", "simulate_sse",
    # gene-content-dependent diversification (genes drive the tree — coevolve genes:species)
    "GeneDiversification", "GeneDiversificationResult", "simulate_gene_diversification",
    # cladogenetic genome evolution (speciation drives gene content — coevolve species:genes)
    "CladogeneticGenome", "CladogeneticGenomeResult", "simulate_cladogenetic_genome",
    # gene-conditioned trait (gene content drives a trait — coevolve genes:traits)
    "GeneConditionedTrait", "GeneConditionedTraitResult", "simulate_gene_conditioned_trait",
    # parallelism
    "run_replicates",
    # Rust engine diagnostic (the built-in model requires the compiled extension)
    "rust_available",
]
