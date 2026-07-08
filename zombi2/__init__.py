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

__version__ = "0.2.0"

# Low-level primitives (kept at top level only; no scikit-style namespace).
from zombi2.genomes.events import EventType, GeneOp, EventRecord, Selection, Region, TargetParams
from zombi2._rust import available as rust_available

# scikit-learn-style namespaces are the single source of truth for the
# public API: importing from them here guarantees that, e.g.,
# ``zombi2.BirthDeath is zombi2.species.BirthDeath`` (same object). The
# namespace modules are thin re-exports over the implementation modules and do
# not redefine anything; see zombi2/species.py, zombi2/genomes.py, etc.
from zombi2.species import (
    Tree, TreeNode, read_newick, prune,
    BirthDeath, Yule, EpisodicBirthDeath, ClaDS, DiversityDependent,
    CladeShiftBirthDeath, simulate_species_tree, add_ghost_lineages,
)
from zombi2.genomes import (
    Gene, Genome, UnorderedGenome, OrderedGene, OrderedGenome,
    NucleotideGenome, Segment, simulate_nucleotide_genomes, NucleotideResult, Block,
    read_gff, GffGenome,
    RateModel, SharedRates, PerGenomeRates, FamilySampledRates, BranchRates,
    EventWeight, TransferModel, ConversionModel, read_family_rates, read_branch_rates,
    GenomeSimulator, GenomeResult, ProfileMatrix,
    simulate_genomes, Genomes, GenomeTrace, read_events_trace,
    build_gene_trees, run_replicates,
)
from zombi2.distributions import (
    Distribution, Fixed, Exponential, Gamma, LogNormal, Uniform, as_distribution,
)
from zombi2.sequences import (
    # relaxed molecular clocks (chronogram -> phylogram; the shared lineage clock family)
    Clock, RateScaledTree, StrictClock, UncorrelatedLogNormalClock,
    UncorrelatedGammaClock, WhiteNoiseClock, AutocorrelatedLogNormalClock,
    CIRClock, RateVariation,
    # substitution models + the gene x lineage substitution clock
    SequenceEvolution, GenePhylograms,
    SubstitutionModel, GammaRates, jc69, k80, hky85, gtr,
    poisson, lg, wag, jtt, dayhoff, make_model, is_protein_model,
    DNA_MODELS, PROTEIN_MODELS, AMINO_ACIDS,
    evolve_on_tree, read_fasta, write_fasta,
)
from zombi2.traits import (
    BrownianMotion, OrnsteinUhlenbeck, MultivariateBrownian, MultivariateOU,
    MultiOptimumOU, ThresholdModel, EarlyBurst, Mk, CorrelatedBinary,
    CorrelatedBinaryK, HiddenStateMk, Cladogenesis, simulate_traits,
    replicate_traits, TraitResult,
    pagel_lambda, pagel_delta, pagel_kappa,
    DEC, simulate_biogeography,
)
from zombi2.coevolve import (
    TraitGeneCoupling, TraitTrajectory, TraitLinkedRates, TraitLinkedResult,
    simulate_trait_linked_genomes,
    GeneDiversification, GeneDiversificationResult, simulate_gene_diversification,
    simulate_co_diversification,
    CladogeneticGenome, CladogeneticGenomeResult, simulate_cladogenetic_genome,
    GeneConditionedTrait, GeneConditionedTraitResult, simulate_gene_conditioned_trait,
    TraitGeneFeedback, TraitGeneFeedbackResult, simulate_trait_gene_feedback,
    BiSSE, MuSSE, HiSSE, QuaSSE, CID, simulate_sse,
)
# NOTE: ABC profile-matching inference has moved out of the core to
# ``ZOMBI2_FUTURE/abc-inference/`` — inference is a Phase-3 Extension, not a core simulation
# level. See that folder's README to revive it as an Extension.

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
    "BranchRates", "EventWeight", "TransferModel", "ConversionModel",
    "read_family_rates", "read_branch_rates",
    # distributions
    "Distribution", "Fixed", "Exponential", "Gamma", "LogNormal", "Uniform", "as_distribution",
    # simulation
    "GenomeSimulator", "GenomeResult", "ProfileMatrix", "simulate_genomes", "Genomes", "GenomeTrace",
    "read_events_trace", "build_gene_trees",
    # trait-conditioned gene families (trait <-> gene-family coupling)
    "TraitGeneCoupling", "TraitTrajectory", "TraitLinkedRates", "TraitLinkedResult",
    "simulate_trait_linked_genomes",
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
    "CorrelatedBinaryK", "HiddenStateMk", "Cladogenesis", "simulate_traits",
    "replicate_traits", "TraitResult",
    "pagel_lambda", "pagel_delta", "pagel_kappa",
    # historical biogeography (range evolution)
    "DEC", "simulate_biogeography",
    # state-dependent diversification (trait drives the tree — coevolve traits:species)
    "BiSSE", "MuSSE", "HiSSE", "QuaSSE", "CID", "simulate_sse",
    # gene-content-dependent diversification (genes drive the tree — coevolve genes:species,
    # + species:genes = the co-diversification joint model)
    "GeneDiversification", "GeneDiversificationResult", "simulate_gene_diversification",
    "simulate_co_diversification",
    # cladogenetic genome evolution (speciation drives gene content — coevolve species:genes)
    "CladogeneticGenome", "CladogeneticGenomeResult", "simulate_cladogenetic_genome",
    # gene-conditioned trait (gene content drives a trait — coevolve genes:traits)
    "GeneConditionedTrait", "GeneConditionedTraitResult", "simulate_gene_conditioned_trait",
    # trait<->genes joint model (trait-gene feedback — coevolve traits:genes + genes:traits)
    "TraitGeneFeedback", "TraitGeneFeedbackResult", "simulate_trait_gene_feedback",
    # parallelism
    "run_replicates",
    # Rust engine diagnostic (the built-in model requires the compiled extension)
    "rust_available",
]
