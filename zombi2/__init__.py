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
from .tree import Tree, TreeNode
from .species_model import BirthDeath, Yule, EpisodicBirthDeath
from .species_sim import simulate_species_tree
from .genome import Gene, Genome, UnorderedGenome, OrderedGene, OrderedGenome
from .rates import (
    RateModel,
    UniformRates,
    GenomeWiseRates,
    FamilySampledRates,
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
from .simulation import simulate_genomes, Genomes

__all__ = [
    "__version__",
    # events
    "EventType", "GeneOp", "EventRecord", "Selection", "Region", "TargetParams",
    # tree
    "Tree", "TreeNode",
    # species tree
    "BirthDeath", "Yule", "EpisodicBirthDeath", "simulate_species_tree",
    # genome
    "Gene", "Genome", "UnorderedGenome", "OrderedGene", "OrderedGenome",
    # rates & transfers
    "RateModel", "UniformRates", "GenomeWiseRates", "FamilySampledRates",
    "EventWeight", "TransferModel",
    # distributions
    "Distribution", "Fixed", "Exponential", "Gamma", "LogNormal", "Uniform", "as_distribution",
    # simulation
    "GenomeSimulator", "GenomeResult", "ProfileMatrix", "simulate_genomes", "Genomes",
    "build_gene_trees",
]
