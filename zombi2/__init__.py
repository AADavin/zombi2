"""ZOMBI2 — simulation of species trees (backward) and gene families (forward).

Public API (see the individual modules for details)::

    import zombi2 as z
    species = z.SpeciesTreeModel(birth=1.0, death=0.3, n_tips=20, age=5.0)
    rates = z.RateModel(z.EventRates(duplication=0.2, transfer=0.1, loss=0.25,
                                     origination=0.5))
    result = z.Simulation(species, rates, seed=42).run()
"""

from __future__ import annotations

__version__ = "0.1.0.dev0"

from .events import EventType, GeneOp, EventRecord, Selection, Region, TargetParams
from .tree import Tree, TreeNode
from .species_model import SpeciesTreeModel
from .species_sim import SpeciesTreeSimulator, simulate_species_tree
from .genome import Gene, Genome, UnorderedGenome
from .rates import EventRates, RateModel, ConstantGeneWiseRates
from .genome_sim import GenomeSimulator
from .profiles import ProfileMatrix
from .simulation import Simulation, SimulationResult

__all__ = [
    "__version__",
    "EventType",
    "GeneOp",
    "EventRecord",
    "Selection",
    "Region",
    "TargetParams",
    "Tree",
    "TreeNode",
    "SpeciesTreeModel",
    "SpeciesTreeSimulator",
    "simulate_species_tree",
    "Gene",
    "Genome",
    "UnorderedGenome",
    "EventRates",
    "RateModel",
    "ConstantGeneWiseRates",
    "GenomeSimulator",
    "ProfileMatrix",
    "Simulation",
    "SimulationResult",
]
