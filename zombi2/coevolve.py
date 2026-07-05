"""Coevolution / coupling namespace (scikit-learn-style).

Re-exports the coupling and state-dependent-coevolution public API so users
can write::

    from zombi2.coevolve import simulate_coupled, BiSSE

Every name here is the *same object* as the corresponding top-level
``zombi2`` attribute -- this module is a thin, additive namespace over the
implementation modules (:mod:`zombi2.coupling`, :mod:`zombi2.trait_coupling`,
:mod:`zombi2.gene_diversification`, :mod:`zombi2.cladogenetic_genome`,
:mod:`zombi2.gene_conditioned_trait`, :mod:`zombi2.sse`); it does not redefine
anything.
"""

from __future__ import annotations

from .coupling import (
    CouplingSpec, PottsRates, pathway_blocks, simulate_coupled, CoupledResult,
)
from .trait_coupling import (
    TraitGeneCoupling, TraitTrajectory, TraitLinkedRates, TraitLinkedResult,
    simulate_trait_linked_genomes,
)
from .gene_diversification import (
    GeneDiversification, GeneDiversificationResult, simulate_gene_diversification,
    simulate_co_diversification,
)
from .cladogenetic_genome import (
    CladogeneticGenome, CladogeneticGenomeResult, simulate_cladogenetic_genome,
)
from .gene_conditioned_trait import (
    GeneConditionedTrait, GeneConditionedTraitResult,
    simulate_gene_conditioned_trait,
)
from .trait_gene_feedback import (
    TraitGeneFeedback, TraitGeneFeedbackResult, simulate_trait_gene_feedback,
)
from .sse import BiSSE, MuSSE, HiSSE, QuaSSE, simulate_sse

__all__ = [
    "CouplingSpec", "PottsRates", "pathway_blocks", "simulate_coupled",
    "CoupledResult",
    "TraitGeneCoupling", "TraitTrajectory", "TraitLinkedRates",
    "TraitLinkedResult", "simulate_trait_linked_genomes",
    "GeneDiversification", "GeneDiversificationResult",
    "simulate_gene_diversification", "simulate_co_diversification",
    "CladogeneticGenome", "CladogeneticGenomeResult",
    "simulate_cladogenetic_genome",
    "GeneConditionedTrait", "GeneConditionedTraitResult",
    "simulate_gene_conditioned_trait",
    "TraitGeneFeedback", "TraitGeneFeedbackResult", "simulate_trait_gene_feedback",
    "BiSSE", "MuSSE", "HiSSE", "QuaSSE", "simulate_sse",
]
