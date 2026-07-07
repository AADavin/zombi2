"""Coevolution / coupling namespace (scikit-learn-style).

Re-exports the coupling and state-dependent-coevolution public API so users
can write::

    from zombi2.coevolve import simulate_trait_linked_genomes, BiSSE

Every name here is the *same object* as the corresponding top-level
``zombi2`` attribute -- this module is a thin, additive namespace over the
implementation modules (:mod:`zombi2.coevolve.trait_coupling`,
:mod:`zombi2.coevolve.gene_diversification`,
:mod:`zombi2.coevolve.cladogenetic_genome`,
:mod:`zombi2.coevolve.gene_conditioned_trait`,
:mod:`zombi2.coevolve.trait_gene_feedback`, :mod:`zombi2.coevolve.sse`);
it does not redefine anything.
"""

from __future__ import annotations

from zombi2.coevolve.trait_coupling import (
    TraitGeneCoupling, TraitTrajectory, TraitLinkedRates, TraitLinkedResult,
    simulate_trait_linked_genomes,
)
from zombi2.coevolve.gene_diversification import (
    GeneDiversification, GeneDiversificationResult, simulate_gene_diversification,
    simulate_co_diversification,
)
from zombi2.coevolve.cladogenetic_genome import (
    CladogeneticGenome, CladogeneticGenomeResult, simulate_cladogenetic_genome,
)
from zombi2.coevolve.gene_conditioned_trait import (
    GeneConditionedTrait, GeneConditionedTraitResult,
    simulate_gene_conditioned_trait,
)
from zombi2.coevolve.trait_gene_feedback import (
    TraitGeneFeedback, TraitGeneFeedbackResult, simulate_trait_gene_feedback,
)
from zombi2.coevolve.sse import BiSSE, MuSSE, HiSSE, QuaSSE, simulate_sse

__all__ = [
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
