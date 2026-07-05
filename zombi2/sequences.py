"""Sequence-evolution namespace (scikit-learn-style).

Re-exports the sequence / substitution public API so users can write::

    from zombi2.sequences import lg, gtr, evolve_on_tree

Every name here is the *same object* as the corresponding top-level
``zombi2`` attribute -- this module is a thin, additive namespace over the
implementation modules (:mod:`zombi2.sequence_evolution`,
:mod:`zombi2.sequence_sim`); it does not redefine anything.

The relaxed molecular clocks (``RateVariation`` and the rest of the
:class:`~zombi2.rate_variation.Clock` family) that feed the shared lineage clock now live in
their own :mod:`zombi2.clocks` namespace.
"""

from __future__ import annotations

from .sequence_evolution import SequenceEvolution, GenePhylograms
from .sequence_sim import (
    SubstitutionModel, GammaRates, jc69, k80, hky85, gtr,
    poisson, lg, wag, jtt, dayhoff, make_model, is_protein_model,
    DNA_MODELS, PROTEIN_MODELS, AMINO_ACIDS,
    evolve_on_tree, read_fasta, write_fasta,
)

__all__ = [
    "SequenceEvolution", "GenePhylograms",
    "SubstitutionModel", "GammaRates", "jc69", "k80", "hky85", "gtr",
    "poisson", "lg", "wag", "jtt", "dayhoff", "make_model", "is_protein_model",
    "DNA_MODELS", "PROTEIN_MODELS", "AMINO_ACIDS",
    "evolve_on_tree", "read_fasta", "write_fasta",
]
