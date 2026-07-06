"""Sequence-evolution namespace (scikit-learn-style).

Re-exports the sequence / substitution public API so users can write::

    from zombi2.sequences import lg, gtr, evolve_on_tree

Every name here is the *same object* as the corresponding top-level
``zombi2`` attribute -- this module is a thin, additive namespace over the
implementation modules (:mod:`zombi2.sequences.evolution`,
:mod:`zombi2.sequences.models`, :mod:`zombi2.sequences.clocks`); it does not
redefine anything.

Sequence evolution has two parts, both here: how *fast* substitutions accumulate
on each branch -- the relaxed molecular **clocks** (:class:`Clock` family, which
turn a chronogram into a phylogram) -- and *what* changes -- the substitution
**models** (JC/HKY/GTR + empirical amino-acid models).
"""

from __future__ import annotations

from zombi2.sequences.clocks import (
    Clock, RateScaledTree, StrictClock, UncorrelatedLogNormalClock,
    UncorrelatedGammaClock, WhiteNoiseClock, AutocorrelatedLogNormalClock,
    CIRClock, RateVariation,
)
from zombi2.sequences.evolution import SequenceEvolution, GenePhylograms
from zombi2.sequences.models import (
    SubstitutionModel, GammaRates, jc69, k80, hky85, gtr,
    poisson, lg, wag, jtt, dayhoff, make_model, is_protein_model,
    DNA_MODELS, PROTEIN_MODELS, AMINO_ACIDS,
    evolve_on_tree, read_fasta, write_fasta,
)

__all__ = [
    "Clock", "RateScaledTree", "StrictClock", "UncorrelatedLogNormalClock",
    "UncorrelatedGammaClock", "WhiteNoiseClock", "AutocorrelatedLogNormalClock",
    "CIRClock", "RateVariation",
    "SequenceEvolution", "GenePhylograms",
    "SubstitutionModel", "GammaRates", "jc69", "k80", "hky85", "gtr",
    "poisson", "lg", "wag", "jtt", "dayhoff", "make_model", "is_protein_model",
    "DNA_MODELS", "PROTEIN_MODELS", "AMINO_ACIDS",
    "evolve_on_tree", "read_fasta", "write_fasta",
]
