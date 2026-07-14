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
    poisson, lg, wag, jtt, dayhoff, make_model, is_protein_model, is_codon_model,
    DNA_MODELS, PROTEIN_MODELS, CODON_MODELS, AMINO_ACIDS,
    evolve_on_tree, read_fasta, write_fasta,
)
from zombi2.sequences.codon_models import gy94, mg94, make_codon_model

# NOTE: the codon *utilities* (translate, GENETIC_CODE, SENSE_CODONS, STOP_CODONS, expected_dnds)
# are intentionally not re-exported here — import them from ``zombi2.sequences.codon_models``. This
# keeps the ``zombi2.sequences`` namespace to the model API and avoids clashing with the separate
# ``zombi2.experimental.codon_selection.translate``.

__all__ = [
    "Clock", "RateScaledTree", "StrictClock", "UncorrelatedLogNormalClock",
    "UncorrelatedGammaClock", "WhiteNoiseClock", "AutocorrelatedLogNormalClock",
    "CIRClock", "RateVariation",
    "SequenceEvolution", "GenePhylograms",
    "SubstitutionModel", "GammaRates", "jc69", "k80", "hky85", "gtr",
    "poisson", "lg", "wag", "jtt", "dayhoff", "make_model", "is_protein_model",
    "is_codon_model", "DNA_MODELS", "PROTEIN_MODELS", "CODON_MODELS", "AMINO_ACIDS",
    "evolve_on_tree", "read_fasta", "write_fasta",
    "gy94", "mg94", "make_codon_model",
]
