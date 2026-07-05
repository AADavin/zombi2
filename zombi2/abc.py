"""Approximate Bayesian Computation namespace (scikit-learn-style).

Re-exports the profile-matching / ABC public API so users can write::

    from zombi2.abc import match_profiles, cooccurrence_summary

Every name here is the *same object* as the corresponding top-level
``zombi2`` attribute -- this module is a thin, additive namespace over the
implementation module :mod:`zombi2.matching`; it does not redefine anything.
"""

from __future__ import annotations

from .matching import (
    match_profiles, match_profiles_smc, match_coupled, ABCFit,
    default_summary, default_gene_tree_summary, cooccurrence_summary,
    cooccurrence_features, event_count_summary, frequency_spectrum,
    genome_sizes, copy_number_spectrum,
)

__all__ = [
    "match_profiles", "match_profiles_smc", "match_coupled", "ABCFit",
    "default_summary", "default_gene_tree_summary", "cooccurrence_summary",
    "cooccurrence_features", "event_count_summary", "frequency_spectrum",
    "genome_sizes", "copy_number_spectrum",
]
