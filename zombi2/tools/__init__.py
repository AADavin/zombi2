"""ZOMBI2 tools — the analysis / interop complement to the simulator.

``zombi2`` proper is a **simulation** core: it generates species trees, gene families,
traits, and sequences. ``zombi2.tools`` is the adjacent layer of **bounded, validated
computations on ZOMBI2 outputs** — likelihoods, scores, statistics, distances, and format
conversions. The bar for a tool is deliberately narrow: a computation with a *right answer*
that is **bounded, built, and validated**. Open-ended methodologies (ABC / MCMC inference
frameworks) are out of scope here and remain deferred to a future inference release.

Unlike :mod:`zombi2.experimental`, everything in ``zombi2.tools`` is **stable**: it has
cleared validation and its API will not change without notice. But like ``experimental`` it
is a **distinct, labelled surface** — nothing here is re-exported into the top-level
``zombi2`` namespace. Import tools explicitly::

    from zombi2.tools import reconciliation_likelihood

The tools so far are **ALElite** (:mod:`zombi2.tools.reconciliation`): the ALE-style marginal
reconciliation likelihood ``P(gene tree | species tree, DTL rates)`` of a simulated gene
family, validated against closed-form oracles; and **RED** (:mod:`zombi2.tools.red`): the
Relative Evolutionary Divergence of every node of a tree (Parks et al. 2018), GTDB's rate-
normalised relative-age scale. See ``docs/tools/`` for the layer's scope and roadmap.

ALElite's interop complement is **reconparser** (:mod:`zombi2.tools.reconparser`): parsers that
*read* the output of the established reconciliation programs (ALE, AleRax) into ``ete3`` trees
and ``pandas`` DataFrames — the bridge for comparing a real reconciliation against a ZOMBI2
simulation. It needs the optional ``reconparser`` extra (``pip install 'zombi2[reconparser]'``)
and, like the others, is imported from its own submodule rather than re-exported here::

    from zombi2.tools.reconparser import ALEParser, AleRaxRun, AleRaxFamily
"""

from __future__ import annotations

from .reconciliation import (
    reconciliation_likelihood,
    ReconciliationLikelihood,
    SpeciesTree,
    GeneTree,
    UndatedDTL,
    DatedDTL,
    undated_loglik,
    reldated_loglik,
    dated_loglik,
    dated_extinction,
    undated_joint_loglik,
    dated_joint_loglik,
    FamilyScore,
    score_reconciliations,
    write_scores_tsv,
)
from .red import relative_evolutionary_divergence

__all__ = [
    "reconciliation_likelihood",
    "ReconciliationLikelihood",
    "SpeciesTree",
    "GeneTree",
    "UndatedDTL",
    "DatedDTL",
    "undated_loglik",
    "reldated_loglik",
    "dated_loglik",
    "dated_extinction",
    "undated_joint_loglik",
    "dated_joint_loglik",
    "FamilyScore",
    "score_reconciliations",
    "write_scores_tsv",
    "relative_evolutionary_divergence",
]
