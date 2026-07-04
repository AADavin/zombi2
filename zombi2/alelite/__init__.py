"""ALElite — ALE-style reconciliation likelihoods for ZOMBI2 scenarios.

ALElite computes the *marginal* reconciliation likelihood ``P(gene tree | species tree,
DTL rates)`` — the quantity ALE (Amalgamated Likelihood Estimation) reports — for a
**single, perfectly known** gene tree. In ZOMBI2 the gene trees and species tree are exact
(no CCP / no gene-tree uncertainty), so the "amalgamation" over a sample of gene trees
collapses to one tree at probability 1, and the likelihood is a plain dynamic-programming
sum over every reconciliation of that one tree against the species tree.

Two model engines:

* :class:`UndatedDTL` — the ALEml_undated / GeneRax ``UndatedDTL`` model: per-branch event
  *odds* (no dates), transfers to any branch. Standard and fast; not a faithful match to a
  dated simulator, but it is what people actually run on real data.
* ``DatedDTL`` (next) — the Szöllősi-2013 dated model: the species tree is time-sliced and
  transfers go only to lineages alive at the same instant, exactly matching ZOMBI2's
  ``_choose_recipient``. This is the faithful likelihood for a ZOMBI2 scenario.

Design: the package is deliberately self-contained (its only ZOMBI2 seam is
:meth:`SpeciesTree.from_tree`), so the numeric core can be lifted out — or reimplemented in
Rust — without dragging in the simulator.
"""

from __future__ import annotations

from .species import SpeciesTree
from .genetree import GeneTree
from .undated import UndatedDTL, undated_loglik, reldated_loglik, undated_joint_loglik
from .dated import DatedDTL, dated_loglik, dated_extinction, dated_joint_loglik

__all__ = [
    "SpeciesTree", "GeneTree",
    "UndatedDTL", "undated_loglik", "reldated_loglik", "undated_joint_loglik",
    "DatedDTL", "dated_loglik", "dated_extinction", "dated_joint_loglik",
]
