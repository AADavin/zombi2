"""ALElite — ALE-style reconciliation likelihoods for ZOMBI2 scenarios.

ALElite computes the *marginal* reconciliation likelihood ``P(gene tree | species tree,
DTL rates)`` — the quantity ALE (Amalgamated Likelihood Estimation) reports — for a
**single, perfectly known** gene tree. In ZOMBI2 the gene trees and species tree are exact
(no CCP / no gene-tree uncertainty), so the "amalgamation" over a sample of gene trees
collapses to one tree at probability 1, and the likelihood is a plain dynamic-programming
sum over every reconciliation of that one tree against the species tree.

Three model engines:

* :class:`UndatedDTL` — the ALEml_undated / GeneRax ``UndatedDTL`` model: per-branch event
  *odds* (no dates), transfers to any branch. Standard and fast; not a faithful match to a
  dated simulator, but it is what people actually run on real data.
* **reldated** (:func:`reldated_loglik`) — the undated per-branch model, but a transfer may
  only land on a branch that overlaps the donor **in time** (using the species-tree dates),
  forbidding time-inconsistent transfers without full slicing.
* :class:`DatedDTL` — the Szöllősi-2013 dated model: the species tree is time-sliced and
  transfers go only to lineages alive at the same instant, exactly matching ZOMBI2's
  ``_choose_recipient``. This is the faithful likelihood for a ZOMBI2 scenario; its rates are
  the per-unit-time δ/τ/λ a simulation was generated under.

The convenience wrapper :func:`reconciliation_likelihood` is the one-call entry point: hand it
a ZOMBI2 :class:`~zombi2.genomes.reconciliation.Reconciliation` (or an explicit gene
tree + species tree) and DTL rates and it returns the log-likelihood of the family's **extant**
gene tree under the chosen model.

Design: the package is deliberately self-contained (its only ZOMBI2 seams are
:func:`zombi2.tree.read_newick` and the species-tree node API used by
:meth:`SpeciesTree.from_tree`), so the numeric core can be lifted out — or reimplemented in
Rust — without dragging in the simulator.
"""

from __future__ import annotations

from .species import SpeciesTree
from .genetree import GeneTree
from .undated import UndatedDTL, undated_loglik, reldated_loglik, undated_joint_loglik
from .dated import DatedDTL, dated_loglik, dated_extinction, dated_joint_loglik
from .scoring import FamilyScore, score_reconciliations, write_scores_tsv

__all__ = [
    "reconciliation_likelihood",
    "ReconciliationLikelihood",
    "SpeciesTree", "GeneTree",
    "UndatedDTL", "undated_loglik", "reldated_loglik", "undated_joint_loglik",
    "DatedDTL", "dated_loglik", "dated_extinction", "dated_joint_loglik",
    "FamilyScore", "score_reconciliations", "write_scores_tsv",
]

_MODELS = ("dated", "undated", "reldated")


def reconciliation_likelihood(reconciliation=None, species_tree=None, *,
                              gene_tree=None,
                              duplication: float = 0.0, transfer: float = 0.0, loss: float = 0.0,
                              model: str = "dated", origination: str = "root",
                              n_steps: int = 100) -> float:
    """Marginal reconciliation log-likelihood ``P(gene tree | species tree, DTL rates)``.

    The single-call entry point to ALElite. Two ways to supply the gene tree:

    * pass a ZOMBI2 ``reconciliation`` (a
      :class:`~zombi2.genomes.reconciliation.Reconciliation`, e.g. one value of
      ``Genomes.reconciliations()``) — its **extant** (survivors-only) tree is scored; or
    * pass ``gene_tree`` explicitly, either a :class:`GeneTree` or a reconciled-extant Newick
      string (leaf labels ``"<species>|<gid>"``).

    ``species_tree`` is a :class:`zombi2.tree.Tree` (the extant/reconstructed species tree) or a
    prebuilt :class:`SpeciesTree`. ``model`` selects the engine:

    * ``"dated"`` (default) — the faithful time-sliced likelihood; ``duplication``/``transfer``/
      ``loss`` are per-unit-time δ/τ/λ (ZOMBI2's native units, directly comparable to a
      simulation's rates). ``n_steps`` sets the per-slice integration resolution.
    * ``"undated"`` — the ALEml_undated / GeneRax model; rates are per-branch odds.
    * ``"reldated"`` — undated with transfers restricted to time-overlapping branches.

    ``origination`` is ``"root"`` (the family is present on the root branch — exact for ZOMBI2's
    root-seeded families) or ``"uniform"`` (root gene node averaged over all branches). Returns a
    finite ``float`` (``<= 0``), or ``-inf`` for a reconciliation impossible under the rates.
    """
    if model not in _MODELS:
        raise ValueError(f"unknown model {model!r}; choose from {_MODELS}")
    if (reconciliation is None) == (gene_tree is None):
        raise ValueError("pass exactly one of `reconciliation` or `gene_tree`")
    if species_tree is None:
        raise ValueError("species_tree is required")

    sp = species_tree if isinstance(species_tree, SpeciesTree) else SpeciesTree.from_tree(species_tree)

    if reconciliation is not None:
        gt = GeneTree.from_reconciliation(reconciliation)
    elif isinstance(gene_tree, GeneTree):
        gt = gene_tree
    else:
        gt = GeneTree.from_newick(gene_tree)

    if model == "dated":
        return dated_loglik(gt, sp, DatedDTL(duplication, transfer, loss),
                            origination=origination, n_steps=n_steps)
    transfers = "dated" if model == "reldated" else "global"
    return undated_loglik(gt, sp, UndatedDTL(duplication, transfer, loss),
                          origination=origination, transfers=transfers)


class ReconciliationLikelihood:
    """A reusable ALElite scorer: fix a species tree, DTL rates, and model once, then score
    many gene families against the same background.

    ``species_tree`` is a :class:`zombi2.tree.Tree` (or a prebuilt :class:`SpeciesTree`);
    ``model`` is ``"dated"`` (rates are per-unit-time δ/τ/λ), ``"undated"``, or ``"reldated"``
    (rates are per-branch odds). See :func:`reconciliation_likelihood` for the parameter
    meanings — this class carries the same settings but avoids rebuilding the species-tree index
    on every call.

    Example::

        scorer = ReconciliationLikelihood(genomes.species_tree,
                                          duplication=0.1, transfer=0.05, loss=0.15)
        for recon in genomes.reconciliations().values():
            if recon.extant is not None:
                print(scorer.score(recon))
    """

    def __init__(self, species_tree, *, duplication: float = 0.0, transfer: float = 0.0,
                 loss: float = 0.0, model: str = "dated", origination: str = "root",
                 n_steps: int = 100):
        if model not in _MODELS:
            raise ValueError(f"unknown model {model!r}; choose from {_MODELS}")
        self.sp = (species_tree if isinstance(species_tree, SpeciesTree)
                   else SpeciesTree.from_tree(species_tree))
        self.duplication = duplication
        self.transfer = transfer
        self.loss = loss
        self.model = model
        self.origination = origination
        self.n_steps = n_steps

    def score(self, reconciliation=None, *, gene_tree=None) -> float:
        """Log-likelihood of one family's extant gene tree. Pass a ``reconciliation`` or an
        explicit ``gene_tree`` (a :class:`GeneTree` or reconciled-extant Newick)."""
        return reconciliation_likelihood(
            reconciliation, self.sp, gene_tree=gene_tree,
            duplication=self.duplication, transfer=self.transfer, loss=self.loss,
            model=self.model, origination=self.origination, n_steps=self.n_steps,
        )

    def score_all(self, reconciliations) -> list["FamilyScore"]:
        """Score every extant family in a ``{family: Reconciliation}`` mapping (e.g.
        ``Genomes.reconciliations()``) under this scorer's model, returning
        :class:`FamilyScore` rows."""
        return score_reconciliations(
            self.sp, reconciliations, self.duplication, self.transfer, self.loss,
            models=(self.model,), origination=self.origination, n_steps=self.n_steps,
        )

    def __repr__(self) -> str:
        return (f"ReconciliationLikelihood(model={self.model!r}, "
                f"duplication={self.duplication}, transfer={self.transfer}, loss={self.loss})")
