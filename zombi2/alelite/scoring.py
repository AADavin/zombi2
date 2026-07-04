"""Per-family ALE reconciliation likelihoods for a simulated gene-family set.

Turns a whole gene-family simulation into a table: for every family with a surviving copy,
the log-likelihood its **extant** reconciled gene tree gets under each ALE model (``dated``,
``undated``, ``reldated``), at given DTL rates. Extinction is built once per model and every
family scored against it — the Rust per-family kernels when the extension is built, the pure
Python engines otherwise.

The rates mean different things per model (they are the *same numbers* passed to each): for
``dated`` they are per-unit-time δ/τ/λ (ZOMBI2's native units, directly comparable to what a
simulation was generated under); for ``undated``/``reldated`` they are per-branch odds.
"""

from __future__ import annotations

from collections import namedtuple
from pathlib import Path

from . import _rust
from .dated import DatedDTL, _DatedEngine
from .genetree import GeneTree
from .species import SpeciesTree
from .undated import UndatedDTL, _transfer_neighbors, _tree_loglik, extinction

#: One row of the score table. ``logliks`` maps model name -> log-likelihood.
FamilyScore = namedtuple("FamilyScore", ["family", "extant_tips", "logliks"])

MODELS = ("dated", "undated", "reldated")


def _score_model(model: str, trees, sp: SpeciesTree, dup, transfer, loss,
                 origination: str, n_steps: int, backend: str) -> list[float]:
    if not trees:
        return []
    use_rust = backend in ("auto", "rust") and _rust.available_family()
    if backend == "rust" and not use_rust:
        raise RuntimeError("backend='rust' requested but the zombi2_core extension is not built")

    if model == "dated":
        if use_rust:
            return list(_rust.dated_family_loglik(trees, sp, dup, transfer, loss, origination, n_steps))
        eng = _DatedEngine(sp, DatedDTL(dup, transfer, loss), n_steps)
        return [eng.gene_loglik(gt, origination) for gt in trees]

    transfers = "dated" if model == "reldated" else "global"
    if use_rust:
        return list(_rust.undated_family_loglik(trees, sp, dup, transfer, loss, origination, transfers))
    m = UndatedDTL(dup, transfer, loss)
    pD, pT, pL, pS = m.probs()
    nb = _transfer_neighbors(sp, transfers)
    E = extinction(sp, m, nb)
    return [_tree_loglik(gt, sp, pD, pT, pS, E, nb, origination) for gt in trees]


def score_reconciliations(species_tree, reconciliations, dup: float, transfer: float, loss: float,
                          *, models=("dated", "undated"), origination: str = "root",
                          n_steps: int = 100, backend: str = "auto") -> list[FamilyScore]:
    """Score every extant family's reconciled gene tree under each ALE ``model``.

    ``species_tree`` is a :class:`zombi2.tree.Tree` (the extant/reconstructed species tree);
    ``reconciliations`` is ``{family: Reconciliation}`` (from ``Genomes.reconciliations()``).
    Families with no surviving copy are skipped (no observable gene tree). Returns a list of
    :class:`FamilyScore` in reconciliation order.
    """
    bad = [m for m in models if m not in MODELS]
    if bad:
        raise ValueError(f"unknown model(s) {bad}; choose from {MODELS}")
    sp = SpeciesTree.from_tree(species_tree)
    fams, trees = [], []
    for fam, recon in reconciliations.items():
        if recon.extant is None:
            continue
        fams.append(fam)
        trees.append(GeneTree.from_reconciliation(recon))
    tips = [sum(nd.is_leaf for nd in gt.nodes) for gt in trees]
    cols = {m: _score_model(m, trees, sp, dup, transfer, loss, origination, n_steps, backend)
            for m in models}
    return [FamilyScore(fams[i], tips[i], {m: cols[m][i] for m in models})
            for i in range(len(fams))]


def write_scores_tsv(rows, path, models=("dated", "undated")) -> None:
    """Write a score table to ``path`` — one row per family, one ``<model>_loglik`` column each."""
    header = "family\textant_copies\t" + "\t".join(f"{m}_loglik" for m in models)
    lines = [header]
    for r in rows:
        vals = "\t".join(f"{r.logliks[m]:.6f}" for m in models)
        lines.append(f"{r.family}\t{r.extant_tips}\t{vals}")
    Path(path).write_text("\n".join(lines) + "\n")
