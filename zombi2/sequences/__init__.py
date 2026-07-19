"""Sequences — level 3: a sequence evolving inside a gene, along its gene tree.

A sequence lives **inside a gene**, so it sees the species tree only through its gene tree
(``SPEC §1``): :func:`simulate_sequences` takes the ``{family: GeneTree}`` a genome run produced
(``GenomesResult.gene_trees``) and evolves one sequence down each family's *complete* gene tree under
a substitution **model** (the menu — :func:`jc69` · :func:`k80` · :func:`hky85` · :func:`gtr`) and a
substitution **rate** (``scope(base) × modifiers``; ``SPEC §5``). Sequences are **target-only** in
v1 — nothing drives *out* of a sequence yet (``SPEC §10``).

``substitution`` is a per-site rate (a bare number, default ``1.0``: a gene-tree branch of ``Δt`` time
gets ``substitution · Δt`` substitutions/site — the **strict clock**), optionally times a **lineage
clock**: ``substitution = 1.0 * mod.ByLineage(spread=)`` is the uncorrelated ("relaxed") clock, one
i.i.d. rate multiplier drawn per **species lineage** and shared by every gene passing through it
(``SPEC §5``, the by-lineage rate modifier). The other clocks (``FromParent`` drift, ``Markov`` hops),
the per-family ``ByFamily`` speed, across-site ``+Γ``, protein/codon models, real-genome-at-root, the
``record=`` memory dial, and the CLI are named later slices; each is a pure addition.

The result is a :class:`SequencesResult` bundle mirroring the other levels (``result-api.md``):
``.alignments`` (the observable sequence at every **extant** tip), ``.ancestral`` (the reconstructed
sequence at every **internal** node), and ``.seed``. Genuine substitution ``.events`` are the deferred
opt-in ``record=`` slice, not the default spine (a substitution log is not compact the way the
speciation / D-T-L-O logs are).
"""

from __future__ import annotations

import pathlib
from dataclasses import dataclass

import numpy as np

from ..genomes import GenomesResult
from ..rates.modifiers import ByLineage
from ..rates.rate import as_rate
from ..rates.scope import PerSite
from .evolution import evolve_gene_tree
from .substitution_models import SubstitutionModel, decode

_WRITE_OUTPUTS = ("alignments", "ancestral")  # the write vocabulary


@dataclass
class SequencesResult:
    """What :func:`simulate_sequences` returns.

    - ``alignments`` — ``{family: {tip_label: sequence}}``: the observable gene alignment, one entry
      per **extant** gene-tree tip (label ``g<copy>_n<species>``, matching the gene tree's Newick
      leaves). Empty for a family with no surviving copy.
    - ``ancestral`` — ``{family: {node_label: sequence}}``: the reconstructed sequence at every
      **internal** gene-tree node (label ``<kind>_n<species>_i<preorder-index>``, unique within the
      family), including the founding ``origination`` node's root sequence.
    - ``seed`` — the run's seed.
    """

    alignments: dict[int, dict[str, str]]
    ancestral: dict[int, dict[str, str]]
    seed: int | None

    def write(self, directory, outputs=("alignments",)) -> None:
        """Write chosen ``outputs`` to ``directory`` (created if needed), one FASTA per family:

        - ``"alignments"`` → ``sequences_alignment_fam<family>.fasta`` (skipped for empty families).
        - ``"ancestral"`` → ``sequences_ancestral_fam<family>.fasta``.
        """
        unknown = [o for o in outputs if o not in _WRITE_OUTPUTS]
        if unknown:
            raise ValueError(f"unknown write outputs {unknown}; choose from {list(_WRITE_OUTPUTS)}")
        d = pathlib.Path(directory)
        d.mkdir(parents=True, exist_ok=True)
        if "alignments" in outputs:
            for fam, aln in self.alignments.items():
                if aln:
                    (d / f"sequences_alignment_fam{fam}.fasta").write_text(_fasta(aln))
        if "ancestral" in outputs:
            for fam, anc in self.ancestral.items():
                if anc:
                    (d / f"sequences_ancestral_fam{fam}.fasta").write_text(_fasta(anc))


def _fasta(records: dict[str, str], width: int = 70) -> str:
    """Serialise ``{name: sequence}`` to FASTA text (sequences wrapped at ``width`` columns)."""
    lines: list[str] = []
    for name, seq in records.items():
        lines.append(f">{name}")
        lines.extend(seq[i:i + width] for i in range(0, len(seq), width))
    return "\n".join(lines) + "\n"


def _split(gene_tree, states_by_id: dict[int, np.ndarray],
           model: SubstitutionModel) -> tuple[dict[str, str], dict[str, str]]:
    """Label one family's evolved nodes and split them into the extant-tip alignment and the
    internal-node ancestral set. Pre-order over the complete tree; the per-node index makes the
    internal labels unique. Extant-tip labels ``g<copy>_n<species>`` are unique by construction
    (a copy id sits in one lineage at the present); a duplicate would be a bug, so it is asserted."""
    alignment: dict[str, str] = {}
    ancestral: dict[str, str] = {}
    idx = 0
    stack = [gene_tree.complete]
    while stack:
        node = stack.pop()
        if node.is_leaf:
            if node.kind == "extant":
                label = f"g{node.copy}_n{node.species}"
                if label in alignment:
                    raise AssertionError(f"duplicate extant tip label {label!r} in family {gene_tree.family}")
                alignment[label] = decode(states_by_id[id(node)], model.alphabet)
        else:
            ancestral[f"{node.kind}_n{node.species}_i{idx}"] = decode(states_by_id[id(node)], model.alphabet)
        idx += 1
        stack.extend(reversed(node.children))
    return alignment, ancestral


def _all_species(gene_trees) -> list[int]:
    """The sorted set of species-branch ids the gene trees touch — the lineages the clock is drawn
    over. Collected from the gene trees (so it works with a bare ``{family: GeneTree}`` too); every
    branch that needs a clock value has its species branch present as some node's ``species``."""
    ids: set[int] = set()
    for gt in gene_trees.values():
        stack = [gt.complete]
        while stack:
            n = stack.pop()
            ids.add(n.species)
            stack.extend(n.children)
    return sorted(ids)


def simulate_sequences(gene_trees, *, model: SubstitutionModel, length: int,
                       substitution=1.0, seed=None) -> SequencesResult:
    """Evolve one sequence down each family's gene tree under a substitution ``model``.

    ``gene_trees`` is a ``{family: GeneTree}`` mapping (a genome run's ``GenomesResult.gene_trees``),
    or a :class:`~zombi2.genomes.GenomesResult` directly (its ``gene_trees`` are used). Each family's
    *complete* gene tree is evolved, so the true history is complete and ancestral sequences exist for
    extinct/lost lineages too; the observable ``alignments`` are the extant tips.

    ``model`` is a substitution model from the menu (:func:`jc69` · :func:`k80` · :func:`hky85` ·
    :func:`gtr`). ``length`` is the number of sites. ``substitution`` is the per-site substitution
    rate (default ``1.0``): a branch of ``Δt`` time accrues ``substitution · Δt`` substitutions/site.
    The root sequence of each family is drawn from the model's stationary frequencies. Deterministic
    given ``seed``.

    ``substitution`` may carry a **lineage clock**: ``1.0 * mod.ByLineage(spread=)`` draws one i.i.d.
    rate multiplier per species lineage (**shared across families**, drawn once before evolving) and
    rescales each gene-tree branch by the clock of the species branch it sits on. Any other modifier
    (the ``FromParent`` / ``Markov`` clocks, the ``ByFamily`` per-family speed, ``+Γ``) or a
    non-``PerSite`` scope is a later slice and raises.
    """
    if isinstance(gene_trees, GenomesResult):
        gene_trees = gene_trees.gene_trees
    if not isinstance(model, SubstitutionModel):
        raise TypeError(f"model must be a SubstitutionModel (e.g. hky85(kappa=2.0)), got {model!r}")
    if isinstance(length, bool) or not isinstance(length, int) or length < 1:
        raise ValueError(f"length must be a positive integer, got {length!r}")
    rate = as_rate(substitution, default_scope=PerSite)
    if not isinstance(rate.scope, PerSite):
        raise ValueError(
            f"substitution has a {type(rate.scope).__name__} scope, but the sequence engine wires only "
            f"PerSite (the default) this slice — drop the scope wrapper or use PerSite(...)."
        )
    clock_mod = None
    if rate.modifiers:
        if len(rate.modifiers) == 1 and isinstance(rate.modifiers[0], ByLineage):
            clock_mod = rate.modifiers[0]
        else:
            offenders = ", ".join(sorted({type(m).__name__ for m in rate.modifiers
                                          if not isinstance(m, ByLineage)}) or ["a second ByLineage"])
            raise ValueError(
                f"substitution carries {offenders}, but this slice wires only a single ByLineage clock "
                "(the uncorrelated lineage clock) — the FromParent / Markov clocks, the ByFamily "
                "per-family speed, and +Γ across-site heterogeneity are later slices."
            )
    rate_base = rate.base

    rng = np.random.default_rng(seed)
    # the lineage clock: one i.i.d. draw per species branch, drawn once here and shared by every
    # family (so a hot species runs hot for all its genes). None ⇒ the strict clock (factor 1).
    clock = None
    if clock_mod is not None:
        clock = {sid: clock_mod.draw(rng) for sid in _all_species(gene_trees)}
    alignments: dict[int, dict[str, str]] = {}
    ancestral: dict[int, dict[str, str]] = {}
    for family in sorted(gene_trees):  # sorted for reproducibility given the seed
        gt = gene_trees[family]
        states = evolve_gene_tree(gt.complete, model, length, rate_base, clock, rng)
        alignments[family], ancestral[family] = _split(gt, states, model)
    return SequencesResult(alignments, ancestral, seed)


# The substitution-model menu is reached through its own module — the one canonical path,
# like `zombi2.rates.scope` / `zombi2.rates.modifiers` — never re-exported here:
#     from zombi2.sequences import substitution_models as sm;  sm.hky85(2.0)
__all__ = ["simulate_sequences", "SequencesResult"]
