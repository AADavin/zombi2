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
sequence at every **internal** node), ``.phylograms`` (each gene tree with branch lengths in
substitutions/site — the ground-truth tree behind each alignment), ``.species_phylogram`` (the species
tree scaled the same way — the molecular clock made visible), and ``.seed``. Genuine substitution
``.events`` are the deferred opt-in ``record=`` slice, not the default spine (a substitution log is not
compact the way the speciation / D-T-L-O logs are).
"""

from __future__ import annotations

import pathlib
from dataclasses import dataclass

import numpy as np

from ..genomes import GenomesResult
from ..genomes.gene_trees import GeneNode, GeneTree
from ..rates.modifiers import ByLineage
from ..rates.rate import as_rate
from ..rates.scope import PerSite
from ..species import Node, Tree, prune
from .evolution import evolve_gene_tree
from .substitution_models import SubstitutionModel, decode

_WRITE_OUTPUTS = ("alignments", "ancestral", "phylograms", "species_phylogram")  # the write vocabulary


@dataclass
class SequencesResult:
    """What :func:`simulate_sequences` returns.

    - ``alignments`` — ``{family: {g<copy>: sequence}}``: the observable gene alignment, one entry per
      **extant** gene-tree tip, keyed by its (unique, per-segment) gene id — the same labels as the
      gene tree's / phylogram's Newick leaves. Empty for a family with no surviving copy.
    - ``ancestral`` — ``{family: {g<copy>: sequence}}``: the reconstructed sequence at every
      **internal** gene-tree node, keyed by its gene id (which includes the family's root gene).
    - ``phylograms`` — ``{family: {"complete": newick, "extant": newick | None}}``: each gene tree with
      branch lengths in **substitutions/site** (``base × lineage-clock × Δt``) — the ground-truth tree
      behind each alignment. **Every** node is labelled by its gene id ``g<copy>``, so the tips match
      the ``alignments`` keys and the internal nodes match the ``ancestral`` keys (the phylogram pairs
      one-to-one with the sequences). ``"extant"`` is ``None`` for a family with no survivor.
    - ``species_phylogram`` — ``{"complete": newick, "extant": newick | None}``, or ``None`` when the
      run was given bare gene trees rather than a ``GenomesResult``: the **species tree** with branch
      lengths in substitutions/site — the molecular clock made visible (which lineages ran hot / cold).
    - ``seed`` — the run's seed.
    """

    alignments: dict[int, dict[str, str]]
    ancestral: dict[int, dict[str, str]]
    phylograms: dict[int, dict[str, str | None]]
    species_phylogram: dict[str, str | None] | None
    seed: int | None

    def write(self, directory, outputs=("alignments", "phylograms")) -> None:
        """Write chosen ``outputs`` to ``directory`` (created if needed):

        - ``"alignments"`` → ``sequences_alignment_fam<family>.fasta`` (skipped for empty families).
        - ``"ancestral"`` → ``sequences_ancestral_fam<family>.fasta``.
        - ``"phylograms"`` → ``sequences_phylogram_fam<family>_{complete,extant}.nwk`` (subs/site).
        - ``"species_phylogram"`` → ``sequences_species_phylogram_{complete,extant}.nwk`` (subs/site;
          nothing written when the run was given bare gene trees).
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
        if "phylograms" in outputs:
            for fam, ph in self.phylograms.items():
                (d / f"sequences_phylogram_fam{fam}_complete.nwk").write_text(ph["complete"] + "\n")
                if ph["extant"] is not None:
                    (d / f"sequences_phylogram_fam{fam}_extant.nwk").write_text(ph["extant"] + "\n")
        if "species_phylogram" in outputs and self.species_phylogram is not None:
            sp = self.species_phylogram
            (d / "sequences_species_phylogram_complete.nwk").write_text(sp["complete"] + "\n")
            if sp["extant"] is not None:
                (d / "sequences_species_phylogram_extant.nwk").write_text(sp["extant"] + "\n")


def _fasta(records: dict[str, str], width: int = 70) -> str:
    """Serialise ``{name: sequence}`` to FASTA text (sequences wrapped at ``width`` columns)."""
    lines: list[str] = []
    for name, seq in records.items():
        lines.append(f">{name}")
        lines.extend(seq[i:i + width] for i in range(0, len(seq), width))
    return "\n".join(lines) + "\n"


def _split(gene_tree, states_by_id: dict[int, np.ndarray],
           model: SubstitutionModel) -> tuple[dict[str, str], dict[str, str]]:
    """Label one family's evolved nodes by their **gene id** and split them into the extant-tip
    alignment and the internal-node ancestral set. Gene ids are per-segment (each node has a unique
    ``copy``), so ``g<copy>`` uniquely names every node and — for the tips — matches the gene tree's
    and phylogram's Newick leaves, pairing the alignment with its tree."""
    alignment: dict[str, str] = {}
    ancestral: dict[str, str] = {}
    stack = [gene_tree.complete]
    while stack:
        node = stack.pop()
        seq = decode(states_by_id[id(node)], model.alphabet)
        if node.is_leaf:
            if node.kind == "extant":
                alignment[f"g{node.copy}"] = seq
        else:
            ancestral[f"g{node.copy}"] = seq
        stack.extend(node.children)
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


def _clock_factor(clock, species: int) -> float:
    """The lineage clock on a species branch — 1.0 under the strict clock (``clock is None``) or for a
    branch no gene passed through (so none was drawn for it)."""
    return 1.0 if clock is None else clock.get(species, 1.0)


def _scaled_gene_tree(gt: GeneTree, rate_base: float, clock) -> GeneTree:
    """A copy of the gene tree whose node ``time`` holds the cumulative **substitutions/site** from the
    root (``base × clock[species] × Δt`` summed along the path). Feeding it to ``GeneTree.to_newick``
    then emits a *phylogram* (branch lengths in subs/site); and because its prune-to-extant merges
    branches by that same cumulative measure, a suppressed branch spanning several species branches gets
    the **sum** of its pieces for free — the exact trick the chronogram uses with time."""
    root = gt.complete
    scaled_root = GeneNode(root.kind, root.species, 0.0, root.copy)
    stack = [(root, scaled_root)]
    while stack:
        onode, snode = stack.pop()
        for ochild in onode.children:
            blen = rate_base * _clock_factor(clock, ochild.species) * (ochild.time - onode.time)
            schild = GeneNode(ochild.kind, ochild.species, snode.time + blen, ochild.copy)
            snode.children.append(schild)
            stack.append((ochild, schild))
    return GeneTree(gt.family, scaled_root)


def _gene_newick(root: GeneNode) -> str:
    """Newick of a (scaled) gene tree labelling **every** node — leaf and internal — by its gene id
    ``g<copy>``, so the tips match the ``alignments`` keys and the internal nodes match the
    ``ancestral`` keys (both keyed ``g<copy>``): the phylogram pairs one-to-one with the sequences.
    Branch lengths are node-``time`` differences (substitutions/site on a scaled tree). Iterative —
    gene trees run past CPython's recursion guard, so recursion would crash on deep trees."""
    stack: list[list] = [[root, None, 0, []]]      # [node, parent_time, next_child, child_strings]
    result = ""
    while stack:
        frame = stack[-1]
        node, parent_time, ci, parts = frame
        if ci < len(node.children):
            frame[2] = ci + 1
            stack.append([node.children[ci], node.time, 0, []])
            continue
        bl = "" if parent_time is None else f":{node.time - parent_time:.6g}"
        s = f"g{node.copy}{bl}" if node.is_leaf else f"({','.join(parts)})g{node.copy}{bl}"
        stack.pop()
        if stack:
            stack[-1][3].append(s)
        else:
            result = s
    return result + ";"


def _scaled_species_tree(tree: Tree, rate_base: float, clock) -> Tree:
    """A copy of the species tree whose branch lengths are **substitutions/site** (``base ×
    clock[branch] × Δt``). Node times become the cumulative subs/site from the root, so
    ``Tree.to_newick`` / ``prune`` emit and merge the phylogram exactly as they do a dated tree."""
    scaled: dict[int, Node] = {}
    scaled_end: dict[int, float] = {}
    order: list[int] = []
    stack = [tree.root]
    while stack:  # pre-order: a parent is visited before its children
        i = stack.pop()
        order.append(i)
        if tree.nodes[i].children is not None:
            stack.extend(tree.nodes[i].children)
    for i in order:
        nd = tree.nodes[i]
        blen = rate_base * _clock_factor(clock, i) * (nd.end_time - nd.birth_time)
        start = 0.0 if nd.parent is None else scaled_end[nd.parent]
        scaled_end[i] = start + blen
        scaled[i] = Node(i, nd.parent, start, start + blen, nd.children, nd.fate)
    return Tree(scaled, tree.root)


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

    The result carries the **phylograms** the sequences were drawn along — each gene tree, and (when a
    ``GenomesResult`` is passed) the species tree, with branch lengths converted from time to
    substitutions/site by the same ``base × clock × Δt``.
    """
    species_tree = None
    if isinstance(gene_trees, GenomesResult):
        species_tree = gene_trees.complete_tree
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
    phylograms: dict[int, dict[str, str | None]] = {}
    for family in sorted(gene_trees):  # sorted for reproducibility given the seed
        gt = gene_trees[family]
        states = evolve_gene_tree(gt.complete, model, length, rate_base, clock, rng)
        alignments[family], ancestral[family] = _split(gt, states, model)
        scaled = _scaled_gene_tree(gt, rate_base, clock)  # branch lengths in subs/site
        ext = scaled.extant
        phylograms[family] = {"complete": _gene_newick(scaled.complete),
                              "extant": _gene_newick(ext) if ext is not None else None}

    species_phylogram = None
    if species_tree is not None:
        sp_scaled = _scaled_species_tree(species_tree, rate_base, clock)
        sp_extant = prune(sp_scaled, keep="extant")
        species_phylogram = {"complete": sp_scaled.to_newick(),
                             "extant": sp_extant.to_newick() if sp_extant is not None else None}

    return SequencesResult(alignments, ancestral, phylograms, species_phylogram, seed)


# The substitution-model menu is reached through its own module — the one canonical path,
# like `zombi2.rates.scope` / `zombi2.rates.modifiers` — never re-exported here:
#     from zombi2.sequences import substitution_models as sm;  sm.hky85(2.0)
__all__ = ["simulate_sequences", "SequencesResult"]
