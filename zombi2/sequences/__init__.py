"""Sequences — level 3: a sequence evolving inside a gene, along its gene tree.

A sequence lives **inside a gene**, so it sees the species tree only through its gene tree
(``SPEC §1``): :func:`simulate_sequences` takes a **genome run** (a
:class:`~zombi2.genomes.GenomesResult`) and evolves one sequence down each family's *complete* gene
tree under a substitution **model** (the menu — nucleotide ``jc69`` · ``k80`` · ``hky85`` · ``gtr``,
or protein ``poisson`` · ``jtt`` · ``dayhoff`` · ``wag`` · ``lg``; :mod:`.substitution_models`) and a
substitution **rate** (``scope(base) × modifiers``; ``SPEC §5``). Sequences are **target-only** in
v1 — nothing drives *out* of a sequence yet (``SPEC §10``).

The whole genome run is required, not just its gene trees, because a level below reads the level
above: the **species tree** is what the lineage clock rides (one rate per species branch, shared by
every family passing through it — ``SPEC §5``) and what the ``species_phylogram`` is drawn on. Bare
gene trees would run, but silently without either, so they are rejected.

``substitution`` is a per-site rate (a bare number, default ``1.0``: a gene-tree branch of ``Δt`` time
gets ``substitution · Δt`` substitutions/site — the **strict clock**), optionally times a **lineage
clock**: ``substitution = 1.0 * mod.ByLineage(spread=)`` is the uncorrelated ("relaxed") clock, one
i.i.d. rate multiplier drawn per **species lineage** and shared by every gene passing through it
(``SPEC §5``, the by-lineage rate modifier). The other clocks (``FromParent`` drift, ``Markov`` hops),
the per-family ``ByFamily`` speed, across-site ``+Γ``, codon models, real-genome-at-root, and the
``record=`` memory dial are named later slices; each is a pure addition.

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
from ..progress import progress_bar
from .evolution import evolve_gene_tree
from .substitution_models import SubstitutionModel, decode

_WRITE_OUTPUTS = ("alignments", "ancestral", "founding", "phylograms", "species_phylogram")

#: The rate grammar this level wires (SPEC §5) — read by the engine gate in :func:`simulate_sequences`
#: and by the CLI's help, so a modifier is never advertised without being implemented. ``ByLineage``
#: on the substitution rate *is* the uncorrelated ("relaxed") lineage clock.
WIRED_MODIFIERS = (ByLineage,)


@dataclass
class SequencesResult:
    """What :func:`simulate_sequences` returns.

    - ``alignments`` — ``{family: {g<copy>: sequence}}``: the observable gene alignment, one entry per
      **extant** gene-tree tip, keyed by its (unique, per-segment) gene id — the same labels as the
      gene tree's / phylogram's Newick leaves. Empty for a family with no surviving copy.
    - ``ancestral`` — ``{family: {g<copy>: sequence}}``: the reconstructed sequence at every
      **internal** gene-tree node, keyed by its gene id (which includes the family's root gene).
    - ``founding`` — ``{family: sequence}``: the sequence the family started with, at its
      **origination** — the state the phylogram's root branch leads *from*. It is drawn from the
      model's stationary frequencies and then evolves across the stem into the root gene's sequence,
      so it is not the same string as ``ancestral[family]["g<root copy>"]`` unless the stem is empty.
      Kept out of ``ancestral`` on purpose: those keys pair one-to-one with phylogram nodes, and the
      origination is a point on a branch, not a node.
    - ``phylograms`` — ``{family: {"complete": newick, "extant": newick | None}}``: each gene tree with
      branch lengths in **substitutions/site** (``base × lineage-clock × Δt``) — the ground-truth tree
      behind each alignment. **Every** node is labelled by its gene id ``g<copy>``, so the tips match
      the ``alignments`` keys and the internal nodes match the ``ancestral`` keys (the phylogram pairs
      one-to-one with the sequences). ``"extant"`` is ``None`` for a family with no survivor.
    - ``species_phylogram`` — ``{"complete": newick, "extant": newick | None}``: the **species tree**
      with branch lengths in substitutions/site — the molecular clock made visible (which lineages ran
      hot / cold). Always present: a run always comes from a genome run, which carries its tree.
    - ``seed`` — the run's seed.
    """

    alignments: dict[int, dict[str, str]]
    ancestral: dict[int, dict[str, str]]
    founding: dict[int, str]
    phylograms: dict[int, dict[str, str | None]]
    species_phylogram: dict[str, str | None]
    seed: int | None

    def write(self, directory, outputs=("alignments", "phylograms")) -> None:
        """Write chosen ``outputs`` to ``directory`` (created if needed):

        - ``"alignments"`` → ``fam<family>.fasta`` (skipped for empty families).
        - ``"ancestral"`` → ``sequences_ancestral_fam<family>.fasta``.
        - ``"founding"`` → ``sequences_founding.fasta``, one record ``fam<family>`` per family: the
          sequence each family originated with, before its stem.
        - ``"phylograms"`` → ``phylogram_fam<family>_{complete,extant}.nwk`` (subs/site).
        - ``"species_phylogram"`` → ``sequences_species_phylogram_{complete,extant}.nwk`` (subs/site).
        """
        unknown = [o for o in outputs if o not in _WRITE_OUTPUTS]
        if unknown:
            raise ValueError(f"unknown write outputs {unknown}; choose from {list(_WRITE_OUTPUTS)}")
        d = pathlib.Path(directory)
        d.mkdir(parents=True, exist_ok=True)
        if "alignments" in outputs:
            for fam, aln in self.alignments.items():
                if aln:
                    (d / f"fam{fam}.fasta").write_text(_fasta(aln))
        if "ancestral" in outputs:
            for fam, anc in self.ancestral.items():
                if anc:
                    (d / f"sequences_ancestral_fam{fam}.fasta").write_text(_fasta(anc))
        if "founding" in outputs and self.founding:
            (d / "sequences_founding.fasta").write_text(
                _fasta({f"fam{fam}": seq for fam, seq in sorted(self.founding.items())}))
        if "phylograms" in outputs:
            for fam, ph in self.phylograms.items():
                (d / f"phylogram_fam{fam}_complete.nwk").write_text(ph["complete"] + "\n")
                if ph["extant"] is not None:
                    (d / f"phylogram_fam{fam}_extant.nwk").write_text(ph["extant"] + "\n")
        if "species_phylogram" in outputs:
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
    over. Collected from the gene trees rather than the species tree, so the draws depend only on the
    branches genes actually pass through; every branch that needs a clock value has its species
    branch present as some node's ``species`` (a branch no gene crossed keeps the factor 1.0)."""
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
    family's **origination** (``base × clock[species] × Δt`` summed along the path). Feeding it to
    ``GeneTree.to_newick`` then emits a *phylogram* (branch lengths in subs/site); and because its
    prune-to-extant merges branches by that same cumulative measure, a suppressed branch spanning
    several species branches gets the **sum** of its pieces for free — the exact trick the chronogram
    uses with time.

    Counting from origination rather than from the root is what gives the root its own branch: the
    founding gene evolves across the stem, so the scaled root sits ``base × clock × stem`` in, not
    at zero."""
    root = gt.complete
    stem = rate_base * _clock_factor(clock, root.species) * (root.time - gt.origination)
    scaled_root = GeneNode(root.kind, root.species, stem, root.copy)
    stack = [(root, scaled_root)]
    while stack:
        onode, snode = stack.pop()
        for ochild in onode.children:
            blen = rate_base * _clock_factor(clock, ochild.species) * (ochild.time - onode.time)
            schild = GeneNode(ochild.kind, ochild.species, snode.time + blen, ochild.copy)
            snode.children.append(schild)
            stack.append((ochild, schild))
    return GeneTree(gt.family, scaled_root, 0.0)   # origination is the zero of the scaled measure


def _gene_newick(root: GeneNode) -> str:
    """Newick of a (scaled) gene tree labelling **every** node — leaf and internal — by its gene id
    ``g<copy>``, so the tips match the ``alignments`` keys and the internal nodes match the
    ``ancestral`` keys (both keyed ``g<copy>``): the phylogram pairs one-to-one with the sequences.
    Branch lengths are node-``time`` differences (substitutions/site on a scaled tree). The root's
    parent measure is 0 — the family's origination — so it carries the stem like every other branch.
    Iterative — gene trees run past CPython's recursion guard, so recursion would crash on deep trees."""
    stack: list[list] = [[root, 0.0, 0, []]]       # [node, parent_time, next_child, child_strings]
    result = ""
    while stack:
        frame = stack[-1]
        node, parent_time, ci, parts = frame
        if ci < len(node.children):
            frame[2] = ci + 1
            stack.append([node.children[ci], node.time, 0, []])
            continue
        bl = f":{node.time - parent_time:.6g}"
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


def simulate_sequences(genomes, *, model: SubstitutionModel, length: int,
                       substitution=1.0, seed=None, progress=False) -> SequencesResult:
    """Evolve one sequence down each family's gene tree under a substitution ``model``.

    ``genomes`` is a **genome run** — the :class:`~zombi2.genomes.GenomesResult` that
    ``genomes.simulate_genomes_unordered(...)`` returned. Its ``gene_trees`` are what the sequences
    evolve along and its ``complete_tree`` is the species tree the lineage clock rides; bare gene
    trees are rejected (they would run, but with no clock and no species phylogram — a silent
    degradation). Each family's *complete* gene tree is evolved, so the true history is complete and
    ancestral sequences exist for extinct/lost lineages too; the observable ``alignments`` are the
    extant tips.

    ``model`` is a substitution model from the menu (:mod:`.substitution_models`) — nucleotide
    ``jc69`` · ``k80`` · ``hky85`` · ``gtr``, or protein ``poisson`` · ``jtt`` · ``dayhoff`` ·
    ``wag`` · ``lg``; its alphabet is what the sequences are written in (``ACGT`` or the 20 amino
    acids). ``length`` is the number of sites. ``substitution`` is the per-site substitution
    rate (default ``1.0``): a branch of ``Δt`` time accrues ``substitution · Δt`` substitutions/site.
    The root sequence of each family is drawn from the model's stationary frequencies. Deterministic
    given ``seed``.

    ``substitution`` may carry a **lineage clock**: ``1.0 * mod.ByLineage(spread=)`` draws one i.i.d.
    rate multiplier per species lineage (**shared across families**, drawn once before evolving) and
    rescales each gene-tree branch by the clock of the species branch it sits on. Any other modifier
    (the ``FromParent`` / ``Markov`` clocks, the ``ByFamily`` per-family speed, ``+Γ``) or a
    non-``PerSite`` scope is a later slice and raises.

    The result carries the **phylograms** the sequences were drawn along — each gene tree and the
    species tree, with branch lengths converted from time to substitutions/site by the same
    ``base × clock × Δt``.
    """
    if not isinstance(genomes, GenomesResult):
        raise TypeError(
            f"the sequence level runs on a genome run, got {type(genomes).__name__} — pass the "
            "GenomesResult that genomes.simulate_genomes_unordered(...) returned: the whole run, "
            "not its .gene_trees (the ordered / nucleotide results are a later slice). A sequence "
            "lives inside a gene, but its clock rides the *species* branch that gene sits on — one "
            "draw per lineage, shared by every family — so the run needs the species tree too."
        )
    species_tree = genomes.complete_tree
    gene_trees = genomes.gene_trees
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
    founding: dict[int, str] = {}
    phylograms: dict[int, dict[str, str | None]] = {}
    bar = progress_bar(len(gene_trees), "sequences", unit="family", enabled=progress)
    for family in sorted(gene_trees):  # sorted for reproducibility given the seed
        bar.update()
        gt = gene_trees[family]
        states, founding_states = evolve_gene_tree(gt.complete, model, length, rate_base, clock, rng,
                                                   gt.origination)
        alignments[family], ancestral[family] = _split(gt, states, model)
        founding[family] = decode(founding_states, model.alphabet)
        scaled = _scaled_gene_tree(gt, rate_base, clock)  # branch lengths in subs/site
        ext = scaled.extant
        phylograms[family] = {"complete": _gene_newick(scaled.complete),
                              "extant": _gene_newick(ext) if ext is not None else None}

    bar.close()

    sp_scaled = _scaled_species_tree(species_tree, rate_base, clock)   # the clock made visible
    sp_extant = prune(sp_scaled, keep="extant")
    species_phylogram = {"complete": sp_scaled.to_newick(),
                         "extant": sp_extant.to_newick() if sp_extant is not None else None}

    return SequencesResult(alignments, ancestral, founding, phylograms, species_phylogram, seed)


# The substitution-model menu is reached through its own module — the one canonical path,
# like `zombi2.rates.scope` / `zombi2.rates.modifiers` — never re-exported here:
#     from zombi2.sequences import substitution_models as sm;  sm.hky85(2.0)
__all__ = ["simulate_sequences", "SequencesResult"]
