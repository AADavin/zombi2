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
i.i.d. rate multiplier drawn per **species lineage** and shared by every gene passing through it, and
``substitution = 1.0 * mod.FromParent(spread=)`` is the **autocorrelated** clock, where the rate drifts
parent→child down the species tree so close relatives run at similar rates (``SPEC §5``). The rest
(``Markov`` hops, the per-family ``ByFamily`` speed, across-site ``+Γ``, codon models, and the
``record=`` memory dial) are named later slices; each is a pure addition.

The result is a :class:`SequencesResult` bundle mirroring the other levels:
``.alignments`` (the observable sequence at every **extant** tip), ``.ancestral`` (the reconstructed
sequence at every **internal** node), ``.phylograms`` (each gene tree with branch lengths in
substitutions/site — the ground-truth tree behind each alignment), ``.species_phylogram`` (the species
tree scaled the same way — the molecular clock made visible), ``.genomes`` and ``.initial_genome``
(every node's whole genome, assembled, and the one the run started with — a **nucleotide** run only),
and ``.seed``. Genuine substitution
``.events`` are the deferred opt-in ``record=`` slice, not the default spine (a substitution log is not
compact the way the speciation / D-T-L-O logs are).
"""

from __future__ import annotations

import pathlib
from collections.abc import Mapping
from dataclasses import dataclass, field

import numpy as np

from ..genomes import GenomesResult
from ..genomes.events import node_label
from ..genomes.gene_trees import GeneNode, GeneTree
from ..rates.modifiers import ByLineage, FromParent
from ..rates.rate import as_rate
from ..rates.scope import PerSite
from ..tree import Node, Tree, prune
from ..progress import progress_bar
from .evolution import evolve_gene_tree
from .substitution_models import BASES, SubstitutionModel, decode, encode, jc69

_WRITE_OUTPUTS = ("alignments", "ancestral", "founding", "phylograms", "species_phylogram",
                  "genomes", "initial_genome")

#: complement of each base, for reading a block laid down on the reverse strand
_COMPLEMENT = str.maketrans("ACGT", "TGCA")

#: The rate grammar this level wires (SPEC §5) — read by the engine gate in :func:`simulate_sequences`
#: and by the CLI's help, so a modifier is never advertised without being implemented. On the
#: substitution rate these are the two lineage clocks: ``ByLineage`` the uncorrelated ("relaxed")
#: clock, ``FromParent`` the autocorrelated clock (the rate drifts parent→child down the species tree).
WIRED_MODIFIERS = (ByLineage, FromParent)


@dataclass
class SequencesResult:
    """What :func:`simulate_sequences` returns.

    - ``alignments`` — ``{family: {g<copy>: sequence}}``: the observable gene alignment, one entry per
      **extant** gene-tree tip, keyed by its (unique, per-segment) gene id — the same labels as the
      gene tree's / phylogram's Newick leaves. Empty for a family with no surviving copy.
    - ``ancestral`` — ``{family: {g<copy>: sequence}}``: the true sequence at every node that is **not**
      an extant tip — internal nodes (the family's root gene included) and the dead tips, where a copy
      was lost or its species went extinct. With ``alignments`` it accounts for every node of the tree
      exactly once, so every label in the complete phylogram names a sequence.
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
    - ``genomes`` — ``{lineage: {chromosome id: sequence}}``: **every** node's assembled genome, its
      blocks concatenated in physical order (reverse-complemented where the genome carries them
      inverted) — extant tips, ancestors and the lineages that went extinct alike. The same coverage
      as the genome level's own ``genomes``, which is keyed the same way: the observed ones are the
      extant tips, ``{node_label(n.id) for n in complete_tree.extant()}``. Only a **nucleotide** genome
      run has any — an unordered or ordered run has gene families, not coordinates, so there is no
      genome to lay out and this is empty.
    - ``initial_genome`` — ``{chromosome id: sequence}``: the genome the run **started** with, at the
      root lineage's origination. Not in ``genomes``, because it belongs to no node: the root branch is
      real simulated time, so the root *node*'s genome is this one plus whatever happened along the
      stem. It stands to ``genomes`` as ``founding`` stands to ``ancestral``.
    - ``seed`` — the run's seed.
    - ``unit`` — what the integer key of ``alignments`` / ``ancestral`` / ``founding`` / ``phylograms``
      **names**: ``"family"`` (a gene family id) on an unordered or ordered run, ``"block"`` (an index
      into the genome run's ``root_blocks``) on a nucleotide one, where every block evolves and spacer
      has no family. They are different numbering schemes over the same ints, so a gene family id is
      **not** a key here on a nucleotide run — go through
      :meth:`~zombi2.genomes.NucleotideGenomesResult.block_of`. It is also what the filenames say.
    """

    alignments: dict[int, dict[str, str]]
    ancestral: dict[int, dict[str, str]]
    founding: dict[int, str]
    phylograms: dict[int, dict[str, str | None]]
    species_phylogram: dict[str, str | None]
    seed: int | None
    # A nucleotide run's genomes are assembled lazily (see :class:`_AssembledGenomes`) so they do not
    # all sit in memory at once; the shape a caller sees is unchanged — ``{lineage: {chromosome: seq}}``.
    genomes: "Mapping[str, dict[int, str]]" = field(default_factory=dict)
    initial_genome: dict[int, str] = field(default_factory=dict)
    unit: str = "family"

    @property
    def _stem(self) -> str:
        """The filename stem for a per-unit output, so a file never claims to be a family when it is
        a block: ``fam<n>.fasta`` against ``block<n>.fasta``."""
        return {"family": "fam", "block": "block"}[self.unit]

    def write(self, directory,
              outputs=("alignments", "phylograms", "species_phylogram", "genomes",
                       "initial_genome")) -> None:
        """Write chosen ``outputs`` to ``directory`` (created if needed). ``<u>`` below is
        ``fam<family>`` on an unordered or ordered run and ``block<index>`` on a nucleotide one — the
        integer keys mean different things, so the files say which (see :attr:`unit`):

        - ``"alignments"`` → ``<u>.fasta`` (skipped for empty families).
        - ``"ancestral"`` → ``sequences_ancestral_<u>.fasta``.
        - ``"founding"`` → ``sequences_founding.fasta``, one record ``<u>`` apiece: the sequence each
          family originated with, before its stem.
        - ``"phylograms"`` → ``phylogram_<u>_{complete,extant}.nwk`` (subs/site).
        - ``"species_phylogram"`` → ``clock_species_tree_{complete,extant}.nwk``: the species tree
          with its branches in substitutions/site — the molecular clock made visible.
        - ``"genomes"`` → ``genome_<lineage>.fasta``, one file per node — extant, extinct and
          ancestral alike — with one record per chromosome. Nucleotide runs only; nothing is written
          otherwise. The big one: a real genome times every node in the tree.
        - ``"initial_genome"`` → ``genome_initial.fasta``, the genome the run started with.
        """
        unknown = [o for o in outputs if o not in _WRITE_OUTPUTS]
        if unknown:
            raise ValueError(f"unknown write outputs {unknown}; choose from {list(_WRITE_OUTPUTS)}")
        d = pathlib.Path(directory)
        d.mkdir(parents=True, exist_ok=True)
        u = self._stem
        if "alignments" in outputs:
            for fam, aln in self.alignments.items():
                if aln:
                    (d / f"{u}{fam}.fasta").write_text(_fasta(aln))
        if "ancestral" in outputs:
            for fam, anc in self.ancestral.items():
                if anc:
                    (d / f"sequences_ancestral_{u}{fam}.fasta").write_text(_fasta(anc))
        if "founding" in outputs and self.founding:
            (d / "sequences_founding.fasta").write_text(
                _fasta({f"{u}{fam}": seq for fam, seq in sorted(self.founding.items())}))
        if "phylograms" in outputs:
            for fam, ph in self.phylograms.items():
                (d / f"phylogram_{u}{fam}_complete.nwk").write_text(ph["complete"] + "\n")
                if ph["extant"] is not None:
                    (d / f"phylogram_{u}{fam}_extant.nwk").write_text(ph["extant"] + "\n")
        if "species_phylogram" in outputs:
            sp = self.species_phylogram
            (d / "clock_species_tree_complete.nwk").write_text(sp["complete"] + "\n")
            if sp["extant"] is not None:
                (d / "clock_species_tree_extant.nwk").write_text(sp["extant"] + "\n")
        # every genome is written the same way and named by whose it is — a node label, or "initial"
        for token, genomes in (("genomes", self.genomes),
                               ("initial_genome",
                                {"initial": self.initial_genome} if self.initial_genome else {})):
            if token in outputs:
                for lineage, chroms in genomes.items():
                    (d / f"genome_{lineage}.fasta").write_text(
                        _fasta({f"{lineage}_chr{cid}": seq for cid, seq in chroms.items()}))


class _AssembledGenomes(Mapping):
    """Every node's genome, assembled **on demand** rather than all held at once.

    A nucleotide run reconstructs a whole genome for every node of the tree — hundreds of megabases
    across a real genome times a real tree. Materialising them all at once would roughly double the
    run's peak memory, on top of the per-block ``alignments``/``ancestral`` where the very same letters
    already live. So this keeps only the cheap **layout** per node —
    ``{chromosome id: [(block, gene, strand), …]}`` from
    :meth:`~zombi2.genomes.NucleotideGenomesResult.assembly` — and concatenates a node's blocks into
    its genome string only when that node is asked for. Iterating (as :meth:`SequencesResult.write`
    does) then builds one node's genome, writes it, and lets it go before the next, so the assembled
    genomes never all coexist.

    The genome level says *what* to concatenate; this puts the letters in, reading each block from the
    ``alignments`` of an extant tip or the ``ancestral`` set of every other node — reverse-complemented
    where the genome carries it inverted. It is a read-only mapping of exactly the documented shape
    ``{lineage: {chromosome id: sequence}}``: indexing, ``.items()``, ``len`` and ``in`` all behave as
    a dict does — only *when* each string is built has changed. Get the order or the strand wrong and
    the genome still looks like a genome, which is why the tests check it nucleotide by nucleotide
    against the run's own trace-back."""

    __slots__ = ("_layouts", "_alignments", "_ancestral", "_extant")

    def __init__(self, layouts: dict[str, dict[int, list]], alignments, ancestral,
                 extant_labels: set[str]) -> None:
        self._layouts = layouts                 # {label: {cid: [(block, gene, strand), …]}}
        self._alignments = alignments
        self._ancestral = ancestral
        self._extant = extant_labels            # labels whose blocks read from `alignments`

    def __getitem__(self, label: str) -> dict[int, str]:
        pieces_by_cid = self._layouts[label]    # KeyError on an unknown label, exactly like a dict
        src = self._alignments if label in self._extant else self._ancestral
        chroms: dict[int, str] = {}
        for cid, pieces in pieces_by_cid.items():
            parts = []
            for (block, gene, strand) in pieces:
                seq = src[block][f"g{gene}"]
                parts.append(seq if strand == 1 else seq.translate(_COMPLEMENT)[::-1])
            chroms[cid] = "".join(parts)
        return chroms

    def __iter__(self):
        return iter(self._layouts)

    def __len__(self) -> int:
        return len(self._layouts)


def _fasta(records: dict[str, str], width: int = 70) -> str:
    """Serialise ``{name: sequence}`` to FASTA text (sequences wrapped at ``width`` columns)."""
    lines: list[str] = []
    for name, seq in records.items():
        lines.append(f">{name}")
        lines.extend(seq[i:i + width] for i in range(0, len(seq), width))
    return "\n".join(lines) + "\n"


def _split(gene_tree, states_by_id: dict[int, np.ndarray],
           model: SubstitutionModel) -> tuple[dict[str, str], dict[str, str]]:
    """Label one family's evolved nodes by their **gene id** and split them into the **observable**
    half — the extant tips — and everything else. Gene ids are per-segment (each node has a unique
    ``copy``), so ``g<copy>`` uniquely names every node and matches the gene tree's and phylogram's
    Newick labels, pairing the sequences with their tree.

    Everything else is internal nodes *and* the dead tips: a copy that a loss ended, and one whose
    species went extinct. Both are nodes of the tree with a sequence at them, so leaving them out
    would give a phylogram whose tips name sequences that exist nowhere — and would make an extinct
    lineage's genome unreconstructable."""
    alignment: dict[str, str] = {}
    ancestral: dict[str, str] = {}
    stack = [gene_tree.complete]
    while stack:
        node = stack.pop()
        seq = decode(states_by_id[id(node)], model.alphabet)
        observable = node.is_leaf and node.kind == "extant"
        (alignment if observable else ancestral)[f"g{node.copy}"] = seq
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


def _preorder(tree) -> list[int]:
    """Species-tree node ids, parent before child — the order the autocorrelated clock descends."""
    order: list[int] = []
    stack = [tree.root]
    while stack:
        i = stack.pop()
        order.append(i)
        kids = tree.nodes[i].children
        if kids is not None:
            stack.extend(kids)
    return order


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


def simulate_sequences(genomes, *, model: SubstitutionModel, length: int | None = None,
                       intergene_model: SubstitutionModel | None = None, intergene_speed=3.0,
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
    The founding sequence of each family is drawn from the model's stationary frequencies. Deterministic
    given ``seed``.

    ``substitution`` may carry a **lineage clock** — one factor per species branch, shared across
    families, computed once before evolving, rescaling each gene-tree branch by the clock of the species
    branch it sits on: ``1.0 * mod.ByLineage(spread=)`` is the uncorrelated clock (each branch drawn
    i.i.d.), and ``1.0 * mod.FromParent(spread=)`` is the autocorrelated clock (the factor drifts
    parent→child down the species tree). Any other modifier (the ``Markov`` clock, the ``ByFamily``
    per-family speed, ``+Γ``) or a non-``PerSite`` scope is a later slice and raises.

    On a **nucleotide** genome run every root block is evolved — spacer as well as genes — each at its
    own length in bp, so ``length`` does not apply and is rejected. ``model`` evolves the genes and
    ``intergene_model`` (default ``jc69``) the spacer, at ``intergene_speed`` times the rate (default
    ``3.0``). Because the whole genome is covered, the run also **puts the genomes back together**:
    ``.genomes`` holds every node's chromosomes, blocks concatenated in physical order — the complete
    tree, reconstructed — and ``.initial_genome`` the one the run started with.

    The result carries the **phylograms** the sequences were drawn along — each gene tree and the
    species tree, with branch lengths converted from time to substitutions/site by the same
    ``base × clock × Δt``.
    """
    from ..genomes import NucleotideGenomesResult

    nucleotide = isinstance(genomes, NucleotideGenomesResult)
    if not nucleotide and not isinstance(genomes, GenomesResult):
        raise TypeError(
            f"the sequence level runs on a genome run, got {type(genomes).__name__} — pass the "
            "GenomesResult that genomes.simulate_genomes_unordered(...) returned, or the "
            "NucleotideGenomesResult from simulate_genomes_nucleotide(...): the whole run, not its "
            ".gene_trees. A sequence lives inside a gene, but its clock rides the *species* branch "
            "that gene sits on — one draw per lineage, shared by every family — so the run needs "
            "the species tree too."
        )
    species_tree = genomes.complete_tree
    if not isinstance(model, SubstitutionModel):
        raise TypeError(f"model must be a SubstitutionModel (e.g. hky85(kappa=2.0)), got {model!r}")
    if intergene_model is not None and not isinstance(intergene_model, SubstitutionModel):
        raise TypeError(f"intergene_model must be a SubstitutionModel, got {intergene_model!r}")

    if nucleotide:
        # Every recovered root block evolves — spacer as well as genes — so the run reconstructs the
        # whole genome rather than the declared loci. Each block brings its own length in bp, which
        # is why a single `length` would contradict the coordinates the genome recorded.
        if length is not None:
            raise ValueError(
                "length does not apply to a nucleotide genome run: every block carries its own "
                "length in bp, so one number here would contradict the coordinates the genomes run "
                "wrote. Drop it — the genome sets the lengths.")
        for name, m in (("model", model), ("intergene_model", intergene_model)):
            if m is not None and m.alphabet != BASES:
                raise ValueError(
                    f"{name}={m.name} is a protein model, but a nucleotide genome is measured in base "
                    "pairs and its blocks are read on either strand — amino acids have no complement "
                    "to read back. Use a nucleotide model (jc69 / k80 / hky85 / gtr).")
        if intergene_model is None:
            intergene_model = jc69()          # flat and parameterless: the null for unconstrained DNA
        if isinstance(intergene_speed, bool) or not isinstance(intergene_speed, (int, float)) \
                or intergene_speed <= 0:
            raise ValueError(f"intergene_speed must be a positive number, got {intergene_speed!r}")
        gene_trees = genomes.block_trees
        blocks = genomes.root_blocks
        genic = {span: fam for fam, span in genomes.gene_spans.items()}
        # per block: its length, whether it is genic, the model it evolves under and its speed
        per_block = {}
        for i, (src, a, b) in enumerate(blocks):
            is_gene = (src, a, b) in genic
            per_block[i] = (b - a, model if is_gene else intergene_model,
                            1.0 if is_gene else float(intergene_speed))
        # Founded from a real FASTA: a block's founding sequence is the supplied DNA at its own root
        # coordinates, encoded to states, rather than a stationary draw. A de-novo source is not in
        # initial_sequence (it arose mid-run), so its blocks still draw from the model. `None` per block
        # ⇒ draw, exactly as before, so a run without one is unchanged.
        founding_seed: dict[int, "np.ndarray | None"] = {}
        for i, (src, a, b) in enumerate(blocks):
            root = genomes.initial_sequence.get(src)
            if root is None:
                founding_seed[i] = None
                continue
            f_model = model if (src, a, b) in genic else intergene_model
            if f_model.alphabet != BASES:
                raise ValueError(
                    f"the run was founded from a FASTA (DNA), but {f_model.name} is a protein model — "
                    "a nucleotide sequence cannot found an amino-acid alignment")
            founding_seed[i] = encode(root[a:b], f_model.alphabet)
    else:
        gene_trees = genomes.gene_trees
        if length is None:
            raise ValueError("length is required: the number of sites each family evolves")
        if isinstance(length, bool) or not isinstance(length, int) or length < 1:
            raise ValueError(f"length must be a positive integer, got {length!r}")
        if intergene_model is not None:
            raise ValueError(
                "intergene_model applies to a nucleotide genome run, where blocks are genes or "
                "spacer. An unordered or ordered run has gene families only, so there is nothing "
                "for a second model to evolve.")
        per_block = None
    rate = as_rate(substitution, default_scope=PerSite)
    if not isinstance(rate.scope, PerSite):
        raise ValueError(
            f"substitution has a {type(rate.scope).__name__} scope, but the sequence engine wires only "
            f"PerSite (the default) this slice — drop the scope wrapper or use PerSite(...)."
        )
    clock_mod = None
    if rate.modifiers:
        if len(rate.modifiers) == 1 and isinstance(rate.modifiers[0], (ByLineage, FromParent)):
            clock_mod = rate.modifiers[0]
        else:
            offenders = ", ".join(sorted({type(m).__name__ for m in rate.modifiers
                                          if not isinstance(m, (ByLineage, FromParent))})
                                  or ["a second clock"])
            raise ValueError(
                f"substitution carries {offenders}, but this slice wires a single lineage clock — one "
                "ByLineage (uncorrelated) or one FromParent (autocorrelated). The Markov clock, the "
                "ByFamily per-family speed, and +Γ across-site heterogeneity are later slices."
            )
    rate_base = rate.base

    rng = np.random.default_rng(seed)
    # the lineage clock: one factor per species branch, computed once here and shared by every family
    # (so a hot species runs hot for all its genes). None ⇒ the strict clock (factor 1). ByLineage draws
    # each branch i.i.d.; FromParent drifts the factor parent→child down the species tree (autocorrelated).
    clock = None
    if isinstance(clock_mod, FromParent):
        clock = {}
        for i in _preorder(species_tree):                       # parent before child
            p = species_tree.nodes[i].parent
            clock[i] = clock_mod.initial() if p is None else clock_mod.descend(clock[p], rng)
    elif clock_mod is not None:
        clock = {sid: clock_mod.draw(rng) for sid in _all_species(gene_trees)}
    alignments: dict[int, dict[str, str]] = {}
    ancestral: dict[int, dict[str, str]] = {}
    founding: dict[int, str] = {}
    phylograms: dict[int, dict[str, str | None]] = {}
    # One transition-CDF cache per model, shared across every block that model evolves. Branch lengths
    # recur across blocks (a block passing straight through a species branch reuses its length), so a
    # run-wide cache builds a few hundred matrices where a per-block cache rebuilt tens of thousands.
    # Keyed by model identity — genes and spacer are different models and must not share a cache.
    cdf_caches: dict[int, dict[float, np.ndarray]] = {}
    bar = progress_bar(len(gene_trees), "sequences", unit="family", enabled=progress)
    for family in sorted(gene_trees):  # sorted for reproducibility given the seed
        bar.update()
        gt = gene_trees[family]
        if per_block is None:
            f_len, f_model, f_rate = length, model, rate_base
        else:                       # a nucleotide block: its own length, and spacer runs faster
            f_len, f_model, speed = per_block[family]
            f_rate = rate_base * speed
        seed_states = None if per_block is None else founding_seed[family]
        cache = cdf_caches.setdefault(id(f_model), {})
        states, founding_states = evolve_gene_tree(gt.complete, f_model, f_len, f_rate, clock, rng,
                                                   gt.origination, founding=seed_states,
                                                   cdf_cache=cache)
        alignments[family], ancestral[family] = _split(gt, states, f_model)
        founding[family] = decode(founding_states, f_model.alphabet)
        scaled = _scaled_gene_tree(gt, f_rate, clock)  # branch lengths in subs/site
        ext = scaled.extant
        phylograms[family] = {"complete": _gene_newick(scaled.complete),
                              "extant": _gene_newick(ext) if ext is not None else None}

    bar.close()

    sp_scaled = _scaled_species_tree(species_tree, rate_base, clock)   # the clock made visible
    sp_extant = prune(sp_scaled, keep="extant")
    species_phylogram = {"complete": sp_scaled.to_newick(),
                         "extant": sp_extant.to_newick() if sp_extant is not None else None}
    # A nucleotide run evolved every block, so **every** node's genome can be put back together —
    # one map, as at the genome level. Which sequences each node reads is the split the level already
    # makes: an extant tip's genes are tips of their block trees, everything else's are not. The
    # concatenation is deferred to read-time (see :class:`_AssembledGenomes`): only the cheap per-node
    # layout is captured now, so hundreds of megabases of genome do not all sit in memory beside the
    # per-block sequences they are built from.
    assembled: "Mapping[str, dict[int, str]]" = {}
    initial_genome: dict[int, str] = {}
    if nucleotide:
        # Capture the layouts in the same order the eager build used — extant nodes (read from
        # `alignments`) sorted first, then the rest (read from `ancestral`) — so the map iterates,
        # and `write` emits its files, in exactly the previous order.
        extant_ids = sorted(n.id for n in species_tree.extant())
        extant_id_set = set(extant_ids)
        extant_labels = {node_label(i) for i in extant_ids}
        ordered_ids = extant_ids + [i for i in sorted(species_tree.nodes) if i not in extant_id_set]
        layouts = {node_label(i): genomes.assembly(i) for i in ordered_ids}
        assembled = _AssembledGenomes(layouts, alignments, ancestral, extant_labels)
        # The genome the run started with. Its blocks were all laid down at the start, so each one's
        # sequence there is its `founding` draw — the state the stem leads *from*. It is not a node,
        # so it is in neither map above; the same reason `founding` is not in `ancestral`.
        for cid, pieces in genomes.initial_assembly().items():
            initial_genome[cid] = "".join(
                founding[block] if strand == 1 else founding[block].translate(_COMPLEMENT)[::-1]
                for (block, strand) in pieces)

    return SequencesResult(alignments, ancestral, founding, phylograms, species_phylogram, seed,
                           assembled, initial_genome, "block" if nucleotide else "family")


# The substitution-model menu is reached through its own module — the one canonical path,
# like `zombi2.rates.scope` / `zombi2.rates.modifiers` — never re-exported here:
#     from zombi2.sequences import substitution_models as sm;  sm.hky85(2.0)
__all__ = ["simulate_sequences", "SequencesResult"]
