"""Sequences — level 3: a sequence evolving inside a gene, along its gene tree.

A sequence lives **inside a gene**, so it sees the species tree only through its gene tree
(``SPEC §1``): `simulate_sequences()` takes a **genome run** (a
`FamilyGenomesResult`) and evolves one sequence down each family's *complete* gene
tree under a substitution **model** (the menu — nucleotide ``jc69`` · ``k80`` · ``hky85`` · ``gtr``,
or protein ``poisson`` · ``jtt`` · ``dayhoff`` · ``wag`` · ``lg``; `substitution_models`) and a
substitution **rate** (``scope(base) × modifiers``; ``SPEC §5``). Sequences are a **target** here — a
trait can drive the substitution rate (``SPEC §3``, Traits–Sequences, conditioned) — and, through
``result.gc()``, a **driver** as well: a finished run's GC content per lineage drives what comes after
it, a trait on the same species tree or a further sequence run (`composition`). Never the genome,
whose gene trees the sequences were grown along; the pair cannot be joined either.

The whole genome run is required, not just its gene trees, because a level below reads the level
above: the **species tree** is what the lineage clock rides (one rate per species branch, shared by
every family passing through it — ``SPEC §5``) and what the ``species_phylogram`` is drawn on. Bare
gene trees would run, but silently without either, so they are rejected.

``substitution`` is a per-site rate (a bare number, default ``1.0``: a gene-tree branch of ``Δt`` time
gets ``substitution · Δt`` substitutions/site — the **strict clock**), optionally times a **lineage
clock**: ``substitution = PerSite(1.0).varying_among('lineages', LogNormal(0.0, σ))`` is the
uncorrelated ("relaxed") clock, one
i.i.d. rate multiplier drawn per **species lineage** and shared by every gene passing through it, and
``substitution = PerSite(1.0).varying_among('lineages', Drift(LogNormal(0.0, σ)))`` is the
**autocorrelated** clock, where the rate drifts
parent→child down the species tree so close relatives run at similar rates (``SPEC §5``). It may also
carry a ``scaled_by(trait, {...})``, which reads a **trait grown first** and lets a lineage's state
set how fast its sequences evolve; a clock and a driver compose (modifiers multiply), and a driver that
switches mid-branch is **integrated** across the switch rather than sampled once for the branch
(`clock`). Any other modifier — ``Markov`` hops, a draw among families — raises.

Rate variation **across sites** is not a modifier and does not go in ``substitution``: it belongs to
the model, where the field puts it. ``model=hky85(2.0).across_sites(gamma_shape=0.5, invariant=0.1)``
is ``HKY85+I+G4``, and the classes are normalised to mean 1, so a branch length stays the mean
substitutions per site (`substitution_models`).

A family's sites may also be split into **partitions**, each under its own model —
``partitions=[(hky85(kappa=2.0), 600), (jc69(), 400)]`` in place of ``model=`` and ``length=`` — on a
family or ordered run. They share one alphabet and one substitution rate, so the family keeps its one
phylogram. Experimental (``SPEC §9``): Python API, no CLI flag yet.

The result is a `SequencesResult` bundle mirroring the other levels:
``.alignments`` (the observable sequence at every **extant** tip), ``.ancestral`` (the reconstructed
sequence at every **internal** node), ``.phylograms`` (each gene tree with branch lengths in
substitutions/site — the ground-truth tree behind each alignment), ``.species_phylogram`` (the species
tree scaled the same way — the molecular clock made visible), ``.genomes`` / ``.node_genomes`` and ``.initial_genome``
(every node's whole genome, assembled, and the one the run started with — a **nucleotide** run only),
and ``.seed``. There is no per-substitution event log: it would not be
compact the way the speciation and D/T/L/O logs are.
"""

from __future__ import annotations

import math
import os
import pathlib
from collections.abc import Mapping
from dataclasses import dataclass, field

import numpy as np

from ..genomes import FamilyGenomesResult
from ..genomes.events import gene_label
from ..genomes.gene_trees import GeneNode, GeneTree
from ..params.conditioned import check_mapping_fires, driven_mods, names_a_live_level, resolve_driver
from ..rng import resolve_seed, seed_sequence, stream
from ..params.mapping import Between
from ..params.evaluate import (DRAWN, INHERITED, Modifier, check_one_memory, describe, matches_declared)
from ..params.connection import Driven
from ..params.parameter import Rate, as_rate
from ..params.scope import PerSite
from ..tree import Node, Tree, prune
from .._runtime.outputs import fresh_dirs, grouped_dir
from .._runtime.progress import progress_bar
from .._runtime.summary import write_summary
from .clock import Clock, resolve_clock
from .evolution import evolve_gene_tree
from .lineage_models import Models
from .substitution_models import (BASES, SubstitutionModel, _with_frequencies, dayhoff, decode,
                                  encode, gtr, hky85, jc69, jtt, k80, lg, poisson, wag)

_WRITE_OUTPUTS = ("summary", "alignments", "ancestral", "founding", "phylograms", "species_phylogram",
                  "genomes", "initial_genome")

#: complement of each base, for reading a block laid down on the reverse strand
_COMPLEMENT = str.maketrans("ACGT", "TGCA")

#: The rate grammar this level supports (SPEC §5) — read by the engine gate in `simulate_sequences()`
#: and by the CLI's help, so a modifier is never advertised without being implemented. On the
#: substitution rate these are the two lineage clocks — a draw among lineages the uncorrelated
#: ("relaxed") clock, an inherited value the autocorrelated clock (the rate drifts parent→child down
#: the species tree) — and ``scaled_by``, the conditioned driver a trait grown first supplies (SPEC §3:
#: Traits→Sequences can be conditioned). A clock and a driver compose: modifiers multiply.
IMPLEMENTED_MODIFIERS = ((DRAWN, "lineages"), (INHERITED, "lineages"), Driven)


@dataclass
class SequencesResult:
    """What `simulate_sequences()` returns.

    - ``alignments`` — ``{family: {n<species>_g<copy>: sequence}}``: the observable gene alignment, one entry per
      **extant** gene-tree tip, keyed by its (unique, per-segment) gene id — the same labels as the
      gene tree's / phylogram's Newick leaves. Empty for a family with no surviving copy.
    - ``ancestral`` — ``{family: {n<species>_g<copy>: sequence}}``: the true sequence at every node that is **not**
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
      behind each alignment. **Every** node is labelled ``n<species>_g<copy>``, so the tips match
      the ``alignments`` keys and the internal nodes match the ``ancestral`` keys (the phylogram pairs
      one-to-one with the sequences). ``"extant"`` is ``None`` for a family with no survivor.
    - ``species_phylogram`` — ``{"complete": newick, "extant": newick | None}``: the **species tree**
      with branch lengths in substitutions/site — the molecular clock made visible (which lineages ran
      hot / cold). Always present: a run always comes from a genome run, which carries its tree.
    - ``node_genomes`` — ``{lineage: {chromosome id: sequence}}``: **every** node's assembled genome,
      its blocks concatenated in physical order (reverse-complemented where the genome carries them
      inverted) — extant tips, ancestors and the lineages that went extinct alike. The same coverage,
      and the same name, as the genome level's ``node_genomes``. ``genomes`` is the extant tips alone
      — the observed genomes — exactly as it is one level down. Only a **nucleotide** genome run has
      any: a family or ordered run has gene families, not coordinates, so there is no genome to lay
      out and both are empty.
    - ``initial_genome`` — ``{chromosome id: sequence}``: the genome the run **started** with, at the
      root lineage's origination. Not in ``genomes``, because it belongs to no node: the root branch is
      real simulated time, so the root *node*'s genome is this one plus whatever happened along the
      stem. It stands to ``genomes`` as ``founding`` stands to ``ancestral``.
    - ``seed`` — the run's seed.
    - ``unit`` — what the integer key of ``alignments`` / ``ancestral`` / ``founding`` / ``phylograms``
      **names**: ``"family"`` (a gene family id) on a family or ordered run, ``"block"`` (an index
      into the genome run's ``root_blocks``) on a nucleotide one, where every block evolves and spacer
      has no family. They are different numbering schemes over the same ints, so a gene family id is
      **not** a key here on a nucleotide run — go through
      `block_of()`. It is also what the filenames say.
    - ``alphabet`` — what the sequences are **written in**: ``"ACGT"``, or the 20 amino acids. One
      per run (every partition of a family shares it), and what `composition` checks the letters it is
      asked to count against.
    """

    #: ``{block: {record: sequence}}`` — read as a plain mapping. On a nucleotide run with
    #: insertions it is a `_SplicedAlignments` view rather than a dict, so a block's rows arrive
    #: with the runs inserted into it opened as real columns.
    alignments: "Mapping[int, dict[str, str]]"
    ancestral: "Mapping[int, dict[str, str]]"
    founding: dict[int, str]
    phylograms: dict[int, dict[str, str | None]]
    species_phylogram: dict[str, str | None]
    seed: int | None
    # A nucleotide run's genomes are assembled lazily (see `_AssembledGenomes`) so they do not
    # all sit in memory at once; the shape a caller sees is unchanged — ``{lineage: {chromosome: seq}}``.
    node_genomes: "Mapping[str, dict[int, str]]" = field(default_factory=dict)
    initial_genome: dict[int, str] = field(default_factory=dict)
    unit: str = "family"
    alphabet: str = ""
    #: the tip labels of the **extant** species — what `genomes` filters ``node_genomes`` down to.
    #: Kept because a sequences result carries its trees as Newick text, not as a `Tree`, so there is
    #: otherwise nothing here that says which of the labels is a survivor.
    extant_tips: tuple[str, ...] = ()
    #: Per host block, which inserted runs sit inside it and which record carries which — the plan
    #: `alignments` splices by. Kept because it is the *plan*, not the rows: the rows are the
    #: alignment over again, and a whole-genome run cannot afford a second copy of that.
    #: Empty unless the run came from a nucleotide genome that used ``insertion``.
    _insertions: dict = field(default_factory=dict, repr=False)
    #: The per-block rows underneath the spliced view — what a genome is assembled from, and what a
    #: test checks the splice against. Same objects `alignments` / `ancestral` read through.
    _raw_rows: dict = field(default_factory=dict, repr=False)
    _raw_ancestral: dict = field(default_factory=dict, repr=False)

    def __repr__(self) -> str:
        n = sum(len(a) for a in self.alignments.values())
        return (f"SequencesResult({n} sequences across {len(self.alignments)} {self.unit} "
                f"alignments, seed={self.seed})")

    @property
    def genomes(self) -> dict[str, dict[int, str]]:
        """The observed genomes — the assembled genome of each **extant** tip, keyed by tip name
        (``n5``): what you would have sequenced. Empty unless the run came from a **nucleotide**
        genome run, which is the only kind that has coordinates to lay out.

        `node_genomes` is the whole record — every node, ancestors and extinct lineages included.
        The pair is the genome level's `~zombi2.genomes.FamilyGenomesResult.genomes` /
        ``node_genomes``, and it reads here the way ``alignments`` reads against ``ancestral``:
        what you observe, and the truth behind it."""
        return {t: self.node_genomes[t] for t in self.extant_tips if t in self.node_genomes}

    def composition(self, letters: str):
        """The share of a lineage's sequence that is one of ``letters``, at every instant, as a
        conditioning driver (`~zombi2.sequences.composition.Composition`)::

            proteins = simulate_sequences(g, model=lg(), length=300, seed=1)
            simulate_discrete(tree, states=["mesophile", "thermophile"], start="mesophile", seed=2,
                              switch=PerLineage(0.2).scaled_by(proteins.composition("KR"),
                                                               Curve(lambda x: 40.0 ** (x - 0.1))))

        This is how an **amino-acid frequency** is asked for: one residue (``"K"``) or a set of them
        (``"KR"``, ``"AVLIMFWP"``). The letters must be in this run's `alphabet`; `gc` is the same
        driver over ``"GC"``, named because it is the one people ask for by name.

        A number, so it takes a `~zombi2.params.mapping.Curve` or a `~zombi2.params.mapping.Scalar`, and
        it drives what comes **after** a sequence — a trait, or a further sequence run — never the
        genome the gene trees came from."""
        from .composition import Composition
        if not isinstance(letters, str):
            raise TypeError(
                f"composition() takes the letters to count as a string — composition('KR'), not "
                f"{letters!r}.")
        return Composition(self, letters.upper())

    def gc(self, family: object = None):
        """This run's **GC content** as a conditioning driver: `composition` over ``"GC"``, the
        fraction of a lineage's DNA that is G or C, pooled over every family the run evolved::

            seqs = simulate_sequences(g, model=hky85(2.0), length=300, seed=1)
            simulate_continuous(tree, rate=PerLineage(1.0).scaled_by(seqs.gc(),
                                                                     Curve(lambda x: 4.0 * x)),
                                seed=2)

        Nucleotide runs only, because G and C are also glycine and cysteine: on a protein run the
        call is ambiguous rather than wrong, so it is refused and `composition` asked for instead.
        ``family`` is here only to refuse one."""
        if family is not None:
            raise ValueError(
                f"gc() is pooled over every family this run evolved, so it takes no family — got "
                f"{family!r}. One family's GC is not offered: it is undefined wherever that family is "
                f"absent, and a driver has to answer for every branch the target walks. Evolve that "
                f"family in a run of its own if its GC alone should drive the rate.")
        if set(self.alphabet) != set(BASES):
            raise ValueError(
                f"gc() is GC content, so it needs DNA; this run's alphabet is {self.alphabet!r}. A "
                f"protein run's G and C are glycine and cysteine — ask for those by name with "
                f"composition('GC') if that is what you meant.")
        return self.composition("GC")

    @property
    def _stem(self) -> str:
        """The filename stem for a per-unit output, so a file never claims to be a family when it is
        a block: ``fam<n>.fasta`` against ``block<n>.fasta``."""
        return {"family": "fam", "block": "block"}[self.unit]

    def summary(self) -> dict:
        """What this run produced, as a plain dict — the payload of ``sequences_summary.json``.

        The identity is the number worth having: it is what says whether the alignments carry signal,
        and it depends on the height of the tree the run went down, which no flag shows you. The run
        prints it and warns when it is near the floor; this is the machine-readable copy of that."""
        aligned = {k: v for k, v in self.alignments.items() if v}
        sites = sorted({len(s) for aln in aligned.values() for s in aln.values()})
        return {
            "level": "sequences",
            "seed": self.seed,
            "unit": self.unit,
            {"family": "families", "block": "blocks"}[self.unit] + "_with_sequences": len(aligned),
            "sequences": sum(len(a) for a in aligned.values()),
            "ancestral_sequences": sum(len(a) for a in self.ancestral.values()),
            # one length on a family run; a nucleotide run gives every block its own, so report the span
            "sites": {"min": sites[0], "max": sites[-1]} if sites else {"min": None, "max": None},
            "mean_pairwise_identity": mean_pairwise_identity(aligned),
            "assembled_genomes": len(self.node_genomes),
        }

    def write(self, directory, outputs=("alignments", "phylograms", "species_phylogram", "genomes",
                                        "initial_genome", "summary"), *, flat: bool = False) -> None:
        """Write chosen ``outputs`` to ``directory`` (created if needed). ``<u>`` below is
        ``fam<family>`` on a family or ordered run and ``block<index>`` on a nucleotide one — the
        integer keys mean different things, so the files say which (see `unit`):

        - ``"summary"`` → ``sequences_summary.json``: what came out (`summary`); the ancestral count
          is dropped when ``"ancestral"`` was not asked for.
        - ``"alignments"`` → ``<u>.fasta`` under ``alignments/`` (skipped for empty families).
        - ``"ancestral"`` → ``sequences_ancestral_<u>.fasta`` under ``ancestral/``.
        - ``"founding"`` → ``sequences_founding.fasta``, one record ``<u>`` apiece: the sequence each
          family originated with, before its stem.
        - ``"phylograms"`` → ``phylogram_<u>_{complete,extant}.nwk`` (subs/site) under ``phylograms/``.
        - ``"species_phylogram"`` → ``clock_species_tree_{complete,extant}.nwk``: the species tree
          with its branches in substitutions/site — the molecular clock made visible.
        - ``"genomes"`` → ``genome_<lineage>.fasta`` under ``genomes/``, one file per node — extant,
          extinct and ancestral alike — with one record per chromosome. Nucleotide runs only; nothing
          is written otherwise. The big one: a real genome times every node in the tree.
        - ``"initial_genome"`` → ``genome_initial.fasta``, in ``genomes/`` with the rest: it is a
          whole-genome FASTA like they are, and it belongs beside them.

        Everything that is one file per family or per node gets a subdirectory, or the two trees and
        the one founding FASTA would be lost among thousands; ``flat=True`` writes everything into
        ``directory`` instead. Nothing is created for an output this run has none of, so a family run
        leaves no empty ``genomes/`` behind.
        """
        unknown = [o for o in outputs if o not in _WRITE_OUTPUTS]
        if unknown:
            raise ValueError(f"unknown write outputs {unknown}; choose from {list(_WRITE_OUTPUTS)}")
        d = pathlib.Path(directory)
        d.mkdir(parents=True, exist_ok=True)
        # a run's directory describes that run: clear the per-unit directories this write is
        # about to fill, so nothing from a previous run survives inside them (see fresh_dirs)
        fresh_dirs(d, ("alignments", "ancestral", "phylograms", "genomes"), flat)
        u = self._stem
        if "alignments" in outputs and any(self.alignments.values()):
            into = grouped_dir(d, "alignments", flat)
            for fam, aln in self.alignments.items():
                if aln:
                    _write_fasta(into / f"{u}{fam}.fasta", aln)
        if "ancestral" in outputs and any(self.ancestral.values()):
            into = grouped_dir(d, "ancestral", flat)
            for fam, anc in self.ancestral.items():
                if anc:
                    _write_fasta(into / f"sequences_ancestral_{u}{fam}.fasta", anc)
        if "founding" in outputs and self.founding:
            _write_fasta(d / "sequences_founding.fasta",
                         {f"{u}{fam}": seq for fam, seq in sorted(self.founding.items())})
        if "phylograms" in outputs and self.phylograms:
            into = grouped_dir(d, "phylograms", flat)
            for fam, ph in self.phylograms.items():
                complete = ph["complete"]
                assert complete is not None        # only the extant member of a pair can be absent
                (into / f"phylogram_{u}{fam}_complete.nwk").write_text(complete + "\n", encoding="utf-8")
                if ph["extant"] is not None:
                    (into / f"phylogram_{u}{fam}_extant.nwk").write_text(ph["extant"] + "\n", encoding="utf-8")
        if "summary" in outputs:
            # The written summary describes the run *as written*, which is not quite what `summary()`
            # describes. Ancestral sequences are reconstructed in memory either way but only land on
            # disk when asked for, so a default run reported a count of ancestral sequences beside a
            # directory that had none — and someone parsing the JSON, which is the point of shipping
            # JSON, concluded the dataset held reconstructions and then could not find them. When you
            # inherit a folder you cannot tell "never written" from "lost in transfer". Dropping the
            # key says there are none here, which is true; a 0 would claim none were reconstructed,
            # which is not.
            written = self.summary()
            if "ancestral" not in outputs:
                written.pop("ancestral_sequences", None)
            write_summary(d / "sequences_summary.json", written)
        if "species_phylogram" in outputs:
            sp = self.species_phylogram
            complete = sp["complete"]
            assert complete is not None
            (d / "clock_species_tree_complete.nwk").write_text(complete + "\n", encoding="utf-8")
            if sp["extant"] is not None:
                (d / "clock_species_tree_extant.nwk").write_text(sp["extant"] + "\n", encoding="utf-8")
        # every genome is written the same way and named by whose it is — a node label, or "initial"
        for token, genomes in (("genomes", self.node_genomes),
                               ("initial_genome",
                                {"initial": self.initial_genome} if self.initial_genome else {})):
            if token in outputs and genomes:   # both land in genomes/: an assembled genome either way
                into = grouped_dir(d, "genomes", flat)
                for lineage, chroms in genomes.items():
                    _write_fasta(into / f"genome_{lineage}.fasta",
                                 {f"{lineage}_chr{cid}": seq for cid, seq in chroms.items()})


def _gap_what_is_not_carried(layouts, alignments, ancestral, root_blocks) -> None:
    """Write ``-`` into every position of a block's sequence that its own lineage does not carry.

    A block is evolved over its whole ancestral extent, because that is the coordinate space its
    tree lives in; an indel then leaves a lineage holding only part of it. Without this the alignment
    hands back bases that lineage does not have — a 120 bp row for a copy carrying 1 bp of the block —
    and those rows are what the per-family FASTA files are written from. The gapped row is the true
    alignment, and it costs nothing to say: the sub-ranges in ``layouts`` already record where the
    gaps go.

    Safe to do **in place**, though the assembly reads these same strings: a gap only ever falls
    outside a carried sub-range, so slicing one out still yields unbroken sequence.

    A record no layout mentions is left alone — a copy that died mid-branch is in no node's genome,
    and blanking it would be a change to every run, indels or not. Without indels every piece is a
    whole block, so nothing here is gapped and the output is what it always was."""
    carried: dict[tuple[int, str], list[tuple[int, int]]] = {}
    for label, by_cid in layouts.items():
        for pieces in by_cid.values():
            for (block, gene, _strand, lo, hi) in pieces:
                carried.setdefault((block, f"{label}_{gene_label(gene)}"), []).append((lo, hi))
    for (block, record), spans in carried.items():
        width = root_blocks[block][2] - root_blocks[block][1]
        spans.sort()
        if len(spans) == 1 and spans[0] == (0, width):
            continue                                    # the whole block: nothing to gap
        table = alignments if record in alignments.get(block, ()) else ancestral
        seq = table[block][record]
        out, at = [], 0
        for lo, hi in spans:
            if lo > at:
                out.append("-" * (lo - at))
            out.append(seq[lo:hi])
            at = hi
        if at < width:
            out.append("-" * (width - at))
        table[block][record] = "".join(out)


def _insertion_plan(genomes, layouts, names) -> dict:
    """``{host block: (runs, carried)}`` — where each inserted run sits inside its host block, and
    which record carries which. The plan `SequencesResult.alignments` splices by.

    ``runs`` is ``[(host offset, inserted block, width), …]`` ordered by ``(offset, source)``: the
    source id is a counter, so two runs that landed at the same offset in different lineages get a
    stable order and stay two runs rather than one. ``carried`` is ``{host record: {inserted block:
    that node's record in it}}``.

    The host of a run is read off the genome's **own layout** rather than from the insertion event:
    an inserted block sitting physically between two pieces of the same host block is inside it, and
    doing it this way gets a duplicated gene right for free, since each copy has its own pieces. On a
    host carried **reversed** the pieces come out in descending coordinate order, so the two that
    meet are the left piece's ``lo`` and the right piece's ``hi`` — the offset itself is the same
    number either way, because it is an ancestral position and inversion does not move those. The
    splice needs no orientation either: an alignment is in ancestral orientation throughout, and
    strand is applied only when a genome is assembled.

    A run between two *different* blocks belongs to neither — it landed on a block boundary, or a
    rearrangement carried it off — and is left alone."""
    from ..genomes.nucleotide import Origination

    inserted = {i for i, (src, _a, _b) in enumerate(genomes.root_blocks)
                if src in {e.source for e in genomes.events
                           if isinstance(e, Origination) and e.kind == "insertion"}}
    if not inserted:
        return {}
    width = {i: b - a for i, (_s, a, b) in enumerate(genomes.root_blocks)}
    source = {i: s for i, (s, _a, _b) in enumerate(genomes.root_blocks)}
    at_offset: dict[int, dict[int, int]] = {}                    # host -> {ins block: offset}
    carried: dict[int, dict[str, dict[int, str]]] = {}           # host -> {record: {ins: record}}
    for label, by_chrom in layouts.items():
        for pieces in by_chrom.values():
            for k in range(1, len(pieces) - 1):
                block, gene, _strand, _lo, _hi = pieces[k]
                if block not in inserted:
                    continue
                left, right = pieces[k - 1], pieces[k + 1]
                if left[0] != right[0] or left[0] == block:
                    continue                                     # no single host to belong to
                meet = (left[4], right[3]) if left[2] == 1 else (left[3], right[4])
                if meet[0] != meet[1]:
                    continue                                     # not adjacent: material between them went
                host = left[0]
                at_offset.setdefault(host, {})[block] = meet[0]
                carried.setdefault(host, {}).setdefault(
                    f"{label}_{gene_label(left[1])}", {})[block] = f"{label}_{gene_label(gene)}"
    return {host: (sorted(((off, b, width[b]) for b, off in blocks.items()),
                          key=lambda r: (r[0], source[r[1]])),
                   carried.get(host, {}))
            for host, blocks in at_offset.items()}


class _SplicedAlignments(Mapping):
    """One block's alignment **as a locus**: its own columns, with the runs inserted into it opened
    at the offsets they landed on. What ``alignments`` and ``ancestral`` hand back.

    An insertion has no position in the block it lands in — novel DNA descends from nothing, so it
    arrives on a source of its own and is a block of its own, with its own tree. Left at that, a
    gene's alignment would show every deletion as a gap and no insertion at all: correct about what
    evolved, quiet about half of what happened. So the columns are put back where they belong, and a
    lineage that lacks a run shows gaps there exactly as a deletion reads.

    A block that *is* such a run is hidden from here, because its letters are now columns of its
    host: each base appears in exactly one alignment. It keeps its entry in `phylograms` and
    `founding`, which is not an inconsistency but the honest shape — an inserted run has its own
    history, starting at the insertion, and the columns of a spliced alignment therefore do **not**
    all share one tree. That is true of any alignment with indels in it.

    Spliced on access rather than stored: the plan is small and the rows are the whole alignment
    again. The per-block rows underneath are what a genome is assembled from, and they stay as they
    are — the assembly slices a block's own coordinates, which the spliced view no longer counts in."""

    __slots__ = ("_raw", "_plan", "_width", "_hidden", "_span")

    def __init__(self, raw: dict, plan: dict, width: dict, hidden: frozenset) -> None:
        self._raw, self._plan, self._width, self._hidden = raw, plan, width, hidden
        self._span: dict[int, int] = {}

    def _columns(self, block: int) -> int:
        """How many columns ``block`` has once its runs are opened — its own width plus theirs,
        and theirs may hold runs of their own: an insertion can land inside one that arrived
        earlier. Memoised, and it terminates because a run is always younger than what it landed in."""
        if block not in self._span:
            runs = self._plan.get(block, ((), {}))[0]
            self._span[block] = self._width[block] + sum(self._columns(b) for (_o, b, _w) in runs)
        return self._span[block]

    def _rows(self, block: int) -> dict[str, str]:
        """``block``'s rows with its runs spliced in, hidden or not — what `__getitem__` guards."""
        rows = self._raw.get(block, {})
        plan = self._plan.get(block)
        if not plan:
            return rows
        runs, carried = plan
        inner = {b: self._rows(b) for (_o, b, _w) in runs}     # spliced, so nesting composes
        out: dict[str, str] = {}
        for record, seq in rows.items():
            mine, parts, at = carried.get(record, {}), [], 0
            for (offset, ins_block, _w) in runs:
                parts.append(seq[at:offset])
                theirs = mine.get(ins_block)
                got = None if theirs is None else inner[ins_block].get(theirs)
                parts.append(got if got is not None else "-" * self._columns(ins_block))
                at = offset
            parts.append(seq[at:])
            out[record] = "".join(parts)
        return out

    def __getitem__(self, block: int) -> dict[str, str]:
        if block in self._hidden or block not in self._raw:
            raise KeyError(block)
        return self._rows(block)

    def __iter__(self):
        return (b for b in self._raw if b not in self._hidden)

    def __len__(self) -> int:
        return sum(1 for _ in self)


class _AssembledGenomes(Mapping):
    """Every node's genome, assembled **on demand** rather than all held at once.

    A nucleotide run reconstructs a whole genome for every node of the tree — hundreds of megabases
    across a real genome times a real tree. Materialising them all at once would roughly double the
    run's peak memory, on top of the per-block ``alignments``/``ancestral`` where the very same letters
    already live. So this keeps only the cheap **layout** per node —
    ``{chromosome id: [(block, gene, strand), …]}`` from
    `assembly()` — and concatenates a node's blocks into
    its genome string only when that node is asked for. Iterating (as `SequencesResult.write()`
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
            for (block, gene, strand, lo, hi) in pieces:
                # the copy a node carries is a copy *of that node*, so its record is named for this
                # label — which is already `n<id>`, the first half of the record name.
                # `[lo, hi)` is the part of the block this node actually carries: the whole of it
                # unless an indel took a stretch out of the middle or opened a gap inside it.
                seq = src[block][f"{label}_{gene_label(gene)}"][lo:hi]
                parts.append(seq if strand == 1 else seq.translate(_COMPLEMENT)[::-1])
            chroms[cid] = "".join(parts)
        return chroms

    def __iter__(self):
        return iter(self._layouts)

    def __len__(self) -> int:
        return len(self._layouts)


def _identity_counts(seqs) -> tuple[int, int]:
    """``(matching positions, positions compared)`` over **every** within-family pair of ``seqs``.

    Counted a column at a time, not a pair at a time: at one column, a residue shared by ``c`` of the
    sequences contributes ``c(c-1)/2`` matching pairs, so a family of ``k`` sequences over ``n`` sites
    costs one pass per residue of the alphabet rather than one per pair. The work is linear in the
    alignment — ``k·n`` — while the pairs it accounts for are quadratic, which is what makes counting
    all of them affordable. Sequences of unequal length are compared over their shared prefix, as a
    pair-at-a-time count did."""
    if len(seqs) < 2:
        return 0, 0
    n = min(len(s) for s in seqs)
    if not n:
        return 0, 0
    k = len(seqs)
    arr = np.frombuffer("".join(s[:n] for s in seqs).encode(), dtype=np.uint8).reshape(k, n)
    # A gap is not a residue and two of them are not a match — `np.unique` walks the bytes that are
    # actually there, so without this exclusion a pair of lineages that both deleted the same stretch
    # would read as identical across it. Pairs are counted per column over the sequences that have a
    # residue there, which on an alignment with no gaps is every pair at every column, exactly as
    # this counted before indels could put one in.
    gap = ord("-")
    real = np.count_nonzero(arr != gap, axis=0).astype(np.int64)
    compared = int((real * (real - 1) // 2).sum())
    matched = 0
    for residue in np.unique(arr):
        if residue == gap:
            continue
        c = np.count_nonzero(arr == residue, axis=0).astype(np.int64)
        matched += int((c * (c - 1) // 2).sum())
    return matched, compared


def mean_pairwise_identity(alignments) -> float | None:
    """Mean identity over **every** within-family sequence pair, or ``None`` when no family holds two
    sequences to compare.

    Exhaustive, and the same quantity however the run was written: it used to be a bounded random
    sample, and the streamed path sampled differently, so one flag that only chooses how much memory
    a run uses moved a number in the report while leaving every alignment byte-identical. Nothing
    said the number was an estimate, so it read as a fingerprint two matching runs could disagree on.
    `_identity_counts` makes counting them all linear in the alignment, so there is no longer a reason
    to estimate."""
    matched = compared = 0
    for a in alignments.values():
        m, c = _identity_counts(list(a.values()))
        matched += m
        compared += c
    return matched / compared if compared else None


def _write_fasta(path, records: dict[str, str], width: int = 70) -> None:
    """Write ``{name: sequence}`` to ``path`` as FASTA (sequences wrapped at ``width`` columns),
    streaming one record at a time straight to the file so a whole-genome sequence is never first
    copied into one big string (nor a list of every wrapped line). Byte-for-byte what building the
    text and writing it produced — including the lone ``"\\n"`` an empty record set wrote."""
    with open(path, "w", encoding="utf-8") as f:
        if not records:
            f.write("\n")                       # the degenerate case "\n".join([]) + "\n" produced
            return
        for name, seq in records.items():
            f.write(f">{name}\n")
            # writelines drains the generator in C — as fast as a joined string but without ever
            # holding the whole wrapped genome (nor a list of its lines) in memory at once.
            f.writelines(f"{seq[i:i + width]}\n" for i in range(0, len(seq), width))


def _split(gene_tree, states_by_id: dict[int, np.ndarray], names,
           model: SubstitutionModel) -> tuple[dict[str, str], dict[str, str]]:
    """Label one family's evolved nodes by their **gene id** and split them into the **observable**
    half — the extant tips — and everything else. Gene ids are per-segment (each node has a unique
    ``copy``), so ``g<copy>`` uniquely names every node and matches the gene tree's and phylogram's
    Newick labels, pairing the sequences with their tree.

    Everything else is internal nodes *and* the dead tips: a copy that a loss ended, and one whose
    species went extinct. Both are nodes of the tree with a sequence at them, so leaving them out
    would give a phylogram whose tips name sequences that exist nowhere — and would make an extinct
    lineage's genome unreconstructable."""
    nodes = []
    stack = [gene_tree.complete]
    while stack:
        node = stack.pop()
        nodes.append(node)
        stack.extend(node.children)
    # Decode every node in one gather + one ASCII decode rather than one call per node: all of a
    # block's nodes share its length, so their states stack into a single (n_nodes, length) array
    # whose flat decode is the node strings back to back — sliced out below. Byte-identical to
    # decoding each node on its own, but it pays the numpy/ASCII-decode overhead once, not per node.
    stacked = np.stack([states_by_id[id(n)] for n in nodes])
    flat = decode(stacked, model.alphabet)
    length = stacked.shape[1]
    alignment: dict[str, str] = {}
    ancestral: dict[str, str] = {}
    for i, node in enumerate(nodes):
        seq = flat[i * length:(i + 1) * length]
        observable = node.is_leaf and node.kind == "extant"
        (alignment if observable else ancestral)[f"{names[node.species]}_{gene_label(node.copy)}"] = seq
    return alignment, ancestral


def _calibrate(substitution, divergence: float, tree: Tree) -> Rate:
    """Turn a wanted **divergence** into the rate that produces it, and keep whatever clock shape
    ``substitution`` was carrying.

    ``divergence`` is substitutions per site along a root-to-tip path, so the base is simply
    ``divergence / height``: the rate is per unit time, and the path is ``height`` long. It is the one
    number a user can judge — 0.2 is a readable alignment, 2.0 is noise — where a rate cannot be
    judged without knowing the height of a tree they have not measured. The same rate on a tree ten
    times taller means something completely different, which is why no default rate can be right.

    The two arguments say different things and compose: ``substitution`` says what *kind* of clock
    (strict, or relaxed by a modifier), ``divergence`` says how far it drifts. A base given alongside
    is refused rather than overridden — silently replacing a number someone typed is how a run comes
    to differ from what its command line says.

    So what ``substitution`` may be here is exactly the shapes that name **no base**: nothing at all
    (the strict clock), a bare ``Random``, which is a value rather than a parameter and says how the
    rate varies without saying how fast it is, or a scope with no number in front of it and the
    verbs chained on — ``PerSite().varying_among('lineages', LogNormal(0.0, 0.3))``. The last is the
    one the written form can carry, since a bare ``Random`` is not an expression a flag accepts.

    A **driven** rate is refused here for a modelling reason, not a coding one: ``divergence / height``
    is the base only when the modifiers average to 1 along a root-to-tip path, which is what
    mean-correcting a draw among lineages and an inherited value buys (SPEC §5) and what a driver
    deliberately does
    not promise — its factor is whatever the trait's state says. Solving as though it did would
    produce a run whose realised divergence is off by the driver's mean factor while the log claims
    the number that was asked for. Set the base yourself alongside the driver.
    """
    if not isinstance(divergence, (int, float)) or isinstance(divergence, bool):
        raise ValueError(f"divergence must be a number of substitutions per site, got {divergence!r}")
    if not (divergence > 0) or divergence == float("inf"):
        raise ValueError(f"divergence must be positive and finite, got {divergence!r}")
    # Judged on what the caller wrote, not on a coerced rate: `as_rate` fills in a default base as
    # readily as it fills in a default scope, so a rate that has been through it cannot tell "shape
    # only" from "a base I chose". A `Rate` whose `base` is None is the written shape-only form.
    if substitution is None:
        mods: tuple[Modifier, ...] = ()
    elif isinstance(substitution, Modifier):
        mods = (substitution,)
    elif isinstance(substitution, Rate) and substitution.base is None:
        mods = substitution.modifiers
    else:
        raise ValueError(
            f"substitution names a base, and divergence={divergence} would override it — the base is "
            f"what divergence solves for. Give the clock's shape alone "
            f"(substitution=PerSite().varying_among('lineages', LogNormal(0.0, …))) to calibrate a "
            f"relaxed clock, or drop divergence and set the base yourself.")
    if any(isinstance(m, Driven) for m in mods):
        raise ValueError(
            f"substitution is driven by another level, and divergence={divergence} cannot solve for "
            "its base: divergence / height is the base only when the modifiers average to 1 along a "
            "root-to-tip path, which the two lineage clocks are mean-corrected to do and a driver is "
            "not — its factor is whatever the driver's state says. The realised divergence would be "
            "off by the driver's mean factor while this claimed the number you asked for. Write the "
            "base yourself: substitution=PerSite(0.01).scaled_by(driver, {…}).")
    height = max(n.end_time for n in tree.nodes.values()) - min(n.birth_time for n in tree.nodes.values())
    if height <= 0:
        raise ValueError("divergence needs a tree with height to divide by; this one has none")
    # A written shape keeps whichever scope it was written from, so the engine's own scope check gets
    # to report it; the bare-modifier and strict forms name none, and take this level's default.
    written = substitution.scope if isinstance(substitution, Rate) else None
    return Rate(base=divergence / height, scope=written or PerSite, modifiers=mods)


def _scaled_gene_tree(gt: GeneTree, rate_base: float, clock: "Clock | None") -> GeneTree:
    """A copy of the gene tree whose node ``time`` holds the cumulative **substitutions/site** from the
    family's **origination** (`Clock.branch_length` summed along the path). Feeding it to
    ``GeneTree.to_newick`` then emits a *phylogram* (branch lengths in subs/site); and because its
    prune-to-extant merges branches by that same cumulative measure, a suppressed branch spanning
    several species branches gets the **sum** of its pieces for free — the exact trick the chronogram
    uses with time.

    Counting from origination rather than from the root is what gives the root its own branch: the
    founding gene evolves across the stem, so the scaled root sits one stem's worth in, not at zero.

    The lengths come from the same `Clock` the sampler reads, and by branch *endpoints* rather than by
    a per-branch factor, which is what keeps a **driven** rate consistent: a trait that switches
    mid-branch makes the driver's contribution an integral over the branch, and a phylogram that
    scaled by a single sample would not be the tree its own alignment was drawn along."""
    root = gt.complete
    if clock is None:
        stem = rate_base * (root.time - gt.origination)
    else:
        stem = clock.branch_length(rate_base, root.species, gt.origination, root.time)
    scaled_root = GeneNode(root.kind, root.species, stem, root.copy)
    stack = [(root, scaled_root)]
    while stack:
        onode, snode = stack.pop()
        for ochild in onode.children:
            if clock is None:
                blen = rate_base * (ochild.time - onode.time)
            else:
                blen = clock.branch_length(rate_base, ochild.species, onode.time, ochild.time)
            schild = GeneNode(ochild.kind, ochild.species, snode.time + blen, ochild.copy)
            snode.children.append(schild)
            stack.append((ochild, schild))
    return GeneTree(gt.family, scaled_root, 0.0)   # origination is the zero of the scaled measure


def _gene_newick(root: GeneNode, names) -> str:
    """Newick of a (scaled) gene tree labelling **every** node — leaf and internal — by its gene id
    ``n<species>_g<copy>``, so the tips match the ``alignments`` keys and the internal nodes match the
    ``ancestral`` keys (both keyed ``n<species>_g<copy>``): the phylogram pairs one-to-one with the sequences.
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
        bl = f":{node.time - parent_time:.7g}"
        tag = f"{names[node.species]}_{gene_label(node.copy)}"
        s = f"{tag}{bl}" if node.is_leaf else f"({','.join(parts)}){tag}{bl}"
        stack.pop()
        if stack:
            stack[-1][3].append(s)
        else:
            result = s
    return result + ";"


def _scaled_species_tree(tree: Tree, rate_base: float, clock: "Clock | None") -> Tree:
    """A copy of the species tree whose branch lengths are **substitutions/site**
    (`Clock.branch_length` over each whole branch). Node times become the cumulative subs/site from
    the root, so ``Tree.to_newick`` / ``prune`` emit and merge the phylogram exactly as they do a
    dated tree.

    This is where the clock is made visible, so it must show *all* of it: a branch whose driver
    switched partway along gets the integral across the switch, the same number the gene phylograms
    and the sampler used."""
    scaled: dict[int, Node] = {}
    scaled_end: dict[int, float] = {}
    order: list[int] = []
    stack = [tree.root]
    while stack:  # pre-order: a parent is visited before its children
        i = stack.pop()
        order.append(i)
        stack.extend(tree.nodes[i].children)
    for i in order:
        nd = tree.nodes[i]
        blen = (rate_base * (nd.end_time - nd.birth_time) if clock is None
                else clock.branch_length(rate_base, i, nd.birth_time, nd.end_time))
        start = 0.0 if nd.parent is None else scaled_end[nd.parent]
        scaled_end[i] = start + blen
        scaled[i] = Node(i, nd.parent, start, start + blen, nd.children, nd.fate)
    return Tree(scaled, tree.root)


@dataclass(frozen=True)
class StreamedSequences:
    """A sequence run written **straight to disk**, family by family — what ``stream_to=`` returns.

    Thin on purpose: the outputs *are* the files, so this carries where they are and how big the run
    was, not the run itself. The sequence level is where a run's memory actually goes — every family's
    alignment plus its ancestral sequences, all live at once — so this is the handle for the size at
    which a `SequencesResult` would not fit."""

    directory: str
    seed: int | None
    n_families: int
    n_sequences: int
    outputs: tuple
    #: mean pairwise identity over every within-family pair — what the CLI reports and warns on.
    #: ``None`` when no family held two sequences to compare.
    identity: "float | None" = None
    #: the site count every alignment carries (one length on a family run)
    sites: "int | None" = None
    #: how many ancestral sequences went by — the sink sees them whether or not they were written
    n_ancestral: int = 0
    unit: str = "family"

    def __repr__(self) -> str:
        return (f"StreamedSequences({self.n_sequences} sequences across {self.n_families} "
                f"{self.unit} alignments, streamed to {self.directory!r}, seed={self.seed})")

    def summary(self) -> dict:
        """The same payload `SequencesResult.summary()` builds, from what the sink counted as it wrote.

        A streamed run has to produce the same files as an in-memory one at the same seed — that is
        the whole contract — so it produces this one too. Two fields it cannot know: the ancestral
        count and the assembled genomes, neither of which a streamed family run writes."""
        return {
            "level": "sequences",
            "seed": self.seed,
            "unit": self.unit,
            {"family": "families", "block": "blocks"}[self.unit] + "_with_sequences":
                self.n_families,
            "sequences": self.n_sequences,
            "ancestral_sequences": self.n_ancestral,
            "sites": {"min": self.sites, "max": self.sites},
            "mean_pairwise_identity": self.identity,
            "assembled_genomes": 0,
        }


class _Sink:
    """Writes one family's sequences the moment they exist, then forgets them.

    Same files, same names and same contents as `SequencesResult.write` — it has to be, or a streamed
    run would be a different dataset from an in-memory one at the same seed. The founding sequences
    are the one output that is not per family: they share a single FASTA, so the handle stays open and
    each family appends a record, which keeps this flat in memory too."""

    def __init__(self, directory, outputs: tuple, unit: str, flat: bool) -> None:
        self.dir = pathlib.Path(directory)
        self.dir.mkdir(parents=True, exist_ok=True)
        self.outputs, self.unit, self.flat = outputs, unit, flat
        # once, here — a streamed run resolves its directories per family, so clearing on each would
        # leave only the last family written
        fresh_dirs(self.dir, ("alignments", "ancestral", "phylograms", "genomes"), flat)
        self.stem = {"family": "fam", "block": "block"}[unit]
        self.n_families = self.n_sequences = 0
        self._founding = (open(self.dir / "sequences_founding.fasta", "w", encoding="utf-8")
                          if "founding" in outputs else None)
        # Mean pairwise identity, accumulated family by family as each one is written. The CLI
        # reports it and warns when a run has saturated, which is the single most useful thing it
        # says about a sequence run — and computing it the in-memory way would need every alignment
        # at once, which is exactly what is not being kept. Counting every pair *within* a family
        # needs only that family, so the running total here is the same number the in-memory path
        # reaches: a run reports one identity whichever way it was written.
        self._matched = self._compared = 0
        self.sites: "int | None" = None
        self.n_ancestral = 0

    def _into(self, name: str) -> pathlib.Path:
        return grouped_dir(self.dir, name, self.flat)

    @property
    def identity(self) -> "float | None":
        """Mean identity over every within-family pair, or ``None`` if no family held two sequences."""
        return self._matched / self._compared if self._compared else None

    def _count_pairs(self, aln: dict) -> None:
        matched, compared = _identity_counts(list(aln.values()))
        self._matched += matched
        self._compared += compared

    def family(self, fam: int, aln: dict, anc: dict, fnd: str, phylo: dict) -> None:
        u = f"{self.stem}{fam}"
        # families with a SURVIVING copy, matching what an in-memory run reports: a family that left
        # no extant gene has an empty alignment and writes no FASTA, so counting it here would make
        # a streamed run claim more families than the same run in memory
        self.n_families += 1 if aln else 0
        self.n_sequences += len(aln)
        self.n_ancestral += len(anc)
        if self.sites is None and aln:
            self.sites = len(next(iter(aln.values())))
        self._count_pairs(aln)
        if "alignments" in self.outputs and aln:
            _write_fasta(self._into("alignments") / f"{u}.fasta", aln)
        if "ancestral" in self.outputs and anc:
            _write_fasta(self._into("ancestral") / f"sequences_ancestral_{u}.fasta", anc)
        if self._founding is not None:
            self._founding.write(f">{u}\n")
            for i in range(0, len(fnd), 70):
                self._founding.write(fnd[i:i + 70] + "\n")
        if "phylograms" in self.outputs:
            into = self._into("phylograms")
            (into / f"phylogram_{u}_complete.nwk").write_text(phylo["complete"] + "\n",
                                                              encoding="utf-8")
            if phylo["extant"] is not None:
                (into / f"phylogram_{u}_extant.nwk").write_text(phylo["extant"] + "\n",
                                                                encoding="utf-8")

    def finish(self, species_phylogram: dict) -> None:
        """The run-sized outputs, once every family has gone by."""
        if self._founding is not None:
            self._founding.close()
        if "species_phylogram" in self.outputs:
            (self.dir / "clock_species_tree_complete.nwk").write_text(
                species_phylogram["complete"] + "\n", encoding="utf-8")
            if species_phylogram["extant"] is not None:
                (self.dir / "clock_species_tree_extant.nwk").write_text(
                    species_phylogram["extant"] + "\n", encoding="utf-8")


#: what a streamed run writes when ``outputs`` is not given — the same set `SequencesResult.write`
#: defaults to, minus the two a family run never has anyway (a genome is a nucleotide-run output, and
#: nucleotide runs cannot stream: assembling a genome needs every block at once).
_DEFAULT_STREAM_OUTPUTS = ("alignments", "phylograms", "species_phylogram", "summary")


def _resolve_partitions(model, partitions,
                        length) -> "tuple[tuple[SubstitutionModel | Models, int], ...]":
    """The **site blocks** one family evolves, as ``((model, sites), …)``.

    One entry for the ordinary ``model=`` + ``length=`` run, several for a partitioned one — so the
    engine below has a single shape to walk and the un-partitioned case is literally the
    one-partition case, not a branch beside it.

    Everything else this owns is a refusal. Giving ``model=`` or ``length=`` *alongside*
    ``partitions`` is a second answer to a question the partitions already answer, and two answers
    are worse than none. Mixing alphabets across partitions is the one refusal that is a statement
    about the **model** rather than about the call: it is meaningless in ``SPEC §5``'s sense, not
    unimplemented, because the partitions concatenate into one sequence per gene copy and there is
    no string that is half DNA and half protein.
    """
    if partitions is None:
        if length is None:
            raise ValueError("length is required: the number of sites each family evolves")
        if isinstance(length, bool) or not isinstance(length, int) or length < 1:
            raise ValueError(f"length must be a positive integer, got {length!r}")
        return ((model, length),)

    if model is not None:
        raise ValueError(
            f"model={model.name if isinstance(model, (SubstitutionModel, Models)) else model!r} was given "
            "alongside partitions, and each partition already carries its own model — so this would "
            "be a second answer to the same question, and nothing here could say which one a site "
            "should follow. Drop one: partitions=[(model, sites), …] for a split sequence, or "
            "model=… with length=… for one model over all of it.")
    if length is not None:
        raise ValueError(
            f"length={length!r} was given alongside partitions, but a family's length is the sum of "
            "the partitions' site counts, so one number here would contradict them. Drop it — the "
            "partitions set the length.")

    listed = list(partitions)
    if not listed:
        raise ValueError(
            "partitions is empty, so a family would evolve no sites at all. Give at least one "
            "(model, sites) pair, or drop partitions and use model=… with length=….")
    parts: "list[tuple[SubstitutionModel | Models, int]]" = []
    for i, item in enumerate(listed):
        pair = tuple(item) if isinstance(item, (tuple, list)) else ()
        if len(pair) != 2:
            # a bare model is the likely slip, and a model's repr is a whole rate matrix — name it
            if isinstance(item, (SubstitutionModel, Models)):
                shown = f"the model {item.name} on its own"
            elif pair and isinstance(pair[0], SubstitutionModel):
                shown = f"{len(pair)} values starting with the model {pair[0].name}"
            else:
                shown = repr(item)
            raise ValueError(
                f"partition {i} is {shown}: every partition is a (model, sites) pair — which "
                "substitution model that stretch of the sequence evolves under, and how many sites "
                "it covers. For example partitions=[(hky85(kappa=2.0), 600), (jc69(), 400)].")
        m, n = pair
        if not isinstance(m, (SubstitutionModel, Models)):
            raise ValueError(
                f"partition {i}'s model is {m!r}, which is not a SubstitutionModel — take one from "
                "the menu (jc69(), hky85(kappa=2.0), lg(), …) or build your own with "
                "substitution_models.reversible(S, frequencies).")
        if isinstance(n, bool) or not isinstance(n, int) or n < 1:
            raise ValueError(
                f"partition {i} ({m.name}) covers {n!r} sites: a partition's site count must be a "
                "positive whole number of sites. A partition of zero sites is not a partition — "
                "leave it out.")
        parts.append((m, n))

    first = parts[0][0]
    if any(m.alphabet != first.alphabet for m, _ in parts):
        j, mj = next((i, m) for i, (m, _) in enumerate(parts) if m.alphabet != first.alphabet)
        raise ValueError(
            f"partition 0 is {first.name} over {first.alphabet!r} and partition {j} is {mj.name} "
            f"over {mj.alphabet!r}, but the partitions are concatenated into one sequence per "
            "gene copy and a sequence has one alphabet. That is meaningless rather than "
            "unimplemented: there is no string a half-nucleotide, half-protein gene could be "
            "written in. Use models over the same alphabet — all nucleotide, or all protein. To "
            "evolve a DNA gene and a protein gene, run the level twice.")
    return tuple(parts)


def _resolve_profiles(profiles, model, length) -> dict:
    """``{family: ((model, 1), (model, 1), …)}`` — one single-site model per position, built from a
    site profile.

    A **profile** is an ``(L, K)`` array whose row *i* is the equilibrium frequencies at position
    *i*, in the model's own ``alphabet`` order. Each row becomes a model of its own over the base
    model's exchangeabilities (`_with_frequencies`), so what changes down the sequence is which
    residues belong where, not the chemistry relating them.

    **Built once, for the whole run, on purpose.** The transition-matrix cache is keyed by
    ``id(model)``, so a model created and dropped inside the family loop could have its id reused by
    a later one and collide with a cache entry computed for a different matrix. Holding every
    per-site model for the run's lifetime is what keeps those ids stable. Do not move this into the
    loop.
    """
    if not isinstance(profiles, dict):
        raise TypeError(
            f"profiles must be a dict of {{family: array}} — one (L, K) array per family you have a "
            f"profile for, not {type(profiles).__name__}. Families you leave out evolve under `model` "
            f"as usual.")
    k = len(model.alphabet)
    out = {}
    for family, profile in profiles.items():
        arr = np.asarray(profile, dtype=float)
        if arr.ndim != 2 or arr.shape[1] != k:
            raise ValueError(
                f"the profile for family {family!r} is {arr.shape}, but it must be (L, {k}) — one row "
                f"per site and one column per state of the model's alphabet {model.alphabet!r}, in "
                f"that order.")
        if length is not None and arr.shape[0] != length:
            raise ValueError(
                f"the profile for family {family!r} has {arr.shape[0]} rows but length={length} — a "
                f"profile carries one row per site, so the two have to agree.")
        if not np.isfinite(arr).all() or (arr < 0).any():
            raise ValueError(f"the profile for family {family!r} must be finite and non-negative.")
        totals = arr.sum(axis=1)
        if (totals <= 0).any():
            raise ValueError(
                f"row {int(np.argmin(totals))} of the profile for family {family!r} sums to zero, so "
                f"that site has no state it could be in.")
        arr = arr / totals[:, None]
        if (arr <= 0).any():
            raise ValueError(
                f"the profile for family {family!r} gives a state a frequency of exactly zero, which "
                f"makes that site's rate matrix degenerate. Add a pseudocount — a real profile says a "
                f"residue is unlikely at a position, not that it is impossible.")
        out[family] = tuple((_with_frequencies(model, row, name=f"{model.name}+profile[{i}]"), 1)
                            for i, row in enumerate(arr))
    return out


def _evolve_partitions(gt, parts, rate, clock, rng, cdf_caches, names, founding=None):
    """Evolve one gene tree partition by partition and hand back the family's whole sequences:
    ``(alignment, ancestral, founding_string)``, each sequence the partitions concatenated in order.

    Every partition is evolved down the **same** tree at the **same** rate, so the family has one
    phylogram and that phylogram is exact for all of it: every model is normalised to one expected
    substitution per site per unit branch length (`substitution_models._reversible_model`) and every
    set of across-site rate classes to a mean of 1, so a branch of ``Δt`` accrues ``rate · Δt``
    substitutions per site in each partition alike. Give the partitions relative *speeds* and that
    stops being true — one phylogram would then be a weighted average of the trees the partitions
    were actually drawn along — which is why there is no per-partition rate here.

    The random draws are consumed **partition by partition**, in the order given, so a run is
    reproducible from its seed; and with a single partition the sequence of draws is exactly what an
    un-partitioned run always took, which is what keeps the default path byte-identical. That case
    also returns its dicts untouched rather than rebuilding them, so the common run does no
    concatenation work at all.

    ``founding`` (a nucleotide run's supplied DNA) is sliced per partition at the same offsets, so
    each block founds from its own stretch.
    """
    pieces = []
    at = 0
    for model, n in parts:
        # One CDF cache for the whole run: it is keyed by (model, branch length), so several models
        # share it safely and branch lengths — which recur massively across families — are still
        # computed once per model.
        if isinstance(model, Models):
            # a per-lineage model set: the family founds under the model of the species branch its
            # ORIGINATION sits on (the stem lies on the root gene's species branch, which is the same
            # branch `evolve_gene_tree` charges the stem to), and each branch below picks its own
            per_species, base = model.per_species, model.at(gt.complete.species)
        else:
            per_species, base = None, model
        states, founding_states = evolve_gene_tree(
            gt.complete, base, n, rate, clock, rng, gt.origination,
            founding=None if founding is None else founding[at:at + n],
            cdf_cache=cdf_caches, models=per_species)
        at += n
        aln, anc = _split(gt, states, names, base)
        pieces.append((aln, anc, decode(founding_states, base.alphabet)))
    if len(pieces) == 1:
        return pieces[0]
    # The node keys are the same in every partition — they come from the same tree, walked the same
    # way — so partition 0's keys are the whole set, and joining in partition order is the sequence.
    alignment = {key: "".join(p[0][key] for p in pieces) for key in pieces[0][0]}
    ancestral = {key: "".join(p[1][key] for p in pieces) for key in pieces[0][1]}
    return alignment, ancestral, "".join(p[2] for p in pieces)


def simulate_sequences(genomes, *, model: SubstitutionModel | None = None,
                       length: int | None = None, partitions=None, profiles=None,
                       intergene_model: SubstitutionModel | None = None, intergene_speed=3.0,
                       substitution=None, divergence=None, seed=None, parallel=False,
                       stream_to=None, outputs=None, flat: bool = False,
                       progress=False) -> "SequencesResult | StreamedSequences":
    """Evolve one sequence down each family's gene tree under a substitution ``model``.

    ``genomes`` is a **genome run** — the `FamilyGenomesResult` that
    ``genomes.simulate_genomes_family(...)`` returned. Its ``gene_trees`` are what the sequences
    evolve along and its ``complete_tree`` is the species tree the lineage clock rides; bare gene
    trees are rejected (they would run, but with no clock and no species phylogram — a silent
    degradation). Each family's *complete* gene tree is evolved, so the true history is complete and
    ancestral sequences exist for extinct/lost lineages too; the observable ``alignments`` are the
    extant tips.

    ``model`` is a substitution model from the menu (`substitution_models`) — nucleotide
    ``jc69`` · ``k80`` · ``hky85`` · ``gtr``, or protein ``poisson`` · ``jtt`` · ``dayhoff`` ·
    ``wag`` · ``lg``; its alphabet is what the sequences are written in (``ACGT`` or the 20 amino
    acids). ``length`` is the number of sites. ``substitution`` is the per-site substitution
    rate (default ``1.0``): a branch of ``Δt`` time accrues ``substitution · Δt`` substitutions/site.
    The founding sequence of each family is drawn from the model's stationary frequencies. Deterministic
    given ``seed``.

    ``substitution`` may carry a **lineage clock** — one factor per species branch, shared across
    families, computed once before evolving, rescaling each gene-tree branch by the clock of the species
    branch it sits on: ``PerSite(1.0).varying_among('lineages', LogNormal(0.0, σ))`` is the
    uncorrelated clock (each branch drawn
    i.i.d.), and ``PerSite(1.0).varying_among('lineages', Drift(LogNormal(0.0, σ)))`` is the
    autocorrelated clock (the factor drifts
    parent→child down the species tree).

    It may also carry a **driver** — ``PerSite(1.0).scaled_by(habitat, {"cave": 0.5, "surface": 1.0})``,
    where ``habitat`` is a trait grown first (the `~zombi2.traits.TraitsResult`, or the path to the
    ``trait_events.tsv`` it wrote). That is conditioning, not a joint run: SPEC §3 allows the pair
    Traits–Sequences to be conditioned and never joined, so naming a live level (``"trait"``) raises
    rather than starting one. A clock and a driver **compose** — the factors multiply (SPEC §5), so a
    lineage's dealt tempo and its state both count — and several drivers on one rate multiply too.
    A discrete driver switches *mid-branch*, and the branch length is the driver **integrated** across
    the branch rather than one sample of it (`clock`), so the phylograms are the trees the alignments
    were actually drawn along. Any other modifier (the ``Markov`` clock, a draw among
    families), a second lineage clock, or a non-``PerSite`` scope raises.

    Rate variation **across sites** rides on ``model``, not on ``substitution``:
    ``model=hky85(2.0).across_sites(gamma_shape=0.5, invariant=0.1)`` sorts the sites into a
    discretised-Gamma set of rate classes plus a class that never changes. The two axes are
    orthogonal and compose — the clock says which *lineages* run fast, the model which *sites* do.

    ``partitions`` splits a family's sites into blocks each under **its own** model, in place of
    ``model=`` and ``length=``::

        partitions=[(hky85(kappa=2.0), 600), (jc69(), 400)]

    — a 1000-site gene whose first 600 sites evolve under HKY85 and whose last 400 evolve under
    JC69, concatenated in that order into one sequence per gene copy. Each partition's model may
    carry its own ``across_sites`` classes, which is how the field usually spells a codon-position
    split. Giving ``model=`` or ``length=`` alongside is refused rather than merged: the partitions
    already answer both, and the length is their site counts summed. Every partition must be over the
    **same alphabet**, because they concatenate into one sequence.

    All the partitions share the run's one substitution rate, and so the family keeps its one
    phylogram — exactly, not approximately: every model is normalised to one expected substitution
    per site per unit branch length, and every set of rate classes to a mean of 1, so a branch of
    ``Δt`` accrues ``rate · Δt`` substitutions per site in each partition alike. There is
    deliberately no per-partition speed; it would make one phylogram a weighted average of the trees
    the partitions were really drawn along. Family and ordered runs only — a **nucleotide** run's
    blocks already carry their own lengths and their own gene/spacer models, so it refuses
    ``partitions``. Experimental (SPEC §9): Python API, no CLI flag yet.

    ``profiles`` gives chosen families a **site profile** — an ``(L, K)`` array whose row *i* is the
    equilibrium frequencies at position *i*, in the model's own ``alphabet`` order::

        profiles={4: array_of_shape_300_by_20}

    Every model on the menu gives a gene one set of frequencies shared by all its sites; a profile
    says which residues belong at *each* position instead, which is what a buried hydrophobic site
    and the loop beside it differ by. The base model's **exchangeabilities are kept** — which pairs
    interchange easily is chemistry, and the profile's business is only where each residue belongs —
    so the site's matrix is the run's model over the site's own frequencies. Families you leave out
    evolve under ``model`` untouched, and a family's profile changes nothing about any other family.

    A row of exact zeros is refused: it makes that site's matrix degenerate, and a real profile says
    a residue is unlikely at a position rather than impossible — add a pseudocount. A **flat**
    profile, every row the model's own frequencies, is the model without one; it is statistically
    identical rather than byte-identical, because L single-site models consume the random stream
    differently from one L-site model. Profiles compose with ``across_sites`` — a profile says *which*
    residues, ``+Γ`` says *how fast* — and are refused alongside ``partitions`` (both decide a
    family's per-site models) and alongside ``parallel`` (which ships one shared partition set to
    every worker).

    An **amino-acid** profile needs a protein model and so belongs to a family or ordered run: a
    nucleotide genome is measured in base pairs and read on either strand, so it refuses protein
    models outright. Profiles still apply there, over the four bases — a row per base pair — and the
    row count must equal that block's length, which the genome run already fixed. Experimental
    (SPEC §9): Python API, no CLI flag yet.

    On a **nucleotide** genome run every root block is evolved — spacer as well as genes — each at its
    own length in bp, so ``length`` does not apply and is rejected. ``model`` evolves the genes and
    ``intergene_model`` (default ``jc69``) the spacer, at ``intergene_speed`` times the rate (default
    ``3.0``). Each carries **its own** across-site variation: decorating ``model`` with
    ``across_sites`` does not reach the spacer, whose default ``jc69()`` stays flat — the spacer's
    job is to be the unconstrained null, and silently giving it the genes' Gamma would make it
    something else. Give ``intergene_model`` a decorated model to vary the spacer too. Because the whole genome is covered, the run also **puts the genomes back together**:
    ``.node_genomes`` holds every node's chromosomes, blocks concatenated in physical order — the complete
    tree, reconstructed — and ``.initial_genome`` the one the run started with.

    The result carries the **phylograms** the sequences were drawn along — each gene tree and the
    species tree, with branch lengths converted from time to substitutions/site by the same
    ``base × clock × Δt``.

    ``parallel`` opts into evolving the gene trees **concurrently** — one gene tree per worker
    process, which is where a run's time goes when the sequences are long or a nucleotide genome has
    thousands of blocks. ``False`` (the default) runs the serial engine above. ``True`` uses every
    core; a positive ``int`` sets the worker count. It is a **separate engine**: each family draws
    from its own RNG stream (spawned from ``seed``), so every worker count returns the *same* bytes,
    but that realisation differs from the serial one for a given seed — parallel is a speed choice,
    made once, not a drop-in for the default. Threads would not help here (numpy releases the GIL too
    little for these array sizes), so this is process-backed and only pays off above a work threshold.
    Because it spawns processes, a script that calls it with ``parallel`` set must guard its entry with
    ``if __name__ == "__main__":`` (the standard multiprocessing requirement); the ``zombi2`` CLI
    already does, so ``--parallel`` there needs nothing extra.

    ``stream_to=DIR`` writes each family's files as it finishes and keeps nothing, returning a
    `StreamedSequences` handle instead of a `SequencesResult`. This is the level where a run's memory
    actually goes — every alignment and every ancestral sequence live at once — so it is the dial for
    a run whose result would not fit. The files are the same ones ``.write(DIR)`` would leave, so a
    streamed run and an in-memory one at the same seed are the same dataset; ``outputs`` picks which,
    exactly as ``.write`` does, and ``flat`` is passed through the same way. It composes with
    ``parallel``. A **nucleotide** run cannot stream: it puts whole genomes back together, and that
    needs every block's sequence at once, which is the opposite of keeping nothing.
    """
    from ..genomes import NucleotideGenomesResult, OrderedGenomesResult, StreamedRun, read_run

    # A written run is a genome run too. `zombi2 sequences --from DIR` has always reopened one; from
    # Python the same handoff was a dead end, which mattered most for `stream_to=` — the feature
    # whose whole point is that the run does not fit in memory, and whose handle then could not be
    # passed on. A path or a StreamedRun reads back here, so both front doors take the same step.
    if isinstance(genomes, (str, os.PathLike, StreamedRun)):
        genomes = read_run(genomes)

    nucleotide = isinstance(genomes, NucleotideGenomesResult)
    # An ordered run is admitted here as a family one: this level reads a genome run's `gene_trees`
    # and its `complete_tree`, and the ordered result carries both — the coordinates it adds are
    # simply not something a sequence needs. It used to be refused, because `OrderedGenomesResult`
    # is not a subclass of `FamilyGenomesResult` and the gate tested the class rather than what it
    # had to supply. That contradicted the documentation everywhere ("a family or ordered run", in
    # this module's own docstring four times, in Chapter 7 and in Appendix B), and the CLI never hit
    # it because its directory handoff rebuilds a FamilyGenomesResult from `genome_events.tsv` on
    # the way in — so `zombi2 genomes --resolution ordered` into `zombi2 sequences` worked while the
    # same two calls in Python did not.
    if not nucleotide and not isinstance(genomes, (FamilyGenomesResult, OrderedGenomesResult)):
        raise TypeError(
            f"the sequence level runs on a genome run, got {type(genomes).__name__} — pass the "
            "result that genomes.simulate_genomes_family(...), simulate_genomes_ordered(...) or "
            "simulate_genomes_nucleotide(...) returned: the whole run, not its .gene_trees. A "
            "sequence lives inside a gene, but its clock rides the *species* branch that gene sits "
            "on — one draw per lineage, shared by every family — so the run needs the species tree "
            "too."
        )
    if stream_to is not None and nucleotide:
        raise ValueError(
            "a nucleotide run cannot stream: it reassembles every node's genome, which needs every "
            "block's sequence in memory at once — the opposite of what streaming does. Run it "
            "in memory, or use --resolution family at the genome level.")
    if outputs is not None and stream_to is None:
        raise ValueError(
            "outputs applies to a streamed run (stream_to=DIR), which writes the files itself; for "
            "an in-memory run choose them when you call result.write(outputs=...).")
    species_tree = genomes.complete_tree
    # With partitions the models arrive inside them, and `_resolve_partitions` checks each one; this
    # is the plain path, where `model` is the whole answer and the common mistake is worth naming
    # exactly as it always was.
    if partitions is None and not isinstance(model, (SubstitutionModel, Models)):
        if model is None:
            raise ValueError(
                "no model: give model=… (one substitution model for every site of every family) "
                "together with length=…, or partitions=[(model, sites), …] to split a family's "
                "sites into blocks each under its own model, or "
                "Models().set_by(Clade({…}), {…}) to give each clade its own.")
        raise TypeError(f"model must be a SubstitutionModel (e.g. hky85(kappa=2.0)) or a "
                        f"per-lineage set (Models().set_by(Clade({{…}}), {{…}})), got {model!r}")
    if intergene_model is not None and not isinstance(intergene_model, (SubstitutionModel, Models)):
        raise TypeError(f"intergene_model must be a SubstitutionModel or a per-lineage set "
                        f"(Models().set_by(Clade({{…}}), {{…}})), got {intergene_model!r}")

    if profiles is not None and partitions is not None:
        raise ValueError(
            "profiles and partitions both decide which model each site of a family evolves under, so "
            "giving both leaves no rule for which wins. A profile already gives every site its own "
            "model — if you want blocks of sites to differ in something a profile cannot say, that is "
            "partitions; if you want per-site frequencies, that is profiles. Not both.")
    if profiles is not None and parallel:
        raise ValueError(
            "profiles are not implemented for the parallel engine yet: it ships one shared partition set to "
            "every worker, and a profile is per family. Run this level serially (drop `parallel`), or "
            "drop `profiles`.")
    site_profiles = {} if profiles is None else _resolve_profiles(profiles, model, length)

    if nucleotide:
        # Every recovered root block evolves — spacer as well as genes — so the run reconstructs the
        # whole genome rather than the declared loci. Each block brings its own length in bp, which
        # is why a single `length` would contradict the coordinates the genome recorded.
        if partitions is not None:
            raise ValueError(
                "partitions do not apply to a nucleotide genome run: every block already carries "
                "its own length in bp and its own model — `model` for a gene, `intergene_model` for "
                "the spacer — so the genome has already said which stretch takes which. Drop "
                "partitions; the genome sets both the lengths and the split. Partitions are for a "
                "family or ordered run, where a family is one undivided sequence until you divide "
                "it.")
        # `partitions is None` here, so the guard above already refused a missing or non-model
        # `model` — from this point it is the gene model, not an optional one.
        for label, m in (("model", model), ("intergene_model", intergene_model)):
            if isinstance(m, Models):
                raise ValueError(
                    f"{label} is a per-lineage model set, which a nucleotide genome run does not "
                    f"read yet: this path evolves blocks of a whole genome rather than a family's "
                    f"gene trees, and a block's stretch of a species branch is chosen by the "
                    f"genome's coordinates rather than by the gene tree walk the set is applied on. "
                    f"Use one model here, or run the family / ordered resolution.")
        assert isinstance(model, SubstitutionModel)
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
        # flat and parameterless: the null for unconstrained DNA. Named separately so the
        # argument's optionality ends here rather than trailing through the loops below.
        spacer: SubstitutionModel = intergene_model if intergene_model is not None else jc69()
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
            per_block[i] = (b - a, model if is_gene else spacer,
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
            f_model = model if (src, a, b) in genic else spacer
            if f_model.alphabet != BASES:
                raise ValueError(
                    f"the run was founded from a FASTA (DNA), but {f_model.name} is a protein model — "
                    "a nucleotide sequence cannot found an amino-acid alignment")
            founding_seed[i] = encode(root[a:b], f_model.alphabet)
        parts = None                # a nucleotide block's model and length come from `per_block`
    else:
        gene_trees = genomes.gene_trees
        parts = _resolve_partitions(model, partitions, length)
        # A per-lineage model set is painted against THIS run's tree, once, before any family is
        # evolved — the same shape as a driven rate resolving its trajectory above. It is checked
        # here too: a label that names no lineage, or a lineage with no model, is a run that would
        # otherwise quietly evolve part of the tree under the wrong matrix.
        parts = tuple((m.resolve(species_tree) if isinstance(m, Models) else m, n) for m, n in parts)
        if intergene_model is not None:
            raise ValueError(
                "intergene_model applies to a nucleotide genome run, where blocks are genes or "
                "spacer. A family or ordered run has gene families only, so there is nothing "
                "for a second model to evolve.")
        per_block = None
    if divergence is not None:
        rate = _calibrate(substitution, divergence, genomes.complete_tree)
    else:
        rate = as_rate(1.0 if substitution is None else substitution, default_scope=PerSite)
    if rate.scope is not PerSite:
        assert rate.scope is not None      # `as_rate` and `_calibrate` both fill the scope in
        raise ValueError(
            f"substitution has a {rate.scope.__name__} scope, but the sequence engine reads the "
            f"substitution rate per site and cannot read it any other way. Write PerSite(...), or "
            f"drop the scope and let the level fill in its own."
        )
    # The rate's modifiers, sorted into the two things this level reads. SPEC §5: modifiers multiply,
    # so a clock and a driver compose — one says which lineages were dealt a fast tempo, the other
    # what their state makes of it — and the gate below rejects only what the level cannot honour.
    clocks = tuple(m for m, _ in rate.carried_modifiers(unit='lineages'))
    drivers = driven_mods(rate)
    # This level is the one that does NOT take a third-party modifier, so the gate is a plain
    # isinstance rather than `is_implemented`. Every other engine evaluates its rate through
    # `Rate.effective`, which multiplies in whatever `factor()` returns; this one reads its two kinds
    # of modifier itself — the clock is *drawn among lineages* before any site evolves, not evaluated
    # at an event — so a modifier declaring itself implemented here would be accepted and then never
    # called. Silently returning the undriven answer is precisely what SPEC §5 forbids, so it is
    # refused by name instead, and `Modifier.implemented_for` documents the omission.
    unimplemented = sorted({describe(m) for m in rate.modifiers
                            if not matches_declared(m, IMPLEMENTED_MODIFIERS)})
    if unimplemented:
        raise ValueError(
            f"substitution carries {', '.join(unimplemented)}, which the sequence engine does not read. It "
            "takes a lineage clock — varying_among('lineages', LogNormal(0.0, 0.3)) (uncorrelated) or "
            "varying_among('lineages', Drift(LogNormal(0.0, 0.3))) (autocorrelated), and "
            "several of one kind compose — and any number of scaled_by drivers, which multiply. "
            "set_by is not read here (a replaced base has nowhere to go: this level draws its clock "
            "among lineages rather than evaluating a rate), and neither is the Markov clock, a draw "
            "among families, or a modifier of your own: this "
            "level reads its modifiers directly rather than through the rate, so one it did not ship "
            "could not be honoured. Rate variation across sites is not a modifier "
            "at all — it belongs to the model: model=hky85(...).across_sites(gamma_shape=0.5), or "
            "--gamma-shape."
        )
    # SPEC §5's one-memory-structure-per-axis rule, in the one place every level calls: a lineage's
    # factor is either drawn afresh or inherited and perturbed, and those are two accounts of the same
    # thing. Several of one kind compose, as any two modifiers do; a Driven is a different axis and
    # composes with either.
    check_one_memory(clocks, label="substitution", unit="lineages")
    for m in drivers:
        if isinstance(m.mapping, Between):
            raise ValueError(
                "substitution carries scaled_by(..., Between(...)), and a donor/recipient kernel is "
                "meaningless in a rate: a rate is read on one lineage, and there is no second lineage "
                "for the pair's first half to name. A Between belongs in the genome level's "
                "transfer_to weight, where the two ends of a transfer exist. Weight the "
                "substitution rate by the lineage's own state instead — scaled_by(driver, {state: "
                "factor})."
            )
        if names_a_live_level(m.driver):
            raise ValueError(
                f"substitution is driven by {m.driver!r}, which names a level growing beside the run "
                "— the joint spelling of a driver (SPEC §5). Traits and Sequences cannot be joined "
                "(SPEC §3): a sequence lives inside a gene and never feeds back into the trait, so "
                "there is nothing for the two to decide together. Grow the trait first and condition "
                "on it — pass the TraitsResult, or the path to the trait_events.tsv it wrote."
            )
    rate_base = rate.base
    # The only rate with no base of its own is one whose base a `set_by` replaces, and the gate above
    # refuses a `set_by` here — this level draws its clock rather than evaluating a rate, so there is
    # nowhere for a replaced base to go.
    assert rate_base is not None

    # Conditioning: resolve each driver ONCE into a DriverTrajectory (value + next-switch lookups,
    # keyed by the shared species node id), before choosing an engine — this is shared validation and
    # shared input, not one engine's business, exactly as the genome level does it. A mapping whose
    # states never occur in the driver would leave every branch at the default factor, so the run
    # would be the undriven model wearing a driven rate; that is refused here, naming the driver.
    # No driven rate ⇒ `driven` is empty and everything below is what it was.
    driven: list = []
    if drivers:
        by_key: dict = {}
        for m in drivers:
            by_key.setdefault(m.key, m)
        trajs = {key: resolve_driver(m.driver, species_tree, step=m.step, level="sequences")
                 for key, m in by_key.items()}
        for m in drivers:
            label = m.driver if isinstance(m.driver, str) else f"<{type(m.driver).__name__}>"
            check_mapping_fires(m.mapping, trajs[m.key].states(), driver_label=label)
            # A scheduled mapping entry makes the factor a function of TIME, and this level does not
            # read one: `IMPLEMENTED_MODIFIERS` leaves out OnTime on purpose, so `changing_at` is
            # refused here and a schedule reaching in through a mapping would be the same model by
            # another door. It has to be refused rather than run, because a rate this level cannot
            # step would silently hold the schedule's opening factor for the whole run — the branch
            # lengths would come out as if the schedule were a plain number, and nothing would say so.
            if m.mapping.next_change(0.0) != math.inf:
                raise ValueError(
                    "a time schedule inside a mapping is not read at the sequences level: this level "
                    "walks each gene tree branch by branch and does not step at a wall-clock time, "
                    "so the schedule's first factor would stand for the whole run. Give this state "
                    "one factor, or put the schedule on the level that grows along the tree — a "
                    "genome or a trait rate — and drive the substitution rate from what that "
                    "produced.")
        driven = [(m, trajs[m.key]) for m in drivers]

    names = species_tree.labels()   # e<id> for a lineage that died; n<id> for the rest
    sink = None
    if stream_to is not None:
        chosen = tuple(outputs) if outputs is not None else _DEFAULT_STREAM_OUTPUTS
        unknown = [o for o in chosen if o not in _WRITE_OUTPUTS]
        if unknown:
            raise ValueError(f"unknown stream outputs {unknown}; choose from {list(_WRITE_OUTPUTS)}")
        sink = _Sink(stream_to, chosen, "family", flat)
    alignments: dict[int, dict[str, str]] = {}
    ancestral: dict[int, dict[str, str]] = {}
    founding: dict[int, str] = {}
    phylograms: dict[int, dict[str, str | None]] = {}
    seed = resolve_seed(seed)      # drawn if none was given, so either engine below records it
    if not parallel:
        # Serial reference engine — the default, left exactly as it was. One shared generator draws the
        # clock, then each family is walked in turn. `parallel` selects a *separate* engine (decision A),
        # so turning it on gives a different-but-valid realisation for a seed; this path never changes.
        rng, _ = stream("sequences", seed)
        clock = resolve_clock(clocks, driven, species_tree, gene_trees, rng)
        # One transition-CDF cache per model, shared across every block that model evolves. Branch lengths
        # recur across blocks (a block passing straight through a species branch reuses its length), so a
        # run-wide cache builds a few hundred matrices where a per-block cache rebuilt tens of thousands.
        # Keyed by model identity — genes and spacer are different models and must not share a cache.
        # A profile keyed to a family this run does not have is a typo that would otherwise apply
        # to nobody and say nothing — the same silence the Driven mapping guard exists to break.
        stray = sorted(set(site_profiles) - set(gene_trees), key=str)
        if stray:
            raise ValueError(
                f"profiles names {len(stray)} family/families this run does not have: "
                f"{', '.join(map(str, stray[:5]))}{' …' if len(stray) > 5 else ''}. The run has "
                f"{len(gene_trees)} of them, keyed {min(gene_trees, key=str)}…{max(gene_trees, key=str)}. "
                f"A profile for a family that is not here applies to nothing.")
        cdf_caches: dict[int, dict[float, np.ndarray]] = {}
        bar = progress_bar(len(gene_trees), "sequences", unit="family", enabled=progress)
        for family in sorted(gene_trees):  # sorted for reproducibility given the seed
            bar.update()
            gt = gene_trees[family]
            if per_block is None:
                f_parts, f_rate = parts, rate_base
            else:                       # a nucleotide block: its own length, and spacer runs faster
                f_len, f_model, speed = per_block[family]
                f_parts, f_rate = ((f_model, f_len),), rate_base * speed
            if family in site_profiles:          # a profile replaces the family's model, site by site
                f_parts = site_profiles[family]
                # `length` is rejected on a nucleotide run — every block carries its own — so the
                # row-count check in `_resolve_profiles` had nothing to compare against and this is
                # where it lands. Without it a short profile silently shortens the sequence, and the
                # alignment stops agreeing with the coordinates the genome run wrote.
                if per_block is not None and len(f_parts) != per_block[family][0]:
                    raise ValueError(
                        f"the profile for block {family!r} has {len(f_parts)} rows but that block is "
                        f"{per_block[family][0]} bp. A profile carries one row per site, and on a "
                        f"nucleotide run the genome already fixed the length.")
            seed_states = None if per_block is None else founding_seed[family]
            aln, anc, fnd = _evolve_partitions(gt, f_parts, f_rate, clock, rng, cdf_caches, names,
                                               founding=seed_states)
            scaled = _scaled_gene_tree(gt, f_rate, clock)  # branch lengths in subs/site
            ext = scaled.extant
            phylo = {"complete": _gene_newick(scaled.complete, names),
                     "extant": _gene_newick(ext, names) if ext is not None else None}
            if sink is None:
                alignments[family], ancestral[family] = aln, anc
                founding[family], phylograms[family] = fnd, phylo
            else:
                sink.family(family, aln, anc, fnd, phylo)   # written and forgotten
        bar.close()
    else:
        # Parallel engine (opt-in): one gene tree per process, each under its own spawned RNG stream, so
        # any worker count is bit-identical to any other. The clock is a shared read-only draw, so it is
        # taken once here from a reserved stream (index 0) and shipped to every worker.
        from .._runtime.parallel import guard_pool_workers, resolve_workers
        from ._pergenetree import evolve_families
        workers = guard_pool_workers(resolve_workers(parallel))
        spawned = seed_sequence("sequences", seed)[0].spawn(1 + len(gene_trees))
        clock = resolve_clock(clocks, driven, species_tree, gene_trees,
                              np.random.default_rng(spawned[0]))
        alignments, ancestral, founding, phylograms = evolve_families(
            gene_trees, per_block, model, intergene_model, length, rate_base, clock,
            founding_seed if nucleotide else None, spawned[1:], workers, progress, names,
            sink=None if sink is None else sink.family, partitions=parts)

    sp_scaled = _scaled_species_tree(species_tree, rate_base, clock)   # the clock made visible
    sp_extant = prune(sp_scaled, keep="extant")
    species_phylogram = {"complete": sp_scaled.to_newick(),
                         "extant": sp_extant.to_newick() if sp_extant is not None else None}
    # A nucleotide run evolved every block, so **every** node's genome can be put back together —
    # one map, as at the genome level. Which sequences each node reads is the split the level already
    # makes: an extant tip's genes are tips of their block trees, everything else's are not. The
    # concatenation is deferred to read-time (see `_AssembledGenomes`): only the cheap per-node
    # layout is captured now, so hundreds of megabases of genome do not all sit in memory beside the
    # per-block sequences they are built from.
    assembled: "Mapping[str, dict[int, str]]" = {}
    initial_genome: dict[int, str] = {}
    insertion_plan: dict = {}                  # only a nucleotide run with insertions has one
    if nucleotide:
        # Capture the layouts in the same order the eager build used — extant nodes (read from
        # `alignments`) sorted first, then the rest (read from `ancestral`) — so the map iterates,
        # and `write` emits its files, in exactly the previous order.
        extant_ids = sorted(species_tree.extant_leaves())
        extant_id_set = set(extant_ids)
        extant_labels = {names[i] for i in extant_ids}
        ordered_ids = extant_ids + [i for i in sorted(species_tree.nodes) if i not in extant_id_set]
        layouts = {names[i]: genomes.assembly(i) for i in ordered_ids}
        _gap_what_is_not_carried(layouts, alignments, ancestral, genomes.root_blocks)
        insertion_plan = _insertion_plan(genomes, layouts, names)
        # `alignments` / `ancestral` become the spliced view; the per-block rows underneath stay as
        # they are, because that is what a genome is assembled from — the assembly slices a block's
        # OWN coordinates, and the spliced columns are not in that frame.
        block_rows, block_ancestral = alignments, ancestral
        widths = {i: b - a for i, (_s, a, b) in enumerate(genomes.root_blocks)}
        folded = frozenset(b for _host, (runs, _c) in insertion_plan.items()
                           for (_off, b, _w) in runs)
        spliced: "Mapping[int, dict[str, str]]" = _SplicedAlignments(
            block_rows, insertion_plan, widths, folded)
        spliced_anc: "Mapping[int, dict[str, str]]" = _SplicedAlignments(
            block_ancestral, insertion_plan, widths, folded)
        alignments, ancestral = spliced, spliced_anc            # type: ignore[assignment]
        assembled = _AssembledGenomes(layouts, block_rows, block_ancestral, extant_labels)
        # The genome the run started with. Its blocks were all laid down at the start, so each one's
        # sequence there is its `founding` draw — the state the stem leads *from*. It is not a node,
        # so it is in neither map above; the same reason `founding` is not in `ancestral`.
        for cid, pieces in genomes.initial_assembly().items():
            initial_genome[cid] = "".join(
                piece if strand == 1 else piece.translate(_COMPLEMENT)[::-1]
                for (block, strand, lo, hi) in pieces
                for piece in (founding[block][lo:hi],))

    if sink is not None:
        sink.finish(species_phylogram)
        handle = StreamedSequences(str(stream_to), seed, sink.n_families, sink.n_sequences,
                                   sink.outputs, sink.identity, sink.sites, sink.n_ancestral)
        if "summary" in sink.outputs:
            # same rule as the in-memory write: the file describes the directory it sits in, so it
            # counts ancestral sequences only where they were actually written. A streamed run and an
            # in-memory one at the same seed are the same dataset, summary included.
            written = handle.summary()
            if "ancestral" not in sink.outputs:
                written.pop("ancestral_sequences", None)
            write_summary(pathlib.Path(stream_to) / "sequences_summary.json", written)
        return handle
    return SequencesResult(alignments, ancestral, founding, phylograms, species_phylogram, seed,
                           assembled, initial_genome, "block" if nucleotide else "family",
                           # a nucleotide run's models are all forced to DNA above (and its `parts` is
                           # None — each block brings its own); elsewhere every partition shares one
                           # alphabet, so the first one speaks for the run
                           BASES if parts is None else parts[0][0].alphabet,
                           tuple(names[i] for i in sorted(species_tree.extant_leaves())),
                           insertion_plan,
                           alignments._raw if isinstance(alignments, _SplicedAlignments) else alignments,
                           ancestral._raw if isinstance(ancestral, _SplicedAlignments) else ancestral)


__all__ = ["simulate_sequences", "SequencesResult", "StreamedSequences",
           "mean_pairwise_identity",
           # the substitution-model menu, re-exported: the TypeError raised for a bad
           # `model=` names these symbols, so they must be importable from the module it names
           "jc69", "k80", "hky85", "gtr", "poisson", "jtt", "dayhoff", "wag", "lg",
           "SubstitutionModel",
           # the per-lineage model set: the TypeError for a bad `model=` names it too
           "Models"]
