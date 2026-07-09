"""Block-based genome selection (**P2.6**): language-model codon selection on the GENE blocks of a
simulated nucleotide genome, along each block's own gene tree, with intergene blocks left neutral.

The nucleotide genome model (:func:`~zombi2.simulate_nucleotide_genomes`) evolves a real annotated
genome forward under structural events (inversion, duplication, loss, transfer, ...), then decomposes
the outcome into **blocks** -- maximal never-cut ancestral intervals, each with its own gene tree.
Because *Design S* guarantees a gene is never split by a breakpoint, every gene is **exactly one
block** spanning the whole coding interval with a single genealogy. That is what makes per-block
codon selection well-posed: a whole CDS evolves as one unit down one tree.

This module overlays :class:`~zombi2.experimental.codon_selection.CodonSelection` on those gene blocks:

* a **gene block** matched to a supplied :class:`~zombi2.experimental.genome_selection.CDS` (so its
  coding **strand** and reading frame are known) evolves as coding DNA -- mutation on DNA, selection
  on the encoded protein -- **along that block's gene tree**, so duplications / transfers / speciations
  all inherit the ancestral protein's constraint and dN/dS emerges as an output;
* every other block drifts **neutrally** under the same nucleotide model, exactly as the core pipeline
  does: intergene blocks, and any gene block we cannot codon-evolve (a novel *originated* gene that is
  not on the seed chromosome, a frame / stop / ambiguity problem, or a run with no root FASTA).

Evolved block sequences are cached back onto the ``NucleotideResult`` (``_block_seqs`` /
``_seq_species``), so :meth:`~zombi2.genomes.nucleotide_sim.NucleotideResult.node_sequence` and
:meth:`~zombi2.genomes.nucleotide_sim.NucleotideResult.gene_alignments` reassemble the DNA at every
node **unchanged** -- reassembly is agnostic to how each block was evolved.

This is the block-based counterpart of P2.5's
:class:`~zombi2.experimental.genome_selection.GenomeSelection` (one contig, the species tree, no
gene-family dynamics): here selection runs on the *full* structural simulation, each block on its own
DTL gene tree -- the "specify a real genome at the root and let it evolve" scenario.

Notes / v1 scope:

* ``gamma`` (across-site rate heterogeneity) applies to the **neutral** (intergene / fallback) blocks
  only; coding-block rate heterogeneity is *emergent* from selection (conserved sites evolve slower),
  so a nucleotide Gamma is neither needed nor applied there.
* Codon and intergene blocks share one branch-length scale (neutral substitutions per nucleotide
  site); the codon model is normalised per nucleotide site to match, so a gene block at ``beta == 0``
  diverges at the same neutral rate as an intergene block.
* **Limitation:** pseudogenization keeps a block classified ``"gene"`` but flips a lineage to
  non-functional at a ``G`` state-change edge in its tree; this version keeps codon selection on the
  whole gene-block tree, so post-pseudogenization lineages remain under selection (a per-lineage
  switch to neutral at the ``G`` edge is a later refinement).

Like the rest of the feature, nothing here imports torch/esm at module load; *evolving* needs scipy
(part of ``zombi2[selection]``).
"""
from __future__ import annotations

import warnings
from dataclasses import dataclass, field

import numpy as np

from zombi2.experimental import ExperimentalWarning, warn_experimental
from zombi2.experimental.codon_selection import CodonSelection, translate
from zombi2.experimental.genome_selection import CDS
from zombi2.experimental.selection import Critic
from zombi2.sequences.models import SubstitutionModel, evolve_on_tree, reverse_complement

__all__ = ["BlockSelectionReport", "NucleotideGenomeSelection", "simulate_nucleotide_selection"]

_ACGT = frozenset("ACGT")


@dataclass
class BlockSelectionReport:
    """What happened to each block when overlaying selection.

    ``n_selected`` gene blocks evolved under codon selection; ``n_neutral_fallback`` gene blocks fell
    back to neutral (each with a ``(block_id, gene_id, reason)`` in ``fallbacks``); intergene blocks are
    always neutral. ``n_empty`` blocks had no gene tree (nothing originated for them). Under a
    zero-divergence run (no clock, ``subst_rate <= 0``) selection is skipped entirely and every block
    just propagates the root, so ``n_selected`` is 0 and no fallbacks are recorded.
    """

    n_blocks: int = 0
    n_gene_blocks: int = 0
    n_selected: int = 0
    n_intergene: int = 0
    n_neutral_fallback: int = 0
    n_empty: int = 0
    fallbacks: list = field(default_factory=list)  # [(block_id, gene_id, reason)]

    def _note_fallback(self, block_id: int, gene_id, reason: str) -> None:
        self.n_neutral_fallback += 1
        self.fallbacks.append((block_id, gene_id, reason))


class NucleotideGenomeSelection:
    """Overlay language-model codon selection on the gene blocks of a nucleotide-genome simulation.

    ``critic`` / ``beta`` / ``nuc_model`` configure an inner
    :class:`~zombi2.experimental.codon_selection.CodonSelection`; the **same** ``nuc_model`` drives the
    neutral intergene evolution, so coding and non-coding share one mutation process.
    """

    def __init__(self, critic: Critic, *, beta: float = 1.0,
                 nuc_model: SubstitutionModel | None = None):
        warn_experimental("NucleotideGenomeSelection")
        # the inner CodonSelection is an implementation detail; don't leak its experimental warning
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", ExperimentalWarning)
            self.codon = CodonSelection(critic, beta=beta, nuc_model=nuc_model)
        self.nuc = self.codon.nuc

    def evolve_blocks(self, result, cds, *, root_fasta: str | None = None, gamma=None,
                      subst_rate: float = 1.0, clock=None, rng: np.random.Generator | None = None,
                      seed: int | None = None) -> BlockSelectionReport:
        """Evolve every block of ``result`` and cache the sequences back onto it.

        Mirrors :meth:`~zombi2.genomes.nucleotide_sim.NucleotideResult.simulate_sequences` (same gene
        trees, branch lengths and reassembly contract) but sends **gene blocks matched to a ``CDS``**
        through codon selection. ``cds`` is the coding annotation of the *root* genome (a ``list[CDS]``
        from :func:`~zombi2.experimental.genome_selection.read_cds_gff`, carrying strand + phase); a
        gene block is matched to its CDS by its ancestral ``[start, end)`` span (exact under Design S).
        ``root_fasta`` is the root genome (upper-cased here); gene blocks need it (their coding DNA is a
        real slice). Returns a :class:`BlockSelectionReport`; afterwards ``result.node_sequence(node)``
        and ``result.gene_alignments()`` work.
        """
        from zombi2.genomes.reconciliation import _node_tree
        from zombi2.sequences.evolution import SequenceEvolution, _annotate

        if rng is None:
            rng = np.random.default_rng(seed)
        if root_fasta is not None:
            root_fasta = root_fasta.upper()                # match P2.5: soft-masking is not preserved
            if len(root_fasta) != result.root_length:
                raise ValueError(f"root_fasta length {len(root_fasta)} != root_length "
                                 f"{result.root_length}")

        # by-span CDS lookup: Design S makes a gene block's ancestral span identical to its CDS span
        by_span: dict[tuple[int, int], CDS] = {}
        for c in cds:
            key = (c.start, c.end)
            if key in by_span:
                raise ValueError(f"two CDS share the span {key}: {by_span[key].name!r} and {c.name!r}")
            by_span[key] = c

        zero = clock is None and subst_rate <= 0.0     # no substitutions: root propagates unchanged
        se = clock or (None if zero else SequenceEvolution(root_rate=subst_rate))
        segments = None if zero else se._lineage_segments(result.species_tree, rng)[0]
        total_age = result.species_tree.total_age
        records_by_block, species_by_block = result._block_records()
        seed_source = result._seed_source()

        rep = BlockSelectionReport()
        block_seqs: dict[int, dict] = {}
        for a in result.blocks:
            rep.n_blocks += 1
            is_gene = a.kind == "gene"
            if is_gene:
                rep.n_gene_blocks += 1
            else:
                rep.n_intergene += 1

            root_node = _node_tree(records_by_block[a.block_id], species_by_block[a.block_id], total_age)
            if root_node is None:
                block_seqs[a.block_id] = {}
                rep.n_empty += 1
                continue
            subst: dict = {}
            if not zero:
                _annotate(root_node, segments, max(0.0, se.family_speed.sample(rng)), subst)
            root_seq = (root_fasta[a.start:a.end]
                        if (root_fasta is not None and a.source == seed_source) else None)

            # a gene block on the seed chromosome, with a real coding slice and a matching CDS -> select
            # (skipped when zero: there is no divergence, so don't pay the critic for a no-op)
            if is_gene and not zero and root_seq is not None:
                c = by_span.get((a.start, a.end))
                if c is None:
                    rep._note_fallback(a.block_id, a.gene_id, "no matching CDS for this gene block")
                else:
                    try:
                        block_seqs[a.block_id] = self._evolve_gene_block(root_node, subst, root_seq, c, rng)
                        rep.n_selected += 1
                        continue
                    except ValueError as exc:               # frame / stop / ambiguity -> neutral, reported
                        rep._note_fallback(a.block_id, a.gene_id, str(exc))
            elif is_gene and not zero:
                rep._note_fallback(a.block_id, a.gene_id,
                                   "gene block not on the seed chromosome (novel origination?)"
                                   if root_fasta is not None else "no root_fasta for gene block")

            block_seqs[a.block_id] = self._evolve_neutral(root_node, subst, root_seq, a.length, gamma, rng)

        result._block_seqs = block_seqs
        result._seq_species = species_by_block
        return rep

    def _evolve_gene_block(self, root_node, subst, root_seq: str, c: CDS, rng) -> dict:
        """Codon-evolve one gene block along its gene tree, honouring the CDS strand. ``root_seq`` is the
        forward (source) slice; for a - strand gene it is reverse-complemented to the 5'->3' coding
        sequence before evolving, a terminal stop codon is split off and kept fixed, then each node's
        result is reverse-complemented back to forward orientation -- so the stored block stays
        source-forward, which is what ``node_sequence`` expects."""
        if c.phase != 0:
            raise ValueError(f"CDS {c.name!r} has phase {c.phase} != 0 (phased/partial not supported)")
        coding = root_seq if c.strand == 1 else reverse_complement(root_seq)
        bad = next((i for i, ch in enumerate(coding) if ch not in _ACGT), None)
        if bad is not None:
            raise ValueError(f"CDS {c.name!r} coding position {bad} is not ACGT ({coding[bad]!r})")
        stop = ""
        if len(coding) >= 3 and translate(coding[-3:]) == "*":     # keep a terminal stop fixed
            stop, coding = coding[-3:], coding[:-3]
        if not coding:
            raise ValueError(f"CDS {c.name!r} has no sense codons to evolve")
        ev = self.codon.evolve_coding_family(root_node, subst, coding, rng=rng)
        ev = {gid: seq + stop for gid, seq in ev.items()}
        if c.strand == -1:
            ev = {gid: reverse_complement(seq) for gid, seq in ev.items()}
        return ev

    def _evolve_neutral(self, root_node, subst, root_seq, length: int, gamma, rng) -> dict:
        """Neutral nucleotide evolution for an intergene / fallback block. Keeps the seed slice (with
        ``gamma`` across-site heterogeneity) when we have one, freezing any non-ACGT base (N / gap /
        IUPAC -- it carries no signal); otherwise draws a random root of ``length`` from the model."""
        ev = evolve_on_tree(root_node, subst, self.nuc, rng, root_seq=root_seq, length=length, gamma=gamma)
        if root_seq is not None:
            fixed = {i for i, ch in enumerate(root_seq) if ch not in _ACGT}
            if fixed:
                ev = {g: "".join(root_seq[i] if i in fixed else s[i] for i in range(len(s)))
                      for g, s in ev.items()}
        return ev


def simulate_nucleotide_selection(species_tree, genome: str, cds, *, critic: Critic,
                                  beta: float = 1.0, nuc_model: SubstitutionModel | None = None,
                                  gamma=None, subst_rate: float = 1.0, clock=None,
                                  seed: int | None = None, **sim_kwargs):
    """Turnkey P2.6: evolve a real annotated genome down ``species_tree`` with ESM2 codon selection.

    Derives ``gene_intervals`` from ``cds`` (so every gene block maps 1:1 to a CDS by span), runs the
    structural nucleotide simulation, then overlays selection on the gene blocks. Extra structural
    rates (``inversion=``, ``duplication=``, ``loss=``, ``transfer=`` ...) pass through ``sim_kwargs``
    to :func:`~zombi2.simulate_nucleotide_genomes`. Returns ``(result, report)``; then call
    ``result.node_sequence(node)`` for the evolved genome at any node.
    """
    from zombi2.genomes.nucleotide_sim import simulate_nucleotide_genomes

    cds = list(cds)
    gene_intervals = [(c.start, c.end, c.name) for c in cds]
    result = simulate_nucleotide_genomes(species_tree, gene_intervals=gene_intervals,
                                         root_length=len(genome), retain_internal=True,
                                         seed=seed, **sim_kwargs)
    sel = NucleotideGenomeSelection(critic, beta=beta, nuc_model=nuc_model)
    report = sel.evolve_blocks(result, cds, root_fasta=genome, gamma=gamma,
                               subst_rate=subst_rate, clock=clock, seed=seed)
    return result, report
