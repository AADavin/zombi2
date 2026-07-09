"""P2.6: block-based genome selection. ESM2 codon selection on the GENE blocks of a nucleotide-genome
simulation, along each block's own gene tree; intergene blocks neutral; reassembled at every node.

Contract tests are torch-free (a length-agnostic FixedProfile-style critic stands in for ESM2). They
exercise the real structural simulation + the strand/frame handling + the reassembly contract."""
from __future__ import annotations

import ast
import inspect

import numpy as np
import pytest

from zombi2 import BirthDeath, simulate_species_tree
import zombi2.experimental as ex
from zombi2.experimental import (
    CDS, NucleotideGenomeSelection, simulate_nucleotide_selection,
)
from zombi2.experimental.codon_selection import SENSE_CODONS, translate
from zombi2.experimental.selection import Critic
from zombi2.sequences.models import AMINO_ACIDS, reverse_complement

_AA_INDEX = {a: i for i, a in enumerate(AMINO_ACIDS)}


class _NativeCritic(Critic):
    """A torch-free critic that prefers each site's *native* residue, for any protein length: a
    profile peaked (``hi``) on the passed protein's own amino acids. beta then makes it purifying."""

    def __init__(self, hi: float = 0.9):
        self.hi = hi

    def profile(self, protein: str) -> np.ndarray:
        lo = (1.0 - self.hi) / 19.0
        P = np.full((len(protein), 20), lo)
        for i, a in enumerate(protein):
            P[i, _AA_INDEX[a]] = self.hi           # protein is stop-free here, so every residue maps
        return P


def _stop_free(n_codons: int, rng, *, terminal_stop: bool = False) -> str:
    """A random 5'->3' coding DNA of ``n_codons`` sense codons (no internal stops)."""
    codons = [SENSE_CODONS[i] for i in rng.integers(len(SENSE_CODONS), size=n_codons)]
    if terminal_stop:
        codons[-1] = "TAA"
    return "".join(codons)


def _genome_and_cds(rng, length: int = 300):
    """A root genome + a mixed-strand CDS annotation. Coding regions are stop-free on their coding
    strand (a - strand gene stores the reverse-complement in the forward genome); a terminal stop is
    included in one gene to exercise stop handling. Intergenic DNA is random ACGT."""
    cds = [
        CDS(21, 60, +1, "gA"),                     # 39 nt = 13 codons, + strand
        CDS(90, 132, -1, "gB"),                    # 42 nt = 14 codons, - strand
        CDS(162, 201, +1, "gC"),                   # 39 nt, + strand, terminal stop
        CDS(231, 282, -1, "gD"),                   # 51 nt, - strand
    ]
    g = list("".join("ACGT"[i] for i in rng.integers(4, size=length)))   # random ACGT background
    for c in cds:
        n = (c.end - c.start) // 3
        coding = _stop_free(n, rng, terminal_stop=(c.name == "gC"))
        placed = coding if c.strand == 1 else reverse_complement(coding)
        g[c.start:c.end] = list(placed)
    return "".join(g), cds


def _tree(n_tips=6, seed=2):
    return simulate_species_tree(BirthDeath(1.0, 0.3), n_tips=n_tips, age=1.0, seed=seed)


# --------------------------------------------------------------------------- #
# reassembly + strand correctness
# --------------------------------------------------------------------------- #
def test_root_reproduces_input_genome_under_selection():
    # the root has no incoming branch, so every block's root state is its input slice -> node_sequence
    # at the root must reproduce the input genome EXACTLY. This exercises both strands (rc in / rc out)
    # and the terminal-stop split, all under active selection (beta > 0).
    rng = np.random.default_rng(0)
    genome, cds = _genome_and_cds(rng)
    tree = _tree()
    result, report = simulate_nucleotide_selection(
        tree, genome, cds, critic=_NativeCritic(0.9), beta=2.0,
        inversion=0.01, subst_rate=0.5, seed=7)
    assert result.node_sequence(tree.root) == genome
    assert report.n_selected == 4 and report.n_neutral_fallback == 0


def test_all_node_sequences_have_the_right_length():
    rng = np.random.default_rng(1)
    genome, cds = _genome_and_cds(rng)
    tree = _tree()
    result, _ = simulate_nucleotide_selection(
        tree, genome, cds, critic=_NativeCritic(0.9), beta=2.0,
        inversion=0.01, duplication=0.004, loss=0.004, subst_rate=0.4, seed=3)
    for n in tree.nodes_preorder():
        assert len(result.node_sequence(n)) == result.node_genomes[n].size()


# --------------------------------------------------------------------------- #
# report accounting
# --------------------------------------------------------------------------- #
def test_report_counts_gene_and_intergene_blocks():
    rng = np.random.default_rng(2)
    genome, cds = _genome_and_cds(rng)
    result, report = simulate_nucleotide_selection(
        tree := _tree(), genome, cds, critic=_NativeCritic(0.9), beta=1.5,
        inversion=0.01, subst_rate=0.3, seed=5)
    assert report.n_gene_blocks == 4                      # one block per CDS (Design S)
    assert report.n_selected == 4                         # all framed + fasta present
    assert report.n_intergene >= 4                        # gaps between/around the genes
    assert report.n_blocks == report.n_gene_blocks + report.n_intergene


# --------------------------------------------------------------------------- #
# the actual selection signal (the point of P2.6)
# --------------------------------------------------------------------------- #
def test_selection_suppresses_protein_divergence_in_the_block_pipeline():
    # one + strand gene, high divergence: strong purifying selection (large beta, native critic) must
    # leave leaf proteins far closer to the root than neutral codon evolution (beta=0) does.
    rng = np.random.default_rng(4)
    length = 150
    gene = CDS(3, 93, +1, "g")                            # 90 nt = 30 codons, + strand
    g = list("".join("ACGT"[i] for i in rng.integers(4, size=length)))
    g[gene.start:gene.end] = list(_stop_free(30, rng))
    genome = "".join(g)
    root_protein = translate(genome[gene.start:gene.end])

    def mean_leaf_protein_div(beta):
        tree = _tree(n_tips=8, seed=11)
        result, _ = simulate_nucleotide_selection(
            tree, genome, [gene], critic=_NativeCritic(0.95), beta=beta,
            inversion=0.004, subst_rate=1.0, seed=9)
        aln = result.gene_alignments().get("g", {})
        assert aln                                        # the gene survived in some leaves
        tot = n = 0
        for s in aln.values():
            p = translate(s)
            tot += sum(x != y for x, y in zip(p, root_protein))
            n += len(root_protein)
        return tot / n

    neutral = mean_leaf_protein_div(0.0)
    purifying = mean_leaf_protein_div(8.0)
    assert neutral > 0.1                                  # there IS divergence to suppress
    assert purifying < 0.5 * neutral                      # selection roughly halves it (or better)


def test_beta_zero_is_neutral_and_beta_positive_is_purifying():
    # the configured codon model: dN/dS ~ 1 at beta=0, well below 1 under selection
    protein = translate(_stop_free(40, np.random.default_rng(6)))
    crit = _NativeCritic(0.9)
    assert NucleotideGenomeSelection(crit, beta=0.0).codon.dnds(protein) == pytest.approx(1.0, abs=1e-6)
    assert NucleotideGenomeSelection(crit, beta=4.0).codon.dnds(protein) < 0.6


# --------------------------------------------------------------------------- #
# fallbacks (reported, never crash)
# --------------------------------------------------------------------------- #
def test_no_fasta_falls_back_to_neutral_for_gene_blocks():
    rng = np.random.default_rng(7)
    genome, cds = _genome_and_cds(rng)
    tree = _tree()
    gene_intervals = [(c.start, c.end, c.name) for c in cds]
    from zombi2.genomes.nucleotide_sim import simulate_nucleotide_genomes
    result = simulate_nucleotide_genomes(tree, gene_intervals=gene_intervals,
                                         root_length=len(genome), retain_internal=True,
                                         inversion=0.01, seed=4)
    sel = NucleotideGenomeSelection(_NativeCritic(0.9), beta=2.0)
    report = sel.evolve_blocks(result, cds, root_fasta=None, subst_rate=0.3, seed=1)
    assert report.n_selected == 0
    assert report.n_neutral_fallback == report.n_gene_blocks - report.n_empty
    assert all("no root_fasta" in reason for _, _, reason in report.fallbacks)
    # reassembly still works (everything went neutral)
    for n in tree.nodes_preorder():
        assert len(result.node_sequence(n)) == result.node_genomes[n].size()


def test_out_of_frame_gene_falls_back_but_others_are_selected():
    rng = np.random.default_rng(8)
    genome, cds = _genome_and_cds(rng)
    # break the frame of gB by declaring a length not divisible by 3 (41 nt); it must fall back
    cds = [c if c.name != "gB" else CDS(c.start, c.start + 41, c.strand, c.name) for c in cds]
    result, report = simulate_nucleotide_selection(
        _tree(), genome, cds, critic=_NativeCritic(0.9), beta=2.0,
        inversion=0.01, subst_rate=0.3, seed=2)
    assert report.n_selected == 3                         # gA, gC, gD selected
    assert any(gid == "gB" and "multiple of 3" in reason for _, gid, reason in report.fallbacks)


# --------------------------------------------------------------------------- #
# experimental hygiene
# --------------------------------------------------------------------------- #
def test_module_has_no_top_level_ml_imports():
    from zombi2.experimental import nucleotide_selection
    tree = ast.parse(inspect.getsource(nucleotide_selection))
    top: set[str] = set()
    for node in tree.body:
        if isinstance(node, ast.Import):
            top.update(a.name.split(".")[0] for a in node.names)
        elif isinstance(node, ast.ImportFrom) and node.module:
            top.add(node.module.split(".")[0])
    assert "torch" not in top and "esm" not in top, top


def test_exports_stay_in_the_experimental_namespace():
    import zombi2
    for name in ("NucleotideGenomeSelection", "simulate_nucleotide_selection", "BlockSelectionReport"):
        assert name in ex.__all__ and hasattr(ex, name)
        assert not hasattr(zombi2, name)
