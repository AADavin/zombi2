"""Forward simulation of signed gene order under the ZOMBI2 **nucleotide** genome model.

The genome is a set of genes with real base-pair spans on a fixed karyotype (``n_chrom`` linear
chromosomes, matching the real clade). Inversions act at nucleotide coordinates — an inversion
reverses an arc of a chromosome and flips the strand of every gene it fully covers; a breakpoint
that falls inside a gene cuts it, exactly as a real inversion would. No DNA is simulated — only the
gene layout descends the tree — so a run is fast and its output is the observable synteny needs:
each extant genome's signed gene order.

Rate convention. ``inversion`` is a rate **per gene per unit tree time** (per Myr on a dated tree);
the genome-wide nucleotide rate handed to the engine is ``inversion × n_total``. ``mean_length`` is
the mean inversion length **in genes**, converted to base pairs through the fixed per-gene spacing.
"""
from __future__ import annotations

from ete3 import Tree as ETree

from zombi2 import genomes
from zombi2.species import read_newick

GENE_LENGTH = 1000       # bp per gene (an arbitrary unit; only gene order and length-in-genes matter)
SPACING = 1.4            # chromosome bp per gene of coding — ~70% coding, leaving intergenic room


def load_dated_tree(nwk_path: str):
    """Load a dated species tree (branch lengths in Myr). Returns ``(tree, {node id: species})``.

    Real dated trees are ultrametric only up to rounding, so every tip is declared extant rather
    than left for ZOMBI to infer from depth (which it refuses to guess)."""
    nwk = open(nwk_path).read()
    species = [leaf.name for leaf in ETree(nwk, format=1).get_leaves()]
    tree, namemap = read_newick(nwk, tip_fates={s: "extant" for s in species})
    return tree, namemap


def simulate_signed_order(tree, namemap, *, inversion: float, mean_length: float,
                          n_total: int, n_chrom: int, seed: int) -> dict[str, list[tuple]]:
    """Evolve ``n_total`` genes on ``n_chrom`` chromosomes down ``tree`` under inversions only.

    ``inversion`` — inversions per gene per Myr; ``mean_length`` — mean inversion length in genes.
    Returns ``{species: [(chromosome, family, strand ∈ {+1,-1}), …]}`` in genomic order.
    """
    genes_per_chrom = n_total // n_chrom
    n_total = genes_per_chrom * n_chrom
    root_length = int(genes_per_chrom * GENE_LENGTH * SPACING)
    bp_per_gene = root_length / genes_per_chrom
    inv_genome = inversion * n_total                        # per-gene rate -> genome-wide rate
    inv_len_bp = max(1, int(round(mean_length * bp_per_gene)))
    res = genomes.simulate_genomes_nucleotide(
        tree, inversion=inv_genome, inversion_length=inv_len_bp,
        genes=genes_per_chrom, gene_length=GENE_LENGTH, chromosomes=n_chrom,
        root_length=root_length, topology="linear", seed=seed)
    out: dict[str, list[tuple]] = {}
    for node in res.complete_tree.extant():
        genome = res.genomes[node.id]
        out[namemap[node.id]] = [
            (chrom.id, b.gene, res.gene_strands.get(b.gene, 1) * b.strand)
            for chrom in genome.chromosomes for b in chrom.blocks if b.is_gene]
    return out
