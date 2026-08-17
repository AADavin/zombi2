"""The dated tree, and the forward simulation of signed gene order down it.

Both the data this study treats as observed and every candidate the fit proposes come from the same
generator, ZOMBI2's nucleotide genome model. That is the point of the design: the truth is known,
so what the fit recovers can be checked against it.

**All genomes are simulated, down the whole tree.** One run evolves the genome along every branch
and yields a genome at every extant tip, which is what the statistics need, because they are
computed over pairs of tips.

Genes have real base-pair spans on linear chromosomes. An inversion reverses an arc at nucleotide
coordinates and flips the strand of the genes it covers. A translocation moves a segment between
chromosomes. No DNA is simulated, only the gene layout descends the tree, so a run takes about a
second.

Two properties of the nucleotide model matter for reading the results. A breakpoint may never fall
strictly inside a gene, so genes are never split and every genome keeps all of them, one to one.
And the requested extent is snapped to the nearest legal breakpoint, so the extent handed to the
engine is a request rather than the arc that is realised.

Rate convention. ``inversion`` and ``translocation`` are given here as rates **per gene per unit of
tree time**; the engine takes genome-wide rates, so each is multiplied by the gene count before it
is handed over. ``inversion_extent`` is a mean extent **in genes**, converted to base pairs through
the fixed per-gene spacing.
"""
from __future__ import annotations

import re

from zombi2 import genomes, species, tree as treemod

GENE_LENGTH = 1000       # bp of coding per gene
SPACING = 1.4            # chromosome bp per bp of coding, so about 70% coding

_BRANCH = re.compile(r":(\d+(?:\.\d+)?(?:[eE][-+]?\d+)?)")


def dated_tree(*, seed: int, n_tips: int, crown_depth: float, birth: float, death: float):
    """A birth-death tree of ``n_tips`` extant tips, rescaled to ``crown_depth`` time units.

    The tree is ultrametric, so every tip sits at the same depth and the patristic distance
    between two tips runs from 0 to ``2 * crown_depth``. Rescaling only fixes the units the rates
    are quoted in; it does not change the tree's shape.
    """
    result = species.simulate_species_tree(birth, death, n_extant=n_tips, seed=seed)
    raw = result.extant_tree
    depth = max(raw.nodes[i].end_time for i in raw.extant_leaves())
    factor = crown_depth / depth
    newick = _BRANCH.sub(lambda m: f":{float(m.group(1)) * factor:.10g}", raw.to_newick())
    tree, _ = treemod.read_newick(newick, assume_extant=True)
    # A ZOMBI tree's labels *are* its node ids, so read_newick returns an empty name map and
    # `labels()` is the id -> name map the statistics key their genomes on.
    return tree, tree.labels(), newick


def simulate_signed_order(tree, namemap, *, inversion: float, inversion_extent: float,
                          translocation: float = 0.0, n_genes: int, n_chromosomes: int,
                          seed: int) -> dict[str, list[tuple]]:
    """Evolve ``n_genes`` genes on ``n_chromosomes`` linear chromosomes down ``tree``.

    ``inversion`` and ``translocation`` are per gene per unit time. ``inversion_extent`` is the
    mean inversion extent in genes. Returns ``{species: [(chromosome, gene, strand)]}`` in
    genomic order, one entry per extant tip.
    """
    per_chromosome = n_genes // n_chromosomes
    n_genes = per_chromosome * n_chromosomes
    root_length = int(per_chromosome * GENE_LENGTH * SPACING)
    bp_per_gene = root_length / per_chromosome
    extent_bp = max(1.0, inversion_extent * bp_per_gene)
    result = genomes.simulate_genomes_nucleotide(
        tree,
        inversion=inversion * n_genes,
        inversion_extent=extent_bp,
        translocation=translocation * n_genes,
        translocation_extent=extent_bp,
        genes=per_chromosome, gene_length=GENE_LENGTH, chromosomes=n_chromosomes,
        root_length=root_length, topology="linear", seed=seed)
    # `result.genomes` holds one genome per extant tip, keyed by the tip's label — the same label
    # `pairwise_divergence` keys its pairs on, so the two line up without a lookup.
    strands = result.gene_strands
    return {
        label: [(chromosome.id, block.gene, strands.get(block.gene, 1) * block.strand)
                for chromosome in genome.chromosomes
                for block in chromosome.blocks if block.is_gene]
        for label, genome in result.genomes.items()}
