"""Figure: one circular chromosome as a signed gene order.

Chapter 5 (ordered genomes) reads a leaf as ``[ 0+ 0+ 1+ 3+ 4− ]`` — a run of gene tokens, each a
family with a strand. This draws that same chromosome as the ring it is: five genes spaced evenly
around a circle, each an arrow that points **the way its strand reads**. Colour is by family, so the
two copies of family ``0`` — a tandem duplication — share a colour and sit adjacent, and family ``4``,
left backwards by an inversion, is the one arrow pointing against the others.

The genome is the leaf ``n4`` of the chapter's first example (``seed=231``); the figure is generated
from that run, not drawn by hand, so it cannot drift from what the code produces. Change the
chapter's seed and this has to be re-run — the two constants below are the chapter's.

**Drawn by Phylustrator**, which is the answer to a question the manual review asked of this figure:
the chunky arrow bent along a ring is `plot(genome, layout="circular")` with `gene_style="arrow"`,
the same call the gallery's inversion example makes. It used to be a hand-rolled ``drawsvg`` ring —
about eighty lines of trigonometry for arcs, arrowheads and label placement — which meant this one
figure drew genes by rules nothing else in the project shared. The family numbers around the ring
are `gene_labels`, added to Phylustrator for it.

House style: no title inside the figure (the manual captions it). The colours are the categorical
four the gallery paints traits and states with, so a *family* reads as an identity here and there.

Run:  python figures/scripts/fig_ordered_chromosome.py
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import phylustrator as ph

from zombi2 import species
from zombi2.genomes import simulate_genomes_ordered
from zombi_style import save, FS_LABEL

#: Families are told apart by hue, not by shade. Two genes of one family share a colour, so a tandem
#: duplication reads as two neighbours of one colour; four hues are enough for this chromosome and
#: they are the gallery's, which keeps one categorical vocabulary across the project.
PALETTE = ["#4C9AA6", "#E08A3C", "#8B6B9E", "#3C8D6E"]

SEED = 231                    # the chapter's first example
LEAF = "n4"                   # the leaf it reads

W = 760                       # square: a ring has no long side


def _leaf_chromosome(node_label: str = LEAF):
    """The chromosome of one extant leaf of the chapter's first example."""
    tree = species.simulate_species_tree(birth=1.0, death=0.1, n_extant=4, seed=SEED)
    g = simulate_genomes_ordered(tree, duplication=0.3, loss=0.2, origination=0.15, inversion=0.5,
                                 chromosomes=1, initial_families=5, seed=SEED)
    return g.node_genomes[int(node_label[1:])][0]


def render() -> None:
    chrom = _leaf_chromosome(LEAF)
    genes = [ph.genomes.Gene(family=str(gene.family), strand=gene.strand, position=i)
             for i, gene in enumerate(chrom.genes)]
    genome = ph.genomes.Genome(LEAF, [ph.genomes.Chromosome(str(chrom.id), genes,
                                                            topology=chrom.topology)])
    families = sorted({g.family for g in genes}, key=int)
    palette = {fam: PALETTE[i % len(PALETTE)] for i, fam in enumerate(families)}

    # ring_gene_frac slims the body so the flared head reads as a head rather than a bulge; the
    # margin is what the family numbers outside the ring sit in.
    style = ph.Style(width=W, height=W, margin=int(W * 0.17), gene_stroke_width=1.0,
                     ring_gene_frac=0.20, font_size=FS_LABEL)
    figure = (ph.genomes.plot(genome, layout="circular", style=style)
              + ph.genomes.genes(by="family", palette=palette)
              + ph.genomes.gene_labels(pad=0.13))
    save(figure.as_svg(), "ordered_chromosome")
    print(f"  {LEAF}: [ " + " ".join(f"{g.family}{'+' if g.strand > 0 else '-'}" for g in genes)
          + " ]  " + ", ".join(f"{f}={palette[f]}" for f in families))


if __name__ == "__main__":
    render()
