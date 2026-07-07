# Genomes

A genome in ZOMBI2 evolves **along a species tree**, and you can model it at three levels of
resolution — pick the coarsest one that answers your question.

<figure markdown="span">
  ![The three genome levels: unordered gene families, ordered chromosomes, and nucleotide genomes.](../img/genome_models.svg){ width="560" }
</figure>

## Gene families (unordered)

The default: a genome is a **bag of gene families** with copy numbers, evolving by
**duplication, transfer, loss and origination** (DTL) — no positional structure. Fast, and
enough for phylogenetic profiles, reconciliation and transfer studies. Rates can be shared,
per-family-sampled (ZOMBI1 style), genome-wise, or per-branch. See
[gene families & rates](gene-families.md), [transfers](transfers.md), and
[gene trees & output](gene-trees-and-output.md).

## Ordered chromosomes

Genes sit on an **ordered, circular chromosome** with a strand, and the genome also undergoes
**inversions and transpositions** — so gene *order* and orientation carry signal, not just
presence and absence. See [ordered genomes](ordered-genomes.md).

## Nucleotide genomes

The finest level: a **nucleotide-resolution** genome of root-anchored segments, with
variable-length structural events (inversions, transpositions, indels), an explicit
gene/intergene structure, homologous replacement, and GFF import to start from a real genome.
Every atom carries its own gene tree, and ancestral sequences can be reconstructed at every
node. See [nucleotide genomes](nucleotide-genomes.md).

---

All three run the same way — `simulate_genomes(tree, ...)` in Python, or `zombi2 genomes` on the
command line (`--genome-model unordered` (default) or `nucleotide`; ordered chromosomes are
selected with the `genome_factory` argument). Growth can be bounded with a hard `max_family_size`
cap or a soft `carrying_capacity` — see [bounding growth](growth.md).
