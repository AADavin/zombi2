# Genomes II: ordered

The previous chapter put genes on the tree as a *bag of families*, how many copies of each, and nothing more. This chapter gives them **structure**. A genome becomes one or more **chromosomes**, each an ordered row of genes, and each gene carries a strand, the direction it reads in. This is the **ordered** resolution.

![An example of an ordered genome. Genes acquire a relative position to each other in a circular chromosome. The orientation of each gene is also registered.](figures/ordered_chromosome.pdf){width=58%}

## Chromosomes split, merge, appear and die

A genome has a **karyotype**: `chromosomes=N` chromosomes, each with a `topology`: `"circular"` (the default) or `"linear"`, or a per-chromosome list like `["circular", "linear"]` for a mixed set. The founding `initial_families` genes are distributed across the chromosomes in turn. The number of chromosomes can also evolve through different events.

The karyotype itself evolves: four events change the **number** of chromosomes ([Ge10](https://aadavin.github.io/zombi2/gallery.html#genomes)<!--gallery:genome_karyotype-->):

- **`fission`** *(per chromosome)*. A chromosome splits in two at a uniform cut, both halves keeping its topology; a single-gene chromosome cannot split.
- **`fusion`** *(per chromosome)*. Two chromosomes of a genome merge into one. Only two of the **same topology**: a ring and a molecule with two ends cannot become one molecule, so a circular chromosome never fuses with a linear one, and a genome holding one of each never fuses at all.
- **`chromosome_origination`** *(per lineage)*. A de-novo replicon appears: a new chromosome, empty and circular, a plasmid.
- **`chromosome_loss`** *(per chromosome)*. A whole chromosome dies, and every gene on it is recorded as a loss. A lineage never loses its *last* chromosome this way.

## The chromosome network

Chromosomes are tracked. A chromosome gets a new id at every event that passes it on or reshapes it: a speciation, a fission, a fusion. Each of those edges is recorded. So the run leaves behind not just the chromosomes at the tips but the *genealogy* that connects them: the **chromosome network**. It is the middle of three genealogies that nest, the species tree containing the chromosome network, which contains the gene trees:

```
species tree  ⊃  chromosome network  ⊃  gene trees
```

The containment is at each instant: every gene sits on one chromosome, and every chromosome in one lineage. That holds even though a gene's lineage can change chromosome by translocation and species by transfer.

It is a **network** and not a tree because of one event: **fusion joins two chromosome lineages into one**, two parents and one child. Fission and speciation are ordinary splits (one parent, two children); `initial` and `origination` are roots: the chromosomes the run began with, and the new replicons `chromosome_origination` creates, kept apart so you can tell which is which. Loss is a leaf. It is a directed graph, and it is written the way graphs are, as an **edge list**: `chromosome_events.tsv` on disk, `.chromosome_events` in Python, one row per event. A run with fissions and fusions gives:

```
  time   kind          parents -> children
  0.00   initial            -  -> 0          a chromosome the run started with
  0.00   initial            -  -> 1          a chromosome the run started with
  0.75   fission            0  -> 2, 3       a bifurcation
  1.52   fusion          2, 1  -> 4          a reticulation (two parents)
  1.52   fission            4  -> 5, 6
  2.02   fusion          3, 5  -> 7
  2.11   loss               7  -> -          chromosome 7 (and its genes) dies
  2.19   speciation         6  -> 8, 9
  2.49   fission            9  -> 10, 11
  3.27   speciation         8  -> 12, 13
  3.30   fusion        10, 11  -> 14
```

The rows carry no lineage column; which species holds a chromosome is read off `gene_order.tsv`, where every chromosome id sits beside its node.

## Events act on segments

Once genes have neighbours, **a gene-level event acts on a segment**, a stretch of consecutive genes, not on one gene. A duplication copies a segment, a loss removes one, a transfer sends one sideways. That produces the signature of real genome evolution: neighbouring genes sharing a history because they were copied, moved or lost *together*.

How many genes are affected simultaneously by one event? That is its **extent**. By default, extents are computed from a geometric distribution, but this can be changed to other probability distributions. A bare number is the mean of the draw, so `duplication_extent=3` copies about three adjacent genes. The default is a single gene, so every event takes one gene unless you set an extent.

Topology decides where a segment stops. On a circular chromosome a segment that reaches the last gene continues from the first, wrapping position 0, so only the whole chromosome bounds it; on a linear one it stops at the last gene, and a draw that would overrun the end is cut short there.

One consequence of events affecting multiple genes at the same time is that **a rate counts starts, not hits.** `duplication=0.2` with `duplication_extent=3` does *not* mean each gene is duplicated 0.2 times per unit time. A gene is copied whenever any duplication's segment covers it, which happens about three times as often as the starts alone, so roughly `0.2 × 3 = 0.6`. If you want genes duplicated at a known rate, divide that rate by the mean extent.

The same reading holds at the nucleotide resolution of Chapter 5, where the extent is in base pairs rather than genes.

## Rearrangements: inversion, transposition, translocation

Three more events reshape the order without creating or destroying genes:

- **Inversion.** Reverse a segment in place, flipping every gene's strand: `2+ 3+ 4+` becomes `4− 3− 2−` ([Ge9](https://aadavin.github.io/zombi2/gallery.html#genomes)<!--gallery:genome_inversion-->).
- **Transposition.** Cut a segment out and reinsert it **elsewhere on the same chromosome**, at a spot drawn uniformly over what is left after the cut.
- **Translocation.** Move a segment to a **different chromosome** of the same genome, landing at a random position. A no-op if the genome has only one chromosome.

All three are counted **per copy**: every gene is a potential start, so a bigger genome rearranges more often; each also takes `PerLineage`, a fixed budget independent of genome size (Appendix A). A moved segment, transposed or translocated, lands **inverted** with probability `inversion_probability` (default `0`, so it keeps its orientation).

## On the command line

The ordered resolution is `--resolution ordered`. It adds two sets of flags to the Chapter 3 events, each still a plain number: the chromosome flags, `--chromosomes N` and `--topology` taking `circular` (the default), `linear`, or a comma-separated list of one per chromosome; and the **extent** flags, `--inversion-extent`, `--duplication-extent` and the rest, each the mean number of genes an event takes. Leave an extent out and every event takes a single gene, which for an inversion means flipping one gene's strand and shuffling nothing:

```bash
# chromosomes split and merge along the tree
zombi2 genomes out/ --resolution ordered \
    --origination 0.5 --fission 0.05 --fusion 0.05 --chromosomes 2 --seed 1

# segmental duplications, losses and inversions on three chromosomes
zombi2 genomes out/ --resolution ordered \
    --duplication 0.2 --loss 0.2 --origination 0.5 \
    --inversion 0.3 --chromosomes 3 --seed 1

# segments relocate and move between chromosomes, sometimes inverting
zombi2 genomes out/ --resolution ordered \
    --origination 0.5 --transposition 0.2 --translocation 0.1 \
    --inversion-probability 0.5 --chromosomes 2 --seed 1
```
