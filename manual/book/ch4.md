# Genomes II: ordered

The previous chapter put genes on the tree as a *bag of families*, how many copies of each, and nothing more. This chapter gives them **structure**. A genome becomes one or more **chromosomes**, each an ordered row of genes, and each gene carries a strand, the direction it reads in. This is the **ordered** resolution.

```python
from zombi2 import species, genomes

tree = species.simulate_species_tree(birth=1.0, death=0.1, n_extant=4, seed=231)
g = genomes.simulate_genomes_ordered(
    tree, duplication=0.3, loss=0.2, origination=0.15, inversion=0.5,
    chromosomes=1, initial_families=5, seed=231)
```

Reading one extant leaf: `.genomes` holds each tip's chromosomes as `Chromosome` objects (an `id`, a `topology`, and ordered genes each carrying `family` and `strand`), `g.gene_order(node)` gives the same layout as plain tuples, and `gene_order.tsv` writes it for every node, ancestors included (Appendix B):

```
leaf n4, chromosome 4 (circular):  [ 0+ 0+ 1+ 3+ 4− ]
```

Each gene is written as its family with the strand as `+` or `−` (the strand is the integer `+1` or `−1`). This leaf has one chromosome of five genes, in which family `0` sits in a pair of tandem copies and family `4` points backwards, left that way by an inversion. The chromosome is numbered 4 on a run that started with one, because a chromosome gets a new id at every event: the chromosome network section below explains. The gene tree of a family is unchanged from Chapter 3, the true genealogy read off the same event log:

```python
g.gene_trees[0].to_newick("extant")
# (((n5_g25:0.02994724,n6_g30:0.02994724)speciation_n3:0.3253988,(n4_g21:0.227248,n4_g22:0.227248)duplication_n4:0.128098)speciation_n1:0.1181814,n2_g9:0.4735275)speciation_n0:1.090775;
```

Internal nodes are labelled `<event>_n<species>`: the event that ended that gene copy, and the species branch it happened on, so `duplication_n4` is a duplication on branch `n4`.

![Leaf `n4`, the same chromosome `[ 0+ 0+ 1+ 3+ 4− ]` drawn as the ring it is. Each gene is an arrow that points the way its strand reads, and its colour marks its family. The two copies of family `0` are a tandem duplication, one colour side by side; family `4`, left backward by an inversion, is the one arrow pointing against the flow.](figures/ordered_chromosome.pdf){width=58%}

## The karyotype

A genome has a **karyotype**: `chromosomes=N` chromosomes, each with a `topology`: `"circular"` (the default) or `"linear"`, or a per-chromosome list like `["circular", "linear"]` for a mixed set. The founding `initial_families` genes are dealt round-robin across them. Topology is not just a label: it decides where a segmental event stops, which the next section takes up.

## Events act on segments

Once genes have neighbours, **a gene-level event acts on a segment**, a stretch of consecutive genes, not on one gene. A duplication copies a segment, a loss removes one, a transfer sends one sideways. That produces the signature of real genome evolution: neighbouring genes sharing a history because they were copied, moved or lost *together*.

How much does an event take? That is its **extent**, set per event type as `<event>_extent`. A bare number is the mean of a **geometric** draw, so `duplication_extent=3` copies about three adjacent genes: often two or three, sometimes one, occasionally many more. Write `Fixed(3)`, from `zombi2.params.distributions` beside `Geometric`, for exactly three every time, or name any other distribution. The default is a single gene, so out of the box every event touches one gene and you recover the simplest behaviour.

Topology decides where a segment stops, which is the promise the karyotype section made. On a circular chromosome a segment that reaches the last gene continues from the first, wrapping position 0, so only the whole chromosome bounds it; on a linear one it stops at the last gene, and a draw that would overrun the end is cut short there.

```python
from zombi2.params.distributions import Geometric

tree = species.simulate_species_tree(birth=1.0, death=0.1, n_extant=3, seed=332)
g = genomes.simulate_genomes_ordered(
    tree, duplication=0.35, loss=0.3,
    duplication_extent=Geometric(mean=3),      # ~3 adjacent genes copied at once
    chromosomes=1, initial_families=5, seed=332)
```

```
leaf n3:  [ 4+ 1+ 4+ 1+ 2+ 3+ ]
            └────┘ └────┘
            the segment 4 1, duplicated as a unit and landed in tandem
```

The segment `4 1` appears twice: a single segmental duplication copied those two adjacent genes together. (Family `0` is absent because it was lost earlier, and the chromosome is a ring, so the loss closed the gap between `4` and `1`; a printed ring starts at an arbitrary gene.) A duplication puts its copy **in tandem**, immediately after the original segment; a transferred segment arrives as one block, at a random position on a random chromosome of the recipient. **Origination is the exception**: a family is born once, as a single new gene, so it has no extent.

### How often, where, how much

A segmental event answers three questions. Two of them are numbers you set; the third is drawn by the run:

- **how often** it starts, the rate;
- **where** it starts, a gene drawn from the genome;
- **how much** it takes, the extent.

In Chapter 3 the third answer was always "one gene", so the rate was the whole story. Here it is not, and the two numbers **multiply**.

One consequence catches people out: **a rate counts starts, not hits.** `duplication=0.2` with `duplication_extent=3` does *not* mean each gene is duplicated 0.2 times per unit time. A gene is copied whenever any duplication's segment covers it, which happens about three times as often as the starts alone, so roughly `0.2 × 3 = 0.6`. If you want genes duplicated at a known rate, divide that rate by the mean extent.

The same reading holds at the nucleotide resolution of Chapter 5, where the extent is in base pairs rather than genes.

## Rearrangements: inversion, transposition, translocation

Three more events reshape the order without creating or destroying genes:

- **Inversion.** Reverse a segment in place, flipping every gene's strand: `2+ 3+ 4+` becomes `4− 3− 2−`. On a circular chromosome the segment may span the origin ([Ge9](https://aadavin.github.io/zombi2/gallery.html#genomes)<!--gallery:genome_inversion-->).
- **Transposition.** Cut a segment out and reinsert it **elsewhere on the same chromosome**, at a spot drawn uniformly over what is left after the cut.
- **Translocation.** Move a segment to a **different chromosome** of the same genome, landing at a random position. A no-op if the genome has only one chromosome.

All three are counted **per copy**: every gene is a potential start, so a bigger genome rearranges more often; each also takes `PerLineage`, a fixed budget however large the genome (Appendix A). A moved segment, transposed or translocated, lands **inverted** with probability `inversion_probability` (default `0`, so it keeps its orientation).

## Chromosomes split, merge, appear and die

The karyotype itself evolves: four events change the **number** of chromosomes ([Ge10](https://aadavin.github.io/zombi2/gallery.html#genomes)<!--gallery:genome_karyotype-->):

- **`fission`** *(per chromosome)*. A chromosome splits in two at a uniform cut, both halves keeping its topology; a single-gene chromosome cannot split.
- **`fusion`** *(per chromosome)*. Two chromosomes of a genome merge into one. Only two of the **same topology**: a ring and a molecule with two ends cannot become one molecule, so a circular chromosome never fuses with a linear one, and a genome holding one of each never fuses at all.
- **`chromosome_origination`** *(per lineage)*. A de-novo replicon appears: a new chromosome, empty and circular, a plasmid.
- **`chromosome_loss`** *(per chromosome)*. A whole chromosome dies, and every gene on it is recorded as a loss. A lineage never loses its *last* chromosome this way.

```python
tree = species.simulate_species_tree(birth=1.0, death=0.1, n_extant=3, seed=42)
g = genomes.simulate_genomes_ordered(
    tree, duplication=0.15, loss=0.1, origination=0.25,
    chromosomes=2, fission=0.25, fusion=0.25,
    chromosome_origination=0.03, chromosome_loss=0.03,
    initial_families=5, seed=42)
```

## The chromosome network

Chromosomes are tracked. A chromosome gets a new id at every event that hands it on or reshapes it: a speciation, a fission, a fusion. Each of those edges is recorded. So the run leaves behind not just the chromosomes at the tips but the *genealogy* that connects them: the **chromosome network**. It is the middle of three genealogies that nest, the species tree containing the chromosome network, which contains the gene trees:

```
species tree  ⊃  chromosome network  ⊃  gene trees
```

The containment is at each instant: every gene sits on one chromosome, and every chromosome in one lineage. That holds even though a gene's lineage can change chromosome by translocation and species by transfer.

It is a **network** and not a tree because of one event: **fusion joins two chromosome lineages into one**, two parents and one child. Fission and speciation are ordinary splits (one parent, two children); `initial` and `origination` are roots: the chromosomes the run began with, and the new replicons `chromosome_origination` creates, kept apart so you can tell which is which. Loss is a leaf. It is a directed graph, and it is written the way graphs are, as an **edge list**: `chromosome_events.tsv` on disk, `.chromosome_events` in Python, one row per event. The run above gives:

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

## What carries over from Chapter 3

Everything Chapter 3 built runs here unchanged: the four events, the rates and their verbs, the family cap, the named families. Three pieces meet the segment, and each is worth a word.

### Families that differ, once events cover several at once

`varying_among('families', ...)` works here as it does in Chapter 3, but with one difference that matters. A segment covers several families at once, so the weight is applied to **the segment, averaged over the genes it covers**, not to the gene the event happened to start on.

Weighting the starting gene is the obvious implementation and the wrong one: a fast family's rate would then apply to whatever sat beside it, so you would be describing the *neighbourhood* of a fast family rather than the family, and the neighbourhood is reshuffled by every inversion and translocation, so the parameter would not name a stable thing. Averaging over the segment keeps what you wrote true: a segment of heavily-weighted genes is favoured, a mixed one sits between, an ordinary one is unweighted.

With no weights set every segment averages to one, so a run without weights is unchanged.

### Who receives a transfer

The recipient rule is Chapter 3's, unchanged: `transfer_to` takes `"uniform"`, `"distance"` / `Distance(decay=)`, a `Clades(...)` kernel over named clades, or `Recipients().weighted_by(...)` read off a trait (Chapter 8). What is ordered about an ordered transfer is the block that moves; who receives it is the same question and the same answer as at the family resolution.

```python
tree = species.simulate_species_tree(birth=1.0, death=0.3, n_extant=16, seed=1)
flows = genomes.Between({("A", "B"): 1.0, ("B", "A"): 1.0}, default=0.0)
g = genomes.simulate_genomes_ordered(
    tree, transfer=1.0, transfer_extent=3, initial_families=20, seed=2,
    transfer_to=genomes.Clades({"A": ["n27", "n28"], "B": ["n21", "n26"]}, flows))
```

Every transferred block now crosses between the two named clades and never lands inside either. The numbers are weights, normalised over the lineages alive when a transfer occurs, so they redistribute transfers without changing how many happen. A weight of 0 means "cannot receive", so a transfer whose every candidate weighs 0 does not happen at all, leaving the donor exactly as it was.

One thing to watch when you combine a restrictive rule with a tight `max_family_size`: the two thin transfers independently. A block is refused when it would take any family it carries past the cap, and refused again when the kernel forbids the pair, so the realised amount of transfer can sit well below the rate you declared. Raise the cap while you are measuring the weights.

### A rate or an extent can be driven by a trait

Every rate here also takes `scaled_by`, so a habitat can decide how often a lineage rearranges its gene order, and every extent takes it too, so the same habitat can decide how long the rearranged segments are. The mechanism is Chapter 8's and is not repeated here. Every parameter in this chapter is an argument of the one call, so a rate and its extent can each read the same trait, on separate axes. (`Extent(4)` is a bare `4` in wrapper form, there for the verb to chain onto.)

```python
from zombi2 import species, genomes, traits
from zombi2.params import Extent, PerCopy

tree = species.simulate_species_tree(birth=1.0, death=0.2, n_extant=30, seed=1)

# a host-restricted lineage inverts four times as often, and each inversion
# covers three times as many genes
habitat = traits.simulate_discrete(tree, states=["host", "free"], switch=0.5, seed=1)
g = genomes.simulate_genomes_ordered(
    tree,
    inversion=PerCopy(0.3).scaled_by(habitat, {"host": 4.0, "free": 1.0}),
    inversion_extent=Extent(4).scaled_by(habitat, {"host": 3.0, "free": 1.0}),
    initial_families=10, seed=1)
```

What belongs here is why the per-family draw above and a trait driver sit apart. A trait driver attaches to the **lineage**: at any instant it is one factor for that lineage's whole genome, so it composes with any extent unchanged and the run is drawn exactly as it would be without it. `varying_among('families', ...)` attaches to the **contents**, so it has to weight each segment by what it covers. The two therefore cannot be set on the same **rate** yet, and a rate given both is refused with the reason: combining them would mean weighting by the product of a lineage factor and a segment factor, which is neither model on its own. A driven extent, or a driven `transfer_to`, is a different axis and runs alongside a per-family draw unchanged.

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
