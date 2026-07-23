# Genomes II: ordered

The previous chapter put genes on the tree as a *bag of families* — how many copies of each, and nothing more. This chapter gives them **structure**. A genome becomes one or more **chromosomes**, each an ordered run of genes, and each gene knows which way it points. This is the **ordered** resolution.

```python
from zombi2 import species
from zombi2.genomes import simulate_genomes_ordered

tree = species.simulate_species_tree(birth=1.0, death=0.1, n_extant=4, seed=2)
g = simulate_genomes_ordered(
    tree, duplication=0.3, loss=0.2, origination=0.15, inversion=0.5,
    chromosomes=1, initial_families=5, seed=2)
```

Reading one extant leaf:

```
leaf n2, chromosome 2 (circular):  [ 0+ 1+ 3+ 3+ 4− ]
```

Each gene is written as its family with the strand as `+` or `−` (the strand is the integer `+1` or `−1`). This leaf has one chromosome of five genes, in which family `3` sits in a run of two tandem copies and family `4` points backwards, left that way by an inversion. The gene tree of a family is unchanged from Chapter 4 — the true genealogy, read off the same event log:

```python
g.gene_trees[0].to_newick("extant")
# (((g24:0.0396561,g28:0.0396561)speciation_n3:0.405992,g19:0.445648)speciation_n1:0.2334,g9:0.679048)speciation_n0:0.118056;
```

![Leaf `n2`, the same chromosome `[ 0+ 1+ 3+ 3+ 4− ]` drawn as the ring it is. Each gene is an arrow that points the way its strand reads, and its shade marks its family. The two copies of family `3` are a tandem duplication — one shade, side by side; family `4`, left backward by an inversion, is the one arrow pointing against the flow.](figures/ordered_chromosome.pdf){width=58%}

## The karyotype

A genome has a **karyotype**: `chromosomes=N` chromosomes, each with a `topology` — `"circular"` (the default) or `"linear"`, or a per-chromosome list like `["circular", "linear"]` for a mixed set. The founding `initial_families` genes are dealt round-robin across them. Topology is not just a label: it decides where a segmental event stops, which the section on segments below takes up.

## Chromosomes split, merge, appear and die

On top of the karyotype, four events change the **number** of chromosomes:

- **`fission`** *(per chromosome)* — a chromosome splits in two.
- **`fusion`** *(per chromosome)* — two chromosomes of a genome merge into one.
- **`chromosome_origination`** *(per lineage)* — a de-novo replicon appears: a new chromosome, empty and circular, a plasmid.
- **`chromosome_loss`** *(per chromosome)* — a whole chromosome dies, and every gene on it is recorded as a loss. A lineage never loses its *last* chromosome this way.

```python
tree = species.simulate_species_tree(birth=1.0, death=0.1, n_extant=3, seed=42)
g = simulate_genomes_ordered(
    tree, duplication=0.15, loss=0.1, origination=0.25,
    chromosomes=2, fission=0.25, fusion=0.25,
    chromosome_origination=0.03, chromosome_loss=0.03,
    initial_families=5, seed=42)
```

## The chromosome network

Chromosomes are tracked. A chromosome id is re-minted at every event that reshapes it — a speciation, a fission, a fusion — and each of those edges is recorded. So the run leaves behind not just the chromosomes at the tips but the *genealogy* that connects them: the **chromosome network**. It is the middle of three tiers that nest, the species tree containing the chromosome network, which contains the gene trees:

```
species tree  ⊃  chromosome network  ⊃  gene trees
```

It is a **network** and not a tree because of one event: **fusion joins two chromosome lineages into one**, two parents and one child. Fission and speciation are ordinary splits (one parent, two children); origination is a root; loss is a leaf. The whole thing is a directed graph, and it is recorded the way graphs are, as an **edge list** — `chromosome_events`, one row per event. The run above gives:

```
  time   kind          parents -> children
  0.00   origination        -  -> 0          an initial chromosome
  0.00   origination        -  -> 1          an initial chromosome
  0.97   loss               1  -> -          chromosome 1 (and its genes) dies
  2.19   speciation         0  -> 2, 3
  3.27   speciation         2  -> 4, 5
  3.35   fission            4  -> 6, 7        a bifurcation
  3.67   fusion          6, 7  -> 8           a reticulation (two parents)
```

## Events act on segments

Once genes have neighbours, **a gene-level event acts on a segment** — a run of consecutive genes — not on a single gene. A duplication copies a segment; a loss removes a segment; a transfer sends a segment sideways. This is what produces the signature of real genome evolution: neighbouring genes that share a history because they were copied, moved, or lost *together*.

How long a segment? Its length — the **extension** — is drawn per event from a distribution you set per event type, `<event>_extension`. The default is `Geometric(mean=1)`, which is degenerate at one gene, so out of the box every event touches a single gene and you recover the simplest behaviour. Raise the mean to make segments longer.

```python
from zombi2.rates.distributions import Geometric

tree = species.simulate_species_tree(birth=1.0, death=0.1, n_extant=3, seed=27)
g = simulate_genomes_ordered(
    tree, duplication=0.35, loss=0.3,
    duplication_extension=Geometric(mean=3),      # duplications copy ~3 adjacent genes at once
    chromosomes=1, initial_families=5, seed=27)
```

```
leaf n2:  [ 4+ 1+ 4+ 1+ 2+ 3+ ]
            └────┘ └────┘
            the segment 4 1, duplicated as a unit and landed in tandem
```

The segment `4 1` appears twice: a single segmental duplication copied those two adjacent genes together. (Family `0` is absent — it was lost earlier, which is what left `4` and `1` next to each other.) A duplication puts its copy **in tandem**, immediately after the original run; a transferred segment arrives together on the recipient. **Origination is the exception**: a family is born once, as a single new gene, so it has no extension.

## Rearrangements: inversion, transposition, translocation

Three further events act on a segment and reshape the order without creating or destroying genes:

- **Inversion** — reverse a segment in place, flipping the strand of every gene in it. The classic signed-permutation move: `+2 +3 +4` becomes `−4 −3 −2`. On a circular chromosome the segment may span the origin; reversal on a ring is well defined.
- **Transposition** — cut a segment out and reinsert it **elsewhere on the same chromosome**.
- **Translocation** — move a segment to a **different chromosome** of the same genome. A no-op if the genome has only one chromosome.

A moved segment — transposed or translocated — lands **inverted** with probability `inversion_probability` (default `0`, so it keeps its orientation).

```python
tree = species.simulate_species_tree(birth=1.0, death=0.1, n_extant=3, seed=0)
g = simulate_genomes_ordered(
    tree, duplication=0.15, loss=0.15, origination=0.1,
    inversion=0.3, transposition=0.25, translocation=0.2,
    transposition_extension=Geometric(mean=2), inversion_probability=0.5,
    chromosomes=2, initial_families=6, seed=0)
```



## The `OrderedGenomesResult` object

I THINK THIS NEEDS AN UPDATE. DOES ORDERED GENOMES STILL STORE REARRANGEMENTS? IF YES, DISCUSS WITH ME

`simulate_genomes_ordered` returns an **`OrderedGenomesResult`**, the ordered counterpart of `GenomesResult` — the same spine, with the structured extras:

- `.complete_tree` — the species tree the genomes ran on, extinct lineages included.
- `.genomes` — a dict from node id to that node's genome, now a tuple of **`Chromosome`** objects. Each `Chromosome` has an `id`, a `topology`, and an ordered list of **`Gene`** objects (`id`, `family`, `strand`).
- `.initial_genome` — the genome the run **started** with, at the root lineage's origination. It is not `.genomes[root]`: a node sits at the **end** of its branch, and the root branch is real simulated time, so events happen along it. Written to its own `initial_genome.tsv`, with no `lineage` column, because it belongs to no node.
- `.events` — the gene-genealogy log, exactly as in Chapter 4 (origination, duplication, transfer, loss, speciation), from which `.gene_trees` and `.profiles` are derived unchanged. Position and orientation are *not* here; they live in the genomes and in the two logs below.
- `.rearrangements` — the inversion / transposition / translocation log.
- `.chromosome_events` — the chromosome network, as an edge list.
- `.gene_trees`, `.profiles`, `.seed` — as before.

with the methods `.family_counts(node_id)` (the multiset view), `.gene_order(node_id)` (the layout — `(chromosome, position, strand, family, gene id)` per gene), and `.write(dir, outputs=[...])`.

```python
g.genomes[2]                     # the chromosomes of node n2
g.gene_order(2)                  # its layout, gene by gene
g.chromosome_events              # the chromosome network, as an edge list
g.gene_trees[0].to_newick()      # a family's gene tree — unchanged from unordered
```

## Usage from Python

```python
from zombi2 import species
from zombi2.genomes import simulate_genomes_ordered
from zombi2.rates.distributions import Geometric
from zombi2.rates import modifiers as mod

tree = species.simulate_species_tree(birth=1.0, death=0.2, n_extant=30, seed=1)

# several chromosomes that split, merge, and gain a plasmid
g = simulate_genomes_ordered(
    tree, chromosomes=6, topology="linear",
    fission=0.05, fusion=0.05, chromosome_origination=0.02, chromosome_loss=0.02,
    origination=0.4, initial_families=20, seed=1)

# ordered genome: the Chapter 4 events, plus inversions, on a single chromosome
g = simulate_genomes_ordered(
    tree, duplication=0.2, loss=0.2, origination=0.3, inversion=0.3,
    chromosomes=1, initial_families=20, seed=1)

# segmental everything: duplications, losses and inversions act on segments of genes
g = simulate_genomes_ordered(
    tree, duplication=0.2, loss=0.25, inversion=0.3,
    duplication_extension=Geometric(mean=4), loss_extension=Geometric(mean=3),
    inversion_extension=Geometric(mean=5), initial_families=15, seed=1)

# rearrangements: relocate and move segments between chromosomes, sometimes inverting them
g = simulate_genomes_ordered(
    tree, duplication=0.2, transposition=0.2, translocation=0.2,
    inversion_probability=0.5, chromosomes=3, initial_families=15, seed=1)

# rates can still depend on time (the skyline), as at every level
g = simulate_genomes_ordered(
    tree, inversion=1.0 * mod.OnTime({0: 1.0, 2: 0.2}), initial_families=10, seed=1)

# the outputs
g.genomes                             # every node's chromosomes
g.gene_order(next(iter(g.genomes)))   # a node's layout
g.chromosome_events                   # the chromosome network
g.rearrangements                      # the inversion/transposition/translocation log
g.gene_trees                          # one gene tree per family, as in Chapter 4
```

## Usage from the CLI

The ordered resolution is `--resolution ordered`. It adds the chromosome and segment flags to the Chapter 4 events, each still a plain number:

```bash
# chromosomes split and merge along the tree
zombi2 genomes out/ --resolution ordered \
    --origination 0.5 --fission 0.05 --fusion 0.05 --chromosomes 2 --seed 1

# segmental duplications, losses and inversions on three chromosomes
zombi2 genomes out/ --resolution ordered \
    --duplication 0.2 --loss 0.2 --origination 0.5 --inversion 0.3 --chromosomes 3 --seed 1

# segments relocate and move between chromosomes, sometimes inverting
zombi2 genomes out/ --resolution ordered \
    --origination 0.5 --transposition 0.2 --translocation 0.1 \
    --inversion-probability 0.5 --chromosomes 2 --seed 1
```



## Outputs

`.write(dir, outputs=[...])` materialises the chosen products to disk:

```python
g.write("out/", outputs=("events", "profiles", "gene_order",
                         "rearrangements", "chromosome_events"))
```

```
out/genome_events.tsv        the whole history: the genealogy, where each event
                             happened, and the rearrangements — in time order
out/profiles.tsv             family × extant-species copy counts
out/gene_order.tsv           every node's layout, one row per gene
out/chromosome_events.tsv    the chromosome network (edge list)
```

chromosome_events.tsv` is the network's ground truth, its columns the edge list above — `time · kind · lineage · parents · children`, one row per event.

The full list of files lives in Appendix B.
