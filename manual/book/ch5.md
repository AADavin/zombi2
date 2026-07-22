# Genomes II: ordered

The previous chapter put genes on the tree as a *bag of families* — how many copies of each, and nothing more. This chapter gives them **structure**. A genome becomes one or more **chromosomes**, each an ordered run of genes, and each gene knows which way it points. This is the **ordered** resolution.

Everything from Chapter 4 is still here and still works the same way: the four events, the rate grammar, transfers, the event log, the gene trees, the profiles. Position is orthogonal to genealogy — two genes being neighbours says nothing about how they are related — so for the same events the gene trees and profiles are exactly the ones the unordered engine would produce. What is new is *which chromosome a gene sits on*, *where along it*, and *which way it points*, plus the events that change those things.

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
leaf n2, chromosome 2 (circular):  [ 0+ 1+ 2+ 2+ 2+ 3+ 4+ ]
```

Each gene is written as its family with the strand as `+` or `−` (the strand is the integer `+1` or `−1`). This leaf has one chromosome of seven genes, in which family `2` sits in a run of three tandem copies. The gene tree of a family is unchanged from Chapter 4 — the true genealogy, read off the same event log:

```python
g.gene_trees[0].to_newick("extant")
# (g27:0.679048,g13:0.679048)speciation_n0;
```

## The karyotype

A genome is seeded with a **karyotype**: `chromosomes=N` chromosomes, each with a `topology` — `"circular"` (the default) or `"linear"`, or a per-chromosome list like `["circular", "linear"]` for a mixed set. The founding `initial_families` genes are dealt round-robin across them. Topology is not just a label: it decides where a segmental event stops, which the section on segments below takes up.

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

Chromosomes carry a genuine **identity**. A chromosome id is re-minted at every event that reshapes it — a speciation, a fission, a fusion — and each of those edges is recorded. So the run leaves behind not just the chromosomes at the tips but the *genealogy* that connects them: the **chromosome network**. It is the middle of three tiers that nest, the species tree containing the chromosome network, which contains the gene trees:

```
species tree  ⊃  chromosome network  ⊃  gene trees
```

It is a **network** and not a tree because of one event: **fusion joins two chromosome lineages into one**, two parents and one child. Fission and speciation are ordinary splits (one parent, two children); origination is a root; loss is a leaf. The whole thing is a directed graph, and it is recorded the way graphs are, as an **edge list** — `chromosome_events`, one row per event. The run above gives:

```
  time   kind          parents -> children
  0.00   origination        -  -> 0          a seed chromosome
  0.00   origination        -  -> 1          a seed chromosome
  0.97   loss               1  -> -          chromosome 1 (and its genes) dies
  2.19   speciation         0  -> 2, 3
  3.27   speciation         2  -> 4, 5
  3.35   fission            4  -> 6, 7        a bifurcation
  3.67   fusion          6, 7  -> 8           a reticulation (two parents)
```

You can read the karyotype's whole history off this table. Two seed chromosomes; one dies early; the other descends through the species splits; then chromosome `4` fissions into `6` and `7`, which promptly fuse back into `8`. The leaf below that fusion holds a single chromosome, `8`, carrying `[ 0+ 0+ 5+ 2+ ]`. No `eNewick` string is involved: a multi-rooted, reticulating graph is not a tree, so it is not forced into Newick.

## Events act on segments

Chromosomes are half of what structure buys you. The other half is this: once genes have neighbours, **a gene-level event acts on a segment** — a run of consecutive genes — not on a single gene. A duplication copies a *block*; a loss removes a *block*; a transfer sends a *block* sideways. This is what produces the signature of real genome evolution: neighbouring genes that share a history because they were copied, moved, or lost *together*.

How long a segment? Its length — the **extension** — is drawn per event from a distribution you set per event type, `<event>_extension`. The default is `Geometric(mean=1)`, which is degenerate at one gene, so out of the box every event touches a single gene and you recover the simplest behaviour. Raise the mean to make blocks longer.

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
            the block 4 1, duplicated as a unit and landed in tandem
```

The block `4 1` appears twice: a single segmental duplication copied those two adjacent genes together. (Family `0` is absent — it was lost earlier, which is what left `4` and `1` next to each other.) A duplication puts its copy **in tandem**, immediately after the original run; a transferred block arrives together on the recipient. **Origination is the exception**: a family is born once, as a single new gene, so it has no extension.

## Circular chromosomes have no ends

A segment runs rightwards from the gene it starts on. Where it stops is set by the chromosome's `topology`. On a **linear** chromosome the run stops at the last gene. On a **circular** one there is no last gene, so the run carries on past the first: on a ring the first gene and the last are neighbours.

That is what the example above shows. Family `0` was lost early, leaving the chromosome `1 2 3 4 4`; the duplicated block was the trailing `4` together with the leading `1`, a run across the origin.

A run that crosses the origin re-anchors the chromosome: position 0 moves so that the run sits at the front. Nothing biological changes, because on a circle position 0 is an index and not a feature of the molecule. That is why the leaf above starts at family `4` rather than at family `1`.

A run is never longer than the chromosome. If the extension distribution asks for more genes than there are, the run is the whole chromosome. A loss then empties it. The empty chromosome stays in the karyotype (only `chromosome_loss` removes one), exactly as a de-novo replicon starts out empty. An inversion reverses the whole ring, which is the same molecule read the other way round.

The distinction is not cosmetic. If runs on a circular chromosome stopped at position 0, every run that started near the end would be cut short. The genes around the origin would be duplicated, lost and moved less often than the rest, and blocks would come out shorter than the extension you asked for. Wrapping removes both effects. On a linear chromosome both are kept, because a linear replicon really does have ends.

## Rearrangements: inversion, transposition, translocation

Three further events act on a segment and reshape the order without creating or destroying genes. They are **identity-preserving** — a gene keeps its id, so nothing is written to the gene genealogy. They only reorder, and they are logged separately.

- **Inversion** *(per chromosome)* — reverse a segment in place, flipping the strand of every gene in it. The classic signed-permutation move: `+2 +3 +4` becomes `−4 −3 −2`. On a circular chromosome the segment may span the origin; reversal on a ring is well defined.
- **Transposition** *(per chromosome)* — cut a segment out and reinsert it **elsewhere on the same chromosome**.
- **Translocation** *(per gene copy)* — move a segment to a **different chromosome** of the same genome. A no-op if the genome has only one chromosome.

A moved block — transposed or translocated — lands **inverted** with probability `inversion_probability` (default `0`, so it keeps its orientation).

```python
tree = species.simulate_species_tree(birth=1.0, death=0.1, n_extant=3, seed=0)
g = simulate_genomes_ordered(
    tree, duplication=0.15, loss=0.15, origination=0.1,
    inversion=0.3, transposition=0.25, translocation=0.2,
    transposition_extension=Geometric(mean=2), inversion_probability=0.5,
    chromosomes=2, initial_families=6, seed=0)
```

All three land in one `rearrangements` log:

```
  translocation  chrom 0 [0:1]  -> chrom 1, pos 2      flipped=False
  translocation  chrom 8 [2:3]  -> chrom 6, pos 1      flipped=False
  inversion      chrom 5 [0:1]
  transposition  chrom 7 [0:1]  -> pos 0 (same chrom)  flipped=False
  transposition  chrom 6 [1:2]  -> pos 1 (same chrom)  flipped=True
```

Every row names its run the same way: `start` is the position the run began at, `length` is how many genes it covered, counted rightwards from `start` and wrapping past the origin on a circular chromosome. So `start + length` larger than the chromosome's gene count means the run crossed the origin. The destination of a transposition or a translocation is an index into what was left after the run was cut out, so it can never fall inside the run itself.

Translocation is the one that carries a gene lineage **across** to another chromosome, which is why it is counted per gene copy, like transfer. But note what it does *not* do: the chromosomes themselves come through unchanged, so a translocation writes no edge into the chromosome network. Only a gene has moved between them.

## The events, and their older names

Readers from the rearrangement and comparative-genomics literature already have names for these events; as in every chapter, the names live in one table and organise nothing.

| What it does | Here | From the literature |
|---|---|---|
| split one chromosome, or merge two | `fission=…`, `fusion=…` | Chromosome fission / fusion |
| a de-novo chromosome appears | `chromosome_origination=…` | Plasmid / replicon gain |
| the chromosome genealogy as a whole | fission + fusion → `chromosome_events` | Karyotype / chromosome-number evolution |
| a block of adjacent genes copied together | `duplication` + `duplication_extension=` | Segmental duplication |
| reverse a gene segment, strands flipped | `inversion=…` | Inversion / reversal (GRIMM, MGR) |
| cut a segment, paste it elsewhere on the same chromosome | `transposition=…` | Transposition |
| a segment moves to another chromosome | `translocation=…` | Reciprocal translocation |

## The `OrderedGenomesResult` object

`simulate_genomes_ordered` returns an **`OrderedGenomesResult`**, the ordered counterpart of `GenomesResult` — the same spine, with the structured extras:

- `.complete_tree` — the species tree the genomes ran on, extinct lineages included.
- `.genomes` — a dict from node id to that node's genome, now a tuple of **`Chromosome`** objects. Each `Chromosome` has an `id`, a `topology`, and an ordered list of **`Gene`** objects (`id`, `family`, `strand`).
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

# segmental everything: duplications, losses and inversions act on blocks of genes
g = simulate_genomes_ordered(
    tree, duplication=0.2, loss=0.25, inversion=0.3,
    duplication_extension=Geometric(mean=4), loss_extension=Geometric(mean=3),
    inversion_extension=Geometric(mean=5), initial_families=15, seed=1)

# rearrangements: relocate and move blocks between chromosomes, sometimes inverting them
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
zombi2 genomes -t out/species_complete.nwk --resolution ordered \
    --origination 0.5 --fission 0.05 --fusion 0.05 --chromosomes 2 --seed 1 -o out/

# segmental duplications, losses and inversions on three chromosomes
zombi2 genomes -t out/species_complete.nwk --resolution ordered \
    --duplication 0.2 --loss 0.2 --origination 0.5 --inversion 0.3 --chromosomes 3 --seed 1 -o out/

# blocks relocate and move between chromosomes, sometimes inverting
zombi2 genomes -t out/species_complete.nwk --resolution ordered \
    --origination 0.5 --transposition 0.2 --translocation 0.1 \
    --inversion-probability 0.5 --chromosomes 2 --seed 1 -o out/
```

The segment-length distributions (how long a duplicated or inverted run is) are structured objects, so like the rate modifiers they stay in Python; the CLI takes the rates as numbers. The same tree-reading rules apply as in Chapter 4, and `--write` selects among the richer ordered outputs.

## Outputs

`.write(dir, outputs=[...])` materialises the chosen products to disk:

```python
g.write("out/", outputs=("events", "profiles", "gene_order",
                         "rearrangements", "chromosome_events"))
```

```
out/genome_events.tsv        the gene genealogy (the source of truth)
out/profiles.tsv             family × extant-species copy counts
out/gene_order.tsv           the observed genomes' layout, one row per gene
out/rearrangements.tsv       inversions, transpositions, translocations
out/chromosome_events.tsv    the chromosome network (edge list)
```

The first three are written by default; the two logs are opt-in. `gene_order.tsv` is the ordered genome's headline output — the signed gene order of every observed leaf:

```
species  chromosome  position  strand  family  gene
21       21          8         1       8       218
21       21          9         -1      9       219
27       27          0         1       0       270
```

and `chromosome_events.tsv` is the network's ground truth, its columns the edge list above — `time · kind · lineage · parents · children`, one row per event. The full list of files lives in Appendix B.
