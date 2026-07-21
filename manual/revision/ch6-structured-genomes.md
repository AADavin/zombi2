# Genomes II: ordered

The previous chapter put genes on the tree as a *bag of families* — how many copies of each, and nothing more. This chapter gives them **structure**: an order along a chromosome and an orientation, so a genome becomes a signed sequence of genes rather than a multiset. This is the **ordered** resolution, and it is the second of the three: it keeps everything from Chapter 5 and adds *position*. The third resolution — nucleotide, where every gene and gap gets a length in base pairs — is the next chapter.

Everything below layers on the unordered core. The four events, the rate grammar, transfers, the event log, the gene trees, the profiles — all still here, all working the same way. Position is orthogonal to genealogy: two genes being neighbours says nothing about how they are related, and indeed the gene trees and profiles a run produces are exactly the ones the unordered engine would produce for the same events. What is new is *where each gene sits* and *which way it points*, and the events that rearrange that.

```python
from zombi2 import species
from zombi2.genomes import simulate_genomes_ordered

tree = species.simulate_species_tree(birth=1.0, death=0.1, n_extant=4, seed=2)
g = simulate_genomes_ordered(
    tree, duplication=0.3, loss=0.2, origination=0.15, inversion=0.5,
    chromosomes=1, initial_families=5, seed=2)
```

A genome is now a list of **chromosomes**, each an ordered run of genes; a gene knows its `family` and its `strand`. Reading one extant leaf:

```
leaf n2:  [ 0+ 1+ 2+ 2+ 2+ 3+ 4+ ]
```

Each gene is written as its family with the strand as `+` / `−` — the strand is the integer `+1` / `−1`, and in these listings `+` means `+1`. This leaf has one chromosome of seven genes; family `2` sits in a run of three tandem copies, and the whole thing is one linear order. The gene tree of a family is unchanged from Chapter 5 — the true genealogy, read off the same event log:

```python
g.gene_trees[0].to_newick("extant")
# (g27:0.679048,g13:0.679048)speciation_n0
```

## Events act on segments

The one idea that makes ordered evolution different from unordered is this: **every gene-level event acts on a segment** — a run of consecutive genes — not on a single gene. When a duplication fires it copies a *block*; when a loss fires it removes a *block*; a transfer sends a *block* sideways. This is what produces the signature of real genome evolution: neighbouring genes that share a history because they were copied, moved, or lost *together*.

How long a segment? Its length — the **extension** — is drawn per event from a distribution you can set per event type, `<event>_extension`. The default is `Geometric(mean=1)`, which is degenerate at one gene, so out of the box every event touches a single gene and you recover the simplest behaviour; raise the mean to make blocks longer.

```python
from zombi2.rates.distributions import Geometric

tree = species.simulate_species_tree(birth=1.0, death=0.1, n_extant=3, seed=10)
g = simulate_genomes_ordered(
    tree, duplication=0.35, loss=0.3,
    duplication_extension=Geometric(mean=3),      # duplications copy ~3 adjacent genes at once
    chromosomes=1, initial_families=5, seed=10)
```

```
leaf n1:  [ 0+ 1+ 3+ 4+ 1+ 3+ 4+ ]
                └───────┘ └───────┘
                the block 1 3 4, duplicated as a unit and landed in tandem
```

The block `1 3 4` appears twice: a single segmental duplication copied those three adjacent genes together. (Family `2` is absent — it was lost along the way.) A duplication puts its copy **in tandem**, immediately after the original run; a transferred block arrives together on the recipient. **Origination is the exception**: a family is born once, as a single new gene, so it has no extension.

## Rearrangements: inversion, transposition, translocation

Three events reshape the order without creating or destroying genes. They are **identity-preserving** — a gene keeps its id, so nothing is written to the gene genealogy; they only reorder, and are logged separately. Each also acts on an extension.

- **Inversion** *(per chromosome)* — reverse a segment in place, flipping the strand of every gene in it. The classic signed-permutation move: `+2 +3 +4` becomes `−4 −3 −2`.
- **Transposition** *(per chromosome)* — cut a segment out and reinsert it **elsewhere on the same chromosome**.
- **Translocation** *(per gene copy)* — move a segment to a **different chromosome** of the same genome (a no-op if the genome has only one chromosome).

A moved block — transposed or translocated — lands **inverted** with probability `inversion_probability` (default `0`, i.e. it keeps its orientation).

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

Transposition stays within a chromosome; translocation is the one that carries a gene lineage **across** to another chromosome — which is why it is counted per gene copy, like transfer, and why (unlike fission or fusion below) it is *not* an edge in the chromosome network: the chromosomes themselves are unchanged, only a gene has hopped between them.

## Chromosomes and the tier

A genome is seeded with a **karyotype**: `chromosomes=N` chromosomes, each with a `topology` — `"circular"` (the default) or `"linear"`, or a per-chromosome list like `["circular", "linear"]` for a mixed set. The founding `initial_families` genes are dealt round-robin across them. Topology is a label: it does not change how inversions or the tier behave.

On top of the karyotype, four events change the **number** of chromosomes — the *chromosome tier*:

- **`fission`** *(per chromosome)* — a chromosome splits in two.
- **`fusion`** *(per chromosome)* — two chromosomes of a genome merge into one.
- **`chromosome_origination`** *(per lineage)* — a de-novo replicon appears (an empty new chromosome — a plasmid).
- **`chromosome_loss`** *(per chromosome)* — a whole chromosome and all its genes die (its genes are recorded as losses). A lineage never loses its *last* chromosome this way.

```python
tree = species.simulate_species_tree(birth=1.0, death=0.1, n_extant=3, seed=42)
g = simulate_genomes_ordered(
    tree, duplication=0.15, loss=0.1, origination=0.25,
    chromosomes=2, fission=0.25, fusion=0.25,
    chromosome_origination=0.03, chromosome_loss=0.03,
    initial_families=5, seed=42)
```

## The chromosome network

Chromosomes carry a genuine **identity**. A chromosome id is re-minted at every event that reshapes it — a speciation, a fission, a fusion — and each of those edges is recorded, so the run leaves behind not just the chromosomes at the tips but the *genealogy* that connects them: the **chromosome network**. It is the middle tier of three that nest — the species tree contains the chromosome network, which contains the gene trees:

```
species tree  ⊃  chromosome network  ⊃  gene trees
```

It is a genuine **network**, not a tree, because of one event: **fusion joins two chromosome lineages into one** (two parents, one child — a reticulation). Fission and speciation are ordinary splits (one parent, two children); origination is a root; loss is a leaf. The whole thing is a directed graph, and — because a graph is just a graph — it is recorded the way graphs are, as an **edge list**: `chromosome_events`, one row per event.

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

You can read the karyotype's whole history off this table: two seed chromosomes; one dies early; the other descends through the species splits; then chromosome `4` fissions into `6` and `7`, which promptly fuse back into `8`. The genome at the leaf below has a single chromosome — `[ 0+ 5+ 2+ ]` — the survivor of that history. (This edge list is the network's ground truth; a portable graph file — GraphML, DOT — is a later convenience, and no `eNewick` string is involved: a multi-rooted, reticulating graph is not a tree, so it is not forced into Newick.)

## The events, and their older names

Readers from the rearrangement and comparative-genomics literature already have names for these events; as in every chapter, the names live in one table and organise nothing.

| What it does | Here | From the literature |
|---|---|---|
| reverse a gene segment, strands flipped | `inversion=…` | Inversion / reversal (GRIMM, MGR) |
| cut a segment, paste it elsewhere on the same chromosome | `transposition=…` | Transposition |
| a segment moves to another chromosome | `translocation=…` | Reciprocal translocation |
| a block of adjacent genes copied together | `duplication` + `duplication_extension=` | Segmental duplication |
| split one chromosome, or merge two | `fission=…`, `fusion=…` | Chromosome fission / fusion |
| a de-novo chromosome appears | `chromosome_origination=…` | Plasmid / replicon gain |
| the chromosome network as a whole | fission + fusion → `chromosome_events` | Karyotype / chromosome-number evolution |

## The `OrderedGenomesResult` object

`simulate_genomes_ordered` returns an **`OrderedGenomesResult`**, the ordered counterpart of `GenomesResult` — the same spine, with the structured extras:

- `.complete_tree` — the species tree the genomes ran on, extinct lineages included.
- `.genomes` — a dict from node id to that node's genome, now a tuple of **`Chromosome`** objects; each `Chromosome` has an `id`, a `topology`, and an ordered list of **`Gene`** objects (`id`, `family`, `strand`).
- `.events` — the gene-genealogy log, exactly as in Chapter 5 (origination, duplication, transfer, loss, speciation), from which `.gene_trees` and `.profiles` are derived unchanged. Position and orientation are *not* here — they live in the genomes and the two logs below.
- `.rearrangements` — the inversion / transposition / translocation log.
- `.chromosome_events` — the chromosome network, as the edge list above.
- `.gene_trees`, `.profiles`, `.seed` — as before.

with the methods `.family_counts(node_id)` (the multiset view), `.gene_order(node_id)` (the layout — `(chromosome, position, strand, family, gene id)` per gene), and `.write(dir, outputs=[...])`.

```python
g.genomes[2]                     # the chromosomes of node n2
g.gene_order(2)                  # its layout, gene by gene
g.chromosome_events              # the chromosome network (edge list)
g.gene_trees[0].to_newick()      # a family's gene tree — unchanged from unordered
```

## Usage from Python

```python
from zombi2 import species
from zombi2.genomes import simulate_genomes_ordered
from zombi2.rates.distributions import Geometric
from zombi2.rates import modifiers as mod

tree = species.simulate_species_tree(birth=1.0, death=0.2, n_extant=30, seed=1)

# ordered genome: the Chapter 5 events, plus inversions, on a single chromosome
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

# several chromosomes that split, merge, and gain a plasmid — the chromosome network
g = simulate_genomes_ordered(
    tree, chromosomes=6, topology="linear",
    fission=0.05, fusion=0.05, chromosome_origination=0.02, chromosome_loss=0.02,
    origination=0.4, initial_families=20, seed=1)

# rates can still depend on time (the skyline), as at every level
g = simulate_genomes_ordered(
    tree, inversion=1.0 * mod.OnTime({0: 1.0, 2: 0.2}), initial_families=10, seed=1)

# the outputs
g.genomes                             # every node's chromosomes
g.gene_order(next(iter(g.genomes)))   # a node's layout
g.rearrangements                      # the inversion/transposition/translocation log
g.chromosome_events                   # the chromosome network, as an edge list
g.gene_trees                          # one gene tree per family, as in Chapter 5
```

## Usage from the CLI

The ordered resolution is `--resolution ordered`. It adds the segmental and chromosome flags to the Chapter 5 events, each still a plain number:

```bash
# segmental duplications, losses and inversions on three chromosomes
zombi2 genomes -t out/species_complete.nwk --resolution ordered --duplication 0.2 --loss 0.2 --origination 0.5 --inversion 0.3 --chromosomes 3 --seed 1 -o out/

# blocks relocate and move between chromosomes (sometimes inverting), and chromosomes split and merge
zombi2 genomes -t out/species_complete.nwk --resolution ordered --origination 0.5 --transposition 0.2 --translocation 0.1 --inversion-probability 0.5 --fission 0.05 --fusion 0.05 --chromosomes 2 --seed 1 -o out/
```

The segment-length distributions (how long a duplicated or inverted run is) are structured objects, so like the rate modifiers they stay in Python; the CLI takes the rates as numbers. The same tree-reading rules apply as in Chapter 5, and `--write` selects among the richer ordered outputs below.

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

`gene_order.tsv` is the ordered genome's headline output — the signed gene order of every observed leaf:

```
species  chromosome  position  strand  family  gene
2        3           0         1       0       9
2        3           1         1       5       10
2        3           2         1       2       11
```

