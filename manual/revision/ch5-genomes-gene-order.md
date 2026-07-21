# Genomes II: gene order

The previous chapter put genes on the tree as a *bag of families* — how many copies of each, and nothing more. This chapter gives them **structure**: an order along a chromosome and an orientation, so a genome becomes a signed sequence of genes rather than a multiset. This is the **ordered** resolution, and it keeps everything from Chapter 4 while adding *position*.

Everything below layers on the unordered core. The four events, the rate grammar, transfers, the event log, the gene trees, the profiles — all still here, all working the same way. Position is orthogonal to genealogy: two genes being neighbours says nothing about how they are related, and the gene trees and profiles a run produces are exactly the ones the unordered engine would produce for the same events. What is new is *where each gene sits* and *which way it points*, and the events that rearrange that.

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

Each gene is written as its family with the strand as `+` / `−` — the strand is the integer `+1` / `−1`, and in these listings `+` means `+1`. This leaf has one chromosome of seven genes; family `2` sits in a run of three tandem copies, and the whole thing is one linear order. The gene tree of a family is unchanged from Chapter 4 — the true genealogy, read off the same event log:

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

Transposition stays within a chromosome; translocation is the one that carries a gene lineage **across** to another chromosome. That is why it is counted per gene copy, like transfer, and why it does not touch the chromosomes themselves — only a gene has hopped between them. The events that change the chromosomes are the subject of the next chapter.

## The events, and their older names

Readers from the rearrangement and comparative-genomics literature already have names for these events; as in every chapter, the names live in one table and organise nothing.

| What it does | Here | From the literature |
|---|---|---|
| reverse a gene segment, strands flipped | `inversion=…` | Inversion / reversal (GRIMM, MGR) |
| cut a segment, paste it elsewhere on the same chromosome | `transposition=…` | Transposition |
| a segment moves to another chromosome | `translocation=…` | Reciprocal translocation |
| a block of adjacent genes copied together | `duplication` + `duplication_extension=` | Segmental duplication |

## The `OrderedGenomesResult` object

`simulate_genomes_ordered` returns an **`OrderedGenomesResult`**, the ordered counterpart of `GenomesResult` — the same spine, with the structured extras:

- `.complete_tree` — the species tree the genomes ran on, extinct lineages included.
- `.genomes` — a dict from node id to that node's genome, now a tuple of **`Chromosome`** objects; each `Chromosome` has an `id`, a `topology`, and an ordered list of **`Gene`** objects (`id`, `family`, `strand`).
- `.events` — the gene-genealogy log, exactly as in Chapter 4 (origination, duplication, transfer, loss, speciation), from which `.gene_trees` and `.profiles` are derived unchanged. Position and orientation are *not* here — they live in the genomes and in the logs.
- `.rearrangements` — the inversion / transposition / translocation log.
- `.chromosome_events` — the chromosome network, covered in the next chapter.
- `.gene_trees`, `.profiles`, `.seed` — as before.

with the methods `.family_counts(node_id)` (the multiset view), `.gene_order(node_id)` (the layout — `(chromosome, position, strand, family, gene id)` per gene), and `.write(dir, outputs=[...])`.

```python
g.genomes[2]                     # the chromosomes of node n2
g.gene_order(2)                  # its layout, gene by gene
g.gene_trees[0].to_newick()      # a family's gene tree — unchanged from unordered
```

## Usage from Python

```python
from zombi2 import species
from zombi2.genomes import simulate_genomes_ordered
from zombi2.rates.distributions import Geometric
from zombi2.rates import modifiers as mod

tree = species.simulate_species_tree(birth=1.0, death=0.2, n_extant=30, seed=1)

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
g.rearrangements                      # the inversion/transposition/translocation log
g.gene_trees                          # one gene tree per family, as in Chapter 4
```

## Usage from the CLI

The ordered resolution is `--resolution ordered`. It adds the segmental flags to the Chapter 4 events, each still a plain number:

```bash
# segmental duplications, losses and inversions on three chromosomes
zombi2 genomes -t out/species_complete.nwk --resolution ordered --duplication 0.2 --loss 0.2 --origination 0.5 --inversion 0.3 --chromosomes 3 --seed 1 -o out/

# blocks relocate and move between chromosomes, sometimes inverting
zombi2 genomes -t out/species_complete.nwk --resolution ordered --origination 0.5 --transposition 0.2 --translocation 0.1 --inversion-probability 0.5 --chromosomes 2 --seed 1 -o out/
```

The segment-length distributions (how long a duplicated or inverted run is) are structured objects, so like the rate modifiers they stay in Python; the CLI takes the rates as numbers. The same tree-reading rules apply as in Chapter 4, and `--write` selects among the richer ordered outputs.

## Outputs

`.write(dir, outputs=[...])` materialises the chosen products to disk:

```python
g.write("out/", outputs=("events", "profiles", "gene_order", "rearrangements"))
```

```
out/genome_events.tsv        the gene genealogy (the source of truth)
out/profiles.tsv             family × extant-species copy counts
out/gene_order.tsv           every node's layout, one row per gene
out/rearrangements.tsv       inversions, transpositions, translocations
```

`gene_order.tsv` is the ordered genome's headline output: the signed gene order of every node, one row per gene. Ancestors are included, not just the observed leaves, so node 0 below is the root and node 1 an internal branch:

```
species  chromosome  position  strand  family  gene
0        0           0         1       1       1
0        0           1         -1      2       5
1        1           0         1       1       7
1        1           1         -1      2       8
2        2           0         1       1       11
```

Ancestral rows are what make the rearrangement log usable. `rearrangements.tsv` gives each inversion a start and a length on a branch; to check what it did, or to replay it, you need the genome the branch started from — that is its parent's rows here. Without them the log can only be read at the tips.

### Replaying a run

`genome_events.tsv` records which gene copy each event created or ended, but not where on the chromosome it happened. That is deliberate: an event is about identity and descent, which is the same at every resolution, so the log is shared with the unordered core (Chapter 5) unchanged. The `event_positions` output adds the coordinates alongside it:

```python
g.write("out/", outputs=("gene_order", "rearrangements", "event_positions"))
```

```
time   kind                lineage  chromosome  start  length  family  donor  recipient  dest_position
0.0    origination         0        0           0      1       0
0.209  transfer_donor      2        2           0      1               2      1
0.209  transfer_recipient  1        1           0      1               2      1
0.345  loss                4        4           3      1
0.644  duplication         4        4           3      1                                 4
```

Every row belongs to one branch. `lineage` names it, and `chromosome`, `start` and `length` are coordinates in that branch's genome as it stood just before the event. So you can pull out the rows for a single branch and know everything that happened to it. `dest_position` says where a duplication's copy block landed. Origination carries its `family`, because it is the only event whose gene does not come from a genome you already have.

A transfer spans two branches, so it writes two rows — one on each — and both name the whole edge in `donor` and `recipient`. The `transfer_donor` row says what left; the `transfer_recipient` row says where it arrived. Pair them on time, donor and recipient.

Those three files together are enough to reconstruct the whole run: start from a node's parent in `gene_order.tsv`, apply that branch's rows from `rearrangements.tsv` and `genome_event_positions.tsv` in time order, and you get the node's own rows back. Rows sharing a timestamp apply in the order written.

The full list of files lives in Appendix B.
