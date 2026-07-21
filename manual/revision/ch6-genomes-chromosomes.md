# Genomes III: chromosomes

The previous chapter moved genes around *within* chromosomes. This one changes the chromosomes themselves: how many there are, how they split and merge, and how they are related to one another. It is the same `simulate_genomes_ordered` run, with a second set of events layered on top.

The distinction is worth keeping straight. A translocation moves a gene from one chromosome to another, but the chromosomes are unchanged by it. The events here do change them — a chromosome becomes two, two become one, a new one appears, one dies.

## Chromosomes and the tier

A genome is seeded with a **karyotype**: `chromosomes=N` chromosomes, each with a `topology` — `"circular"` (the default) or `"linear"`, or a per-chromosome list like `["circular", "linear"]` for a mixed set. The founding `initial_families` genes are dealt round-robin across them. Topology is a label: it does not change how inversions or the tier behave.

On top of the karyotype, four events change the **number** of chromosomes — the *chromosome tier*:

- **`fission`** *(per chromosome)* — a chromosome splits in two.
- **`fusion`** *(per chromosome)* — two chromosomes of a genome merge into one.
- **`chromosome_origination`** *(per lineage)* — a de-novo replicon appears (an empty new chromosome — a plasmid).
- **`chromosome_loss`** *(per chromosome)* — a whole chromosome and all its genes die (its genes are recorded as losses). A lineage never loses its *last* chromosome this way.

```python
from zombi2 import species
from zombi2.genomes import simulate_genomes_ordered

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

You can read the karyotype's whole history off this table: two seed chromosomes; one dies early; the other descends through the species splits; then chromosome `4` fissions into `6` and `7`, which promptly fuse back into `8`. The genome at the leaf below has a single chromosome — `[ 0+ 5+ 2+ ]` — the survivor of that history. No `eNewick` string is involved: a multi-rooted, reticulating graph is not a tree, so it is not forced into Newick.

## The events, and their older names

| What it does | Here | From the literature |
|---|---|---|
| split one chromosome, or merge two | `fission=…`, `fusion=…` | Chromosome fission / fusion |
| a de-novo chromosome appears | `chromosome_origination=…` | Plasmid / replicon gain |
| the chromosome network as a whole | fission + fusion → `chromosome_events` | Karyotype / chromosome-number evolution |

## Usage from Python

The result object is the `OrderedGenomesResult` of the previous chapter; `.chromosome_events` is the accessor this one adds.

```python
from zombi2 import species
from zombi2.genomes import simulate_genomes_ordered

tree = species.simulate_species_tree(birth=1.0, death=0.2, n_extant=30, seed=1)

# several chromosomes that split, merge, and gain a plasmid
g = simulate_genomes_ordered(
    tree, chromosomes=6, topology="linear",
    fission=0.05, fusion=0.05, chromosome_origination=0.02, chromosome_loss=0.02,
    origination=0.4, initial_families=20, seed=1)

g.chromosome_events              # the chromosome network, as an edge list
g.genomes[2]                     # the chromosomes of node n2
```

## Usage from the CLI

The chromosome flags join the ordered resolution's others:

```bash
# chromosomes split and merge along the tree
zombi2 genomes -t out/species_complete.nwk --resolution ordered \
    --origination 0.5 --fission 0.05 --fusion 0.05 --chromosomes 2 --seed 1 -o out/
```

## Outputs

The chromosome tier adds one file to those of the previous chapter:

```python
g.write("out/", outputs=("events", "profiles", "gene_order", "chromosome_events"))
```

```
out/chromosome_events.tsv    the chromosome network (edge list)
```

Its columns are the edge list above — `time · kind · lineage · parents · children` — one row per event, which is the network's ground truth. The full list of files lives in Appendix B.
