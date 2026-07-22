# Genomes I: unordered

Genomes live inside the species tree, and they can be simulated at different levels of **resolution**. The simplest case is a single gene family evolving along the species tree. The most complex is a genome in which every nucleotide is tracked across several chromosomes. ZOMBI2 offers three resolutions, one per chapter: **unordered** here, **ordered** in Chapter 5, **nucleotide** in Chapter 6.

The **unordered** resolution is genomes made of gene families and nothing more: no position along a chromosome, no DNA sequence. Genes are copied, lost, born from nothing, and passed sideways between lineages.

## The four events

An unordered genome evolves by four kinds of event, applied to every lineage as it runs down the species tree:

- **Duplication** — a gene copy is copied, so its family gains a member in that lineage.
- **Transfer** — a copy jumps from one lineage to another that is alive at the same moment. This is the only event that crosses lineages, and it gets its own section below.
- **Loss** — a gene copy is deleted; a family that loses its last copy is gone from that lineage.
- **Origination** — a brand-new family appears in a lineage, with one copy.

You give ZOMBI2 a rate for each, and it plays the events out along the tree, starting from the root genome and letting speciation hand a lineage's genome down to both its children. Out comes the genome of *every* lineage in the tree together with the event log that produced it.

```python
from zombi2 import species
from zombi2.genomes import simulate_genomes_unordered

tree = species.simulate_species_tree(birth=1.0, death=0.3, n_extant=20, seed=1)
g = simulate_genomes_unordered(
    tree, duplication=0.2, loss=0.25, origination=0.5, initial_families=20, seed=1)
```

The root starts with `initial_families` families of one copy each, recorded as originations at the crown; from there the four rates drive everything.

The two ways of stocking a genome are worth separating. `initial_families` puts families there at the start, and `origination` adds new ones as the run goes. From Python `initial_families` defaults to 0, so a caller says what it wants; from the command line `--initial-families` defaults to 100, so a bare run hands back a genome rather than a hundred empty ones. `--origination` is 0 by default either way: nothing arrives that you did not ask for. Every run's `.log` records the value it used, so a run is never ambiguous about which.

## What the rate depends on

The rates follow the **same grammar as the species level** (`base` optionally wrapped in a scope, optionally multiplied by modifiers). The scope answers *per what*, and the default is the natural one for each event. Duplication, transfer, and loss are counted **per copy**: a family with ten copies is ten times as likely to duplicate or lose one as a family with a single copy, which is what you want — more genes, more chances. Origination is counted **per lineage**: acquiring a wholly new family is a property of the lineage, not of any gene it already has.

Rates can also depend on **time**. Multiplying a base rate by an `OnTime` modifier makes it change at set moments — the skyline, or episodic, genome: fast early and slow later, or any schedule you give.

```python
from zombi2.rates import modifiers as mod
# lots of new families early, then origination shuts off after time 2
g = simulate_genomes_unordered(tree, origination=1.0 * mod.OnTime({0: 1.0, 2: 0.0}), seed=1)
```

## Lateral gene transfers

Transfer is the one event that couples lineages, and it is what makes the unordered resolution more than four independent birth–death processes. When a transfer fires, a copy is picked from the whole pool of live genes, and it is delivered to another lineage that is **alive at that same instant**.

Three arguments shape what a transfer does:

- **`transfer_to`** — who receives. `"uniform"` (the default) picks any other contemporaneous lineage with equal chance; `"distance"` makes closer relatives likelier, weighting recipients by how far they sit from the donor on the tree. The distance version is *scale-free*: its strength means the same whether your tree is measured in years or in millions of them. A third rule weights recipients by a trait instead, so that only competent lineages take DNA up; it is a coupling, and Chapter 9 covers it.
- **`replacement`** — what happens on arrival. By default the incoming copy is **additive**: the recipient simply gains a copy. With `replacement=True` it **overwrites** a copy of the same family already present, and falls back to additive when the recipient has none.
- **`self_transfer`** — whether a lineage may donate to itself. Off by default. With additive arrival the lineage gains a copy, so the gene content changes as it would under a duplication, but the event is recorded as a transfer. Combined with `replacement=True` it is not a duplication at all: the arriving copy overwrites a paralogue in the same genome, which is gene conversion.

```python
tree = species.simulate_species_tree(birth=1.0, death=0.4, n_extant=30, seed=7)
# horizontal transfer biased toward close relatives, overwriting resident copies
g = simulate_genomes_unordered(
    tree, transfer=0.5, transfer_to="distance", replacement=True,
    origination=0.4, initial_families=10, seed=3)
```

One consequence is worth stating plainly: a transfer can arrive **from a lineage that later goes extinct**. A genome run happens on the complete tree, dead branches included, so a gene can enter a survivor from a donor that leaves no other trace.

## The `GenomesResult` object

`simulate_genomes_unordered` returns a **GenomesResult** which carries:

- `.complete_tree` — the species tree the genomes ran on, extinct lineages and all.
- `.genomes` — a dict from node to that node's genome.
- `.events` — the event log: every gene event with its time and lineage — origination, duplication, transfer, loss, and the *speciations* at a split.
- `.profiles` — the family × extant-species copy-count table.
- `.gene_trees` — one `GeneTree` per family.
- `.seed` — the seed, so the run reproduces.

and two methods:

- `.family_counts(node_id)` — a `Counter` collapsing a node's genome to `family → number of copies`, when you want the multiset rather than the individual copies.
- `.write(dir)` — materialise the outputs to disk: the event log (`genome_events.tsv`) and the profiles (`profiles.tsv`).

```python
n5 = g.genomes[5]                    # the gene copies in node n5
counts = g.family_counts(5)          # {family: copies} for the same node
g.write("out/")                      # genome_events.tsv + profiles.tsv
```

## Profiles and gene trees

Two products are usually what you came for, and ZOMBI2 derives both from the run's recorded history.

**Profiles** are the classic comparative-genomics view: how many copies of each gene family sit in each extant species. They are read straight off the observed genomes, so the run stays lean and you materialise them on access.

```python
g.profiles.matrix        # families × extant-species copy counts, a NumPy array
g.profiles.presence      # the same as 0/1 presence/absence
g.profiles.to_tsv()      # the table as text
```

**Gene trees** are the deeper output. Every family has two trees: the `.complete` tree with every gene lineage, and the `.extant` tree pruned to the genes that survive.

```python
gt = g.gene_trees[7]                 # the gene tree of family 7
gt.to_newick("extant")               # the surviving copies as Newick ...
gt.to_newick("complete")             # ... or the whole genealogy
gt.origination                       # when the family began
```

The root of a gene tree carries a branch length, as the species tree's does. A family starts at its origination and the founding gene then lives for a while before its first duplication, transfer or speciation, and that wait is the root's branch: the family's stem. It matters most for a family that originates late, whose whole history is short. A gene that originated and never split at all is a one-node tree, written as its own lifespan — `g55:0.263097;`.

## Usage from Python

The whole range is one function call:

```python
from zombi2 import species
from zombi2.rates import modifiers as mod
from zombi2.genomes import simulate_genomes_unordered

tree = species.simulate_species_tree(birth=1.0, death=0.3, n_extant=30, seed=1)

# a plain duplication–loss–origination run
g = simulate_genomes_unordered(
    tree, duplication=0.2, loss=0.25, origination=0.5, initial_families=20, seed=1)

# origination only — every family stays a single copy, none is ever duplicated
g = simulate_genomes_unordered(tree, origination=0.6, seed=1)

# a skyline: new families pour in early, then origination shuts off after time 2
g = simulate_genomes_unordered(tree, origination=1.0 * mod.OnTime({0: 1.0, 2: 0.0}), seed=1)

# horizontal transfer, biased toward close relatives, overwriting resident copies
g = simulate_genomes_unordered(
    tree, transfer=0.5, transfer_to="distance", replacement=True,
    origination=0.4, initial_families=10, seed=3)

# the genomes you observe are the extant tips
observed = {n.id: g.genomes[n.id] for n in g.complete_tree.extant()}

# and the outputs, derived from that history
g.profiles.matrix                                # family × extant-species copy counts
some_family = next(iter(g.gene_trees))
g.gene_trees[some_family].to_newick("extant")    # that family's surviving gene tree
g.write("out/")                                  # the event log + profiles, on disk
```

## Usage from the CLI

`zombi2 genomes` evolves gene families along a species tree read from a Newick file. The unordered resolution is the default, and each rate is a plain number:

```bash
# duplication–loss–origination along a species tree
zombi2 genomes -t out/species/species_complete.nwk --duplication 0.2 --loss 0.25 --origination 0.5 --seed 1 -o out/

# horizontal transfer biased toward close relatives, overwriting resident copies
zombi2 genomes -t out/species/species_complete.nwk --transfer 0.5 --transfer-to distance --replacement --origination 0.4 --seed 3 -o out/
```

## Outputs

By default a run writes the event log and the profiles:

```
out/genome_events.tsv    the gene genealogy (the source of truth)
out/profiles.tsv         family × extant-species copy counts
```

Two more are written on request, because both are large and most runs do not need them:

```
out/genomes.tsv                    every node's genes, ancestors included
out/gene_tree_fam<f>_*.nwk         each family's genealogy, complete and extant
```

`genomes.tsv` is one row per gene copy, so a lineage carrying six genes has six rows:

```
lineage  family  copy
n0       0       0
n0       1       1
n0       2       2
n1       0       3
n1       0       8
```

Read a block of rows sharing a `lineage` and you have that lineage's genome. `family` says which gene family a copy belongs to, so the two `n1` rows above are two copies of family `0` — one of them a duplicate. `copy` is the individual gene, and it is the same identifier the event log uses, so any gene here can be traced back to the event that made it. `profiles.tsv` is this same information counted rather than listed, and only for the extant tips; `genomes.tsv` keeps the ancestors, the root included, which is what you want if you are scoring a reconstruction of ancestral gene content.

`--write gene_trees` writes a Newick per family.

`--write` picks which of the four to keep, and a `names.tsv` joins the `n<id>` columns back to your labels whenever you brought your own tree. Everything is derived from the event log, so writing none of them loses nothing that cannot be replayed.

They land in `out/genomes/`, and the gene trees — one file per family per view — in `out/genomes/gene_trees/` below it, so a run of a few hundred families does not bury the rest. `--flat` writes everything straight into `out/` instead. The full list of files lives in Appendix B.
