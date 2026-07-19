# Genomes I: unordered

The species tree is the backbone; this chapter puts genes on it. At the **unordered** level a lineage carries a *multiset of gene families* — how many copies of each family it has, and nothing more: no position along a chromosome, no DNA sequence. That comes in the next chapter. Here a genome is a bag of families, and it changes as the tree grows: genes are copied, lost, born from nothing, and passed sideways between lineages. This is the first and simplest of the three genome resolutions, and the one the other two are built on.

## The four events

An unordered genome evolves by four kinds of event, applied to every lineage as it runs down the species tree:

- **Duplication** — a gene copy is copied, so its family gains a member in that lineage.
- **Loss** — a gene copy is deleted; a family that loses its last copy is gone from that lineage.
- **Origination** — a brand-new family appears in a lineage, with one copy.
- **Transfer** — a copy jumps from one lineage to another that is alive at the same moment. This is the only event that crosses lineages, and it gets its own section below.

You give ZOMBI2 a rate for each, and it plays the events out along the tree, starting from the root genome and letting speciation hand a lineage's genome down to both its children. Out comes the genome of *every* lineage in the tree together with the event log that produced it.

```python
from zombi2 import species
from zombi2.genomes import simulate_genomes_unordered

tree = species.simulate_species_tree(birth=1.0, death=0.3, n_extant=20, seed=1)
g = simulate_genomes_unordered(
    tree, duplication=0.2, loss=0.25, origination=0.5, initial_families=20, seed=1)
```

The root starts with `initial_families` families of one copy each, recorded as originations at the crown; from there the four rates drive everything.

## What the rate depends on

The rates follow the **same grammar as the species level** (`base` optionally wrapped in a scope, optionally multiplied by modifiers), so nothing new has to be learned — only what the defaults mean here.

The scope answers *per what*, and the default is the natural one for each event. Duplication, transfer, and loss are counted **per copy**: a family with ten copies is ten times as likely to duplicate or lose one as a family with a single copy, which is what you want — more genes, more chances. Origination is counted **per lineage**: acquiring a wholly new family is a property of the lineage, not of any gene it already has.

Rates can also depend on **time**. Multiplying a base rate by a `OnTime` modifier makes it change at set moments — the skyline, or episodic, genome, fast early and slow later, or any schedule you give:

```python
from zombi2.rates import modifiers as mod
# lots of new families early, then origination shuts off after time 2
g = simulate_genomes_unordered(tree, origination=1.0 * mod.OnTime({0: 1.0, 2: 0.0}), seed=1)
```

Two richer dials are part of the design but land in a later release, and this version rejects them loudly rather than pretending: **per-family heterogeneity** — letting each family carry its own rate, or a single per-family *speed* that scales all of a family's rates together, so some families churn and others sit still — and **scope overrides** on the genome rates. For now every family shares the given rate, and the defaults (per copy, per lineage) are fixed.

## Transfer — the sideways event

Transfer is the one event that couples lineages, and it is what makes the unordered level more than four independent birth–death processes. When a transfer fires, a copy is picked from the whole pool of live genes, and it is delivered to another lineage that is **alive at that same instant**. Because it links contemporaries, the engine cannot walk one branch at a time the way the other three events allow; it runs a single timeline over the whole tree at once, exactly as the species level grows all living lineages together.

Three arguments shape what a transfer does:

- **`transfer_to`** — who receives. `"uniform"` (the default) picks any other contemporaneous lineage with equal chance; `"distance"` makes closer relatives likelier, weighting recipients by how far they sit from the donor on the tree. The distance version is *scale-free* — its strength means the same whether your tree is measured in years or in millions of them — and `transfer_to="distance"` is shorthand for `Distance(decay=1.0)`, whose `decay` you can turn up to concentrate transfers among near neighbours.
- **`replacement`** — what happens on arrival. By default the incoming copy is **additive**: the recipient simply gains a copy. With `replacement=True` it **overwrites** a copy of the same family already present, and falls back to additive when the recipient has none — the acquired gene displacing its resident homolog.
- **`self_transfer`** — whether a lineage may donate to itself. Off by default; on, a lineage can be its own recipient.

```python
tree = species.simulate_species_tree(birth=1.0, death=0.4, n_extant=30, seed=7)
# horizontal transfer biased toward close relatives, overwriting resident copies
g = simulate_genomes_unordered(
    tree, transfer=0.5, transfer_to="distance", replacement=True,
    origination=0.4, initial_families=10, seed=3)
```

One consequence is worth stating plainly: a transfer can arrive **from a lineage that later goes extinct**. Genomes evolve on the whole tree, dead branches included (next section), so a gene can enter a survivor from a donor that leaves no other trace — the sideways echo of a lost lineage.

## The complete tree

Just as a birth–death tree is really two trees, a genome run happens on the **complete** tree — every lineage that ever lived, including the ones that went extinct — not only on the survivors. This is deliberate: it is what lets a transfer come from the dead, and it is what makes the true gene-family history complete rather than a trace of it. What you *observe* is the genomes at the **extant tips**; the rest are the hidden history behind them.

```python
# the genomes you actually observe: the extant tips of the complete tree
observed = {n.id: g.genomes[n.id] for n in g.complete_tree.extant()}
```

## The `GenomesResult` object

`simulate_genomes_unordered` returns a **`GenomesResult`**, not a bare genome — a run produces the genome of every node *and* the event log that generated them, and both matter. It shares the shape of every level's result (`SpeciesResult`, `SequencesResult`, `TraitsResult`), so the four levels stay symmetric.

A `GenomesResult` carries:

- `.complete_tree` — the species tree the genomes ran on, extinct lineages and all.
- `.genomes` — a dict from node to that node's genome, a tuple of `GeneCopy` objects. The key is the node's integer id, the one that prints as `n<id>` in the Newick (node `5` is `n5`); each `GeneCopy` knows its own `id` and its `family`.
- `.events` — the event log: every gene event with its time and lineage — origination, duplication, transfer, loss, and the *speciations* that re-id each gene at a split. Gene ids are **per branch** (the ZOMBI1 model: every event ends a gene and starts fresh ids for its descendants), so this is the complete source of truth the profiles and gene trees are derived from.
- `.profiles` — the family × extant-species copy-count table (the next section).
- `.gene_trees` — one `GeneTree` per family (the next section).
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

**Profiles** are the classic comparative-genomics view: how many copies of each gene family sit in each extant species — families down the rows, species across the columns, zero where a family is absent. They are read straight off the observed genomes, so the run stays lean and you materialise them on access.

```python
g.profiles.matrix        # families × extant-species copy counts, a NumPy array
g.profiles.presence      # the same as 0/1 presence/absence
g.profiles.to_tsv()      # the table as text
```

**Gene trees** are the deeper output. Every family has its own gene tree — the true genealogy of its copies, growing *inside* the complete species tree. As in ZOMBI1, a gene lives on a single species branch and carries that branch's id; every event — a duplication, a transfer, or a speciation — ends it and starts fresh ids for its descendants, until it finally ends by being lost, by its species going extinct, or at an extant tip. So every tip is a distinct gene. `.gene_trees` gives one `GeneTree` per family, and — exactly like the species result — each carries **two trees**: the `.complete` tree with every gene lineage, and the `.extant` tree pruned to the genes that survive. Every node knows the species branch it sits on and the event that ended it.

```python
gt = g.gene_trees[7]                 # the gene tree of family 7
gt.to_newick("extant")               # the surviving copies as Newick ...
gt.to_newick("complete")             # ... or the whole genealogy
```

This is a recorded *history*, not a reconciliation: reconciliation is the inference problem of fitting an *observed* gene tree back onto a species tree, and it belongs to the analysis tools — here we simulated the embedding, so we simply write it down.

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

*(A `zombi2 genomes` command mirroring these calls arrives with the genome CLI; today this is the Python engine, and `.write` already puts the outputs on disk.)*
