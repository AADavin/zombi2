# A first simulation

This chapter walks you through your first complete ZOMBI2 run, start to finish. You will build
a species tree, evolve gene families along it, look at everything the simulator writes to disk,
and finally load two of those results back into Python. By the end you will have a full
simulated dataset in a folder called `out/` and a clear picture of what each file means.

This is the simplest use of ZOMBI2: a *hierarchical* run over the first two levels from the
previous chapter — a species tree, then gene families evolved along it. Each stage is seeded, so
the same `seed` always gives byte-identical results, and your run will match the numbers you see
here.

## Simulating a species tree

Start with the species tree. We condition a birth–death process on 20 extant tips and a crown
age of 5.

```python
import zombi2 as z

tree = z.simulate_species_tree(z.BirthDeath(birth=1.0, death=0.3),
                               n_tips=20, age=5.0, seed=1)

print(tree.to_newick())   # timed Newick
```

The same run from the command line writes the tree straight to a folder:

```bash
zombi2 species --birth 1 --death 0.3 --tips 20 --age 5 --seed 1 -o out/
```

The result is a timed tree with 20 surviving tips. `BirthDeath` is the model; `birth` and
`death` are the per-lineage speciation and extinction rates.

![A simulated species tree with the duplication, loss and transfer events of one gene family marked on its branches.](figures/species_tree_events.pdf){width=100%}

## Evolving gene families along the tree

Now evolve gene families forward along that tree. Each family starts as a single copy at its
point of origination and then duplicates, transfers between branches, and is lost, at the rates
you supply.

```python
genomes = z.simulate_genomes(tree, duplication=0.2, transfer=0.1, loss=0.25,
                             origination=0.5, initial_families=40, seed=42)

print(genomes.profiles.matrix.shape)   # (n_families, n_species) copy numbers
```

Equivalently, from the command line, using the tree you just wrote:

```bash
zombi2 genomes --tree out/species_tree.nwk --dup 0.2 --trans 0.1 \
               --loss 0.25 --orig 0.5 -o out/
```

Here `initial_families=40` seeds the root genome with 40 families; `origination=0.5` is the rate at
which new families appear along branches. The four core rates — duplication, transfer, loss,
origination — drive the whole forward process. The built-in model runs on the Rust engine
automatically.

## The output folder

Write everything to disk with a single call:

```python
genomes.write("out/")
```

That single call writes the full layout below. On the command line the default `zombi2 genomes`
command writes only the common subset (`Profiles.tsv`, `Presence.tsv` and `gene_trees/`); add
`--write all` for everything shown here — Chapter 9 covers `--write` in full. Here is what you
will find:

| File / folder | Contents |
|---|---|
| `species_tree.nwk` | the timed species tree (Newick) |
| `species_nodes.tsv` | node name, time, `is_leaf`, `is_extant` |
| `gene_family_events/<fid>_events.tsv` | per-family event log (O/D/T/S/L with lineage ids) |
| `gene_trees/<fid>_complete.nwk`, `_extant.nwk` | reconstructed gene trees (with / without losses) |
| `Transfers.tsv` | one row per transfer: time, family, donor, recipient, ids |
| `Gene_family_summary.tsv` | per family: origin, event counts, extant copies, species present |
| `Profiles.tsv`, `Presence.tsv` | families × species copy-number and presence matrices |

Two things are worth understanding before you open these files.

The **event log** is a full genealogy. Every event re-mints the gene-lineage ids it touches, and
speciations are logged too, so each family's `_events.tsv` records a complete parent-to-children
history. Reconstructing a gene tree from it is then pure post-processing.

The **gene trees** come in two flavours. The *complete* tree keeps every lineage, including
losses (leaves labelled `LOSS_<id>`). The *extant* tree is pruned to lineages with a surviving
copy, with unifurcations suppressed; its leaves are labelled `<species>_<geneid>`. The number of
leaves in an extant tree equals the family's total copy count across all species.

![A gene family's *complete* gene tree: every copy implied by its duplication, transfer, and loss events — including a copy that was later lost (the dashed lineage ending in a cross). The *extant* tree is this tree pruned to the surviving copies.](figures/gene_tree.pdf){width=100%}

## Reading the results back in Python

### A gene tree

`gene_trees()` returns, for each family id, a pair of Newick strings — complete and extant.

```python
trees = genomes.gene_trees()          # {fid: (complete, extant) newick}
complete, extant = trees["7"]
```

`extant` is `None` when a family has no survivors. Otherwise it is the tree you would reconcile
against the species tree, implied by that family's duplication, transfer and loss events.

### The profile matrix

The profile matrix is the key object for phylogenetic-profiling analyses: gene families as rows,
extant species as columns.

```python
P = genomes.profiles
P.matrix        # integer copy numbers
P.presence()    # binary presence / absence
P.families      # row labels
P.species       # column labels
P.to_tsv()      # tab-separated text (as written to Profiles.tsv)
```

`P.matrix` is exactly what was written to `Profiles.tsv`, and `P.presence()` matches
`Presence.tsv`. Each entry counts how many copies of a family survive in a species.

If you want the raw events instead of reconstructed trees, iterate the event log directly:

```python
from zombi2 import EventType

for r in genomes.event_log:
    if r.event is EventType.TRANSFER:
        print(r.time, r.donor, "->", r.recipient, [op.gid for op in r.genes])
```

Each `EventRecord` carries a `time`, an `event` (an `EventType`), a `branch`, optional `donor`
and `recipient`, and `genes` — the per-gene rows, the first being the source lineage and the
rest its descendants.

::: note
`Transfers.tsv` gives you the same transfer events in flat tabular form, one row each, if you
would rather not walk the log in code.
:::

You now have a complete simulated dataset: a timed species tree, gene families evolved along it,
reconstructed gene trees, and the profile matrix — all reproducible from a pair of seeds. The
following chapters unpack each stage: species-tree models and conditioning, shared versus
per-family rates, transfers, and bounding genome growth — before the later parts add the
remaining levels (traits and sequences) and the joint coevolution models.
