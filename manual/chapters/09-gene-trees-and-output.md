# Gene trees and output

A gene-family simulation produces one artefact that everything else is derived from: the *event
log*. This chapter treats that log as the ground truth it is, shows how the complete and extant
gene trees are read back out of it, documents the profile matrix and its sparse form, and specifies
every file the output folder can contain — including the compact event trace that lets large runs
skip the per-family reconstruction entirely.

## The event log is a full genealogy

Every event re-mints the gene-lineage ids it touches, and speciations are recorded alongside the
origination, duplication, transfer and loss events, so each family's slice of the log is a complete
parent-to-children genealogy. Nothing is inferred after the fact: the simulator writes the true
species branch of every event, so gene-tree reconstruction and reconciliation are pure
post-processing (as in ZOMBI 1), not likelihood or parsimony inference.

You can walk the log directly. Each record is an `EventRecord` carrying `time`, `event` (an
`EventType`), `branch`, optional `donor` and `recipient`, and `genes` — the per-gene rows, the first
of which is the source lineage and the rest its descendants. Each gene row exposes `.gid` (the
gene-lineage id) and `.role`.

```python
from zombi2 import EventType

for r in genomes.event_log:
    if r.event is EventType.TRANSFER:
        print(r.time, r.donor, "->", r.recipient, [op.gid for op in r.genes])
```

`genomes.event_log.by_family()` (also reachable as `genomes.gene_families`) groups the flat log into
`{family_id: [EventRecord, ...]}`, the per-family lists that every reconstruction consumes.

## Reconstructing the gene trees

`gene_trees()` turns the per-family records into Newick strings, one pair per family:

```python
trees = genomes.gene_trees()          # {family_id: (complete, extant)}
complete, extant = trees["7"]
```

- **complete** — every lineage the family ever had, including the ones that were lost. A lost tip is
  labelled `LOSS_<gid>`.
- **extant** — pruned to the lineages with a surviving copy, with the resulting degree-two nodes
  suppressed. An extant tip is labelled `<species>_<gid>`. `None` when the family has no survivors.

The number of leaves in the extant tree equals the family's total copy count across all species —
the same number you get by summing that family's row of the profile matrix. Passing
`annotate_species=True` additionally labels each internal gene node `<gid>|<species-branch>`, which
is what lets the tree be drawn against the species tree.

![The two trees `gene_trees()` returns for one family, side by side. **A**, the *complete* tree: every copy the family's duplication, transfer and loss events imply — including a lineage that was later lost (dashed, ending in a cross). **B**, the *extant* tree: the same tree pruned to the copies that survive to the present, with the resulting degree-two nodes suppressed. The loss is exactly the difference between the two. (The fully annotated *reconciliation* — either tree tagged with its species mapping — is the separate `reconciliations()` output described below.)](figures/gene_tree_panels.pdf){width=100%}

A fully annotated reconciliation is available separately. `reconciliations()` returns
`{family: Reconciliation(complete, extant, events)}`, where `complete` and `extant` are the two
trees annotated with their species mapping and `events` is the S/D/T/L event list. Reconciliation
here is exact annotation, not LCA inference, because the generating species branch of every event is
already known.

## The profile matrix

The `ProfileMatrix` is the key object for phylogenetic-profiling analyses: gene families (rows) by
extant species (columns), read straight off the final genome at each leaf. It is exposed as
`genomes.profiles`.

```python
P = genomes.profiles
P.families        # row labels (family ids)
P.species         # column labels (species names)
P.matrix          # dense integer copy-number array
P.presence()      # dense binary presence/absence (int8)
P.to_tsv()        # dense TSV text (to_tsv(presence=True) for the 0/1 form)
```

A copy-number profile is overwhelmingly zero — most families live in a handful of species — so the
dense `families × species` array is $O(N^2)$ in tip count and becomes the memory wall long before
the simulation itself does. `ProfileMatrix` therefore stores its data in **COO** form: three
parallel arrays holding, one entry per *present* cell, the family index, the species index and the
copy number. This is $O(\text{non-zeros})$, roughly $O(N)$.

The dense `matrix` and `presence()` views densify only when you touch them. The summaries below
never do — they run straight off the COO arrays and stay cheap at a million tips:

```python
P.nnz                     # number of present (non-zero) cells
P.coo                     # the raw (rows, cols, data) arrays
P.presence_per_family()   # species count per family
P.copies_per_family()     # total copies per family
P.copies_per_species()    # genome size per species
```

For on-disk use, `to_coo_tsv()` writes the sparse long table — one `family<TAB>species<TAB>copies`
row per present cell — with two leading `#species` and `#families` header lines that record the full
label sets in order, so the round trip through `from_coo_tsv()` is lossless even for all-absent rows
or columns. `from_tsv()` / `from_coo_tsv()` load either format back, and `from_leaf_genomes()` builds
a matrix directly from a `{node: Genome}` map without ever allocating the dense array.

## The output folder

`genomes.write("out/")` writes the ZOMBI-1-style output folder; the CLI's `--write` flag (whose
parts map onto the same `include=` set) selects the same components by name. `species_tree.nwk` and
`species_nodes.tsv` are always written; every other component is opt-in, and omitted components do
no work — in particular `trees` is what drives the expensive gene-tree reconstruction.

```python
genomes.write("out/", include=["profiles", "trees", "transfers"])
```

```bash
zombi2 genomes -t species_tree.nwk --dup 0.1 --trans 0.05 --loss 0.1 \
    --write profiles trees transfers summary -o out/
```

The selectable components are exactly `Genomes.WRITE_PARTS`: `profiles`, `trace`, `trees`, `events`,
`transfers`, `summary` (the CLI also accepts `all`). Each maps to these paths:

| Component | Path(s) written |
|---|---|
| *(always)* | `species_tree.nwk`, `species_nodes.tsv` |
| `profiles` | `Profiles.tsv` + `Presence.tsv` (dense), or `Profiles_sparse.tsv` (with `--sparse`) |
| `trees` | `gene_trees/<fid>_complete.nwk`, `gene_trees/<fid>_extant.nwk` |
| `events` | `gene_family_events/<fid>_events.tsv` — one file per family |
| `transfers` | `Transfers.tsv` — one row per transfer |
| `summary` | `Gene_family_summary.tsv` — per-family origin, event counts, copies |
| `trace` | `Events_trace.tsv` — the whole genealogy in one file |

`species_nodes.tsv` carries `name`, `time`, `is_leaf`, `is_extant` per node. A per-family events file
has columns `time`, `event`, `branch`, `donor`, `recipient`, `nodes`, where `nodes` is a
`;`-separated list of `role=gid` pairs. `Transfers.tsv` gives `time`, `family`, `donor_branch`,
`recipient_branch`, `parent_id`, `donor_copy_id`, `transfer_id`. `Gene_family_summary.tsv` reports,
per family, the origin time and branch, the D / T / L / speciation counts, and the extant copy and
species-present totals — the last two taken from the sparse profile so no dense matrix is built.

### Sparse profiles

`--sparse` (Python `sparse=True`) replaces the dense `Profiles.tsv` / `Presence.tsv` pair with a
single `Profiles_sparse.tsv` in the COO long format described above. The dense pair is the one output
that does not scale — it is `families × species` — so at large tip counts the sparse form is the only
viable profile output. It only affects the profile component; the CLI rejects `--sparse` unless
`profiles` is among the requested parts.

### The event trace and its fast path

`events` writes one small file per family; `trace` writes the entire log as a single
`Events_trace.tsv`. Its header is:

```
time    event   branch  donor   recipient       family  parent  child1  child2
```

One row per event, with the touched gene-lineage ids as `parent`, `child1`, `child2` (blank when the
event has no children). This is the scalable substitute for the `gene_family_events/` directory, and
it is a complete record: the gene trees can be rebuilt from it later, on demand, with no need to keep
the simulation object around.

The trace also unlocks a genuinely faster run. `simulate_genomes(..., output="trace")` returns a
`GenomeTrace` rather than a full `Genomes`: it holds the genealogy in its cheapest form and builds
the per-event Python objects and gene trees only if you ask.

```python
trace = z.simulate_genomes(tree, duplication=0.1, transfer=0.05, loss=0.1, output="trace")
trace.write("out/", include=["trace", "profiles"])   # the cheap components
genomes = trace.genomes()                             # promote on demand
```

`trace.write(...)` defaults to `include=("trace", "profiles")`, the two components that cost almost
nothing; writing `Events_trace.tsv` from the engine's raw columns never materialises a single
`EventRecord`. Requesting any heavier component (`trees`, `events`, `transfers`, `summary`)
transparently promotes the trace to a full `Genomes` first, so the fast path is exactly the subset
`{trace, profiles}`. On the command line the same fast path fires automatically when `--write` is
`trace`, or `trace profiles`.

::: note
The compact trace carries **no speciation rows** — a lineage keeps its id across speciations, so they
are redundant and omitted. When the trace is promoted (via `trace.genomes()`) or reloaded with a
species tree, it is *replayed* against that tree to re-insert the implied speciations and re-mint
per-instance gene ids, yielding the full O/S/D/T/L genealogy that ordinary reconstruction expects.
:::
