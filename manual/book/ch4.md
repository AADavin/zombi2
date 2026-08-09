# Genomes I: gene families

Genomes live inside the species tree, and they can be simulated at three **resolutions**, one per chapter: **family** here, **ordered** in Chapter 5, **nucleotide** in Chapter 6. The simplest is a gene family evolving along the tree; the most detailed tracks every nucleotide across several chromosomes.

The **family** resolution is genomes made of gene families and nothing more: no position along a chromosome, no DNA sequence. Genes are copied, lost, born from nothing, and passed sideways between lineages.

## The four events

A genome at the family resolution evolves by four kinds of event, applied to every lineage as it runs down the species tree:

- **Duplication.** A gene copy is copied, so its family gains a member in that lineage.
- **Transfer.** A copy jumps from one lineage to another that is alive at the same moment.
- **Loss.** A gene copy is deleted; a family that loses its last copy is gone from that lineage.
- **Origination.** A brand-new family appears in a lineage, with one copy.

![The four events, one to a panel. The pale band is a **lineage** of the species tree and the dark line inside it a **gene copy**; time runs left to right. A duplication splits one copy into two inside the one lineage. A transfer sends a copy into another lineage alive at that same moment, where it appears as a new copy. A loss ends a copy. An origination starts one, and a family with it, from nothing. Speciation is not here: it is the species level's event, and what it does to a genome is hand the whole of it to both daughters.](figures/four_events_print.png){width=100%}

You give a rate for each, and the events play out along the tree from the initial genome, speciation handing a lineage's genome down to both children. Out comes the genome of *every* lineage together with the event log that produced it.

```python
from zombi2 import species, genomes
from zombi2.params import PerCopy, PerLineage

tree = species.simulate_species_tree(birth=1.0, death=0.3, n_extant=20, seed=1)
g = genomes.simulate_genomes_family(
    tree, duplication=0.2, loss=0.25, origination=0.5, initial_families=20, seed=1)
```

The initial genome, at the beginning of the stem, starts with `initial_families` families of one copy each; from there the four rates drive everything.

### Four things to know

- **Families need not evolve at the same pace.** `varying_among('families', ...)` draws a factor per family (`--loss "PerCopy(0.25).varying_among('families', LogNormal(0.0, 0.6))"`). Whether you write one draw or two decides what varies together: two separately written draws vary their rates independently, while **one** `Random` object read by several rates gives that family a single tempo across all of them — `speed = Random('families', LogNormal(0.0, 0.5))`, then `duplication=PerCopy(0.2).varying_among(speed), loss=PerCopy(0.25).varying_among(speed)`.
- **A family's copies within one genome are capped**, at `max_family_size=10` by default, because a duplication rate above the loss rate grows without bound. Set `max_family_size=None` when you are measuring rates: a cap that binds discards events and pulls the realised rates below the ones you declared.
- **There is no floor.** Loss is counted per copy and the last copy is a copy like any other, so a high loss rate can leave a lineage with nothing. `genome_summary.json` reports `empty_genomes` and the command warns, because an empty genome is otherwise invisible.
- **The chromosome-based resolutions do have a floor**: a loss never takes a chromosome below its last gene. That is what a chromosome is, not a bound on genome size.

## What the rate depends on

The rates follow the **same grammar as the species level** — a `base` in its scope, with verbs chained onto it. By default, duplication, transfer, and loss are counted **per copy**: a family with ten copies is ten times as likely to duplicate or lose one as a family with a single copy, which is usually what you want — more genes, more chances. Origination is counted **per lineage** (i.e. per branch of the species tree): acquiring a wholly new family is a property of the lineage, not of any gene it already has, and it takes no other scope.

Duplication, transfer and loss do take another, and it is worth knowing because it is a different model rather than a different spelling:

```python
from zombi2.params import PerCopy, PerLineage

loss = PerCopy(0.25)      # every copy independently at risk — a big genome loses often
loss = PerLineage(0.25)   # a deletion budget — the lineage loses at 0.25 whatever it holds
```

The second is how you write a lineage whose losses are set by its own biology rather than by how many genes it happens to carry: deletion-biased genomes lose at their own pace, and shrinking does not slow them down. The same `0.25` means something a hundred times different in a genome of a hundred genes, so this is a choice to make deliberately rather than a default to leave alone.

Rates can also depend on **time**. Chaining `changing_at` onto a rate makes it change at set moments. That is the skyline, or episodic, genome: fast early and slow later, or any schedule you give.

```python
# lots of new families early, then origination shuts off after time 2
g = genomes.simulate_genomes_family(
    tree, origination=PerLineage(1.0).changing_at({0: 1.0, 2: 0.0}), seed=1)
```

Rates can also depend on **where in the tree** a lineage sits. `Clade` names a group by a few of its tips — the clade is the subtree below their most recent common ancestor — and the rate reads which group a lineage is in, exactly as it would read an evolved trait (Chapter 9). A lineage in no named clade falls in the group `"rest"`. It is the rate-side counterpart of the `Clades` transfer rule below, reading the same fact off the same tree.

```python
from zombi2.params import Clade

# the clade below the MRCA of n27 and n51 loses genes three times as fast
g = genomes.simulate_genomes_family(
    tree, loss=PerCopy(0.25).scaled_by(Clade({"fast": ["n27", "n51"]}), {"fast": 3.0}),
    duplication=0.2, origination=0.5, seed=1)
```

## Lateral gene transfers

When a transfer fires, a copy is picked from the whole pool of live genes, and it is delivered to another lineage that is **alive at that same instant**.

Three arguments shape what a transfer does:

- **`transfer_to`**, who receives. `"uniform"` (the default) picks any other contemporaneous lineage with equal chance; `"distance"` makes closer relatives likelier, weighting recipients by how far they sit from the donor on the tree. The distance version is *scale-free*. A third rule, `Clades(...)`, weights recipients by **named clades of the tree**, described below.
- **`replacement`**, what happens on arrival. By default the incoming copy is **additive**: the recipient simply gains a copy. With `replacement=True` it **overwrites** a copy of the same family already present, and falls back to additive when the recipient has none.
- **`self_transfer`**, whether a lineage may donate to itself. Off by default. With additive arrival the lineage gains a copy, so the gene content changes as it would under a duplication, but the event is recorded as a transfer.

```python
tree = species.simulate_species_tree(birth=1.0, death=0.4, n_extant=30, seed=7)
# horizontal transfer biased toward close relatives, overwriting resident copies
g = genomes.simulate_genomes_family(
    tree, transfer=0.5, transfer_to="distance", replacement=True,
    origination=0.4, initial_families=10, seed=3)
```

Transfers can arrive **from a lineage that later goes extinct** [@szollosi2013lgtdead]. A genome run happens on the complete tree, dead branches included, so a gene can enter a survivor from a donor that leaves no other trace. That is the feature the software was originally named for.

### Transfer between named clades

`"distance"` biases transfer by relatedness, but sometimes you want to name the groups yourself: "let genes flow between these two clades, and nowhere else." `Clades` does that. You name each clade, by a few of its tips (the clade is the subtree below their MRCA) or by a node id, and give a `Between` kernel of weights, one for each ordered **(donor clade, recipient clade)** pair.

```python
from zombi2 import species, genomes

sp = species.simulate_species_tree(birth=1.0, death=0.3, n_extant=16, seed=1)
# genes flow only BETWEEN clade A and clade B, never within either, never to the rest
flows = genomes.Between({("A", "B"): 1.0, ("B", "A"): 1.0}, default=0.0)
g = genomes.simulate_genomes_family(
    sp, transfer=1.0, initial_families=20, seed=2,
    transfer_to=genomes.Clades({"A": ["n27", "n28"], "B": ["n21", "n26"]}, flows))
```

There is one group you do not have to name: **`"rest"`** is every lineage outside the clades you did name, so `Between({("rest", "A"): 8.0})` makes clade A a transfer hotspot that the whole rest of the tree donates into, without having to enumerate the rest of the tree.

Each entry is a weight, read the same way `"distance"`'s weights are: normalised over the lineages alive at the instant a transfer fires. Naming only `("A", "B")` and `("B", "A")` and setting `default=0.0` means every other pairing weighs 0: a clade-A donor can reach clade B but not another clade-A lineage, and the rest of the tree neither sends nor receives. Drop the `default=0.0` and unlisted pairs return to weight 1 (baseline), so `Between({("A", "B"): 5.0})` *enriches* A→B fivefold while leaving everything else to happen normally. A weight of 0 means "cannot receive", exactly as in Chapter 9: when a donor's every candidate weighs 0, the transfer has nowhere to land and does not fire.

`--transfer-to` takes the same four rules on the command line, in the same written form: `uniform`, `distance`, `Distance(decay=3.0)`, a `Clades(...)` expression, or a `Recipients().weighted_by(...)` recipient weight (Chapter 9). The bare `distance` is `Distance(decay=1.0)`, so write the constructor when you want another decay.

`transfer_to` chooses who receives, and takes one rule. All four rules (`"uniform"`, `"distance"` / `Distance(decay=)`, `Clades(...)` and `Recipients().weighted_by(driver, {...})`) work unchanged at the ordered and nucleotide resolutions. What differs between the resolutions is *what moves*: one gene copy here, a block of consecutive genes in Chapter 5, an arc of DNA in Chapter 6. Who receives it is chosen the same way in all three, so the rules are described once, here.

## The `FamilyGenomesResult` object

`simulate_genomes_family` returns a **FamilyGenomesResult** which carries:

- `.complete_tree`, the species tree the genomes ran on, extinct lineages and all.
- `.genomes`, a dict from node to that node's genome.
- `.initial_genome`, the genome the run **started** with, at the root lineage's origination. It is not `.genomes[root]`: a node sits at the **end** of its branch, and the root branch is real simulated time, so events happen along it. Written to its own `initial_genome.tsv`, with no `lineage` column, because it belongs to no node.
- `.events`, the event log: every gene event with its time, its family, and the gene copies it ended and began, covering origination, duplication, transfer, loss, and the *speciations* at a split. A transfer is written under the kind that says what it did on arrival: `transfer_additive` when the arriving copy is an addition, `transfer_replacing` when it overwrote a homolog. So filter for both, or on the prefix: there is no bare `transfer` row.
- `.profiles`, the family × extant-species copy-count table.
- `.gene_trees`, one `GeneTree` per family.
- `.seed`, so the run reproduces.

and its methods:

- `.family_counts(node_id)`, a `Counter` collapsing a node's genome to `family → number of copies`, when you want the multiset rather than the individual copies.
- `.has_family(node_id, name)`, whether a family you named with `family_names=` has at least one copy in that node's genome.
- `.presence(name)` and `.completion(name)`, a named family's presence and a module's completion along every lineage at every instant rather than at one node — the driver views a conditioned rate reads (Chapter 9). They need `family_names=` or `modules=` to have been declared.
- `.summary()`, what the run produced, as a dict — the payload of `genome_summary.json`.
- `.write(dir, outputs=[...])`, materialise the outputs to disk, listed under *Outputs* below.

```python
n5 = g.genomes[5]                    # the gene copies in node n5
counts = g.family_counts(5)          # {family: copies} for the same node
g.write("out/")                      # the run's files, on disk
```

## Profiles and gene trees

Two products are derived from the run's recorded history.

**Profiles** are the classic comparative-genomics view [@pellegrini1999profiles]: how many copies of each gene family sit in each extant species. They are read off the observed genomes on access, so the run itself stays lean.

```python
g.profiles.matrix        # families × extant-species copy counts, a NumPy array
g.profiles.presence      # the same as 0/1 presence/absence
g.profiles.to_tsv()      # the table as text
```

**Gene trees** are the deeper output. Every family has a `.complete` tree with every gene lineage, and, when at least one copy survives, an `.extant` tree pruned to the genes that do. A family that died out has no extant tree: `.extant` is `None`, `to_newick("extant")` returns `None`, and the run writes no `_extant` file for it.

```python
gt = g.gene_trees[7]                 # the gene tree of family 7
gt.to_newick("extant")               # the surviving copies as Newick ...
gt.to_newick("complete")             # ... or the whole genealogy
gt.origination                       # when the family began
```

The root of a gene tree carries a branch length, as the species tree's does. A family starts at its origination, and the founding gene lives a while before its first duplication, transfer or speciation; that wait is the root's branch, the family's stem. A gene that originated and never split at all is a one-node tree, written as its own lifespan, `n26_g545:0.6614353;`.

## Evolving families in parallel

Because families are independent, and no event ever mixes two families, a run can evolve them **concurrently**, one family per worker process. It is off by default; `parallel` turns it on.

```python
g = genomes.simulate_genomes_family(
    tree, duplication=0.2, loss=0.25, origination=0.5,
    initial_families=1000, seed=1, parallel=8)      # 8 workers
```

`parallel=True` uses every core and an integer sets the worker count; on the command line it is `--parallel` for all cores or `--parallel 8` for eight. It is a **separate engine**, not a faster path through the default one: each family draws from its own random stream, so the result is identical for any worker count, but it differs from a serial run of the same seed. Both are valid draws of the same process.

From a script, put the call under a `__main__` guard, the standard requirement for anything that starts worker processes, since they re-import your script and would otherwise run it again from the top:

```python
if __name__ == "__main__":
    g = genomes.simulate_genomes_family(tree, duplication=0.2, seed=1, parallel=8)
```

Leave it out and the run stops with a message saying so. A notebook or `python -c` has no script to re-import, so neither needs the guard; the library notices and runs single-process there instead.

A **conditioned** rate (Chapter 9) runs here too. Conditioning does not tie families to one another: the driver was grown before this run and is an input to it, so a lineage's factor at a given moment is the same number whichever family is asking, and no family can reach another through it. Each worker reads the driver as a lookup, and the decomposition is untouched.

For very large runs of hundreds of thousands of families, or a million, the difficulty stops being speed and becomes memory: the finished result itself no longer fits. `stream_to` writes each family to a directory the moment it is done, and hands back a light handle holding a path rather than a `FamilyGenomesResult` holding everything. Memory then stays flat however many families you run, so a run that would have held 2 GB streams in about 40 MB, and the sequence level reads the families back off the disk afterwards. Choose which files to write with `outputs=`, exactly as `.write` takes them. On the command line this is `--stream`.

```python
run = genomes.simulate_genomes_family(
    tree, origination=2.0, initial_families=5000, seed=1,
    parallel=8, stream_to="out/", outputs=("events", "profiles", "species_tree"))
run.path("events")            # out/genome_events.tsv, the log, ready to replay
```

The handle goes straight into the next level, and so does the directory: `sequences.simulate_sequences(run, …)` and `sequences.simulate_sequences("out/", …)` both reopen it. `genomes.read_run("out/")` gives the run itself back, the event log and the gene trees derived from it, for anything you would rather do in Python than on the command line.

## On the command line

`zombi2 genomes` evolves gene families along the species tree already in the run directory, or one given with `--from` (a Newick file, or another run). The family resolution is the default, and a rate is a plain number or a quoted expression in the written form of Chapter 2.

```bash
# duplication–loss–origination along a species tree
zombi2 genomes out/ --duplication 0.2 --loss 0.25 --origination 0.5 --seed 1

# horizontal transfer biased toward close relatives, overwriting resident copies
zombi2 genomes out/ --transfer 0.5 --transfer-to distance --replacement \
    --origination 0.4 --seed 3
```

## Outputs

| File | What it holds |
|---|---|
| `genome_events.tsv` | the gene genealogy: every event with its time, its family, and the copies it ended and began |
| `profiles.tsv` | family × extant-species copy counts |
| `genomes.tsv` | every node's gene content, ancestors included, one row per copy |
| `initial_genome.tsv` | the genome the run started with |
| `gene_trees/` | each family's genealogy, complete and extant |
| `species_complete.nwk` | the species tree the run evolved along, so the run stands alone; every other file is indexed by its node labels. Written by `result.write()` and by `--stream`, not by an ordinary `zombi2 genomes`, which leaves the copy under `species/` to serve |
| `genome_summary.json` | what came out: events and families counted, `empty_genomes`, and whether the family-size cap bound |

A gene copy is written `n<species>_g<copy>` — `n19_g358` — the token the gene-tree leaf, the
alignment header and the event log all carry. `genomes.tsv` gives the two halves in their own
columns, `n19` and `g358`, so join a row to a tree, a sequence or an event by putting them back
together with the `_`. The copy id alone is unique across the run, so splitting the composite the
other way works just as well. Appendix B gives the columns and the formats.
