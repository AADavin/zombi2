# Genomes I: gene families

Genomes live inside the species tree, and they can be simulated at three **resolutions**, one per chapter: **family** here, **ordered** in Chapter 5, **nucleotide** in Chapter 6. The simplest is a gene family evolving along the tree; the most detailed tracks every nucleotide across several chromosomes.

The **family** resolution is genomes made of gene families and nothing more: no position along a chromosome, no DNA sequence. Genes are copied, lost, born from nothing, and passed sideways between lineages.

## The four events

A genome at the family resolution evolves by four kinds of event, applied to every lineage as it runs down the species tree:

- **Duplication** — a gene copy is copied, so its family gains a member in that lineage.
- **Transfer** — a copy jumps from one lineage to another that is alive at the same moment. 
- **Loss** — a gene copy is deleted; a family that loses its last copy is gone from that lineage.
- **Origination** — a brand-new family appears in a lineage, with one copy.

You give a rate for each, and the events play out along the tree from the initial genome, speciation handing a lineage's genome down to both children. Out comes the genome of *every* lineage together with the event log that produced it.

```python
from zombi2 import species, genomes
from zombi2.rates import modifiers as mod

tree = species.simulate_species_tree(birth=1.0, death=0.3, n_extant=20, seed=1)
g = genomes.simulate_genomes_family(
    tree, duplication=0.2, loss=0.25, origination=0.5, initial_families=20, seed=1)
```

The initial genome, at the begining of the stem, starts with `initial_families` families of one copy each; from there the four rates drive everything.

### Families can evolve at different paces

By default, ZOMBI2 simulates gene families with the same rates. In reality, some gene families are more prone to be transferred than others (think of antibiotic resistance); some families are lost more easily (accesory genes); some families are rarely lost (core gene families like ribosomal proteins). There are multiple ways to simulate families that evolve at different paces. Two ways are to vary **one rate** family by family, and to give each family **one tempo** that scales all of its rates.

The first says that families differ in a particular respect. A resistance gene is transferred more than most, and that says nothing about how often it is duplicated, so `ByFamily` goes on `transfer` alone. Put it on `loss` instead and you separate the accessory families from the core ones, leaving gain untouched.

The second says that families differ overall — some are simply more volatile than others, in every way at once. That is `family_speed`, one draw per family multiplying every rate that family has.

Either way each family draws one multiplier and keeps it for the whole run, and the draw is mean-corrected so `E[factor] = 1`: widening the spread spreads families further apart without moving the average family off the base rate you typed. So where you put it decides what varies *together*:

```python
# each rate varies by family on its own — losing fast is not duplicating fast
g = genomes.simulate_genomes_family(
    tree,
    duplication = 0.2  * mod.ByFamily(spread=0.6),
    transfer    = 0.1  * mod.ByFamily(spread=0.6),
    loss        = 0.25 * mod.ByFamily(spread=0.6),
    initial_families = 100, seed = 42)

# one tempo per family, scaling every rate it has — a fast family is fast at everything
g = genomes.simulate_genomes_family(
    tree, duplication=0.2, transfer=0.1, loss=0.25,
    family_speed = mod.ByFamily(spread=0.5),
    initial_families = 100, seed = 42)
```

The two compose: `family_speed` for a family's overall tempo, and a `ByFamily` on one rate for extra variation particular to it. On the command line the rate keeps its written form, `--loss "0.25 * ByFamily(spread=0.6)"`.

### How large a family may get

Growth compounds: a duplication rate above the loss rate multiplies without bound, and with `ByFamily` some families draw a rate well above the one you typed. So a family's copies **within one genome** are capped, and the cap is on by default.

The cap is a plain count of copies **in one genome** — there is no "per what?" left to ask, because it is compared against a single genome's own copies:

```python
from zombi2 import species, genomes

tree = species.simulate_species_tree(birth=1.0, n_extant=8, seed=1)

# ten copies of one family in one genome — the default
g = genomes.simulate_genomes_family(tree, duplication=0.5, loss=0.1,
                                    initial_families=10, max_family_size=10, seed=1)

# no ceiling: what you want when you are measuring rates, since a cap that bites
# discards events and pulls the realised rates below the ones you declared
g = genomes.simulate_genomes_family(tree, duplication=0.5, loss=0.1,
                                    initial_families=10, max_family_size=None, seed=1)
```

At the cap the family stops duplicating, and the ceiling holds for arrivals too, so a transfer cannot push it past sideways.

It used to be written with a scope, and `scope.PerLineage(n)` multiplied that number by the size of the *species tree* — so the shipped default was over a thousand copies on a fifty-tip tree, and the cap you got moved when you added species. Both scope spellings are refused now, with the arithmetic each used to imply, so an old script fails loudly rather than running at a cap different from the one it reads.

### When a genome empties

There is a ceiling but no floor. Loss is counted per copy, and the last copy is a copy like any other, so a loss rate well above the duplication and origination rates strips a lineage of every gene it has:

```python
g = genomes.simulate_genomes_family(tree, loss=5.0, initial_families=10, seed=3)
g.summary()["empty_genomes"]            # 8 — every extant genome came out with none
```

That is a real outcome of the model, not a failure. The lineage carries an empty genome to the tip: `profiles.tsv` has no row for it, and there is no gene tree for a sequence to run down.

It is easy to miss, because an empty genome shows up as an absence. So the run says so. `genome_summary.json` reports `empty_genomes`, the number of extant genomes that came out with no genes at all, and `zombi2 genomes` prints a line on standard error when it happens. Lower `loss`, raise `origination`, or start with more `initial_families` if that is not the model you meant.

The chromosome-based resolutions do have a floor: a loss never takes a chromosome below its last gene. That is a statement about what a chromosome is, not a bound on genome size — see Chapter 5.

## What the rate depends on

The rates follow the **same grammar as the species level** (`base` optionally wrapped in a scope, optionally multiplied by modifiers). The scope answers *per what*, and the default is the natural one for each event. Duplication, transfer, and loss are counted **per copy**: a family with ten copies is ten times as likely to duplicate or lose one as a family with a single copy, which is what you want: more genes, more chances. Origination is counted **per lineage** (i.e. per branch of the species tree): acquiring a wholly new family is a property of the lineage, not of any gene it already has.

Rates can also depend on **time**. Multiplying a base rate by an `OnTime` modifier makes it change at set moments — the skyline, or episodic, genome: fast early and slow later, or any schedule you give.

```python
# lots of new families early, then origination shuts off after time 2
g = genomes.simulate_genomes_family(
    tree, origination=1.0 * mod.OnTime({0: 1.0, 2: 0.0}), seed=1)
```

## Lateral gene transfers

When a transfer fires, a copy is picked from the whole pool of live genes, and it is delivered to another lineage that is **alive at that same instant**.

Three arguments shape what a transfer does:

- **`transfer_to`** — who receives. `"uniform"` (the default) picks any other contemporaneous lineage with equal chance; `"distance"` makes closer relatives likelier, weighting recipients by how far they sit from the donor on the tree. The distance version is *scale-free*. A third rule, `Clades(...)`, weights recipients by **named clades of the tree** — see below.
- **`replacement`** — what happens on arrival. By default the incoming copy is **additive**: the recipient simply gains a copy. With `replacement=True` it **overwrites** a copy of the same family already present, and falls back to additive when the recipient has none.
- **`self_transfer`** — whether a lineage may donate to itself. Off by default. With additive arrival the lineage gains a copy, so the gene content changes as it would under a duplication, but the event is recorded as a transfer. 

```python
tree = species.simulate_species_tree(birth=1.0, death=0.4, n_extant=30, seed=7)
# horizontal transfer biased toward close relatives, overwriting resident copies
g = genomes.simulate_genomes_family(
    tree, transfer=0.5, transfer_to="distance", replacement=True,
    origination=0.4, initial_families=10, seed=3)
```

One consequence is worth stating plainly: a transfer can arrive **from a lineage that later goes extinct** [@szollosi2013lgtdead]. A genome run happens on the complete tree, dead branches included, so a gene can enter a survivor from a donor that leaves no other trace. This was in fact the feature that gave originally the name to this software.

### Transfer between named clades

`"distance"` biases transfer by relatedness, but sometimes you want to name the groups yourself — "let genes flow between these two clades, and nowhere else." `Clades` does that. You name each clade — by a few of its tips (the clade is the subtree below their MRCA) or by a node id — and give a `Between` table of weights, one for each ordered **(donor clade, recipient clade)** pair.

```python
from zombi2 import species, genomes

sp = species.simulate_species_tree(birth=1.0, death=0.3, n_extant=16, seed=1)
# genes flow only BETWEEN clade A and clade B — never within either, never to the rest
flows = genomes.Between({("A", "B"): 1.0, ("B", "A"): 1.0}, default=0.0)
g = genomes.simulate_genomes_family(
    sp, transfer=1.0, initial_families=20, seed=2,
    transfer_to=genomes.Clades({"A": ["n27", "n28"], "B": ["n21", "n26"]}, flows))
```

Each entry is a weight, read the same way `"distance"`'s weights are: normalised over the lineages alive at the instant a transfer fires. Naming only `("A", "B")` and `("B", "A")` and setting `default=0.0` means every other pairing weighs 0 — a clade-A donor can reach clade B but not another clade-A lineage, and the rest of the tree neither sends nor receives. Drop the `default=0.0` and unlisted pairs return to weight 1 (baseline), so `Between({("A", "B"): 5.0})` *enriches* A→B fivefold while leaving everything else to happen normally. A weight of 0 means "cannot receive", exactly as in Chapter 9: when a donor's every candidate weighs 0, the transfer has nowhere to land and does not fire.

`Clades` is written in Python. On the command line `--transfer-to` takes `uniform`, `distance`, or a `DrivenBy` recipient weight (Chapter 9).

## The `FamilyGenomesResult` object

`simulate_genomes_family` returns a **FamilyGenomesResult** which carries:

- `.complete_tree` — the species tree the genomes ran on, extinct lineages and all.
- `.genomes` — a dict from node to that node's genome.
- `.initial_genome` — the genome the run **started** with, at the root lineage's origination. It is not `.genomes[root]`: a node sits at the **end** of its branch, and the root branch is real simulated time, so events happen along it. Written to its own `initial_genome.tsv`, with no `lineage` column, because it belongs to no node.
- `.events` — the event log: every gene event with its time and lineage — origination, duplication, transfer, loss, and the *speciations* at a split.
- `.profiles` — the family × extant-species copy-count table.
- `.gene_trees` — one `GeneTree` per family.
- `.seed` — the seed, so the run reproduces.

and two methods:

- `.family_counts(node_id)` — a `Counter` collapsing a node's genome to `family → number of copies`, when you want the multiset rather than the individual copies.
- `.write(dir, outputs=[...])` — materialise the outputs to disk, listed under *Outputs* below.

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

**Gene trees** are the deeper output. Every family has two trees: the `.complete` tree with every gene lineage, and the `.extant` tree pruned to the genes that survive.

```python
gt = g.gene_trees[7]                 # the gene tree of family 7
gt.to_newick("extant")               # the surviving copies as Newick ...
gt.to_newick("complete")             # ... or the whole genealogy
gt.origination                       # when the family began
```

The root of a gene tree carries a branch length, as the species tree's does. A family starts at its origination, and the founding gene lives a while before its first duplication, transfer or speciation; that wait is the root's branch, the family's stem. A gene that originated and never split at all is a one-node tree, written as its own lifespan — `g55:0.263097;`.

## Usage from Python

The whole range is one function call:

```python
from zombi2 import species, genomes
from zombi2.rates import modifiers as mod

tree = species.simulate_species_tree(birth=1.0, death=0.3, n_extant=30, seed=1)

# a plain duplication–loss–origination run
g = genomes.simulate_genomes_family(
    tree, duplication=0.2, loss=0.25, origination=0.5, initial_families=20, seed=1)

# origination only — every family stays a single copy, none is ever duplicated
g = genomes.simulate_genomes_family(tree, origination=0.6, seed=1)

# a skyline: new families pour in early, then origination shuts off after time 2
g = genomes.simulate_genomes_family(
    tree, origination=1.0 * mod.OnTime({0: 1.0, 2: 0.0}), seed=1)

# horizontal transfer, biased toward close relatives, overwriting resident copies
g = genomes.simulate_genomes_family(
    tree, transfer=0.5, transfer_to="distance", replacement=True,
    origination=0.4, initial_families=10, seed=3)

# the genomes you observe are the extant tips
observed = {n.id: g.genomes[n.id] for n in g.complete_tree.extant_leaves()}

# and the outputs, derived from that history
g.profiles.matrix                                # family × extant-species copy counts
some_family = next(iter(g.gene_trees))
g.gene_trees[some_family].to_newick("extant")    # that family's surviving gene tree
g.write("out/")                                  # the run's files, on disk
```

## Usage from the CLI

`zombi2 genomes` evolves gene families along a species tree read from a Newick file. The family resolution is the default, and each rate is a plain number:

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
| `genome_events.tsv` | the gene genealogy — every event with its time and lineage |
| `profiles.tsv` | family × extant-species copy counts |
| `genomes.tsv` | every node's gene content, ancestors included, one row per copy |
| `initial_genome.tsv` | the genome the run started with |
| `gene_trees/` | one Newick per family, complete and extant |

A gene copy is written `g<id>`, the same token the gene-tree leaf and the alignment header carry, so a
row of `genomes.tsv` joins straight onto a tree or a sequence, and onto the event that made it.
Appendix B gives the columns and the formats.

## Evolving families in parallel

Because families are independent — a transfer moves a copy between lineages, but no event ever mixes two families — a run can evolve them **concurrently**, one family per worker process. It is off by default; `parallel` turns it on.

```python
g = genomes.simulate_genomes_family(
    tree, duplication=0.2, loss=0.25, origination=0.5,
    initial_families=1000, seed=1, parallel=8)      # 8 workers
```

`parallel=True` uses every core and an integer sets the worker count; on the command line it is `--parallel` for all cores or `--parallel 8` for eight. It is a **separate engine**, not a faster path through the default one: each family draws from its own random stream, so the result is identical for any worker count, but it differs from a serial run of the same seed — both are valid draws of the same process.

A **conditioned** rate (Chapter 9) runs here too. Conditioning does not couple families: the driver was grown before this run and is an input to it, so a lineage's factor at a given moment is the same number whichever family is asking, and no family can reach another through it. Each worker reads the driver as a lookup, and the decomposition is untouched.

For very large runs — hundreds of thousands of families, or a million — the difficulty stops being speed and becomes memory: the finished result itself no longer fits. `stream_to` writes each family to a directory the moment it is done, and hands back a light handle holding a path rather than a `FamilyGenomesResult` holding everything. Memory then stays flat however many families you run — a run that would have held 2 GB streams in about 40 MB — and the sequence level reads the families back off the disk afterwards. Choose which files to write with `outputs=`, exactly as `.write` takes them. On the command line this is `--stream`.

```python
run = genomes.simulate_genomes_family(
    tree, origination=2.0, initial_families=5000, seed=1,
    parallel=8, stream_to="out/", outputs=("events", "profiles"))
run.path("events")            # out/genome_events.tsv — the log, ready to replay
```

