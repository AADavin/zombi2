# Conditioning

The book has so far run the levels one at a time: species tree, then genomes, then sequences. They already depend on the tree, but sometimes you want more than that background dependency. There are two ways — conditioning, in this chapter, and joining, in the next.

**Conditioning** takes a value that has already evolved and makes it drive a rate in the run that comes next. Three examples:

- **Cave fish lose their eyes.** A habitat trait, cave or surface, has already evolved down the species tree. Wherever a lineage sits in the dark, its genes are lost four times faster than on the surface.
- **Endosymbionts shed their genomes.** A lifestyle trait, free-living or host-restricted, drives gene loss across the board, so the lineages that moved inside a host are the ones that end up with the small genomes.
- **Competent bacteria pick genes up.** A trait for natural competence raises the rate at which new gene families appear in a lineage, so the trait leaves its mark on gene *gain* rather than on loss.

In each case one evolved value is read by another rate. The value read is the **driver**; the rate reading it is the **target**. They are not interchangeable, and the asymmetry is the point: a driver is a *value* that already varies from lineage to lineage — a habitat state, a gene count — while a target is a *rate*, a "how often", multiplied by a factor the driver's value picks out. The arrow runs one way, and nothing flows back.


![The shape of a conditioned run. The **driver** is a level already simulated — a habitat trait here, its two states shown below it. The **target** is a rate in the run that comes next. The **modifier** is what joins them: `DrivenBy` carries one multiplier per state of the driver, so a branch's habitat sets that branch's loss rate. The driver is finished and written to a file before the second run starts, which is what lets this be two ordinary commands.](figures/conditioning_print.png){width=95%}

That one-wayness is what makes conditioning cheap. The habitat is unaffected by how many genes a lineage has, so the trait can be grown on its own and written to a file before the genome run that reads it starts. Two ordinary commands, in order.

## The mechanism

`mod.DrivenBy` takes a `source` and a `mapping`. The `source` says which driver is read, and it is either a finished result object or the path to the event log that result wrote — the same model either way, since the file is only how a second process reads what the first one grew. The `mapping` says what the driver's value does to the rate, in one of three shapes: a **Table**, one factor per state of a discrete driver (`{"cave": 4.0, "surface": 1.0}`, and any state you leave out keeps its rate unchanged); a **Curve**, a function of a numeric driver (`lambda n: math.exp(0.05 * n)`, with `bound=` to cap it); or a **Scalar**, the log-link `exp(strength · value)` for a driver that is already an indicator or a single covariate. Whatever the shape, the factor is dimensionless and non-negative, because it multiplies a rate.

The canonical case is a trait driving gene loss: lineages in the dark lose genes faster than lineages in the light. It is two ordinary runs, in order.

```python
from zombi2 import species, traits, genomes
from zombi2.rates import modifiers as mod

tree = species.simulate_species_tree(birth=1.0, death=0.3, n_extant=20, seed=1)
# 1. grow the driver: a habitat trait down the species tree
habitat = traits.simulate_discrete(tree, states=["cave", "surface"], switch=0.1, seed=1)

# 2. grow the genomes, with loss reading the habitat on each lineage
genomes.simulate_genomes_family(tree,
    loss = 0.25 * mod.DrivenBy(habitat, {"cave": 4.0, "surface": 1.0}),
    duplication=0.2, origination=0.5, seed=2)
```

Passing `habitat` hands the result over in memory; passing `"out/trait_events.tsv"` after `habitat.write("out/", outputs=("events",))` hands over the same thing on disk. Either way the driver is read **wherever it changes**, not once per branch: the event log gives each branch's constant stretches, and the genome engine steps its Gillespie at every switch, so a lineage that changes habitat halfway down a branch loses genes at one rate before the switch and another after it. A discrete driver is therefore followed exactly. A continuous one has no switches to step at, so it is approximated: the branch is cut into stretches of at most `step` time units (`DrivenBy(..., step=)`, by default 1% of the tree's height) and read at each stretch's midpoint.

What can drive is a trait, discrete or continuous; a named gene family, present or absent; a module's completion, a number between 0 and 1; or a finished run's GC content. The gene and the module are below, under *A gene as the driver*, and GC under *A sequence as the driver*. What a trait can drive fits in six rows:

| The driver | What it drives | Written like this | Mapping |
|---|---|---|---|
| a discrete trait | `loss`, `duplication`, `origination`, `transfer` — the gene-family rates, at the **family** resolution | `loss = 0.25 * mod.DrivenBy(source, {…})` | Table |
| a discrete trait | **every rate** of an **ordered** run — those four, plus `inversion`, `transposition`, `translocation` and the chromosome tier | `inversion = 0.3 * mod.DrivenBy(source, {…})` | Table |
| a discrete trait | **every rate** of a **nucleotide** run — the same list, on a genome measured in base pairs | `inversion = 1.5 * mod.DrivenBy(source, {…})` | Table |
| a discrete trait | **every extent** of an **ordered** or **nucleotide** run — how much an event takes, rather than how often one happens | `loss_extent = 150 * mod.DrivenBy(source, {…})` | Table |
| a discrete trait | `transfer_to` — which lineage a transfer lands on, at **every** resolution | `transfer_to = mod.DrivenBy(source, {…})` | Table, or Between (reads the donor too) |
| a discrete or continuous trait | `substitution` — how fast the sequences inside a gene evolve, at the **sequences** level | `substitution = 0.05 * mod.DrivenBy(source, {…})` | Table, Curve or Scalar |

The last row is the one level further down, and Chapter 7 covers it beside the clocks it sits with: the factor multiplies the substitution rate on each species branch, so a lineage's state sets how fast its genes' sequences change. It composes with either lineage clock.

Driving a **rate** makes a lineage delete more often; driving an **extent** makes each deletion bigger. Both say "this lineage sheds more", they are different processes, and set together they multiply. An extent's unit is set by the resolution, genes at ordered and base pairs at nucleotide, and a driven extent is written in Python: the `--*-extent` flags take a plain number.

## Driving transfer

A transfer joins two lineages, a donor and a recipient. A trait can drive either end, and the two are different models.

Driving the **rate** drives the donor. `transfer = 0.1 * mod.DrivenBy(competence, {"competent": 3.0, "normal": 1.0})` makes a competent lineage donate three times as often as a normal one, so the run has more transfer in it.

Driving **`transfer_to`** drives the recipient. `transfer_to = mod.DrivenBy(competence, {"competent": 3.0, "normal": 1.0})` makes a competent lineage three times likelier to be the one a transfer lands on. No rate changes: the same transfers happen and go somewhere else.

The two expressions look alike, but their numbers mean different things. In a rate the number multiplies the rate of the lineage it is read on. In `transfer_to` it is a **weight**: the engine reads it on every lineage alive at that instant and draws the recipient in proportion, so five candidates at weight 1 and five at weight 2 send two thirds of the transfers to the weight-2 group. Weights are normalised, so doubling them all changes nothing — and that is why `transfer_to` takes the modifier on its own, with no base in front of it. `transfer_to = 1.0 * mod.DrivenBy(...)` is an error. A weight of 0 means the lineage cannot receive, and when every candidate weighs 0 the transfer has nowhere to land, so it does not happen.

Driving the rate and driving `transfer_to` are independent, and a run may use either or both:

```python
competence = traits.simulate_discrete(tree, states=["competent", "normal"],
                                      switch=0.3, seed=1)

genomes.simulate_genomes_family(tree,
    transfer    = 0.1 * mod.DrivenBy(competence, {"competent": 3.0, "normal": 1.0}),
    transfer_to =       mod.DrivenBy(competence, {"competent": 3.0, "normal": 1.0}),
    initial_families=10, seed=2)
```

### Reading the donor as well

The `transfer_to` above reads the driver on the recipient only, so it can say "competent lineages take DNA up more often" but not "genes move between two habitats and not within them". A **`Between`** kernel closes that gap: in place of one factor per recipient state, it gives a weight per ordered **(donor state, recipient state)** pair, and the engine reads the driver on the donor as well as on the candidate.

```python
from zombi2.rates.mapping import Between

habitat = traits.simulate_discrete(tree, states=["marine", "soil"], switch=0.3, seed=1)
within = Between({("marine", "marine"): 1.0, ("soil", "soil"): 1.0}, default=0.0)

genomes.simulate_genomes_family(tree,
    transfer    = 0.5,
    transfer_to = mod.DrivenBy(habitat, within),
    initial_families=10, seed=2)
```

Every transfer now stays within one habitat, because every cross-habitat pair weighs 0. `default` (1.0) is the weight of any pair you leave out, so `Between({("marine", "soil"): 3.0})` enriches one direction against a baseline instead of forbidding the rest, and `default=0.0` is the "only the flows I name" idiom. A `Between` is always a recipient weight: driving a rate with one is refused, because a rate has no donor to read it on.

This is the trait-driven twin of Chapter 4's `Clades`. There the groups come from the tree, here from an evolved trait; the kernel and the steering are the same. `transfer_to` takes one rule, so a driven `transfer_to` cannot also be `"distance"`.

## A trait driving another trait

Both ends of the arrow can sit at the same level, and nothing changes: a trait grown along a fixed tree can be finished first, so this is conditioning, written with the same `DrivenBy`. The target is `rate` for a continuous trait, its variance-rate σ², and `switch` for a discrete one.

```python
depth = traits.simulate_discrete(tree, states=["deep", "shallow"], switch=0.3, seed=3)

# a body size that diffuses four times faster in the deep
size = traits.simulate_continuous(tree, start=0.0,
    rate = 1.0 * mod.DrivenBy(depth, {"deep": 4.0, "shallow": 1.0}), seed=4)

# a diet that switches five times faster there
diet = traits.simulate_discrete(tree, states=["carnivore", "herbivore"],
    switch = 0.2 * mod.DrivenBy(depth, {"deep": 5.0, "shallow": 1.0}), seed=5)
```

A `switch` written per transition drives only the transitions you name — `{"carnivore->herbivore": 0.2 * mod.DrivenBy(depth, {…}), "herbivore->carnivore": 0.2}`.

## What can be driven, and what cannot yet

A genome run takes a driver at all three resolutions. At the **family** resolution that is the four gene-family rates in the table above. At the **ordered** resolution it is every rate the engine has, so a trait can make a lineage rearrange its gene order or reshape its karyotype more often, and a driven extent makes each rearranged block longer. At the **nucleotide** resolution it is the same list on a genome measured in base pairs, which is where a driven `loss` becomes genome reduction as it is usually meant. At the **sequences** level it is `substitution`, the one rate that level has. `transfer_to` is not a rate, but it takes a driver at all three resolutions too: the same rules, kernel and weights everywhere, because only what a transfer *moves* differs between them.

Two things are refused, and they are refused for different reasons. **`ByFamily` and `DrivenBy` cannot be set in the same run** at the family or ordered resolution, because one weights lineages by a driver and the other weights the segment by what it covers, so combining them means weighting by the product (`family_speed` counts as a `ByFamily` here). That is a limit of the code. **A genome cannot be driven by a sequence**, and that one is a limit of the model: a sequence lives inside a gene, so it was grown along the gene trees the genome run produced, and a genome reading it back would be conditioning a run on its own output. Genomes and Sequences can only be *joined*, and no joint engine for the pair exists. GC content therefore drives what comes after a sequence — a trait on the same species tree, or a further sequence run. Either way a driven rate an engine cannot honour raises, naming the reason, rather than being silently dropped.

## On the command line

Conditioning folds into the target level's own command. There is no conditioning command and no object to build: you grow the driver, then make an ordinary genome run whose `loss` happens to be `DrivenBy` instead of a bare number, with the rate in its written form.

```bash
# 1. a species tree
zombi2 species out/ --birth 1 --death 0.3 --n-extant 20 --seed 1

# 2. the driver: a habitat trait, writing its event log
zombi2 traits out/ --kind discrete \
    --states cave,surface --switch 0.1 --seed 1 --write values tree events

# 3. the target: genomes whose loss reads that trait
zombi2 genomes out/ --duplication 0.2 --origination 0.5 --seed 2 \
    --loss "0.25 * DrivenBy('out/traits/trait_events.tsv', {'cave': 4.0, 'surface': 1.0})"

# 4. one level further down: sequences whose substitution rate reads the same trait
zombi2 sequences out/ --model hky85 --length 1000 --seed 3 \
    --substitution "0.05 * DrivenBy('out/traits/trait_events.tsv', {'cave': 0.5, 'surface': 1.0})"
```

Each conditioned run writes a `conditioned_on` file naming the levels it read, so re-running the trait afterwards refuses rather than leaving the runs beneath it silently mismatched. Pass `--force` to re-run it and clear them. The guard works **within one run directory**: it is the directory that holds the record, so a driver and a target written to two different directories (with `--from`) are not linked, and re-running the driver there will not warn you.

Both halves of transfer take that same text: the rate with a base number in front of it, the recipient weight without one.

```bash
# the driver: a competence trait, into its own directory
zombi2 traits comp/ --kind discrete --from out/ \
    --states competent,normal --switch 0.3 --seed 1 --write events

driver=comp/traits/trait_events.tsv
zombi2 genomes comp_genomes/ --from out/ --initial-families 10 --seed 2 \
    --transfer    "0.1 * DrivenBy('$driver', {'competent': 3.0, 'normal': 1.0})" \
    --transfer-to "DrivenBy('$driver', {'competent': 3.0, 'normal': 1.0})"
```

A trait target reads that same file the same way, and every rate flag it has takes the expression: `--rate` for a continuous trait's variance-rate, `--switch` for a discrete one's switching rate. Use `--name` to keep two traits in one run directory, since a trait driving a trait needs somewhere for both to sit.

```bash
zombi2 traits size/ --from out/ --kind continuous --start 0.0 --seed 4 \
    --rate "1.0 * DrivenBy('out/traits/trait_events.tsv', {'cave': 4.0, 'surface': 1.0})"

zombi2 traits out/ --kind discrete --name diet --states plant,fish --seed 5 \
    --switch "0.2 * DrivenBy('out/traits/trait_events.tsv', {'cave': 5.0, 'surface': 1.0})"
```

`--switch` also takes the two shapes its keyword does — a `{'a->b': rate}` dict and a `k x k` matrix — and each entry of those is a rate in the same written form, so only the transitions you name need be driven.

## A gene as the driver

Everything so far has a trait doing the driving. A **gene family** can do it too, and it is the same
relation the other way round: a trait can make gene loss faster, and this makes a trait's rate depend
on whether a lineage carries a gene.

Name the family when you grow the genomes, then ask the run for its presence:

```python
tree = species.simulate_species_tree(birth=1.0, death=0.2, n_extant=40, seed=4)

# 1. grow the driver: genomes in which one family is named, so it can be referred to
g = genomes.simulate_genomes_family(tree, initial_families=20, family_names=["tox"],
                                    duplication=0.1, loss=0.15, seed=9)

# 2. grow the trait, with its switch rate reading whether that family is there
pathogenicity = traits.simulate_discrete(
    tree, states=["harmless", "pathogenic"], start="harmless", seed=2,
    switch = 0.1 * mod.DrivenBy(g.presence("tox"), {"present": 8.0, "absent": 1.0}))
```

`g.presence("tox")` is a driver like any other, so the mapping is an ordinary table over its two
states, `present` and `absent`. Only families you **named** with `family_names=` can be asked for: a
family that arose during the run has an id but nothing stable to call it by.

The signal is exact and it changes **during** a branch, not only at the nodes. A lineage that loses
its last copy of the family halfway along a branch is `present` before that instant and `absent`
after it, and the driven rate switches there — the same treatment a trait's mid-branch switch gets.
Presence is read off the family's gene tree, so a lineage that gains the family by transfer has it
from the moment the transfer landed.

This is conditioning, so the genome is finished before the trait starts and does not react to it.
For a gene whose presence shapes the tree that the genome is itself evolving on, see the next
chapter: `joint` grows both at once, and there the target is speciation rather than a trait.

### A module rather than a single gene

Genes rarely act alone. A **module** is a named group of families — a pathway, a complex, an operon
— and what a lineage has of it is a matter of degree: all six flagellar genes, or four of them, or
none. `modules=` declares the grouping and `completion` reads it, as a number between 0 and 1:

```python
from zombi2.rates.mapping import Curve

flg = [f"flg{i}" for i in range(6)]

g = genomes.simulate_genomes_family(tree, initial_families=20, family_names=flg,
                                    modules={"flagellum": flg},
                                    duplication=0.05, loss=0.2, seed=9)

motility = traits.simulate_discrete(
    tree, states=["sessile", "motile"], start="sessile", seed=2,
    switch = 0.05 * mod.DrivenBy(g.completion("flagellum"),
                                 Curve(lambda f: 0.05 + 30.0 * f ** 4)))
```

Completion is a **continuous** driver, so it takes a `Curve` rather than a table — and that is where
a threshold belongs, alongside every other response shape. `lambda f: 8.0 if f > 0.8 else 1.0` reads
"eight times faster once four fifths of it is there"; the `f ** 4` above is a softer version of the
same idea, where the last genes matter most.

**Why a fraction and not a yes-or-no.** Under independent loss, the chance that *every* family of a
module survives falls off geometrically with its size. Measured on a 200-tip tree: a module of three
was complete at 189 tips, one of six at none of them. A complete/incomplete driver would therefore be
a constant for anything but the smallest modules, and a run driven by a constant has told you
nothing. The fraction is always informative, and it degrades gracefully — a module of one family is
exactly that family's presence, 1 or 0.

Members must be families you named with `family_names=`. An anonymous family's id comes from the
order events happened to fire in, so a module built on one would mean something else at another seed.

### At the nucleotide resolution

All three resolutions answer, and `presence` and `completion` mean the same thing at each. What
differs is how a family is named and how it can be taken away. At the nucleotide resolution a family
is a **declared gene**, so its name comes from the GFF that declared it (its `ID` or `Name`) — the
evenly-spaced `genes=` layout lays its genes down unnamed, and there is nothing there to ask for. And
what removes it is an arc of DNA rather than a whole copy, so a gene is gone once a deletion has taken
it: the driver still reads `present` and `absent`, off the same recovered gene tree.

```python
# a GFF is a path, or the lines themselves — here two genes on one 3 kb replicon
annotation = ["##sequence-region c1 1 3000",
              "c1\tzombi2\tgene\t201\t800\t.\t+\t.\tID=dnaA",
              "c1\tzombi2\tgene\t1401\t2000\t.\t+\t.\tID=dnaN"]

nt = genomes.simulate_genomes_nucleotide(tree, gff=annotation, loss=1.0, loss_extent=300,
                                         modules={"operon": ["dnaA", "dnaN"]}, seed=9)

traits.simulate_discrete(tree, states=["harmless", "pathogenic"], start="harmless", seed=2,
                         switch = 0.1 * mod.DrivenBy(nt.presence("dnaA"),
                                                     {"present": 8.0, "absent": 1.0}))
```

## A sequence as the driver

The level furthest down can drive too, through the one number a finished sequence says that another
level can use: its **GC content**. `result.gc()` gives the fraction of a lineage's DNA that is G or C
at every instant, and it is read like any other continuous driver.

```python
from zombi2.sequences import simulate_sequences
from zombi2.sequences.substitution_models import hky85, lg

seqs = simulate_sequences(g, model=hky85(2.0), length=300, seed=1)

traits.simulate_discrete(tree, states=["mesophile", "thermophile"], start="mesophile", seed=2,
                         switch = 0.2 * mod.DrivenBy(seqs.gc(),
                                                     Curve(lambda x: 20.0 ** (x - 0.5))))
```

**GC is pooled over the whole lineage**, across every family the run evolved, because that is the
quantity the field measures and the one a lineage's mutational bias acts on. One family's GC is not
offered: it is undefined wherever that family is absent, and a driver has to answer for every branch
the target walks. Grow that family in a run of its own if its GC alone is what should drive the rate.

A run gives one sequence per gene-tree node, so pooling gives one GC per species node — the value at
the end of its branch — and the path between two nodes is the straight line between them, cut into
stretches of at most `step`, exactly as a continuous trait's is. It is a number, so it takes a `Curve`
or a `Scalar`, never a `{state: factor}` table.

### Any letters, so any amino-acid frequency

`gc()` is the named door onto a general one. `composition(letters)` counts whichever letters of the
run's own alphabet you ask for, so an amino-acid frequency is the same driver with different letters —
one residue, or a set of them:

```python
proteins = simulate_sequences(g, model=lg(), length=300, seed=1)

traits.simulate_discrete(tree, states=["mesophile", "thermophile"], start="mesophile", seed=2,
                         switch = 0.2 * mod.DrivenBy(proteins.composition("KR"),
                                                     Curve(lambda x: 40.0 ** (x - 0.1))))
```

`gc()` *is* `composition("GC")`, with one extra check: it refuses a protein run, where G and C are
glycine and cysteine and the call is ambiguous rather than wrong. Letters outside the run's alphabet
are refused too — they occur nowhere, so the driver would read 0.0 on every lineage and the run would
be the undriven model wearing a driven rate.

What this does **not** reach is a pattern — a motif, a site. A composition is a count over a lineage's
whole complement; a motif is a property of one copy, and a sequence exists only at the nodes, so a
motif driver would need both a pooling rule and an interpolated switch point rather than the exact
one a gene's presence has.

**It drives forwards only.** A sequence lives inside a gene, so it was grown along the gene trees a
genome run produced. Handing GC back to a genome rate asks a run to condition on its own output, and
that pair — Genomes and Sequences — can only be joined, never conditioned; the genome level refuses it
by name. What GC can drive is what comes after it: a trait on the same species tree, or a further
sequence run.

## Outputs

Conditioning adds no format. A conditioned run writes the target level's own files, plus a
`conditioned_on` record naming the levels its rates read; the driver's own files sit in the same
run directory, so the pairing that produced the pattern is kept on disk.

| File | What it holds |
|---|---|
| the target level's usual outputs | whatever that level writes — see its own chapter |
| `conditioned_on` | the levels this run's rates read, one per line |
| `trait_values.tsv` · `trait_events.tsv` · `trait_tree.nwk` | the driver, when it was a trait |

Appendix B gives the columns and the formats.
