# Conditioning

The book has so far run the levels one at a time: species tree, then genomes, then sequences. They already depend on the tree, but sometimes you want more than that background dependency. There are two ways — conditioning, in this chapter, and joining, in the next.

**Conditioning** takes a value that has already evolved and makes it drive a rate in the run that comes next. Three examples:

- **Cave fish lose their eyes.** A habitat trait, cave or surface, has already evolved down the species tree. Wherever a lineage sits in the dark, its genes are lost four times faster than on the surface.
- **Endosymbionts shed their genomes.** A lifestyle trait, free-living or host-restricted, drives gene loss across the board, so the lineages that moved inside a host are the ones that end up with the small genomes.
- **Competent bacteria pick genes up.** A trait for natural competence raises the rate at which new gene families appear in a lineage, so the trait leaves its mark on gene *gain* rather than on loss.

In each case one evolved value is read by another rate. The value read is the **driver**; the rate reading it is the **target**. They are not interchangeable, and the asymmetry is the point: a driver is a *value* that already varies from lineage to lineage — a habitat state, a gene count — while a target is a *rate*, a "how often", multiplied by a factor the driver's value picks out. The arrow runs one way, and nothing flows back.


![The shape of a conditioned run. The **driver** is a level already simulated — a habitat trait here, its two states shown below it. The **target** is a rate in the run that comes next. The **modifier** is what joins them: `DrivenBy` carries one multiplier per state of the driver, so a branch's habitat sets that branch's loss rate. The driver is finished and written to a file before the second run starts, which is what lets this be two ordinary commands.](figures/conditioning_print.png){width=95%}

That one-wayness is what makes conditioning cheap. The habitat is unaffected by how many genes a lineage has, so the trait can be grown on its own and written to a file before the genome run that reads it starts. Two ordinary commands, in order.

**Joining** simulates two levels **at once**, because that ordering is no longer available. Suppose a trait drives **speciation** itself: large-bodied lineages split twice as fast as small ones. You cannot grow the trait first, because a trait is grown *along a tree* and the tree is what this trait shapes. Nor the tree first, because its branching rate needs a trait value that does not exist yet. Neither can be finished before the other starts, so neither can be a file handed on. They are grown together, in one run whose Gillespie races speciation, extinction and trait change against one another, each event reading the other level's current state.

Both of these chapters turn on one question:

> **Can the driver be grown first, on its own, and handed over?**

If yes, it is **conditioning**, and this chapter. If no, it is **joining**, and the next one. Underneath, both are the same single mechanism — a modifier, `mod.DrivenBy` — and only the `source` differs.

```python
from zombi2.rates import modifiers as mod

# conditioned
loss  = 0.25 * mod.DrivenBy("trait_events.tsv", {"cave": 4.0, "surface": 1.0})

# joint
birth = 1.0  * mod.DrivenBy("trait",            {"small": 1.0, "large": 2.0})
```

## The mechanism

`mod.DrivenBy` takes two things: a `source` and a `mapping`. The `source` we have just split into file versus live level. The `mapping` is the other half, and it answers a separate question: once you know the driver's value on a lineage, what factor does the rate get multiplied by? It comes in three shapes:

- **Table** — a discrete driver becomes a dict, one factor per state: `{"cave": 4.0, "surface": 1.0}`. Any state you leave out keeps its rate unchanged, and states are matched by their written form, so an integer-labelled trait still finds its entry.
- **Curve** — a numeric driver becomes a function: `lambda n: math.exp(0.05 * n)`, a rate that rises smoothly with a lineage's gene count. Pass `bound=` to cap the factor when the driver has no ceiling of its own.
- **Scalar** — a single log-link coefficient, `exp(strength · value)`, the natural response when the driver is already a 0/1 indicator or one continuous covariate. `Scalar(0.0)` is the null: factor 1 everywhere.

The `source` says where the driver lives; the `mapping` says how its value is read. Whatever the shape, the factor it returns is dimensionless and non-negative, because it is going to multiply a rate.

The canonical case is **a trait driving gene gain or loss** (Traits → Genomes) — lineages in the dark lose genes faster than lineages in the light — and it takes two runs:

```python
from zombi2 import species, traits, genomes

tree = species.simulate_species_tree(birth=1.0, death=0.3, n_extant=20, seed=1)
# 1. grow the driver: a habitat trait down the species tree
habitat = traits.simulate_discrete(tree, states=["cave", "surface"], switch=0.1, seed=1)

# 2. grow the genomes, with loss reading the habitat on each lineage
genomes.simulate_genomes_family(tree,
    loss = 0.25 * mod.DrivenBy(habitat, {"cave": 4.0, "surface": 1.0}),
    duplication=0.2, origination=0.5, seed=2)
```

The `source` here is the grown `TraitsResult` itself. That is the in-memory shortcut for the file: it is still conditioning, still two runs in order, but with no `write` and re-read in between. Hand it a filename instead and nothing else changes:

```python
habitat.write("out/", outputs=("events",))          # writes out/trait_events.tsv (a bare
                                                    # write() puts it where you point it)
genomes.simulate_genomes_family(tree,
    loss = 0.25 * mod.DrivenBy("out/trait_events.tsv", {"cave": 4.0, "surface": 1.0}),
    duplication=0.2, origination=0.5, seed=2)
```

What a trait can drive fits in six rows:

| The driver | What it drives | Written like this | Mapping |
|---|---|---|---|
| a discrete trait | `loss`, `duplication`, `origination`, `transfer` — the gene-family rates, at the **family** resolution | `loss = 0.25 * mod.DrivenBy(source, {…})` | Table |
| a discrete trait | **every rate** of an **ordered** run — those four, plus `inversion`, `transposition`, `translocation` and the chromosome tier | `inversion = 0.3 * mod.DrivenBy(source, {…})` | Table |
| a discrete trait | **every rate** of a **nucleotide** run — the same list, on a genome measured in base pairs | `inversion = 1.5 * mod.DrivenBy(source, {…})` | Table |
| a discrete trait | **every extent** of an **ordered** or **nucleotide** run — how much an event takes, rather than how often one happens | `loss_extent = 150 * mod.DrivenBy(source, {…})` | Table |
| a discrete trait | `transfer_to` — which lineage a transfer lands on, at **every** resolution | `transfer_to = mod.DrivenBy(source, {…})` | Table, or Between (reads the donor too) |
| a discrete or continuous trait | `substitution` — how fast the sequences inside a gene evolve, at the **sequences** level | `substitution = 0.05 * mod.DrivenBy(source, {…})` | Table, Curve or Scalar |

The last row is the one level further down, and Chapter 7 covers it beside the clocks it sits with: the factor multiplies the substitution rate on each species branch, so a lineage's state sets how fast its genes' sequences change. It composes with either lineage clock.

The rate rows and the extent row are the pair worth holding apart. Driving a **rate** makes a lineage delete more often; driving an **extent** makes each deletion bigger. Both are ways of saying "this lineage sheds more", they are different processes, and set together they multiply.

An extent's unit is set by the resolution: genes at ordered, base pairs at nucleotide. Write a driven extent in Python — the `--*-extent` flags take a plain number.

`source` in both rows is the grown `TraitsResult`, or the path to the `trait_events.tsv` it wrote — the trait event log, which a driven run replays against the shared tree.

The driver file is the trait event log from Chapter 8: an `initial` row for the state at t=0, then every switch with its time. Replayed against the shared species tree, that rebuilds each branch's constant stretches, which is what lets the genome engine step its Gillespie at every switch — so a lineage that changes habitat halfway down a branch loses genes at one rate before the switch and another after it. For a **discrete** driver the driven rate follows the trait exactly, not as a per-branch average. A **continuous** driver has no switches to step at, so it is approximated instead: the branch is cut into stretches of at most `step` time units (`DrivenBy(..., step=)`; by default 1% of the tree's height) and read at each stretch's midpoint.

### Two ways a trait can drive transfer

A transfer joins two lineages, a donor and a recipient. A trait can drive either end, and the two are different models. The second row of the table is the recipient end.

Driving the **rate** drives the donor. `transfer = 0.1 * mod.DrivenBy(competence, {"competent": 3.0, "normal": 1.0})` makes a competent lineage donate three times as often as a normal one. That changes how much horizontal transfer happens in the run.

Driving `transfer_to` drives the recipient. `transfer_to = mod.DrivenBy(competence, {"competent": 3.0, "normal": 1.0})` makes a competent lineage three times likelier than a normal one to be the lineage a transfer lands on. That changes no rate at all. The same transfers happen; they go somewhere else.

The two expressions look alike, but their numbers mean different things. In a rate, the number is a multiplier: it multiplies the rate of the lineage it is read on. In `transfer_to` it is a **weight**: the engine reads it on every lineage alive at that instant, and the recipient is drawn in proportion. Five candidates at weight 1 and five at weight 2 send two thirds of the transfers to the weight-2 group, because ten of the fifteen units of weight are theirs. Weights are normalised, so doubling all of them changes nothing.

That is why `transfer_to` takes the modifier on its own, with no number in front of it. A rate has a base, `0.1` per copy per unit time; a weight does not. Writing `transfer_to = 1.0 * mod.DrivenBy(...)` is an error.

A weight of 0 means the lineage cannot receive, which is often the point: only a competent lineage takes DNA up. That has one consequence worth stating plainly. If at some instant every candidate weighs 0, the transfer has nowhere to land, so it does not happen. While no eligible recipient is alive, the run's transfer rate is 0.

Driving the rate and driving `transfer_to` are independent, and a run may use either or both:

```python
competence = traits.simulate_discrete(tree, states=["competent", "normal"],
                                      switch=0.3, seed=1)

genomes.simulate_genomes_family(tree,
    transfer    = 0.1 * mod.DrivenBy(competence, {"competent": 3.0, "normal": 1.0}),
    transfer_to =       mod.DrivenBy(competence, {"competent": 3.0, "normal": 1.0}),
    initial_families=10, seed=2)
```

#### The recipient weight can read the donor too

The `transfer_to` above reads the driver on the *recipient* only. So it can say "competent lineages take DNA up more often", but not "genes move between two states and not within them" — the weight on a candidate does not know who is donating. A **`Between`** kernel closes that gap. In place of one factor per recipient state, it gives a weight per ordered **(donor state, recipient state)** pair, and the engine reads the driver on the donor as well as the candidate:

```python
from zombi2.rates.mapping import Between

habitat = traits.simulate_discrete(tree, states=["marine", "soil"], switch=0.3, seed=1)
within = Between({("marine", "marine"): 1.0, ("soil", "soil"): 1.0}, default=0.0)

genomes.simulate_genomes_family(tree,
    transfer    = 0.5,
    transfer_to = mod.DrivenBy(habitat, within),
    initial_families=10, seed=2)
```

Every transfer now stays within one habitat: a marine donor can only reach a marine recipient, a soil donor only a soil one, because every cross-habitat pair weighs 0. The numbers are weights, read exactly as before — normalised over the candidates, so they redistribute transfers without changing how many happen. `default` (1.0) is the weight of any pair you leave out, so `Between({("marine", "soil"): 3.0})` *enriches* one direction against a baseline rather than forbidding the rest; `default=0.0` is the "only the flows I name" idiom. A `Between` is a recipient weight, never a rate: driving a rate with one is refused, because a rate has no donor to condition on.

This is the trait-driven twin of Chapter 4's `Clades`. There the groups are clades — a fact about the tree — so no other level is read; here they are an evolved trait, so one is. Same kernel, same steering; only where the groups come from differs.

Combining a driven `transfer_to` with the `"distance"` rule of Chapter 4 is not supported: `transfer_to` takes one rule.

Everything outside those rows raises an error. That is the discipline everywhere in ZOMBI2 — a modifier a level cannot honour raises an error rather than being silently dropped.

### A trait can drive another trait

Both ends of the arrow can sit at the same level, and nothing changes. The question is the same one — can the driver be grown first? — and a trait grown along a fixed tree can, so this is conditioning, written with the same `DrivenBy`. The target is `rate` for a continuous trait, its variance-rate σ², and `switch` for a discrete one:

```python
depth = traits.simulate_discrete(tree, states=["deep", "shallow"], switch=0.3, seed=3)

# a body size that diffuses four times faster in the deep
size = traits.simulate_continuous(tree, start=0.0,
    rate = 1.0 * mod.DrivenBy(depth, {"deep": 4.0, "shallow": 1.0}), seed=4)

# a diet that switches five times faster there
diet = traits.simulate_discrete(tree, states=["carnivore", "herbivore"],
    switch = 0.2 * mod.DrivenBy(depth, {"deep": 5.0, "shallow": 1.0}), seed=5)
```

A `switch` written per transition drives only the transitions you name — `{"carnivore->herbivore": 0.2 * mod.DrivenBy(depth, {…}), "herbivore->carnivore": 0.2}`. Either way the driver is read wherever it changes rather than once per branch, so a lineage that moves into the deep halfway down a branch evolves at one rate before the move and another after it.

### What can be conditioned, and what cannot yet

**What can drive**: a trait, discrete or continuous; a named gene family, present or absent; and a
module's completion, a number between 0 and 1. The last two are below, under *A gene as the driver* —
they are ordinary drivers, spelled the way a trait is, and everything in this section applies to them
equally.

**What can be driven** is the rest of this section. A trait's own rate takes a driver, as above. In a genome run conditioning works at all three resolutions. At the **family** resolution it covers the four gene-family rates in the table above. At the **ordered** resolution it covers every rate the engine has, so a trait can make a lineage rearrange its gene order more often — inversion, transposition, translocation — or reshape its karyotype more often, and a driven extent makes each rearranged block longer. At the **nucleotide** resolution it covers the same list on a genome measured in base pairs, which is where a driven `loss` becomes genome reduction as it is usually meant: a lifestyle trait shedding DNA rather than dropping family tokens. At the **sequences** level it covers `substitution`, the one rate that level has.

`transfer_to` sits outside that per-resolution list because it is not a rate. It is the choice slot, and it works at all three resolutions: the same four rules, the same kernel, the same weights. What a transfer moves differs — one copy, a block of genes, an arc of DNA — and who receives it does not.

What is not implemented yet:

- **`ByFamily` and `DrivenBy` cannot be set in the same run** at the family or ordered resolution. One weights lineages by a driver, the other weights the segment by what it covers, so combining them means weighting by the product. `family_speed` counts as a `ByFamily` here. Use one or the other.
- **A sequence cannot drive anything.** The arrow runs one way: a trait or a gene drives `substitution`, but nothing reads a sequence back out.
- **A nucleotide run cannot be the driver.** A family or ordered run offers `presence` and `completion`; the nucleotide one names its genes differently and has neither yet. It can still be the *target*.

These are limits of the implementation, not of the model — the rate grammar (`SPEC §5`) is the same everywhere, and each engine gains a modifier when its own code learns to read it. Until then, a driven rate an engine cannot honour raises rather than being silently dropped.

Notice too that conditioning **folds into the target level's own command**. There is no conditioning command and no object to build; you grow the driver, then make an ordinary genome run whose `loss` happens to be `DrivenBy` instead of a bare number. That holds on the command line as well, where the rate keeps its written form:

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

A trait target reads that same file the same way. Only the continuous `--rate` takes an expression here: `--switch` is a bare number, so a driven switch rate is written in Python.

```bash
zombi2 traits size/ --from out/ --kind continuous --start 0.0 --seed 4 \
    --rate "1.0 * DrivenBy('out/traits/trait_events.tsv', {'cave': 4.0, 'surface': 1.0})"
```

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
