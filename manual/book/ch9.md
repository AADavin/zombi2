# Conditioning and joining

The four levels of the book have so far discussed how to run things individually: first species tree, then genomes, then sequences. The levels are not independent in the sense that all depend on the species tree. However, in some cases we want to crank up that dependency in a way that the species tree is more than a simple background. There are two main ways to do this: conditioning and joining.

When we condition, we take a value of one level and make it drive a rate in a different level. This is more clearly seen with examples:

- **Cave fish lose their eyes.** A habitat trait, cave or surface, has already evolved down the species tree. Wherever a lineage sits in the dark, its genes are lost four times faster than on the surface.
- **Endosymbionts shed their genomes.** A lifestyle trait, free-living or host-restricted, drives gene loss across the board, so the lineages that moved inside a host are the ones that end up with the small genomes.
- **Competent bacteria pick genes up.** A trait for natural competence raises the rate at which new gene families appear in a lineage, so the trait leaves its mark on gene *gain* rather than on loss.

In each case one level's value has been read by another level's rate. The value doing the reading-from is the **driver**; the rate doing the reading is the **target**. The two are not interchangeable, and the asymmetry is the whole point: a driver is a *value* that already varies from lineage to lineage — a habitat state, a gene count — while a target is a *rate*, a "how often", which gets multiplied by a factor that the driver's value picks out. The arrow runs one way, from the driver's value to the target's rate, and the driver never notices it is being read.

That one-wayness is what makes conditioning cheap. Because the habitat does not care how many genes a lineage has, the trait can be grown completely, on its own, and written to a file before the genome run that reads it ever starts. Two ordinary commands, in order.

When we join, we simulate **simultaneously** two levels, because that ordering is no longer available to us. Suppose a trait drives not gene loss but **speciation** itself: large-bodied lineages split twice as fast as small ones. Now you cannot grow the trait first, because a trait is grown *along a tree* and the tree is precisely what this trait is busy shaping. Nor can you grow the tree first, because its branching rate needs a trait value that does not exist yet. Neither level can be finished before the other starts, so neither can be a file handed to the next command. They have to be grown together, in a single run whose Gillespie races speciation, extinction and trait change against one another, each event reading the other level's current state.

So the whole chapter turns on one question:

> **Can the driver be grown first, on its own, and handed over?**

If yes, it is **conditioning**: two runs, and the coupling's `source` is a file. If no, it is **joining**: one run, and the `source` is the name of a level growing beside it. Underneath, both are the same single mechanism — a modifier, `mod.DrivenBy` — and only the `source` differs.

```python
loss  = 0.25 * mod.DrivenBy("trait_driver.tsv", {"cave": 4.0, "surface": 1.0})   # conditioned
birth = 1.0  * mod.DrivenBy("trait",            {"small": 1.0, "large": 2.0})    # joint
```

## Conditioning

`mod.DrivenBy` takes two things: a `source` and a `mapping`. The `source` we have just split into file versus live level. The `mapping` is the other half, and it answers a separate question: once you know the driver's value on a lineage, what factor does the rate get multiplied by? It comes in three shapes:

- **Table** — a discrete driver becomes a dict, one factor per state: `{"cave": 4.0, "surface": 1.0}`. Any state you leave out keeps its rate unchanged, and states are matched by their written form, so an integer-labelled trait still finds its entry.
- **Curve** — a numeric driver becomes a function: `lambda n: math.exp(0.05 * n)`, a rate that rises smoothly with a lineage's gene count. Pass `bound=` to cap the factor when the driver has no ceiling of its own.
- **Scalar** — a single log-link coefficient, `exp(strength · value)`, the natural response when the driver is already a 0/1 indicator or one continuous covariate. `Scalar(0.0)` is the null: factor 1 everywhere.

The `source` says where the driver lives; the `mapping` says how its value is read. Whatever the shape, the factor it returns is dimensionless and non-negative, because it is going to multiply a rate.

Conditioning goes in exactly one direction today: **a trait drives gene gain or loss** (Traits → Genomes). The cave example is the canonical one — lineages in the dark lose genes faster than lineages in the light — and it takes two runs:

```python
# 1. grow the driver: a habitat trait down the species tree
habitat = traits.simulate_discrete(tree, states=["cave", "surface"], switch=0.1, seed=1)

# 2. grow the genomes, with loss reading the habitat on each lineage
genomes.simulate_genomes_unordered(tree,
    loss = 0.25 * mod.DrivenBy(habitat, {"cave": 4.0, "surface": 1.0}),
    duplication=0.2, origination=0.5, seed=2)
```

The `source` here is the grown `TraitsResult` itself. That is the in-memory shortcut for the file: it is still conditioning, still two runs in order, but with no `write` and re-read in between. Hand it a filename instead and nothing else changes:

```python
habitat.write("out/", outputs=("driver",))          # writes out/trait_driver.tsv
genomes.simulate_genomes_unordered(tree,
    loss = 0.25 * mod.DrivenBy("out/trait_driver.tsv", {"cave": 4.0, "surface": 1.0}),
    duplication=0.2, origination=0.5, seed=2)
```

That is the whole of conditioning today, and it fits in two rows:

| The driver | What it drives | Written like this | Mapping |
|---|---|---|---|
| a discrete trait | `loss`, `duplication`, `origination`, `transfer` — the rates of an unordered genome run | `loss = 0.25 * mod.DrivenBy(source, {…})` | Table |
| a discrete trait | `transfer_to` — which lineage a transfer lands on | `transfer_to = mod.DrivenBy(source, {…})` | Table |

`source` in both rows is the grown `TraitsResult`, or the path to the `trait_driver.tsv` it wrote.

The driver file is the one from Chapter 8: a discrete trait's stochastic character map, cut into the constant stretches of each branch. That segmentation is what lets the genome engine step its Gillespie at every switch, so a lineage that changes habitat halfway down a branch loses genes at one rate before the switch and another after it. The coupling is exact, not a per-branch average.

### Two ways a trait can drive transfer

A transfer joins two lineages, a donor and a recipient. A trait can drive either end, and the two are different models. The second row of the table is the recipient end.

Driving the **rate** drives the donor. `transfer = 0.1 * mod.DrivenBy(competence, {"competent": 3.0, "normal": 1.0})` makes a competent lineage donate three times as often as a normal one. That changes how much horizontal transfer happens in the run.

Driving `transfer_to` drives the recipient. `transfer_to = mod.DrivenBy(competence, {"competent": 3.0, "normal": 1.0})` makes a competent lineage three times likelier than a normal one to be the lineage a transfer lands on. That changes no rate at all. The same transfers happen; they go somewhere else.

The two expressions look alike, but their numbers mean different things. In a rate, the number is a multiplier: it multiplies the rate of the lineage it is read on. In `transfer_to` it is a **weight**: the engine reads it on every lineage alive at that instant, and the recipient is drawn in proportion. Five candidates at weight 1 and five at weight 2 send two thirds of the transfers to the weight-2 group, because ten of the fifteen units of weight are theirs. Weights are normalised, so doubling all of them changes nothing.

That is why `transfer_to` takes the modifier on its own, with no number in front of it. A rate has a base, `0.1` per copy per unit time; a weight does not. Writing `transfer_to = 1.0 * mod.DrivenBy(...)` is an error.

A weight of 0 means the lineage cannot receive, which is often the point: only a competent lineage takes DNA up. That has one consequence worth stating plainly. If at some instant every candidate weighs 0, the transfer has nowhere to land, so it does not happen. While no eligible recipient is alive, the run's transfer rate is 0.

The two couplings are independent, and a run may use either or both:

```python
competence = traits.simulate_discrete(tree, states=["competent", "normal"],
                                      switch=0.3, seed=1)

genomes.simulate_genomes_unordered(tree,
    transfer    = 0.1 * mod.DrivenBy(competence, {"competent": 3.0, "normal": 1.0}),
    transfer_to =       mod.DrivenBy(competence, {"competent": 3.0, "normal": 1.0}),
    initial_families=10, seed=2)
```

Combining a driven `transfer_to` with the `"distance"` rule of Chapter 4 is not supported: `transfer_to` takes one rule.

Everything outside those two rows raises an error. Every rate of the ordered resolution of Chapter 5 is refused, and so are the rates of sequence evolution and of a trait run. The ordered resolution refuses a driven `transfer_to` as well. The driver has to be a discrete trait: only a discrete trait carries the character map that cuts a branch into constant segments, so a continuous trait is refused as a driver. That is the discipline everywhere in ZOMBI2 — a modifier a level cannot honour raises an error rather than being silently dropped.

Notice too that the coupling **folds into the target level's own command**. There is no separate coupling step and no coupling object to build; you grow the driver, then make an ordinary genome run whose `loss` happens to be `DrivenBy` instead of a bare number. That holds on the command line as well, where the rate keeps its written form:

```bash
# 1. a species tree
zombi2 species --birth 1 --death 0.3 --n-extant 20 --seed 1 -o out/

# 2. the driver: a habitat trait, writing the driver file
zombi2 traits --kind discrete -t out/species/species_complete.nwk \
    --states cave,surface --switch 0.1 --seed 1 -o out/ --write values tree driver

# 3. the target: genomes whose loss reads that trait
zombi2 genomes -t out/species/species_complete.nwk \
    --loss "0.25 * DrivenBy('out/trait_driver.tsv', {'cave': 4.0, 'surface': 1.0})" \
    --duplication 0.2 --origination 0.5 --seed 2 -o out/
```

Both halves of transfer take that same text: the rate with a base number in front of it, the recipient weight without one.

```bash
# the driver: a competence trait, into its own directory
zombi2 traits --kind discrete -t out/species/species_complete.nwk \
    --states competent,normal --switch 0.3 --seed 1 -o comp/ --write driver

zombi2 genomes -t out/species/species_complete.nwk --initial-families 10 \
    --transfer    "0.1 * DrivenBy('comp/trait_driver.tsv', {'competent': 3.0, 'normal': 1.0})" \
    --transfer-to "DrivenBy('comp/trait_driver.tsv', {'competent': 3.0, 'normal': 1.0})" \
    --seed 2 -o comp_genomes/
```

## Joining

When the driver cannot be grown first, there is nothing to write to a file, because the file would have to be written onto a tree that does not exist yet. Instead both levels are grown in one Gillespie that races every kind of event against every other: a speciation event reads the current trait state to set its rate, a trait-change event evolves the trait on the tree as it has grown so far, and a speciation hands the parent's state down to both daughters. Out come both levels at once. Because these drivers only change *at* events, the rate is constant between them and the race is exact — no thinning, no approximation.

Version 1 ships two joint pairs, and in both a level reaches back into the species tree, so the tree itself is an output rather than an input:

**A trait drives speciation** — `P(Species, Traits)`. A body-size state makes large lineages speciate twice as fast as small ones. The trait enters as a **process spec**, `traits.discrete(...)`, rather than as a finished run: it is a description of a process to be grown with the tree, not a result:

```python
from zombi2 import joint, traits, genomes
from zombi2.rates import modifiers as mod

joint.simulate_joint(
    birth = 1.0 * mod.DrivenBy("trait", {"small": 1.0, "large": 2.0}),
    death = 0.2,
    trait = traits.discrete(states=["small", "large"], switch=0.1),
    n_extant = 100, seed = 1)
```

Driving `death` with the same modifier makes extinction state-dependent too, and a birth *and* death that both read the trait is the model the literature calls BiSSE — with more than two states, MuSSE.

**Gene content drives speciation** — `P(Species, Genomes)`. The presence of a key gene, a toxin or a transporter, lifts a lineage's speciation rate, and the genome and the tree grow together. The genome enters as a process spec in the same way, and the family whose presence does the driving has to be seeded at the root:

```python
joint.simulate_joint(
    birth  = 1.0 * mod.DrivenBy("genomes:toxin", {"present": 1.8, "absent": 1.0}),
    death  = 0.2,
    genome = genomes.unordered(duplication=0.2, loss=0.25, origination=0.5,
                               families=["toxin"]),
    n_extant = 100, seed = 1)
```

The live `source` names either the family, `"genomes:toxin"`, or the lineage's total gene count, `"genomes:count"` — and a count is a number, so that is where a **Curve** earns its place: `mod.DrivenBy("genomes:count", lambda n: math.exp(0.02 * n))` makes gene-rich lineages speciate faster along a smooth response rather than a lookup table.

The three live sources in full:

| The driver | The rates it can drive | `source` | Mapping |
|---|---|---|---|
| a discrete trait | `birth`, `death` | `"trait"` | Table |
| a named family's presence | `birth`, `death` | `"genomes:<family>"` | Table on `present`/`absent` |
| a lineage's gene count | `birth`, `death` | `"genomes:count"` | Curve, or Scalar |

Give one driver per run, `trait=` or `genome=`, not both, and drive `birth`, `death`, or both with it.

Stop the run at a size with `n_extant=` or at an age with `total_time=`, exactly as in Chapter 3, and give one or the other. What comes back is a `JointResult` carrying **both** grown levels: `.species` always, and then either `.trait` or `.genome`, the same result objects the standalone commands return. They share one `complete_tree`, because there was only ever one tree — the one they grew between them. Joint runs are Python-only for now; the conditioned half of the chapter is the part that has a command line.

### Not everything that looks like a connection is one

One distinction keeps this chapter from swallowing material that belongs elsewhere. A trait that jumps at a speciation event, or a genome that changes only at splits, looks like a coupling to the species level, but it is not one. It is the level reading the tree it *already lives on*, which every level does for free. A cladogenetic trait shift and a punctuational burst of gene change are options of the trait's and the genome's own models, and they stay in Chapters 8 and 6 respectively. A coupling, in this chapter's sense, always reads a *different* level. When in doubt, ask whether the rate reads a value that some *other* level produced; if it only reads the tree, it is not here.

## Literature

The state-dependent models arrive under a wall of acronyms, and a reader who wants "a BiSSE model" should be able to find the door. The names live here, in one table, and organise nothing else in the chapter.

| What it does | ZOMBI2 | From the literature |
|---|---|---|
| a binary trait drives speciation (and extinction) | `simulate_joint`: `birth`, `death` `= … * mod.DrivenBy("trait", {…})` | BiSSE |
| a multi-state trait drives speciation | the same, with more states in the Table | MuSSE |

## Outputs

A conditioned run writes what any ordinary level run writes — its genome or trait output — plus the **driver file** that fed it, so the pairing that produced the pattern is kept on disk alongside the result. A joint run writes **both** levels from one call: the grown species tree (`species_complete.nwk`, `species_extant.nwk`, `species_events.tsv`) together with the trait it grew (`trait_values.tsv`, `trait_changes.tsv`, `trait_tree.nwk`) or the genomes it grew (`genome_events.tsv`, `profiles.tsv`), each in the format it would have had from its own command. Because a joint run grows the tree, the tree it writes is a *complete* tree in the sense of Chapter 3, with the extinct lineages that shaped the trait or gene distribution still in place — which matters here more than anywhere, since those are exactly the lineages whose fate the coupling decided. The full list of files lives in Appendix B.
