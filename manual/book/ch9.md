# Conditioning

The book has so far run the levels one at a time: species tree, then genomes, then sequences. They already depend on the tree, but sometimes you want more than that background dependency. There are two ways, conditioning in this chapter and joining in the next.

**Conditioning** takes a value that has already evolved and lets it drive the run that comes next. Three examples:

- **Cave fish lose their eyes.** A habitat trait, cave or surface, has already evolved down the species tree. Wherever a lineage sits in the dark, its genes are lost four times faster than on the surface.
- **Endosymbionts shed their genomes.** A lifestyle trait, free-living or host-restricted, drives gene loss across the board, so the lineages that moved inside a host are the ones that end up with the small genomes.
- **Competent bacteria pick genes up.** A trait for natural competence raises the rate at which new gene families appear in a lineage, so the trait leaves its mark on gene *gain* rather than on loss.

In each case one evolved value is read by the run that comes next. The value read is the **driver**; what reads it there is the **target**. They are not interchangeable, and the asymmetry is the point: a driver is a *value* that already varies from lineage to lineage, a habitat state or a gene count, while a target is whatever that value is attached to, most often a rate, multiplied by the factor the driver's value picks out. The arrow runs one way, and nothing flows back.

![The shape of a conditioned run, and the four words the rest of this chapter uses. The **driver** is a level already simulated, a habitat trait here, its two states shown below it. The **target** is what the factor is attached to, a rate here, in the run that comes next. The **modifier** is what joins them, `DrivenBy`, and what it carries is the **mapping**: one multiplier per state of the driver, so a branch's habitat sets that branch's loss rate. The driver is finished and written to a file before the second run starts, which is what lets this be two ordinary commands.](figures/conditioning_print.png){width=95%}

That one-wayness is what makes conditioning cheap. The habitat is unaffected by how many genes a lineage has, so the trait can be grown on its own and written to a file before the genome run that reads it starts. Two ordinary commands, in order.

## The mechanism

A written rate shows all four parts of Figure 9.1 at once. In

```
loss = 0.25 * mod.DrivenBy(habitat, {"cave": 4.0, "surface": 1.0})
```

`loss` is the target, `0.25` its base rate, `DrivenBy` the modifier, `habitat` the driver and the dict the mapping. Those are the argument names too: `DrivenBy(driver, mapping)`.

The driver is either a finished result object or the path to the file that result wrote, the same model either way, since the file is only how a second process reads what the first one grew. The mapping says what the driver's value does, and which shape you need follows from the driver rather than from what it drives.

![The mappings. A driver with named states takes a `Table`, `{"cave": 4.0, "surface": 1.0}`, where any state you leave out keeps its rate unchanged. A driver whose value is a number takes a `Curve`, any function you write, with `bound=` to cap it, or a `Scalar`, the log-link `exp(strength · value)` for a driver that is already an indicator or a single covariate. `Between` is the odd one, reading the driver at both ends of a transfer, and only `transfer_to`, the third kind of target below, can use it.](figures/conditioning_mappings_print.png){width=92%}

A bare dict is read as a Table and a bare function as a Curve, so neither of the common two needs a constructor. Whatever the shape, the factor is dimensionless and non-negative. The canonical case is a trait driving gene loss: lineages in the dark lose genes faster than lineages in the light. It is two ordinary runs, in order.

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

Passing `habitat` hands the result over in memory; passing `"out/trait_events.tsv"` after `habitat.write("out/", outputs=("events",))` hands over the same thing on disk.

Either way the driver is read **wherever it changes**, not once per branch. The event log gives each branch's constant stretches, and the genome engine steps its Gillespie at every switch, so a lineage that changes habitat halfway down a branch loses genes at one rate before the switch and another after it. A discrete driver is therefore followed exactly. A continuous one has no switches to step at, so it is approximated: the branch is cut into stretches of at most `step` time units (`DrivenBy(..., step=)`, by default 1% of the tree's height) and read at each stretch's midpoint.

What a factor can be attached to comes in three kinds, and only the first is a rate. **How often**: driving `loss` makes a lineage delete more often. **How much**: driving `loss_extent` makes each deletion bigger. **A choice**: driving `transfer_to` changes which lineage a transfer lands on. That triple is not Chapter 5's three questions about a segmental event: *where* an event starts is drawn by the engine and takes no modifier. The first two say "this lineage sheds more" by different means, and set together they multiply. An extent's unit is set by the resolution, genes at ordered and base pairs at nucleotide, and a driven extent is written in Python, since the `--*-extent` flags take a plain number.

### Driving who receives a transfer

`transfer_to` sits outside everything above, and it is worth saying how. It is the one target whose number is not a multiplier: it does not say how often or how much anything happens, only which lineage a transfer that has already fired lands on. Its numbers therefore mean something else, and it takes a mapping of its own that no rate can use.

In a rate the factor multiplies the rate of the lineage it is read on. In `transfer_to` it is a **weight**: the engine reads it on every lineage alive at that instant and draws the recipient in proportion, so five candidates at weight 1 and five at weight 2 send two thirds of the transfers to the weight-2 group. Weights are normalised, so doubling them all changes nothing, which is why `transfer_to` takes the modifier on its own with no base in front of it. `transfer_to = 1.0 * mod.DrivenBy(...)` is an error. A weight of 0 means the lineage cannot receive, and when every candidate weighs 0 the transfer has nowhere to land, so it does not happen.

A plain table reads the driver on the recipient only, so it can say "competent lineages take DNA up more often" but not "genes move between two habitats and not within them". A **`Between`** kernel closes that gap, giving a weight per ordered (donor state, recipient state) pair and reading the driver on the donor as well as on the candidate.

```python
from zombi2.rates.mapping import Between

habitat = traits.simulate_discrete(tree, states=["marine", "soil"], switch=0.3, seed=1)
within = Between({("marine", "marine"): 1.0, ("soil", "soil"): 1.0}, default=0.0)

genomes.simulate_genomes_family(tree, transfer=0.5, initial_families=10, seed=2,
                                transfer_to=mod.DrivenBy(habitat, within))
```

Every transfer now stays within one habitat, because every cross-habitat pair weighs 0. `default` (1.0) is the weight of any pair you leave out, so `Between({("marine", "soil"): 3.0})` enriches one direction against a baseline instead of forbidding the rest. This is the trait-driven twin of Chapter 4's `Clades`: there the groups come from the tree, here from an evolved trait, and the kernel and the steering are the same. Driving a rate or an extent with a `Between` is refused, because both are read on one lineage and have no donor to read it on, and `transfer_to` takes one rule, so a driven `transfer_to` cannot also be `"distance"`.

## What can condition what

A conditioned pair is a **driver** and a **target**, and every level offers both. As a driver a level offers a value that varies from lineage to lineage: a family's presence, a GC content, a character's state. As a target it offers the things that value can be attached to, which come in three kinds and are not all rates.

![What each level offers, as a driver and as a target. Targets come in three kinds and only the first is a rate: **how often** an event fires, **how much** it takes, and **which one** receives a transfer. The last two are narrower than the first for reasons of the model, not of the code: a gene-family event moves exactly one gene and a substitution changes exactly one site, so there is nothing to size; and only a transfer has a recipient to choose. The two genome columns share a cell wherever the resolutions agree.](figures/conditioning_parts_print.png){width=88%}

Any level can be a driver and any level can be a target, so the question is which pairs are models. Eleven are.

![What can condition what. Rows drive, columns are driven, and the eleven numbered pairs are the conditioned models listed below. The five shaded cells are the pairs that are not models. Four of them would need two genomes for one lineage. The fifth is a sequence driving a gene family: a sequence lives inside a gene, so a genome reading one back would condition a run on its own output. Gene families are different, because a genome is a *set* of families, and two runs are two parts of one genome. The inset box marks a level conditioning itself, which is ordinary conditioning done twice on one level.](figures/conditioning_map_print.png){width=95%}

**1. A gene family drives a gene family.** Two family runs, the second reading the first: a mobile element makes transfer likelier for the rest of the genome, written `transfer = 0.1 * mod.DrivenBy(first.presence("IS1"), {"present": 5.0, "absent": 1.0})`.

**2. A gene family drives a sequence.** Presence or completion setting the substitution rate: lose the repair gene and evolve faster. Chapter 7 covers this beside the clocks it composes with.

**3. A gene family drives a trait.** Carry the toxin family and turn pathogenic: `switch = 0.1 * mod.DrivenBy(g.presence("tox"), {"present": 8.0, "absent": 1.0})`. Only families you **named** with `family_names=` can be asked for, because a family that arose during the run has an id but nothing stable to call it by. Presence is exact and changes during a branch, so a lineage that loses its last copy halfway along one is `present` before that instant and `absent` after it.

**4. An ordered or nucleotide genome drives a sequence.** As **2**, with the genome at a resolution that also has coordinates. The driver reads the same; only what the genome run itself tracked differs.

**5. An ordered or nucleotide genome drives a trait.** As **3**, and the driver may also be a **module**, a named group of families such as a pathway or an operon, declared with `modules=` and read as `g.completion("flagellum")`. Completion is a fraction rather than a state, so it takes a Curve, and that is where a threshold belongs: `lambda f: 8.0 if f > 0.8 else 1.0` reads "eight times faster once four fifths of it is there". A fraction rather than a yes-or-no for a measured reason: under independent loss the chance that *every* family survives falls off geometrically with the module's size, and on a 200-tip tree a module of three was complete at 189 tips where one of six was complete at none, so a complete/incomplete driver would be a constant for all but the smallest modules.

**6. A sequence drives a sequence.** One gene's sequence driving another's substitution rate, which is compensatory evolution.

**7. A sequence drives a trait.** The same driver on a trait's rate.

**8. A trait drives a gene family.** The cave-fish run above. All four rates take a driver, `transfer` included, so a competence trait can make a lineage donate three times as often, and `transfer_to` can make it three times likelier to receive.

**9. A trait drives an ordered or nucleotide genome.** The same over eleven rates, and over the extents besides, so a habitat can decide both how often a lineage rearranges its gene order and how much each rearrangement takes.

**10. A trait drives a sequence.** The substitution rate, covered in Chapter 7 beside the clocks it sits with: the factor multiplies the rate on each species branch and composes with either lineage clock.

**11. A trait drives a trait.** One character driving another's `rate` or `switch`, and the plainest case of a level conditioning itself. A body size that diffuses four times faster in the deep is `rate = 1.0 * mod.DrivenBy(depth, {"deep": 4.0, "shallow": 1.0})`; a `switch` written per transition drives only the transitions you name.

Two things are refused, for different reasons. **`Drawn(per='family')` and `DrivenBy` cannot be set in the same run** at the family or ordered resolution, because one weights lineages by a driver and the other weights the segment by what it covers, so combining them means weighting by the product. That is a limit of the code. **A genome cannot be driven by a sequence**, and that is a limit of the model: a sequence was grown along the gene trees the genome run produced, so a genome reading it back would condition a run on its own output. Those two cells are shaded on the map for that reason, and the pair can only be joined (Chapter 10). A sequence therefore drives what comes *after* it, a trait on the same species tree or a further sequence run. Either way a driven rate an engine cannot honour raises, naming the reason, rather than being silently dropped.

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

Each conditioned run writes a `conditioned_on` file naming the levels it read, so re-running the trait afterwards refuses rather than leaving the runs beneath it silently mismatched. Pass `--force` to re-run it and clear them. The guard works **within one run directory**: it is the directory that holds the record, so a driver and a target written to two different directories (with `--from`) are not linked, and re-running the driver there will not warn you. It also works **between** levels and not within one: a trait driving a second trait sits under `traits/` on both sides, so the record is written but re-running the driver does not invalidate the run that read it.

Both halves of transfer take that same text, the rate with a base number in front of it and the recipient weight without one.

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

Only a trait driver has a written form, a path to its event log, so every other cell is Python only. That is the asterisk on the map.

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
from zombi2.rates.mapping import Curve
from zombi2.sequences import simulate_sequences
from zombi2.sequences.substitution_models import hky85, lg

g = genomes.simulate_genomes_family(tree, initial_families=20, duplication=0.1, loss=0.15, seed=9)
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

Conditioning adds no format. A conditioned run writes the target level's own files, listed in that level's own chapter, and adds one record: `conditioned_on`, naming the levels this run reads through `DrivenBy`, whether on a rate, an extent or `transfer_to`, one per line. When the driver was grown in the same run directory its own files sit beside them, so the pairing that produced the pattern is kept on disk. Appendix B gives the columns and the formats.
