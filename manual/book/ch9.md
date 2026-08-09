# Conditioning

The book has so far run the levels one at a time: species tree, then genomes, then sequences, then traits. Each already depends on the tree. Sometimes you want more than that — you want a value that has *evolved* to steer the run that comes next.

Take olfactory genes. A habitat trait, aquatic or terrestrial, evolves down the species tree. Wherever a lineage is aquatic, it loses those genes four times faster. Two more of the same shape:

- **Endosymbionts shed their genomes.** A lifestyle trait, free-living or host-restricted, drives gene loss across the board, so the lineages that moved inside a host end up with the small genomes.
- **Competent bacteria pick genes up.** A trait for natural competence raises the rate at which new families arrive in a lineage, so the trait leaves its mark on gene *gain* rather than on loss.

What makes these **conditioning** is that the driver can be finished first. A lineage's habitat does not depend on how many genes it has, so the trait can be grown on its own and written to a file before the genome run starts. Two ordinary commands, in order. The next chapter is for the cases where that does not work.

## The four words

![The shape of a conditioned run, and the four words the rest of this chapter uses. The **driver** is a level already simulated, a habitat trait here, its two states shown below it. The **target** is what the factor is attached to, a rate here, in the run that comes next. The **verb** is what joins them, `scaled_by` here, and what it carries is the **mapping**: one multiplier per state of the driver, so a branch's habitat sets that branch's loss rate. The driver is finished and written to a file before the second run starts, which is what lets this be two ordinary commands.](figures/conditioning_print.png){width=95%}

The **driver** is the value that is read: a habitat state, a gene's presence, a GC content. The **target** is what reads it — usually a rate, in the run that comes next. The **verb** is what joins them, and the **mapping** is what the verb carries: it turns the driver's value into a number.

Driver and target are not interchangeable. A driver is a value that already varies from lineage to lineage; a target is a parameter of the next run. The arrow goes one way and nothing comes back.

## Writing it

All four words appear in one line:

```
loss = PerCopy(0.25).scaled_by(habitat, {"aquatic": 4.0, "terrestrial": 1.0})
```

`loss` is the target, `PerCopy(0.25)` its base rate in its scope, `scaled_by` the verb, `habitat` the driver, and the dict the mapping. Those last two are the argument names as well: `scaled_by(driver, mapping)`.

Written out, the canonical case is two ordinary runs:

```python
from zombi2 import species, traits, genomes
from zombi2.params import PerCopy, PerLineage, Recipients

tree = species.simulate_species_tree(birth=1.0, death=0.3, n_extant=20, seed=1)
# 1. grow the driver: a habitat trait down the species tree
habitat = traits.simulate_discrete(tree, states=["aquatic", "terrestrial"], switch=0.1, seed=1)

# 2. grow the genomes, with loss reading the habitat on each lineage
genomes.simulate_genomes_family(tree,
    loss = PerCopy(0.25).scaled_by(habitat, {"aquatic": 4.0, "terrestrial": 1.0}),
    duplication=0.2, origination=0.5, seed=2)
```

The driver is either the finished result object, as here, or the path to the file that result wrote. It is the same model either way: the file is only how a second process reads what the first one grew. Passing `habitat` hands it over in memory; passing `"out/trait_events.tsv"` after `habitat.write("out/", outputs=("events",))` hands over the same thing on disk.

### Setting a rate instead of scaling it

A factor is a multiple of a base you had to invent. The literature often states the rate itself — *the loss rate is 1.0 in the water* — and `set_by` says that directly:

```
loss = PerCopy().set_by(habitat, {"aquatic": 1.0, "terrestrial": 0.25})
```

There is no base in front, because the driver supplies the whole number, in the rate's own units. Writing one anyway raises rather than silently discarding it. The scope still applies, so a per-copy rate set to 1.0 is 1.0 per copy. A rate carries one `set_by` and any number of ordinary factors. Three targets take it: the family and ordered genome resolutions, and a continuous trait's rate.

### A clade as the driver

One driver needs nothing grown first. `Clade` reads which named group of the tree a lineage sits in, which is a fact about the tree the run is already walking:

```
loss = PerCopy(0.2).scaled_by(Clade({"fast": ["n12", "n27"]}), {"fast": 3.0})
```

Name a clade either by a list of tips, giving the subtree below their most recent common ancestor, or by a single node id. Clades must not overlap, and a lineage in none of them is in the group `"rest"`, which a mapping may name like any other state. Membership never changes along a branch, so this is the cheapest driver there is. Chapter 4's `Clades`, plural, is a different thing: it is the `transfer_to` rule that weighs the donor's clade against the recipient's.

## The driver is read wherever it changes

A driven rate is not read once per branch. The event log gives each branch's constant stretches, and the engine steps its Gillespie at every switch, so a lineage that changes habitat halfway down a branch loses genes at one rate before the switch and another after it. A discrete driver is therefore followed exactly.

A continuous driver has no switches to step at, so it is approximated. The branch is cut into stretches of at most `step` time units — `scaled_by(..., step=)`, by default a hundredth of the tree's height — and read at each stretch's midpoint.

## Choosing the mapping

![The mappings. A driver with named states takes a `Table`, `{"aquatic": 4.0, "terrestrial": 1.0}`, where any state you leave out keeps its rate unchanged. A driver whose value is a number takes a `Curve`, any function you write, with `bound=` to cap it, or a `Scalar`, the log-link `exp(strength · value)` for a driver that is already a single measured number. `Between` is the odd one, reading the driver at both ends of a transfer, and only `transfer_to`, the third kind of target below, can use it.](figures/conditioning_mappings_print.png){width=92%}

Which mapping you need follows from the **driver**, not from what it drives. If the driver has named states, give a name to each: `{"aquatic": 4.0, "terrestrial": 1.0}`. If its value is a number, give a function of that number: `lambda x: 2.0 ** x`. A bare dict is read as a `Table` and a bare function as a `Curve`, so neither of the common two needs a constructor.

`Scalar(strength)` is a third shape for a number, giving `exp(strength · value)`. It is the form the literature usually writes when a single measured quantity enters a rate, and it is worth knowing that a `Scalar` is a `Curve` with the exponential already chosen.

Whatever the shape, the number that comes out is non-negative and has no units.

## What a driver can be attached to

![What each level offers, as a driver and as a target. Targets come in three kinds and only the first is a rate: **how often** an event fires, **how much** it takes, and **which one** receives a transfer. The last two are narrower than the first for reasons of the model, not of the code: a gene-family event moves exactly one gene and a substitution changes exactly one site, so there is nothing to size; and only a transfer has a recipient to choose. The two genome columns share a cell wherever the resolutions agree.](figures/conditioning_parts_print.png){width=88%}

A factor can be attached to three kinds of thing, and only the first is a rate.

**How often.** Driving `loss` makes a lineage delete more often.

**How much.** Driving `loss_extent` makes each deletion bigger. An extent's unit comes from the resolution — genes at ordered, base pairs at nucleotide — and a driven extent is written in Python, because the `--*-extent` flags take a plain number.

**Which one.** Driving `transfer_to` changes which lineage a transfer lands on, and is the subject of the next section.

The first two say "this lineage sheds more" by different means, and set together they multiply. This is not Chapter 5's three questions about a segmental event: *where* an event starts is drawn by the engine and takes no modifier.

## Choosing who receives a transfer

`transfer_to` sits outside everything above. Its number is not a multiplier: it does not say how often or how much anything happens, only which lineage a transfer that has already fired lands on.

So the number means something else. It is a **weight**. The engine reads it on every lineage alive at that instant and draws the recipient in proportion, so five candidates at weight 1 and five at weight 2 send two thirds of the transfers to the weight-2 group. Weights are normalised, so doubling them all changes nothing — which is why `transfer_to` is written from `Recipients()`, with no base in front of it. `transfer_to = PerCopy(1.0).weighted_by(...)` is an error. A weight of 0 means the lineage cannot receive, and when every candidate weighs 0 the transfer has nowhere to land, so it does not happen.

A plain table reads the driver on the recipient only. That says "competent lineages take DNA up more often", but it cannot say "genes move between two habitats and not within them". `Between` closes that gap: it gives a weight per ordered (donor state, recipient state) pair, so it reads the driver on the donor as well as on the candidate.

```python
from zombi2.params import Between, Recipients

habitat = traits.simulate_discrete(tree, states=["marine", "soil"], switch=0.3, seed=1)
within = Between({("marine", "marine"): 1.0, ("soil", "soil"): 1.0}, default=0.0)

genomes.simulate_genomes_family(tree, transfer=0.5, initial_families=10, seed=2,
                                transfer_to=Recipients().weighted_by(habitat, within))
```

Every transfer now stays within one habitat, because every cross-habitat pair weighs 0. `default` (1.0) is the weight of any pair you leave out, so `Between({("marine", "soil"): 3.0})` enriches one direction against a baseline instead of forbidding the rest. This is the trait-driven twin of Chapter 4's `Clades`: there the groups come from the tree, here from an evolved trait, and the steering is the same. A `Between` on a rate or an extent is refused, because both are read on one lineage and have no donor to read the other end on. And `transfer_to` takes one rule, so a driven `transfer_to` cannot also be `"distance"`.

## What can drive what

![What can condition what. Rows drive, columns are driven, and the eleven numbered pairs are the conditioned models listed below. The five shaded cells are the pairs that are not models. Three of them would need two genomes for one lineage. The other two are a sequence driving a genome: a sequence lives inside a gene, so a genome reading one back would condition a run on its own output. Gene families are different, because a genome is a *set* of families, and two runs are two parts of one genome. The inset box marks a level conditioning itself, which is ordinary conditioning done twice on one level.](figures/conditioning_map_print.png){width=95%}

Every level offers a driver and every level offers a target, so the question is which pairs are models. Eleven are, and they are numbered on the map.

| | Driver | Target | What it says |
|---|---|---|---|
| **1** | a gene family | a gene family | a mobile element makes transfer likelier for the rest of the genome |
| **2** | a gene family | a sequence | lose the repair gene and evolve faster (Chapter 7) |
| **3** | a gene family | a trait | carry the toxin family and turn pathogenic |
| **4** | an ordered or nucleotide genome | a sequence | as **2**, with coordinates in the genome run |
| **5** | an ordered or nucleotide genome | a trait | as **3**, with coordinates in the genome run |
| **6** | a sequence | a sequence | compensatory evolution: one gene's sequence sets another's rate |
| **7** | a sequence | a trait | the same driver on a trait's rate |
| **8** | a trait | a gene family | the habitat run above; all four rates take a driver |
| **9** | a trait | an ordered or nucleotide genome | eleven rates, and the extents besides |
| **10** | a trait | a sequence | the substitution rate (Chapter 7) |
| **11** | a trait | a trait | one character sets another's `rate` or `switch` |

Three of those need a word more.

**A gene family as a driver** is `g.presence("IS1")`, and only families you **named** with `family_names=` can be asked for, because a family that arose during the run has an id but nothing stable to call it by. Presence is exact and changes during a branch, so a lineage that loses its last copy halfway along one is `present` before that instant and `absent` after it. The driver may also be a **module**, a named group of families such as a pathway or an operon, declared with `modules=` and read as `g.completion("flagellum")` — at every resolution, since a module is a set of declared families. Completion is a fraction rather than a state, so it takes a curve, and that is where a threshold belongs: `lambda f: 8.0 if f > 0.8 else 1.0` reads "eight times faster once four fifths of it is there". It is a fraction for a measured reason. Under independent loss, the chance that *every* family survives falls off geometrically with the module's size — on a 200-tip tree a module of three was complete at 189 tips, where one of six was complete at none. A complete-or-not driver would therefore be a constant for all but the smallest modules.

**A trait driving a trait** (**11**) is the plainest case of a level conditioning itself. A body size that diffuses four times faster in the deep is `rate = PerLineage(1.0).scaled_by(depth, {"deep": 4.0, "shallow": 1.0})`; a `switch` written per transition drives only the transitions you name.

**A sequence as a driver** (**6**, **7**) reads its composition, and has a section of its own below.

Two things are refused, for different reasons.

**`varying_among('families', ...)` and a driven rate cannot be set in the same run**, at the family or ordered resolution. One weights the lineages by a driver; the other weights by family. Setting both on the rates means weighting by the product, which is not what either says. A driven *extent*, or a driven `transfer_to`, is a different axis and runs alongside a per-family draw unchanged. That is a limit of the code.

**A genome cannot be driven by a sequence**, and that is a limit of the model. A sequence was grown along the gene trees the genome run produced, so a genome reading it back would condition a run on its own output. Those cells are shaded on the map, and the pair can only be joined (Chapter 10).

Either way, a driven rate an engine cannot honour raises and names the reason, rather than being silently dropped.

### At the nucleotide resolution

All three resolutions answer, and `presence` and `completion` mean the same thing at each. What differs is how a family is named and how it can be taken away. At the nucleotide resolution a family is a **declared gene**, so its name comes from the GFF that declared it (its `ID` or `Name`) — the evenly-spaced `genes=` layout lays its genes down unnamed, and there is nothing there to ask for. And what removes it is an arc of DNA rather than a whole copy, so a gene is gone once a deletion has taken it: the driver still reads `present` and `absent`, off the same recovered gene tree.

```python
# a GFF is a path, or the lines themselves — here two genes on one 3 kb replicon
annotation = ["##sequence-region c1 1 3000",
              "c1\tzombi2\tgene\t201\t800\t.\t+\t.\tID=dnaA",
              "c1\tzombi2\tgene\t1401\t2000\t.\t+\t.\tID=dnaN"]

nt = genomes.simulate_genomes_nucleotide(tree, gff=annotation, loss=1.0, loss_extent=300,
                                         modules={"operon": ["dnaA", "dnaN"]}, seed=9)

traits.simulate_discrete(tree, states=["harmless", "pathogenic"], start="harmless", seed=2,
                         switch = PerLineage(0.1).scaled_by(nt.presence("dnaA"),
                                                            {"present": 8.0, "absent": 1.0}))
```

## On the command line

Conditioning folds into the target level's own command. There is no conditioning command and no object to build: you grow the driver, then make an ordinary genome run whose `loss` happens to carry a `scaled_by` instead of being a bare number, with the rate in its written form.

```bash
# 1. a species tree
zombi2 species out/ --birth 1 --death 0.3 --n-extant 20 --seed 1

# 2. the driver: a habitat trait, writing its event log
zombi2 traits out/ --kind discrete \
    --states aquatic,terrestrial --switch 0.1 --seed 1 --write values tree events

# 3. the target: genomes whose loss reads that trait
zombi2 genomes out/ --duplication 0.2 --origination 0.5 --seed 2 \
    --loss "PerCopy(0.25).scaled_by('out/traits/trait_events.tsv', {'aquatic': 4.0, 'terrestrial': 1.0})"

# 4. one level further down: sequences whose substitution rate reads the same trait
zombi2 sequences out/ --model hky85 --length 1000 --seed 3 \
    --substitution "PerSite(0.05).scaled_by('out/traits/trait_events.tsv', {'aquatic': 0.5, 'terrestrial': 1.0})"
```

A conditioned run records what it read, so that re-running the driver afterwards cannot leave the runs beneath it quietly out of step. The record is a `conditioned_on` file naming the levels this run reads, and re-running the driver refuses unless you pass `--force`, which re-runs it and clears them.

The guard works **within one run directory**, because it is the directory that holds the record. A driver and a target written to two different directories, as the `--from` examples below do, are not linked: no marker is written there and re-running the driver will not warn you. It also works **between** levels rather than within one, so a trait driving a second trait writes the record but re-running the first does not invalidate the second.

Both halves of transfer take that same text, the rate with a base number in front of it and the recipient weight without one.

```bash
# the driver: a competence trait, into its own directory
zombi2 traits comp/ --kind discrete --from out/ \
    --states competent,normal --switch 0.3 --seed 1 --write events

driver=comp/traits/trait_events.tsv
zombi2 genomes comp_genomes/ --from out/ --initial-families 10 --seed 2 \
    --transfer    "PerCopy(0.1).scaled_by('$driver', {'competent': 3.0, 'normal': 1.0})" \
    --transfer-to "Recipients().weighted_by('$driver', {'competent': 3.0, 'normal': 1.0})"
```

A trait target reads that same file the same way: `--rate` for a continuous trait's variance-rate and `--switch` for a discrete one's switching rate both take the expression. `--liability`, the threshold model's variance-rate, takes no modifier yet — a bare number. Use `--name` to keep two traits in one run directory, since a trait driving a trait needs somewhere for both to sit.

```bash
zombi2 traits size/ --from out/ --kind continuous --start 0.0 --seed 4 \
    --rate "PerLineage(1.0).scaled_by('out/traits/trait_events.tsv', {'aquatic': 4.0, 'terrestrial': 1.0})"

zombi2 traits size/ --from out/ --kind discrete --name diet --states plant,fish --seed 5 \
    --switch "PerLineage(0.2).scaled_by('out/traits/trait_events.tsv', {'aquatic': 5.0, 'terrestrial': 1.0})"
```

Only a trait driver and a clade have a written form — a path to the file the trait wrote, and the clade's own literal — so every other cell of the map is Python only. That is the asterisk on the map.

## A sequence as the driver

The level furthest down can drive too, through the one number a finished sequence says that another level can use: its **GC content**. `result.gc()` gives the fraction of a lineage's DNA that is G or C at every instant, and it is read like any other continuous driver.

```python
from zombi2.params import Curve, PerLineage
from zombi2.sequences import simulate_sequences
from zombi2.sequences.substitution_models import hky85, lg

g = genomes.simulate_genomes_family(tree, initial_families=20, duplication=0.1, loss=0.15, seed=9)
seqs = simulate_sequences(g, model=hky85(2.0), length=300, seed=1)

traits.simulate_discrete(tree, states=["mesophile", "thermophile"], start="mesophile", seed=2,
                         switch = PerLineage(0.2).scaled_by(seqs.gc(),
                                                            Curve(lambda x: 20.0 ** (x - 0.5))))
```

**GC is pooled over the whole lineage**, across every family the run evolved, because that is the quantity the field measures and the one a lineage's mutational bias acts on. One family's GC is not offered: it is undefined wherever that family is absent, and a driver has to answer for every branch the target walks. Grow that family in a run of its own if its GC alone is what should drive the rate.

A run gives one sequence per gene-tree node, so pooling gives one GC per species node: the value at the end of its branch. Between two nodes the path is the straight line between them, cut into stretches of at most `step`, exactly as a continuous trait's is. GC is a number, so it takes a `Curve` or a `Scalar`, never a `{state: factor}` table.

### Any letters, so any amino-acid frequency

`gc()` is the named door onto a general one. `composition(letters)` counts whichever letters of the run's own alphabet you ask for, so an amino-acid frequency is the same driver with different letters — one residue, or a set of them:

```python
proteins = simulate_sequences(g, model=lg(), length=300, seed=1)

traits.simulate_discrete(tree, states=["mesophile", "thermophile"], start="mesophile", seed=2,
                         switch = PerLineage(0.2).scaled_by(proteins.composition("KR"),
                                                            Curve(lambda x: 40.0 ** (x - 0.1))))
```

`gc()` *is* `composition("GC")`, with one extra check: it refuses a protein run, where G and C are glycine and cysteine and the call is ambiguous rather than wrong. Letters outside the run's alphabet are refused too. They occur nowhere, so the driver would read 0.0 on every lineage, and the run would be the undriven model while its log recorded a driven rate.

What this does **not** reach is a pattern — a motif, a site. A composition is a count over a lineage's whole complement, while a motif is a property of one copy. A sequence also exists only at the nodes, so a motif driver would need a pooling rule and an interpolated switch point, where a gene's presence has an exact one.

**It drives forwards only.** A sequence lives inside a gene, so it was grown along the gene trees a genome run produced. Handing GC back to a genome rate asks a run to condition on its own output, and that pair — genomes and sequences — can only be joined, never conditioned; the genome level refuses it by name. What GC can drive is what comes after it: a trait on the same species tree, or a further sequence run.

## Outputs

Conditioning adds no format. A conditioned run writes the target level's own files, listed in that level's own chapter, and adds one record: `conditioned_on`, naming the levels this run reads, on a rate or on `transfer_to`, one per line. When the driver was grown in the same run directory its own files sit beside them, so the pairing that produced the pattern is kept on disk. Appendix B gives the columns and the formats.
