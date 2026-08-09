# Conditioning

The book has so far run the levels one at a time: species tree, then genomes, then sequences, then traits. Each already depends on the tree. Sometimes you want more than that — you want a value that has *evolved* to steer the run that comes next.

Take olfactory genes. A habitat trait, aquatic or terrestrial, evolves down the species tree. Wherever a lineage is aquatic, it loses those genes four times faster. Two more of the same shape:

- **Endosymbionts shed their genomes.** A lifestyle trait, free-living or host-restricted, drives gene loss across the board, so the lineages that moved inside a host end up with the small genomes.
- **Competent bacteria pick genes up.** A trait for natural competence raises the rate at which new families arrive in a lineage, so the trait leaves its mark on gene *gain* rather than on loss.

What makes these **conditioning** is that the value doing the steering can be finished first. A lineage's habitat does not depend on how many genes it has, so the trait can be grown on its own and written to a file before the genome run starts. Two ordinary runs, in order. The next chapter is for the cases where that does not work.

## The three parts

![A conditioned run. The **driver** is a level already simulated, a habitat trait here, with the two states each lineage switches between shown below it. The **target** is what it controls, a rate here, in the run that comes next. The **connection** is the arrow: which multiplier each state of the driver hands over, so that a branch's habitat sets that branch's loss rate.](figures/conditioning_print.png){width=95%}

Chapter 2 named the three, and this chapter takes them one at a time. The **driver** is the value that is read. The **target** is what reads it. The **connection** is what joins them: it decides what number arrives, and what that number does when it gets there.

Driver and target are not interchangeable. A driver is a value that already varies from lineage to lineage; a target is a parameter of the run that comes next. The arrow goes one way and nothing comes back.

Here is the habitat run, whole:

```python
from zombi2 import species, traits, genomes
from zombi2.params import PerCopy, PerLineage, Recipients

tree = species.simulate_species_tree(birth=1.0, death=0.3, n_extant=20, seed=1)
# 1. grow the driver: a habitat trait down the species tree
habitat = traits.simulate_discrete(tree, states=["aquatic", "terrestrial"], switch=0.1, seed=1)

# 2. grow the target: genomes whose loss reads the habitat on each lineage
genomes.simulate_genomes_family(tree,
    loss = PerCopy(0.25).scaled_by(habitat, {"aquatic": 4.0, "terrestrial": 1.0}),
    duplication=0.2, origination=0.5, seed=2)
```

and the same run on the command line, where each level is its own command and the driver reaches the second one as the file the first one wrote:

```bash
zombi2 species out/ --birth 1 --death 0.3 --n-extant 20 --seed 1

zombi2 traits out/ --kind discrete \
    --states aquatic,terrestrial --switch 0.1 --seed 1 --write values tree events

zombi2 genomes out/ --duplication 0.2 --origination 0.5 --seed 2 \
    --loss "PerCopy(0.25).scaled_by('out/traits/trait_events.tsv', {'aquatic': 4.0, 'terrestrial': 1.0})"
```

There is no conditioning command and no object to build. It is an ordinary genome run whose `loss` happens to carry a `scaled_by` instead of a bare number.

### Across levels, and within one

![What can condition what. Rows drive, columns are driven, and the eleven numbered pairs are the models. The five shaded cells are the pairs that are not: three would need two genomes for one lineage, and two are a sequence driving a genome, which would condition a run on its own output. The four boxed cells are on the diagonal — a level conditioning itself.](figures/conditioning_map_print.png){width=95%}

The habitat run crosses two levels: a trait drives a genome. That is the common case but not the only one. A level can also condition **itself** — one trait driving a second trait, one gene family driving the rest of the genome, one gene's sequence driving another's. Those are the boxed cells on the diagonal of the map, and nothing about them is special: they are still two runs, the first finished before the second starts.

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

The shaded cells are refused, and for a reason of the model rather than of the code. **A genome cannot be driven by a sequence**: a sequence was grown along the gene trees the genome run produced, so a genome reading it back would be conditioned on its own output. That pair can only be joined (Chapter 10).

One combination is refused for a reason of the code. **`varying_among('families', ...)` and a driven rate cannot be set in the same run**, at the family or ordered resolution: one weights the lineages by a driver, the other weights by family, and setting both on the rates means weighting by the product, which is not what either says. A driven *extent*, or a driven `transfer_to`, is a different axis and runs alongside a per-family draw unchanged.

Either way, a driven rate an engine cannot honour raises and names the reason, rather than being silently dropped.

## The driver

A driver is a finished value that varies along the tree. It is grown by an ordinary run, and it reaches the next run in one of two ways: as the result object, in Python, or as the path to the file that result wrote. It is the same model either way — the file is only how a second process reads what the first one grew.

```python
habitat.write("out/", outputs=("events",))     # then pass "out/trait_events.tsv"
```

**The driver is read wherever it changes**, not once per branch. Its event log gives each branch's constant stretches, and the engine steps its Gillespie at every switch, so a lineage that changes habitat halfway down a branch loses genes at one rate before the switch and another after it. A discrete driver is followed exactly.

A continuous driver has no switches to step at, so it is approximated. The branch is cut into stretches of at most `step` time units — `scaled_by(..., step=)`, by default a hundredth of the tree's height — and read at each stretch's midpoint.

Because a driver is a file on disk, re-running it after the fact would leave the runs beneath it quietly out of step. So a conditioned run writes a `conditioned_on` record naming the levels it reads, and re-running one of those refuses unless you pass `--force`. The guard works **within one run directory**, since it is the directory that holds the record, and **between** levels rather than within one: a trait driving a second trait writes the record, but re-running the first does not invalidate the second.

### What can be a driver

**A trait**, the case every example so far has used: a discrete trait offers its state, a continuous one its value.

**A gene family**, `g.presence("IS1")`. Only families you **named** with `family_names=` can be asked for, because a family that arose during the run has an id but nothing stable to call it by. Presence is exact and changes during a branch, so a lineage that loses its last copy halfway along one is `present` before that instant and `absent` after it.

**A module**, a named group of families — a pathway, an operon — declared with `modules=` and read as `g.completion("flagellum")`. It is a **fraction**, not a state, so it takes a curve, and that is where a threshold belongs: `lambda f: 8.0 if f > 0.8 else 1.0` reads "eight times faster once four fifths of it is there". It is a fraction for a measured reason: under independent loss the chance that *every* family survives falls off geometrically with the module's size, so a complete-or-not driver would be a constant for all but the smallest. On a 200-tip tree a module of three was complete at 189 tips, one of six at none.

**A clade**, the one driver that needs nothing grown first, being a fact about the tree the run is already walking. Chapter 4 writes it: `Clade({"fast": ["n12", "n27"]})`, read like any other set of named states. Membership never changes along a branch, which makes it the cheapest driver there is.

**A sequence's composition**, `seqs.gc()` — the fraction of a lineage's DNA that is G or C — or `seqs.composition(letters)` for whichever letters of the run's own alphabet you ask for, so that an amino-acid frequency is the same driver with different letters.

```python
from zombi2.params import Curve
from zombi2.sequences import simulate_sequences
from zombi2.sequences.substitution_models import hky85

g = genomes.simulate_genomes_family(tree, initial_families=20, duplication=0.1, loss=0.15, seed=9)
seqs = simulate_sequences(g, model=hky85(2.0), length=300, seed=1)

traits.simulate_discrete(tree, states=["mesophile", "thermophile"], start="mesophile", seed=2,
                         switch = PerLineage(0.2).scaled_by(seqs.gc(),
                                                            Curve(lambda x: 20.0 ** (x - 0.5))))
```

Composition is **pooled over the whole lineage**, across every family the run evolved, because that is the quantity the field measures and the one a mutational bias acts on. One family's GC is not offered: it is undefined wherever that family is absent, and a driver has to answer for every branch the target walks. A run gives one sequence per gene-tree node, so pooling gives one value per species node — the value at the end of its branch — and between two nodes the path is the straight line between them, cut into stretches of at most `step`. There is no motif driver, for the same reason: a composition is a count over a lineage's entire complement, a motif a property of one copy.

`gc()` *is* `composition("GC")`, with one extra check: it refuses a protein run, where G and C are glycine and cysteine and the call is ambiguous rather than wrong. Letters outside the run's alphabet are refused too — they occur nowhere, so the driver would read 0.0 on every lineage and the run would be the undriven model while its log recorded a driven rate.

Only a trait and a clade have a **written form** — a path to the file the trait wrote, and the clade's own literal — so every other driver here is Python only. That is the asterisk on the map.

### A note on the resolutions

`presence` and `completion` mean the same thing at all three genome resolutions. What differs is how a family is **named**: at the nucleotide resolution a family is a declared gene, so its name comes from the GFF that declared it (its `ID` or `Name`), and the evenly-spaced `genes=` layout lays its genes down unnamed, with nothing there to ask for. What takes a gene away is an arc of DNA rather than a whole copy, but the driver still reads `present` and `absent`, off the same recovered gene tree.

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

## The target

A target is a parameter of the run that comes next — the thing you write the connection *on*. In Python it is the keyword you would otherwise give a number to; on the command line it is the flag, taking the same expression as text.

Every target belongs to a level, and a level takes only the targets it has. A run refuses a driven parameter it cannot honour rather than accepting and ignoring it.

### What can be a target

![What each level offers, as a driver and as a target. Targets come in three kinds and only the first is a rate: **how often** an event fires, **how much** it takes, and **which one** receives a transfer. The last two are narrower than the first for reasons of the model, not of the code: a gene-family event moves exactly one gene and a substitution changes exactly one site, so there is nothing to size; and only a transfer has a recipient to choose. The two genome columns share a cell wherever the resolutions agree.](figures/conditioning_parts_print.png){width=88%}

**How often**, a rate. Driving `loss` makes a lineage delete more often. Every level has these: four rates at the family resolution, eleven at the structured ones, the substitution rate for sequences, and `rate` or `switch` for traits.

**How much**, an extent. Driving `loss_extent` makes each deletion bigger, in genes at the ordered resolution and base pairs at the nucleotide one. A driven extent is Python only, because the `--*-extent` flags take a plain number. Set together with the rate, the two multiply: one says the lineage sheds more often, the other that it sheds more each time.

**Which one**, a choice. `transfer_to` decides which lineage a transfer that has already fired lands on. It is the odd target, because its number is not a multiplier but a **weight**: the engine reads it on every lineage alive at that instant and draws the recipient in proportion, so five candidates at weight 1 and five at weight 2 send two thirds of the transfers to the weight-2 group. Weights are normalised, so doubling them all changes nothing — which is why `transfer_to` is written from `Recipients()`, with no base in front of it, and why `PerCopy(1.0).weighted_by(...)` is an error. A weight of 0 means the lineage cannot receive, and when every candidate weighs 0 the transfer has nowhere to land.

This is not Chapter 5's three questions about a segmental event: *where* an event starts is drawn by the engine and takes no modifier.

On the command line, a trait target reads a driver the same way a genome does — `--rate` for a continuous trait's variance-rate, `--switch` for a discrete one's switching rate — and `--name` keeps two traits in one directory, since a trait driving a trait needs somewhere for both to sit. `--liability`, the threshold model's variance-rate, takes no modifier yet: a bare number.

```bash
zombi2 traits size/ --from out/ --kind discrete --name diet --states plant,fish --seed 5 \
    --switch "PerLineage(0.2).scaled_by('out/traits/trait_events.tsv', {'aquatic': 5.0, 'terrestrial': 1.0})"
```

A conditioned run writes the target level's own files and adds one record, `conditioned_on`, naming the levels it read. Appendix B gives the columns and the formats.

## The connection

The connection is written as two things: a **verb**, which says what the number does when it arrives, and a **mapping**, which says what the number is.

```
loss = PerCopy(0.25).scaled_by(habitat, {"aquatic": 4.0, "terrestrial": 1.0})
```

`scaled_by` is the verb and the dict is the mapping. They are the argument names too: `scaled_by(driver, mapping)`. The verb is recorded on the object, so a run's log says back what you typed, and a verb the target cannot honour is refused by name rather than quietly reinterpreted.

### Ways of connecting

Three verbs, because a number arriving at a target can do three things.

**`scaled_by` multiplies.** The mapping gives a factor and the factor multiplies the base in front of it. This is the common one.

**`set_by` replaces.** A factor is a multiple of a base you had to invent. The literature often states the rate itself — *the loss rate is 1.0 in the water* — and `set_by` says that directly:

```
loss = PerCopy().set_by(habitat, {"aquatic": 1.0, "terrestrial": 0.25})
```

There is no base in front, because the driver supplies the whole number in the rate's own units; writing one anyway raises rather than silently discarding it. The scope still applies, so a per-copy rate set to 1.0 is 1.0 per copy. A rate carries one `set_by` and any number of ordinary factors. Three targets take it: the family and ordered genome resolutions, and a continuous trait's rate.

**`weighted_by` compares.** Its target is `transfer_to`, and its number is the weight above rather than a multiplier.

![The mappings. A driver with named states takes a `Table`, `{"aquatic": 4.0, "terrestrial": 1.0}`, where any state you leave out keeps its rate unchanged. A driver whose value is a number takes a `Curve`, any function you write, with `bound=` to cap it, or a `Scalar`, the log-link `exp(strength · value)` for a driver that is already a single measured number. `Between` is the odd one, reading the driver at both ends of a transfer, and only `transfer_to` can use it.](figures/conditioning_mappings_print.png){width=92%}

Which mapping you need follows from the **driver**, not from what it drives. If the driver has named states, give a name to each: `{"aquatic": 4.0, "terrestrial": 1.0}`. If its value is a number, give a function of that number: `lambda x: 2.0 ** x`. A bare dict is read as a `Table` and a bare function as a `Curve`, so neither of the common two needs a constructor. `Scalar(strength)` is a third shape for a number, giving `exp(strength · value)`; it is the form the literature usually writes when a single measured quantity enters a rate, and it is a `Curve` with the exponential already chosen. Whatever the shape, the number that comes out is non-negative and has no units.

`Between` is the fourth, and the only one that reads the driver at **both** ends. A plain table reads it on the recipient only, which says "competent lineages take DNA up more often" but cannot say "genes move between two habitats and not within them". `Between` gives a weight per ordered (donor state, recipient state) pair:

```python
from zombi2.params import Between

habitat = traits.simulate_discrete(tree, states=["marine", "soil"], switch=0.3, seed=1)
within = Between({("marine", "marine"): 1.0, ("soil", "soil"): 1.0}, default=0.0)

genomes.simulate_genomes_family(tree, transfer=0.5, initial_families=10, seed=2,
                                transfer_to=Recipients().weighted_by(habitat, within))
```

Every transfer now stays within one habitat, because every cross-habitat pair weighs 0. `default` (1.0) is the weight of any pair you leave out, so `Between({("marine", "soil"): 3.0})` enriches one direction against a baseline instead of forbidding the rest. This is the trait-driven twin of Chapter 4's `Clades`: there the groups come from the tree, here from an evolved trait, and the steering is the same. A `Between` on a rate or an extent is refused, since both are read on one lineage and have no donor to read the other end on. And `transfer_to` takes one rule, so a driven `transfer_to` cannot also be `"distance"`.

Both halves of transfer take the same text on the command line, the rate with a base number in front of it and the recipient weight without one:

```bash
# the driver: a competence trait, grown into its own directory off the same tree
zombi2 traits comp/ --kind discrete --from out/ \
    --states competent,normal --switch 0.3 --seed 1 --write events

driver=comp/traits/trait_events.tsv
zombi2 genomes comp_genomes/ --from out/ --initial-families 10 --seed 2 \
    --transfer    "PerCopy(0.1).scaled_by('$driver', {'competent': 3.0, 'normal': 1.0})" \
    --transfer-to "Recipients().weighted_by('$driver', {'competent': 3.0, 'normal': 1.0})"
```

A driver and a target in two different directories, as those two commands are, are not linked by the `conditioned_on` guard: it is one directory that holds the record.
