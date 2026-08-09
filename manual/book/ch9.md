# Conditioning

The book has so far run the levels one at a time: species tree, then genomes, then sequences, then traits. Each already depends on the tree. Sometimes you want to simulate more complex scenarios with dependencies, for example:

- A trait controls how quickly duplications occur.
- A gene controls how quickly a different gene is lost.
- The GC content of an organism controls some trait.

As long as the thing doing the controlling can be simulated first and independently — even within one level, if the two entities are separable — we talk about **conditioning**. If instead the two have to grow at the same time, we are dealing with **joining**, which is the next chapter.

## The three parts

![A conditioned run. The **driver** is a level already simulated, a habitat trait here, with the two states each lineage switches between shown below it. The **target** is what it controls, a rate here, in the run that comes next. The **connection** is the arrow: which multiplier each state of the driver hands over, so that a branch's habitat sets that branch's loss rate.](figures/conditioning_print.png){width=95%}

Chapter 2 named the three, and this chapter takes them one at a time. The **driver** is the value that is read. The **target** is what reads it. The **connection** is what joins them: it decides what number arrives, and what that number does when it gets there.

### Across levels, and within one

![What can condition what. Rows drive, columns are driven, and the eleven numbered pairs are the models. The five shaded cells are the pairs that are not: three would need two genomes for one lineage, and two are a sequence driving a genome, which would condition a run on its own output. The four boxed cells are on the diagonal — a level conditioning itself.](figures/conditioning_map_print.png){width=95%}

Not everything can act as a driver and not everything can be a target. A sequence cannot control the gene it grows inside, by construction — although it can control a different one, if the two are simulated in order. The map has every pair that is a model, and the eleven are these:

| | Driver | Target | What it says |
|---|---|---|---|
| **1** | a gene family | a gene family | a mobile element makes transfer likelier for the rest of the genome |
| **2** | a gene family | a sequence | lose the repair gene and evolve faster |
| **3** | a gene family | a trait | carry the toxin family and turn pathogenic |
| **4** | an ordered or nucleotide genome | a sequence | as **2**, with coordinates in the genome run |
| **5** | an ordered or nucleotide genome | a trait | as **3**, with coordinates in the genome run |
| **6** | a sequence | a sequence | compensatory evolution: one gene's sequence sets another's rate |
| **7** | a sequence | a trait | GC content sets how fast a trait changes |
| **8** | a trait | a gene family | habitat sets the loss rate; all four rates take a driver |
| **9** | a trait | an ordered or nucleotide genome | eleven rates, and the extents besides |
| **10** | a trait | a sequence | habitat sets the substitution rate |
| **11** | a trait | a trait | one character sets another's `rate` or `switch` |

## The driver

A driver is the input to the function that controls the target. The easiest one to think about is a trait that controls a rate at the genome level, but a trait is not the only thing that can drive. The presence of a gene family works too, and so does how complete a declared group of families is.

### What can be a driver

| Driver | What it offers | How you write it |
|---|---|---|
| a discrete trait | one of its states | the trait's result, or the path to the event log it wrote |
| a continuous trait | a number | the same, read every `step` of time — a hundredth of the tree's height by default |
| a gene family | `present` or `absent` | `g.presence("IS1")`, for families named with `family_names=` |
| a module | a fraction, 0 to 1 | `g.completion("flagellum")`, for a group of families declared with `modules=` |
| a sequence's composition | a number, 0 to 1 | `seqs.gc()`, or `seqs.composition("KR")` for any letters of the alphabet |
| a clade | one of the named groups | `Clade({"fast": ["n12", "n27"]})` — a fact about the tree, so nothing is grown first |

A driver is read wherever it changes, not once per branch: a lineage that switches habitat halfway down one loses genes at one rate before the switch and another after it. Only a trait and a clade have a written form, a path and a literal, so every other driver here is Python only.

## The target

A target is the parameter the connection is written on — what would otherwise be a plain number. In Python it is the keyword, on the command line the flag. A level takes only the targets it has, and refuses a driven parameter it cannot honour rather than accepting and ignoring it.

### What can be a target

| Target | Kind | Level |
|---|---|---|
| `duplication`, `transfer`, `loss`, `origination` | how often | genomes, every resolution |
| `inversion`, `transposition`, `translocation`, `fission`, `fusion`, `chromosome_origination`, `chromosome_loss` | how often | genomes, ordered and nucleotide |
| `substitution` | how often | sequences |
| `rate` (continuous), `switch` (discrete) | how often | traits |
| every event's extent | how much | genomes, ordered and nucleotide — Python only |
| `transfer_to` | which one | genomes, every resolution |

Only the first kind is a rate. An **extent** says how much an event takes, in genes at the ordered resolution and base pairs at the nucleotide one; driven together with its rate, the two multiply. **`transfer_to`** is the odd one: its number is not a multiplier but a weight, read on every candidate lineage and drawn from in proportion, so a weight of 0 means the lineage cannot receive at all.

## The connection

The connection is written as two things: a **verb**, which says what the number does when it arrives, and a **mapping**, which says what the number is.

```
loss = PerCopy(0.25).scaled_by(habitat, {"aquatic": 4.0, "terrestrial": 1.0})
```

`scaled_by` is the verb and the dict is the mapping. Those are the argument names too: `scaled_by(driver, mapping)`.

### Ways of connecting

Three verbs, because a number arriving at a target can do three things.

| Verb | What the number does | Written on |
|---|---|---|
| `scaled_by` | multiplies the base in front of it | a rate, an extent |
| `set_by` | replaces the base, in the rate's own units, so nothing is written in front | a rate |
| `weighted_by` | weighs the candidates against each other | `transfer_to`, from `Recipients()` |

`set_by` is for when the literature states the rate itself — *the loss rate is 1.0 in the water* — rather than a multiple of a base you had to invent. `weighted_by` needs no base because weights are normalised, so doubling them all changes nothing.

Four mappings, and which one you need follows from the **driver**, not from what it drives.

| The driver gives | Mapping | What you give it |
|---|---|---|
| named states | `Table` | one factor per state; a bare dict is read as one, and a state left out is unchanged |
| a number | `Curve` | any function of it; a bare function is read as one, and `bound=` caps it |
| a number | `Scalar` | a strength, giving `exp(strength × value)` |
| a pair of states | `Between` | one weight per ordered (donor, recipient) pair, `default=` for the rest — `transfer_to` only |

Whatever the shape, the number that comes out is non-negative and has no units.

### Writing one

A conditioned run is two ordinary runs, the driver first:

```python
from zombi2 import species, traits, genomes
from zombi2.params import PerCopy

tree = species.simulate_species_tree(birth=1.0, death=0.3, n_extant=20, seed=1)
habitat = traits.simulate_discrete(tree, states=["aquatic", "terrestrial"], switch=0.1, seed=1)

genomes.simulate_genomes_family(tree,
    loss = PerCopy(0.25).scaled_by(habitat, {"aquatic": 4.0, "terrestrial": 1.0}),
    duplication=0.2, origination=0.5, seed=2)
```

and the same thing on the command line, where the driver reaches the second run as the file the first one wrote:

```bash
zombi2 species out/ --birth 1 --death 0.3 --n-extant 20 --seed 1

zombi2 traits out/ --kind discrete \
    --states aquatic,terrestrial --switch 0.1 --seed 1 --write values tree events

zombi2 genomes out/ --duplication 0.2 --origination 0.5 --seed 2 \
    --loss "PerCopy(0.25).scaled_by('out/traits/trait_events.tsv', {'aquatic': 4.0, 'terrestrial': 1.0})"
```

There is no conditioning command and no object to build: the second run is an ordinary genome run whose `loss` happens to carry a `scaled_by` instead of a bare number. It records what it read in a `conditioned_on` file, and re-running the driver afterwards refuses unless you pass `--force`, so that the runs beneath it cannot be left quietly out of step.
