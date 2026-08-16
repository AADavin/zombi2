# Conditioning

The book has so far run the levels one at a time: species tree, then genomes, then sequences, then traits. Each already depends on the tree. Sometimes you want to simulate more complex scenarios with dependencies, for example:

- A trait controls how quickly duplications occur.
- A gene controls how quickly a different gene is lost.
- The GC content of an organism controls some trait.

As long as the thing doing the controlling can be simulated first and independently — even within one level, when the two are separable: one gene family driving another, one trait driving another — we talk about **conditioning**. If instead the two have to grow at the same time, we are dealing with **joining**, which is the next chapter.

## The three parts

![A conditioned run. The **driver** is a level already simulated, a habitat trait here, with the two states each lineage switches between, and how fast, shown below it. The **target** is what it controls: a parameter of the run that comes next, named with what kind it is — a rate here, counted per copy. The **connection** is the arrow, carrying the verb that joins them and the multiplier each state hands over, so that a branch's habitat sets that branch's loss rate.](figures/conditioning_print.png){width=95%}

Chapter 1 named the three, and this chapter takes them one at a time. The **driver** is the value that is read. The **target** is what reads it. The **connection** is what joins them: it decides what number arrives, and what that number does when it gets there.

### Across levels, and within one

![What can condition what. Rows drive, columns are driven. The five shaded cells are the pairs that are not: three would need two genomes for one lineage, and two are a sequence driving a genome, which would condition a run on its own output. The three boxed cells are on the diagonal — a level conditioning itself.](figures/conditioning_map_print.png){width=95%}

Not everything can act as a driver and not everything can be a target. The map's rows and columns are the same four kinds — a trait, a gene family, an ordered or nucleotide genome, a sequence — rows driving, columns driven. Five cells cannot be connected. Three would need a second genome run for the same lineage, and one genome per lineage is the model — row **5**'s mobile element drives *within* the one genome run it sits in, which is why that cell works. The other two are a sequence driving a genome, which would condition a run on its own output: a sequence cannot control the gene it grows inside, by construction, though it can control a different gene's *sequence* run if the two are simulated in order (row **10**).

| # | Driver | Target | What it says | Gallery |
|---|---|---|---|---|
| **1** | a trait | a gene family | habitat sets the loss rate; all four rates take a driver | [Co1–Co6](https://aadavin.github.io/zombi2/gallery.html#conditioning)<!--gallery:genome_reduction--><!--gallery:genome_expansion--><!--gallery:hgt_uptake--><!--gallery:continuous_conditioning--><!--gallery:curve_saturating--><!--gallery:curve_optimum--> |
| **2** | a trait | an ordered or nucleotide genome | eleven rates, and the extents besides | [Co7](https://aadavin.github.io/zombi2/gallery.html#conditioning)<!--gallery:climate_inversions--> |
| **3** | a trait | a sequence | habitat sets the substitution rate | [Co8](https://aadavin.github.io/zombi2/gallery.html#conditioning)<!--gallery:climate_substitution--> |
| **4** | a trait | a trait | one character sets another's `rate` or `switch` | [Co9–Co10](https://aadavin.github.io/zombi2/gallery.html#conditioning)<!--gallery:driven--><!--gallery:trait_drives_trait--> |
| **5** | a gene family | a gene family | a mobile element makes transfer likelier for the rest of the genome | [Co11](https://aadavin.github.io/zombi2/gallery.html#conditioning)<!--gallery:mobile_element--> |
| **6** | a gene family | a sequence | lose the repair gene and evolve faster | [Co12](https://aadavin.github.io/zombi2/gallery.html#conditioning)<!--gallery:repair_gene--> |
| **7** | a gene family | a trait | carry the toxin family and turn pathogenic | [Co13](https://aadavin.github.io/zombi2/gallery.html#conditioning)<!--gallery:gene_drives_trait--> |
| **8** | an ordered or nucleotide genome | a sequence | as **6**, with coordinates in the genome run | [Co14](https://aadavin.github.io/zombi2/gallery.html#conditioning)<!--gallery:operon_substitution--> |
| **9** | an ordered or nucleotide genome | a trait | as **7**, with coordinates in the genome run | [Co15–Co16](https://aadavin.github.io/zombi2/gallery.html#conditioning)<!--gallery:module_drives_metabolism--><!--gallery:operon_trait--> |
| **10** | a sequence | a sequence | one gene's composition indexes something about the lineage, and that sets another gene's rate | [Co17–Co18](https://aadavin.github.io/zombi2/gallery.html#conditioning)<!--gallery:gc_drives_sequence--><!--gallery:named_family_drives_sequence--> |
| **11** | a sequence | a trait | GC content sets how fast a trait changes | [Co19](https://aadavin.github.io/zombi2/gallery.html#conditioning)<!--gallery:gc_drives_trait--> |

## The driver

A driver is the input to the function that controls the target. The easiest one to think about is a trait that controls a rate at the genome level, but a trait is not the only thing that can drive. The presence of a gene family works too, and so does how complete a declared group of families is.

### What can be a driver

| Driver | What it offers | How you write it |
|---|---|---|
| a discrete trait | one of its states | the trait's result, or the path to the event log it wrote |
| a continuous trait | a number | the same, read every `step` of time — a hundredth of the tree's height by default |
| a gene family | `present` or `absent` | `g.presence("IS1")`, for families declared with `families=[family("IS1")]` |
| a module | a fraction, 0 to 1 | `g.completion("flagellum")`, for a group of families declared with `modules=` |
| a sequence's composition | a number, 0 to 1 | `seqs.gc()`, or `seqs.composition("KR")` for any letters of the alphabet |
| **one family's** composition | a number, 0 to 1 | the same, on a run restricted to it with `families=["marker"]`, plus an `absent=` |

Every row is a level **grown first** and then read. Where a lineage sits in the tree, and when it is alive, need nothing grown: `Clade` and `changing_at` (Appendix A) read facts the run already has, so they need no driver and work at every level whether or not anything is being conditioned.

A driver is read wherever it changes, not once per branch: a lineage that switches habitat halfway down one loses genes at one rate before the switch and another after it ([Co4](https://aadavin.github.io/zombi2/gallery.html#conditioning)<!--gallery:continuous_conditioning-->). That works because a discrete driver changes at moments the run can step to exactly. A continuous one never stops changing, so there are no such moments: it is sampled every `step` instead — the connection's own argument, `scaled_by(size, Curve(...), step=0.01)`, a hundredth of the tree's height unless you set it — and the rate holds between samples. On the sequence level, where a driven rate becomes a branch length, the factor is **integrated** along the branch across the driver's switches, not read once at the top.

### One family's composition

A composition is pooled over whatever the sequence run evolved, so on a whole genome it is the lineage's whole complement and belongs to no family in particular. To read **one** family's, restrict the run to it:

```python
from zombi2.genomes import family, simulate_genomes_family
from zombi2.params import Curve, PerSite
from zombi2.sequences import lg, simulate_sequences
from zombi2.species import simulate_species_tree

ct = simulate_species_tree(birth=1.0, n_extant=25, seed=1).complete_tree
g = simulate_genomes_family(ct, initial_families=8, duplication=0.05, loss=0.25,
                            origination=0.1, seed=2,
                            families=[family("marker"), family("ribosomal")])

marker = simulate_sequences(g, families=["marker"], model=lg(), length=300, seed=3)
hotness = marker.composition("IVYWREL", absent=0.40)

ribosomal = simulate_sequences(
    g, families=["ribosomal"], model=lg(), length=300, seed=4,
    substitution=PerSite(1.0).scaled_by(
        hotness, Curve(lambda s: 4.0 ** ((s - 0.32) / 0.20))))  # 4x per +0.20 of share, =1 at 0.32
```

A word on what a composition is for. It is not a mechanism: no count of residues reaches out and sets another gene's rate. What it is good for is standing in for a property of the **lineage**. The share of a protein that is I, V, Y, W, R, E or L rises with the temperature its owner lives at, which is the compositional signature of thermophily that holds across the tree of life [@zeldovich2007ivywrel]; a genome's GC content indexes its mutational regime much the same way. Written like that the driver says something a reader can check. Written as a mechanism it does not.

`families=` names families the genome run declared, and both sequence runs read the **same** genome run, so they share one genome history — the same lineages carrying the same families — and the driver is read against branches the target's run knows. Two separate genome runs would have put the driver on a history the target never met.

What that costs is the branches where the family is absent. There is no sequence there and so nothing to count, and a driver has to answer for every branch the target walks. `absent=` is that answer, and it is required rather than guessed: without it the run raises, because carrying the parent's value forward would drive those branches as though the family were still there, which is a different model ([Co18](https://aadavin.github.io/zombi2/gallery.html#conditioning)<!--gallery:named_family_drives_sequence-->).

This one is Python only, and so is every genome and sequence driver — `presence`, `completion`, a composition — because each is an object in memory. On the command line the only driver is a **file**: the trait event log a rate expression names, as in Writing one below. Reaching further would mean writing a per-node table for each of these and teaching the driver loader to read it.

## The target

A target is the parameter the connection is written on — what would otherwise be a plain setting: a number, or a recipient rule. In Python it is the keyword, on the command line the flag. A level takes only the targets it has, and refuses a driven parameter it cannot honour rather than accepting and ignoring it.

### What can be a target

| Target | Kind | Level |
|---|---|---|
| `duplication`, `transfer`, `loss`, `origination` | how often | genomes, every resolution |
| `inversion`, `transposition`, `translocation`, `fission`, `fusion`, `chromosome_origination`, `chromosome_loss` | how often | genomes, ordered and nucleotide |
| `substitution` | how often | sequences |
| `rate` (continuous), `switch` (discrete) | how often | traits |
| every event's extent | how much | genomes, ordered and nucleotide — Python only |
| `transfer_to` | which one | genomes, every resolution |

## The connection

A driver does not speak in the units of a target. A habitat is `aquatic` or `terrestrial`; a loss rate is a number of losses per copy per unit time. The connection is what gets from one to the other, so it has to say two things — what number this driver value is worth, which is the **mapping**, and what that number then does to the parameter, which is the **verb**.

```
loss = PerCopy(0.25).scaled_by(habitat, {"aquatic": 4.0, "terrestrial": 1.0})
```

That line reads as a sentence: loss is scaled by habitat, four-fold in water and unchanged on land. `scaled_by` is the verb and the dict is the mapping, and those are the argument names too: `scaled_by(driver, mapping)`.

### Ways of connecting

Three verbs, because a number arriving at a target can do three things.

| Verb | What the number does | Written on |
|---|---|---|
| `scaled_by` | multiplies the base in front of it | a rate, an extent |
| `set_by` | replaces the base, in the rate's own units, so nothing is written in front | a rate |
| `weighted_by` | weighs the candidates against each other | `transfer_to`, from `Recipients()` |

`set_by` is for when the literature states the rate itself — *the loss rate is 1.0 in the water* — rather than a multiple of a base you had to invent. `weighted_by` needs no base because weights are normalised, so doubling them all changes nothing. Written out:

```python
from zombi2 import traits
from zombi2.params import PerCopy, Recipients

habitat    = traits.simulate_discrete(ct, states=["aquatic", "terrestrial"], switch=0.1, seed=5)
competence = traits.simulate_discrete(ct, states=["competent", "closed"], switch=0.3, seed=6)

# set_by: the mapping IS the rate, so the scope is written empty
loss = PerCopy().set_by(habitat, {"aquatic": 1.0, "terrestrial": 0.25})

# weighted_by: a competence trait decides who takes DNA up
transfer_to = Recipients().weighted_by(competence, {"competent": 5.0, "closed": 1.0})
```

`Recipients()` is `transfer_to`'s empty slot, as the bare scope is a rate's: it exists so the verb has something to chain onto. The weights are per candidate lineage, normalised over whoever is alive when a transfer fires — Chapter 3's rule — so they choose the recipient and never change how many transfers happen.

Four mappings, and which one you need follows from the **driver**, not from what it drives.

| The driver gives | Mapping | What you give it |
|---|---|---|
| named states | `Table` | one factor per state; a bare dict is read as one, and a state left out is unchanged |
| a number | `Curve` | any function of it; a bare function is read as one, and `bound=` puts a ceiling on the factor, which an unbounded function needs |
| a number | `Scalar` | a strength, giving `exp(strength × value)` |
| a pair of states | `Between` | one weight per ordered (donor, recipient) pair, `default=` for the rest — `transfer_to` only |

Whatever the shape, the number that comes out is non-negative and has no units — it multiplies, or weighs, what it lands on. `set_by` is the one exception: its mapping gives the rate itself, in the rate's own units.

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

There is no conditioning command and no object to build: the second run is an ordinary genome run whose `loss` happens to carry a `scaled_by` instead of a bare number. It records what it read in a `conditioned_on` file (Appendix B), and re-running the driver afterwards refuses — the `--force` goes on that re-run — so the runs beneath it cannot be left quietly out of step.
