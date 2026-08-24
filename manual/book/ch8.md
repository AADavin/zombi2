# Dependent runs

The book has so far run the levels one at a time: a species tree first, then genomes on it, then sequences, then traits. This order exists because each level needs the previous one: a genome evolves along the species tree, and a sequence evolves inside a gene of the genome. In ZOMBI2 this order of simulation is called the hierarchy.

![The four levels and the order they are simulated in. Each level needs the one before it.](figures/fig-2-1-four-levels_print.png){width=60%}

Sometimes you want to create more complex evolutionary scenarios, for example:

- A trait controls how quickly genes are duplicated.
- A gene controls how quickly a different gene is lost.
- The GC content of an organism controls one of its traits.
- A trait controls how quickly lineages speciate, and the trait evolves along the tree that it shapes.
- The composition of a gene controls how quickly its genome loses genes, and the genome controls which copies of that gene exist.

In all these scenarios, one simulation depends on another. This chapter first distinguishes the three kinds of run, then explains the connection, which is how every dependency is written, and finally presents the joint models.

## Independent, conditioned and joint runs

Two levels can be simulated in three ways.

![The three kinds of run. Left: two independent runs, one after the other. Middle: a conditioned run, where the second run depends on the results of the first. Right: a joint run, where both levels are simulated at the same time.](figures/execution_print.png){width=95%}

In an **independent** run the levels are simulated one after the other. The second simulation runs along the output of the first, a genome along the species tree for example, but its parameters do not depend on what happened in the first run. Every run in the previous chapters was of this kind.

In a **conditioned** run the levels are also simulated one after the other, but the parameters of the second run depend on the results of the first. A habitat trait is simulated first, and the loss rate of the genome run then depends on the habitat of each lineage. In Python, the result of the first run is passed to the second as an object. On the command line, the second run reads the file that the first run wrote. The second run records what it depended on in a `conditioned_on` file (Appendix B), and re-running the driver afterwards requires `--force`, so the runs that depend on it cannot be left out of date silently.

In a **joint** run the two levels are simulated at the same time, in a single run. This is necessary when neither level can be simulated first, and there are two ways it happens. The first is a cycle of two connections: a trait controls the substitution rate of a gene, and the composition of that gene controls how the trait changes, so each simulation would need the other one to be finished. The second is a single connection that points against the hierarchy: the composition of a gene controls the loss rate of its genome, but the genome must be simulated to know which copies of the gene exist. In both cases there is no possible order, so one run simulates both levels.

A dependency is always created by a connection. A joint run without any connection is refused: without a dependency the result would be identical to two independent runs, and the error message says to run them independently instead.

## Connections

Every dependency in this chapter is written in the same way, whatever the levels involved and whether the run is conditioned or joint.

![A connection. The driver is the value the connection takes as input. The target is the parameter that depends on it. The link states how each value of the driver becomes a factor on the target.](figures/connection_abstract_print.png){width=80%}

A connection has three parts. The **driver** is the value the connection takes as input: the state of a trait, the presence of a gene family, the GC content of a gene. The **target** is the parameter that depends on it: a loss rate, a substitution rate, a speciation rate. The **link** is what joins them: it states what number each value of the driver produces, and what that number does to the target.

![A connection in full. The driver is a habitat trait with two states; the diagram under the box shows the states and how fast lineages switch between them. The target is the loss rate of the genome level. The link is the verb `scaled_by` with a table: an aquatic lineage loses genes four times faster.](figures/conditioning_print.png){width=95%}

### The driver

A driver can be the state of a discrete trait, the value of a continuous trait, the presence of a gene family, the completion of a module, or the composition of a sequence. Appendix C lists them all, with what each one offers and a gallery example.

A driver is evaluated wherever it changes, not once per branch. A lineage that switches habitat halfway along a branch loses genes at one rate before the switch and at another after it ([Co4](https://aadavin.github.io/zombi2/gallery.html#conditioning)<!--gallery:continuous_conditioning-->). A discrete driver changes at specific moments, and the run steps to them exactly. A continuous driver changes all the time, so it is sampled every `step` of time instead, an argument written on the link: `scaled_by(size, Curve(...), step=0.01)`. On the sequence level, where a driven rate becomes a branch length, the factor is integrated along the branch across the driver's changes.

In a conditioned run, a trait driver works both from Python and from the command line, where it is read from the trait's event file. The genome and sequence drivers (a presence, a completion, a composition) are objects in memory and work from Python only. A composition driver must also declare `absent=`, the value used on lineages that carry no copy of the family. Without it the run raises an error, because using the parent's last value would treat the family as still present, and that is a different model ([Co20](https://aadavin.github.io/zombi2/gallery.html#conditioning)<!--gallery:named_family_drives_sequence-->).

### The target

A target is a parameter that would otherwise be a plain number: a rate, an extent, the recipient rule of transfers, or the substitution model of a sequence. In Python it is a keyword argument, on the command line a flag. Each level accepts only its own targets, and refuses a driven parameter it does not have. The targets are the rates of every level: the four rates of a family genome, the eleven rates of an ordered genome, the thirteen of a nucleotide one, and their extents, the substitution rate of a sequence, the rate or switch of a trait, and, in joint runs only, the birth and death rates of the species tree. Two targets are not rates: `transfer_to`, the choice of a transfer recipient, and the substitution model of a sequence run. The full table is in Appendix C.

### The link

A driver and a target have different units. A habitat is `aquatic` or `terrestrial`; a loss rate is a number of losses per copy per unit of time. The link converts one into the other, and it has two parts: the **mapping** states what number each value of the driver is worth, and the **verb** states what that number does to the target.

```
loss = PerCopy(0.25).scaled_by(habitat, {"aquatic": 4.0, "terrestrial": 1.0})
```

This line defines the loss rate of a genome run: a base rate of 0.25 per copy, multiplied by 4 in aquatic lineages. `scaled_by` is the verb and the dictionary is the mapping, and these are also the names of the arguments: `scaled_by(driver, mapping)`.

There are three verbs. `scaled_by` multiplies the base rate and also applies to extents. `set_by` replaces the base rate, and its mapping then gives the rate in its own units. `weighted_by` applies to `transfer_to` and weighs the candidate recipient lineages against each other.

There are four mappings, and the right one depends on the driver. A `Table` gives one factor per state of a discrete driver. A `Curve` is a function of a numeric driver. A `Scalar` gives the factor `exp(strength × value)`. A `Between` gives one weight per pair of donor and recipient states, for transfers only. The factor that comes out of a mapping is a number without units and cannot be negative. Appendix C gives both tables, with their rules and the gallery examples.

### What can connect to what

Not every pair can be connected. A genome cannot be driven by a second genome run of the same lineage, because one genome per lineage is the model. A sequence can drive the genome that carries it, but only in a joint run: the genome must be simulated to know which copies of the sequence exist. The species tree can only be a target of joint runs, because no level can be simulated before the tree it evolves along. The full catalog is in Appendix C: a map of the pairs, and the table of every connection with the gallery examples that run it conditioned or joint.

<!-- --8<-- [start:joint] -->

## Joint runs

![What can be joined. Every arrow has two heads, because neither level can be simulated first. A loop marks a pair whose two parts belong to the same level. The pair with no arrow, species and sequences, cannot be joined.](figures/joining_map_print.png){width=88%}

Levels can also be simulated jointly, as the figure shows. Two pairs reach the species tree: a trait or the gene content drives speciation, and the tree is then a result of the run. Three pairs join two different levels along a tree that is handed to the run: a trait and a genome, a trait and a sequence, a genome and a sequence. In the remaining three cases the two joined parts belong to the same level: a gene family changes the rates of the genome that carries it, one trait changes how another trait evolves, and each of two genes changes how fast the other one evolves. A joint run holds exactly two parts, and asking for a third is refused. The species tree and a sequence cannot be joined: a sequence evolves along a gene tree, which evolves along the species tree, so joining the two would require simulating the genome as well.

### Naming the driver

The connections of a joint run are written with the same verbs and mappings as everywhere else. The difference is the driver. It cannot be passed as a finished result, because it does not exist yet, so it is named instead: `"trait"`, `"genomes:toxin"`, `"sequences:rpoB"`. The name refers to another level of the same run, and the run grows both together. A run with a single trait can call it plainly `"trait"`; names exist so that a run can hold two. Appendix C lists the named drivers and what each one offers.

A family whose presence drives must be declared in the genome spec, with `families=[family("toxin")]`. A gene whose composition drives a rate must declare it with `offers=composition("KR", absent=0.02)`: which letters are counted, and the value used on lineages that carry no copy. A gene that drives a rate usually also declares `start=`, a second substitution model whose stationary frequencies are used only to found the gene. The gene starts with that composition and then evolves under its own model, so its composition moves toward its own equilibrium. Without `start=` a gene is founded at its equilibrium, its composition barely changes, and the rate that depends on it stays constant ([Jo2](https://aadavin.github.io/zombi2/gallery.html#joining)<!--gallery:trait_and_sequence-->).

A joint run is written from the specs of the levels it simulates, for example `joint.simulate(species.birth_death(...), traits.discrete(...))`. If the species spec is present, the tree is simulated, and the stop condition (`n_extant=`, `total_time=`) is written on that spec; if it is not, the tree is passed with `tree=` and the run ends where the tree ends. The result is a `JointResult` that carries both simulated levels, the same result objects the ordinary functions return, and both levels write their usual files. When the two joined parts belong to the same level, the run stays on that level's own function, with `joint=True`: `simulate_traits`, `simulate_genomes_family` or `simulate_sequences`. The flag is checked both ways: asking for it without a live driver is an error, and depending on a live driver without it is an error. One more restriction: a transfer needs the set of lineages alive at an instant, which is unknown while the tree is still being simulated, so a genome run that grows the tree refuses transfer. The gallery presents every joint model in full, with its code ([Jo1](https://aadavin.github.io/zombi2/gallery.html#joining)<!--gallery:genome_and_sequence--> to [Jo12](https://aadavin.github.io/zombi2/gallery.html#joining)<!--gallery:quasse-->).

### Exact or sliced

A joint run advances by racing every possible event against every other, with the same Gillespie algorithm as everywhere in ZOMBI2. A discrete driver changes only at events, so between two events every rate is constant and the race is exact. When a speciation fires, it uses the trait state of the lineage it happens on; when a trait switches, the tree is untouched; both daughters of a split start from the parent's state.

Two kinds of driver change between events: a continuous trait diffuses at every instant, and a composition moves with every substitution. There is no interval over which the driven rate holds still, so there is nothing exact to race. These models slice time instead. Time is cut into steps of `step`, and the driver is held fixed within each slice, so the rates are constant inside it. At the slice boundary the driver moves: a diffusing trait by its exact transition law, a composition by the substitutions that accumulated. The approximation is only in the timing: the target uses the value the driver had at the start of the slice.

| Model | Tree | Advances by |
|---|---|---|
| a trait drives speciation | simulated | events, exact |
| gene content drives speciation | simulated | events, exact |
| a diffusing trait drives speciation | simulated | slices of `step` |
| a trait and a genome | handed | events, exact |
| two traits, each driving the other | handed | events, exact |
| a trait and a gene's sequence | handed | slices; the composition is taken at the start of the slice |
| a genome and a gene's sequence | handed | slices; the composition is held within a slice |
| a gene family drives its own genome | handed | events, exact |
| two genes, each depending on the other | handed | slices; both compositions held within a slice |

`step` is written on the link, because the right size depends on how the value is taken: a steep curve needs a finer step than a flat one. In a joint run it has no default, and leaving it out is an error, because the timescale belongs to a model that does not exist yet. The practical check is to halve `step`, rerun, and confirm that the numbers you report move by less than their seed-to-seed spread ([Jo12](https://aadavin.github.io/zombi2/gallery.html#joining)<!--gallery:quasse-->).

### On the command line

`zombi2 joint` covers the two models that drive speciation from a discrete driver. The rate is written exactly as in Python, and the remaining flags are the ones `zombi2 traits` or `zombi2 genomes` would take. Giving flags from both levels at once is an error.

```bash
zombi2 joint out/ --death 0.2 --states small,large --switch 0.3 \
    --n-extant 100 --seed 1 \
    --birth "PerLineage(1.0).scaled_by('trait', {'small': 1.0, 'large': 3.0})"

zombi2 joint out/ --origination 0.2 --loss 0.1 --family-names toxin \
    --n-extant 60 --seed 1 \
    --birth "PerLineage(1.0).scaled_by('genomes:toxin', {'present': 3.0, 'absent': 1.0})"
```

The other models are Python only, deliberately: each needs something no flag can carry, such as a curve, a rate for a single named family, a list of specs, or a start model built as a matrix.

### Literature

The state-dependent diversification models are usually known by their acronyms. This table gives the correspondence.

| What it does | From the literature | Gallery |
|-------------------|--------------------------------|---|
| a binary trait drives speciation (and extinction) | BiSSE [@maddison2007bisse] | [Jo8](https://aadavin.github.io/zombi2/gallery.html#joining)<!--gallery:bisse--> |
| a multi-state trait drives speciation | MuSSE [@fitzjohn2012diversitree] | [Jo10](https://aadavin.github.io/zombi2/gallery.html#joining)<!--gallery:musse--> |
| a trait drives speciation **and** jumps at the split | ClaSSE [@goldberg2012classe] | [Jo11](https://aadavin.github.io/zombi2/gallery.html#joining)<!--gallery:classe--> |
| a **continuous** trait drives speciation | QuaSSE [@fitzjohn2010quasse] | [Jo12](https://aadavin.github.io/zombi2/gallery.html#joining)<!--gallery:quasse--> |

<!-- --8<-- [end:joint] -->
