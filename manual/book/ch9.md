# Joining

Conditioning works because the driver can be finished first. Sometimes it cannot:

- Body size decides how fast a lineage splits, and body size evolves along the very tree it is shaping.
- Carrying a key gene decides which lineages diversify, and gene content evolves along that same tree.

Growing the trait first would need a tree that does not exist yet; growing the tree first would need the trait. Neither can go first, so **one run simulates both**. That is what joining is. The test is one question: *can the driver be grown on its own and handed over?* If it can, condition. If it cannot, **join** ([Jo7, Jo10](https://aadavin.github.io/zombi2/gallery.html#joining)<!--gallery:key_innovation--><!--gallery:musse-->).

In both examples above the species tree is one of the two things being simulated, so it comes out of the run rather than going into it. That is a property of these two models rather than of joining: a trait and a gene family can drive each other on a tree you hand the run, and the tree is an input there like anywhere else.

![A joint run, drawn as chapter 8 draws a conditioned one — because it is the same relationship, and the only difference is the arrow coming back. Body size scales the speciation rate; the tree that rate builds is the tree body size evolves along. Neither box is only a driver or only a target, so neither is labelled one: each is the other's, which is precisely why no order exists to simulate them in.](figures/joining_print.png){width=95%}

## The three parts, without the order

A joint run still has a driver, a target and a connection, written exactly as in Chapter 8. What changes is that the driver cannot be handed over, because it does not exist yet. So it is **named** instead, and the run grows it alongside the level that reads it.

Where the tree is one of the two, the driver cannot reach the run as a file either — Chapter 8's command-line handoff would need the driver written onto a tree that does not exist yet. Both levels are simulated by one Gillespie, racing every kind of event against every other. A speciation reads the trait state on the lineage it happens on, which is what makes its rate depend on the trait; an extinction does the same; a trait change moves one lineage to another state and leaves the tree alone. When a speciation occurs, both daughters start from the parent's state. No approximation is needed: these drivers change only *at* events, so between two events every rate holds steady, which is what a Gillespie step already assumes. A driver that is a *diffusing number* is the exception, and the first of four models in this chapter that slice; it comes later in its own section.

### What can be joined

![What can be joined. Every arrow has two heads, because neither level can be finished before the other starts: a joint pair still has a driver and a target, but no order to simulate them in. Compare the figure at the head of this chapter, and the conditioning map in Chapter 8, where the arrows run one way. Solid arrows are the eight joint models ZOMBI2 has, which is every pair that can be joined at all; a loop is a level joined to itself. One pair has no arrow at all.](figures/joining_map_print.png){width=88%}

Every joinable pair is built, and they split in two. In the two that reach **Species**, the tree is one of the levels being simulated, so it comes out of the run. In the other six the tree is handed to the run like any other input: **Genomes with Traits**, where gene content and a character set each other's rates; **Traits with Sequences**, where a character sets a gene's substitution rate and that gene's composition sets how fast the character switches; **Genomes with Sequences**, where the genome decides which sequences exist and their composition decides how fast it changes; **Genomes with itself**, where one named family sets the rates of the genome it sits in; **Traits with itself**, two characters each reading the other; and **Sequences with itself**, two genes each reading how much of the other is a given set of letters. Joining is **pairwise**: one run holds two levels, and asking for a third is refused.

One pair has no arrow, and it is not waiting to be built. The species tree and a sequence are not a joinable pair, because a sequence lives on a gene tree, which lives on the species tree, so simulating those two together would mean simulating the genome as well.

A trait and a sequence are a pair, and the reason is that both directions exist. A trait drives a gene's substitution rate, and that gene's composition drives how fast a trait changes ([Co8, Co19](https://aadavin.github.io/zombi2/gallery.html#conditioning)<!--gallery:climate_substitution--><!--gallery:gc_drives_trait-->). Either one alone is conditioning, because the driver can be finished first. Write both at once and there is a cycle, and a cycle has no order to simulate in — which is the arrow between them.

## The driver

A joint driver is a level being simulated beside this one, named as a string rather than handed over. That name is what separates joining from conditioning: a **finished result** was produced before this run and makes it conditioned, while a **name** is something this run is still producing.

| Driver | What it reads | Mapping |
|---|---|---|
| `"traits:<name>"`, discrete | that trait's current state | a table over the states |
| `"traits:<name>"`, continuous | that trait's current value | a curve, or a `Scalar` |
| `"genomes:<family>"` | whether that family is there | a table over `present` / `absent` |
| `"genomes:count"` | how many genes the lineage has | a curve, or a `Scalar` |
| `"sequences:<name>"` | how much of that gene is a given set of letters | a curve, or a `Scalar` |

A run holding one trait may call it plainly `"trait"`; the name is what lets a run hold two.

The level itself is given as a **process spec** — `traits.discrete(...)`, `traits.continuous(...)`, `genomes.genome(...)` or `sequences.gene(...)` — a description of a process still to be grown rather than a finished result. A spec keeps all of its own options, and one of them changes the model: `at_speciation=` gives the trait a chance of jumping *at* each split, so the daughters diverge at the moment they are born while how fast a lineage splits still depends on the state it is in. For a discrete trait it is a **probability** in [0, 1], the chance the daughters shift state at the split — a different reading from the continuous trait's jump variance of Chapter 7. The trait's event log tells the two apart, `on_speciation` against `on_branch`.

A family whose presence does the driving has to be declared in the genome spec, with `families=[family("toxin")]`.

## The target

What can be driven depends on which levels are in the run. `birth` and `death` are the targets when Species is one of them, and driving those is what makes the tree an output: the driver decides which lineages split and which die, so it decides the shape of the tree it is evolving along. On a tree handed to the run the targets are the ordinary ones each level already has — the genome's four rates, a trait's `switch`.

What comes back is a `JointResult` carrying **both** simulated levels. `.species` is always there, and then `.trait`, `.genome`, `.sequences`, or two of them. They are the same result objects the standalone functions return, they share one `complete_tree` because there was only ever one tree, and the run writes both levels exactly as their own commands would.

## Writing one: the participants

A joint run is written as its **participants** — a process spec for each thing being simulated — plus, as ordinary inputs, whatever is not. Write `species.birth_death(...)` among the participants and the tree is one of the things being simulated, coming out of the run; leave it out and pass `tree=` instead, and the tree goes in like any other input.

The stop condition rides on the species spec, because stopping is a fact about growing a tree, and everything else in the run rides the tree that spec produces. On a handed tree there is no species spec and no stop condition: the run ends where the tree does.

One rule decides the function: **`joint.simulate` is needed only when more than one level is in the run.** A level reading *itself* folds into that level's own function — exactly as a conditioned run folds into the driven level's — and the examples below say which is which as they come.

### A trait drives speciation

Body size makes large lineages speciate twice as fast as small ones, and go extinct half as often:

```python
from zombi2 import joint, species, traits
from zombi2.params import PerLineage

joint.simulate(
    species.birth_death(
        birth = PerLineage(1.0).scaled_by("trait", {"small": 1.0, "large": 2.0}),
        death = PerLineage(0.2).scaled_by("trait", {"small": 2.0, "large": 1.0}),
        n_extant = 100),
    traits.discrete(states=["small", "large"], switch=0.1),
    seed = 1)
```

### Gene content drives speciation

```python
from zombi2 import genomes

joint.simulate(
    species.birth_death(
        birth = PerLineage(1.0).scaled_by("genomes:toxin", {"present": 1.8, "absent": 1.0}),
        death = 0.2, n_extant = 100),
    genomes.genome(duplication=0.2, loss=0.25, origination=0.5,
                   families=[genomes.family("toxin")]),
    seed = 1)
```

### A diffusing trait drives speciation

The two above race exactly. This one does not; nor do the three later models that read a composition, and for the same reason.

Body size is a number that diffuses, and a lineage splits at a rate that rises with it. The driver is a number rather than a state, so the mapping is a `Curve` — a function from the value to a factor — and never a `{state: factor}` table:

```python
import math
from zombi2.params import Curve

joint.simulate(
    species.birth_death(
        birth = PerLineage(0.4).scaled_by("trait", Curve(lambda x: math.exp(0.5 * x)),
                                          step = 0.05),
        death = 0.05, n_extant = 70),
    traits.continuous(start=0.0, rate=1.0),
    seed = 11)
```

A discrete driver changes only at events, and an event ends the Gillespie step, so between events the rate is a constant and the race is exact. A diffusion changes at every instant. There is no interval over which this birth rate holds still, and so nothing exact to draw the next event against.

So the run **slices**. Time is cut into steps of `step`, and inside a step every lineage's size is held where it was. The rates are then constant and the race inside the slice is the ordinary one. At the boundary each lineage's size moves by the exact transition law of its own diffusion. The trait is therefore exact, and what is approximated is its grip on speciation: a lineage splits at the rate its size had at the top of the slice rather than at that instant.

`step` is written on the connection rather than on the run, because it belongs to that reading: a steep curve needs a finer step than a flat one. In a conditioned run `step` has a default, because the driver is finished and a hundredth of its tree's height is a scale the run can see; here the timescale belongs to a model still being grown, so there is no default, and leaving it out is refused with the reason. The check is to halve `step`, rerun, and see whether the numbers you report move by more than their seed-to-seed spread ([Jo12](https://aadavin.github.io/zombi2/gallery.html#joining)<!--gallery:quasse-->).

## On a tree you hand the run

The three models above simulate the species tree. The six that follow take one, and are written the same way — participants, minus the species spec, plus the tree, handed as `tree=` or inside the genome run that carries it.

### A trait and a genome, each other's driver

A cave costs a lineage its genes; losing the eye commits it further to the cave. Neither can be simulated first, because to grow the habitat you need to know whether the eye is still there, and whether the eye survives depends on the habitat ([Jo5](https://aadavin.github.io/zombi2/gallery.html#joining)<!--gallery:cave_genomes-->).

```python
from zombi2 import genomes, joint, species, traits
from zombi2.genomes import family
from zombi2.params import PerCopy, PerLineage

tree = species.simulate_species_tree(birth=1.0, n_extant=40, seed=4).complete_tree

joint.simulate(
    genomes.genome(duplication=0.05, origination=12.0, initial_families=60,
                   loss=PerCopy(0.30).scaled_by("trait", {"cave": 5.0, "surface": 1.0}),
                   families=[family("eye")]),
    traits.discrete(states=["surface", "cave"], start="surface",
                    switch={"surface->cave": PerLineage(0.02).scaled_by(
                                "genomes:eye", {"present": 1.0, "absent": 25.0}),
                            "cave->surface": 0.10}),
    tree = tree, seed = 1)
```

**Transfer needs the set of lineages alive at an instant.** On a handed tree that set is known, so transfer works here; where the tree is being simulated it is still forming, so transfer is refused there.

The result carries both levels, `.trait` and `.genome`, and both write the files their own commands write.

### Two traits, each reading the other

One trait grown first and then read is conditioning, and Chapter 8 covers it. This is the case with no order — body size sets how readily a lineage goes underground, and living underground sets how readily it grows ([Jo4](https://aadavin.github.io/zombi2/gallery.html#joining)<!--gallery:trait_loop-->). One level, so — by the participants rule — it stays on the trait level's own function: `simulate_traits` takes a list of specs, and `joint=True` says the two read each other live.

```python
from zombi2 import traits
from zombi2.params import PerLineage

traits.simulate_traits(tree, [
    traits.discrete(name="habitat", states=["surface", "cave"], start="surface",
                    switch={"surface->cave": PerLineage(0.05).scaled_by(
                                "traits:size", {"small": 1.0, "large": 8.0}),
                            "cave->surface": 0.1}),
    traits.discrete(name="size", states=["small", "large"], start="small",
                    switch={"small->large": PerLineage(0.05).scaled_by(
                                "traits:habitat", {"surface": 1.0, "cave": 6.0}),
                            "large->small": 0.1})],
    joint = True, seed = 1)
```

This one is **exact**, and for a reason worth knowing. Two traits that read each other are one Markov chain over the pairs of their states: from *(surface, small)* the only moves are to *(cave, small)* at the habitat's rate read with size sitting at small, or to *(surface, large)* at size's rate read with habitat sitting at surface. Nothing moves both at once, because two switches never coincide. So the pair has an ordinary generator, and the same branch walk a single trait takes runs it with nothing thinned and nothing approximated.

What comes back is one result per trait, keyed by name — each exactly what `simulate_discrete` returns, so everything that reads one reads these.


### A trait and a gene's sequence

The two ends of the map, joined. A character sets how fast a gene evolves; how much of that gene is a given set of residues sets how readily the character switches:

```python
import numpy as np
from zombi2 import genomes, joint, sequences, species, traits
from zombi2.genomes import family
from zombi2.params import Curve, PerLineage, PerSite
from zombi2.sequences.substitution_models import AMINO_ACIDS, reversible

ct = species.simulate_species_tree(birth=1.0, n_extant=20, seed=1).complete_tree
g = genomes.simulate_genomes_family(ct, initial_families=3, duplication=0.0, loss=0.0, seed=2,
                                    families=[family("rpoB")])

m = sequences.lg()      # m.Q is the rate matrix, m.stationary its equilibrium frequencies
S = (m.Q / m.stationary[None, :] + (m.Q / m.stationary[None, :]).T) / 2.0
np.fill_diagonal(S, 0.0)
pi = m.stationary.copy()
pi[[AMINO_ACIDS.index(c) for c in "KR"]] *= 0.12
kr_poor = reversible(S, pi / pi.sum(), name="KR-poor LG", alphabet=AMINO_ACIDS)

r = joint.simulate(
    traits.discrete(name="habitat", states=["cold", "hot"], start="cold",
                    switch={"cold->hot": PerLineage(0.5).scaled_by(
                                "sequences:rpoB", Curve(lambda x: 0.05 + 30.0 * x), step=0.05),
                            "hot->cold": 0.3}),
    sequences.gene(name="rpoB", model=sequences.lg(), length=250, start=kr_poor,
                   offers=sequences.composition("KR", absent=0.02),
                   substitution=PerSite(0.6).scaled_by("trait", {"hot": 4.0, "cold": 1.0})),
    genomes=g, seed=5)
```

Two keywords appear here for the first time. `offers=` is what the gene **publishes** for the others to read — here its lysine-plus-arginine share, with `absent=` the value a lineage carrying no copy reads instead, the same declaration Chapter 8 asks for. `start=` is a second model whose **stationary frequencies alone** are used: the gene is founded from them and then evolves under its own model, so it arrives with a foreign composition and **ameliorates** toward its own — which is what gives the driver somewhere to move. Without it a gene founds at its own equilibrium and sits there, and a rate reading it reads a constant. The `kr_poor` block above builds such a start: LG's chemistry over lysine- and arginine-depleted frequencies.

The **genome run** is what you hand over here rather than a bare tree, and for the reason the pair is on the map at all: a sequence lives on a gene tree, which lives on the species tree, so the genome run is the one thing that carries both.

Species time is sliced. Inside a slice the trait takes its own Gillespie, so its switches mid-slice are exact; the gene's branch length is then the trait's factor **integrated** across those switches and drawn once, which is exact as well — the matrix does not change when the trait switches, only the length does. What is approximated is the composition the trait reads: it belongs to the top of the slice rather than to each instant ([Jo2](https://aadavin.github.io/zombi2/gallery.html#joining)<!--gallery:trait_and_sequence-->).

### A genome and a gene's sequence

The last cell, and the one whose two halves are most tangled. The genome decides which sequences exist — when a gene is copied, when a copy moves lineage, when a copy dies. The sequences decide how fast the genome changes, because a genome rate reads their composition.

```python
from zombi2 import genomes, joint, sequences, species
from zombi2.params import Curve, PerCopy

ct = species.simulate_species_tree(birth=1.0, n_extant=20, seed=1).complete_tree
at_rich = sequences.hky85(2.0, frequencies=(0.40, 0.10, 0.10, 0.40))

r = joint.simulate(
    genomes.genome(duplication=0.15, origination=0.05, initial_families=25,
                   families=[genomes.family("hisA")],
                   loss=PerCopy(0.15).scaled_by(
                       "sequences:hisA", Curve(lambda gc: 30.0 ** ((0.35 - gc) / 0.2)),
                       step=0.05)),
    sequences.gene(name="hisA", model=sequences.hky85(2.0), length=250, start=at_rich,
                   substitution=0.8,
                   offers=sequences.composition("GC", absent=0.35)),
    tree=ct, seed=3)
```

Both levels are participants, so both come out. The **tree** is the one thing handed over, and the gene trees are not: they come out of the genome participant. That is what separates this from every conditioned run at this level, where the gene trees go in. `start=` works as above — founded AT-rich, the gene ameliorates toward its model's own composition, and the loss rate reads the moving GC.

Species time is sliced, for the reason every sliced model has: a composition moves with every substitution, so a genome rate reading it is never constant. Inside a slice the composition is held where it was. At an event the picked copy's sequence is carried to that instant before it is cloned, or before it ends — otherwise two copies would each redraw the stretch they actually shared ([Jo1](https://aadavin.github.io/zombi2/gallery.html#joining)<!--gallery:genome_and_sequence-->).

This is the family resolution only. At the nucleotide resolution an event takes a run of base pairs and can fall inside a gene, so a copy's sequence stops being one string; carrying that through a live race is different work, and it is refused by name.

### A level joined to itself

One gene family setting the rates of the genome it sits in is one level, one engine and one result, so — by the participants rule again — it stays on `simulate_genomes_family`, and `joint=True` is how the run says what it is ([Jo6](https://aadavin.github.io/zombi2/gallery.html#joining)<!--gallery:mobile_element_joint-->).

```python
genomes.simulate_genomes_family(
    tree, initial_families=25, duplication=0.05, loss=0.12, joint=True, seed=7, max_family_size=8,
    families=[family("IS1", origin=("n11", None), transfer=PerCopy(0.30), loss=0.08)],
    transfer=PerCopy(0.025).scaled_by("genomes:IS1", {"present": 30.0, "absent": 1.0}))
```

An insertion sequence makes the genome it is in donate genes thirty times as often, and it moves itself, so it spreads into lineages that never inherited it and makes them donors too. Inside `family(...)`, rates are that family's own and **replace** the run's for it — `IS1` transfers at 0.30 and is lost at 0.08 whatever the genome-wide numbers say — and `origin=("n11", None)` plants it at the start of branch `n11`. `joint=True` is checked both ways: asking for it when nothing reads a live driver is an error, and reading one without it is an error too.

### Two genes, each reading the other's composition

The sequence level's own loop, and the last of the three self-joins. Two genes are named, each one's substitution rate reads how much of the *other* is a given set of letters, and one run walks both:

```python
import numpy as np
from zombi2.genomes import family, simulate_genomes_family
from zombi2.sequences import composition, gene, lg, simulate_sequences
from zombi2.sequences.substitution_models import AMINO_ACIDS, reversible
from zombi2.params import Curve, PerSite
from zombi2.species import simulate_species_tree

ct = simulate_species_tree(birth=1.0, n_extant=20, seed=1).complete_tree
g = simulate_genomes_family(ct, initial_families=4, duplication=0.0, loss=0.0, seed=2,
                            families=[family("hisA"), family("hisF")])

# where the two genes start from: the same KR-poor LG the trait-and-sequence model built
m = lg()
S = (m.Q / m.stationary[None, :] + (m.Q / m.stationary[None, :]).T) / 2.0
np.fill_diagonal(S, 0.0)
pi = m.stationary.copy()
pi[[AMINO_ACIDS.index(c) for c in "KR"]] *= 0.12
kr_poor = reversible(S, pi / pi.sum(), name="KR-poor LG", alphabet=AMINO_ACIDS)

r = simulate_sequences(g, joint=True, seed=3, genes=[
    gene(name=name, model=lg(), length=250, start=kr_poor,
         offers=composition("KR", absent=0.02),
         substitution=PerSite(0.5).scaled_by(f"sequences:{other}",
                                             Curve(lambda x: 0.10 + 25.0 * x), step=0.05))
    for name, other in (("hisA", "hisF"), ("hisF", "hisA"))])
```

`offers=` and `start=` are the trait-and-sequence section's, doing the same work: each gene publishes its KR share and founds away from its own equilibrium, so each rate has something moving to read.

This one slices, for the reason the diffusing trait did: a composition moves with every substitution, so there is no interval over which either rate holds still. The run holds both compositions fixed across `step` of species time and releases them at the boundary; inside a slice the transition matrix is the ordinary one. The walk is by time rather than by family — every living copy of every gene advances together — which is what a cycle requires and what an ordinary run, which finishes one family before starting the next, cannot do ([Jo3](https://aadavin.github.io/zombi2/gallery.html#joining)<!--gallery:sequence_loop-->).

A joint **sequence** run has no `species_phylogram`. That file is the clock made visible, and here each gene runs at a rate the other sets, so no one set of branch lengths belongs to the run rather than to a gene. Each gene's own phylogram is written as always, with branch lengths accumulated slice by slice.

## The nine models, side by side

Slicing enters exactly where a driver moves *between* events — a diffusion, or a composition that shifts with every substitution. Everything discrete races exactly:

| Model | Tree | Advances by |
|---|---|---|
| a trait drives speciation | simulated | events — exact |
| gene content drives speciation | simulated | events — exact |
| a diffusing trait drives speciation | simulated | slices of `step` |
| a trait and a genome | handed | events — exact |
| two traits | handed | events — exact |
| a trait and a gene's sequence | handed | slices — the composition is read at the slice top |
| a genome and a gene's sequence | handed | slices — the composition is held within a slice |
| a genome joined to itself | handed | events — exact |
| two genes' compositions | handed | slices — both compositions held within a slice |

## On the command line

`zombi2 joint` covers the chapter's first two models. The rate is written exactly as in Python, and the flags that build the second level are the ones `zombi2 traits` or `zombi2 genomes` would take. Giving flags from both is an error rather than a silent choice between them.

```bash
zombi2 joint out/ --death 0.2 --states small,large --switch 0.3 \
    --n-extant 100 --seed 1 \
    --birth "PerLineage(1.0).scaled_by('trait', {'small': 1.0, 'large': 3.0})"

zombi2 joint out/ --origination 0.2 --loss 0.1 --family-names toxin \
    --n-extant 60 --seed 1 \
    --birth "PerLineage(1.0).scaled_by('genomes:toxin', {'present': 3.0, 'absent': 1.0})"
```

**These two are the whole of the command.** Every other model in this chapter is Python only, and deliberately: each needs at least one thing no flag can carry — a curve, a rate written for a single named family, a list of process specs, a start model built as a matrix. A command line is for the runs you write often enough to want short, and these are not those yet.

## Literature

The state-dependent models arrive under a wall of acronyms, and a reader who wants "a BiSSE model" should be able to find the door.

| What it does | From the literature | Gallery |
|-------------------|--------------------------------|---|
| a binary trait drives speciation (and extinction) | BiSSE [@maddison2007bisse] | [Jo8](https://aadavin.github.io/zombi2/gallery.html#joining)<!--gallery:bisse--> |
| a multi-state trait drives speciation | MuSSE [@fitzjohn2012diversitree] | [Jo10](https://aadavin.github.io/zombi2/gallery.html#joining)<!--gallery:musse--> |
| a trait drives speciation **and** jumps at the split | ClaSSE [@goldberg2012classe] | [Jo11](https://aadavin.github.io/zombi2/gallery.html#joining)<!--gallery:classe--> |
| a **continuous** trait drives speciation | QuaSSE [@fitzjohn2010quasse] | [Jo12](https://aadavin.github.io/zombi2/gallery.html#joining)<!--gallery:quasse--> |
