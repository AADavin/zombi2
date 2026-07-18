# Species trees

The species tree is the backbone every other level runs on, so it is where almost every ZOMBI2 workflow begins. This chapter is about making one: the birth–death process behind it, the ways its rates can vary, what you choose to observe, and what becomes of the lineages that die.

## The birth–death process

A species tree grows by two kinds of event: a lineage **speciates**, splitting in two, or it **goes extinct** and stops. You give ZOMBI2 a **speciation rate** and an **extinction rate**, and it plays the birth–death process out: every lineage alive at a given moment has the same constant chance per unit time of splitting or dying, independently of the rest. Out comes a dated, bifurcating tree of the lineages that survive to the present.

The two rates set the tempo. Their difference fixes how fast diversity builds up. Their ratio fixes how much of the history is hidden, because a lineage that goes extinct takes its part of the tree with it. With extinction set to zero nothing is ever lost, and the tree you get is the whole tree that grew: this is the classic **Yule** (pure-birth) process, which in ZOMBI2 is just birth–death with the death rate at zero, not a separate model. As extinction rises, the tree of survivors becomes a thinner and thinner trace of the one that actually grew; what became of the lineages that died is taken up later, under *Extinct lineages*.

You also say when to stop: grow the tree to a fixed **age**, or until it reaches a fixed **number of tips**. Both work.

```python
import zombi2 as z
# a birth–death tree of 20 surviving lineages, crown age 5
tree = z.simulate_species_tree(birth=1.0, death=0.3, n_tips=20, seed=1)
```

By default each rate is **per lineage**: every branch alive is an independent chance for the event to fire. You can instead make a rate **global**, one budget shared by the whole tree, with `per="global"`.

## What the rate depends on

So far the rates have been constant, but a birth or death rate need not be. It can depend on **time**, on **how crowded the tree is**, or on a lineage's **ancestry**. You express each the same way: multiply the base rate by a **modifier** that names what it depends on.

- **Time** — the rate changes at set moments, fast early and slow later, or any schedule you give. This is the skyline, or episodic, tree. `birth = 1.0 * Time({0: 1.0, 3: 0.3})` runs at full rate until age 3, then a third of it.
- **Diversity** — the rate slows as the tree fills up, so diversity levels off toward a carrying capacity instead of growing without bound: `birth = 1.0 * Diversity(cap=100)`.
- **Ancestry** — each lineage inherits its parent's rate, nudged at every split, so rates wander across the tree and close relatives resemble each other: `birth = 1.0 * Inherited(spread=0.2)`.

The modifiers live in `zombi2.modifiers`. Each is a dimensionless factor on the base rate, and you can stack them with `*` (a rate that changes in time *and* saturates). Birth and death are bent independently.

| From the literature | What it does | Here |
|---|---|---|
| skyline / episodic birth–death | rates change at set times | `1.0 * Time({…})` |
| diversity-dependent diversification | rate slows as the tree fills | `1.0 * Diversity(cap=…)` |
| ClaDS | rates drift, inherited at each split | `1.0 * Inherited(spread=…)` |
| mass extinction | a fraction culled at an instant | `mass_extinctions=[(t, f)]` |

A **mass extinction** belongs here too, as the extinction rate spiking at a single instant: a fraction of the living culled at a chosen time. Because it is a pulse and not a steady rate, it is its own argument. `mass_extinctions=[(3.0, 0.75)]` kills three-quarters of the lineages alive at age 3.

Whether a rate depends only on time, or on the tree as it grows, quietly decides how the tree gets simulated — which we return to once we look at what you actually observe.

## Sampling

Two more choices decide not how the tree grows but how much of it you get to see.

By default you see every surviving species, but real datasets are incomplete. **`sampling`** keeps only a fraction of the extant tips, chosen at random, so `sampling=0.5` gives you half. This is the standard incomplete-sampling correction, and because it only thins a tree that already grew, it costs nothing.

**`fossils`** does the opposite: it recovers some of the lineages that died. Fossils are picked up along the branches of the complete tree at a rate you set, so `fossils=0.1` scatters fossil observations through the tree's history. This is the fossilised birth–death process, and unlike extant sampling it needs the extinct lineages to exist, so it grows the tree forward.

```python
# see only half the survivors
tree = z.simulate_species_tree(birth=1.0, death=0.3, n_tips=20, sampling=0.5, seed=1)

# recover fossils of extinct lineages along the branches
tree = z.simulate_species_tree(birth=1.0, death=0.3, age=6.0, fossils=0.1, seed=1)
```

## Extinct lineages

*[Draft — the concepts here are settled, but the extinct-lineages API is still to be designed with you; the syntax below is a placeholder.]*

Every birth–death tree is really two trees. The **complete** tree contains every lineage that ever lived, including the ones that went extinct. The **reconstructed** tree keeps only the survivors, the extant species, and it is what you get by default, because it is almost always what you want.

When the rates are simple enough, ZOMBI2 never grows the extinct lineages at all: it samples the reconstructed tree directly from the distribution the process implies, working backward from the present. This is fast and exact, and it is why you never chose "forward" or "backward" anywhere above — the engine takes that shortcut whenever the rates allow, and grows the tree forward only when something (diversity, ancestry, a mass extinction, fossils) needs the extinct lineages to be there.

Sometimes you want them there anyway. Keeping the complete tree hands you the extinct lineages in full. And sometimes you want a reconstructed tree with the dead grafted back on *approximately*, without simulating each one in detail: **ghost lineages** add extinct tips to a reconstructed tree after the fact. Ghosts are the natural tool for stress-testing a method that has to cope with extinction, without paying to grow the whole complete tree.

## The `Tree` object

*[Draft — depends on the final `Tree` API.]*

`simulate_species_tree` returns a `Tree`. It carries the topology and the dated branch lengths, and lets you ask for its tips, its internal nodes, and which tips are extant versus extinct. You can write it to Newick, walk it, or hand it straight to the next level as the tree that genomes, sequences, or traits will evolve along.

## Usage from Python

The whole range is one function call:

```python
import zombi2 as z
from zombi2 import modifiers as mod

# constant-rate birth–death
tree = z.simulate_species_tree(birth=1.0, death=0.3, n_tips=20, seed=1)

# Yule (pure birth)
tree = z.simulate_species_tree(birth=1.0, n_tips=50, seed=1)

# skyline birth that also slows with diversity, global rates
tree = z.simulate_species_tree(
    birth = 1.0 * mod.Time({0: 1.0, 3: 0.5}) * mod.Diversity(cap=100),
    death = 0.3, per="global", age=8.0, seed=1)

# a mass extinction and incomplete sampling
tree = z.simulate_species_tree(
    birth=1.0, death=0.3, mass_extinctions=[(3.0, 0.75)],
    sampling=0.5, age=5.0, seed=1)
```

## Usage from the CLI

*[Draft — the CLI needs to be re-fitted to this API; the modifier syntax on the command line is still to be designed.]*

```bash
# constant-rate birth–death, 20 tips
zombi2 species --birth 1.0 --death 0.3 --tips 20 --seed 1 -o my_tree
```

## Outputs

*[Draft — to finalise once the outputs are settled.]*

A run writes the tree in Newick, alongside the record of events (the speciations and extinctions) behind it. When the tree distinguishes survivors from the dead, both the complete and the reconstructed tree are written, so downstream levels can run on whichever you need. The full list of files lives in Appendix B.
