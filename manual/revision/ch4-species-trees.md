# Species trees

The species tree is the backbone every other level runs on, so it is where almost every ZOMBI2 workflow begins. This chapter is about making one: the birth–death process behind it, the ways its rates can vary, what you choose to observe, and what becomes of the lineages that die.

## The birth–death process

A species tree grows by two kinds of event: a lineage **speciates**, splitting in two, or it **goes extinct** and stops. You give ZOMBI2 a **speciation rate** and an **extinction rate**, and it plays the birth–death process out: every lineage alive at a given moment has the same constant chance per unit time of splitting or dying, independently of the rest. Out comes a dated, bifurcating tree of the lineages that survive to the present.

The two rates set the tempo. Their difference fixes how fast diversity builds up. Their ratio fixes how much of the history is hidden, because a lineage that goes extinct takes its part of the tree with it. With extinction set to zero nothing is ever lost, and the tree you get is the whole tree that grew: this is the classic **Yule** (pure-birth) process. As extinction rises, the tree of survivors becomes a thinner and thinner trace of the one that actually grew; what became of the lineages that died is taken up later, under *Extinct lineages*.

You also say when to stop: grow the tree to a fixed **total time** (`total_time`), or until it reaches a fixed **number of surviving lineages** (`n_extant`). Both work.

```python
from zombi2 import species
# a birth–death tree of 20 surviving lineages
result = species.simulate_species_tree(birth=1.0, death=0.3, n_extant=20, seed=1)
```

## Simulating species trees with variable rates

So far the rates have been constant, but a birth or death rate need not be. It can depend on **time**, on **how crowded the tree is**, or on a lineage's **ancestry**. You express each the same way: multiply the base rate by a **modifier** that names what it depends on.

- **On time** — the rates change at set points in time. This is the skyline, or episodic, tree. `birth = 1.0 * mod.OnTime({0: 1.0, 3: 0.3})` runs at full rate until time 3, then at a third of it.
- **On total diversity** — the rate slows as the tree fills up, so diversity levels off at a carrying capacity instead of growing without bound: `birth = 1.0 * mod.OnTotalDiversity(cap=100)`.
- **On the parent's rate** — each lineage inherits its parent's rate, nudged at every split, so rates wander across the tree and close relatives resemble each other: `birth = 1.0 * mod.FromParent(spread=0.2)`.

The modifiers live in `zombi2.rates.modifiers`. Each is a dimensionless factor on the base rate, and you can stack them with `*` to get a rate that changes in time *and* saturates.

Birth and death are modified independently. Give both a `FromParent` and each lineage draws its own speciation factor and its own extinction factor at every split, so the two rates drift without correlation.

## Other models

A few models do not fit the modifier framework. ZOMBI2 provides one: a **mass extinction**, where at one instant only a fraction of the living lineages survive. `mass_extinctions=[(3.0, 0.75)]` kills three-quarters of the lineages alive at time 3. It is a pulse rather than a steady rate, so it is its own argument.

## A summary of models

| What it does | Here | From the literature |
|---|---|---|
| rates change at set times | `1.0 * mod.OnTime({…})` | skyline / episodic birth–death |
| rate slows as the tree fills | `1.0 * mod.OnTotalDiversity(cap=…)` | diversity-dependent diversification |
| rates drift, inherited at each split | `1.0 * mod.FromParent(spread=…)` | ClaDS |
| a fraction culled at an instant | `mass_extinctions=[(t, f)]` | mass extinction |

## Sampling

Real datasets are incomplete, but by default you see every surviving species. **`sampling`** keeps a fraction of the extant tips, chosen at random, so `sampling=0.5` gives you half. It thins a tree that has already grown, so it costs nothing.

```python
# see only half the survivors
result = species.simulate_species_tree(birth=1.0, death=0.3, n_extant=20, sampling=0.5, seed=1)
```

## Extinct lineages

The **complete** tree contains every lineage that ever lived, including the ones that went extinct. The **extant** tree keeps only the survivors, and it is what you get by default because it is almost always what you want. ZOMBI2 grows the tree forward in time and records the complete tree. The extant tree is the survivors pruned out of it, with internal nodes that lose all their descendants suppressed, so it stays dated and bifurcating.

## The `SpeciesResult` object

`simulate_species_tree` returns a **`SpeciesResult`**, which carries:

- `.extant_tree` — the survivors' tree, dated and bifurcating; this is what you get by default and hand to the next level.
- `.complete_tree` — the whole tree that grew, with the extinct lineages still on it.
- `.events` — the event log: every speciation and extinction with its time, the source of truth the run exists to record.

As at every level, the bundle also carries `.seed` and `.write(dir, outputs=[...])` to write the chosen outputs to disk. Each tree carries its topology and dated branch lengths, and lets you ask for its tips, its internal nodes, and which tips are extant. Hand `.extant_tree` to the next level as the tree that genomes, sequences, or traits will evolve along.

## Usage from Python

The whole range is one function call:

```python
from zombi2 import species
from zombi2.rates import scope, modifiers as mod   # the rate grammar: scopes + modifiers

# constant-rate birth–death (per lineage, the default)
result = species.simulate_species_tree(birth=1.0, death=0.3, n_extant=20, seed=1)

# Yule (pure birth) — death defaults to 0
result = species.simulate_species_tree(birth=1.0, n_extant=50, seed=1)

# skyline birth that also slows with diversity, with a global death rate
result = species.simulate_species_tree(
    birth = 1.0 * mod.OnTime({0: 1.0, 3: 0.5}) * mod.OnTotalDiversity(cap=100),
    death = scope.Global(0.3), total_time=8.0, seed=1)

# a mass extinction and incomplete sampling
result = species.simulate_species_tree(
    birth=1.0, death=0.3, mass_extinctions=[(3.0, 0.75)],
    sampling=0.5, total_time=5.0, seed=1)
```

## Usage from the CLI

The command mirrors the Python call. The base rates, the stop condition and the sampling fraction each have a flag:

```bash
# a birth–death tree of 20 surviving lineages
zombi2 species --birth 1.0 --death 0.3 --n-extant 20 --seed 1 -o out/

# grow to time 5, with a mass extinction at time 3 and half the survivors sampled
zombi2 species --birth 1.0 --death 0.4 --total-time 5 --mass-extinction 3 0.75 --sampling 0.5 --seed 2 -o out/
```

## Outputs

A run writes two Newick trees by default: the **extant** tree of survivors (`species_extant.nwk`) and the **complete** tree, which also carries the extinct and unsampled lineages (`species_complete.nwk`). The survivors are tips of both trees; the dead and unsampled are tips of the complete tree only.

The **event log** (`species_events.tsv`) is always written: every speciation and extinction with its time. It is the ground truth the simulator exists to record. The full list of files lives in Appendix B.
