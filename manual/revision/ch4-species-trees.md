# Species trees

The species tree is the backbone every other level runs on, so it is where almost every ZOMBI2 workflow begins. This chapter is about making one: the birth–death process behind it, the ways its rates can vary, what you choose to observe, and what becomes of the lineages that die.

## The birth–death process

A species tree grows by two kinds of event: a lineage **speciates**, splitting in two, or it **goes extinct** and stops. You give ZOMBI2 a **speciation rate** and an **extinction rate**, and it plays the birth–death process out: every lineage alive at a given moment has the same constant chance per unit time of splitting or dying, independently of the rest. Out comes a dated, bifurcating tree of the lineages that survive to the present.

The two rates set the tempo. Their difference fixes how fast diversity builds up. Their ratio fixes how much of the history is hidden, because a lineage that goes extinct takes its part of the tree with it. With extinction set to zero nothing is ever lost, and the tree you get is the whole tree that grew: this is the classic **Yule** (pure-birth) process, which in ZOMBI2 is just birth–death with the death rate at zero, not a separate model. As extinction rises, the tree of survivors becomes a thinner and thinner trace of the one that actually grew; what became of the lineages that died is taken up later, under *Extinct lineages*.

You also say when to stop: grow the tree to a fixed **total time**, or until it reaches a fixed **number of tips**. Both work.

```python
from zombi2 import species
# a birth–death tree of 20 surviving lineages
result = species.simulate_species_tree(birth=1.0, death=0.3, n_tips=20, seed=1)
```

By default each rate is counted **per lineage**: every branch alive is an independent chance for the event to fire. To make a rate a single shared clock for the whole tree instead, wrap it: `birth = scope.Global(1.0)`. The scope wrappers live in `zombi2.scope` (`scope.Global`, `scope.PerLineage`, …), and `Global` is capitalised because `global` is a reserved word in Python.

## What the rate depends on

So far the rates have been constant, but a birth or death rate need not be. It can depend on **time**, on **how crowded the tree is**, or on a lineage's **ancestry**. You express each the same way: multiply the base rate by a **modifier** that names what it depends on.

- **Time** — the rate changes at set moments, fast early and slow later, or any schedule you give. This is the skyline, or episodic, tree. `birth = 1.0 * mod.Time({0: 1.0, 3: 0.3})` runs at full rate until time 3, then a third of it.
- **Diversity** — the rate slows as the tree fills up, so diversity levels off toward a carrying capacity instead of growing without bound: `birth = 1.0 * mod.Diversity(cap=100)`.
- **Ancestry** — each lineage inherits its parent's rate, nudged at every split, so rates wander across the tree and close relatives resemble each other: `birth = 1.0 * mod.Inherited(spread=0.2)`.

The modifiers live in `zombi2.modifiers`. Each is a dimensionless factor on the base rate, and you can stack them with `*` (a rate that changes in time *and* saturates). Birth and death are bent independently. Note the two ways of shaping a rate: you *wrap* it to set the scope (`scope.Global`), and you *multiply* it to bend it (`* mod.Diversity`).

| From the literature | What it does | Here |
|---|---|---|
| skyline / episodic birth–death | rates change at set times | `1.0 * mod.Time({…})` |
| diversity-dependent diversification | rate slows as the tree fills | `1.0 * mod.Diversity(cap=…)` |
| ClaDS | rates drift, inherited at each split | `1.0 * mod.Inherited(spread=…)` |
| mass extinction | a fraction culled at an instant | `mass_extinctions=[(t, f)]` |

A **mass extinction** belongs here too, as the extinction rate spiking at a single instant: a fraction of the living culled at a chosen time (measured forward from the crown, like every time in ZOMBI2). Because it is a pulse and not a steady rate, it is its own argument. `mass_extinctions=[(3.0, 0.75)]` kills three-quarters of the lineages alive at time 3.

Whether a rate depends only on time, or on the tree as it grows, quietly decides how the tree gets simulated — which we return to once we look at what you actually observe.

## Sampling

Two more choices decide not how the tree grows but how much of it you get to see.

By default you see every surviving species, but real datasets are incomplete. **`sampling`** keeps only a fraction of the extant tips, chosen at random, so `sampling=0.5` gives you half. This is the standard incomplete-sampling correction, and because it only thins a tree that already grew, it costs nothing.

**`fossils`** does the opposite: it recovers some of the lineages that died. Fossils are picked up along the branches of the complete tree at a rate you set, so `fossils=0.1` scatters fossil observations through the tree's history. Because it needs the extinct lineages to exist, it grows the tree forward. In v1 the fossils are a **side output**: the sampled lineages and their ages, reported alongside the trees. The fossil does not remove its lineage and does not appear in the extant tree. The full fossilised birth–death process — fossils placed as dated sampled-ancestors *inside* the tree — is deferred to a future `tools` command.

```python
# see only half the survivors
result = species.simulate_species_tree(birth=1.0, death=0.3, n_tips=20, sampling=0.5, seed=1)

# recover fossils of extinct lineages along the branches
result = species.simulate_species_tree(birth=1.0, death=0.3, total_time=6.0, fossils=0.1, seed=1)
```

## Extinct lineages

*[Draft — the complete-vs-extant behaviour and backward sampling are settled for v1; ghost lineages are deferred to a later release.]*

Every birth–death tree is really two trees. The **complete** tree contains every lineage that ever lived, including the ones that went extinct. The **extant** tree keeps only the survivors, the extant species, and it is what you get by default, because it is almost always what you want.

When the rates are simple enough, ZOMBI2 never grows the extinct lineages at all: it samples the extant tree directly from the distribution the process implies, working backward from the present. This is fast and exact, and it is why you never chose "forward" or "backward" anywhere above — the engine takes that shortcut whenever the rates allow, and grows the tree forward only when something (diversity, ancestry, a mass extinction, fossils) needs the extinct lineages to be there.

Sometimes you want them there anyway, and keeping the complete tree hands you the extinct lineages in full. A third option — a extant tree with the dead grafted back on *approximately*, without simulating each one in detail — is **ghost lineages**: extinct tips added to a extant tree after the fact. Ghosts are the natural tool for stress-testing a method that has to cope with extinction, without paying to grow the whole complete tree. They run on a different paradigm and are **set aside for v1** — planned, but not in the first release.

## The `SpeciesResult` object

*[Draft — depends on the final result API.]*

`simulate_species_tree` returns a **`SpeciesResult`**, not a bare tree — a birth–death run produces *two* trees plus the event log, and no single tree object can hold all three. Every level returns a bundle of this shape (`GenomesResult`, `SequencesResult`, `TraitsResult`), so the four levels stay symmetric.

A `SpeciesResult` carries:

- `.extant_tree` — the survivors' tree, dated and bifurcating; this is what you get by default and hand to the next level.
- `.complete_tree` — the whole tree that grew, with the extinct lineages still on it.
- `.fossils` — the sampled fossil lineages and their ages, present only when you asked for `fossils`.
- `.events` — the event log, every speciation and extinction with its time: the compact source of truth the run exists to record.

The bundle also shares the common spine of every result — `.events`, `.seed`, and `.write(dir, include=[...])` to materialise the chosen outputs to disk. Each tree carries its topology and dated branch lengths and lets you ask for its tips, its internal nodes, and which tips are extant versus extinct. Hand `.extant_tree` straight to the next level as the tree that genomes, sequences, or traits will evolve along.

## Usage from Python

The whole range is one function call:

```python
from zombi2 import species, modifiers as mod
from zombi2 import scope              # scope wrappers: Global, PerLineage, …

# constant-rate birth–death (per lineage, the default)
result = species.simulate_species_tree(birth=1.0, death=0.3, n_tips=20, seed=1)

# Yule (pure birth) — death defaults to 0
result = species.simulate_species_tree(birth=1.0, n_tips=50, seed=1)

# skyline birth that also slows with diversity, with a global death rate
result = species.simulate_species_tree(
    birth = 1.0 * mod.Time({0: 1.0, 3: 0.5}) * mod.Diversity(cap=100),
    death = scope.Global(0.3), total_time=8.0, seed=1)

# a mass extinction and incomplete sampling
result = species.simulate_species_tree(
    birth=1.0, death=0.3, mass_extinctions=[(3.0, 0.75)],
    sampling=0.5, total_time=5.0, seed=1)
```

## Usage from the CLI

*[Draft — the CLI needs to be re-fitted to this API; the modifier syntax on the command line is still to be designed.]*

```bash
# constant-rate birth–death, 20 tips
zombi2 species --birth 1.0 --death 0.3 --tips 20 --seed 1 -o my_tree
```

## Outputs

A run writes two Newick trees by default: the **extant** tree of survivors and the **complete** tree carrying the extinct lineages too, under the `_extant` and `_complete` names, with tips labelled *extant*, *extinct*, or *unsampled* so the three are told apart. The **event log** — every speciation and extinction with its time — is always written; it is the ground truth the simulator exists to record. And if you asked for fossils, the sampled fossil lineages are written too. The full list of files lives in Appendix B.
