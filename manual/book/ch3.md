# Species trees

The species tree is the backbone every other level runs on, so it is where almost every ZOMBI2 workflow begins. This chapter is about making one: the birth–death process behind it, the ways its rates can vary, what you choose to observe, and what becomes of the lineages that die.

## The birth–death process

A species tree grows by two kinds of event: a lineage **speciates**, splitting in two, or it **goes extinct** and stops. You give ZOMBI2 a **speciation rate** and an **extinction rate**, and it plays the birth–death process out: every lineage alive at a given moment has the same constant chance per unit time of splitting or dying, independently of the rest. 

The two rates set the tempo. Their difference fixes how fast diversity builds up. Their ratio fixes how much of the history is hidden, because a lineage that goes extinct takes its part of the tree with it. With extinction set to zero nothing is ever lost, and the tree you get is the whole tree that grew: this is the classic **Yule** (pure-birth) process. As extinction rises, the tree of survivors becomes a thinner and thinner trace of the one that actually grew; what became of the lineages that died is taken up later, under *Extinct lineages*.

![A species tree grown by the birth–death process. Every lineage alive at a given moment has the same chance per unit time of splitting or of dying. The lineages that died are drawn dashed and stop where they died; the survivors reach the present. Both are in the complete tree, and only the solid ones are in the extant tree.](figures/species_tree.pdf){width=100%}

You also say when to stop: grow the tree to a fixed **total time** (`total_time`), or until it reaches a fixed **number of surviving lineages** (`n_extant`). Both work.

The two differ in one practical way. `n_extant` bounds the run by construction; `total_time` does not, because standing diversity grows like exp((birth − death) · t), so a rate slightly too high or a time slightly too long is the difference between a thousand lineages and ten million. A run that passes **100,000 standing lineages** therefore stops with an error rather than filling memory. Raise `max_lineages` if that is the size you want, or set it to `None` to lift the guard. It never truncates: a tree cut off at a size is no longer a sample from the process you asked for, so handing one back would be worse than not running at all.

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

![Three ways a rate can vary, one tree apiece — all three stopped at the same 25 surviving lineages, so what differs is how they got there. **A** `OnTime`: the rate drops at time 2, so an early burst gives way to a long slow tail. **B** `OnTotalDiversity`: the rate falls as the tree fills toward its cap, and splits thin out near the present. **C** `FromParent`: each lineage inherits its parent's rate, so one clade radiates late while its sister stays sparse. Solid lineages survive to the present and dashed ones died, as in the previous figure.](figures/variable_rates.pdf){width=100%}

The modifiers live in `zombi2.rates.modifiers`. Each is a dimensionless factor on the base rate, and you can stack them with `*` to get a rate that changes in time *and* saturates.

Birth and death are modified independently. Give both a `FromParent` and each lineage draws its own speciation factor and its own extinction factor at every split, so the two rates drift without correlation.

## Other models

A few models do not fit the modifier framework. ZOMBI2 provides one: a **mass extinction**, where at one instant only a fraction of the living lineages survive. `mass_extinctions=[(3.0, 0.75)]` kills three-quarters of the lineages alive at time 3. It is a pulse rather than a steady rate, so it is its own argument.

![A mass extinction as a survival pulse. The tree grows under a constant birth–death process until, at one instant, a fraction of the standing lineages die together — the cohort of dots along the vertical wall. Survivors are solid and extinct lineages dashed. The lineages-through-time curve below shares the time axis and shows the diversity crash at the pulse and the recovery after it. This tree was grown with `mass_extinctions=[(2.5, 0.75)]`.](figures/mass_extinction.pdf){width=100%}

## A summary of models

| What it does | Here | From the literature |
|---|---|---|
| rates change at set times | `1.0 * mod.OnTime({…})` | skyline / episodic birth–death |
| rate slows as the tree fills | `1.0 * mod.OnTotalDiversity(cap=…)` | diversity-dependent diversification |
| rates drift, inherited at each split | `1.0 * mod.FromParent(spread=…)` | ClaDS |
| a fraction culled at an instant | `mass_extinctions=[(t, f)]` | mass extinction |

## Sampling

Two more choices decide not how the tree grows, but how much of it you get to see.

Real datasets are incomplete, but by default you see every surviving species. **`sampling`** keeps a fraction of the extant tips, chosen at random, so `sampling=0.5` gives you half. It thins a tree that has already grown, so it costs nothing.

**`fossils`** does the opposite: it recovers lineages from the past. Fossils are picked up along **every** branch of the complete tree at a rate you set — a surviving lineage's branch as readily as an extinct one — so `fossils=0.1` scatters dated samples through the tree's history. They are a side output — the sampled lineages and their ages, reported alongside the trees. A fossil does not remove its lineage, and it does not appear in the extant tree.

![Sampling and fossils, the two ways a dataset falls short of the whole tree. A single complete tree shows every lineage's fate. Sampled species reach the present as solid lines and are the data you keep. Lineages alive today but not sampled reach the present as dashed lines ending in an open ring. Lineages that went extinct are dashed and stop where they died. Fossils are dated samples recovered along any branch of the complete tree, a surviving lineage's branch as readily as an extinct one, shown as black diamonds. The data is the solid tips together with the diamonds; the dashed lineages are never observed. This tree was grown with `sampling=0.6, fossils=0.15`.](figures/sampling_fossils.pdf){width=100%}

```python
# see only half the survivors
result = species.simulate_species_tree(birth=1.0, death=0.3, n_extant=20, sampling=0.5, seed=1)

# recover fossils of extinct lineages along the branches
result = species.simulate_species_tree(birth=1.0, death=0.3, total_time=6.0, fossils=0.1, seed=1)
```



## The `SpeciesResult` object

`simulate_species_tree` returns a **`SpeciesResult`**, which carries:

- `.extant_tree` — the survivors' tree, dated and bifurcating; this is what you get by default and hand to the next level.
- `.complete_tree` — the whole tree that grew, with the extinct lineages still on it.
- `.fossils` — the sampled fossil lineages and their ages, present when you asked for `fossils`.
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

The command mirrors the Python call. The base rates, the stop condition, and the sampling and fossil knobs each have a flag:

```bash
# a birth–death tree of 20 surviving lineages
zombi2 species out/ --birth 1.0 --death 0.3 --n-extant 20 --seed 1

# grow to time 5, with a mass extinction at time 3 and half the survivors sampled
zombi2 species out/ --birth 1.0 --death 0.4 --total-time 5 --mass-extinction 3 0.75 --sampling 0.5 --seed 2
```

## Outputs

A run writes two Newick trees by default: the **extant** tree of survivors (`species_extant.nwk`) and the **complete** tree, which also carries the extinct and unsampled lineages (`species_complete.nwk`). The survivors are tips of both trees; the dead and unsampled are tips of the complete tree only.

Both trees give the root a branch length, which many simulators leave off. A run begins with a single lineage and that lineage lives for a while before it first splits, so the root's branch is the **stem**: the time from the origin to the crown. It is ordinary simulated time — genes are gained and lost along it, traits drift along it — and a tree written without it would start at the crown and lose that history. In the complete tree the stem runs from the origin to the first speciation; in the extant tree it runs from the origin to the most recent common ancestor of the survivors, absorbing whatever branches were pruned away above it.

The **event log** (`species_events.tsv`) is always written: every speciation and extinction with its time. It is the ground truth the simulator exists to record. If you asked for fossils, the sampled fossil lineages are written too.

All three land in `out/species/`, the command grouping its files by level; `--flat` writes them straight into `out/` instead. The full list of files lives in Appendix B.
