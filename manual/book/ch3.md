# Species trees

The species tree is the backbone every other level runs on, so it is where almost every workflow begins. 

## The birth–death process

A species tree grows by two kinds of event: a lineage **speciates**, splitting in two, or it **goes extinct** and stops. You give a **speciation rate** and an **extinction rate**, and every lineage alive at a given moment has the same constant chance per unit time of splitting or dying, independently of the rest.

The two rates set the tempo. Their difference fixes how fast diversity builds up. Their ratio fixes how much of the history is hidden, because a lineage that goes extinct takes its part of the tree with it. With extinction set to zero nothing is ever lost, and the tree you get is the whole tree that grew: this is the classic **Yule** (pure-birth) process. The lineages that died are kept: they are in the complete tree, and *Outputs* below says where.

![A species tree grown by the birth–death process. Every lineage alive at a given moment has the same chance per unit time of splitting or of dying. The lineages that died are drawn dashed and stop where they died; the survivors reach the present. Both are in the complete tree, and only the solid ones are in the extant tree.](figures/species_tree.pdf){width=100%}

You also say when to stop: grow to a fixed **total time** (`total_time`), or to a fixed **number of surviving lineages** (`n_extant`).[^stopping]

`n_extant` bounds the run by construction; `total_time` does not, because standing diversity grows like exp((birth − death) · t), so a rate slightly too high or a time slightly too long is the difference between a thousand lineages and ten million. A run that passes **100,000 standing lineages** stops with an error rather than filling memory. Raise `max_lineages` if that is the size you want, or set it to `None` to lift the guard. It never truncates: a tree cut off at a size is no longer a sample from the process you asked for.

`total_time` is not conditioned on survival. A run that dies out raises an error rather than handing back a tree with no present — so if you loop over seeds and skip the failures, you are conditioning on survival yourself, and everything downstream inherits that.

[^stopping]: The two rules also give different tree shapes, which matters if you are going to estimate rates from the trees. `n_extant` stops the first moment that many lineages are alive together, then draws one more waiting time and puts the present where the next event would have fired, so the two newest tips have a real branch length instead of a zero-length one. With `death=0` that is exactly the standard way to sample a tree of a given size [@hartmann2010sampling]. With extinction it is not: the run stops the *first* time it touches the target, so the trees come out shallower than a birth–death process conditioned on that many tips — nothing measurable at `death=0`, about a tenth of the tree height at 10 tips with `death` at 0.4 of `birth`, a third to a half at 10 tips with `death` at 0.8 of `birth`, and back to nothing by 50 tips at moderate turnover. Use `total_time` if you need the conditioned distribution exactly; otherwise say in your methods which rule made the trees.

```python
from zombi2 import species
# a birth–death tree of 20 surviving lineages
result = species.simulate_species_tree(birth=1.0, death=0.3, n_extant=20, seed=1)
```

## Simulating species trees with variable rates

A birth or death rate need not be constant. It can depend on **time**, on **how crowded the tree is**, or on a lineage's **ancestry**, and you express each the same way: multiply the base rate by a **modifier** naming what it depends on.

- **On time** — the rates change at set points in time. This is the skyline, or episodic, tree. `birth = 1.0 * mod.OnTime({0: 1.0, 3: 0.3})` runs at full rate until time 3, then at a third of it. Each entry holds from its own time up to the next, and the earliest entry also applies *backwards* to the origin — so `OnTime({3: 0.3})`, with no entry at 0, runs at a third of the rate for the whole tree rather than only after time 3. Start the schedule at 0 whenever you mean "full rate until".
- **On total diversity** — the rate slows as the tree fills up, so diversity levels off at a carrying capacity instead of growing without bound: `birth = 1.0 * mod.OnTotalDiversity(cap=100)`.
- **On the parent's rate** — each lineage inherits its parent's rate, nudged at every split, so rates wander across the tree and close relatives resemble each other: `birth = 1.0 * mod.FromParent(spread=0.2)`.
- **By lineage** — each lineage draws its own rate independently, with no memory of its parent: `birth = 1.0 * mod.ByLineage(spread=0.2)`. Same spread of rates as `FromParent`, none of it inherited.

The last two are the two answers to one question — where does a lineage's rate come from? — and they are worth holding side by side, because they differ in what they do to the *shape* of the tree and not to the spread of rates. Inherited variation lets a fast clade stay fast, so it hoards the tips and the tree comes out lopsided; independent variation reshuffles at every split, so imbalance stays near what a constant rate gives. `FromParent` is the model to fit when you believe diversification is a heritable property of a clade; `ByLineage` is its null, and the honest thing to compare against. A rate carrying both is refused: there is no model in which a lineage's rate is inherited from its parent and independent of it at once.

![Three ways a rate can vary, one tree apiece — all three stopped at the same 25 surviving lineages, so what differs is how they got there. **A** `OnTime`: the rate drops at time 2, so an early burst gives way to a long slow tail. **B** `OnTotalDiversity`: the rate falls as the tree fills toward its cap, and splits thin out near the present. **C** `FromParent`: each lineage inherits its parent's rate, so one clade radiates late while its sister stays sparse. Solid lineages survive to the present and dashed ones died, as in the previous figure.](figures/variable_rates.pdf){width=100%}

The modifiers live in `zombi2.rates.modifiers`. Each is a dimensionless factor on the base rate, and you can stack them with `*` to get a rate that changes in time *and* saturates.

Birth and death are modified independently. Give both a `FromParent` and each lineage draws its own speciation factor and its own extinction factor at every split, so the two rates drift without correlation; the same holds for `ByLineage`.

Both draws are **mean-corrected**, so widening `spread` spreads the lineages out without moving the average one off the base rate you typed. And under either, the lineage that speciates or dies is drawn in proportion to its own rate rather than uniformly — which is the point: a fast lineage is likelier to be the one that splits.

## Other models

One model does not fit the modifier framework: a **mass extinction**, where at one instant only a fraction of the living lineages survive. `mass_extinctions=[(3.0, 0.75)]` kills three-quarters of those alive at time 3.

![A mass extinction as a survival pulse. The tree grows under a constant birth–death process until, at one instant, a fraction of the standing lineages die together — the cohort of dots along the vertical wall. Survivors are solid and extinct lineages dashed. The lineages-through-time curve below shares the time axis and shows the diversity crash at the pulse and the recovery after it. This tree was grown with `mass_extinctions=[(2.5, 0.75)]`.](figures/mass_extinction.pdf){width=100%}

## A summary of models

| What it does | Here | From the literature |
|---|---|---|
| rates change at set times | `1.0 * mod.OnTime({…})` | skyline / episodic birth–death [@stadler2011mammalian] |
| rate slows as the tree fills | `1.0 * mod.OnTotalDiversity(cap=…)` | diversity-dependent diversification [@rabosky2008densitydependent; @etienne2012diversitydependence] |
| rates drift, inherited at each split | `1.0 * mod.FromParent(spread=…)` | ClaDS [@maliet2019clads] |
| rates vary, drawn afresh per lineage | `1.0 * mod.ByLineage(spread=…)` | uncorrelated ("relaxed") rates |
| a fraction culled at an instant | `mass_extinctions=[(t, f)]` | mass extinction |

## Sampling

Two more choices decide not how the tree grows, but how much of it you get to see.

By default you see every surviving species. **`sampling`** keeps a random fraction of the extant tips, so `sampling=0.5` gives you half [@stadler2009incomplete]. It thins a tree that has already grown, so it costs nothing.

`n_extant` counts survivors, and sampling happens afterwards, so the two compose rather than cancel: `n_extant=20, sampling=0.5` grows to 20 survivors and then shows you about 10 of them. If you want 20 tips in hand, ask for 40. The rest are not gone — they stay in the complete tree with the fate `unsampled`, which is why the run's summary counts them separately from the extinct.

**`fossils`** does the opposite: it recovers lineages from the past [@heath2014fossilized; @gavryushkina2014sampledancestor]. Fossils are picked up along **every** branch of the complete tree at a rate you set — a surviving lineage's branch as readily as an extinct one — so `fossils=0.1` scatters dated samples through its history. They are a side output, reported alongside the trees; a fossil does not remove its lineage and does not appear in the extant tree.

![Sampling and fossils, the two ways a dataset falls short of the whole tree. A single complete tree shows every lineage's fate. Sampled species reach the present as solid lines and are the data you keep. Lineages alive today but not sampled reach the present as dashed lines ending in an open ring. Lineages that went extinct are dashed and stop where they died. Fossils are dated samples recovered along any branch of the complete tree, a surviving lineage's branch as readily as an extinct one, shown as black diamonds. The data is the solid tips together with the diamonds; the dashed lineages are never observed. This tree was grown with `sampling=0.6, fossils=0.15`.](figures/sampling_fossils.pdf){width=100%}

```python
# see only half the survivors
result = species.simulate_species_tree(
    birth=1.0, death=0.3, n_extant=20, sampling=0.5, seed=1)

# recover fossils of extinct lineages along the branches
result = species.simulate_species_tree(
    birth=1.0, death=0.3, total_time=6.0, fossils=0.1, seed=2)
```

## The `SpeciesResult` object

`simulate_species_tree` returns a **`SpeciesResult`**, which carries:

- `.complete_tree` — the whole tree that grew, with the extinct lineages still on it. **This is what the next level runs along**, which is what lets a gene be transferred out of a lineage that later dies.
- `.extant_tree` — the survivors' tree, dated and bifurcating: the tree an analysis would be handed, and a projection of the run rather than its substrate.
- `.fossils` — the sampled fossil lineages and their ages, present when you asked for `fossils`.
- `.events` — the event log: every speciation and extinction with its time, the source of truth the run exists to record.

As at every level it also carries `.seed` and `.write(dir, outputs=[...])`. Each tree carries its topology and dated branch lengths, and lets you ask for its tips, its internal nodes, and which are extant.

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
    death = scope.Global(0.3), total_time=8.0, seed=2)

# a mass extinction and incomplete sampling
result = species.simulate_species_tree(
    birth=1.0, death=0.3, mass_extinctions=[(3.0, 0.75)],
    sampling=0.5, total_time=5.0, seed=2)
```

## Usage from the CLI

The command mirrors the Python call — the base rates, the stop condition, and the sampling and fossil knobs each have a flag:

```bash
# a birth–death tree of 20 surviving lineages
zombi2 species out/ --birth 1.0 --death 0.3 --n-extant 20 --seed 1

# grow to time 5, with a mass extinction at time 3 and half the survivors sampled
zombi2 species out/ --birth 1.0 --death 0.4 --total-time 5 \
    --mass-extinction 3 0.75 --sampling 0.5 --seed 2
```

## Outputs

| File | What it holds |
|---|---|
| `species_complete.nwk` | the whole tree that grew, with the extinct and unsampled lineages still on it — **what the next level reads** |
| `species_extant.nwk` | the tree of sampled survivors: what an analysis would be handed |
| `species_events.tsv` | every speciation and extinction, with its time |
| `species_fates.tsv` | each tip's fate: `extant`, `extinct` or `unsampled` |
| `species_fossils.tsv` | the sampled fossil lineages and their ages, when you asked for `fossils` |
| `species_summary.json` | what the run produced: counts, tree height, stem, total branch length, realised rates |

Both trees give the root a branch length, the stem of Chapter 2: in the complete tree it runs from the origin to the first speciation, in the extant tree from the origin to the most recent common ancestor of the survivors, absorbing whatever branches were pruned away above it. They land in `out/species/`, or straight into `out/` with `--flat`; Appendix B lists every file with its format and its default.
