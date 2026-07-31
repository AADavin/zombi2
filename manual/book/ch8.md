# Trait evolution

The trait level evolves **phenotypes**: a body size, a habitat, the presence or absence of a structure. A trait evolves along the species tree like everything else here. There are two kinds, continuous and discrete:

```python
from zombi2 import species, traits
from zombi2.rates import modifiers as mod

tree = species.simulate_species_tree(birth=1.0, death=0.3, n_extant=20, seed=1)
result = traits.simulate_continuous(tree, start=0.0, rate=1.0, seed=1)   # a real value
result = traits.simulate_discrete(tree, states=["marine", "terrestrial"],
                                  switch=0.1, seed=1)   # a discrete state
```

## Continuous traits

A continuous trait **diffuses**. Give it a starting value and a rate, and it wanders down every branch, its variance growing in proportion to elapsed time. On its own that is **Brownian motion**, the null model of continuous trait evolution:

```python
# BM — a body size diffusing from 0 at variance-rate σ² = 1.0
traits.simulate_continuous(tree, start=0.0, rate=1.0, seed=1)
```

Here `rate` is the Brownian variance-rate σ², the trait level's reading of "how fast", and like every rate in ZOMBI2 it accepts modifiers. For example, two variations on this are:

```python
# OU — the same diffusion, pulled toward an optimum value
traits.simulate_continuous(tree, start=0.0, rate=1.0,
                           reverts_to=2.0, pull=0.5, seed=1)

# early burst — the diffusion rate itself decays through time
traits.simulate_continuous(tree, start=0.0,
                           rate=1.0 * mod.OnTime({0: 1.0, 5: 0.2}), seed=1)
```

The **Ornstein–Uhlenbeck** process is Brownian motion with a rubber band: `reverts_to` is the optimum it is pulled back toward, and `pull` is how hard. **Early burst** (or ACDC) is a diffusion rate that decays as the tree ages, so most of the divergence happens near the root; it is written with the same `mod.OnTime` that gives the species tree its skyline.

The rest of the modifier vocabulary applies to `rate` unchanged, each with a name in the comparative-methods literature: `mod.FromParent(spread=…)` makes σ² drift from parent to daughter (variable-rates BM, the trait twin of ClaDS), `mod.OnTotalDiversity(cap=…)` slows σ² as the clade fills. Two knobs sit alongside `rate`: `regimes=` paints a multi-optimum OU, where clades pull toward different optima (a discrete trait supplies the painting and `reverts_to` becomes one optimum per regime), and `at_speciation=` adds a jump *at* each split rather than along the branches. These are not alternatives to each other. `regimes=` and `at_speciation=` combine, and both combine with the modifiers on `rate`, so nothing here is a knob you have to give up to use another. The one exception is stated where it bites: `regimes=` takes a plain σ², not a modified one.

## Discrete traits

A discrete trait takes a finite set of states and switches between them along the branches — a continuous-time Markov chain, the field's **Mk model** [@lewis2001mk]:

```python
# Mk — habitat flips between two states at rate 0.1
traits.simulate_discrete(tree, states=["marine", "terrestrial"],
                         switch=0.1, start="marine", seed=1)
```

When the flips are not symmetric, replace the single rate with a small matrix of directed rates:

```python
# asymmetric — gains are commoner than losses
traits.simulate_discrete(tree, states=["absent", "present"],
                         switch={"absent->present": 0.2, "present->absent": 0.05},
                         seed=1)
```

Either shape may carry a modifier, so a switch rate can be **driven by another trait** grown first on the same tree: `switch = 0.2 * mod.DrivenBy(habitat, {"marine": 5.0, "terrestrial": 1.0})`. That is the conditioning of Chapter 9.

A **threshold** trait is the third case, and it is a bridge back to the continuous world. An observed discrete state can be driven by an underlying continuous **liability** that itself does Brownian motion; the state you see is which side of a threshold the liability currently sits on:

```python
# threshold — a discrete state read off an underlying continuous liability
traits.simulate_discrete(tree, states=["absent", "present"],
                         liability=1.0, threshold=0.0, seed=1)
```

## Correlated traits

Two traits that evolve independently are two separate calls, in either order. Two traits that drift *together* cannot be simulated one before the other, because each is entangled with the other as it unfolds. Correlation is specified as per-trait rates plus a correlation overlay:

```python
traits.simulate_continuous(tree,
    start={"size": 0.0, "limb": 0.0},
    rate={"size": 1.0, "limb": 0.8},        # one variance-rate per trait
    correlation={("size", "limb"): 0.6},    # the overlay, ∈ [−1, 1]
    seed=1)
```

The overlay is a dimensionless number in `[−1, 1]`, not a covariance matrix. Under `correlation=` the per-trait rates are plain numbers.

The overlay carries the reversion and the speciation jumps as well. Add `reverts_to` and `pull` and each trait reverts to its own optimum at its own strength; add `at_speciation` and the jump at each split is drawn under the same correlation the diffusion uses. Each argument takes one value shared across the traits, or a dict giving one per trait:

```python
traits.simulate_continuous(tree,
    start={"size": 0.0, "limb": 0.0},
    rate={"size": 1.0, "limb": 0.8},
    correlation={("size", "limb"): 0.6},
    reverts_to={"size": 3.0, "limb": -1.0},   # each trait its own optimum
    pull={"size": 1.5, "limb": 0.4},          # and its own strength
    at_speciation=0.5, seed=1)
```

The restriction is worth stating, because the literature's multivariate OU is larger. Here the reversion is **diagonal**: a trait's deviation pulls that trait back and no other, and everything the traits share lives in the diffusion. A drift matrix, where size being above its optimum pulls limb along with it, is a different model; passing one is refused by name rather than read as its diagonal.

The same overlay handles *discrete* correlation with no extra machinery, through the threshold model: give each trait a liability, correlate the liabilities, and put the thresholds on top. Correlated presence/absence characters, the setting Pagel's method was built for, are then one call:

```python
traits.simulate_discrete(tree, states=["absent", "present"],
    liability={"wings": 1.0, "flight": 1.0},
    correlation={("wings", "flight"): 0.7}, threshold=0.0, seed=1)
```

## What the modifiers reach

Brownian motion is where this chapter starts, not where it stops, and the way past it is not a longer menu of models. The diffusion rate is a rate like any other in ZOMBI2, so it takes the modifiers of Chapter 2, and most of the named alternatives to BM are one modifier on that rate:

```python
rate = 1.0 * mod.OnTime({0: 4.0, 1: 1.0})       # fast early, then settling — early burst
rate = 1.0 * mod.FromParent(spread=0.3)         # each clade inherits and drifts in tempo
rate = 1.0 * mod.OnTotalDiversity(cap=100)      # the rate eases off as the clade fills

# σ² reads a second trait, grown first on the same tree
habitat = traits.simulate_discrete(tree, states=["marine", "terrestrial"], switch=0.3, seed=1)
size = traits.simulate_continuous(tree, start=0.0, seed=2,
    rate=1.0 * mod.DrivenBy(habitat, {"marine": 4.0, "terrestrial": 1.0}))
```

The last is one trait driving another, which is the conditioning of Chapter 9: the driver is grown first, and the run that reads it comes second.

Two arguments of `simulate_continuous` go further in the same spirit. `reverts_to` and `pull` make the walk revert toward an optimum instead of wandering freely, and `at_speciation` adds a jump at every split, so change concentrates at branching rather than accumulating along branches.

None of these is a separate model with its own function and its own parameters, which is why they combine: a trait that bursts early *and* reverts to an optimum is one rate with one modifier and two arguments, not a model somebody had to implement.

```python
# a burst that decays, under a pull toward an optimum of 2
traits.simulate_continuous(tree, start=0.0, reverts_to=2.0, pull=0.5, seed=1,
                           rate=1.0 * mod.OnTime({0: 4.0, 1: 1.0}))
```

The combination is exact, not an approximation of one. A branch running from `t₀` to `t₁` ends normally distributed around `θ + (x−θ)·e^{−α(t₁−t₀)}` with variance `∫ e^{−2α(t₁−s)}·σ²(s) ds`, integrated across every point where the schedule, the standing diversity or a driver changes σ². The weight inside that integral is what mean reversion does to old variance: what the trait accrued early has been pulled back toward the optimum by the time the branch ends, so it counts for less than what it accrued late. Dropping the weight and using Brownian motion's `∫ σ²(s) ds` would overstate the variance by an order of magnitude on a typical branch. The table at the end of the chapter gives each combination its usual name in the literature.

## Models from the literature

Trait models arrive under a thicket of names, and a reader who wants "an OU model" or "a threshold model" should be able to find it. The names live here, in one table, and organise nothing else in the chapter.

| What it does | ZOMBI2 | From the literature |
|---|---|---|
| a value diffusing | `simulate_continuous(rate=…)` | Brownian motion (BM) [@felsenstein1985comparative] |
| diffusion pulled to an optimum | `simulate_continuous(rate=…, reverts_to=…, pull=…)` | Ornstein–Uhlenbeck (OU) [@hansen1997stabilizing; @butler2004phylogenetic] |
| diffusion rate decays through time | `simulate_continuous(rate=1.0 * mod.OnTime({…}))` | Early burst (EB / ACDC) [@harmon2010earlyburst] |
| diffusion rate drifts between lineages | `simulate_continuous(rate=1.0 * mod.FromParent(spread=…))` | Variable-rates BM [@maliet2019clads] |
| diffusion rate slows as the clade fills | `simulate_continuous(rate=1.0 * mod.OnTotalDiversity(cap=…))` | Diversity-dependent / ecological limits [@etienne2012diversitydependence] |
| the optimum differs between painted clades | `simulate_continuous(regimes=…, reverts_to={…}, pull=…)` | Multi-optimum OU (OUM) [@beaulieu2012ouwie] |
| the value jumps at each split | `at_speciation=…` (either kind) | Cladogenetic / punctuational change |
| traits evolving together | one `simulate_continuous(rate={…}, correlation={…})` call | Multivariate BM |
| traits reverting together, each to its own optimum | `simulate_continuous(rate={…}, correlation={…}, reverts_to={…}, pull=…)` | Multivariate OU, diagonal drift [@clavel2015mvmorph] |
| a discrete state switching | `simulate_discrete(states=…, switch=…)` | Mk (k-state Markov) |
| discrete driven by continuous liability | `simulate_discrete(liability=…, threshold=…)` | Threshold / liability (Wright–Felsenstein) [@felsenstein2012threshold] |
| discrete traits evolving together | `simulate_discrete(liability={…}, correlation={…})` | Correlated binary / Pagel [@pagel1994correlated] |

## The objects

A run returns a **`TraitsResult`** bundle:

- `.values` — the observable vector: the trait's value at each **extant tip**. This is the comparative-data matrix a method would be handed.
- `.node_values` — the value at **every** node (extant, extinct, and internal alike), the true ancestors at each split, from the same process that produced the tips.
- `.events` — the timestamped event log, the same shape as the genome level's: each entry is a change on a lineage at a time, from one state to another, and its `kind` is `on_branch` (a switch along a branch) or `on_speciation` (a jump at a split). For a discrete trait this log is the source of truth. A continuous trait diffuses with no along-branch events, so its log holds only the `at_speciation` jumps and is empty without them.
- `.history` — for a **discrete** trait, the per-branch stochastic character map derived from that log: the ordered list of `(state, duration)` segments each branch passed through. It is `None` for a continuous trait, which has no map, and for a threshold trait, whose liability crossings are un-timed.

For discrete traits the stored values are the state labels you gave (not integer indices), so `.values` and `.node_values` already read back in your own vocabulary.

## Usage from Python

```python
from zombi2 import species, traits
from zombi2.rates import modifiers as mod

# a species tree from the previous chapters, then a trait riding along it
# (the complete tree — a trait evolves on extinct lineages too)
tree = species.simulate_species_tree(
    birth=1.0, death=0.3, n_extant=30, seed=1).complete_tree

# continuous: body size under Brownian motion
size = traits.simulate_continuous(tree, start=0.0, rate=1.0, seed=1)
size.values                 # {extant tip: float}

# continuous: an Ornstein–Uhlenbeck trait pulled toward an optimum of 2
temp = traits.simulate_continuous(tree, start=0.0, rate=1.0,
                                  reverts_to=2.0, pull=0.5, seed=1)

# discrete: habitat flipping between two states
habitat = traits.simulate_discrete(tree, states=["marine", "terrestrial"],
                                   switch=0.1, start="marine", seed=1)
habitat.values              # {extant tip: "marine" | "terrestrial"}
habitat.events              # the realized flips, in time order

# two continuous traits that drift together — one joint call
bodyplan = traits.simulate_continuous(tree,
    start={"size": 0.0, "limb": 0.0},
    rate={"size": 1.0, "limb": 0.8},
    correlation={("size", "limb"): 0.6}, seed=1)
```

## Usage from the CLI

The state space is `--kind`, and it is required. It decides which of the other flags apply — `--rate` and the OU knobs belong to a continuous trait, `--states` and `--switch` to a discrete one — so there is no default that would not silently pick a model for you:

```bash
# a continuous (Brownian) trait along a species tree
zombi2 traits out/ --kind continuous \
    --start 0.0 --rate 1.0 --seed 1

# a discrete two-state trait
zombi2 traits out/ --kind discrete \
    --states marine,terrestrial --switch 0.1 --seed 1

# the same, also writing the driver file a conditioned genome or trait run reads (Chapter 9)
zombi2 traits out/ --kind discrete \
    --states cave,surface --switch 0.1 --seed 1 \
    --write values events tree
```

The trait evolves on the **complete** tree, extinct lineages included, so `species_complete.nwk` is the file to hand it. An external tree works too; if it is not ultrametric you must declare each tip's fate with `--tip-fates`, because ZOMBI will not guess which early-ending tips are extinct.

## Outputs

| File | What it holds |
|---|---|
| `trait_values.tsv` | the value at every node — extant tips, extinct lineages and internal nodes |
| `trait_events.tsv` | a `root` row for the state at t=0, then every realized switch with its time (discrete traits) |
| `trait_tree.nwk` | the complete tree with every node annotated `[&trait=…]`, which opens in FigTree or iTOL |

Because the value at every node comes from the same process that produced the tips, these carry the
*exact* ancestral states, not a reconstruction. `trait_events.tsv` is also the driver file a
conditioned genome or trait run reads (Chapter 9) — given the shared tree, the root state plus the switches
rebuild the trait on every lineage at every instant; Appendix B gives the columns and the formats.
