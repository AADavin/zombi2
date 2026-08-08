# Traits

The trait level evolves **phenotypes**: a body size, a habitat, the presence or absence of a structure. A trait evolves along the species tree like everything else here. There are two kinds, continuous and discrete.

![The two kinds, on the same tree. A **continuous** trait starts at one value and drifts along every branch, so each tip ends at a number of its own and close relatives end up close. A **discrete** trait sits in one state and switches to another now and then; the dots mark the two switches, and every branch below a switch is in the new state until it switches again. The chips at the tips are the state each one ended in. The two are different functions, `simulate_continuous` and `simulate_discrete`, because they answer different questions.](figures/trait_kinds_print.png){width=100%}

```python
from zombi2 import species, traits
from zombi2.rates import ScaledBy
from zombi2.rates import modifiers as mod

tree = species.simulate_species_tree(birth=1.0, death=0.3, n_extant=20, seed=1)
result = traits.simulate_continuous(tree, start=0.0, rate=1.0, seed=1)   # a real value
result = traits.simulate_discrete(tree, states=["marine", "terrestrial"],
                                  switch=0.1, seed=1)   # a discrete state
```

## Continuous traits

A continuous trait **diffuses**. Give it a starting value and a rate, and it wanders down every branch, its variance growing in proportion to elapsed time. On its own that is **Brownian motion**, the null model of continuous trait evolution:

```python
# BM: a body size diffusing from 0 at variance-rate σ² = 1.0
traits.simulate_continuous(tree, start=0.0, rate=1.0, seed=1)
```

Here `rate` is the Brownian variance-rate σ², the trait level's reading of "how fast", and like every rate in ZOMBI2 it accepts modifiers. For example, two variations on this are:

```python
# OU: the same diffusion, pulled toward an optimum value
traits.simulate_continuous(tree, start=0.0, rate=1.0,
                           reverts_to=2.0, pull=0.5, seed=1)

# early burst: the diffusion rate itself decays through time
traits.simulate_continuous(tree, start=0.0,
                           rate=1.0 * mod.OnTime({0: 1.0, 5: 0.2}), seed=1)
```

The **Ornstein–Uhlenbeck** process is Brownian motion with a rubber band: `reverts_to` is the optimum it is pulled back toward, and `pull` is how hard. **Early burst** (or ACDC) is a diffusion rate that decays as the tree ages, so most of the divergence happens near the root; it is written with the same `mod.OnTime` that gives the species tree its skyline.

That is the pattern for everything past Brownian motion here, and it is not a longer menu of models. The diffusion rate is a rate like any other in ZOMBI2, so the modifiers this level reads apply to it unchanged, each with a name in the comparative-methods literature:

```python
rate = 1.0 * mod.OnTime({0: 4.0, 1: 1.0})       # fast early, then settling: early burst
rate = 1.0 * mod.Inherited(per='lineage', spread=0.3)         # each clade inherits and drifts in tempo
rate = 1.0 * mod.OnTotalDiversity(cap=100)      # the rate eases off as the clade fills

# σ² reads a second trait, grown first on the same tree
habitat = traits.simulate_discrete(tree, states=["marine", "terrestrial"], switch=0.3, seed=1)
size = traits.simulate_continuous(tree, start=0.0, seed=2,
    rate=1.0 * ScaledBy(habitat, {"marine": 4.0, "terrestrial": 1.0}))
```

The last is one trait driving another, which is the conditioning of Chapter 9: the driver is grown first, and the run that reads it comes second.

`Drawn(per='lineage')`, Chapter 3's null for inherited tempo, is not among them: the level refuses it by name. Appendix A lists what each level accepts.

Two more knobs sit alongside `rate`. `regimes=` paints a multi-optimum OU, where clades pull toward different optima (a discrete trait supplies the painting and `reverts_to` becomes one optimum per regime), and `at_speciation=` adds a jump *at* each split rather than along the branches, so change concentrates at branching. The value is the jump variance, so `at_speciation=0.5` gives a jump of width √0.5. None of these is a separate model with its own function and its own parameters, which is why they combine: a trait that bursts early *and* reverts to an optimum is one rate with one modifier and two arguments.

```python
# a burst that decays, under a pull toward an optimum of 2
traits.simulate_continuous(tree, start=0.0, reverts_to=2.0, pull=0.5, seed=1,
                           rate=1.0 * mod.OnTime({0: 4.0, 1: 1.0}))
```

The combination is exact, not an approximation of one. A branch running from `t₀` to `t₁` ends normally distributed around `θ + (x−θ)·e^{−α(t₁−t₀)}` with variance `∫ e^{−2α(t₁−s)}·σ²(s) ds`, integrated across every point where the schedule, the standing diversity or a driver changes σ². The weight inside that integral is what mean reversion does to old variance: what the trait accrued early has been pulled back toward the optimum by the time the branch ends, so it counts for less than what it accrued late. Dropping the weight and using Brownian motion's `∫ σ²(s) ds` would overstate the variance — by about a quarter on an average branch of this tree at `pull=0.5`, and by more the longer the branch or the stronger the pull.

`regimes=` is the one knob that asks you to give things up, and it says so rather than ignoring them: it takes a plain σ² (not a modified one), one jump variance shared across regimes (not one per regime), and one trait (so not `correlation=`).

## Discrete traits

A discrete trait takes a finite set of states and switches between them along the branches: a continuous-time Markov chain, the field's **Mk model** [@lewis2001mk].

```python
# Mk: habitat flips between two states at rate 0.1
traits.simulate_discrete(tree, states=["marine", "terrestrial"],
                         switch=0.1, start="marine", seed=1)
```

When the flips are not symmetric, replace the single rate with a small matrix of directed rates:

```python
# asymmetric: gains are commoner than losses
traits.simulate_discrete(tree, states=["absent", "present"],
                         switch={"absent->present": 0.2, "present->absent": 0.05},
                         seed=1)
```

Either of those two shapes may carry a modifier, a bare rate or a `{'from->to': rate}` entry, though not the third spelling, a `k×k` matrix, whose entries are numbers by construction. So a switch rate can be **driven by another trait** grown first on the same tree: `switch = 0.2 * ScaledBy(habitat, {"marine": 5.0, "terrestrial": 1.0})`. That is the conditioning of Chapter 9.

`at_speciation=` works here too, and means something different: it is the probability, in `[0, 1]`, that a daughter hops at the split, to a state drawn uniformly from the others.

A **threshold** trait is the third case, and it is a bridge back to the continuous world. An observed discrete state can be driven by an underlying continuous **liability** that itself does Brownian motion; the state you see is which side of a threshold the liability currently sits on:

```python
# threshold: a discrete state read off an underlying continuous liability
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

## Literature

Trait models arrive under a thicket of names, and a reader who wants "an OU model" or "a threshold model" should be able to find it. The names live here, in one table, and organise nothing else in the chapter.

| What it does | ZOMBI2 | From the literature |
|---|---|---|
| a value diffusing | `simulate_continuous(rate=…)` | Brownian motion (BM) [@felsenstein1985comparative] |
| diffusion pulled to an optimum | `simulate_continuous(rate=…, reverts_to=…, pull=…)` | Ornstein–Uhlenbeck (OU) [@hansen1997stabilizing; @butler2004phylogenetic] |
| diffusion rate decays through time | `simulate_continuous(rate=1.0 * mod.OnTime({…}))` | Early burst (EB / ACDC) [@harmon2010earlyburst] |
| diffusion rate drifts between lineages | `simulate_continuous(rate=1.0 * mod.Inherited(per='lineage', spread=…))` | Variable-rates BM [@maliet2019clads] |
| diffusion rate slows as the clade fills | `simulate_continuous(rate=1.0 * mod.OnTotalDiversity(cap=…))` | Diversity-dependent / ecological limits [@etienne2012diversitydependence] |
| the optimum differs between painted clades | `simulate_continuous(regimes=…, reverts_to={…}, pull=…)` | Multi-optimum OU (OUM) [@beaulieu2012ouwie] |
| the value jumps at each split | `at_speciation=…` (continuous, and Mk `switch=` traits, not threshold ones) | Cladogenetic / punctuational change |
| traits evolving together | one `simulate_continuous(rate={…}, correlation={…})` call | Multivariate BM |
| traits reverting together, each to its own optimum | `simulate_continuous(rate={…}, correlation={…}, reverts_to={…}, pull=…)` | Multivariate OU, diagonal drift [@clavel2015mvmorph] |
| a discrete state switching | `simulate_discrete(states=…, switch=…)` | Mk (k-state Markov) |
| discrete driven by continuous liability | `simulate_discrete(liability=…, threshold=…)` | Threshold / liability (Wright–Felsenstein) [@felsenstein2012threshold] |
| discrete traits evolving together | `simulate_discrete(liability={…}, correlation={…})` | Correlated binary / Pagel [@pagel1994correlated] |

## The objects

A run returns a **`TraitsResult`** bundle:

- `.values`, the observable vector: the trait's value at each **extant tip**, keyed by the tip's name (`n5`), the same names the Newick and `trait_values.tsv` use, so the dataset joins the tree it came from. `.values_by_id` is the same thing keyed by bare node id, for joining against `.node_values`. This is the comparative-data matrix a method would be handed.
- `.node_values`, the value at **every** node (extant, extinct, and internal alike), the true ancestors at each split, from the same process that produced the tips.
- `.events`, the timestamped event log — one row per change, like the genome level's, but with its own columns: each entry is a change on a lineage at a time, from one state to another, and its `kind` is `initial` (the value at t=0), `on_branch` (a switch along a branch) or `on_speciation` (a jump at a split). For an Mk (`switch=`) trait this log is the source of truth. A continuous trait diffuses with no along-branch events, so its log holds the `initial` row and any `at_speciation` jumps. The `initial` row is written for a continuous and an Mk trait, so a plain Brownian run's log has exactly one entry. A threshold trait keeps no log at all: its state is read off a continuous liability, so there are no timed crossings to record and `.events` is empty — its `.node_values` is the record.
- `.history`, for a **discrete** trait, the per-branch stochastic character map derived from that log: the ordered list of `(state, duration)` segments each branch passed through. It is `None` for a continuous trait, which has no map, and for a threshold trait, whose liability crossings are un-timed.

For discrete traits the stored values are the state labels you gave (not integer indices), so `.values` and `.node_values` already read back in your own vocabulary.

A trait evolves on the **complete** tree, extinct lineages included, so that is the tree to hand it:

```python
tree = species.simulate_species_tree(
    birth=1.0, death=0.3, n_extant=30, seed=1).complete_tree
habitat = traits.simulate_discrete(tree, states=["marine", "terrestrial"],
                                   switch=0.1, start="marine", seed=1)
habitat.values              # {"n5": "marine" | "terrestrial", …}, keyed like the tree's tips
habitat.events              # the realized flips, in time order
```

## On the command line

The state space is `--kind`, and it is required. It decides which of the other flags apply, since `--rate` and the OU knobs belong to a continuous trait while `--states` and `--switch` belong to a discrete one, so there is no default that would not silently pick a model for you:

```bash
# a continuous (Brownian) trait along a species tree
zombi2 traits out/ --kind continuous \
    --start 0.0 --rate 1.0 --seed 1

# a discrete two-state trait
zombi2 traits out/ --kind discrete \
    --states marine,terrestrial --switch 0.1 --seed 1

# the states the next command drives from — a discrete run writes trait_events.tsv by default,
# the driver file a conditioned genome, sequence or trait run reads (Chapter 9)
zombi2 traits out/ --kind discrete \
    --states cave,surface --switch 0.1 --seed 1

# a second trait, in its own directory, whose switch rate reads the first (Chapter 9)
zombi2 traits out/ --kind discrete --name diet --states plant,fish --seed 2 \
    --switch "0.2 * ScaledBy('out/traits/trait_events.tsv', {'cave': 5.0, 'surface': 1.0})"
```

Every rate flag takes a rate in its written form, `--switch` as much as `--rate`, so the expression above is the same text the Python API takes. `--switch` reads the other two shapes its keyword does as well: a `{'a->b': rate}` dict and a `k x k` matrix.

Two keywords have no flag. `correlation=` grows several traits in one call, and the command line grows one trait per run. `regimes=` takes a discrete result object, the painting, handed to the call. Both stay in the Python API.

The trait evolves on the **complete** tree, extinct lineages included, so `species_complete.nwk` is the file to hand it. An external tree works too; if it is not ultrametric you must declare each tip's fate with `--tip-fates`, because ZOMBI will not guess which early-ending tips are extinct.

## Outputs

| File | What it holds |
|---|---|
| `trait_values.tsv` | the value at every node: extant tips, extinct lineages and internal nodes |
| `trait_events.tsv` | an `initial` row for the state at t=0, then every realized switch with its time (discrete traits) |
| `trait_tree.nwk` | the complete tree with every node annotated `[&trait=…]`, which opens in FigTree or iTOL |

Because the value at every node comes from the same process that produced the tips, these carry the
*exact* ancestral states, not a reconstruction. `trait_events.tsv` is also the driver file a
conditioned genome, sequence or trait run reads (Chapter 9): given the shared tree, the initial
state plus the switches rebuild the trait on every lineage at every instant. Appendix B gives the
columns and the formats.
