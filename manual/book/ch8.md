# Traits

The trait level evolves **phenotypes**: a body size, a habitat, the presence or absence of a structure. A trait evolves along the species tree like everything else here. There are two kinds, continuous and discrete.

![The two kinds, on the same tree. A **continuous** trait starts at one value and drifts along every branch, so every node ends at a number of its own — each branch is painted with the value at its far end — and close relatives end up close. A **discrete** trait sits in one state and switches to another now and then; the dots mark the two switches, and every branch below a switch is in the new state until it switches again. The chips at the tips are the state each one ended in. The two are different functions, `simulate_continuous` and `simulate_discrete`, because they answer different questions.](figures/trait_kinds_print.png){width=100%}

## Continuous traits

A continuous trait **diffuses**. Give it a starting value and a rate, and it wanders down every branch, its variance growing in proportion to elapsed time. On its own that is **Brownian motion**, the null model of continuous trait evolution:

```python
from zombi2 import species, traits
from zombi2.params import PerLineage

tree = species.simulate_species_tree(birth=1.0, death=0.3, n_extant=20, seed=1)
```

```python
# BM: a body size diffusing from 0 at variance-rate σ² = 1.0
traits.simulate_continuous(tree, start=0.0, rate=1.0, seed=1)
```

Here `rate` is the Brownian variance-rate σ², the trait level's reading of "how fast", and like every rate in ZOMBI2 it takes the verbs. For example, two variations on this are:

```python
# OU: the same diffusion, pulled toward an optimum value
traits.simulate_continuous(tree, start=0.0, rate=1.0,
                           reverts_to=2.0, pull=0.5, seed=1)

# early burst: the diffusion rate itself decays through time
traits.simulate_continuous(tree, start=0.0,
                           rate=PerLineage(1.0).changing_at({0: 1.0, 5: 0.2}), seed=1)
```

The **Ornstein–Uhlenbeck** process is Brownian motion with a rubber band: `reverts_to` is the optimum it is pulled back toward, and `pull` is how hard. **Early burst** (or ACDC) is a diffusion rate that decays as the tree ages, so most of the divergence happens near the root; it is written with the same `changing_at` that gives the species tree its skyline.

Two more arguments sit alongside `rate`. `regimes=` paints a multi-optimum OU, where clades pull toward different optima (a discrete trait supplies the painting and `reverts_to` becomes one optimum per regime), and `at_speciation=` adds a jump *at* each split rather than along the branches, so change concentrates at branching. The value is the jump variance, so `at_speciation=0.5` gives a jump of width √0.5. None of these is a separate model with its own function and its own parameters, which is why they combine: a trait that bursts early *and* reverts to an optimum is one rate with one verb and two arguments.

`regimes=` is the one argument that asks you to give things up, and it says so rather than ignoring them: it takes a plain σ² (not a modified one), one jump variance shared across regimes (not one per regime), and one trait (so not `correlation=`).

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

## Literature

Trait models arrive under a thicket of names, and a reader who wants "an OU model" or "a threshold model" should be able to find it. The names live here, in one table, and organise nothing else in the chapter.

| What it does | ZOMBI2 | From the literature |
|---|---|---|
| a value diffusing | `simulate_continuous(rate=…)` | Brownian motion (BM) [@felsenstein1985comparative] |
| diffusion pulled to an optimum | `simulate_continuous(rate=…, reverts_to=…, pull=…)` | Ornstein–Uhlenbeck (OU) [@hansen1997stabilizing; @butler2004phylogenetic] |
| diffusion rate decays through time | `simulate_continuous(rate=PerLineage(1.0).changing_at({…}))` | Early burst (EB / ACDC) [@harmon2010earlyburst] |
| diffusion rate drifts between lineages | `simulate_continuous(rate=PerLineage(1.0).varying_among('lineages', Drift(LogNormal(0.0, …))))` | Variable-rates BM [@maliet2019clads] |
| diffusion rate slows as the clade fills | `simulate_continuous(rate=PerLineage(1.0).scaled_by(TotalDiversity(cap=…)))` | Diversity-dependent / ecological limits [@etienne2012diversitydependence] |
| the optimum differs between painted clades | `simulate_continuous(regimes=…, reverts_to={…}, pull=…)` | Multi-optimum OU (OUM) [@beaulieu2012ouwie] |
| the value jumps at each split | `at_speciation=…` (continuous, and Mk `switch=` traits, not threshold ones) | Cladogenetic / punctuational change |
| traits evolving together | one `simulate_continuous(rate={…}, correlation={…})` call | Multivariate BM |
| traits reverting together, each to its own optimum | `simulate_continuous(rate={…}, correlation={…}, reverts_to={…}, pull=…)` | Multivariate OU, diagonal drift [@clavel2015mvmorph] |
| a discrete state switching | `simulate_discrete(states=…, switch=…)` | Mk (k-state Markov) |
| discrete driven by continuous liability | `simulate_discrete(liability=…, threshold=…)` | Threshold / liability (Wright–Felsenstein) [@felsenstein2012threshold] |
| discrete traits evolving together | `simulate_discrete(liability={…}, correlation={…})` | Correlated binary / Pagel [@pagel1994correlated] |

## What a run gives back

A run returns a `TraitsResult`. **`.values` is the observable vector** — the trait at each *extant
tip*, keyed by the tip's name (`n5`), the same names the Newick and `trait_values.tsv` use, so the
dataset joins the tree it came from. **`.node_values` is every node**, extant, extinct and internal
alike: the true ancestors at each split, from the same process that produced the tips. For a discrete
trait both read back in the state labels you gave, not integer indices. Appendix B lists the rest.

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

The state space is `--kind`, and it is required. It decides which of the other flags apply, since `--rate` and the OU flags belong to a continuous trait while `--states` and `--switch` belong to a discrete one, so there is no default that would not silently pick a model for you:

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
    --switch "PerLineage(0.2).scaled_by('out/traits/trait_events.tsv', {'cave': 5.0, 'surface': 1.0})"
```

Every rate flag takes a rate in its written form, `--switch` as much as `--rate`, so the expression above is the same text the Python API takes. `--switch` reads the other two shapes its keyword does as well: a `{'a->b': rate}` dict and a `k x k` matrix.

Two keywords have no flag. `correlation=` grows several traits in one call, and the command line grows one trait per run. `regimes=` takes a discrete result object, the painting, handed to the call. Both stay in the Python API.

The trait evolves on the **complete** tree, extinct lineages included, so `species_complete.nwk` is the file to hand it. An external tree works too; if it is not ultrametric you must declare each tip's fate with `--tip-fates`, because ZOMBI will not guess which early-ending tips are extinct.
