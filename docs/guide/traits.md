# Trait evolution

> **Reference:** see the [continuous-trait](../models/continuous-traits.md),
> [discrete-trait](../models/discrete-traits.md) and
> [biogeography](../models/biogeography.md) catalog pages.

Once you have a tree, you can evolve **traits** on it — a body size, an expression level, a
discrete character such as habitat or the presence of a structure. ZOMBI2 simulates the classic
phylogenetic-comparative models (Felsenstein 1985 and successors) with one function,
`simulate_traits`, and a family of model objects.

```python
from zombi2.species import simulate_species_tree, BirthDeath
from zombi2.traits import simulate_traits, BrownianMotion

tree = simulate_species_tree(BirthDeath(1.0, 0.3), n_tips=30, age=5.0, seed=1)
result = simulate_traits(tree, BrownianMotion(sigma2=0.5), seed=1)

result.values                 # {extant leaf: value} — the observable tip data
result.ancestral_states()     # {internal node: value} — exact, not inferred
```

All models share one engine: the trait starts at the root and is evolved branch by branch in
pre-order, each node inheriting its parent's end state. **Continuous** models draw the exact
branch endpoint; **discrete** models simulate the Markov jumps exactly, so the realized history
along each branch — a *stochastic character map* — comes for free. It works on any `Tree`: a
simulated species tree, or a gene tree loaded with [`read_newick`](../reference/api.md#tree).

## The result

`simulate_traits` returns a `TraitResult`:

| Access | Meaning |
| --- | --- |
| `result.values` | Values at the **extant** tips (the observable comparative data) |
| `result.labeled_values()` | Same, with discrete state indices decoded to their labels |
| `result.ancestral_states()` | Values at every internal node (exact ancestral states) |
| `result.node_values` | Every node (root, internal, tips) |
| `result.history` | Per-branch `(state, duration)` segments — the stochastic map (discrete only) |
| `result.changes()` | Realized transitions `(node, time, from, to)` (discrete only) |
| `result.to_tsv()` / `result.to_newick()` | Tip table / annotated Newick (`[&trait=…]`) |

## Continuous traits

### Brownian motion

A random walk: `dX = trend·dt + σ·dW`. The tips are jointly multivariate-normal with mean
`x0 + trend·depth` and covariance `sigma2 · C`, where `C` is the shared-path-length matrix.

```python
from zombi2.traits import simulate_traits, BrownianMotion

simulate_traits(tree, BrownianMotion(sigma2=0.5, x0=0.0, trend=0.0), seed=1)
```

<figure markdown="span">
![A continuous trait evolving by Brownian motion along a tree](../img/trait_bm.svg)
<figcaption>A continuous trait wandering down the branches by Brownian motion — the value at
each tip is the endpoint of its root-to-tip random walk.</figcaption>
</figure>

### Ornstein–Uhlenbeck

Stabilizing selection: the trait is pulled toward an optimum `theta` with strength `alpha`
while it diffuses (`dX = alpha·(theta − X)·dt + σ·dW`). It starts at the optimum unless you
set `x0`.

```python
from zombi2.traits import simulate_traits, OrnsteinUhlenbeck

simulate_traits(tree, OrnsteinUhlenbeck(sigma2=0.4, alpha=2.0, theta=10.0), seed=1)
```

<figure markdown="span">
![An Ornstein–Uhlenbeck trait pulled toward an optimum](../img/trait_ou.svg)
<figcaption>Ornstein–Uhlenbeck: the trait diffuses but is pulled back toward an optimum, so
lineages cluster around it instead of wandering freely as under Brownian motion.</figcaption>
</figure>

### Early burst / ACDC

Brownian motion whose rate changes exponentially through time, `σ²(t) = sigma2·e^{rate·t}`.
`rate < 0` is an *early burst* (most divergence happens early — an adaptive radiation),
`rate > 0` accelerates, `rate = 0` is plain Brownian motion.

```python
from zombi2.traits import simulate_traits, EarlyBurst

simulate_traits(tree, EarlyBurst(sigma2=1.0, rate=-0.8), seed=1)
```

## Correlated continuous traits

A vector-valued trait with a rate (covariance) matrix `R` couples its dimensions — the model
of correlated evolution (`mvMORPH`). Each node's value is a length-`k` array.

```python
from zombi2.traits import simulate_traits, MultivariateBrownian, MultivariateOU

R = [[1.0, 0.9],
     [0.9, 1.0]]                       # strong positive correlation between the two dimensions
simulate_traits(tree, MultivariateBrownian(R), seed=1)

# multivariate OU: pull each dimension toward an optimum (alpha may be a scalar, vector, or matrix)
simulate_traits(tree, MultivariateOU(R, alpha=1.5, theta=[0.0, 5.0]), seed=1)
```

<figure markdown="span">
![Two correlated continuous traits evolving jointly](../img/trait_multivariate.svg)
<figcaption>A multivariate trait: the rate matrix <em>R</em> correlates the dimensions, so the
two traits tend to move together down the tree.</figcaption>
</figure>

## Discrete traits

### The Mk model

A continuous-time Markov chain over `k` states with a rate matrix `Q` (Lewis 2001). Convenience
constructors cover the standard sub-models; the raw constructor takes any `Q` (all-rates-different,
or an ordered/meristic character).

```python
from zombi2.traits import simulate_traits, Mk

mk = simulate_traits(tree, Mk.equal_rates(3, 0.4,                  # ER: one shared rate
                                          states=["marine", "brackish", "fresh"]),
                     seed=2)

mk.labeled_values()                    # {extant leaf: "marine" | "brackish" | "fresh"}
mk.history[node]                       # [(state, duration), ...] — the stochastic map
mk.changes()                           # the transition events

Mk.symmetric([[0, 2, 1], [2, 0, 3], [1, 3, 0]])   # SYM: symmetric Q
Mk.ordered(4, 0.5)                                 # ordered: adjacent-only steps (i <-> i±1)
Mk([[0, 1, 2], [3, 0, 1], [1, 1, 0]])             # ARD: any user-supplied rate matrix
```

<figure markdown="span">
![A discrete character changing state along a tree under an Mk model](../img/trait_mk.svg)
<figcaption>A discrete character under an Mk model: the state jumps along the branches (the
exact stochastic character map), and the tips inherit whatever state they end in.</figcaption>
</figure>

The transition structure is entirely in the `Q` you pass: `equal_rates` is all-to-all at one
rate, `symmetric` makes `i→j` and `j→i` equal, `ordered` is the tridiagonal nearest-neighbour
chain (the character-state analogue of [`RateVariation`](sequences.md)), and the raw
constructor takes an arbitrary Markov chain. From the CLI, `--model mk` is equal-rates by
default, `--ordered` gives the adjacent-only chain, and `--q-matrix FILE` reads an arbitrary `Q`.

Every `Mk` also exposes its analytic quantities:

```python
from zombi2.traits import Mk

m = Mk.equal_rates(3, 0.4)
m.transition_matrix(1.0)               # P(t) = exp(Q·t)
m.stationary_distribution()            # π with πQ = 0
```

The root state is `"uniform"` by default; pass `root="stationary"`, an index, or a probability
vector.

### Correlated binary characters (Pagel 1994)

Two binary characters **X** and **Y** evolving jointly, one changing at a time. Each trait's
gain/loss rate may depend on the other trait's state — that dependence *is* correlated evolution.
`CorrelatedBinary.independent(...)` builds the null model where they evolve independently.

```python
from zombi2.traits import simulate_traits, CorrelatedBinary

# Y tracks X: Y is gained quickly when X = 1 and lost quickly when X = 0
m = CorrelatedBinary(x_gain_y0=0.5, x_gain_y1=0.5, x_loss_y0=0.5, x_loss_y1=0.5,
                     y_gain_x0=0.05, y_gain_x1=2.0, y_loss_x0=2.0, y_loss_x1=0.05)
res = simulate_traits(tree, m, seed=1)
res.labeled_values()                   # {extant leaf: (X, Y)}

null = CorrelatedBinary.independent(x_gain=0.5, x_loss=0.5, y_gain=0.5, y_loss=0.5)
```

### Hidden rate classes (corHMM)

`HiddenStateMk` gives the observed character **hidden rate classes**: its transition rates
depend on an unobserved class that itself switches along the tree, capturing rate heterogeneity
a plain `Mk` cannot (Beaulieu et al. 2013). Tips report the observed state; `full_label` and
`changes()` expose the hidden dimension.

```python
from zombi2.traits import simulate_traits, HiddenStateMk

slow = [[0, 0.1], [0.1, 0]]
fast = [[0, 3.0], [3.0, 0]]
hmm = HiddenStateMk(observed_rates=[slow, fast], hidden_rate=0.5,
                    observed_states=[0, 1], hidden_states=["slow", "fast"])
res = simulate_traits(tree, hmm, seed=1)
res.labeled_values()                   # observed 0/1 (hidden class collapsed)
res.full_label(res.node_values[tree.extant_leaves()[0]])   # (observed, hidden), e.g. (1, 'fast')
```

### The threshold model

An unobserved continuous **liability** evolves by Brownian motion; the observed discrete state
is the interval it falls in, cut by an ordered set of thresholds (Felsenstein 2012). The
evolving value is the liability; the observed state comes from `labeled_values()`.

```python
from zombi2.traits import simulate_traits, ThresholdModel

th = simulate_traits(tree, ThresholdModel(thresholds=[0.0]), seed=1)   # binary
th.values                              # liabilities (continuous, latent)
th.labeled_values()                    # observed 0/1 states
```

<figure markdown="span">
![A threshold model: a latent liability crossing a threshold flips the observed state](../img/trait_threshold.svg)
<figcaption>The threshold model: a latent liability evolves by Brownian motion (top), and the
observed discrete state is whichever interval it falls in at each point (bottom).</figcaption>
</figure>

## Adaptation to different regimes (multi-optimum OU)

Different parts of the tree can adapt toward different optima (`OUwie`). The **regimes** come
from a discrete stochastic map — simulate them with `Mk` on the *same* tree, then run an OU with
one `theta` per regime. `alpha` and `sigma2` may also vary by regime.

```python
from zombi2.traits import simulate_traits, Mk, MultiOptimumOU

regimes = simulate_traits(tree, Mk.equal_rates(2, 0.4), seed=1)            # paint 2 regimes
mou = MultiOptimumOU(regimes, theta=[-5.0, 5.0], alpha=4.0, sigma2=0.4)
simulate_traits(tree, mou, seed=2)                                         # tips track their optimum
```

## Pagel's tree transforms

Pagel's (1999) λ, κ and δ transform the tree's branch/node lengths; run any trait model on the
result to model departures from a strict clock.

```python
from zombi2.traits import simulate_traits, BrownianMotion, pagel_lambda, pagel_delta, pagel_kappa

simulate_traits(pagel_lambda(tree, 0.5), BrownianMotion(0.5), seed=1)      # scale signal
pagel_delta(tree, 2.0)       # node depths ^delta (>1 late, <1 early change)
pagel_kappa(tree, 0.0)       # branch lengths ^kappa (0 = speciational: unit branches)
```

- **λ** scales internal (shared) depths while holding tip depths fixed: `1` is the original
  tree (full phylogenetic signal), `0` a star tree (independent tips).
- **δ** raises node depths to a power (root and tips fixed).
- **κ** raises each branch length to a power; `0` gives a speciational model (change per
  speciation event, not per unit time).

## Historical biogeography (DEC)

A species' "trait" can be its **geographic range** — a subset of discrete areas. The
Dispersal–Extinction–Cladogenesis model (Ree & Smith 2008) evolves the range by *dispersal*
(gaining an area) and local *extinction* (losing one) along branches, plus a **cladogenetic**
split of the ancestral range between daughters at each speciation (narrow sympatry, subset
sympatry, or vicariance). Because of the node process it has its own driver,
`simulate_biogeography`.

```python
from zombi2.traits import DEC, simulate_biogeography

dec = DEC(areas=["A", "B", "C"], dispersal=0.1, extinction=0.1, max_range_size=3)
res = simulate_biogeography(tree, dec, root_state={"A"}, seed=1)

res.labeled_values()                   # {extant leaf: ('A', 'B') ...} — the observed ranges
res.ancestral_states()                 # ancestral ranges at every internal node
res.changes()                          # anagenetic dispersal / extinction events along branches
```

<figure markdown="span">
![Geographic range evolution under the DEC model](../img/dec.svg)
<figcaption>Historical biogeography (DEC): a lineage's range gains and loses discrete areas by
dispersal and local extinction, and the ancestral range is split between daughters at each
speciation.</figcaption>
</figure>

`dispersal` may be a scalar or an area-by-area matrix, `extinction` a scalar or per-area vector,
and `max_range_size` caps how many areas a range may span.

## Output

```python
from zombi2.traits import simulate_traits, BrownianMotion

res = simulate_traits(tree, BrownianMotion(0.5), seed=1)
print(res.to_tsv())                    # node<TAB>trait, one row per extant tip
print(res.to_newick())                 # Newick with [&trait=…] on every node
res.to_tsv(nodes="all")                # include internal (ancestral) nodes too
```

Continuous traits are written as their value, multivariate traits as `{a,b,c}`, and discrete
traits as their state label.
