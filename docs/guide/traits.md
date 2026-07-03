# Trait evolution

Once you have a tree, you can evolve **traits** on it — a body size, an expression level, a
discrete character such as habitat or the presence of a structure. ZOMBI2 simulates the classic
phylogenetic-comparative models (Felsenstein 1985 and successors) with one function,
`simulate_traits`, and a family of model objects.

```python
import zombi2 as z

tree = z.simulate_species_tree(z.BirthDeath(1.0, 0.3), n_tips=30, age=5.0, seed=1)
result = z.simulate_traits(tree, z.BrownianMotion(sigma2=0.5), seed=1)

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
z.simulate_traits(tree, z.BrownianMotion(sigma2=0.5, x0=0.0, trend=0.0), seed=1)
```

### Ornstein–Uhlenbeck

Stabilizing selection: the trait is pulled toward an optimum `theta` with strength `alpha`
while it diffuses (`dX = alpha·(theta − X)·dt + σ·dW`). It starts at the optimum unless you
set `x0`.

```python
z.simulate_traits(tree, z.OrnsteinUhlenbeck(sigma2=0.4, alpha=2.0, theta=10.0), seed=1)
```

### Early burst / ACDC

Brownian motion whose rate changes exponentially through time, `σ²(t) = sigma2·e^{rate·t}`.
`rate < 0` is an *early burst* (most divergence happens early — an adaptive radiation),
`rate > 0` accelerates, `rate = 0` is plain Brownian motion.

```python
z.simulate_traits(tree, z.EarlyBurst(sigma2=1.0, rate=-0.8), seed=1)
```

## Correlated continuous traits

A vector-valued trait with a rate (covariance) matrix `R` couples its dimensions — the model
of correlated evolution (`mvMORPH`). Each node's value is a length-`k` array.

```python
R = [[1.0, 0.9],
     [0.9, 1.0]]                       # strong positive correlation between the two dimensions
z.simulate_traits(tree, z.MultivariateBrownian(R), seed=1)

# multivariate OU: pull each dimension toward an optimum (alpha may be a scalar, vector, or matrix)
z.simulate_traits(tree, z.MultivariateOU(R, alpha=1.5, theta=[0.0, 5.0]), seed=1)
```

## Discrete traits

### The Mk model

A continuous-time Markov chain over `k` states with a rate matrix `Q` (Lewis 2001). Convenience
constructors cover the standard sub-models; the raw constructor takes any `Q` (all-rates-different,
or an ordered/meristic character).

```python
mk = z.simulate_traits(tree, z.Mk.equal_rates(3, 0.4,               # ER: one shared rate
                                               states=["marine", "brackish", "fresh"]),
                       seed=2)

mk.labeled_values()                    # {extant leaf: "marine" | "brackish" | "fresh"}
mk.history[node]                       # [(state, duration), ...] — the stochastic map
mk.changes()                           # the transition events

z.Mk.symmetric([[0, 2, 1], [2, 0, 3], [1, 3, 0]])   # SYM: symmetric Q
```

Every `Mk` also exposes its analytic quantities:

```python
m = z.Mk.equal_rates(3, 0.4)
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
# Y tracks X: Y is gained quickly when X = 1 and lost quickly when X = 0
m = z.CorrelatedBinary(x_gain_y0=0.5, x_gain_y1=0.5, x_loss_y0=0.5, x_loss_y1=0.5,
                       y_gain_x0=0.05, y_gain_x1=2.0, y_loss_x0=2.0, y_loss_x1=0.05)
res = z.simulate_traits(tree, m, seed=1)
res.labeled_values()                   # {extant leaf: (X, Y)}

null = z.CorrelatedBinary.independent(x_gain=0.5, x_loss=0.5, y_gain=0.5, y_loss=0.5)
```

### Hidden rate classes (corHMM)

`HiddenStateMk` gives the observed character **hidden rate classes**: its transition rates
depend on an unobserved class that itself switches along the tree, capturing rate heterogeneity
a plain `Mk` cannot (Beaulieu et al. 2013). Tips report the observed state; `full_label` and
`changes()` expose the hidden dimension.

```python
slow = [[0, 0.1], [0.1, 0]]
fast = [[0, 3.0], [3.0, 0]]
hmm = z.HiddenStateMk(observed_rates=[slow, fast], hidden_rate=0.5,
                      observed_states=[0, 1], hidden_states=["slow", "fast"])
res = z.simulate_traits(tree, hmm, seed=1)
res.labeled_values()                   # observed 0/1 (hidden class collapsed)
res.full_label(res.node_values[tree.extant_leaves()[0]])   # (observed, hidden), e.g. (1, 'fast')
```

### The threshold model

An unobserved continuous **liability** evolves by Brownian motion; the observed discrete state
is the interval it falls in, cut by an ordered set of thresholds (Felsenstein 2012). The
evolving value is the liability; the observed state comes from `labeled_values()`.

```python
th = z.simulate_traits(tree, z.ThresholdModel(thresholds=[0.0]), seed=1)   # binary
th.values                              # liabilities (continuous, latent)
th.labeled_values()                    # observed 0/1 states
```

## Adaptation to different regimes (multi-optimum OU)

Different parts of the tree can adapt toward different optima (`OUwie`). The **regimes** come
from a discrete stochastic map — simulate them with `Mk` on the *same* tree, then run an OU with
one `theta` per regime. `alpha` and `sigma2` may also vary by regime.

```python
regimes = z.simulate_traits(tree, z.Mk.equal_rates(2, 0.4), seed=1)        # paint 2 regimes
mou = z.MultiOptimumOU(regimes, theta=[-5.0, 5.0], alpha=4.0, sigma2=0.4)
z.simulate_traits(tree, mou, seed=2)                                       # tips track their optimum
```

## Pagel's tree transforms

Pagel's (1999) λ, κ and δ transform the tree's branch/node lengths; run any trait model on the
result to model departures from a strict clock.

```python
z.simulate_traits(z.pagel_lambda(tree, 0.5), z.BrownianMotion(0.5), seed=1)   # scale signal
z.pagel_delta(tree, 2.0)     # node depths ^delta (>1 late, <1 early change)
z.pagel_kappa(tree, 0.0)     # branch lengths ^kappa (0 = speciational: unit branches)
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
dec = z.DEC(areas=["A", "B", "C"], dispersal=0.1, extinction=0.1, max_range_size=3)
res = z.simulate_biogeography(tree, dec, root_state={"A"}, seed=1)

res.labeled_values()                   # {extant leaf: ('A', 'B') ...} — the observed ranges
res.ancestral_states()                 # ancestral ranges at every internal node
res.changes()                          # anagenetic dispersal / extinction events along branches
```

`dispersal` may be a scalar or an area-by-area matrix, `extinction` a scalar or per-area vector,
and `max_range_size` caps how many areas a range may span.

## Output

```python
res = z.simulate_traits(tree, z.BrownianMotion(0.5), seed=1)
print(res.to_tsv())                    # node<TAB>trait, one row per extant tip
print(res.to_newick())                 # Newick with [&trait=…] on every node
res.to_tsv(nodes="all")                # include internal (ancestral) nodes too
```

Continuous traits are written as their value, multivariate traits as `{a,b,c}`, and discrete
traits as their state label.
