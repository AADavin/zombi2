# Trait evolution

Once you have a tree, you can evolve *traits* on it: a body size, an expression level, a discrete
character such as habitat or the presence of a structure. ZOMBI2 simulates the classic phylogenetic
comparative models with a single driver, `simulate_traits`, and a family of model objects that plug
into it.

```python
from zombi2.species import BirthDeath, simulate_species_tree
from zombi2.traits import (
    simulate_traits, BrownianMotion, OrnsteinUhlenbeck, EarlyBurst,
    MultivariateBrownian, MultivariateOU, Mk, CorrelatedBinary,
    CorrelatedBinaryK, HiddenStateMk, ThresholdModel, MultiOptimumOU,
    DEC, simulate_biogeography,
    pagel_lambda, pagel_delta, pagel_kappa,
)

tree = simulate_species_tree(BirthDeath(1.0, 0.3), n_tips=30, age=5.0, seed=1)
result = simulate_traits(tree, BrownianMotion(sigma2=0.5), seed=1)

result.values                 # {extant leaf: value} — the observable tip data
result.ancestral_states()     # {internal node: value} — exact, not inferred
```

The trait works on any `Tree`: a simulated species tree, or a gene tree loaded from disk. Because
ZOMBI2 simulates the process forward, the ancestral states it returns are the true history, not a
reconstruction.

## The shared overlay engine

Every model, continuous or discrete, runs on one engine. A value is dropped at the root and evolved
branch by branch in pre-order traversal, each node inheriting its parent's end state. The two model
families differ only in how a branch is traversed:

- *Continuous* models draw the exact branch endpoint. Along a branch of length $t$, the change is a
  single normal draw whose mean and variance the model prescribes; nothing between the endpoints is
  simulated because nothing needs to be.
- *Discrete* models simulate the Markov jumps exactly. The engine steps through the continuous-time
  chain, drawing waiting times and transitions, so the realized history along each branch, a
  *stochastic character map*, comes for free rather than being sampled after the fact.

![Two ways a branch is traversed. **A**, continuous models draw the branch endpoint in one step — the change over a branch of length $t$ is a single Normal draw, and the path in between is never simulated. **B**, discrete models simulate the Markov jumps exactly, so the realized history along the branch (a stochastic character map) comes for free.](figures/trait_overlay.pdf){width=100%}

`simulate_traits` returns a `TraitResult`. Its accessors expose the process at whatever granularity
you need:

| Access | Meaning |
| --- | --- |
| `result.values` | Values at the extant tips (the observable comparative data) |
| `result.labeled_values()` | Same, with discrete state indices decoded to their labels |
| `result.ancestral_states()` | Values at every internal node (exact ancestral states) |
| `result.node_values` | Every node: root, internal, and tips |
| `result.history` | Per-branch `(state, duration)` segments, the stochastic map (discrete only) |
| `result.changes()` | Realized transitions `(node, time, from, to)` (discrete only) |
| `result.to_tsv()` / `result.to_newick()` | Tip table / annotated Newick (`[&trait=…]`) |

## Trait models at a glance

| Model | Kind | What it captures |
|---|---|---|
| `BrownianMotion` | continuous | neutral random walk |
| `OrnsteinUhlenbeck` | continuous | stabilizing selection toward an optimum |
| `EarlyBurst` | continuous | a rate that changes through time (adaptive radiation) |
| `MultivariateBrownian` / `MultivariateOU` | continuous | correlated multi-trait evolution |
| `Mk` | discrete | a $k$-state Markov character |
| `CorrelatedBinary` | discrete | two binary characters evolving jointly (Pagel) |
| `CorrelatedBinaryK` | discrete | $k$ binary characters evolving jointly (Pagel, $k \ge 2$) |
| `HiddenStateMk` | discrete | hidden rate classes (corHMM) |
| `ThresholdModel` | discrete | a discrete state from a latent continuous liability |
| `MultiOptimumOU` | continuous | OU with a different optimum per regime |
| `DEC` | range | geographic-range evolution (dispersal / extinction / cladogenesis) |

## Continuous traits

### Brownian motion

Brownian motion is the reference model of quantitative-character evolution [@felsenstein1985comparative]:
a random walk, $dX = \mathrm{trend}\cdot dt + \sigma\, dW$. The tips are jointly multivariate-normal
with mean $x_0 + \mathrm{trend}\cdot\mathrm{depth}$ and covariance $\mathrm{sigma2}\cdot C$, where $C$
is the shared-path-length matrix, so more closely related tips are more strongly correlated.

```python
simulate_traits(tree, BrownianMotion(sigma2=0.5, x0=0.0, trend=0.0), seed=1)
```

![A continuous trait wanders down the branches by Brownian motion; each tip value is the endpoint of its root-to-tip random walk. Branches are painted by the trait value (viridis).](figures/trait_bm.pdf){width=100%}

`sigma2` is the diffusion rate, `x0` the root value, and `trend` a directional drift added per unit
time.

### Ornstein–Uhlenbeck

The Ornstein–Uhlenbeck process adds stabilizing selection [@hansen1997stabilizing; @butler2004phylogenetic]:
the trait diffuses but is pulled toward an optimum `theta` with strength `alpha`,
$dX = \alpha\,(\theta - X)\,dt + \sigma\, dW$. It starts at the optimum unless you set `x0`. The pull
bounds the variance, so lineages cluster around `theta` instead of wandering freely.

```python
ou = OrnsteinUhlenbeck(sigma2=0.4, alpha=2.0, theta=10.0)
simulate_traits(tree, ou, seed=1)
```

![Ornstein–Uhlenbeck: the trait diffuses but is pulled back toward an optimum, so lineages cluster around it.](figures/trait_ou.pdf){width=100%}

### Early burst / ACDC

`EarlyBurst` is Brownian motion whose rate changes exponentially through time,
$\sigma^2(t) = \mathrm{sigma2}\cdot e^{\,\mathrm{rate}\cdot t}$ [@harmon2010earlyburst]. A negative
`rate` is an *early burst*, in which most divergence happens early, the signature of an adaptive
radiation; a positive `rate` accelerates through time; `rate = 0` recovers plain Brownian motion.

```python
simulate_traits(tree, EarlyBurst(sigma2=1.0, rate=-0.8), seed=1)
```

![A trait under early burst (ACDC). The rate decays through time (grey strip: fast early, slow late), so most change happens early — the deep branches carry big divergences (clades split apart) while young clades stay similar. Branches are painted by the trait value (viridis).](figures/trait_earlyburst.pdf){width=100%}

### Correlated continuous traits

A vector-valued trait with a rate (covariance) matrix `R` couples its dimensions, the standard model
of correlated multivariate evolution [@clavel2015mvmorph]. Each node's value is a length-$k$ array,
and off-diagonal entries of `R` make the dimensions evolve together.

```python
R = [[1.0, 0.9],
     [0.9, 1.0]]      # strong positive correlation between the dimensions
simulate_traits(tree, MultivariateBrownian(R), seed=1)

# multivariate OU: pull each dimension toward an optimum
# (alpha may be a scalar, vector, or matrix)
simulate_traits(tree, MultivariateOU(R, alpha=1.5, theta=[0.0, 5.0]), seed=1)
```

![Multivariate Brownian motion. Two traits diffuse jointly under a rate matrix $R$; each dot is one simulated endpoint, the ring is the 95% covariance ellipse, and the lines are sample paths. **Left**, a correlated $R$ tilts the cloud so the traits move together; **right**, an independent $R$ gives an upright circle.](figures/trait_multivariate.pdf){width=100%}

## Discrete traits

### The Mk model

A discrete character evolves as a continuous-time Markov chain over $k$ states with a rate matrix `Q`
[@lewis2001mk]. Convenience constructors cover the standard sub-models, and the raw constructor takes
any `Q`, so you can specify an all-rates-different or an ordered/meristic character.

```python
# equal-rates Mk over 3 labeled states (one shared rate)
mk = simulate_traits(
    tree, Mk.equal_rates(3, 0.4, states=["marine", "brackish", "fresh"]), seed=2)

mk.labeled_values()   # {extant leaf: "marine" | "brackish" | "fresh"}
mk.history[node]      # [(state, duration), ...] — the stochastic map
mk.changes()          # the transition events

Mk.symmetric([[0, 2, 1], [2, 0, 3], [1, 3, 0]])   # SYM: symmetric Q
Mk.ordered(4, 0.5)   # ordered: adjacent-only chain
Mk([[0, 1, 2], [3, 0, 1], [1, 1, 0]])   # ARD: any user-supplied Q
```

![A discrete character under the Mk model. **A**, the model — an equal-rates chain over three states (all the transition structure lives in the rate matrix $Q$). **B**, a realization: the exact stochastic character map, each branch painted in the state it is in and switching at the instant of a transition; the tip chips give each lineage's observed state.](figures/trait_mk.pdf){width=100%}

The transition structure lives entirely in the `Q` you pass. `equal_rates` is all-to-all at one
rate, `symmetric` makes $i\to j$ and $j\to i$ equal, `ordered` is the tridiagonal nearest-neighbour
chain, and the raw constructor takes an arbitrary Markov chain. Every `Mk` also exposes its analytic
quantities:

```python
m = Mk.equal_rates(3, 0.4)
m.transition_matrix(1.0)               # P(t) = exp(Q·t)
m.stationary_distribution()            # pi with pi Q = 0
```

The root state is `"uniform"` by default; pass `root="stationary"`, an index, or a probability
vector.

### Correlated binary characters (Pagel)

`CorrelatedBinary` evolves two binary characters **X** and **Y** jointly, one changing at a time
[@pagel1994correlated]. Each trait's gain and loss rate may depend on the *other* trait's state, and
that dependence is exactly what correlated evolution means. The `independent` constructor builds the
null model in which the two evolve separately, giving the comparison that tests for correlation.

```python
# Y tracks X: Y is gained quickly when X = 1 and lost quickly when X = 0
m = CorrelatedBinary(
    x_gain_y0=0.5, x_gain_y1=0.5, x_loss_y0=0.5, x_loss_y1=0.5,
    y_gain_x0=0.05, y_gain_x1=2.0, y_loss_x0=2.0, y_loss_x1=0.05)
res = simulate_traits(tree, m, seed=1)
res.labeled_values()                   # {extant leaf: (X, Y)}

null = CorrelatedBinary.independent(
    x_gain=0.5, x_loss=0.5, y_gain=0.5, y_loss=0.5)
```

![Correlated binary characters. **A**, the joint state space $\{00, 01, 10, 11\}$: horizontal edges flip Y, vertical edges flip X, and each arrow's width is its rate. With these rates Y is gained fast only when X = 1 and lost fast when X = 0, so Y is driven toward X. **B**, one realization on a tree: branches are drawn heavy where X = 1 and light where X = 0, and the two tip columns give each leaf's $(X, Y)$ — the clades where X switches on are exactly the clades where Y is present.](figures/trait_pagel.pdf){width=100%}

`CorrelatedBinaryK` generalises this to **$k \ge 2$** binary characters evolving jointly over the
$2^k$ configurations of $\{0,1\}^k$, still changing one trait at a time (every transition that would
flip two or more bits at once has rate zero — Pagel's single-change constraint). Each single-bit flip
may depend on the *other* traits' states, and that dependence is the correlated evolution. It is a
plain $2^k$-state `Mk`, so it simulates like any discrete trait, with each tip labelled by its
$k$-tuple. Build it with the compact convenience constructors — `independent` (the Pagel **null**,
each trait on its own), `equal_rates` (one shared gain/loss), `partner_coupling` (each trait modulated
by one designated partner), or `from_table` (the full rate table):

```python
# 3 binary traits; trait 1 "tracks" trait 0 as its partner: when
# trait 0 is present, trait 1 is gained 3x faster and lost half as
# fast, so the two come to co-occur (traits 0 and 2 evolve freely)
m = CorrelatedBinaryK.partner_coupling(
    gains=[0.5, 0.3, 0.5], losses=[0.5, 0.5, 0.5],
    partners=[None, 0, None], boost_gain=[1, 3, 1], boost_loss=[1, 0.5, 1])
res = simulate_traits(tree, m, seed=1)
res.labeled_values()                   # {extant leaf: (t0, t1, t2)}

# the k = 3 Pagel null
null = CorrelatedBinaryK.independent(gains=0.5, losses=[0.5, 0.5, 0.5])
```

### Hidden rate classes (corHMM)

`HiddenStateMk` gives the observed character *hidden rate classes*: its transition rates depend on an
unobserved class that itself switches along the tree, capturing rate heterogeneity that a plain `Mk`
cannot represent [@beaulieu2013hidden]. Tips report the observed state; `full_label` and `changes()`
expose the hidden dimension.

```python
slow = [[0, 0.1], [0.1, 0]]
fast = [[0, 3.0], [3.0, 0]]
hmm = HiddenStateMk(observed_rates=[slow, fast], hidden_rate=0.5,
                      observed_states=[0, 1], hidden_states=["slow", "fast"])
res = simulate_traits(tree, hmm, seed=1)
res.labeled_values()                   # observed 0/1 (hidden class collapsed)
res.full_label(res.node_values[tree.extant_leaves()[0]])   # (observed, hidden)
```

![Hidden rate classes. **A**, the four joint states, (observed 0/1) $\times$ (hidden slow/fast): the observed character flips slowly in the slow class and quickly in the fast class (arrow width = rate), while the hidden class itself switches along the tree. **B**, a realization — branches are heavy where the hidden class is fast and light where it is slow, small circles mark observed-state changes, and the tip chips give the observed state. The observed changes cluster on the fast (heavy) branches, the signature of hidden rate classes.](figures/trait_hiddenmk.pdf){width=100%}

### The threshold model

In the threshold model an unobserved continuous *liability* evolves by Brownian motion, and the
observed discrete state is the interval the liability falls in, cut by an ordered set of thresholds
[@felsenstein2012threshold]. The evolving value is the liability; the observed state comes from
`labeled_values()`. Because the underlying process is continuous, the model links quantitative and
discrete evolution.

```python
th = simulate_traits(tree, ThresholdModel(thresholds=[0.0]), seed=1)   # binary
th.values                              # liabilities (continuous, latent)
th.labeled_values()                    # observed 0/1 states
```

![The threshold model. A latent liability evolves by Brownian motion — the branches are painted by its value (viridis) — and the observed discrete state at each tip is which side of the threshold the liability ends on: open square = below (state 0), filled = above (state 1). The threshold is marked on the colour bar, so colours to its left are below and to its right above. The continuous gradient flowing into a binary tip pattern is the model's point.](figures/trait_threshold.pdf){width=100%}

Discrete simulation, whether Mk, correlated-binary, hidden-rate, or the DEC ranges below, always
yields a full stochastic character map [@nielsen2002mapping; @huelsenbeck2003stochastic]; `history`
and `changes()` read it back.

## Adaptation to different regimes (multi-optimum OU)

Different parts of the tree can adapt toward different optima, the multi-optimum OU of `OUwie`
[@beaulieu2012ouwie]. The *regimes* come from a discrete stochastic map: simulate them with `Mk` on
the *same* tree, then run an OU with one `theta` per regime. `alpha` and `sigma2` may also vary by
regime.

```python
regimes = simulate_traits(tree, Mk.equal_rates(2, 0.4), seed=1)   # 2 regimes
mou = MultiOptimumOU(regimes, theta=[-5.0, 5.0], alpha=4.0, sigma2=0.4)
simulate_traits(tree, mou, seed=2)   # tips track their optimum
```

This is the natural way to model a shift in selective regime: a discrete character paints where the
optimum changes, and the continuous trait chases whichever optimum its lineage currently sits under.

![Multi-optimum OU. A discrete character paints two regimes on the tree (the tip bars: regime A dark, regime B light); the continuous trait then evolves under OU with a different optimum per regime. Branches are painted by the trait value (viridis) — regime-A lineages are pulled toward the low optimum, regime-B lineages toward the high one, so the tips split into two clusters.](figures/trait_multioptimum.pdf){width=100%}

## Pagel's tree transforms

Pagel's (1999) $\lambda$, $\kappa$ and $\delta$ rescale a tree's branch and node lengths *before* you
run a trait model on it [@pagel1999inferring]. Each warps the **tempo of trait evolution** — how much
expected trait change accrues along each branch — so a Brownian or OU trait can depart from a
constant, strictly time-proportional rate. (This is the 1999 tree-transform family, distinct from the
1994 correlated-character model earlier in this chapter, which also carries Pagel's name.) Each is a
*deterministic*, one-parameter map `Tree` $\to$ `Tree`; the stochasticity still comes from whichever
trait model you then evolve on the transformed tree.

```python
# reshape trait tempo, then run any model on the transformed tree
simulate_traits(pagel_lambda(tree, 0.5), BrownianMotion(0.5), seed=1)
pagel_delta(tree, 2.0)     # node depths ^delta (>1 late, <1 early change)
pagel_kappa(tree, 0.0)     # branch lengths ^kappa (0 = speciational)
```

- **$\lambda$** scales internal (shared) depths while holding tip depths fixed: `1` is the original
  tree (full phylogenetic signal), `0` a star tree (independent tips). It tunes how much of a
  Brownian trait's covariance is shared versus tip-specific — something no molecular clock does.
- **$\delta$** raises node depths to a power, root and tips fixed: values above `1` concentrate
  change late (toward the tips), below `1` early (toward the root).
- **$\kappa$** raises each branch length to a power; `0` gives a speciational model — change accrues
  per speciation event rather than per unit time.

![Pagel's 1999 tree transforms, each shown on one small tree at three parameter values. **lambda** scales the internal (shared) depths while pinning the tips, interpolating from the original tree to a star; **delta** bends node depths toward the root (early change) or the present (late change); **kappa** raises branch lengths to a power, collapsing toward equal, speciational branches at `0`. The transform is applied to the tree; the trait model is then run on the result.](figures/pagel_transforms.pdf){width=100%}

::: note
These are **not** the molecular clocks of Chapter 15. Pagel's transforms reshape the timetree that a
*continuous- or discrete-trait* model (BM/OU/Mk) evolves over — a deterministic warp of *trait
tempo*. The relaxed clocks of Chapter 15 instead turn a timetree into a *substitution* phylogram
(expected substitutions per site) to drive *sequence* simulation, drawing random per-branch rates.
Different object transformed, different model consuming it, different units — neither can express the
other.
:::

## Historical biogeography (DEC)

A species' "trait" can be its *geographic range*, a subset of discrete areas. The
Dispersal–Extinction–Cladogenesis model evolves the range by *dispersal* (gaining an area) and local
*extinction* (losing one) along branches, plus a *cladogenetic* split of the ancestral range between
daughters at each speciation, by narrow sympatry, subset sympatry, or vicariance [@ree2008dec].
Because of that node process, DEC has its own driver, `simulate_biogeography`.

```python
dec = DEC(areas=["A", "B", "C"], dispersal=0.1, extinction=0.1, max_range_size=3)
res = simulate_biogeography(tree, dec, root_state={"A"}, seed=1)

res.labeled_values()      # {extant leaf: ('A', 'B') ...} — tip ranges
res.ancestral_states()    # ancestral ranges at every internal node
res.changes()             # anagenetic dispersal / extinction events
```

![The DEC model of geographic-range evolution. **A**, along a branch a range gains areas by dispersal and loses them by local extinction — a Markov chain over the lattice of ranges (subsets of areas). **B**, at each speciation the ancestral range is split between the two daughters by narrow sympatry, subset sympatry, or vicariance.](figures/dec.pdf){width=100%}

`dispersal` may be a scalar or an area-by-area matrix, `extinction` a scalar or a per-area vector, and
`max_range_size` caps how many areas a range may span.

![A DEC history simulated on a tree. The root begins in area A; along the branches the range gains areas by dispersal (green triangles) and loses them by local extinction (red crosses), and at each speciation the ancestral range is passed to the daughters. Each tip's realized range is drawn as a row of area cells.](figures/dec_tree.pdf){width=100%}

## From the command line

The `zombi2 traits` command runs a single model on a tree you provide and writes the tip and ancestral
values:

```bash
T=species_tree.nwk

zombi2 traits -t $T --model bm --sigma2 0.5 --seed 1 -o out/
zombi2 traits -t $T --model ou --alpha 2.0 --theta 10.0 -o out/
zombi2 traits -t $T --model eb --rate -0.8 -o out/
zombi2 traits -t $T --model mk --states 3 -o out/         # equal-rates Mk
zombi2 traits -t $T --model mk --ordered -o out/          # adjacent-only chain
zombi2 traits -t $T --model mk --q-matrix Q.tsv -o out/   # arbitrary Q
zombi2 traits -t $T --model threshold -o out/
zombi2 traits -t $T --model dec --areas A,B,C --dispersal 0.3 -o out/
```

For the Mk model, `--model mk` is equal-rates by default, `--ordered` gives the adjacent-only
(tridiagonal) chain, and `--q-matrix FILE` reads an arbitrary `Q`. `--replicates N` repeats the
simulation on the same tree.

::: tip
To simulate a multi-optimum OU or any regime-dependent model, first paint the regimes with a discrete
`Mk` run on the tree, then pass that `TraitResult` to `MultiOptimumOU`. The Python API is the way to
compose the two stages; the CLI runs one model at a time.
:::

## Output

Continuous traits are written as their value, multivariate traits as `{a,b,c}`, and discrete traits
as their state label. `to_tsv()` writes a `node<TAB>trait` table, and `to_newick()` a Newick with
`[&trait=…]` on every node; pass `nodes="all"` to include the internal (ancestral) nodes.
