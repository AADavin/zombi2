# Relaxed molecular clocks

Every tree ZOMBI2 produces is a **timetree**: branch lengths are amounts of *time*. Sequence
evolution, though, does not accumulate time — it accumulates *substitutions*, and the number of
substitutions along a branch is its duration multiplied by an evolutionary **rate**. If that rate
were the same everywhere the two would be interchangeable, but real rates vary from branch to
branch: some lineages evolve fast, some slow, and the pattern drifts over the tree. A **relaxed
molecular clock** is the model of that variation. Applying one rescales every branch from time into
expected substitutions per site and turns a **chronogram** (branch lengths in time) into a
**phylogram** (branch lengths in substitutions).

Chapter 15 met one relaxed clock already — the shared lineage clock inside `SequenceEvolution`. This
chapter steps back and presents the whole **family** of clocks as first-class models with one common
interface, so you can rescale *any* tree (a species tree, a gene tree, anything read from Newick),
compare models, and plug whichever you like into a sequence-evolution run.

## The strict clock and what relaxing it means

The baseline is the **strict clock**: one rate on every branch. It multiplies the whole tree by a
constant, so the phylogram is just the chronogram uniformly stretched — the relative branch
proportions are untouched. `StrictClock(rate=1.0)` therefore reproduces the input tree exactly:

```python
import zombi2 as z

tree = z.simulate_species_tree(z.BirthDeath(1.0, 0.3), n_tips=30, age=8.0, seed=1)

z.StrictClock(rate=1.0).scale(tree, seed=1).to_newick() == tree.to_newick()   # True
```

`scale` is the shared entry point of every clock. It takes a `Tree` and a seed and returns a
`RateScaledTree` — the phylogram, plus the rates it drew. A **relaxed** clock is anything that lets
that per-branch rate vary. The clocks differ only in *how* the rate is drawn, and they split into two
kinds:

- **Uncorrelated** clocks draw an **independent** multiplier for each branch. A branch's rate tells
  you nothing about its neighbours' — rate is a property of the branch alone.
- **Autocorrelated** clocks anchor each branch's rate to its parent's, so the rate evolves *along*
  the tree and nearby lineages have similar rates.

## Drawing a phylogram

Every clock is used the same way. Construct it, call `scale`, and read the result:

```python
clock = z.UncorrelatedLogNormalClock(sigma=0.5)   # each branch an i.i.d. multiplier
phylogram = clock.scale(tree, seed=2)

phylogram.to_newick()        # the rescaled tree; branch lengths in substitutions/site
phylogram.branch_rate        # {node: the rate multiplier applied to that node's branch}
phylogram.branch_lengths     # {node: the substitution length of that node's branch}
```

The result is a valid tree in its own right; `to_newick()` writes it out, and it round-trips through
`read_newick` like any other. Because a clock is seeded, the same clock and seed always give the same
phylogram:

```python
a = clock.scale(tree, seed=7).to_newick()
b = clock.scale(tree, seed=7).to_newick()
assert a == b
```

A clock works on a **gene** tree just as well as a species tree — anything with time-valued branch
lengths. Load a reconciled gene tree with `read_newick` and scale it directly.

## The uncorrelated clocks

These are the PhyloBayes-style relaxed clocks [@drummond2006relaxed]. Each branch draws its rate
independently, from a distribution centred (in the first two cases) on mean rate 1, so on average the
tree's total length is preserved.

**Uncorrelated lognormal** — the workhorse. Each branch draws its rate as
$\exp\!\big(\mathcal{N}(-\sigma^2/2,\ \sigma)\big)$, a lognormal with mean 1. The single parameter
`sigma` is the amount of heterogeneity: `sigma = 0` is the strict clock, larger values spread the
rates wider.

```python
z.UncorrelatedLogNormalClock(sigma=0.5).scale(tree, seed=2)
```

**Uncorrelated gamma** — the same idea with a gamma distribution: `rate ~ Gamma(shape, 1/shape)`, mean
1, variance `1/shape`. Here a *large* `shape` means *less* variation (rates concentrate near 1); a
small shape spreads them out.

```python
z.UncorrelatedGammaClock(shape=3.0).scale(tree, seed=2)
```

**White noise** — an uncorrelated clock whose variance depends on branch *length*. The rate is the
integral of a white-noise process over the branch, so it is gamma-distributed with mean 1 and
variance $\sigma^2/\Delta t$ inversely proportional to the branch duration $\Delta t$. Long branches
average the noise away and sit near rate 1; short branches are highly variable. This branch-length
dependence is what sets it apart from the plain gamma clock above.

```python
z.WhiteNoiseClock(sigma=0.5).scale(tree, seed=2)
```

::: note
The uncorrelated clocks are memoryless across branches. Two branches meeting at a node can have very
different rates — there is no smoothing. If you expect rates to change gradually down the tree,
reach for an autocorrelated clock instead.
:::

## The autocorrelated clocks

Here a branch inherits its rate from its parent, so the rate is smooth along the tree.

**Autocorrelated lognormal** — a geometric random walk down the tree,
$R_{\text{child}} = R_{\text{parent}} \cdot \exp\!\big(\mathcal{N}(0,\ \sigma\sqrt{\ell})\big)$, where
$\ell$ is the branch length in time [@thorne1998autocorrelated]. A child's rate is centred on its
parent's, and the variance grows with elapsed time. `sigma = 0` freezes the walk into a strict clock
at `root_rate`. This is exactly the clock `SequenceEvolution`'s `branch_sigma` selects (see
Chapter 15), now available as a model in its own right.

```python
z.AutocorrelatedLogNormalClock(sigma=0.3).scale(tree, seed=2)
```

**Cox–Ingersoll–Ross (CIR)** — a mean-reverting diffusion of the rate itself
[@lepage2007general]. The instantaneous rate follows

$$dr = \theta\,(\mu - r)\,dt + \sigma\sqrt{r}\;dW,$$

which stays strictly positive and is pulled back toward the long-run mean $\mu$ (the `mean` argument)
at speed $\theta$ (`theta`); $\sigma$ (`sigma`) is the volatility. Unlike the lognormal walk the rate
also varies *within* a branch — the path is simulated by sub-stepping, so a long branch is split into
several rate segments — and, because it reverts to the mean, the tree's total length stays close to
the mean times the total time rather than drifting.

```python
z.CIRClock(theta=1.0, sigma=0.4, mean=1.0).scale(tree, seed=2)
```

**Discrete-bin (GTDB)** — the model of Chapter 15's `RateVariation`, also an autocorrelated clock: an
ordered set of rate bins with a nearest-neighbour Markov walk between them along the tree. It is
included here for completeness and shares the same interface.

```python
z.RateVariation(bins=[0.25, 0.5, 1.0, 2.0, 4.0], switch_rate=1.0).scale(tree, seed=2)
```

The whole family, at a glance:

| clock | class | kind | main parameter |
|-------|-------|------|----------------|
| strict | `StrictClock` | — | `rate` |
| uncorrelated lognormal | `UncorrelatedLogNormalClock` | uncorrelated | `sigma` |
| uncorrelated gamma | `UncorrelatedGammaClock` | uncorrelated | `shape` |
| white noise | `WhiteNoiseClock` | uncorrelated | `sigma` |
| autocorrelated lognormal | `AutocorrelatedLogNormalClock` | autocorrelated | `sigma` |
| Cox–Ingersoll–Ross | `CIRClock` | autocorrelated | `theta`, `sigma` |
| discrete-bin (GTDB) | `RateVariation` | autocorrelated | `bins`, `switch_rate` |

They all live in the `zombi2.clocks` namespace (and at the top level as `z.StrictClock`, and so on),
sharing the `Clock` interface: `scale(tree, seed=...)` returns a `RateScaledTree`.

## Feeding a sequence-evolution run

The clocks are the **shared lineage clock** of Chapter 15's gene × lineage model. `SequenceEvolution`
took `branch_sigma` (the autocorrelated lognormal) or a `RateVariation`; it now takes *any* clock via
`lineage=`, so you can drive a genome simulation's gene trees with any model in the family:

```python
genomes = z.simulate_genomes(tree, duplication=0.2, transfer=0.1, loss=0.2,
                             origination=0.5, seed=1)

se = z.SequenceEvolution(lineage=z.UncorrelatedGammaClock(shape=3.0),  # shared lineage clock
                         family_speed=z.LogNormal(0.0, 0.4))           # per-family speed
phylo = se.scale(genomes, seed=2)
```

Everything Chapter 15 says about the lineage clock still holds: it is drawn once and shared by every
family, then multiplied by each family's own speed. Passing `branch_sigma` remains a shorthand for
`lineage=AutocorrelatedLogNormalClock` at the same drift $\sigma$.

## Usage from the CLI

The `zombi2 sequence` command grows a `--clock` selector for the family. As in Chapter 15 it runs on
a prior `genomes` result whose event trace was written (`--write trace`):

```bash
zombi2 genomes -t species_tree.nwk --dup 0.2 --trans 0.1 --loss 0.2 --orig 0.5 \
    --write trace profiles -o run/

# an uncorrelated lognormal lineage clock, plus per-family speed
zombi2 sequence --genomes run/ --clock uncorrelated-lognormal --clock-sigma 0.5 \
    --family-speed 0.5 -o run/

# a mean-reverting CIR clock instead
zombi2 sequence --genomes run/ --clock cir --clock-theta 1.0 --clock-sigma 0.4 -o run/
```

`--clock` chooses the model — `strict`, `autocorrelated-lognormal`, `uncorrelated-lognormal`,
`uncorrelated-gamma`, `white-noise`, `cir`, or `discrete-bin`. Its parameter is supplied by
`--clock-sigma` (the spread of the lognormal, white-noise and CIR clocks), `--clock-shape` (the gamma
shape), `--clock-theta` (the CIR mean-reversion speed), and `--clock-mean` (the target/strict/root
rate, default 1); the discrete-bin clock reads its bins from `--branch-bins`. The output is exactly as
in Chapter 15 — one phylogram per family under `run/gene_trees/`, plus `gene_family_speeds.tsv` and
`branch_rates.tsv` recording the drawn rates.

The historical flags still work: `--branch-speed SIGMA` is the autocorrelated lognormal clock and
`--branch-bins R1,R2,...` the discrete-bin one, so old command lines are unchanged.

::: tip
Because `sequence` only replays the event trace on disk, one expensive `genomes` run feeds any number
of clock experiments. Sweep `--clock` and its parameters — strict, then a couple of uncorrelated
models, then CIR — reusing the same gene content, and compare the phylograms.
:::

::: note
The rescaled branch lengths already carry the clock, so adding `--subst-model` (Chapter 15) to
simulate an alignment needs no extra rate: one unit of the rescaled branch length is one expected
substitution per site by construction.
:::
