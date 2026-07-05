# Matching empirical profiles (ABC)

The rest of ZOMBI2 runs *forward*: given rates, it produces gene families. `match_profiles`
runs the other way — given an **empirical copy-number profile** (families × extant species)
and the species tree it was observed on, it searches for the duplication/transfer/loss/
origination rates under which the forward model reproduces that profile. It is a plain
**Approximate Bayesian Computation** (ABC) rejection sampler:

1. draw a rate set from your priors,
2. simulate a profile under it,
3. reduce the simulated and the empirical matrix to a vector of **summary statistics**, and
4. keep the draws whose summaries land closest to the empirical one.

## Match summaries, not the table

You cannot match a profile *table* directly — a simulation produces a different set of
families, in a different order, with different labels. So matching is done on
**permutation-invariant** summaries of the whole matrix. The default (see
[`default_summary`](../reference/api.md#profile-matching-abc)) is three blocks:

| Statistic | What it captures |
|---|---|
| **Gene frequency spectrum** | how many families are present in exactly 1, 2, … *S* species (the pangenome core/shell/cloud curve — the most informative single statistic) |
| **Genome sizes** | total gene copies per species |
| **Copy-number spectrum** | how many present cells hold 1, 2, 3, ≥4 copies (separates duplication from transfer) |

Each component is scaled by its spread across the simulated batch, so the distance is
scale-free without hand-tuned weights.

## Basic usage

```python
import zombi2 as z

tree = z.simulate_species_tree(z.BirthDeath(1.0, 0.3), n_tips=40, age=5.0, seed=1)

empirical = z.ProfileMatrix.from_tsv("profiles.tsv")   # or any ProfileMatrix

fit = z.match_profiles(
    tree, empirical,
    priors={"duplication": (0, 1.0), "transfer": (0, 0.5),
            "loss": (0, 1.5), "origination": (0, 3.0)},
    n_sims=4000, accept=0.02, seed=1,
)
print(fit)                 # median [95% CI] per rate
fit.summary()              # {rate: {mean, median, lo95, hi95}}
```

A prior is given per rate as a `(low, high)` tuple (uniform), a bare float (fixed), a
[`Distribution`](../reference/api.md#distributions), a scipy frozen distribution, or a
callable `rng -> float`. Rates you omit are held at 0. The `empirical` argument is a
`ProfileMatrix`, or a path / TSV string that `ProfileMatrix.from_tsv` can read.

!!! note "The tree sets the rate scale"
    Rates are per unit branch length, so the input tree must be the one the profiles were
    observed on, with meaningful branch lengths.

## Reading the result

`match_profiles` returns an [`ABCFit`](../reference/api.md#profile-matching-abc). Because
the rates are only partly identifiable from a profile (see below), the result is a
**posterior sample, not a point estimate**:

```python
fit.posterior              # {rate: array of accepted values}
fit.summary()              # per-rate mean / median / 95% credible interval
fit.best                   # the single closest-matching draw (one point on a ridge)
fit.plot_spectra()         # posterior-predictive check (needs matplotlib)
```

`plot_spectra()` overlays the empirical frequency spectrum on the spectra of the accepted
simulations. If the empirical curve sits inside the accepted band, the model can reproduce
your data; if it falls outside, no setting of these rates fits — a sign the model is
missing something. `spectra_data()` returns the same arrays if you would rather plot them
yourself.

## What is — and isn't — identifiable

!!! warning "Loss sits on a ridge"
    From copy number alone, the **gain-side rates (duplication, transfer, origination) are
    well identified**, but **loss is not**: a fully lost lineage leaves little observable
    trace, so a range of (origination, loss) pairs fit a profile equally well. Read the
    `summary()` credible intervals, not `best` alone, and expect the loss interval to be
    wide. More data (a bigger tree) sharpens the gain-side rates but does **not** remove the
    loss bias — it is systematic, not noise.

The [`profile_matching_experiment.py`](https://github.com/AADavin/zombi2/blob/main/examples/profile_matching_experiment.py)
and [`profile_matching_scaling.py`](https://github.com/AADavin/zombi2/blob/main/examples/profile_matching_scaling.py)
examples demonstrate the recovery and this identifiability limit end to end.

## Sharper estimates

Three tools push against the limits above. Their honest scope: they sharpen the
identifiable rates; none fully resolves the loss ridge.

### Regression adjustment

A post-hoc correction (Beaumont 2002) that regresses the accepted parameters on their
summary residuals and shifts them toward the empirical target — **no extra simulations**:

```python
fit.regression_adjust()            # {rate: adjusted accepted values}
fit.summary(adjusted=True)         # summary of the adjusted posterior
```

### Gene-tree information

If you have gene trees, pass an empirical [`Genomes`](../reference/api.md#simulation-driver)
and set `gene_trees=True`. The summary then also uses duplication/transfer/loss **event
counts**, which pin the gain-side rates sharply. Simulations then return the full genealogy
(`output="genomes"`) rather than the counts-only path:

```python
emp = z.simulate_genomes(tree, duplication=0.3, transfer=0.1, loss=0.5,
                         origination=1.5, initial_families=15, seed=1)   # a Genomes
fit = z.match_profiles(tree, emp, priors=priors, gene_trees=True, n_sims=2000, seed=1)
```

The three event-count columns would be drowned by the many profile columns, so
`default_gene_tree_summary` automatically up-weights them (via the general `feature_weights`
option on the distance).

### Sequential Monte Carlo

`match_profiles_smc` evolves a population of particles over rounds with a shrinking
tolerance (Toni et al. 2009), reaching a given tolerance with fewer simulations than plain
rejection. It returns a **weighted** `ABCFit` (`summary()` reports weighted intervals):

```python
fit = z.match_profiles_smc(
    tree, empirical,
    priors={"duplication": (0, 1.0), "transfer": (0, 0.5),
            "loss": (0, 1.5), "origination": (0, 3.0)},
    rounds=5, n_particles=200, seed=1,
)
fit.n_simulations          # total simulations used
```

SMC needs **uniform (or fixed) priors**. For copy-number profiles it typically matches
rejection (the posterior width is set by the data, not the tolerance); it pays off when the
model is sharply identified.

## Choosing the model

By default `match_profiles` fits the four scalar rates shared by every family (the built-in
model, which runs on the Rust engine). You can instead fit richer models:

```python
# per-family heterogeneous rates: fit the *means* of FamilySampledRates distributions
z.match_profiles(tree, emp, priors=priors, model="family", family_shape=2.0,
                 max_family_size=25, n_sims=2000, seed=1)

# any custom builder params -> RateModel
z.match_profiles(tree, emp, priors={"d": (0, 1), "l": (0, 1)},
                 model=lambda p: z.SharedRates(duplication=p["d"], loss=p["l"]), seed=1)
```

!!! tip "Cap growth on the Python engine"
    `model="family"` and custom models run on the Python engine, which builds an object per
    gene copy — a high-duplication draw can blow up. Pass `max_family_size=` to bound it.
    The default uniform model on the Rust engine is immune (it tracks integer counts).

## Speed

The inner loop is embarrassingly parallel and, for the default uniform model, runs on the
[Rust engine](rust-engine.md) automatically (the engine is chosen by the model — there is no
`engine` argument). `processes > 1` distributes the simulations across worker processes —
results are identical regardless of the count:

```python
fit = z.match_profiles(tree, empirical, priors=priors, n_sims=20_000,
                       processes=8, seed=1)     # call from a __main__ guard
```

## Custom summaries

Pass your own `statistics` — any callable `ProfileMatrix -> 1-D array` (or, with
`gene_trees=True`, `Genomes -> 1-D array`). It may expose a `feature_weights` attribute, or
you can pass `feature_weights=` directly, to weight components in the distance. The
[`frequency_spectrum`](../reference/api.md#profile-matching-abc), `genome_sizes`, and
`copy_number_spectrum` building blocks are exported for reuse.
