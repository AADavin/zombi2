# Gene families & rates

Along the fixed species tree, ZOMBI2 runs one forward continuous-time (Gillespie) process
over all co-existing branches, with four core events:

| Event | Effect |
|---|---|
| **Origination (O)** | a brand-new family appears (one copy) on a branch |
| **Duplication (D)** | a gene copy splits into two |
| **Transfer (T)** | a copy is gained by another lineage alive at that time |
| **Loss (L)** | a copy is removed |

Speciation is implicit: at each species-tree node a branch's genome is inherited by both
children.

<figure markdown="span">
![A species tree with duplication, loss and transfer events on its branches](../img/species_tree_events.svg)
<figcaption>One gene family evolving along a species tree: a duplication (□), a loss (○) and a
transfer (curved arrow) placed on the branches where the Gillespie process fires them.</figcaption>
</figure>

## Rate models

Rates are supplied by a **rate model**. Two ship today; both are subclasses of `RateModel`.

### Uniform rates — every family the same

```python
import zombi2 as z

rates = z.SharedRates(duplication=0.2, transfer=0.1, loss=0.25, origination=0.5)
genomes = z.simulate_genomes(tree, rates, initial_families=40, seed=42)
```

D/T/L are **per gene copy** (the family-level rate scales with copy number); origination is
**per branch**. There is a shorthand that builds `SharedRates` for you:

```python
genomes = z.simulate_genomes(tree, duplication=0.2, transfer=0.1, loss=0.25,
                             origination=0.5, initial_families=40, seed=42)
```

### Per-family sampled rates — ZOMBI-1 style

Give each family its **own** D/T/L, drawn from distributions the first time the family
appears and kept for its lifetime:

```python
rates = z.FamilySampledRates(
    duplication=z.Gamma(2, 0.06),      # built-in distribution
    transfer=z.Exponential(0.08),
    loss=z.Gamma(2, 0.07),
    origination=0.5,                   # per-branch (a single rate)
)
genomes = z.simulate_genomes(tree, rates, initial_families=40, seed=42)
```

Distribution arguments accept:

- a **built-in**: `z.Gamma(shape, scale)`, `z.Exponential(mean)`, `z.LogNormal(mu, sigma)`,
  `z.Uniform(low, high)`, `z.Fixed(value)`;
- any **scipy.stats frozen distribution** (e.g. `scipy.stats.gamma(2, scale=0.1)`);
- a **callable** `rng -> float`.

Negative draws are clipped to 0.

### Genome-wise rates

`SharedRates` is *gene-wise*: the total duplication/transfer/loss rate scales with the
number of gene copies, so a family's size follows an exponential birth–death. `PerGenomeRates`
instead fires each event at a **constant per-genome rate**, independent of genome size (a
target copy is then chosen uniformly):

```python
genomes = z.simulate_genomes(tree, z.PerGenomeRates(duplication=1.0, transfer=0.3,
                                                     loss=0.5, origination=0.4),
                             initial_families=20, seed=1)
```

A useful consequence: family sizes grow *linearly* rather than exponentially, so
genome-wise models are far less prone to runaway growth.

### Branch-wise rates

`BranchRates` makes rates vary **per species-tree branch** by scaling any base rate model
with a per-branch factor (one scalar per branch, scaling duplication/transfer/loss
together; origination is left unscaled). It composes with the base model, so branch and
family heterogeneity combine. Choose one factor source:

```python
base = z.SharedRates(duplication=0.2, transfer=0.1, loss=0.2, origination=0.4)

# 1. autocorrelated (relaxed clock): related lineages have similar rates
z.simulate_genomes(tree, z.BranchRates(base, autocorr_sigma=0.5), seed=1)

# 2. i.i.d. per branch, drawn from a distribution
z.simulate_genomes(tree, z.BranchRates(base, per_branch=z.LogNormal(0.0, 0.5)), seed=1)

# 3. an explicit {branch_name: factor} map (branches not listed keep root_rate)
z.simulate_genomes(tree, z.BranchRates(base, factors={"i3": 10.0}), seed=1)
```

For the relaxed clock, `factor(child) = factor(parent) · exp(N(0, σ·√branch_length))`, so
the drift accumulates with time and `σ = 0` recovers the base model.

## Seeding the root genome

`initial_families` sets how many gene families the root genome starts with (each originated at
time 0). Additional families appear over time at the origination rate.

## The result

`simulate_genomes` returns a `Genomes` object:

```python
genomes.species_tree      # the input tree
genomes.profiles          # ProfileMatrix (families × extant species)
genomes.event_log         # full chronological event log
genomes.gene_families     # {family_id: [EventRecord, ...]}
genomes.gene_trees()      # {family_id: (complete_newick, extant_newick)}
genomes.write("out/")
```

See [Transfers](transfers.md) for transfer mechanics, [Bounding growth](growth.md) for
caps, and [Gene trees & output](gene-trees-and-output.md) for what comes out.
