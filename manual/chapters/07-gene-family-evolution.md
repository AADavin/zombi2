# Gene-family evolution

ZOMBI2 offers three models of genome evolution, differing in how much structure a genome
carries — from an unordered set of families to a real nucleotide sequence.

- **Independent gene families** — a genome is an unordered *set* of gene copies. Families have
  no position, no neighbours, no length; each evolves on its own. This is the fastest model, and
  the right one when what you care about is gene content and copy number: presence/absence
  profiles, reconciliations, family sizes.
- **Ordered gene families** — genes sit on an *ordered chromosome*. Position now matters:
  rearrangements shuffle neighbourhoods, and a transfer can land beside a syntenic locus. But
  distance is counted in *genes*, not nucleotides — length stays abstract.
- **Nucleotide genomes** — the chromosome is a real sequence of base pairs. Events have real
  lengths, genes and intergenes have real coordinates, and structural events act at nucleotide
  resolution.

![The three models of genome evolution, in increasing structure: an unordered set of gene families; genes on an ordered chromosome (order matters, length does not); and a real nucleotide sequence with genes and intergenes at true coordinates.](figures/genome_models.pdf)

This chapter covers the first and simplest — gene families as an unordered set, evolving by
duplication, transfer, loss and origination. Ordered genomes are the subject of Chapter 10, and
nucleotide genomes of Chapter 11.

Once the species tree is fixed, ZOMBI2 populates it with genes: a single forward continuous-time
(Gillespie) process runs over every branch alive at a given moment, firing discrete events that
create, move, and remove gene copies. The result is a set of gene families, each with its own
history, threaded through the species tree.

## The DTL(O) process

The genome of a lineage is a collection of gene copies. As simulated time advances, four kinds of
event can fire on any living branch:

| Event | Effect |
|---|---|
| **Origination (O)** | a brand-new family appears (one copy) on a branch |
| **Duplication (D)** | a gene copy splits into two |
| **Transfer (T)** | a copy is gained by another lineage alive at that time |
| **Loss (L)** | a copy is removed |

Speciation is implicit: at each species-tree node the branch's genome is inherited, intact, by
both children. No gene event is attached to the node itself; the DTL(O) events fire only along the
branches between nodes.

The process is a Gillespie simulation over all branches alive at once. Duplication, transfer, and
loss act **per gene copy**, so a family's total event rate scales with its current copy number and
its size follows an exponential birth–death. Origination acts **per branch**, seeding new families
over time independently of what is already present.

![One gene family evolving along a species tree: a duplication, a loss, and a transfer, each placed on the branch where the Gillespie process fired it.](figures/species_tree_events.pdf)

## Rate models

Rates are supplied by a *rate model*, a subclass of `RateModel`. The choice of rate model is where
per-copy, per-genome, per-family, and per-branch heterogeneity enters.

### Uniform rates

`UniformRates` gives every family the same per-copy D/T/L rates and a shared per-branch origination
rate:

```python
from zombi2.genomes import (
    UniformRates, GenomeWiseRates, FamilySampledRates, TransferModel, simulate_genomes,
)
from zombi2.distributions import Gamma, Exponential, LogNormal, Uniform, Fixed

rates = UniformRates(duplication=0.2, transfer=0.1, loss=0.25, origination=0.5)
genomes = simulate_genomes(tree, rates, initial_families=40, seed=42)
```

Because the family-level rate scales with copy number, this is a gene-wise model: bigger families
experience more events. There is a shorthand that builds `UniformRates` for you from the same
keywords:

```python
genomes = simulate_genomes(tree, duplication=0.2, transfer=0.1, loss=0.25,
                             origination=0.5, initial_families=40, seed=42)
```

### Genome-wise rates

`GenomeWiseRates` fires each event at a **constant per-genome rate**, independent of genome size; a
target copy is then chosen uniformly among the copies present:

```python
genomes = simulate_genomes(tree, GenomeWiseRates(duplication=1.0, transfer=0.3,
                                                     loss=0.5, origination=0.4),
                             initial_families=20, seed=1)
```

Because the rate no longer scales with copy number, family sizes grow *linearly* rather than
exponentially, which makes genome-wise models far less prone to runaway growth than gene-wise ones.

### Per-family sampled rates

`FamilySampledRates` gives each family its **own** D/T/L rates, drawn from distributions the first
time the family appears and kept fixed for its lifetime:

```python
rates = FamilySampledRates(
    duplication=Gamma(2, 0.06),      # built-in distribution
    transfer=Exponential(0.08),
    loss=Gamma(2, 0.07),
    origination=0.5,                   # per-branch, a single rate
)
genomes = simulate_genomes(tree, rates, initial_families=40, seed=42)
```

Each distribution argument accepts a built-in (`Gamma(shape, scale)`, `Exponential(mean)`,
`LogNormal(mu, sigma)`, `Uniform(low, high)`, `Fixed(value)`), any `scipy.stats` frozen
distribution, or a callable `rng -> float`. Negative draws are clipped to 0.

::: note
`initial_families` sets how many families the root genome starts with, each originated at time 0.
Additional families continue to appear over the run at the origination rate.
:::

### From the command line

The same four rates drive the `genomes` command:

```bash
zombi2 genomes --tree species_tree.nwk \
    --dup 0.2 --trans 0.1 --loss 0.25 --orig 0.5 \
    --initial-families 40 --seed 42 -o out/
```

`simulate_genomes` returns a `Genomes` object exposing the input `species_tree`, the `profiles`
matrix (families $\times$ extant species), the chronological `event_log`, the per-family
`gene_families` records, and `gene_trees()`; `genomes.write("out/")` serialises them.

### Rate models at a glance

| Rate model | Rate applies… | Family size grows | Use when |
|---|---|---|---|
| `UniformRates` | per copy, same for every family | exponentially (gene-wise) | one set of D/T/L/O rates for the whole genome |
| `GenomeWiseRates` | per genome, size-independent | linearly | you want bounded growth that ignores copy number |
| `FamilySampledRates` | per copy, drawn per family | exponentially, family-specific | families should differ — each draws its own D/T/L |

All three take the same four rates (duplication, transfer, loss per copy or per genome;
origination per branch) and compose with the growth caps below.

## Transfers

The rate model decides *how often* a transfer fires; a `TransferModel` decides *what a transfer
does* — who receives the copy, and whether it adds to or overwrites the recipient's genome:

```python
genomes = simulate_genomes(
    tree, transfer=0.3,
    transfers=TransferModel(
        replacement=0.2,      # additive vs replacement
        distance_decay=2.0,   # recipient choice by phylogenetic distance
        allow_self=False,     # self-transfer would act as a duplication
    ),
)
```

The default `TransferModel()` is additive, uniform-recipient, and forbids self-transfer.

### Recipient choice

By default the recipient is drawn **uniformly** among the lineages alive at the transfer time.
Setting `distance_decay=`$\lambda$ favours phylogenetically close recipients: a candidate $r$ is
weighted by $\exp(-\lambda d)$, where $d = 2(t - t_{\mathrm{MRCA}})$ is the patristic distance
between donor and candidate at the transfer time $t$. Larger $\lambda$ means more local transfers;
distant transfers are damped but never forbidden.

### Additive versus replacement

- **Additive** (`replacement=0`): the recipient simply gains a copy, a net $+1$.
- **Replacement** (`replacement=p`): with probability $p$ the transfer also removes one
  pre-existing copy of that family in the recipient — a net-zero swap. Replacement is only possible
  when the recipient already carries the family; otherwise the transfer is additive. `replacement=1`
  makes every eligible transfer a replacement. In the event log a replacement appears as a transfer
  plus a compensating loss.

### Self-transfer

With `allow_self=True` the donor is itself an eligible recipient. A self-transfer creates a second
copy in the same genome, which is mechanically a duplication. This lets you drop explicit
duplications and run a transfer/loss-only model:

```python
simulate_genomes(tree, transfer=1.0, duplication=0.0,
                   transfers=TransferModel(allow_self=True),
                   max_family_size=0.5, seed=1)
```

::: warning
Self-transfers grow families exactly as duplications do, so always pair `allow_self=True` with a
growth cap.
:::

## Bounding growth

A family's copy number is a birth–death process in disguise: with duplication rate above loss rate,
its expected size grows like $e^{(d-l)t}$ without bound. Both duplication and transfer create
copies, so both must be reined in. ZOMBI2 offers a hard cap and a soft cap, which compose.

### Hard cap: `max_family_size`

A single ceiling on family size, enforced across **all** copy-creating events:

```python
simulate_genomes(tree, duplication=0.5, transfer=0.2, loss=0.1, origination=0.3,
                   max_family_size=0.5)     # cap = round(0.5 * N_species)
```

An **integer** is an absolute cap (`max_family_size=20`); a **float** is a fraction or multiple of
the number of species (`0.5` gives half the tips, `2.0` twice the tips). At the cap, duplication is
rate-suppressed and an additive transfer that would overflow is turned into a replacement, so family
size never exceeds the ceiling.

### Soft cap: `carrying_capacity`

A logistic, duplication-only alternative that lives on the rate model: the per-copy duplication rate
is scaled by $\max(0,\, 1 - n/K)$, so family size settles *around* $K$ with a proper stationary
distribution:

```python
UniformRates(duplication=0.5, loss=0.1, origination=0.3, carrying_capacity=20)
```

::: note
`carrying_capacity` shapes duplication smoothly but does **not** bound transfers; use
`max_family_size` when you need a firm ceiling that also accounts for transfer-driven growth. You
can set both.
:::

If a run still diverges — for example `allow_self=True` with no cap — the simulator raises a clear
error pointing you to `max_family_size` rather than hanging.
