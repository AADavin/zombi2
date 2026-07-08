# Unordered genomes

The **unordered** level is the simplest and fastest of the three genome models (Chapter 7): a genome
is an unordered *set* of gene families, each present in some copy number, with no position, no
neighbours and no length. It is the right model whenever what you care about is gene *content* —
presence/absence profiles, copy number, reconciliations, family sizes.

Once the species tree is fixed, ZOMBI2 populates it with genes: a single forward continuous-time
(Gillespie) process runs over every branch alive at a given moment, firing discrete events that
create, move, and remove gene copies. The result is a set of gene families, each with its own
history, threaded through the species tree.

## The four unordered events

The genome of a lineage is a collection of gene copies. As simulated time advances, four kinds of
event can fire on any living branch — the four events of the unordered level (the ordered and
nucleotide levels add more; see Chapter 7 and their own chapters):

| Event | Effect |
|---|---|
| **Origination (O)** | a brand-new family appears (one copy) on a branch |
| **Duplication (D)** | a gene copy splits into two |
| **Transfer (T)** | a copy is gained by another lineage alive at that time |
| **Loss (L)** | a copy is removed |

A fifth, optional event — **gene conversion (C)**, one copy of a family overwriting another copy of
the same family — is covered in its own section below; it is off unless a `conversion` rate is set.

Speciation is implicit: at each species-tree node the branch's genome is inherited, intact, by
both children. No gene event is attached to the node itself; the four events fire only along the
branches between nodes.

The process is a Gillespie simulation over all branches alive at once. Duplication, transfer, and
loss act **per gene copy**, so a family's total event rate scales with its current copy number and
its size follows an exponential birth–death. Origination acts **per branch**, seeding new families
over time independently of what is already present.

![Three gene families evolving along one species tree, each drawn in its own colour. Every family originates at the root and then accrues its own duplication, transfer and loss, placed on the branch where the shared Gillespie process fired it — one process, run over every branch alive at each moment, threads all three histories through the same tree at once.](figures/species_tree_events_multi.pdf){width=95%}

## Rates

Within the unordered level, `--rate-model` — or a `RateModel` object in Python — chooses **how the
four rates vary across gene families**. Chapter 7 summarises the rate models in a table; this
section works through them. `shared` (the default), `per-genome` and `family` are available on the
command line; explicit per-family and per-branch rate tables are supplied from files with
`--family-rates` / `--branch-rates` (below).

### Shared rates

`SharedRates` gives every family the same per-copy D/T/L rates and a shared per-branch origination
rate:

```python
from zombi2.genomes import (
    SharedRates, PerGenomeRates, FamilySampledRates, TransferModel, simulate_genomes,
)
from zombi2.distributions import Gamma, Exponential, LogNormal, Uniform, Fixed

rates = SharedRates(duplication=0.2, transfer=0.1, loss=0.25, origination=0.5)
genomes = simulate_genomes(tree, rates, initial_families=40, seed=42)
```

Because the family-level rate scales with copy number, this is a gene-wise model: bigger families
experience more events. There is a shorthand that builds `SharedRates` for you from the same
keywords:

```python
genomes = simulate_genomes(tree, duplication=0.2, transfer=0.1, loss=0.25,
                             origination=0.5, initial_families=40, seed=42)
```

### Per-genome rates

`PerGenomeRates` fires each event at a **constant per-genome rate**, independent of genome size; a
target copy is then chosen uniformly among the copies present:

```python
genomes = simulate_genomes(tree, PerGenomeRates(duplication=1.0, transfer=0.3,
                                                     loss=0.5, origination=0.4),
                             initial_families=20, seed=1)
```

Because the rate no longer scales with copy number, family sizes grow *linearly* rather than
exponentially, which makes per-genome models far less prone to runaway growth than gene-wise ones.

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

To instead **fix** particular families' rates by hand, pass a `rates` table — a
`{family_id: (dup, transfer, loss)}` map. Listed families use exactly those rates; families not
listed fall back to the distributions above (so with the default rates of `0`, only the tabulated
families are active):

```python
rates = FamilySampledRates(rates={"1": (3, 2, 1), "2": (4, 0, 1)})   # families 1 and 2, rest inert
```

The same table is read from a TSV on the command line with `--family-rates` (below).

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

Custom rate tables are supplied from TSV files. `--family-rates` gives named families their own
rates (columns `family duplication transfer loss`; unlisted families fall back to `--dup/--trans/
--loss`), and `--branch-rates` gives named branches a transfer **emission** factor and/or
**receptivity** weight (columns `branch emission receptivity`, either optional):

```bash
zombi2 genomes --tree species_tree.nwk \
    --family-rates families.tsv --branch-rates branches.tsv \
    --initial-families 40 --seed 42 -o out/
```

A `--branch-rates` file with only a `receptivity` column stays on the fast Rust engine; one with an
`emission` column (or `--family-rates`) runs on the Python engine.

`simulate_genomes` returns a `Genomes` object exposing the input `species_tree`, the `profiles`
matrix (families $\times$ extant species), the chronological `event_log`, the per-family
`gene_families` records, and `gene_trees()`; `genomes.write("out/")` serialises them.

All the rate models take the same four rates (duplication, transfer, loss per copy or per genome;
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

### Transfer emission and receptivity

A branch can be made more prone to transfer in two independent ways — as a **donor** or as a
**recipient**:

- **Emission** (how often a branch *donates*) is a *rate*, so it is scaled per branch with
  `BranchRates`, restricted to transfer with `events=("transfer",)`:

  ```python
  # branch i3 donates transfers 5x as often (its duplication/loss are untouched)
  rates = BranchRates(SharedRates(transfer=0.2), factors={"i3": 5.0}, events=("transfer",))
  ```

- **Receptivity** (how likely a branch is to *receive*) is a *selection weight* on recipient choice,
  set on the `TransferModel`. Each candidate's weight is multiplied by its receptivity (branches not
  listed default to `1`), composing with `distance_decay`:

  ```python
  # branch i7 receives transfers 10x as readily as an equidistant unlisted branch
  simulate_genomes(tree, transfer=0.3, transfers=TransferModel(receptivity={"i7": 10.0}))
  ```

  A weight of `0` makes a branch never receive; receptivity is honoured by both the Python and the
  Rust engines. On the command line both dials come from one `--branch-rates` file (below).

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

## Gene conversion

**Intra-genome gene conversion** is the intra-genome analogue of a transfer: within a *single*
genome, one copy of a family overwrites ("converts") another copy of the **same family**. It is
**non-reciprocal** (the donor/template copy is unchanged, the recipient is overwritten) and a
**replacement, not a duplication** — the copy number does not change. Repeated, it homogenises a
family's copies (**concerted evolution**) and, because the overwritten copy's own history is erased,
it pulls the reconstructed within-family coalescences toward the present.

![Intra-genome gene conversion. **A** — the mechanism: a duplication makes two copies; a later
conversion has the template copy overwrite the other, so copy number is unchanged but the two
surviving copies now coalesce at the conversion time, not at the (older) duplication. **B** — a
complete gene tree from a run with conversion, with the conversion (filled circle) shown against the
duplications (squares) and losses (crosses).](figures/gene_conversion.pdf){width=100%}

The rate is a per-copy `conversion` rate on `SharedRates` (or the `conversion=` shorthand); a
`ConversionModel` sets *what a conversion does* once it fires — the donor directionality — exactly as
a `TransferModel` does for transfers:

```python
from zombi2 import simulate_genomes, SharedRates, ConversionModel

genomes = simulate_genomes(
    tree,
    SharedRates(duplication=0.4, loss=0.1, conversion=1.0),
    conversions=ConversionModel(bias=0.0),   # 0 = unbiased donor; 1 = toward the founder
    initial_families=40, seed=42,
)
```

A conversion needs both a donor and a recipient, so it fires only on families holding **two or more
copies**: a family with $n$ copies is converted at total rate $c\,n$. On a long lineage that keeps a
family at two copies, conversions arrive as a Poisson process of rate $2c$ and each resets the pair's
coalescence, so the expected within-family coalescence depth is $1/(2c)$ — the oracle the model is
validated against.

### Donor directionality

The **recipient** (overwritten) copy is always chosen uniformly; `ConversionModel(bias=…)` controls
only the **donor**. `bias=0` (the default) draws the donor uniformly among the family's other copies.
`bias=1` always picks the family's oldest lineage (the founder), homogenising the other copies
*toward* it. In between, with probability `bias` the oldest copy donates. Bias is inert when a family
holds exactly two copies (there is only one possible donor).

### From the command line

```bash
zombi2 genomes --tree species_tree.nwk \
    --dup 0.4 --loss 0.1 --conversion 1.0 --conversion-bias 1.0 \
    --initial-families 40 --seed 42 -o out/
```

::: note
Conversion reshapes gene *trees* rather than copy-number profiles, so a run with `--conversion` takes
the full (Python) engine rather than the Rust counts-only fast path, and is supported on unordered
genomes with `--rate-model shared`.
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
SharedRates(duplication=0.5, loss=0.1, origination=0.3, carrying_capacity=20)
```

::: note
`carrying_capacity` shapes duplication smoothly but does **not** bound transfers; use
`max_family_size` when you need a firm ceiling that also accounts for transfer-driven growth. You
can set both.
:::

If a run still diverges — for example `allow_self=True` with no cap — the simulator raises a clear
error pointing you to `max_family_size` rather than hanging.
