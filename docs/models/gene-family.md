# Gene-family evolution (DTL)

Gene families evolve **forward** along a fixed species tree by a **duplication–transfer–loss (DTL)**
process, with new families entering by **origination**. A family is born as a single copy on a branch,
copied by **duplication** (a copy splits in two) or **transfer** (a copy is gained by another lineage
alive at that time), and removed by **loss** — one continuous-time (Gillespie) process over all
co-existing branches. Speciation is implicit: at each species-tree node a branch's genome is inherited
by both children. The result is a **gene tree per family** (reconciled against the species tree) and a
**presence/absence profile** — the copy number of every family in every extant species. The models below
are variants of this process that differ only in **how the D/T/L rates are set** — shared across every
family, drawn per family, held constant per genome, or scaled per branch.

| Model | Rates | Reach for it when |
| --- | --- | --- |
| **SharedRates** | one per-copy D/T/L for every family (Rust engine) | the default DTL backbone, fast and uniform |
| **PerGenomeRates** | constant per-genome totals; families grow linearly | you want size-independent rates and no runaway growth |
| **FamilySampledRates** | each family draws its own D/T/L (ZOMBI-1 style) | families should differ in their evolutionary rates |
| **BranchRates** | a per-branch factor scaling any base model | rates vary across the species tree (relaxed clock) |

D/T/L are **per gene copy** unless a model says otherwise, so a family's copy number follows an
exponential birth–death and the genome-wide rate scales with the number of copies; **origination** is
always **per branch** (independent of genome size). All rates are per unit of tree time — see
[Conventions § Rates and time](../contributing/conventions.md#rates-and-time). See
[Transfers](../guide/transfers.md) for how a transfer recipient is chosen.

## The models

### SharedRates

Every gene family shares the same per-copy `duplication`, `transfer`, and `loss` rates, plus a per-branch
`origination` rate. The default DTL model and the one run by the compiled **Rust engine** on the
counts-only path. Optional `carrying_capacity` (K) adds logistic density dependence, damping duplication
as a family grows; the simulator's `max_family_size` (CLI `--max-family-size`) is a hard ceiling. This is the model built for you when you pass bare
`--dup/--trans/--loss/--orig` (or the `simulate_genomes` rate shorthand).

### PerGenomeRates

Each event type fires at a **constant per-genome rate**, independent of how many copies the genome holds —
in contrast to `SharedRates`, whose totals scale with genome size. A useful consequence: family sizes grow
**linearly** rather than exponentially, so these models are far less prone to runaway growth. When an event
fires the target copy is chosen uniformly; origination is per branch. Selected on the CLI with
`--rate-model per-genome` (Python engine).

### FamilySampledRates

Each gene family draws its **own** `duplication`/`transfer`/`loss` rates from distributions the first time
it is seen, then keeps them for life (the ZOMBI-1 style of rate heterogeneity). Pass a float for a fixed
rate or a distribution (e.g. `Exponential`, `Gamma`) to sample per family. Origination stays a single
per-branch rate. Python API only. Also honours an optional `carrying_capacity`.

### BranchRates

Wraps **any** base rate model and multiplies its D/T/L weights on each species-tree branch by a per-branch
factor (origination is left unscaled), so branch heterogeneity composes with whatever family/uniform rates
the base uses. Provide exactly one factor source: `autocorr_sigma` (a **relaxed clock** — the factor drifts
lognormally down the tree, so relatives have similar rates), `per_branch` (a distribution drawn i.i.d. per
branch), or an explicit `factors` `{branch_name: factor}` map. Python API only.

## Command line

Bare `--dup/--trans/--loss/--orig` selects `SharedRates` (the default `--rate-model shared`, Rust engine);
`--rate-model per-genome` selects `PerGenomeRates` (Python). `FamilySampledRates` and `BranchRates` are
Python-API only. `--write` selects the output parts (default `profiles trees`); `species_tree.nwk` is
always copied through, and the run writes a `genomes.log` manifest.

```bash
# a species tree to evolve genomes along
zombi2 species --birth 1 --death 0.3 --tips 20 --age 5 --seed 1 -o run/

# SharedRates DTL (default), full output
zombi2 genomes -t run/species_tree.nwk --dup 0.2 --trans 0.1 --loss 0.25 --orig 0.5 \
  --initial-families 40 --write all --seed 42 -o run/

# PerGenomeRates: constant per-genome totals, linear growth
zombi2 genomes -t run/species_tree.nwk --rate-model per-genome \
  --dup 0.5 --trans 0.2 --loss 0.5 --orig 0.4 --write profiles --seed 42 -o run/
```

## Python

Models live in `zombi2.genomes` (and re-export at the top level, so `zombi2.SharedRates` also works):

```python
from zombi2.species import BirthDeath, simulate_species_tree
from zombi2.genomes import (SharedRates, PerGenomeRates, FamilySampledRates,
                            BranchRates, simulate_genomes)
from zombi2.distributions import Exponential

tree = simulate_species_tree(BirthDeath(1.0, 0.3), n_tips=20, age=5.0, seed=1)

# SharedRates (default): one per-copy D/T/L for every family
rates = SharedRates(duplication=0.2, transfer=0.1, loss=0.25, origination=0.5)
genomes = simulate_genomes(tree, rates, initial_families=40, seed=42)

# PerGenomeRates: constant per-genome totals, linear growth
genomes = simulate_genomes(
    tree, PerGenomeRates(duplication=0.5, transfer=0.2, loss=0.5, origination=0.4),
    initial_families=40, seed=42)

# FamilySampledRates: each family draws its own D/T/L
fs = FamilySampledRates(duplication=Exponential(0.2), transfer=Exponential(0.1),
                        loss=Exponential(0.25), origination=0.5)
genomes = simulate_genomes(tree, fs, initial_families=40, seed=42)

# BranchRates: a per-branch factor (relaxed clock) scaling any base model
br = BranchRates(SharedRates(0.2, 0.1, 0.25, 0.5), autocorr_sigma=0.5)
genomes = simulate_genomes(tree, br, initial_families=40, seed=42)
```

There is a shorthand that builds `SharedRates` for you (pass the rate model **or** the shorthand, not
both):

```python
genomes = simulate_genomes(tree, duplication=0.2, transfer=0.1, loss=0.25,
                           origination=0.5, initial_families=40, seed=42)
```

## Output

The returned `Genomes` object (and `--write`) exposes:

- **Profile matrix** — `genomes.profiles`; `Profiles.tsv` (copy counts, families × extant species) and
  `Presence.tsv` (its 0/1 binarization). `--sparse` writes `Profiles_sparse.tsv` instead.
- **Gene trees** — `genomes.gene_trees()`; `gene_trees/<family>_complete.nwk` (all lineages) and
  `_extant.nwk` (survivors only), each **reconciled** with the species tree.
- **Reconciliations & per-family events** — `gene_family_events/<family>_events.tsv` records where each
  family's D/T/L/S events map onto the species tree.
- **Events trace** — `Events_trace.tsv`, one compact chronological log of every event, from which gene
  trees can be reconstructed on demand. Also `Transfers.tsv` and `Gene_family_summary.tsv`.

Node and family names follow the [standard naming](../contributing/conventions.md#naming) (`g*` gene
lineages, plain integers for families); event codes are **O**rigination, **D**uplication, **T**ransfer,
**L**oss, **S**peciation — see [Conventions § Outputs](../contributing/conventions.md#outputs).

## Validation

- **SharedRates.** With duplication and loss only, the mean copy number at a leaf matches
  `exp((duplication − loss) · age)` over many replicates
  (`test_genome_dtl.py::test_dl_mean_copy_number`); the compiled Rust engine and the pure-Python
  reference engine agree on the mean family count within Monte-Carlo error over a shared model
  (`test_rust.py::test_rust_matches_python_engine`).
- **PerGenomeRates.** Under a size-independent per-genome model the realized duplication/loss event
  count equals `rate × total_branch_length` — a Poisson oracle checked to Monte-Carlo error for
  duplication-only and loss-only runs, with tripling the rate scaling the mean count by the same
  factor (`test_extensibility.py::test_per_genome_rates_event_counts_match_poisson_oracle`).
- **FamilySampledRates.** On an ultrametric tree each family's realized presence across the extant
  leaves matches its own closed-form survival probability `exp(−loss × T)`: binning families by their
  sampled loss rate, the observed per-leaf presence fraction tracks the per-bin oracle
  `mean(exp(−loss × T))`, confirming the simulator honours each family's cached rate
  (`test_extensibility.py::test_family_sampled_loss_calibrates_to_per_family_rate`).
- **BranchRates.** A branch scaled by a factor `f` yields about `f` times the events of the unscaled
  branch — in the low-rate (linear) regime each factor's mean loss count matches its own closed-form
  pure-death expectation `M·(1 − exp(−loss·f·age))`, and the event-count ratio calibrates to `f`
  rather than merely being "more" (`test_branch_rates.py::test_branch_factor_scales_event_count_proportionally`).

## References

- Davín, A. A., Tricou, T., Tannier, E., de Vienne, D. M. & Szöllősi, G. J. (2020). Zombi: a phylogenetic
  simulator of trees, genomes and sequences that accounts for dead lineages. *Bioinformatics* 36(4):
  1286–1288.
