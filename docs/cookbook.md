# Cookbook

Task-oriented recipes. Every snippet is self-contained and runnable as-is; each imports
the symbols it uses from ZOMBI2's namespaced API.

Where a recipe needs a species tree, it uses this one:

```python
from zombi2.species import simulate_species_tree, BirthDeath

tree = simulate_species_tree(BirthDeath(birth=1.0, death=0.3), n_tips=100, age=5.0, seed=1)
```

For the concepts behind these, see the [user guide](guide/species-trees.md) and the
[FAQ](faq.md); for the `zombi2` command, see the [CLI page](cli.md).

## Species trees

### A pure-birth (Yule) tree

```python
from zombi2.species import simulate_species_tree, Yule

tree = simulate_species_tree(Yule(birth=1.0), n_tips=100, age=5.0, seed=1)
```

`Yule(b)` is exactly `BirthDeath(b, death=0)`.

### A birth–death tree conditioned on the number of tips

```python
from zombi2.species import simulate_species_tree, BirthDeath

tree = simulate_species_tree(BirthDeath(birth=1.0, death=0.3), n_tips=5000, age=5.0, seed=1)
```

You get *exactly* `n_tips` extant species at the given `age` (the tree is a reconstructed,
conditioned birth–death process).

### Interpret the age as stem instead of crown

```python
from zombi2.species import simulate_species_tree, BirthDeath

tree = simulate_species_tree(BirthDeath(1.0, 0.3), n_tips=100, age=5.0,
                             age_type="stem", seed=1)   # default is "crown"
```

### Episodic (skyline) rates with incomplete sampling

Rates are piecewise-constant through time, ordered from the present backward; `shifts` are
the epoch boundaries as ages. `sampling_fraction=ρ` models incomplete extant sampling.

```python
from zombi2.species import simulate_species_tree, EpisodicBirthDeath

model = EpisodicBirthDeath(birth=[1.0, 2.0], death=[0.3, 0.1], shifts=[2.0],
                           sampling_fraction=0.5)
tree = simulate_species_tree(model, n_tips=100, age=5.0, seed=1)
```

### Use your own species tree (from Newick)

```python
from zombi2 import read_newick
from zombi2.genomes import simulate_genomes

with open("my_tree.nwk") as f:
    tree = read_newick(f.read())

genomes = simulate_genomes(tree, duplication=0.2, loss=0.2, origination=0.5, seed=1)
```

Branch lengths are read as durations. (This is what the CLI's `genomes` subcommand does.)

### Add ghost (extinct) lineages

Un-prune the reconstructed tree — graft back the lineages that went extinct (or were not
sampled). Pass the **same model** you built the tree with; ghosts are added in place:

```python
from zombi2.species import simulate_species_tree, BirthDeath, add_ghost_lineages

model = BirthDeath(birth=1.0, death=0.5)
tree = simulate_species_tree(model, n_tips=100, age=5.0, seed=1)

add_ghost_lineages(tree, model, seed=7)   # method="htransform" for a rejection-free sampler

ghosts = [n for n in tree.leaves() if not n.is_extant]   # extinct (e*) tips; extant tips untouched
```

Only birth–death with extinction (or `sampling_fraction < 1`) produces ghosts — a pure-birth
tree has none. See [ghost lineages](guide/ghost-lineages.md).

### A forward (complete) tree with extinct lineages

By default `simulate_species_tree` samples the *reconstructed* tree (extant tips only). Pass
`direction="forward"` to grow the *complete* tree forward in time instead, keeping extinct
lineages. Give **exactly one** of `age` (grow for that long; survivor count is random) or
`n_tips` (condition on that many extant tips):

```python
from zombi2.species import simulate_species_tree, BirthDeath

# grow forward for a fixed age -> complete tree, random number of survivors
tree = simulate_species_tree(BirthDeath(1.0, 0.5), age=5.0, direction="forward", seed=1)

# ...or condition on N extant tips (extinct lineages still included)
tree = simulate_species_tree(BirthDeath(1.0, 0.5), n_tips=50, direction="forward", seed=1)

extinct = [n for n in tree.leaves() if not n.is_extant]
```

This is the native alternative to un-pruning a backward tree
([ghost lineages](guide/ghost-lineages.md)), and it also supports fossilized birth–death
(dated/fossil tips). See [species trees](guide/species-trees.md).

## Gene families: rate models

### The same rates for every family (shared)

The keyword shorthand builds a `SharedRates` for you:

```python
from zombi2.genomes import simulate_genomes

genomes = simulate_genomes(tree, duplication=0.2, transfer=0.1, loss=0.25,
                           origination=0.5, initial_families=40, seed=42)
```

### Per-family rates drawn from distributions (ZOMBI1 style)

Each family samples its own D/T/L once, from the given distributions:

```python
from zombi2 import Gamma, Exponential
from zombi2.genomes import simulate_genomes, FamilySampledRates

genomes = simulate_genomes(tree, FamilySampledRates(
    duplication=Gamma(2, 0.06), transfer=Exponential(0.08),
    loss=Gamma(2, 0.07), origination=0.5), initial_families=40, seed=42)
```

### Genome-wise (per-genome) rates

Each event type fires at a constant per-genome rate (independent of copy number), so family
sizes grow linearly rather than exponentially:

```python
from zombi2.genomes import simulate_genomes, PerGenomeRates

genomes = simulate_genomes(tree, PerGenomeRates(
    duplication=0.5, transfer=0.3, loss=0.4, origination=0.5), seed=42)
```

### Custom distributions

Any distribution argument accepts a built-in, a `scipy.stats` frozen distribution, or a
plain `rng -> float` callable:

```python
import scipy.stats as st
from zombi2 import Gamma
from zombi2.genomes import FamilySampledRates

FamilySampledRates(duplication=Gamma(2, 0.06), loss=0.2, origination=0.5)   # built-in
FamilySampledRates(duplication=st.gamma(2, scale=0.06), loss=0.2, origination=0.5)  # scipy
FamilySampledRates(duplication=lambda rng: rng.gamma(2, 0.06), loss=0.2, origination=0.5)  # callable
```

Built-ins: `Gamma`, `Exponential`, `LogNormal`, `Uniform`, `Fixed` (all from `zombi2`).

## Transfers

Pass a `TransferModel` (from `zombi2.genomes`) to control what a transfer does (see [transfers](guide/transfers.md)).

### Additive vs replacement

```python
from zombi2.genomes import simulate_genomes, TransferModel

# additive (default): recipient gains a copy (net +1)
simulate_genomes(tree, transfer=0.3, loss=0.2, origination=0.5, seed=1)

# replacement: with prob 0.2 the transfer also removes a pre-existing copy (net 0)
simulate_genomes(tree, transfer=0.3, loss=0.2, origination=0.5,
                 transfers=TransferModel(replacement=0.2), seed=1)
```

### Prefer phylogenetically nearby recipients

```python
from zombi2.genomes import simulate_genomes, TransferModel

simulate_genomes(tree, transfer=0.3, loss=0.2, origination=0.5,
                 transfers=TransferModel(distance_decay=2.0), seed=1)
```

`distance_decay=None` (default) picks recipients uniformly; larger values favour close
relatives.

### A transfer/loss-only model (self-transfer as duplication)

With `allow_self=True` a lineage can transfer to itself — mechanically a duplication — so
you can drop explicit duplications entirely. Pair it with a growth cap:

```python
from zombi2.genomes import simulate_genomes, TransferModel

simulate_genomes(tree, transfer=1.0, duplication=0.0, loss=0.3, origination=0.5,
                 transfers=TransferModel(allow_self=True),
                 max_family_size=0.5, seed=1)
```

## Controlling family growth

Both duplication and transfer create copies, so families can grow like `e^{(d−l)t}`. Two
controls (see [bounding growth](guide/growth.md)).

### A hard cap

`max_family_size` is an absolute integer, or a float read as a fraction of the number of
species:

```python
from zombi2.genomes import simulate_genomes

simulate_genomes(tree, duplication=0.5, transfer=0.2, loss=0.1, origination=0.5,
                 max_family_size=0.5, seed=1)   # cap = round(0.5 * N_species)
```

### A soft cap

`carrying_capacity` on the rate model applies logistic (duplication-only) density
dependence; family size settles around `K`:

```python
from zombi2.genomes import simulate_genomes, SharedRates

simulate_genomes(tree, SharedRates(0.5, 0.0, 0.1, 0.5, carrying_capacity=10), seed=1)
```

## Gene order & rearrangements

Put genes on an ordered chromosome and enable inversions/transpositions. The rearrangement
rates live on `SharedRates`; the ordered genome is selected via `genome_factory` (see
[ordered genomes](guide/ordered-genomes.md)):

```python
from zombi2.genomes import simulate_genomes, SharedRates, OrderedGenome

genomes = simulate_genomes(
    tree,
    SharedRates(duplication=0.2, transfer=0.1, loss=0.2, origination=0.5,
                inversion=0.3, transposition=0.3),
    genome_factory=lambda ids: OrderedGenome(ids, extension=0.5),
    seed=1)
```

`extension` is the segment-length knob (probability of extending a rearrangement to the
next gene).

## Getting results out

`simulate_genomes(...)` returns a `Genomes` object:

```python
from zombi2.genomes import simulate_genomes

genomes = simulate_genomes(tree, duplication=0.2, transfer=0.1, loss=0.25,
                           origination=0.5, seed=42)
```

### The presence/copy-number matrix

```python
genomes.profiles.matrix                    # families × extant-species copy numbers (ndarray)
genomes.profiles.families                  # family ids (row order)
print(genomes.profiles.to_tsv())           # copy-number table
print(genomes.profiles.to_tsv(presence=True))   # 0/1 presence table
```

### Reconstructed gene trees

```python
family = genomes.profiles.families[0]
complete, extant = genomes.gene_trees()[family]   # full history, and pruned-to-survivors
```

### The event log

```python
for record in genomes.event_log.records[:10]:
    print(record.time, record.event, record.family)
```

### Write the full ZOMBI1 output

```python
genomes.write("out/")   # species tree, event tables, gene trees, transfers, profiles
```

See [gene trees & output](guide/gene-trees-and-output.md) for the file layout.

## Reproducibility

Pass a `seed` for a reproducible run, or hand in your own NumPy generator with `rng`:

```python
from zombi2.genomes import simulate_genomes

simulate_genomes(tree, duplication=0.2, loss=0.2, origination=0.5, seed=42)

import numpy as np
simulate_genomes(tree, duplication=0.2, loss=0.2, origination=0.5,
                 rng=np.random.default_rng(42))
```

The same seed and inputs reproduce a run within one engine. (Python and Rust engines use
different RNG streams — see the [FAQ](faq.md).)

## Scaling up

### Many replicates in parallel

`run_replicates` runs independent replicates across CPU cores, writing each to
`outdir/replicate_<i>/` and returning a summary per replicate (see
[running in parallel](guide/parallel.md)):

```python
from zombi2.species import BirthDeath
from zombi2.genomes import run_replicates

summaries = run_replicates(
    100, "batch/", BirthDeath(1.0, 0.3),
    n_tips=200, age=5.0,
    duplication=0.2, transfer=0.1, loss=0.25, origination=0.5,
    seed=0)                 # processes=None uses all cores; processes=1 runs serially
```

### Profiles only (the fast path)

The built-in model already runs on Rust automatically. When you only need the
copy-number/presence matrix — not the event log or gene trees — pass `output="profiles"` to
skip the genealogy entirely (much faster; the right path for large datasets):

```python
from zombi2.genomes import simulate_genomes

pm = simulate_genomes(tree, duplication=0.05, transfer=0.03, loss=0.1,
                      origination=0.5, max_family_size=0.3, seed=42,
                      output="profiles")       # -> ProfileMatrix, no event log / gene trees

pm.matrix.shape                                # (n_families, n_extant_species)
```

The default `output="genomes"` returns the full `Genomes`. Both require the compiled
`zombi2_core` extension for the built-in model — see [the Rust engine](guide/rust-engine.md).

## From the command line

The same common cases are available without writing Python — see the
[command-line interface](cli.md):

```bash
zombi2 species --tips 5000 --age 5 -o out/
zombi2 genomes --tree out/species_tree.nwk --dup 0.2 --trans 0.1 --loss 0.25 --orig 0.5 -o out/
```
