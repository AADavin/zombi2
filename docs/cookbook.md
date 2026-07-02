# Cookbook

Task-oriented recipes. Every snippet is self-contained and runnable as-is; each assumes:

```python
import zombi2 as z
```

Where a recipe needs a species tree, it uses this one:

```python
tree = z.simulate_species_tree(z.BirthDeath(birth=1.0, death=0.3), n_tips=100, age=5.0, seed=1)
```

For the concepts behind these, see the [user guide](guide/species-trees.md) and the
[FAQ](faq.md); for the `zombi2` command, see the [CLI page](cli.md).

## Species trees

### A pure-birth (Yule) tree

```python
tree = z.simulate_species_tree(z.Yule(birth=1.0), n_tips=100, age=5.0, seed=1)
```

`z.Yule(b)` is exactly `z.BirthDeath(b, death=0)`.

### A birth–death tree conditioned on the number of tips

```python
tree = z.simulate_species_tree(z.BirthDeath(birth=1.0, death=0.3), n_tips=5000, age=5.0, seed=1)
```

You get *exactly* `n_tips` extant species at the given `age` (the tree is a reconstructed,
conditioned birth–death process).

### Interpret the age as stem instead of crown

```python
tree = z.simulate_species_tree(z.BirthDeath(1.0, 0.3), n_tips=100, age=5.0,
                               age_type="stem", seed=1)   # default is "crown"
```

### Episodic (skyline) rates with incomplete sampling

Rates are piecewise-constant through time, ordered from the present backward; `shifts` are
the epoch boundaries as ages. `sampling_fraction=ρ` models incomplete extant sampling.

```python
model = z.EpisodicBirthDeath(birth=[1.0, 2.0], death=[0.3, 0.1], shifts=[2.0],
                             sampling_fraction=0.5)
tree = z.simulate_species_tree(model, n_tips=100, age=5.0, seed=1)
```

### Use your own species tree (from Newick)

```python
with open("my_tree.nwk") as f:
    tree = z.read_newick(f.read())

genomes = z.simulate_genomes(tree, duplication=0.2, loss=0.2, origination=0.5, seed=1)
```

Branch lengths are read as durations. (This is what the CLI's `genomes` subcommand does.)

### Add ghost (extinct) lineages

Un-prune the reconstructed tree — graft back the lineages that went extinct (or were not
sampled). Pass the **same model** you built the tree with; ghosts are added in place:

```python
model = z.BirthDeath(birth=1.0, death=0.5)
tree = z.simulate_species_tree(model, n_tips=100, age=5.0, seed=1)

z.add_ghost_lineages(tree, model, seed=7)   # method="htransform" for a rejection-free sampler

ghosts = [n for n in tree.leaves() if not n.is_extant]   # ghost_* tips; extant tips untouched
```

Only birth–death with extinction (or `sampling_fraction < 1`) produces ghosts — a pure-birth
tree has none. See [ghost lineages](guide/ghost-lineages.md).

## Gene families: rate models

### The same rates for every family (uniform)

The keyword shorthand builds a `UniformRates` for you:

```python
genomes = z.simulate_genomes(tree, duplication=0.2, transfer=0.1, loss=0.25,
                             origination=0.5, initial_size=40, seed=42)
```

### Per-family rates drawn from distributions (ZOMBI-1 style)

Each family samples its own D/T/L once, from the given distributions:

```python
genomes = z.simulate_genomes(tree, z.FamilySampledRates(
    duplication=z.Gamma(2, 0.06), transfer=z.Exponential(0.08),
    loss=z.Gamma(2, 0.07), origination=0.5), initial_size=40, seed=42)
```

### Genome-wise (per-genome) rates

Each event type fires at a constant per-genome rate (independent of copy number), so family
sizes grow linearly rather than exponentially:

```python
genomes = z.simulate_genomes(tree, z.GenomeWiseRates(
    duplication=0.5, transfer=0.3, loss=0.4, origination=0.5), seed=42)
```

### Custom distributions

Any distribution argument accepts a built-in, a `scipy.stats` frozen distribution, or a
plain `rng -> float` callable:

```python
import scipy.stats as st

z.FamilySampledRates(duplication=z.Gamma(2, 0.06), loss=0.2, origination=0.5)   # built-in
z.FamilySampledRates(duplication=st.gamma(2, scale=0.06), loss=0.2, origination=0.5)  # scipy
z.FamilySampledRates(duplication=lambda rng: rng.gamma(2, 0.06), loss=0.2, origination=0.5)  # callable
```

Built-ins: `z.Gamma`, `z.Exponential`, `z.LogNormal`, `z.Uniform`, `z.Fixed`.

## Transfers

Pass a `z.TransferModel` to control what a transfer does (see [transfers](guide/transfers.md)).

### Additive vs replacement

```python
# additive (default): recipient gains a copy (net +1)
z.simulate_genomes(tree, transfer=0.3, loss=0.2, origination=0.5, seed=1)

# replacement: with prob 0.2 the transfer also removes a pre-existing copy (net 0)
z.simulate_genomes(tree, transfer=0.3, loss=0.2, origination=0.5,
                   transfers=z.TransferModel(replacement=0.2), seed=1)
```

### Prefer phylogenetically nearby recipients

```python
z.simulate_genomes(tree, transfer=0.3, loss=0.2, origination=0.5,
                   transfers=z.TransferModel(distance_decay=2.0), seed=1)
```

`distance_decay=None` (default) picks recipients uniformly; larger values favour close
relatives.

### A transfer/loss-only model (self-transfer as duplication)

With `allow_self=True` a lineage can transfer to itself — mechanically a duplication — so
you can drop explicit duplications entirely. Pair it with a growth cap:

```python
z.simulate_genomes(tree, transfer=1.0, duplication=0.0, loss=0.3, origination=0.5,
                   transfers=z.TransferModel(allow_self=True),
                   max_family_size=0.5, seed=1)
```

## Controlling family growth

Both duplication and transfer create copies, so families can grow like `e^{(d−l)t}`. Two
controls (see [bounding growth](guide/growth.md)).

### A hard cap

`max_family_size` is an absolute integer, or a float read as a fraction of the number of
species:

```python
z.simulate_genomes(tree, duplication=0.5, transfer=0.2, loss=0.1, origination=0.5,
                   max_family_size=0.5, seed=1)   # cap = round(0.5 * N_species)
```

### A soft cap

`carrying_capacity` on the rate model applies logistic (duplication-only) density
dependence; family size settles around `K`:

```python
z.simulate_genomes(tree, z.UniformRates(0.5, 0.0, 0.1, 0.5, carrying_capacity=10), seed=1)
```

## Gene order & rearrangements

Put genes on an ordered chromosome and enable inversions/transpositions. The rearrangement
rates live on `UniformRates`; the ordered genome is selected via `genome_factory` (see
[ordered genomes](guide/ordered-genomes.md)):

```python
genomes = z.simulate_genomes(
    tree,
    z.UniformRates(duplication=0.2, transfer=0.1, loss=0.2, origination=0.5,
                   inversion=0.3, transposition=0.3),
    genome_factory=lambda ids: z.OrderedGenome(ids, extension=0.5),
    seed=1)
```

`extension` is the segment-length knob (probability of extending a rearrangement to the
next gene).

## Getting results out

`z.simulate_genomes(...)` returns a `Genomes` object:

```python
genomes = z.simulate_genomes(tree, duplication=0.2, transfer=0.1, loss=0.25,
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
for record in genomes.event_log[:10]:
    print(record.time, record.event, record.family)
```

### Write the full ZOMBI-1 output

```python
genomes.write("out/")   # species tree, event tables, gene trees, transfers, profiles
```

See [gene trees & output](guide/gene-trees-and-output.md) for the file layout.

## Reproducibility

Pass a `seed` for a reproducible run, or hand in your own NumPy generator with `rng`:

```python
z.simulate_genomes(tree, duplication=0.2, loss=0.2, origination=0.5, seed=42)

import numpy as np
z.simulate_genomes(tree, duplication=0.2, loss=0.2, origination=0.5,
                   rng=np.random.default_rng(42))
```

The same seed and inputs reproduce a run within one engine. (Python and Rust engines use
different RNG streams — see the [FAQ](faq.md).)

## Scaling up

### Many replicates in parallel

`z.run_replicates` runs independent replicates across CPU cores, writing each to
`outdir/replicate_<i>/` and returning a summary per replicate (see
[running in parallel](guide/parallel.md)):

```python
summaries = z.run_replicates(
    100, "batch/", z.BirthDeath(1.0, 0.3),
    n_tips=200, age=5.0,
    duplication=0.2, transfer=0.1, loss=0.25, origination=0.5,
    seed=0)                 # processes=None uses all cores; processes=1 runs serially
```

### The Rust fast paths

Optional native engines for large runs (require the compiled `zombi2_core` extension —
see [the Rust fast path](guide/rust-fast-path.md)). Guard on availability:

```python
if z.rust_available():
    # profiles only (fastest) -> a ProfileMatrix
    profiles = z.simulate_profiles_fast(tree, duplication=0.05, transfer=0.03, loss=0.1,
                                        origination=0.5, max_family_size=0.3, seed=42)

    # a full Genomes (event log + gene trees), just faster — a drop-in for simulate_genomes
    genomes = z.simulate_genomes_fast(tree, duplication=0.05, loss=0.1,
                                      origination=0.5, seed=42)

    # simulate, reconstruct, and write the whole ZOMBI-1 output in Rust
    z.simulate_and_write_fast(tree, "out/", duplication=0.05, loss=0.1,
                              origination=0.5, seed=42)
```

## From the command line

The same common cases are available without writing Python — see the
[command-line interface](cli.md):

```bash
zombi2 species --birth 1 --death 0.3 --tips 5000 --age 5 --seed 1 -o out/
zombi2 genomes --tree out/species_tree.nwk --dup 0.2 --trans 0.1 --loss 0.25 --orig 0.5 -o out/
zombi2 all --birth 1 --death 0.2 --tips 100 --age 5 --dup 0.2 --loss 0.25 --orig 0.5 --seed 42 -o out/
```
