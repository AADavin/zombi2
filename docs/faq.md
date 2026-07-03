# FAQ

### Can I simulate gene families on my own species tree?

Yes. Use the `genomes` subcommand with `--tree` pointing at any Newick file, or
`z.read_newick(text)` + `z.simulate_genomes(tree, ...)` in Python. Branch lengths are read
as durations. See the [command-line interface](cli.md).

### Is the CLI using the Rust engine?

Yes — the built-in model runs on Rust automatically, so `zombi2 genomes` / `zombi2 all` use
it (there is no `--fast` flag any more). Add `--profiles-only` for the reduced, counts-only
output. See the [command-line interface](cli.md).

### How do I get just the profile matrix, fast?

Call `simulate_genomes` with `output="profiles"` — it runs the Rust engine over per-family
counts only (no event log or gene trees), which is the right path for large presence/absence
datasets and for ABC:

```python
pm = z.simulate_genomes(tree, duplication=0.05, transfer=0.03, loss=0.1,
                        origination=0.5, output="profiles")   # -> ProfileMatrix
```

The default `output="genomes"` returns the full `Genomes`. There is no longer a separate
`simulate_*_fast` family of functions — Rust is selected automatically. See
[the Rust engine](guide/rust-engine.md).

### Why don't the Rust and Python engines give identical results?

They use different random-number streams, so results are **statistically** equivalent, not
bit-identical. A given `seed` is reproducible **within** one engine.

### Do I need Rust to use ZOMBI2?

For the **built-in** model (`UnorderedGenome` + `UniformRates`, plus the full
`TransferModel`), yes — it runs on Rust with no silent Python fallback, so `simulate_genomes`
raises a clear "build the extension" error if it isn't compiled. **Flexible** models
(`FamilySampledRates`, `GenomeWiseRates`, `BranchRates`, soft `carrying_capacity`, ordered or
nucleotide genomes, a custom `sampler`) run on pure Python automatically and need no Rust. See
[installation](installation.md) and [the Rust engine](guide/rust-engine.md).

### How do I stop a gene family from growing without bound?

Both duplication and transfer add copies, so families can grow exponentially. Use
`max_family_size` for a **hard** cap (an integer, or a float fraction of the number of
species like `0.5`), or `carrying_capacity=K` on the rate model for a **soft** logistic
limit. See [bounding growth](guide/growth.md).

### `Yule` vs `BirthDeath`?

`z.Yule(birth)` is pure birth (no extinction) — the same as `z.BirthDeath(birth, death=0)`.

### Crown age vs stem age?

`--age-type crown` (default) treats `--age` as the age of the tips' most-recent common
ancestor (the root of the reconstructed tree); `--age-type stem` treats it as the stem
age. In the Python API this is the `age_type` argument to `simulate_species_tree`. See
[species trees](guide/species-trees.md).

### How do I give every gene family its own rates (ZOMBI-1 style)?

Pass a `z.FamilySampledRates(...)` rate model with distributions instead of the scalar
shorthand:

```python
z.simulate_genomes(tree, z.FamilySampledRates(
    duplication=z.Gamma(2, 0.06), transfer=z.Exponential(0.08),
    loss=z.Gamma(2, 0.07), origination=0.5), seed=42)
```

See [gene families & rates](guide/gene-families.md).

### How do I run many replicates at once?

`z.run_replicates(...)` distributes independent replicates across CPU cores. See
[running in parallel](guide/parallel.md).

### Where is the full API reference?

In the [API reference](reference/api.md), auto-generated from the docstrings.
