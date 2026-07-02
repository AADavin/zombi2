# FAQ

### Can I simulate gene families on my own species tree?

Yes. Use the `genomes` subcommand with `--tree` pointing at any Newick file, or
`z.read_newick(text)` + `z.simulate_genomes(tree, ...)` in Python. Branch lengths are read
as durations. See the [command-line interface](cli.md).

### Is the CLI using the Rust engine?

Only when you pass `--fast`. By default `zombi2 genomes` and `zombi2 all` run the
pure-Python simulator. With `--fast` the output is the same; add `--profiles-only` for the
reduced (profiles-only) output. See the [command-line interface](cli.md).

### What are the Rust fast paths, and when should I use each?

All three cover the built-in `UnorderedGenome` + `UniformRates` model and require the
compiled `zombi2_core` extension:

| Function | Returns | Use when |
| --- | --- | --- |
| `z.simulate_profiles_fast(tree, ...)` | a `ProfileMatrix` (counts only) | you need just the presence/copy-number matrix at scale |
| `z.simulate_genomes_fast(tree, ...)` | a full `Genomes` (with `.event_log`, `.gene_trees()`, `.write()`) | you want the complete result, just faster — a drop-in for `simulate_genomes` |
| `z.simulate_and_write_fast(tree, "out/", ...)` | writes the whole output in Rust | you're generating large datasets straight to disk |

The pure-Python `z.simulate_genomes` stays the default. See the
[Rust fast path](guide/rust-fast-path.md) guide.

### Why don't the Rust and Python engines give identical results?

They use different random-number streams, so results are **statistically** equivalent, not
bit-identical. A given `seed` is reproducible **within** one engine.

### Do I need Rust to use ZOMBI2?

No. Rust is entirely optional acceleration. Everything works on pure Python; the Rust
functions and `--fast` simply become unavailable (with a clear message) if the extension
isn't built. See [installation](installation.md).

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
