# Command-line interface

Installing the package puts a `zombi2` command on your PATH — a thin wrapper over the library
with two subcommands that mirror the two-step design: build a species tree, then evolve gene
families along it.

```bash
# 1. a species tree -> out/species_tree.nwk  (runs with defaults)
zombi2 species -o out/

# 2. gene families along that tree (or any Newick tree)
zombi2 genomes --tree out/species_tree.nwk \
    --dup 0.2 --trans 0.1 --loss 0.25 --orig 0.5 --max-family-size 0.5 --seed 42 -o out/
```

`species` writes only `species_tree.nwk`; `genomes` reads a tree from `--tree` and writes the
full output (see [gene trees & output](guide/gene-trees-and-output.md)). Run
`zombi2 <command> -h` for a command's options.

## `species` — the species tree

`species` runs a birth–death process in one of two models, chosen by `--model`:

* **`backward`** (default) — the *reconstructed* tree: only surviving (extant) lineages,
  **conditioned on the number of tips**. `--tips` sets how many extant species you get; `--age`
  sets the **timescale** — the crown or stem age, in the same time units as the DTL rates.
  They are orthogonal, not alternatives, so you give both (both have defaults).
* **`forward`** — the *complete* tree grown forward in time, keeping extinct lineages. Here you
  give **exactly one** stopping condition: `--age` (grow for that long; the number of survivors
  is random) **or** `--tips` (grow until that many extant species; the age is random).

Everything has a default, so a bare run works:

```bash
zombi2 species -o out/                                            # backward, 50 tips, age 1.0
zombi2 species --tips 5000 --age 5 --birth 1 --death 0.3 -o out/  # backward, custom
zombi2 species --model forward --age 5 -o out/                   # complete tree, grown 5 time units
zombi2 species --model forward --tips 50 -o out/                 # complete tree, until 50 extant
```

## Options

### `species`

| Option | Meaning |
| --- | --- |
| `--model {backward,forward}` | reconstructed (default) or complete forward tree |
| `--birth` / `--death` | speciation / extinction rate (default `1.0` / `0.3`) |
| `--tips` | number of extant species (backward default `50`; forward: `--tips` **or** `--age`) |
| `--age` | tree age / timescale (backward default `1.0`; forward: `--tips` **or** `--age`) |
| `--age-type {crown,stem}` | interpret `--age` as crown (default) or stem age [backward] |
| `--max-attempts` | [forward] retries before giving up when the process goes extinct (default 10000) |
| `--max-lineages` | [forward] abort a run exceeding this many live lineages (default 1000000) |
| `--seed` / `-o` / `--out` | RNG seed / output directory |

### `genomes`

| Option | Meaning |
| --- | --- |
| `--tree` / `-t` | input species tree in Newick format |
| `--dup` `--trans` `--loss` `--orig` | per-copy duplication / transfer / loss / origination rates |
| `--initial-size` | number of gene families seeded at the root (default 20) |
| `--max-family-size` | growth cap — integer = absolute, decimal = fraction of N (e.g. `0.5`) |
| `--profiles-only` | write only the profile matrices (skip event log / gene trees) — see below |
| `--seed` / `-o` / `--out` | RNG seed / output directory |

Run `zombi2 <command> -h` for the authoritative list.

## The Rust engine and `--profiles-only`

The built-in model runs on the native **Rust** engine automatically — there is no flag to turn
it on, and no pure-Python fallback for it (which keeps results reproducible against one engine).
The `genomes` command therefore needs the compiled `zombi2_core` extension; it is **not** built
by `pip install`, so build it once (see [the Rust engine](guide/rust-engine.md)) or the command
exits with a build hint.

By default `genomes` writes the full output — event tables, complete and extant gene trees,
transfers, summary, and the profile matrices. Add `--profiles-only` to skip the event log and
gene trees and write just `species_tree.nwk` + `Profiles.tsv` / `Presence.tsv` — the fastest
path, for when the copy-number/presence matrix is all you need:

```bash
zombi2 genomes --tree out/species_tree.nwk --dup 0.2 --loss 0.25 --orig 0.5 --profiles-only -o out/
```

## Scope

The CLI covers the common **uniform-rate** case. For family-sampled or genome-wise rates,
custom transfer mechanics, ordered genomes, or replicate parallelism, use the Python API
(see [gene families & rates](guide/gene-families.md) and the other guides).
