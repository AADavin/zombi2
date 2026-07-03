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

### Fancier models

**Episodic (skyline) rates.** Give several `--birth`/`--death` values plus the shift ages
between them (`K-1` ages for `K` rate epochs, ordered present → past):

```bash
zombi2 species --birth 1 2 --death 0.3 0.1 --shifts 2 --age 5 --tips 40 -o out/
```

**Fossilized birth–death & incomplete sampling** (forward only). `--fossilization` (ψ) adds
dated fossil tips through time, `--sampling-fraction` (ρ) samples only a fraction of extant
species, and `--removal` (`r < 1`) keeps sampled ancestors:

```bash
zombi2 species --model forward --age 6 --fossilization 0.3 -o out/       # fossils
zombi2 species --model forward --age 6 --sampling-fraction 0.5 -o out/   # 50% sampled
```

**Ghost lineages** (backward only). `--ghosts` un-prunes the reconstructed tree, grafting the
extinct/unsampled lineages back on (they get `e*` names, like any extinct tip);
`--ghost-method` picks the
sampler. See [ghost lineages](guide/ghost-lineages.md):

```bash
zombi2 species --tips 50 --death 0.6 --ghosts -o out/
```

**Reproducible runs.** Every run **always** writes its parameters to `<out>/species_tree.log`
(or `<out>/genomes.log`) — version, command line, seed, and the parameters used.
`--log-level {low,medium,high}` sets the detail: `low` is the bare minimum to reproduce,
`medium` (default) adds the core scientific parameters, `high` adds a timestamp and every
argument.

## Options

### `species`

| Option | Meaning |
| --- | --- |
| `--model {backward,forward}` | reconstructed (default) or complete forward tree |
| `--birth` / `--death` | speciation / extinction rate(s) (default `1.0` / `0.3`); several values + `--shifts` = episodic |
| `--shifts` | episodic rate-shift ages (`K-1` for `K` rate values), present → past |
| `--tips` | number of extant species (backward default `50`; forward: `--tips` **or** `--age`) |
| `--age` | tree age / timescale (backward default `1.0`; forward: `--tips` **or** `--age`) |
| `--age-type {crown,stem}` | interpret `--age` as crown (default) or stem age [backward] |
| `--sampling-fraction` | [forward] fraction of extant species sampled, ρ (default `1.0`) |
| `--fossilization` | [forward] fossil sampling rate ψ — fossilized birth–death (default `0`) |
| `--removal` | [forward] removal probability `r` on sampling (`r<1` → sampled ancestors; default `1.0`) |
| `--ghosts` / `--ghost-method {rejection,htransform}` | [backward] un-prune extinct/unsampled ghost lineages |
| `--max-attempts` / `--max-lineages` | [forward] extinction-retry cap / live-lineage cap |
| `--log-level {low,medium,high}` | detail of the always-written `<out>/species_tree.log` (default medium) |
| `--seed` / `-o` / `--out` | RNG seed / output directory |

### `genomes`

| Option | Meaning |
| --- | --- |
| `--tree` / `-t` | input species tree in Newick format |
| `--rate-model {uniform,genome-wise}` | `uniform` (default, Rust): same per-copy rates for all families; `genome-wise` (Python): constant per-genome rates, linear growth |
| `--dup` `--trans` `--loss` `--orig` | duplication / transfer / loss / origination rates |
| `--initial-size` | number of gene families seeded at the root (default 20) |
| `--max-family-size` | growth cap — integer = absolute, decimal = fraction of N (e.g. `0.5`) |
| `--profiles-only` | write only the profile matrices (skip event log / gene trees) — see below |
| `--log-level {low,medium,high}` | detail of the always-written `<out>/genomes.log` (default medium) |
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
