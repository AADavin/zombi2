# Command-line interface

Installing the package puts a `zombi2` command on your PATH. It is a thin wrapper over the
library — it only builds model objects and calls the simulate functions — with three
subcommands that mirror the two-step design.

```bash
# 1. species tree only (backward birth–death) -> out/species_tree.nwk
zombi2 species --birth 1 --death 0.3 --tips 5000 --age 5 --seed 1 -o out/

# 2. gene families along a supplied Newick tree (your own, or one from `species`)
zombi2 genomes --tree out/species_tree.nwk \
    --dup 0.2 --trans 0.1 --loss 0.25 --orig 0.5 --max-family-size 0.5 --seed 42 -o out/

# 3. species tree, then gene families along it, in one run
zombi2 all --birth 1 --death 0.2 --tips 20 --age 5 \
    --dup 0.2 --trans 0.1 --loss 0.25 --orig 0.5 --seed 42 -o out/
```

`species` and `all` build the tree; `genomes` reads one from `--tree` (any Newick file).
So `zombi2 species … -o out/` followed by `zombi2 genomes --tree out/species_tree.nwk …`
is exactly the split form of `zombi2 all`. `species` writes only `species_tree.nwk`;
`genomes` and `all` write the full output (see
[Gene trees & output](guide/gene-trees-and-output.md)).

## Options

| Option | Commands | Meaning |
| --- | --- | --- |
| `--birth` / `--death` | `species`, `all` | speciation / extinction rate (`--death` defaults to 0 = Yule) |
| `--tips` / `--age` | `species`, `all` | number of extant species N / tree age |
| `--age-type {crown,stem}` | `species`, `all` | interpret `--age` as crown (default) or stem age |
| `--tree` / `-t` | `genomes` | input species tree in Newick format |
| `--dup` `--trans` `--loss` `--orig` | `genomes`, `all` | per-copy duplication / transfer / loss / origination rates |
| `--initial-size` | `genomes`, `all` | number of gene families seeded at the root (default 20) |
| `--max-family-size` | `genomes`, `all` | growth cap — integer = absolute, decimal = fraction of N (e.g. `0.5`) |
| `--profiles-only` | `genomes`, `all` | write only the profile matrices (skip event log / gene trees) — see below |
| `--seed` | all | RNG seed for reproducibility |
| `-o` / `--out` | all | output directory |

Run `zombi2 <command> --help` for the authoritative list.

## The Rust engine and `--profiles-only`

The built-in model runs on the native **Rust** engine automatically — there is no flag to
turn it on, and no pure-Python fallback for it (which keeps results reproducible against one
engine). The CLI therefore needs the compiled `zombi2_core` extension; it is **not** built by
`pip install`, so build it once (see [the Rust engine](guide/rust-engine.md)) or the command
exits with a build hint.

By default `genomes`/`all` write the full output — event tables, complete and extant gene
trees, transfers, summary, and the profile matrices. Add `--profiles-only` to skip the event
log and gene trees and write just `species_tree.nwk` + `Profiles.tsv` / `Presence.tsv` — the
fastest path, for when the copy-number/presence matrix is all you need:

```bash
zombi2 genomes --tree out/species_tree.nwk --dup 0.2 --loss 0.25 --orig 0.5 --profiles-only -o out/
```

## Scope

The CLI covers the common **uniform-rate** case. For family-sampled or genome-wise rates,
custom transfer mechanics, ordered genomes, or replicate parallelism, use the Python API
(see [Gene families & rates](guide/gene-families.md) and the other guides).
