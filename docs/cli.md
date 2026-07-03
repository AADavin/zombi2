# Command-line interface

Installing the package puts a `zombi2` command on your PATH — a thin wrapper over the library.
Everything starts from a species tree; you then evolve gene families and/or a phenotypic trait
along it.

```bash
# 1. a species tree -> out/species_tree.nwk  (runs with defaults)
zombi2 species -o out/

# 2. gene families along that tree (or any Newick tree)
zombi2 genomes --tree out/species_tree.nwk \
    --dup 0.2 --trans 0.1 --loss 0.25 --orig 0.5 --max-family-size 0.5 --seed 42 -o out/

# 3. a phenotypic trait along that tree
zombi2 trait --tree out/species_tree.nwk --model ou --alpha 2 --theta 5 --seed 1 -o out/
```

`species` writes only `species_tree.nwk`; `genomes` reads a tree from `--tree` and writes the
full output (see [gene trees & output](guide/gene-trees-and-output.md)); `trait` reads a tree and
writes the tip and ancestral trait values (see [trait evolution](guide/traits.md)). Run
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

## `trait` — a phenotypic trait

`trait` evolves one trait **along a species tree you provide** (`--tree`, like `genomes`) and
writes the result to `-o`. It always outputs both the tip values and the **ancestral** node
values. Pick a model with `--model`:

```bash
T=out/species_tree.nwk

zombi2 trait -t $T --model bm --sigma2 0.5 -o out/                     # Brownian motion
zombi2 trait -t $T --model ou --alpha 2 --theta 5 -o out/             # pulled toward an optimum
zombi2 trait -t $T --model eb --sigma2 1 --rate -1.5 -o out/          # early burst (rate < 0)
zombi2 trait -t $T --model mk --states 3 --rate 0.5 -o out/           # discrete 3-state
zombi2 trait -t $T --model threshold --thresholds 0 -o out/          # binary from a liability
zombi2 trait -t $T --model dec --areas A,B,C --dispersal 0.3 -o out/ # geographic ranges
```

It writes two files: **`traits.tsv`** (a `node`/`trait` table over *every* node — tips named
`n*`, ancestral nodes `i*`/`root`) and **`trait_tree.nwk`** (the tree with each value annotated
as `[&trait=…]`).

**Replicates.** `--replicates N` simulates the trait `N` times with the same parameters and
writes a **wide** `traits.tsv` — one column per replicate (`rep_1 … rep_N`), one row per node —
instead of the annotated tree. Handy for an empirical distribution or method-testing dataset:

```bash
zombi2 trait -t $T --model bm --sigma2 0.5 --replicates 100 --seed 1 -o out/
```

**The Mk chain.** For `--model mk` the transition structure is up to you: the default is
equal-rates (every state ↔ every state), `--ordered` restricts to adjacent-only steps
(`i ↔ i±1`, a meristic character), and `--q-matrix FILE` reads an arbitrary rate matrix — a
whitespace/comma-separated `k×k` grid (rows = *from*-state, columns = *to*-state; the diagonal is
ignored, blank/`#` lines skipped), which overrides `--states`/`--rate`/`--ordered`.

**DEC (geographic ranges).** `--model dec` evolves a range over discrete areas by dispersal /
extinction along branches plus cladogenetic range splits at speciations. Set the areas with
`--areas` (a count like `3`, or labels like `A,B,C`), the `--dispersal` / `--extinction` rates,
optionally cap the range with `--max-range-size`, and pin the root range with `--root-range`
(e.g. `A`). Ranges are written as `{A,B}`.

## Options

### `species`

| Option | Meaning |
| --- | --- |
| `--model {backward,forward}` | reconstructed (default) or complete forward tree |
| `--birth` / `--death` | speciation / extinction rate (default `1.0` / `0.3`) |
| `--tips` | number of extant species (backward default `50`; forward: `--tips` **or** `--age`) |
| `--age` | tree age / timescale (backward default `1.0`; forward: `--tips` **or** `--age`) |
| `--age-type {crown,stem}` | interpret `--age` as crown (default) or stem age [backward] |
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

### `trait`

| Option | Meaning |
| --- | --- |
| `--tree` / `-t` | input species tree in Newick format (required) |
| `--model {bm,ou,eb,mk,threshold,dec}` | trait model (default `bm`) |
| `--sigma2` | diffusion rate [bm/ou/eb/threshold] |
| `--x0` / `--trend` | root value / directional drift [bm/eb/threshold; OU root defaults to `--theta`] |
| `--alpha` / `--theta` | OU mean-reversion strength / optimum [ou] |
| `--rate` | EB rate-of-change (negative = early burst) [eb], or the per-transition rate [mk] |
| `--states` / `--ordered` / `--q-matrix` | mk: number of states / adjacent-only chain / arbitrary Q from a file |
| `--thresholds` | comma-separated liability cut points [threshold] |
| `--areas` `--dispersal` `--extinction` `--max-range-size` `--root-range` | DEC range-evolution parameters |
| `--replicates` | simulate this many times → wide one-column-per-replicate table |
| `--seed` / `-o` / `--out` | RNG seed / output directory (required) |

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
