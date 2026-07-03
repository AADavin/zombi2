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

**Reproducible runs.** Every run **always** writes the full set of parameters to
`<out>/species_tree.log` (or `<out>/genomes.log`) — version, timestamp, command line, seed, and
every option used — so any run can be reproduced.

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
| `--seed` / `-o` / `--out` | RNG seed / output directory |

### `genomes`

| Option | Meaning |
| --- | --- |
| `--tree` / `-t` | input species tree in Newick format |
| `--rate-model {uniform,genome-wise}` | `uniform` (default, Rust): same per-copy rates for all families; `genome-wise` (Python): constant per-genome rates, linear growth |
| `--dup` `--trans` `--loss` `--orig` | duplication / transfer / loss / origination rates |
| `--initial-size` | number of gene families seeded at the root (default 20) |
| `--max-family-size` | growth cap — integer = absolute, decimal = fraction of N (e.g. `0.5`) |
| `--output {profiles,trees,events,transfers,summary,all}` | which files to write (one or more; default `profiles trees`); `profiles` alone → the fast path — see below |
| `--sparse` | write the profile as a sparse `Profiles_sparse.tsv` instead of the dense matrix (needs `profiles` in `--output`) |
| `--annotate-species` | label internal gene-tree nodes `<gid>\|<species-branch>` (e.g. `g570\|i5`) |
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

## Choosing the output, and the Rust engine

The built-in model runs on the native **Rust** engine automatically — there is no flag to turn
it on, and no pure-Python fallback for it (which keeps results reproducible against one engine).
The `genomes` command therefore needs the compiled `zombi2_core` extension; it is **not** built
by `pip install`, so build it once (see [the Rust engine](guide/rust-engine.md)) or the command
exits with a build hint.

`--output` selects which files to write — any of `profiles`, `trees`, `events`, `transfers`,
`summary`, or `all` (default: `profiles trees`); `species_tree.nwk` is always written, and a
component you don't ask for does no work (e.g. omitting `trees` skips the gene-tree
reconstruction). Asking for **only** `profiles` takes the Rust counts-only fast path (no
genealogy); `--sparse` then writes the profile as a scalable `Profiles_sparse.tsv` long table
instead of the dense matrix:

```bash
zombi2 genomes --tree out/species_tree.nwk --dup 0.2 --loss 0.25 --orig 0.5 --output all -o out/
zombi2 genomes --tree out/species_tree.nwk --dup 0.2 --loss 0.25 --orig 0.5 --output profiles --sparse -o out/
```

## Scope

The CLI covers the common **uniform-rate** case. For family-sampled or genome-wise rates,
custom transfer mechanics, ordered genomes, or replicate parallelism, use the Python API
(see [gene families & rates](guide/gene-families.md) and the other guides).
