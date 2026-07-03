# Command-line interface

Installing the package puts a `zombi2` command on your PATH — a thin wrapper over the library.
Everything starts from a species tree; you then evolve gene families and/or a phenotypic trait
along it — or run the inverse, fitting rates to an empirical profile.

```bash
# 1. a species tree -> out/species_tree.nwk  (runs with defaults)
zombi2 species -o out/

# 2. gene families along that tree (or any Newick tree)
zombi2 genomes --tree out/species_tree.nwk \
    --dup 0.2 --trans 0.1 --loss 0.25 --orig 0.5 --max-family-size 0.5 --seed 42 -o out/

# 3. a phenotypic trait along that tree
zombi2 trait --tree out/species_tree.nwk --model ou --alpha 2 --theta 5 --seed 1 -o out/

# 4. fit gene-family rates to an empirical profile (ABC inference)
zombi2 abc --tree out/species_tree.nwk --profiles empirical_Profiles.tsv \
    --dup 0 1 --loss 0 1.5 --orig 0 4 --n-sims 1000 --seed 1 -o out/
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

**Reproducible runs.** Every run **always** writes the full set of parameters to a log in
`<out>/` — `species_tree.log`, `genomes.log`, or `trait.log` — with the version, timestamp,
command line, seed, and every option used, so any run can be reproduced.

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

## `abc` — fit rates to an empirical profile

`abc` runs the inverse of `genomes`: given a species tree and an **empirical copy-number profile**
(a `families × species` TSV, like the `Profiles.tsv` that `genomes` writes), it fits the D/T/L/O
rates by **Approximate Bayesian Computation** — simulate many profiles from the priors, keep the
draws whose summary statistics land closest to the data, and report the posterior.

Give each rate you want to fit a **prior**: two values `LOW HIGH` (a uniform prior) or one value
(fixed); omit a rate to hold it at 0. At least one rate must be a range (there must be something
to fit).

```bash
zombi2 abc -t out/species_tree.nwk --profiles empirical_Profiles.tsv \
    --dup 0 1 --trans 0 0.5 --loss 0 1.5 --orig 0 4 --n-sims 1000 --seed 1 -o out/
```

It writes **`summary.tsv`** (per-parameter posterior mean / median / 95% CI), **`posterior.tsv`**
(the accepted draws, one column per fitted rate), **`spectra.tsv`** (a posterior-predictive check:
the empirical gene-frequency spectrum against the accepted simulations), and `abc.log`.

**Models.** `--model uniform` (default) fits one shared scalar rate per type on the fast Rust
engine; `--model family` fits each rate as the **mean** of a per-family distribution (ZOMBI-1
style, Python engine — pass `--max-family-size` to bound runaway growth). **ABC-SMC.** `--smc`
switches from rejection to a sequential sampler (`--rounds`, `--particles`, `--quantile`) that
shrinks the tolerance across rounds. **`--regression-adjust`** adds a bias-corrected posterior
(Beaumont 2002) to `summary.tsv`.

Copy-number profiles identify the **gain-side** rates (duplication, origination) well, but **loss**
and **transfer** sit on an identifiability ridge — expect wide intervals there. The `spectra.tsv`
check tells you whether *any* rates reproduce the data.

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
| `--rate-model {uniform,genome-wise,nucleotide}` | `uniform` (default, Rust): same per-copy rates for all families; `genome-wise` (Python): constant per-genome rates, linear growth; `nucleotide`: nucleotide-resolution genomes with variable-length structural events ([see below](#nucleotide-genomes-rate-model-nucleotide)) |
| `--dup` `--trans` `--loss` `--orig` | duplication / transfer / loss / origination rates (per copy; **per nucleotide** for `--rate-model nucleotide`) |
| `--initial-size` | genomes seeded at the root (default: 20 gene families; `1` root chromosome for `nucleotide`) |
| `--max-family-size` | growth cap — integer = absolute, decimal = fraction of N (e.g. `0.5`) [not used by `nucleotide`] |
| `--inversion` `--transposition` | [nucleotide] per-nucleotide inversion / transposition rates |
| `--root-length` `--extension` | [nucleotide] root chromosome length (nt) / geometric event-length parameter (mean `1/(1-extension)`) |
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

### `abc`

| Option | Meaning |
| --- | --- |
| `--tree` / `-t` | species tree the empirical data evolved along (required) |
| `--profiles` | empirical copy-number profile TSV, `families × species` (required) |
| `--dup` `--trans` `--loss` `--orig` | prior per rate: two values `LOW HIGH` (uniform) or one (fixed); omitted → held at 0 |
| `--model {uniform,family}` | `uniform` (default, Rust): shared scalar rates; `family` (Python): fit per-family rate means |
| `--family-shape` | [family] Gamma shape for per-family dispersion (default `2.0`) |
| `--n-sims` / `--accept` | [rejection] number of simulations (default `1000`) / accepted fraction (default `0.05`) |
| `--processes` | [rejection] parallel worker processes (default: serial) |
| `--smc` `--rounds` `--particles` `--quantile` | run ABC-SMC instead of rejection, with these controls |
| `--regression-adjust` | also write the regression-adjusted posterior (Beaumont 2002) |
| `--initial-size` / `--max-family-size` | families seeded per sim / growth cap (advised for `--model family`) |
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

## Nucleotide genomes (`--rate-model nucleotide`)

`--rate-model nucleotide` switches `genomes` to a **nucleotide-resolution** model: each genome is
a circular sequence that evolves by **variable-length structural events** — inversion, deletion,
tandem duplication, transposition, HGT transfer and origination. Genes are not predefined; they
**emerge** as *atoms* — maximal intervals of ancestral sequence with a single shared history — so
you still get a phylogenetic profile and per-atom gene trees, but derived from real sequence
structure rather than an atomic gene-family model.

The shared `--dup` / `--trans` / `--loss` / `--orig` flags become **per-nucleotide** rates here
(so use small values, on the order of `1e-3`); `--inversion` and `--transposition` are the extra
structural events, `--root-length` sets the starting chromosome length, and `--extension` the mean
event length (`1/(1-extension)` nt). `--initial-size` seeds root chromosomes (default `1`).

```bash
zombi2 genomes -t out/species_tree.nwk --rate-model nucleotide \
    --inversion 0.001 --dup 0.0006 --loss 0.0006 --root-length 1000 --seed 1 -o out/
```

`--output profiles` writes the emergent atom profile (`Profiles.tsv` / `Presence.tsv`, atoms ×
species) plus `atoms.tsv` and the per-leaf `Mosaics.tsv`, taking the fast Rust path; the default
`profiles trees` also writes the per-atom `gene_trees/` and their reconciliations
(`Reconciled_complete.nwk` / `Reconciled_extant.nwk` / `Reconciliation_events.tsv`). `--sparse`
applies to the profile as usual. (The family-model `events` / `transfers` / `summary` outputs do
not apply to this model.)

## Scope

The CLI covers the common **uniform-rate** case. For family-sampled or genome-wise rates,
custom transfer mechanics, ordered genomes, or replicate parallelism, use the Python API
(see [gene families & rates](guide/gene-families.md) and the other guides).
