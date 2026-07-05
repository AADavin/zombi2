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

# 5. gene families whose loss/gain is conditioned on a trait
zombi2 coevolve --couple traits:genes --tree out/species_tree.nwk \
    --trait-model mk --states 2 --trait-center --responsive 0.3 --effect-loss 3 --seed 1 -o out/
```

`species` writes only `species_tree.nwk`; `genomes` reads a tree from `--tree` and writes the
full output (see [gene trees & output](guide/gene-trees-and-output.md)); `trait` reads a tree and
writes the tip and ancestral trait values (see [trait evolution](guide/traits.md)).

Run `zombi2 <command> -h` for a command's options — each command's help groups its flags into
labelled model sections (the general options, then one section per model), so the parameters of
different models are kept clearly apart. `zombi2 --version` prints the version.

## `species` — the species tree

`species` runs a birth–death process in one of two modes, chosen by `--mode`:

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
zombi2 species --mode forward --age 5 -o out/                   # complete tree, grown 5 time units
zombi2 species --mode forward --tips 50 -o out/                 # complete tree, until 50 extant
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
zombi2 species --mode forward --age 6 --fossilization 0.3 -o out/       # fossils
zombi2 species --mode forward --age 6 --sampling-fraction 0.5 -o out/   # 50% sampled
```

**Mass extinctions** (forward only). `--mass-extinction AGE FRACTION` fires an instantaneous,
tree-wide pulse: at `AGE` before the present, every lineage then alive dies with probability
`FRACTION`. Repeat the flag for several pulses. Unlike raising `--death` over an episodic epoch
(which spreads extra extinction across a window), a pulse wipes out a fraction of the standing
diversity in a single instant. It needs `--age` (the pulse time is an age before a fixed present):

```bash
zombi2 species --mode forward --age 5 --mass-extinction 2.5 0.75 -o out/        # 75% die at age 2.5
zombi2 species --mode forward --age 5 \
    --mass-extinction 1.0 0.5 --mass-extinction 2.5 0.75 -o out/                 # two pulses
```

The pulse's victims become extinct `e*` leaves, so a downstream `zombi2 genomes` run sees the
mass extinction's genomic aftermath (families lost with the dead clades, transfers from the dead).

**Per-lineage & diversity-dependent rates** (forward only). `--diversification` picks a
heterogeneous-rate process instead of constant-rate birth–death. `clads` gives every lineage its
own speciation rate that shifts lognormally at each split (`--clads-alpha`, `--clads-sigma`,
`--turnover`); `diversity-dependent` makes speciation decline toward a carrying capacity
`--carrying-capacity/-K`:

```bash
# ClaDS: per-lineage rates, α<1 slows speciation toward the present, ε=μ/λ turnover
zombi2 species --mode forward --diversification clads \
    --birth 1.0 --clads-alpha 0.9 --clads-sigma 0.2 --turnover 0.1 --age 5 -o out/

# diversity-dependent: radiates then saturates near K=50
zombi2 species --mode forward --diversification diversity-dependent \
    --birth 2 --death 0.2 -K 50 --age 15 -o out/
```

Both take `--age` or `--tips` (diversity-dependent needs `--tips ≤ K`) and compose with
`--mass-extinction` and `--sampling-fraction`.

**Clade-specific rate shifts** (forward only). `--clade-shift AGE BIRTH DEATH` schedules a
diversification shift: at `AGE` before the present, one lineage then alive (chosen at random,
since contemporaneous lineages are exchangeable) and *all of its descendants* switch to speciation
`BIRTH` / extinction `DEATH`. Repeat for several radiating or collapsing clades — the discrete,
hand-specified version of clade rate heterogeneity. Needs `--age`:

```bash
# one clade starts diversifying fast at age 3, on a slow background
zombi2 species --mode forward --birth 0.6 --death 0.4 --age 5 \
    --clade-shift 3.0 2.0 0.1 -o out/
# several clades shifting to different regimes
zombi2 species --mode forward --age 6 \
    --clade-shift 4.0 2.0 0.1 --clade-shift 2.0 0.3 0.5 -o out/
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

**Models.** `--rate-model uniform` (default) fits one shared scalar rate per type on the fast Rust
engine; `--rate-model family` fits each rate as the **mean** of a per-family distribution (ZOMBI-1
style, Python engine — pass `--max-family-size` to bound runaway growth). **ABC-SMC.** `--smc`
switches from rejection to a sequential sampler (`--rounds`, `--particles`, `--quantile`) that
shrinks the tolerance across rounds. **`--regression-adjust`** adds a bias-corrected posterior
(Beaumont 2002) to `summary.tsv`.

Copy-number profiles identify the **gain-side** rates (duplication, origination) well, but **loss**
and **transfer** sit on an identifiability ridge — expect wide intervals there. The `spectra.tsv`
check tells you whether *any* rates reproduce the data.

## `sequence` — substitution branch lengths

`genomes` produces gene trees in **time** units. `sequence` rescales them into
**substitutions/site** under a gene × lineage clock, as a separate step — so you can retune the
rate model without re-simulating gene content. It reads a prior `genomes` run's `species_tree.nwk`
and `Events_trace.tsv` (so run `genomes` with `trace` in `--write`), replays the genealogy, and
writes the phylograms:

```bash
zombi2 genomes  -t out/species_tree.nwk --dup 0.2 --trans 0.1 --loss 0.2 --orig 0.5 \
    --write trace profiles -o run/
zombi2 sequence --genomes run/ --branch-speed 0.4 --family-speed 0.5 -o run/
```

The model is `rate(family g, species branch b) = R_b · s_g`: a **shared lineage clock** `R_b`
(either `--branch-speed` autocorrelated lognormal, or `--branch-bins` for the discrete-bin GTDB
model) times a **per-family speed** `s_g ~ LogNormal(0, --family-speed)`. It writes
`gene_trees/<family>_extant_subst.nwk` (+ `_complete_subst.nwk`) and, for reproducibility,
`gene_family_speeds.tsv` / `branch_rates.tsv`. See
[Rate variation](guide/rate-variation.md#family-sequence-evolution) for the model.

## `coevolve --couple traits:genes` — trait-conditioned gene families

The **`traits:genes`** edge of [`coevolve`](coevolution_models.md) links the two halves of
the toolkit: it evolves a phenotypic trait along the tree, then evolves a **panel** of gene
families whose loss and gain **depend on the local trait value**, so the resulting profile carries
a known, trait-linked signal (the forward generator behind reading gene content as a record of a
trait's history — e.g. dating the tree from the Great Oxidation Event). It simulates the trait
with any [`trait`](#trait-a-phenotypic-trait) model (`--trait-model`), builds the coupling, and
writes the gene-family output alongside the trait and a coupling manifest. (This was the standalone
`coevolve-genetrait` command before it was folded into `coevolve`.)

```bash
T=out/species_tree.nwk

# a binary aerobic(1)/anaerobic(0) trait; 30% of a 40-family panel respond to it
zombi2 coevolve --couple traits:genes -t $T \
    --trait-model mk --states 2 --rate 0.3 --trait-center \
    --panel 40 --responsive 0.3 --weight 1 --effect-loss 3 \
    --loss 0.4 --trans 1.0 --write all --seed 7 -o out/
```

A *responsive* family is retained where the trait favours it (loss scaled by
`exp(-effect_loss · weight · trait)`) and purged where it does not; gain is field-blind
horizontal transfer, so the **net** gene content of a lineage tracks its trait. `--responsive`
chooses which families respond — a count (`8`), a fraction (`0.3`), an id/index list
(`F3,F7,12`), or `@file` of ids — and `--signed` randomises the weight signs so some families
co-occur with a high trait value and others with a low one. `--trait-center` centers a discrete
trait's states (recommended for a binary character, giving a symmetric two-sided coupling), and
`--trait-steps K` sets the within-branch resolution for a continuous trait (discrete traits use
their exact stochastic map). `--effect-gain` optionally scales a lineage's transfer activity by
the trait too (off by default).

It writes the gene-family files selected by `--write` (as [`genomes`](#choosing-the-output-and-the-rust-engine)),
and always adds **`traits.tsv`** / **`trait_tree.nwk`** (the trait at every node) and
**`coupling.tsv`** (the per-family weights and effect sizes — the trait↔gene linkage on record
for downstream inference). Reuse a precomputed trait instead of simulating one with
`--trait-file traits.tsv` (a `node`/`value` table over **every** node — tips and ancestors —
with numeric values, as `zombi2 trait` writes). See
[Trait-linked gene families](guide/trait-linked-genomes.md) for the model.

## Options

### `species`

| Option | Meaning |
| --- | --- |
| `--mode {backward,forward}` | reconstructed (default) or complete forward tree |
| `--diversification {constant,clads,diversity-dependent}` | [forward] rate process: constant-rate BD (default), ClaDS per-lineage rates, or diversity-dependent |
| `--birth` / `--death` | speciation / extinction rate(s) (default `1.0` / `0.3`); several values + `--shifts` = episodic; single λ₀/μ for clads/diversity-dependent |
| `--shifts` | episodic rate-shift ages (`K-1` for `K` rate values), present → past |
| `--clads-alpha` / `--clads-sigma` / `--turnover` | [clads] per-branch rate trend α, jump spread σ, and turnover ε=μ/λ |
| `--carrying-capacity` / `-K` | [diversity-dependent] carrying capacity K in λ(n)=λ₀(1−n/K) |
| `--clade-shift AGE BIRTH DEATH` | [forward] at AGE, a random lineage + descendants switch to (BIRTH, DEATH) (repeatable) |
| `--tips` | number of extant species (backward default `50`; forward: `--tips` **or** `--age`) |
| `--age` | tree age / timescale (backward default `1.0`; forward: `--tips` **or** `--age`) |
| `--age-type {crown,stem}` | interpret `--age` as crown (default) or stem age [backward] |
| `--sampling-fraction` | [forward] fraction of extant species sampled, ρ (default `1.0`) |
| `--fossilization` | [forward] fossil sampling rate ψ — fossilized birth–death (default `0`) |
| `--removal` | [forward] removal probability `r` on sampling (`r<1` → sampled ancestors; default `1.0`) |
| `--mass-extinction AGE FRACTION` | [forward] instantaneous pulse: fraction `FRACTION` of lineages die at age `AGE` before the present (repeatable) |
| `--ghosts` / `--ghost-method {rejection,htransform}` | [backward] un-prune extinct/unsampled ghost lineages |
| `--max-attempts` / `--max-lineages` | [forward] extinction-retry cap / live-lineage cap |
| `--seed` / `-o` / `--out` | RNG seed / output directory |

### `genomes`

| Option | Meaning |
| --- | --- |
| `--tree` / `-t` | input species tree in Newick format |
| `--rate-model {uniform,genome-wise,nucleotide}` | `uniform` (default, Rust): same per-copy rates for all families; `genome-wise` (Python): constant per-genome rates, linear growth; `nucleotide`: nucleotide-resolution genomes with variable-length structural events ([see below](#nucleotide-genomes-rate-model-nucleotide)) |
| `--dup` `--trans` `--loss` `--orig` | duplication / transfer / loss / origination rates (per copy; **per nucleotide** for `--rate-model nucleotide`) |
| `--initial-families` | number of gene families seeded at the root (default: 20) [`uniform`/`genome-wise`] |
| `--max-family-size` | growth cap — integer = absolute, decimal = fraction of N (e.g. `0.5`) [not used by `nucleotide`] |
| `--inversion` `--transposition` | [nucleotide] per-nucleotide inversion / transposition rates |
| `--initial-chromosomes` | [nucleotide] number of root chromosomes seeded at the root (default: 1) |
| `--root-length` `--extension` | [nucleotide] root chromosome length (nt) / geometric event-length parameter (mean `1/(1-extension)`) |
| `--gff FILE` | [nucleotide] a GFF3 annotation (optionally `.gz`) — copies the chromosome length + gene coordinates (overlaps trimmed) to start genic mode from a real genome; supersedes `--genes`/`--root-length`. `--gff-seqid ID` picks a sequence in a multi-record file |
| `--genes FILE` | [nucleotide] BED/TSV of gene intervals (`start end [name]`) on the root chromosome — enables *genic mode* (genes are never split; genes & intergenes recovered as separate tree sets) |
| `--pseudogenization` `--replacement` | [nucleotide, genic] probability a loss demotes a gene to intergene (sequence retained) / a transfer is a homologous replacement |
| `--write {profiles,trace,trees,events,transfers,summary,all}` | which files to write (one or more; default `profiles trees`); `profiles` alone → the counts-only fast path; `trace` (± `profiles`) → the compact `Events_trace.tsv` fast path — see below |
| `--sparse` | write the profile as a sparse `Profiles_sparse.tsv` instead of the dense matrix (needs `profiles` in `--write`) |
| `--annotate-species` | label internal gene-tree nodes `<gid>\|<species-branch>` (e.g. `g570\|i5`) |
| `--seed` / `-o` / `--out` | RNG seed / output directory |

Substitution branch lengths (sequence evolution) are a **separate step** — run
[`zombi2 sequence`](#sequence) on a genomes run, rather than a flag here.

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
| `--rate-model {uniform,family}` | `uniform` (default, Rust): shared scalar rates; `family` (Python): fit per-family rate means |
| `--family-shape` | [family] Gamma shape for per-family dispersion (default `2.0`) |
| `--n-sims` / `--accept` | [rejection] number of simulations (default `1000`) / accepted fraction (default `0.05`) |
| `--processes` | [rejection] parallel worker processes (default: serial) |
| `--smc` `--rounds` `--particles` `--quantile` | run ABC-SMC instead of rejection, with these controls |
| `--regression-adjust` | also write the regression-adjusted posterior (Beaumont 2002) |
| `--initial-families` / `--max-family-size` | families seeded per sim / growth cap (advised for `--rate-model family`) |
| `--seed` / `-o` / `--out` | RNG seed / output directory (required) |

### `sequence`

| Option | Meaning |
| --- | --- |
| `--genomes DIR` | a prior `genomes` output directory — reads its `species_tree.nwk` + `Events_trace.tsv` (run `genomes` with `trace` in `--write`) (required) |
| `--family-speed SIGMA` | per-family intrinsic substitution speed `~ LogNormal(0, SIGMA)`, constant per family (`0` = every family the same) |
| `--branch-speed SIGMA` | shared lineage clock — autocorrelated lognormal relaxed clock, drift `SIGMA` per `√time` (`0` = strict). Exclusive with `--branch-bins` |
| `--branch-bins R1,R2,...` | alternative lineage clock — the discrete-bin GTDB model: ordered rate multipliers, a Markov walk between adjacent bins (`--branch-switch-rate`, `--branch-up-bias`) |
| `--seed` / `-o` / `--out` | RNG seed / output directory (required) |

### `coevolve --couple traits:genes`

| Option | Meaning |
| --- | --- |
| `--couple traits:genes` | select the trait-conditioned-genes edge (required) |
| `--tree` / `-t` | input species tree in Newick format (required) |
| `--trait-model {bm,ou,eb,mk,threshold}` | trait to evolve then couple to gene families (default `bm`); its parameters are the [`trait`](#trait-a-phenotypic-trait) flags (`--sigma2`, `--alpha`/`--theta`, `--rate`, `--states`/`--ordered`/`--q-matrix`, `--thresholds`, …) |
| `--trait-file TSV` | reuse a precomputed trait instead — a numeric `node`/`value` table over **every** node (as `zombi2 trait` writes); overrides `--trait-model` |
| `--trait-center` | [discrete] center the state values around their mean (two-sided coupling; recommended for a binary trait) |
| `--trait-steps K` | [continuous] within-branch resolution — sub-segment each branch into K pieces (default `16`; ignored for discrete traits) |
| `--panel` | number of gene families in the panel (default `50`) |
| `--loss` `--trans` `--dup` `--orig` | panel base rates — baseline per-copy loss (default `0.5`), transfer/HGT gain (default `1.0`), duplication, origination |
| `--responsive SPEC` | which families respond: a count, a fraction (e.g. `0.3`), an id/index list (`F3,F7,12`), or `@FILE` (default `0.3`) |
| `--weight` / `--signed` | coupling weight of each responsive family (default `1.0`) / randomise its sign |
| `--effect-loss` | retention coupling strength: loss scales by `exp(-effect_loss · weight · trait)` (default `2.0`; `0` = uncoupled) |
| `--effect-gain` | optional donor-side HGT-activity coupling: transfer scales by `exp(effect_gain · trait)` (default `0`) |
| `--write {profiles,trace,trees,events,transfers,summary,all}` | which gene-family files to write (default `profiles trees`); `traits.tsv` / `trait_tree.nwk` / `coupling.tsv` are always written too |
| `--sparse` / `--annotate-species` | sparse profile table / label internal gene-tree nodes (as in `genomes`) |
| `--seed` / `-o` / `--out` | RNG seed / output directory (required) |

Run `zombi2 <command> -h` for the authoritative list.

## Choosing the output, and the Rust engine

The built-in model runs on the native **Rust** engine automatically — there is no flag to turn
it on, and no pure-Python fallback for it (which keeps results reproducible against one engine).
The `genomes` command therefore needs the compiled `zombi2_core` extension; it is **not** built
by `pip install`, so build it once (see [the Rust engine](guide/rust-engine.md)) or the command
exits with a build hint.

`--write` selects which files to write — any of `profiles`, `trace`, `trees`, `events`,
`transfers`, `summary`, or `all` (default: `profiles trees`); `species_tree.nwk` is always
written, and a component you don't ask for does no work (e.g. omitting `trees` skips the gene-tree
reconstruction). Asking for **only** `profiles` takes the Rust counts-only fast path (no
genealogy); `--sparse` then writes the profile as a scalable `Profiles_sparse.tsv` long table
instead of the dense matrix:

```bash
zombi2 genomes --tree out/species_tree.nwk --dup 0.2 --loss 0.25 --orig 0.5 --write all -o out/
zombi2 genomes --tree out/species_tree.nwk --dup 0.2 --loss 0.25 --orig 0.5 --write profiles --sparse -o out/
```

### The event trace — a fast intermediate for very large datasets

Between the counts-only `profiles` and the full genealogy sits **`--write trace`**. It writes a
single compact file, `Events_trace.tsv` (one row per event: origination / duplication / transfer
/ loss / speciation, with the gene-lineage ids), instead of the per-family `gene_family_events/`
directory — one file rather than one-per-family, so it scales to millions of families without an
inode explosion. Crucially it **skips gene-tree reconstruction and never materialises the
per-event objects**, which is the real bottleneck on large trees, so it runs at close to
counts-only speed while the gene trees remain **reconstructable later on demand** from the trace.

```bash
# fast: the event trace + the sparse profile, no gene trees
zombi2 genomes --tree out/species_tree.nwk --dup 0.2 --loss 0.25 --orig 0.5 \
    --write trace profiles --sparse -o out/

# reconstruct the gene trees too (pays for the genealogy, still writes the trace file)
zombi2 genomes --tree out/species_tree.nwk --dup 0.2 --loss 0.25 --orig 0.5 \
    --write trace trees -o out/
```

From the Python API the same thing is `simulate_genomes(tree, ..., output="trace")`, which returns
a `GenomeTrace` — a lazy handle whose `.profiles`, `.gene_trees()` and `.reconciliations()` build
only what you ask for.

## Nucleotide genomes (`--rate-model nucleotide`)

`--rate-model nucleotide` switches `genomes` to a **nucleotide-resolution** model: each genome is
a circular sequence that evolves by **variable-length structural events** — inversion, deletion,
tandem duplication, transposition, HGT transfer and origination. Genes are not predefined; they
**emerge** as *blocks* — maximal intervals of ancestral sequence with a single shared history — so
you still get a phylogenetic profile and per-block gene trees, but derived from real sequence
structure rather than an atomic gene-family model.

The shared `--dup` / `--trans` / `--loss` / `--orig` flags become **per-nucleotide** rates here
(so use small values, on the order of `1e-3`); `--inversion` and `--transposition` are the extra
structural events, `--root-length` sets the starting chromosome length, and `--extension` the mean
event length (`1/(1-extension)` nt). `--initial-chromosomes` seeds root chromosomes (default `1`).

```bash
zombi2 genomes -t out/species_tree.nwk --rate-model nucleotide \
    --inversion 0.001 --dup 0.0006 --loss 0.0006 --root-length 1000 --seed 1 -o out/
```

`--write profiles` writes the emergent block profile (`Profiles.tsv` / `Presence.tsv`, blocks ×
species) plus `blocks.tsv` and the per-leaf `Mosaics.tsv`, taking the fast Rust path; the default
`profiles trees` also writes the per-block `gene_trees/` and their reconciliations
(`Reconciled_complete.nwk` / `Reconciled_extant.nwk` / `Reconciliation_events.tsv`). `--sparse`
applies to the profile as usual. (The family-model `events` / `transfers` / `summary` outputs do
not apply to this model.)

**Genes & intergenes.** Pass `--genes genes.tsv` (a BED/TSV of `start end [name]` intervals on the
root chromosome) to declare genes explicitly. Event breakpoints then fall only in intergene
positions, so **genes are never split** — each gene is one block, each intergene stretch fragments
into intergene blocks. `--pseudogenization P` makes a loss that hits a gene demote it to intergene
with probability `P` (sequence retained, a state change in the gene's tree); `--replacement P`
makes a transfer a homologous replacement (the copy replaces the recipient's syntenic locus, found
via flanking genes; additive when there is no homolog). Genic mode runs on the Python engine and
adds `genes.tsv` (the annotation), `Gene_trees/` and `Intergene_trees/` (the two tree sets), a
`kind`/`gene_id` column in `blocks.tsv`, and `Pseudogenizations.tsv`.

```bash
zombi2 genomes -t out/species_tree.nwk --rate-model nucleotide \
    --genes genes.tsv --pseudogenization 0.3 --replacement 0.4 \
    --inversion 0.001 --loss 0.0008 --write profiles trees -o out/
```

**Starting from a real genome.** `--gff FILE` copies a real annotation's chromosome length and
gene coordinates (the intergenes are the gaps), so a simulation can start from an actual bacterium.
Overlapping genes — common in bacteria — are trimmed to be disjoint (each gene's start clipped to
the previous gene's end; a swallowed gene dropped), and the count is reported. The file may be
gzipped; for a chromosome-plus-plasmids file the most-annotated sequence is used unless
`--gff-seqid ID` selects another. Genes keep their annotation names (locus tag / `Name`).

```bash
# evolve the E. coli K-12 chromosome (from its RefSeq GFF) along a tree
zombi2 genomes -t out/species_tree.nwk --rate-model nucleotide \
    --gff ecoli.gff --inversion 2e-6 --loss 1.5e-6 --pseudogenization 0.3 \
    --write profiles trees -o out/
```

**Sequences and ancestral genomes.** `--write ancestral` simulates the DNA and reconstructs the
genome at **every** node of the tree (the root reproduces the input genome). It writes
`Architecture/<node>.tsv` (each node's oriented gene/intergene mosaic), gzipped
`Genomes/<node>.fasta.gz` (the full assembled DNA of every node), and `Gene_alignments/<gene>.fasta`
(the extant per-gene alignments). Pick the substitution model with `--subst-model {jc69,k80,hky85,
gtr}` (`--kappa`, `--base-freqs`, `--gtr-rates`, `--gamma-shape`, `--subst-rate`). `--genome-fasta
FILE` seeds the root from the real genome DNA (so the reconstructed root is byte-identical to the
input); without it, root sequences are drawn at random.

```bash
zombi2 genomes -t out/species_tree.nwk --rate-model nucleotide \
    --gff ecoli.gff --genome-fasta ecoli.fna --subst-model hky85 --subst-rate 0.05 \
    --write ancestral -o out/
```

## Scope

The CLI covers the common **uniform-rate** case. For family-sampled or genome-wise rates,
custom transfer mechanics, ordered genomes, or replicate parallelism, use the Python API
(see [gene families & rates](guide/gene-families.md) and the other guides).
