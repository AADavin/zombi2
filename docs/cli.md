# Command-line interface

Installing the package puts a `zombi2` command on your PATH — a thin wrapper over the library.

```bash
# 1. a species tree -> out/species_tree.nwk  (runs with defaults)
zombi2 species -o out/

# 2. gene families along that tree (or any Newick tree)
zombi2 genomes --tree out/species_tree.nwk \
    --dup 0.2 --trans 0.1 --loss 0.25 --orig 0.5 --seed 42 -o out/

# 3. a phenotypic trait along that tree
zombi2 trait --tree out/species_tree.nwk --model ou --alpha 2 --theta 5 --seed 1 -o out/

# 4. DNA sequences along the gene trees (from a genomes run written with --write trace)
zombi2 sequence --genomes out/ --subst-model hky85 --kappa 4 --seed 7 -o out/seq/
```

`species` writes `species_tree.nwk` (plus `species_nodes.tsv` and a run log; forward mode also
writes the extant-only tree); `genomes` reads a tree from `--tree` and writes the
full output (see [gene trees & output](guide/genomes.md#gene-trees-output)); `trait` reads a tree and
writes the tip and ancestral trait values (see [trait evolution](guide/traits.md)). To **couple**
levels so one drives another, see [`coevolve`](guide/coevolution.md).

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
tree-wide survival pulse — at `AGE` before the present, each live lineage dies with probability
`FRACTION` — and is repeatable for several pulses (needs `--age`). See the
[advanced diversification](guide/species-trees.md#mass-extinctions-instantaneous-pulses) catalog for details.

```bash
zombi2 species --mode forward --age 5 --mass-extinction 2.5 0.75 -o out/        # 75% die at age 2.5
```

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
sampler. See [ghost lineages](guide/species-trees.md):

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

## `sequence` — substitution branch lengths and DNA/protein alignments

`genomes` produces gene trees in **time** units. `sequence` is a separate step that (1) **rescales**
them into **substitutions/site** under a gene × lineage clock, and (2) optionally **simulates a
sequence alignment** along each rescaled tree — so you can retune the rate model or draw sequences
without re-simulating gene content. It reads a prior `genomes` run's `species_tree.nwk` and
`Events_trace.tsv` (so run `genomes` with `trace` in `--write`), replays the genealogy, and writes
the phylograms:

```bash
zombi2 genomes  -t out/species_tree.nwk --dup 0.2 --trans 0.1 --loss 0.2 --orig 0.5 \
    --write trace profiles -o run/
zombi2 sequence --genomes run/ --branch-speed 0.4 --family-speed 0.5 -o run/
```

The rescaling model is `rate(family g, species branch b) = R_b · s_g`: a **shared lineage clock**
`R_b` (either `--branch-speed` autocorrelated lognormal, or `--branch-bins` for the discrete-bin
GTDB model) times a **per-family speed** `s_g ~ LogNormal(0, --family-speed)`. It writes
`gene_trees/<family>_extant_subst.nwk` (+ `_complete_subst.nwk`) and, for reproducibility,
`gene_family_speeds.tsv` / `branch_rates.tsv`. See
[Sequences](guide/sequences.md) for the model.

**Simulating sequences.** Add `--subst-model` and `sequence` also evolves a DNA or protein
alignment down each rescaled gene tree, writing `alignments/<family>.fasta` alongside the
phylograms. Pick a DNA model (`jc69`, `k80`, `hky85`, `gtr`) or a protein one (`poisson`, `lg`,
`wag`, `jtt`, `dayhoff`) — DNA vs protein is auto-detected from the name. `--seq-length` sets the
alignment length (default `300`), `--gamma-shape` adds across-site rate heterogeneity, and the
DNA-model knobs are `--kappa` (k80/hky85 transition/transversion ratio), `--base-freqs A C G T`
(hky85/gtr) and `--gtr-rates AC AG AT CG CT GT` (gtr). Seed each family's root from real sequences
with `--root-fasta FILE` (a FASTA keyed by family id, whose per-family length overrides
`--seq-length`) instead of a random draw:

```bash
zombi2 sequence --genomes run/ --subst-model hky85 --kappa 4 --seed 7 -o seq/
```

This writes one `alignments/<family>.fasta` per gene family (plus the rescaled `gene_trees/`).

## `coevolve` — coupled models

Everything above is a **pipeline**: build a tree, then overlay a trait or gene families on it,
each stage independent. `coevolve` is for the case where one level **drives** another. A coupling
is a **directed edge** `driver:target` — the driver's state modulates the target's rates — selected
with the repeatable flag `coevolve --couple driver:target`. There are six edges among the three
levels **species** / **traits** / **genes** (e.g. `traits:species` = a trait sets speciation, the
SSE models; `traits:genes` = a trait sets gene loss/gain; and their reverses). Edges that point
*into* species grow the tree as an output (forward-only, no `-t`); the rest are overlays on a given
tree.

```bash
# a trait sets speciation/extinction and the tree is grown jointly with it (BiSSE)
zombi2 coevolve --couple traits:species --sse-model bisse \
    --lambda0 1 --lambda1 3 --q01 0.1 --q10 0.1 --tips 200 --seed 1 -o out/
```

Add **`--null {neutral,cid,timing}`** to any edge to generate its matched *decoupled* null instead
of the coupled model — the arrow cut but the target's variance kept, for calibrating a detector's
false-positive rate (writes a `null_manifest.tsv` recording what was cut). See **[null models of
coevolution](guide/coevolution_nulls.md)**.

See **[coevolution models](guide/coevolution.md)** for the full reference — all six edges, their
parameters, the joint (both-arrow) models, and the CLI options.

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
| `--genome-model {unordered,ordered,nucleotide}` | genome level: `unordered` (default) evolves gene families with no positional structure; `ordered` places genes on a chromosome where order matters (adds inversion/transposition on gene segments); `nucleotide`: nucleotide-resolution genomes with variable-length structural events ([see below](#nucleotide-genomes-genome-model-nucleotide)) |
| `--rate-model {shared,per-genome,family}` | rate heterogeneity for the unordered level: `shared` (default, Rust): same per-copy rates for all families; `per-genome` (Python): constant per-genome rates, linear growth; `family` (Python): each family its own rates, from `--family-rates` |
| `--dup` `--trans` `--loss` `--orig` | duplication / transfer / loss / origination rates (per copy; **per nucleotide** for `--genome-model nucleotide`) |
| `--conversion` `--conversion-bias` | intra-genome gene-conversion rate (per copy; one copy overwrites another of the same family — concerted evolution) and its donor directionality in `[0,1]` (0 = uniform donor, 1 = the oldest copy). Unordered genomes, `--rate-model shared`; runs on the Python engine |
| `--family-rates FILE` | TSV of explicit per-family rates (`family duplication transfer loss`); selects `--rate-model family`; unlisted families fall back to `--dup/--trans/--loss` [unordered, Python] |
| `--branch-rates FILE` | TSV of per-branch transfer `emission` (donation-rate factor) and/or `receptivity` (absorption weight) (`branch emission receptivity`, either optional) [unordered; receptivity-only stays on Rust] |
| `--initial-families` | number of gene families seeded at the root (default: 20) [`--genome-model unordered`] |
| `--max-family-size` | growth cap — integer = absolute, decimal = fraction of N (e.g. `0.5`) [not used by `--genome-model nucleotide`] |
| `--inversion` `--transposition` | [ordered/nucleotide] inversion / transposition (a segment moved elsewhere in the genome) rates — per gene copy for `ordered`, per nucleotide for `nucleotide` |
| `--n-chromosomes N` | [ordered/nucleotide] number of chromosomes seeded at the root (default 1). [ordered] the root's initial families are spread across them; [nucleotide] each is an independent full-length copy of the root chromosome (its own source). Gene rearrangements stay within a chromosome, while a transfer may land a copy on any chromosome. With `N > 1` the run also writes the layout (`Gene_order.tsv` for ordered / `Chromosomes.tsv` for nucleotide) |
| `--linear-chromosomes` | [ordered] chromosomes are linear — segments never wrap the origin (default: circular, as for bacteria). Nucleotide chromosomes are always circular |
| `--fission` `--fusion` | [ordered/nucleotide] chromosome-tier rates (per chromosome, default 0): a chromosome splits in two (linear: one breakpoint; circular: two) / two chromosomes merge into one |
| `--chromosome-origination` `--chromosome-loss` | [ordered/nucleotide] chromosome-tier rates (default 0): a de-novo replicon (a *plasmid*) appears (per genome) / a whole chromosome and its genes are lost (per chromosome). Any chromosome-tier rate also writes `Karyotype_trace.tsv` (the fission/fusion/origination/loss genealogy) |
| `--initial-chromosomes` | [nucleotide] **deprecated** alias for `--n-chromosomes` |
| `--root-length` `--mean-length` | [ordered/nucleotide] root chromosome length (nt) / mean inversion–transposition segment length (geometric; genes for ordered, nt for nucleotide) |
| `--gff FILE` | [nucleotide] a GFF3 annotation (optionally `.gz`) — copies each sequence's length + gene coordinates (overlaps trimmed) to start genic mode from a real genome; supersedes `--genes`/`--root-length`. A **multi-sequence** GFF seeds **one chromosome per sequence** (a chromosome + its plasmids); `--gff-seqid ID` instead picks a single sequence. (`--write ancestral`/`bed` + `--genome-fasta` need a single sequence.) |
| `--genes FILE` | [nucleotide] BED/TSV of gene intervals (`start end [name]`) on the root chromosome — enables *genic mode* (genes are never split; genes & intergenes recovered as separate tree sets) |
| `--pseudogenization` `--replacement` | [nucleotide, genic] probability a loss demotes a gene to intergene (sequence retained) / a transfer is a homologous replacement |
| `--write {profiles,trace,trees,events,transfers,summary,branch_events,layout,karyotype,ancestral,bed,all}` | which files to write (one or more; default `profiles trees`); `profiles` alone → the counts-only fast path; `trace` (± `profiles`) → the compact `Events_trace.tsv` fast path; `branch_events` → per-species-branch event counts (`Branch_events.tsv`, with an `is_extant` flag); `layout`/`karyotype` are ordered-only (`Gene_order.tsv` / `Karyotype_trace.tsv`, added automatically for a multi-chromosome or fission/fusion run); the nucleotide model writes the equivalent `Chromosomes.tsv` / `Karyotype_trace.tsv` automatically for a multi-chromosome or chromosome-tier run; `ancestral`/`bed` are nucleotide-only ([see below](#nucleotide-genomes-genome-model-nucleotide)) — see below |
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

### `sequence`

| Option | Meaning |
| --- | --- |
| `--genomes DIR` | a prior `genomes` output directory — reads its `species_tree.nwk` + `Events_trace.tsv` (run `genomes` with `trace` in `--write`) (required) |
| `--family-speed SIGMA` | per-family intrinsic substitution speed `~ LogNormal(0, SIGMA)`, constant per family (`0` = every family the same) |
| `--branch-speed SIGMA` | shared lineage clock — autocorrelated lognormal relaxed clock, drift `SIGMA` per `√time` (`0` = strict). Exclusive with `--branch-bins` |
| `--branch-bins R1,R2,...` | alternative lineage clock — the discrete-bin GTDB model: ordered rate multipliers, a Markov walk between adjacent bins (`--branch-switch-rate`, `--branch-up-bias`) |
| `--subst-model MODEL` | simulate an alignment per family: DNA (`jc69`, `k80`, `hky85`, `gtr`) or protein (`poisson`, `lg`, `wag`, `jtt`, `dayhoff`); auto-detected. Omit to only rescale the trees (no sequences) |
| `--seq-length N` | alignment length in sites (default `300`); ignored where `--root-fasta` seeds a family's root |
| `--root-fasta FILE` | FASTA (optionally `.gz`) of per-family root sequences keyed by family id — seeds each family's root instead of a random draw; its length overrides `--seq-length` per family |
| `--gamma-shape ALPHA` | discrete-Gamma across-site rate heterogeneity shape (default: none) |
| `--kappa K` | [DNA k80/hky85] transition/transversion ratio (default `2.0`) |
| `--base-freqs A C G T` | [DNA hky85/gtr] equilibrium base frequencies (default equal) |
| `--gtr-rates AC AG AT CG CT GT` | [DNA gtr] the 6 exchangeabilities (default all `1`) |
| `--seed` / `-o` / `--out` | RNG seed / output directory (required) |

Run `zombi2 <command> -h` for the authoritative list. For the `coevolve` options, see
[coevolution models](guide/coevolution.md#command-line).

## Choosing the output, and the Rust engine

The built-in model runs on the native **Rust** engine automatically — there is no flag to turn
it on, and no pure-Python fallback for it (which keeps results reproducible against one engine).
The `genomes` command therefore needs the compiled `zombi2_core` extension, but `pip install
zombi2` pulls a prebuilt engine wheel (Linux/macOS/Windows, CPython 3.10+), so nothing extra is
needed. From a source checkout you build it once from `rust/` (see
[the Rust engine](guide/rust-engine.md)); if it is missing the command exits with a build hint.

`--write` selects which files to write — any of `profiles`, `trace`, `trees`, `events`,
`transfers`, `summary`, `branch_events`, or `all` (default: `profiles trees`); `species_tree.nwk`
is always written, and a component you don't ask for does no work (e.g. omitting `trees` skips the
gene-tree reconstruction). `branch_events` writes `Branch_events.tsv` — one row per species-tree
branch with the count of each event that fired on it (D/T/L/O, plus inversion/transposition for
ordered genomes; transfers split into `transfer_out`/`transfer_in`) and an `is_extant` flag, so the
extant-tree view is a filter. The nucleotide model adds two more parts, `ancestral` (simulate DNA +
reconstruct every node's genome) and `bed` (BED gene annotations), described in the nucleotide
section. Asking for **only** `profiles` takes the Rust counts-only fast path (no
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

## Nucleotide genomes (`--genome-model nucleotide`)

`--genome-model nucleotide` switches `genomes` to a **nucleotide-resolution** model — genomes evolve
by variable-length structural events (inversion, deletion, tandem duplication, transposition,
transfer, origination), genes emerge as blocks, and the model can start from a real GFF genome and
even simulate ancestral DNA. It is a substantial model with its own flags; see
[nucleotide genomes](guide/genomes.md) for the full walkthrough.

## Scope

The CLI covers the common **shared-rate** case (and per-genome rates via `--rate-model
per-genome`). For family-sampled rates, custom transfer mechanics, ordered genomes, or replicate
parallelism, use the Python API (see [gene families & rates](guide/genomes.md) and the
other guides).
