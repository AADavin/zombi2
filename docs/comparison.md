# Comparison with ZOMBI1

ZOMBI2 is a ground-up redesign of the original
[ZOMBI](https://github.com/AADavin/Zombi). It keeps the ideas that made ZOMBI useful —
duplication/transfer/loss/origination along a species tree, per-family rate variation,
the output format — while rebuilding the internals around clean interfaces.

## What carries over

- **The DTLO model** — duplication, transfer, loss, and origination of gene families
  along a species tree.
- **Per-family sampled rates** — the ZOMBI1 style where each family draws its own D/T/L
  from distributions is available as `FamilySampledRates(...)`.
- **Output format** — `Genomes.write("out/")` produces the familiar files: the species
  tree, per-family event tables, reconstructed complete and extant gene trees, a transfer
  log, a per-family summary, and presence/copy-number matrices. See
  [gene trees & output](guide/genomes.md#gene-trees-output).

## What's new or different

- **Interface-first architecture.** Species-tree models, rate models, genome
  representations, and transfer mechanics are all pluggable interfaces. New behaviour
  arrives as a *subclass* the simulator already knows how to drive, rather than as edits to
  a monolithic script. See [extending ZOMBI2](contributing/adding-a-model.md).
- **Species trees, backward *or* forward.** Sample the *reconstructed* tree directly,
  conditioned on exactly `N` extant tips (cheap and exact), **or** simulate the *complete*
  tree **forward** in time — extinct lineages included — with fossilized birth–death
  (dated/fossil tips) and incomplete sampling. See [species trees](guide/species-trees.md).
- **Trait levels.** Evolve phenotypic traits along the tree — Brownian motion,
  Ornstein–Uhlenbeck, early burst, Mk, threshold, and DEC biogeography. See
  [trait evolution](guide/traits.md).
- **Nucleotide genomes.** Nucleotide-resolution genomes with variable-length structural
  events — inversions, transpositions, indels — genes and intergenes, and GFF import. See
  [nucleotide genomes](guide/genomes.md).
- **Coevolution modes.** Couple any two levels along a directed edge (`coevolve --couple
  driver:target`) — state-dependent diversification (SSE), cladogenetic change, key
  innovations, and trait-linked gene families. See [coevolution](guide/coevolution.md).
- **Principled growth control.** A hard `max_family_size` cap (absolute or a fraction of
  the number of species) and a soft logistic `carrying_capacity`. See
  [bounding growth](guide/genomes.md).
- **Ghost lineages.** `add_ghost_lineages` un-prunes the reconstructed tree, grafting back
  the extinct/unsampled lineages the backward process leaves out. See
  [ghost lineages](guide/species-trees.md).
- **Scale.** The built-in model runs on a native **Rust** engine, making large trees
  (thousands to tens of thousands of tips) practical; `run_replicates` parallelises
  independent replicates across cores. See [the Rust engine](guide/rust-engine.md)
  and [running in parallel](guide/parallel.md).
- **Distributions.** Rate distributions accept the built-ins (`Gamma`, `Exponential`,
  `LogNormal`, `Uniform`, `Fixed`), any `scipy.stats` frozen distribution, or a plain
  `rng -> float` callable.
