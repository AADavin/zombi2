# Comparison with ZOMBI-1

ZOMBI2 is a ground-up redesign of the original
[ZOMBI](https://github.com/AADavin/Zombi). It keeps the ideas that made ZOMBI useful —
duplication/transfer/loss/origination along a species tree, per-family rate variation,
the output format — while rebuilding the internals around clean interfaces.

## What carries over

- **The DTLO model** — duplication, transfer, loss, and origination of gene families
  along a species tree.
- **Per-family sampled rates** — the ZOMBI-1 style where each family draws its own D/T/L
  from distributions is available as `z.FamilySampledRates(...)`.
- **Output format** — `Genomes.write("out/")` produces the familiar files: the species
  tree, per-family event tables, reconstructed complete and extant gene trees, a transfer
  log, a per-family summary, and presence/copy-number matrices. See
  [gene trees & output](guide/gene-trees-and-output.md).

## What's new or different

- **Interface-first architecture.** Species-tree models, rate models, genome
  representations, and transfer mechanics are all pluggable interfaces. New behaviour
  arrives as a *subclass* the simulator already knows how to drive, rather than as edits to
  a monolithic script. See [extending ZOMBI2](guide/extending.md).
- **Backward, conditioned species trees.** The species tree is a reconstructed birth–death
  process **conditioned on the number of extant tips**, so you get exactly `N` species at a
  chosen age — cheaply and exactly — instead of running a forward process until it happens
  to land on `N`. See [species trees](guide/species-trees.md).
- **More rate models out of the box.** Beyond uniform and family-sampled rates, there are
  genome-wise rates and rate variation across the branches of the tree. See
  [rate variation](guide/rate-variation.md).
- **Richer transfers.** A `z.TransferModel` controls additive vs replacement transfers,
  phylogenetic-distance-weighted recipient choice, and self-transfers. See
  [transfers](guide/transfers.md).
- **Gene order and rearrangements.** An optional `z.OrderedGenome` puts genes on an ordered
  chromosome and adds inversions and transpositions. See
  [ordered genomes](guide/ordered-genomes.md).
- **Principled growth control.** A hard `max_family_size` cap (absolute or a fraction of
  the number of species) and a soft logistic `carrying_capacity`. See
  [bounding growth](guide/growth.md).
- **Scale.** A modern Gillespie core plus an optional **Rust** engine make large trees
  (thousands to tens of thousands of tips) practical; `z.run_replicates` parallelises
  independent replicates across cores. See [the Rust fast path](guide/rust-fast-path.md)
  and [running in parallel](guide/parallel.md).
- **Distributions.** Rate distributions accept the built-ins (`z.Gamma`, `z.Exponential`,
  `z.LogNormal`, `z.Uniform`, `z.Fixed`), any `scipy.stats` frozen distribution, or a plain
  `rng -> float` callable.

## On the roadmap

The interface-first design is meant to make the next models drop-in rather than rewrites —
for example further species-tree models (episodic/skyline is already here;
diversity-dependent and fossilized birth–death are planned) and gene-family **coupling**
(a Potts-style model of non-independence between families). See the
[species-tree models roadmap](species_tree_models.md).
