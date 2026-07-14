# Changelog

All notable changes to ZOMBI2 are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and this
project aims to follow [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- **Codon substitution models** — `zombi2 sequence --subst-model gy94`/`mg94` evolve in-frame coding
  DNA over the 61 sense codons with `dN/dS` set directly by `--omega` (`<1` purifying, `1` neutral,
  `>1` positive selection) and a ti/tv bias `--kappa`. GY94 (Goldman & Yang 1994) weights by the
  target-codon frequency, MG94 (Muse & Gaut 1994) by the introduced-nucleotide frequency; codon
  frequencies come from `F1×4`/`F3×4`/`F61`. Stop codons are never produced, and `--seq-length`
  counts codons. On the API as `zombi2.gy94` / `zombi2.mg94` (and `zombi2.sequences.codon_models`,
  with `translate` and `expected_dnds`); validated against detailed balance, exact `omega` recovery,
  the synonymous ti/tv ratio, and equilibrium base composition.
- **Codon site models (dN/dS across sites)** — `--omega-model {m1a,m2a,m3,m7,m8}` lets `dN/dS` vary
  among codon sites (the Nielsen–Yang site models: M1a/M7 purifying/neutral nulls, M2a/M8 add a
  positive-selection class, M3 discrete, M7/M8 a discretised `Beta(p,q)`). Each `ω` class is a codon
  matrix sharing one mutation process and stationary distribution, normalised on a single shared
  scale so purifying classes evolve slower and the gene's genome-wide `dN/dS` is the
  proportion-weighted mean `ω`. On the API as `zombi2.m1a`…`zombi2.m8` and `zombi2.CodonSiteModel`;
  mutually exclusive with `--gamma-shape`.
- **Per-branch event tables** — `zombi2 genomes --write branch_events` writes `Branch_events.tsv`,
  one row per species-tree branch with the count of each event that fired on it (origination,
  duplication, `transfer_out`/`transfer_in`, loss, and inversion/transposition for ordered
  genomes) plus a `total`. An `is_extant` column, derived from node times, makes the extant-tree
  view a simple filter. Exposed on the API as `zombi2.genomes.simulation.branch_events_table`.
- **BED gene annotations** — `zombi2 genomes --genome-model nucleotide --write bed` (genic mode)
  writes `genes.bed` for the root genome and `BED/<node>.bed` for every node's genome after
  rearrangements, in standard BED6 (0-based half-open), ready for genome browsers and `bedtools`.
  Per-node BED chromosomes match `Genomes/<node>.fasta.gz` from `--write ancestral`.

### Changed

- Documented that `--transposition` (a segment moved elsewhere in the genome) applies to both the
  `ordered` and `nucleotide` genome models, not only `nucleotide`.
- **`--genome-resolution` replaces `--genome-model`** (naming consolidation C4). The genome
  representation axis (unordered / ordered / nucleotide) is a *resolution*, not a "level" or a
  stochastic "model"; `--genome-model` is kept as an accepted (deprecated) alias. The metavar is now
  `RESOLUTION`, and "level" is reserved for the four simulation domains (species, genomes, traits,
  sequences).
- **Coevolve nodes standardise on `genomes`** (C6): `--couple traits:genomes`, `genomes:species`,
  etc. The old node spelling `genes` is still accepted but warns — one domain word everywhere.
- **Rate vocabulary clarified in the docs and manual** (naming consolidation C2/C3). The guide,
  the manual, and `conventions.md` now state a single two-word rule — a **rate** is a quantity per
  unit time, a **modifier** is a dimensionless multiplier on a rate — and reserve **odds** for the
  undated-ALE tools alone (`tools reconcile --model undated/reldated`, `tools simulate`), where the
  D/T/L parameters are per-branch odds, not rates. Retired class/function names (`PerGenomeRates`,
  `BranchRates`, `BranchModifier`, `read_branch_rates`, the `shared` rate model) and per-branch
  rate/modifier phrasing were updated to the `per-lineage` / `PerLineageRates` / `LineageRates` /
  `read_lineage_rates` vocabulary throughout the prose. (The sequence-clock `--branch-*` /
  `--family-speed` flags are σ/spread parameters and are left unchanged pending a naming decision;
  see `docs/design/naming-consolidation.md` §C3.)

### Deprecated

- **Renamed rate names retired from the public API surface** (naming consolidation, see
  `docs/design/naming-consolidation.md`). The five backwards-compatible aliases
  `SharedRates`→`PerCopyRates`, `PerGenomeRates`→`PerLineageRates`, `BranchRates`→`LineageRates`,
  `BranchModifier`→`LineageModifier`, and `read_branch_rates`→`read_lineage_rates` still work but now
  emit a `DeprecationWarning` and no longer appear in `zombi2.__all__`, `dir(zombi2)`, or the API
  reference (they resolve via a PEP-562 `__getattr__`; the deep-module spelling
  `zombi2.genomes.rates.SharedRates` stays silent). The deprecated CLI flags `--rate-model` and
  `--initial-chromosomes` are hidden from `--help` (still accepted, with a warning). All are
  scheduled for removal in **0.4.0**.

### Removed

- **Experimental protein-language-model selection** — `zombi2 experimental selection`,
  `zombi2.experimental.selection` (+ `codon_selection` / `genome_selection` /
  `nucleotide_selection` / `realism`), and the `zombi2[selection]` extra with its `torch` /
  `fair-esm` dependencies have been removed (naming consolidation, `docs/design/naming-consolidation.md`
  §C9). It was the ESM2-critic realism feature; the active protein-realism work continues off-repo.
  `experimental ils` is unaffected and is now the sole experimental model. "Selection" now
  unambiguously refers to the core codon `dN/dS` models (`--subst-model gy94/mg94`, `--omega`).

### Fixed

- `--initial-chromosomes` (deprecated) no longer silently overrides the canonical `--n-chromosomes`;
  passing both with conflicting values is now an error.

## [0.2.0] - 2026-07-07

First public release. ZOMBI2 is a ground-up redesign of
[ZOMBI](https://github.com/AADavin/Zombi).

### Added

- **Species trees** — backward (reconstructed) and forward (complete) birth–death and Yule,
  episodic/skyline rate shifts, fossilized birth–death and incomplete sampling, ClaDS /
  diversity-dependent / clade-shift diversification, mass extinctions, and ghost lineages.
  Backward trees scale to millions of tips on a laptop via a native Rust engine.
- **Gene families** — duplication, transfer, loss and origination along a species tree, with
  uniform, family-sampled and genome-wise rate models; ordered chromosomes with inversions and
  transpositions; and nucleotide-resolution genomes. Output as full event logs, compact event
  traces, or counts-only sparse profiles, with parallel replicates.
- **Traits** — Brownian motion, Ornstein–Uhlenbeck, early burst, Mk, threshold,
  correlated-binary, hidden-state and multi-optimum models, Pagel transforms, and DEC
  biogeography, evolved over a phylogeny.
- **Sequences** — a gene × lineage relaxed-clock family (strict, UCLN, UGAM, white-noise,
  autocorrelated-lognormal, CIR) that rescales gene trees from time into substitutions/site,
  plus nucleotide substitution models (JC, K80, HKY, GTR + Gamma).
- **Coevolution** — a unified `coevolve` command that couples species, traits and genes along
  six directed edges via `--couple driver:target`.
- **Command-line interface** — `zombi2 species | genomes | trait | sequence | coevolve | tools`,
  with grouped, sectioned help and a per-run reproducibility manifest (version, seed, full command
  line and resolved parameters).
- **Packaging** — MIT license, `CITATION.cff`, and PyPI metadata.

### Notes

- ABC inference of DTL rates from an empirical copy-number profile (`zombi2.matching`) is **not**
  part of this release: it has been moved to a separate future-extensions archive and is neither
  importable from `zombi2` nor exposed on the command line.
- The default (built-in) gene-family engine is a compiled Rust extension; build it once with
  maturin (see [installation](https://aadavin.github.io/zombi2/docs/installation/)).
