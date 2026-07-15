# Changelog

All notable changes to ZOMBI2 are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and this
project aims to follow [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- **Opportunity knob for species diversification (`per=`)** — `BirthDeath(birth, death, per="shared")`
  (CLI `zombi2 species --per shared`, alias `--rate-per`) grows a species tree under one *shared*
  diversification clock: the total speciation rate is a fixed `birth` regardless of how many lineages
  stand, so diversity grows **linearly** rather than exponentially (a diversity-dependent process in
  disguise, `λ(N)=birth/N`). The default `per="lineage"` is the usual per-lineage birth–death. This is
  the first step (Part 3, phase A) of making the rate *opportunity* — the unit the clock rides on —
  a selectable knob at every level; see `docs/design/opportunity-knob.md`. Forward-only; same seed →
  byte-identical to the previous `SharedBirthDeath`.
- **Opportunity knob for gene-family rates (`Rates(per=…)`)** — the gene D/T/L rate models unify under
  one `zombi2.Rates(duplication=…, …, per="copy"|"lineage")` (Part 3, phase B). `per="copy"` (default)
  scales the total rate by copy number → exponential families (the built-in Rust model); `per="lineage"`
  is a constant rate per genome → linear families. `FamilySampledRates` gains the same `per=`, so
  per-family heterogeneity and the opportunity are now **independently selectable** (per-family × per-
  lineage no longer needs the modifier route). CLI `--rate-per` gains a `--per` alias. Same seed →
  byte-identical (Rust and Python paths).
- **Per-event opportunity mixing (`Per(unit, rate)`)** — each of a gene family's duplication / loss /
  transfer rates can carry its **own** opportunity, overriding the model-level `per` (Part 3, phase D).
  `Rates(duplication=Per("shared", 0.5), loss=Per("copy", 0.3))` is a **self-limiting** family: a shared
  (tree-wide) duplication clock with per-copy loss, so births are capped while deaths grow with copy
  number and the family stays bounded near `copies ≈ dup/loss` instead of exploding. Exposed as
  `zombi2.Per`; API-only for now. Rearrangements / `carrying_capacity` require every event per-copy, and
  a shared *transfer* clock is not yet supported. Byte-identical for every non-mixed model.
- **Shared gene-family clock (`Rates(per="shared")`)** — the gene analogue of `SharedBirthDeath`
  (Part 3, phase C). Duplication and loss become **one tree-wide clock per family**: the total rate is
  a constant regardless of how many lineages carry the family or how many copies it holds, so the
  family's size grows **linearly and pooled** (`#events ≈ base × time`, independent of the tree's
  size) — a fire is localised to a copy chosen uniformly across the whole family. CLI
  `zombi2 genomes --per shared`. Unordered genomes only for now; duplication/loss (transfer,
  rearrangements, and rate *modifiers* — `LineageRates` / `--lineage-rates` — are rejected, since a
  modifier would silently bypass the shared clock); origination stays per lineage. Implemented as a
  "shared pool" beside the per-branch Gillespie, inert for every other model → those stay byte-identical.
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
- **Sequence-clock flags folded into `--clock`** (naming consolidation). The discrete-bin knobs are
  now `--clock-bins` / `--clock-switch-rate` / `--clock-up-bias` (the old `--branch-bins` /
  `--branch-switch-rate` / `--branch-up-bias` remain as accepted, deprecated aliases), and
  `--branch-speed SIGMA` is deprecated in favour of `--clock autocorrelated-lognormal --clock-sigma SIGMA`
  (byte-identical). This retires the last "branch" residue on the CLI and the redundant parallel
  clock interface — the `--clock` model selector plus its `--clock-*` parameters is now the one way.
  The old spellings still work but warn.
- **Coevolve `traits:genomes` classes use one stem** (C8): `TraitLinkedRates`, `TraitLinkedResult`,
  and `simulate_trait_linked_genomes` were renamed to `TraitGeneRates`, `TraitGeneResult`, and
  `simulate_trait_conditioned_genomes` — matching the edge's config (`TraitGeneCoupling`) and its
  joint model (`TraitGeneFeedback`), so the whole edge reads with one `TraitGene*` stem. The old
  names still work but warn (removed in 0.4.0).
- **CLI commands are plural nouns** (C5): `zombi2 traits` and `zombi2 sequences` (matching `species`
  / `genomes`, the packages, and the guide). The singular `trait` / `sequence` still work but warn.
  The run-manifest filename follows the command, so these runs now write `traits.log` / `sequences.log`
  (were `trait.log` / `sequence.log`).
- **Output filenames are `lowercase_snake`** (naming consolidation C7). Every file and directory a
  run writes is lowercased: `Profiles.tsv`→`profiles.tsv`, `Presence.tsv`→`presence.tsv`,
  `Events_trace.tsv`→`events_trace.tsv`, `Transfers.tsv`→`transfers.tsv`,
  `Reconciled_complete.nwk`/`Reconciled_extant.nwk`→`reconciled_*.nwk`,
  `Reconciliation_*.tsv`→`reconciliation_*.tsv`, `Tree_distances.tsv`→`tree_distances.tsv`,
  `RED.tsv`→`red.tsv`, `Karyotype_trace.tsv`, `Chromosomes.tsv`, `Mosaics.tsv`, `Genes.gff`→`genes.gff`,
  and the `Gene_trees/`, `Intergene_trees/`, `Gene_alignments/`, `Architecture/`, `Genomes/`, `BED/`
  directories. `species_tree.log`→`species.log` (matching the `<command>.log` convention);
  `Branch_events.tsv`→`branch_events.tsv` (kept "branch" — it is a genuine per-tree-edge table); and
  the plain-ILS single-file output `gene_trees.nwk`→`ils_gene_trees.nwk` so it no longer shadows the
  `gene_trees/` directory. **Breaking change with no compatibility shim** — update downstream scripts;
  file *contents* are byte-identical, only the names changed. (`species_tree.nwk` / `species_nodes.tsv`
  were already lowercase.)
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

- **`SharedBirthDeath` → `BirthDeath(per="shared")`**, and the CLI token **`--diversification shared`
  → `--per shared`**. The shared-clock birth–death is now an *opportunity* setting on `BirthDeath`, so
  the standalone name is redundant. The old spellings still work but warn (the class on construction,
  the flag on use) and `SharedBirthDeath` has left `zombi2.__all__`; both are removed in **0.4.0**.
- **`PerCopyRates` / `PerLineageRates` → `Rates(per="copy"|"lineage")`.** The gene D/T/L opportunity is
  now a knob on `Rates`, so the two named classes are redundant presets: they still work but warn on
  construction and have left `zombi2.__all__` (along with `SharedRates`, the older alias for
  `PerCopyRates`). The `--rate-per genome` token remains the deprecated spelling of `--rate-per lineage`.
  All removed in **0.4.0**.
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
