# Changelog

All notable changes to ZOMBI2 are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and this
project aims to follow [Semantic Versioning](https://semver.org/spec/v2.0.0.html). Until the
first tagged release the changelog tracks the development line under **Unreleased**.

## [Unreleased]

First public development line (`0.2.0.dev0`). ZOMBI2 is a ground-up redesign of
[ZOMBI](https://github.com/AADavin/Zombi).

### Added

- **Species trees** — backward (reconstructed) and forward (complete) birth–death and Yule,
  episodic/skyline rate shifts, fossilized birth–death and incomplete sampling, ClaDS /
  diversity-dependent / clade-shift diversification, mass extinctions, and ghost lineages.
  Backward trees scale to millions of tips on a laptop via a native Rust engine.
- **Gene families** — duplication, transfer, loss and origination along a species tree, with
  uniform, family-sampled and genome-wise rate models; ordered chromosomes with inversions and
  transpositions; nucleotide-resolution genomes; and a gene-family coupling model
  (non-independence). Output as full event logs, compact event traces, or counts-only sparse
  profiles, with parallel replicates.
- **Traits** — Brownian motion, Ornstein–Uhlenbeck, early burst, Mk, threshold,
  correlated-binary, hidden-state and multi-optimum models, Pagel transforms, and DEC
  biogeography, evolved over a phylogeny.
- **Sequences** — a gene × lineage relaxed-clock family (strict, UCLN, UGAM, white-noise,
  autocorrelated-lognormal, CIR) that rescales gene trees from time into substitutions/site,
  plus nucleotide substitution models (JC, K80, HKY, GTR + Gamma).
- **Coevolution** — a unified `coevolve` command that couples species, traits and genes along
  six directed edges via `--couple driver:target`.
- **Command-line interface** — `zombi2 species | genomes | trait | sequence | coevolve`, with
  grouped, sectioned help and a per-run reproducibility manifest (version, seed, full command
  line and resolved parameters).
- **Packaging** — GPL-3.0-or-later license, `CITATION.cff`, and PyPI metadata.

### Notes

- ABC inference of DTL rates from an empirical copy-number profile (`zombi2.matching`) is
  available as an experimental Python API but is withheld from the command line in this release.
- The default (built-in) gene-family engine is a compiled Rust extension; build it once with
  maturin (see [installation](docs/installation.md)).
