# Changelog

All notable changes to ZOMBI2 are recorded here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/). While the project is pre-1.0, a **minor**
bump (0.x.0) carries new features or breaking changes and a **patch** bump (0.x.y) carries fixes.

A release is cut with `scripts/release.sh patch|minor|major` (the version is computed, not typed),
which moves the entries below from `[Unreleased]` into a dated version section.

## [Unreleased]

### Added
- `analyses/` — a home for self-contained validation studies, each regenerating from fixed seeds:
  **RED** node-age validation (does Relative Evolutionary Divergence recover node ages under a
  realistic clock?) and the **yeast inversion-rate** study (inferring the genome inversion rate from
  synteny with the nucleotide genome model). (#227)
- `zombi2.tree` — one home for the `Tree` datatype and its toolkit, all free functions: `with_stem`,
  `make_ultrametric`, `rescale`, `relative_evolutionary_divergence`, `red_scaled`, `distance`
  (Robinson–Foulds / branch-score), and `read_newick(assume_extant=)`. (#228)
- CLI: `zombi2 tools tree` (prune / round / stem / rescale / RED) and `zombi2 tools treedist`. (#228)
- The **autocorrelated** molecular clock (`FromParent`) at the sequence level, alongside the
  uncorrelated `ByLineage`. (#228)
- Release tooling: this `CHANGELOG.md` and `scripts/release.sh`.

### Changed
- `Tree`, `Node`, `read_newick`, `to_newick`, `prune` moved from `zombi2.species` to `zombi2.tree`
  — import them from `zombi2.tree` now (`zombi2.species` keeps the simulator and `Event`). (#228)
- `docs/design/MAP.md` correctness pass: the quarantine now points at the sibling `ZOMBI2_LEGACY/`,
  the `joint` CLI is listed, sequences/traits are marked built, and stale signatures/sections fixed.
  (#228)

### Fixed
- `zombi2 tools treedist` matches tips by label, not by parse-order node ids, so two external trees
  are compared by taxon rather than by structural position. (#228)

## [0.3.0] - 2026-07-23

### Added
- The clean, pure-Python core grown from `docs/design/SPEC.md`: species trees, genomes at three
  resolutions (unordered ⊂ ordered ⊂ nucleotide), sequences, traits, and the joint engine — with the
  cross-level rate grammar (`zombi2.rates`) and a CLI. First release of the rewrite; `pip install`
  needs no build step.

### Removed
- The Rust engine (`zombi2_core`) and the old codebase, retired to the sibling `ZOMBI2_LEGACY/`.
  (0.2.0 was yanked on PyPI, superseded by this clean core.)
