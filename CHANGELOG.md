# Changelog

All notable changes to ZOMBI2 are recorded here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/). While the project is pre-1.0, a **minor**
bump (0.x.0) carries new features or breaking changes and a **patch** bump (0.x.y) carries fixes.

A release is cut with `scripts/release.sh patch|minor|major` (the version is computed, not typed),
which moves the entries below from `[Unreleased]` into a dated version section.

## [Unreleased]

### Added
- **Cross-level staleness guard.** A level refuses to re-run in place when a later level built from it
  is already in the run directory — re-running would leave that downstream output silently mismatched.
  This covers both the pipeline chain (re-running `genomes` would orphan the `sequences` under it) and
  `DrivenBy` **conditioning** (re-running a `traits` run that a `genomes` rate was conditioned on orphans
  that genomes, and the sequences beneath it — the dependency is recorded in a small `conditioned_on`
  marker). `--force` re-runs anyway and removes the now-stale downstream, so a run's levels can never
  quietly disagree. The forward pipeline is unaffected (each level is run once); applies to the default
  grouped layout — `--flat` commingles the levels and is left to the user. (#243)
- **Onboarding nudges.** After a level runs, a one-line hint points at the key output file and the next
  command in the pipeline (e.g. `next: zombi2 genomes out/`) — suppressed by `--quiet`, so a scripted
  batch stays quiet. The top-level `--help` now leads with the levels and the plain quickstart, with the
  `DrivenBy` / joint coupling note moved below them. (#244)

### Fixed
- A non-nucleotide `sequences` run no longer leaves an empty `sequences/genomes/` directory behind (the
  assembled-genome FASTAs are a nucleotide-run output). (#244)

## [0.6.0] - 2026-07-26

### Added
- `--quiet` on `zombi2 tools format`, matching the level commands — it suppresses the summary line,
  for a log file or a batch of runs. (#242)

### Changed
- **Gene-copy ids are now written `g<id>` in every genome table** — `genomes.tsv`, `genome_events.tsv`,
  `gene_order.tsv`, `blocks.tsv` and `initial_genome.tsv`, in both the `copy` and `parent` columns.
  This is the same token the gene-tree Newick leaves, the alignment FASTA headers and the homology
  tables already use, so a gene copy now joins across every file of a run with no translation step
  (previously the tables held a bare integer while the trees/alignments held `g<id>`). The column
  names are unchanged; this is a breaking change to the *values* those columns hold. (#242)
- **A joint (`zombi2 joint`) driver mapping that names a state outside the trait's declared alphabet is
  now refused.** A key that can never occur — a typo such as `{"caev": 4.0}` for a cave/surface trait —
  is caught up front instead of silently applying to nothing. The conditioned (file-driven) path is
  unchanged: it still refuses only a mapping that matches *none* of the driver's observed states. (#242)

### Fixed
- The run-completion line no longer doubles a trailing slash: a directory given as `out/` (as the
  quickstart shows) now echoes as `wrote out/`, not `wrote out//`. (#242)
- A conditioned rate whose `DrivenBy` points at a missing file now reports it as a missing driver file
  (with the path and how to fix it), instead of a bare `[Errno 2] No such file or directory`. (#242)
- `parallel=` no longer crashes with a raw `BrokenProcessPool` when called from a notebook, `python -c`,
  or a stdin heredoc (where worker processes cannot re-import your program). It falls back to
  single-process with a one-line note; a `.py` script or the CLI still uses every core. (#242)

## [0.5.0] - 2026-07-24

### Added
- **Opt-in parallelism** for the unordered-genome and sequence engines (`parallel=` in Python,
  `--parallel` on the CLI): a separate, worker-count-invariant engine that evolves independent units
  (gene families, gene trees) across processes. It covers D/T/L/O with `uniform`/`distance` transfer;
  a driven or clade-based `transfer_to` falls back to the serial engine (loudly). Also `stream_to=`
  to write each family straight to disk for the many-families regime. The serial path is the default
  and is unchanged. (#233)
- **`Clades` transfer_to rule + the `Between` kernel** — weight horizontal transfers by the donor's
  and the recipient's **named clade**, so transfer can be steered to run *between* two clades rather
  than within them (a topological rule, sibling of `"distance"`). Name a clade by a few of its tips or
  a node id. The trait-driven form is `mod.DrivenBy(trait, Between({...}))`. (#235)

### Changed
- **Faster nucleotide sequence runs at lower memory**: a run-wide CDF cache and lazy genome assembly,
  plus batched decoding and a streaming FASTA writer. (#234)

## [0.4.0] - 2026-07-24

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
