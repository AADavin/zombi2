# ZOMBI2 — working notes for Claude

## Read this first

**Before changing any code, docs, CLI, or manual chapter, read [`docs/design/SPEC.md`](docs/design/SPEC.md).**

It is the authoritative specification for how ZOMBI2 is organised: the four levels and their exact
words, the independent / conditioned / joint framework, which couplings are allowed, how rates are
defined, the level ontology (process / sampling / mechanics), the canonical vocabulary, and the figure
conventions. The codebase still contains **fossils** from an earlier lexicon (a `coevolve` command,
"diamond", "propensity", "opportunity", single-letter levels, dead code in `grammar.py`). When code
disagrees with the SPEC, the code is wrong and must be aligned to it — see SPEC.md §13 for the known
gaps. Principle: **concepts → code → chapter** (fix the code before documenting new behaviour).

Do not reintroduce the old lexicon. If a convention genuinely needs to change, change SPEC.md **first**,
then propagate.

## Project

ZOMBI2 is a phylogenetic simulator (Rust engine + Python library + CLI) that simulates four levels of
evolution — species trees, genomes, sequences, traits — independently, conditioned, or jointly, and
records the true history behind every dataset. Author: Adrián Davín.

## Run environment

- Python: `/Users/aadria/miniconda3/bin/python` (3.12). Bare `python`/`zombi2` are not on PATH.
- After any `rust/` change, rebuild: `maturin build --release -m rust/Cargo.toml` then force-reinstall.
- PDF/manual toolchain: `xelatex` at `/Library/TeX/texbin`, `pandoc` + `rsvg-convert` on PATH.

## Manual revision (in progress)

A chapter-by-chapter rewrite of the manual is underway with Adrián, tracked in `manual/revision/`:
- `dashboard.yml` — every round-1 comment + the ratified decisions (source of truth for the revision).
- `proposed-index.md` — the agreed 11-chapter index with section detail.
- `SPEC.md` (in `docs/design/`) is the distilled constitution; the manual is its exposition.

Work happens on a branch in an isolated worktree, not the shared main checkout.
