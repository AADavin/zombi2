# ZOMBI2 — working notes for Claude

## Read this first

**Before changing any code, docs, CLI, or manual chapter, read both:**
- **[`docs/design/SPEC.md`](docs/design/SPEC.md)** — the **model and the words** (the four levels, the
  independent / conditioned / joint framework, the rate grammar, the canonical vocabulary, figure conventions).
- **[`docs/design/MAP.md`](docs/design/MAP.md)** — the **files and the names** (every module, every public
  name, its one canonical home; the clean core vs the `legacy/` quarantine; the rebuild order).

ZOMBI2 is being rebuilt as a **clean core grown from SPEC**, not a migration. The old codebase is being
**quarantined** in `legacy/` at the repo root (read-only, not importable); features are ported out of it
deliberately, one level at a time, renamed to MAP. The active tree is just the clean core. When code
disagrees with SPEC or MAP, the code is the **fossil** and must be aligned. Principle: **concepts → code →
chapter** (fix the code before documenting new behaviour).

Do not reintroduce the old lexicon. If a convention genuinely needs to change, change SPEC.md (words) or
MAP.md (shape) **first**, then propagate.

## Output files — keep the cheatsheet in sync

Every file ZOMBI2 writes is catalogued in one table per level in
[`manual/revision/appendix-b-output-files.md`](manual/revision/appendix-b-output-files.md) (Appendix B;
the docs site single-sources it at `docs/reference/output-files.md`, so editing the appendix updates both
book and site). **Whenever you add or change an output** — a new `.write()` token, a new filename or
format, a changed default, new columns, or a new Python-only accessor like `.gene_trees` — **update that
table in the same change.** The table's columns are Output · File · Format · Default (yes / no /
Python) · Contents.

## Project

ZOMBI2 is a phylogenetic simulator (Python library + CLI) that simulates four levels of evolution —
species trees, genomes, sequences, traits — independently, conditioned, or jointly, and records the
true history behind every dataset. Author: Adrián Davín.

The clean core is **pure Python**: nothing in `zombi2/` imports the Rust extension. A Rust engine is
stage 2 of the rebuild; today only `legacy/` binds to it.

## Run environment

- Python: `/Users/aadria/miniconda3/bin/python` (3.12). Bare `python`/`zombi2` are not on PATH.
- The clean core needs **no build step**. If `rust/` ever gets wired into `zombi2/`, rebuild after any
  `rust/` change: `maturin build --release -m rust/Cargo.toml` then force-reinstall.
- PDF/manual toolchain: `xelatex` at `/Library/TeX/texbin`, `pandoc` + `rsvg-convert` on PATH.

## Manual revision (in progress)

A chapter-by-chapter rewrite of the manual is underway with Adrián, tracked in `manual/revision/`:
- `dashboard.yml` — every round-1 comment + the ratified decisions (source of truth for the revision).
- `proposed-index.md` — the agreed 9-chapter index with section detail.
- `SPEC.md` (in `docs/design/`) is the distilled constitution; the manual is its exposition.

Work happens on a branch in an isolated worktree, not the shared main checkout.
