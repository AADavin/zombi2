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
[`manual/book/appendix-b.md`](manual/book/appendix-b.md) (Appendix B;
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

## The manual

The book lives in `manual/book/`, one file per chapter: `ch1.md` … `ch9.md`, plus `appendix-a.md`
(the Gillespie algorithm) and `appendix-b.md` (output files). SPEC is the constitution; the manual is
its exposition.

Ch4–Ch6 are the genome **resolution** ladder — unordered ⊂ ordered ⊂ nucleotide — one chapter per rung.
Every level chapter closes the same way: **The objects → Usage from Python → Usage from the CLI →
Outputs.** Ch1, Ch2 and Ch9 are essays and are exempt.

Two rules learned the hard way:
- **Run every example before trusting a chapter.** They drift behind the code in both directions —
  flags get invented, and shipped features get described as forthcoming.
- **Document only what ships.** Cut an unbuilt model rather than marking it "not yet". Where a level is
  genuinely partial, say so plainly (Ch6 does this: the nucleotide resolution has no CLI and no
  `.write()`, so it documents Python accessors only).

Only Ch3 and Ch4 are published to the docs site so far, via snippet includes in `docs/guide/`; Appendix
B is included at `docs/reference/output-files.md`. Renaming a chapter file breaks those includes and
CI's `mkdocs --strict` will catch it.

Work happens on a branch in an isolated worktree, not the shared main checkout.
