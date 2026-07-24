# ZOMBI2 — working notes for Claude

## Read this first

**Before changing any code, docs, CLI, or manual chapter, read:**

- **[`docs/design/SPEC.md`](docs/design/SPEC.md)** — the **model and the words** (the four levels, the
  independent / conditioned / joint framework, the rate grammar, the canonical vocabulary).

**Any change to the code must stay consistent with what is already written** — with the structure and
decisions derivable from the code and the manual themselves. When in doubt, defer to Adrián A. Davín,
who has the last word.

## Output files — keep the cheatsheet in sync

Every file ZOMBI2 writes is catalogued in one table per level in
[`manual/book/appendix-b.md`](manual/book/appendix-b.md) (Appendix B;
the docs site single-sources it at `docs/reference/output-files.md`, so editing the appendix updates both
book and site). **Whenever you add or change an output** — a new `.write()` token, a new filename or
format, a changed default, new columns, or a new Python-only accessor like `.gene_trees` — **update that
table in the same change.** The table's columns are Output · File · Format · Default (yes / no /
Python) · Contents.

## Run environment

- Python: `/Users/aadria/miniconda3/bin/python` (3.12). Bare `python`/`zombi2` are not on PATH.
- The core needs **no build step** (pure Python).
- PDF/manual toolchain: `xelatex` at `/Library/TeX/texbin`, `pandoc` + `rsvg-convert` on PATH.

## Analyses

`analyses/<study>/` holds **self-contained validation studies** built on the shipped API — each its own
scripts, data, figures, and `REPORT.md`, regenerating deterministically from fixed seeds (see
[`analyses/README.md`](analyses/README.md)). 

## The manual

The book lives in `manual/book/`, one file per chapter: `ch1.md` … `ch9.md`, plus `appendix-a.md`
(the Gillespie algorithm), `appendix-b.md` (output files) and `appendix-c.md` (the `zombi2 tools`
read-back commands). 

Two rules learned the hard way:
- **Run every example before trusting a chapter.** 
- **Document only what ships.** Cut an unbuilt model rather than marking it "not yet". Where a level is
  genuinely partial, say so plainly, and delete the caveat the moment it stops being true.

Every chapter is published to the docs site as a snippet include — Ch1–Ch9 under `docs/guide/`, the
appendices under `docs/reference/`. Renaming a chapter file breaks those includes and CI's
`mkdocs --strict` will catch it.
