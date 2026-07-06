# The ZOMBI2 manual (build system)

The PDF manual is authored as **one Markdown file per chapter** in `chapters/`, built with
pandoc → XeLaTeX. Each chapter compiles to its own PDF for isolated preview, and a single
`make manual` assembles the book. The assembly is a mechanical pandoc pass over the chapter
files **on disk** — so you can write or edit *one* chapter without ever loading the whole
manual into context. That is the whole point of this layout.

## Layout

```
manual/
  chapters/NN-title.md   the chapters (order = numeric filename prefix)
  metadata.yaml          title, author, document class, link colours
  preamble.tex           LaTeX preamble (callout boxes, figure defaults, captions)
  callouts.lua           maps ::: note / ::: warning / ::: tip to boxes
  Makefile               figures + per-chapter PDFs + the merged book
  figures/               SVG→PDF, auto-built from ../docs/img (git-ignored)
  build/                 output PDFs (git-ignored)
```

Bibliography: `../docs/references.bib` (the verified database). Figures: the SVGs in
`../docs/img/` (converted to PDF automatically).

## Build

```bash
cd manual
make build/01-introduction.pdf   # one chapter (fast; what you use while writing)
make chapters                    # every chapter as its own PDF
make manual                      # the merged book -> build/zombi2-manual.pdf
make clean
```

## Writing a chapter (for a future editor)

To work on a chapter you only need: **that one chapter file**, this README, and
`../docs/references.md` (for citation keys). You never need the other chapters.

- **Headings.** Start the file with `# Chapter Title` (one level-1 heading → one `\chapter`);
  use `##`/`###` for sections.
- **Figures.** `![Caption text.](figures/NAME.pdf)` where `NAME.svg` exists in `../docs/img/`
  (the Makefile converts it). The alt text becomes the printed caption. Add new figures by
  dropping a colour SVG into `../docs/img/` (see `../docs/guide` conventions).
- **Citations.** `[@bibkey]`, e.g. `[@davin2020zombi]`. Keys are listed in
  `../docs/references.md` / `../docs/references.bib`. A References section is generated
  automatically (per chapter in previews; once, at the end, in the book).
- **Callouts.** Fenced divs:
  ```
  ::: note
  A short aside.
  :::
  ```
  Classes: `note`, `warning`, `tip`.
- **Code.** Fenced code blocks with a language tag (```` ```python ````) are syntax-highlighted.
- **Cross-references.** Link to another chapter/section by its heading id, e.g.
  `[unordered genomes](#unordered-genomes)` (resolves in the merged book).

## Scope

The manual is **Concepts + Tutorial** (Parts I–VI). The CLI and Python API reference stay in
the online docs. Planned parts: I Getting started · II Species trees · III Gene families ·
IV Traits · V Coevolution · VI Sequence evolution. Add chapters incrementally as
`chapters/NN-title.md`.

## Requirements

pandoc, a LaTeX engine (XeLaTeX / TeX Live), and `rsvg-convert` (or edit the Makefile to use
`cairosvg`). All were present at setup time.
