# The ZOMBI2 manual (build system)

The PDF manual is authored as **one Markdown file per chapter** in `book/`, built with
pandoc → XeLaTeX. Each chapter compiles to its own PDF for isolated preview, and a single
`make manual` assembles the book. The assembly is a mechanical pandoc pass over the chapter
files **on disk** — so you can write or edit *one* chapter without ever loading the whole
manual into context. That is the whole point of this layout.

## Layout

```
manual/
  book/ch1.md … ch8.md   the chapters
  book/appendix-a.md     the appendices
  book/appendix-b.md
  book/appendix-c.md
  book/appendix-d.md
  book/figures/          chapter figures, as SVGs (see Figures below)
  metadata.yaml          title, author, document class, link colours
  preamble.tex           LaTeX preamble (callout boxes, glyphs, figures, captions)
  callouts.lua           maps ::: note / ::: warning / ::: tip to boxes
  Makefile               figures + per-chapter PDFs + the merged book
  figures/               SVG→PDF, auto-built from ../figures/svg/ (git-ignored)
  build/                 output PDFs (git-ignored)
```

The chapter order is an explicit list, not a wildcard: a wildcard sorts the appendices in
front of `ch1`, and expands to nothing — building an empty book instead of failing — if the
directory is ever renamed. `book/README.md` maps each file to its chapter title.

Bibliography: `../docs/references.bib` (the verified database).

## Build

```bash
cd manual
make build/ch1.pdf               # one chapter (fast; what you use while writing)
make chapters                    # every chapter as its own PDF
make manual                      # the merged book -> build/zombi2-manual.pdf
make clean
```

## Writing a chapter (for a future editor)

To work on a chapter you only need: **that one chapter file**, this README, and
`../docs/references.bib` (for citation keys). You never need the other chapters.

- **Headings.** Start the file with `# Chapter Title` (one level-1 heading → one `\chapter`);
  use `##`/`###` for sections.
- **Figures.** Drop the SVG in `book/figures/` and write
  `![Caption text.](figures/NAME_print.png)`. The alt text becomes the printed caption. The
  `_print` file is **generated**: the chapter SVGs paint with literal colours, so the Makefile
  only rasterises them to PNG. The site's figures are also reachable: any `NAME.svg` in
  `../figures/svg/` is auto-converted, and a chapter reaches it as
  `![Caption.](figures/NAME.pdf)`.
- **Citations.** `[@bibkey]`, e.g. `[@davin2020zombi]`. Keys are listed in
  `../docs/references.bib`. A References section is generated
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
  `[gene families](#genomes-i-gene-families)` (resolves in the merged book).

## The tools manual (retired)

There was a second book, `zombi2-tools-manual.pdf`, generated from `../docs/tools/` by
`tools_to_chapters.py`. Neither the script nor the second book exists any more. The tools are
documented as Appendix D of this manual, and the pages under `../docs/tools/` are generated
from that appendix by snippet include.

## Scope

The manual is **Concepts + Tutorial**: eight chapters and four appendices. The CLI and Python
API reference stay in the online docs. The chapter order lives in two places, `book/README.md`
and the Makefile's `CHAP_MDS`: add or reorder chapters in both.

## Requirements

pandoc, a LaTeX engine (XeLaTeX / TeX Live), and `rsvg-convert` (or edit the Makefile to use
`cairosvg`). All were present at setup time.
