# ZOMBI2 manual — the Markdown chapters in book/ -> per-chapter PDFs and one merged book.
#
#   make figures                  convert ../figures/svg/*.svg -> figures/*.pdf
#   make build/ch1.pdf            build ONE chapter (fast preview, isolated)
#   make chapters                 build every chapter as its own PDF
#   make manual                   build the merged book, build/zombi2-manual.pdf
#   make clean
#
# Each chapter is an independent source file; the book is assembled by a single
# pandoc pass over the chapter files ON DISK, so nothing needs the whole manual in memory.

PANDOC ?= pandoc
ENGINE ?= xelatex
BIB    := ../docs/references.bib
IMG    := ../figures/svg
DATE   := $(shell date '+%B %Y')

# The book, in reading order. Spelled out rather than wildcarded: a wildcard sorts
# appendix-a/appendix-b in FRONT of ch1, and — worse — silently expands to nothing if the
# directory is ever renamed, which builds an empty book instead of failing. With the list
# explicit, a missing or renamed file stops the build.
BOOK      := book
CHAP_MDS  := ch1.md ch2.md ch3.md ch4.md ch5.md ch6.md ch7.md ch8.md ch9.md \
             appendix-a.md appendix-b.md
CHAPTERS  := $(addprefix $(BOOK)/,$(CHAP_MDS))
CHAP_PDFS := $(patsubst $(BOOK)/%.md,build/%.pdf,$(CHAPTERS))
SVGS      := $(wildcard $(IMG)/*.svg)
FIGS      := $(patsubst $(IMG)/%.svg,figures/%.pdf,$(SVGS))

# Chapter figures live beside the chapters, in book/figures/, and are authored for the web:
# they paint with var(--ink)/var(--paper) and carry a prefers-color-scheme block. librsvg
# resolves neither, so converting one straight gives black-on-black silhouettes — the custom
# properties have to be flattened to their light-theme literals first. Which figures to build
# is read out of the chapters, so adding one needs no edit here.
BOOK_FIGDIR := $(BOOK)/figures
BOOK_FIGS   := $(sort $(addprefix $(BOOK_FIGDIR)/,$(notdir \
                 $(shell grep -ohE 'figures/[A-Za-z0-9_-]+_print\.png' $(CHAPTERS)))))

COMMON := --from=markdown --pdf-engine=$(ENGINE) --top-level-division=chapter \
          --number-sections --citeproc --bibliography=$(BIB) \
          --lua-filter=callouts.lua --include-in-header=preamble.tex \
          --syntax-highlighting=tango --resource-path=.:figures:$(BOOK) \
          --metadata-file=metadata.yaml --metadata=date:"$(DATE)"

DEPS := preamble.tex callouts.lua metadata.yaml $(BIB) $(FIGS)

.PHONY: all manual chapters figures clean
all: manual

# --- figures: SVG -> PDF (vector) ---
figures: $(FIGS)
figures/%.pdf: $(IMG)/%.svg
	@mkdir -p figures
	rsvg-convert -f pdf -o $@ $<

# --- chapter figures: web SVG -> flattened SVG -> print PNG ---
$(BOOK_FIGDIR)/%_print.svg: $(BOOK_FIGDIR)/%.svg
	sed -e 's/var(--ink)/#1a1a1a/g' -e 's/var(--paper)/#ffffff/g' $< > $@
$(BOOK_FIGDIR)/%_print.png: $(BOOK_FIGDIR)/%_print.svg
	rsvg-convert -f png -z 3 -o $@ $<

# --- one PDF per chapter (isolated preview) ---
chapters: $(CHAP_PDFS)
build/%.pdf: $(BOOK)/%.md $(DEPS) $(BOOK_FIGS)
	@mkdir -p build
	$(PANDOC) $(COMMON) --toc --toc-depth=1 $< -o $@

# --- the merged book (single pass over all chapters) ---
manual: build/zombi2-manual.pdf
build/zombi2-manual.pdf: $(CHAPTERS) $(DEPS) $(BOOK_FIGS)
	@mkdir -p build
	$(PANDOC) $(COMMON) --toc --toc-depth=2 $(CHAPTERS) -o $@

clean:
	rm -rf build figures
	rm -f $(BOOK_FIGDIR)/*_print.svg $(BOOK_FIGDIR)/*_print.png
