# ZOMBI2 manual — Markdown chapters -> per-chapter PDFs and one merged book.
#
#   make figures                  convert docs/img/*.svg -> figures/*.pdf
#   make build/01-introduction.pdf   build ONE chapter (fast preview, isolated)
#   make chapters                 build every chapter as its own PDF
#   make manual                   build the merged book, build/zombi2-manual.pdf
#   make clean
#
# Each chapter is an independent source file; the book is assembled by a single
# pandoc pass over the chapter files ON DISK, so nothing needs the whole manual in memory.

PANDOC ?= pandoc
ENGINE ?= xelatex
BIB    := ../docs/references.bib
IMG    := ../docs/img
DATE   := $(shell date '+%B %Y')

CHAPTERS  := $(sort $(wildcard chapters/*.md))
CHAP_PDFS := $(patsubst chapters/%.md,build/%.pdf,$(CHAPTERS))
SVGS      := $(wildcard $(IMG)/*.svg)
FIGS      := $(patsubst $(IMG)/%.svg,figures/%.pdf,$(SVGS))

COMMON := --from=markdown --pdf-engine=$(ENGINE) --top-level-division=chapter \
          --number-sections --citeproc --bibliography=$(BIB) \
          --lua-filter=callouts.lua --include-in-header=preamble.tex \
          --syntax-highlighting=tango --resource-path=.:figures \
          --metadata-file=metadata.yaml --metadata=date:"$(DATE)"

DEPS := preamble.tex callouts.lua metadata.yaml $(BIB) $(FIGS)

.PHONY: all manual chapters figures clean
all: manual

# --- figures: SVG -> PDF (vector) ---
figures: $(FIGS)
figures/%.pdf: $(IMG)/%.svg
	@mkdir -p figures
	rsvg-convert -f pdf -o $@ $<

# --- one PDF per chapter (isolated preview) ---
chapters: $(CHAP_PDFS)
build/%.pdf: chapters/%.md $(DEPS)
	@mkdir -p build
	$(PANDOC) $(COMMON) --toc --toc-depth=1 $< -o $@

# --- the merged book (single pass over all chapters) ---
manual: build/zombi2-manual.pdf
build/zombi2-manual.pdf: $(CHAPTERS) $(DEPS)
	@mkdir -p build
	$(PANDOC) $(COMMON) --toc --toc-depth=2 $(CHAPTERS) -o $@

clean:
	rm -rf build figures
