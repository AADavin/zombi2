# ZOMBI2 manual — Markdown chapters -> per-chapter PDFs and one merged book.
#
#   make figures                  convert docs/img/*.svg -> figures/*.pdf
#   make build/01-introduction.pdf   build ONE chapter (fast preview, isolated)
#   make chapters                 build every chapter as its own PDF
#   make manual                   build the merged book, build/zombi2-manual.pdf
#   make tools-manual             build the tools-only book, build/zombi2-tools-manual.pdf
#   make supplements              build standalone supplements (NOT part of the book)
#   make clean
#
# Each chapter is an independent source file; the book is assembled by a single
# pandoc pass over the chapter files ON DISK, so nothing needs the whole manual in memory.

PANDOC ?= pandoc
PYTHON ?= python3
ENGINE ?= xelatex
BIB    := ../docs/references.bib
IMG    := ../docs/img
DATE   := $(shell date '+%B %Y')

CHAPTERS  := $(sort $(wildcard chapters/*.md))
CHAP_PDFS := $(patsubst chapters/%.md,build/%.pdf,$(CHAPTERS))
# standalone supplements: independent chapters, deliberately NOT part of the merged book
SUPP      := $(sort $(wildcard supplements/*.md))
SUPP_PDFS := $(patsubst supplements/%.md,build/supp-%.pdf,$(SUPP))
SVGS      := $(wildcard $(IMG)/*.svg)
FIGS      := $(patsubst $(IMG)/%.svg,figures/%.pdf,$(SVGS))

# The tools book is GENERATED from docs/tools/ (the single source of truth) by
# tools_to_chapters.py, so there is no hand-maintained copy to drift. Listed in
# nav order (see mkdocs.yml); the overview page leads.
TOOLS_ORDER := index reconciliation-likelihood undated-simulation red treedist \
               recon-accuracy reconciliation-parsing geneorder-export
TOOLS_SRC   := $(addprefix ../docs/tools/,$(addsuffix .md,$(TOOLS_ORDER)))
# only the figures those pages actually reference (not the whole docs/img set)
TOOLS_FIGS  := $(addprefix figures/,$(addsuffix .pdf,$(basename \
               $(notdir $(shell grep -ohE '\.\./img/[A-Za-z0-9_-]+\.svg' $(TOOLS_SRC))))))

COMMON := --from=markdown --pdf-engine=$(ENGINE) --top-level-division=chapter \
          --number-sections --citeproc --bibliography=$(BIB) \
          --lua-filter=callouts.lua --include-in-header=preamble.tex \
          --syntax-highlighting=tango --resource-path=.:figures \
          --metadata-file=metadata.yaml --metadata=date:"$(DATE)"

DEPS_BASE := preamble.tex callouts.lua metadata.yaml $(BIB)
DEPS      := $(DEPS_BASE) $(FIGS)

.PHONY: all manual tools-manual chapters supplements figures clean
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

# --- the tools-only book (generated from docs/tools/, then one pandoc pass) ---
tools-manual: build/zombi2-tools-manual.pdf
build/zombi2-tools-manual.pdf: $(TOOLS_SRC) tools_to_chapters.py preamble-tools.tex $(DEPS_BASE) $(TOOLS_FIGS)
	@mkdir -p build/tools-chapters
	$(PYTHON) tools_to_chapters.py --out build/tools-chapters $(TOOLS_SRC)
	$(PANDOC) $(COMMON) --include-in-header=preamble-tools.tex --toc --toc-depth=2 \
	  --metadata=title:"ZOMBI2 Tools" \
	  --metadata=subtitle:"The zombi2.tools analysis & interoperability layer" \
	  build/tools-chapters/*.md -o $@

# --- standalone supplements (independent PDFs, each with its own title page) ---
supplements: $(SUPP_PDFS)
build/supp-%.pdf: supplements/%.md $(DEPS)
	@mkdir -p build
	$(PANDOC) $(COMMON) --toc --toc-depth=2 \
	  --metadata=title:"ZOMBI2 — Supplement" \
	  --metadata=subtitle:"An experimental supplement (zombi2.experimental)" $< -o $@

# per-supplement title (an explicit rule overrides the pattern rule above for this one target)
build/supp-incomplete-lineage-sorting.pdf: supplements/incomplete-lineage-sorting.md $(DEPS)
	@mkdir -p build
	$(PANDOC) $(COMMON) --toc --toc-depth=2 \
	  --metadata=title:"ZOMBI2 — Incomplete Lineage Sorting" \
	  --metadata=subtitle:"An experimental supplement (zombi2.experimental)" $< -o $@

clean:
	rm -rf build figures
