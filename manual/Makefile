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

# The ZOMBI2 release this manual documents, shown on the title page. Read from the package rather
# than typed here, so the book can never claim a version the code does not have. Empty means the
# read failed (a moved file, a changed line) — stop rather than print "Version  —".
VERSION := $(shell sed -n 's/^__version__ = "\(.*\)"/\1/p' ../zombi2/__init__.py)
ifeq ($(strip $(VERSION)),)
$(error could not read __version__ from ../zombi2/__init__.py)
endif

# The book, in reading order. Spelled out rather than wildcarded: a wildcard sorts
# appendix-a/appendix-b in FRONT of ch1, and — worse — silently expands to nothing if the
# directory is ever renamed, which builds an empty book instead of failing. With the list
# explicit, a missing or renamed file stops the build.
BOOK      := book
CHAP_MDS  := ch1.md ch2.md ch3.md ch4.md ch5.md ch6.md ch7.md ch8.md ch9.md ch10.md \
             appendix-a.md appendix-b.md appendix-c.md
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
          --metadata-file=metadata.yaml --metadata=date:"Version $(VERSION) — $(DATE)"

DEPS := preamble.tex callouts.lua metadata.yaml $(BIB) $(FIGS) ../zombi2/__init__.py

.PHONY: all manual chapters figures docx clean
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

# --- the book as .docx, for reading and commenting ---
# Word cannot show a PDF image, and the chapters name one of two artefacts: `figures/NAME.pdf` for a
# generated figure and `figures/NAME_print.png` for a hand-authored chapter diagram. So the sources
# are staged with every image as PNG and the .pdf references rewritten, and pandoc is pointed at the
# staging directory. No `--pdf-engine`, no LaTeX filter: a docx takes neither.
docx: build/zombi2-manual.docx
build/zombi2-manual.docx: $(CHAPTERS) $(BIB) $(BOOK_FIGS) $(FIGS)
	@mkdir -p build/docx/figures
	@cp ../figures/png/*.png build/docx/figures/ 2>/dev/null || true
	@cp $(BOOK_FIGDIR)/*_print.png build/docx/figures/ 2>/dev/null || true
	@for f in $(CHAP_MDS); do \
		sed -e 's|(figures/\([A-Za-z0-9_-]*\)\.pdf)|(figures/\1.png)|g' $(BOOK)/$$f > build/docx/$$f; \
	done
	$(PANDOC) --from=markdown --number-sections --citeproc --bibliography=$(BIB) \
	          --resource-path=build/docx --metadata-file=metadata.yaml \
	          --metadata=date:"Version $(VERSION) — $(DATE)" --toc --toc-depth=2 \
	          $(addprefix build/docx/,$(CHAP_MDS)) -o $@
	@$(PYTHON) fix_docx_images.py $@

PYTHON ?= python3

clean:
	rm -rf build figures
	rm -f $(BOOK_FIGDIR)/*_print.svg $(BOOK_FIGDIR)/*_print.png
