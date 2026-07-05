# ZOMBI2 figure house style

The agreed style for every manual figure. The reference implementations are
`fig_diversity_dependent.py`, `fig_species_tree_events.py`, `fig_gene_tree.py`,
`fig_clads.py`, `fig_mass_extinction.py`, `fig_trait_ou.py`. Shared constants live in
`figures/scripts/zombi_style.py`.

## Canvas & rendering
- `drawsvg` canvas ~1000–1180 px wide. Make the canvas tall/wide enough that **no label
  is clipped** (check the bottom axis label and right-edge labels).
- Figures are converted SVG → PDF (`rsvg-convert`) and placed in the manual, where they
  are centered on the page (captioned `figure` floats). Information-dense figures render
  wider (up to full text width); simple trees can be smaller.

## Typography
- Font: `FONT` (Helvetica). Sizes come from the shared scale in `zombi_style.py`:
  `FS_TITLE=32`, `FS_LABEL=22`, `FS_ANNOT=20`, `FS_TICK=18`. Never hard-code small sizes.
- **Title**: one short bold line, **horizontally centered** at the top
  (`text_anchor="middle"`, `x = W/2`). **No subtitle.**
- **ASCII only** in text — `rsvg`/`cairosvg` render `->`... no: do NOT use `→ λ μ θ × ≤`
  etc. (they become tofu boxes). Use words ("to", "lambda", "theta") or *draw* the symbol.

## Legends
- Always **clear of the tree/plot** — never overlapping the data.
- Prefer **top-left**. Multi-entry legends are a single vertical column.

## Colour — used only where it carries information
- Default is black & white: near-black `INK` branches/curves, grey dashed for
  extinct/secondary lineages.
- **Lineage-through-time / skyline curves are BLACK** (no accent colour).
- The one place colour is used: a **continuous quantity painted along branches** →
  **viridis** gradient (e.g. the OU trait). ClaDS per-lineage rates stay grayscale.

## Symbol conventions (gene-family event figures)
- Event markers are **solid black**: filled square = duplication, filled triangle (or a
  black arc + arrowhead on the species tree) = transfer, **cross (×, two crossing black
  segments) = loss**. Do not colour event markers.
- A small symbol→event legend, as a top-left column, in the order **Duplication,
  Transfer, Loss**.

## Axes
- A visible time axis: baseline line + tick marks + numeric labels + a
  "time (root to present)" label, consistent across figures.
