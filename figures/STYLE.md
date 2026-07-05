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
- The main place colour is used: a **continuous quantity painted along branches** →
  **viridis** gradient (e.g. the OU trait). ClaDS per-lineage rates stay grayscale.
- **Categorical exception (Adrián's call):** where a small set of *categorical identities* must
  be told apart and grey genuinely fails, use a distinct colour-blind-safe hue per category
  (Paul Tol 'bright': `#4477AA` / `#EE6677` / `#228833`). Applied to the **DEC** areas (and the
  gain=green / loss=red dispersal-extinction arrows). Do NOT do this for binary states (open vs
  filled reads fine) or where 2–3 greys suffice (e.g. the Mk 3-state map).

## Symbol conventions (gene-family event figures)
- Event markers are **solid black**: filled square = duplication, filled triangle (or a
  black arc + arrowhead on the species tree) = transfer, **cross (×, two crossing black
  segments) = loss**. Do not colour event markers.
- A small symbol→event legend, as a top-left column, in the order **Duplication,
  Transfer, Loss**.

## Axes
- A visible time axis: baseline line + tick marks + numeric labels + a
  "time (root to present)" label, consistent across figures.

## Tree aesthetics (Adrián's cross-cutting rule)
Trees should look like one another across the whole manual — **wide and dense**, like the
trait trees of Chapters 12–13 (`fig_trait_ou.py` is the reference). Never tall-and-narrow.

- **Aspect ratio: landscape, ~1.5–1.7:1** (wider than tall). If a `species_style` height grows
  with the tip count, **cap it** (`height=min(820, k·n_leaves + c)`) so more tips make the tree
  *denser*, not *taller*. A portrait tree is a bug — widen it.
- **Canvas width ~1150–1300 px** for a standalone tree; height then follows from the aspect.
- **Target tip count ~20–30** for an illustrative tree (the old ~10-tip trees read as sparse).
  Simulate more species / pick a seed with more leaves rather than settling for a small tree.
- **Tip labels are optional and usually dropped** on a dense tree — a coloured tip square or a
  bare leaf reads better than 25 single letters. Label only the tips a caption refers to (e.g.
  the extinct `e1…` lineages), and leave the rest unlabelled.
- **Fonts** stay on the shared `FS_*` scale (`FS_TITLE/FS_LABEL/FS_TICK`); `species_style(...,
  font_size=FS_TICK)` for the leaf/axis text. A multi-panel figure that composes several trees
  side by side is *very* wide once on the page, so its fonts must be set correspondingly larger
  (see `fig_model_ghosts.py`).
