"""Shared 'house style' for all ZOMBI2 publication figures.

Every figure imports from here so the whole set stays visually consistent:
same font, same palette, same stroke weights, same canvas conventions. Tweak a
value here once and re-run the figure scripts to restyle the entire set.
"""

from __future__ import annotations

import phylustrator as ph

# --- Typography -----------------------------------------------------------
FONT = "Helvetica"          # falls back gracefully to Arial on systems without it

# Figure font sizes, in drawsvg user units, tuned so text stays readable once a
# ~1000-1180px-wide canvas is scaled down to roughly full text width in the
# manual (~9-11pt on the page). Bump these once to rescale text across every
# drawsvg figure that imports them. Everything below the title shares one size
# so ticks, inline annotations and legend/axis labels read consistently.
FS_TITLE = 32               # figure title (bold)
FS_LABEL = 22               # axis titles, legend entries, curve labels
FS_ANNOT = 22               # inline annotations
FS_TICK  = 22               # tick numbers and other small numeric labels

# --- Core palette ---------------------------------------------------------
# A restrained, print-friendly set. Branches are near-black (not pure #000,
# which reads harsh at small sizes); accents come from ColorBrewer Set1, the
# same family Phylustrator/ZOMBI use for gene-family events.
INK        = "#1a1a1a"      # branches, axes, primary text
MUTED      = "#8a8a8a"      # secondary text (internal-node labels, captions)
PANEL      = "#ffffff"      # background

ACCENT = {
    "origination": "#984EA3",   # purple
    "duplication": "#377EB8",   # blue
    "transfer":    "#4DAF4A",   # green
    "loss":        "#E41A1C",   # red
    "speciation":  "#999999",   # grey
    "highlight":   "#FDBF6F",   # warm sand, for clade shading
}

# --- Binary-state figures: a tasteful pale two-tone -----------------------
# Several coevolution figures (SSE, key-innovation, punctuational genome, HiSSE,
# cladogenetic) encode a BINARY lineage state as heavy-vs-light branches. Their
# default is B&W (near-black heavy / grey light — kept as a `<name>_bw.svg`), but
# a restrained COLOUR version reads more clearly on the page. The pair below is
# deliberately pale/neutral and harmonises with the viridis ramp used elsewhere:
# a muted teal for the "active" state (present / state-1 / carrier / fast-hidden)
# and a warm pale taupe for the "inactive" one. Colour-blind-safe (teal vs taupe
# differ in both hue and lightness); event markers stay solid INK per house style.
STATE_ON  = "#2f7d84"      # active state: heavy branch / filled chip  (muted teal)
STATE_OFF = "#b9b0a4"      # inactive state: light branch / open chip  (warm taupe)

# --- Coupled / Potts figures ---------------------------------------------
# Categorical palette for the pathway MODULES of the coupled (Potts) chapter, where a
# few module identities must be told apart (the STYLE.md "categorical exception").
# Colour-blind-safe and deliberately distinct from the green/red used below for the
# co-occur/avoid (and protected/lost) *semantics*, so the two never get confused.
MODULE_COLORS = ["#4477AA", "#E08A3C", "#7B5EA7", "#4C9AA6"]   # blue, orange, purple, teal
COOCCUR = "#2f8f4e"         # J>0 / partners present / kept / protected  (green)
AVOID   = "#cc4b3c"         # J<0 / partners absent / purged / fast loss (red)

# --- Stroke weights (px) --------------------------------------------------
BRANCH_W = 2.6


def species_style(width: int = 820, height: int = 680, **overrides) -> ph.TreeStyle:
    """House style for a time-calibrated species tree.

    Clean rectangular cladogram: no tip/node dots by default (labels carry the
    identity), near-black branches, generous margins for labels and a time axis.
    Pass any :class:`phylustrator.TreeStyle` field as a keyword to override.
    """
    params = dict(
        width=width,
        height=height,
        margin=82.0,
        root_stub_length=14.0,
        branch_stroke_width=BRANCH_W,
        branch_color=INK,
        leaf_r=0.0,          # no dot at the tip; the label is enough
        node_r=0.0,          # clean bifurcations, no internal-node dots
        font_size=17,
        font_family=FONT,
    )
    params.update(overrides)
    return ph.TreeStyle(**params)
