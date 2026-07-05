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
# drawsvg figure that imports them.
FS_TITLE = 32               # figure title (bold)
FS_LABEL = 22               # axis titles, legend entries, curve labels
FS_ANNOT = 20               # inline annotations
FS_TICK  = 18               # tick numbers and other small numeric labels

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
