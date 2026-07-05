"""Figure: a 10-species time-calibrated species tree simulated with ZOMBI2.

Reads the Newick that ``zombi2 species`` wrote, relabels the tips with clean
display names, and renders a consistent-style vertical tree to SVG + PNG.

Run:  python figures/scripts/fig_species_tree.py
"""

from __future__ import annotations

import string
from pathlib import Path

import drawsvg as draw
import phylustrator as ph
from phylustrator.io import read_newick

from zombi_style import INK, MUTED, species_style, FS_TITLE, FS_TICK

# --- Paths ----------------------------------------------------------------
FIG_DIR = Path(__file__).resolve().parent.parent          # .../figures
TREE_NWK = FIG_DIR / "species_tree_10" / "species_tree.nwk"
OUT_STEM = FIG_DIR / "species_tree_10" / "species_tree"   # -> .svg / .png


def clean_leaf_labels(tree) -> None:
    """Relabel tips top-to-bottom as A, B, C, ... for a tidy illustration.

    ``get_leaves()`` returns tips in the tree's drawing (top-to-bottom) order,
    so the labels come out alphabetical down the figure.
    """
    for leaf, letter in zip(tree.get_leaves(), string.ascii_uppercase):
        leaf.name = letter


def nice_ticks(depth: float, n: int = 4) -> list[float]:
    """``n+1`` evenly spaced ticks from 0 to ``depth`` (rounded to the depth)."""
    return [round(depth * i / n, 6) for i in range(n + 1)]


def main() -> None:
    tree = read_newick(TREE_NWK)
    clean_leaf_labels(tree)

    # generous top margin leaves a clean band for the centered title above the tree
    style = species_style(height=740, margin=110, font_size=FS_TICK)
    d = ph.VerticalTreeDrawer(tree, style=style)
    d.draw()
    d.add_leaf_names(color=INK, padding=12)

    depth = d.total_tree_depth
    d.add_time_axis(
        ticks=nice_ticks(depth, 4),
        label="Time (root to present)",
        tick_size=6.0,
        padding=14.0,
        stroke_width=1.6,
    )

    # title: one short bold line, horizontally centered at the top
    d.drawing.append(draw.Text("A birth-death species tree", FS_TITLE, 0,
                               -style.height / 2 + 44, font_weight="bold",
                               font_family=style.font_family, text_anchor="middle",
                               dominant_baseline="central", fill=INK))

    d.save_svg(f"{OUT_STEM}.svg")
    d.save_png(f"{OUT_STEM}.png", dpi=300)
    print(f"wrote {OUT_STEM}.svg and {OUT_STEM}.png  ({len(tree.get_leaves())} tips)")


if __name__ == "__main__":
    main()
