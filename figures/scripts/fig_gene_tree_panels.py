"""Figure: the COMPLETE gene tree beside the EXTANT (pruned) gene tree.

Two panels of the *same* gene family, side by side:

  * A - COMPLETE: every copy the family ever had, including the lineage that was
        later lost (drawn dashed, ending in a cross). This is what the simulator
        actually produced.
  * B - EXTANT: the same tree pruned to the copies that survive to the present,
        with the resulting degree-two nodes suppressed. This is the tree you get
        back for a real, sampled family.

Reading the two together is the whole point: the extant tree is what you would
ever observe, the complete tree is the history that generated it, and the loss is
exactly the difference between them.

The trees are REAL ZOMBI2 output (not hand-placed). A 10-tip backward birth-death
tree (seed 24) is run through ``simulate_genomes`` (seed 5); family 4 has a clean
signature -- one duplication, one transfer, one lost lineage -- with the
duplication early (t~0.25) and the transfer late (t~0.92) so the two are widely
separated on the tree. Its ``(complete, extant)`` pair is drawn directly from
``gene_trees()``.

Encoding (same house style as fig_gene_tree / the species-tree figures)
  * filled square      -> duplication (two sister copies in the same species)
  * filled triangle    -> transfer  (donor-kept copy vs transferred copy)
  * cross (x)          -> loss (a dead-end lineage, drawn dashed)
  * solid vs dashed    -> surviving vs lost lineage

Run:  python figures/scripts/fig_gene_tree_panels.py
"""

from __future__ import annotations

import re
from pathlib import Path

import cairosvg
import drawsvg as draw
import phylustrator as ph
from phylustrator.io import read_newick

import zombi2 as z
from fig_species_tree_events import draw_cross as _draw_cross_ink
from fig_species_tree_extinct import annotate_depths, draw_skeleton, mark_survival
from zombi_style import INK, MUTED, species_style

FIG_DIR = Path(__file__).resolve().parent.parent
OUT_DIR = FIG_DIR / "gene_tree_panels"
OUT_STEM = OUT_DIR / "gene_tree_panels"
COMPLETE_NWK = OUT_DIR / "complete.nwk"
EXTANT_NWK = OUT_DIR / "extant.nwk"

# --- scenario (reproducible) ------------------------------------------------
TREE_SEED, N_TIPS = 24, 10
GENOME_SEED = 5
FAMILY = "4"
# species (ZOMBI node) -> display letter, from the species tree's top-to-bottom
# leaf order (see fig_gene_tree_panels header / detail_ch9_v2.py).
SPECIES_LETTER = {"n4": "A", "n5": "B", "n10": "C", "n2": "D", "n6": "E",
                  "n1": "F", "n7": "G", "n8": "H", "n3": "I", "n9": "J"}

# reconciliation annotation for family 4 (from its event log):
#   duplication at gene node g104 (two copies in species F), early (t~0.25)
#   transfer    at gene node g340: donor-kept copy g693, transferred copy g694,
#               late (t~0.92) -- so dup and transfer sit far apart on the tree
DUP_NODE = "g104"
TRANSFER_NODE = "g340"
TRANSFER_KEPT = "g693"          # child that stays in the donor lineage (-> G, H)
TRANSFER_MOVED = "n9_g694"      # child that was transferred away (leaf -> J);
                               # a tip, so its node name still carries the species

# per-panel geometry. The two panels compose into one wide (~2:1) figure, so on
# the page labels scale down; fonts are therefore set MUCH larger than a
# single-panel figure would use (matching fig_model_ghosts).
PW, PH = 900, 900
FS_TITLE_P = 46
FS_PANEL = 40
FS_LABEL_P = 34
FS_ANNOT_P = 30
FS_TICK_P = 30
MARKER_R = 13.0
FOOTER = 130


def regenerate_trees():
    """Run the scenario and (re)write the complete/extant Newick fixtures."""
    model = z.BirthDeath(birth=1.0, death=0.4)
    tree = z.simulate_species_tree(model, n_tips=N_TIPS, age=1.0,
                                   direction="backward", seed=TREE_SEED)
    g = z.simulate_genomes(tree, duplication=0.25, transfer=0.2, loss=0.3,
                           origination=0.5, initial_families=30, seed=GENOME_SEED)
    complete, extant = g.gene_trees()[FAMILY]
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    COMPLETE_NWK.write_text(complete.strip() + "\n")
    EXTANT_NWK.write_text(extant.strip() + "\n")


def species_of(tip_name: str) -> str:
    if tip_name.startswith("LOSS_"):
        return ""
    return SPECIES_LETTER.get(tip_name.split("_")[0], tip_name.split("_")[0])


def label_tips(tree):
    """Rename tips <species>_<copy> (so duplicates read C_1 / C_2); blank the loss
    tip (marked only by the dead-end + cross). Returns the list of loss tips."""
    loss_tips, copies = [], {}
    for lf in tree.get_leaves():
        if lf.name.startswith("LOSS_"):
            loss_tips.append(lf)
            lf.name = ""
        else:
            sp = species_of(lf.name)
            copies[sp] = copies.get(sp, 0) + 1
            lf.name = f"{sp}_{copies[sp]}"
    return loss_tips


def panel(nwk_path, label, subtitle, *, complete):
    tree = read_newick(str(nwk_path))
    present = annotate_depths(tree)
    mark_survival(tree, present, tol=1e-6 * present)
    name2node = {n.name: n for n in tree.traverse()}
    loss_tips = label_tips(tree)

    style = species_style(width=PW, height=PH, margin=132, font_size=FS_TICK_P)
    d = ph.VerticalTreeDrawer(tree, style=style)
    d._calculate_layout()
    draw_skeleton(d, tree)                       # solid survivors, dashed lost lineage
    d.add_leaf_names(color=INK, padding=14)

    # duplication square
    if DUP_NODE in name2node:
        d._draw_shape_at(*name2node[DUP_NODE].coordinates, "square", INK, r=MARKER_R)
    # transfer triangle + donor/transferred copy labels. Place the labels toward
    # whichever side has room: to the LEFT (anchored end) when the transfer sits
    # in the right third of the tree (so the text can't run off the panel or hit
    # the tip labels), otherwise to the RIGHT. This keeps the figure robust to
    # where the transfer happens to fall.
    if TRANSFER_NODE in name2node:
        tn = name2node[TRANSFER_NODE]
        d._draw_shape_at(*tn.coordinates, "triangle", INK, r=MARKER_R)
        xn = tn.coordinates[0]
        right_edge = d.root_x + present * d.sf
        put_left = xn > d.root_x + 0.62 * (right_edge - d.root_x)
        anchor = "end" if put_left else "start"
        dx = -22 if put_left else 22
        for role, child in (("donor-kept copy", TRANSFER_KEPT),
                            ("transferred copy", TRANSFER_MOVED)):
            if child in name2node:
                cy = name2node[child].coordinates[1]
                d.add_text(role, xn + dx, cy - 16, font_size=FS_ANNOT_P,
                           color=MUTED, text_anchor=anchor)
    # loss crosses (complete panel only; extant has none)
    for lf in loss_tips:
        _draw_cross_ink(d, *lf.coordinates, MARKER_R, stroke_width=4.0)

    ticks = [round(present * i / 2, 6) for i in range(3)]
    d.add_time_axis(ticks=ticks, tick_labels=[f"{t:.1f}" for t in ticks],
                    label="Time (root to present)", tick_size=8.0, padding=16.0,
                    stroke_width=1.8)

    d.add_text(label, x=-PW / 2 + 26, y=-PH / 2 + 58, font_size=FS_TITLE_P,
               color=INK, weight="bold")
    d.add_text(subtitle, x=-PW / 2 + 90, y=-PH / 2 + 58, font_size=FS_PANEL, color=INK)
    return d, len(tree.get_leaves())


def _footer_legend(cx, y):
    """Raw-SVG key: the three event glyphs + solid/dashed lineage swatch,
    one centred row. Kept colour-neutral (INK), matching STYLE.md."""
    fam = "Helvetica"
    out = []
    x = cx - 720
    r = 15

    def text(tx, s):
        out.append(f'<text x="{tx}" y="{y}" font-size="{FS_LABEL_P}" font-family="{fam}" '
                   f'text-anchor="start" dominant-baseline="middle" fill="{INK}">{s}</text>')

    # duplication: filled square
    out.append(f'<rect x="{x - r}" y="{y - r}" width="{2*r}" height="{2*r}" fill="{INK}" />')
    text(x + r + 14, "Duplication"); x += 330
    # transfer: filled triangle
    out.append(f'<polygon points="{x},{y - r} {x + r},{y + r} {x - r},{y + r}" fill="{INK}" />')
    text(x + r + 14, "Transfer"); x += 300
    # loss: cross
    out.append(f'<line x1="{x - r}" y1="{y - r}" x2="{x + r}" y2="{y + r}" stroke="{INK}" '
               f'stroke-width="4.5" stroke-linecap="round" />')
    out.append(f'<line x1="{x - r}" y1="{y + r}" x2="{x + r}" y2="{y - r}" stroke="{INK}" '
               f'stroke-width="4.5" stroke-linecap="round" />')
    text(x + r + 14, "Loss"); x += 230
    # solid swatch
    out.append(f'<line x1="{x}" y1="{y}" x2="{x + 54}" y2="{y}" stroke="{INK}" '
               f'stroke-width="4.5" stroke-linecap="round" />')
    text(x + 66, "surviving copy"); x += 430
    # dashed swatch
    out.append(f'<line x1="{x}" y1="{y}" x2="{x + 54}" y2="{y}" stroke="{INK}" '
               f'stroke-width="4.5" stroke-dasharray="8,6" stroke-linecap="butt" />')
    text(x + 66, "lost lineage")
    return "\n".join(out)


def compose(drawers, out_stem):
    n = len(drawers)
    W, H = PW * n, PH + FOOTER
    parts = [f'<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" '
             f'width="{W}" height="{H}" viewBox="0 0 {W} {H}">',
             f'<rect x="0" y="0" width="{W}" height="{H}" fill="white" />']
    for i, d in enumerate(drawers):
        svg = re.sub(r"<\?xml[^>]*\?>\s*", "", d.drawing.as_svg())
        svg = svg.replace("<svg ", f'<svg x="{PW * i}" y="0" ', 1)
        parts.append(svg)
    parts.append(_footer_legend(W / 2, PH + FOOTER / 2))
    parts.append("</svg>")
    parent = "\n".join(parts)
    out_stem = Path(out_stem)
    out_stem.parent.mkdir(parents=True, exist_ok=True)
    out_stem.with_suffix(".svg").write_text(parent, encoding="utf-8")
    cairosvg.svg2png(bytestring=parent.encode("utf-8"),
                     write_to=str(out_stem.with_suffix(".png")), scale=300 / 72.0)


def main():
    regenerate_trees()
    a, na = panel(COMPLETE_NWK, "A", "Complete gene tree", complete=True)
    b, nb = panel(EXTANT_NWK, "B", "Extant gene tree", complete=False)
    compose([a, b], OUT_STEM)
    print(f"wrote {OUT_STEM}.svg / .png  (complete={na} tips, extant={nb} tips)")


if __name__ == "__main__":
    main()
