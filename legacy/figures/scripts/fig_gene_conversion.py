"""Figure: intra-genome gene conversion -- mechanism (A) and a complete gene tree (B).

Two panels, same house style as fig_gene_tree / fig_gene_tree_panels:

  * A - MECHANISM (schematic). One gene family inside one lineage. A duplication makes two
        copies; later a conversion has one copy (the template/donor) overwrite the other. Copy
        number is unchanged, but the overwritten copy's own history is erased, so the two
        surviving copies now coalesce at the CONVERSION time (recent), not at the DUPLICATION
        (old). That downward pull of the coalescence is the whole point of the model.

  * B - A COMPLETE GENE TREE (real ZOMBI2 output). A small family simulated with
        SharedRates(conversion=...) that happens to carry exactly one conversion; the complete
        gene tree is drawn with the conversion marked, so the concerted-evolution signature is
        visible directly against the duplications and losses.

Glyphs (solid black, per STYLE.md): filled square = duplication, filled CIRCLE = conversion,
cross (x) = loss; solid vs dashed = surviving vs overwritten/lost lineage.

Run:  python figures/scripts/fig_gene_conversion.py
"""

from __future__ import annotations

import re
from pathlib import Path

import cairosvg
import drawsvg as draw
import phylustrator as ph
from phylustrator.io import read_newick

import zombi2 as z
from zombi2.genomes.events import EventType
from fig_species_tree_events import draw_cross as _draw_cross_ink
from fig_species_tree_extinct import annotate_depths, draw_skeleton, mark_survival
from zombi_style import INK, MUTED, species_style

FIG_DIR = Path(__file__).resolve().parent.parent
OUT_DIR = FIG_DIR / "gene_conversion"
OUT_STEM = OUT_DIR / "gene_conversion"

# per-panel geometry (composed into one wide ~2:1 figure, so fonts are set large)
PW, PH = 900, 900
FS_TITLE_P = 46
FS_PANEL = 40
FS_LABEL_P = 34
FS_ANNOT_P = 30
FS_TICK_P = 30
MARKER_R = 14.0
FOOTER = 130


# =============================================================================== Panel A (schematic)

def _t(x, y, s, size, *, anchor="start", weight="normal", color=INK, italic=False):
    style = ' font-style="italic"' if italic else ""
    return (f'<text x="{x:.1f}" y="{y:.1f}" font-size="{size}" font-family="Helvetica" '
            f'text-anchor="{anchor}" dominant-baseline="middle" font-weight="{weight}"{style} '
            f'fill="{color}">{s}</text>')


def _line(x1, y1, x2, y2, *, w=4.5, dash=None, color=INK, cap="round"):
    d = f' stroke-dasharray="{dash}"' if dash else ""
    return (f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" stroke="{color}" '
            f'stroke-width="{w}" stroke-linecap="{cap}"{d} />')


def _square(x, y, r):
    return f'<rect x="{x-r:.1f}" y="{y-r:.1f}" width="{2*r:.1f}" height="{2*r:.1f}" fill="{INK}" />'


def _circle(x, y, r):
    return f'<circle cx="{x:.1f}" cy="{y:.1f}" r="{r:.1f}" fill="{INK}" />'


def _cross(x, y, r, w=4.5):
    return (_line(x - r, y - r, x + r, y + r, w=w) + _line(x - r, y + r, x + r, y - r, w=w))


def _conv(x, y, r):
    """Conversion glyph: an OPEN (white-filled) circle with a centred cross -- distinct from the
    solid dot used for origination and from the bare cross used for loss."""
    cr = r * 0.58
    return (f'<circle cx="{x:.1f}" cy="{y:.1f}" r="{r:.1f}" fill="white" stroke="{INK}" '
            f'stroke-width="3.2" />' + _cross(x, y, cr, w=3.2))


def _draw_conv_glyph(d, x, y, r):
    """Same glyph, appended to a phylustrator drawer (Panel B)."""
    cr = r * 0.58
    d.drawing.append(draw.Circle(x, y, r, fill="white", stroke=INK, stroke_width=3.2))
    for x1, y1, x2, y2 in ((x - cr, y - cr, x + cr, y + cr), (x - cr, y + cr, x + cr, y - cr)):
        d.drawing.append(draw.Line(x1, y1, x2, y2, stroke=INK, stroke_width=3.2,
                                   stroke_linecap="round"))


def panel_a() -> str:
    """The mechanism, as a raw SVG group of size (PW, PH). Time runs top (past) -> bottom (present)."""
    o = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{PW}" height="{PH}" '
         f'viewBox="0 0 {PW} {PH}">',
         f'<rect x="0" y="0" width="{PW}" height="{PH}" fill="white" />']

    # panel label
    o.append(_t(40, 58, "A", FS_TITLE_P, weight="bold"))
    o.append(_t(96, 58, "The mechanism", FS_PANEL))

    # time axis (down = toward the present), words not numbers (schematic)
    ax = 96
    y0, y_pres = 170, 812
    o.append(_line(ax, y0, ax, y_pres, w=2.4, color=MUTED, cap="butt"))
    o.append(f'<polygon points="{ax},{y_pres+16} {ax-8},{y_pres-2} {ax+8},{y_pres-2}" '
             f'fill="{MUTED}" />')
    o.append(_t(ax - 16, (y0 + y_pres) / 2, "time", FS_ANNOT_P, anchor="middle", color=MUTED)
             .replace('<text ', '<text transform="rotate(-90 %.0f %.0f)" ' % (ax - 16, (y0 + y_pres) / 2)))

    # event times and copy positions
    y_dup, y_conv = 300, 585
    xX, xY = 360, 560          # copy X (donor/template) and copy Y (recipient)
    xroot = (xX + xY) / 2

    # family origin -> duplication (single ancestral lineage)
    o.append(_line(xroot, 200, xroot, y_dup))
    # duplication: horizontal split + square glyph (glyph carries the identity; no label)
    o.append(_line(xX, y_dup, xY, y_dup))
    o.append(_square(xroot, y_dup, MARKER_R))

    # copy X: solid all the way down (the template that survives)
    o.append(_line(xX, y_dup, xX, y_pres))
    # copy Y BEFORE conversion: its own history, later erased -> dashed, capped by a loss cross
    o.append(_line(xY, y_dup, xY, y_conv, dash="9,7", cap="butt"))

    # conversion: X donates a copy that overwrites Y. Connector to Y drawn first; Y's position
    # continues solid below (now occupied by the converted copy of X); loss cross on the old Y;
    # the conversion glyph sits on X's lineage last, on top.
    o.append(_line(xX, y_conv, xY, y_conv))
    o.append(_line(xY, y_conv, xY, y_pres))
    o.append(_cross(xY, y_conv, MARKER_R - 1))
    o.append(_conv(xX, y_conv, MARKER_R))
    o.append(_t(xX - 26, y_conv, "conversion", FS_ANNOT_P, anchor="end"))

    # tips at the present
    o.append(_t(xX, y_pres + 34, "copy 1", FS_ANNOT_P, anchor="middle"))
    o.append(_t(xY, y_pres + 34, "copy 2", FS_ANNOT_P, anchor="middle"))

    # bracket: the two copies now coalesce at the conversion (recent), not the duplication (old)
    bx = 648
    o.append(_line(bx, y_conv, bx, y_pres, w=2.4, color=INK, cap="butt"))
    o.append(_line(bx - 8, y_conv, bx, y_conv, w=2.4, color=INK, cap="butt"))
    o.append(_line(bx - 8, y_pres, bx, y_pres, w=2.4, color=INK, cap="butt"))
    o.append(_t(bx + 14, (y_conv + y_pres) / 2 - 18, "copies coalesce", FS_ANNOT_P - 4))
    o.append(_t(bx + 14, (y_conv + y_pres) / 2 + 16, "at the conversion,", FS_ANNOT_P - 4))
    o.append(_t(bx + 14, (y_conv + y_pres) / 2 + 48, "not the duplication", FS_ANNOT_P - 4,
                color=MUTED))

    # copy-number invariant note
    o.append(_t(xroot, 128, "copy number is unchanged: 2 copies before and after",
                FS_ANNOT_P - 4, anchor="middle", color=MUTED))

    o.append("</svg>")
    return "\n".join(o)


# ============================================================================= Panel B (real output)

def build_scenario():
    """Find a small, legible family with EXACTLY one duplication and one conversion (and as few
    extra losses as possible), by a reproducible seed search."""
    model = z.BirthDeath(birth=1.0, death=0.3)
    for seed in range(1, 1200):
        tree = z.simulate_species_tree(model, n_tips=6, age=1.0, direction="backward", seed=seed)
        g = z.simulate_genomes(
            tree, z.SharedRates(duplication=0.5, conversion=2.2, loss=0.06, origination=0.5),
            conversions=z.ConversionModel(bias=1.0), initial_families=12, seed=seed)
        trees = g.gene_trees()
        for fam, records in g.gene_families.items():
            dups = [r for r in records if r.event is EventType.DUPLICATION]
            convs = [r for r in records if r.event is EventType.CONVERSION]
            losses = [r for r in records if r.event is EventType.LOSS]
            if len(dups) != 1 or len(convs) != 1 or fam not in trees:
                continue
            complete, _extant = trees[fam]
            n_tips = complete.count(",") + 1
            # the conversion itself contributes one loss (the overwritten copy); allow at most one
            # more so the tree stays legible
            if 4 <= n_tips <= 7 and len(losses) <= 2:
                return tree, g, fam, records, convs[0], complete, seed
    raise RuntimeError("no clean single-duplication single-conversion family found")


def _species_letters(tree):
    leaves = [n for n in tree.nodes_preorder() if not n.children]
    return {lf.name: chr(ord("A") + i) for i, lf in enumerate(leaves)}


def panel_b(complete_nwk: str, records, conv_record, species_letter):
    """Draw the complete gene tree; mark duplications (square), the conversion (circle) and
    losses (cross)."""
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    nwk_path = OUT_DIR / "complete.nwk"
    nwk_path.write_text(complete_nwk.strip() + "\n")
    tree = read_newick(str(nwk_path))
    present = annotate_depths(tree)
    mark_survival(tree, present, tol=1e-6 * present)
    # gid -> node, resolving internal nodes (named by gid) AND tips (named "<species>_<gid>"), so a
    # conversion/duplication child gid can be located even when it is a tip
    gid2node = {(n.name.split("_")[-1] if "_" in n.name else n.name): n for n in tree.traverse()}

    # tips like <species-letter>_<copy>; blank + collect the loss tips
    loss_tips, copies = [], {}
    for lf in tree.get_leaves():
        if lf.name.startswith("LOSS_"):
            loss_tips.append(lf)
            lf.name = ""
        else:
            sp = species_letter.get(lf.name.split("_")[0], lf.name.split("_")[0])
            copies[sp] = copies.get(sp, 0) + 1
            lf.name = f"{sp}_{copies[sp]}"

    style = species_style(width=PW, height=PH, margin=140, font_size=FS_TICK_P)
    d = ph.VerticalTreeDrawer(tree, style=style)
    d._calculate_layout()
    draw_skeleton(d, tree)
    d.add_leaf_names(color=INK, padding=14)

    dup_nodes = {r.genes[0].gid for r in records if r.event is EventType.DUPLICATION}
    conv_node = conv_record.genes[0].gid
    for gid in dup_nodes:
        if gid in gid2node and gid != conv_node:
            d._draw_shape_at(*gid2node[gid].coordinates, "square", INK, r=MARKER_R)
    if conv_node in gid2node:
        cx, cy = gid2node[conv_node].coordinates
        # name the two products of the conversion; place the labels to the LEFT of the node (into
        # the open incoming-branch space) so they never collide with the tip labels at the edge —
        # donor above the branch, converted below, matching the two children's vertical order
        d.add_text("donor copy", cx - 22, cy - 28, font_size=FS_TICK_P - 4, color=MUTED,
                   text_anchor="end")
        d.add_text("converted copy", cx - 22, cy + 28, font_size=FS_TICK_P - 4, color=MUTED,
                   text_anchor="end")
        _draw_conv_glyph(d, cx, cy, MARKER_R)   # glyph on top, after the labels
    for lf in loss_tips:
        _draw_cross_ink(d, *lf.coordinates, MARKER_R, stroke_width=4.0)

    ticks = [round(present * i / 2, 6) for i in range(3)]
    d.add_time_axis(ticks=ticks, tick_labels=[f"{t:.1f}" for t in ticks],
                    label="Time (root to present)", tick_size=8.0, padding=16.0, stroke_width=1.8)
    d.add_text("B", x=-PW / 2 + 26, y=-PH / 2 + 58, font_size=FS_TITLE_P, color=INK, weight="bold")
    d.add_text("A complete gene tree", x=-PW / 2 + 90, y=-PH / 2 + 58,
               font_size=FS_PANEL, color=INK)
    return d, len(tree.get_leaves())


# ===================================================================================== compose + main

def _footer_legend(cx, y):
    out, r = [], 15
    x = cx - 620

    def text(tx, s):
        out.append(_t(tx, y, s, FS_LABEL_P))

    out.append(_square(x, y, r)); text(x + r + 14, "Duplication"); x += 300
    out.append(_conv(x, y, r)); text(x + r + 14, "Conversion"); x += 300
    out.append(_cross(x, y, r)); text(x + r + 14, "Loss"); x += 190
    out.append(_line(x, y, x + 54, y)); text(x + 66, "surviving copy"); x += 370
    out.append(_line(x, y, x + 54, y, dash="8,6", cap="butt")); text(x + 66, "overwritten / lost")
    return "\n".join(out)


def compose(panel_a_svg, drawer_b, out_stem):
    W, H = PW * 2, PH + FOOTER
    parts = [f'<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" '
             f'width="{W}" height="{H}" viewBox="0 0 {W} {H}">',
             f'<rect x="0" y="0" width="{W}" height="{H}" fill="white" />']
    a = re.sub(r"<\?xml[^>]*\?>\s*", "", panel_a_svg)
    a = a.replace("<svg ", '<svg x="0" y="0" ', 1)
    parts.append(a)
    b = re.sub(r"<\?xml[^>]*\?>\s*", "", drawer_b.drawing.as_svg())
    b = b.replace("<svg ", f'<svg x="{PW}" y="0" ', 1)
    parts.append(b)
    parts.append(_footer_legend(W / 2, PH + FOOTER / 2))
    parts.append("</svg>")
    parent = "\n".join(parts)
    out_stem = Path(out_stem)
    out_stem.parent.mkdir(parents=True, exist_ok=True)
    out_stem.with_suffix(".svg").write_text(parent, encoding="utf-8")
    cairosvg.svg2png(bytestring=parent.encode("utf-8"),
                     write_to=str(out_stem.with_suffix(".png")), scale=300 / 72.0)


def main():
    tree, g, fam, records, conv_record, complete, seed = build_scenario()
    letters = _species_letters(tree)
    b, nb = panel_b(complete, records, conv_record, letters)
    compose(panel_a(), b, OUT_STEM)
    print(f"wrote {OUT_STEM}.svg / .png  (seed {seed}, family {fam}, complete={nb} tips)")


if __name__ == "__main__":
    main()
