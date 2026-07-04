"""Figure: the emergent signal of the Potts model — co-occurrence modules in the profiles.

Coupling makes co-functional families present/absent TOGETHER across genomes. Simulated with
``simulate_coupled`` over a species tree (pathway_blocks: 3 modules of 4 families). Left: the
phylogenetic-profile matrix (families x genomes), families grouped by module — whole modules
switch on/off together across genomes. Right: the family x family co-occurrence matrix —
featureless when families are independent (J=0), block-diagonal once they are coupled. That
block structure is exactly the signal inverse-Potts / DCA methods recover; ZOMBI2 generates
it with a known ground-truth J.

Colour, didactic.  Run:  python figures/scripts/fig_potts_modules.py
"""

from __future__ import annotations

import warnings
from pathlib import Path

import cairosvg
import drawsvg as draw
import numpy as np

from zombi2 import BirthDeath, simulate_coupled, simulate_species_tree
from zombi2.coupling import pathway_blocks

from zombi_style import FONT, INK

warnings.filterwarnings("ignore")
OUT = Path(__file__).resolve().parent.parent / "potts_modules"

W, H = 1240, 720
MOD = ["#4477AA", "#E08A3C", "#2E8B7A"]           # 3 modules
SIZES = [4, 4, 4]
N = sum(SIZES)


def simulate(within):
    tree = simulate_species_tree(BirthDeath(1.0, 0.2), n_tips=34, age=5.0, seed=7)
    spec = pathway_blocks(SIZES, within=within, between=0.0, base_loss=1.0, transfer=0.4,
                          beta=1.4, h=[-2.0] * N)
    A = (np.asarray(simulate_coupled(tree, spec, seed=3).profiles.matrix) > 0).astype(int)
    return A


def cooc(A):
    C = np.corrcoef(A)
    return np.nan_to_num(C, nan=0.0)


def module_of(i):
    return next(m for m, s in enumerate([sum(SIZES[:k + 1]) for k in range(3)]) if i < s)


def heat(v):
    v = max(0.0, min(1.0, v))                     # co-occurrence 0..1 -> white..teal
    r = int(round(255 - 200 * v)); g = int(round(255 - 120 * v)); b = int(round(255 - 130 * v))
    return "#%02x%02x%02x" % (r, g, b)


def main():
    A_c = simulate(1.2)
    A_i = simulate(0.0)
    # order genomes (columns) so similar module-content groups together — reveals the blocks
    modsum = np.vstack([A_c[m * 4:(m + 1) * 4].sum(0) for m in range(3)])
    order = sorted(range(A_c.shape[1]), key=lambda j: tuple(-modsum[m, j] for m in range(3)))

    d = draw.Drawing(W, H, origin=(0, 0))
    d.append(draw.Rectangle(0, 0, W, H, fill="white"))
    d.append(draw.Text("The emergent signal — coupling writes co-occurrence modules into the profiles",
                       20, 40, 44, font_family=FONT, text_anchor="start", font_weight="bold", fill=INK))
    d.append(draw.Text("co-functional families share presence/absence profiles (Pellegrini 1999); "
                       "ZOMBI2 injects a known J and the modules re-emerge", 13.5, 40, 70,
                       font_family=FONT, text_anchor="start", fill="#777"))

    # ---- left: profile matrix ----
    d.append(draw.Text("Phylogenetic profile  (families x genomes)", 15, 96, 118, font_family=FONT,
                       text_anchor="start", font_weight="bold", fill=INK))
    x0, y0, cw, ch = 150, 150, 16.5, 20
    ns = A_c.shape[1]
    for i in range(N):
        for jj, j in enumerate(order):
            x, y = x0 + jj * cw, y0 + i * ch
            if A_c[i, j]:
                d.append(draw.Rectangle(x, y, cw + 0.5, ch + 0.5, fill=MOD[module_of(i)]))
            else:
                d.append(draw.Rectangle(x, y, cw + 0.5, ch + 0.5, fill="#f3f4f6"))
    d.append(draw.Rectangle(x0, y0, ns * cw, N * ch, fill="none", stroke=INK, stroke_width=1.2))
    for m in range(1, 3):                                            # module separators
        d.append(draw.Line(x0, y0 + m * 4 * ch, x0 + ns * cw, y0 + m * 4 * ch, stroke="white", stroke_width=2))
    for m in range(3):                                              # module brackets + labels
        ym = y0 + m * 4 * ch
        d.append(draw.Rectangle(x0 - 16, ym, 8, 4 * ch, fill=MOD[m]))
        d.append(draw.Text(f"module {chr(65 + m)}", 11.5, x0 - 22, ym + 2 * ch, font_family=FONT,
                           text_anchor="end", dominant_baseline="central", fill=MOD[m], font_weight="bold"))
    d.append(draw.Text(f"{ns} genomes  (columns, ordered by content)", 12, x0 + ns * cw / 2, y0 + N * ch + 20,
                       font_family=FONT, text_anchor="middle", fill="#777"))
    d.append(draw.Text("whole modules are present or absent together across genomes", 12.5, 150,
                       y0 + N * ch + 44, font_family=FONT, text_anchor="start", fill="#555", font_style="italic"))

    # ---- right: co-occurrence matrices ----
    def matrix(A, mx, my, title, sub):
        C = cooc(A)
        s = 13.5
        d.append(draw.Text(title, 14, mx, my - 26, font_family=FONT, text_anchor="start",
                           font_weight="bold", fill=INK))
        d.append(draw.Text(sub, 11.5, mx, my - 10, font_family=FONT, text_anchor="start", fill="#777"))
        for i in range(N):
            for j in range(N):
                d.append(draw.Rectangle(mx + j * s, my + i * s, s + 0.4, s + 0.4, fill=heat(C[i, j])))
        d.append(draw.Rectangle(mx, my, N * s, N * s, fill="none", stroke=INK, stroke_width=1.1))
        for m in range(1, 3):
            d.append(draw.Line(mx + m * 4 * s, my, mx + m * 4 * s, my + N * s, stroke="#c8c8c8", stroke_width=0.8))
            d.append(draw.Line(mx, my + m * 4 * s, mx + N * s, my + m * 4 * s, stroke="#c8c8c8", stroke_width=0.8))
        for m in range(3):
            d.append(draw.Rectangle(mx + m * 4 * s, my - 6, 4 * s, 4, fill=MOD[m]))

    mx = 820
    matrix(A_i, mx, 168, "Independent  (J = 0)", "each family on its own — no blocks")
    matrix(A_c, mx, 468, "Coupled  (Potts)", "block-diagonal — the modules emerge")

    # colourbar (between the two matrices, clear of the labels)
    cbx, cby, cbw = mx + 262, 372, 16
    for k in range(66):
        d.append(draw.Rectangle(cbx, cby + k * 3, cbw, 3.4, fill=heat(1 - k / 65)))
    d.append(draw.Rectangle(cbx, cby, cbw, 198, fill="none", stroke=INK, stroke_width=0.8))
    d.append(draw.Text("co-occurrence", 11.5, cbx + cbw / 2, cby - 12, font_family=FONT,
                       text_anchor="middle", fill="#777"))
    d.append(draw.Text("high", 11, cbx + cbw + 5, cby + 6, font_family=FONT, text_anchor="start",
                       dominant_baseline="central", fill="#777"))
    d.append(draw.Text("0", 11, cbx + cbw + 5, cby + 192, font_family=FONT, text_anchor="start",
                       dominant_baseline="central", fill="#777"))

    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "potts_modules.svg").write_text(d.as_svg(), encoding="utf-8")
    cairosvg.svg2png(bytestring=d.as_svg().encode(), write_to=str(OUT / "potts_modules.png"), scale=300 / 72.0)
    print("wrote potts_modules")


if __name__ == "__main__":
    main()
