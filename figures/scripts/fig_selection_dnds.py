"""Figure (experimental / language-model selection): emergent dN/dS from codon mutation-selection.

Left: the model splits a codon substitution into MUTATION (a nucleotide change, from the neutral
model) times SELECTION (a fixation factor from the language model's preference for the encoded amino
acid). A synonymous change keeps the amino acid, so its factor is 1 (neutral); a non-synonymous change
is scrutinised by the critic. dN/dS is therefore an OUTPUT of the model, not a parameter.

Right: the model's expected genome-wide dN/dS (omega) as the selection strength beta rises -- computed
from the real codon model with a preference peaked on a natural protein (ubiquitin). omega is exactly 1
at beta = 0 (neutral) and decreases monotonically toward 0 under stronger purifying selection.

House style: one centered bold title, ASCII text, viridis-free (a single accent curve).
Run:  PYTHONPATH=. python figures/scripts/fig_selection_dnds.py
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import cairosvg
import drawsvg as draw
import numpy as np

from zombi_style import FONT, INK, MUTED, FS_TITLE, FS_LABEL, FS_TICK

from zombi2.experimental.codon_selection import CodonSelection
from zombi2.experimental.selection import FixedProfileCritic
from zombi2.sequences.models import AMINO_ACIDS

OUT = Path(__file__).resolve().parent.parent.parent / "docs" / "img"
W, H = 1180, 560
ACCENT = "#377EB8"     # the omega(beta) curve
GOOD = "#4DAF4A"       # synonymous / neutral
SCRUT = "#E41A1C"      # non-synonymous / scrutinised
UBIQUITIN = "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG"


def _peaked(protein, hi=0.9):
    lo = (1.0 - hi) / 19.0
    idx = {a: i for i, a in enumerate(AMINO_ACIDS)}
    P = np.full((len(protein), 20), lo)
    for i, a in enumerate(protein):
        P[i, idx[a]] = hi
    return FixedProfileCritic(P)


def _omega_curve():
    critic = _peaked(UBIQUITIN, hi=0.9)
    betas = np.array([0.0, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 9.0, 13.0, 18.0])
    omegas = np.array([CodonSelection(critic, beta=float(b)).dnds(UBIQUITIN) for b in betas])
    return betas, omegas


CBW, CBH = 132, 62


def _codon_box(d, x, y, codon, aa, tag, fill):
    d.append(draw.Rectangle(x, y, CBW, CBH, rx=8, fill="white", stroke=fill, stroke_width=2.4))
    d.append(draw.Text(codon, FS_LABEL, x + CBW / 2, y + 26, font_family=FONT, text_anchor="middle",
                       fill=INK, font_weight="bold"))
    d.append(draw.Text(f"{aa} ({tag})", FS_TICK, x + CBW / 2, y + 48, font_family=FONT,
                       text_anchor="middle", fill=fill))


def _arrow(d, x0, y0, x1, y1, label, color):
    d.append(draw.Line(x0, y0, x1, y1, stroke=color, stroke_width=2.4))
    ang = np.arctan2(y1 - y0, x1 - x0)
    for da in (0.5, -0.5):
        d.append(draw.Line(x1, y1, x1 - 13 * np.cos(ang + da), y1 - 13 * np.sin(ang + da),
                           stroke=color, stroke_width=2.4))
    d.append(draw.Text(label, FS_TICK, (x0 + x1) / 2, (y0 + y1) / 2 - 10, font_family=FONT,
                       text_anchor="middle", fill=color))


def left_panel(d, ox, oy, pw, ph):
    d.append(draw.Text("Codon substitution = mutation x selection", FS_LABEL, ox + pw / 2, oy - 8,
                       font_family=FONT, text_anchor="middle", font_weight="bold", fill=INK))
    cx, cy = ox, oy + (ph - CBH) / 2 - 10                       # current codon, left, vertically centred
    tx = ox + pw - CBW                                          # target column, right
    syn_y, non_y = oy + 24, oy + ph - CBH - 40
    _codon_box(d, cx, cy, "ACA", "Thr", "current", MUTED)
    _codon_box(d, tx, syn_y, "ACG", "Thr", "same aa", GOOD)
    _codon_box(d, tx, non_y, "GCA", "Ala", "new aa", SCRUT)
    _arrow(d, cx + CBW + 4, cy + 12, tx - 6, syn_y + CBH / 2, "1 nt change", GOOD)
    _arrow(d, cx + CBW + 4, cy + CBH - 12, tx - 6, non_y + CBH / 2, "1 nt change", SCRUT)
    d.append(draw.Text("synonymous: h = 1 (neutral)", FS_TICK, tx + CBW / 2, syn_y - 12,
                       font_family=FONT, text_anchor="middle", fill=GOOD))
    d.append(draw.Text("non-synonymous: h(beta) from ESM2", FS_TICK, tx + CBW / 2, non_y + CBH + 22,
                       font_family=FONT, text_anchor="middle", fill=SCRUT))
    d.append(draw.Text("mutation: neutral DNA model    selection: ESM2 on the amino acid",
                       FS_TICK, ox, oy + ph + 20, font_family=FONT, text_anchor="start",
                       fill=MUTED, font_style="italic"))


def right_panel(d, ox, oy, pw, ph, betas, omegas):
    d.append(draw.Text("Emergent dN/dS as selection strengthens", FS_LABEL, ox + pw / 2, oy - 8,
                       font_family=FONT, text_anchor="middle", font_weight="bold", fill=INK))
    bmax = float(betas.max())
    x_at = lambda b: ox + (b / bmax) * pw            # noqa: E731
    y_at = lambda w: oy + ph - w * ph                # noqa: E731 (omega in [0,1])
    # axes
    d.append(draw.Line(ox, oy, ox, oy + ph, stroke=INK, stroke_width=1.5))
    d.append(draw.Line(ox, oy + ph, ox + pw, oy + ph, stroke=INK, stroke_width=1.5))
    for w in (0.0, 0.25, 0.5, 0.75, 1.0):
        yy = y_at(w)
        d.append(draw.Line(ox - 5, yy, ox, yy, stroke=INK, stroke_width=1.5))
        d.append(draw.Text(f"{w:.2f}", FS_TICK, ox - 10, yy + 6, font_family=FONT,
                           text_anchor="end", fill=MUTED))
    for b in (0, 5, 10, 15):
        xx = x_at(b)
        d.append(draw.Line(xx, oy + ph, xx, oy + ph + 5, stroke=INK, stroke_width=1.5))
        d.append(draw.Text(str(b), FS_TICK, xx, oy + ph + 22, font_family=FONT,
                           text_anchor="middle", fill=MUTED))
    d.append(draw.Text("beta (selection strength)", FS_LABEL, ox + pw / 2, oy + ph + 48,
                       font_family=FONT, text_anchor="middle", fill=INK))
    d.append(draw.Text("dN/dS (omega)", FS_LABEL, ox - 62, oy + ph / 2, font_family=FONT,
                       text_anchor="middle", fill=INK, transform=f"rotate(-90,{ox - 62},{oy + ph / 2})"))
    # neutral reference at omega = 1
    d.append(draw.Line(ox, y_at(1.0), ox + pw, y_at(1.0), stroke=MUTED, stroke_width=1.2,
                       stroke_dasharray="4,4"))
    d.append(draw.Text("omega = 1 (neutral)", FS_TICK, ox + pw - 6, y_at(1.0) + 18, font_family=FONT,
                       text_anchor="end", fill=MUTED, font_style="italic"))
    # the curve
    pts = [(x_at(b), y_at(w)) for b, w in zip(betas, omegas)]
    d.append(draw.Lines(*[c for xy in pts for c in xy], close=False, fill="none",
                        stroke=ACCENT, stroke_width=3.2))
    for (px, py) in pts:
        d.append(draw.Circle(px, py, 3.4, fill=ACCENT, stroke="white", stroke_width=0.8))
    # annotate beta=1
    i1 = int(np.argmin(np.abs(betas - 1.0)))
    px, py = x_at(betas[i1]), y_at(omegas[i1])
    d.append(draw.Circle(px, py, 5.5, fill="none", stroke=INK, stroke_width=1.6))
    d.append(draw.Text(f"beta = 1: omega ~ {omegas[i1]:.2f}", FS_TICK, px + 12, py - 10,
                       font_family=FONT, text_anchor="start", fill=INK))


def main():
    betas, omegas = _omega_curve()
    d = draw.Drawing(W, H, origin=(0, 0))
    d.append(draw.Rectangle(0, 0, W, H, fill="white"))
    d.append(draw.Text("Selection on coding DNA yields emergent dN/dS", FS_TITLE, W / 2, 44,
                       font_family=FONT, text_anchor="middle", font_weight="bold", fill=INK))
    left_panel(d, 80, 130, 470, 330)
    right_panel(d, 720, 130, 380, 300, betas, omegas)
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "selection_dnds.svg").write_text(d.as_svg(), encoding="utf-8")
    cairosvg.svg2png(bytestring=d.as_svg().encode(), write_to=str(OUT / "selection_dnds.png"),
                     scale=200 / 72.0)
    print(f"wrote {OUT}/selection_dnds.svg / .png ; omega range {omegas.max():.3f}..{omegas.min():.3f}")


if __name__ == "__main__":
    main()
