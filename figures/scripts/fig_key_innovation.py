"""Figure: genes:species -- key-innovation diversification.

A binary "driver" gene family sets each lineage's speciation rate: carriers (D+) branch faster than
non-carriers (D-). The driver is gained by origination and, crucially, by frequency-dependent
transfer (donated more often when more lineages already carry it), and lost at some rate -- so the
gene content and the tree must grow together. This edge produces the tree.

  * Panel A (the model): two states, D- and D+; the fork under each is its speciation rate (D+
    heavier). The gain arrow is transfer + origination, the return arrow is loss.
  * Panel B (a realization): branches heavy where the lineage carries the driver, light where it
    does not; a small + marks where a lineage gains the driver and a - where it loses it. The
    carrier (heavy) clades are the speciose ones -- a genomic cause of a diversification-rate shift.

House style: B&W, one centered title, ASCII text.

Run:  /Users/aadria/miniconda3/bin/python figures/scripts/fig_key_innovation.py
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import cairosvg
import drawsvg as draw

from zombi2.coevolve import GeneDiversification, simulate_gene_diversification

from fig_sse import spec_fork
from fig_trait_pagel import curved_arrow, rate_width, state_node, _layout
from model_common import zombi_to_ete3
from zombi_style import (FONT, INK, MUTED, STATE_ON, STATE_OFF,
                         FS_TITLE, FS_LABEL, FS_ANNOT, FS_TICK)

OUT_DIR = Path(__file__).resolve().parent.parent

W, H = 1200, 700
GREY = "#9a9a9a"

# Colour version (default -> key_innovation.svg, embedded) + preserved B&W (key_innovation_bw.svg).
# ON = carries the driver (heavy branch / filled chip), OFF = no driver; swapped by render(bw=...).
ON_COL, OFF_COL = STATE_ON, STATE_OFF


def chip(d, cx, cy, on, s=13):
    """Tip chip: filled = driver present, open = absent.  Colour-aware."""
    d.append(draw.Rectangle(cx - s, cy - s, 2 * s, 2 * s,
                            fill=ON_COL if on else "white", stroke=ON_COL, stroke_width=1.6))


def event_marker(d, cx, cy, gain, r=10.0):
    """A gain (+) / loss (-) marker: an open circle with the symbol drawn as strokes so it is
    geometrically centred (SVG text glyphs for '+'/'-' do not centre reliably)."""
    d.append(draw.Circle(cx, cy, r, fill="white", stroke=INK, stroke_width=1.8))
    a = r * 0.52                                            # half-length of each symbol arm
    d.append(draw.Line(cx - a, cy, cx + a, cy, stroke=INK, stroke_width=2.2))   # horizontal (both)
    if gain:
        d.append(draw.Line(cx, cy - a, cx, cy + a, stroke=INK, stroke_width=2.2))  # vertical (+)


N_TIPS = 16

L0, DRIVER_SPEC, MU, TRANSFER, LOSS, ORIG = 0.8, 1.6, 0.2, 1.2, 0.18, 0.09


# --------------------------------------------------------------------------- panel A: the model
def panel_model(d, cx0, cy0):
    gx = 250
    d.append(draw.Text("A   the model", FS_LABEL, cx0 - 60, cy0 - 150, font_family=FONT,
                       text_anchor="start", fill=INK, font_weight="bold"))
    P = {"D-": (cx0, cy0), "D+": (cx0 + gx, cy0)}
    curved_arrow(d, P["D-"], P["D+"], +1, 22, rate_width(TRANSFER), "gain")   # transfer + origination
    curved_arrow(d, P["D+"], P["D-"], -1, 22, rate_width(LOSS), "loss")
    for lab, (x, y) in P.items():
        state_node(d, x, y, lab)
    spec_fork(d, cx0, cy0 + 46, L0)
    spec_fork(d, cx0 + gx, cy0 + 46, L0 * 3.0)          # carriers branch faster
    d.append(draw.Text("slow speciation", FS_TICK, cx0, cy0 + 120, font_family=FONT,
                       text_anchor="middle", fill=INK))
    d.append(draw.Text("fast speciation", FS_TICK, cx0 + gx, cy0 + 120, font_family=FONT,
                       text_anchor="middle", fill=INK))
    d.append(draw.Text("gain = transfer (proportional to prevalence) + origination", FS_TICK,
                       cx0 + gx / 2, cy0 - 96, font_family=FONT, text_anchor="middle", fill=MUTED,
                       font_style="italic"))
    d.append(draw.Text("carriers speciate faster, so the driver's clades grow", FS_ANNOT,
                       cx0 + gx / 2, cy0 + 158, font_family=FONT, text_anchor="middle", fill=INK,
                       font_style="italic"))


# --------------------------------------------------------------------------- panel B: a realization
def _gains(ete, has):
    """Branches where a lineage acquired the driver (parent absent, child present)."""
    return sum(1 for n in ete.traverse()
               if not n.is_root() and has[n.name] and not has[n.up.name])


def _pick(seed_range):
    # driver absent at the root: it must originate, then spread by transfer -- so a gain is on show
    model = GeneDiversification(1, lambda0=L0, mu0=MU, driver_speciation=DRIVER_SPEC,
                                transfer=TRANSFER, loss=LOSS, origination=ORIG, root_drivers=0)
    best = None
    for s in seed_range:
        res = simulate_gene_diversification(model, n_tips=N_TIPS, seed=s)
        ete = zombi_to_ete3(res.tree)
        leaves = ete.get_leaves()
        if len(leaves) > 30:
            continue
        has = {node.name: (0 in fs) for node, fs in res.node_drivers.items()}
        extant = [n for n in leaves if n.is_extant]
        prev = sum(1 for n in extant if has[n.name]) / max(1, len(extant))
        if not (0.45 <= prev <= 0.82 and _gains(ete, has) == 1):   # exactly one key innovation
            continue
        tfo, present, _, _ = _layout(ete)
        gain_t = min(tfo[n.name] for n in ete.traverse()
                     if not n.is_root() and has[n.name] and not has[n.up.name])
        # one gain, prevalence near 0.62, and a mild nudge toward an earlier gain (room to radiate)
        score = -abs(prev - 0.62) - 0.3 * (gain_t / present)
        if best is None or score > best[0]:
            best = (score, s, res, ete, has, prev)
    if best is None:                                    # fall back: whatever spread most
        for s in seed_range:
            res = simulate_gene_diversification(model, n_tips=N_TIPS, seed=s)
            ete = zombi_to_ete3(res.tree)
            if len(ete.get_leaves()) <= 28:
                has = {node.name: (0 in fs) for node, fs in res.node_drivers.items()}
                extant = [n for n in ete.get_leaves() if n.is_extant]
                prev = sum(1 for n in extant if has[n.name]) / max(1, len(extant))
                return (0.0, s, res, ete, has, prev)
    return best


def panel_realization(d, ox, oy, pw, ph):
    _, seed, res, ete, has, prev = _pick(range(1, 260))

    tfo, present, ys, nleaf = _layout(ete)
    x_at = lambda t: ox + 40 + (t / present) * (pw - 150)      # noqa: E731
    y_at = lambda k: oy + 40 + (k / max(1, nleaf - 1)) * (ph - 96)   # noqa: E731

    d.append(draw.Text("B   a simulated realization", FS_LABEL, ox, oy - 6, font_family=FONT,
                       text_anchor="start", fill=INK, font_weight="bold"))

    for n in ete.traverse():
        if n.is_root():
            continue
        y = y_at(ys[n.name])
        on = has[n.name]
        d.append(draw.Line(x_at(tfo[n.up.name]), y, x_at(tfo[n.name]), y,
                           stroke=ON_COL if on else OFF_COL, stroke_width=5.2 if on else 2.4,
                           stroke_linecap="butt"))
        if has[n.name] != has[n.up.name]:                       # gain / loss at this branch
            mx = x_at(tfo[n.up.name])
            event_marker(d, mx + 13, y, has[n.name])
        if n.is_leaf() and not n.is_extant:
            xx, ah = x_at(tfo[n.name]), 5.5
            d.append(draw.Line(xx - ah, y - ah, xx + ah, y + ah, stroke=INK, stroke_width=2.0))
            d.append(draw.Line(xx - ah, y + ah, xx + ah, y - ah, stroke=INK, stroke_width=2.0))
    for n in ete.traverse("postorder"):
        if not n.is_leaf():
            x = x_at(tfo[n.name])
            yy = [y_at(ys[c.name]) for c in n.children]
            d.append(draw.Line(x, min(yy), x, max(yy), stroke=INK, stroke_width=2.4))

    colc = ox + pw - 30
    d.append(draw.Text("driver", FS_TICK, colc, oy + 8, font_family=FONT, text_anchor="middle",
                       fill=INK, font_weight="bold"))
    for n in ete.get_leaves():
        if n.is_extant:
            chip(d, colc, y_at(ys[n.name]), has[n.name], s=9)
    return seed


def render(bw=False):
    global ON_COL, OFF_COL
    ON_COL, OFF_COL = (INK, GREY) if bw else (STATE_ON, STATE_OFF)

    d = draw.Drawing(W, H, origin=(0, 0))
    d.append(draw.Rectangle(0, 0, W, H, fill="white"))
    d.append(draw.Text("Key-innovation diversification (genes to species)",
                       FS_TITLE, W / 2, 46, font_family=FONT, text_anchor="middle",
                       font_weight="bold", fill=INK))
    ly = 82
    d.append(draw.Line(W / 2 - 300, ly, W / 2 - 268, ly, stroke=ON_COL, stroke_width=5.2))
    d.append(draw.Text("carries the driver", FS_TICK, W / 2 - 260, ly, font_family=FONT,
                       text_anchor="start", dominant_baseline="central", fill=INK))
    d.append(draw.Line(W / 2 - 90, ly, W / 2 - 58, ly, stroke=OFF_COL, stroke_width=2.4))
    d.append(draw.Text("no driver", FS_TICK, W / 2 - 50, ly, font_family=FONT,
                       text_anchor="start", dominant_baseline="central", fill=INK))
    event_marker(d, W / 2 + 78, ly, True)
    d.append(draw.Text("gain / loss", FS_TICK, W / 2 + 96, ly, font_family=FONT,
                       text_anchor="start", dominant_baseline="central", fill=INK))

    panel_model(d, 210, 320)
    seed = panel_realization(d, 560, 150, 600, 470)

    name = "key_innovation"
    suffix = "_bw" if bw else ""
    out = OUT_DIR / name
    out.mkdir(parents=True, exist_ok=True)
    (out / f"{name}{suffix}.svg").write_text(d.as_svg(), encoding="utf-8")
    cairosvg.svg2png(bytestring=d.as_svg().encode(), write_to=str(out / f"{name}{suffix}.png"),
                     scale=300 / 72.0)
    print(f"wrote {out}/{name}{suffix}.svg / .png  (tree seed {seed})")


if __name__ == "__main__":
    render(bw=False)   # colour -> key_innovation.svg (embedded)
    render(bw=True)    # preserved B&W -> key_innovation_bw.svg
