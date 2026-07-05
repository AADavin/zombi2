"""Figure: HiSSE -- hidden-state diversification, the honest null.

Diversification is driven by an *unobserved* class (a slow regime and a fast one), NOT by the
observed binary trait -- within each hidden class the two observed states speciate at the SAME rate.
So the tree has real fast/slow clades, but the observed character is spread across them and cannot
explain the diversity. This is the null a raw BiSSE fit would wrongly read as a trait effect.

  * Panel A (the model): the four joint states (observed 0/1) x (hidden slow/fast). The fork under
    each state is its speciation rate -- heavy in the fast row, light in the slow row, and identical
    across the observed columns (the observed trait is neutral for diversification). Horizontal
    arrows flip the observed trait; vertical arrows switch the hidden regime.
  * Panel B (a realization): branches heavy where the hidden class is fast (speciose, bushy clades)
    and light where it is slow; the tip chips give the OBSERVED state. The bushy clades are the
    fast-hidden ones -- the observed chips do not track them.

House style: B&W, one centered title, ASCII text.

Run:  /Users/aadria/miniconda3/bin/python figures/scripts/fig_sse_hisse.py
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import cairosvg
import drawsvg as draw

from zombi2.coevolve import BiSSE, HiSSE, simulate_sse

from fig_sse import spec_fork
from fig_trait_pagel import curved_arrow, rate_width, state_node, chip, _layout
from model_common import zombi_to_ete3
from zombi_style import FONT, INK, MUTED, FS_TITLE, FS_LABEL, FS_ANNOT, FS_TICK

OUT_DIR = Path(__file__).resolve().parent.parent

W, H = 1200, 700
GREY = "#9a9a9a"
N_TIPS = 16

L_SLOW, L_FAST, MU, QOBS, QHID = 0.7, 2.4, 0.25, 0.5, 0.25
SLOW = BiSSE(lambda0=L_SLOW, lambda1=L_SLOW, mu0=MU, mu1=MU, q01=QOBS, q10=QOBS)
FAST = BiSSE(lambda0=L_FAST, lambda1=L_FAST, mu0=MU, mu1=MU, q01=QOBS, q10=QOBS)


# --------------------------------------------------------------------------- panel A: the model
def panel_model(d, cx0, cy0):
    gx = gy = 176
    d.append(draw.Text("A   the model", FS_LABEL, cx0 - 60, cy0 - 150, font_family=FONT,
                       text_anchor="start", fill=INK, font_weight="bold"))
    P = {"0s": (cx0, cy0), "1s": (cx0 + gx, cy0),
         "0f": (cx0, cy0 + gy), "1f": (cx0 + gx, cy0 + gy)}
    B = 15
    curved_arrow(d, P["0s"], P["1s"], +1, B, rate_width(QOBS), f"{QOBS:g}")     # obs flip (slow row)
    curved_arrow(d, P["1s"], P["0s"], -1, B, rate_width(QOBS), "")
    curved_arrow(d, P["0f"], P["1f"], -1, B, rate_width(QOBS), f"{QOBS:g}")     # obs flip (fast row)
    curved_arrow(d, P["1f"], P["0f"], +1, B, rate_width(QOBS), "")
    curved_arrow(d, P["0s"], P["0f"], +1, B, rate_width(QHID), f"{QHID:g}")     # hidden switch
    curved_arrow(d, P["0f"], P["0s"], -1, B, rate_width(QHID), "")
    curved_arrow(d, P["1s"], P["1f"], -1, B, rate_width(QHID), f"{QHID:g}")
    curved_arrow(d, P["1f"], P["1s"], +1, B, rate_width(QHID), "")
    for key, (x, y) in P.items():
        state_node(d, x, y, key[0])                        # observed digit
        spec_fork(d, x, y + 34, L_SLOW if key[1] == "s" else L_FAST)   # fork width = speciation rate
    d.append(draw.Text("slow", FS_TICK, cx0 - 116, cy0, font_family=FONT, text_anchor="middle",
                       dominant_baseline="central", fill=MUTED, font_style="italic"))
    d.append(draw.Text("fast", FS_TICK, cx0 - 116, cy0 + gy, font_family=FONT, text_anchor="middle",
                       dominant_baseline="central", fill=MUTED, font_style="italic"))
    d.append(draw.Text("flip observed", FS_TICK, cx0 + gx / 2, cy0 - 96, font_family=FONT,
                       text_anchor="middle", fill=MUTED, font_style="italic"))
    d.append(draw.Text("node = observed;  fork width = speciation rate", FS_TICK, cx0 + gx / 2,
                       cy0 + gy + 108, font_family=FONT, text_anchor="middle", fill=MUTED))
    d.append(draw.Text("same speciation rate in both observed columns", FS_ANNOT, cx0 + gx / 2,
                       cy0 + gy + 138, font_family=FONT, text_anchor="middle", fill=INK,
                       font_style="italic"))


# --------------------------------------------------------------------------- panel B: a realization
def _fullmap(res):
    """{node name -> (observed 0/1, hidden 0=slow/1=fast)}."""
    return {n.name: res.full_label(i) for n, i in res.node_values.items()}


def _pick(seed_range):
    model = HiSSE(classes=[SLOW, FAST], hidden_transition=QHID)
    best = None
    for s in seed_range:
        res = simulate_sse(model, n_tips=N_TIPS, seed=s)
        fm = _fullmap(res)
        ete = zombi_to_ete3(res.tree)
        leaves = ete.get_leaves()
        if len(leaves) > 30:
            continue
        extant = [n for n in leaves if n.is_extant]
        frac_fast = sum(1 for n in extant if fm[n.name][1] == 1) / max(1, len(extant))
        frac_obs1 = sum(1 for n in extant if fm[n.name][0] == 1) / max(1, len(extant))
        # want fast to dominate the tips (bushy) but observed roughly balanced (uncorrelated)
        if not (0.55 <= frac_fast <= 0.85 and 0.3 <= frac_obs1 <= 0.7):
            continue
        score = -abs(frac_obs1 - 0.5)
        if best is None or score > best[0]:
            best = (score, s, res, ete, fm)
    return best


def panel_realization(d, ox, oy, pw, ph):
    picked = _pick(range(1, 260))
    _, seed, res, ete, fm = picked
    hid_chg = {}
    for node, t, frm, to in res.changes():
        if frm[1] != to[1]:
            hid_chg.setdefault(node.name, []).append((t, to[1]))

    tfo, present, ys, nleaf = _layout(ete)
    x_at = lambda t: ox + 40 + (t / present) * (pw - 150)      # noqa: E731
    y_at = lambda k: oy + 40 + (k / max(1, nleaf - 1)) * (ph - 96)   # noqa: E731

    d.append(draw.Text("B   a simulated realization", FS_LABEL, ox, oy - 6, font_family=FONT,
                       text_anchor="start", fill=INK, font_weight="bold"))

    def seg(x1, x2, y, fast):
        d.append(draw.Line(x1, y, x2, y, stroke=INK if fast else GREY,
                           stroke_width=5.2 if fast else 2.4, stroke_linecap="butt"))

    for n in ete.traverse():
        if n.is_root():
            continue
        y = y_at(ys[n.name])
        t0, cur = tfo[n.up.name], fm[n.up.name][1]
        for tt, to in sorted(hid_chg.get(n.name, [])):
            seg(x_at(t0), x_at(tt), y, cur == 1)
            t0, cur = tt, to
        seg(x_at(t0), x_at(tfo[n.name]), y, cur == 1)
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
    d.append(draw.Text("obs", FS_TICK, colc, oy + 8, font_family=FONT, text_anchor="middle",
                       fill=INK, font_weight="bold"))
    for n in ete.get_leaves():
        if n.is_extant:
            chip(d, colc, y_at(ys[n.name]), fm[n.name][0] == 1, s=9)
    d.append(draw.Text("bushy clades are the fast-hidden ones; the observed trait does not track them",
                       FS_TICK, ox + pw / 2, oy + ph + 6, font_family=FONT, text_anchor="middle",
                       fill=INK, font_style="italic"))
    return seed


def render():
    d = draw.Drawing(W, H, origin=(0, 0))
    d.append(draw.Rectangle(0, 0, W, H, fill="white"))
    d.append(draw.Text("Hidden-state diversification (HiSSE) -- the null",
                       FS_TITLE, W / 2, 46, font_family=FONT, text_anchor="middle",
                       font_weight="bold", fill=INK))
    ly = 82
    d.append(draw.Line(W / 2 - 320, ly, W / 2 - 288, ly, stroke=INK, stroke_width=5.2))
    d.append(draw.Text("hidden: fast", FS_TICK, W / 2 - 280, ly, font_family=FONT,
                       text_anchor="start", dominant_baseline="central", fill=INK))
    d.append(draw.Line(W / 2 - 150, ly, W / 2 - 118, ly, stroke=GREY, stroke_width=2.4))
    d.append(draw.Text("hidden: slow", FS_TICK, W / 2 - 110, ly, font_family=FONT,
                       text_anchor="start", dominant_baseline="central", fill=INK))
    d.append(draw.Rectangle(W / 2 + 40, ly - 9, 18, 18, fill=INK, stroke=INK, stroke_width=1.2))
    d.append(draw.Text("observed 1", FS_TICK, W / 2 + 66, ly, font_family=FONT,
                       text_anchor="start", dominant_baseline="central", fill=INK))

    panel_model(d, 210, 300)
    seed = panel_realization(d, 560, 150, 600, 470)

    name = "sse_hisse"
    out = OUT_DIR / name
    out.mkdir(parents=True, exist_ok=True)
    (out / f"{name}.svg").write_text(d.as_svg(), encoding="utf-8")
    cairosvg.svg2png(bytestring=d.as_svg().encode(), write_to=str(out / f"{name}.png"),
                     scale=300 / 72.0)
    print(f"wrote {out}/{name}.svg / .png  (tree seed {seed})")


if __name__ == "__main__":
    render()
