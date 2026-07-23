"""Figure: sampling and fossils — the two ways your data is less than the whole tree.

Chapter 3's sampling section adds two knobs that do not change how the tree grows, only how much
of it you see. ``sampling`` keeps a fraction of the living tips; the rest are alive today but absent
from your data. ``fossils`` does the opposite, recovering dated samples scattered along **every**
branch of the complete tree (a survivor's branch as readily as an extinct one). This figure grows
one complete tree with the current engine and shows every lineage's fate at once:

  * sampled species     solid, reaches the present            -> in your data
  * unsampled (extant)   dashed, reaches the present, open ring -> alive today, not in your data
  * extinct              dashed, stops short of the present    -> died, never seen
  * fossil               solid black diamond on a branch       -> a dated past sample, on any branch

So the solid tree plus the diamonds is the data; the dashed lineages are the dark tree behind it —
the whole point of ZOMBI2's forward engine, which keeps what a backward simulation never sees. The
tree is **simulated here**, and the fates and fossils are read straight off the run.

House style: B&W, ASCII text, no title inside the figure — the manual captions it.

Run:  python figures/scripts/fig_sampling_fossils.py
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import drawsvg as draw

from zombi2 import species
from zombi_style import save, FONT, INK, MUTED, FS_LABEL, FS_TICK

# One run with both knobs on: sampling 60% of the survivors, a steady fossil rate. The seed gives a
# ~30-tip tree with a healthy mix of all three fates and a handful of fossils on different branches.
BIRTH, DEATH = 1.0, 0.35
SAMPLING, FOSSILS = 0.6, 0.15
TOTAL_TIME, SEED = 4.0, 31

W = 1400
XL, XR = 116, 1300
TREE_TOP = 210
DY = 24.0                 # vertical spacing per tip (spacious: ~30 tips)
DASH = "6,5"


def _has_sampled(nodes: dict) -> dict:
    """Map id -> True if the node or a descendant is a sampled survivor (fate 'extant').
    That subtree is the *sampled* tree, drawn solid; everything else is the dashed dark tree."""
    keep: dict[int, bool] = {}

    def rec(i: int) -> bool:
        n = nodes[i]
        # list-then-any, never a short-circuiting generator (it would skip later children).
        s = (n.fate == "extant") if not n.children else any([rec(c) for c in n.children])
        keep[i] = s
        return s

    rec(next(i for i, n in nodes.items() if n.parent is None))
    return keep


def _order(nodes: dict):
    """Ladderised leaf order; internal nodes take the mean y of their children."""
    y: dict[int, float] = {}
    counter = [0]

    def size(i: int) -> int:
        n = nodes[i]
        return 1 if not n.children else sum(size(c) for c in n.children)

    def rec(i: int) -> None:
        n = nodes[i]
        if not n.children:
            y[i] = counter[0]
            counter[0] += 1
        else:
            for c in sorted(n.children, key=size):
                rec(c)
            y[i] = sum(y[c] for c in n.children) / len(n.children)

    rec(next(i for i, n in nodes.items() if n.parent is None))
    return y, counter[0]


def _diamond(d, cx, cy, r):
    d.append(draw.Lines(cx, cy - r, cx + r, cy, cx, cy + r, cx - r, cy, close=True,
                        fill=INK, stroke="white", stroke_width=1.1))


def render() -> None:
    result = species.simulate_species_tree(birth=BIRTH, death=DEATH, sampling=SAMPLING,
                                           fossils=FOSSILS, total_time=TOTAL_TIME, seed=SEED)
    nodes = result.complete_tree.nodes
    present = max(n.end_time for n in nodes.values())
    keep = _has_sampled(nodes)
    yrank, nleaf = _order(nodes)

    H = int(TREE_TOP + DY * (nleaf - 1) + 150)
    top = TREE_TOP

    def X(t):
        return XL + t / present * (XR - XL)

    def Y(i):
        return top + yrank[i] * DY

    d = draw.Drawing(W, H, origin=(0, 0))
    d.append(draw.Rectangle(0, 0, W, H, fill="white"))

    # faint present reference behind the tree
    ys = [Y(i) for i, n in nodes.items() if not n.children]
    d.append(draw.Line(X(present), min(ys) - 16, X(present), max(ys) + 16,
                       stroke="#cccccc", stroke_width=1.0, stroke_dasharray="2,4"))

    def seg(x1, y1, x2, y2, dashed):
        if dashed:
            d.append(draw.Line(x1, y1, x2, y2, stroke=INK, stroke_width=2.0,
                               stroke_dasharray=DASH, stroke_linecap="butt"))
        else:
            d.append(draw.Line(x1, y1, x2, y2, stroke=INK, stroke_width=2.6, stroke_linecap="round"))

    for i, n in nodes.items():
        px = X(0) - 14 if n.parent is None else X(nodes[n.parent].end_time)
        seg(px, Y(i), X(n.end_time), Y(i), not keep[i])
        for c in (n.children or ()):
            seg(X(n.end_time), Y(i), X(n.end_time), Y(c), not keep[c])
        if n.fate == "unsampled":                    # alive today but not in the sample: open ring
            d.append(draw.Circle(X(present), Y(i), 4.6, fill="white", stroke=INK, stroke_width=2.0))

    for nid, t in result.fossils:                    # fossils, scattered along the branches
        _diamond(d, X(t), Y(nid), 6.2)

    # legend: a single vertical column, bottom-left, clear of the tree
    row = FS_LABEL * 1.85
    lx = 40
    ly = max(ys) - 3 * row          # anchor the four-row block at the bottom-left

    def label(y, txt):
        d.append(draw.Text(txt, FS_LABEL, lx + 58, y, font_family=FONT, text_anchor="start",
                           dominant_baseline="central", fill=INK))

    d.append(draw.Line(lx, ly, lx + 40, ly, stroke=INK, stroke_width=2.6, stroke_linecap="round"))
    label(ly, "sampled species (in your data)")
    y2 = ly + row
    d.append(draw.Line(lx, y2, lx + 40, y2, stroke=INK, stroke_width=2.0, stroke_dasharray=DASH,
                       stroke_linecap="butt"))
    d.append(draw.Circle(lx + 40, y2, 4.6, fill="white", stroke=INK, stroke_width=2.0))
    label(y2, "unsampled - alive today, not sampled")
    y3 = y2 + row
    d.append(draw.Line(lx, y3, lx + 40, y3, stroke=INK, stroke_width=2.0, stroke_dasharray=DASH,
                       stroke_linecap="butt"))
    label(y3, "extinct - died before the present")
    y4 = y3 + row
    _diamond(d, lx + 20, y4, 6.2)
    label(y4, "fossil - a dated sample on a branch")

    # time axis
    ya = max(ys) + 44
    for k in range(5):
        t = present * k / 4
        d.append(draw.Line(X(t), ya, X(t), ya + 6, stroke=INK, stroke_width=1.2))
        d.append(draw.Text(f"{t:.0f}" if t == int(t) else f"{t:.1f}", FS_TICK, X(t), ya + 26,
                           font_family=FONT, text_anchor="middle", fill=MUTED))
    d.append(draw.Text("time (origin to present)", FS_LABEL, (XL + XR) / 2, ya + 56, font_family=FONT,
                       text_anchor="middle", fill=MUTED))

    save(d, "sampling_fossils")
    from collections import Counter
    f = Counter(n.fate for n in nodes.values() if not n.children)
    print(f"  ({nleaf} tips: {f.get('extant',0)} sampled, {f.get('unsampled',0)} unsampled, "
          f"{f.get('extinct',0)} extinct; {len(result.fossils)} fossils)")


if __name__ == "__main__":
    render()
