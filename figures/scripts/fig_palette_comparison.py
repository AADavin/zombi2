"""Reference sheet: the same Brownian-motion tree rendered in 8 candidate colormaps,
each with its colour bar, so a palette can be chosen for the continuous trait figures.

Run:  python figures/scripts/fig_palette_comparison.py
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))   # local: colormaps, zombi_style

import drawsvg as draw
from PIL import Image

import phylustrator as ph
from zombi2 import BirthDeath, BrownianMotion, simulate_species_tree, simulate_traits

from colormaps import CMAPS, cmap, hexc
from model_common import zombi_to_ete3
from zombi_style import species_style

OUT = Path(__file__).resolve().parent.parent / "palette_comparison"
PALS = ["viridis", "magma", "plasma", "cividis", "mako", "rocket", "coolwarm", "RdBu"]


def main():
    ztree = simulate_species_tree(BirthDeath(1.0, 0.3), n_tips=14, age=1.0, direction="backward", seed=3)
    res = simulate_traits(ztree, BrownianMotion(sigma2=2.0), seed=7)
    name2val = {n.name: float(v) for n, v in res.node_values.items()}
    vmin, vmax = min(name2val.values()), max(name2val.values())
    norm = lambda v: (v - vmin) / (vmax - vmin) if vmax > vmin else 0.5   # noqa: E731

    OUT.mkdir(parents=True, exist_ok=True)
    imgs = []
    for name in PALS:
        cm = cmap(name)
        node_to_rgb = {k: cm(norm(v)) for k, v in name2val.items()}
        tree = zombi_to_ete3(ztree)
        style = species_style(width=520, height=430, margin=38, font_size=12)
        d = ph.VerticalTreeDrawer(tree, style=style)
        d._calculate_layout()
        d.plot_continuous_variable(node_to_rgb, stroke_width=4.5)
        for i in range(60):
            d.drawing.append(draw.Rectangle(-150 + i * 5, -style.height / 2 + 34, 6, 14,
                                            fill=hexc(cm(i / 59.0))))
        d.add_title(name, font_size=22)
        p = OUT / f"_{name}.png"
        d.save_png(str(p), dpi=105)
        imgs.append(Image.open(p))

    w, h = imgs[0].size
    sheet = Image.new("RGB", (w * 2, h * 4), "white")
    for i, im in enumerate(imgs):
        sheet.paste(im, ((i % 2) * w, (i // 2) * h))
    sheet.save(OUT / "palette_comparison.png")
    for name in PALS:
        (OUT / f"_{name}.png").unlink()
    print(f"wrote {OUT / 'palette_comparison.png'}  ({', '.join(PALS)})")


if __name__ == "__main__":
    main()
