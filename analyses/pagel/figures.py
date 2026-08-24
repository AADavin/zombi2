#!/usr/bin/env python3
"""The paper figure: four panels, composed. Reads ``data/`` and ``results.json`` beside it,
resimulates one replicate for panel A, and writes ``figures/casestudy_pagel.png``.

    python figures.py

Panel A is one replicate of the feedback experiment (seed 8), drawn radially with
Phylustrator: branches painted by habitat, arrowheads at the habitat switches, an inner
ring for the presence of family A and an outer ring for genome size. Panel B is the two
connections in the gallery's diagram language. Panel C is genome size by habitat across
the four experiments, from the ``genes`` column of ``data/states``. Panel D is the share
of replicates on which Pagel's test is significant, from ``results.json``.

Needs Phylustrator (the plotting library built to read ZOMBI2's output) and, for panel B,
the repo's ``gallery/helpers.py``; run it from a repo checkout.
"""
from __future__ import annotations

import csv
import json
import pathlib
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = pathlib.Path(__file__).parent
FIG = HERE / "figures"
sys.path.insert(0, str(HERE.parent.parent / "gallery"))

PAL = {"aerobic": "#C44E52", "anaerobic": "#4C72B0"}
INK, MUTED, FAINT = "#1a1a1a", "#8a8a8a", "#c9c9c9"
EXPERIMENTS = [("feedback", "both\nconnections"), ("trait2gen", "connection 1"),
               ("gen2trait", "connection 2"), ("null", "no\nconnections")]
SEED_A = 8      # the replicate panel A shows

plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False,
                     "savefig.bbox": "tight"})


def panel_a(out: pathlib.Path) -> None:
    """One replicate of the feedback experiment, radial, with the two tip rings."""
    import phylustrator as ph
    from zombi2 import joint, traits
    from zombi2.genomes import family, genome as genome_spec
    from zombi2.params import PerCopy, PerLineage
    from zombi2.species import simulate_species_tree

    from experiment import ARMS, LOSS_BASE, N_EXTANT, SWITCH_BACK, SWITCH_BASE

    L, S = ARMS["feedback"]
    ct = simulate_species_tree(birth=1.0, n_extant=N_EXTANT, seed=SEED_A).complete_tree
    r = joint.simulate(
        genome_spec(duplication=0.05, origination=8.0, initial_families=40,
                    loss=PerCopy(LOSS_BASE).scaled_by("trait", {"anaerobic": L, "aerobic": 1.0}),
                    families=[family("A")]),
        traits.discrete(states=["aerobic", "anaerobic"], start="aerobic",
                        switch={"aerobic->anaerobic": PerLineage(SWITCH_BASE).scaled_by(
                                    "genomes:A", {"present": 1.0, "absent": S}),
                                "anaerobic->aerobic": SWITCH_BACK}),
        tree=ct, seed=SEED_A)
    lab, tips = ct.labels(), sorted(ct.extant_leaves())
    A = r.genome.family_names["A"]
    history = {lab[i]: segs for i, segs in r.trait.history.items()}
    switches = [{"kind": f"to {e.to_state}", "node": lab[e.lineage], "x": float(e.time)}
                for e in r.trait.events if e.kind != "initial" and lab.get(e.lineage)]
    presence = {lab[n]: ("present" if r.genome.family_counts(n)[A] > 0 else "absent")
                for n in tips}
    sizes = {lab[n]: len(r.genome.genomes[lab[n]]) for n in tips}

    fig = (ph.trees.plot(ph.trees.loads(ct.to_newick()), layout="radial", skeleton=False,
                         style=ph.Style(width=1250, height=1250, margin=88, branch_width=3.9))
           + ph.trees.color_history(history, palette=PAL)
           + ph.trees.branch_events(switches,
                                    styles={f"to {s}": ("triangle_right", c)
                                            for s, c in PAL.items()},
                                    size=7.5, legend=False)
           + ph.trees.ring(presence, palette={"present": "#111111", "absent": "#ffffff"},
                           gap=14, width=16, edge="#666666", edge_width=1.8)
           + ph.trees.ring(sizes, cmap="cividis", gap=38, width=16)
           + ph.trees.colorbar("number of genes", loc="top-left", width=196, height=14,
                               size=20, inset=24)
           + ph.trees.legend(entries={"aerobic": PAL["aerobic"], "anaerobic": PAL["anaerobic"]},
                             size=20, dy=86, inset=24)
           + ph.trees.legend(entries={"A present": "#111111", "A absent": "#ffffff"},
                             size=20, dy=158, inset=24))
    fig.save(str(out))


def panel_b(out: pathlib.Path) -> None:
    """The two connections, numbered the way panels C and D name them."""
    import helpers as h

    h.joint_png(str(out), [
        (("traits", "habitat", []), ("genomes", "loss rate", []),
         "(1) in the anaerobic habitat, genes are lost 5x faster"),
        (("genomes", "family A", []), ("traits", "switch rate", []),
         "(2) with no copy of A, turns anaerobic 12x faster"),
    ], frame=None, keys=False, sentence_size=15.5, wrap_width=230, font_scale=1.18, pitch=152)


def panel_c(out: pathlib.Path) -> None:
    """Genome size at the tips by habitat, per experiment, from the states tables."""
    sizes = {arm: {"aerobic": [], "anaerobic": []} for arm, _ in EXPERIMENTS}
    for arm, _ in EXPERIMENTS:
        for f in sorted((HERE / "data" / "states").glob(f"{arm}_r*.tsv")):
            for row in csv.DictReader(open(f), delimiter="\t"):
                sizes[arm][row["habitat"]].append(int(row["genes"]))
    fig, ax = plt.subplots(figsize=(7.2, 3.6))
    x = np.arange(len(EXPERIMENTS)) * 1.15
    for off, hab in ((-0.19, "aerobic"), (0.19, "anaerobic")):
        data = [sizes[a][hab] for a, _ in EXPERIMENTS]
        bp = ax.boxplot(data, positions=x + off, widths=0.3, patch_artist=True,
                        showfliers=False, medianprops={"color": "#111111", "lw": 1.4})
        for b in bp["boxes"]:
            b.set_facecolor(PAL[hab]); b.set_alpha(0.85); b.set_edgecolor("#333333")
        for w in bp["whiskers"] + bp["caps"]:
            w.set_color("#555555")
    ax.set_xticks(x); ax.set_xticklabels([]); ax.tick_params(axis="x", length=0)
    ax.set_ylabel("gene number")
    handles = [plt.Rectangle((0, 0), 1, 1, facecolor=PAL[h], alpha=0.85, edgecolor="#333333")
               for h in ("aerobic", "anaerobic")]
    ax.legend(handles, ["aerobic tips", "anaerobic tips"], frameon=False, loc="lower right")
    ax.set_xlim(-0.6, x[-1] + 0.6)
    fig.savefig(out, dpi=300)
    plt.close(fig)


def panel_d(out: pathlib.Path) -> None:
    """The share of replicates on which Pagel's test is significant, with Wilson intervals."""
    d = json.loads((HERE / "results.json").read_text())
    fig, ax = plt.subplots(figsize=(7.2, 3.6))
    x = np.arange(len(EXPERIMENTS)) * 1.15
    for off, pair, color, label in ((-0.19, "A", INK, "habitat & family A"),
                                    (0.19, "ctrl", FAINT, "habitat & control")):
        rates = [d[a][pair]["rate"] for a, _ in EXPERIMENTS]
        los = [d[a][pair]["rate"] - d[a][pair]["wilson95"][0] for a, _ in EXPERIMENTS]
        his = [d[a][pair]["wilson95"][1] - d[a][pair]["rate"] for a, _ in EXPERIMENTS]
        ax.bar(x + off, rates, width=0.36, color=color, label=label,
               yerr=np.array([los, his]), error_kw={"ecolor": MUTED, "capsize": 3, "lw": 1.1})
    ax.axhline(0.05, color=MUTED, lw=0.9, ls=(0, (3, 3)))
    ax.set_xticks(x)
    ax.set_xticklabels([lab for _, lab in EXPERIMENTS], fontsize=9.5)
    ax.set_ylabel("Pagel's test significant\n(share of replicates)")
    ax.set_ylim(0, 1.0)
    ax.legend(frameon=False, loc="upper right")
    ax.set_xlim(-0.6, x[-1] + 0.6)
    fig.savefig(out, dpi=300)
    plt.close(fig)


def compose(out: pathlib.Path, a, b, c, d) -> None:
    """Tree on the left; diagram, sizes and test stacked on the right, letters in reading order."""
    from PIL import Image, ImageChops, ImageDraw, ImageFont

    def trim(im, bg=(255, 255, 255)):
        diff = ImageChops.difference(im.convert("RGB"), Image.new("RGB", im.size, bg))
        return im.crop(diff.getbbox())

    tree, diag, cS, dS = (trim(Image.open(p)) for p in (a, b, c, d))
    W, T = 2500, 1380
    tree_x, tree_y = 100, 96
    right_x = tree_x + T + 90
    rw = W - right_x - 30
    diag = diag.resize((rw, int(diag.height * rw / diag.width)))
    cS = cS.resize((rw, int(cS.height * rw / cS.width)))
    dS = dS.resize((rw, int(dS.height * rw / dS.width)))
    tree = tree.resize((T, T))
    gap = max((T - diag.height - cS.height - dS.height) // 2, 16)
    H = tree_y + max(T, diag.height + cS.height + dS.height + 2 * gap) + 40
    canvas = Image.new("RGB", (W, H), "white")
    canvas.paste(tree, (tree_x, tree_y))
    y, ys = tree_y, []
    for panel in (diag, cS, dS):
        canvas.paste(panel, (right_x, y))
        ys.append(y)
        y += panel.height + gap
    draw = ImageDraw.Draw(canvas)
    try:
        font = ImageFont.truetype("/System/Library/Fonts/Helvetica.ttc", 60)
    except OSError:
        font = ImageFont.load_default()
    draw.text((16, 16), "A", fill="#111111", font=font)
    for letter, yy in zip("BCD", ys):
        draw.text((right_x - 66, yy + (6 if letter == "B" else 0)), letter,
                  fill="#111111", font=font)
    canvas.save(out)


def main() -> int:
    FIG.mkdir(exist_ok=True)
    a, b = FIG / "panel_a.png", FIG / "panel_b.png"
    c, d = FIG / "panel_c.png", FIG / "panel_d.png"
    panel_a(a)
    panel_b(b)
    panel_c(c)
    panel_d(d)
    out = FIG / "casestudy_pagel.png"
    compose(out, a, b, c, d)
    print(" ", out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
