"""The README's overview figure — one tree, four levels.

    python figures/scripts/fig_overview.py        # -> assets/overview.png


Every panel is the SAME run: one 30-tip birth-death tree, its genomes, its sequences and a trait on
it. Panel A shows the complete tree (extinct lineages included — the thing ZOMBI2 keeps and an
empirical tree cannot); B, C and D show the extant tree, because that is what an analysis is handed.
Reading across the four is then reading one dataset at four levels.

  A  species  — the complete history: survivors solid, extinct lineages dashed
  B  genomes  — each tip's gene order, homologues ribboned, coloured by ancestral position
  C  sequences— one family's alignment
  D  traits   — two characters drifting together, and the correlation they leave
"""
from __future__ import annotations

import os
import subprocess
import sys

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, ROOT)
# The panels ARE gallery cards — the extinct-lineages tree, the synteny tracks, the alignment beside
# a tree, the dependent-traits pair — so this borrows their helpers rather than restating the
# recipes, which would let the README and the gallery drift apart.
sys.path.insert(0, os.path.join(ROOT, "gallery"))

import matplotlib
matplotlib.use("Agg")
import matplotlib.image as mpimg
import matplotlib.pyplot as plt

import phylustrator as ph
import helpers as h
from zombi2.traits import simulate_continuous
from zombi2.tree import read_newick

OUT = os.path.join(ROOT, "figures", "_overview")      # panels + the cached run, not committed
ASSETS = os.path.join(ROOT, "assets")                 # where the README reads the figure from
PY = "/Users/aadria/miniconda3/bin/python"
RUN = os.path.join(OUT, "_data", "overview")

W, H = 1240, 840                       # one shape for every panel, so the grid has no letterboxing
PANEL = ph.Style(width=W, height=H, margin=72, branch_width=2.4)
#: D carries its signal IN the branches, so they are drawn twice as thick — a colour
#: gradient needs area to be read, where a topology only needs a line.
PANEL_TRAIT = ph.Style(width=W, height=H, margin=72, branch_width=4.8)


def _build_run():
    """One run, four levels — the dataset every panel reads."""
    if os.path.isdir(os.path.join(RUN, "sequences", "alignments")):
        return

    def z(*args):
        subprocess.run([PY, "-m", "zombi2.cli.main", *[str(a) for a in args], "--quiet"],
                       check=True, cwd=ROOT, stdout=subprocess.DEVNULL)

    # death well below birth: enough extinct lineages for panel A to be worth drawing, and 30
    # survivors to carry the other three
    # seed 56: the most balanced 30-tip tree in the first 60 seeds (Colless 29), and 16
    # extinct lineages — enough for panel A to be worth drawing
    z("species", RUN, "--birth", 1.0, "--death", 0.45, "--n-extant", 30, "--seed", 56)
    # Gentle rearrangement on purpose: a synteny panel says something only when the blocks are
    # mostly collinear and a few have moved. Shuffled hard, every row is a fresh permutation and the
    # ribbons cross into noise, which reads as decoration rather than as data.
    z("genomes", RUN, "--resolution", "ordered", "--initial-families", 14,
      "--duplication", 0.015, "--loss", 0.015, "--transfer", 0.0, "--origination", 0.0,
      "--inversion", 0.10, "--inversion-extent", 3, "--chromosomes", 1, "--seed", 5)
    # --divergence keeps the alignment from saturating (else every column varies — no signal)
    z("sequences", RUN, "--model", "jc69", "--length", 60, "--divergence", 0.35, "--seed", 7)


def _extant():
    return ph.trees.read(os.path.join(RUN, "species", "species_extant.nwk"))


def _complete():
    with open(os.path.join(RUN, "species", "species_complete.nwk"), encoding="utf-8") as f:
        text = f.read()
    ct, _ = read_newick(text)
    return ct, text


# ── A. species: the whole history, not just the survivors ──────────────────────────────────────
def panel_species(out):
    ct, text = _complete()
    tree = ph.trees.loads(text)
    (ph.trees.plot(tree, dashed=h.dashed_extinct(tree, ct), style=PANEL)
     + ph.trees.time_axis("time", tick_size=20, label_size=26)).save(out)


# ── B. genomes: gene order at every tip ────────────────────────────────────────────────────────
def _reference_order():
    """The gene order the run started with — the ancestral arrangement every tip is read against."""
    with open(os.path.join(RUN, "genomes", "initial_genome.tsv"), encoding="utf-8") as f:
        ix = {n: k for k, n in enumerate(f.readline().rstrip("\n").split("\t"))}
        return [line.rstrip("\n").split("\t")[ix["family"]] for line in f if line.strip()]


def panel_genomes(out):
    fig = ph.trees.plot(_extant(), style=PANEL)
    genomes = ph.zombi.read_genomes(RUN)
    panel = ph.genomes.tracks(list(genomes.values()), reference=_reference_order(),
                              cmap="viridis", opacity=0.38)
    ph.beside(fig, panel, width=W, height=H, tree_fraction=0.27, gap=22, pad=34).save(out)


# ── C. sequences: the alignment the tree produced ──────────────────────────────────────────────
def _best_single_copy(M) -> str:
    """The family single-copy in the most genomes — the cleanest one-row-per-tip alignment."""
    best, best_n = M.cols[0], -1
    for j, c in enumerate(M.cols):
        n = sum(1 for row in M.values if row[j] == 1)
        if n > best_n:
            best, best_n = c, n
    return best


def panel_sequences(out):
    fam = _best_single_copy(ph.zombi.read_profiles(RUN))
    aln = ph.zombi.read_alignment(RUN, fam)
    fig = ph.trees.plot(_extant(), style=PANEL)
    ph.beside(fig, ph.genomes.alignment(aln, letters=False, legend=True),
              width=W, height=H, tree_fraction=0.27, gap=22, pad=34, footer=62).save(out)


# ── D. traits: two characters that drift together ──────────────────────────────────────────────
def panel_traits(out):
    """Two correlated continuous traits — the gallery's dependent-continuous card.

    A single diffusing character on a tree is a tree with coloured branches, which panel A already
    is; two traits that drift *together* is a different picture and a different claim. The two trees
    stack on the left, one per trait, and the tip scatter on the right is the correlation itself —
    the only element in the whole figure that is not a phylogeny.

    Composited here rather than through ``h.composite_two_trees_scatter``: that builds its own figure
    at its own aspect, and a panel of this grid has to fill the same box as the other three or it
    reads as a small picture in a large frame.
    """
    ct, _ = _complete()
    res = simulate_continuous(ct, start={"x": 0.0, "y": 0.0}, rate={"x": 1.0, "y": 1.0},
                              correlation={("x", "y"): 0.9}, seed=1)
    extant = ph.trees.read(os.path.join(RUN, "species", "species_extant.nwk"))

    pngs = []
    for key in ("x", "y"):
        png = out.replace(".png", f"_{key}.png")
        values = {f"n{i}": v[key] for i, v in res.node_values.items()}
        h.render_tree_for_composite(extant, png, values, width=1000, height=560, branch_width=4.5)
        pngs.append(png)

    tips = list(ct.extant_leaves())
    xs = [res.node_values[n.id]["x"] for n in tips]
    ys = [res.node_values[n.id]["y"] for n in tips]

    fig = plt.figure(figsize=(W / 100, H / 100), dpi=100)
    grid = fig.add_gridspec(2, 2, width_ratios=[2.15, 1.25], hspace=0.06, wspace=0.10,
                            left=0.015, right=0.965, top=0.97, bottom=0.10)
    for row, (png, name) in enumerate(zip(pngs, ("trait x", "trait y"))):
        ax = fig.add_subplot(grid[row, 0])
        ax.imshow(mpimg.imread(png))
        ax.set_axis_off()
        ax.text(0.085, 0.99, name, transform=ax.transAxes, fontsize=19,
                va="bottom", ha="left", color="#1a1a1a")
    ax = fig.add_subplot(grid[:, 1])
    ax.scatter(xs, ys, c=xs, cmap="viridis", s=64, edgecolors="white", linewidths=0.7)
    ax.set_xlabel("trait x", fontsize=18)
    ax.set_ylabel("trait y", fontsize=18)
    ax.tick_params(labelsize=14)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    fig.savefig(out, dpi=170, facecolor="white")
    plt.close(fig)


PANELS = [
    ("A", "Species trees", panel_species),
    ("B", "Genomes", panel_genomes),
    ("C", "Sequences", panel_sequences),
    ("D", "Traits", panel_traits),
]


def main(border: bool = False, name: str = "readme_overview.png"):
    os.makedirs(os.path.join(OUT, "_data"), exist_ok=True)
    os.makedirs(ASSETS, exist_ok=True)
    _build_run()
    paths = []
    for letter, title, fn in PANELS:
        p = os.path.join(OUT, f"panel_{letter}.png")
        if not os.path.isfile(p):
            fn(p)
            print("rendered", letter, title)
        paths.append(p)

    fig, axes = plt.subplots(2, 2, figsize=(15.4, 10.6))
    for ax, path, (letter, title, _fn) in zip(axes.ravel(), paths, PANELS):
        ax.imshow(mpimg.imread(path))
        if border:
            # a hairline, not a frame: enough to say where a panel ends when two of them are white
            # to the edge, and light enough not to become part of the picture
            ax.set_xticks([]); ax.set_yticks([])
            for spine in ax.spines.values():
                spine.set_visible(True)
                spine.set_linewidth(0.6)
                spine.set_color("#c9d2ce")
        else:
            ax.set_axis_off()
        ax.text(0.0, 1.012, letter, transform=ax.transAxes, fontsize=21, fontweight="bold",
                va="bottom", ha="left", color="#1a1a1a")
        ax.text(0.5, 1.012, title, transform=ax.transAxes, fontsize=19, va="bottom", ha="center",
                color="#1a1a1a")
    fig.subplots_adjust(left=0.012, right=0.988, top=0.955, bottom=0.012, wspace=0.05, hspace=0.11)
    grid = os.path.join(ASSETS, name)
    fig.savefig(grid, dpi=170, facecolor="white")
    plt.close(fig)
    print("wrote", grid)


if __name__ == "__main__":
    os.makedirs(ASSETS, exist_ok=True)
    main(border=False, name="overview.png")
