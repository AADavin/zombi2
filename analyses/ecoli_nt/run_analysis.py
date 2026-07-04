#!/usr/bin/env python3
"""Evolve the E. coli K-12 chromosome along a simulated 12-species tree (nucleotide genic model),
then compute summary statistics and render the report figures.

Outputs (in this directory):
  figs/*.pdf   — the report figures
  stats.tex    — LaTeX \\newcommand macros with every number the report cites
"""
from __future__ import annotations

import json
import os
from collections import Counter

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Patch

import zombi2 as z
from zombi2.events import EventType

HERE = os.path.dirname(os.path.abspath(__file__))
FIGS = os.path.join(HERE, "figs")
GFF = os.path.join(HERE, "ecoli.gff")

# ---------------------------------------------------------------- house style
INK, MUTED = "#1a1a1a", "#8a8a8a"
ACCENT = {"origination": "#984EA3", "duplication": "#377EB8", "transfer": "#4DAF4A",
          "loss": "#E41A1C", "speciation": "#999999", "highlight": "#FDBF6F",
          "gene": "#377EB8", "intergene": "#cfcfcf", "pseudo": "#E41A1C"}
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
    "font.size": 10, "axes.edgecolor": INK, "axes.linewidth": 0.8,
    "axes.titlesize": 11, "axes.titleweight": "bold", "xtick.color": INK,
    "ytick.color": INK, "text.color": INK, "axes.labelcolor": INK,
    "figure.dpi": 150, "savefig.bbox": "tight",
})

N_TIPS, AGE, TREE_SEED, SIM_SEED = 12, 1.0, 42, 7
RATES = dict(inversion=4e-6, loss=7e-6, duplication=3e-6, transfer=3e-6,
             transposition=1.5e-6, origination=2.0, extension=0.9996,
             pseudogenization=0.4, replacement=0.4)


# ---------------------------------------------------------------- newick draw
def parse_newick(s):
    s = s.strip().rstrip(";")
    pos = 0

    def clade():
        nonlocal pos
        node = {"name": "", "length": 0.0, "children": []}
        if pos < len(s) and s[pos] == "(":
            pos += 1
            while True:
                node["children"].append(clade())
                if s[pos] == ",":
                    pos += 1
                    continue
                if s[pos] == ")":
                    pos += 1
                    break
        start = pos
        while pos < len(s) and s[pos] not in ",();":
            pos += 1
        label = s[start:pos]
        if ":" in label:
            nm, ln = label.split(":", 1)
            node["name"] = nm
            try:
                node["length"] = float(ln)
            except ValueError:
                node["length"] = 0.0
        else:
            node["name"] = label
        return node

    return clade()


def _layout(node, x0=0.0, counter=None):
    counter = counter or [0]
    x = x0 + node["length"]
    if not node["children"]:
        node["_x"], node["_y"] = x, counter[0]
        counter[0] += 1
        return node["_y"], node["_y"]
    ys = [_layout(c, x, counter) for c in node["children"]]
    node["_x"] = x
    node["_y"] = np.mean([(lo + hi) / 2 for lo, hi in ys])
    return min(lo for lo, _ in ys), max(hi for _, hi in ys)


def draw_tree(ax, newick, *, tip_label=lambda n: n, color=INK, lw=1.6,
              tip_color=None, tip_fontsize=8, edge_event=None):
    root = parse_newick(newick)
    _layout(root)

    def walk(node):
        for c in node["children"]:
            col = color
            if edge_event:
                col = edge_event(c) or color
            ax.plot([node["_x"], c["_x"]], [c["_y"], c["_y"]], color=col, lw=lw,
                    solid_capstyle="round", zorder=2)
            walk(c)
        if node["children"]:
            ys = [c["_y"] for c in node["children"]]
            ax.plot([node["_x"], node["_x"]], [min(ys), max(ys)], color=color, lw=lw,
                    solid_capstyle="round", zorder=2)
        else:
            tc = tip_color(node["name"]) if tip_color else INK
            ax.plot(node["_x"], node["_y"], "o", ms=3.2, color=tc, zorder=3)
            ax.text(node["_x"] + 0.012 * _xmax[0], node["_y"], tip_label(node["name"]),
                    va="center", ha="left", fontsize=tip_fontsize, color=INK)
    global _xmax
    _xmax = [0]

    def maxx(n):
        _xmax[0] = max(_xmax[0], n["_x"])
        for c in n["children"]:
            maxx(c)
    maxx(root)
    walk(root)
    ax.set_xlim(-0.02 * _xmax[0], _xmax[0] * 1.28)
    ax.axis("off")
    return root


# ---------------------------------------------------------------- simulate
def main():
    os.makedirs(FIGS, exist_ok=True)
    gff = z.read_gff(GFF)
    tree = z.simulate_species_tree(z.BirthDeath(birth=1.0, death=0.35),
                                   n_tips=N_TIPS, age=AGE, seed=TREE_SEED)
    res = z.simulate_nucleotide_genomes(
        tree, root_length=gff.length, gene_intervals=gff.genes,
        seed=SIM_SEED, **RATES)

    leaves = sorted(res.leaf_genomes, key=lambda n: n.name)
    names = [n.name for n in leaves]

    # --- per-species genome composition ---
    comp = {}
    for n in leaves:
        segs = res.leaf_genomes[n]._segments
        func = [s for s in segs if s.gene_id is not None and s.is_gene]
        pseud = [s for s in segs if s.gene_id is not None and not s.is_gene]
        comp[n.name] = dict(
            size=res.leaf_genomes[n].size(),
            gene_copies=len(func),
            distinct_genes=len({s.gene_id for s in func}),
            pseudogenes=len(pseud),
        )

    # --- profile (copy number of each atom per species) + pangenome ---
    atom_ids, sp_names, M = res.profile_matrix()
    kind = {a.atom_id: a.kind for a in res.atoms}
    gid = {a.atom_id: a.gene_id for a in res.atoms}
    order = [sp_names.index(nm) for nm in names]          # column order = leaf order
    M = M[:, order]
    gene_rows = np.array([i for i, aid in enumerate(atom_ids) if kind[aid] == "gene"])
    Mg = M[gene_rows]
    prevalence = (Mg > 0).sum(axis=1)                     # in how many of 12 species present
    ancestral_genes = len(gff.genes)
    n_gene_families = Mg.shape[0]
    core = int((prevalence == N_TIPS).sum())
    unique = int((prevalence == 1).sum())
    accessory = int(((prevalence >= 1) & (prevalence < N_TIPS)).sum())
    duplicated = int((Mg.max(axis=1) > 1).sum())          # paralogs somewhere

    # --- event tally from the log ---
    ev = Counter()
    for r in res.event_log:
        e = r.event
        if e is EventType.LOSS:
            ev["pseudogenization" if (len(r.genes) == 2 and r.genes[1].role == "pseudogenized")
                else "loss"] += 1
        elif e is EventType.ORIGINATION:
            ev["origination"] += 1
        elif e is EventType.DUPLICATION:
            ev["duplication"] += 1
        elif e is EventType.TRANSFER:
            ev["transfer"] += 1
        elif e is EventType.INVERSION:
            ev["inversion"] += 1
        elif e is EventType.TRANSPOSITION:
            ev["transposition"] += 1
    pseudos = res.pseudogenizations()

    # ============================================================ FIGURES
    # Fig 1 — species tree
    fig, ax = plt.subplots(figsize=(5.0, 3.4))
    draw_tree(ax, tree.to_newick(), color=INK, lw=2.0, tip_fontsize=9)
    ax.set_title("Simulated 12-species tree", loc="left")
    fig.savefig(os.path.join(FIGS, "fig1_species_tree.pdf"))
    plt.close(fig)

    # per-species genome dynamics (relative to the ancestral E. coli gene set)
    anc_names = {nm for _a, _b, nm in gff.genes}
    anc_rows = np.array([i for i in range(len(gene_rows)) if gid[atom_ids[gene_rows[i]]] in anc_names])
    nov_rows = np.array([i for i in range(len(gene_rows)) if gid[atom_ids[gene_rows[i]]] not in anc_names])
    lost = np.array([(Mg[anc_rows, j] == 0).sum() for j in range(N_TIPS)])
    gained = np.array([(Mg[nov_rows, j] > 0).sum() if len(nov_rows) else 0 for j in range(N_TIPS)])
    dup_sp = np.array([(Mg[:, j] > 1).sum() for j in range(N_TIPS)])
    pseud_sp = np.array([comp[nm]["pseudogenes"] for nm in names])

    # Fig 2 — genome dynamics per lineage: erosion (down) vs innovation (up)
    fig, ax = plt.subplots(figsize=(6.6, 3.6))
    x = np.arange(N_TIPS)
    ax.bar(x, -lost, color=ACCENT["loss"], label="genes lost")
    ax.bar(x, -pseud_sp, bottom=-lost, color=ACCENT["highlight"], label="pseudogenised")
    ax.bar(x, gained, color=ACCENT["transfer"], label="novel genes gained")
    ax.bar(x, dup_sp, bottom=gained, color=ACCENT["duplication"], label="duplicated families")
    ax.axhline(0, color=INK, lw=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(names, rotation=45, ha="right", fontsize=8)
    ax.set_ylabel("gene families  (loss ↓ / gain ↑)")
    ax.set_title("Genome dynamics per lineage", loc="left")
    ax.legend(frameon=False, fontsize=7.5, ncol=2, loc="upper center",
              bbox_to_anchor=(0.5, 1.02))
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    fig.savefig(os.path.join(FIGS, "fig2_gene_content.pdf"))
    plt.close(fig)

    # Fig 3 — pangenome prevalence spectrum
    fig, ax = plt.subplots(figsize=(5.4, 3.3))
    counts = np.bincount(prevalence, minlength=N_TIPS + 1)[1:]
    xs = np.arange(1, N_TIPS + 1)
    cols = [ACCENT["highlight"] if k < N_TIPS else ACCENT["duplication"] for k in xs]
    ax.bar(xs, counts, color=cols)
    ax.set_yscale("log")
    ax.set_xlabel("number of species a gene family is present in")
    ax.set_ylabel("gene families (log)")
    ax.set_title("Pan-genome prevalence spectrum", loc="left")
    ax.set_xticks(xs)
    ax.legend(handles=[Patch(color=ACCENT["duplication"], label=f"core ({core})"),
                       Patch(color=ACCENT["highlight"], label=f"accessory ({accessory})")],
              frameon=False, fontsize=8, loc="upper center")
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    fig.savefig(os.path.join(FIGS, "fig3_pangenome.pdf"))
    plt.close(fig)

    # Fig 4 — copy-number heatmap of the most variable gene families
    var = Mg.var(axis=1)
    top = np.argsort(var)[::-1][:60]
    H = Mg[top]
    fig, ax = plt.subplots(figsize=(6.2, 4.2))
    vmax = max(2, int(H.max()))
    im = ax.imshow(H, aspect="auto", cmap="magma_r", vmin=0, vmax=vmax, interpolation="nearest")
    ax.set_xticks(np.arange(len(names)))
    ax.set_xticklabels(names, rotation=45, ha="right", fontsize=7)
    ax.set_yticks([])
    ax.set_ylabel(f"{len(top)} most variable gene families")
    ax.set_title("Gene copy number (gains, losses, duplications)", loc="left")
    cb = fig.colorbar(im, ax=ax, fraction=0.035, pad=0.02, ticks=range(vmax + 1))
    cb.set_label("copies", fontsize=8)
    fig.savefig(os.path.join(FIGS, "fig4_heatmap.pdf"))
    plt.close(fig)

    # Fig 5 — event tally
    fig, ax = plt.subplots(figsize=(5.4, 3.2))
    keys = ["duplication", "transfer", "loss", "pseudogenization", "inversion",
            "transposition", "origination"]
    vals = [ev.get(k, 0) for k in keys]
    ev_col = {"pseudogenization": ACCENT["highlight"], "inversion": MUTED,
              "transposition": MUTED}
    cols = [ev_col.get(k, ACCENT.get(k, MUTED)) for k in keys]
    ax.barh(range(len(keys))[::-1], vals, color=cols)
    ax.set_yticks(range(len(keys))[::-1])
    ax.set_yticklabels([k.capitalize() for k in keys], fontsize=8)
    for i, v in enumerate(vals):
        ax.text(v, len(keys) - 1 - i, f" {v}", va="center", fontsize=8, color=INK)
    ax.set_xlabel("events over the whole tree")
    ax.set_title("Structural events fired", loc="left")
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    fig.savefig(os.path.join(FIGS, "fig5_events.pdf"))
    plt.close(fig)

    # Fig 6 — an example gene tree discordant with the species tree (a duplicated family)
    gene_trees = res.gene_trees()
    # pick a gene family with paralogs (max copy >=2) and present in several species
    cand = [(atom_ids[gene_rows[i]], prevalence[i], Mg[i].max())
            for i in range(len(gene_rows))
            if Mg[i].max() >= 2 and 4 <= prevalence[i] <= N_TIPS]
    example = None
    if cand:
        cand.sort(key=lambda t: (-t[2], -t[1]))
        aid = cand[0][0]
        complete, extant = gene_trees.get(aid, (None, None))
        if extant:
            example = (aid, gid[aid], extant)
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.6))
    draw_tree(axes[0], tree.to_newick(), color=INK, lw=1.8, tip_fontsize=8)
    axes[0].set_title("Species tree", loc="left")
    if example:
        aid, gname, extant = example
        draw_tree(axes[1], extant, color=ACCENT["duplication"], lw=1.6,
                  tip_label=lambda s: s.split("|")[0].split("_")[0], tip_fontsize=7,
                  tip_color=lambda s: INK)
        axes[1].set_title(f"Gene tree — {gname} (with duplication)", loc="left")
    else:
        axes[1].axis("off")
        axes[1].text(0.5, 0.5, "no duplicated family", ha="center", va="center")
    fig.savefig(os.path.join(FIGS, "fig6_genetree.pdf"))
    plt.close(fig)

    # ============================================================ STATS -> LaTeX
    tot_len = tree.total_length if hasattr(tree, "total_length") else None
    mean_genes = np.mean([comp[nm]["distinct_genes"] for nm in names])
    mean_pseud = np.mean([comp[nm]["pseudogenes"] for nm in names])
    mean_size = np.mean([comp[nm]["size"] for nm in names])
    retention = 100 * mean_genes / ancestral_genes
    genic_frac = 100 * sum(b - a for a, b, _ in gff.genes) / gff.length

    macros = {
        "EcNtips": N_TIPS, "EcAge": f"{AGE:g}", "EcTreeSeed": TREE_SEED, "EcSimSeed": SIM_SEED,
        "EcChromLen": f"{gff.length:,}", "EcRawGenes": gff.n_features,
        "EcGenes": len(gff.genes), "EcTrimmed": gff.n_trimmed, "EcDropped": gff.n_dropped,
        "EcGenicFrac": f"{genic_frac:.1f}",
        "EcAtoms": len(res.atoms), "EcGeneAtoms": len(res.gene_atoms()),
        "EcInterAtoms": len(res.intergene_atoms()),
        "EcNfamilies": n_gene_families, "EcCore": core, "EcAccessory": accessory,
        "EcUnique": unique, "EcDuplicated": duplicated,
        "EcMeanGenes": f"{mean_genes:.0f}", "EcMeanPseud": f"{mean_pseud:.0f}",
        "EcRetention": f"{retention:.1f}", "EcMeanSize": f"{mean_size/1e6:.3f}",
        "EcMinGenes": min(comp[nm]["distinct_genes"] for nm in names),
        "EcMaxGenes": max(comp[nm]["distinct_genes"] for nm in names),
        "EcDup": ev.get("duplication", 0), "EcTrans": ev.get("transfer", 0),
        "EcLoss": ev.get("loss", 0), "EcPseudo": ev.get("pseudogenization", 0),
        "EcInv": ev.get("inversion", 0), "EcTransp": ev.get("transposition", 0),
        "EcOrig": ev.get("origination", 0),
        "EcInvR": f"{RATES['inversion']:.0e}", "EcLossR": f"{RATES['loss']:.0e}",
        "EcDupR": f"{RATES['duplication']:.0e}", "EcTransR": f"{RATES['transfer']:.0e}",
        "EcTranspR": f"{RATES['transposition']:.0e}", "EcOrigR": f"{RATES['origination']:g}",
        "EcExt": f"{RATES['extension']:g}", "EcMeanEvLen": f"{1/(1-RATES['extension']):,.0f}",
        "EcPseudoP": f"{RATES['pseudogenization']:g}", "EcReplP": f"{RATES['replacement']:g}",
        "EcExampleGene": example[1] if example else "n/a",
    }
    with open(os.path.join(HERE, "stats.tex"), "w") as f:
        for k, v in macros.items():
            f.write(f"\\newcommand{{\\{k}}}{{{v}}}\n")
    with open(os.path.join(HERE, "stats.json"), "w") as f:
        json.dump({k: str(v) for k, v in macros.items()}, f, indent=2)

    print("=== SUMMARY ===")
    for k, v in macros.items():
        print(f"{k:16s} {v}")


if __name__ == "__main__":
    main()
