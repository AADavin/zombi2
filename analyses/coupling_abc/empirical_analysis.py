"""Empirical analysis: gene-family co-occurrence in a real 43-genome eggNOG dataset (ZOMBI2_DATA),
testing whether genes with a shared function tend to appear together across genomes.

For each COG functional category we ask whether its orthologous groups co-occur more than a
*random same-size* set of groups -- a permutation null that keeps every group's real presence
pattern (so any tree-driven co-occurrence is present in both observed and null) and only randomises
the grouping. We set aside the non-functional buckets from the start ("function unknown" S,
"general prediction" R), so every category shown is a real function.

    python analyses/coupling_abc/empirical_analysis.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))   # use this repo's own zombi2
from eggnog import load_eggnog

HERE = Path(__file__).resolve().parent
FIG = HERE / "figures"; FIG.mkdir(exist_ok=True)
DATA = Path("/Users/aadria/Desktop/CLAUDE/ZOMBI2_DATA")
EGGNOG_FILE = "eggnog_annotations_mini_43 (1).tsv"     # the full 43-genome annotation

INK, MUTED, BLUE, ORANGE, GREY = "#1A1A1A", "#6A6A6A", "#2B6CB0", "#C05621", "#B0B7BE"
plt.rcParams.update({
    "figure.dpi": 120, "savefig.dpi": 200, "font.size": 10.5, "font.family": "serif",
    "axes.edgecolor": INK, "axes.labelcolor": INK, "text.color": INK,
    "xtick.color": INK, "ytick.color": INK, "axes.linewidth": 0.8,
    "axes.spines.top": False, "axes.spines.right": False, "axes.titlelocation": "left",
    "axes.titleweight": "bold", "axes.titlepad": 8, "legend.frameon": False,
})


def save(fig, name):
    fig.savefig(FIG / f"{name}.pdf", bbox_inches="tight")
    fig.savefig(FIG / f"{name}.png", bbox_inches="tight")
    plt.close(fig)


COG_NAME = {
    "J": "Translation", "K": "Transcription", "L": "Replication/repair",
    "D": "Cell cycle", "V": "Defense", "T": "Signal transduction", "M": "Cell wall/membrane",
    "N": "Cell motility", "U": "Trafficking/secretion", "O": "PTM/chaperones",
    "C": "Energy", "G": "Carbohydrate", "E": "Amino acid", "F": "Nucleotide",
    "H": "Coenzyme", "I": "Lipid", "P": "Inorganic ion", "Q": "Secondary metabolites",
}
DROP = {"S", "R", "?"}          # non-functional buckets, set aside from the start

# --- load the full 43-genome profile ----------------------------------------------
pm, category = load_eggnog(DATA / EGGNOG_FILE, genomes_tsv=DATA / "02_miniset2_genomes.tsv")
G = len(pm.species)
P = (pm.matrix > 0).astype(float)
cat = np.array([category[og] for og in pm.families])
freq = np.bincount(P.sum(1).astype(int), minlength=G + 1)[1:G + 1]     # pangenome (all groups)
print(f"empirical: {P.shape[0]} orthologous groups x {G} genomes")

# functional groups that vary across genomes (the ones that carry a co-occurrence signal)
keep = (P.std(1) > 0) & ~np.isin(cat, list(DROP))
Pk, ck = P[keep], cat[keep]
C = np.corrcoef(Pk)
rng = np.random.default_rng(0)
print(f"functional variable groups used: {keep.sum()}")

# --- per-category cohesion vs random same-size groups -----------------------------
percat = {}
for c in sorted(set(ck)):
    idx = np.where(ck == c)[0]; n = len(idx)
    if n < 15:
        continue
    internal = float(C[np.ix_(idx, idx)][np.triu_indices(n, 1)].mean())
    null = np.empty(300)
    for i in range(300):
        r = rng.choice(len(ck), n, replace=False)
        null[i] = C[np.ix_(r, r)][np.triu_indices(n, 1)].mean()
    percat[c] = dict(n=int(n), prev=float(Pk[idx].mean()), internal=internal,
                     excess=float(internal - null.mean()),
                     p=float((np.sum(null >= internal) + 1) / (len(null) + 1)))
n_sig = sum(1 for v in percat.values() if v["p"] < 0.05)
print(f"per-category: {n_sig}/{len(percat)} functional categories significantly cohesive (p<0.05)")

# --- same-function vs different-function co-occurrence, permutation null -----------
iu = np.triu_indices(len(ck), 1)
corr = C[iu]; same = ck[iu[0]] == ck[iu[1]]
obs = float(corr[same].mean() - corr[~same].mean())
null = np.empty(500)
for i in range(500):
    cp = rng.permutation(ck); sm = cp[iu[0]] == cp[iu[1]]
    null[i] = corr[sm].mean() - corr[~sm].mean()
z = float((obs - null.mean()) / null.std())
p = float((np.sum(null >= obs) + 1) / (len(null) + 1))
print(f"same- vs different-function co-occurrence gap = {obs:+.4f}  (z={z:.1f}, p={p:.4f})")

results = dict(n_ogs=int(P.shape[0]), n_genomes=G, n_functional_variable=int(keep.sum()),
               gap=obs, z=z, p=p, n_sig=n_sig, n_categories=len(percat),
               frequency_spectrum=freq.tolist(), per_category=percat)
(HERE / "empirical_results.json").write_text(json.dumps(results, indent=2, default=float))


# ==================================================================================
# Figure 5 — how gene presence is distributed, and the same-function co-occurrence test
# ==================================================================================
def fig_empirical():
    fig, (a, b) = plt.subplots(1, 2, figsize=(11.4, 4.5))
    fig.subplots_adjust(wspace=0.30)
    a.bar(range(1, G + 1), freq, color=MUTED, width=0.85)
    a.set_xlabel("number of genomes a gene group is present in")
    a.set_ylabel("gene groups")
    a.set_title(f"a  How gene presence is distributed ({P.shape[0]} groups, {G} genomes)")

    b.hist(null, bins=30, color="#C7D0DA", label="if function didn't matter (shuffled)")
    b.axvline(obs, color=ORANGE, lw=2.2, label=f"observed  (p = {p:.3f})")
    b.axvline(0, color=MUTED, lw=0.8, ls=(0, (3, 3)))
    b.set_xlabel("extra co-occurrence among same-function genes")
    b.set_ylabel("shuffled replicates")
    b.set_title("b  Same-function genes appear together more than by chance")
    b.legend(loc="upper left", fontsize=9)
    save(fig, "fig5_empirical")


# ==================================================================================
# Figure 6 — which functions form co-occurrence modules
# ==================================================================================
def fig_cohesion():
    from matplotlib.patches import Patch
    items = sorted(percat.items(), key=lambda kv: kv[1]["excess"], reverse=True)
    labels = [f"{c}  {COG_NAME.get(c, c)}" for c, _ in items]
    vals = [v["excess"] for _, v in items]
    ns = [v["n"] for _, v in items]
    cols = [BLUE if v["p"] < 0.05 else GREY for _, v in items]
    sig = ["*" if v["p"] < 0.05 else "" for _, v in items]
    fig, ax = plt.subplots(figsize=(7.7, 0.44 * len(items) + 1.2))
    y = range(len(items))
    ax.barh(y, vals, color=cols, height=0.7)
    ax.set_yticks(y); ax.set_yticklabels(labels, fontsize=9.5); ax.invert_yaxis()
    ax.axvline(0, color=MUTED, lw=0.8)
    ax.set_xlabel("how much more genes of this function appear together\n(vs. a random group of genes the same size)")
    ax.set_title("Which functions travel together across genomes")
    for i, (v, n, s) in enumerate(zip(vals, ns, sig)):
        ax.text(v + 0.004 if v >= 0 else v - 0.004, i, f"n={n}{(' ' + s) if s else ''}",
                va="center", ha="left" if v >= 0 else "right", fontsize=7.8, color=MUTED)
    ax.legend(handles=[Patch(color=BLUE, label="more than chance (p<0.05)"),
                       Patch(color=GREY, label="not significant")],
              loc="lower right", fontsize=9)
    ax.margins(x=0.15)
    save(fig, "fig6_cohesion")


fig_empirical()
fig_cohesion()
print("wrote fig5_empirical, fig6_cohesion, empirical_results.json")
print("top functional modules:",
      [(c, round(v["excess"], 3)) for c, v in sorted(percat.items(), key=lambda kv: -kv[1]["excess"])][:6])
