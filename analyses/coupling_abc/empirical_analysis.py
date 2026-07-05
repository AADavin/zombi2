"""Empirical analysis: gene-family co-occurrence in a real eggNOG dataset (ZOMBI2_DATA),
and a test of the COG-category-as-coupling idea, on the full 43-genome annotation.

Two tests, both against a permutation null that keeps every group's real presence pattern (so
the phylogenetic confound is present in observed and null) and only randomises the category:

1. Per-category cohesion -- is each COG category more internally co-occurring than a *random
   same-size* group of orthologous groups? This is the clean test of "is this category a module".
2. A functional-category aggregate (within- vs between-category), excluding the catch-all
   "Function unknown" (S) and "General prediction" (R): do functional groups co-occur by category?

A positive result supports using COG *functional* categories to scaffold the coupling matrix J.

    python analyses/coupling_abc/empirical_analysis.py
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

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
    "R": "General prediction", "S": "Function unknown",
}
NONFUNC = {"S", "R"}          # catch-alls: not functional modules, excluded from the scaffold

# --- load the full 43-genome profile ----------------------------------------------
pm, category = load_eggnog(DATA / EGGNOG_FILE, genomes_tsv=DATA / "02_miniset2_genomes.tsv")
G = len(pm.species)
P = (pm.matrix > 0).astype(float)
cat = np.array([category[og] for og in pm.families])
print(f"empirical: {P.shape[0]} orthologous groups x {G} genomes")

# shell OGs (variable presence) with a known category carry the co-occurrence signal
keep = (P.std(1) > 0) & (cat != "?")
Pk, ck = P[keep], cat[keep]
C = np.corrcoef(Pk)
rng = np.random.default_rng(0)

# --- (1) per-category cohesion vs random same-size groups -------------------------
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
n_sig = sum(1 for c, v in percat.items() if v["p"] < 0.05 and c not in NONFUNC)
n_func = sum(1 for c in percat if c not in NONFUNC)
print(f"per-category: {n_sig}/{n_func} functional categories significantly cohesive (p<0.05)")

# --- (2) functional-only aggregate (exclude S/R), label-permutation null ----------
fi = np.where(~np.isin(ck, list(NONFUNC)))[0]
Cf, cf = C[np.ix_(fi, fi)], ck[fi]
iuf = np.triu_indices(len(fi), 1)
corrf = Cf[iuf]; samef = cf[iuf[0]] == cf[iuf[1]]
obs = float(corrf[samef].mean() - corrf[~samef].mean())
null = np.empty(500)
for i in range(500):
    cp = rng.permutation(cf); sm = cp[iuf[0]] == cp[iuf[1]]
    null[i] = corrf[sm].mean() - corrf[~sm].mean()
z = float((obs - null.mean()) / null.std())
p = float((np.sum(null >= obs) + 1) / (len(null) + 1))
print(f"functional-only within-vs-between gap = {obs:+.4f}  (z={z:.1f}, p={p:.4f}); "
      f"{len(fi)} functional OGs")

freq = np.bincount(P.sum(1).astype(int), minlength=G + 1)[1:G + 1]     # pangenome spectrum
results = dict(n_ogs=int(P.shape[0]), n_genomes=G, n_shell_used=int(keep.sum()),
               functional_gap=obs, functional_z=z, functional_p=p,
               n_functional_ogs=int(len(fi)), n_sig=n_sig, n_functional=n_func,
               frequency_spectrum=freq.tolist(), per_category=percat)
(HERE / "empirical_results.json").write_text(json.dumps(results, indent=2, default=float))


# ==================================================================================
# Figure 5 — pangenome spectrum + the functional-category co-occurrence test
# ==================================================================================
def fig_empirical():
    fig, (a, b) = plt.subplots(1, 2, figsize=(11.4, 4.5))
    fig.subplots_adjust(wspace=0.30)
    a.bar(range(1, G + 1), freq, color=MUTED, width=0.85)
    a.set_xlabel("number of genomes an OG is present in")
    a.set_ylabel("orthologous groups")
    a.set_title(f"a  Pangenome spectrum ({P.shape[0]} OGs, {G} genomes)")

    b.hist(null, bins=30, color="#C7D0DA", label="permutation null")
    b.axvline(obs, color=ORANGE, lw=2.2, label=f"observed  (p = {p:.3f})")
    b.axvline(0, color=MUTED, lw=0.8, ls=(0, (3, 3)))
    b.set_xlabel("within $-$ between category co-occurrence")
    b.set_ylabel("null replicates")
    b.set_title("b  Functional-category OGs co-occur more (S, R excluded)")
    b.legend(loc="upper left", fontsize=9)
    save(fig, "fig5_empirical")


# ==================================================================================
# Figure 6 — per-category cohesion + significance (the coupling scaffold)
# ==================================================================================
def fig_cohesion():
    from matplotlib.patches import Patch
    items = sorted(percat.items(), key=lambda kv: kv[1]["excess"], reverse=True)
    labels, vals, cols, ns, sig = [], [], [], [], []
    for c, v in items:
        labels.append(f"{c}  {COG_NAME.get(c, c)}")
        vals.append(v["excess"]); ns.append(v["n"])
        catchall = c in NONFUNC
        signif = v["p"] < 0.05 and not catchall
        cols.append(BLUE if signif else (ORANGE if catchall else GREY))
        sig.append("*" if v["p"] < 0.05 else "")
    fig, ax = plt.subplots(figsize=(7.8, 0.42 * len(items) + 1.3))
    y = range(len(items))
    ax.barh(y, vals, color=cols, height=0.7)
    ax.set_yticks(y); ax.set_yticklabels(labels, fontsize=9); ax.invert_yaxis()
    ax.axvline(0, color=MUTED, lw=0.8)
    ax.set_xlabel("internal cohesion  (within-category $-$ random same-size groups)")
    ax.set_title("Which COG categories are co-occurrence modules")
    for i, (v, n, s) in enumerate(zip(vals, ns, sig)):
        ax.text(v + 0.004 if v >= 0 else v - 0.004, i, f"n={n}{(' ' + s) if s else ''}",
                va="center", ha="left" if v >= 0 else "right", fontsize=7.5, color=MUTED)
    ax.legend(handles=[Patch(color=BLUE, label="functional module (p<0.05)"),
                       Patch(color=GREY, label="not significant"),
                       Patch(color=ORANGE, label="catch-all (excluded)")],
              loc="lower right", fontsize=8.5)
    ax.margins(x=0.16)
    save(fig, "fig6_cohesion")


fig_empirical()
fig_cohesion()
print("wrote fig5_empirical, fig6_cohesion, empirical_results.json")
print("top functional modules:",
      [(c, round(v["excess"], 3)) for c, v in sorted(percat.items(), key=lambda kv: -kv[1]["excess"])
       if c not in NONFUNC][:6])
