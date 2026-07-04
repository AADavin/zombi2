"""Empirical analysis: gene-family co-occurrence in a real eggNOG dataset (ZOMBI2_DATA),
and a test of the COG-category-as-coupling idea.

Tests whether orthologous groups sharing a COG functional category co-occur more than groups
in different categories -- against a label-permutation null that preserves each group's actual
presence pattern (so the phylogenetic confound is present in both observed and null, and only the
category assignment is randomized). A positive result supports using COG categories to structure
the coupling matrix J of the Potts model.

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

INK, MUTED, BLUE, ORANGE = "#1A1A1A", "#6A6A6A", "#2B6CB0", "#C05621"
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


COG_NAME = {  # standard COG functional categories (abbreviated)
    "J": "Translation", "K": "Transcription", "L": "Replication/repair",
    "D": "Cell cycle", "V": "Defense", "T": "Signal transduction", "M": "Cell wall/membrane",
    "N": "Cell motility", "U": "Trafficking/secretion", "O": "PTM/chaperones",
    "C": "Energy", "G": "Carbohydrate", "E": "Amino acid", "F": "Nucleotide",
    "H": "Coenzyme", "I": "Lipid", "P": "Inorganic ion", "Q": "Secondary metabolites",
    "R": "General prediction", "S": "Function unknown",
}

# --- load the empirical profile ---------------------------------------------------
pm, category = load_eggnog(DATA / "eggnog_annotations_mini_43.tsv",
                           genomes_tsv=DATA / "02_miniset2_genomes.tsv")
G = len(pm.species)
P = (pm.matrix > 0).astype(float)
cat = np.array([category[og] for og in pm.families])
print(f"empirical: {P.shape[0]} orthologous groups x {G} genomes")

# shell OGs (variable presence) with a known category carry the co-occurrence signal
keep = (P.std(1) > 0) & (cat != "?")
Pk, ck = P[keep], cat[keep]
C = np.corrcoef(Pk)
iu = np.triu_indices(len(ck), k=1)
corr = C[iu]
same = ck[iu[0]] == ck[iu[1]]
obs = corr[same].mean() - corr[~same].mean()

# label-permutation null (phylogeny-controlled: presence patterns fixed, categories shuffled)
rng = np.random.default_rng(0)
null = np.empty(500)
for i in range(500):
    cp = rng.permutation(ck)
    sm = cp[iu[0]] == cp[iu[1]]
    null[i] = corr[sm].mean() - corr[~sm].mean()
z = float((obs - null.mean()) / null.std())
p = float((np.sum(null >= obs) + 1) / (len(null) + 1))
print(f"within-vs-between category co-occurrence gap = {obs:+.4f}  (z={z:.2f}, p={p:.4f})")

# per-category internal cohesion (mean within-category correlation minus global between baseline)
between = corr[~same].mean()
cohesion = {}
for c in sorted(set(ck)):
    idx = np.where(ck == c)[0]
    if len(idx) < 20:
        continue
    sub = C[np.ix_(idx, idx)]
    m = sub[np.triu_indices(len(idx), k=1)].mean()
    cohesion[c] = (m - between, len(idx))

freq = np.bincount(P.sum(1).astype(int), minlength=G + 1)[1:G + 1]     # pangenome spectrum

results = dict(n_ogs=int(P.shape[0]), n_genomes=G, n_shell_used=int(keep.sum()),
               gap=obs, z=z, p=p, null_mean=float(null.mean()), null_sd=float(null.std()),
               frequency_spectrum=freq.tolist(),
               cohesion={c: dict(excess=v[0], n=v[1]) for c, v in cohesion.items()})
(HERE / "empirical_results.json").write_text(json.dumps(results, indent=2, default=float))


# ==================================================================================
# Figure 5 — pangenome spectrum + the COG-category co-occurrence test
# ==================================================================================
def fig_empirical():
    fig, (a, b) = plt.subplots(1, 2, figsize=(11.2, 4.5))
    fig.subplots_adjust(wspace=0.32)
    a.bar(range(1, G + 1), freq, color=MUTED, width=0.7)
    a.set_xlabel("number of genomes an OG is present in")
    a.set_ylabel("orthologous groups")
    a.set_title(f"a  Pangenome spectrum ({P.shape[0]} OGs, {G} genomes)")
    a.set_xticks(range(1, G + 1))

    b.hist(null, bins=30, color="#B8C4D0", label="permutation null")
    b.axvline(obs, color=ORANGE, lw=2.2, label=f"observed  (p = {p:.3f})")
    b.axvline(0, color=MUTED, lw=0.8, ls=(0, (3, 3)))
    b.set_xlabel("within $-$ between category co-occurrence")
    b.set_ylabel("null replicates")
    b.set_title("b  Same-category OGs co-occur more (vs shuffled labels)")
    b.legend(loc="upper left", fontsize=9)
    save(fig, "fig5_empirical")


# ==================================================================================
# Figure 6 — which COG categories are most internally cohesive (informing J structure)
# ==================================================================================
def fig_cohesion():
    items = sorted(cohesion.items(), key=lambda kv: kv[1][0], reverse=True)
    labels = [f"{c}  {COG_NAME.get(c, c)}" for c, _ in items]
    vals = [v[0] for _, v in items]
    ns = [v[1] for _, v in items]
    fig, ax = plt.subplots(figsize=(7.4, 0.42 * len(items) + 1.2))
    y = range(len(items))
    ax.barh(y, vals, color=[BLUE if v > 0 else MUTED for v in vals], height=0.68)
    ax.set_yticks(y); ax.set_yticklabels(labels, fontsize=9)
    ax.invert_yaxis()
    ax.axvline(0, color=MUTED, lw=0.8)
    ax.set_xlabel("internal cohesion  (within-category $-$ global-between co-occurrence)")
    ax.set_title("Co-occurrence cohesion by COG functional category")
    for i, (v, n) in enumerate(zip(vals, ns)):
        ax.text(v + (0.0005 if v >= 0 else -0.0005), i, f"n={n}",
                va="center", ha="left" if v >= 0 else "right", fontsize=7.5, color=MUTED)
    save(fig, "fig6_cohesion")


fig_empirical()
fig_cohesion()
print("wrote fig5_empirical, fig6_cohesion, empirical_results.json")
print("top cohesive categories:",
      [(c, round(v[0], 4)) for c, v in sorted(cohesion.items(), key=lambda kv: -kv[1][0])[:6]])
