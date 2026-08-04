"""Species-tree examples: forward birth–death trees, with the diversification model made visible.

The last example is not a single tree: it runs the level thousands of times and measures the trees,
which is the other way to use this level.
"""

from __future__ import annotations

import helpers as h
from helpers import Example

import phylustrator as ph
from zombi2.species import simulate_species_tree
from zombi2.rates import modifiers as mod
from zombi2.tree import gamma_statistic


def yule(out):
    sp = simulate_species_tree(birth=1.0, n_extant=100, seed=7)   # pure birth, no extinction
    style = ph.Style(width=h.TREE_W, height=820, margin=88, branch_width=2.4)
    (ph.trees.plot(ph.trees.loads(sp.complete_tree.to_newick()), style=style)
     + ph.trees.time_axis("time", tick_size=22, label_size=28)).save(out)


def extinct_lineages(out):
    sp = simulate_species_tree(birth=1.0, death=0.6, n_extant=50, seed=3)
    ct = sp.complete_tree
    tree = ph.trees.loads(ct.to_newick())
    dashed = h.dashed_extinct(tree, ct)
    style = ph.Style(width=1400, height=760, margin=88, branch_width=2.0)   # ~40% wider than tall
    (ph.trees.plot(tree, dashed=dashed, style=style)
     + ph.trees.time_axis("time", tick_size=22, label_size=28)).save(out)


def mass_extinction(out):
    sp = simulate_species_tree(birth=1.5, death=0.15, total_time=3.8,
                               mass_extinctions=[(3.0, 0.75)], seed=2)   # 85 extant + 56 extinct
    ct = sp.complete_tree
    tree = ph.trees.loads(ct.to_newick())
    dashed = h.dashed_extinct(tree, ct)
    present = max(n.end_time for n in ct.nodes.values())
    tree_png = out.replace(".png", "_tree.png")
    TW, TM = 1320, 70
    style = ph.Style(width=TW, height=560, margin=TM, branch_width=1.3)
    (ph.trees.plot(tree, dashed=dashed, style=style)
     + ph.trees.time_marker(3.0, label="mass extinction", color="#111111", label_size=20)).save(tree_png)
    times, counts = h.lineages_through_time(ct)

    def panel(ax):
        ax.step(times, counts, where="post", color="#333333", lw=1.8)
        ax.fill_between(times, counts, step="post", color="#999999", alpha=0.12)
        ax.axvline(3.0, color="#111111", lw=1.4, ls="--")     # the extinction pulse

    h.composite_below(tree_png, present, out, panel, "lineages",
                      tree_w=TW, margin=TM, figsize=(14, 9), axis_fontsize=17)


def rate_shift(out):
    rate = 1.0 * mod.OnTime({0: 0.8, 2.0: 1.9, 3.5: 0.8})   # slow, then fast, then slow
    sp = simulate_species_tree(birth=rate, total_time=4.5, seed=3)     # 155 extant
    tree = ph.trees.loads(sp.complete_tree.to_newick())
    style = ph.Style(width=1400, height=900, margin=88, branch_width=1.5)
    (ph.trees.plot(tree, style=style)
     + ph.trees.time_marker(2.0, label="faster", label_size=18)        # slow -> fast
     + ph.trees.time_marker(3.5, label="slower", label_size=18)        # fast -> slow
     + ph.trees.time_axis("time", tick_size=22, label_size=28)).save(out)


def diversity_dependent(out):
    sp = simulate_species_tree(birth=1.0 * mod.OnTotalDiversity(cap=100), total_time=10.0, seed=6)
    ct = sp.complete_tree
    present = max(n.end_time for n in ct.nodes.values())
    tree_png = out.replace(".png", "_tree.png")
    TW, TM = 1100, 55
    ph.trees.plot(ph.trees.loads(ct.to_newick()),
            style=ph.Style(width=TW, height=640, margin=TM, branch_width=1.6)).save(tree_png)
    times, counts = h.lineages_through_time(ct)

    def panel(ax):
        ax.step(times, counts, where="post", color="#333333", lw=1.8)
        ax.fill_between(times, counts, step="post", color="#999999", alpha=0.12)
        ax.axhline(100, color="#111111", lw=1.2, ls="--")     # the cap

    h.composite_below(tree_png, present, out, panel, "lineages",
                      tree_w=TW, margin=TM, figsize=(12, 9), axis_fontsize=16)


# --- many trees, measured (local to this module: helpers.py is shared) ------

_CONST, _SLOW = "#3A7CA5", "#C25A3C"      # constant rate · diversity-dependent
_REPS, _TIPS, _CAP, _SHOWN = 2000, 100, 110, 60


def shape_statistics(out):
    """Simulate _REPS trees under each of two processes and compare their γ distributions."""
    import numpy as np
    import matplotlib.pyplot as plt          # helpers.py already selected the Agg backend

    processes = (("constant rate", 1.0, _CONST),
                 ("diversity-dependent", 1.0 * mod.OnTotalDiversity(cap=_CAP), _SLOW))
    gammas, curves = {}, {}
    for name, birth, _ in processes:
        gs, cs = [], []
        for s in range(1, _REPS + 1):
            sp = simulate_species_tree(birth=birth, n_extant=_TIPS, seed=s)
            gs.append(gamma_statistic(sp.extant_tree))
            if s <= _SHOWN:                  # a subsample, drawn as lineage curves
                cs.append(h.lineages_through_time(sp.complete_tree))
        gammas[name], curves[name] = gs, cs

    fig, (axl, axr) = plt.subplots(1, 2, figsize=(13.2, 5.2),
                                   gridspec_kw={"width_ratios": [1.0, 1.35], "wspace": 0.20})
    for name, _, col in processes:
        for times, counts in curves[name]:
            # lined up at the present, the way a slowdown reads: constant rate climbs straight to
            # the last moment, diversity-dependence flattens as the cap fills
            axl.step([t - times[-1] for t in times], counts, where="post",
                     color=col, lw=0.9, alpha=0.25)
    axl.set_yscale("log")
    axl.set_yticks([1, 10, 100])
    axl.set_yticklabels(["1", "10", "100"])
    axl.set_xlabel("time before the present", fontsize=13)
    axl.set_ylabel("lineages", fontsize=13)
    axl.set_title(f"{_SHOWN} of the trees, lined up at the present", fontsize=12.5, loc="left",
                  color="#555555")
    handles = [plt.Line2D([], [], color=col, lw=2.6) for _, _, col in processes]
    axl.legend(handles, ["constant rate", f"diversity-dependent (cap {_CAP})"], loc="upper left",
               frameon=False, fontsize=12.5, handlelength=1.4, borderpad=0.1)

    edges = np.arange(-11.0, 4.01, 0.25)
    for name, _, col in processes:
        axr.hist(gammas[name], bins=edges, color=col, alpha=0.62, edgecolor="white", lw=0.3)
    top = max(np.histogram(g, bins=edges)[0].max() for g in gammas.values())
    axr.set_ylim(0, top * 1.30)
    for name, _, col in processes:
        m, sd = float(np.mean(gammas[name])), float(np.std(gammas[name]))
        axr.text(m, top * 1.06, f"{name}\nγ = {m:+.2f} ± {sd:.2f}", ha="center", va="bottom",
                 color=col, fontsize=12.5, linespacing=1.4)
    axr.set_xlabel("γ  (Pybus & Harvey)", fontsize=13)
    axr.set_ylabel("trees", fontsize=13)
    axr.set_title(f"{_REPS} trees per process, {_TIPS} tips each", fontsize=12.5, loc="left",
                  color="#555555")
    for ax in (axl, axr):
        for spine in ("top", "right"):
            ax.spines[spine].set_visible(False)
        ax.tick_params(labelsize=11.5)
    fig.savefig(out, dpi=125, bbox_inches="tight")
    plt.close(fig)


# --- copy-paste-reproducible snippets shown on the detail view --------------
_DASH = '''
tree = ph.trees.loads(ct.to_newick())
lab = ct.labels()                           # {id: 'n<id>'} — or 'e<id>' where the lineage went extinct
extinct = {lab[n.id] for n in ct.extinct_leaves()}
dashed = set()                              # dash a branch whose whole subtree is extinct
for node in tree.walk("postorder"):
    if node.is_leaf and node.name in extinct:
        dashed.add(node.name)
    elif node.children and all(c.name in dashed for c in node.children):
        dashed.add(node.name)'''

_C_BASIC = '''\
### simulate  —  ZOMBI2 (Python API)
from zombi2.species import simulate_species_tree

# pure birth (Yule): birth only, death defaults to 0 — no lineage goes extinct
sp = simulate_species_tree(birth=1.0, n_extant=100, seed=7)

### plot  —  Phylustrator
import phylustrator as ph

tree = ph.trees.loads(sp.complete_tree.to_newick())
(ph.trees.plot(tree, style=ph.Style(branch_width=2.4))
 + ph.trees.time_axis("time", tick_size=22, label_size=28)).save("yule.png")'''

_C_EXTINCT = '''\
### simulate
from zombi2.species import simulate_species_tree

sp = simulate_species_tree(birth=1.0, death=0.6, n_extant=50, seed=3)
ct = sp.complete_tree

### plot  —  extinct lineages dashed
import phylustrator as ph
''' + _DASH + '''
(ph.trees.plot(tree, dashed=dashed, style=ph.Style(width=1400, height=760, branch_width=2.0))
 + ph.trees.time_axis("time", tick_size=22, label_size=28)).save("extinct.png")'''

_C_MASSEXT = '''\
### simulate  —  a pulse at t = 3 culls 75% of lineages (time, fraction_lost)
from zombi2.species import simulate_species_tree

sp = simulate_species_tree(birth=1.5, death=0.15, total_time=3.8,
                           mass_extinctions=[(3.0, 0.75)], seed=2)
ct = sp.complete_tree

### plot  —  tree (pulse marked) + a diversity skyline, on a shared time axis
import phylustrator as ph
import helpers as h                         # gallery helper: composites tree + panel
''' + _DASH + '''
(ph.trees.plot(tree, dashed=dashed, style=ph.Style(branch_width=1.3))
 + ph.trees.time_marker(3.0, label="mass extinction", color="#111", label_size=20)).save("tree.png")

times, counts = h.lineages_through_time(ct)  # standing diversity over time
def panel(ax):
    ax.step(times, counts, where="post", color="#333", lw=1.8)
    ax.fill_between(times, counts, step="post", color="#999", alpha=0.12)
    ax.axvline(3.0, color="#111", ls="--")         # the pulse

present = max(n.end_time for n in ct.nodes.values())
h.composite_below("tree.png", present, "massext.png", panel, "lineages")'''

_C_RATESHIFT = '''\
### simulate  —  speciation slow (0.8), then fast (1.9), then slow again, by time
from zombi2.species import simulate_species_tree
from zombi2.rates import modifiers as mod

rate = 1.0 * mod.OnTime({0: 0.8, 2.0: 1.9, 3.5: 0.8})
sp = simulate_species_tree(birth=rate, total_time=4.5, seed=3)

### plot  —  dashed lines mark the regime changes
import phylustrator as ph

tree = ph.trees.loads(sp.complete_tree.to_newick())
(ph.trees.plot(tree, style=ph.Style(branch_width=1.5))
 + ph.trees.time_marker(2.0, label="faster", label_size=18)          # slow -> fast
 + ph.trees.time_marker(3.5, label="slower", label_size=18)          # fast -> slow
 + ph.trees.time_axis("time", tick_size=22, label_size=28)).save("rateshift.png")'''

_C_DIVERSITY = '''\
### simulate  —  speciation slows as diversity approaches a cap of 100
from zombi2.species import simulate_species_tree
from zombi2.rates import modifiers as mod

sp = simulate_species_tree(birth=1.0 * mod.OnTotalDiversity(cap=100), total_time=10.0, seed=6)
ct = sp.complete_tree

### plot  —  tree + a diversity skyline, on a shared time axis
import phylustrator as ph
import helpers as h

ph.trees.plot(ph.trees.loads(ct.to_newick()), style=ph.Style(branch_width=1.6)).save("tree.png")
times, counts = h.lineages_through_time(ct)

def panel(ax):
    ax.step(times, counts, where="post", color="#333", lw=1.8)
    ax.fill_between(times, counts, step="post", color="#999", alpha=0.12)
    ax.axhline(100, color="#111", ls="--")       # the cap

present = max(n.end_time for n in ct.nodes.values())
h.composite_below("tree.png", present, "diversity.png", panel, "lineages")'''

_C_SHAPE = '''\
### measure  —  gamma is a zombi2 tool, so one tree needs no code at all
zombi2 tools tree out/species/species_extant.nwk --gamma        # gamma  -6.31

### simulate  —  2000 trees of 100 tips under each process (a few seconds)
from zombi2.species import simulate_species_tree
from zombi2.rates import modifiers as mod
from zombi2.tree import gamma_statistic     # 0 on average under constant rates

gammas = {}
for name, birth in (("constant rate", 1.0),
                    ("diversity-dependent", 1.0 * mod.OnTotalDiversity(cap=110))):
    gammas[name] = [gamma_statistic(simulate_species_tree(birth=birth, n_extant=100,
                                                          seed=s).extant_tree)
                    for s in range(1, 2001)]
# constant rate         gamma = +0.05 +- 0.99      (the constant-rate null: mean 0, sd 1)
# diversity-dependent   gamma = -6.27 +- 1.08
# Colless imbalance does not separate them at all: diversity-dependence rescales the waiting
# times and leaves the topology alone, so under the same seed it returns the identical value.

### plot  —  the two distributions
import matplotlib.pyplot as plt
import numpy as np

edges = np.arange(-11.0, 4.01, 0.25)
for name, colour in (("constant rate", "#3A7CA5"), ("diversity-dependent", "#C25A3C")):
    plt.hist(gammas[name], bins=edges, color=colour, alpha=0.62, label=name)
plt.xlabel("gamma")
plt.ylabel("trees")
plt.legend()
plt.savefig("shape.png", dpi=125, bbox_inches="tight")'''


EXAMPLES = [
    Example("basic", "Yule tree", "Pure birth, no extinction — a forward tree of 100 lineages.",
            "pure birth (Yule)", yule, code=_C_BASIC),
    Example("extinct", "Extinct lineages",
            "The full history behind 50 survivors — their branches solid, extinct lineages dashed.",
            "birth–death", extinct_lineages, code=_C_EXTINCT),
    Example("massext", "Mass extinction",
            "A pulse at t&nbsp;=&nbsp;3 culls 75% of lineages — the skyline drops sharply at the dashed "
            "line, then recovers.",
            "mass extinction · +&nbsp;skyline", mass_extinction, code=_C_MASSEXT),
    Example("rateshift", "Rate shifts",
            "Speciation runs slow, then fast, then slow — the burst packs branches between the two "
            "dashed regime lines.",
            "time-varying birth", rate_shift, code=_C_RATESHIFT),
    Example("diversity", "Diversity-dependent",
            "Speciation slows as diversity fills up; the skyline rises and plateaus at the cap of 100.",
            "birth–death · +&nbsp;skyline", diversity_dependent, code=_C_DIVERSITY),
    Example("shape", "Shape statistics over many trees",
            "Two thousand trees under each of two processes. "
            "<code>zombi2 tools tree --gamma</code> separates them almost perfectly.",
            "simulation study · 4000&nbsp;trees", shape_statistics, code=_C_SHAPE),
]
