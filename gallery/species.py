"""Species-tree examples: forward birth–death trees, with the diversification model made visible."""

from __future__ import annotations

import helpers as h
from helpers import Example

import phylustrator as ph
from zombi2.species import simulate_species_tree
from zombi2.rates import modifiers as mod


def yule(out):
    sp = simulate_species_tree(birth=1.0, n_extant=100, seed=7)   # pure birth, no extinction
    style = ph.Style(width=h.TREE_W, height=820, margin=88, branch_width=2.4)
    (ph.trees.plot(ph.trees.loads(sp.complete_tree.to_newick()), style=style)
     + ph.trees.time_axis("time", tick_size=22, label_size=28)).save(out)


def extinct_lineages(out):
    sp = simulate_species_tree(birth=1.0, death=0.6, n_extant=50, seed=3)
    ct = sp.complete_tree
    tree = ph.trees.loads(ct.to_newick())
    dashed = h.dashed_extinct(tree, {f"n{n.id}" for n in ct.extinct()})
    style = ph.Style(width=1400, height=760, margin=88, branch_width=2.0)   # ~40% wider than tall
    (ph.trees.plot(tree, dashed=dashed, style=style)
     + ph.trees.time_axis("time", tick_size=22, label_size=28)).save(out)


def mass_extinction(out):
    sp = simulate_species_tree(birth=1.5, death=0.15, total_time=3.8,
                               mass_extinctions=[(3.0, 0.75)], seed=2)   # 85 extant + 56 extinct
    ct = sp.complete_tree
    tree = ph.trees.loads(ct.to_newick())
    dashed = h.dashed_extinct(tree, {f"n{n.id}" for n in ct.extinct()})
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


# --- copy-paste-reproducible snippets shown on the detail view --------------
_DASH = '''
tree = ph.trees.loads(ct.to_newick())
extinct = {f"n{n.id}" for n in ct.extinct()}
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


EXAMPLES = [
    Example("basic", "Yule tree", "Pure birth, no extinction — a forward tree of 100 lineages.",
            "pure birth (Yule)", yule, code=_C_BASIC),
    Example("extinct", "Extinct lineages",
            "The full history behind 100 survivors — their branches solid, extinct lineages dashed.",
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
]
