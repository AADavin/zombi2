"""Trait-evolution examples: a trait evolving down the tree, branches coloured by its value."""

from __future__ import annotations

import helpers as h
from helpers import Example

import phylustrator as ph
from zombi2.params import Drift, LogNormal, PerLineage, TotalDiversity
from zombi2.species import simulate_species_tree
from zombi2.traits import simulate_continuous, simulate_discrete

SEED = 42
N = 100

#: the carrying capacity the diversity-dependent card's tree and trait rate share
_DD_CAP = 80

# two binary characters X, Y as compound Mk states; only single-bit flips are named (all else 0).
# the dependence: gaining Y is slow while X is absent (00->01) but fast once X is present (10->11).
_DEP_STATES = ("00", "01", "10", "11")
_DEP_RATES = {"00->01": 0.05, "01->00": 2.0, "10->11": 2.0, "11->10": 0.05,   # flip Y (rate ~ X)
              "00->10": 0.5, "10->00": 0.5, "01->11": 0.5, "11->01": 0.5}      # flip X
# one binary palette per trait (shared "absent"); X and Y are shown as two lanes on each branch.
_ABSENT = "#c2cac8"
_XPAL = {"0": _ABSENT, "1": "#3C8D6E"}      # X present = green
_YPAL = {"0": _ABSENT, "1": "#8B6B9E"}      # Y present = muted purple
# the driven example: a discrete habitat (the driver) and the body size it speeds up (the target)
_HABITAT = {"stable": "#3A7CA5", "fluctuating": "#E4572E"}
_SWITCH, _BASE, _FACTOR = 0.4, 0.25, 20.0


def _tree():
    sp = simulate_species_tree(birth=1.0, n_extant=N, seed=SEED)
    return sp.complete_tree, ph.trees.loads(sp.complete_tree.to_newick())


def _style():
    # headroom: these are 100-tip trees whose first branch reaches the top-left corner, where the
    # colorbar and the legend are drawn. Without a row of its own the key shares one with a branch.
    return ph.Style(width=1200, height=1070, margin=80, headroom=70, branch_width=1.4)


def brownian_motion(out):
    ct, tree = _tree()
    res = simulate_continuous(ct, start=0.0, rate=1.0, seed=SEED)
    (ph.trees.plot(tree, style=_style())
     + ph.trees.color_branches(h.node_values(res))
     + ph.trees.colorbar("trait value", width=240, height=16, size=20)
     + ph.trees.time_axis("time", tick_size=22, label_size=28)).save(out)


def ornstein_uhlenbeck(out):
    ct, tree = _tree()
    res = simulate_continuous(ct, start=6.0, rate=1.0, reverts_to=0.0, pull=2.5, seed=SEED)
    (ph.trees.plot(tree, style=_style())
     + ph.trees.color_branches(h.node_values(res))
     + ph.trees.colorbar("trait value  —  OU: optimum 0, pull 2.5, σ 1, start 6",
                   width=240, height=16, size=20)
     + ph.trees.time_axis("time", tick_size=22, label_size=28)).save(out)


def discrete_states(out):
    ct, tree = _tree()
    res = simulate_discrete(ct, states=("marine", "terrestrial"), switch=0.5, start="marine", seed=5)
    history = {f"n{i}": segs for i, segs in res.history.items()}   # per-branch (state, duration)
    (ph.trees.plot(tree, style=_style())
     + ph.trees.color_history(history, palette={"marine": "#4C9AA6", "terrestrial": "#E08A3C"})
     + ph.trees.legend("habitat")
     + ph.trees.time_axis("time", tick_size=22, label_size=28)).save(out)


# --- the rest of chapter 8's continuous menu: the same diffusion with one argument added ---------
#
# Each is drawn on the SAME tree and seed as `bm`, so a reader can put any of them beside it and the
# only difference is the argument named on the card. The chapter introduces them in this order.


def early_burst(out):
    """The diffusion rate itself decays with time, so the spread is made near the root."""
    ct, tree = _tree()
    res = simulate_continuous(ct, start=0.0, seed=SEED,
                              rate=PerLineage(2.0).changing_at({0: 1.0, 4.5: 0.01}))
    (ph.trees.plot(tree, style=_style())
     + ph.trees.color_branches(h.node_values(res))
     + ph.trees.colorbar("trait value  —  early burst: the diffusion rate falls 100-fold at t = 4.5",
                         width=240, height=16, size=20)
     + ph.trees.time_marker(4.5, color="#111111")   # unlabelled: the key above names the moment
     + ph.trees.time_axis("time", tick_size=22, label_size=28)).save(out)


def regime_optima(out):
    """Two clades pulled to two different optima — the OU with a painting, `regimes=`."""
    ct, tree = _tree()
    painting = simulate_discrete(ct, states=("upland", "lowland"), switch=0.12,
                                 start="upland", seed=3)   # 8 switches: few regimes, each large
    res = simulate_continuous(ct, start=0.0, rate=0.6, regimes=painting,
                              reverts_to={"upland": -4.0, "lowland": 4.0}, pull=3.0, seed=SEED)
    (ph.trees.plot(tree, style=_style())
     + ph.trees.color_branches(h.node_values(res))
     + ph.trees.colorbar("trait value  —  optima: upland -4, lowland 4, pull 3",
                         width=240, height=16, size=20)
     + ph.trees.time_axis("time", tick_size=22, label_size=28)).save(out)


# --- change at the splits: the one trait model a painted tree cannot show on its own -------------
#
# A branch is painted with the value at its far end, so a tree alone cannot say *when* along the
# branch the change happened: jumps at the splits and diffusion along them paint the same picture.
# What separates them is the field's own test. Take the pairs of tips that are each other's closest
# relative and ask how far apart the two have grown. Under diffusion that difference keeps growing
# with the time since they split; under jumps the split itself did the work and the difference does
# not care how long ago it was.
#
# Its own tree, larger than the other cards': the test is one point per pair of sister tips, and a
# hundred-tip tree has only thirty-five such pairs.

_JUMP_TIPS = 300


def speciation_jumps(out):
    """The jump model beside the test that tells it from diffusion."""
    import statistics

    sp = simulate_species_tree(birth=1.0, n_extant=_JUMP_TIPS, seed=SEED)
    ct = sp.complete_tree
    tree = ph.trees.loads(ct.to_newick())
    extant = set(ct.extant_leaves())
    present = max(ct.nodes[i].end_time for i in extant)
    pairs = [(n.end_time, n.children) for n in ct.nodes.values()
             if len(n.children) == 2 and all(c in extant for c in n.children)]

    jumps = simulate_continuous(ct, start=0.0, rate=0.0, at_speciation=1.0, seed=SEED)
    diffusion = simulate_continuous(ct, start=0.0, rate=1.0, seed=SEED)

    png = out.replace(".png", "_tree.png")
    h.render_tree_for_composite(tree, png, h.node_values(jumps), width=1150, height=1150,
                                branch_width=1.0)

    def panel(ax):
        for res, color, name in ((jumps, "#3C8D6E", "jumps at the splits"),
                                 (diffusion, "#B0B8B5", "diffusion along them")):
            xs = [present - t for t, _ in pairs]
            ys = [abs(res.node_values[a] - res.node_values[b]) for _, (a, b) in pairs]
            ax.scatter(xs, ys, s=15, color=color, alpha=0.7, linewidths=0, label=name, zorder=2)
            bins: dict[int, list] = {}
            for x, y in zip(xs, ys):
                bins.setdefault(int(x / 0.2), []).append(y)
            # a bin with a handful of pairs in it says nothing; the tail of the axis is sparse
            pts = sorted((0.2 * k + 0.1, statistics.mean(v)) for k, v in bins.items() if len(v) > 5)
            ax.plot([x for x, _ in pts], [m for _, m in pts], color=color, lw=2.6, zorder=3)
        ax.set_xlabel("time since the two split", fontsize=13)
        ax.set_ylabel("difference between the two", fontsize=13)
        ax.set_ylim(bottom=0.0)
        ax.set_xlim(-0.03, 1.15)                 # beyond this there are only a few pairs left
        ax.legend(fontsize=11, frameon=False, loc="upper left")
        ax.tick_params(labelsize=11)

    h.composite_beside(png, out, panel, figsize=(13, 6.6), ratios=(2.5, 1.6))


# --- the rest of the Literature table: a rate that varies between lineages, one that answers to the
# --- tree's own size, two traits reverting each to its own optimum, and a state read off a liability


def varying_rate_trait(out):
    """The diffusion rate itself varies between lineages, inherited and nudged at each split.

    Same drift law the species tree can put on its birth rate, and the same signature: whole clades
    wander far while their sisters barely move, instead of every branch contributing alike."""
    ct, tree = _tree()
    res = simulate_continuous(
        ct, start=0.0, seed=SEED,
        rate=PerLineage(1.0).varying_among("lineages", Drift(LogNormal(0.0, 0.6))))
    (ph.trees.plot(tree, style=_style())
     + ph.trees.color_branches(h.node_values(res))
     + ph.trees.colorbar("trait value  —  the rate drifts, inherited at each split",
                         width=240, height=16, size=20)
     + ph.trees.time_axis("time", tick_size=22, label_size=28)).save(out)


def diversity_dependent_trait(out):
    """The diffusion rate answers to how full the tree is, so change stops as diversity saturates.

    Its own tree, and a diversity-dependent one: the species tree's birth rate is scaled by the same
    `TotalDiversity` the trait's is, so the tree fills up and levels off instead of growing without
    bound. The curve below shares the time axis and is the number of lineages alive — the quantity
    both rates are divided by. Where it flattens, the colour stops changing."""
    sp = simulate_species_tree(birth=PerLineage(1.4).scaled_by(TotalDiversity(cap=_DD_CAP)),
                               death=0.05, total_time=9.0, seed=5)
    ct = sp.complete_tree
    tree = ph.trees.loads(ct.to_newick())
    res = simulate_continuous(ct, start=0.0, seed=SEED,
                              rate=PerLineage(6.0).scaled_by(TotalDiversity(cap=_DD_CAP)))
    png = out.replace(".png", "_tree.png")
    present = max(ct.nodes[i].end_time for i in ct.extant_leaves())
    (ph.trees.plot(tree, dashed=h.dashed_extinct(tree, ct),
                   style=ph.Style(width=1250, height=830, margin=80, headroom=70, branch_width=1.4))
     + ph.trees.color_branches(h.node_values(res))
     + ph.trees.colorbar("trait value", width=240, height=16, size=20)).save(png)
    times, counts = h.lineages_through_time(ct)

    def panel(ax):
        ax.step(times, counts, where="post", color="#333333", lw=1.8)
        ax.set_ylim(0, _DD_CAP * 1.15)

    h.composite_below(png, present, out, panel, "lineages alive")


def multivariate_ou(out):
    """Two traits drifting together, each pulled to its own optimum at its own strength."""
    ct, _ = _tree()
    res = simulate_continuous(ct, start={"x": 0.0, "y": 0.0}, rate={"x": 1.0, "y": 1.0},
                              correlation={("x", "y"): 0.8},
                              reverts_to={"x": 4.0, "y": -4.0}, pull={"x": 2.0, "y": 0.6},
                              seed=SEED)
    vx = {f"n{i}": v["x"] for i, v in res.node_values.items()}
    vy = {f"n{i}": v["y"] for i, v in res.node_values.items()}
    px, py = out.replace(".png", "_x.png"), out.replace(".png", "_y.png")
    for png, vals in ((px, vx), (py, vy)):
        h.render_tree_for_composite(ph.trees.loads(ct.to_newick()), png, vals,
                                    width=1000, height=520, branch_width=1.2)
    tips = list(ct.extant_leaves())
    xs = [res.node_values[n]["x"] for n in tips]
    ys = [res.node_values[n]["y"] for n in tips]
    h.composite_two_trees_scatter(px, py, xs, ys, out)


#: the threshold model's palette: the state is read off a liability nobody sees
_THRESHOLD = {"small": "#B7C4CC", "large": "#2E5F6E"}


def threshold_trait(out):
    """A discrete state read off a continuous liability: which side of the threshold it sits on.

    The crossings carry no times — the liability is what evolves, and the state is a reading of it —
    so this is the state at each node rather than a painted history, and a branch takes the state its
    far end is in."""
    ct, tree = _tree()
    res = simulate_discrete(ct, states=("small", "large"), liability=1.0, threshold=0.0, seed=9)
    (ph.trees.plot(tree, style=_style())
     + ph.trees.color_branches({f"n{i}": v for i, v in res.node_values.items()},
                               palette=_THRESHOLD)
     + ph.trees.legend("state")
     + ph.trees.time_axis("time", tick_size=22, label_size=28)).save(out)


def asymmetric_switch(out):
    """One direction commoner than the other: a state that is easy to gain and hard to lose."""
    ct, tree = _tree()
    res = simulate_discrete(ct, states=("absent", "present"), start="absent", seed=7,
                            switch={"absent->present": 0.35, "present->absent": 0.03})
    history = {f"n{i}": segs for i, segs in res.history.items()}
    (ph.trees.plot(tree, style=_style())
     + ph.trees.color_history(history, palette={"absent": _ABSENT, "present": "#3C8D6E"})
     + ph.trees.legend("structure")
     + ph.trees.time_axis("time", tick_size=22, label_size=28)).save(out)


def correlated(out):
    ct, _ = _tree()
    res = simulate_continuous(ct, start={"x": 0.0, "y": 0.0}, rate={"x": 1.0, "y": 1.0},
                              correlation={("x", "y"): 0.9}, seed=SEED)
    vx = {f"n{i}": v["x"] for i, v in res.node_values.items()}
    vy = {f"n{i}": v["y"] for i, v in res.node_values.items()}
    px, py = out.replace(".png", "_x.png"), out.replace(".png", "_y.png")
    for png, vals in ((px, vx), (py, vy)):
        h.render_tree_for_composite(ph.trees.loads(ct.to_newick()), png, vals,
                                    width=1000, height=520, branch_width=1.2)
    tips = list(ct.extant_leaves())
    xs = [res.node_values[n]["x"] for n in tips]
    ys = [res.node_values[n]["y"] for n in tips]
    h.composite_two_trees_scatter(px, py, xs, ys, out)


def dependent_characters(out):
    sp = simulate_species_tree(birth=1.0, n_extant=55, seed=7)             # Yule, all lineages extant
    ct = sp.complete_tree
    res = simulate_discrete(ct, states=_DEP_STATES, switch=_DEP_RATES, start="00", seed=3)
    raw = res.history                                                      # {id: [(compound, dur), …]}
    xh = {f"n{i}": [(s[0], d) for s, d in segs] for i, segs in raw.items()}   # project onto X …
    yh = {f"n{i}": [(s[1], d) for s, d in segs] for i, segs in raw.items()}   # … and onto Y
    # two trees, coloured by each character's history (like the continuous "correlated" figure) — you
    # can see Y (purple) tends to be present where X (green) is; the model panel says why.
    px, py = out.replace(".png", "_x.png"), out.replace(".png", "_y.png")
    for png, hist, pal in ((px, xh, _XPAL), (py, yh, _YPAL)):
        (ph.trees.plot(ph.trees.loads(ct.to_newick()),
                       style=ph.Style(width=1000, height=520, margin=45, branch_width=2.6),
                       skeleton=False)
         + ph.trees.color_history(hist, palette=pal)).save(png)
    h.composite_two_trees_panel(px, py, lambda ax: h.draw_grid_markov(ax, _DEP_RATES, _XPAL, _YPAL),
                                out, x_label="character X", y_label="character Y")


def driven_trait(out):
    """A trait driving a trait, with the same ``scaled_by`` a trait uses to drive a genome rate. The
    habitat is grown first and held fixed; body size then diffuses down the same tree, twenty times
    faster on the stretches of branch where the habitat fluctuates.

    The figure is the tree painted by the habitat, with the driver·mapping·target diagram above it.
    What the driving *did* to body size is the Conditioning section's "One trait drives another",
    which paints the same tree twice; repeating it here only made this figure taller."""
    ct = simulate_species_tree(birth=1.0, n_extant=50, seed=7).complete_tree
    hab = simulate_discrete(ct, states=["stable", "fluctuating"], switch=_SWITCH, start="stable", seed=5)
    simulate_continuous(ct, start=0.0, seed=9,           # the target the diagram names
            rate=PerLineage(_BASE).scaled_by(hab, {"fluctuating": _FACTOR, "stable": 1.0}))
    lab = ct.labels()                      # {id: 'n<id>'} — the tree's own names, never built by hand
    tree_png = out.replace(".png", "_tree.png")
    (ph.trees.plot(ph.trees.loads(ct.to_newick()), skeleton=False,
                   # margin clears the time axis and its label: tighter and "time" is cut in half
                   style=ph.Style(width=1320, height=640, margin=96, branch_width=2.6))
     + ph.trees.color_history({lab[i]: segs for i, segs in hab.history.items()},
                              palette=_HABITAT)
     + ph.trees.time_axis("time", tick_size=20, label_size=26, bold=False)).save(tree_png)
    diag = h.conditioning_png(
        out.replace(".png", "_diag.png"),
        driver=("traits", "habitat", "two states"),
        connection=("scaled_by", "table"),
        target_level="traits",
        targets=[("body size", "rate · per lineage", f"fluctuating × {_FACTOR:g}    stable × 1")],
        chain=(("stable", "fluctuating"),
               [(f"{_SWITCH:g}", f"{_SWITCH:g}")],
               (_HABITAT["stable"], _HABITAT["fluctuating"])))
    h.composite_under_diagram(out, diag, [(tree_png, "habitat", _HABITAT)])


_C_BM = '''\
### simulate  —  Brownian motion down a 100-tip tree
from zombi2.species import simulate_species_tree
from zombi2.traits import simulate_continuous

sp = simulate_species_tree(birth=1.0, n_extant=100, seed=42)
res = simulate_continuous(sp.complete_tree, start=0.0, rate=1.0, seed=42)

### plot  —  branches coloured by trait value, with a colour key
import phylustrator as ph

tree = ph.trees.loads(sp.complete_tree.to_newick())
values = {f"n{i}": v for i, v in res.node_values.items()}
(ph.trees.plot(tree)
 + ph.trees.color_branches(values)
 + ph.trees.colorbar("trait value", width=240, height=16, size=20)
 + ph.trees.time_axis("time", tick_size=22, label_size=28)).save("bm.png")'''

_C_OU = '''\
### simulate  —  Ornstein–Uhlenbeck (pulled to an optimum), 100 tips
from zombi2.species import simulate_species_tree
from zombi2.traits import simulate_continuous

sp = simulate_species_tree(birth=1.0, n_extant=100, seed=42)
res = simulate_continuous(sp.complete_tree, start=6.0, rate=1.0, reverts_to=0.0, pull=2.5, seed=42)

### plot  —  the colour key names the OU parameters
import phylustrator as ph

tree = ph.trees.loads(sp.complete_tree.to_newick())
values = {f"n{i}": v for i, v in res.node_values.items()}
(ph.trees.plot(tree)
 + ph.trees.color_branches(values)
 + ph.trees.colorbar("trait value  —  OU: optimum 0, pull 2.5, σ 1, start 6", width=240, height=16, size=20)
 + ph.trees.time_axis("time", tick_size=22, label_size=28)).save("ou.png")'''

_C_DISCRETE = '''\
### simulate  —  a two-state trait hopping between habitats, 100 tips
from zombi2.species import simulate_species_tree
from zombi2.traits import simulate_discrete

sp = simulate_species_tree(birth=1.0, n_extant=100, seed=42)
res = simulate_discrete(sp.complete_tree, states=("marine", "terrestrial"),
                        switch=0.5, start="marine", seed=5)

### plot  —  color_history paints each branch as segments (states change *along* a branch)
import phylustrator as ph

tree = ph.trees.loads(sp.complete_tree.to_newick())
history = {f"n{i}": segs for i, segs in res.history.items()}   # [(state, duration), ...] per branch
(ph.trees.plot(tree)
 + ph.trees.color_history(history, palette={"marine": "#4C9AA6", "terrestrial": "#E08A3C"})
 + ph.trees.legend("habitat")
 + ph.trees.time_axis("time", tick_size=22, label_size=28)).save("discrete.png")'''

_C_CORRELATED = '''\
### simulate  —  two traits evolving together (correlation 0.9)
from zombi2.species import simulate_species_tree
from zombi2.traits import simulate_continuous

sp = simulate_species_tree(birth=1.0, n_extant=100, seed=42)
ct = sp.complete_tree
res = simulate_continuous(ct, start={"x": 0.0, "y": 0.0}, rate={"x": 1.0, "y": 1.0},
                          correlation={("x", "y"): 0.9}, seed=42)

### plot  —  two trees (coloured by x, then y) stacked left; a tip x-vs-y scatter right
import phylustrator as ph
import helpers as h                          # gallery helper: two trees + scatter

vx = {f"n{i}": v["x"] for i, v in res.node_values.items()}
vy = {f"n{i}": v["y"] for i, v in res.node_values.items()}
h.render_tree_for_composite(ph.trees.loads(ct.to_newick()), "tree_x.png", vx)
h.render_tree_for_composite(ph.trees.loads(ct.to_newick()), "tree_y.png", vy)
tips = list(ct.extant_leaves())
xs = [res.node_values[n]["x"] for n in tips]
ys = [res.node_values[n]["y"] for n in tips]
h.composite_two_trees_scatter("tree_x.png", "tree_y.png", xs, ys, "correlated.png")'''


_C_DEPENDENT = '''\
### simulate  —  two binary characters X, Y where Y's flip rate depends on X's state
from zombi2.species import simulate_species_tree
from zombi2.traits import simulate_discrete

# compound Mk states 00/01/10/11; only single-bit flips are named (everything else stays 0).
# gaining Y is slow while X is absent (00->01) but fast once X is present (10->11): the dependence.
rates = {"00->01": 0.05, "01->00": 2.0, "10->11": 2.0, "11->10": 0.05,   # flip Y (rate ~ X)
         "00->10": 0.5,  "10->00": 0.5, "01->11": 0.5, "11->01": 0.5}     # flip X
sp = simulate_species_tree(birth=1.0, n_extant=55, seed=7)
ct = sp.complete_tree
res = simulate_discrete(ct, states=("00", "01", "10", "11"), switch=rates, start="00", seed=3)

### plot  —  two trees (coloured by each character's history) + the 2x2 model chain, like "correlated"
import phylustrator as ph
import helpers as h                           # gallery helper: two trees + a side panel

xpal = {"0": "#c2cac8", "1": "#3C8D6E"}       # X: absent / present (green)
ypal = {"0": "#c2cac8", "1": "#8B6B9E"}       # Y: absent / present (muted purple)
raw = res.history                             # {id: [(compound_state, dur), …]}
xh = {f"n{i}": [(s[0], d) for s, d in segs] for i, segs in raw.items()}   # project onto X
yh = {f"n{i}": [(s[1], d) for s, d in segs] for i, segs in raw.items()}   # project onto Y
for png, hist, pal in (("tree_x.png", xh, xpal), ("tree_y.png", yh, ypal)):
    (ph.trees.plot(ph.trees.loads(ct.to_newick()), skeleton=False)
     + ph.trees.color_history(hist, palette=pal)).save(png)
# right panel = the 2x2 CTMC (the model), kept as the key; Y present tends to track X present
h.composite_two_trees_panel("tree_x.png", "tree_y.png",
                            lambda ax: h.draw_grid_markov(ax, rates, xpal, ypal),
                            "dependent.png", x_label="character X", y_label="character Y")'''


_C_DRIVEN = '''\
### simulate  —  a discrete trait drives a continuous one: the driver first, then the target
from zombi2.species import simulate_species_tree
from zombi2.traits import simulate_continuous, simulate_discrete
from zombi2.params import Drift, LogNormal, PerLineage, TotalDiversity

ct = simulate_species_tree(birth=1.0, n_extant=50, seed=7).complete_tree
hab = simulate_discrete(ct, states=["stable", "fluctuating"], switch=0.4, start="stable", seed=5)
# the same scaled_by that lets a trait drive a genome rate lets it drive another trait: body size
# diffuses 20x faster on the stretches of branch where the habitat fluctuates
size = simulate_continuous(ct, start=0.0, seed=9,
        rate=PerLineage(0.25).scaled_by(hab, {"fluctuating": 20.0, "stable": 1.0}))

### plot  —  the tree coloured by the habitat, the driver that sets body size's diffusion rate
import phylustrator as ph

pal = {"stable": "#3A7CA5", "fluctuating": "#E4572E"}
lab = ct.labels()                                # {id: 'n<id>'}
(ph.trees.plot(ph.trees.loads(ct.to_newick()), skeleton=False)
 + ph.trees.color_history({lab[i]: segs for i, segs in hab.history.items()}, palette=pal)
 + ph.trees.time_axis("time", bold=False)).save("tree.png")
# the figure then composites the driver->mapping->target diagram (habitat -> body size) on top'''


_C_EARLY_BURST = '''\
### simulate  —  the diffusion rate falls 100-fold partway down the tree
from zombi2.species import simulate_species_tree
from zombi2.traits import simulate_continuous
from zombi2.params import Drift, LogNormal, PerLineage, TotalDiversity

sp = simulate_species_tree(birth=1.0, n_extant=100, seed=42)
res = simulate_continuous(sp.complete_tree, start=0.0, seed=42,
                          rate=PerLineage(2.0).changing_at({0: 1.0, 4.5: 0.01}))

### plot  —  the moment the rate drops, marked
import phylustrator as ph

tree = ph.trees.loads(sp.complete_tree.to_newick())
(ph.trees.plot(tree)
 + ph.trees.color_branches({f"n{i}": v for i, v in res.node_values.items()})
 + ph.trees.colorbar("trait value", width=240, height=16, size=20)
 + ph.trees.time_marker(4.5)
 + ph.trees.time_axis("time")).save("early_burst.png")'''

_C_REGIMES = '''\
### simulate  —  a discrete trait paints the tree, and each painted clade has its own optimum
from zombi2.species import simulate_species_tree
from zombi2.traits import simulate_continuous, simulate_discrete

sp = simulate_species_tree(birth=1.0, n_extant=100, seed=42)
painting = simulate_discrete(sp.complete_tree, states=("upland", "lowland"),
                             switch=0.12, start="upland", seed=3)
res = simulate_continuous(sp.complete_tree, start=0.0, rate=0.6, regimes=painting,
                          reverts_to={"upland": -4.0, "lowland": 4.0}, pull=3.0, seed=42)

### plot  —  branches coloured by the trait, so the two optima show as blocks
import phylustrator as ph

tree = ph.trees.loads(sp.complete_tree.to_newick())
(ph.trees.plot(tree)
 + ph.trees.color_branches({f"n{i}": v for i, v in res.node_values.items()})
 + ph.trees.colorbar("trait value", width=240, height=16, size=20)
 + ph.trees.time_axis("time")).save("regimes.png")'''

_C_ASYMMETRIC = '''\
### simulate  —  a structure that is easy to gain and hard to lose
from zombi2.species import simulate_species_tree
from zombi2.traits import simulate_discrete

sp = simulate_species_tree(birth=1.0, n_extant=100, seed=42)
res = simulate_discrete(sp.complete_tree, states=("absent", "present"), start="absent", seed=7,
                        switch={"absent->present": 0.35, "present->absent": 0.03})

### plot  —  each branch painted by the state it was in
import phylustrator as ph

tree = ph.trees.loads(sp.complete_tree.to_newick())
(ph.trees.plot(tree)
 + ph.trees.color_history({f"n{i}": segs for i, segs in res.history.items()},
                          palette={"absent": "#c2cac8", "present": "#3C8D6E"})
 + ph.trees.legend("structure")
 + ph.trees.time_axis("time")).save("asymmetric.png")'''


_C_JUMPS = '''\
### simulate  —  the same tree twice: all the change at the splits, or all of it along the branches
from zombi2.species import simulate_species_tree
from zombi2.traits import simulate_continuous

sp = simulate_species_tree(birth=1.0, n_extant=300, seed=42)
jumps     = simulate_continuous(sp.complete_tree, start=0.0, rate=0.0, at_speciation=1.0, seed=42)
diffusion = simulate_continuous(sp.complete_tree, start=0.0, rate=1.0, seed=42)

### the test  —  for each pair of tips that are each other's closest relative, how far apart they
### have grown against how long ago they split. Diffusion keeps separating them; a jump does not.
ct = sp.complete_tree
extant = set(ct.extant_leaves())
present = max(ct.nodes[i].end_time for i in extant)
pairs = [(n.end_time, n.children) for n in ct.nodes.values()
         if len(n.children) == 2 and all(c in extant for c in n.children)]

for res in (jumps, diffusion):
    x = [present - t for t, _ in pairs]
    y = [abs(res.node_values[a] - res.node_values[b]) for _, (a, b) in pairs]'''


_C_VARYING_RATE = '''\
### simulate  —  the diffusion rate itself varies between lineages, inherited and nudged
from zombi2.species import simulate_species_tree
from zombi2.traits import simulate_continuous
from zombi2.params import Drift, LogNormal, PerLineage

sp = simulate_species_tree(birth=1.0, n_extant=100, seed=42)
res = simulate_continuous(sp.complete_tree, start=0.0, seed=42,
                          rate=PerLineage(1.0).varying_among("lineages", Drift(LogNormal(0.0, 0.6))))
# the same drift law a species tree can put on its birth rate — see the Species section'''

_C_DIVERSITY = '''\
### simulate  —  both rates answer to how full the tree is: the tree's own, and the trait's
from zombi2.species import simulate_species_tree
from zombi2.traits import simulate_continuous
from zombi2.params import PerLineage, TotalDiversity

sp = simulate_species_tree(birth=PerLineage(1.4).scaled_by(TotalDiversity(cap=80)),
                           death=0.05, total_time=9.0, seed=5)
res = simulate_continuous(sp.complete_tree, start=0.0, seed=42,
                          rate=PerLineage(6.0).scaled_by(TotalDiversity(cap=80)))
# the panel below the tree is the lineages-through-time curve: the quantity both rates read'''

_C_MV_OU = '''\
### simulate  —  two traits drifting together, each pulled to its own optimum at its own strength
from zombi2.species import simulate_species_tree
from zombi2.traits import simulate_continuous

sp = simulate_species_tree(birth=1.0, n_extant=100, seed=42)
res = simulate_continuous(sp.complete_tree,
                          start={"x": 0.0, "y": 0.0}, rate={"x": 1.0, "y": 1.0},
                          correlation={("x", "y"): 0.8},
                          reverts_to={"x": 4.0, "y": -4.0},     # each its own optimum
                          pull={"x": 2.0, "y": 0.6},            # and its own strength
                          seed=42)'''

_C_THRESHOLD = '''\
### simulate  —  a discrete state read off a continuous liability
from zombi2.species import simulate_species_tree
from zombi2.traits import simulate_discrete

sp = simulate_species_tree(birth=1.0, n_extant=100, seed=42)
res = simulate_discrete(sp.complete_tree, states=("small", "large"),
                        liability=1.0, threshold=0.0, seed=9)
# the liability is what evolves; the state is which side of `threshold` it is on. The crossings
# carry no times, so this run has no event log — `res.node_values` is the state at each node'''


EXAMPLES = [
    Example("bm", "Brownian motion", "Free diffusion — sister lineages drift apart with time.",
            "rate", brownian_motion, code=_C_BM),
    Example("ou", "Ornstein–Uhlenbeck", "Pulled to an optimum: a high start (yellow) converges to blue.",
            "rate · pull · reverts_to", ornstein_uhlenbeck, code=_C_OU),
    Example("early_burst", "Early burst",
            "The diffusion rate falls a hundred-fold at <code>t&nbsp;=&nbsp;4.5</code>, so the deep "
            "clades separate and then every one of them freezes. "
            "<code>rate&nbsp;=&nbsp;PerLineage(2.0).changing_at({0:&nbsp;1.0,&nbsp;4.5:&nbsp;0.01})</code>.",
            "rate · changing_at", early_burst, code=_C_EARLY_BURST),
    Example("varying_rate", "A rate that varies between lineages",
            "The diffusion rate is inherited at each split and nudged, so whole clades wander while "
            "their sisters barely move — the trait-level reading of the same "
            "<code>Drift</code> a species tree can put on its birth rate.",
            "rate · Drift", varying_rate_trait, code=_C_VARYING_RATE),
    Example("diversity_dependent", "A rate that answers to diversity",
            "The diffusion rate is divided by how full the tree is. The curve below shares the time "
            "axis: where the lineages level off, the colour stops changing.",
            "rate · TotalDiversity", diversity_dependent_trait, code=_C_DIVERSITY),
    Example("regimes", "Two optima, one tree",
            "A discrete trait paints the tree and each painted clade reverts to its own optimum: "
            "upland to &minus;4, lowland to 4. <code>regimes=</code> with a "
            "<code>reverts_to</code> per state.",
            "regimes · reverts_to", regime_optima, code=_C_REGIMES),
    Example("jumps", "Change at the splits",
            "All the change happens <b>at</b> the splits, none along the branches. A painted tree "
            "cannot show that, so the panel runs the test: sister tips that split long ago are no "
            "more different than sister tips that split recently, which is what diffusion cannot do.",
            "at_speciation", speciation_jumps, code=_C_JUMPS),
    Example("discrete", "Discrete states",
            "A two-state trait hops between habitats; each branch is painted by its state history.",
            "switch", discrete_states, code=_C_DISCRETE),
    Example("asymmetric", "Gains and losses at different rates",
            "One direction commoner than the other, written as a matrix of directed rates: the "
            "structure is gained ten times more readily than it is lost, so it spreads and only "
            "rarely goes back.",
            "switch · directed", asymmetric_switch, code=_C_ASYMMETRIC),
    Example("threshold", "A state read off a liability",
            "The state is which side of a threshold a continuous liability sits on, so relatives "
            "flip back and forth near the boundary. The crossings carry no times — the liability is "
            "what evolves — so a branch takes the state its far end is in.",
            "liability · threshold", threshold_trait, code=_C_THRESHOLD),
    Example("correlated", "Dependent continuous traits",
            "Two traits evolve together (r&nbsp;=&nbsp;0.9) — two trees, coloured by each trait, and the "
            "tip scatter.",
            "rate · correlated", correlated, code=_C_CORRELATED),
    Example("mv_ou", "Two traits, two optima",
            "The correlation carries the reversion as well: both traits drift together, and each is "
            "pulled to its own optimum at its own strength — x hard to 4, y gently to &minus;4.",
            "correlation · reverts_to", multivariate_ou, code=_C_MV_OU),
    Example("dependent", "Dependent discrete traits",
            "Two binary characters where one's flip rate depends on the other's state. X in green, Y in "
            "purple, so you can see Y is present where X is; the 2×2 chain is the model.",
            "switch · dependent", dependent_characters, code=_C_DEPENDENT),
]

#: A trait driving a trait is **conditioning**, not a trait model: one trait is grown first and held
#: fixed, and the second reads it through the same `scaled_by` a genome rate would. The code lives
#: here, beside the trait models it is built from; the card belongs in the Conditioning section, and
#: `build.py` puts it there.
CONDITIONING = [
    Example("driven", "A trait driving a trait",
            "A habitat trait sets how fast body size diffuses: twenty times faster where the habitat "
            "fluctuates. <code>rate&nbsp;=&nbsp;PerLineage(0.25).scaled_by(habitat,&nbsp;{…})</code>.",
            "trait → trait", driven_trait, code=_C_DRIVEN),
]
