"""The figure for the transfers example: world, visible truth, evidence, verdict.

Reads ``data/`` (the truth tables from experiment.py) and ``recs/`` (ALE's reconciliations,
pulled from the cluster); every number in the figure is derived from them at render time.
Writes ``figures/transfers.png``.

    python figures.py

A  Complete tree      the truth with its extinct lineages
B  Extant tree        the true history on the sampled tree, switching mid-branch
C  Vote tree          the evidence: habitat read from transfer partners (+ flow matrix inset D)
E  ALE inference      green/red against the exact truth, accuracy-through-time aligned below
"""
import csv
import pathlib
import re
import sys
from collections import Counter, defaultdict

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import phylustrator as ph
from phylustrator import color as ph_color
from PIL import Image, ImageChops, ImageDraw, ImageFont

HERE = pathlib.Path(__file__).parent
sys.path.insert(0, str(HERE))
from experiment import build_world, Projection  # noqa: E402
from score import parse_s_tree  # noqa: E402

S = HERE / "figures"
S.mkdir(exist_ok=True)
PAL = {"A": "#1F6E8C", "B": "#E3A857"}   # petrol / amber
OK, BAD = "#2E6B45", "#C2453C"
INK, MUTED = "#1a1a1a", "#8a8a8a"
LEG = dict(size=24, inset=56, dy=20)
ph_color._COLORMAPS["habvote"] = [(227, 168, 87), (238, 238, 238), (31, 110, 140)]

sp, ct, tr, g = build_world()
proj = Projection(ct, tr)
et = sp.extant_tree
et_lab = et.labels()
lab = ct.labels()
nw_nostem = re.sub(r"\)(n\d+):[0-9.eE+-]+;$", r")\1;", et.to_newick())

branches = list(csv.DictReader(open(HERE / "data/branches.tsv"), delimiter="\t"))
tip_state = {b["branch"]: b["tip_state"] for b in branches if b["kind"] == "terminal"}
truth_hab = {b["branch"]: b["truth_habitat"] for b in branches}
internal = [b["branch"] for b in branches if b["kind"] == "internal"]
by_clade = {frozenset(b["clade"].split(";")): b["branch"] for b in branches}
one_rec = sorted((HERE / "recs").glob("*.uml_rec"))[0]
ale_to_ours = {n: by_clade[c] for n, c in parse_s_tree(one_rec).items() if c in by_clade}

votes = defaultdict(lambda: [0.0, 0.0])
flow = defaultdict(float)
partners = defaultdict(list)
for path in (HERE / "recs").glob("*.uTs"):
    for line in open(path):
        if line.startswith("#") or not line.strip():
            continue
        src, dst, freq = line.split()
        a, b = ale_to_ours.get(src), ale_to_ours.get(dst)
        if a is None or b is None:
            continue
        q = float(freq)
        flow[(truth_hab[a], truth_hab[b])] += q
        partners[a].append((b, q))
        partners[b].append((a, q))
        for me, other in ((a, b), (b, a)):
            if other in tip_state:
                votes[me][0 if tip_state[other] == "A" else 1] += q

ale_state = dict(tip_state)
for rep in internal:
    v = Counter()
    for p, q in partners.get(rep, ()):
        if p in tip_state:
            v[tip_state[p]] += q
    if v:
        (top, n), *rest = v.most_common()
        if not rest or n > rest[0][1]:
            ale_state[rep] = top

hist_true, spans = {}, {}
for nid in et.nodes:
    chain = proj.chain(nid)
    t0 = ct.nodes[chain[-1]].birth_time
    spans[nid] = (t0, ct.nodes[nid].end_time)
    segs = []
    for n in reversed(chain):
        node = ct.nodes[n]
        for state, local_end in tr.history[n]:
            end_abs = min(node.birth_time + local_end, node.end_time)
            if end_abs <= t0:
                continue
            if segs and segs[-1][0] == state:
                segs[-1] = (state, end_abs - t0)
            else:
                segs.append((state, end_abs - t0))
    hist_true[et_lab[nid]] = segs

# ---------- panel renders ----------
hist_ct = {lab[i]: segs for i, segs in tr.history.items()}
(ph.trees.plot(ph.trees.loads(ct.to_newick()), skeleton=False,
               style=ph.Style(width=1500, height=625, branch_width=3.1))
 + ph.trees.color_history(hist_ct, palette=PAL)
 + ph.trees.legend(entries={"habitat A": PAL["A"], "habitat B": PAL["B"]}, **LEG)
 ).save(S / "f_a.png")

(ph.trees.plot(ph.trees.loads(et.to_newick()), skeleton=False,
               style=ph.Style(width=1000, height=450, branch_width=2.6))
 + ph.trees.color_history(hist_true, palette=PAL)).save(S / "f_b.png")

vals = {}
for nid in et.nodes:
    v = votes.get(str(nid))
    if v and sum(v) > 0:
        vals[et_lab[nid]] = v[0] / (v[0] + v[1])
(ph.trees.plot(ph.trees.loads(nw_nostem), skeleton=False,
               style=ph.Style(width=1000, height=450, branch_width=2.6))
 + ph.trees.color_branches(vals, cmap="habvote", limits=(0.0, 1.0))
 + ph.trees.colorbar(size=24, width=260, inset=40,
                     labels=("all B", "all A"))).save(S / "f_c.png")

hist_ok = {}
for nid in et.nodes:
    key, name = str(nid), et_lab[nid]
    if key not in ale_state or name not in hist_true:
        continue
    segs = []
    for state, end in hist_true[name]:
        verdict = "ok" if state == ale_state[key] else "bad"
        if segs and segs[-1][0] == verdict:
            segs[-1] = (verdict, end)
        else:
            segs.append((verdict, end))
    hist_ok[name] = segs
(ph.trees.plot(ph.trees.loads(nw_nostem), skeleton=False,
               style=ph.Style(width=1500, height=625, branch_width=3.1))
 + ph.trees.color_history(hist_ok, palette={"ok": OK, "bad": BAD})
 + ph.trees.legend(entries={"matches the truth": OK, "differs": BAD}, **LEG)
 ).save(S / "f_e.png")

# the flow matrices: realized transfers and ALE's inferred mass, shared scale, Blues
true_M = np.zeros((2, 2))
for t in csv.DictReader(open(HERE / "data/transfers.tsv"), delimiter="\t"):
    true_M[0 if t["donor_habitat"] == "A" else 1,
           0 if t["recipient_habitat"] == "A" else 1] += 1
ale_M = np.array([[flow[("A", "A")], flow[("A", "B")]],
                  [flow[("B", "A")], flow[("B", "B")]]])
fmax = max((true_M / true_M.sum()).max(), (ale_M / ale_M.sum()).max())
for fname, title, M in (("f_d1.png", "realized transfers", true_M),
                        ("f_d2.png", "inferred by ALE", ale_M)):
    fig, ax = plt.subplots(figsize=(3.0, 2.9))
    F = M / M.sum()
    ax.imshow(F, cmap="Blues", vmin=0, vmax=fmax * 1.5)
    for i in range(2):
        for j in range(2):
            ax.text(j, i, f"{M[i, j]:.0f}\n({F[i, j]:.0%})", ha="center",
                    va="center", fontsize=12, color=INK)
    ax.set_xticks([0, 1]); ax.set_yticks([0, 1])
    ax.set_xticklabels(["A", "B"], fontsize=13)
    ax.set_yticklabels(["A", "B"], fontsize=13)
    for t, c in zip(ax.get_xticklabels(), [PAL["A"], PAL["B"]]):
        t.set_color(c); t.set_fontweight("bold")
    for t, c in zip(ax.get_yticklabels(), [PAL["A"], PAL["B"]]):
        t.set_color(c); t.set_fontweight("bold")
    ax.set_xlabel("recipient habitat", fontsize=11)
    ax.set_ylabel("donor habitat", fontsize=11)
    ax.set_title(title, fontsize=13)
    ax.spines[:].set_visible(False)
    ax.tick_params(length=0)
    fig.tight_layout()
    fig.savefig(S / fname, dpi=250)
    plt.close(fig)

# the accuracy curve, in the tree's exact pixel frame (margins 50 style = 100 px at 2x)
root_nid = next(nid for nid in et.nodes
                if (et.nodes[nid].parent if isinstance(et.nodes[nid].parent, (int, type(None)))
                    else et.nodes[nid].parent.id) not in et.nodes)
t_start, t_end = ct.nodes[root_nid].end_time, max(e for _, e in spans.values())
ts = np.linspace(t_start, t_end - 1e-9, 400)
share, alive = [], []
for t in ts:
    ok = tot = 0
    for nid, (t0, t1) in spans.items():
        if t0 <= t <= t1 and str(nid) in ale_state:
            tot += 1
            ok += proj.habitat_of_branch_at(nid, t) == ale_state[str(nid)]
    share.append(ok / tot if tot else np.nan)
    alive.append(tot)
fig = plt.figure(figsize=(15.0, 2.5), dpi=200)
ax = fig.add_axes([100 / 3000, 0.30, 2800 / 3000, 0.62])
ax.plot(ts, np.array(share) * 100, color=INK, linewidth=1.8)
ax2 = ax.twinx()
ax2.plot(ts, alive, color=MUTED, linewidth=1.0, linestyle=(0, (1, 2)))
for a in (ax, ax2):
    a.set_xlim(t_start, t_end)
    a.spines["top"].set_visible(False)
ax.set_ylim(0, 105)
ax.axhline(100, color=MUTED, linewidth=0.7, linestyle=(0, (2, 2)))
ax.set_xlabel("Time", fontsize=10)
ax.tick_params(labelsize=8)
ax2.tick_params(axis="y", colors=MUTED, labelsize=8)
ax.annotate("lineages correctly labelled (%)", (0.012, 0.84), xycoords="axes fraction",
            fontsize=9.5, color=INK)
ax.annotate("lineages alive", (0.012, 0.70), xycoords="axes fraction",
            fontsize=9.5, color=MUTED)
fig.savefig(S / "f_curve.png")
plt.close(fig)

# ---------- composition ----------
def trim(im, bg=(255, 255, 255)):
    d = ImageChops.difference(im.convert("RGB"), Image.new("RGB", im.size, bg))
    return im.crop(d.getbbox())

fa, fb = trim(Image.open(S / "f_a.png")), trim(Image.open(S / "f_b.png"))
fc = trim(Image.open(S / "f_c.png"))
fd1, fd2 = trim(Image.open(S / "f_d1.png")), trim(Image.open(S / "f_d2.png"))
fe, fl = Image.open(S / "f_e.png"), Image.open(S / "f_curve.png")

gap = 70
W = 3000
left_w = 1930                       # B and C stacked here; D takes the right third
fa = fa.resize((W, int(fa.height * W / fa.width)))
fb = fb.resize((left_w, int(fb.height * left_w / fb.width)))
fc = fc.resize((left_w, int(fc.height * left_w / fc.width)))
fd_w = W - left_w - gap - 240
fd1 = fd1.resize((fd_w, int(fd1.height * fd_w / fd1.width)))
fd2 = fd2.resize((fd_w, int(fd2.height * fd_w / fd2.width)))

f_letter = ImageFont.truetype("/System/Library/Fonts/Helvetica.ttc", 56)
f_title = ImageFont.truetype("/System/Library/Fonts/Helvetica.ttc", 46)
title_h = 78
row2_h = title_h + fb.height + 50 + title_h + fc.height
H = (title_h + fa.height) + 60 + row2_h + 60 + (title_h + fe.height + fl.height) + 40
canvas = Image.new("RGB", (3040, H), "white")
draw = ImageDraw.Draw(canvas)

y = 0
draw.text((40, y + 10), "A", fill=INK, font=f_letter)
draw.text((132, y + 18), "Complete tree", fill=INK, font=f_title)
y += title_h
canvas.paste(fa, (20, y))
y += fa.height + 60

row2_top = y
draw.text((40, y + 4), "B", fill=INK, font=f_letter)
draw.text((132, y + 12), "Extant tree", fill=INK, font=f_title)
y += title_h
canvas.paste(fb, (20, y))
y += fb.height + 50
draw.text((40, y + 4), "C", fill=INK, font=f_letter)
draw.text((132, y + 12), "Evidence from transfers", fill=INK, font=f_title)
y += title_h
canvas.paste(fc, (20, y))
y += fc.height + 60

xd = 20 + left_w + gap + 120
stack_h = fd1.height + 50 + fd2.height
yd = row2_top + (row2_h - stack_h) // 2
canvas.paste(fd1, (xd, yd))
canvas.paste(fd2, (xd, yd + fd1.height + 50))
draw.text((xd - 75, yd + 4), "D", fill=INK, font=f_letter)

draw.text((40, y + 4), "E", fill=INK, font=f_letter)
draw.text((132, y + 12), "ALE inference", fill=INK, font=f_title)
y += title_h
canvas.paste(fe, (20, y))
canvas.paste(fl, (20, y + fe.height))

canvas.save(S / "transfers.png")
print("saved", canvas.size)
