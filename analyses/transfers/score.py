"""Score ALE's reconciliations against the truth tables.

Reads `recs/fam_*.nwk.ale.uTs` (ALE's transfer calls, `from  to  freq`) and one
`.uml_rec` (its S: line carries ALE's species tree, whose internal ids are mapped back to
our branch names by clade content against `data/branches.tsv`). Reports:

1. **Transfer recovery.** For every detectable true transfer, the truth says which pair of
   extant-tree branches a perfect method would name (the projections). Recall: how many of
   those pairs ALE calls; precision: how many of ALE's calls are true pairs. Both at several
   frequency thresholds, since a `.uTs` row is a sampled frequency, not a yes/no.

2. **The habitat network.** The same-habitat share of ALE's inferred network (endpoints
   scored by their branch's majority habitat) next to the oracle's projected ceiling, and
   the vote classifier (tips-only, then propagation, votes weighted by ALE's frequencies)
   scored against the true habitat of every internal branch.

    python score.py
"""
from __future__ import annotations

import csv
import pathlib
import re
from collections import Counter, defaultdict

HERE = pathlib.Path(__file__).parent
DATA, RECS = HERE / "data", HERE / "recs"


# ---------- a tiny newick parser for ALE's S: line ----------
def parse_s_tree(uml_rec_path):
    """{ALE node name: frozenset of tip names} for every node of ALE's species tree."""
    s_line = next(line for line in open(uml_rec_path) if line.startswith("S:"))
    nwk = s_line.split("\t", 1)[1].strip().rstrip(";")
    clades, pos = {}, 0

    def take_name():
        nonlocal pos
        m = re.match(r"[^,():;]*", nwk[pos:])
        name = m.group(0)
        pos += len(name)
        m = re.match(r":[0-9.eE+-]+", nwk[pos:])   # the branch length, if any
        if m:
            pos += len(m.group(0))
        return name

    def parse():
        nonlocal pos
        if nwk[pos] == "(":
            pos += 1
            below = parse()
            while nwk[pos] == ",":
                pos += 1
                below |= parse()
            assert nwk[pos] == ")"
            pos += 1
            name = take_name()
        else:
            name = take_name()
            below = {name}
        clades[name] = frozenset(below)
        return set(below)

    parse()
    return clades


def main() -> int:
    branches = list(csv.DictReader(open(DATA / "branches.tsv"), delimiter="\t"))
    by_clade = {frozenset(b["clade"].split(";")): b["branch"] for b in branches}
    truth_hab = {b["branch"]: b["truth_habitat"] for b in branches}
    tip_state = {b["branch"]: b["tip_state"] for b in branches if b["kind"] == "terminal"}
    internal = [b["branch"] for b in branches if b["kind"] == "internal"]

    one_rec = sorted(RECS.glob("*.uml_rec"))[0]
    ale_to_ours = {}
    for name, clade in parse_s_tree(one_rec).items():
        ours = by_clade.get(clade)
        if ours is not None:
            ale_to_ours[name] = ours
    n_unmapped = sum(1 for b in branches if b["branch"] not in set(ale_to_ours.values()))
    print(f"ALE species-tree nodes mapped to our branches: {len(ale_to_ours)} "
          f"(unmapped on our side: {n_unmapped})")

    # --- ALE's calls: {family: {(from, to): freq}} in OUR branch names ---
    calls = {}
    for path in sorted(RECS.glob("*.uTs")):
        fam = int(re.search(r"fam_(\d+)", path.name).group(1))
        rows = {}
        for line in open(path):
            if line.startswith("#") or not line.strip():
                continue
            src, dst, freq = line.split()
            a, b = ale_to_ours.get(src), ale_to_ours.get(dst)
            if a is not None and b is not None:
                rows[(a, b)] = rows.get((a, b), 0.0) + float(freq)
        calls[fam] = rows
    print(f"families with reconciliations: {len(calls)}")

    # --- 1. transfer recovery against the projected truth ---
    transfers = list(csv.DictReader(open(DATA / "transfers.tsv"), delimiter="\t"))
    truth_pairs = defaultdict(set)
    for t in transfers:
        if t["detectable"] == "1" and t["donor_proj"] and t["recipient_proj"] \
                and int(t["family"]) in calls:
            truth_pairs[int(t["family"])].add((t["donor_proj"], t["recipient_proj"]))
    n_truth = sum(len(v) for v in truth_pairs.values())
    print(f"detectable true (donor,recipient) pairs in scored families: {n_truth}")

    for thr in (0.1, 0.3, 0.5, 0.8):
        called = {(f, p) for f, rows in calls.items() for p, q in rows.items() if q >= thr}
        true_set = {(f, p) for f, ps in truth_pairs.items() for p in ps}
        tp = len(called & true_set)
        # a call with donor and recipient swapped still finds the event, count separately
        swapped = len({(f, (b, a)) for f, (a, b) in called} & true_set - {
            (f, p) for f, p in called if (f, p) in true_set})
        prec = tp / len(called) if called else 0.0
        rec = tp / len(true_set)
        print(f"freq >= {thr}: calls {len(called):5d}  precision {prec:.1%}  "
              f"recall {rec:.1%}  (+{swapped} more as swapped-direction hits)")

    # --- 2. the habitat network ALE sees ---
    weight = defaultdict(float)
    partners = defaultdict(list)          # branch -> [(partner, freq)]
    for fam, rows in calls.items():
        for (a, b), q in rows.items():
            weight[(a, b)] += q
            partners[a].append((b, q))
            partners[b].append((a, q))
    same = sum(q for (a, b), q in weight.items() if truth_hab[a] == truth_hab[b])
    total = sum(weight.values())
    print(f"\nALE network: {total:.0f} transfers (frequency-weighted), "
          f"same-habitat share {same / total:.1%}")

    def classify(known):
        out = {}
        for rep in internal:
            votes = Counter()
            for p, q in partners.get(rep, ()):
                if p in known:
                    votes[known[p]] += q
            if votes:
                (top, n), *rest = votes.most_common()
                if not rest or n > rest[0][1]:
                    out[rep] = top
        return out

    v1 = classify(dict(tip_state))
    known = dict(tip_state)
    for _ in range(50):
        merged = dict(tip_state)
        merged.update(classify(known))
        if merged == known:
            break
        known = merged
    v2 = {r: known[r] for r in internal if r in known}

    for name, pred in (("tips-only votes ", v1), ("with propagation", v2)):
        hits = sum(1 for r, h in pred.items() if h == truth_hab[r])
        print(f"{name}: classified {len(pred)}/{len(internal)} internal branches "
              f"({len(pred) / len(internal):.0%}), accuracy among classified "
              f"{hits}/{len(pred)} ({hits / len(pred):.0%})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
