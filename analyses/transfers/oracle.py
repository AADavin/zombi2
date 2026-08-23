"""What a perfect reconciliation could recover, before any inference tool runs.

Reads the truth tables `experiment.py` wrote and reports two ceilings:

1. **Visibility.** How many true transfers left any extant trace, how many are detectable in
   principle (the gene tree holds both the transferred clade and its donor context), and how
   the same-habitat signal survives extinction plus projection onto the extant tree.

2. **Classification.** The habitat of every internal extant-tree branch, inferred from the
   projected transfer network alone: terminal branches carry their observable tip state, and
   every transfer is one vote (a partner in habitat X votes X). Two classifiers: tips-only
   votes, and label propagation, where classified branches vote in later rounds. Scored
   against the duration-majority truth per branch.

Any reconciliation method scored later sits somewhere below these numbers.

    python oracle.py
"""
from __future__ import annotations

import csv
import pathlib
from collections import Counter, defaultdict

HERE = pathlib.Path(__file__).parent
DATA = HERE / "data"


def read(name):
    with open(DATA / name) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def frac(part, whole):
    return f"{part}/{whole} ({part / whole:.1%})"


def classify(internal, partners, known):
    out = {}
    for rep in internal:
        votes = Counter(known[p] for p in partners.get(rep, ()) if p in known)
        if votes:
            (top, n), *rest = votes.most_common()
            if not rest or n > rest[0][1]:
                out[rep] = top
    return out


def main() -> int:
    transfers = read("transfers.tsv")
    branches = read("branches.tsv")

    # --- visibility ---
    n = len(transfers)
    trace = [t for t in transfers if t["trace"] == "1"]
    det = [t for t in transfers if t["detectable"] == "1"]
    print(f"true transfers:              {n}")
    print(f"trace survives:              {frac(len(trace), n)}")
    print(f"detectable in principle:     {frac(len(det), n)}")
    dead_donor = sum(1 for t in det if t["donor_alive"] == "0")
    dead_recip = sum(1 for t in det if t["recipient_alive"] == "0")
    print(f"dead donor among detectable:     {frac(dead_donor, len(det))}")
    print(f"dead recipient among detectable: {frac(dead_recip, len(det))}")

    def same(rows, a, b):
        pairs = [(t[a], t[b]) for t in rows if t[a] and t[b]]
        return sum(1 for x, y in pairs if x == y) / len(pairs)

    print(f"same-habitat, all true transfers:       {same(transfers, 'donor_habitat', 'recipient_habitat'):.1%}")
    print(f"same-habitat, detectable (true states): {same(det, 'donor_habitat', 'recipient_habitat'):.1%}")
    print(f"same-habitat, detectable (projected):   {same(det, 'donor_proj_habitat', 'recipient_proj_habitat'):.1%}")

    # --- the classification ceiling ---
    truth = {b["branch"]: b["truth_habitat"] for b in branches}
    switched = {b["branch"] for b in branches if b["switched_mid_branch"] == "1"}
    terminal = {b["branch"]: b["tip_state"] for b in branches if b["kind"] == "terminal"}
    internal = [b["branch"] for b in branches if b["kind"] == "internal"]

    partners = defaultdict(list)
    for t in det:
        pd, pr = t["donor_proj"], t["recipient_proj"]
        if pd and pr and pd != pr:
            partners[pd].append(pr)
            partners[pr].append(pd)

    print(f"\ninternal branches: {len(internal)} "
          f"(switched habitat mid-branch: {sum(1 for r in internal if r in switched)})")
    zero = sum(1 for r in internal if not partners.get(r))
    print(f"branches untouched by any detectable transfer: {zero}")

    v1 = classify(internal, partners, dict(terminal))
    known = dict(terminal)
    for _ in range(50):
        merged = dict(terminal)
        merged.update(classify(internal, partners, known))
        if merged == known:
            break
        known = merged
    v2 = {r: known[r] for r in internal if r in known}

    for name, pred in (("tips-only votes ", v1), ("with propagation", v2)):
        hits = sum(1 for r, h in pred.items() if h == truth[r])
        print(f"{name}: classified {frac(len(pred), len(internal))}, "
              f"accuracy among classified {frac(hits, len(pred))}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
