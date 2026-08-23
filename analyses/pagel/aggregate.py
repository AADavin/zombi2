#!/usr/bin/env python3
"""Rejection rates with Wilson intervals, per arm and per character pair, from fits.tsv.

    python aggregate.py          # writes results.json
"""
from __future__ import annotations

import csv
import json
import math
import os

HERE = os.path.dirname(os.path.abspath(__file__))
ALPHA = 0.05


def wilson(k: int, n: int, z: float = 1.959964) -> tuple[float, float]:
    if n == 0:
        return (0.0, 1.0)
    p = k / n
    d = 1 + z * z / n
    c = p + z * z / (2 * n)
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n))
    return ((c - h) / d, (c + h) / d)


def main() -> int:
    rows = list(csv.DictReader(open(os.path.join(HERE, "fits.tsv")), delimiter="\t"))
    out = {}
    for arm in ("feedback", "trait2gen", "gen2trait", "null"):
        out[arm] = {}
        for pair in ("cox", "ctrl"):
            sub = [r for r in rows if r["arm"] == arm and r["pair"] == pair]
            fitted = [r for r in sub if r["note"] == ""]
            k = sum(1 for r in fitted if float(r["P"]) < ALPHA)
            n = len(fitted)
            lo, hi = wilson(k, n)
            out[arm][pair] = {"fitted": n, "skipped": len(sub) - n,
                              "rejections": k, "rate": round(k / n, 4) if n else None,
                              "wilson95": [round(lo, 4), round(hi, 4)]}
    with open(os.path.join(HERE, "results.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    for arm, pairs in out.items():
        for pair, s in pairs.items():
            print(f"{arm:10s} {pair:4s} n={s['fitted']:3d} skip={s['skipped']:3d} "
                  f"rate={s['rate']} CI={s['wilson95']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
