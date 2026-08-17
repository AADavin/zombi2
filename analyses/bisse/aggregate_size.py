#!/usr/bin/env python3
"""Aggregate the tree-size BiSSE fits into power-vs-size verdicts.

Reads fits_n500.tsv and fits_n1000.tsv (whichever exist), plus verdicts.json for
the 150-tip points, writes verdicts_size.json and prints the table.

    /Users/aadria/miniconda3/bin/python aggregate_size.py
"""

from __future__ import annotations

import json
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent
ALPHA = 0.05


def wilson_ci(k: int, n: int, z: float = 1.959964) -> tuple[float, float]:
    if n == 0:
        return (float("nan"), float("nan"))
    p = k / n
    denom = 1 + z * z / n
    centre = (p + z * z / (2 * n)) / denom
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / denom
    return (max(0.0, centre - half), min(1.0, centre + half))


def read_fits(path: Path) -> list[dict]:
    rows, header = [], None
    for line in path.read_text().splitlines():
        parts = line.rstrip("\n").split("\t")
        if header is None:
            header = parts
            continue
        r = dict(zip(header, parts))
        # tag: n500_f3.0_r012
        bits = r["tag"].split("_")
        r["size"] = int(bits[0][1:])
        r["factor"] = float(bits[1][1:])
        r["p"] = float(r["p"]) if r["p"] not in ("NA", "") else float("nan")
        r["lambda0"] = float(r["lambda0"]) if r["lambda0"] not in ("NA", "") else float("nan")
        r["lambda1"] = float(r["lambda1"]) if r["lambda1"] not in ("NA", "") else float("nan")
        rows.append(r)
    return rows


def main() -> int:
    rows = []
    for name in ("fits_n500.tsv", "fits_n1000.tsv"):
        p = HERE / name
        if p.exists():
            rows.extend(read_fits(p))
    if not rows:
        raise SystemExit("no size fits found")

    v150 = json.loads((HERE / "verdicts.json").read_text())
    verdicts = {"alpha": ALPHA, "cells": {}}
    print(f"{'n_extant':>8} {'factor':>7} {'char':>8} {'fit':>4} {'sig':>4} "
          f"{'rate':>6} {'95% CI':>15}")
    # the 150-tip points come from the existing verdicts
    for f in ("3.0", "5.0"):
        for char in ("driver", "control"):
            cell = v150["per_factor"][f][char]
            verdicts["cells"][f"n150_f{f}_{char}"] = cell
            lo, hi = cell["wilson_ci95"]
            print(f"{150:>8} {f:>7} {char:>8} {cell['n_fit']:>4} "
                  f"{cell['n_significant']:>4} {cell['positive_rate_among_fit']:>6.3f} "
                  f"[{lo:.3f}, {hi:.3f}]")
    for size in (500, 1000):
        for f in (3.0, 5.0):
            for char in ("driver", "control"):
                sel = [r for r in rows if r["size"] == size and r["factor"] == f
                       and r["character"] == char]
                if not sel:
                    continue
                ok = [r for r in sel if r["status"] == "ok"]
                sig = [r for r in ok if r["p"] < ALPHA]
                right = [r for r in sig if r["lambda1"] > r["lambda0"]]
                rate = len(sig) / len(ok) if ok else float("nan")
                lo, hi = wilson_ci(len(sig), len(ok))
                verdicts["cells"][f"n{size}_f{f}_{char}"] = {
                    "n_replicates": len(sel), "n_fit": len(ok),
                    "n_skipped_invariant":
                        sum(r["status"] == "skipped_invariant" for r in sel),
                    "n_error": sum(r["status"].startswith("error") for r in sel),
                    "n_significant": len(sig),
                    "positive_rate_among_fit": rate,
                    "wilson_ci95": [lo, hi],
                    "n_significant_right_direction": len(right)}
                print(f"{size:>8} {f:>7} {char:>8} {len(ok):>4} {len(sig):>4} "
                      f"{rate:>6.3f} [{lo:.3f}, {hi:.3f}]")

    (HERE / "verdicts_size.json").write_text(json.dumps(verdicts, indent=1))
    print(f"\nwrote {HERE / 'verdicts_size.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
