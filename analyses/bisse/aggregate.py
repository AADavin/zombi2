#!/usr/bin/env python3
"""Aggregate the BiSSE fits into the verdicts the paper reports.

Reads fits.tsv (one row per replicate x character) and results.json (for the design),
writes verdicts.json and prints the table.

    /Users/aadria/miniconda3/bin/python aggregate.py
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ALPHA = 0.05


def wilson_ci(k: int, n: int, z: float = 1.959964) -> tuple[float, float]:
    """Wilson score 95% CI for a binomial proportion."""
    if n == 0:
        return (float("nan"), float("nan"))
    p = k / n
    denom = 1 + z * z / n
    centre = (p + z * z / (2 * n)) / denom
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / denom
    return (max(0.0, centre - half), min(1.0, centre + half))


def main() -> int:
    rows = []
    header = None
    for line in (HERE / "fits.tsv").read_text().splitlines():
        parts = line.rstrip("\n").split("\t")
        if header is None:
            header = parts
            continue
        rows.append(dict(zip(header, parts)))
    for r in rows:
        r["factor"] = float(r["tag"].split("_")[0][1:])   # f3.0_r017 -> 3.0
        for k in ("n_tips", "n_present"):
            r[k] = int(r[k])
        for k in ("lnL_full", "lnL_null", "lambda0", "lambda1", "mu0", "mu1",
                  "q01", "q10", "chisq", "p"):
            r[k] = float(r[k]) if r[k] not in ("NA", "") else float("nan")

    factors = sorted({r["factor"] for r in rows})
    verdicts: dict = {"alpha": ALPHA, "per_factor": {}}
    print(f"{'factor':>7} {'char':>8} {'fit':>4} {'skip':>5} {'err':>4} "
          f"{'sig':>4} {'rate':>6} {'95% CI':>15} {'right-dir':>9}")
    for f in factors:
        verdicts["per_factor"][repr(f)] = {}
        for char in ("driver", "control"):
            sel = [r for r in rows if r["factor"] == f and r["character"] == char]
            ok = [r for r in sel if r["status"] == "ok"]
            skipped = [r for r in sel if r["status"] == "skipped_invariant"]
            errors = [r for r in sel if r["status"].startswith("error")]
            sig = [r for r in ok if r["p"] < ALPHA]
            # among the significant fits, how many put the higher speciation rate on
            # the state that actually carries the family (lambda1 > lambda0)?
            right = [r for r in sig if r["lambda1"] > r["lambda0"]]
            n_ok, n_sig = len(ok), len(sig)
            rate = n_sig / n_ok if n_ok else float("nan")
            lo, hi = wilson_ci(n_sig, n_ok)
            verdicts["per_factor"][repr(f)][char] = {
                "n_replicates": len(sel), "n_fit": n_ok,
                "n_skipped_invariant": len(skipped), "n_error": len(errors),
                "n_significant": n_sig,
                "positive_rate_among_fit": rate,
                "wilson_ci95": [lo, hi],
                "n_significant_right_direction": len(right),
                "median_p": float(np.median([r["p"] for r in ok])) if ok else None,
                "median_lambda1_minus_lambda0":
                    float(np.median([r["lambda1"] - r["lambda0"] for r in ok]))
                    if ok else None,
            }
            print(f"{f:>7} {char:>8} {n_ok:>4} {len(skipped):>5} {len(errors):>4} "
                  f"{n_sig:>4} {rate:>6.3f} [{lo:.3f}, {hi:.3f}] "
                  f"{len(right):>4}/{n_sig}")

    # headline numbers for the prose
    null_drv = verdicts["per_factor"]["1.0"]["driver"]
    ctl_all_ok = [r for r in rows if r["character"] == "control" and r["status"] == "ok"]
    ctl_all_sig = [r for r in ctl_all_ok if r["p"] < ALPHA]
    lo, hi = wilson_ci(len(ctl_all_sig), len(ctl_all_ok))
    verdicts["headline"] = {
        "fp_driver_at_null": null_drv,
        "fp_control_pooled_all_factors": {
            "n_fit": len(ctl_all_ok), "n_significant": len(ctl_all_sig),
            "positive_rate": len(ctl_all_sig) / len(ctl_all_ok) if ctl_all_ok else None,
            "wilson_ci95": [lo, hi]},
    }
    ctl_coupled_ok = [r for r in ctl_all_ok if r["factor"] > 1.0]
    ctl_coupled_sig = [r for r in ctl_coupled_ok if r["p"] < ALPHA]
    lo, hi = wilson_ci(len(ctl_coupled_sig), len(ctl_coupled_ok))
    verdicts["headline"]["fp_control_pooled_coupled_only"] = {
        "n_fit": len(ctl_coupled_ok), "n_significant": len(ctl_coupled_sig),
        "positive_rate": len(ctl_coupled_sig) / len(ctl_coupled_ok)
        if ctl_coupled_ok else None,
        "wilson_ci95": [lo, hi]}

    (HERE / "verdicts.json").write_text(json.dumps(verdicts, indent=1))
    print(f"\nwrote {HERE / 'verdicts.json'}")
    h = verdicts["headline"]
    print(f"driver at f=1 (true null): "
          f"{h['fp_driver_at_null']['n_significant']}/{h['fp_driver_at_null']['n_fit']} "
          f"= {h['fp_driver_at_null']['positive_rate_among_fit']:.3f}")
    p = h["fp_control_pooled_coupled_only"]
    print(f"control on coupled trees (f>1): {p['n_significant']}/{p['n_fit']} "
          f"= {p['positive_rate']:.3f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
