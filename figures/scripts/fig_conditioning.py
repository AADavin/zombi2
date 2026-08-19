"""Figure: the three parts of a connection — driver, link, target.

Chapter 8's connection figure, also used in chapter 1, drawn by the **same generator the gallery
uses** (``gallery/helpers.conditioning_png``). It was a hand-authored SVG, which meant the book and
the gallery drew the same diagram twice, from two sources, in two styles — and a convention agreed
for one of them held in the other only by memory. There is one now, so a change to the standard
reaches both.

Run:  python figures/scripts/fig_conditioning.py
"""

from __future__ import annotations

import pathlib
import sys

ROOT = pathlib.Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(ROOT / "gallery"))

import helpers as h                                             # noqa: E402

_HAB = {"aquatic": "#3A7CA5", "terrestrial": "#C9A227"}


def main() -> None:
    out = ROOT / "figures" / "svg" / "conditioning.svg"
    out.parent.mkdir(parents=True, exist_ok=True)
    h.conditioning_png(
        str(out), transparent=False,
        driver=("traits", "habitat", "two states"),
        link=("scaled_by", "table"),
        target_level="genomes",
        targets=[("loss rate", "rate · per copy", "aquatic × 4    terrestrial × 1")],
        chain=(("terrestrial", "aquatic"), [("0.30", "0.30")],
               (_HAB["terrestrial"], _HAB["aquatic"])))
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
