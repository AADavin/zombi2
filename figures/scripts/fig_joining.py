"""Figure: a joint run — the conditioning diagram with one more arrow.

Chapter 10's opening figure, drawn by the **same generator** as chapter 9's
(``gallery/helpers.conditioning_png``), because that is the argument. A joint run is not a different
kind of relationship from a conditioned one; it is the same relationship when the driver cannot be
grown first, because what it drives is what it grows along. Everything about the picture is the same
except the arrow coming back, and that is exactly the difference.

Two things drop out when the loop closes. The DRIVER / TARGET headers go, because neither box is only
one of them — each is the other's driver — and the mapping moves under the returning arrow, where
there is room for it.

The retired hand-drawn `joining.svg` claimed the tree was an output in a way that was not true of
every joint model, which is why it went. This says less and is true of all of them.

Run:  python figures/scripts/fig_joining.py
"""

from __future__ import annotations

import pathlib
import shutil
import sys

ROOT = pathlib.Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(ROOT / "gallery"))

import helpers as h                                             # noqa: E402


def main() -> None:
    out = ROOT / "figures" / "svg" / "joining.svg"
    out.parent.mkdir(parents=True, exist_ok=True)
    h.conditioning_png(
        str(out),
        driver=("traits", "body size", "two states"),
        connection=("scaled_by", "table"),
        target_level="species",
        targets=[("birth", "rate · per lineage", "large × 3    small × 1")],
        returns=True)
    # the book and the README both read it from manual/book/figures/, as they do the conditioning one
    book = ROOT / "manual" / "book" / "figures" / "joining.svg"
    shutil.copyfile(out, book)
    print(f"wrote {out}\n      {book}")


if __name__ == "__main__":
    main()
