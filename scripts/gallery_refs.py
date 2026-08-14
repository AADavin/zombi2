#!/usr/bin/env python
"""Keep the manual's gallery citations — `(Ge7)` — pointing at the right example.

A citation names an example by its **id**, in an HTML comment beside the number:

    …the segment may span the origin ([Ge7](…/gallery.html#genomes)<!--gallery:genome_inversion-->).

The number is not the identity; the id is. The gallery derives `Ge7` from the example's position at
build time, so inserting an example anywhere renumbers everything after it — and this rewrites the
citations to match, from the gallery's own sources rather than from a list kept in step by hand.

    python scripts/gallery_refs.py            # rewrite the citations
    python scripts/gallery_refs.py --check    # fail if any is stale (what CI runs)

The comment is invisible in both outputs: pandoc drops raw HTML on the way to LaTeX, and MkDocs
drops it as a comment. Deleting it by hand only breaks the renumbering, never the page.
"""
from __future__ import annotations

import argparse
import pathlib
import re
import sys

ROOT = pathlib.Path(__file__).resolve().parent.parent
BOOK = ROOT / "manual" / "book"
URL = "https://aadavin.github.io/zombi2/gallery.html"
SECTION = {"Sp": "species", "Ge": "genomes", "Sq": "sequences",
           "Tr": "traits", "Co": "conditioning", "Jo": "joining"}

#: `[Ge7](url#genomes)<!--gallery:the_id-->` or a range, `[Sq1–Sq4](…)<!--…--><!--…-->`
CITATION = re.compile(r"\[(?P<shown>[A-Za-z]{2}\d+(?:[–,] ?[A-Za-z]{2}\d+)*)\]\((?P<url>[^)]*)\)"
                      r"(?P<tags>(?:<!--gallery:[a-z0-9_]+-->)+)")
_TAG = re.compile(r"<!--gallery:([a-z0-9_]+)-->")


def numbers(build=None) -> dict[str, str]:
    """`{example id: number}`, from the gallery's own sources — the same walk `build.py` makes.

    ``build`` is the already-imported module when a caller has one. The test suite does: it imports
    `build` with the drawing packages stubbed, because neither PIL nor Phylustrator is installed in
    the test job and the numbering needs neither."""
    if build is None:
        sys.path.insert(0, str(ROOT / "gallery"))
        import build                                # noqa: PLC0415 — needs the path set first
    out = {}
    for slug, _, _, examples in build.LEVELS:
        for i, ex in enumerate(examples, start=1):
            out[ex.id] = f"{build.PREFIX[slug]}{i}"
    return out


def _shown(got: list[str]) -> str:
    """`Co4`, `Co1–Co6` when they run consecutively, `Co9, Co15` when they do not."""
    if len(got) == 1:
        return got[0]
    prefix = got[0][:2]
    ns = [int(g[2:]) for g in got]
    if all(g[:2] == prefix for g in got) and ns == list(range(ns[0], ns[0] + len(ns))):
        return f"{got[0]}–{got[-1]}"
    return ", ".join(got)


def fix(text: str, nums: dict[str, str]) -> tuple[str, list[str]]:
    """The text with every citation renumbered, and what changed."""
    notes: list[str] = []

    def one(m: re.Match) -> str:
        ids = _TAG.findall(m["tags"])
        missing = [i for i in ids if i not in nums]
        if missing:
            notes.append(f"cites {', '.join(missing)}, which the gallery no longer has")
            return m.group(0)
        got = [nums[i] for i in ids]
        shown = _shown(got)
        anchor = SECTION[got[0][:2]]
        if shown != m["shown"]:
            notes.append(f"{m['shown']} is now {shown}")
        return f"[{shown}]({URL}#{anchor}){m['tags']}"

    return CITATION.sub(one, text), notes


def review(nums: dict[str, str], write: bool) -> tuple[int, list[str]]:
    """Renumber (or just report) every citation in the book. Returns ``(seen, stale)``."""
    stale, seen = [], 0
    for path in sorted(BOOK.glob("*.md")):
        text = path.read_text(encoding="utf-8")
        seen += len(CITATION.findall(text))
        fixed, notes = fix(text, nums)
        if notes:
            stale += [f"{path.name}: {n}" for n in notes]
        if fixed != text and write:
            path.write_text(fixed, encoding="utf-8")
    return seen, stale


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--check", action="store_true", help="report staleness instead of fixing it")
    args = ap.parse_args()

    seen, stale = review(numbers(), write=not args.check)
    if not seen:
        print("no gallery citations found — has the comment been stripped?", file=sys.stderr)
        return 1
    if stale:
        print("\n".join(stale), file=sys.stderr)
        if args.check:
            print(f"\n{len(stale)} stale citation(s). Run: python scripts/gallery_refs.py",
                  file=sys.stderr)
            return 1
        print(f"rewrote {len(stale)} citation(s) of {seen}")
        return 0
    print(f"{seen} gallery citations, all current")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
