#!/usr/bin/env python3
"""Turn the MkDocs `docs/tools/` pages into manual chapters for the tools PDF.

`docs/tools/*.md` is the single source of truth for the tools reference. It is
authored in MkDocs-Material flavour; the PDF manual is pandoc + XeLaTeX. This
script bridges the two so the tools PDF regenerates from the docs with no
hand-maintained copy to drift:

  * MkDocs admonitions  `!!! note "Title"` (+ warning/tip)  ->  pandoc fenced
    divs `::: note` (mapped to callout boxes by callouts.lua), with the custom
    title kept as a bold lead-in.
  * figures  `![alt](../img/NAME.svg)`  ->  `![alt](figures/NAME.pdf)` (the
    Makefile rasterises the SVG next to the manual's other figures).
  * cross-links: a link to a sibling tools page `(page.md)` becomes an
    intra-PDF jump `#<that-chapter-anchor>`; a link that escapes the tools
    section `(../guide/x.md)` becomes an absolute URL on the live docs site.

Everything outside those three rewrites is passed through verbatim, and nothing
inside fenced code blocks is touched.

Usage:
    tools_to_chapters.py --out DIR  src1.md src2.md ...

The output files are `NN-<basename>.md`, numbered in the order given so a plain
`DIR/*.md` glob assembles them in nav order.
"""

from __future__ import annotations

import argparse
import os
import re
import sys

# Base URL of the live documentation (mkdocs.yml `site_url`). Links that leave
# the tools section point here so they stay clickable in the PDF.
SITE = "https://aadavin.github.io/zombi2/docs"

_FENCE = re.compile(r"^(```|~~~)")
_ADMON = re.compile(r'^(\s*)!!!\s+(note|warning|tip)\b(?:\s+"([^"]*)")?\s*$')
_LINK = re.compile(r"(!?)\[([^\]]*)\]\(([^)]+)\)")
_HEADING = re.compile(r"^(#{1,6})\s+(.*?)\s*$")


def pandoc_id(text: str) -> str:
    """Reproduce pandoc's `auto_identifiers` for a heading's text.

    Keep letters/digits/`_-.`, turn whitespace runs into single hyphens,
    lowercase, then drop everything before the first letter.
    """
    kept = "".join(ch for ch in text if ch.isalnum() or ch in "_-. \t")
    kept = re.sub(r"[\s]+", "-", kept).lower()
    m = re.search(r"[a-z]", kept)
    kept = kept[m.start():] if m else ""
    return kept or "section"


def first_heading(path: str) -> str:
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("# "):
                return line[2:].strip()
    raise SystemExit(f"{path}: no level-1 heading found")


def transform_url(url: str, anchors: dict[str, str]) -> str:
    """Rewrite one link/image target for the PDF (see module docstring)."""
    if re.match(r"^(https?:|mailto:|#)", url):
        return url

    path, frag = (url.split("#", 1) + [None])[:2] if "#" in url else (url, None)

    # figures: ../img/NAME.svg -> figures/NAME.pdf (defensive; images also hit this)
    m = re.match(r"^\.\./img/(.+)\.svg$", path)
    if m:
        return "figures/" + os.path.basename(m.group(1)) + ".pdf"

    # a link that leaves tools/ (../guide/x.md, ../contributing/y.md) -> live docs URL
    m = re.match(r"^\.\./(.+)\.md$", path)
    if m:
        return f"{SITE}/{m.group(1)}/" + (f"#{frag}" if frag else "")

    # a sibling tools page (page.md) -> intra-PDF anchor
    m = re.match(r"^([A-Za-z0-9_.-]+)\.md$", path)
    if m:
        if frag:  # explicit section anchor wins
            return f"#{frag}"
        base = m.group(1)
        if base in anchors:
            return f"#{anchors[base]}"
    return url


def rewrite_links(line: str, anchors: dict[str, str], remap: dict[str, str]) -> str:
    def repl(m: re.Match) -> str:
        bang, text, url = m.groups()
        new = transform_url(url.strip(), anchors)
        if new.startswith("#") and new[1:] in remap:  # same-page link -> unique id
            new = "#" + remap[new[1:]]
        return f"{bang}[{text}]({new})"

    return _LINK.sub(repl, line)


def same_page_remap(text: str, stem: str) -> dict[str, str]:
    """Same-page `#frag` links collide once chapters merge (many pages share a
    "Validation" section). For each such fragment that names a heading in *this*
    file, mint a file-unique id `<stem>-<frag>` used for both the heading and the
    link, so the jump stays within the chapter.
    """
    in_fence, frags, ids = False, set(), set()
    for line in text.split("\n"):
        if _FENCE.match(line):
            in_fence = not in_fence
            continue
        if in_fence:
            continue
        h = _HEADING.match(line)
        if h:
            ids.add(pandoc_id(h.group(2)))
        for m in _LINK.finditer(line):
            url = m.group(3).strip()
            if url.startswith("#"):
                frags.add(url[1:])
    return {f: f"{stem}-{f}" for f in frags & ids}


def convert(text: str, anchors: dict[str, str], stem: str = "") -> str:
    remap = same_page_remap(text, stem)
    lines = text.split("\n")
    out: list[str] = []
    i, in_fence = 0, False
    while i < len(lines):
        line = lines[i]
        if _FENCE.match(line):
            in_fence = not in_fence
            out.append(line)
            i += 1
            continue
        if not in_fence:
            heading = _HEADING.match(line)
            if heading and pandoc_id(heading.group(2)) in remap:
                line = f"{line.rstrip()} {{#{remap[pandoc_id(heading.group(2))]}}}"
                out.append(line)
                i += 1
                continue
            m = _ADMON.match(line)
            if m:
                indent, kind, title = m.group(1), m.group(2), m.group(3)
                base = len(indent)
                body: list[str] = []
                i += 1
                while i < len(lines):
                    bl = lines[i]
                    if bl.strip() == "":
                        body.append("")
                        i += 1
                        continue
                    if (len(bl) - len(bl.lstrip(" "))) > base:  # still indented => in block
                        body.append(bl[base + 4:] if len(bl) >= base + 4 else bl.lstrip(" "))
                        i += 1
                    else:
                        break
                while body and body[-1] == "":
                    body.pop()
                out.append(f"::: {kind}")
                if title:
                    out.append(f"**{title}**")
                    out.append("")
                out.extend(rewrite_links(b, anchors, remap) for b in body)
                out.append(":::")
                out.append("")
                continue
            line = rewrite_links(line, anchors, remap)
        out.append(line)
        i += 1
    return "\n".join(out)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", required=True, help="output directory for the chapter files")
    ap.add_argument("sources", nargs="+", help="docs/tools/*.md in nav order")
    args = ap.parse_args()

    os.makedirs(args.out, exist_ok=True)
    # basename (no extension) -> the pandoc anchor of that page's H1
    anchors = {
        os.path.splitext(os.path.basename(p))[0]: pandoc_id(first_heading(p))
        for p in args.sources
    }

    for n, src in enumerate(args.sources):
        stem = os.path.splitext(os.path.basename(src))[0]
        with open(src, encoding="utf-8") as fh:
            converted = convert(fh.read(), anchors, stem)
        dst = os.path.join(args.out, f"{n:02d}-{stem}.md")
        with open(dst, "w", encoding="utf-8") as fh:
            fh.write(converted)
        print(f"  {src} -> {dst}", file=sys.stderr)


if __name__ == "__main__":
    main()
