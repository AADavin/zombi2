#!/usr/bin/env python3
"""Generate the GitHub wiki from the canonical MkDocs pages under ``docs/``.

The wiki is a *mirror*, not a second source: edit the pages in ``docs/`` and run this
script manually to regenerate the wiki (it is currently dormant — there is no CI Action
for it, and the MkDocs site already renders the same ``docs/``). It:

* flattens ``docs/``'s nested layout to the wiki's flat page names,
* rewrites intra-doc ``*.md`` links to wiki page names,
* converts MkDocs ``!!! admonition`` blocks to Markdown blockquotes (GitHub-flavoured
  Markdown has no admonitions),
* replaces the mkdocstrings API page (which only renders under MkDocs) with a stub, and
* writes ``Home.md``, ``_Sidebar.md`` and ``_Footer.md``.

Usage::

    python tools/sync_wiki.py <output-dir>      # e.g. a checkout of the *.wiki.git repo
"""

from __future__ import annotations

import os
import re
import sys

# Each entry: (source path relative to docs/, wiki page name, sidebar label, section).
# Order defines the sidebar order. ``section`` groups pages under a heading (None = top).
NAV: list[tuple[str, str, str, str | None]] = [
    ("index.md",                        "Home",                       "Home",                    None),
    ("installation.md",                 "Installation",               "Installation",            None),
    ("quickstart.md",                   "Quickstart",                 "Quickstart",              None),
    ("cli.md",                          "Command-Line-Interface",     "Command-line interface",  None),
    ("guide/species-trees.md",          "Species-Trees",              "Species trees",           "User guide"),
    ("guide/ghost-lineages.md",         "Ghost-Lineages",             "Ghost lineages",          "User guide"),
    ("guide/gene-families.md",          "Gene-Families-and-Rates",    "Gene families & rates",   "User guide"),
    ("guide/transfers.md",              "Transfers",                  "Transfers",               "User guide"),
    ("guide/growth.md",                 "Bounding-Growth",            "Bounding growth",         "User guide"),
    ("guide/gene-trees-and-output.md",  "Gene-Trees-and-Output",      "Gene trees & output",     "User guide"),
    ("guide/coupling.md",               "Gene-Family-Coupling",       "Gene-family coupling",    "User guide"),
    ("guide/ordered-genomes.md",        "Ordered-Genomes",            "Ordered genomes",         "User guide"),
    ("guide/nucleotide-genomes.md",     "Nucleotide-Genomes",         "Nucleotide genomes",      "User guide"),
    ("guide/rate-variation.md",         "Rate-Variation",             "Rate variation",          "User guide"),
    ("guide/traits.md",                 "Trait-Evolution",            "Trait evolution",         "User guide"),
    ("guide/trait-linked-genomes.md",   "Trait-Linked-Gene-Families", "Trait-linked gene families", "User guide"),
    ("guide/parallel.md",               "Running-in-Parallel",        "Running in parallel",     "User guide"),
    ("guide/rust-engine.md",            "Rust-Engine",                "Rust engine",             "User guide"),
    ("guide/extending.md",              "Extending-ZOMBI2",           "Extending ZOMBI2",        "User guide"),
    ("cookbook.md",                     "Cookbook",                   "Cookbook",                None),
    ("faq.md",                          "FAQ",                        "FAQ",                     None),
    ("comparison.md",                   "Comparison-with-ZOMBI-1",    "Comparison with ZOMBI-1", None),
    ("contributing.md",                 "Contributing",               "Contributing",            None),
    ("reference/api.md",                "API-Reference",              "API reference",           None),
    ("species_tree_models.md",          "Species-Tree-Models-Roadmap","Species-tree models",     "Roadmap"),
    ("models/coevolution.md",           "Coevolution-Coupled-Models", "Coevolution (coupled models)", "Roadmap"),
]

# Pages we cannot faithfully mirror (mkdocstrings autodoc renders only under MkDocs).
# We emit a stub for these instead of copying the source.
STUB_PAGES = {"reference/api.md"}

REPO = "https://github.com/AADavin/zombi2"

# docs-relative posix path -> wiki page name, for link rewriting.
PATH_TO_WIKI = {src: wiki for src, wiki, _label, _section in NAV}


def convert_admonitions(text: str) -> str:
    """Turn ``!!! note "Title"`` / ``??? note`` blocks into Markdown blockquotes."""
    lines = text.splitlines()
    out: list[str] = []
    i = 0
    header = re.compile(r'^(\s*)(?:!!!|\?\?\?\+?)\s+(\S+)(?:\s+"([^"]*)")?\s*$')
    while i < len(lines):
        m = header.match(lines[i])
        if not m:
            out.append(lines[i])
            i += 1
            continue
        indent, kind, title = m.group(1), m.group(2), m.group(3)
        label = title if title else kind.capitalize()
        out.append(f"{indent}> **{label}**")
        i += 1
        # Consume the indented body (indented deeper than the marker, or blank).
        body_indent = len(indent) + 4
        block: list[str] = []
        while i < len(lines):
            ln = lines[i]
            if ln.strip() == "":
                block.append("")
                i += 1
                continue
            if len(ln) - len(ln.lstrip()) >= body_indent:
                block.append(ln[body_indent:])
                i += 1
                continue
            break
        while block and block[-1] == "":
            block.pop()
        if block:
            out.append(f"{indent}>")
            for b in block:
                out.append(f"{indent}> {b}".rstrip())
        out.append("")
    return "\n".join(out) + ("\n" if text.endswith("\n") else "")


def rewrite_links(text: str, src: str) -> tuple[str, list[str]]:
    """Rewrite ``[t](target.md#anchor)`` links to flat wiki page names. Returns the new
    text and a list of unresolved ``.md`` targets (for warnings)."""
    src_dir = os.path.dirname(src)
    warnings: list[str] = []

    def repl(match: re.Match) -> str:
        label, target = match.group(1), match.group(2)
        if target.startswith(("http://", "https://", "#", "mailto:")):
            return match.group(0)
        path, _, anchor = target.partition("#")
        if not path.endswith(".md"):
            return match.group(0)  # image/pdf/other — leave as-is
        resolved = os.path.normpath(os.path.join(src_dir, path)).replace(os.sep, "/")
        wiki = PATH_TO_WIKI.get(resolved)
        if wiki is None:
            warnings.append(f"{src}: unresolved link -> {target}")
            return match.group(0)
        return f"[{label}]({wiki}{('#' + anchor) if anchor else ''})"

    new = re.sub(r'\[([^\]]+)\]\(([^)]+)\)', repl, text)
    return new, warnings


def api_stub() -> str:
    return (
        "# API reference\n\n"
        "The full API reference is generated from the source docstrings and is only\n"
        "available in the MkDocs documentation site (it uses `mkdocstrings`, which does\n"
        "not render on the GitHub wiki).\n\n"
        f"Build it locally from the [repository]({REPO}):\n\n"
        "```bash\n"
        'pip install -e ".[docs]"\n'
        "mkdocs serve   # then open the \"API reference\" page\n"
        "```\n"
    )


def build_sidebar() -> str:
    lines = ["### ZOMBI2", ""]
    current_section: str | None = "__top__"
    for _src, wiki, label, section in NAV:
        if section != current_section:
            if section is not None:
                lines.append("")
                lines.append(f"**{section}**")
            current_section = section
        indent = "  " if section is not None else ""
        lines.append(f"{indent}- [{label}]({wiki})")
    lines += ["", "---", "",
              "_This wiki is generated from `docs/`; edit there, not here._"]
    return "\n".join(lines) + "\n"


def build_footer() -> str:
    return (f"ZOMBI2 — simulation of species trees and gene families · "
            f"[repository]({REPO}) · "
            f"successor to [ZOMBI](https://github.com/AADavin/Zombi)\n")


def main(argv: list[str]) -> int:
    if len(argv) != 2:
        print(__doc__)
        return 2
    out_dir = argv[1]
    docs_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "docs")
    os.makedirs(out_dir, exist_ok=True)

    all_warnings: list[str] = []
    written = 0
    for src, wiki, _label, _section in NAV:
        if src in STUB_PAGES:
            content = api_stub()
        else:
            path = os.path.join(docs_dir, src)
            if not os.path.exists(path):
                all_warnings.append(f"missing source page: {src}")
                continue
            with open(path, encoding="utf-8") as f:
                content = f.read()
            content = convert_admonitions(content)
            content, warns = rewrite_links(content, src)
            all_warnings.extend(warns)
        with open(os.path.join(out_dir, f"{wiki}.md"), "w", encoding="utf-8") as f:
            f.write(content)
        written += 1

    with open(os.path.join(out_dir, "_Sidebar.md"), "w", encoding="utf-8") as f:
        f.write(build_sidebar())
    with open(os.path.join(out_dir, "_Footer.md"), "w", encoding="utf-8") as f:
        f.write(build_footer())

    print(f"wrote {written} pages + _Sidebar + _Footer to {out_dir}/")
    for w in all_warnings:
        print(f"  WARN {w}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
