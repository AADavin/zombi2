"""Render the manual's pandoc citations as readable links on the docs site.

A chapter cites the literature in pandoc's form, `[@maliet2019clads]`, because the book is built with
`--citeproc` and that is what turns those keys into citations and a reference list. The site includes
the same chapter verbatim and runs no citeproc, so without this the reader sees the raw key —
`[@maliet2019clads]` — which is exactly what the one pre-existing citation looked like on the website.

Rather than duplicate the references or drop them from the site, this rewrites each key on the way in,
using `docs/references.bib` as the single source both outputs read:

    [@maliet2019clads]                 ->  (Maliet et al. 2019)
    [@gillespie1976; @gillespie1977]   ->  (Gillespie 1976; Gillespie 1977)

Each name links to its DOI, so a citation on the site is a click rather than a dead end.

**An unknown key raises.** A citation that silently renders as literal text is the failure this hook
exists to prevent, so a mistyped key stops the build instead — and CI runs `mkdocs --strict`, which
means a bad reference cannot reach the site.

Like `manual_callouts` and `manual_figures`, this runs as a Markdown *preprocessor* below
`pymdownx.snippets` (priority 32): before the include expands, a chapter page is still just its
`--8<--` line, so there is nothing to rewrite yet.
"""

from __future__ import annotations

import pathlib
import re

from markdown.extensions import Extension
from markdown.preprocessors import Preprocessor

BIB = pathlib.Path(__file__).resolve().parent.parent / "docs" / "references.bib"

#: `[@key]`, or several separated by `;` — pandoc's form, and the only one the chapters use.
CITE = re.compile(r"\[@([^\]]+)\]")
_ENTRY = re.compile(r"@\w+\{([^,]+),(.*?)\n\}", re.S)


def _field(body: str, name: str) -> str:
    # The `(?:\n|$)` matters: the entry pattern captures up to `\n}`, so the LAST field of an entry
    # has no trailing newline. Requiring one silently returned "" for it — and `doi` is conventionally
    # last, so every citation came out unlinked while the build stayed green.
    m = re.search(rf"{name}\s*=\s*[{{\"](.+?)[}}\"],?\s*(?:\n|$)", body, re.S)
    return " ".join(m.group(1).replace("{", "").replace("}", "").split()) if m else ""


def _short(authors: str, year: str) -> str:
    """The author-year label: one name, two joined by `&`, three or more as `et al.`."""
    names = [a.strip() for a in authors.split(" and ") if a.strip()]
    surnames = [n.split(",")[0].strip() if "," in n else n.split()[-1] for n in names]
    if not surnames:
        return year
    if len(surnames) == 1:
        who = surnames[0]
    elif len(surnames) == 2:
        who = f"{surnames[0]} & {surnames[1]}"
    else:
        who = f"{surnames[0]} et al."
    return f"{who} {year}".strip()


def _load() -> dict[str, tuple[str, str]]:
    """`{key: (label, doi url)}` from the shared bibliography."""
    out: dict[str, tuple[str, str]] = {}
    text = BIB.read_text(encoding="utf-8")
    for m in _ENTRY.finditer(text):
        key, body = m.group(1).strip(), m.group(2)
        doi = _field(body, "doi")
        out[key] = (_short(_field(body, "author"), _field(body, "year")),
                    f"https://doi.org/{doi}" if doi else "")
    return out


class _Citations(Preprocessor):
    def __init__(self, md, entries):
        super().__init__(md)
        self.entries = entries

    def run(self, lines: list[str]) -> list[str]:
        # Skip fenced code. A preprocessor runs before the block parser, so a fence is still just
        # lines here; without this, a `[@…]` inside an example would be rewritten as a citation.
        out, fence = [], None
        for line in lines:
            marker = re.match(r"\s*(`{3,}|~{3,})", line)
            if marker:
                token = marker.group(1)[0]
                fence = None if fence and token == fence else (fence or token)
            out.append(line if fence else CITE.sub(self._one, line))
        return out

    def _one(self, m: re.Match) -> str:
        keys = [k.strip().lstrip("@") for k in m.group(1).split(";") if k.strip()]
        parts = []
        for k in keys:
            if k not in self.entries:
                raise KeyError(
                    f"citation [@{k}] is not in {BIB.name}. Add the entry, or fix the key — a key "
                    f"that is not there would render on the site as literal text.")
            label, url = self.entries[k]
            parts.append(f"[{label}]({url})" if url else label)
        return "(" + "; ".join(parts) + ")"


class ManualCitationExtension(Extension):
    def extendMarkdown(self, md):
        md.preprocessors.register(_Citations(md, _load()), "manual_citations", 32)


def on_config(config):
    config.markdown_extensions.append(ManualCitationExtension())
    return config
