"""Point the manual's chapter figures at their web originals on the docs site.

A chapter figure is authored once, as an SVG in `manual/book/figures/`, painted with
`var(--ink)`/`var(--paper)` so it follows the reader's colour scheme. The book cannot use that file:
librsvg resolves neither the custom properties nor the `prefers-color-scheme` block, so the Makefile
flattens it and rasterises it to `NAME_print.png`, and the chapter names *that*. The print PNG is a
build artefact — never committed — so the site, which includes the chapter verbatim, would link to a
file that is not there.

Rather than commit the raster or duplicate the figure, this rewrites the reference on the way in:
`figures/NAME_print.png` becomes `../img/NAME.svg`, the symlink in `docs/img/` that points back at
the one authored file. The site then gets the vector, with dark mode working as authored.

Like `manual_callouts`, this has to run as a Markdown *preprocessor* below `pymdownx.snippets`
(priority 32): before the include expands, a chapter page is still just its `--8<--` line.
"""

from __future__ import annotations

import re

from markdown.extensions import Extension
from markdown.preprocessors import Preprocessor

# `![alt](figures/fig-2-1-four-levels_print.png){width=40%}`
#   -> `![alt](../img/fig-2-1-four-levels.svg){ width="40%" }`
# Pandoc's `{width=40%}` is respelled rather than dropped, so the chapter keeps deciding how wide its
# figure is: unquoted, python-markdown's attr_list leaves the braces as literal text on the page.
_FIGURE = re.compile(
    r"\]\(figures/([A-Za-z0-9_-]+)_print\.png\)"  # the print raster the book names
    r"(?:\{ *width *= *([0-9.]+%?) *\})?"  # pandoc's optional width
)

# Just below pymdownx.snippets (32), so chapters are inlined before this runs.
_PRIORITY = 31


def _rewrite(m: re.Match[str]) -> str:
    width = f'{{ width="{m.group(2)}" }}' if m.group(2) else ""
    return f"](../img/{m.group(1)}.svg){width}"


class _FigurePreprocessor(Preprocessor):
    def run(self, lines: list[str]) -> list[str]:
        return [_FIGURE.sub(_rewrite, line) for line in lines]


class _FigureExtension(Extension):
    def extendMarkdown(self, md) -> None:
        md.preprocessors.register(_FigurePreprocessor(md), "manual_figures", _PRIORITY)


def on_config(config):
    config.markdown_extensions.append(_FigureExtension())
    return config
