"""Render the manual's callouts on the docs site.

The chapters are written for pandoc, where a callout is a fenced div — `::: note` … `:::` — which
`manual/callouts.lua` maps to a LaTeX box. MkDocs has no such syntax: its admonitions are `!!! note`
with an indented body. There is no spelling both toolchains accept, so the single source stays
pandoc-native and this bridges it on the way in.

It has to run as a Markdown *preprocessor* rather than a `on_page_markdown` hook: that hook fires
before the extensions do, when a chapter page is still just its `--8<--` line, so the chapter text
does not exist yet to convert. Registering below `pymdownx.snippets` (priority 32) means the include
has been expanded by the time this sees the page.

Fixing the collision matters beyond looks. `mkdocstrings` claims *every* `::: identifier` line as an
autodoc block, so an unconverted `::: note` makes it try to import a module named `note` and abort
the build. Only the classes `callouts.lua` knows are rewritten; anything else starting with `:::` is
passed through untouched, so a genuine `::: zombi2.species` block still reaches its handler.
"""

from __future__ import annotations

import re

from markdown.extensions import Extension
from markdown.preprocessors import Preprocessor

# Kept in step with the `map` in manual/callouts.lua — the callout classes the book actually uses.
CALLOUT_CLASSES = ("note", "warning", "tip")

_OPEN = re.compile(r"^::: *(" + "|".join(CALLOUT_CLASSES) + r") *$")
_CLOSE = re.compile(r"^::: *$")

# Just below pymdownx.snippets (32), so chapters are inlined before this runs.
_PRIORITY = 30


class _CalloutPreprocessor(Preprocessor):
    def run(self, lines: list[str]) -> list[str]:
        out: list[str] = []
        in_callout = False
        for line in lines:
            if not in_callout:
                opened = _OPEN.match(line)
                if opened:
                    out.append(f"!!! {opened.group(1)}")
                    out.append("")
                    in_callout = True
                else:
                    out.append(line)
                continue
            if _CLOSE.match(line):
                out.append("")
                in_callout = False
                continue
            # Blank lines stay blank; an admonition body is indented by four spaces.
            out.append(f"    {line}" if line.strip() else line)
        return out


class ManualCalloutExtension(Extension):
    def extendMarkdown(self, md) -> None:
        md.preprocessors.register(_CalloutPreprocessor(md), "zombi2_manual_callouts", _PRIORITY)


def on_config(config):
    """Register the bridge as a Markdown extension instance (python-markdown accepts those, so the
    extension does not have to be an installed, importable package)."""
    config.markdown_extensions.append(ManualCalloutExtension())
    return config
