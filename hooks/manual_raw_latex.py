"""Drop the manual's raw-LaTeX blocks on the way into the docs site.

A chapter may carry a block pandoc passes straight to LaTeX::

    ```{=latex}
    \\appendix
    ```

which is how Appendix A tells the book that the appendices start here. It is typesetting, not
content, and it means nothing to the site: python-markdown reads it as an ordinary fenced code
block and prints ``\\appendix`` to the reader in a grey box.

There is no spelling both toolchains accept, so — as with `manual_callouts` and `manual_figures` —
the single source stays pandoc-native and this drops the block on the way in. Only ``{=FORMAT}``
fences are touched, which is pandoc's raw-attribute syntax and cannot be a language name, so a
normal code block is never at risk.

Like the others, it has to run as a Markdown *preprocessor* below `pymdownx.snippets` (priority
32): before the include expands, a chapter page is still just its ``--8<--`` line.
"""

from __future__ import annotations

import re

from markdown.extensions import Extension
from markdown.preprocessors import Preprocessor

_FENCE = re.compile(r"^(`{3,}|~{3,}) *\{=[A-Za-z]+\} *$")

# Just below pymdownx.snippets (32), so chapters are inlined before this runs.
_PRIORITY = 30


class _RawLatexPreprocessor(Preprocessor):
    def run(self, lines: list[str]) -> list[str]:
        out: list[str] = []
        closing: str | None = None
        for line in lines:
            if closing is not None:                       # inside a raw block: drop until it ends
                if line.strip() == closing:
                    closing = None
                continue
            match = _FENCE.match(line.strip())
            if match:
                closing = match.group(1)
                continue
            out.append(line)
        return out


class _RawLatexExtension(Extension):
    def extendMarkdown(self, md) -> None:
        md.preprocessors.register(_RawLatexPreprocessor(md), "manual_raw_latex", _PRIORITY)


def on_config(config):
    config.markdown_extensions.append(_RawLatexExtension())
    return config
