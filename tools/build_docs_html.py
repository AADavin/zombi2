#!/usr/bin/env python3
"""Build the whole documentation into ONE self-contained, offline HTML file.

Reads the same ``docs/`` pages as the wiki sync (single source of truth), converts each
with Python-Markdown + Pygments, and writes a single ``zombi2-docs.html`` with an inline
sidebar, styles, and syntax highlighting — no server, no assets, just a file you open in a
browser. Re-run it whenever the docs change.

Usage::

    python tools/build_docs_html.py [output.html]      # default: ./zombi2-docs.html
"""

from __future__ import annotations

import html as _html
import os
import re
import sys

import markdown
from pygments.formatters import HtmlFormatter

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from sync_wiki import NAV, REPO, STUB_PAGES  # reuse the page list + order

DOCS = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "docs")


def anchor(src: str) -> str:
    return src[:-3].replace("/", "-").replace("_", "-") if src.endswith(".md") else src


PATH_ANCHOR = {src: anchor(src) for src, _w, _l, _s in NAV}


def rewrite_links(body: str, src: str) -> str:
    """Point intra-doc ``*.md`` links at the corresponding in-page ``#anchor``."""
    src_dir = os.path.dirname(src)

    def repl(m: re.Match) -> str:
        href = m.group(2)
        if href.startswith(("http://", "https://", "#", "mailto:")):
            return m.group(0)
        path = href.split("#", 1)[0]
        if not path.endswith(".md"):
            return m.group(0)
        resolved = os.path.normpath(os.path.join(src_dir, path)).replace(os.sep, "/")
        anc = PATH_ANCHOR.get(resolved)
        return f'{m.group(1)}#{anc}{m.group(3)}' if anc else m.group(0)

    return re.sub(r'(href=")([^"]+)(")', repl, body)


def render_page(src: str) -> str:
    if src in STUB_PAGES:
        return (
            "<h1>API reference</h1><p>The full API reference is generated from the source "
            "docstrings by <code>mkdocstrings</code> and only renders in the MkDocs site "
            f"(<code>mkdocs serve</code>). See the <a href=\"{REPO}\">repository</a>.</p>"
        )
    md = markdown.Markdown(
        extensions=["fenced_code", "tables", "admonition", "toc", "codehilite", "attr_list", "sane_lists"],
        extension_configs={"codehilite": {"guess_lang": False}},
    )
    with open(os.path.join(DOCS, src), encoding="utf-8") as f:
        return rewrite_links(md.convert(f.read()), src)


def build_sidebar() -> str:
    out, section = ['<nav id="side"><div class="brand">ZOMBI2</div>'], "__top__"
    for src, _w, label, sec in NAV:
        if src not in STUB_PAGES and not os.path.exists(os.path.join(DOCS, src)):
            continue
        if sec != section:
            if sec is not None:
                out.append(f'<div class="sec">{_html.escape(sec)}</div>')
            section = sec
        cls = "lnk sub" if sec is not None else "lnk"
        out.append(f'<a class="{cls}" href="#{anchor(src)}">{_html.escape(label)}</a>')
    out.append("</nav>")
    return "".join(out)


CSS = """
*{box-sizing:border-box}
:root{--bg:#fff;--fg:#1b1c1d;--muted:#5c6166;--acc:#0b6b63;--border:#e3e6e8;--card:#f6f8fa;--side:#fafbfc}
body{margin:0;font:16px/1.7 -apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,Helvetica,Arial,sans-serif;color:var(--fg);background:var(--bg)}
#wrap{display:flex;align-items:flex-start;max-width:1180px;margin:0 auto}
#side{position:sticky;top:0;align-self:flex-start;width:250px;height:100vh;overflow:auto;padding:20px 14px;border-right:1px solid var(--border);background:var(--side);font-size:14px}
#side .brand{font-weight:600;font-size:16px;margin:0 8px 12px}
#side .sec{margin:14px 8px 4px;font-size:12px;text-transform:uppercase;letter-spacing:.04em;color:var(--muted)}
#side a.lnk{display:block;padding:5px 8px;border-radius:6px;color:var(--fg);text-decoration:none}
#side a.lnk.sub{padding-left:18px;color:var(--muted)}
#side a.lnk:hover{background:var(--card);color:var(--acc)}
#side a.lnk.active{background:#e7f4f2;color:var(--acc);font-weight:500}
main{flex:1;min-width:0;padding:32px 40px 120px;max-width:860px}
main h1{font-size:30px;margin:.2em 0 .5em;padding-top:12px}
main h2{font-size:23px;margin:1.6em 0 .5em;padding-bottom:.2em;border-bottom:1px solid var(--border)}
main h3{font-size:18px;margin:1.3em 0 .4em}
main p,main li{color:var(--fg)}
main a{color:var(--acc);text-decoration:none}
main a:hover{text-decoration:underline}
code{font:13.5px/1.5 "SF Mono",SFMono-Regular,Menlo,Consolas,monospace;background:var(--card);padding:1.5px 5px;border-radius:5px}
pre{background:var(--card);border:1px solid var(--border);border-radius:10px;padding:14px 16px;overflow:auto}
pre code{background:none;padding:0}
table{border-collapse:collapse;margin:1em 0;display:block;overflow:auto}
th,td{border:1px solid var(--border);padding:7px 12px;text-align:left}
th{background:var(--card);font-weight:600}
blockquote,.admonition{margin:1em 0;padding:.4em 1em;border-left:4px solid var(--acc);background:var(--card);border-radius:0 8px 8px 0}
.admonition-title{font-weight:600;margin:.2em 0}
.doc-page{scroll-margin-top:16px;border-bottom:1px dashed var(--border);margin-bottom:24px}
.topnote{color:var(--muted);font-size:13px;margin:0 0 18px}
@media(max-width:820px){#side{display:none}main{padding:20px}}
"""


def main(argv: list[str]) -> int:
    out_path = argv[1] if len(argv) > 1 else "zombi2-docs.html"
    pyg = HtmlFormatter(style="default").get_style_defs(".codehilite")

    sections = []
    for src, _w, _l, _s in NAV:
        if src not in STUB_PAGES and not os.path.exists(os.path.join(DOCS, src)):
            continue
        sections.append(f'<section id="{anchor(src)}" class="doc-page">{render_page(src)}</section>')

    doc = (
        "<!doctype html><html lang=en><head><meta charset=utf-8>"
        "<meta name=viewport content='width=device-width,initial-scale=1'>"
        "<title>ZOMBI2 documentation</title>"
        f"<style>{CSS}\n{pyg}</style></head><body><div id=wrap>"
        f"{build_sidebar()}<main>"
        "<p class=topnote>Offline build of the ZOMBI2 docs — regenerate with "
        "<code>python tools/build_docs_html.py</code>.</p>"
        f"{''.join(sections)}</main></div>"
        "<script>"
        "var links=[].slice.call(document.querySelectorAll('#side a.lnk'));"
        "var secs=links.map(function(a){return document.getElementById(a.getAttribute('href').slice(1));});"
        "function onScroll(){var y=window.scrollY+90,i=0;for(var k=0;k<secs.length;k++){if(secs[k]&&secs[k].offsetTop<=y)i=k;}"
        "links.forEach(function(a,j){a.classList.toggle('active',j===i);});}"
        "document.addEventListener('scroll',onScroll);onScroll();"
        "</script></body></html>"
    )
    with open(out_path, "w", encoding="utf-8") as f:
        f.write(doc)
    print(f"wrote {out_path}  ({len(sections)} pages, {len(doc)//1024} KB)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
