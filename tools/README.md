# tools/

Developer tooling for the documentation. Neither script is part of the `zombi2` package
or needed to use it.

## `build_docs_html.py`

Bundles the whole MkDocs site into a single, self-contained **offline HTML file** (styles
and images inlined), for sharing docs without hosting them:

```bash
python tools/build_docs_html.py [out.html]   # default: ./zombi2-docs.html
```

The output (`zombi2-docs.html`) is a **build artifact and is not tracked in git** (see
`.gitignore`) — regenerate it on demand, or publish the MkDocs site instead
(`mkdocs build`).

## `sync_wiki.py`

Mirrors the canonical `docs/` pages into a GitHub **wiki** checkout (flattening the nav,
rewriting intra-doc links, converting MkDocs admonitions to plain Markdown). Run it
manually against a local wiki clone.

It is **currently dormant**: there is no CI Action wired up to run it, and the MkDocs site
already renders the same `docs/`. The canonical, published documentation is the MkDocs
site — treat the wiki as an optional mirror.
