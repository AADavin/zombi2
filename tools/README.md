# tools/

Developer tooling for the documentation. It is not part of the `zombi2` package or needed
to use it.

## `build_docs_html.py`

Bundles the whole MkDocs site into a single, self-contained **offline HTML file** (styles
and images inlined), for sharing docs without hosting them:

```bash
python tools/build_docs_html.py [out.html]   # default: ./zombi2-docs.html
```

The page list and its order come from `mkdocs.yml`'s `nav`, so the bundle and the published
site always contain the same pages. A page the nav names but that is missing from `docs/` is
an error, not a skip.

The output (`zombi2-docs.html`) is a **build artifact and is not tracked in git** (see
`.gitignore`) — regenerate it on demand, or publish the MkDocs site instead
(`mkdocs build`).
