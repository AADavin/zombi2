# Installation

ZOMBI2 needs Python ≥ 3.10 and depends only on **numpy**.

## From source

```bash
git clone https://github.com/AADavin/zombi2.git
cd zombi2
pip install -e .
```

The editable install also exposes the `zombi2` command-line tool.

## Development extras

To run the test suite and the statistical checks:

```bash
pip install -e ".[dev]"   # adds pytest and scipy
pytest
```

## Optional dependencies

- **scipy** — only needed if you pass a `scipy.stats` frozen distribution to
  `FamilySampledRates`; the built-in distributions (`Gamma`, `Exponential`, …) need
  nothing beyond numpy.
- **matplotlib / pandas** — used only by the demo notebook
  (`examples/zombi2_demo.ipynb`), not by the library.

## Building these docs

```bash
pip install mkdocs-material "mkdocstrings[python]"
mkdocs serve      # live preview at http://127.0.0.1:8000
mkdocs build      # static site in site/
```
