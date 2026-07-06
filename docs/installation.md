# Installation

ZOMBI2 needs Python ≥ 3.10. The library depends only on **numpy**, but the built-in
(default) gene-family engine is a compiled **Rust** extension that you build once with
[maturin](https://www.maturin.rs/).

## From source

```bash
git clone https://github.com/AADavin/zombi2.git
cd zombi2
pip install -e . maturin

# build the compiled gene-family engine (required for the default `genomes` model)
cd rust && maturin build --release -i python3
pip install --force-reinstall target/wheels/*.whl
```

The editable install exposes the `zombi2` command-line tool. Species trees, traits,
sequences, and the flexible-rate genome models run in pure Python — but the **default
`genomes` engine needs the compiled extension**, and without it `zombi2 genomes` exits
with a build hint that points back to this step.

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
