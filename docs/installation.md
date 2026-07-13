# Installation

ZOMBI2 needs Python ≥ 3.10.

## With pip

```bash
pip install zombi2
```

This installs the pure-Python `zombi2` package together with its compiled engine
`zombi2_core`, for which prebuilt binary wheels are published for Linux, macOS and
Windows (CPython 3.10 and newer). No Rust toolchain is needed. If no wheel matches
your platform, pip builds `zombi2_core` from source, which requires a
[Rust toolchain](https://rustup.rs).

## From source (development)

The compiled engine (`zombi2_core`, in `rust/`) is a **separate package that
`zombi2` depends on**. Build and install it *first*, then install `zombi2`
editable — otherwise pip tries to fetch `zombi2_core` from an index before it is
published:

```bash
git clone https://github.com/AADavin/zombi2.git
cd zombi2

pip install ./rust          # compiled gene-family engine (needs a Rust toolchain)
pip install -e ".[dev]"     # the library + dev tools (pytest, scipy)
pytest
```

The editable install exposes the `zombi2` command-line tool. Species trees, traits,
sequences, and the flexible-rate genome models run in pure Python — but the **default
`genomes` engine needs the compiled extension**, and without it `zombi2 genomes` exits
with a build hint that points back to this step.

## Optional dependencies

- **scipy** — only needed if you pass a `scipy.stats` frozen distribution to
  `FamilySampledRates`; the built-in distributions (`Gamma`, `Exponential`, …) need
  nothing beyond numpy.
- **`[selection]`** — the experimental protein-language-model selection models
  (`zombi2.experimental.selection`) need PyTorch and ESM-2. Install them with
  `pip install 'zombi2[selection]'`, which pulls in `torch`, `fair-esm` and `scipy`. Heavy
  and opt-in — the default install stays torch-free.

## Building these docs

```bash
pip install mkdocs-material "mkdocstrings[python]"
mkdocs serve      # live preview at http://127.0.0.1:8000
mkdocs build      # static site in site/
```
