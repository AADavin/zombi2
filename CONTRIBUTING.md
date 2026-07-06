# Contributing to ZOMBI2

Thanks for your interest in ZOMBI2! Bug reports, feature ideas and pull requests are all
welcome.

## Development setup

ZOMBI2 needs Python ≥ 3.10. The library itself is pure Python (only `numpy`), but the default
gene-family engine is a compiled Rust extension that you build once with
[maturin](https://www.maturin.rs/).

```bash
git clone https://github.com/AADavin/zombi2.git
cd zombi2

# the compiled gene-family engine FIRST (needed for the default `genomes` model
# and its tests) -- zombi2 depends on zombi2_core, which isn't on an index during
# development, so building it here satisfies that pin locally
pip install ./rust

# the library + dev tools (pytest, scipy)
pip install -e ".[dev]"
```

Species trees, traits, sequences and the flexible-rate genome models run in pure Python, but
the built-in `genomes` engine — and the tests that cover it — need the compiled extension.

## Running the tests

```bash
pytest
```

The suite is deterministic (seeded) and hermetic — no network, no external data. If
`zombi2.rust_available()` is `False`, the Rust-backed tests **skip** rather than fail, so make
sure you built the engine (above) before concluding the suite passed. CI builds it and treats a
missing engine as an error, so the compiled path is always exercised there.

## Documentation

```bash
pip install -e ".[docs]"
mkdocs serve            # live preview at http://127.0.0.1:8000
mkdocs build --strict   # what CI runs
```

The book-style manual lives in `manual/` (Pandoc → XeLaTeX); see `manual/README.md` for how to
build it.

## Submitting changes

1. Branch off `main`.
2. Keep each pull request focused. Add or update tests for any behaviour change, and prefer
   asserting real invariants — seeded determinism, conservation laws, analytic/oracle values —
   over "it runs without error".
3. Match the surrounding style: public modules, classes and functions carry docstrings.
4. Open a pull request. CI must pass — `pytest` on Python 3.10–3.12 with the Rust engine built,
   plus a strict documentation build.

## License

ZOMBI2 is released under the **GNU General Public License v3.0 or later** (see
[`LICENSE`](LICENSE)). By contributing, you agree that your contributions are licensed under the
same terms.
