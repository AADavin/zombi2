# Installation

ZOMBI2 needs Python 3.10 or newer. Its only runtime dependency is **numpy**; the built-in
gene-family engine additionally requires a compiled Rust extension, built at install time
(see below).

## From source

Clone the repository and install the package in editable mode:

```bash
git clone https://github.com/AADavin/zombi2.git
cd zombi2
pip install -e . maturin
```

The editable install exposes the `zombi2` command-line tool and the Python package.

## Building the compiled engine

ZOMBI2's gene-family simulator is backed by a Rust engine, compiled to a native Python
extension with maturin. From the repository root:

```bash
cd rust && maturin build --release -i python3 && pip install --force-reinstall target/wheels/*.whl
```

This compiles the engine in release mode and installs the resulting wheel over the editable
package, wiring the compiled backend into `zombi2`.

## Development extras

To run the test suite and the statistical checks, install the optional `dev` dependencies:

```bash
pip install -e ".[dev]"   # adds pytest and scipy
pytest
```

## Optional dependencies

- **scipy** — needed only when passing a `scipy.stats` frozen distribution to
  `FamilySampledRates`; the built-in distributions (`Gamma`, `Exponential`, and so on) require
  nothing beyond numpy.
- **matplotlib / pandas** — used only by the demo notebook, not by the library itself.

## Verifying the installation

Confirm that the command-line tool is on your path:

```bash
zombi2 --help
```

This lists the available subcommands (`species`, `genomes`, `trait`, and others). Then confirm
the Python package imports cleanly:

```python
import zombi2
```

If both succeed, the core installation is complete. To exercise the compiled gene-family engine
end to end, simulate a small species tree and evolve genomes along it:

```bash
zombi2 species --birth 1 --death 0.3 --tips 50 --age 5 --seed 1 -o out/
zombi2 genomes --tree out/species_tree.nwk \
    --dup 0.2 --trans 0.1 --loss 0.25 --orig 0.5 --seed 42 -o out/
```

A successful `genomes` run confirms that the compiled backend is correctly installed.
