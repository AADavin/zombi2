# Contributing

## Dev setup

```bash
git clone https://github.com/AADavin/zombi2
cd zombi2
pip install -e ".[dev]"
pytest
```

The test suite lives in `tests/` and runs in a couple of seconds. A `conftest.py` at the
repo root puts the project on `sys.path`, so `import zombi2` works even without the
editable install.

## Package layout

The package mirrors the **four levels of evolution** — each is a subpackage that re-exports
its public API from an `__init__.py` (the scikit-learn-style namespace), plus a `coevolve/`
layer for couplings among them and a shared kernel at the top level:

| Subpackage | Responsibility |
| --- | --- |
| `zombi2/species/` | species-tree models + simulation (`model.py`, `sim.py`, `forward.py`, `ghosts.py`) |
| `zombi2/genomes/` | gene families + genome structure (`rates.py`, `genome.py`, `genome_sim.py`, `simulation.py`, `reconciliation.py`, `profiles.py`, `events.py`, `transfers.py`, `nucleotide_*.py`, `gff.py`) |
| `zombi2/traits/` | trait models (`models.py`) + biogeography (`biogeography.py`, the DEC model) |
| `zombi2/sequences/` | substitution models (`models.py`), the gene×lineage clock (`evolution.py`), and the relaxed molecular clocks (`clocks.py`) |
| `zombi2/coevolve/` | the couplings among the four levels (`sse.py`, `trait_coupling.py`, `gene_diversification.py`, …) |

Shared kernel at the top level:

| Module | Responsibility |
| --- | --- |
| `tree.py` | `Tree` / `TreeNode`, Newick read/write |
| `distributions.py` | distribution helpers (`Gamma`, `Exponential`, `as_distribution`, …) |
| `_sampling.py` | the Gillespie event sampler (`EventSampler`) |
| `_rust.py` | native Rust engine bindings; the built-in model runs here automatically |
| `parallel.py` | `run_replicates` (replicate-level parallelism) |
| `cli.py` | the `zombi2` command-line wrapper |
| `rust/` | the Rust `zombi2_core` extension (built with maturin) |

## Adding a model

ZOMBI2 grows by **adding models, not editing the engine** — most additions are new subclasses.
Before you start, read the contract:

- **[Adding a model](contributing/adding-a-model.md)** — the interface to implement for each
  level, and the end-to-end checklist (implement → export → CLI → validate → document).
- **[Conventions](contributing/conventions.md)** — the names, outputs, seeding, and CLI grammar
  every model follows.
- **[Extending ZOMBI2](guide/extending.md)** — a worked example of the gene-family seams.

**The hard rule:** *no model enters the core without an oracle or a statistical test* — a check
that only asserts "it runs without error" is not validation. See [Validation](validation.md).

## Building the docs

```bash
pip install -e ".[docs]"
mkdocs serve      # live-reload preview at http://127.0.0.1:8000
mkdocs build      # static site into site/
```

Docstrings use the NumPy style and are pulled into the [API reference](reference/api.md) by
`mkdocstrings`.

### The wiki is generated

The GitHub wiki is **not** edited by hand — it is a mirror of the narrative pages in
`docs/`, produced by `tools/sync_wiki.py` and pushed by the `sync-wiki` GitHub Action on
every change to `docs/`. Edit the pages under `docs/`; the wiki updates itself. (The
[API reference](reference/api.md) is MkDocs-only, since it is generated from docstrings.)

## Building the Rust engine

Only needed if you're working on the Rust engine in `rust/` or want to benchmark it:

```bash
pip install maturin
cd rust
maturin build --release -i python3
pip install --force-reinstall target/wheels/*.whl
```

The Rust engine must stay behaviourally consistent with the pure-Python reference (they agree
statistically, not bit-for-bit, since the RNG streams differ). The built-in model runs on
Rust; flexible models run on Python. See [the Rust engine](guide/rust-engine.md).

## Tests and style

- Add or update tests in `tests/` for any behaviour change; keep the suite fast.
- Match the surrounding code's style, naming, and comment density.
- The public API is what `zombi2/__init__.py` exports (`__all__`); update it when adding
  user-facing symbols.
