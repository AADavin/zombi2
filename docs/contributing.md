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

The library is organised so each subsystem is a small, focused module:

| Module | Responsibility |
| --- | --- |
| `species_model.py` | species-tree models (`BirthDeath`, `Yule`, `EpisodicBirthDeath`) |
| `species_sim.py` | backward, conditioned species-tree simulation |
| `tree.py` | `Tree` / `TreeNode`, Newick read/write |
| `rates.py` | rate models (`UniformRates`, `FamilySampledRates`, `GenomeWiseRates`, `BranchRates`) |
| `rate_variation.py` | rate variation across the tree |
| `distributions.py` | distribution helpers (`Gamma`, `Exponential`, `as_distribution`, …) |
| `genome.py` | genome representations (`UnorderedGenome`, `OrderedGenome`) |
| `transfers.py` | `TransferModel` (additive/replacement, distance decay, self-transfer) |
| `genome_sim.py` | the forward Gillespie gene-family simulator |
| `events.py` | event records and types |
| `reconciliation.py` | reconstructing gene trees from the event log |
| `profiles.py` | `ProfileMatrix` (copy-number / presence) |
| `simulation.py` | `simulate_genomes` + the `Genomes` result and `.write()` |
| `parallel.py` | `run_replicates` (replicate-level parallelism) |
| `fast.py` | optional Rust fast paths (`simulate_profiles_fast`, `simulate_genomes_fast`, `simulate_and_write_fast`) |
| `cli.py` | the `zombi2` command-line wrapper |
| `rust/` | the Rust `zombi2_core` extension (built with maturin) |

## Adding a model

Because the simulator only talks to interfaces, most additions are new subclasses that need
no engine changes:

- **A new species-tree model** — subclass the species model in `species_model.py`.
- **A new rate model** — subclass `z.RateModel` (see `rates.py`). Implement how per-family
  or per-copy rates are produced; the forward simulator will drive it unchanged.
- **A new genome representation** — subclass the genome interface in `genome.py`.

See [extending ZOMBI2](guide/extending.md) for a worked example.

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

Only needed if you're working on the fast paths in `rust/` or want to benchmark them:

```bash
pip install maturin
cd rust
maturin build --release -i python3
pip install --force-reinstall target/wheels/*.whl
```

The Rust paths must stay behaviourally consistent with the pure-Python engine (they agree
statistically, not bit-for-bit, since the RNG streams differ). Keep pure Python the
default; Rust is an optional accelerator. See [the Rust fast path](guide/rust-fast-path.md).

## Tests and style

- Add or update tests in `tests/` for any behaviour change; keep the suite fast.
- Match the surrounding code's style, naming, and comment density.
- The public API is what `zombi2/__init__.py` exports (`__all__`); update it when adding
  user-facing symbols.
