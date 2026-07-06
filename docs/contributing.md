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
| `rates.py` | rate models (`SharedRates`, `FamilySampledRates`, `PerGenomeRates`, `BranchRates`) |
| `rate_variation.py` | substitution-rate variation across one tree (relaxed clock) |
| `sequence_evolution.py` | gene × lineage substitution clock over reconciled gene trees (`SequenceEvolution`) |
| `distributions.py` | distribution helpers (`Gamma`, `Exponential`, `as_distribution`, …) |
| `genome.py` | genome representations (`UnorderedGenome`, `OrderedGenome`) |
| `transfers.py` | `TransferModel` (additive/replacement, distance decay, self-transfer) |
| `genome_sim.py` | the forward Gillespie gene-family simulator |
| `events.py` | event records and types |
| `reconciliation.py` | reconstructing gene trees from the event log |
| `profiles.py` | `ProfileMatrix` (copy-number / presence) |
| `simulation.py` | `simulate_genomes` + the `Genomes` result and `.write()` |
| `parallel.py` | `run_replicates` (replicate-level parallelism) |
| `_rust.py` | native Rust engine bindings; the built-in model runs here automatically |
| `ghosts.py` | ghost lineages (`add_ghost_lineages`, un-pruning) |
| `matching.py` | ABC profile matching (`match_profiles`, `match_profiles_smc`) |
| `species_forward.py` | forward-in-time / fossilized birth–death trees |
| `nucleotide_genome.py` / `nucleotide_sim.py` | the nucleotide genome model |
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
