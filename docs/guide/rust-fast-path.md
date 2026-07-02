# The Rust fast path (optional)

ZOMBI2 ships an optional native engine that runs the forward Gillespie in Rust, in three
flavours:

* **`simulate_profiles_fast`** — genomes are per-family *counts* only; returns just the
  presence/copy-number **profile matrix** (the σ dataset). Fastest: the 10 000-tip matrix in
  **~0.4 s vs ~21 s** for the pure-Python sim (≈50×).
* **`simulate_genomes_fast`** — tracks individual gene lineages and emits the **full event
  genealogy**, returning a complete [`Genomes`](gene-families.md) with `.event_log`,
  `.gene_trees()` and `.write()` — a drop-in for `simulate_genomes` (~3× at 10k tips; limited
  by materialising ~2M Python objects).
* **`simulate_and_write_fast`** — simulates, reconstructs the gene trees, **and writes the
  full ZOMBI-1 output to disk entirely in Rust**, returning only a small summary. This is the
  scale path for "simulate-and-write-everything" — nothing large crosses back into Python
  (**~10× vs Python simulate + write** at 10k tips: ~6.5 s vs ~66 s).

The pure-Python [`simulate_genomes`](gene-families.md) remains the default. All three Rust
engines cover exactly the built-in `UnorderedGenome` + `UniformRates` model (per-copy
duplication / transfer / loss, per-branch origination, additive uniform-recipient
transfers, optional hard `max_family_size`) and the default `TransferModel`.

## Building the extension

The extension isn't built by `pip install -e .`; build it once with
[maturin](https://www.maturin.rs/):

```bash
pip install maturin
cd rust && maturin build --release -i python3
pip install --force-reinstall rust/target/wheels/*.whl
```

Check it's available:

```python
import zombi2 as z
z.rust_available()   # True once the wheel is installed
```

## Usage

```python
tree = z.simulate_species_tree(z.BirthDeath(1.0, 0.4), n_tips=10_000, age=10.0, seed=1)

profiles = z.simulate_profiles_fast(
    tree, duplication=0.05, transfer=0.03, loss=0.1, origination=0.5,
    initial_size=200, max_family_size=0.3, seed=42,
)
print(profiles.matrix.shape)     # (n_families, n_extant_species)
profiles.to_tsv()                # same ProfileMatrix as the Python engine
```

It accepts the keyword shorthand or a `UniformRates` object; other rate models
(`GenomeWiseRates`, `FamilySampledRates`, `BranchRates`), soft `carrying_capacity`, and
rearrangements raise a clear error — use `simulate_genomes` for those.

## Full event log + gene trees

`simulate_genomes_fast` returns a full `Genomes`, so everything the Python engine offers
downstream works unchanged:

```python
g = z.simulate_genomes_fast(tree, duplication=0.1, transfer=0.05, loss=0.15,
                            origination=0.5, initial_size=100, max_family_size=0.5, seed=1)

len(g.event_log)                 # every O/D/T/L/S record, as EventRecord objects
complete, extant = g.gene_trees()["1"]   # reconstructed per-family gene trees
g.write("out/")                  # full ZOMBI-1 output (event tables, trees, transfers, ...)
```

The Rust engine generates the genealogy; Python materialises the `EventLog` and runs the
existing reconciliation/writers. Gene ids are integers (not `g`-prefixed) and the RNG differs
from the Python engine, so results are statistically — not bit — identical. The reconciliation
invariant holds: each family's extant gene-tree leaf count equals its extant copy number.

## Writing everything from Rust (scale path)

To generate a large dataset and get files on disk, `simulate_and_write_fast` does the whole
pipeline — simulate, reconstruct gene trees, write every output — in Rust:

```python
summary = z.simulate_and_write_fast(
    tree, "out/", duplication=0.05, transfer=0.03, loss=0.1, origination=0.5,
    initial_size=200, max_family_size=0.3, seed=42,
)
# {'path': 'out/', 'n_families': ..., 'n_events': ..., 'n_species': ...}
```

It writes the identical file set as `Genomes.write` (`gene_family_events/`, `gene_trees/`,
`Transfers.tsv`, `Gene_family_summary.tsv`, `Profiles.tsv`, `Presence.tsv`; Python writes only
the two tiny tree-only files). Because it never builds the multi-million-object Python log, it
is far faster than `simulate_genomes(...).write(...)` at scale.

## Reproducibility & correctness

The Rust engine uses its own PRNG, so results are **statistically** equivalent to the
Python engine, not bit-identical. A given `seed` is reproducible *within* this engine. The
test suite checks that mean copy number and family counts agree with the Python engine
within Monte-Carlo error.

## Performance notes

At 10 000 tips (~2.1M events):

| Path | Time | vs Python |
| --- | --- | --- |
| `simulate_profiles_fast` (profiles only) | ~0.4 s | ~50× |
| `simulate_genomes_fast` (full `Genomes`, materialised) | ~5.5 s | ~3× |
| `simulate_and_write_fast` (simulate + trees + write) | ~6.5 s | ~10× |
| `simulate_genomes` + `.write()` (pure Python) | ~66 s | 1× |

`simulate_genomes_fast` is "only" ~3× because Python still materialises ~2M `EventRecord`
objects; `simulate_and_write_fast` avoids that entirely (everything stays in Rust) and is the
right choice when you want files on disk at scale.

## What's next

Gene-tree reconstruction and the writers are now in Rust. The remaining items are the richer
`TransferModel` mechanics (replacement / distance / self) in the Rust engines, and
within-simulation per-family threading (rayon).
