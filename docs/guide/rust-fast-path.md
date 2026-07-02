# The Rust fast path (optional)

ZOMBI2 ships an optional native engine that runs the forward Gillespie in Rust, in two
flavours:

* **`simulate_profiles_fast`** — genomes are per-family *counts* only; returns just the
  presence/copy-number **profile matrix** (the σ dataset). Fastest: the 10 000-tip matrix in
  **~0.4 s vs ~21 s** for the pure-Python sim (≈50×).
* **`simulate_genomes_fast`** — tracks individual gene lineages and emits the **full event
  genealogy**, returning a complete [`Genomes`](gene-families.md) with `.event_log`,
  `.gene_trees()` and `.write()` — a drop-in for `simulate_genomes` at large scale.

The pure-Python [`simulate_genomes`](gene-families.md) remains the default. Both Rust engines
cover exactly the built-in `UnorderedGenome` + `UniformRates` model (per-copy
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

## Reproducibility & correctness

The Rust engine uses its own PRNG, so results are **statistically** equivalent to the
Python engine, not bit-identical. A given `seed` is reproducible *within* this engine. The
test suite checks that mean copy number and family counts agree with the Python engine
within Monte-Carlo error.

## Performance notes

At 10 000 tips (~2.1M events): `simulate_genomes_fast` builds the full `Genomes` in **~5.5 s
vs ~17 s** for `simulate_genomes` (≈3×). The Rust side generates the genealogy; the ~3× (not
50×) is because Python still materialises ~2M `EventRecord` objects. And the *downstream*
`gene_trees()` (~18 s) and `write()` (~33 s) are still pure Python — so for the
"simulate-and-write-everything at scale" workflow, those now dominate.

## What's next

The genealogy is in Rust; the remaining large wins are to move **gene-tree reconstruction and
the output writers into Rust** (so nothing large crosses back into Python), and to add the
richer `TransferModel` mechanics (replacement / distance / self) — which also opens the door
to within-simulation per-family threading.
