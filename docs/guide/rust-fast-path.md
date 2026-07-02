# The Rust fast path (optional)

For large-scale studies where you only need the **presence/copy-number profile matrix**
(the σ dataset) — not the full event log or gene trees — ZOMBI2 ships an optional native
engine that runs the forward Gillespie in Rust over per-family *counts*. On a 10 000-tip
tree it produces the profile matrix in **~0.4 s vs ~21 s** for the pure-Python full
simulation (≈50× faster).

The pure-Python [`simulate_genomes`](gene-families.md) remains the default and the **only**
path that produces the event log, transfers, and reconstructed gene trees. The Rust engine
covers exactly the built-in `UnorderedGenome` + `UniformRates` model (per-copy
duplication / transfer / loss, per-branch origination, additive uniform-recipient
transfers, optional hard `max_family_size`).

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

## Reproducibility & correctness

The Rust engine uses its own PRNG, so results are **statistically** equivalent to the
Python engine, not bit-identical. A given `seed` is reproducible *within* this engine. The
test suite checks that mean copy number and family counts agree with the Python engine
within Monte-Carlo error.

## What's next

This is the first increment of the Rust core. It's profiles-only by design; porting the
full event log / gene-tree genealogy and the richer transfer mechanics
(`TransferModel`) to Rust — which also opens the door to within-simulation per-family
threading — is planned follow-up work.
