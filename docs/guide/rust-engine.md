# The Rust engine

The built-in gene-family model runs on a native **Rust** engine (`zombi2_core`). There is no
separate "fast" function and no engine switch: [`simulate_genomes`](genomes.md) routes
the built-in model to Rust automatically, and flexible models to Python. The engine is chosen
by the model, so a given `seed` is always reproducible against a single engine.

* **Built-in model → Rust (required).** The default `UnorderedGenome` with a plain
  `PerCopyRates` (per-copy duplication / transfer / loss, per-branch origination, optional
  hard `max_family_size`) — plus the **full `TransferModel`** (replacement, distance-weighted
  recipients, self-transfers) — runs on Rust. If the extension isn't built, `simulate_genomes`
  raises a clear error telling you to build it (there is no silent Python fallback for the
  built-in model, which is what keeps results reproducible against one engine).
* **Flexible models → Python.** `FamilySampledRates`, `PerGenomeRates`, `BranchRates`, soft
  `carrying_capacity`, ordered genomes, rearrangements, or a custom `sampler` run on the
  pure-Python engine automatically.

## Getting the extension

`pip install zombi2` installs the engine automatically: `zombi2_core` ships as
prebuilt binary wheels (Linux, macOS, Windows; CPython 3.10+), so no Rust
toolchain is needed.

From a source checkout it is a separate package you install from `rust/` — before
an editable install of `zombi2`, since `zombi2` depends on it:

```bash
pip install ./rust      # needs a Rust toolchain (rustup); wraps maturin
pip install -e .
```

Check it's available:

```python
from zombi2 import rust_available
rust_available()   # True once zombi2_core is installed
```

## Two outputs, one function

`simulate_genomes` returns either a full genealogy or just the profile matrix, selected by
`output`:

```python
from zombi2.species import BirthDeath, simulate_species_tree
from zombi2.genomes import simulate_genomes

tree = simulate_species_tree(BirthDeath(1.0, 0.4), n_tips=10_000, age=10.0, seed=1)

# output="genomes" (default): the full Genomes — event log, gene trees, write()
g = simulate_genomes(tree, duplication=0.1, transfer=0.05, loss=0.15,
                     origination=0.5, initial_families=100, max_family_size=0.5, seed=1)
len(g.event_log)                        # every O/D/T/L/S record, as EventRecord objects
complete, extant = g.gene_trees()["1"]  # reconstructed per-family gene trees
g.write("out/")                         # full ZOMBI1 output (tables, trees, transfers, ...)

# output="profiles": just the copy-number matrix — the fast counts-only path (the σ dataset)
pm = simulate_genomes(tree, duplication=0.05, transfer=0.03, loss=0.1, origination=0.5,
                      initial_families=200, max_family_size=0.3, seed=42, output="profiles")
pm.matrix.shape                         # (n_families, n_extant_species)
```

`output="profiles"` skips gene-ids, the event log and gene-tree reconstruction, so it is much
faster and is the right path for generating large presence/absence datasets.
`output="genomes"` tracks the full genealogy.

## Reproducibility & correctness

The Rust engine uses its own PRNG, so results are **statistically** equivalent to the
Python reference engine, not bit-identical. A given `seed` is reproducible *within* Rust. Gene
ids keep the `g`-prefixed form (`g0`, `g1`, …), so the output schema matches the Python engine.
The reconciliation invariant holds: each family's extant gene-tree leaf count equals its extant
copy number. The test suite checks that mean copy number and family counts agree with the
Python engine within Monte-Carlo error.

## Capability boundaries (handled automatically)

* **Sampled-ancestor trees.** The Rust engine assumes a strictly **binary** species tree.
  Trees with degree-two nodes — FBD *sampled ancestors*, from forward simulation with
  `removal < 1` — are routed to the Python engine (which passes genomes through such nodes).
  You do not need to do anything; the correct engine is selected per tree.
* **The nucleotide model.** [`simulate_nucleotide_genomes`](genomes.md) runs on
  Python by default (it emits the full event log and per-block gene trees). Pass
  `output="profiles"` for the Rust path over leaf segments — much faster, enough for
  `profile_matrix()` / `leaf_mosaic()` / `trace_back()`, but with no event log.

## Performance notes

At 10 000 tips (~2.1M events), relative to the pure-Python reference engine:

| Call | Time | vs Python |
| --- | --- | --- |
| `simulate_genomes(..., output="profiles")` | ~0.4 s | ~50× |
| `simulate_genomes(...)` (full `Genomes`, materialised) | ~5.5 s | ~3× |

The full-`Genomes` path is "only" ~3× because Python still materialises ~2M `EventRecord`
objects from the Rust genealogy; the counts-only path avoids that entirely. The forward
Gillespie is sequential; the remaining large win is a per-family decomposition of the
*simulation* (families are independent given the tree).
