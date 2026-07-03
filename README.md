# ZOMBI2

**A phylogenetic simulator of species trees and gene families.**

ZOMBI2 simulates evolution in two steps: build a **species tree**, then evolve
**gene families** along it under duplication, transfer, loss, and origination (DTL).
It is a ground-up redesign of [ZOMBI](https://github.com/AADavin/Zombi), with a fast
Rust engine, a simple command-line interface, and a composable Python library.

Use it to generate benchmark datasets for phylogenetic methods — gene trees,
reconciliations, and copy-number profiles — with fully reproducible, seeded runs.

---

## Installation

ZOMBI2 needs Python ≥ 3.10 and depends only on **numpy**.

```bash
git clone https://github.com/AADavin/zombi2.git
cd zombi2
pip install -e .
```

This puts a `zombi2` command on your PATH.

The gene-family engine is compiled Rust and is built once, separately:

```bash
pip install maturin
cd rust && maturin build --release -i python3
pip install --force-reinstall rust/target/wheels/*.whl
```

---

## Command line

Everything starts from a species tree; you then evolve gene families along it. Every
option has a sensible default, and every run writes a `.log` recording the version,
seed, and full command line so it can be reproduced exactly.

### 1. Simulate a species tree

```bash
zombi2 species --birth 1 --death 0.3 --tips 50 --age 5 --seed 1 -o out/
```

Writes `out/species_tree.nwk` (timed Newick). By default this is the *reconstructed*
tree — extant lineages only, conditioned on `--tips`. `--tips` sets the number of
extant species and `--age` sets the timescale; they are independent, so give both.

### 2. Evolve gene families along it

```bash
zombi2 genomes --tree out/species_tree.nwk \
    --dup 0.2 --trans 0.1 --loss 0.25 --orig 0.5 --seed 42 -o out/
```

Writes the full ZOMBI-style output to `out/`:

| File / folder | Contents |
|---|---|
| `Profiles.tsv` / `Presence.tsv` | families × species copy-number / presence matrix |
| `gene_trees/` | per-family reconstructed gene trees (complete and extant) |
| `gene_family_events/` | per-family event log (origination, duplication, transfer, loss) |
| `Transfers.tsv` | every transfer (donor, recipient) |
| `Gene_family_summary.tsv` | per-family event counts and extant copies |

Pass `--output profiles` for the counts-only fast path (no gene trees or event log) —
a 10,000-tip profile matrix in under a second.

### Common options

| Option | Command | Meaning |
|---|---|---|
| `--birth` / `--death` | `species` | speciation / extinction rate (`--death 0` = Yule) |
| `--tips` / `--age` | `species` | number of extant species / tree timescale |
| `--tree` / `-t` | `genomes` | input species tree (any Newick file) |
| `--dup` `--trans` `--loss` `--orig` | `genomes` | per-copy DTL and origination rates |
| `--max-family-size` | `genomes` | growth cap — integer, or decimal fraction of N (e.g. `0.5`) |
| `--output` | `genomes` | which files to write (`profiles trees` by default) |
| `--seed` | both | RNG seed for reproducibility |
| `-o` / `--out` | both | output directory |

Run `zombi2 <command> -h` for the full list. Two further commands round out the tool:
`zombi2 trait` evolves a phenotypic trait along a tree (BM, OU, Mk, threshold, DEC),
and `zombi2 abc` fits DTL rates to an empirical profile by Approximate Bayesian
Computation. See [`docs/cli.md`](docs/cli.md) for all four commands.

---

## Python library

The same two steps from Python, where models are first-class objects you can compose:

```python
import zombi2 as z

# 1. species tree (backward). Yule(birth) == BirthDeath(birth, death=0).
tree = z.simulate_species_tree(z.BirthDeath(birth=1.0, death=0.3),
                               n_tips=20, age=5.0, seed=1)

# 2. gene families along it
genomes = z.simulate_genomes(tree, duplication=0.2, transfer=0.1, loss=0.25,
                             origination=0.5, initial_size=40, seed=42)

print(tree.to_newick())
print(genomes.profiles.matrix)          # families × species copy numbers
complete, extant = genomes.gene_trees()["1"]
genomes.write("out/")                   # trees, event tables, transfers, profiles
```

The library exposes the flexibility the CLI leaves out:

- **Rate models** beyond uniform — `z.FamilySampledRates` (each family draws its own
  DTL rates from distributions, ZOMBI-1 style) and `z.GenomeWiseRates`.
- **Species-tree models** — episodic/skyline rates, fossilized birth–death, incomplete
  sampling, and forward (complete) trees with extinct lineages.
- **Transfer mechanics** — additive vs. replacing transfers, distance-decayed recipients
  (`z.TransferModel`).
- **Genome representations** — an ordered-chromosome model with inversions and
  transpositions (`z.OrderedGenome`).

```python
# every family draws its own D/T/L from distributions
genomes = z.simulate_genomes(tree, z.FamilySampledRates(
    duplication=z.Gamma(2, 0.06), transfer=z.Exponential(0.08),
    loss=z.Gamma(2, 0.07), origination=0.5), initial_size=40, seed=42)
```

The architecture is interface-first: new rate models, genome representations, and
species-tree processes arrive as subclasses that the simulator uses unchanged.

---

## Documentation

Full guides and API reference live in `docs/` (build with
`pip install -e ".[docs]" && mkdocs serve`). Start with
[`docs/quickstart.md`](docs/quickstart.md), the
[command-line reference](docs/cli.md), and the user guide under `docs/guide/`.

## Development

```bash
pip install -e ".[dev]"   # adds pytest and scipy
pytest
```
