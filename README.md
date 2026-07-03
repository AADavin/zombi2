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

ZOMBI2 needs Python ≥ 3.10.

```bash
git clone https://github.com/AADavin/zombi2.git
cd zombi2
pip install -e . maturin
cd rust && maturin build --release -i python3 && pip install --force-reinstall target/wheels/*.whl
```

This installs the `zombi2` command and the compiled gene-family engine.

---

## Command line

### 1. Simulate a species tree

```bash
zombi2 species --birth 1 --death 0.3 --tips 50 --age 5 --seed 1 -o out/
```

### 2. Evolve gene families along it

```bash
zombi2 genomes --tree out/species_tree.nwk \
    --dup 0.2 --trans 0.1 --loss 0.25 --orig 0.5 --seed 42 -o out/
```

This writes the full output to `out/`:

| File / folder | Contents |
|---|---|
| `Profiles.tsv` / `Presence.tsv` | families × species copy-number / presence matrix |
| `gene_trees/` | per-family reconstructed gene trees |
| `gene_family_events/` | per-family event log (origination, duplication, transfer, loss) |
| `Transfers.tsv` | every transfer (donor, recipient) |
| `Gene_family_summary.tsv` | per-family event counts and extant copies |

Two further commands round out the tool: `zombi2 trait` evolves a phenotypic trait
along a tree, and `zombi2 abc` fits DTL rates to an empirical profile by Approximate
Bayesian Computation. Run `zombi2 <command> -h` for options, or see
[`docs/cli.md`](docs/cli.md).

---

## Performance

The core models run on a native Rust engine and scale to millions of tips on a laptop.
A backward species tree of 1M tips builds in ~19 s; gene families over a 100k-tip tree
in ~19 s; and the counts-only profile path handles 100k tips in ~1.4 s.

![ZOMBI2 performance](performance_analysis/figures/overview.png)

---

## Models

ZOMBI2 ships a broad range of models, all reachable from the Python API.

**Species-tree models**

- Backward (reconstructed) and forward (complete) birth–death
- Episodic / skyline rate shifts
- Fossilized birth–death and incomplete sampling
- Ghost lineages

**Genome models**

- Uniform DTL rates (default, Rust engine)
- Family-sampled rates — each family draws its own DTL from distributions (ZOMBI-1 style)
- Genome-wise rates
- Ordered chromosomes with inversions and transpositions
- Nucleotide-resolution genomes, where genes emerge as *atoms* from structural events
- Gene-family coupling (a Potts model of non-independence)

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
genomes.write("out/")                   # trees, event tables, transfers, profiles
```

---

## Documentation

Full guides and API reference live in `docs/` (build with
`pip install -e ".[docs]" && mkdocs serve`). Start with
[`docs/quickstart.md`](docs/quickstart.md) and the
[command-line reference](docs/cli.md).

## Development

```bash
pip install -e ".[dev]"   # adds pytest and scipy
pytest
```
