# ZOMBI2

[![CI](https://github.com/AADavin/zombi2/actions/workflows/ci.yml/badge.svg)](https://github.com/AADavin/zombi2/actions/workflows/ci.yml)
[![License: GPL-3.0-or-later](https://img.shields.io/badge/license-GPL--3.0--or--later-blue.svg)](LICENSE)
![Python 3.10+](https://img.shields.io/badge/python-3.10%2B-blue)

**A simulator suite for genome evolution.**

ZOMBI2 simulates how genomes evolve along a phylogeny — **species trees**, the **gene
families** that live on them, **genome structure**, **phenotypic traits**, and molecular
**sequences**, plus their **coevolution** — as one composable, seeded, fully reproducible
suite. Use it to generate benchmark datasets with known ground truth for phylogenetic and
comparative methods.

---

## Install

```bash
pip install zombi2
```

Prebuilt wheels are published for Linux, macOS, and Windows (CPython 3.10+), including the
native engine — no toolchain required. Building from source is covered in the
[installation guide](docs/installation.md).

---

## Quickstart

Each level is its own subcommand — run whichever you need. Here a dated species tree, then
gene families evolving along it under duplication, transfer, loss, and origination:

```bash
zombi2 species --birth 1 --death 0.3 --tips 50 --age 5 --seed 1  -o run/
zombi2 genomes --tree run/species_tree.nwk \
    --dup 0.2 --trans 0.1 --loss 0.25 --orig 0.5 --seed 42       -o run/
```

`zombi2 <command> -h` documents each of `species`, `genomes`, `trait`, `sequence`, and
`coevolve`; see the [quickstart](docs/quickstart.md) and [CLI reference](docs/cli.md).

From Python, every model is a first-class object you can compose:

```python
import zombi2 as z

tree = z.simulate_species_tree(z.BirthDeath(birth=1.0, death=0.3), n_tips=20, age=5.0, seed=1)
genomes = z.simulate_genomes(tree, duplication=0.2, transfer=0.1, loss=0.25,
                             origination=0.5, initial_families=40, seed=42)

genomes.write("run/")           # gene trees, event tables, transfers, copy-number profiles
```

---

## Models

ZOMBI2 is a library of models organized by the level of evolution they act on. Each links to
its guide:

- **[Species trees](docs/guide/species-trees.md)** — birth–death (backward and forward),
  episodic/skyline shifts, fossilized birth–death, incomplete sampling, diversity-dependent
  and per-lineage (ClaDS) diversification, mass extinctions, and ghost lineages
  ([full list](docs/species_tree_models.md)).
- **[Gene families](docs/guide/gene-families.md)** — duplication, transfer, loss, and
  origination (DTL), with shared, family-sampled, genome-wise, or per-branch rates and a full
  [transfer model](docs/guide/transfers.md).
- **[Genome structure](docs/guide/ordered-genomes.md)** — ordered chromosomes with inversions
  and transpositions, and [nucleotide-resolution genomes](docs/guide/nucleotide-genomes.md)
  where genes emerge as blocks from structural events.
- **[Traits](docs/guide/traits.md)** — Brownian motion, Ornstein–Uhlenbeck, and early burst
  (continuous); Mk and threshold (discrete); DEC biogeography.
- **[Sequences](docs/cli.md#sequence--substitution-branch-lengths)** — substitution models
  (JC/K80/HKY/GTR + Gamma, and empirical amino-acid models) simulated along the gene trees.
- **[Molecular clocks](docs/guide/rate-variation.md)** — a relaxed-clock family (strict,
  uncorrelated lognormal/gamma, white-noise, autocorrelated, Cox–Ingersoll–Ross) that
  rescales timetrees into substitutions per site.
- **[Coevolution](docs/coevolution_models.md)** — couple any two levels with
  `coevolve --couple driver:target`: state-dependent diversification (SSE), cladogenetic
  change, key innovations, and [trait-linked gene families](docs/guide/trait-linked-genomes.md).

---

## Performance

The built-in models run on a native **Rust** engine and scale to millions of tips on a
laptop: a backward species tree of 1M tips builds in ~6 s, and gene families over a 100k-tip
tree simulate in ~1 s as copy-number profiles. On the same task, ZOMBI2 runs gene-family
simulations **over 1000× faster than ZOMBI 1** (48 s → 38 ms at 1,000 tips). Details in the
[Rust engine guide](docs/guide/rust-engine.md).

![ZOMBI2 performance](analyses/performance/figures/overview.png)

---

## Documentation

Guides, a command-line reference, and the full API live in [`docs/`](docs/) (build locally
with `pip install -e ".[docs]" && mkdocs serve`); a book-style [manual](manual/) is built with
Pandoc. Start with the [quickstart](docs/quickstart.md).

## Citation

If you use ZOMBI2, please cite it via [`CITATION.cff`](CITATION.cff) (GitHub's *Cite this
repository* button). A dedicated ZOMBI2 paper is in preparation; until then, cite the original
[ZOMBI](https://github.com/AADavin/Zombi).

## License

ZOMBI2 is released under the [GNU General Public License v3.0 or later](LICENSE).
