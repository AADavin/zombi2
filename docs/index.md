# ZOMBI2

**ZOMBI2** is a Python library for simulating **species-tree** and **gene-family**
evolution in phylogenetics — a ground-up redesign of
[ZOMBI](https://github.com/AADavin/Zombi).

The workflow is two symmetric steps:

1. Simulate a **species tree** *backward* in time — a reconstructed birth–death process
   conditioned on the number of extant species.
2. Simulate **gene families** *forward* in time along that fixed tree — duplication,
   transfer, loss and origination, with optional gene order and rearrangements.

```python
import zombi2 as z

tree = z.simulate_species_tree(z.BirthDeath(birth=1.0, death=0.3), n_tips=20, age=5.0, seed=1)
genomes = z.simulate_genomes(tree, duplication=0.2, transfer=0.1, loss=0.25,
                             origination=0.5, seed=42)

print(genomes.profiles.matrix)          # gene families × extant species (copy numbers)
complete, extant = genomes.gene_trees()["1"]
genomes.write("out/")                   # trees, event tables, transfers, profiles
```

## What's in the box

- **Species trees** — constant-rate birth–death and Yule, conditioned on the number of
  extant tips; a [roadmap](species_tree_models.md) of further models (episodic/skyline,
  diversity-dependent, fossilized birth–death).
- **Gene families** — duplication / transfer / loss / origination along the tree, with
  two rate models: the same rates for every family (`UniformRates`) or **per-family rates
  drawn from distributions** (`FamilySampledRates`, ZOMBI-1 style).
- **Transfers** — additive or **replacement**, **phylogenetic-distance-weighted**
  recipient choice, and optional **self-transfer** (= duplication).
- **Gene order** — an ordered-chromosome genome (`OrderedGenome`) with segment-based
  events plus **inversions** and **transpositions**.
- **Growth control** — a hard family-size cap (`max_family_size`, absolute or a fraction
  of the number of species) and a soft logistic `carrying_capacity`.
- **Outputs** — reconstructed **complete and pruned gene trees**, per-family event
  tables, a transfers table, a gene-family summary, and the **presence/absence profile
  matrix** (the dataset for phylogenetic-profiling analyses).

## Design philosophy

ZOMBI2 is **interface-first**: one Gillespie simulator programs only against `Genome`,
`RateModel` and `EventSampler` protocols, so new genome representations, rate models and
event types drop in as subclasses without touching the engine. See
[Extending ZOMBI2](guide/extending.md).

## Where next

- New to ZOMBI2? Start with [Installation](installation.md) and the
  [Quickstart](quickstart.md).
- Then work through the **User guide** in the navigation.
- The headline scientific goal — the **non-independence of gene families** (a Potts-model
  extension) — is written up in [this background note](background.md).
