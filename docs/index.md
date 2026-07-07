# ZOMBI2

**ZOMBI2** simulates phylogenetic evolution end to end — **species trees**, then
**gene families**, **phenotypic traits** and molecular **sequences** along them, plus
**coupled (coevolving)** processes. It is a ground-up redesign of
[ZOMBI](https://github.com/AADavin/Zombi), with a fast Rust engine, a composable Python
library and a command-line interface.

The workflow is two symmetric steps:

1. Simulate a **species tree** *backward* in time — a reconstructed birth–death process
   conditioned on the number of extant species.
2. Simulate **gene families** *forward* in time along that fixed tree — duplication,
   transfer, loss and origination, with optional gene order and rearrangements.

```python
from zombi2.species import BirthDeath, simulate_species_tree
from zombi2.genomes import simulate_genomes

tree = simulate_species_tree(BirthDeath(birth=1.0, death=0.3), n_tips=20, age=5.0, seed=1)
genomes = simulate_genomes(tree, duplication=0.2, transfer=0.1, loss=0.25,
                           origination=0.5, seed=42)

print(genomes.profiles.matrix)          # gene families × extant species (copy numbers)
complete, extant = genomes.gene_trees()["1"]
genomes.write("out/")                   # trees, event tables, transfers, profiles
```

## What's in the box

- **Species trees** — backward (reconstructed) and forward (complete) birth–death and
  Yule, episodic/skyline rate shifts, fossilized birth–death with incomplete sampling,
  heterogeneous-rate diversification (ClaDS, diversity-dependent, clade shifts), mass
  extinctions and ghost lineages. The Rust engine scales to millions of tips.
- **Gene families** — duplication / transfer / loss / origination along the tree, with
  uniform (`SharedRates`), per-family-sampled (`FamilySampledRates`, ZOMBI-1 style) and
  genome-wise rate models; transfers additive or **replacement** with distance-weighted
  recipients; ordered chromosomes with **inversions** and **transpositions**;
  nucleotide-resolution genomes; and a **Potts model of gene-family coupling**
  (non-independence). Output as full event logs, compact event traces, or counts-only
  sparse **profile matrices**.
- **Traits** — Brownian motion, Ornstein–Uhlenbeck, early burst, Mk, threshold and
  related models, plus **DEC biogeography**, evolved along a phylogeny.
- **Sequences** — a gene × lineage **relaxed-clock** family that rescales gene trees from
  time into substitutions/site, plus nucleotide substitution models (JC / K80 / HKY / GTR
  + Γ).
- **Coevolution** — couple species, traits and genes along six directed edges with
  `coevolve --couple driver:target`.
- **Growth control** — a hard family-size cap (`max_family_size`, absolute or a fraction
  of the number of species) and a soft logistic `carrying_capacity`.

## Design philosophy

ZOMBI2 is **interface-first**: one Gillespie simulator programs only against `Genome`,
`RateModel` and `EventSampler` protocols, so new genome representations, rate models and
event types drop in as subclasses without touching the engine. See
[Extending ZOMBI2](guide/extending.md).

## Where next

- New to ZOMBI2? Start with [Installation](installation.md) and the
  [Quickstart](quickstart.md).
- Then work through the **User guide** in the navigation.
- Heterogeneous-rate diversification, mass extinctions and the **Potts model of
  gene-family coupling** (the non-independence of gene families) are all implemented — the
  [command-line reference](cli.md) and the user guide cover them.
