# ZOMBI2

A Python library for simulating **species-tree** and **gene-family** evolution in
phylogenetics. A ground-up redesign of [ZOMBI](https://github.com/AADavin/Zombi).

- The **species tree** is simulated *backward* in time (reconstructed birth–death,
  conditioned on the number of extant species).
- **Gene families** are simulated *forward* in time along that fixed species tree, with
  duplication, transfer, loss and origination events.

This is **v1**: an order-free, family-level model with independent gene families. The
architecture is deliberately interface-first so that later additions — gene order and
rearrangements, gene length, genome-wise rates, gene-family coupling (Potts model),
ghost lineages — arrive as new subclasses rather than a rewrite.

## Quick start

Two symmetric steps — simulate the species tree, then the gene families along it.
Models are objects; simple cases have a keyword shorthand.

```python
import zombi2 as z

# 1. species tree (backward). Yule(birth) == BirthDeath(birth, death=0).
tree = z.simulate_species_tree(z.BirthDeath(birth=1.0, death=0.3),
                               n_tips=20, age=5.0, seed=1)

# 2a. gene families where every family shares the same D/T/L (shorthand):
genomes = z.simulate_genomes(tree, duplication=0.2, transfer=0.1, loss=0.25,
                             origination=0.5, initial_size=40, seed=42)

# 2b. ...or every family draws its OWN D/T/L from distributions (ZOMBI-1 style):
genomes = z.simulate_genomes(tree, z.FamilySampledRates(
    duplication=z.Gamma(2, 0.06), transfer=z.Exponential(0.08),
    loss=z.Gamma(2, 0.07), origination=0.5), initial_size=40, seed=42)

print(tree.to_newick())
print(genomes.profiles.matrix)          # families x extant-species copy numbers
complete, extant = genomes.gene_trees()["1"]   # per-family reconstructed gene trees
genomes.write("out/")                   # trees, event tables, transfers, profiles
```

`genomes.write("out/")` produces ZOMBI-1-style output: `species_tree.nwk`,
per-family event tables (`gene_family_events/`), reconstructed **complete** and
**extant** gene trees (`gene_trees/`), `Transfers.tsv`, `Gene_family_summary.tsv`, and
the presence/copy-number matrices (`Profiles.tsv` / `Presence.tsv`).

**Bounding gene-family growth.** With duplication > loss a family's size grows like
`e^{(d−l)t}`. Pass `carrying_capacity=K` (logistic density dependence — family size
settles around `K`) and/or `max_copies=N` (a hard cap) to any rate model:

```python
z.UniformRates(duplication=0.5, loss=0.1, origination=0.3, carrying_capacity=20)
```

Rate models are the flexible knob. `z.UniformRates(...)` and `z.FamilySampledRates(...)`
ship today; new models (genome-wise rates, gene-family coupling / Potts) are subclasses
of `z.RateModel` that need no change to the simulator. Distribution arguments accept a
built-in (`z.Gamma`, `z.Exponential`, `z.LogNormal`, `z.Uniform`, `z.Fixed`), any
`scipy.stats` frozen distribution, or a callable `rng -> float`.

### CLI

```bash
zombi2 species --birth 1 --death 0.3 --tips 20 --age 5 --seed 1 -o out/
zombi2 all --birth 1 --death 0.2 --tips 20 --age 5 \
    --dup 0.2 --trans 0.1 --loss 0.25 --orig 0.5 --seed 42 -o out/
```

## Development

```bash
uv pip install -e ".[dev]"   # or: pip install -e ".[dev]"
pytest
```
