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

```python
import zombi2 as z

species = z.SpeciesTreeModel(birth=1.0, death=0.3, n_tips=20, age=5.0)
rates = z.RateModel(z.EventRates(duplication=0.2, transfer=0.1, loss=0.25, origination=0.5))

result = z.Simulation(species, rates, seed=42).run()
print(result.species_tree.to_newick())
print(result.profiles.matrix)          # families x extant-species copy numbers
result.write("out/")
```

## Development

```bash
uv pip install -e ".[dev]"   # or: pip install -e ".[dev]"
pytest
```
